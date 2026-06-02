//! Metadata-based normalization. Per-sample factor = user-selected numeric
//! metadata column (e.g. dry weight).
//!
//! Samples whose metadata value is missing (`None`) are **dropped** —
//! their entire column in the output matrix is set to NaN so DAM's
//! NaN-aware machinery naturally excludes them from per-feature stats.
//! Samples with a non-positive value still error (`MetadataValueNonPositive`)
//! because zero/negative metadata is a data-entry problem, not absence.
//! The caller (`start_dam`) is responsible for a *preflight* that ensures
//! dropping these samples doesn't leave either DAM group below 2; see
//! `NormalizationError::InsufficientSamplesAfterDrop`.

use ndarray::Array2;
use tracing::warn;

use crate::data::GroupMapping;
use crate::normalize::types::NormalizationError;

/// Returns `(out, median_factor, nan_cells_in, nan_cells_out)`. Samples whose
/// metadata value is `None` are NaN-marked in `out`; their NaN cells are
/// counted in `nan_cells_out` but NOT in `nan_cells_in` (they weren't NaN
/// in `raw`).
///
/// `sample_cols` is the per-mode list of sample column names from the
/// `IonModeTable` being normalized (length == `raw.ncols()`). In dual-mode
/// `mapping.sample_names` is the UNION of all modes' samples, so we MUST
/// look up metadata values by name via `mapping.metadata_value_of`, not by
/// positional index against the union slice (Finding #1 in the 2026-05-25
/// audit — release builds previously applied POS metadata values to NEG
/// sample columns because the union-indexed slice was indexed with per-mode
/// `j`; the `debug_assert_eq!(col_values.len(), n_samples)` panicked in
/// debug but was elided in release).
pub fn apply_metadata(
    raw: &Array2<f64>,
    mapping: &GroupMapping,
    sample_cols: &[String],
    column: &str,
) -> Result<(Array2<f64>, f64, usize, usize), NormalizationError> {
    // Existence check: distinguish "column missing entirely" from "column
    // present but this sample has no row in the mapping".
    if mapping.metadata_values(column).is_none() {
        return Err(NormalizationError::MetadataColumnMissing {
            column: column.to_string(),
        });
    }
    let n_samples = raw.ncols();
    // See apply_sum for the debug_assert vs assert rationale (PR-M).
    debug_assert_eq!(
        sample_cols.len(),
        n_samples,
        "sample_cols must be aligned to raw column axis (got {} vs {})",
        sample_cols.len(),
        n_samples
    );
    // Per-sample factor: `Some(v)` for samples with a value (must be > 0),
    // `None` for samples we will drop (NaN-out the entire column).
    let mut factor: Vec<Option<f64>> = Vec::with_capacity(n_samples);
    let mut dropped: Vec<String> = Vec::new();
    for name in sample_cols.iter() {
        // PR-O (Finding #2 from PR-H/J review): metadata_value_of returns
        // None when (a) the column is missing OR (b) the sample name is
        // absent from mapping.sample_names. The upstream column-existence
        // check above already excluded (a), so a None here can only mean
        // (b). By construction this is unreachable: mapping.sample_names
        // is built as the union of every ion mode's `sample_cols`
        // (src/ui/stage1_input.rs:638 `union_sample_cols`), so every name
        // in any per-mode sample_cols IS in the mapping. Pre-PR-O this
        // surfaced as `MetadataColumnMissing` — a wrong-variant message
        // that would mislead a future maintainer if the invariant ever
        // broke. By construction this is unreachable for a consistent
        // (mapping, tables) pair, but a stale pair after a Stage 1 file
        // re-pick can produce it — so rather than panic, return
        // `SampleNotInMapping` and let it surface as a recoverable Stage 2
        // banner error (`convert-defensive-panics-to-errors`).
        let Some(opt) = mapping.metadata_value_of(name, column) else {
            return Err(NormalizationError::SampleNotInMapping {
                sample: name.clone(),
                column: column.to_string(),
            });
        };
        match opt {
            Some(v) if v > 0.0 => factor.push(Some(v)),
            Some(v) => {
                return Err(NormalizationError::MetadataValueNonPositive {
                    sample: name.clone(),
                    column: column.to_string(),
                    value: v,
                });
            }
            None => {
                factor.push(None);
                dropped.push(name.clone());
            }
        }
    }
    if !dropped.is_empty() {
        warn!(
            column,
            dropped_count = dropped.len(),
            samples = ?dropped,
            "metadata normalization: samples without a value in this column will be dropped (NaN-marked)"
        );
    }
    // Write-back via the one shared loop: `None` factors NaN-out their column
    // (the metadata drop case — `nan_in` only where `raw` was already NaN,
    // `nan_out` always), `Some(f)` rescales by `raw / f * median_factor`. The
    // `median_factor` is the NaN-aware median over the `Some` factors, exactly
    // as the prior `median_of(&finite_factors)`.
    Ok(crate::normalize::apply_factors_and_count(raw, &factor))
}

/// Return the sample names that the Metadata-column normalization would
/// drop given the current mapping. Used by Stage 2 UI to render a yellow
/// pre-DAM warning and by `start_dam` for the per-group preflight.
/// Returns an empty `Vec` when the column is missing or fully populated.
///
/// Unlike `apply_metadata`, this helper iterates the union-aligned
/// `col_values` directly (the UI hint covers every sample in the metadata
/// CSV regardless of which ion mode each one belongs to), so the per-mode
/// `sample_cols` API isn't applicable here — we resolve names via
/// `mapping.sample_name(j)` against the union slice on purpose.
pub fn dropped_samples(mapping: &GroupMapping, column: &str) -> Vec<String> {
    let Some(col_values) = mapping.metadata_values(column) else {
        return Vec::new();
    };
    col_values
        .iter()
        .enumerate()
        .filter_map(|(j, opt)| {
            if opt.is_none() {
                Some(
                    mapping
                        .sample_name(j)
                        .map(|s| s.to_string())
                        .unwrap_or_else(|| format!("(sample {j})")),
                )
            } else {
                None
            }
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::test_support::{cols, mapping_with_column};
    use crate::normalize::types::NormalizationError;
    use ndarray::array;

    #[test]
    fn metadata_sample_missing_from_mapping_errors() {
        // A stale (mapping, per-mode sample_cols) pair: `Z` appears in the
        // columns passed to apply_metadata but is absent from the mapping.
        // Pre-change this panicked; now it returns SampleNotInMapping so the
        // skew surfaces as a recoverable Stage 2 banner error.
        let raw = array![[1.0, 2.0]];
        let mapping = mapping_with_column(&["A", "B"], "dw", &["10", "20"]);
        let sc = cols(&["A", "Z"]);
        let err = apply_metadata(&raw, &mapping, &sc, "dw").expect_err("Z absent must error");
        match err {
            NormalizationError::SampleNotInMapping { sample, column } => {
                assert_eq!(sample, "Z");
                assert_eq!(column, "dw");
            }
            other => panic!("expected SampleNotInMapping, got {other:?}"),
        }
    }

    #[test]
    fn metadata_hand_computed() {
        // Column [10, 20, 40]; median = 20; out[i, j] = raw[i, j] / col[j] * 20
        let raw = array![[5.0, 10.0, 20.0], [50.0, 100.0, 200.0]];
        let mapping = mapping_with_column(&["A", "B", "C"], "dw", &["10", "20", "40"]);
        let sc = cols(&["A", "B", "C"]);
        let (out, scale, _, _) = apply_metadata(&raw, &mapping, &sc, "dw").unwrap();
        assert!((scale - 20.0).abs() < 1e-9);
        // A: 5/10*20 = 10
        assert!((out[[0, 0]] - 10.0).abs() < 1e-9);
        // B: 10/20*20 = 10 (identity since column == median)
        assert!((out[[0, 1]] - 10.0).abs() < 1e-9);
        // C: 20/40*20 = 10
        assert!((out[[0, 2]] - 10.0).abs() < 1e-9);
    }

    #[test]
    fn metadata_missing_value_drops_sample() {
        // Sample B has an empty value -> its column is NaN-marked in the
        // output. Sample A is normalized as usual.
        let raw = array![[1.0, 2.0], [3.0, 4.0]];
        let mapping = mapping_with_column(&["A", "B"], "dw", &["10", ""]);
        let sc = cols(&["A", "B"]);
        let (out, scale, nan_in, nan_out) =
            apply_metadata(&raw, &mapping, &sc, "dw").expect("must NOT error on missing");
        // Only one valid factor (10), so the median is 10. Sample A: identity.
        assert!((scale - 10.0).abs() < 1e-9);
        assert!((out[[0, 0]] - 1.0).abs() < 1e-9);
        assert!((out[[1, 0]] - 3.0).abs() < 1e-9);
        // Sample B's entire column is NaN.
        assert!(out[[0, 1]].is_nan());
        assert!(out[[1, 1]].is_nan());
        // NaN counts: 0 in (raw had no NaN), 2 out (sample B's two cells).
        assert_eq!(nan_in, 0);
        assert_eq!(nan_out, 2);
    }

    #[test]
    fn metadata_dropped_samples_helper_lists_missing() {
        let mapping = mapping_with_column(&["A", "B", "C"], "dw", &["10", "", "20"]);
        let dropped = dropped_samples(&mapping, "dw");
        assert_eq!(dropped, vec!["B".to_string()]);
    }

    #[test]
    fn metadata_dropped_samples_helper_empty_for_unknown_column() {
        let mapping = mapping_with_column(&["A"], "dw", &["10"]);
        assert!(dropped_samples(&mapping, "nope").is_empty());
    }

    #[test]
    fn metadata_negative_value_errors() {
        let raw = array![[1.0, 2.0]];
        let mapping = mapping_with_column(&["A", "B"], "dw", &["10", "-5"]);
        let sc = cols(&["A", "B"]);
        let err = apply_metadata(&raw, &mapping, &sc, "dw").expect_err("negative must error");
        match err {
            NormalizationError::MetadataValueNonPositive {
                sample,
                column,
                value,
            } => {
                assert_eq!(sample, "B");
                assert_eq!(column, "dw");
                assert!((value - -5.0).abs() < 1e-9);
            }
            other => panic!("expected MetadataValueNonPositive, got {other:?}"),
        }
    }

    #[test]
    fn metadata_unknown_column_errors() {
        let raw = array![[1.0, 2.0]];
        let mapping = mapping_with_column(&["A", "B"], "dw", &["10", "20"]);
        let sc = cols(&["A", "B"]);
        let err =
            apply_metadata(&raw, &mapping, &sc, "nope").expect_err("unknown column must error");
        assert!(matches!(
            err,
            NormalizationError::MetadataColumnMissing { ref column } if column == "nope"
        ));
    }

    /// Regression test for Finding #1: dual-mode metadata normalization
    /// must use each ion mode's own dry_weight values (looked up by sample
    /// NAME), NOT the union-positional values which mapped POS dw to NEG
    /// columns. Asymmetric ratios make the divergence numerically visible:
    /// POS uses ratio 1:1.2, NEG uses ratio 1:2.
    #[test]
    fn metadata_dual_mode_neg_table_uses_neg_metadata_values() {
        // Union mapping: 4 samples with distinct dry_weight values.
        let mapping = mapping_with_column(
            &["POS_S01", "POS_S02", "NEG_S01", "NEG_S02"],
            "dw",
            &["10", "12", "50", "100"], // POS: 10, 12 (ratio 1.2); NEG: 50, 100 (ratio 2)
        );

        // POS table: 1 feature, 2 samples.
        let pos_raw = array![[100.0, 100.0]];
        let pos_cols = cols(&["POS_S01", "POS_S02"]);
        let (pos_out, pos_scale, _, _) =
            apply_metadata(&pos_raw, &mapping, &pos_cols, "dw").expect("POS ok");
        // Per-mode factors [10, 12]; median = 11.
        assert!((pos_scale - 11.0).abs() < 1e-9, "POS scale = {pos_scale}");
        assert!((pos_out[[0, 0]] - 110.0).abs() < 1e-9); // 100/10*11
        assert!((pos_out[[0, 1]] - 100.0 / 12.0 * 11.0).abs() < 1e-6);

        // NEG table: 1 feature, 2 samples — identical raw values to POS but
        // with VERY different metadata, so the output magnitudes must differ.
        let neg_raw = array![[100.0, 100.0]];
        let neg_cols = cols(&["NEG_S01", "NEG_S02"]);
        let (neg_out, neg_scale, _, _) =
            apply_metadata(&neg_raw, &mapping, &neg_cols, "dw").expect("NEG ok");
        // Per-mode factors [50, 100]; median = 75. Post-fix:
        //   NEG_S01 = 100 / 50 * 75 = 150
        //   NEG_S02 = 100 / 100 * 75 = 75
        // Pre-fix (buggy): factor array was union-length [10,12,50,100],
        //   first 2 entries used for NEG → NEG_S01 = 100/10*(median 31) =
        //   310, NEG_S02 = 100/12*31 ≈ 258. The post-fix 150 and 75 cannot
        //   be reached by any single union-positional combination, so this
        //   test bracket genuinely guards against the regression.
        assert!((neg_scale - 75.0).abs() < 1e-9, "NEG scale = {neg_scale}");
        assert!(
            (neg_out[[0, 0]] - 150.0).abs() < 1e-9,
            "NEG_S01 expected 150 (uses NEG dw=50), got {}",
            neg_out[[0, 0]]
        );
        assert!(
            (neg_out[[0, 1]] - 75.0).abs() < 1e-9,
            "NEG_S02 expected 75 (uses NEG dw=100), got {}",
            neg_out[[0, 1]]
        );
    }
}
