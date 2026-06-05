//! Sample-axis (column-wise) normalization for Stage 2.
//!
//! Six methods (`None`, `Sum`, `Median`, `Metadata`, `Quantile`, `PQN`)
//! operate on the `MetabolomicsTable.intensity_raw` matrix and produce a
//! working `intensity` matrix that the DAM runner stores back. `None` is
//! the default; the other five preserve magnitude by scaling to the median
//! per-sample factor so the Welch path's downstream arcsinh stays in a
//! useful range.

pub mod median;
pub mod metadata;
pub mod pqn;
pub mod quantile;
pub mod sum;
pub mod types;

pub use types::{NormalizationConfig, NormalizationError, NormalizationMethod, PqnReference};

use ndarray::Array2;
use tracing::info;

use crate::data::GroupMapping;

/// Apply `config.method` to `raw` and return the working matrix. Thin logging
/// wrapper over [`run`]: emits exactly one INFO `tracing` event per non-`None`
/// invocation in a greppable single-line format (`None` emits nothing). Pure —
/// no mutation of inputs. See [`run`] for the per-mode `sample_cols` contract,
/// and [`validate`] for the non-logging sibling used by the Stage 2 preflight.
pub fn apply(
    config: &NormalizationConfig,
    raw: &Array2<f64>,
    mapping: &GroupMapping,
    sample_cols: &[String],
) -> Result<Array2<f64>, NormalizationError> {
    run(config, raw, mapping, sample_cols, true)
}

/// Validate that `config.method` applies cleanly to `raw` WITHOUT emitting the
/// `normalize:` INFO summary line and WITHOUT returning the matrix: `Ok(())` on
/// success, or the IDENTICAL `Err` [`apply`] would return for the same inputs.
/// Used by the Stage 2 `start_dam` preflight so a dual-mode run logs each mode's
/// normalization once — from the real DAM worker — instead of twice. Lower-level
/// per-method diagnostics (e.g. the Metadata missing-value WARN) still fire;
/// only the summary INFO line is suppressed.
pub fn validate(
    config: &NormalizationConfig,
    raw: &Array2<f64>,
    mapping: &GroupMapping,
    sample_cols: &[String],
) -> Result<(), NormalizationError> {
    run(config, raw, mapping, sample_cols, false).map(|_| ())
}

/// Dispatch `config.method` over `raw`, emitting the `normalize:` INFO summary
/// line only when `log_stats`. The single shared body behind [`apply`] (logs)
/// and [`validate`] (silent) — one dispatch, one per-method implementation.
/// Pure — no mutation of inputs.
///
/// `sample_cols` is the per-mode sample column names for `raw` (i.e.
/// `IonModeTable.sample_cols`). In dual-mode the `mapping` is built from the
/// UNION of all modes' sample names, so per-sample lookups inside
/// `apply_metadata` and `apply_pqn` MUST use these per-mode names (NOT the
/// union-indexed positional accessor on `mapping`) to avoid cross-mode
/// index confusion (Findings #1, #2, #3 in the 2026-05-25 audit).
fn run(
    config: &NormalizationConfig,
    raw: &Array2<f64>,
    mapping: &GroupMapping,
    sample_cols: &[String],
    log_stats: bool,
) -> Result<Array2<f64>, NormalizationError> {
    if raw.nrows() == 0 || raw.ncols() == 0 {
        return Err(NormalizationError::EmptyMatrix);
    }

    let n_features = raw.nrows();
    let n_samples = raw.ncols();

    match &config.method {
        NormalizationMethod::None => Ok(raw.clone()),
        NormalizationMethod::Sum => {
            let (out, scale, nan_in, nan_out) = sum::apply_sum(raw, mapping, sample_cols)?;
            log_norm_stats(
                log_stats,
                "Sum",
                n_samples,
                n_features,
                "",
                "",
                &NormStats {
                    scale,
                    nan_in,
                    nan_out,
                },
            );
            Ok(out)
        }
        NormalizationMethod::Median => {
            let (out, scale, nan_in, nan_out) = median::apply_median(raw, mapping, sample_cols)?;
            log_norm_stats(
                log_stats,
                "Median",
                n_samples,
                n_features,
                "",
                "",
                &NormStats {
                    scale,
                    nan_in,
                    nan_out,
                },
            );
            Ok(out)
        }
        NormalizationMethod::Metadata { column } => {
            let (out, scale, nan_in, nan_out) =
                metadata::apply_metadata(raw, mapping, sample_cols, column)?;
            log_norm_stats(
                log_stats,
                "Metadata",
                n_samples,
                n_features,
                &format!("column={column} "),
                "",
                &NormStats {
                    scale,
                    nan_in,
                    nan_out,
                },
            );
            Ok(out)
        }
        NormalizationMethod::Quantile => {
            let (out, scale, nan_in, nan_out) = quantile::apply_quantile(raw, mapping)?;
            log_norm_stats(
                log_stats,
                "Quantile",
                n_samples,
                n_features,
                "",
                "",
                &NormStats {
                    scale,
                    nan_in,
                    nan_out,
                },
            );
            Ok(out)
        }
        NormalizationMethod::Pqn { reference } => {
            let (out, scale, nan_in, nan_out, ref_used) =
                pqn::apply_pqn(raw, mapping, sample_cols, reference)?;
            let ref_label = match reference {
                PqnReference::AllSamples => "AllSamples".to_string(),
                PqnReference::Group(name) => format!("Group({name})"),
            };
            log_norm_stats(
                log_stats,
                "PQN",
                n_samples,
                n_features,
                &format!("reference={ref_label} "),
                &format!("reference_features_used={ref_used} "),
                &NormStats {
                    scale,
                    nan_in,
                    nan_out,
                },
            );
            Ok(out)
        }
    }
}

/// Per-run normalization stats fed to [`log_norm_stats`].
pub(crate) struct NormStats {
    pub scale: f64,
    pub nan_in: usize,
    pub nan_out: usize,
}

/// Emit the single greppable `normalize: …` INFO line shared by every method.
/// A no-op when `enabled` is `false` (the `validate` / Stage 2 preflight path).
/// `pre_samples` is inserted between `method=` and `samples=` (Metadata's
/// `column=…`, PQN's `reference=…`); `pre_scale` between `features=` and
/// `scaling_to_median_factor=` (PQN's `reference_features_used=…`). Both empty
/// for Sum / Median / Quantile. The two slots reproduce each method's prior
/// field order byte-for-byte.
fn log_norm_stats(
    enabled: bool,
    method: &str,
    n_samples: usize,
    n_features: usize,
    pre_samples: &str,
    pre_scale: &str,
    s: &NormStats,
) {
    if !enabled {
        return;
    }
    info!(
        "normalize: method={method} {pre_samples}samples={n_samples} features={n_features} {pre_scale}scaling_to_median_factor={scale:.6e} nan_cells_in={nan_in} nan_cells_out={nan_out}",
        scale = s.scale,
        nan_in = s.nan_in,
        nan_out = s.nan_out,
    );
}

/// The single factor→write-back loop shared by all four divide-and-rescale
/// methods. `factor[j] == None` NaN-outs the entire column (Metadata's drop
/// case: `nan_in` counts only cells already NaN in `src`, `nan_out` counts
/// every cell); `Some(f)` rescales finite cells by `src / f * median_factor`
/// and passes NaN through. `median_factor` is the NaN-aware median over the
/// `Some` factors. Returns `(out, median_factor, nan_in, nan_out)`.
pub(crate) fn apply_factors_and_count(
    src: &Array2<f64>,
    factor: &[Option<f64>],
) -> (Array2<f64>, f64, usize, usize) {
    let (n_features, n_samples) = (src.nrows(), src.ncols());
    let finite: Vec<f64> = factor.iter().filter_map(|f| *f).collect();
    let median_factor = sum::median_of(&finite);
    let mut out = Array2::<f64>::zeros((n_features, n_samples));
    let (mut nan_in, mut nan_out) = (0usize, 0usize);
    for j in 0..n_samples {
        match factor[j] {
            None => {
                for i in 0..n_features {
                    if src[[i, j]].is_nan() {
                        nan_in += 1;
                    }
                    out[[i, j]] = f64::NAN;
                    nan_out += 1;
                }
            }
            Some(f) => {
                for i in 0..n_features {
                    let v = src[[i, j]];
                    if v.is_nan() {
                        nan_in += 1;
                        nan_out += 1;
                        out[[i, j]] = f64::NAN;
                    } else {
                        out[[i, j]] = v / f * median_factor;
                    }
                }
            }
        }
    }
    (out, median_factor, nan_in, nan_out)
}

/// Driver for the per-sample-scalar methods (Sum, Median): compute each
/// sample's factor via `factor_fn(column)`, raise the same `NanFactor` /
/// `ZeroFactor` errors (per-mode `sample_label` attribution) on a NaN / zero
/// factor, then delegate the write-back to [`apply_factors_and_count`].
/// `_mapping` is kept for signature uniformity with the metadata/PQN methods.
pub(crate) fn apply_per_sample_factor(
    raw: &Array2<f64>,
    _mapping: &GroupMapping,
    sample_cols: &[String],
    method: &'static str,
    factor_fn: impl Fn(&[f64]) -> f64,
) -> Result<(Array2<f64>, f64, usize, usize), NormalizationError> {
    let (n_features, n_samples) = (raw.nrows(), raw.ncols());
    // See the apply_sum doc comment for the debug_assert vs assert rationale.
    debug_assert_eq!(
        sample_cols.len(),
        n_samples,
        "sample_cols must be aligned to raw column axis (got {} vs {})",
        sample_cols.len(),
        n_samples
    );
    let mut factor: Vec<Option<f64>> = Vec::with_capacity(n_samples);
    for j in 0..n_samples {
        let col: Vec<f64> = (0..n_features).map(|i| raw[[i, j]]).collect();
        let f = factor_fn(&col);
        if f.is_nan() {
            return Err(NormalizationError::NanFactor {
                sample: sum::sample_label(sample_cols, j),
                method,
            });
        }
        if f == 0.0 {
            return Err(NormalizationError::ZeroFactor {
                sample: sum::sample_label(sample_cols, j),
                method,
            });
        }
        factor.push(Some(f));
    }
    Ok(apply_factors_and_count(raw, &factor))
}

/// Shared `#[cfg(test)]` fixture builders for the normalization unit tests,
/// consolidated here so `sum`/`median`/`metadata`/`pqn`/`quantile`/`mod` no
/// longer each carry their own copy. Bit-neutral: every helper produces the
/// same `GroupMapping` the per-module copies did.
#[cfg(test)]
pub(crate) mod test_support {
    use crate::data::{GroupMapping, load_group_mapping};
    use std::io::Write;

    /// Build a `GroupMapping` from parallel `samples` / `groups` slices via a
    /// temp `sample,group` CSV. Single-group fixtures pass `&["G1"; n]`.
    pub(crate) fn mapping_for(samples: &[&str], groups: &[&str]) -> GroupMapping {
        let mut content = String::from("sample,group\n");
        for (s, g) in samples.iter().zip(groups.iter()) {
            content.push_str(&format!("{s},{g}\n"));
        }
        let mut f = tempfile::Builder::new().suffix(".csv").tempfile().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        let cols: Vec<String> = samples.iter().map(|s| s.to_string()).collect();
        load_group_mapping(f.path(), &cols).expect("mapping")
    }

    /// Build a `GroupMapping` with one extra column `col` (group fixed to `G1`),
    /// used by the metadata-normalization tests.
    pub(crate) fn mapping_with_column(
        samples: &[&str],
        col: &str,
        values: &[&str],
    ) -> GroupMapping {
        let mut content = format!("sample,group,{col}\n");
        for (s, v) in samples.iter().zip(values.iter()) {
            content.push_str(&format!("{s},G1,{v}\n"));
        }
        let mut f = tempfile::Builder::new().suffix(".csv").tempfile().unwrap();
        f.write_all(content.as_bytes()).unwrap();
        let cols: Vec<String> = samples.iter().map(|s| s.to_string()).collect();
        load_group_mapping(f.path(), &cols).expect("mapping")
    }

    pub(crate) fn cols(names: &[&str]) -> Vec<String> {
        names.iter().map(|s| s.to_string()).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::test_support::{cols, mapping_for};
    use super::*;
    use ndarray::array;

    #[test]
    fn none_returns_clone() {
        let raw = array![[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]];
        let mapping = mapping_for(&["A", "B"], &["G1", "G2"]);
        let sc = cols(&["A", "B"]);
        let out = apply(&NormalizationConfig::default(), &raw, &mapping, &sc).unwrap();
        assert_eq!(out, raw);
    }

    #[test]
    fn empty_matrix_errors() {
        let mapping = mapping_for(&["A"], &["G1"]);
        let sc = cols(&["A"]);
        let zero_rows = Array2::<f64>::zeros((0, 1));
        let zero_cols = Array2::<f64>::zeros((3, 0));
        assert!(matches!(
            apply(&NormalizationConfig::default(), &zero_rows, &mapping, &sc),
            Err(NormalizationError::EmptyMatrix)
        ));
        let sc_empty: Vec<String> = vec![];
        assert!(matches!(
            apply(
                &NormalizationConfig::default(),
                &zero_cols,
                &mapping,
                &sc_empty
            ),
            Err(NormalizationError::EmptyMatrix)
        ));
    }

    #[test]
    fn sum_preserves_magnitude() {
        // Sample A total = 1000; sample B total = 2000; sample C total = 1500
        // median_factor = 1500; each column post-norm sum == 1500
        let raw = array![
            [100.0, 200.0, 150.0],
            [400.0, 800.0, 600.0],
            [500.0, 1000.0, 750.0]
        ];
        let mapping = mapping_for(&["A", "B", "C"], &["G1", "G1", "G1"]);
        let sc = cols(&["A", "B", "C"]);
        let out = apply(
            &NormalizationConfig {
                method: NormalizationMethod::Sum,
            },
            &raw,
            &mapping,
            &sc,
        )
        .unwrap();
        for j in 0..3 {
            let col_sum: f64 = (0..3).map(|i| out[[i, j]]).sum();
            assert!((col_sum - 1500.0).abs() < 1e-9, "col {j} sum = {col_sum}");
        }
    }

    #[test]
    fn shape_preserved_for_every_method() {
        let raw = array![[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
        let mapping = mapping_for(&["A", "B", "C"], &["G1", "G2", "G1"]);
        let sc = cols(&["A", "B", "C"]);
        let methods = [
            NormalizationMethod::None,
            NormalizationMethod::Sum,
            NormalizationMethod::Median,
            NormalizationMethod::Quantile,
            NormalizationMethod::Pqn {
                reference: PqnReference::AllSamples,
            },
        ];
        for m in methods {
            let out = apply(
                &NormalizationConfig { method: m.clone() },
                &raw,
                &mapping,
                &sc,
            )
            .unwrap();
            assert_eq!(out.shape(), raw.shape(), "method {m:?} must preserve shape");
        }
    }

    #[test]
    fn validate_ok_when_apply_ok() {
        // A non-`None` method `apply` accepts; `validate` must mirror it —
        // return `Ok(())` and discard the matrix.
        let raw = array![
            [100.0, 200.0, 150.0],
            [400.0, 800.0, 600.0],
            [500.0, 1000.0, 750.0]
        ];
        let mapping = mapping_for(&["A", "B", "C"], &["G1", "G1", "G1"]);
        let sc = cols(&["A", "B", "C"]);
        let config = NormalizationConfig {
            method: NormalizationMethod::Sum,
        };
        assert!(apply(&config, &raw, &mapping, &sc).is_ok());
        assert!(validate(&config, &raw, &mapping, &sc).is_ok());
    }

    #[test]
    fn validate_returns_same_err_as_apply() {
        // Degenerate-PQN fixture (mirrors `pqn::tests`): sample D's quotient
        // median is 0 → `PqnDegenerateSamples`. `validate` must surface the
        // identical error `apply` would, proving the Setup-gate check still
        // fires when logging is suppressed.
        let raw = array![
            [10.0, 10.0, 10.0, 0.0],
            [20.0, 20.0, 20.0, 0.0],
            [30.0, 30.0, 30.0, 0.0],
            [40.0, 40.0, 40.0, 100.0],
        ];
        let mapping = mapping_for(&["A", "B", "C", "D"], &["G1", "G1", "G1", "G1"]);
        let sc = cols(&["A", "B", "C", "D"]);
        let config = NormalizationConfig {
            method: NormalizationMethod::Pqn {
                reference: PqnReference::AllSamples,
            },
        };
        assert!(matches!(
            apply(&config, &raw, &mapping, &sc),
            Err(NormalizationError::PqnDegenerateSamples { .. })
        ));
        assert!(matches!(
            validate(&config, &raw, &mapping, &sc),
            Err(NormalizationError::PqnDegenerateSamples { .. })
        ));
    }

    #[test]
    fn validate_none_returns_ok() {
        // Parity with `apply`'s `None` passthrough: `validate(None)` is `Ok(())`.
        let raw = array![[1.0, 2.0], [3.0, 4.0]];
        let mapping = mapping_for(&["A", "B"], &["G1", "G2"]);
        let sc = cols(&["A", "B"]);
        assert!(validate(&NormalizationConfig::default(), &raw, &mapping, &sc).is_ok());
    }
}
