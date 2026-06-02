//! Sum normalization. Per-sample factor = NaN-aware sum; magnitude is
//! preserved by multiplying the result by the median sample sum.

use ndarray::Array2;

use crate::data::GroupMapping;
use crate::normalize::types::NormalizationError;

/// Returns `(out, median_factor, nan_cells_in, nan_cells_out)`.
///
/// `sample_cols` is the per-mode sample column names (length == raw.ncols()).
/// In dual-mode `mapping.sample_names` is the UNION across modes; using
/// `mapping.sample_name(j)` for per-mode `j` would resolve a NEG column index
/// to a POS sample name in error messages (PR-K, the leftover from the
/// PR-H/J review). `sample_label` now keys on `sample_cols` so error
/// attribution is per-mode-correct. `_mapping` is kept on the signature for
/// API uniformity with apply_metadata/apply_pqn even though Sum doesn't need
/// it for the math.
pub fn apply_sum(
    raw: &Array2<f64>,
    mapping: &GroupMapping,
    sample_cols: &[String],
) -> Result<(Array2<f64>, f64, usize, usize), NormalizationError> {
    // The per-sample factor is the NaN-aware running sum (all-NaN column →
    // NaN → NanFactor; finite sum of 0.0 → ZeroFactor). The accumulation
    // order is i = 0..n_features, identical to the pre-refactor loop, so the
    // running sum is floating-point bit-equal. The driver owns the per-column
    // NaN/zero dispatch and the write-back.
    crate::normalize::apply_per_sample_factor(raw, mapping, sample_cols, "Sum", |col| {
        let mut sum = 0.0;
        let mut n_finite = 0usize;
        for &v in col {
            if !v.is_nan() {
                sum += v;
                n_finite += 1;
            }
        }
        if n_finite == 0 { f64::NAN } else { sum }
    })
}

/// Resolve sample `j`'s display label from the per-mode `sample_cols`
/// slice (real name if `j < sample_cols.len()`, else `(sample j)`). Keyed
/// on per-mode names so dual-mode error messages attribute the failure to
/// the correct ion mode's sample (PR-K, closing the leftover from the
/// PR-H/J review where apply_sum/apply_median used the union-indexed
/// `mapping.sample_name(j)` and could name a POS sample on a NEG-mode
/// failure).
pub(crate) fn sample_label(sample_cols: &[String], j: usize) -> String {
    sample_cols
        .get(j)
        .cloned()
        .unwrap_or_else(|| format!("(sample {j})"))
}

/// NaN-aware median of a 1D slice. Returns NaN for an all-NaN / empty input.
pub(crate) fn median_of(values: &[f64]) -> f64 {
    let mut v: Vec<f64> = values.iter().copied().filter(|x| !x.is_nan()).collect();
    if v.is_empty() {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).expect("non-NaN compare"));
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        (v[n / 2 - 1] + v[n / 2]) / 2.0
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::test_support::{cols, mapping_for};
    use crate::normalize::types::NormalizationError;
    use ndarray::array;

    #[test]
    fn sum_hand_computed_3x3() {
        // Sums: A=6, B=15, C=24; median_factor = 15
        // out[i,j] = raw[i,j] / sum_j * 15
        let raw = array![[1.0, 4.0, 7.0], [2.0, 5.0, 8.0], [3.0, 6.0, 9.0]];
        let mapping = mapping_for(&["A", "B", "C"], &["G1"; 3]);
        let sc = cols(&["A", "B", "C"]);
        let (out, scale, nan_in, nan_out) = apply_sum(&raw, &mapping, &sc).unwrap();
        assert!((scale - 15.0).abs() < 1e-9);
        assert_eq!(nan_in, 0);
        assert_eq!(nan_out, 0);
        // Column A: 1*15/6, 2*15/6, 3*15/6 = 2.5, 5.0, 7.5
        assert!((out[[0, 0]] - 2.5).abs() < 1e-9);
        assert!((out[[1, 0]] - 5.0).abs() < 1e-9);
        assert!((out[[2, 0]] - 7.5).abs() < 1e-9);
        // Column B: identity (sum already == median)
        assert!((out[[0, 1]] - 4.0).abs() < 1e-9);
        assert!((out[[1, 1]] - 5.0).abs() < 1e-9);
        assert!((out[[2, 1]] - 6.0).abs() < 1e-9);
    }

    #[test]
    fn sum_nan_passthrough() {
        let raw = array![[f64::NAN, 4.0], [2.0, f64::NAN], [3.0, 6.0]];
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let sc = cols(&["A", "B"]);
        let (out, _, nan_in, nan_out) = apply_sum(&raw, &mapping, &sc).unwrap();
        assert_eq!(nan_in, 2);
        assert_eq!(nan_out, 2);
        assert!(out[[0, 0]].is_nan());
        assert!(out[[1, 1]].is_nan());
        assert!(!out[[2, 0]].is_nan());
    }

    #[test]
    fn sum_all_nan_sample_errors_with_nan_factor() {
        let raw = array![[1.0, f64::NAN], [2.0, f64::NAN]];
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let sc = cols(&["A", "B"]);
        let err = apply_sum(&raw, &mapping, &sc).expect_err("all-NaN sample must error");
        match err {
            NormalizationError::NanFactor { sample, method } => {
                assert_eq!(sample, "B");
                assert_eq!(method, "Sum");
            }
            other => panic!("expected NanFactor, got {other:?}"),
        }
    }

    #[test]
    fn sum_zero_sum_sample_errors_with_zero_factor() {
        let raw = array![[0.0, 4.0], [0.0, 5.0]];
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let sc = cols(&["A", "B"]);
        let err = apply_sum(&raw, &mapping, &sc).expect_err("zero-sum sample must error");
        match err {
            NormalizationError::ZeroFactor { sample, method } => {
                assert_eq!(sample, "A");
                assert_eq!(method, "Sum");
            }
            other => panic!("expected ZeroFactor, got {other:?}"),
        }
    }

    /// PR-K regression: in dual-mode, an apply_sum failure on the NEG
    /// table must name the offending NEG sample (NOT the POS sample at
    /// the same per-mode column index). Pre-PR-K `sample_label(mapping, j)`
    /// resolved per-mode `j` against the union-aligned mapping, so NEG_S01's
    /// failure was labelled `POS_S01` — a debugging trap. The fix keys
    /// `sample_label` on the per-mode `sample_cols` slice.
    #[test]
    fn sum_dual_mode_neg_failure_names_neg_sample_not_pos() {
        // Union mapping has 4 samples: POS_S01, POS_S02, NEG_S01, NEG_S02.
        let mapping = mapping_for(&["POS_S01", "POS_S02", "NEG_S01", "NEG_S02"], &["G1"; 4]);
        // NEG table: 1 feature × 2 samples. NEG_S01 has sum=0.
        let neg_raw = array![[0.0, 5.0]];
        let neg_cols = cols(&["NEG_S01", "NEG_S02"]);
        let err =
            apply_sum(&neg_raw, &mapping, &neg_cols).expect_err("NEG_S01 zero sum must error");
        match err {
            NormalizationError::ZeroFactor { sample, method } => {
                assert_eq!(
                    sample, "NEG_S01",
                    "must name the NEG sample, not the POS sample at the same per-mode index"
                );
                assert_eq!(method, "Sum");
            }
            other => panic!("expected ZeroFactor, got {other:?}"),
        }
    }
}
