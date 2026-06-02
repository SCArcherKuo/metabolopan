//! Median normalization. Per-sample factor = NaN-aware median; magnitude is
//! preserved by multiplying by the median of per-sample medians.

use ndarray::Array2;

use crate::data::GroupMapping;
use crate::normalize::sum::median_of;
use crate::normalize::types::NormalizationError;

/// Returns `(out, median_factor, nan_cells_in, nan_cells_out)`.
///
/// `sample_cols` is the per-mode sample column names; required (with
/// `_mapping` kept for API uniformity) so error messages report the correct
/// per-mode sample name in dual-mode runs. See `apply_sum` for the rationale
/// (PR-K, leftover from the PR-H/J review).
pub fn apply_median(
    raw: &Array2<f64>,
    mapping: &GroupMapping,
    sample_cols: &[String],
) -> Result<(Array2<f64>, f64, usize, usize), NormalizationError> {
    // Per-sample factor is the NaN-aware median; the driver owns the per-column
    // NaN/zero dispatch (NanFactor / ZeroFactor with "Median") and the
    // write-back. Output is bit-identical to the prior hand-rolled loop.
    crate::normalize::apply_per_sample_factor(raw, mapping, sample_cols, "Median", median_of)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::test_support::{cols, mapping_for};
    use crate::normalize::types::NormalizationError;
    use ndarray::array;

    #[test]
    fn median_hand_computed() {
        // Per-sample medians: A=2 (1,2,3), B=20 (10,20,30), C=200 (100,200,300)
        // median_factor = 20
        let raw = array![[1.0, 10.0, 100.0], [2.0, 20.0, 200.0], [3.0, 30.0, 300.0]];
        let mapping = mapping_for(&["A", "B", "C"], &["G1"; 3]);
        let sc = cols(&["A", "B", "C"]);
        let (out, scale, _, _) = apply_median(&raw, &mapping, &sc).unwrap();
        assert!((scale - 20.0).abs() < 1e-9);
        // Sample A: divide by 2, multiply by 20 -> 10, 20, 30
        assert!((out[[0, 0]] - 10.0).abs() < 1e-9);
        assert!((out[[1, 0]] - 20.0).abs() < 1e-9);
        assert!((out[[2, 0]] - 30.0).abs() < 1e-9);
        // Sample B already at median; identity
        assert!((out[[1, 1]] - 20.0).abs() < 1e-9);
    }

    #[test]
    fn median_nan_aware_per_sample() {
        // Sample A median over non-NaN: median(10, 30, 50) = 30
        let raw = array![
            [10.0, 2.0],
            [f64::NAN, 4.0],
            [30.0, 6.0],
            [f64::NAN, 8.0],
            [50.0, 10.0]
        ];
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let sc = cols(&["A", "B"]);
        let (out, _, nan_in, nan_out) = apply_median(&raw, &mapping, &sc).unwrap();
        assert_eq!(nan_in, 2);
        assert_eq!(nan_out, 2);
        // Sample A's factor = 30 (NaNs dropped); 30 / 30 == 1 then * median_factor
        assert!(out[[2, 0]] > 0.0);
        assert!(out[[1, 0]].is_nan());
    }

    #[test]
    fn median_zero_median_errors() {
        // Sample A has odd count with median value 0
        let raw = array![[0.0, 1.0], [0.0, 2.0], [1.0, 3.0]];
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let sc = cols(&["A", "B"]);
        let err = apply_median(&raw, &mapping, &sc).expect_err("zero median must error");
        assert!(matches!(
            err,
            NormalizationError::ZeroFactor {
                ref sample,
                method: "Median"
            } if sample == "A"
        ));
    }

    #[test]
    fn median_all_nan_sample_errors() {
        let raw = array![[f64::NAN, 1.0], [f64::NAN, 2.0]];
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let sc = cols(&["A", "B"]);
        let err = apply_median(&raw, &mapping, &sc).expect_err("all-NaN sample must error");
        assert!(matches!(
            err,
            NormalizationError::NanFactor {
                ref sample,
                method: "Median"
            } if sample == "A"
        ));
    }

    /// PR-K regression: dual-mode median failure on NEG must name the NEG
    /// sample. Mirrors the same fix that apply_sum needed; both functions
    /// previously used the union-aligned `sample_label(mapping, j)` which
    /// resolved per-mode `j` to a POS sample's name.
    #[test]
    fn median_dual_mode_neg_failure_names_neg_sample_not_pos() {
        let mapping = mapping_for(&["POS_S01", "POS_S02", "NEG_S01", "NEG_S02"], &["G1"; 4]);
        // NEG table: NEG_S01 column is all-NaN.
        let neg_raw = array![[f64::NAN, 1.0], [f64::NAN, 2.0], [f64::NAN, 3.0]];
        let neg_cols = cols(&["NEG_S01", "NEG_S02"]);
        let err =
            apply_median(&neg_raw, &mapping, &neg_cols).expect_err("NEG_S01 all-NaN must error");
        match err {
            NormalizationError::NanFactor { sample, method } => {
                assert_eq!(
                    sample, "NEG_S01",
                    "must name the NEG sample, not the POS sample at the same per-mode index"
                );
                assert_eq!(method, "Median");
            }
            other => panic!("expected NanFactor, got {other:?}"),
        }
    }
}
