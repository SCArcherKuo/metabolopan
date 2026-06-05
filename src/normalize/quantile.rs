//! Quantile normalization. Each run of tied values in a sample (sorted
//! positions `[k, k2)`) is assigned `mean(reference[k..k2])` — the mean of
//! the reference (pooled-quantile) values across every rank position the
//! tie spans. This is metabolopan's *literal* reading of Smyth's remark in
//! the Bioconductor support thread #1569 (2003) that tied items should get
//! "the average of the corresponding pooled quantiles"; it is NOT a
//! consensus Bolstad and Smyth reached — the thread's only worked example
//! is a 2-way tie, where every reading coincides, and it never adjudicated
//! `N ≥ 3`. The canonical tools — including Smyth's own
//! `limma::normalizeQuantiles` with `ties=TRUE` (the default: `rank()`
//! average-rank then `approx()` interpolation) and Bolstad's
//! `preprocessCore::normalize.quantiles` — instead return the reference
//! value at the tie's middle rank. The two coincide for `N == 2` ties or
//! locally-linear reference regions and diverge for `N ≥ 3` ties on a
//! non-linear reference (common at the bottom of LOD-imputed metabolomics
//! samples). metabolopan therefore differs from BOTH tools in that case, by
//! design — see USER_MANUAL.md for the worked example. Reference:
//! <https://support.bioconductor.org/p/1569/>.
//!
//! When samples have unequal non-NaN counts (e.g. heterogeneous
//! missingness across the input), we build the reference on a common
//! fractional-rank grid of size `K = max(n_j)` and linearly interpolate
//! each sample's sorted values onto it (matching the
//! `(r − 1) / (n − 1)` ∈ [0, 1] scheme limma uses). Tied assignment then
//! evaluates the mean-of-tied-positions rule on the interpolated reference
//! at each tied position's fractional rank. When all samples share the same
//! non-NaN count, the fractional grid collapses onto the integer rank
//! grid and the implementation is bit-equal to the equal-length-only
//! version we used prior to this change. NaN cells stay NaN.

use ndarray::Array2;

use crate::data::GroupMapping;
use crate::normalize::types::NormalizationError;

/// Returns `(out, scale_marker, nan_cells_in, nan_cells_out)`. `scale_marker`
/// is the mean reference value across all ranks — exposed purely so the
/// dispatcher log line has a single magnitude field consistent with the
/// other methods. It's NaN for the all-NaN-reference case (which errors
/// before returning).
pub fn apply_quantile(
    raw: &Array2<f64>,
    _mapping: &GroupMapping,
) -> Result<(Array2<f64>, f64, usize, usize), NormalizationError> {
    let (n_features, n_samples) = (raw.nrows(), raw.ncols());

    // Step 1: per-sample sorted (feature_idx, value) ascending — non-NaN only.
    let mut sorted_per_sample: Vec<Vec<(usize, f64)>> = Vec::with_capacity(n_samples);
    for j in 0..n_samples {
        let mut indexed: Vec<(usize, f64)> = (0..n_features)
            .map(|i| (i, raw[[i, j]]))
            .filter(|(_, v)| !v.is_nan())
            .collect();
        indexed.sort_by(|a, b| a.1.partial_cmp(&b.1).expect("non-NaN compare"));
        sorted_per_sample.push(indexed);
    }

    // Common rank grid size: K = max(n_j). When all samples share the same
    // non-NaN count, every per-sample fractional rank position lands on an
    // integer grid index, so the interpolation paths below collapse to direct
    // lookups and the output is bit-equal to the equal-length-only version
    // we shipped before this change. When n_j varies (heterogeneous
    // missingness), this grid lets us still place every sample's
    // smallest-to-largest range onto the SAME [0, 1] rank axis instead of
    // letting longer samples dominate the high ranks (the bug this change
    // closes).
    let k_grid = sorted_per_sample.iter().map(|s| s.len()).max().unwrap_or(0);
    if k_grid == 0 {
        return Err(NormalizationError::ReferenceAllNan { method: "Quantile" });
    }

    // Step 2: build the reference on the common grid Q = (0..K)/(K−1).
    // For each Q[k], every non-empty sample contributes its linearly-
    // interpolated value at q = k/(K−1) on its own sorted-rank axis
    // (q = i/(n_j − 1)). This matches limma::normalizeQuantiles's
    // approx-based grid; the difference vs limma is purely in step 3
    // (our mean-of-tied-positions vs limma's average-rank lookup).
    let mut reference: Vec<f64> = vec![0.0; k_grid];
    let mut contrib: Vec<usize> = vec![0; k_grid];
    for sorted in &sorted_per_sample {
        if sorted.is_empty() {
            continue;
        }
        let values: Vec<f64> = sorted.iter().map(|&(_, v)| v).collect();
        for k in 0..k_grid {
            let target_q = grid_q(k, k_grid);
            reference[k] += interp_at(&values, target_q);
            contrib[k] += 1;
        }
    }
    for k in 0..k_grid {
        if contrib[k] == 0 {
            // Unreachable: k_grid > 0 implies at least one non-empty sample,
            // which contributes to every k via interpolation. Kept as a
            // defensive net for refactors.
            return Err(NormalizationError::ReferenceAllNan { method: "Quantile" });
        }
        reference[k] /= contrib[k] as f64;
    }

    // Step 3: every tied value in a sample is assigned the MEAN of the
    // reference values at the fractional-rank positions the tied items
    // collectively span (metabolopan's literal reading of Smyth's #1569
    // remark; diverges from limma `ties=TRUE` / preprocessCore for N ≥ 3,
    // see the module doc-comment):
    //   for a tied block at sorted positions [k, k2) of length N,
    //   output = mean( ref(s / (n_j − 1)) for s in k..k2 )
    // where ref(q) is the linearly-interpolated reference at q ∈ [0, 1].
    // When n_j == k_grid, ref(s / (n_j − 1)) collapses to reference[s], so
    // this reduces to `mean(reference[k..k2])` — exactly the pre-this-change
    // equal-length path. NaN cells pass through untouched.
    let mut out = Array2::<f64>::zeros((n_features, n_samples));
    let mut nan_in = 0usize;
    let mut nan_out = 0usize;
    for j in 0..n_samples {
        for i in 0..n_features {
            if raw[[i, j]].is_nan() {
                nan_in += 1;
                nan_out += 1;
                out[[i, j]] = f64::NAN;
            }
        }
    }
    for (j, indexed) in sorted_per_sample.iter().enumerate() {
        let n_j = indexed.len();
        if n_j == 0 {
            continue;
        }
        let mut k = 0usize;
        while k < n_j {
            let mut k2 = k + 1;
            while k2 < n_j && indexed[k2].1 == indexed[k].1 {
                k2 += 1;
            }
            let mut sum = 0.0;
            for s in k..k2 {
                let target_q = grid_q(s, n_j);
                sum += interp_at(&reference, target_q);
            }
            let tied_ref = sum / (k2 - k) as f64;
            for &(orig_i, _) in &indexed[k..k2] {
                out[[orig_i, j]] = tied_ref;
            }
            k = k2;
        }
    }

    // Magnitude marker: mean of reference values (a single scalar for the log line).
    let scale_marker = reference.iter().sum::<f64>() / reference.len() as f64;
    Ok((out, scale_marker, nan_in, nan_out))
}

/// Map integer rank index `i ∈ [0, n)` to a fractional position `q ∈ [0, 1]`
/// on the `(i − 1) / (n − 1)`-style grid limma uses. Special-cases `n == 1`
/// to `q = 0.0` to avoid 0/0 — for a singleton sample there's only one rank
/// position so the choice of `q` is immaterial.
fn grid_q(i: usize, n: usize) -> f64 {
    if n <= 1 {
        0.0
    } else {
        i as f64 / (n - 1) as f64
    }
}

/// Linear interpolation of a sorted-ascending value array at fractional
/// rank position `q ∈ [0, 1]`. `q = 0` → `values[0]`; `q = 1` → `values[n−1]`;
/// in between → linear interpolation in rank-axis space. Caller guarantees
/// the slice is non-empty and sorted ascending.
fn interp_at(values: &[f64], q: f64) -> f64 {
    let n = values.len();
    debug_assert!(n > 0, "interp_at requires non-empty input");
    if n == 1 {
        return values[0];
    }
    let q = q.clamp(0.0, 1.0);
    let pos = q * (n - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi {
        values[lo]
    } else {
        let frac = pos - lo as f64;
        values[lo] * (1.0 - frac) + values[hi.min(n - 1)] * frac
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::test_support::mapping_for;
    use ndarray::array;

    #[test]
    fn quantile_equalizes_distributions() {
        // No ties, no NaN. Output columns must each be a permutation of the
        // same reference distribution (mean across ranks).
        // Each column has 4 distinct values (no internal ties).
        let raw = array![
            [5.0, 7.0, 11.0],
            [2.0, 1.0, 8.0],
            [3.0, 9.0, 10.0],
            [4.0, 6.0, 12.0]
        ];
        let mapping = mapping_for(&["A", "B", "C"], &["G1"; 3]);
        let (out, _, _, _) = apply_quantile(&raw, &mapping).unwrap();
        let mut col_a: Vec<f64> = (0..4).map(|i| out[[i, 0]]).collect();
        let mut col_b: Vec<f64> = (0..4).map(|i| out[[i, 1]]).collect();
        let mut col_c: Vec<f64> = (0..4).map(|i| out[[i, 2]]).collect();
        for c in [&mut col_a, &mut col_b, &mut col_c] {
            c.sort_by(|a, b| a.partial_cmp(b).unwrap());
        }
        for r in 0..4 {
            assert!(
                (col_a[r] - col_b[r]).abs() < 1e-9,
                "rank {r}: A={} != B={}",
                col_a[r],
                col_b[r]
            );
            assert!(
                (col_a[r] - col_c[r]).abs() < 1e-9,
                "rank {r}: A={} != C={}",
                col_a[r],
                col_c[r]
            );
        }
    }

    #[test]
    fn quantile_nan_passthrough() {
        let raw = array![[1.0, f64::NAN, 3.0], [2.0, 1.0, f64::NAN], [3.0, 2.0, 1.0]];
        let mapping = mapping_for(&["A", "B", "C"], &["G1"; 3]);
        let (out, _, nan_in, nan_out) = apply_quantile(&raw, &mapping).unwrap();
        assert_eq!(nan_in, 2);
        assert_eq!(nan_out, 2);
        assert!(out[[0, 1]].is_nan());
        assert!(out[[1, 2]].is_nan());
    }

    #[test]
    fn quantile_all_nan_matrix_errors() {
        let raw = Array2::from_elem((3, 2), f64::NAN);
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let err = apply_quantile(&raw, &mapping).expect_err("all-NaN matrix must error");
        assert!(matches!(
            err,
            NormalizationError::ReferenceAllNan { method: "Quantile" }
        ));
    }

    #[test]
    fn quantile_two_way_tie_matches_mean_of_two_reference_positions() {
        // Regression guard: 2-way ties produce identical output under the old
        // linear-interpolation code and the mean-of-tied-positions rule
        // (interpolating between ref[0] and ref[1] with frac=0.5 IS the mean).
        // A has a 2-way tie at the bottom; B has no ties.
        //
        // Sample A sorted: [5, 5, 7] → ranks: 5/5 tied at positions 0..2,
        //                                     7 alone at position 2.
        // Sample B sorted: [1, 8, 10].
        // Reference (mean per rank): ref[0]=mean(5,1)=3, ref[1]=mean(5,8)=6.5,
        //                            ref[2]=mean(7,10)=8.5.
        // A's two 5s: tied_ref = mean(ref[0], ref[1]) = (3+6.5)/2 = 4.75.
        // A's 7: ref[2] = 8.5.
        let raw = array![[5.0, 1.0], [5.0, 8.0], [7.0, 10.0]];
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let (out, _, _, _) = apply_quantile(&raw, &mapping).unwrap();
        let expected_tied = (3.0 + 6.5) / 2.0;
        assert!(
            (out[[0, 0]] - expected_tied).abs() < 1e-9,
            "tied A[0] expected {expected_tied}, got {}",
            out[[0, 0]]
        );
        assert!(
            (out[[1, 0]] - expected_tied).abs() < 1e-9,
            "tied A[1] expected {expected_tied}, got {}",
            out[[1, 0]]
        );
        assert!(
            (out[[2, 0]] - 8.5).abs() < 1e-9,
            "non-tied A[2] expected 8.5, got {}",
            out[[2, 0]]
        );
    }

    #[test]
    fn quantile_four_way_tie_uses_mean_of_all_tied_reference_positions() {
        // Bug-revealing case: 4-way tie spanning a NON-LINEAR section of the
        // reference. The pre-2026-05-26 linear-interpolation between ref[1]
        // and ref[2] (frac=0.5) would return 8.0; the mean-of-tied-positions
        // rule returns 7.167. Divergence ≈ 0.83. (limma `ties=TRUE` /
        // preprocessCore use average-rank lookup: this 4-way tie's average
        // rank is 1.5, so they return ref interpolated at 1.5 =
        // (ref[1]+ref[2])/2 = 8.0 — metabolopan's 7.167 diverges from both.)
        //
        // Sample A: [5, 5, 5, 5] — all values tied (k=0, k2=4).
        // Sample B: [1, 10, 11, 12], Sample C: [2, 8, 9, 13].
        // Reference: ref[0]=mean(5,1,2)=8/3, ref[1]=mean(5,10,8)=23/3,
        //            ref[2]=mean(5,11,9)=25/3, ref[3]=mean(5,12,13)=10.
        // mean-of-tied-positions: tied_ref = mean(ref[0..4]) = (8/3+23/3+25/3+10)/4 = 86/12.
        let raw = array![
            [5.0, 1.0, 2.0],
            [5.0, 10.0, 8.0],
            [5.0, 11.0, 9.0],
            [5.0, 12.0, 13.0],
        ];
        let mapping = mapping_for(&["A", "B", "C"], &["G1"; 3]);
        let (out, _, _, _) = apply_quantile(&raw, &mapping).unwrap();
        let expected = 86.0 / 12.0;
        for i in 0..4 {
            assert!(
                (out[[i, 0]] - expected).abs() < 1e-9,
                "mean-over-tied-positions: A[{i}] expected {expected}, got {}",
                out[[i, 0]]
            );
        }
    }

    /// Heterogeneous non-NaN counts: sample A has 3 non-NaN values, sample B
    /// has 5. The fix maps A's largest value to the **reference's largest**
    /// position (via fractional rank `q=1`), not to `reference[2]` as the
    /// pre-fix code did (which would have demoted A's 100th percentile to
    /// the reference's 60th — the bug this change closes). Pins the
    /// limma-canonical `(r − 1)/(n − 1)` grid behaviour.
    #[test]
    fn quantile_unequal_lengths_map_largest_to_reference_top() {
        // A column: [3, NaN, 6, NaN, 8] → 3 non-NaN sorted [3, 6, 8].
        // B column: [1, 2, 4, 7, 9]    → 5 non-NaN.
        let raw = array![
            [3.0, 1.0],
            [f64::NAN, 2.0],
            [6.0, 4.0],
            [f64::NAN, 7.0],
            [8.0, 9.0]
        ];
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let (out, _, _, _) = apply_quantile(&raw, &mapping).unwrap();

        // Reference (K=5):
        //   A interp: q=0→3, 0.25→4.5, 0.5→6, 0.75→7, 1→8
        //   B direct: [1, 2, 4, 7, 9]
        //   ref = [2, 3.25, 5, 7, 8.5]
        // A's 3 sorted positions map at q = 0, 0.5, 1 → ref[0], ref[2], ref[4] = 2, 5, 8.5
        assert!(
            (out[[0, 0]] - 2.0).abs() < 1e-9,
            "A smallest → 2, got {}",
            out[[0, 0]]
        );
        assert!(
            (out[[2, 0]] - 5.0).abs() < 1e-9,
            "A middle  → 5, got {}",
            out[[2, 0]]
        );
        assert!(
            (out[[4, 0]] - 8.5).abs() < 1e-9,
            "A largest → 8.5 (reference top), got {} — pre-fix value was 6 \
             (reference middle), which incorrectly demoted A's 100th percentile",
            out[[4, 0]]
        );
        // NaN cells stay NaN.
        assert!(out[[1, 0]].is_nan());
        assert!(out[[3, 0]].is_nan());
        // B's 5 sorted positions land exactly on reference (K=n_j so no interp).
        for (i, expected) in [(0, 2.0), (1, 3.25), (2, 5.0), (3, 7.0), (4, 8.5)] {
            assert!(
                (out[[i, 1]] - expected).abs() < 1e-9,
                "B[{i}] expected {expected}, got {}",
                out[[i, 1]]
            );
        }
    }

    /// Heterogeneous non-NaN counts AND a tie inside the shorter sample.
    /// Exercises step 3's `sum += interp_at(reference, grid_q(s, n_j))` loop
    /// where each tied position's fractional rank `q = s/(n_j−1)` lands on a
    /// non-integer point of the K-grid. Catches regressions where the tie
    /// average uses the integer-rank reference (the simpler — and wrong —
    /// shortcut).
    #[test]
    fn quantile_unequal_lengths_with_tie_inside_short_sample() {
        // A column: [3, NaN, 3, NaN, 8] → 3 non-NaN sorted [3, 3, 8] (2-way tie).
        // B column: [1, 2, 4, 7, 9]    → 5 non-NaN.
        let raw = array![
            [3.0, 1.0],
            [f64::NAN, 2.0],
            [3.0, 4.0],
            [f64::NAN, 7.0],
            [8.0, 9.0]
        ];
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let (out, _, _, _) = apply_quantile(&raw, &mapping).unwrap();

        // Reference (K=5):
        //   A interp at q=[0, .25, .5, .75, 1]:
        //     q=0   → A[0]=3
        //     q=.25 → mean(A[0], A[1]) = mean(3, 3) = 3
        //     q=.5  → A[1]=3
        //     q=.75 → mean(A[1], A[2]) = mean(3, 8) = 5.5
        //     q=1   → A[2]=8
        //   B direct: [1, 2, 4, 7, 9]
        //   ref = [(3+1)/2, (3+2)/2, (3+4)/2, (5.5+7)/2, (8+9)/2]
        //       = [2, 2.5, 3.5, 6.25, 8.5]
        //
        // A's 2-way tie at sorted positions 0, 1:
        //   q values = 0/(3−1), 1/(3−1) = 0, 0.5
        //   ref(0) = ref[0] = 2; ref(0.5) = ref[2] = 3.5
        //   tied_ref = (2 + 3.5)/2 = 2.75
        // A's sorted position 2: q = 1 → ref[4] = 8.5
        assert!(
            (out[[0, 0]] - 2.75).abs() < 1e-9,
            "tied A[0] → 2.75 (mean of ref[0]=2 and ref[2]=3.5 — mean-of-tied-positions \
             rule on fractional ranks), got {}",
            out[[0, 0]]
        );
        assert!(
            (out[[2, 0]] - 2.75).abs() < 1e-9,
            "tied A[2] → 2.75 (same as tied peer), got {}",
            out[[2, 0]]
        );
        assert!(
            (out[[4, 0]] - 8.5).abs() < 1e-9,
            "untied A[4] → ref[4] = 8.5, got {}",
            out[[4, 0]]
        );
    }

    #[test]
    fn quantile_three_way_tie_linear_reference_unchanged() {
        // Sanity: 3-way tie + locally-linear reference produces the same
        // answer under both old interpolation and the mean-of-tied-positions
        // rule (the integer-mid-rank ref position equals the mean of the three
        // ref positions when they're evenly spaced).
        //
        // Sample A: [5, 5, 5, 10] — 3-way tie at the bottom.
        // Sample B: [1, 2, 3, 4].
        // Reference: ref[0]=3, ref[1]=3.5, ref[2]=4, ref[3]=7.
        // mean(ref[0..3]) = (3 + 3.5 + 4) / 3 = 3.5 == ref[1] (which is what
        // the old code returned). This test guards against accidentally
        // breaking the linear-reference case during the rewrite.
        let raw = array![[5.0, 1.0], [5.0, 2.0], [5.0, 3.0], [10.0, 4.0]];
        let mapping = mapping_for(&["A", "B"], &["G1"; 2]);
        let (out, _, _, _) = apply_quantile(&raw, &mapping).unwrap();
        let expected = 3.5;
        for i in 0..3 {
            assert!(
                (out[[i, 0]] - expected).abs() < 1e-9,
                "linear-reference 3-way tie A[{i}] expected {expected}, got {}",
                out[[i, 0]]
            );
        }
    }
}
