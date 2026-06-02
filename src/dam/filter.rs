//! NaN-aware reductions for DAM statistics.
//!
//! These helpers all treat NaN as "missing": they drop NaN values up front and operate on
//! the remaining values. Empty inputs (all-NaN slices) return either `f64::NAN` (mean,
//! median, variance, IQR) or 0 (`nunique`).

fn drop_nan(values: &[f64]) -> Vec<f64> {
    values.iter().copied().filter(|x| !x.is_nan()).collect()
}

pub fn nan_aware_mean(values: &[f64]) -> f64 {
    let clean = drop_nan(values);
    if clean.is_empty() {
        f64::NAN
    } else {
        clean.iter().sum::<f64>() / clean.len() as f64
    }
}

/// Median of non-NaN values. For even-count input, returns the mean of the two middle
/// values (linear interpolation, same as numpy's default).
pub fn nan_aware_median(values: &[f64]) -> f64 {
    let mut clean = drop_nan(values);
    if clean.is_empty() {
        return f64::NAN;
    }
    clean.sort_by(|a, b| a.partial_cmp(b).expect("non-NaN compare"));
    let n = clean.len();
    if n % 2 == 1 {
        clean[n / 2]
    } else {
        (clean[n / 2 - 1] + clean[n / 2]) / 2.0
    }
}

/// Sample variance with `ddof` degrees of freedom removed (use `ddof=1` for sample
/// variance, `ddof=0` for population). NaN-aware. Returns NaN when fewer than `ddof+1`
/// non-NaN values remain.
pub fn nan_aware_var(values: &[f64], ddof: usize) -> f64 {
    let clean = drop_nan(values);
    if clean.len() <= ddof {
        return f64::NAN;
    }
    let m = clean.iter().sum::<f64>() / clean.len() as f64;
    let sq: f64 = clean.iter().map(|x| (x - m).powi(2)).sum();
    sq / (clean.len() - ddof) as f64
}

/// Interquartile range with linear interpolation (matches numpy default `linear`).
/// NaN-aware. Returns NaN for empty input.
pub fn nan_aware_iqr(values: &[f64]) -> f64 {
    let mut clean = drop_nan(values);
    if clean.is_empty() {
        return f64::NAN;
    }
    clean.sort_by(|a, b| a.partial_cmp(b).expect("non-NaN compare"));
    quantile(&clean, 0.75) - quantile(&clean, 0.25)
}

fn quantile(sorted: &[f64], q: f64) -> f64 {
    // Caller guarantees non-empty + already sorted.
    let n = sorted.len();
    if n == 1 {
        return sorted[0];
    }
    let pos = q * (n - 1) as f64;
    let lo = pos.floor() as usize;
    let hi = pos.ceil() as usize;
    if lo == hi {
        sorted[lo]
    } else {
        let frac = pos - lo as f64;
        sorted[lo] * (1.0 - frac) + sorted[hi] * frac
    }
}

/// Count unique non-NaN values. Uses bit-pattern comparison after NaN filtering — safe
/// because NaN is removed up front.
pub fn nan_aware_nunique(values: &[f64]) -> usize {
    let mut bits: Vec<u64> = drop_nan(values).into_iter().map(|x| x.to_bits()).collect();
    bits.sort_unstable();
    bits.dedup();
    bits.len()
}

/// Pre-filter for DAM. A feature is **kept** iff ALL of the following hold:
///
/// 1. The numerator group has at least 2 non-NaN values, AND
/// 2. The denominator group has at least 2 non-NaN values, AND
/// 3. The combined `numerator ∪ denominator` non-NaN values have `nunique > 1`, AND
/// 4. The combined IQR is strictly positive.
///
/// Checks 1 + 2 (added 2026-05-29) close the "one group entirely NaN" loophole:
/// previously such features passed the prefilter, surfaced NaN through Welch /
/// Student / BM (because their per-test `n < 2` guard kicked in), then occupied
/// an NS slot in `DamResult.features`. Now they are skipped at the prefilter
/// layer and counted in `DamResult.skipped` instead — the per-method NaN guard
/// stays as defence in depth.
pub fn passes_prefilter(numerator: &[f64], denominator: &[f64]) -> bool {
    let n_num = numerator.iter().filter(|x| !x.is_nan()).count();
    let n_den = denominator.iter().filter(|x| !x.is_nan()).count();
    if n_num < 2 || n_den < 2 {
        return false;
    }
    let combined: Vec<f64> = numerator
        .iter()
        .chain(denominator.iter())
        .copied()
        .collect();
    nan_aware_nunique(&combined) > 1 && nan_aware_iqr(&combined) > 0.0
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn mean_nan_handling() {
        assert!(nan_aware_mean(&[]).is_nan());
        assert!(nan_aware_mean(&[f64::NAN, f64::NAN]).is_nan());
        assert_eq!(nan_aware_mean(&[1.0, 2.0, 3.0]), 2.0);
        assert_eq!(nan_aware_mean(&[1.0, f64::NAN, 3.0]), 2.0);
    }

    #[test]
    fn median_even_and_odd() {
        assert_eq!(nan_aware_median(&[1.0, 2.0, 3.0]), 2.0);
        assert_eq!(nan_aware_median(&[1.0, 2.0, 3.0, 4.0]), 2.5);
        assert_eq!(nan_aware_median(&[1.0, f64::NAN, 3.0]), 2.0);
        assert!(nan_aware_median(&[]).is_nan());
    }

    #[test]
    fn var_sample_default() {
        let v = nan_aware_var(&[1.0, 2.0, 3.0, 4.0, 5.0], 1);
        // mean=3, sum_sq=10, /(5-1)=2.5
        assert!((v - 2.5).abs() < 1e-12);
        assert!(nan_aware_var(&[1.0], 1).is_nan());
    }

    #[test]
    fn iqr_basic() {
        // [1,2,3,4,5] → Q1=2, Q3=4, IQR=2
        let i = nan_aware_iqr(&[1.0, 2.0, 3.0, 4.0, 5.0]);
        assert!((i - 2.0).abs() < 1e-12, "got {i}");
        assert_eq!(nan_aware_iqr(&[5.0, 5.0, 5.0]), 0.0);
    }

    #[test]
    fn nunique_dedup() {
        assert_eq!(nan_aware_nunique(&[1.0, 1.0, 2.0, 2.0, 3.0]), 3);
        assert_eq!(nan_aware_nunique(&[1.0, f64::NAN, 1.0]), 1);
        assert_eq!(nan_aware_nunique(&[]), 0);
    }

    #[test]
    fn prefilter_drops_all_equal() {
        assert!(!passes_prefilter(&[1.0, 1.0, 1.0], &[1.0, 1.0]));
    }

    #[test]
    fn prefilter_drops_iqr_zero() {
        // [1,1,1,1,2] → IQR = 0
        assert!(!passes_prefilter(&[1.0, 1.0, 1.0], &[1.0, 2.0]));
    }

    #[test]
    fn prefilter_keeps_real_variation() {
        assert!(passes_prefilter(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]));
    }

    #[test]
    fn prefilter_drops_all_nan() {
        assert!(!passes_prefilter(&[f64::NAN; 3], &[f64::NAN; 3]));
    }

    /// 2026-05-29 tightening: numerator with < 2 non-NaN values is dropped at
    /// the prefilter layer (previously the feature passed and the t-test /
    /// BM produced NaN that propagated as an NS slot in `DamResult.features`).
    #[test]
    fn prefilter_drops_numerator_below_two_non_nan() {
        // 1 non-NaN in numerator, 3 in denominator.
        assert!(!passes_prefilter(
            &[10.0, f64::NAN, f64::NAN],
            &[5.0, 6.0, 7.0]
        ));
        // 0 non-NaN in numerator.
        assert!(!passes_prefilter(
            &[f64::NAN, f64::NAN, f64::NAN],
            &[5.0, 6.0, 7.0]
        ));
    }

    /// Symmetric counterpart of the numerator check.
    #[test]
    fn prefilter_drops_denominator_below_two_non_nan() {
        assert!(!passes_prefilter(&[1.0, 2.0, 3.0], &[10.0, f64::NAN]));
        assert!(!passes_prefilter(&[1.0, 2.0, 3.0], &[f64::NAN, f64::NAN]));
    }

    /// Boundary: exactly 2 non-NaN per group is the minimum that passes the
    /// per-group check; combined nunique / IQR still apply.
    #[test]
    fn prefilter_keeps_exactly_two_non_nan_per_group_with_variation() {
        assert!(passes_prefilter(
            &[1.0, 2.0, f64::NAN],
            &[3.0, 4.0, f64::NAN]
        ));
    }
}
