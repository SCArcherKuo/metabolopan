//! Welch's t-test, two-tailed, with NaN-aware mean/var and Satterthwaite df.
//!
//! The shared `arcsinh` + Pareto-scaling pre-test transform lives in
//! `crate::dam::transforms` (both Welch and Student use it).

use statrs::distribution::{ContinuousCDF, StudentsT};

use crate::dam::filter::{nan_aware_mean, nan_aware_var};

/// Two-tailed Welch's t-test p value over NaN-aware inputs. Returns NaN when either
/// group has fewer than 2 non-NaN values OR pooled variance is non-positive.
pub fn welch_t_two_tailed(a: &[f64], b: &[f64]) -> f64 {
    let na = a.iter().filter(|x| !x.is_nan()).count();
    let nb = b.iter().filter(|x| !x.is_nan()).count();
    if na < 2 || nb < 2 {
        return f64::NAN;
    }
    let ma = nan_aware_mean(a);
    let mb = nan_aware_mean(b);
    let va = nan_aware_var(a, 1);
    let vb = nan_aware_var(b, 1);
    if va.is_nan() || vb.is_nan() {
        return f64::NAN;
    }
    let se2 = va / na as f64 + vb / nb as f64;
    if !se2.is_finite() || se2 <= 0.0 {
        return f64::NAN;
    }
    let t = (ma - mb) / se2.sqrt();
    let na_f = na as f64;
    let nb_f = nb as f64;
    let df_num = se2.powi(2);
    let df_den = (va / na_f).powi(2) / (na_f - 1.0) + (vb / nb_f).powi(2) / (nb_f - 1.0);
    if df_den <= 0.0 || !df_den.is_finite() {
        return f64::NAN;
    }
    let df = df_num / df_den;
    let dist = match StudentsT::new(0.0, 1.0, df) {
        Ok(d) => d,
        Err(_) => return f64::NAN,
    };
    // Two-tailed: 2 * (1 - cdf(|t|))
    let cdf_abs = dist.cdf(t.abs());
    let p = 2.0 * (1.0 - cdf_abs);
    p.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn well_separated_groups_small_p() {
        // Hand-check: groups separated by ~10×, n=3 each. p must be small.
        let a = [10.0, 12.0, 11.0];
        let b = [1.0, 2.0, 1.5];
        let p = welch_t_two_tailed(&a, &b);
        assert!(p.is_finite() && p < 0.01, "expected very small p, got {p}");
    }

    #[test]
    fn identical_groups_large_p() {
        let a = [5.0, 5.5, 6.0, 5.2];
        let b = [5.1, 5.4, 5.9, 5.3];
        let p = welch_t_two_tailed(&a, &b);
        assert!(p.is_finite() && p > 0.05, "expected p > 0.05, got {p}");
    }

    #[test]
    fn nan_padded_inputs() {
        // Same data as well-separated test but padded with NaN.
        let a = [10.0, f64::NAN, 12.0, 11.0];
        let b = [1.0, 2.0, f64::NAN, 1.5];
        let p = welch_t_two_tailed(&a, &b);
        assert!(
            p.is_finite() && p < 0.01,
            "expected small p with NaN-padding, got {p}"
        );
    }

    #[test]
    fn too_few_samples_returns_nan() {
        assert!(welch_t_two_tailed(&[1.0], &[1.0, 2.0, 3.0]).is_nan());
        assert!(welch_t_two_tailed(&[1.0, 2.0], &[f64::NAN]).is_nan());
    }

    #[test]
    fn zero_variance_returns_nan() {
        // Both groups identical constants → variance 0 in each → se2 = 0 → NaN
        let p = welch_t_two_tailed(&[3.0, 3.0, 3.0], &[3.0, 3.0, 3.0]);
        assert!(p.is_nan(), "got {p}");
    }
}
