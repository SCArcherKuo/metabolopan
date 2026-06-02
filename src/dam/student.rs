//! Student's t-test, two-tailed, classical equal-variance form.
//!
//! Pooled variance, df = na + nb - 2. NaN-aware on the same boundary
//! conditions as Welch. Shares the `arcsinh` + Pareto-scaling pre-test
//! transform with Welch (see `crate::dam::transforms`); only the
//! variance-assumption / df construction differs between the two.

use statrs::distribution::{ContinuousCDF, StudentsT};

use crate::dam::filter::{nan_aware_mean, nan_aware_var};

/// Two-tailed classical Student's t-test p value over NaN-aware inputs.
///
/// Returns `f64::NAN` when:
/// - either group has fewer than 2 non-NaN values,
/// - the pooled variance is non-positive (both groups constant),
/// - any intermediate (`se²`, `se`) is non-finite or zero.
pub fn student_t_two_tailed(a: &[f64], b: &[f64]) -> f64 {
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
    let na_f = na as f64;
    let nb_f = nb as f64;
    let df = na_f + nb_f - 2.0;
    let sp2 = ((na_f - 1.0) * va + (nb_f - 1.0) * vb) / df;
    if !sp2.is_finite() || sp2 <= 0.0 {
        return f64::NAN;
    }
    let se2 = sp2 * (1.0 / na_f + 1.0 / nb_f);
    if !se2.is_finite() || se2 <= 0.0 {
        return f64::NAN;
    }
    let se = se2.sqrt();
    if !se.is_finite() || se == 0.0 {
        return f64::NAN;
    }
    let t = (ma - mb) / se;
    let dist = match StudentsT::new(0.0, 1.0, df) {
        Ok(d) => d,
        Err(_) => return f64::NAN,
    };
    let cdf_abs = dist.cdf(t.abs());
    let p = 2.0 * (1.0 - cdf_abs);
    p.clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn well_separated_groups_small_p() {
        // Same hand-check inputs as Welch — groups separated by ~10×, n=3 each.
        let a = [10.0, 12.0, 11.0];
        let b = [1.0, 2.0, 1.5];
        let p = student_t_two_tailed(&a, &b);
        assert!(p.is_finite() && p < 0.01, "expected very small p, got {p}");
    }

    #[test]
    fn identical_groups_large_p() {
        let a = [5.0, 5.5, 6.0, 5.2];
        let b = [5.1, 5.4, 5.9, 5.3];
        let p = student_t_two_tailed(&a, &b);
        assert!(p.is_finite() && p > 0.05, "expected p > 0.05, got {p}");
    }

    #[test]
    fn nan_padded_inputs() {
        let a = [10.0, f64::NAN, 12.0, 11.0];
        let b = [1.0, 2.0, f64::NAN, 1.5];
        let p = student_t_two_tailed(&a, &b);
        assert!(
            p.is_finite() && p < 0.01,
            "expected small p with NaN-padding, got {p}"
        );
    }

    #[test]
    fn too_few_samples_returns_nan() {
        assert!(student_t_two_tailed(&[1.0], &[1.0, 2.0, 3.0]).is_nan());
        assert!(student_t_two_tailed(&[1.0, 2.0], &[f64::NAN]).is_nan());
    }

    #[test]
    fn zero_pooled_variance_returns_nan() {
        // Both groups identical constants → pooled variance 0 → NaN.
        let p = student_t_two_tailed(&[3.0, 3.0, 3.0], &[3.0, 3.0, 3.0]);
        assert!(p.is_nan(), "got {p}");
    }
}
