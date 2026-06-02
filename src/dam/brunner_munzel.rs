//! Brunner-Munzel test + Cliff's δ effect size. Both NaN-aware.
//!
//! BM behaviour matches SciPy `brunnermunzel(distribution='t')` — mid-ranks, Welch-
//! Satterthwaite-like df, two-tailed p via Student-t.

use statrs::distribution::{ContinuousCDF, StudentsT};

/// Two-tailed Brunner-Munzel p value. Returns NaN when either group has < 2 non-NaN
/// values or when both groups are identical constants (variance = 0).
pub fn brunner_munzel_two_tailed(a: &[f64], b: &[f64]) -> f64 {
    let xa: Vec<f64> = a.iter().copied().filter(|x| !x.is_nan()).collect();
    let xb: Vec<f64> = b.iter().copied().filter(|x| !x.is_nan()).collect();
    let na = xa.len();
    let nb = xb.len();
    if na < 2 || nb < 2 {
        return f64::NAN;
    }

    // Combined sample, with provenance flags.
    let mut combined: Vec<(f64, u8)> = Vec::with_capacity(na + nb);
    for &v in &xa {
        combined.push((v, 0));
    }
    for &v in &xb {
        combined.push((v, 1));
    }
    let mid_ranks_combined = mid_ranks(&combined.iter().map(|&(v, _)| v).collect::<Vec<_>>());

    // Per-group ranks within combined.
    let mut ranks_a_combined: Vec<f64> = Vec::with_capacity(na);
    let mut ranks_b_combined: Vec<f64> = Vec::with_capacity(nb);
    for (i, &(_, g)) in combined.iter().enumerate() {
        if g == 0 {
            ranks_a_combined.push(mid_ranks_combined[i]);
        } else {
            ranks_b_combined.push(mid_ranks_combined[i]);
        }
    }

    // Per-group ranks within group only.
    let ranks_a_within = mid_ranks(&xa);
    let ranks_b_within = mid_ranks(&xb);

    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
    let r1_bar = mean(&ranks_a_combined);
    let r2_bar = mean(&ranks_b_combined);
    // Mean of within-group ranks is trivially (n_k + 1) / 2; the explicit means
    // (`s1_bar`, `s2_bar`) would just duplicate that constant, so we use the closed
    // form directly in the variance formula below.

    let na_f = na as f64;
    let nb_f = nb as f64;
    let nx = na_f + nb_f;

    // s_x^2 = sum (R_ki - S_ki - r_bar_k + (n_k+1)/2)^2 / (n_k - 1)
    let s1_sq = ranks_a_combined
        .iter()
        .zip(ranks_a_within.iter())
        .map(|(&r, &s)| {
            let v = r - s - r1_bar + (na_f + 1.0) / 2.0;
            v * v
        })
        .sum::<f64>()
        / (na_f - 1.0);
    let s2_sq = ranks_b_combined
        .iter()
        .zip(ranks_b_within.iter())
        .map(|(&r, &s)| {
            let v = r - s - r2_bar + (nb_f + 1.0) / 2.0;
            v * v
        })
        .sum::<f64>()
        / (nb_f - 1.0);

    // SciPy `brunnermunzel(distribution='t')` formula: the W denominator
    // inside sqrt is `nx*Sx + ny*Sy` (NOT `(nx+ny) * (Sx/nx + Sy/ny)`,
    // which inflates |W| by sqrt(N/2) for equal n and was the pre-2026-05-26
    // bug). See SciPy v1.14 _stats_py.py::brunnermunzel and lawstat R pkg.
    let s_pool = na_f * s1_sq + nb_f * s2_sq;
    if !s_pool.is_finite() || s_pool <= 0.0 {
        return f64::NAN;
    }
    let w = na_f * nb_f * (r2_bar - r1_bar) / (nx * s_pool.sqrt());

    // Welch-Satterthwaite-like df, matching SciPy's brunnermunzel df formula.
    let df_num = (s1_sq / na_f + s2_sq / nb_f).powi(2);
    let df_den = (s1_sq / na_f).powi(2) / (na_f - 1.0) + (s2_sq / nb_f).powi(2) / (nb_f - 1.0);
    if df_den <= 0.0 || !df_den.is_finite() {
        return f64::NAN;
    }
    let df = df_num / df_den;
    let dist = match StudentsT::new(0.0, 1.0, df) {
        Ok(d) => d,
        Err(_) => return f64::NAN,
    };
    let cdf_abs = dist.cdf(w.abs());
    let p = 2.0 * (1.0 - cdf_abs);
    p.clamp(0.0, 1.0)
}

/// Cliff's δ: `(gt - lt) / (n * m)` where `gt` = strict-greater pairs, `lt` = strict-less.
/// Returns NaN when either group has 0 non-NaN values.
pub fn cliffs_delta(a: &[f64], b: &[f64]) -> f64 {
    let xa: Vec<f64> = a.iter().copied().filter(|x| !x.is_nan()).collect();
    let xb: Vec<f64> = b.iter().copied().filter(|x| !x.is_nan()).collect();
    if xa.is_empty() || xb.is_empty() {
        return f64::NAN;
    }
    let mut gt: i64 = 0;
    let mut lt: i64 = 0;
    for &x in &xa {
        for &y in &xb {
            if x > y {
                gt += 1;
            } else if x < y {
                lt += 1;
            }
        }
    }
    let n = xa.len() as i64;
    let m = xb.len() as i64;
    (gt - lt) as f64 / (n * m) as f64
}

/// Mid-ranks (a.k.a. tied ranks): tied values get the mean of the ranks they would have
/// occupied. Equivalent to scipy.stats.rankdata(values, method='average').
fn mid_ranks(values: &[f64]) -> Vec<f64> {
    let n = values.len();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| values[i].partial_cmp(&values[j]).expect("non-NaN compare"));
    let mut ranks = vec![0.0; n];
    let mut i = 0;
    while i < n {
        let mut j = i + 1;
        while j < n && values[idx[j]] == values[idx[i]] {
            j += 1;
        }
        // Indices [i, j) are tied at the same value.
        let avg_rank = (i + j + 1) as f64 / 2.0;
        for k in i..j {
            ranks[idx[k]] = avg_rank;
        }
        i = j;
    }
    ranks
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bm_well_separated_with_one_overlap_gives_small_p() {
        // [1..5] vs [4.5, 6..9] — single overlap so BM variance is non-zero. Groups
        // are still very well separated; p should be small.
        let a = [1.0, 2.0, 3.0, 4.0, 5.0];
        let b = [4.5, 6.0, 7.0, 8.0, 9.0];
        let p = brunner_munzel_two_tailed(&a, &b);
        assert!(p.is_finite(), "p must be finite, got {p}");
        assert!(p < 0.01, "expected p < 0.01, got {p}");
    }

    #[test]
    fn bm_perfectly_stratified_returns_nan() {
        // [1..5] vs [6..10] — no overlap at all. BM variance estimator is exactly 0
        // (perfectly stratified), so the test is undefined. SciPy returns
        // pvalue=nan for this case (BrunnerMunzelResult(statistic=-inf, pvalue=nan));
        // we match that behavior.
        let a = [1.0, 2.0, 3.0, 4.0, 5.0];
        let b = [6.0, 7.0, 8.0, 9.0, 10.0];
        let p = brunner_munzel_two_tailed(&a, &b);
        assert!(
            p.is_nan(),
            "expected NaN for perfectly stratified data, got {p}"
        );
    }

    #[test]
    fn bm_matches_scipy_brunnermunzel_distribution_t() {
        // [1..5] vs [2..6] — heavy overlap. Hand-computed against the
        // canonical SciPy formula (`nx*Sx + ny*Sy` inside sqrt):
        //   r̄_x = 4.6, r̄_y = 6.4, Sx = Sy = 2.05, df = 8.0
        //   W = 5·5·(6.4 - 4.6) / (10 · sqrt(5·2.05 + 5·2.05))
        //     = 45 / (10 · sqrt(20.5)) ≈ 0.9939
        //   two-tailed p at df=8, |W|=0.9939 ≈ 0.349
        // SciPy `brunnermunzel(distribution='t')` and R `lawstat::
        // brunner.munzel.test` both return p ≈ 0.349. Pre-2026-05-26 the
        // Rust code returned p ≈ 0.155 (inflated significance by sqrt(N/2)).
        let a = [1.0, 2.0, 3.0, 4.0, 5.0];
        let b = [2.0, 3.0, 4.0, 5.0, 6.0];
        let p = brunner_munzel_two_tailed(&a, &b);
        assert!(
            (p - 0.349).abs() < 0.01,
            "expected p ≈ 0.349 (SciPy / R lawstat parity), got {p}"
        );
    }

    #[test]
    fn bm_too_few_samples_nan() {
        assert!(brunner_munzel_two_tailed(&[1.0], &[1.0, 2.0, 3.0]).is_nan());
        assert!(brunner_munzel_two_tailed(&[1.0, f64::NAN], &[1.0, 2.0]).is_nan());
    }

    #[test]
    fn cliffs_delta_extremes() {
        let a = [1.0, 2.0, 3.0];
        let b = [4.0, 5.0, 6.0];
        let d = cliffs_delta(&a, &b);
        assert!((d - (-1.0)).abs() < 1e-12, "got {d}");
        let d2 = cliffs_delta(&b, &a);
        assert!((d2 - 1.0).abs() < 1e-12, "got {d2}");
    }

    #[test]
    fn cliffs_delta_overlap() {
        let a = [1.0, 2.0, 3.0];
        let b = [2.0, 3.0, 4.0];
        // gt: (2>2)=0, (2>3)=0, (2>4)=0 — for x=1: 0; x=2: 0; x=3: 1 (3>2)
        // lt: x=1 vs all: 3; x=2 vs (3,4): 2; x=3 vs (4): 1 → total 6
        // gt=1, lt=6 → (1-6)/9 = -5/9 ≈ -0.5556
        let d = cliffs_delta(&a, &b);
        assert!((d - (-5.0 / 9.0)).abs() < 1e-9, "got {d}");
    }

    #[test]
    fn mid_ranks_with_ties() {
        // [1, 2, 2, 3] → ranks: 1, 2.5, 2.5, 4
        let r = mid_ranks(&[1.0, 2.0, 2.0, 3.0]);
        assert_eq!(r, vec![1.0, 2.5, 2.5, 4.0]);
    }
}
