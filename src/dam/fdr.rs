//! FDR correction for Stage 2 DAM and Stage 3 ORA.
//!
//! Provides two NaN-aware corrections — Benjamini–Hochberg (BH) and
//! Benjamini–Yekutieli (BY) — plus a dispatcher `adjust_pvalues(p, method)`
//! that callers use to stay agnostic of the choice. NaN inputs pass through
//! as NaN and are excluded from the family size `m`.

use serde::{Deserialize, Serialize};

/// FDR correction method.
///
/// `BenjaminiHochberg` is the default — it matches R's `p.adjust` /
/// MetaboAnalyst conventions and is the literature norm for metabolomics
/// DAM. `BenjaminiYekutieli` adds the harmonic factor `c(m) = Σ 1/i` and
/// is the safer choice when entries share dependence (especially Stage 3
/// pathway/module ORA, where compounds are shared across entries).
/// `NoCorrection` returns raw p-values unchanged — exposed ONLY on the
/// Stage 3 ORA setup radio for exploratory runs; the Stage 2 DAM UI does
/// not expose this variant and `apply_snapshot` defensively coerces it
/// to `BenjaminiHochberg` so a hand-crafted snapshot can't push raw p
/// through a ~13 k-feature DAM run. Short label is "None"; serialised
/// to JSON as `"NoCorrection"`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default, Serialize, Deserialize)]
pub enum FdrMethod {
    #[default]
    BenjaminiHochberg,
    BenjaminiYekutieli,
    NoCorrection,
}

impl FdrMethod {
    /// Short label used in CSV `# FDR:` tags and the settings-summary line.
    ///
    /// `NoCorrection` renders as `NoCorrection`, matching the value serde
    /// writes for this variant in a saved session, so an exported CSV and a
    /// settings file name the choice with the same word. It was `None`, which
    /// is the one method word that reads as *absent* rather than as a named
    /// choice — a hazard wherever it is read without the radio in view.
    pub fn short_label(self) -> &'static str {
        match self {
            FdrMethod::BenjaminiHochberg => "BH",
            FdrMethod::BenjaminiYekutieli => "BY",
            FdrMethod::NoCorrection => "NoCorrection",
        }
    }

    /// The noun for the significance QUANTITY this method produces.
    ///
    /// BH and BY produce q-values, which this project names `FDR` throughout
    /// (the term both manuals use almost exclusively). `NoCorrection` produces
    /// no second quantity at all — the value is the raw p-value — so a string
    /// that names it must say `p-value`. Distinct from [`Self::short_label`],
    /// which names the method CHOICE: a control labelled "FDR correction: none"
    /// is true, a column headed `FDR` holding an uncorrected p-value is not.
    pub fn metric_label(self) -> &'static str {
        match self {
            FdrMethod::BenjaminiHochberg | FdrMethod::BenjaminiYekutieli => "FDR",
            FdrMethod::NoCorrection => "p-value",
        }
    }
}

/// Apply the chosen FDR method to `p_values`. Dispatches to
/// [`benjamini_hochberg`] or [`benjamini_yekutieli`]; both share the same
/// NaN-aware contract (NaN ↔ NaN, NaN excluded from `m`). `NoCorrection`
/// returns a copy of `p_values` verbatim (NaNs included) — caller
/// semantics for downstream FDR-keyed fields become "FDR == p".
pub fn adjust_pvalues(p_values: &[f64], method: FdrMethod) -> Vec<f64> {
    match method {
        FdrMethod::BenjaminiHochberg => benjamini_hochberg(p_values),
        FdrMethod::BenjaminiYekutieli => benjamini_yekutieli(p_values),
        FdrMethod::NoCorrection => p_values.to_vec(),
    }
}

/// Apply Benjamini-Hochberg FDR correction. Returns a vector of the same
/// length as `p_values`, with the i-th output being the BH-adjusted
/// q-value for the i-th input (NaN ↔ NaN). The formula is
/// `adj[k] = p[k] * m / rank` followed by monotone-min from the right
/// and cap at 1.0; no harmonic factor (the difference vs BY).
pub fn benjamini_hochberg(p_values: &[f64]) -> Vec<f64> {
    let mut out = vec![f64::NAN; p_values.len()];
    let mut indexed: Vec<(usize, f64)> = p_values
        .iter()
        .enumerate()
        .filter_map(|(i, &p)| if p.is_nan() { None } else { Some((i, p)) })
        .collect();
    if indexed.is_empty() {
        return out;
    }
    indexed.sort_by(|(_, a), (_, b)| a.partial_cmp(b).expect("non-NaN compare"));
    let m = indexed.len();
    let m_f = m as f64;

    let mut raw_adj: Vec<f64> = indexed
        .iter()
        .enumerate()
        .map(|(rank0, &(_, p))| {
            let rank = (rank0 + 1) as f64;
            p * m_f / rank
        })
        .collect();

    for i in (0..raw_adj.len() - 1).rev() {
        if raw_adj[i] > raw_adj[i + 1] {
            raw_adj[i] = raw_adj[i + 1];
        }
    }
    for v in raw_adj.iter_mut() {
        if *v > 1.0 {
            *v = 1.0;
        }
    }

    for (sorted_pos, &(orig_idx, _)) in indexed.iter().enumerate() {
        out[orig_idx] = raw_adj[sorted_pos];
    }
    out
}

/// Apply Benjamini-Yekutieli FDR correction. Returns a vector of the same length as
/// `p_values`, with the i-th output being the BY-adjusted q-value for the i-th input
/// (NaN ↔ NaN). The harmonic factor `c(m) = Σ_{i=1}^{m} 1/i` is computed exactly over
/// the m = count of non-NaN inputs.
pub fn benjamini_yekutieli(p_values: &[f64]) -> Vec<f64> {
    let mut out = vec![f64::NAN; p_values.len()];
    let mut indexed: Vec<(usize, f64)> = p_values
        .iter()
        .enumerate()
        .filter_map(|(i, &p)| if p.is_nan() { None } else { Some((i, p)) })
        .collect();
    if indexed.is_empty() {
        return out;
    }
    indexed.sort_by(|(_, a), (_, b)| a.partial_cmp(b).expect("non-NaN compare"));
    let m = indexed.len();
    let m_f = m as f64;
    let c_m: f64 = (1..=m).map(|i| 1.0 / i as f64).sum();

    // BY adjusted: q_i = p_i * m * c(m) / rank_i; then running-min from the top.
    let mut raw_adj: Vec<f64> = indexed
        .iter()
        .enumerate()
        .map(|(rank0, &(_, p))| {
            let rank = (rank0 + 1) as f64;
            p * m_f * c_m / rank
        })
        .collect();

    // Enforce monotone non-decreasing (in p order) by cumulative-min from the right.
    for i in (0..raw_adj.len() - 1).rev() {
        if raw_adj[i] > raw_adj[i + 1] {
            raw_adj[i] = raw_adj[i + 1];
        }
    }
    // Cap at 1.0.
    for v in raw_adj.iter_mut() {
        if *v > 1.0 {
            *v = 1.0;
        }
    }

    // Scatter back into original positions.
    for (sorted_pos, &(orig_idx, _)) in indexed.iter().enumerate() {
        out[orig_idx] = raw_adj[sorted_pos];
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Golden values from R: `p.adjust(c(0.001, 0.01, 0.03, 0.05, 0.20), method="BY")`
    /// = c(0.011416, 0.057082, 0.114164, 0.142705, 0.456657)
    #[test]
    fn by_matches_r_golden() {
        let p = vec![0.001, 0.01, 0.03, 0.05, 0.20];
        let q = benjamini_yekutieli(&p);
        let expected = [0.011416, 0.057082, 0.114164, 0.142705, 0.456657];
        for (a, e) in q.iter().zip(expected.iter()) {
            assert!((a - e).abs() < 1e-4, "got {a}, expected {e}");
        }
    }

    #[test]
    fn nan_passes_through() {
        let p = vec![0.01, f64::NAN, 0.05];
        let q = benjamini_yekutieli(&p);
        assert!(q[1].is_nan());
        assert!(q[0].is_finite() && q[2].is_finite());
        // Only 2 non-NaN values → c(m=2) = 1 + 0.5 = 1.5
        // sorted: 0.01 (rank 1), 0.05 (rank 2). adj = p * 2 * 1.5 / rank
        // 0.01 → 0.01*3/1 = 0.03; 0.05 → 0.05*3/2 = 0.075
        // monotone enforce (already monotone)
        assert!((q[0] - 0.03).abs() < 1e-9, "got {}", q[0]);
        assert!((q[2] - 0.075).abs() < 1e-9, "got {}", q[2]);
    }

    #[test]
    fn empty_input() {
        assert!(benjamini_yekutieli(&[]).is_empty());
    }

    #[test]
    fn all_nan_input() {
        let q = benjamini_yekutieli(&[f64::NAN; 3]);
        assert!(q.iter().all(|v| v.is_nan()));
    }

    /// Golden values from R: `p.adjust(c(0.001, 0.01, 0.03, 0.05, 0.20), method="BH")`
    /// = c(0.005, 0.025, 0.05, 0.0625, 0.20)
    #[test]
    fn bh_matches_r_golden() {
        let p = vec![0.001, 0.01, 0.03, 0.05, 0.20];
        let q = benjamini_hochberg(&p);
        let expected = [0.005, 0.025, 0.05, 0.0625, 0.20];
        for (a, e) in q.iter().zip(expected.iter()) {
            assert!((a - e).abs() < 1e-9, "got {a}, expected {e}");
        }
    }

    #[test]
    fn bh_nan_passes_through() {
        let p = vec![0.01, f64::NAN, 0.05];
        let q = benjamini_hochberg(&p);
        assert!(q[1].is_nan());
        assert!(q[0].is_finite() && q[2].is_finite());
        // Only 2 non-NaN values → m = 2; no harmonic factor.
        // sorted: 0.01 (rank 1), 0.05 (rank 2). adj = p * 2 / rank
        // 0.01 → 0.02; 0.05 → 0.05. monotone enforce (already monotone).
        assert!((q[0] - 0.02).abs() < 1e-9, "got {}", q[0]);
        assert!((q[2] - 0.05).abs() < 1e-9, "got {}", q[2]);
    }

    #[test]
    fn bh_empty_input() {
        assert!(benjamini_hochberg(&[]).is_empty());
    }

    #[test]
    fn dispatcher_matches_underlying_functions() {
        let p = vec![0.001, 0.01, 0.03, 0.05, 0.20];
        let via_dispatch_bh = adjust_pvalues(&p, FdrMethod::BenjaminiHochberg);
        let direct_bh = benjamini_hochberg(&p);
        let via_dispatch_by = adjust_pvalues(&p, FdrMethod::BenjaminiYekutieli);
        let direct_by = benjamini_yekutieli(&p);
        assert_eq!(via_dispatch_bh, direct_bh);
        assert_eq!(via_dispatch_by, direct_by);
    }

    #[test]
    fn bh_and_by_differ_on_same_input() {
        // Guards against accidental dispatcher aliasing — at least one
        // position must differ by > 1e-12. (BY's c(m) factor inflates the
        // q-values relative to BH on this input.)
        let p = vec![0.001, 0.01, 0.03, 0.05, 0.20];
        let bh = adjust_pvalues(&p, FdrMethod::BenjaminiHochberg);
        let by = adjust_pvalues(&p, FdrMethod::BenjaminiYekutieli);
        let any_differ = bh.iter().zip(by.iter()).any(|(a, b)| (a - b).abs() > 1e-12);
        assert!(
            any_differ,
            "BH and BY must produce different outputs; bh={bh:?} by={by:?}"
        );
    }

    #[test]
    fn fdr_method_default_is_bh() {
        assert_eq!(FdrMethod::default(), FdrMethod::BenjaminiHochberg);
    }

    #[test]
    fn fdr_method_short_labels() {
        assert_eq!(FdrMethod::BenjaminiHochberg.short_label(), "BH");
        assert_eq!(FdrMethod::BenjaminiYekutieli.short_label(), "BY");
        assert_eq!(FdrMethod::NoCorrection.short_label(), "NoCorrection");
    }

    #[test]
    fn fdr_method_metric_labels() {
        // The quantity, not the choice: BH/BY produce q-values (named FDR
        // project-wide), NoCorrection produces a raw p-value.
        assert_eq!(FdrMethod::BenjaminiHochberg.metric_label(), "FDR");
        assert_eq!(FdrMethod::BenjaminiYekutieli.metric_label(), "FDR");
        assert_eq!(FdrMethod::NoCorrection.metric_label(), "p-value");
    }

    #[test]
    fn no_correction_returns_p_values_verbatim_including_nan() {
        let p = vec![0.001, 0.05, 0.5, f64::NAN, 1.0];
        let q = adjust_pvalues(&p, FdrMethod::NoCorrection);
        assert_eq!(q.len(), p.len());
        for (a, b) in q.iter().zip(p.iter()) {
            if b.is_nan() {
                assert!(a.is_nan());
            } else {
                assert_eq!(a, b);
            }
        }
    }
}
