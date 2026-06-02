//! DAM runner: end-to-end pre-filter + per-method statistics + user-selected
//! FDR correction (BH default, BY opt-in).

use anyhow::{Result, anyhow};
use std::collections::HashSet;
use std::sync::mpsc;
use tracing::{debug, info};

use crate::dam::brunner_munzel::{brunner_munzel_two_tailed, cliffs_delta};
use crate::dam::fdr::{FdrMethod, adjust_pvalues};
use crate::dam::filter::{nan_aware_mean, nan_aware_median, nan_aware_var, passes_prefilter};
use crate::dam::student::student_t_two_tailed;
use crate::dam::transforms::arcsinh_in_place;
use crate::dam::types::{DamFeature, DamMethod, DamProgress, DamResult, FcBasis, Trend};
use crate::dam::welch::welch_t_two_tailed;
use crate::data::{GroupMapping, MetabolomicsTable};
use crate::dedup::{DedupReport, run_dedup};
use crate::normalize::NormalizationConfig;

/// Configuration for a single [`run_dam`] invocation: the statistical method,
/// the column-wise normalization, the three pre-test filters / transforms, and
/// the FDR correction. Bundled into one struct so the three adjacent `bool`
/// flags (`drop_unknown` / `dedup_enabled` / `log_transform`) cannot be
/// transposed at the call site — the footgun the prior 11-positional-arg
/// signature carried under `#[allow(clippy::too_many_arguments)]`.
///
/// Field order mirrors the `run_dam` processing sequence and
/// `SessionSettings`'s declaration order. No `Default` impl: the canonical
/// defaults live in `SessionSettings::default()` (single source of truth).
#[derive(Debug, Clone)]
pub struct DamConfig {
    /// Welch (parametric, heteroscedastic) / Student (parametric, homoscedastic)
    /// / Brunner-Munzel (non-parametric).
    pub method: DamMethod,
    /// Column-wise sample normalization applied to `intensity_raw` before
    /// stats. Owned by value (the type is `Clone`, not `Copy`).
    pub normalization: NormalizationConfig,
    /// Drop features whose `inchikey.is_none()` before pre-filter / stats / FDR.
    pub drop_unknown: bool,
    /// Apply the cascade-only InChIKey deduplication mask in the per-feature loop.
    pub dedup_enabled: bool,
    /// Apply the `arcsinh` generalised-log variance-stabilisation step on the
    /// Welch / Student paths (BM ignores it — arcsinh is monotone, so ranks are
    /// invariant under it).
    pub log_transform: bool,
    /// FDR correction applied to the per-feature p vector (BH default, BY opt-in).
    pub fdr_method: FdrMethod,
}

/// Per-feature DAM analysis between `numerator` and `denominator` groups using the
/// chosen method. `drop_unknown == true` removes every feature whose `inchikey.is_none()`
/// before pre-filter / stats / FDR, so the Unknown features do not inflate the
/// multiple-testing burden.
///
/// `log_transform` gates the project's "generalised log" (`arcsinh`)
/// variance-stabilisation step on the Welch / Student paths. `true` =
/// apply `arcsinh_in_place` to the union of numerator + denominator
/// values BEFORE the t-test (matches pre-2026-05-27 behaviour). `false`
/// = skip arcsinh and pass the working matrix to the t-test directly.
/// BM (rank-based) ignores the flag — arcsinh is monotone so ranks are
/// invariant under it. The earlier hard-coded `pareto_scale_in_place`
/// step is gone (no-op for univariate t-statistics — verified empirically
/// at f64 precision; see add-log-transform-and-scaling design D1).
///
/// `fdr_method` selects the FDR correction applied to the per-feature p
/// vector (BH for the literature default and cross-tool reproducibility, BY
/// for the more conservative dependence-tolerant correction). The chosen
/// method is recorded on `DamResult.fdr_method` so the CSV exporter and
/// downstream renderers can surface the choice.
///
/// When `progress_tx` is `Some`, the runner emits a [`DamProgress`] event after every
/// feature so the UI can drive a progress bar. `total` reports the input feature count
/// (i.e. `table.features.len()`); `completed` advances regardless of whether each
/// feature was kept, dropped as Unknown, or skipped by pre-filter, so the bar moves
/// monotonically through the input. The channel is best-effort — `try_send` errors are
/// ignored so a dropped UI receiver never aborts a fetch.
pub async fn run_dam(
    table: &mut MetabolomicsTable,
    mapping: &GroupMapping,
    numerator: &str,
    denominator: &str,
    config: &DamConfig,
    progress_tx: Option<mpsc::Sender<DamProgress>>,
) -> Result<DamResult> {
    // Destructure the config into locals so the body below reads exactly as it
    // did under the prior positional signature — this change is pure
    // re-packaging of the same six values in the same order, so the output is
    // bit-identical. `normalization` stays a borrow (`NormalizationConfig` is
    // `Clone`, not `Copy`); the two enums and three bools are `Copy`.
    let method = config.method;
    let normalization = &config.normalization;
    let drop_unknown = config.drop_unknown;
    let dedup_enabled = config.dedup_enabled;
    let log_transform = config.log_transform;
    let fdr_method = config.fdr_method;

    // Validate group names.
    let groups = mapping.groups();
    if !groups.iter().any(|g| g == numerator) {
        return Err(anyhow!(
            "numerator group '{numerator}' not found; available: {}",
            groups.join(", ")
        ));
    }
    if !groups.iter().any(|g| g == denominator) {
        return Err(anyhow!(
            "denominator group '{denominator}' not found; available: {}",
            groups.join(", ")
        ));
    }
    if numerator == denominator {
        return Err(anyhow!(
            "numerator and denominator must be different groups (got '{numerator}' for both)"
        ));
    }

    // Rebuild the working intensity matrix from raw using the selected
    // normalization. table.intensity_raw is never mutated; only
    // table.intensity is overwritten.
    table.intensity = crate::normalize::apply(
        normalization,
        &table.intensity_raw,
        mapping,
        &table.sample_cols,
    )
    .map_err(|e| anyhow!("{}", e))?;

    // InChIKey deduplication. When `dedup_enabled`, compute the kept-index
    // mask + audit report once; the per-feature loop below consults the
    // mask and skips dup-losers with the same in-loop-skip pattern used by
    // `drop_unknown`. `intensity_raw` / `intensity` are NOT resized — the
    // mask only affects iteration. Per the cascade-only contract,
    // singletons and null-InChIKey features are always kept.
    let (kept_indices_opt, dedup_report_opt): (Option<HashSet<usize>>, Option<DedupReport>) =
        if dedup_enabled {
            let (kept, report) = run_dedup(&table.features);
            (Some(kept), Some(report))
        } else {
            (None, None)
        };

    // Remap mapping-indices to per-table column indices. In single-mode the
    // mapping was built from this table's `sample_cols`, so the remap is the
    // identity. In dual-mode the mapping was built from the UNION of both
    // modes' sample columns — its indices over-index this mode's table by
    // exactly the other mode's sample count. We translate mapping_idx →
    // sample_name → table_idx (and drop names absent from this table) so
    // downstream `table.intensity[[i, j]]` accesses are always in bounds.
    let remap = |group: &str| -> Vec<usize> {
        mapping
            .samples_in(group)
            .iter()
            .filter_map(|&mi| {
                let name = mapping.sample_name(mi)?;
                table.sample_cols.iter().position(|c| c == name)
            })
            .collect()
    };
    let num_idx = remap(numerator);
    let den_idx = remap(denominator);

    // Collect candidate features (post Unknown filter, post pre-filter).
    let mut features: Vec<DamFeature> = Vec::new();
    let mut p_raw: Vec<f64> = Vec::new();
    let mut skipped: usize = 0;
    // Diagnostic: count features where the parametric test sees one group
    // with zero sample variance after arcsinh+Pareto (typically because the
    // raw values are all identical in one group — feature below LOD in every
    // replicate of that condition, or perfectly saturated). Welch's
    // Satterthwaite df collapses to `n - 1` of the non-degenerate group in
    // this case, producing conservative p-values. Behaviour matches R
    // `t.test(var.equal=FALSE)` and SciPy `ttest_ind(equal_var=False)`; we
    // surface a per-run INFO log so users know to consider BM if these
    // features are biologically interesting.
    let mut zero_variance_features: usize = 0;

    let total = table.features.len();
    let report = |i: usize| {
        if let Some(tx) = progress_tx.as_ref() {
            let _ = tx.send(DamProgress {
                completed: i + 1,
                total,
            });
        }
    };

    for (i, feat) in table.features.iter().enumerate() {
        // (3a) Dedup mask skip — dedup-losers are removed BEFORE the
        // existing `drop_unknown` check and pre-filter. Dup-losers do NOT
        // contribute to `DamResult.skipped` (which is reserved for
        // pre-filter drops). Per the `dam-analysis` spec, they appear in
        // `report.dropped` exactly once and never re-enter the pipeline.
        if let Some(ref kept_set) = kept_indices_opt
            && !kept_set.contains(&i)
        {
            report(i);
            continue;
        }
        // (3b) Unknown filter — existing behaviour. Null-InChIKey features
        // are passed through by dedup and then handled here per the
        // pre-existing contract.
        if drop_unknown && feat.inchikey.is_none() {
            // Silent drop; not counted in `skipped` (which is reserved for pre-filter).
            report(i);
            continue;
        }
        let num_vals: Vec<f64> = num_idx.iter().map(|&j| table.intensity[[i, j]]).collect();
        let den_vals: Vec<f64> = den_idx.iter().map(|&j| table.intensity[[i, j]]).collect();

        if !passes_prefilter(&num_vals, &den_vals) {
            skipped += 1;
            report(i);
            continue;
        }

        let numerator_mean = nan_aware_mean(&num_vals);
        let denominator_mean = nan_aware_mean(&den_vals);
        let numerator_median = nan_aware_median(&num_vals);
        let denominator_median = nan_aware_median(&den_vals);

        let (p_value, effect_size, fc_basis, fold_change, log2_fold_change) = match method {
            DamMethod::Welch | DamMethod::Student => {
                // Parametric tests apply an optional arcsinh (project's
                // "generalised log") when `log_transform == true`. The earlier
                // hard-coded `pareto_scale_in_place` step was removed because
                // per-feature linear rescaling is bit-invariant for Welch /
                // Student t-statistics (`raw_t == arcsinh+pareto_t ==
                // arcsinh-only_t == 14.7173367156` empirically verified at
                // f64 precision on the `[10,12,11]` vs `[1,2,1.5]` fixture
                // — see add-log-transform-and-scaling design D1 + dam-analysis
                // Welch/Student spec scenarios).
                let mut all: Vec<f64> = num_vals.iter().chain(den_vals.iter()).copied().collect();
                if log_transform {
                    arcsinh_in_place(&mut all);
                }
                let (scaled_num, scaled_den) = all.split_at(num_vals.len());
                // Pre-flag zero-variance-group features (one group is constant
                // after NaN-aware filtering). arcsinh is strictly monotone, so
                // var(scaled_*) == 0 iff var(raw) == 0 when `log_transform` is
                // true; when false, the slice is raw and the same check applies.
                // `is_effectively_zero_variance` uses a relative tolerance so
                // FP-noise variance at high intensity scale (`va ≈ ε² × c²` for
                // bit-equal-but-not-exact-equal inputs) is caught — pre-
                // 2026-05-29 the diagnostic used exact `va == 0.0` and missed
                // the FP-noise cases that hit the same Welch-Satterthwaite df
                // collapse the counter was designed to surface.
                let va = nan_aware_var(scaled_num, 1);
                let vb = nan_aware_var(scaled_den, 1);
                if is_effectively_zero_variance(scaled_num, va)
                    || is_effectively_zero_variance(scaled_den, vb)
                {
                    zero_variance_features += 1;
                }
                let p = match method {
                    DamMethod::Welch => welch_t_two_tailed(scaled_num, scaled_den),
                    DamMethod::Student => student_t_two_tailed(scaled_num, scaled_den),
                    DamMethod::BrunnerMunzel => unreachable!(),
                };
                // FC scale matches the test scale.
                //
                // log_transform=false: classical raw mean ratio. Backwards-
                // compatible with the pre-2026-05-29 pipeline.
                //
                // log_transform=true: compute log2FC on the arcsinh-transformed
                // scale and back-transform `fc = 2^log2_fc`. Reason: Jensen's
                // inequality on heavy-tailed data lets the raw mean ratio
                // disagree in SIGN with the arcsinh-mean difference that the
                // t-statistic actually tests. Reporting the raw mean ratio
                // alongside an arcsinh-scale p value silently misclassifies
                // outlier-driven features (e.g. num=[0.1×9, 100], den=[5×10]
                // has raw FC=2.02 → "Up" but t-stat=−3.25 → "Down"). With
                // arcsinh-scale log2FC the sign always agrees with the t-stat.
                // For large x, arcsinh(x) ≈ ln(2x), so log2FC asymptotes to
                // log2(GM(num) / GM(den)) — the standard log-FC of limma /
                // DESeq2 et al. For small x where arcsinh(x) ≈ x, log2FC
                // degrades to (mean_num − mean_den)/ln(2), an arithmetic
                // difference in log2 units rather than a true ratio. This is
                // a known consequence of variance-stabilisation; documented in
                // USER_MANUAL.md.
                let (basis, fc, log2_fc) = if log_transform {
                    let mean_a_t = nan_aware_mean(scaled_num);
                    let mean_b_t = nan_aware_mean(scaled_den);
                    let log2_fc = (mean_a_t - mean_b_t) / std::f64::consts::LN_2;
                    (FcBasis::ArcsinhMean, log2_fc.exp2(), log2_fc)
                } else {
                    let fc = numerator_mean / denominator_mean;
                    (FcBasis::Mean, fc, fc.log2())
                };
                (p, None, basis, fc, log2_fc)
            }
            DamMethod::BrunnerMunzel => {
                let p = brunner_munzel_two_tailed(&num_vals, &den_vals);
                // δ is None when either group has < 2 non-NaN values (matches BM NaN guard).
                let na = num_vals.iter().filter(|x| !x.is_nan()).count();
                let nb = den_vals.iter().filter(|x| !x.is_nan()).count();
                let delta = if na < 2 || nb < 2 {
                    None
                } else {
                    Some(cliffs_delta(&num_vals, &den_vals))
                };
                let fc = numerator_median / denominator_median;
                (p, delta, FcBasis::Median, fc, fc.log2())
            }
        };
        p_raw.push(p_value);

        features.push(DamFeature {
            alignment_id: feat.alignment_id.clone(),
            metabolite_name: feat.metabolite_name.clone(),
            inchikey: feat.inchikey.clone(),
            average_rt_min: feat.average_rt_min,
            average_mz: feat.average_mz,
            formula: feat.formula.clone(),
            smiles: feat.smiles.clone(),
            numerator_mean,
            denominator_mean,
            numerator_median,
            denominator_median,
            fold_change,
            log2_fold_change,
            fc_basis,
            p_value,
            p_adjusted: f64::NAN, // placeholder; filled below
            neg_log10_p_adjusted: f64::NAN,
            effect_size,
        });
    }

    // FDR correction over the surviving features' p values (NaN passes through).
    // The caller chooses BH (literature default) or BY (more conservative).
    let p_adj = adjust_pvalues(&p_raw, fdr_method);
    for (feat, &q) in features.iter_mut().zip(p_adj.iter()) {
        feat.p_adjusted = q;
        // q-value of exactly 0.0 happens when the raw p underflows (very
        // large |t|, well-separated groups). Map to +INF rather than NaN so
        // (a) classify_trend reads p_adjusted=0 as "passes any positive FDR
        //     threshold" — feature is correctly classified Up/Down
        // (b) the volcano renderer can dock these points at y_max instead
        //     of silently dropping them via the NaN guard
        // Mirrors the ±INF log2_fold_change saturation handling on the X
        // axis. Reserve NaN for the genuine "p couldn't be computed" cases
        // (BM perfect stratification, Welch n<2 — those keep NaN p_value
        // upstream and propagate NaN through adjust_pvalues).
        feat.neg_log10_p_adjusted = if q.is_nan() {
            f64::NAN
        } else if q <= 0.0 {
            f64::INFINITY
        } else {
            -q.log10()
        };
    }

    // Per-run summary log for parametric features whose Welch-Satterthwaite
    // df collapsed because one group had zero variance. Quiet for normal
    // runs; surfaces a single line when these features exist so the user
    // can decide whether to switch to BM.
    if zero_variance_features > 0 {
        info!(
            zero_variance_features,
            total = features.len(),
            method = ?method,
            "{} feature(s) had zero variance in one group; Welch-Satterthwaite df collapses to n-1 of the non-degenerate group, producing conservative p-values. Consider Brunner-Munzel for these.",
            zero_variance_features
        );
    }

    debug!(
        method = ?method,
        normalization = ?normalization.method,
        fdr_method = ?fdr_method,
        kept = features.len(),
        skipped,
        drop_unknown,
        dedup_enabled,
        dedup_dropped = dedup_report_opt.as_ref().map(|r| r.dropped.len()).unwrap_or(0),
        "run_dam complete"
    );

    Ok(DamResult {
        method,
        numerator: numerator.to_string(),
        denominator: denominator.to_string(),
        features,
        skipped,
        fdr_method,
        dedup_report: dedup_report_opt,
    })
}

/// Returns true when `values`'s sample variance is so small relative to its
/// magnitude that the Welch–Satterthwaite df collapses to `n − 1` of the OTHER
/// group — the diagnostic condition `zero_variance_features` reports.
///
/// `var` is the pre-computed `nan_aware_var(values, 1)` (the caller already
/// has it). Returns `false` when `var` is NaN (`n < 2` case — handled
/// separately by the per-method NaN guards downstream).
///
/// **Why relative tolerance instead of `var == 0.0` exact equality.** When the
/// inputs are mathematically equal but arrive via FP arithmetic at intensity
/// scale `c` (e.g. normalization output where bit-equal pre-norm values get
/// divided by slightly-different per-sample factors), the realised variance is
/// `va ≈ ε² × c²` where `ε ≈ 2.22e−16` is the f64 machine epsilon. At scale
/// `c = 1e6` that's `va ≈ 5e−20` — non-zero by FP standards but still triggers
/// the Welch-Satterthwaite df collapse. An exact `va == 0.0` check would miss
/// these features; the threshold below catches them across every realistic
/// MS-DIAL intensity scale without false-positive on real biological variation
/// (the smallest coefficient of variation in metabolomics is several orders of
/// magnitude above `sqrt(1e−20) = 1e−10`).
fn is_effectively_zero_variance(values: &[f64], var: f64) -> bool {
    /// Variance is treated as effectively zero when it falls below
    /// `(max(|mean|, 1))² × FP_NOISE_RELATIVE`. `1e−20` is well above the
    /// FP-noise ceiling `ε² ≈ 5e−32` (so noise reliably triggers) and well
    /// below any plausible biological signal `var / mean²` (so genuine
    /// variation is never flagged).
    const FP_NOISE_RELATIVE: f64 = 1e-20;
    if !var.is_finite() {
        return false;
    }
    let mean = nan_aware_mean(values);
    let scale_sq = mean.abs().max(1.0).powi(2);
    var < scale_sq * FP_NOISE_RELATIVE
}

/// Classify a single feature against the user's thresholds. Pure function; no side effects.
///
/// - Welch: `Up` iff `p_adjusted < fdr AND log2_fc > log2(fc)`; `Down` iff symmetric.
/// - BM: same as Welch AND `|effect_size| ≥ delta_threshold`. `effect_size == None`
///   classifies as `NotSignificant` regardless of the other fields.
pub fn classify_trend(
    feature: &DamFeature,
    fc_threshold: f64,
    fdr_threshold: f64,
    delta_threshold: f64,
    method: DamMethod,
) -> Trend {
    if feature.p_adjusted.is_nan() || feature.log2_fold_change.is_nan() {
        return Trend::NotSignificant;
    }
    if feature.p_adjusted >= fdr_threshold {
        return Trend::NotSignificant;
    }
    let log2_fc_threshold = fc_threshold.log2();
    let crosses_fc_up = feature.log2_fold_change > log2_fc_threshold;
    let crosses_fc_down = feature.log2_fold_change < -log2_fc_threshold;
    if !crosses_fc_up && !crosses_fc_down {
        return Trend::NotSignificant;
    }
    match method {
        DamMethod::Welch | DamMethod::Student => {
            if crosses_fc_up {
                Trend::Up
            } else {
                Trend::Down
            }
        }
        DamMethod::BrunnerMunzel => match feature.effect_size {
            None => Trend::NotSignificant,
            Some(delta) if delta.is_nan() => Trend::NotSignificant,
            Some(delta) => {
                if delta.abs() < delta_threshold {
                    Trend::NotSignificant
                } else if crosses_fc_up && delta >= delta_threshold {
                    Trend::Up
                } else if crosses_fc_down && delta <= -delta_threshold {
                    Trend::Down
                } else {
                    Trend::NotSignificant
                }
            }
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Exact-zero variance (bit-equal inputs, e.g. below-LOD imputation to the
    /// same constant across replicates) MUST still trigger the diagnostic —
    /// this is the historical case the counter was designed for.
    #[test]
    fn zero_variance_detected_for_bit_equal_inputs() {
        let values = [3.0_f64, 3.0, 3.0];
        let var = nan_aware_var(&values, 1);
        assert_eq!(var, 0.0, "bit-equal inputs must produce var == 0.0 exactly");
        assert!(is_effectively_zero_variance(&values, var));
    }

    /// FP-noise variance at high intensity scale (e.g. post-normalization
    /// arithmetic on bit-equal pre-norm values) used to slip past the
    /// pre-2026-05-29 `va == 0.0` check. Pins the relative-tolerance fix.
    #[test]
    fn zero_variance_detected_for_fp_noise_at_intensity_scale_1e6() {
        // Two nominally-equal values that differ by 1 ULP near 1e6 — variance
        // of order `ε² × c²` ≈ 5e−20 at scale c=1e6, non-zero but the Welch-
        // Satterthwaite df collapses just the same.
        let c: f64 = 1.0e6;
        let values = [c, c + (1.0e-9_f64), c];
        let var = nan_aware_var(&values, 1);
        assert!(
            var > 0.0,
            "FP-noise variance must be strictly positive for this regression test; got {var}"
        );
        assert!(
            is_effectively_zero_variance(&values, var),
            "FP-noise variance ({var}) at intensity scale 1e6 must trigger the \
             diagnostic; the pre-fix exact `var == 0.0` check would have missed it"
        );
    }

    /// Genuine biological signal at coefficient-of-variation ≈ 1 % MUST NOT be
    /// flagged. CV(%) = 100 × sqrt(var)/|mean|; for [1e6, 1.01e6, 0.99e6] it's
    /// ~ 1 % — orders of magnitude above the noise floor.
    #[test]
    fn real_signal_at_cv_one_percent_not_flagged() {
        let values = [1.0e6_f64, 1.01e6, 0.99e6];
        let var = nan_aware_var(&values, 1);
        assert!(
            !is_effectively_zero_variance(&values, var),
            "real biological variation (CV ≈ 1 %, var ≈ 1e8) must not be \
             flagged as zero-variance; got var = {var}"
        );
    }

    /// `var.is_nan()` (n < 2 case) MUST return `false` from the helper — that
    /// path is handled separately by the per-method NaN guard, not by this
    /// diagnostic.
    #[test]
    fn nan_variance_returns_false() {
        let values = [3.0_f64];
        let var = nan_aware_var(&values, 1);
        assert!(var.is_nan(), "single-value var with ddof=1 must be NaN");
        assert!(
            !is_effectively_zero_variance(&values, var),
            "NaN variance must not trigger the diagnostic"
        );
    }

    fn synth_feature(p_adjusted: f64, log2_fc: f64, effect_size: Option<f64>) -> DamFeature {
        DamFeature {
            alignment_id: "x".into(),
            metabolite_name: "X".into(),
            inchikey: Some("XXX".into()),
            average_rt_min: None,
            average_mz: None,
            formula: None,
            smiles: None,
            numerator_mean: 0.0,
            denominator_mean: 0.0,
            numerator_median: 0.0,
            denominator_median: 0.0,
            fold_change: 0.0,
            log2_fold_change: log2_fc,
            fc_basis: FcBasis::Mean,
            p_value: 0.0,
            p_adjusted,
            neg_log10_p_adjusted: 0.0,
            effect_size,
        }
    }

    #[test]
    fn welch_up_when_fc_and_fdr_pass() {
        let f = synth_feature(0.01, 1.5, None);
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.33, DamMethod::Welch),
            Trend::Up
        );
    }

    #[test]
    fn welch_down() {
        let f = synth_feature(0.01, -1.5, None);
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.33, DamMethod::Welch),
            Trend::Down
        );
    }

    #[test]
    fn welch_ns_below_fc() {
        let f = synth_feature(0.01, 0.5, None);
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.33, DamMethod::Welch),
            Trend::NotSignificant
        );
    }

    #[test]
    fn bm_requires_delta_threshold() {
        let f = synth_feature(0.01, 1.5, Some(0.2));
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.33, DamMethod::BrunnerMunzel),
            Trend::NotSignificant
        );
    }

    #[test]
    fn bm_passes_with_delta_above_threshold() {
        let f = synth_feature(0.01, 1.5, Some(0.5));
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.33, DamMethod::BrunnerMunzel),
            Trend::Up
        );
    }

    #[test]
    fn bm_with_none_delta_is_ns() {
        let f = synth_feature(0.01, 1.5, None);
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.33, DamMethod::BrunnerMunzel),
            Trend::NotSignificant
        );
    }

    #[test]
    fn nan_p_is_ns() {
        let f = synth_feature(f64::NAN, 1.5, Some(0.5));
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.33, DamMethod::BrunnerMunzel),
            Trend::NotSignificant
        );
    }

    #[test]
    fn student_up_when_fc_and_fdr_pass() {
        let f = synth_feature(0.01, 1.5, None);
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.33, DamMethod::Student),
            Trend::Up
        );
    }

    #[test]
    fn student_down() {
        let f = synth_feature(0.01, -1.5, None);
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.33, DamMethod::Student),
            Trend::Down
        );
    }

    #[test]
    fn student_ns_below_fc() {
        let f = synth_feature(0.01, 0.5, None);
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.33, DamMethod::Student),
            Trend::NotSignificant
        );
    }

    #[test]
    fn student_ignores_delta_threshold() {
        // Student has no effect size, so a non-zero delta_threshold is irrelevant —
        // passes whenever Welch would pass on the same input.
        let f = synth_feature(0.01, 1.5, None);
        assert_eq!(
            classify_trend(&f, 2.0, 0.05, 0.9, DamMethod::Student),
            Trend::Up
        );
    }
}
