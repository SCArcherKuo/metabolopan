//! Volcano renderer smoke tests via the public API.

use metabolopan::dam::types::{DamFeature, DamMethod, DamResult, FcBasis};
use metabolopan::plot::{VolcanoOpts, render_volcano};

fn feat(log2_fc: f64, neg_log10: f64) -> DamFeature {
    DamFeature {
        alignment_id: "f".into(),
        metabolite_name: "f".into(),
        inchikey: Some("X".into()),
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
        p_value: 0.01,
        p_adjusted: 0.01,
        neg_log10_p_adjusted: neg_log10,
        effect_size: None,
    }
}

#[test]
fn renders_at_preview_resolution() {
    let result = DamResult {
        method: DamMethod::Welch,
        numerator: "A".into(),
        denominator: "B".into(),
        features: vec![
            feat(2.0, 4.0),
            feat(-2.0, 4.0),
            feat(0.0, 0.5),
            feat(f64::INFINITY, 2.0),
            feat(f64::NEG_INFINITY, 2.0),
        ],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        dedup_report: None,
    };
    let opts = VolcanoOpts {
        width_px: 800,
        height_px: 800,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
    };
    let buf = render_volcano(&result, &opts).expect("render");
    assert_eq!(buf.len(), 800 * 800 * 4);
}

#[test]
fn renders_empty_result_without_panicking() {
    let result = DamResult {
        method: DamMethod::Welch,
        numerator: "A".into(),
        denominator: "B".into(),
        features: vec![],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        dedup_report: None,
    };
    let opts = VolcanoOpts {
        width_px: 400,
        height_px: 300,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
    };
    let buf = render_volcano(&result, &opts).expect("empty render");
    assert_eq!(buf.len(), 400 * 300 * 4);
}

#[test]
fn renders_underflowed_q_feature_at_top_without_panic() {
    // q-underflow contract (post-2026-05-26): when BH/BY adjusted q is 0.0
    // (raw p underflowed from very large |t|), neg_log10_p_adjusted is set
    // to +INF (NOT NaN). The volcano renderer must dock these at y_max
    // instead of silently dropping them via the old NaN guard. This test
    // confirms the renderer accepts INF without panicking and the resulting
    // buffer is well-formed.
    let result = DamResult {
        method: DamMethod::Welch,
        numerator: "A".into(),
        denominator: "B".into(),
        features: vec![
            feat(3.0, f64::INFINITY),  // Up + saturated q — dock top
            feat(-3.0, f64::INFINITY), // Down + saturated q — dock top
            feat(2.5, 4.0),            // normal finite point for comparison
        ],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        dedup_report: None,
    };
    let opts = VolcanoOpts {
        width_px: 500,
        height_px: 500,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
    };
    let buf = render_volcano(&result, &opts).expect("INF y must NOT panic the renderer");
    assert_eq!(buf.len(), 500 * 500 * 4);
}

#[test]
fn renders_nan_q_feature_skipped() {
    // NaN neg_log10 is reserved for genuine "p couldn't be computed" cases
    // (BM perfect stratification, Welch n<2). These remain dropped from the
    // plot via the `continue` branch — the NaN-skip behaviour is unchanged.
    let result = DamResult {
        method: DamMethod::Welch,
        numerator: "A".into(),
        denominator: "B".into(),
        features: vec![feat(2.0, f64::NAN), feat(2.5, 4.0)],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        dedup_report: None,
    };
    let opts = VolcanoOpts {
        width_px: 500,
        height_px: 500,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
    };
    let buf = render_volcano(&result, &opts).expect("NaN feature handled by skip");
    assert_eq!(buf.len(), 500 * 500 * 4);
}

#[test]
#[ignore]
fn renders_at_600dpi_export_size() {
    // 6000x6000 = 144 MB RGBA buffer; takes ~2 s. Run with --ignored.
    let result = DamResult {
        method: DamMethod::Welch,
        numerator: "A".into(),
        denominator: "B".into(),
        features: vec![feat(2.0, 4.0)],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        dedup_report: None,
    };
    let opts = VolcanoOpts {
        width_px: 6000,
        height_px: 6000,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
    };
    let buf = render_volcano(&result, &opts).expect("600 DPI render");
    assert_eq!(buf.len(), 6000 * 6000 * 4);
}
