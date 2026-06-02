//! End-to-end DAM pipeline test against `tests/fixtures/msdial_mini.txt`.

use std::io::Write;

use metabolopan::dam::{
    DamConfig, DamMethod, FcBasis, fdr::FdrMethod, run::classify_trend, run_dam, types::Trend,
};
use metabolopan::data::{FeatureMeta, MetabolomicsTable, load_group_mapping, parse_msdial_txt};
use metabolopan::normalize::{NormalizationConfig, NormalizationMethod};
use ndarray::Array2;

/// Default FDR method used by the existing pipeline tests. Pinned to BY so
/// these tests preserve the pre-change bit-for-bit assertions; new BH-vs-BY
/// regression coverage is added at the end of this file via
/// `bh_vs_by_p_value_invariant_and_padjusted_differs`.
const TEST_FDR: FdrMethod = FdrMethod::BenjaminiYekutieli;

/// Base `DamConfig` for these tests: Welch, no normalization, drop-unknown on,
/// dedup off, log-transform on, and the pinned `TEST_FDR` (BY). Pairwise tests
/// derive a variant via struct-update (`DamConfig { fdr_method: …, ..base_cfg() }`)
/// so the single differing field is the visible delta.
fn base_cfg() -> DamConfig {
    DamConfig {
        method: DamMethod::Welch,
        normalization: NormalizationConfig::default(),
        drop_unknown: true,
        dedup_enabled: false,
        log_transform: true,
        fdr_method: TEST_FDR,
    }
}

fn write_tmp_csv(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new()
        .suffix(".csv")
        .tempfile()
        .expect("tempfile");
    f.write_all(content.as_bytes()).expect("write tempfile");
    f
}

fn load_mini_inputs() -> (
    metabolopan::data::MetabolomicsTable,
    metabolopan::data::GroupMapping,
) {
    let table = parse_msdial_txt(std::path::Path::new("tests/fixtures/msdial_mini.txt"))
        .expect("parse mini fixture");
    // Mini fixture has 6 Sample columns: T-1, T-2, T-3, C-1, C-2, C-3.
    let csv = "sample,group\nT-1,T\nT-2,T\nT-3,T\nC-1,C\nC-2,C\nC-3,C\nBk-1,Bk\n";
    let tmp = write_tmp_csv(csv);
    let mapping = load_group_mapping(tmp.path(), &table.sample_cols).expect("mapping");
    (table, mapping)
}

#[tokio::test]
async fn welch_pipeline_drops_unknown_by_default() {
    let (mut table, mapping) = load_mini_inputs();
    let total = table.features.len();
    let unknown = table
        .features
        .iter()
        .filter(|f| f.inchikey.is_none())
        .count();

    let result = run_dam(&mut table, &mapping, "T", "C", &base_cfg(), None)
        .await
        .expect("DAM run");
    assert_eq!(result.method, DamMethod::Welch);
    for f in &result.features {
        assert!(f.inchikey.is_some(), "Unknown feature leaked through");
        // Under log_transform=true the FC moves to the arcsinh scale so its
        // sign agrees with the t-statistic — see the heavy-tail sign-
        // disagreement test below for the motivation.
        assert_eq!(f.fc_basis, FcBasis::ArcsinhMean);
        assert!(f.effect_size.is_none(), "Welch must not produce δ");
    }
    let kept = result.features.len();
    assert!(
        kept <= total - unknown,
        "kept ({kept}) must be ≤ annotated ({})",
        total - unknown
    );
}

#[tokio::test]
async fn bm_pipeline_keeps_unknown_when_flag_off() {
    let (mut table, mapping) = load_mini_inputs();
    let result_with = run_dam(
        &mut table,
        &mapping,
        "T",
        "C",
        &DamConfig {
            method: DamMethod::BrunnerMunzel,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("DAM with drop_unknown=true");
    let result_without = run_dam(
        &mut table,
        &mapping,
        "T",
        "C",
        &DamConfig {
            method: DamMethod::BrunnerMunzel,
            drop_unknown: false,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("DAM with drop_unknown=false");

    assert!(
        result_without.features.len() >= result_with.features.len(),
        "drop_unknown=false should keep at least as many features"
    );

    for f in &result_with.features {
        assert_eq!(f.fc_basis, FcBasis::Median);
    }
}

#[tokio::test]
async fn unknown_group_returns_err() {
    let (mut table, mapping) = load_mini_inputs();
    let err = run_dam(&mut table, &mapping, "not_a_group", "C", &base_cfg(), None)
        .await
        .expect_err("unknown group must err");
    let msg = format!("{err}");
    assert!(
        msg.contains("not_a_group"),
        "error must name the missing group: {msg}"
    );
}

#[tokio::test]
async fn default_normalization_is_bit_equal_passthrough() {
    // Two DAM runs over the same inputs, both with NormalizationConfig::default()
    // (i.e. NormalizationMethod::None), must produce identical p-values,
    // log2 fold changes, and effect sizes — this guards the "None = prior
    // behaviour" promise documented on the dam-analysis spec.
    let (mut t1, mapping) = load_mini_inputs();
    let (mut t2, _) = load_mini_inputs();
    let r1 = run_dam(&mut t1, &mapping, "T", "C", &base_cfg(), None)
        .await
        .expect("run 1");
    let r2 = run_dam(&mut t2, &mapping, "T", "C", &base_cfg(), None)
        .await
        .expect("run 2");
    assert_eq!(r1.features.len(), r2.features.len());
    for (a, b) in r1.features.iter().zip(r2.features.iter()) {
        assert_eq!(a.alignment_id, b.alignment_id);
        // Comparison handles NaN (both), ±inf (both equal), and finite (within 1e-12).
        for (x, y) in [
            (a.p_value, b.p_value),
            (a.log2_fold_change, b.log2_fold_change),
            (a.p_adjusted, b.p_adjusted),
        ] {
            if x.is_nan() {
                assert!(y.is_nan(), "NaN mismatch for {}", a.alignment_id);
            } else if x.is_infinite() {
                assert!(
                    y.is_infinite() && x.signum() == y.signum(),
                    "inf mismatch for {}: {x} vs {y}",
                    a.alignment_id
                );
            } else {
                assert!(
                    (x - y).abs() < 1e-12,
                    "{} mismatch: {x} vs {y}",
                    a.alignment_id
                );
            }
        }
    }
}

#[tokio::test]
async fn sum_normalization_moves_results_vs_none() {
    // Sum normalization on the mini fixture (which has non-uniform sample
    // sums by construction) must produce at least one feature whose p-value
    // or log2FC differs from the None run.
    let (mut t_none, mapping) = load_mini_inputs();
    let (mut t_sum, _) = load_mini_inputs();
    let r_none = run_dam(&mut t_none, &mapping, "T", "C", &base_cfg(), None)
        .await
        .expect("None run");
    let r_sum = run_dam(
        &mut t_sum,
        &mapping,
        "T",
        "C",
        &DamConfig {
            normalization: NormalizationConfig {
                method: NormalizationMethod::Sum,
            },
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("Sum run");
    assert_eq!(r_none.features.len(), r_sum.features.len());
    let any_diff = r_none
        .features
        .iter()
        .zip(r_sum.features.iter())
        .any(|(a, b)| {
            let p_diff = match (a.p_value.is_nan(), b.p_value.is_nan()) {
                (true, true) => false,
                (true, false) | (false, true) => true,
                (false, false) => (a.p_value - b.p_value).abs() > 1e-9,
            };
            let fc_diff = (a.log2_fold_change - b.log2_fold_change).abs() > 1e-9;
            p_diff || fc_diff
        });
    assert!(
        any_diff,
        "Sum normalization should move at least one feature's p or log2FC vs None"
    );
}

#[tokio::test]
async fn classify_trend_round_trip_with_dam_result() {
    let (mut table, mapping) = load_mini_inputs();
    let result = run_dam(&mut table, &mapping, "T", "C", &base_cfg(), None)
        .await
        .expect("DAM run");
    // Just confirm classify_trend never panics across the whole result with default
    // thresholds.
    for feat in &result.features {
        let _ = classify_trend(feat, 2.0, 0.05, 0.33, DamMethod::Welch);
    }
    // At least one feature should classify as something other than ns (since the
    // mini fixture has very well-separated T and C columns by design).
    let any_dam = result
        .features
        .iter()
        .any(|f| classify_trend(f, 1.5, 0.5, 0.33, DamMethod::Welch) != Trend::NotSignificant);
    assert!(
        any_dam,
        "expected at least one DAM feature in the mini fixture"
    );
}

#[tokio::test]
async fn bh_vs_by_p_value_invariant_and_padjusted_differs() {
    // Regression guard: switching only `fdr_method` between BH and BY must
    // (a) keep every per-feature p_value byte-identical (the statistical test
    // is unchanged) and (b) shift at least one p_adjusted by > 1e-12 (BY's
    // c(m) harmonic factor inflates q-values relative to BH). Guards against
    // accidental dispatcher aliasing at the run_dam level.
    let (mut t_bh, mapping) = load_mini_inputs();
    let (mut t_by, _) = load_mini_inputs();
    let r_bh = run_dam(
        &mut t_bh,
        &mapping,
        "T",
        "C",
        &DamConfig {
            fdr_method: FdrMethod::BenjaminiHochberg,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("BH run");
    let r_by = run_dam(
        &mut t_by,
        &mapping,
        "T",
        "C",
        &DamConfig {
            fdr_method: FdrMethod::BenjaminiYekutieli,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("BY run");

    assert_eq!(r_bh.fdr_method, FdrMethod::BenjaminiHochberg);
    assert_eq!(r_by.fdr_method, FdrMethod::BenjaminiYekutieli);
    assert_eq!(r_bh.features.len(), r_by.features.len());

    // (a) Every p_value byte-identical (NaN-aware: both NaN counts as equal).
    for (a, b) in r_bh.features.iter().zip(r_by.features.iter()) {
        assert_eq!(a.alignment_id, b.alignment_id);
        if a.p_value.is_nan() {
            assert!(
                b.p_value.is_nan(),
                "NaN p_value mismatch for {}",
                a.alignment_id
            );
        } else {
            assert!(
                (a.p_value - b.p_value).abs() < 1e-15,
                "p_value differs for {}: {} vs {}",
                a.alignment_id,
                a.p_value,
                b.p_value
            );
        }
    }

    // (b) At least one p_adjusted shifts by > 1e-12.
    let any_padjusted_differs =
        r_bh.features
            .iter()
            .zip(r_by.features.iter())
            .any(
                |(a, b)| match (a.p_adjusted.is_nan(), b.p_adjusted.is_nan()) {
                    (true, true) => false,
                    (true, false) | (false, true) => true,
                    (false, false) => (a.p_adjusted - b.p_adjusted).abs() > 1e-12,
                },
            );
    assert!(
        any_padjusted_differs,
        "BH and BY must produce different p_adjusted on at least one feature"
    );
}

/// Sync-report coverage gap: dam-analysis "Welch with log_transform=false
/// skips arcsinh" scenario asserts the toggle produces an OBSERVABLY
/// DIFFERENT p-value. Pin the difference by running the same Welch
/// pipeline twice on the same fixture and asserting at least one feature's
/// p_value moves by more than 1e-9 across the toggle.
#[tokio::test]
async fn welch_with_log_off_differs_from_log_on() {
    let (mut t_log_on, mapping) = load_mini_inputs();
    let (mut t_log_off, _) = load_mini_inputs();
    let r_on = run_dam(&mut t_log_on, &mapping, "T", "C", &base_cfg(), None)
        .await
        .expect("log=on Welch");
    let r_off = run_dam(
        &mut t_log_off,
        &mapping,
        "T",
        "C",
        &DamConfig {
            log_transform: false,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("log=off Welch");
    assert_eq!(
        r_on.features.len(),
        r_off.features.len(),
        "feature set should be identical across the toggle (same filters)"
    );
    let any_p_differs = r_on
        .features
        .iter()
        .zip(r_off.features.iter())
        .any(
            |(a, b)| match (a.p_value.is_finite(), b.p_value.is_finite()) {
                (true, true) => (a.p_value - b.p_value).abs() > 1e-9,
                _ => false,
            },
        );
    assert!(
        any_p_differs,
        "Welch with log_transform=false MUST produce an observably different p-value \
         on at least one feature than log_transform=true (arcsinh is non-linear and \
         changes the t-statistic)"
    );
}

/// Sync-report coverage gap: dam-analysis "BM bypasses arcsinh regardless of
/// log_transform value" scenario asserts BM p-value + effect_size are
/// bit-equal across the toggle. Pin the property by running BM twice with
/// the toggle flipped and asserting every per-feature value matches.
#[tokio::test]
async fn bm_is_invariant_under_log_transform_toggle() {
    let (mut t_log_on, mapping) = load_mini_inputs();
    let (mut t_log_off, _) = load_mini_inputs();
    let r_on = run_dam(
        &mut t_log_on,
        &mapping,
        "T",
        "C",
        &DamConfig {
            method: DamMethod::BrunnerMunzel,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("BM log=on");
    let r_off = run_dam(
        &mut t_log_off,
        &mapping,
        "T",
        "C",
        &DamConfig {
            method: DamMethod::BrunnerMunzel,
            log_transform: false,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("BM log=off");
    assert_eq!(r_on.features.len(), r_off.features.len());
    for (a, b) in r_on.features.iter().zip(r_off.features.iter()) {
        // Handle NaN equality explicitly — NaN != NaN by f64 contract, but
        // for BM bit-invariance both sides must produce the same NaN-ness.
        match (a.p_value.is_finite(), b.p_value.is_finite()) {
            (true, true) => assert_eq!(
                a.p_value, b.p_value,
                "BM p_value differs across log_transform toggle on feature {}",
                a.alignment_id
            ),
            (false, false) => {} // both NaN — OK
            _ => panic!(
                "BM p_value NaN-ness differs across log_transform toggle on feature {}: \
                 on={}, off={}",
                a.alignment_id, a.p_value, b.p_value
            ),
        }
        assert_eq!(
            a.effect_size, b.effect_size,
            "BM effect_size differs across log_transform toggle on feature {}",
            a.alignment_id
        );
    }
}

/// Heavy-tail Jensen sign-disagreement fixture: raw mean ratio says "Up",
/// arcsinh-scale t-statistic says "Down". This regression test pins the
/// fix from 2026-05-29 — under `log_transform=true` the reported log2_fc
/// MUST match the t-stat sign, so a feature significant in the negative
/// direction (in the arcsinh scale that the t-test actually uses) is
/// classified `Down`, not silently flipped to `Up` by the raw mean ratio.
///
/// Fixture:
///   num = [0.1]×9 + [100.0]   (raw mean = 10.09)
///   den = [4.9, 5.0, 5.1] cycled to 10   (raw mean ≈ 4.99)
/// Raw FC = 10.09 / 4.99 ≈ 2.02 → log2 ≈ +1.01 (looks like UP)
/// arcsinh means: ≈ 0.62 vs ≈ 2.31 → diff ≈ −1.69 / ln2 ≈ −2.44 (DOWN)
/// Welch t at df ≈ 9: ≈ −3.25, two-tailed p ≈ 0.01 (significant)
#[tokio::test]
async fn welch_log_on_arcsinh_log2fc_sign_matches_t_stat_on_heavy_tail() {
    let sample_cols: Vec<String> = (0..20).map(|i| format!("S{i:02}")).collect();
    // One synthetic feature.
    let features = vec![FeatureMeta {
        alignment_id: "F001".into(),
        metabolite_name: "Heavy-tail synthetic".into(),
        inchikey: Some("ZZZZ-HEAVY".into()),
        adduct_type: None,
        average_rt_min: None,
        average_mz: None,
        formula: None,
        smiles: None,
        fill_percent: None,
        ms_ms_matched: None,
        isotope_tracking_weight_number: None,
        total_score: None,
        sn_average: None,
    }];
    // Layout: S00..S09 = num group T (one outlier), S10..S19 = den group C.
    let row: Vec<f64> = vec![
        0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 100.0, // num
        4.9, 5.0, 5.1, 4.9, 5.0, 5.1, 4.9, 5.0, 5.1, 4.9, // den
    ];
    let intensity = Array2::from_shape_vec((1, 20), row.clone()).expect("shape");
    let table = MetabolomicsTable {
        annotated_count: features.iter().filter(|f| f.inchikey.is_some()).count(),
        features,
        sample_cols: sample_cols.clone(),
        intensity_raw: intensity.clone(),
        intensity,
        excluded_cols: vec![],
    };

    // Mapping: S00..S09 → T, S10..S19 → C.
    let mut csv = String::from("sample,group\n");
    for s in sample_cols.iter().take(10) {
        csv.push_str(&format!("{s},T\n"));
    }
    for s in sample_cols.iter().skip(10) {
        csv.push_str(&format!("{s},C\n"));
    }
    let tmp = write_tmp_csv(&csv);
    let mapping = load_group_mapping(tmp.path(), &sample_cols).expect("mapping");

    let mut t_on = MetabolomicsTable {
        annotated_count: table.annotated_count,
        features: table.features.clone(),
        sample_cols: table.sample_cols.clone(),
        intensity_raw: table.intensity_raw.clone(),
        intensity: table.intensity.clone(),
        excluded_cols: vec![],
    };
    let r_on = run_dam(&mut t_on, &mapping, "T", "C", &base_cfg(), None)
        .await
        .expect("Welch log=on");

    assert_eq!(r_on.features.len(), 1, "fixture must survive prefilter");
    let f = &r_on.features[0];
    assert_eq!(f.fc_basis, FcBasis::ArcsinhMean);
    assert!(
        f.p_value.is_finite(),
        "Welch p must be finite, got {}",
        f.p_value
    );
    assert!(
        f.p_value < 0.05,
        "Welch should detect a difference on the transformed scale; got p={}",
        f.p_value
    );
    // The key assertion: arcsinh-scale FC sign MATCHES the t-stat sign.
    // Raw means say UP (+log2FC); arcsinh means say DOWN. Post-fix log2_fc
    // must be NEGATIVE.
    assert!(
        f.log2_fold_change < 0.0,
        "post-fix arcsinh-scale log2_fc must be negative (matches t-stat \
         direction); got {}",
        f.log2_fold_change
    );
    // Sanity: 2^log2_fc == fold_change (back-transform invariant).
    assert!(
        (f.fold_change - f.log2_fold_change.exp2()).abs() < 1e-9,
        "fold_change must equal 2^log2_fold_change; got fc={}, 2^log2={}",
        f.fold_change,
        f.log2_fold_change.exp2()
    );
}

/// `log_transform=false` preserves the raw-mean FC behaviour bit-for-bit:
/// FcBasis::Mean, log2_fc = log2(mean_num / mean_den), fc = mean_num/mean_den.
/// Pins the backward-compatibility guarantee for callers that pin
/// log_transform=false (e.g. published analyses, reproducibility tests).
#[tokio::test]
async fn welch_log_off_uses_raw_mean_fc_backward_compat() {
    let sample_cols: Vec<String> = (0..20).map(|i| format!("S{i:02}")).collect();
    let features = vec![FeatureMeta {
        alignment_id: "F001".into(),
        metabolite_name: "Heavy-tail synthetic".into(),
        inchikey: Some("ZZZZ-HEAVY".into()),
        adduct_type: None,
        average_rt_min: None,
        average_mz: None,
        formula: None,
        smiles: None,
        fill_percent: None,
        ms_ms_matched: None,
        isotope_tracking_weight_number: None,
        total_score: None,
        sn_average: None,
    }];
    // Same fixture as the log_transform=ON case so the contrast is sharp.
    let row: Vec<f64> = vec![
        0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 100.0, 4.9, 5.0, 5.1, 4.9, 5.0, 5.1, 4.9, 5.0,
        5.1, 4.9,
    ];
    let intensity = Array2::from_shape_vec((1, 20), row).expect("shape");
    let mut table = MetabolomicsTable {
        annotated_count: features.iter().filter(|f| f.inchikey.is_some()).count(),
        features,
        sample_cols: sample_cols.clone(),
        intensity_raw: intensity.clone(),
        intensity,
        excluded_cols: vec![],
    };
    let mut csv = String::from("sample,group\n");
    for s in sample_cols.iter().take(10) {
        csv.push_str(&format!("{s},T\n"));
    }
    for s in sample_cols.iter().skip(10) {
        csv.push_str(&format!("{s},C\n"));
    }
    let tmp = write_tmp_csv(&csv);
    let mapping = load_group_mapping(tmp.path(), &sample_cols).expect("mapping");

    let r_off = run_dam(
        &mut table,
        &mapping,
        "T",
        "C",
        &DamConfig {
            log_transform: false,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("Welch log=off");

    let f = &r_off.features[0];
    assert_eq!(f.fc_basis, FcBasis::Mean);
    // Raw mean(num) = (9 * 0.1 + 100) / 10 = 10.09
    // Raw mean(den) ≈ 4.99
    // Raw FC ≈ 2.022 → log2 ≈ +1.016 (POSITIVE = "Up")
    let expected_fc = 10.09 / 4.99;
    assert!(
        (f.fold_change - expected_fc).abs() < 0.01,
        "raw FC: expected ≈ {expected_fc}, got {}",
        f.fold_change
    );
    assert!(
        f.log2_fold_change > 0.0,
        "raw-mean log2_fc must be positive on this heavy-tail fixture (the \
         outlier pulls the numerator mean above the denominator); got {}",
        f.log2_fold_change
    );
    // The whole point of this test: log_transform=false keeps the
    // pre-2026-05-29 behaviour where log2_fc reflects the RAW mean ratio
    // (which on heavy-tail data can disagree with the t-statistic sign —
    // this is the trap that log_transform=true resolves under the new FC
    // computation, see `welch_log_on_arcsinh_log2fc_sign_matches_t_stat_on_heavy_tail`).
}

/// `introduce-dam-config-struct` guard: a fully-specified `DamConfig` (all six
/// fields set, non-default where it matters) drives `run_dam` to a
/// deterministic, well-formed result, and the config round-trips onto the
/// `DamResult`. The whole pre-existing suite already re-runs every prior
/// assertion *through* `DamConfig` (every call site was migrated), so this is
/// the focused pin for the spec's "bit-identical under DamConfig" scenario:
/// two runs of the same explicit config over the same input agree
/// field-for-field.
#[tokio::test]
async fn dam_config_signature_is_bit_identical() {
    let (mut t1, mapping) = load_mini_inputs();
    let (mut t2, _) = load_mini_inputs();
    let config = DamConfig {
        method: DamMethod::BrunnerMunzel,
        normalization: NormalizationConfig::default(),
        drop_unknown: true,
        dedup_enabled: true,
        log_transform: false,
        fdr_method: FdrMethod::BenjaminiHochberg,
    };
    let r1 = run_dam(&mut t1, &mapping, "T", "C", &config, None)
        .await
        .expect("run 1");
    let r2 = run_dam(&mut t2, &mapping, "T", "C", &config, None)
        .await
        .expect("run 2");

    // The config round-trips onto the result.
    assert_eq!(r1.method, DamMethod::BrunnerMunzel);
    assert_eq!(r1.fdr_method, FdrMethod::BenjaminiHochberg);
    assert!(
        r1.dedup_report.is_some(),
        "dedup_enabled = true must attach a DedupReport"
    );

    // Two runs of the same explicit config agree field-for-field.
    assert_eq!(r1.features.len(), r2.features.len());
    for (a, b) in r1.features.iter().zip(r2.features.iter()) {
        assert_eq!(a.alignment_id, b.alignment_id);
        assert_eq!(a.fc_basis, b.fc_basis, "fc_basis for {}", a.alignment_id);
        for (x, y) in [
            (a.p_value, b.p_value),
            (a.p_adjusted, b.p_adjusted),
            (a.log2_fold_change, b.log2_fold_change),
        ] {
            if x.is_nan() {
                assert!(y.is_nan(), "NaN mismatch for {}", a.alignment_id);
            } else {
                assert_eq!(x, y, "value mismatch for {}", a.alignment_id);
            }
        }
    }
}
