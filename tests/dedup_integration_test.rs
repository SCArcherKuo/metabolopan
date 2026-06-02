//! End-to-end integration test for the deduplication pre-pass — invokes
//! `run_dam` with `dedup_enabled = true` / `false` against a synthetic
//! `MetabolomicsTable` containing controlled InChIKey duplicates.
//!
//! See `openspec/changes/add-msdial-duplicate-filter/specs/dam-analysis/spec.md`
//! "System SHALL apply deduplication mask before pre-filter when
//! `dedup_enabled` is true" for the requirement we're validating.

use std::io::Write;

use metabolopan::dam::{DamConfig, DamMethod, fdr::FdrMethod, run_dam};
use metabolopan::data::types::{FeatureMeta, MetabolomicsTable};
use metabolopan::data::{GroupMapping, load_group_mapping};
use metabolopan::dedup::CascadeStep;
use metabolopan::normalize::NormalizationConfig;
use ndarray::Array2;

const TEST_FDR: FdrMethod = FdrMethod::BenjaminiYekutieli;

/// Base `DamConfig` for these dedup tests: Welch, no normalization,
/// drop-unknown OFF, dedup ON, log-transform on, pinned `TEST_FDR`. Variants
/// override one field via struct-update so the delta is visible.
fn base_cfg() -> DamConfig {
    DamConfig {
        method: DamMethod::Welch,
        normalization: NormalizationConfig::default(),
        drop_unknown: false,
        dedup_enabled: true,
        log_transform: true,
        fdr_method: TEST_FDR,
    }
}

/// Build a `FeatureMeta` with the dedup-relevant fields populated. Every
/// non-dedup field is left at a defensible default.
#[allow(clippy::too_many_arguments)]
fn feat(
    alignment_id: &str,
    inchikey: Option<&str>,
    adduct: Option<&str>,
    ms_ms: Option<bool>,
    total: Option<f64>,
    fill: Option<f64>,
    sn: Option<f64>,
) -> FeatureMeta {
    FeatureMeta {
        alignment_id: alignment_id.to_string(),
        metabolite_name: alignment_id.to_string(),
        inchikey: inchikey.map(|s| s.to_string()),
        adduct_type: adduct.map(|s| s.to_string()),
        average_rt_min: Some(1.0),
        average_mz: Some(100.0),
        formula: None,
        smiles: None,
        fill_percent: fill,
        ms_ms_matched: ms_ms,
        isotope_tracking_weight_number: Some(0),
        total_score: total,
        sn_average: sn,
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

/// Synthesise a 10-feature `MetabolomicsTable` with 3 same-InChIKey
/// duplicate groups (2 features each, 6 total dup-participants) plus 3
/// non-duplicate annotated features and 1 null-InChIKey feature.
/// Intensities are constructed so every feature individually passes the
/// pre-filter (`nunique > 1`, `IQR > 0`) — i.e. the only post-pre-filter
/// drops will be the 3 dedup-losers when `dedup_enabled = true`.
fn synth_inputs() -> (MetabolomicsTable, GroupMapping) {
    let features = vec![
        // Group X — A is the cascade winner on Total score.
        feat(
            "X_A",
            Some("XXX"),
            Some("[M+H]+"),
            Some(true),
            Some(0.9),
            Some(80.0),
            Some(50.0),
        ),
        feat(
            "X_B",
            Some("XXX"),
            Some("[M+H]+"),
            Some(true),
            Some(0.8),
            Some(70.0),
            Some(40.0),
        ),
        // Group Y — Q is the winner on Total score.
        feat(
            "Y_P",
            Some("YYY"),
            Some("[M+H]+"),
            Some(true),
            Some(0.85),
            Some(60.0),
            Some(30.0),
        ),
        feat(
            "Y_Q",
            Some("YYY"),
            Some("[M+H]+"),
            Some(true),
            Some(0.95),
            Some(75.0),
            Some(45.0),
        ),
        // Group Z — winner on adduct class (Primary beats Dimer).
        feat(
            "Z_R",
            Some("ZZZ"),
            Some("[M+H]+"),
            Some(true),
            Some(0.85),
            Some(60.0),
            Some(30.0),
        ),
        feat(
            "Z_S",
            Some("ZZZ"),
            Some("[2M+H]+"),
            Some(true),
            Some(0.85),
            Some(99.0),
            Some(99.0),
        ),
        // Three lone annotated features (each unique InChIKey).
        feat(
            "S1",
            Some("S1K"),
            Some("[M+H]+"),
            Some(true),
            Some(0.7),
            Some(50.0),
            Some(20.0),
        ),
        feat(
            "S2",
            Some("S2K"),
            Some("[M+H]+"),
            Some(true),
            Some(0.8),
            Some(55.0),
            Some(25.0),
        ),
        feat(
            "S3",
            Some("S3K"),
            Some("[M+H]+"),
            Some(true),
            Some(0.85),
            Some(65.0),
            Some(35.0),
        ),
        // One null-InChIKey feature (should pass through dedup).
        feat("N1", None, None, None, None, None, None),
    ];
    // Six samples (3 T, 3 C), each feature has well-separated intensities so
    // every row passes nunique>1 and IQR>0.
    let n_features = features.len();
    let sample_cols = vec![
        "T-1".to_string(),
        "T-2".to_string(),
        "T-3".to_string(),
        "C-1".to_string(),
        "C-2".to_string(),
        "C-3".to_string(),
    ];
    // Same intensity pattern per feature, scaled by feature index to ensure
    // distinct row signatures. Pattern: T-* ~ 100..120, C-* ~ 200..220.
    let mut flat: Vec<f64> = Vec::with_capacity(n_features * sample_cols.len());
    for (i, _) in features.iter().enumerate() {
        let base = 100.0 + (i as f64);
        flat.extend_from_slice(&[
            base,
            base + 5.0,
            base + 10.0,
            base * 2.0,
            base * 2.0 + 5.0,
            base * 2.0 + 10.0,
        ]);
    }
    let intensity = Array2::from_shape_vec((n_features, 6), flat).expect("matrix");
    let table = MetabolomicsTable {
        annotated_count: features.iter().filter(|f| f.inchikey.is_some()).count(),
        features,
        sample_cols: sample_cols.clone(),
        intensity_raw: intensity.clone(),
        intensity,
        excluded_cols: vec![],
    };
    let csv = "sample,group\nT-1,T\nT-2,T\nT-3,T\nC-1,C\nC-2,C\nC-3,C\n";
    let tmp = write_tmp_csv(csv);
    let mapping = load_group_mapping(tmp.path(), &table.sample_cols).expect("mapping");
    (table, mapping)
}

#[tokio::test]
async fn dedup_enabled_reduces_feature_count_and_attaches_report() {
    let (mut table, mapping) = synth_inputs();
    let result = run_dam(&mut table, &mapping, "T", "C", &base_cfg(), None)
        .await
        .expect("DAM run");

    let report = result
        .dedup_report
        .as_ref()
        .expect("dedup_report must be Some when dedup_enabled = true");
    // 3 dup-loser groups → 3 dropped
    assert_eq!(report.dropped.len(), 3, "expected 3 dup-losers");
    // 6 dup-participants → 3 winners + 3 annotated singletons + 1 null = 7 kept
    assert_eq!(report.kept_count, 7);
    assert_eq!(report.null_inchikey_passthrough, 1);
    // Post-dedup features.len() == 7 (no further pre-filter drops by design)
    assert_eq!(result.features.len(), 7);
    assert_eq!(result.skipped, 0, "no pre-filter drops in synthetic input");

    // Spot-check one decision: X_A beats X_B on TotalScore.
    let x_loser = report
        .dropped
        .iter()
        .find(|d| d.inchikey == "XXX")
        .expect("X group should have a dropped feature");
    assert_eq!(x_loser.alignment_id, "X_B");
    assert_eq!(x_loser.winner_alignment_id, "X_A");
    assert_eq!(
        x_loser.decided_at,
        CascadeStep::TotalScore,
        "Group X must decide at TotalScore (ms_ms tied, total_score differs); \
         falling to Tiebreak would mean the cascade is not reaching level 1b"
    );
}

#[tokio::test]
async fn dedup_disabled_leaves_features_untouched_and_report_is_none() {
    let (mut table, mapping) = synth_inputs();
    let baseline = run_dam(
        &mut table,
        &mapping,
        "T",
        "C",
        &DamConfig {
            dedup_enabled: false,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("DAM run");

    assert!(
        baseline.dedup_report.is_none(),
        "dedup_report MUST be None when dedup_enabled = false"
    );
    // All 10 features survive (drop_unknown=false keeps the null-InChIKey
    // one too; the synthetic input is designed to pass pre-filter).
    assert_eq!(baseline.features.len(), 10);
    assert_eq!(baseline.skipped, 0);
}

#[tokio::test]
async fn dedup_independent_of_method() {
    // The same input run under Welch / Student / BM with dedup_enabled
    // produces identical DroppedFeature lists (dedup is method-agnostic).
    let methods = [
        DamMethod::Welch,
        DamMethod::Student,
        DamMethod::BrunnerMunzel,
    ];
    let mut prev_dropped: Option<Vec<(String, String)>> = None;
    for method in methods {
        let (mut table, mapping) = synth_inputs();
        let result = run_dam(
            &mut table,
            &mapping,
            "T",
            "C",
            &DamConfig {
                method,
                ..base_cfg()
            },
            None,
        )
        .await
        .expect("DAM run");
        let report = result.dedup_report.expect("dedup_report Some");
        let dropped: Vec<(String, String)> = report
            .dropped
            .iter()
            .map(|d| (d.alignment_id.clone(), d.winner_alignment_id.clone()))
            .collect();
        if let Some(prev) = &prev_dropped {
            assert_eq!(
                &dropped, prev,
                "dedup decision MUST be method-agnostic (method = {method:?})"
            );
        }
        prev_dropped = Some(dropped);
    }
}

/// Dual-mode dedup: parse real POS + NEG fixture files, run `run_dam` with
/// `dedup_enabled = true` on each, and assert both modes produce a populated
/// `DedupReport`. This guards that dedup runs on each ion-mode table
/// independently and that the dedup path is exercised on real MS-DIAL v4 data.
///
/// Skips gracefully when the `data/double-mode/` fixtures are absent (e.g.
/// CI without large data fixtures). Run locally with the full data directory.
#[tokio::test]
async fn dual_mode_dedup_runs_on_real_fixtures() {
    let pos_path = std::path::Path::new("data/double-mode/data-positive.txt");
    let neg_path = std::path::Path::new("data/double-mode/data-negative.txt");
    let meta_path = std::path::Path::new("data/double-mode/metadata.csv");
    if !pos_path.exists() || !neg_path.exists() || !meta_path.exists() {
        eprintln!("skipping dual_mode_dedup_runs_on_real_fixtures: double-mode fixtures absent");
        return;
    }

    use metabolopan::data::{load_group_mapping, parse_msdial_txt};

    let pos_table_full = parse_msdial_txt(pos_path).expect("parse POS fixture");
    let neg_table_full = parse_msdial_txt(neg_path).expect("parse NEG fixture");
    // Build the mapping from the union of both modes' sample columns, matching
    // how start_dam works in the UI (the mapping covers the full dual-mode axis).
    let mut all_sample_cols = pos_table_full.sample_cols.clone();
    for col in &neg_table_full.sample_cols {
        if !all_sample_cols.contains(col) {
            all_sample_cols.push(col.clone());
        }
    }
    let mapping_full = load_group_mapping(meta_path, &all_sample_cols).expect("load mapping");
    let mapping = mapping_full.without_unassigned_samples();
    let mut pos_table = pos_table_full.without_unassigned_samples(&mapping_full);
    let mut neg_table = neg_table_full.without_unassigned_samples(&mapping_full);

    // The double-mode fixture uses Treatment / Control groups.
    let pos_result = run_dam(
        &mut pos_table,
        &mapping,
        "Treatment",
        "Control",
        &DamConfig {
            drop_unknown: true,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("POS DAM run");

    let neg_result = run_dam(
        &mut neg_table,
        &mapping,
        "Treatment",
        "Control",
        &DamConfig {
            drop_unknown: true,
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("NEG DAM run");

    // Both modes must produce a dedup report when dedup_enabled = true.
    let pos_report = pos_result
        .dedup_report
        .expect("POS dedup_report must be Some");
    let neg_report = neg_result
        .dedup_report
        .expect("NEG dedup_report must be Some");

    // kept_count + dropped.len() == features seen by run_dedup (≥ 1 each mode).
    assert!(
        pos_report.kept_count >= 1,
        "POS dedup must keep at least one feature"
    );
    assert!(
        neg_report.kept_count >= 1,
        "NEG dedup must keep at least one feature"
    );

    // The real fixtures have known InChIKey duplicates; at least one should drop.
    // (This is a real-data assertion — if it ever fails, the fixture may have
    //  changed or the dedup is silently not running.)
    assert!(
        pos_report.dropped.len() + neg_report.dropped.len() > 0,
        "real dual-mode fixtures should contain at least one InChIKey duplicate across POS+NEG"
    );
}
