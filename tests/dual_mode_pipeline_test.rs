//! End-to-end dual-mode Stage 3 pipeline test: two DamResults → unioned
//! PubChem (mocked) → KEGG conv (mocked) → conflict-only-strict K +
//! DualModeBreakdown + ORA.
//!
//! Complements `tests/stage3_pipeline_test.rs` (single-mode) and the
//! in-module unit tests in `src/stage3/mod.rs` that exercise the dual-mode
//! math on synthetic fixtures without any network. This test wires the
//! full async orchestrator end-to-end with wiremock'd HTTP.

use metabolopan::dam::types::{DamFeature, FcBasis};
use metabolopan::dam::{DamMethod, DamResult};
use metabolopan::enrichment::EnrichmentDirection;
use metabolopan::kegg::{KeggClient, KeggCompoundSet, SpeciesKegg};
use metabolopan::pubchem::PubchemClient;
use metabolopan::stage3::{Stage3Params, run_stage3};
use reqwest::Url;
use std::sync::OnceLock;
use std::sync::mpsc;
use tempfile::tempdir;
use wiremock::matchers::{method, path, path_regex};
use wiremock::{Mock, MockServer, ResponseTemplate};

async fn serial() -> tokio::sync::MutexGuard<'static, ()> {
    static M: OnceLock<tokio::sync::Mutex<()>> = OnceLock::new();
    M.get_or_init(|| tokio::sync::Mutex::new(())).lock().await
}

fn setup_tmp_root() -> (tempfile::TempDir, tempfile::TempDir) {
    let dp = tempdir().expect("pubchem tempdir");
    let dk = tempdir().expect("kegg tempdir");
    unsafe {
        std::env::set_var("PUBCHEM_CACHE_DIR", dp.path());
        std::env::set_var("KEGG_CACHE_DIR", dk.path());
    }
    (dp, dk)
}

fn feat(inchikey: &str, log2_fc: f64, p_adj: f64) -> DamFeature {
    DamFeature {
        alignment_id: format!("aid-{inchikey}-{log2_fc}"),
        metabolite_name: format!("met-{inchikey}"),
        inchikey: Some(inchikey.to_string()),
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
        p_value: p_adj,
        p_adjusted: p_adj,
        neg_log10_p_adjusted: -p_adj.log10(),
        effect_size: None,
    }
}

/// End-to-end dual-mode run that exercises:
/// - PubChem call sees the UNION of both modes' InChIKeys (one call covers both).
/// - KEGG /conv resolves the unioned CIDs.
/// - DualModeBreakdown partition counts come out correctly.
/// - K contains POS-only, NEG-only, and agree-both cpds; excludes the
///   conflict cpd.
/// - ORA runs on the unioned universe.
#[tokio::test]
async fn end_to_end_dual_mode_pathway_runs_clean() {
    let _g = serial().await;
    let (_dp, _dk) = setup_tmp_root();

    // ── Inchikey → CID layout (PubChem mock) ──
    // K_AGREE → 100 → C_AGREE     (both modes Up → agree_both)
    // K_POSONLY → 101 → C_POSONLY (POS Up, NEG NS / not present)
    // K_NEGONLY → 102 → C_NEGONLY (POS not present, NEG Up)
    // K_CONFL → 103 → C_CONFL     (POS Up, NEG Down → excluded_by_conflict)
    let pubchem_server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path("/compound/inchikey/property/InChIKey/CSV"))
        .respond_with(ResponseTemplate::new(200).set_body_string(
            "\"InChIKey\",\"CID\"\n\
             \"K_AGREE\",\"100\"\n\
             \"K_POSONLY\",\"101\"\n\
             \"K_NEGONLY\",\"102\"\n\
             \"K_CONFL\",\"103\"\n",
        ))
        .mount(&pubchem_server)
        .await;

    let kegg_server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path_regex(r"^/conv/compound/.+"))
        .respond_with(ResponseTemplate::new(200).set_body_string(
            "pubchem:100\tcpd:C_AGREE\n\
             pubchem:101\tcpd:C_POSONLY\n\
             pubchem:102\tcpd:C_NEGONLY\n\
             pubchem:103\tcpd:C_CONFL\n",
        ))
        .mount(&kegg_server)
        .await;

    // POS DamResult: K_AGREE Up, K_POSONLY Up, K_CONFL Up.
    let pos = DamResult {
        method: DamMethod::Welch,
        numerator: "T".into(),
        denominator: "C".into(),
        features: vec![
            feat("K_AGREE", 2.0, 0.001),
            feat("K_POSONLY", 2.5, 0.001),
            feat("K_CONFL", 2.0, 0.001),
        ],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiHochberg,
        dedup_report: None,
    };
    // NEG DamResult: K_AGREE Up, K_NEGONLY Up, K_CONFL Down.
    let neg = DamResult {
        method: DamMethod::Welch,
        numerator: "T".into(),
        denominator: "C".into(),
        features: vec![
            feat("K_AGREE", 2.0, 0.001),
            feat("K_NEGONLY", 2.0, 0.001),
            feat("K_CONFL", -2.0, 0.001),
        ],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiHochberg,
        dedup_report: None,
    };

    // Pathway entries: one pathway covering all four cpds (it will be enriched).
    let species_kegg = SpeciesKegg {
        code: "syn".into(),
        fetched_at: chrono::Utc::now(),
        pathways: vec![KeggCompoundSet {
            id: "syn00001".into(),
            name: "Test enriched".into(),
            compounds: vec![
                "C_AGREE".into(),
                "C_POSONLY".into(),
                "C_NEGONLY".into(),
                "C_CONFL".into(),
                "C_FILLER".into(),
            ],
        }],
    };

    let pubchem = PubchemClient::with_base_url(pubchem_server.uri());
    let kegg = KeggClient::with_base_url(Url::parse(&kegg_server.uri()).unwrap());
    let params = Stage3Params {
        method: DamMethod::Welch,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        direction: EnrichmentDirection::Up,
        min_hit_count: 1,
        min_entry_size: 1,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiHochberg,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };

    let (pub_tx, _pub_rx) = mpsc::channel();
    let (kegg_tx, _kegg_rx) = mpsc::channel();
    let target = metabolopan::app::AnalysisPayload::Pathway { species_kegg };
    let output = run_stage3(
        &pubchem,
        &kegg,
        &[pos, neg],
        &target,
        params,
        pub_tx,
        kegg_tx,
    )
    .await
    .expect("run_stage3 dual-mode ok");

    // Universe N = {C_AGREE, C_POSONLY, C_NEGONLY, C_CONFL} (all four cpds
    // were reachable from at least one mode).
    assert_eq!(output.mapped_universe.len(), 4);
    for c in ["C_AGREE", "C_POSONLY", "C_NEGONLY", "C_CONFL"] {
        assert!(
            output.mapped_universe.contains(c),
            "expected {c} in universe; got: {:?}",
            output.mapped_universe
        );
    }

    // DualModeBreakdown should be Some in dual-mode.
    let b = output
        .dual_mode_breakdown
        .as_ref()
        .expect("expected DualModeBreakdown in dual-mode run");
    // Universe partition:
    //   POS-only: C_POSONLY only      = 1
    //   NEG-only: C_NEGONLY only      = 1
    //   in-both:  C_AGREE + C_CONFL   = 2
    assert_eq!(b.universe_pos_only, 1, "breakdown: {b:?}");
    assert_eq!(b.universe_neg_only, 1, "breakdown: {b:?}");
    assert_eq!(b.universe_in_both, 2, "breakdown: {b:?}");

    // Foreground K (direction = Up):
    //   agree_both: C_AGREE              = 1
    //   pos_only:   C_POSONLY            = 1
    //   neg_only:   C_NEGONLY            = 1
    //   excluded:   C_CONFL (Up vs Down) = 1
    assert_eq!(b.foreground_agree_both, 1, "breakdown: {b:?}");
    assert_eq!(b.foreground_pos_only, 1, "breakdown: {b:?}");
    assert_eq!(b.foreground_neg_only, 1, "breakdown: {b:?}");
    assert_eq!(b.foreground_excluded_conflict, 1, "breakdown: {b:?}");

    // Enriched pathway picks up the 3 non-conflict cpds.
    let row = output
        .enrichment_result
        .rows
        .iter()
        .find(|r| r.entry_id == "syn00001")
        .expect("syn00001 row");
    assert_eq!(row.hits, 3, "expected 3 hits (C_AGREE/POSONLY/NEGONLY)");
}

/// Regression bridge: when dual-mode is invoked with one mode's features
/// fully overlapping the other (effectively single-mode redundantly), the
/// universe / K / hits match what a single-mode call would produce. This
/// guards against drift in the dual-mode plumbing.
#[tokio::test]
async fn dual_mode_with_redundant_modes_matches_single_mode_universe() {
    let _g = serial().await;
    let (_dp, _dk) = setup_tmp_root();

    let pubchem_server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path("/compound/inchikey/property/InChIKey/CSV"))
        .respond_with(
            ResponseTemplate::new(200)
                .set_body_string("\"InChIKey\",\"CID\"\n\"K1\",\"200\"\n\"K2\",\"201\"\n"),
        )
        .mount(&pubchem_server)
        .await;
    let kegg_server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path_regex(r"^/conv/compound/.+"))
        .respond_with(
            ResponseTemplate::new(200)
                .set_body_string("pubchem:200\tcpd:C00001\npubchem:201\tcpd:C00002\n"),
        )
        .mount(&kegg_server)
        .await;

    let dam = DamResult {
        method: DamMethod::Welch,
        numerator: "T".into(),
        denominator: "C".into(),
        features: vec![feat("K1", 2.0, 0.001), feat("K2", 2.0, 0.001)],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiHochberg,
        dedup_report: None,
    };
    let species_kegg = SpeciesKegg {
        code: "syn".into(),
        fetched_at: chrono::Utc::now(),
        pathways: vec![KeggCompoundSet {
            id: "syn00001".into(),
            name: "Test".into(),
            compounds: vec!["C00001".into(), "C00002".into()],
        }],
    };

    let pubchem = PubchemClient::with_base_url(pubchem_server.uri());
    let kegg = KeggClient::with_base_url(Url::parse(&kegg_server.uri()).unwrap());
    let params = Stage3Params {
        method: DamMethod::Welch,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        direction: EnrichmentDirection::Both,
        min_hit_count: 1,
        min_entry_size: 1,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiHochberg,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };
    let (pub_tx, _pub_rx) = mpsc::channel();
    let (kegg_tx, _kegg_rx) = mpsc::channel();
    let target = metabolopan::app::AnalysisPayload::Pathway { species_kegg };
    let dual = run_stage3(
        &pubchem,
        &kegg,
        &[dam.clone(), dam.clone()],
        &target,
        params,
        pub_tx,
        kegg_tx,
    )
    .await
    .expect("dual run");

    // Universe = {C00001, C00002}; same in dual vs single. K = both cpds.
    assert_eq!(dual.mapped_universe.len(), 2);
    let b = dual.dual_mode_breakdown.as_ref().expect("dual mode");
    assert_eq!(b.universe_in_both, 2);
    assert_eq!(b.universe_pos_only, 0);
    assert_eq!(b.universe_neg_only, 0);
    assert_eq!(b.foreground_agree_both, 2);
    assert_eq!(b.foreground_excluded_conflict, 0);
}
