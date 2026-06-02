//! UI-side spawn-path tests for Stage 3. Covers the bug fixed by
//! `fix-stage3-ui-dual-mode-spawn`: a `Vec<DamResult>` of length 2 on
//! `Stage3EnrichSetup` MUST reach the orchestrator intact (not collapsed
//! via `slice::from_ref(&dam_results[0])`).
//!
//! Complements `tests/dual_mode_pipeline_test.rs` (orchestrator math layer)
//! by exercising the UI-layer plumbing — `build_stage3_spawn_inputs`,
//! `Stage3Params` construction, and the orchestrator hand-off with a
//! 2-element slice.

use std::collections::HashSet;
use std::sync::{OnceLock, mpsc};

use metabolopan::app::{AnalysisPayload, AppState, SessionCache, SessionInputs, SessionSettings};
use metabolopan::dam::types::{DamFeature, FcBasis};
use metabolopan::dam::{DamMethod, DamResult};
use metabolopan::kegg::{KeggClient, KeggCompoundSet, SpeciesKegg};
use metabolopan::pubchem::PubchemClient;
use metabolopan::stage3::{Stage3Params, run_stage3};
use metabolopan::ui::stage3_setup::build_stage3_spawn_inputs;
use reqwest::Url;
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

fn dam(numerator: &str, denominator: &str, features: Vec<DamFeature>) -> DamResult {
    DamResult {
        method: DamMethod::Welch,
        numerator: numerator.to_string(),
        denominator: denominator.to_string(),
        features,
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiHochberg,
        dedup_report: None,
    }
}

/// Tiny `SpeciesKegg` with a couple of pathways; enough to satisfy the
/// `build_analysis_payload` Pathway branch and feed `run_stage3`.
fn synth_species_kegg() -> SpeciesKegg {
    SpeciesKegg {
        code: "tst".into(),
        fetched_at: chrono::Utc::now(),
        pathways: vec![
            KeggCompoundSet {
                id: "tst00001".into(),
                name: "Test pathway 1".into(),
                compounds: vec!["C_AGREE".into(), "C_POSONLY".into(), "C_NEGONLY".into()],
            },
            KeggCompoundSet {
                id: "tst00002".into(),
                name: "Test pathway 2".into(),
                compounds: vec!["C_AGREE".into(), "C_NEGONLY".into()],
            },
        ],
    }
}

fn dual_mode_setup_state() -> (AppState, SessionSettings, SessionCache) {
    // POS: K_AGREE Up, K_POSONLY Up.
    let pos = dam(
        "T",
        "C",
        vec![feat("K_AGREE", 2.0, 0.001), feat("K_POSONLY", 2.5, 0.001)],
    );
    // NEG: K_AGREE Up, K_NEGONLY Up.
    let neg = dam(
        "T",
        "C",
        vec![feat("K_AGREE", 2.0, 0.001), feat("K_NEGONLY", 2.5, 0.001)],
    );

    let state = AppState::Stage3EnrichSetup {
        dam_results: vec![pos, neg],
        error: None,
        kegg_fetch: None,
        modules_fetch: None,
    };

    let settings = SessionSettings {
        kegg_species: Some("tst".into()),
        ..SessionSettings::default()
    };

    let cache = SessionCache {
        species_kegg: Some(synth_species_kegg()),
        ..SessionCache::default()
    };

    (state, settings, cache)
}

/// `build_stage3_spawn_inputs` MUST return the full 2-element Vec for dual-mode,
/// the additive `pubchem_total`, and a properly-wired `Stage3Params`. Pre-fix it
/// was impossible to test this — the helper didn't exist; `start_run` did
/// `slice::from_ref(&dam_results[0])` inline.
#[test]
fn build_stage3_spawn_inputs_preserves_dual_mode_vec() {
    let (state, settings, cache) = dual_mode_setup_state();
    let (dam_results_clone, params, target, pubchem_total) =
        build_stage3_spawn_inputs(&state, &settings, &cache)
            .expect("helper should succeed with populated cache");

    assert_eq!(
        dam_results_clone.len(),
        2,
        "dual-mode Vec MUST reach the orchestrator intact (len=2), not collapsed"
    );
    assert_eq!(dam_results_clone[0].features.len(), 2);
    assert_eq!(dam_results_clone[1].features.len(), 2);

    // params reads from dam_results[0].method but reflects the SessionSettings.
    assert_eq!(params.method, DamMethod::Welch);
    assert_eq!(params.min_entry_size, settings.min_entry_size);
    assert_eq!(params.direction, settings.direction);
    assert!(!params.force_refresh_pubchem);
    assert!(!params.force_refresh_kegg_conv);

    // pubchem_total = sum across modes of InChIKey-bearing features (additive,
    // per design D3).
    assert_eq!(
        pubchem_total, 4,
        "pubchem_total must sum across modes: 2 POS + 2 NEG = 4"
    );

    // target is Pathway with the synthesized species.
    match target {
        AnalysisPayload::Pathway { species_kegg } => {
            assert_eq!(species_kegg.code, "tst");
            assert_eq!(species_kegg.pathways.len(), 2);
        }
        AnalysisPayload::Module { .. } => {
            panic!("expected Pathway target, got Module");
        }
    }
}

/// Single-mode path: `build_stage3_spawn_inputs` returns a length-1 Vec and
/// `pubchem_total` reflects only that one mode. Regression baseline for users
/// running single-mode flows.
#[test]
fn build_stage3_spawn_inputs_single_mode_returns_len_one() {
    let pos = dam("T", "C", vec![feat("K_X", 2.0, 0.001)]);
    let state = AppState::Stage3EnrichSetup {
        dam_results: vec![pos],
        error: None,
        kegg_fetch: None,
        modules_fetch: None,
    };
    let settings = SessionSettings {
        kegg_species: Some("tst".into()),
        ..SessionSettings::default()
    };
    let cache = SessionCache {
        species_kegg: Some(synth_species_kegg()),
        ..SessionCache::default()
    };

    let (dam_results_clone, _params, _target, pubchem_total) =
        build_stage3_spawn_inputs(&state, &settings, &cache).expect("helper succeeds");
    assert_eq!(dam_results_clone.len(), 1);
    assert_eq!(pubchem_total, 1);
}

/// Wrong state variant → `None`. Locks the "expects Stage3EnrichSetup" precondition.
#[test]
fn build_stage3_spawn_inputs_wrong_state_returns_none() {
    let state = AppState::Stage1Input {
        slot1_mode: None,
        slot2_revealed: false,
        slot2_mode: None,
        error: None,
    };
    let settings = SessionSettings::default();
    let cache = SessionCache::default();
    assert!(build_stage3_spawn_inputs(&state, &settings, &cache).is_none());
}

/// Missing cache → `None`. Mirrors the start_run error-restoration branch.
#[test]
fn build_stage3_spawn_inputs_missing_cache_returns_none() {
    let (state, settings, _populated) = dual_mode_setup_state();
    let empty_cache = SessionCache::default(); // species_kegg = None
    assert!(build_stage3_spawn_inputs(&state, &settings, &empty_cache).is_none());
}

/// End-to-end: drive `run_stage3` with the helper's outputs against wiremock
/// servers. Asserts that the orchestrator took the dual-mode branch
/// (`dual_mode_breakdown.is_some()`) and that the universe arithmetic holds.
/// This is the durable regression anchor for the bug — pre-fix the helper
/// didn't exist and the inline spawn collapsed to len=1, so `dual_mode_breakdown`
/// was always `None` even with both ion tables loaded.
#[tokio::test]
async fn ui_dual_mode_spawn_reaches_orchestrator_with_breakdown() {
    let _g = serial().await;
    let (_dp, _dk) = setup_tmp_root();

    let pubchem_server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path("/compound/inchikey/property/InChIKey/CSV"))
        .respond_with(ResponseTemplate::new(200).set_body_string(
            "\"InChIKey\",\"CID\"\n\
             \"K_AGREE\",\"100\"\n\
             \"K_POSONLY\",\"101\"\n\
             \"K_NEGONLY\",\"102\"\n",
        ))
        .mount(&pubchem_server)
        .await;

    let kegg_server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path_regex(r"^/conv/compound/.+"))
        .respond_with(ResponseTemplate::new(200).set_body_string(
            "pubchem:100\tcpd:C_AGREE\n\
             pubchem:101\tcpd:C_POSONLY\n\
             pubchem:102\tcpd:C_NEGONLY\n",
        ))
        .mount(&kegg_server)
        .await;

    let (state, settings, cache) = dual_mode_setup_state();
    let (dam_results_clone, params, target, _pubchem_total) =
        build_stage3_spawn_inputs(&state, &settings, &cache).expect("helper succeeds");

    let pubchem = PubchemClient::with_base_url(pubchem_server.uri());
    let kegg = KeggClient::with_base_url(Url::parse(&kegg_server.uri()).unwrap());

    let (pub_tx, _pub_rx) = mpsc::channel();
    let (kegg_tx, _kegg_rx) = mpsc::channel();

    let output = run_stage3(
        &pubchem,
        &kegg,
        &dam_results_clone, // ← full 2-element slice
        &target,
        params,
        pub_tx,
        kegg_tx,
    )
    .await
    .expect("orchestrator ok");

    let breakdown = output
        .dual_mode_breakdown
        .as_ref()
        .expect("dual-mode breakdown MUST be Some(_) when len=2 reaches the orchestrator");

    // Universe partition arithmetic must hold.
    assert_eq!(
        breakdown.universe_pos_only + breakdown.universe_neg_only + breakdown.universe_in_both,
        output.mapped_universe.len(),
        "universe partition sum must equal |N|"
    );

    // K_AGREE was in both modes' DAM features → universe_in_both at least 1.
    assert!(
        breakdown.universe_in_both >= 1,
        "C_AGREE was mapped from both modes; expected universe_in_both >= 1, got {}",
        breakdown.universe_in_both
    );

    // Quick sanity: the PubChem mock saw the UNION of both modes' InChIKeys
    // — i.e. universe contains both C_POSONLY (POS-only) and C_NEGONLY (NEG-only).
    let universe: HashSet<&String> = output.mapped_universe.iter().collect();
    assert!(
        universe.contains(&"C_POSONLY".to_string()),
        "POS-only cpd missing from universe — UI may have dropped POS"
    );
    assert!(
        universe.contains(&"C_NEGONLY".to_string()),
        "NEG-only cpd missing from universe — UI silently dropped NEG (this is the bug)"
    );

    // The `_inputs` ref isn't used; suppress dead-code if rustc flags it.
    let _ = SessionInputs::default();
}

/// `run_stage3` MUST panic in debug builds when handed an oversized slice.
/// Locks the new `debug_assert!(dam_results.len() <= 2)` at the orchestrator entry.
#[tokio::test]
#[should_panic(expected = "1-or-2 mode contract")]
#[cfg(debug_assertions)]
async fn run_stage3_panics_in_debug_when_len_exceeds_two() {
    let _g = serial().await;
    let (_dp, _dk) = setup_tmp_root();
    let three = vec![
        dam("T", "C", vec![]),
        dam("T", "C", vec![]),
        dam("T", "C", vec![]),
    ];
    let target = AnalysisPayload::Pathway {
        species_kegg: synth_species_kegg(),
    };
    let params = Stage3Params {
        method: DamMethod::Welch,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        direction: metabolopan::enrichment::EnrichmentDirection::Both,
        min_hit_count: 1,
        min_entry_size: 1,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiHochberg,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };
    let (pub_tx, _pub_rx) = mpsc::channel();
    let (kegg_tx, _kegg_rx) = mpsc::channel();
    let _ = run_stage3(
        &PubchemClient::new(),
        &KeggClient::new(),
        &three,
        &target,
        params,
        pub_tx,
        kegg_tx,
    )
    .await;
}

/// Empty input returns `Err` (release path) — the `anyhow::bail!` predates the
/// new debug_assert and remains the deterministic behaviour. In debug builds
/// the assert fires first; this test asserts the contract via the `Result`.
#[tokio::test]
async fn run_stage3_empty_slice_returns_err_or_panics() {
    let _g = serial().await;
    let (_dp, _dk) = setup_tmp_root();
    let target = AnalysisPayload::Pathway {
        species_kegg: synth_species_kegg(),
    };
    let params = Stage3Params {
        method: DamMethod::Welch,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        direction: metabolopan::enrichment::EnrichmentDirection::Both,
        min_hit_count: 1,
        min_entry_size: 1,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiHochberg,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };
    let (pub_tx, _pub_rx) = mpsc::channel();
    let (kegg_tx, _kegg_rx) = mpsc::channel();
    // catch_unwind so we accept either the bail Err (release) or the
    // debug_assert panic (debug). Either is a valid contract enforcement.
    let result = std::panic::AssertUnwindSafe(async {
        run_stage3(
            &PubchemClient::new(),
            &KeggClient::new(),
            &[],
            &target,
            params,
            pub_tx,
            kegg_tx,
        )
        .await
    });
    // We're inside #[tokio::test], so just await + assert Err for the
    // release path. Debug builds will panic above the .await with the
    // pre-bail debug_assert; we don't try to catch that here — the
    // separate oversized-slice test covers the assert's message format.
    match result.0.await {
        Err(_) => { /* expected — anyhow::bail */ }
        Ok(_) => panic!("run_stage3 with empty slice should not return Ok"),
    }
}
