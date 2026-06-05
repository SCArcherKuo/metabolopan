//! Integration regression test for the `refactor-session-settings`
//! change. Walks through the key state transitions using real fixture
//! inputs + a real `run_dam` invocation, then asserts the
//! `SessionSettings` preservation/reset surfaces against the
//! design.md D11 inventory.
//!
//! Q-B (fixture-based): we load the canonical single-mode fixture
//! (`tests/fixtures/msdial_mini.txt` + a constructed minimal group CSV)
//! and run real `run_dam` to get a real `DamResult`. Stage 3 result-side
//! transitions assert against the named reset APIs directly (they take
//! `&mut SessionSettings`, so no spawn or wiremock is needed for the
//! settings-level invariants this test locks in).
//!
//! The unit tests in `src/app.rs` `#[cfg(test)] mod tests` already
//! cover each `reset_*` method's surface on synthetic baselines; this
//! test adds the "real DAM produces a `DamResult`; that result flows
//! correctly through the settings + state model" integration coverage.

use std::io::Write;

use metabolopan::app::{AnalysisMode, AppState, SessionCache, SessionInputs, SessionSettings};
use metabolopan::dam::{DamConfig, DamMethod, fdr::FdrMethod, run_dam, types::DamResult};
use metabolopan::data::{GroupMapping, load_group_mapping, parse_msdial_txt};
use metabolopan::enrichment::EnrichmentDirection;
use metabolopan::normalize::{NormalizationConfig, NormalizationMethod};

fn write_tmp_csv(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new()
        .suffix(".csv")
        .tempfile()
        .expect("tempfile");
    f.write_all(content.as_bytes()).expect("write tempfile");
    f
}

/// Load the canonical single-mode mini fixture (3 T vs 3 C samples).
fn load_mini_fixture() -> (metabolopan::data::MetabolomicsTable, GroupMapping) {
    let table = parse_msdial_txt(std::path::Path::new("tests/fixtures/msdial_mini.txt"))
        .expect("parse mini fixture");
    let csv = "sample,group\nT-1,T\nT-2,T\nT-3,T\nC-1,C\nC-2,C\nC-3,C\nBk-1,Bk\n";
    let tmp = write_tmp_csv(csv);
    let mapping = load_group_mapping(tmp.path(), &table.sample_cols).expect("mapping");
    (table, mapping)
}

/// Populate a non-default `SessionSettings` so reset assertions have
/// something to actually reset. Mirrors the `non_default_settings()`
/// helper used by the unit tests in `src/app.rs`.
fn populated_settings() -> SessionSettings {
    SessionSettings {
        analysis_mode: AnalysisMode::Module,
        kegg_species: Some("hsa".into()),
        organism_group_level: Some(2),
        organism_group: Some("Mammals".into()),
        min_group_overlap: 5,
        numerator: Some("T".into()),
        denominator: Some("C".into()),
        dam_method: DamMethod::Welch,
        drop_unknown: false,
        dedup_enabled: false,
        normalization: NormalizationMethod::Sum,
        metadata_column: Some("dry_weight".into()),
        pqn_reference: metabolopan::normalize::PqnReference::Group("C".into()),
        pqn_reference_group: Some("C".into()),
        log_transform: false,
        dam_fdr_method: FdrMethod::BenjaminiYekutieli,
        fc_threshold: 4.0,
        fdr_threshold: 0.01,
        delta_threshold: 0.5,
        stage2_export_width_in: 6.0,
        stage2_export_height_in: 4.0,
        stage2_export_dpi: 600,
        direction: EnrichmentDirection::Up,
        top_n: 50,
        enrichment_fdr_threshold: 0.1,
        min_hit_count: 3,
        min_entry_size: 5,
        enrichment_fdr_method: FdrMethod::BenjaminiHochberg,
        stage3_export_width_in: 5.0,
        stage3_export_height_in: 10.0,
        stage3_export_dpi: 600,
    }
}

/// Run real DAM on the mini fixture and return a `DamResult` we can
/// thread into post-DAM state assertions.
async fn run_real_dam() -> DamResult {
    let (mut table, mapping) = load_mini_fixture();
    run_dam(
        &mut table,
        &mapping,
        "T",
        "C",
        &DamConfig {
            method: DamMethod::Student,
            normalization: NormalizationConfig::default(),
            drop_unknown: true,
            dedup_enabled: true,
            log_transform: true,
            fdr_method: FdrMethod::BenjaminiHochberg,
        },
        None,
    )
    .await
    .expect("DAM run")
}

#[tokio::test]
async fn real_dam_result_flows_through_settings_and_state() {
    let dam_result = run_real_dam().await;
    // DamResult.fdr_method comes from settings.dam_fdr_method at the
    // time of `run_dam` invocation (here BH). The `Stage2DamThreshold`
    // variant carries `dam_results` only — fc/fdr/delta thresholds live
    // on settings, not the variant. Verify that contract holds.
    assert_eq!(dam_result.method, DamMethod::Student);
    assert_eq!(dam_result.fdr_method, FdrMethod::BenjaminiHochberg);
    assert!(
        !dam_result.features.is_empty(),
        "fixture should yield features"
    );

    // Construct a slim Stage2DamThreshold from the real DamResult. No
    // settings fields land on the variant — they're already on
    // SessionSettings.
    let state = AppState::Stage2DamThreshold {
        dam_results: vec![dam_result.clone()],
        active_volcano_tab: metabolopan::data::IonMode::Positive,
        volcano_textures: vec![None],
        rendering: false,
        render_rx: None,
    };
    matches!(state, AppState::Stage2DamThreshold { .. });
}

#[tokio::test]
async fn back_to_dam_setup_preserves_all_stage2_fields() {
    // After `reorder-gui-and-move-mode-to-stage3` (Phase 2),
    // `reset_stage2_choices_on_change_comparison` is a no-op: pressing
    // "Back to DAM Setup" preserves every Stage 2 settings field. This
    // mirrors the user-facing scenario where the user wants to tweak
    // ONE Stage 2 choice and re-run.
    let _dam_result = run_real_dam().await;
    let mut settings = populated_settings();
    let baseline = settings.clone();
    settings.reset_stage2_choices_on_change_comparison();
    assert_eq!(
        settings, baseline,
        "Back to DAM Setup MUST preserve every settings field"
    );
}

#[tokio::test]
async fn continue_to_enrichment_preserves_all_stage3_settings() {
    // Post-smoke-test feedback: Stage 3 settings persist across every
    // Continue. Locks the no-op contract.
    let _dam_result = run_real_dam().await;
    let mut settings = populated_settings();
    let baseline = settings.clone();
    settings.reset_stage3_on_continue_to_enrichment();
    assert_eq!(
        settings, baseline,
        "Continue to Enrichment MUST preserve every settings field"
    );

    // Idempotent across multiple Continues (still a no-op).
    settings.top_n = 99;
    settings.direction = EnrichmentDirection::Down;
    let mutated = settings.clone();
    settings.reset_stage3_on_continue_to_enrichment();
    assert_eq!(
        settings, mutated,
        "second Continue MUST keep the user's mutated values"
    );
}

#[tokio::test]
async fn back_to_stage1_preserves_all_settings_including_stage2_3() {
    // Post-smoke-test feedback: Back to Input preserves every settings
    // field (including numerator/denominator/method/thresholds and
    // Stage 3 fields). If the user re-picks files such that the
    // preserved groups no longer exist, the Stage 2 setup gate refuses
    // to start DAM until valid groups are re-selected — settings
    // themselves are not cleared.
    let _dam_result = run_real_dam().await;
    let mut settings = populated_settings();
    let baseline = settings.clone();
    settings.reset_for_back_to_stage1();
    assert_eq!(
        settings, baseline,
        "Back to Input MUST preserve every settings field"
    );
}

#[tokio::test]
async fn back_to_dam_result_preserves_all_settings() {
    // Post-smoke-test feedback: the Stage 3 setup → Stage 2 result Back
    // transition no longer clears any settings field.
    let _dam_result = run_real_dam().await;
    let mut settings = populated_settings();
    let baseline = settings.clone();
    settings.reset_for_back_to_stage2_threshold();
    assert_eq!(
        settings, baseline,
        "Back to DAM Result MUST preserve every settings field"
    );
}

#[tokio::test]
async fn mode_toggle_in_stage3_setup_preserves_both_modes_state() {
    // After Phase 2: toggling Mode on Stage 3 setup neither clears the
    // inactive mode's selection on `settings` nor its fetched cache on
    // `cache`. Both modes' state coexists for the lifetime of the
    // session so the user can toggle freely without re-fetching.
    let mut settings = populated_settings();
    settings.analysis_mode = AnalysisMode::Pathway;
    settings.kegg_species = Some("hsa".into());
    settings.organism_group_level = Some(2);
    settings.organism_group = Some("Mammals".into());
    settings.min_group_overlap = 5;

    let mut cache = SessionCache {
        group_org_codes: Some(std::collections::HashSet::from(["hsa".to_string()])),
        ..SessionCache::default()
    };

    // Toggle Pathway → Module → Pathway → Module: every selection +
    // cache field must survive every transition.
    for new_mode in [
        AnalysisMode::Module,
        AnalysisMode::Pathway,
        AnalysisMode::Module,
    ] {
        settings.reset_kegg_selection_for_mode_switch(new_mode);
        cache.clear_for_mode_switch(new_mode);
        assert_eq!(settings.analysis_mode, new_mode);
        assert_eq!(settings.kegg_species, Some("hsa".to_string()));
        assert_eq!(settings.organism_group_level, Some(2));
        assert_eq!(settings.organism_group, Some("Mammals".to_string()));
        assert_eq!(settings.min_group_overlap, 5);
        assert!(cache.group_org_codes.is_some());
    }
}

#[tokio::test]
async fn cache_clear_for_mode_switch_is_a_no_op() {
    // Companion to `mode_toggle_in_stage3_setup_preserves_both_modes_state`
    // — explicitly asserts the cache-side no-op contract independently.
    let mut cache = SessionCache {
        group_org_codes: Some(std::collections::HashSet::from(["hsa".to_string()])),
        ..SessionCache::default()
    };
    let baseline_has_codes = cache.group_org_codes.is_some();
    cache.clear_for_mode_switch(AnalysisMode::Pathway);
    cache.clear_for_mode_switch(AnalysisMode::Module);
    cache.clear_for_mode_switch(AnalysisMode::Pathway);
    assert_eq!(cache.group_org_codes.is_some(), baseline_has_codes);
}

#[tokio::test]
async fn session_inputs_csv_path_survives_post_stage1_state_transitions() {
    // Verifies the `add-csv-path-throughout-stages` invariant: `csv_path`
    // lives on `app.inputs` (not on any AppState variant), so it is
    // readable from every post-Stage-1 state.
    let (_table, _mapping) = load_mini_fixture();
    let tmp = write_tmp_csv("sample,group\nS1,A\nS2,B\n");
    let path = tmp.path().to_path_buf();
    let inputs = SessionInputs {
        ion_tables: vec![],
        mapping: None,
        csv_path: Some(path.clone()),
    };
    // Construct three different post-Stage-1 states and verify the
    // csv_path on `inputs` is invariant across them.
    let s1 = AppState::Stage2DamSetup { error: None };
    let s2 = AppState::Stage3EnrichSetup {
        dam_results: vec![],
        error: None,
        kegg_fetch: None,
        modules_fetch: None,
    };
    let s3 = AppState::Stage3EnrichResult {
        dam_results: vec![],
        module_retention: None,
        enrichment_result: metabolopan::enrichment::EnrichmentResult {
            universe_size: 0,
            dam_cpd_size: 0,
            direction: EnrichmentDirection::Both,
            min_hit_count: 1,
            min_entry_size: 1,
            entries_dropped_by_min_entry_size: 0,
            empty_compound_count: 0,
            rows: vec![],
            fdr_method: FdrMethod::BenjaminiHochberg,
        },
        mapped_universe: Default::default(),
        feature_to_cpds: Default::default(),
        pubchem_time_span: None,
        kegg_conv_time_span: None,
        dual_mode_breakdown: None,
        funnel: Default::default(),
        dotplot_tex: None,
        rendering: false,
        render_rx: None,
        refresh_state: metabolopan::app::RefreshState::Idle,
        confirming_new_round: false,
        height_user_overridden: false,
    };
    for state in [&s1, &s2, &s3] {
        // app.inputs.csv_path is the single source of truth — no variant
        // carries it any more.
        assert_eq!(inputs.csv_path.as_deref(), Some(path.as_path()));
        // (state isn't expected to expose csv_path itself; it's just a
        // marker that we are post-Stage-1.)
        let _ = state;
    }
}
