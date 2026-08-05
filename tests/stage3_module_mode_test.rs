//! Integration tests for Track E: Stage 3 orchestrator module-mode branch.
//!
//! Covers `run_stage3` + `AnalysisTarget::Module` behaviour from the
//! `stage3-ui` spec:
//! - Modules are filtered by `|complete_orgs ∩ group_orgs| >= min_overlap`.
//! - `ModuleRetention` summary is populated on `Stage3RunOutput`.
//! - Pathway-mode regression: unchanged behaviour when `AnalysisTarget::Pathway`.
//! - Empty Group → empty result.
//! - CSV export uses the renamed `EntryID,EntryName,...` header (Track A
//!   carryover, verified end-to-end here for module mode).

use metabolopan::app::AnalysisPayload;
use metabolopan::dam::types::{DamFeature, FcBasis};
use metabolopan::dam::{DamMethod, DamResult};
use metabolopan::enrichment::EnrichmentDirection;
use metabolopan::enrichment::export::export_csv;
use metabolopan::kegg::KeggClient;
use metabolopan::kegg::types::{KeggModuleEntry, KeggModulesCache};
use metabolopan::pubchem::PubchemClient;
use metabolopan::stage3::{Stage3Params, run_stage3};

use chrono::Utc;
use reqwest::Url;
use std::collections::HashSet;
use std::sync::OnceLock;
use std::sync::mpsc;
use wiremock::matchers::{method, path, path_regex};
use wiremock::{Mock, MockServer, ResponseTemplate};

async fn serial() -> tokio::sync::MutexGuard<'static, ()> {
    static M: OnceLock<tokio::sync::Mutex<()>> = OnceLock::new();
    M.get_or_init(|| tokio::sync::Mutex::new(())).lock().await
}

fn setup_tmp_cache_dirs() -> (tempfile::TempDir, tempfile::TempDir) {
    let dp = tempfile::tempdir().expect("pubchem tempdir");
    let dk = tempfile::tempdir().expect("kegg tempdir");
    unsafe {
        std::env::set_var("PUBCHEM_CACHE_DIR", dp.path());
        std::env::set_var("KEGG_CACHE_DIR", dk.path());
    }
    (dp, dk)
}

fn dam_feature(inchikey: &str, log2_fc: f64, p_adj: f64) -> DamFeature {
    DamFeature {
        alignment_id: format!("aid-{inchikey}"),
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

fn module_entry(name: &str, compounds: &[&str], orgs: &[&str]) -> KeggModuleEntry {
    KeggModuleEntry {
        name: name.to_string(),
        compounds: compounds.iter().map(|s| s.to_string()).collect(),
        complete_orgs: orgs.iter().map(|s| s.to_string()).collect(),
        fetched_at: Utc::now(),
    }
}

fn build_modules_cache() -> KeggModulesCache {
    let mut cache = KeggModulesCache::default();
    // M00001: Has K1+K2 hits, complete in hsa (in Animals group).
    cache.modules.insert(
        "M00001".to_string(),
        module_entry(
            "Mod 1 (Animals hit)",
            &["C00001", "C00002"],
            &["hsa", "ptr"],
        ),
    );
    // M00002: No K hits, complete in hsa.
    cache.modules.insert(
        "M00002".to_string(),
        module_entry("Mod 2 (Animals no-hit)", &["C00099"], &["hsa"]),
    );
    // M00003: Has K2 hit, complete in ath only (NOT in Animals).
    cache.modules.insert(
        "M00003".to_string(),
        module_entry("Mod 3 (Plants-only)", &["C00002"], &["ath"]),
    );
    // M00004: Empty COMPOUND, complete in hsa.
    cache.modules.insert(
        "M00004".to_string(),
        module_entry("Mod 4 (empty COMPOUND)", &[], &["hsa"]),
    );
    cache
}

#[tokio::test]
async fn module_mode_filters_modules_by_group_overlap() {
    let _g = serial().await;
    let (_dp, _dk) = setup_tmp_cache_dirs();

    // Mock PubChem: K1 → 1001, K2 → 1002.
    let pubchem_server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path("/compound/inchikey/property/InChIKey/CSV"))
        .respond_with(
            ResponseTemplate::new(200)
                .set_body_string("\"InChIKey\",\"CID\"\n\"K1\",\"1001\"\n\"K2\",\"1002\"\n"),
        )
        .mount(&pubchem_server)
        .await;

    // Mock KEGG conv: 1001 → C00001, 1002 → C00002.
    let kegg_server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path_regex(r"^/conv/compound/.+"))
        .respond_with(
            ResponseTemplate::new(200)
                .set_body_string("pubchem:1001\tcpd:C00001\npubchem:1002\tcpd:C00002\n"),
        )
        .mount(&kegg_server)
        .await;

    let dam_result = DamResult {
        method: DamMethod::Welch,
        numerator: "T".into(),
        denominator: "C".into(),
        features: vec![dam_feature("K1", 2.0, 0.001), dam_feature("K2", 2.5, 0.001)],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        dedup_report: None,
    };

    let pubchem = PubchemClient::with_base_url(pubchem_server.uri());
    let kegg = KeggClient::with_base_url(Url::parse(&kegg_server.uri()).unwrap());

    let modules_pack = build_modules_cache();
    let group_org_codes: HashSet<String> = ["hsa", "ptr", "mmu"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let target = AnalysisPayload::Module {
        modules_pack,
        group_level: 2,
        group_name: "Animals".to_string(),
        group_org_codes,
        min_group_overlap: 1,
    };

    let params = Stage3Params {
        method: DamMethod::Welch,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        direction: EnrichmentDirection::Both,
        min_hit_count: 1,
        min_entry_size: 1,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };

    let (pub_tx, _pub_rx) = mpsc::channel();
    let (kegg_tx, _kegg_rx) = mpsc::channel();
    let output = run_stage3(
        &pubchem,
        &kegg,
        std::slice::from_ref(&dam_result),
        &target,
        params,
        pub_tx,
        kegg_tx,
    )
    .await
    .expect("run_stage3 ok");

    // ── Retention summary populated ──
    let retention = output.module_retention.expect("module mode → Some");
    assert_eq!(retention.total_modules, 4);
    // M00001 + M00002 + M00004 are in Animals (hsa); M00003 (ath only) excluded.
    assert_eq!(retention.retained_modules, 3);
    assert_eq!(retention.min_group_overlap, 1);
    assert_eq!(retention.group_level, 2);
    assert_eq!(retention.group_name, "Animals");
    assert_eq!(retention.group_org_count, 3);
    // Time span sanity: oldest ≤ newest.
    assert!(retention.oldest_fetched_at <= retention.newest_fetched_at);

    // ── Empty-COMPOUND count surfaces ──
    // M00004 has empty compounds and was retained → empty_compound_count = 1.
    assert_eq!(output.enrichment_result.empty_compound_count, 1);

    // ── M00001 is the only module with hits ──
    let m1_row = output
        .enrichment_result
        .rows
        .iter()
        .find(|r| r.entry_id == "M00001")
        .expect("M00001 row");
    assert_eq!(m1_row.hits, 2); // C00001 + C00002 both hit
    assert_eq!(m1_row.total, 2); // pathway-restricted-to-universe size

    // ── M00003 (Plants-only) was filtered out before ORA ──
    assert!(
        output
            .enrichment_result
            .rows
            .iter()
            .all(|r| r.entry_id != "M00003"),
        "M00003 should be excluded by Group filter (Plants only, not in Animals)"
    );

    // ── CSV header uses Track-A-renamed columns; leading `# FDR: …` tag is line 0 ──
    let mut buf = Vec::new();
    export_csv(&mut buf, &output.enrichment_result).expect("export csv");
    let s = String::from_utf8(buf).unwrap();
    let lines: Vec<&str> = s.lines().collect();
    assert!(
        lines[0].starts_with("# FDR: "),
        "FDR tag on first line; got {}",
        lines[0]
    );
    assert!(
        lines[1].starts_with("# MinEntrySize: "),
        "MinEntrySize tag on second line; got {}",
        lines[1]
    );
    assert_eq!(
        lines[2],
        "EntryID,EntryName,Hits,Total,Expected,EnrichmentRatio,PValue,FDR,HitKeggIDs"
    );
    // Rows carry module IDs in the EntryID column.
    assert!(s.contains("M00001,"));
}

#[tokio::test]
async fn module_mode_empty_group_yields_empty_result() {
    let _g = serial().await;
    let (_dp, _dk) = setup_tmp_cache_dirs();

    // Trivial PubChem + KEGG conv mocks (universe will be non-empty).
    let pubchem_server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path("/compound/inchikey/property/InChIKey/CSV"))
        .respond_with(
            ResponseTemplate::new(200).set_body_string("\"InChIKey\",\"CID\"\n\"K1\",\"1001\"\n"),
        )
        .mount(&pubchem_server)
        .await;
    let kegg_server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path_regex(r"^/conv/compound/.+"))
        .respond_with(ResponseTemplate::new(200).set_body_string("pubchem:1001\tcpd:C00001\n"))
        .mount(&kegg_server)
        .await;

    let dam_result = DamResult {
        method: DamMethod::Welch,
        numerator: "T".into(),
        denominator: "C".into(),
        features: vec![dam_feature("K1", 2.0, 0.001)],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        dedup_report: None,
    };
    let pubchem = PubchemClient::with_base_url(pubchem_server.uri());
    let kegg = KeggClient::with_base_url(Url::parse(&kegg_server.uri()).unwrap());

    // Empty Group.
    let target = AnalysisPayload::Module {
        modules_pack: build_modules_cache(),
        group_level: 1,
        group_name: "EmptyGroup".to_string(),
        group_org_codes: HashSet::new(),
        min_group_overlap: 1,
    };
    let params = Stage3Params {
        method: DamMethod::Welch,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        direction: EnrichmentDirection::Both,
        min_hit_count: 1,
        min_entry_size: 1,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };
    let (pub_tx, _pub_rx) = mpsc::channel();
    let (kegg_tx, _kegg_rx) = mpsc::channel();
    let output = run_stage3(
        &pubchem,
        &kegg,
        std::slice::from_ref(&dam_result),
        &target,
        params,
        pub_tx,
        kegg_tx,
    )
    .await
    .expect("run_stage3 ok");

    // Universe is still the measurable metabolome (non-empty), but no
    // entries were retained → no rows.
    assert!(output.enrichment_result.rows.is_empty());
    assert_eq!(output.enrichment_result.universe_size, 1);
    let retention = output.module_retention.expect("module mode");
    assert_eq!(retention.retained_modules, 0);
    assert_eq!(retention.group_org_count, 0);
}

#[tokio::test]
async fn module_mode_strict_overlap_threshold_filters_more() {
    let _g = serial().await;
    let (_dp, _dk) = setup_tmp_cache_dirs();

    let pubchem_server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path("/compound/inchikey/property/InChIKey/CSV"))
        .respond_with(
            ResponseTemplate::new(200).set_body_string("\"InChIKey\",\"CID\"\n\"K1\",\"1001\"\n"),
        )
        .mount(&pubchem_server)
        .await;
    let kegg_server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path_regex(r"^/conv/compound/.+"))
        .respond_with(ResponseTemplate::new(200).set_body_string("pubchem:1001\tcpd:C00001\n"))
        .mount(&kegg_server)
        .await;

    let dam_result = DamResult {
        method: DamMethod::Welch,
        numerator: "T".into(),
        denominator: "C".into(),
        features: vec![dam_feature("K1", 2.0, 0.001)],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        dedup_report: None,
    };
    let pubchem = PubchemClient::with_base_url(pubchem_server.uri());
    let kegg = KeggClient::with_base_url(Url::parse(&kegg_server.uri()).unwrap());

    // Group of 3 orgs; with min_overlap=2, only modules with ≥2 of those
    // orgs in COMPLETE survive. M00001 has {hsa, ptr} (both in Group);
    // M00002 has only {hsa} (1 overlap → fails); M00004 has only {hsa}
    // (1 overlap → fails).
    let group_org_codes: HashSet<String> = ["hsa", "ptr", "mmu"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let target = AnalysisPayload::Module {
        modules_pack: build_modules_cache(),
        group_level: 2,
        group_name: "Animals".to_string(),
        group_org_codes,
        min_group_overlap: 2,
    };
    let params = Stage3Params {
        method: DamMethod::Welch,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        direction: EnrichmentDirection::Both,
        min_hit_count: 1,
        min_entry_size: 1,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };
    let (pub_tx, _pub_rx) = mpsc::channel();
    let (kegg_tx, _kegg_rx) = mpsc::channel();
    let output = run_stage3(
        &pubchem,
        &kegg,
        std::slice::from_ref(&dam_result),
        &target,
        params,
        pub_tx,
        kegg_tx,
    )
    .await
    .expect("run_stage3 ok");

    let retention = output.module_retention.expect("module mode");
    // Only M00001 has ≥2 overlap → retained_modules = 1.
    assert_eq!(retention.retained_modules, 1);
    assert_eq!(retention.min_group_overlap, 2);
}

#[tokio::test]
async fn pathway_mode_module_retention_is_none() {
    let _g = serial().await;
    let (_dp, _dk) = setup_tmp_cache_dirs();

    let pubchem_server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path("/compound/inchikey/property/InChIKey/CSV"))
        .respond_with(
            ResponseTemplate::new(200).set_body_string("\"InChIKey\",\"CID\"\n\"K1\",\"1001\"\n"),
        )
        .mount(&pubchem_server)
        .await;
    let kegg_server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path_regex(r"^/conv/compound/.+"))
        .respond_with(ResponseTemplate::new(200).set_body_string("pubchem:1001\tcpd:C00001\n"))
        .mount(&kegg_server)
        .await;

    let dam_result = DamResult {
        method: DamMethod::Welch,
        numerator: "T".into(),
        denominator: "C".into(),
        features: vec![dam_feature("K1", 2.0, 0.001)],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        dedup_report: None,
    };
    let pubchem = PubchemClient::with_base_url(pubchem_server.uri());
    let kegg = KeggClient::with_base_url(Url::parse(&kegg_server.uri()).unwrap());

    let species_kegg = metabolopan::kegg::SpeciesKegg {
        code: "syn".into(),
        fetched_at: Utc::now(),
        pathways: vec![metabolopan::kegg::KeggCompoundSet {
            id: "syn00001".into(),
            name: "Test pathway".into(),
            compounds: vec!["C00001".into()],
        }],
    };
    let target = AnalysisPayload::Pathway { species_kegg };
    let params = Stage3Params {
        method: DamMethod::Welch,
        fc_threshold: 2.0,
        fdr_threshold: 0.05,
        delta_threshold: 0.33,
        direction: EnrichmentDirection::Both,
        min_hit_count: 1,
        min_entry_size: 1,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };
    let (pub_tx, _pub_rx) = mpsc::channel();
    let (kegg_tx, _kegg_rx) = mpsc::channel();
    let output = run_stage3(
        &pubchem,
        &kegg,
        std::slice::from_ref(&dam_result),
        &target,
        params,
        pub_tx,
        kegg_tx,
    )
    .await
    .expect("run_stage3 ok");

    // Pathway mode → no retention summary.
    assert!(output.module_retention.is_none());
    // Pathway mode tested existing pathway with C00001 → one hit.
    let row = output
        .enrichment_result
        .rows
        .iter()
        .find(|r| r.entry_id == "syn00001")
        .expect("syn00001 row");
    assert_eq!(row.hits, 1);
}
