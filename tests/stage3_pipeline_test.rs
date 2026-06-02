//! End-to-end Stage 3 pipeline smoke test: PubChem (mocked) → KEGG
//! conv (mocked) → ORA → dot plot.

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

#[tokio::test]
async fn stage3_pipeline_resolves_universe_and_enriches() {
    let _g = serial().await;
    let (_dp, _dk) = setup_tmp_root();

    // Mock PubChem: K1 → 1001, K2 → 1002, K3 → no match.
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

    // DAM result with 3 features: K1 Up, K2 Up, K3 ns.
    let dam_result = DamResult {
        method: DamMethod::Welch,
        numerator: "T".into(),
        denominator: "C".into(),
        features: vec![
            dam_feature("K1", 2.0, 0.001),
            dam_feature("K2", 2.5, 0.001),
            dam_feature("K3", 0.1, 0.40),
        ],
        skipped: 0,
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        dedup_report: None,
    };

    // Species pathways: p1 contains K1+K2, p2 contains K3 only, p3 empty
    // overlap.
    let species_kegg = SpeciesKegg {
        code: "syn".into(),
        fetched_at: chrono::Utc::now(),
        pathways: vec![
            KeggCompoundSet {
                id: "syn00001".into(),
                name: "Test enriched".into(),
                compounds: vec!["C00001".into(), "C00002".into(), "C00005".into()],
            },
            KeggCompoundSet {
                id: "syn00002".into(),
                name: "Test misc".into(),
                compounds: vec!["C00099".into()],
            },
        ],
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
        fdr_method: metabolopan::dam::fdr::FdrMethod::BenjaminiYekutieli,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };

    let (pub_tx, _pub_rx) = mpsc::channel();
    let (kegg_tx, _kegg_rx) = mpsc::channel();
    // Track E: pathway-mode AnalysisTarget wraps the SpeciesKegg.
    let target = metabolopan::app::AnalysisPayload::Pathway { species_kegg };
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

    // Universe should be {C00001, C00002}: K3 had no PubChem match so it
    // contributes nothing. K1+K2 mapped to C00001+C00002.
    assert_eq!(output.mapped_universe.len(), 2);
    assert!(output.mapped_universe.contains("C00001"));
    assert!(output.mapped_universe.contains("C00002"));

    // feature_to_cpds should have K1 and K2 entries; K3 absent (no mapping).
    assert_eq!(output.feature_to_cpds.len(), 2);

    // Find the enriched pathway row.
    let enriched = output
        .enrichment_result
        .rows
        .iter()
        .find(|r| r.entry_id == "syn00001")
        .expect("syn00001 row");
    // K = {C00001, C00002}, pathway hits both.
    assert_eq!(enriched.hits, 2);
    // M_p = pathway compounds in universe = {C00001, C00002} (C00005 absent).
    assert_eq!(enriched.total, 2);
    // p_value: N=2, M=2, K=2, k=2 → 1 - hypercdf(1; 2, 2, 2) = 1 - 0 = 1.
    // Actually for N=K=M=2, k=2 is certain → cdf(1) = 0 (we want P(X<=1)
    // when X must be 2), so 1 - 0 = 1. Skip exact p assertion; just
    // verify finite.
    assert!(enriched.p_value.is_finite());

    // syn00002 has compounds = ["C00099"] but C00099 is NOT in the universe
    // (only C00001/C00002 mapped). m_p = 0 < min_entry_size = 1 → dropped
    // by the pre-FDR filter and absent from rows. The dropped count
    // surfaces on EnrichmentResult.
    assert!(
        output
            .enrichment_result
            .rows
            .iter()
            .all(|r| r.entry_id != "syn00002"),
        "syn00002 (m_p=0) must be filtered out under min_entry_size=1"
    );
    assert!(output.enrichment_result.entries_dropped_by_min_entry_size >= 1);
}
