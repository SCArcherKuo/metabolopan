use std::time::Instant;
use tokio::sync::mpsc;
use wiremock::matchers::{method, path};
use wiremock::{Mock, MockServer, ResponseTemplate};

use metabolopan::kegg::{KeggClient, KeggProgress};

fn organism_body() -> &'static str {
    // KEGG BRITE "KEGG Organism" hierarchy (`/get/br:br08601`): leading char is
    // the taxonomy level A→E; indentation is significant, so this uses a
    // literal-newline string (no `\`-continuation, which would eat leading
    // whitespace). Group names carry a ` (count)` suffix the parser strips.
    "\
+E\tKEGG Organism
!
AEukaryotes (3)
B  Plants (3)
C    Eudicots (3)
D      Fabales (1)
E        gmx  Glycine max (soybean)
D      Solanales (1)
E        sly  Solanum lycopersicum (tomato)
D      Brassicales (1)
E        ath  Arabidopsis thaliana (thale cress)
#Last updated: June 18, 2026
"
}

fn pathway_list_body() -> &'static str {
    "path:gmx00010\tGlycolysis / Gluconeogenesis - Glycine max (soybean)\n\
     path:gmx00020\tCitrate cycle (TCA cycle) - Glycine max (soybean)\n\
     path:gmx00030\tPentose phosphate pathway - Glycine max (soybean)\n"
}

fn pathway_detail_with_compounds() -> &'static str {
    "\
ENTRY       gmx00010
NAME        Glycolysis / Gluconeogenesis
COMPOUND    C00001 Water
            C00002 ATP
            C00009 Orthophosphate
REFERENCE   ...
///
"
}

fn pathway_detail_no_compounds() -> &'static str {
    "\
ENTRY       gmx00020
NAME        Citrate cycle
REFERENCE   ...
///
"
}

fn build_client(server: &MockServer) -> KeggClient {
    KeggClient::with_base_url(server.uri().parse().expect("base url"))
}

#[tokio::test]
async fn list_organisms_parses_three_records() {
    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/get/br:br08601"))
        .respond_with(ResponseTemplate::new(200).set_body_string(organism_body()))
        .mount(&server)
        .await;

    let client = build_client(&server);
    let organisms = client.list_organisms().await.expect("list organisms");
    assert_eq!(organisms.len(), 3);
    assert_eq!(organisms[0].code, "gmx");
    // BRITE carries no T-numbers; the parser synthesizes `T_{code}`.
    assert_eq!(organisms[0].t_number, "T_gmx");
    assert!(organisms[0].name.contains("Glycine max"));
    assert!(organisms[0].lineage.contains("Fabales"));
}

#[tokio::test]
async fn fetch_species_returns_pathways_and_compounds() {
    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/pathway/gmx"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_list_body()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00010"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_detail_with_compounds()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00020"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_detail_no_compounds()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00030"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_detail_with_compounds()))
        .mount(&server)
        .await;

    let client = build_client(&server);
    let (tx, mut rx) = mpsc::channel::<KeggProgress>(64);

    let species = client
        .fetch_species_pathways("gmx", tx)
        .await
        .expect("fetch species");

    assert_eq!(species.code, "gmx");
    assert_eq!(species.pathways.len(), 3);
    // Pathway IDs come from `path:gmx0001X` lines stripped of the prefix.
    let ids: Vec<&str> = species.pathways.iter().map(|p| p.id.as_str()).collect();
    assert_eq!(ids, vec!["gmx00010", "gmx00020", "gmx00030"]);

    // First pathway has 3 compounds, second has 0, third has 3.
    assert_eq!(
        species.pathways[0].compounds,
        vec!["C00001", "C00002", "C00009"]
    );
    assert!(species.pathways[1].compounds.is_empty());
    assert_eq!(
        species.pathways[2].compounds,
        vec!["C00001", "C00002", "C00009"]
    );

    // Drain progress events: expect 3 (one per pathway).
    let mut events = Vec::new();
    while let Ok(p) = rx.try_recv() {
        events.push(p);
    }
    assert_eq!(
        events.len(),
        3,
        "expected 3 progress events, got {events:?}"
    );
    assert_eq!(
        events.iter().map(|e| e.completed).collect::<Vec<_>>(),
        vec![1, 2, 3]
    );
    assert!(events.iter().all(|e| e.total == 3));
}

#[tokio::test]
async fn fetch_species_propagates_http_error() {
    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/pathway/gmx"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_list_body()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00010"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_detail_with_compounds()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00020"))
        .respond_with(ResponseTemplate::new(500))
        .mount(&server)
        .await;

    let client = build_client(&server);
    let (tx, _rx) = mpsc::channel::<KeggProgress>(64);

    let err = client
        .fetch_species_pathways("gmx", tx)
        .await
        .expect_err("must error on 500");
    let msg = format!("{err:?}");
    assert!(
        msg.contains("gmx00020"),
        "error must mention failing pathway gmx00020: {msg}"
    );
}

#[tokio::test]
async fn fetch_species_throttles_at_least_50ms_between_details() {
    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/pathway/gmx"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_list_body()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00010"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_detail_with_compounds()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00020"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_detail_no_compounds()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00030"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_detail_with_compounds()))
        .mount(&server)
        .await;

    let client = build_client(&server);
    let (tx, _rx) = mpsc::channel::<KeggProgress>(64);

    let start = Instant::now();
    client
        .fetch_species_pathways("gmx", tx)
        .await
        .expect("fetch species");
    let elapsed = start.elapsed();

    // 3 pathways: 2 inter-request sleeps of 50 ms each = ≥ 100 ms.
    assert!(
        elapsed.as_millis() >= 90,
        "expected ≥ ~100 ms total elapsed for 3 pathways with 50 ms throttle, got {elapsed:?}"
    );
}

#[tokio::test]
async fn fetch_species_continues_when_progress_channel_closed() {
    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/pathway/gmx"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_list_body()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00010"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_detail_with_compounds()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00020"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_detail_no_compounds()))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/gmx00030"))
        .respond_with(ResponseTemplate::new(200).set_body_string(pathway_detail_with_compounds()))
        .mount(&server)
        .await;

    let client = build_client(&server);
    let (tx, rx) = mpsc::channel::<KeggProgress>(64);
    drop(rx); // close receiver immediately

    let species = client
        .fetch_species_pathways("gmx", tx)
        .await
        .expect("fetch must still complete when progress channel is closed");
    assert_eq!(species.pathways.len(), 3);
}
