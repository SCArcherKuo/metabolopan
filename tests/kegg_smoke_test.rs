use tokio::sync::mpsc;

use metabolopan::kegg::{KeggClient, KeggProgress};

/// Hits the real `https://rest.kegg.jp` endpoints for `gmx` (Glycine max).
/// Marked `#[ignore]` so it never runs in `cargo test`; run locally with
/// `cargo test -- --ignored kegg_smoke`. Takes ~10s of wall-time (130 pathways
/// × 50 ms throttle + per-request latency).
#[tokio::test]
#[ignore]
async fn smoke_fetch_gmx_against_real_kegg() {
    let client = KeggClient::new();
    let (tx, _rx) = mpsc::channel::<KeggProgress>(256);
    let species = client
        .fetch_species_pathways("gmx", tx)
        .await
        .expect("real KEGG fetch must succeed");
    assert!(
        species.pathways.len() > 50,
        "expected > 50 pathways, got {}",
        species.pathways.len()
    );
    let any_with_c00001 = species
        .pathways
        .iter()
        .any(|p| p.compounds.iter().any(|c| c == "C00001"));
    assert!(
        any_with_c00001,
        "expected at least one pathway to contain C00001 (Water)"
    );
}
