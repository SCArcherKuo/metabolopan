//! Integration tests for the PubChem cache-first resolver.

use metabolopan::pubchem::{PubchemClient, PubchemProgress, resolve_inchikeys_to_cids};
use std::sync::OnceLock;
use std::sync::mpsc;
use tempfile::tempdir;
use wiremock::matchers::{method, path};
use wiremock::{Mock, MockServer, ResponseTemplate};

const ENDPOINT: &str = "/compound/inchikey/property/InChIKey/CSV";

// All tests in this file share the same `set_cache_root_for_tests`
// override (which is a process-wide OnceLock). They also race over the
// PUBCHEM_CACHE_DIR env var. Serialise to keep them deterministic.
async fn serial() -> tokio::sync::MutexGuard<'static, ()> {
    static M: OnceLock<tokio::sync::Mutex<()>> = OnceLock::new();
    M.get_or_init(|| tokio::sync::Mutex::new(())).lock().await
}

fn setup_tmp_root() -> tempfile::TempDir {
    let dir = tempdir().expect("tempdir");
    // Use the env var (read on every cache_dir() call) so each test gets
    // its own dir. We deliberately do NOT call set_cache_root_for_tests
    // because that uses a OnceLock — the first test would pin the path
    // for the rest of the binary's life, and subsequent tests' tempdirs
    // would be ignored.
    unsafe {
        std::env::set_var("PUBCHEM_CACHE_DIR", dir.path());
    }
    dir
}

#[tokio::test]
async fn cold_cache_fetches_all() {
    let _g = serial().await;
    let _tmp = setup_tmp_root();

    let server = MockServer::start().await;
    let body = "\"InChIKey\",\"CID\"\n\"K1\",\"100\"\n\"K2\",\"200\"\n";
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(ResponseTemplate::new(200).set_body_string(body))
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let inputs = vec!["K1".to_string(), "K2".to_string(), "K3".to_string()];
    let result = resolve_inchikeys_to_cids(&client, &inputs, false, None)
        .await
        .expect("ok");

    assert_eq!(result.len(), 3);
    assert_eq!(result["K1"], vec!["100".to_string()]);
    assert_eq!(result["K2"], vec!["200".to_string()]);
    assert!(result["K3"].is_empty(), "K3 should be cached negative");
}

#[tokio::test]
async fn warm_cache_skips_network() {
    let _g = serial().await;
    let _tmp = setup_tmp_root();

    let server = MockServer::start().await;
    let body = "\"InChIKey\",\"CID\"\n\"K1\",\"100\"\n";
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(ResponseTemplate::new(200).set_body_string(body))
        .expect(1) // First call fills the cache.
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let inputs = vec!["K1".to_string()];

    // First call: cold, hits server once.
    resolve_inchikeys_to_cids(&client, &inputs, false, None)
        .await
        .expect("ok cold");

    // Second call: warm, should NOT hit server (mock expects exactly 1).
    let result = resolve_inchikeys_to_cids(&client, &inputs, false, None)
        .await
        .expect("ok warm");
    assert_eq!(result["K1"], vec!["100".to_string()]);
}

#[tokio::test]
async fn mixed_cache_fetches_only_missing() {
    let _g = serial().await;
    let _tmp = setup_tmp_root();

    let server = MockServer::start().await;
    // First fetch: only K1.
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(
            ResponseTemplate::new(200).set_body_string("\"InChIKey\",\"CID\"\n\"K1\",\"100\"\n"),
        )
        .up_to_n_times(1)
        .mount(&server)
        .await;
    // Second fetch: K2 only (K1 cached now).
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(
            ResponseTemplate::new(200).set_body_string("\"InChIKey\",\"CID\"\n\"K2\",\"200\"\n"),
        )
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    resolve_inchikeys_to_cids(&client, &["K1".to_string()], false, None)
        .await
        .expect("ok 1");

    let result =
        resolve_inchikeys_to_cids(&client, &["K1".to_string(), "K2".to_string()], false, None)
            .await
            .expect("ok 2");
    assert_eq!(result["K1"], vec!["100".to_string()]);
    assert_eq!(result["K2"], vec!["200".to_string()]);
}

#[tokio::test]
async fn force_refresh_refetches_everything() {
    let _g = serial().await;
    let _tmp = setup_tmp_root();

    let server = MockServer::start().await;
    // Both calls return the same body — we just want to verify two hits.
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(
            ResponseTemplate::new(200).set_body_string("\"InChIKey\",\"CID\"\n\"K1\",\"100\"\n"),
        )
        .expect(2)
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let inputs = vec!["K1".to_string()];

    // Cold: fetches.
    resolve_inchikeys_to_cids(&client, &inputs, false, None)
        .await
        .expect("ok cold");
    // Force refresh: fetches again.
    resolve_inchikeys_to_cids(&client, &inputs, true, None)
        .await
        .expect("ok refresh");
}

#[tokio::test]
async fn progress_events_carry_total_inputs() {
    let _g = serial().await;
    let _tmp = setup_tmp_root();

    let server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(
            ResponseTemplate::new(200).set_body_string("\"InChIKey\",\"CID\"\n\"K1\",\"100\"\n"),
        )
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let (tx, rx) = mpsc::channel::<PubchemProgress>();
    let inputs = vec!["K1".to_string(), "K2".to_string()];
    resolve_inchikeys_to_cids(&client, &inputs, false, Some(tx))
        .await
        .expect("ok");

    let events: Vec<PubchemProgress> = rx.try_iter().collect();
    assert!(!events.is_empty(), "expected at least one progress event");
    for e in &events {
        assert_eq!(e.total_inputs, 2);
    }
    // Final event must show all inputs accounted for.
    let last = events.last().unwrap();
    assert_eq!(last.from_cache + last.fetched, 2);
}
