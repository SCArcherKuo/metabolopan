//! Wiremock integration tests for the KEGG /conv/compound/pubchem client
//! method and the high-level resolver.

use metabolopan::kegg::{ConvProgress, KeggClient, resolve_cids_to_cpds};
use reqwest::Url;
use std::sync::OnceLock;
use std::sync::mpsc;
use tempfile::tempdir;
use wiremock::matchers::{method, path};
use wiremock::{Mock, MockServer, ResponseTemplate};

async fn serial() -> tokio::sync::MutexGuard<'static, ()> {
    static M: OnceLock<tokio::sync::Mutex<()>> = OnceLock::new();
    M.get_or_init(|| tokio::sync::Mutex::new(())).lock().await
}

fn setup_tmp_root() -> tempfile::TempDir {
    let dir = tempdir().expect("tempdir");
    unsafe {
        std::env::set_var("KEGG_CACHE_DIR", dir.path());
    }
    dir
}

#[tokio::test]
async fn client_returns_partial_matches_with_none_for_missing() {
    let server = MockServer::start().await;
    let body = "pubchem:5793\tcpd:C00031\npubchem:99999\tcpd:C00074\n";
    Mock::given(method("GET"))
        .and(path(
            "/conv/compound/pubchem:5793+pubchem:12345+pubchem:99999",
        ))
        .respond_with(ResponseTemplate::new(200).set_body_string(body))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).unwrap());
    let inputs = vec!["5793".to_string(), "12345".to_string(), "99999".to_string()];
    let result = client.conv_compound_pubchem(&inputs).await.expect("ok");

    assert_eq!(result.len(), 3);
    assert_eq!(result[0].1.as_deref(), Some("C00031"));
    assert_eq!(result[1].1, None);
    assert_eq!(result[2].1.as_deref(), Some("C00074"));
}

#[tokio::test]
async fn client_404_treated_as_all_no_match() {
    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/conv/compound/pubchem:A+pubchem:B"))
        .respond_with(ResponseTemplate::new(404))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).unwrap());
    let inputs = vec!["A".to_string(), "B".to_string()];
    let result = client.conv_compound_pubchem(&inputs).await.expect("404 ok");

    assert_eq!(result.len(), 2);
    assert!(result[0].1.is_none());
    assert!(result[1].1.is_none());
}

#[tokio::test]
async fn client_retries_after_503() {
    // Short backoffs so the test runs quickly. Tests also serialize on
    // these env vars below.
    let _g = serial().await;
    unsafe {
        std::env::set_var("KEGG_CONV_NETWORK_BACKOFF_MS", "10");
        std::env::set_var("KEGG_CONV_403_BACKOFF_MS", "10");
    }

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/conv/compound/pubchem:X"))
        .respond_with(ResponseTemplate::new(503))
        .up_to_n_times(1)
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/conv/compound/pubchem:X"))
        .respond_with(ResponseTemplate::new(200).set_body_string("pubchem:X\tcpd:C12345\n"))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).unwrap());
    let result = client
        .conv_compound_pubchem(&["X".to_string()])
        .await
        .expect("retry succeeded");
    assert_eq!(result[0].1.as_deref(), Some("C12345"));
}

#[tokio::test]
async fn client_retries_after_403_rate_limit() {
    // 403 is KEGG's rate-limit signal — MUST be retried, not surfaced
    // as fatal.
    let _g = serial().await;
    unsafe {
        std::env::set_var("KEGG_CONV_403_BACKOFF_MS", "10");
        std::env::set_var("KEGG_CONV_NETWORK_BACKOFF_MS", "10");
    }

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/conv/compound/pubchem:Y"))
        .respond_with(ResponseTemplate::new(403))
        .up_to_n_times(2)
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/conv/compound/pubchem:Y"))
        .respond_with(ResponseTemplate::new(200).set_body_string("pubchem:Y\tcpd:C00001\n"))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).unwrap());
    let result = client
        .conv_compound_pubchem(&["Y".to_string()])
        .await
        .expect("retry should recover from 403");
    assert_eq!(result[0].1.as_deref(), Some("C00001"));
}

#[tokio::test]
async fn client_persistent_5xx_surfaces_error_after_max_attempts() {
    let _g = serial().await;
    unsafe {
        std::env::set_var("KEGG_CONV_NETWORK_BACKOFF_MS", "10");
        std::env::set_var("KEGG_CONV_403_BACKOFF_MS", "10");
    }

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/conv/compound/pubchem:X"))
        .respond_with(ResponseTemplate::new(503))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).unwrap());
    let err = client
        .conv_compound_pubchem(&["X".to_string()])
        .await
        .expect_err("persistent 503 should fail");
    assert!(format!("{err}").contains("503"));
}

#[tokio::test]
async fn resolver_warm_cache_skips_network() {
    let _g = serial().await;
    let _tmp = setup_tmp_root();

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/conv/compound/pubchem:5793"))
        .respond_with(ResponseTemplate::new(200).set_body_string("pubchem:5793\tcpd:C00031\n"))
        .expect(1)
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).unwrap());
    let inputs = vec!["5793".to_string()];

    // First call hits server, populates cache.
    resolve_cids_to_cpds(&client, &inputs, false, None)
        .await
        .expect("cold ok");
    // Second call: mock expects exactly 1, so this MUST come from cache.
    let result = resolve_cids_to_cpds(&client, &inputs, false, None)
        .await
        .expect("warm ok");
    assert_eq!(result["5793"].as_deref(), Some("C00031"));
}

#[tokio::test]
async fn resolver_progress_carries_total_inputs() {
    let _g = serial().await;
    let _tmp = setup_tmp_root();

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/conv/compound/pubchem:A+pubchem:B"))
        .respond_with(ResponseTemplate::new(200).set_body_string("pubchem:A\tcpd:C00001\n"))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).unwrap());
    let (tx, rx) = mpsc::channel::<ConvProgress>();
    let inputs = vec!["A".to_string(), "B".to_string()];
    resolve_cids_to_cpds(&client, &inputs, false, Some(tx))
        .await
        .expect("ok");

    let events: Vec<ConvProgress> = rx.try_iter().collect();
    assert!(!events.is_empty());
    for e in &events {
        assert_eq!(e.total_inputs, 2);
    }
}
