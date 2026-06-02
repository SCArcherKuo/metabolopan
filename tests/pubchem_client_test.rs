//! Wiremock-driven integration tests for the PubChem PUG REST client.

use metabolopan::pubchem::PubchemClient;
use wiremock::matchers::{header, method, path};
use wiremock::{Mock, MockServer, ResponseTemplate};

const ENDPOINT: &str = "/compound/inchikey/property/InChIKey/CSV";

#[tokio::test]
async fn standard_batch_returns_each_input_once() {
    let server = MockServer::start().await;
    let body = "\"InChIKey\",\"CID\"\n\
                \"KEY-A\",\"5793\"\n\
                \"KEY-B\",\"100\"\n";
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .and(header("Content-Type", "application/x-www-form-urlencoded"))
        .respond_with(ResponseTemplate::new(200).set_body_string(body))
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let inputs = vec![
        "KEY-A".to_string(),
        "KEY-B".to_string(),
        "KEY-C".to_string(),
    ];
    let result = client
        .post_inchikeys_to_cids(&inputs)
        .await
        .expect("ok response");

    assert_eq!(result.len(), 3);
    assert_eq!(result[0].0, "KEY-A");
    assert_eq!(result[0].1, vec!["5793".to_string()]);
    assert_eq!(result[1].0, "KEY-B");
    assert_eq!(result[1].1, vec!["100".to_string()]);
    assert_eq!(result[2].0, "KEY-C");
    assert!(
        result[2].1.is_empty(),
        "KEY-C had no row in response, must be no-match"
    );
}

#[tokio::test]
async fn outbound_request_carries_crate_user_agent() {
    // The PubChem client now builds its reqwest client via the shared
    // `cache_io::http_client`, which sets `User-Agent: <crate>/<version>`
    // (the one intended behavior change of dedup-network-cache-layer).
    let server = MockServer::start().await;
    let expected_ua = concat!(env!("CARGO_PKG_NAME"), "/", env!("CARGO_PKG_VERSION"));
    // This mock matches ONLY when the User-Agent header equals the crate UA.
    // If the client omitted it, no mock matches → wiremock 404 → the client
    // treats it as all-no-match (empty cids), and the assert below fails.
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .and(header("user-agent", expected_ua))
        .respond_with(
            ResponseTemplate::new(200)
                .set_body_string("\"InChIKey\",\"CID\"\n\"KEY-A\",\"5793\"\n"),
        )
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let result = client
        .post_inchikeys_to_cids(&["KEY-A".to_string()])
        .await
        .expect("ok response");

    assert_eq!(result.len(), 1);
    assert_eq!(
        result[0].1,
        vec!["5793".to_string()],
        "request must carry User-Agent = {expected_ua} for the mock to match and return the CID"
    );
}

#[tokio::test]
async fn multi_cid_same_inchikey_is_aggregated() {
    let server = MockServer::start().await;
    let body = "\"InChIKey\",\"CID\"\n\
                \"K1\",\"100\"\n\
                \"K1\",\"200\"\n\
                \"K1\",\"300\"\n";
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(ResponseTemplate::new(200).set_body_string(body))
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let result = client
        .post_inchikeys_to_cids(&["K1".to_string()])
        .await
        .expect("ok");

    assert_eq!(result.len(), 1);
    assert_eq!(
        result[0].1,
        vec!["100".to_string(), "200".to_string(), "300".to_string()]
    );
}

#[tokio::test]
async fn empty_body_means_all_no_match() {
    let server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(ResponseTemplate::new(200).set_body_string(""))
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let inputs = vec!["K1".to_string(), "K2".to_string()];
    let result = client
        .post_inchikeys_to_cids(&inputs)
        .await
        .expect("ok empty");

    assert_eq!(result.len(), 2);
    assert!(result[0].1.is_empty());
    assert!(result[1].1.is_empty());
}

#[tokio::test]
async fn retry_once_after_503_then_succeeds() {
    let server = MockServer::start().await;
    // First call returns 503, second call returns 200 with content.
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(ResponseTemplate::new(503))
        .up_to_n_times(1)
        .mount(&server)
        .await;
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(
            ResponseTemplate::new(200).set_body_string("\"InChIKey\",\"CID\"\n\"K1\",\"5793\"\n"),
        )
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let result = client
        .post_inchikeys_to_cids(&["K1".to_string()])
        .await
        .expect("retry succeeded");

    assert_eq!(result[0].1, vec!["5793".to_string()]);
}

#[tokio::test]
async fn two_consecutive_503s_surface_error() {
    let server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(ResponseTemplate::new(503))
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let err = client
        .post_inchikeys_to_cids(&["K1".to_string()])
        .await
        .expect_err("two 503s should fail");
    let msg = format!("{err}");
    assert!(msg.contains("503"), "error must mention 503: {msg}");
}

#[tokio::test]
async fn empty_input_short_circuits() {
    let server = MockServer::start().await;
    // No mock — function must not hit the network for empty input.

    let client = PubchemClient::with_base_url(server.uri());
    let result = client.post_inchikeys_to_cids(&[]).await.expect("empty ok");
    assert!(result.is_empty());
}

#[tokio::test]
async fn http_404_is_confirmed_all_no_match_not_fatal() {
    // PubChem PUG REST returns 404 PUGREST.NotFound when the entire batch
    // matches zero CIDs. This must be semantically identical to an empty
    // 200 body — every input becomes (key, vec![]) — NOT propagated as
    // an error. Mirrors KEGG /conv's AllNoMatch handling.
    let server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(
            ResponseTemplate::new(404)
                .set_body_string("Status: 404\nCode: PUGREST.NotFound\nMessage: No CID found"),
        )
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let inputs = vec![
        "AAAAAAAAAAAAAA-AAAAAAAAAA-A".to_string(),
        "BBBBBBBBBBBBBB-BBBBBBBBBB-B".to_string(),
    ];
    let result = client
        .post_inchikeys_to_cids(&inputs)
        .await
        .expect("404 must NOT be fatal");
    assert_eq!(result.len(), 2);
    assert_eq!(result[0].0, "AAAAAAAAAAAAAA-AAAAAAAAAA-A");
    assert!(result[0].1.is_empty());
    assert_eq!(result[1].0, "BBBBBBBBBBBBBB-BBBBBBBBBB-B");
    assert!(result[1].1.is_empty());
}

#[tokio::test]
async fn http_400_remains_fatal() {
    // Only 404 is remapped to AllNoMatch; other 4xx (400 Bad Request,
    // 414 URI Too Long, etc.) still surface as Fatal so genuine bugs
    // aren't silenced.
    let server = MockServer::start().await;
    Mock::given(method("POST"))
        .and(path(ENDPOINT))
        .respond_with(ResponseTemplate::new(400))
        .mount(&server)
        .await;

    let client = PubchemClient::with_base_url(server.uri());
    let err = client
        .post_inchikeys_to_cids(&["K1".to_string()])
        .await
        .expect_err("400 must remain fatal");
    let msg = format!("{err}");
    assert!(msg.contains("400"), "error must mention 400: {msg}");
}
