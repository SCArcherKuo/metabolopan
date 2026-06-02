//! Integration tests for Track B: KEGG module REST + cache pipeline.
//!
//! Uses `wiremock` to stub `/list/module` and `/get/<module-id>`,
//! exercising `fetch_modules` end-to-end: cold full fetch, warm cache
//! short-circuit, force-refresh prune semantics, 403 retry on per-module
//! GET, and empty COMPOUND / COMPLETE handling.

use metabolopan::kegg::cache::{self, set_cache_root_for_tests};
use metabolopan::kegg::types::{
    KeggCacheScope, KeggModuleEntry, KeggModulesCache, ModuleFetchProgress,
};
use metabolopan::kegg::{KeggClient, fetch_modules};

use chrono::Utc;
use reqwest::Url;
use std::collections::HashMap;
use std::sync::{Once, OnceLock};
use tokio::sync::mpsc;
use wiremock::matchers::{method, path, path_regex};
use wiremock::{Mock, MockServer, ResponseTemplate};

static INIT: Once = Once::new();

async fn serial() -> tokio::sync::MutexGuard<'static, ()> {
    static M: OnceLock<tokio::sync::Mutex<()>> = OnceLock::new();
    M.get_or_init(|| tokio::sync::Mutex::new(())).lock().await
}

fn ensure_cache_root() {
    INIT.call_once(|| {
        let tmp = tempfile::tempdir().expect("tempdir for cache root");
        let root = tmp.path().to_path_buf();
        std::mem::forget(tmp);
        set_cache_root_for_tests(root);
    });
}

fn reset_modules_cache() {
    let _ = cache::invalidate_cache(KeggCacheScope::Modules);
}

fn list_response(ids: &[(&str, &str)]) -> String {
    ids.iter()
        .map(|(id, name)| format!("{id}\t{name}\n"))
        .collect()
}

fn module_detail(name: &str, compounds: &[&str], orgs: &[&str]) -> String {
    let mut s = format!(
        "ENTRY       FOO            Pathway   Module\nNAME        {name}\nDEFINITION  (K00000)\n"
    );
    if !compounds.is_empty() {
        s.push_str("COMPOUND    ");
        for (i, c) in compounds.iter().enumerate() {
            if i == 0 {
                s.push_str(&format!("{c}  Compound {c}\n"));
            } else {
                s.push_str(&format!("            {c}  Compound {c}\n"));
            }
        }
    }
    if !orgs.is_empty() {
        s.push_str("COMPLETE    ");
        for (i, o) in orgs.iter().enumerate() {
            if i == 0 {
                s.push_str(&format!("{o}  Some species (common)\n"));
            } else {
                s.push_str(&format!("            {o}  Some species (common)\n"));
            }
        }
    }
    s.push_str("REFERENCE   foo\n");
    s
}

#[tokio::test]
async fn cold_full_fetch_populates_cache() {
    let _g = serial().await;
    ensure_cache_root();
    reset_modules_cache();

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/module"))
        .respond_with(ResponseTemplate::new(200).set_body_string(list_response(&[
            ("M00001", "Mod one"),
            ("M00002", "Mod two"),
        ])))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/M00001"))
        .respond_with(ResponseTemplate::new(200).set_body_string(module_detail(
            "Mod one",
            &["C00001", "C00002"],
            &["hsa", "mmu"],
        )))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/M00002"))
        .respond_with(ResponseTemplate::new(200).set_body_string(module_detail(
            "Mod two",
            &["C00003"],
            &["ath"],
        )))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).expect("uri"));
    let (tx, _rx) = mpsc::channel::<ModuleFetchProgress>(16);
    let cache_out = fetch_modules(&client, false, tx).await.expect("fetch ok");

    assert_eq!(cache_out.modules.len(), 2);
    let m1 = &cache_out.modules["M00001"];
    assert_eq!(m1.name, "Mod one");
    assert_eq!(m1.compounds, vec!["C00001", "C00002"]);
    assert!(m1.complete_orgs.contains("hsa"));
    assert!(m1.complete_orgs.contains("mmu"));
    let m2 = &cache_out.modules["M00002"];
    assert_eq!(m2.compounds, vec!["C00003"]);
    assert!(m2.complete_orgs.contains("ath"));

    // Persisted to disk.
    let read = cache::read_modules_cache().expect("read");
    assert_eq!(read.modules.len(), 2);
}

#[tokio::test]
async fn warm_cache_fetches_only_missing() {
    let _g = serial().await;
    ensure_cache_root();
    reset_modules_cache();

    // Pre-populate cache with M00001.
    let mut seed = KeggModulesCache::default();
    seed.modules.insert(
        "M00001".into(),
        KeggModuleEntry {
            name: "Existing".into(),
            compounds: vec!["C00001".into()],
            complete_orgs: ["hsa".to_string()].into_iter().collect(),
            fetched_at: Utc::now(),
        },
    );
    cache::write_modules_cache(&seed).expect("seed");

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/module"))
        .respond_with(ResponseTemplate::new(200).set_body_string(list_response(&[
            ("M00001", "Mod one"),
            ("M00002", "Mod two (new)"),
        ])))
        .mount(&server)
        .await;
    // Only M00002 should be fetched. Register a strict mock that
    // verifies M00001 detail is NOT requested.
    Mock::given(method("GET"))
        .and(path("/get/M00002"))
        .respond_with(ResponseTemplate::new(200).set_body_string(module_detail(
            "Mod two (new)",
            &["C00099"],
            &["mmu"],
        )))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).expect("uri"));
    let (tx, _rx) = mpsc::channel::<ModuleFetchProgress>(16);
    let cache_out = fetch_modules(&client, false, tx).await.expect("fetch ok");

    assert_eq!(cache_out.modules.len(), 2);
    // M00001 untouched (preserved existing data, NOT re-fetched).
    assert_eq!(cache_out.modules["M00001"].name, "Existing");
    // M00002 freshly fetched.
    assert_eq!(cache_out.modules["M00002"].name, "Mod two (new)");

    // Verify M00001 detail was NOT requested.
    let received = server.received_requests().await.expect("requests");
    let m1_hits = received
        .iter()
        .filter(|r| r.url.path() == "/get/M00001")
        .count();
    assert_eq!(m1_hits, 0, "warm cache must NOT re-fetch existing modules");
}

#[tokio::test]
async fn force_refresh_prunes_removed_modules() {
    let _g = serial().await;
    ensure_cache_root();
    reset_modules_cache();

    // Pre-populate with M00001 and M99999 (which is no longer in /list).
    let mut seed = KeggModulesCache::default();
    for id in ["M00001", "M99999"] {
        seed.modules.insert(
            id.into(),
            KeggModuleEntry {
                name: format!("Old {id}"),
                compounds: vec!["C00001".into()],
                complete_orgs: HashMap::<String, ()>::new().keys().cloned().collect(),
                fetched_at: Utc::now(),
            },
        );
    }
    cache::write_modules_cache(&seed).expect("seed");

    let server = MockServer::start().await;
    // /list returns only M00001 — M99999 was retired by KEGG.
    Mock::given(method("GET"))
        .and(path("/list/module"))
        .respond_with(
            ResponseTemplate::new(200).set_body_string(list_response(&[("M00001", "Mod one")])),
        )
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/M00001"))
        .respond_with(ResponseTemplate::new(200).set_body_string(module_detail(
            "Mod one (refreshed)",
            &["C00001"],
            &["hsa"],
        )))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).expect("uri"));
    let (tx, _rx) = mpsc::channel::<ModuleFetchProgress>(16);
    let cache_out = fetch_modules(&client, true, tx).await.expect("fetch ok");

    // M99999 pruned; only M00001 remains.
    assert_eq!(cache_out.modules.len(), 1);
    assert!(cache_out.modules.contains_key("M00001"));
    assert!(!cache_out.modules.contains_key("M99999"));
    assert_eq!(cache_out.modules["M00001"].name, "Mod one (refreshed)");
}

#[tokio::test]
async fn module_get_403_is_retried() {
    // KEGG's 403 = rate-limit; the client retries with backoff. We
    // override the backoff env var to keep the test fast.
    let _g = serial().await;
    ensure_cache_root();
    reset_modules_cache();
    // SAFETY: This test is serialised by the SERIAL Mutex; no concurrent
    // env access.
    unsafe {
        std::env::set_var("KEGG_CONV_403_BACKOFF_MS", "10");
        std::env::set_var("KEGG_CONV_NETWORK_BACKOFF_MS", "10");
    }

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/module"))
        .respond_with(
            ResponseTemplate::new(200).set_body_string(list_response(&[("M00001", "M1")])),
        )
        .mount(&server)
        .await;
    // First request 403; subsequent requests 200.
    Mock::given(method("GET"))
        .and(path("/get/M00001"))
        .respond_with(ResponseTemplate::new(403))
        .up_to_n_times(1)
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path_regex(r"^/get/M00001$"))
        .respond_with(ResponseTemplate::new(200).set_body_string(module_detail(
            "M1",
            &["C00001"],
            &["hsa"],
        )))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).expect("uri"));
    let (tx, _rx) = mpsc::channel::<ModuleFetchProgress>(16);
    let cache_out = fetch_modules(&client, false, tx).await.expect("fetch ok");
    assert_eq!(cache_out.modules.len(), 1);
    assert_eq!(cache_out.modules["M00001"].name, "M1");
}

#[tokio::test]
async fn complete_cache_skips_per_module_fetches_after_listing() {
    // Track G: warm-cache path where /list/module returns the same IDs
    // already present in the cache. fetch_modules should hit /list once,
    // determine missing_ids is empty, emit a single "(cache complete)"
    // initial progress event, and return without any /get/<module>
    // requests.
    let _g = serial().await;
    ensure_cache_root();
    reset_modules_cache();

    let mut seed = KeggModulesCache::default();
    seed.modules.insert(
        "M00001".into(),
        KeggModuleEntry {
            name: "Cached M1".into(),
            compounds: vec!["C00001".into()],
            complete_orgs: ["hsa".to_string()].into_iter().collect(),
            fetched_at: Utc::now(),
        },
    );
    seed.modules.insert(
        "M00002".into(),
        KeggModuleEntry {
            name: "Cached M2".into(),
            compounds: vec!["C00002".into()],
            complete_orgs: ["mmu".to_string()].into_iter().collect(),
            fetched_at: Utc::now(),
        },
    );
    cache::write_modules_cache(&seed).expect("seed");

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/module"))
        .respond_with(
            ResponseTemplate::new(200)
                .set_body_string(list_response(&[("M00001", "M1"), ("M00002", "M2")])),
        )
        .mount(&server)
        .await;
    // No /get mocks registered — if fetch_modules tries to hit /get the
    // mock server will 404 and the assertion below will catch the
    // unexpected request.

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).expect("uri"));
    let (tx, mut rx) = mpsc::channel::<ModuleFetchProgress>(16);
    let cache_out = fetch_modules(&client, false, tx).await.expect("fetch ok");
    assert_eq!(cache_out.modules.len(), 2);

    // First (and only) progress event: total = 0 (no missing), "(cache complete)".
    let first = rx.try_recv().expect("initial progress event");
    assert_eq!(first.total, 0);
    assert_eq!(first.completed, 0);
    assert!(first.current_id.contains("cache complete"));
    // No further events.
    assert!(rx.try_recv().is_err());

    // Verify zero /get/<module> requests were made.
    let received = server.received_requests().await.expect("requests");
    let get_hits = received
        .iter()
        .filter(|r| r.url.path().starts_with("/get/"))
        .count();
    assert_eq!(
        get_hits, 0,
        "warm cache must NOT hit /get/<module> when nothing is missing"
    );
}

#[tokio::test]
async fn initial_progress_event_emits_missing_count() {
    // Track G: when fetch_modules has work to do, the first emitted
    // progress event should carry `total = missing_ids.len()` so the UI
    // shows an accurate progress denominator immediately instead of the
    // misleading "0 / 0" placeholder.
    let _g = serial().await;
    ensure_cache_root();
    reset_modules_cache();

    // Pre-seed cache with M00001 only; /list returns M00001 + M00002.
    let mut seed = KeggModulesCache::default();
    seed.modules.insert(
        "M00001".into(),
        KeggModuleEntry {
            name: "Cached".into(),
            compounds: vec!["C00001".into()],
            complete_orgs: ["hsa".to_string()].into_iter().collect(),
            fetched_at: Utc::now(),
        },
    );
    cache::write_modules_cache(&seed).expect("seed");

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/module"))
        .respond_with(
            ResponseTemplate::new(200)
                .set_body_string(list_response(&[("M00001", "M1"), ("M00002", "M2 new")])),
        )
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/M00002"))
        .respond_with(ResponseTemplate::new(200).set_body_string(module_detail(
            "M2 new",
            &["C00002"],
            &["hsa"],
        )))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).expect("uri"));
    let (tx, mut rx) = mpsc::channel::<ModuleFetchProgress>(16);
    let _cache_out = fetch_modules(&client, false, tx).await.expect("fetch ok");

    // First event = initial "(starting fetch)" with total = 1 (only M00002 missing).
    let first = rx.try_recv().expect("initial progress event");
    assert_eq!(first.total, 1);
    assert_eq!(first.completed, 0);
    assert!(first.current_id.contains("starting fetch"));

    // Subsequent event(s) = per-module from fetch_modules_incremental.
    // At least one more event should have current_id = "M00002".
    let mut saw_m00002 = false;
    while let Ok(ev) = rx.try_recv() {
        if ev.current_id == "M00002" {
            saw_m00002 = true;
            assert_eq!(ev.total, 1);
            assert_eq!(ev.completed, 1);
        }
    }
    assert!(saw_m00002, "expected per-module progress event for M00002");
}

#[tokio::test]
async fn empty_compound_and_complete_blocks_are_cached_normally() {
    let _g = serial().await;
    ensure_cache_root();
    reset_modules_cache();

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/module"))
        .respond_with(ResponseTemplate::new(200).set_body_string(list_response(&[
            ("M00001", "Has both"),
            ("M00002", "No compounds"),
            ("M00003", "No COMPLETE"),
        ])))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/M00001"))
        .respond_with(ResponseTemplate::new(200).set_body_string(module_detail(
            "Has both",
            &["C00001"],
            &["hsa"],
        )))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/M00002"))
        .respond_with(ResponseTemplate::new(200).set_body_string(module_detail(
            "No compounds",
            &[],
            &["hsa"],
        )))
        .mount(&server)
        .await;
    Mock::given(method("GET"))
        .and(path("/get/M00003"))
        .respond_with(ResponseTemplate::new(200).set_body_string(module_detail(
            "No COMPLETE",
            &["C00099"],
            &[],
        )))
        .mount(&server)
        .await;

    let client = KeggClient::with_base_url(Url::parse(&server.uri()).expect("uri"));
    let (tx, _rx) = mpsc::channel::<ModuleFetchProgress>(16);
    let cache_out = fetch_modules(&client, false, tx).await.expect("fetch ok");

    assert_eq!(cache_out.modules.len(), 3);
    assert!(cache_out.modules["M00002"].compounds.is_empty());
    assert!(!cache_out.modules["M00002"].complete_orgs.is_empty());
    assert!(!cache_out.modules["M00003"].compounds.is_empty());
    assert!(cache_out.modules["M00003"].complete_orgs.is_empty());
}
