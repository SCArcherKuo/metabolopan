//! Integration tests for Track C: eager organism load + Group precompute.
//!
//! Covers the `kegg::list_organisms` contract from the `kegg-fetching`
//! spec's MODIFIED + ADDED requirements:
//! - Both cache hit AND fresh fetch paths persist `organism_groups.json`.
//! - The persisted index's `fetched_at` matches the organisms cache's.
//! - The cache-first policy: existing organisms.json shortcuts the REST
//!   fetch regardless of age.

use metabolopan::kegg::cache::{self, set_cache_root_for_tests};
use metabolopan::kegg::types::{KeggOrganism, OrganismsCache};
use metabolopan::kegg::{KeggClient, build_organism_group_index, list_organisms};

use chrono::{Duration, Utc};
use std::sync::{Once, OnceLock};
use wiremock::matchers::{method, path};
use wiremock::{Mock, MockServer, ResponseTemplate};

static INIT: Once = Once::new();

/// All tests in this binary mutate the same on-disk cache files. Cargo
/// runs `#[test]` and `#[tokio::test]` cases in parallel by default, so
/// each test acquires this Mutex on entry to serialise file access.
/// Uses `tokio::sync::Mutex` so the guard can be held across `.await`
/// without tripping `clippy::await_holding_lock`.
fn shared_mutex() -> &'static tokio::sync::Mutex<()> {
    static M: OnceLock<tokio::sync::Mutex<()>> = OnceLock::new();
    M.get_or_init(|| tokio::sync::Mutex::new(()))
}

async fn serial() -> tokio::sync::MutexGuard<'static, ()> {
    shared_mutex().lock().await
}

/// Blocking lock acquisition for plain `#[test]` (non-async) functions
/// that share the same serialisation as the async tests.
/// `tokio::sync::Mutex::blocking_lock` panics if called from inside a
/// tokio runtime, which is the correct behaviour here — `#[test]`
/// functions don't run on a tokio runtime.
fn serial_blocking() -> tokio::sync::MutexGuard<'static, ()> {
    shared_mutex().blocking_lock()
}

/// Wipe both organism-related caches so a test starts from a known clean
/// state. Idempotent (safe when files are absent).
fn reset_caches() {
    use metabolopan::kegg::types::KeggCacheScope;
    let _ = cache::invalidate_cache(KeggCacheScope::Organisms);
    let _ = cache::invalidate_cache(KeggCacheScope::OrganismGroups);
}

fn ensure_cache_root() -> std::path::PathBuf {
    use std::sync::Mutex;
    static DIR: Mutex<Option<std::path::PathBuf>> = Mutex::new(None);
    INIT.call_once(|| {
        let tmp = tempfile::tempdir().expect("tempdir for cache root");
        let root = tmp.path().to_path_buf();
        std::mem::forget(tmp);
        set_cache_root_for_tests(root.clone());
        *DIR.lock().unwrap() = Some(root);
    });
    DIR.lock()
        .unwrap()
        .as_ref()
        .expect("cache root set")
        .clone()
}

fn sample_organisms() -> Vec<KeggOrganism> {
    vec![
        KeggOrganism {
            t_number: "T01001".into(),
            code: "hsa".into(),
            name: "Homo sapiens (human)".into(),
            lineage: "Eukaryotes;Animals;Mammals;Primates".into(),
        },
        KeggOrganism {
            t_number: "T01002".into(),
            code: "mmu".into(),
            name: "Mus musculus (mouse)".into(),
            lineage: "Eukaryotes;Animals;Mammals;Rodents".into(),
        },
        KeggOrganism {
            t_number: "T00041".into(),
            code: "ath".into(),
            name: "Arabidopsis thaliana".into(),
            lineage: "Eukaryotes;Plants;Eudicots;Brassicales".into(),
        },
    ]
}

#[tokio::test]
async fn cache_hit_writes_organism_group_index() {
    let _g = serial().await;
    ensure_cache_root();
    reset_caches();
    // Pre-populate organisms cache.
    let fetched_at = Utc::now() - Duration::days(3);
    let cache_entry = OrganismsCache {
        fetched_at,
        organisms: sample_organisms(),
    };
    cache::write_organisms(&cache_entry).expect("write organisms");

    // list_organisms returns the cached data AND writes the group index.
    let client = KeggClient::new();
    let organisms = list_organisms(&client).await.expect("list ok");
    assert_eq!(organisms.len(), 3);

    let index = cache::read_organism_group_index()
        .expect("read index")
        .expect("index present");
    assert_eq!(index.fetched_at, fetched_at);
    // Level 1: all 3 under Eukaryotes.
    assert_eq!(index.by_level[0]["Eukaryotes"].len(), 3);
    // Level 2: 2 Animals + 1 Plants.
    assert_eq!(index.by_level[1]["Animals"].len(), 2);
    assert_eq!(index.by_level[1]["Plants"].len(), 1);
}

#[tokio::test]
async fn cache_first_stale_cache_short_circuits_rest() {
    let _g = serial().await;
    ensure_cache_root();
    reset_caches();
    // Pre-populate with a 50-day-old cache.
    let stale_fetched_at = Utc::now() - Duration::days(50);
    let cache_entry = OrganismsCache {
        fetched_at: stale_fetched_at,
        organisms: sample_organisms(),
    };
    cache::write_organisms(&cache_entry).expect("write organisms");

    // Mock server expects ZERO requests — cache-first must NOT hit REST.
    let server = MockServer::start().await;
    // Don't register any matchers; any incoming request fails the test.

    let client = KeggClient::with_base_url(server.uri().parse().expect("uri parses"));
    let organisms = list_organisms(&client)
        .await
        .expect("list ok from stale cache");
    assert_eq!(organisms.len(), 3);

    // Verify no requests reached the mock.
    let received = server.received_requests().await.expect("requests");
    assert!(
        received.is_empty(),
        "stale cache must not trigger REST fetch (cache-first policy); received: {received:?}"
    );

    // Group index was also (re)written with the same stale timestamp.
    let index = cache::read_organism_group_index()
        .expect("read")
        .expect("present");
    assert_eq!(index.fetched_at, stale_fetched_at);
}

#[tokio::test]
async fn fresh_fetch_writes_both_caches_with_matching_timestamps() {
    let _g = serial().await;
    ensure_cache_root();
    reset_caches();

    let server = MockServer::start().await;
    Mock::given(method("GET"))
        .and(path("/list/organism"))
        .respond_with(ResponseTemplate::new(200).set_body_string(
            "T01001\thsa\tHomo sapiens (human)\tEukaryotes;Animals;Mammals;Primates\n\
             T00041\tath\tArabidopsis thaliana\tEukaryotes;Plants;Eudicots;Brassicales\n",
        ))
        .mount(&server)
        .await;

    let before = Utc::now();
    let client = KeggClient::with_base_url(server.uri().parse().expect("uri parses"));
    let organisms = list_organisms(&client).await.expect("list ok");
    let after = Utc::now();
    assert_eq!(organisms.len(), 2);

    // Both caches written.
    let orgs_cache = cache::read_organisms().expect("read").expect("present");
    let index = cache::read_organism_group_index()
        .expect("read")
        .expect("present");

    // Timestamps match each other AND lie within the bracketing window.
    assert_eq!(orgs_cache.fetched_at, index.fetched_at);
    assert!(orgs_cache.fetched_at >= before);
    assert!(orgs_cache.fetched_at <= after);

    // Index shape sanity check.
    assert_eq!(index.by_level[0]["Eukaryotes"].len(), 2);
    assert_eq!(index.by_level[1]["Animals"].len(), 1);
    assert_eq!(index.by_level[1]["Plants"].len(), 1);
}

#[test]
fn build_index_round_trips_via_persistent_cache() {
    let _g = serial_blocking();
    ensure_cache_root();
    reset_caches();
    let organisms = sample_organisms();
    let ts = Utc::now();
    let index = build_organism_group_index(&organisms, ts);
    cache::write_organism_group_index(&index).expect("write");
    let read = cache::read_organism_group_index()
        .expect("read")
        .expect("present");
    assert_eq!(read, index);
}
