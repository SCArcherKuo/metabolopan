use chrono::Utc;

use metabolopan::kegg::cache::{self, set_cache_root_for_tests};
use metabolopan::kegg::types::{
    KeggCacheScope, KeggCompoundSet, KeggOrganism, OrganismsCache, SpeciesKegg,
};

// All cache tests in this binary share one OnceLock-guarded cache root, so we
// run them serially using a single fixture path. Different binaries (different
// `tests/*.rs` files) get their own process and thus their own OnceLock.

static INIT: std::sync::Once = std::sync::Once::new();

fn ensure_cache_root() -> std::path::PathBuf {
    use std::sync::Mutex;
    static DIR: Mutex<Option<std::path::PathBuf>> = Mutex::new(None);
    INIT.call_once(|| {
        let tmp = tempfile::tempdir().expect("tempdir for cache root");
        let root = tmp.path().to_path_buf();
        // Leak the TempDir so it lives for the entire test-binary process and is
        // cleaned up by the OS at exit.
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

fn fresh_species(code: &str) -> SpeciesKegg {
    SpeciesKegg {
        code: code.to_string(),
        fetched_at: Utc::now(),
        pathways: vec![
            KeggCompoundSet {
                id: format!("{code}00010"),
                name: "Glycolysis".into(),
                compounds: vec!["C00001".into(), "C00002".into()],
            },
            KeggCompoundSet {
                id: format!("{code}00020"),
                name: "Citrate cycle".into(),
                compounds: vec![],
            },
        ],
    }
}

#[test]
fn species_cache_round_trip() {
    ensure_cache_root();
    let species = fresh_species("rtspc");
    cache::write_species(&species).expect("write");
    let back = cache::read_species("rtspc").expect("read").expect("Some");
    assert_eq!(back, species);
}

#[test]
fn invalidate_removes_existing_and_is_noop_for_missing() {
    ensure_cache_root();
    let species = fresh_species("inv");
    cache::write_species(&species).expect("write");
    assert!(cache::read_species("inv").unwrap().is_some());

    cache::invalidate_cache(KeggCacheScope::Species("inv".into())).expect("invalidate");
    assert!(cache::read_species("inv").unwrap().is_none());

    // Invalidating again is a no-op.
    cache::invalidate_cache(KeggCacheScope::Species("inv".into())).expect("invalidate-2");
}

#[test]
fn organisms_cache_round_trip() {
    ensure_cache_root();
    let orig = OrganismsCache {
        fetched_at: Utc::now(),
        organisms: vec![
            KeggOrganism {
                code: "gmx".into(),
                t_number: "T01710".into(),
                name: "Glycine max (soybean)".into(),
                lineage: "Eukaryota;Plants;Fabales".into(),
            },
            KeggOrganism {
                code: "sly".into(),
                t_number: "T01791".into(),
                name: "Solanum lycopersicum (tomato)".into(),
                lineage: "Eukaryota;Plants;Solanales".into(),
            },
        ],
    };
    cache::write_organisms(&orig).expect("write");
    let back = cache::read_organisms().expect("read").expect("Some");
    assert_eq!(back.organisms.len(), 2);
    assert_eq!(back.organisms[0].code, "gmx");
}
