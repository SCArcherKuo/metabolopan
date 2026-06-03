//! Persistent per-entry-timestamped cache for InChIKey → Vec<CID> lookups,
//! with an advisory lock file guarding writes against concurrent app
//! instances.
//!
//! Cache location resolution mirrors `kegg::cache::cache_dir()`:
//! 1. `set_cache_root_for_tests` override (tests only).
//! 2. `PUBCHEM_CACHE_DIR` environment variable (full path).
//! 3. `<data_dir>/metabolopan/cache/pubchem` (via `dirs::data_dir()`).
//! 4. `./data/cache/pubchem` fallback (when `data_dir()` is unavailable).
//!
//! Negative answers (`cids: []`) are persisted just like positive ones so
//! we never re-query a known no-match. Missing keys in the JSON map mean
//! "not yet queried".

use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::sync::OnceLock;

use crate::cache_io;
use crate::pubchem::types::InchikeyCidsEntry;

static CACHE_ROOT_OVERRIDE: OnceLock<PathBuf> = OnceLock::new();

/// Override the cache root for tests. Idempotent: first call wins.
pub fn set_cache_root_for_tests(path: PathBuf) {
    let _ = CACHE_ROOT_OVERRIDE.set(path);
}

/// Returns the directory where the PubChem cache file lives.
pub fn cache_dir() -> PathBuf {
    cache_io::resolve_cache_dir("PUBCHEM_CACHE_DIR", "pubchem", &CACHE_ROOT_OVERRIDE)
}

fn cache_path() -> PathBuf {
    cache_dir().join("inchikey.json")
}

fn lock_path() -> PathBuf {
    cache_dir().join(".inchikey.lock")
}

fn ensure_cache_dir() -> Result<()> {
    let dir = cache_dir();
    fs::create_dir_all(&dir)
        .with_context(|| format!("failed to create PubChem cache dir at {}", dir.display()))?;
    Ok(())
}

/// Read the full cache map from disk. Returns an empty map if the file
/// does not exist (cold cache).
pub fn read_cache() -> Result<HashMap<String, InchikeyCidsEntry>> {
    cache_io::load_json(&cache_path(), HashMap::new)
}

/// Write the full cache map to disk, guarded by a short-lived advisory lock.
pub fn write_cache(cache: &HashMap<String, InchikeyCidsEntry>) -> Result<()> {
    ensure_cache_dir()?;
    let bytes =
        serde_json::to_vec_pretty(cache).context("failed to serialise PubChem cache to JSON")?;
    cache_io::with_write_lock(&lock_path(), |_| {
        cache_io::atomic_write(&cache_path(), &bytes)
    })
}

/// Remove any orphaned `.inchikey.lock` file. Called at application
/// startup to recover from a crashed prior process. Unconditionally
/// removes the lock; never errors if the lock is absent.
pub fn clear_stale_locks() -> Result<()> {
    let path = lock_path();
    if path.exists() {
        fs::remove_file(&path).with_context(|| {
            format!(
                "failed to remove stale PubChem cache lock {}",
                path.display()
            )
        })?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::Utc;
    use std::sync::Mutex;
    use tempfile::tempdir;

    // Tests share process-wide state (the cache root override and the env
    // var); serialise them.
    static SERIAL: Mutex<()> = Mutex::new(());

    fn setup_tmp_root() -> tempfile::TempDir {
        let dir = tempdir().expect("tempdir");
        unsafe {
            std::env::set_var("PUBCHEM_CACHE_DIR", dir.path());
        }
        dir
    }

    #[test]
    fn roundtrip_positive_and_negative_entries() {
        let _g = SERIAL.lock().unwrap();
        let _tmp = setup_tmp_root();

        let mut cache = HashMap::new();
        cache.insert(
            "INCHIKEY-A".to_string(),
            InchikeyCidsEntry {
                cids: vec!["5793".into()],
                fetched_at: Utc::now(),
            },
        );
        cache.insert(
            "INCHIKEY-B".to_string(),
            InchikeyCidsEntry {
                cids: vec![],
                fetched_at: Utc::now(),
            },
        );
        write_cache(&cache).expect("write");
        let read = read_cache().expect("read");
        assert_eq!(read.len(), 2);
        assert_eq!(
            read.get("INCHIKEY-A").unwrap().cids,
            vec!["5793".to_string()]
        );
        assert!(read.get("INCHIKEY-B").unwrap().cids.is_empty());
    }

    #[test]
    fn missing_key_vs_empty_entry() {
        let _g = SERIAL.lock().unwrap();
        let _tmp = setup_tmp_root();

        let mut cache = HashMap::new();
        cache.insert(
            "INCHIKEY-PRESENT".to_string(),
            InchikeyCidsEntry {
                cids: vec![],
                fetched_at: Utc::now(),
            },
        );
        write_cache(&cache).expect("write");
        let read = read_cache().expect("read");

        // Cache hit, negative answer: present with cids = [].
        let present = read.get("INCHIKEY-PRESENT");
        assert!(present.is_some());
        assert!(present.unwrap().cids.is_empty());

        // Cache miss: absent from map.
        assert!(!read.contains_key("INCHIKEY-MISSING"));
    }

    #[test]
    fn clear_stale_locks_removes_orphan() {
        let _g = SERIAL.lock().unwrap();
        let _tmp = setup_tmp_root();
        ensure_cache_dir().unwrap();

        // Simulate a crashed prior writer.
        let lp = lock_path();
        fs::File::create(&lp).expect("create lock");
        assert!(lp.exists());

        clear_stale_locks().expect("clear");
        assert!(!lp.exists(), "lock should be cleared");
    }

    #[test]
    fn clear_stale_locks_is_noop_when_absent() {
        let _g = SERIAL.lock().unwrap();
        let _tmp = setup_tmp_root();
        ensure_cache_dir().unwrap();
        assert!(!lock_path().exists());
        clear_stale_locks().expect("clear");
    }

    #[test]
    fn write_releases_lock_on_success() {
        let _g = SERIAL.lock().unwrap();
        let _tmp = setup_tmp_root();
        let cache = HashMap::new();
        write_cache(&cache).expect("write");
        assert!(!lock_path().exists(), "lock should be released after write");
    }
}
