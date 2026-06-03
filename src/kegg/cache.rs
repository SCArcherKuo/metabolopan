use anyhow::{Context, Result};
use std::collections::HashMap;
use std::fs;
use std::path::PathBuf;
use std::sync::OnceLock;
use std::thread;
use std::time::{Duration, Instant};
use tracing::warn;

use crate::cache_io;
use crate::kegg::types::{
    CidCpdEntry, KeggCacheScope, KeggModulesCache, OrganismGroupIndex, OrganismsCache, SpeciesKegg,
};

static CACHE_ROOT_OVERRIDE: OnceLock<PathBuf> = OnceLock::new();

/// Override the cache root for tests. Idempotent: first call wins.
pub fn set_cache_root_for_tests(path: PathBuf) {
    let _ = CACHE_ROOT_OVERRIDE.set(path);
}

/// Returns the directory where KEGG cache files are stored.
///
/// Resolution order:
/// 1. `set_cache_root_for_tests` override (tests only).
/// 2. `KEGG_CACHE_DIR` environment variable (full path; the value is used verbatim).
/// 3. `<data_dir>/metabolopan/cache/kegg`, where `<data_dir>` is
///    `dirs::data_dir()` (macOS `~/Library/Application Support`, Linux
///    `~/.local/share`, Windows `%APPDATA%`). Re-anchored off the binary so the
///    app runs from a read-only install (`.app` bundle, `/Applications`).
/// 4. Fallback if `current_exe()` is not resolvable (very rare on supported
///    platforms): `./data/cache/kegg` relative to the current working directory.
pub fn cache_dir() -> PathBuf {
    cache_io::resolve_cache_dir("KEGG_CACHE_DIR", "kegg", &CACHE_ROOT_OVERRIDE)
}

fn organisms_path() -> PathBuf {
    cache_dir().join("organisms.json")
}

fn species_path(code: &str) -> PathBuf {
    cache_dir().join(format!("{code}.json"))
}

fn modules_path() -> PathBuf {
    cache_dir().join("modules.json")
}

fn organism_groups_path() -> PathBuf {
    cache_dir().join("organism_groups.json")
}

fn ensure_cache_dir() -> Result<()> {
    let dir = cache_dir();
    fs::create_dir_all(&dir)
        .with_context(|| format!("failed to create KEGG cache dir at {}", dir.display()))?;
    Ok(())
}

pub fn read_organisms() -> Result<Option<OrganismsCache>> {
    cache_io::load_json(&organisms_path(), || None)
}

pub fn write_organisms(cache: &OrganismsCache) -> Result<()> {
    ensure_cache_dir()?;
    let bytes =
        serde_json::to_vec_pretty(cache).context("failed to serialise OrganismsCache to JSON")?;
    cache_io::atomic_write(&organisms_path(), &bytes)
}

pub fn read_species(code: &str) -> Result<Option<SpeciesKegg>> {
    cache_io::load_json(&species_path(code), || None)
}

pub fn write_species(species: &SpeciesKegg) -> Result<()> {
    ensure_cache_dir()?;
    let bytes =
        serde_json::to_vec_pretty(species).context("failed to serialise SpeciesKegg to JSON")?;
    cache_io::atomic_write(&species_path(&species.code), &bytes)
}

pub fn invalidate_cache(scope: KeggCacheScope) -> Result<()> {
    let path = match scope {
        KeggCacheScope::Organisms => organisms_path(),
        KeggCacheScope::Species(code) => species_path(&code),
        KeggCacheScope::Modules => modules_path(),
        KeggCacheScope::OrganismGroups => organism_groups_path(),
    };
    if path.exists() {
        fs::remove_file(&path)
            .with_context(|| format!("failed to remove cache file {}", path.display()))?;
    }
    Ok(())
}

pub fn read_organism_group_index() -> Result<Option<OrganismGroupIndex>> {
    cache_io::load_json(&organism_groups_path(), || None)
}

pub fn write_organism_group_index(index: &OrganismGroupIndex) -> Result<()> {
    ensure_cache_dir()?;
    let bytes = serde_json::to_vec_pretty(index)
        .context("failed to serialise OrganismGroupIndex to JSON")?;
    cache_io::atomic_write(&organism_groups_path(), &bytes)
}

// ─── CID → cpd cache (Stage 3) ────────────────────────────────────────────
//
// The CID-to-cpd cache lives at `<cache_dir>/cid_to_cpd.json` and grows
// incrementally across Stage 3 sessions. Writes are guarded by an advisory
// `.cid_to_cpd.lock` file; see the kegg-fetching capability spec for the
// contract.

fn cid_to_cpd_path() -> PathBuf {
    cache_dir().join("cid_to_cpd.json")
}

fn cid_to_cpd_lock_path() -> PathBuf {
    cache_dir().join(".cid_to_cpd.lock")
}

pub fn read_cid_to_cpd_cache() -> Result<HashMap<String, CidCpdEntry>> {
    cache_io::load_json(&cid_to_cpd_path(), HashMap::new)
}

pub fn write_cid_to_cpd_cache(cache: &HashMap<String, CidCpdEntry>) -> Result<()> {
    ensure_cache_dir()?;
    let bytes =
        serde_json::to_vec_pretty(cache).context("failed to serialise CID→cpd cache to JSON")?;
    cache_io::with_write_lock(&cid_to_cpd_lock_path(), |_| {
        cache_io::atomic_write(&cid_to_cpd_path(), &bytes)
    })
}

/// Remove any orphaned `.cid_to_cpd.lock` file. Called at application
/// startup to recover from a crashed prior process.
pub fn clear_stale_locks() -> Result<()> {
    for path in [cid_to_cpd_lock_path(), modules_lock_path()] {
        if path.exists() {
            fs::remove_file(&path)
                .with_context(|| format!("failed to remove stale KEGG lock {}", path.display()))?;
        }
    }
    Ok(())
}

// ─── Modules cache + long-running PID lock (Track B) ──────────────────────
//
// Unlike the CID→cpd cache (small writes, brief locks), the modules cache
// is filled by a ~12-min long-running fetch loop. A simple write-time
// lock would let two concurrent app instances both decide a cache is
// missing the same entries and double-fetch — wasting bandwidth and
// risking KEGG's 403 rate-limit. The `.modules.lock` file is held for
// the entire fetch duration, carrying a PID + heartbeat timestamp so a
// crashed holder can be detected (heartbeat >90 s old) and reclaimed.
//
// Cache writes themselves still use atomic temp+rename. The long-running
// lock controls fetch-loop concurrency; atomic write covers final-step
// crash safety.

fn modules_lock_path() -> PathBuf {
    cache_dir().join(".modules.lock")
}

/// Wait for an active fetch holder to clear, up to 30 min.
const MODULES_LOCK_WAIT_TIMEOUT: Duration = Duration::from_secs(30 * 60);
/// Poll interval while waiting for the lock to clear.
const MODULES_LOCK_POLL_INTERVAL: Duration = Duration::from_secs(5);
/// `last_seen_at` older than this is treated as stale (holder crashed).
const MODULES_LOCK_STALE_THRESHOLD: chrono::Duration = chrono::Duration::seconds(90);
/// Heartbeat cadence — rewrite `last_seen_at` no more often than this.
const MODULES_LOCK_HEARTBEAT_INTERVAL: Duration = Duration::from_secs(30);

#[derive(Debug, serde::Serialize, serde::Deserialize)]
struct ModulesLockFile {
    pid: u32,
    last_seen_at: chrono::DateTime<chrono::Utc>,
}

/// RAII guard for the long-running modules fetch lock. Created via
/// `acquire_modules_fetch_lock()`; on drop, removes the lock file if it
/// still belongs to us (PID match). Caller MUST call `heartbeat()`
/// periodically during the fetch loop to keep the lock alive.
pub struct ModulesFetchGuard {
    path: PathBuf,
    pid: u32,
    last_heartbeat_local: Instant,
}

impl ModulesFetchGuard {
    /// Rewrite `last_seen_at` if ≥ heartbeat cadence has elapsed locally.
    /// Cheap when called too often (no-op until the cadence elapses).
    pub fn heartbeat(&mut self) -> Result<()> {
        if self.last_heartbeat_local.elapsed() < MODULES_LOCK_HEARTBEAT_INTERVAL {
            return Ok(());
        }
        let body = ModulesLockFile {
            pid: self.pid,
            last_seen_at: chrono::Utc::now(),
        };
        let bytes =
            serde_json::to_vec(&body).context("failed to serialise modules-lock heartbeat")?;
        cache_io::atomic_write(&self.path, &bytes)?;
        self.last_heartbeat_local = Instant::now();
        Ok(())
    }
}

impl Drop for ModulesFetchGuard {
    fn drop(&mut self) {
        // Read the current lock; only remove if it still belongs to us.
        // (Defensive — if another instance's stale-detection stole the
        // lock from us, we MUST NOT delete its file.)
        match fs::read(&self.path) {
            Ok(bytes) => match serde_json::from_slice::<ModulesLockFile>(&bytes) {
                Ok(body) if body.pid == self.pid => {
                    let _ = fs::remove_file(&self.path);
                }
                Ok(_) => {
                    // PID differs — lock was stolen. Leave it alone.
                }
                Err(_) => {
                    // Malformed — treat as orphan and remove.
                    let _ = fs::remove_file(&self.path);
                }
            },
            Err(_) => {
                // Already gone; nothing to do.
            }
        }
    }
}

/// Acquire the long-running modules fetch lock. Polls with bounded
/// backoff for up to 30 min waiting for any active holder to release.
/// If an existing lock's `last_seen_at` is older than 90 s, treats it as
/// orphaned (holder crashed) and overwrites it.
pub fn acquire_modules_fetch_lock() -> Result<ModulesFetchGuard> {
    ensure_cache_dir()?;
    let path = modules_lock_path();
    let self_pid = std::process::id();
    let start = Instant::now();

    loop {
        // If no lock present, take it.
        if !path.exists() {
            return write_lock(&path, self_pid);
        }
        // Lock exists — inspect its heartbeat.
        match fs::read(&path) {
            Ok(bytes) => match serde_json::from_slice::<ModulesLockFile>(&bytes) {
                Ok(body) => {
                    let age = chrono::Utc::now() - body.last_seen_at;
                    if age > MODULES_LOCK_STALE_THRESHOLD {
                        warn!(
                            existing_pid = body.pid,
                            age_secs = age.num_seconds(),
                            "modules.lock heartbeat is stale; overwriting"
                        );
                        return write_lock(&path, self_pid);
                    }
                    // Still alive — wait and retry.
                }
                Err(_) => {
                    // Malformed lock file — treat as orphan.
                    warn!(path = %path.display(), "modules.lock is malformed; overwriting");
                    return write_lock(&path, self_pid);
                }
            },
            Err(_) => {
                // Race: file disappeared since exists() check. Try take.
                return write_lock(&path, self_pid);
            }
        }
        if start.elapsed() > MODULES_LOCK_WAIT_TIMEOUT {
            anyhow::bail!(
                "modules.lock held by another writer for longer than {:?}",
                MODULES_LOCK_WAIT_TIMEOUT
            );
        }
        thread::sleep(MODULES_LOCK_POLL_INTERVAL);
    }
}

fn write_lock(path: &std::path::Path, pid: u32) -> Result<ModulesFetchGuard> {
    let body = ModulesLockFile {
        pid,
        last_seen_at: chrono::Utc::now(),
    };
    let bytes = serde_json::to_vec(&body).context("failed to serialise modules-lock body")?;
    cache_io::atomic_write(path, &bytes)?;
    Ok(ModulesFetchGuard {
        path: path.to_path_buf(),
        pid,
        last_heartbeat_local: Instant::now(),
    })
}

pub fn read_modules_cache() -> Result<KeggModulesCache> {
    cache_io::load_json(&modules_path(), KeggModulesCache::default)
}

pub fn write_modules_cache(cache: &KeggModulesCache) -> Result<()> {
    ensure_cache_dir()?;
    let bytes =
        serde_json::to_vec_pretty(cache).context("failed to serialise KeggModulesCache to JSON")?;
    cache_io::atomic_write(&modules_path(), &bytes)
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::Utc;
    use std::sync::Mutex;
    use tempfile::tempdir;

    // Tests share process-wide state; serialise them.
    static SERIAL: Mutex<()> = Mutex::new(());

    fn setup_tmp_root() -> tempfile::TempDir {
        let dir = tempdir().expect("tempdir");
        unsafe {
            std::env::set_var("KEGG_CACHE_DIR", dir.path());
        }
        dir
    }

    #[test]
    fn cache_dir_env_override_beats_data_dir_default() {
        // With KEGG_CACHE_DIR set, the env override must win over the new
        // dirs::data_dir() default — precedence is unchanged by the data-dir move.
        let _g = SERIAL.lock().unwrap();
        let tmp = setup_tmp_root();
        assert_eq!(cache_dir(), tmp.path());
    }

    #[test]
    fn cid_to_cpd_roundtrip_positive_and_negative() {
        let _g = SERIAL.lock().unwrap();
        let _tmp = setup_tmp_root();

        let mut cache = HashMap::new();
        cache.insert(
            "5793".to_string(),
            CidCpdEntry {
                cpd: Some("C00031".to_string()),
                fetched_at: Utc::now(),
            },
        );
        cache.insert(
            "12345".to_string(),
            CidCpdEntry {
                cpd: None,
                fetched_at: Utc::now(),
            },
        );
        write_cid_to_cpd_cache(&cache).expect("write");
        let read = read_cid_to_cpd_cache().expect("read");
        assert_eq!(read.len(), 2);
        assert_eq!(read.get("5793").unwrap().cpd.as_deref(), Some("C00031"));
        assert!(read.get("12345").unwrap().cpd.is_none());
    }

    #[test]
    fn cid_to_cpd_missing_key_vs_negative() {
        let _g = SERIAL.lock().unwrap();
        let _tmp = setup_tmp_root();

        let mut cache = HashMap::new();
        cache.insert(
            "12345".to_string(),
            CidCpdEntry {
                cpd: None,
                fetched_at: Utc::now(),
            },
        );
        write_cid_to_cpd_cache(&cache).expect("write");
        let read = read_cid_to_cpd_cache().expect("read");

        // Cache hit, negative.
        assert!(read.get("12345").unwrap().cpd.is_none());
        // Cache miss.
        assert!(!read.contains_key("99999"));
    }

    #[test]
    fn cid_to_cpd_clear_stale_locks_removes_orphan() {
        let _g = SERIAL.lock().unwrap();
        let _tmp = setup_tmp_root();
        ensure_cache_dir().unwrap();

        let lp = cid_to_cpd_lock_path();
        fs::File::create(&lp).expect("create lock");
        assert!(lp.exists());

        clear_stale_locks().expect("clear");
        assert!(!lp.exists());
    }

    #[test]
    fn cid_to_cpd_write_releases_lock_on_success() {
        let _g = SERIAL.lock().unwrap();
        let _tmp = setup_tmp_root();
        let cache = HashMap::new();
        write_cid_to_cpd_cache(&cache).expect("write");
        assert!(!cid_to_cpd_lock_path().exists());
    }
}
