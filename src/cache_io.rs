//! Shared cache-IO + HTTP-client primitives used by both the KEGG and PubChem
//! cache/client modules. These four cache helpers (`atomic_write`,
//! `with_write_lock`, `load_json`, `resolve_cache_dir`) plus the `http_client`
//! builder were duplicated across `kegg/cache.rs` + `pubchem/cache.rs` and the
//! two clients; this module owns the single copy.
//!
//! Deliberately NOT here (per `dedup-network-cache-layer` design D3): the
//! long-running modules-fetch PID/heartbeat lock (`kegg::cache`), which is a
//! different mechanism with no PubChem analog. It keeps its own lifecycle and
//! merely calls `atomic_write` for the heartbeat body.

use anyhow::{Context, Result};
use serde::de::DeserializeOwned;
use std::path::{Path, PathBuf};
use std::sync::OnceLock;
use std::time::{Duration, Instant};

/// Maximum time a short-lived advisory write lock will wait for a prior holder
/// before bailing.
const LOCK_WAIT_TIMEOUT: Duration = Duration::from_secs(30);
const LOCK_POLL_INTERVAL: Duration = Duration::from_millis(100);

/// Atomic write: write to a sibling `.<name>.tmp.<pid>` file, fsync, then rename
/// onto `target`. On the same filesystem, `rename` is atomic on POSIX and on
/// Windows when `MoveFileExA` semantics apply. Lifted verbatim from the
/// (byte-identical) KEGG/PubChem cache copies; the KEGG modules heartbeat lock
/// keeps calling this for its lock-file body (design D3).
pub(crate) fn atomic_write(target: &Path, bytes: &[u8]) -> Result<()> {
    let dir = target.parent().unwrap_or_else(|| Path::new("."));
    let file_name = target.file_name().and_then(|s| s.to_str()).ok_or_else(|| {
        anyhow::anyhow!("cache target {} has no usable file name", target.display())
    })?;
    let tmp_path = dir.join(format!(".{file_name}.tmp.{}", std::process::id()));

    {
        let mut f = std::fs::File::create(&tmp_path)
            .with_context(|| format!("failed to create temp file {}", tmp_path.display()))?;
        std::io::Write::write_all(&mut f, bytes)
            .with_context(|| format!("failed to write temp file {}", tmp_path.display()))?;
        f.sync_all()
            .with_context(|| format!("failed to fsync temp file {}", tmp_path.display()))?;
    }

    std::fs::rename(&tmp_path, target).with_context(|| {
        format!(
            "failed to rename {} -> {}",
            tmp_path.display(),
            target.display()
        )
    })?;
    Ok(())
}

/// Acquire a short-lived advisory lock at `lock_path` (bounded 30 s wait),
/// run `body`, then remove the lock on success OR failure (so a subsequent
/// writer is never permanently blocked). The lock-file's parent dir must
/// already exist (callers run `ensure_cache_dir()` first, as before).
pub(crate) fn with_write_lock<T>(
    lock_path: &Path,
    body: impl FnOnce(&Path) -> Result<T>,
) -> Result<T> {
    let start = Instant::now();
    while lock_path.exists() {
        if start.elapsed() > LOCK_WAIT_TIMEOUT {
            anyhow::bail!(
                "cache lock {} held by another writer for longer than {:?}",
                lock_path.display(),
                LOCK_WAIT_TIMEOUT
            );
        }
        std::thread::sleep(LOCK_POLL_INTERVAL);
    }
    std::fs::File::create(lock_path)
        .with_context(|| format!("failed to create cache lock at {}", lock_path.display()))?;
    let out = body(lock_path);
    // Best-effort release on success OR failure.
    let _ = std::fs::remove_file(lock_path);
    out
}

/// Load-or-default: a missing file yields `on_missing()`; otherwise the file is
/// read and parsed as JSON into `T`. Each caller supplies its existing default
/// flavour (`|| None`, `HashMap::new`, `T::default`) so the read shapes are
/// byte-equivalent to the pre-extraction copies.
pub(crate) fn load_json<T: DeserializeOwned>(
    path: &Path,
    on_missing: impl FnOnce() -> T,
) -> Result<T> {
    if !path.exists() {
        return Ok(on_missing());
    }
    let bytes =
        std::fs::read(path).with_context(|| format!("failed to read {}", path.display()))?;
    serde_json::from_slice(&bytes).with_context(|| format!("failed to parse {}", path.display()))
}

/// Resolve a cache directory: test override cell → env var →
/// `<data_dir>/metabolopan/cache/<leaf>` → `./data/cache/<leaf>` fallback,
/// where `<data_dir>` is `dirs::data_dir()` (macOS `~/Library/Application
/// Support`, Linux `~/.local/share`, Windows `%APPDATA%`). Re-anchored off the
/// binary directory so the app runs from a read-only install (`.app` bundle,
/// `/Applications`, `C:\Program Files`); the `*_CACHE_DIR` env override keeps
/// its precedence ahead of the default. The override cell is passed in by
/// reference so KEGG and PubChem keep SEPARATE `OnceLock` statics (design D2 —
/// a shared static would couple test isolation).
pub(crate) fn resolve_cache_dir(
    env_var: &str,
    leaf: &str,
    override_cell: &OnceLock<PathBuf>,
) -> PathBuf {
    if let Some(p) = override_cell.get() {
        return p.clone();
    }
    if let Ok(env_path) = std::env::var(env_var)
        && !env_path.is_empty()
    {
        return PathBuf::from(env_path);
    }
    match dirs::data_dir() {
        Some(dir) => dir.join("metabolopan").join("cache").join(leaf),
        None => PathBuf::from(format!("data/cache/{leaf}")),
    }
}

/// Build a `reqwest::Client` carrying the crate User-Agent
/// (`<name>/<version>`). Used by BOTH the KEGG and PubChem clients so every
/// outbound request identifies the app (KEGG already sent this exact value;
/// PubChem gains it — the one intended behavior change of this refactor).
pub(crate) fn http_client(timeout: Duration) -> reqwest::Client {
    reqwest::Client::builder()
        .timeout(timeout)
        .user_agent(concat!(
            env!("CARGO_PKG_NAME"),
            "/",
            env!("CARGO_PKG_VERSION")
        ))
        .build()
        .expect("reqwest client builds")
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashMap;
    use tempfile::tempdir;

    #[test]
    fn atomic_write_round_trips() {
        let dir = tempdir().expect("tempdir");
        let target = dir.path().join("out.json");
        atomic_write(&target, b"hello bytes").expect("write");
        assert_eq!(std::fs::read(&target).expect("read"), b"hello bytes");
        // No stray temp file left behind.
        let leftovers: Vec<_> = std::fs::read_dir(dir.path())
            .unwrap()
            .filter_map(|e| e.ok())
            .filter(|e| e.file_name().to_string_lossy().contains(".tmp."))
            .collect();
        assert!(leftovers.is_empty(), "temp file should be renamed away");
    }

    #[test]
    fn load_json_missing_returns_on_missing_default() {
        let dir = tempdir().expect("tempdir");
        let missing = dir.path().join("absent.json");
        // `None` flavour.
        let opt: Option<u32> = load_json(&missing, || None).expect("load");
        assert_eq!(opt, None);
        // `HashMap::new` flavour.
        let map: HashMap<String, u32> = load_json(&missing, HashMap::new).expect("load");
        assert!(map.is_empty());
    }

    #[test]
    fn load_json_present_returns_parsed_value() {
        let dir = tempdir().expect("tempdir");
        let path = dir.path().join("present.json");
        let mut original = HashMap::new();
        original.insert("k".to_string(), 7u32);
        atomic_write(&path, &serde_json::to_vec(&original).unwrap()).expect("write");
        let loaded: HashMap<String, u32> = load_json(&path, HashMap::new).expect("load");
        assert_eq!(loaded.get("k"), Some(&7));
    }

    #[test]
    fn with_write_lock_runs_body_and_releases() {
        let dir = tempdir().expect("tempdir");
        let lock = dir.path().join(".test.lock");
        let target = dir.path().join("data.json");
        with_write_lock(&lock, |_| atomic_write(&target, b"x")).expect("locked write");
        assert!(target.exists(), "body ran");
        assert!(!lock.exists(), "lock released after success");
        // Lock is also released on body failure.
        let err = with_write_lock(&lock, |_| -> Result<()> { anyhow::bail!("boom") });
        assert!(err.is_err());
        assert!(!lock.exists(), "lock released after failure");
    }

    #[test]
    fn resolve_cache_dir_two_cells_stay_independent() {
        // Two distinct override cells → two distinct roots (the test-isolation
        // property that forbids merging the KEGG/PubChem statics, design D2).
        let dir_a = tempdir().expect("tempdir a");
        let dir_b = tempdir().expect("tempdir b");
        let cell_a: OnceLock<PathBuf> = OnceLock::new();
        let cell_b: OnceLock<PathBuf> = OnceLock::new();
        cell_a.set(dir_a.path().to_path_buf()).unwrap();
        cell_b.set(dir_b.path().to_path_buf()).unwrap();
        assert_eq!(
            resolve_cache_dir("KEGG_CACHE_DIR", "kegg", &cell_a),
            dir_a.path()
        );
        assert_eq!(
            resolve_cache_dir("PUBCHEM_CACHE_DIR", "pubchem", &cell_b),
            dir_b.path()
        );
    }

    #[test]
    fn resolve_cache_dir_default_is_under_data_dir() {
        // No override cell + an env var that is never set → the default resolves
        // under `dirs::data_dir()` (re-anchored off the binary directory), with a
        // CWD-relative fallback only when `data_dir()` is unavailable. The
        // override-cell branch (above) and the env-var branch (exercised by the
        // kegg/pubchem cache unit + integration tests) are unchanged.
        let cell: OnceLock<PathBuf> = OnceLock::new();
        let got = resolve_cache_dir("METABOLOPAN_NEVER_SET_CACHE_DIR", "kegg", &cell);
        match dirs::data_dir() {
            Some(base) => assert_eq!(got, base.join("metabolopan").join("cache").join("kegg")),
            None => assert_eq!(got, PathBuf::from("data/cache/kegg")),
        }
    }
}
