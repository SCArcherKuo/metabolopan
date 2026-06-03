use std::fs::{self, File, OpenOptions};
use std::io;
use std::path::{Path, PathBuf};
use std::time::SystemTime;

use chrono::{DateTime, Duration, Utc};
use tracing_subscriber::EnvFilter;
use tracing_subscriber::fmt::{format::Writer, time::FormatTime};

/// Fixed `EnvFilter` directive for the file sink. Independent of
/// `RUST_LOG` by design (D2): crate targets at `info`, HTTP/TLS
/// dependency crates pinned at `warn`. The on-disk session log is
/// bug-report-grade signal, not the live developer feed.
pub const FILE_SINK_DIRECTIVE: &str = "metabolopan=info,\
    h2=warn,hyper=warn,reqwest=warn,rustls=warn,hickory_resolver=warn";

/// Returns `<data_dir>/metabolopan/logs/`, where `<data_dir>` is
/// `dirs::data_dir()` — the same anchor as `kegg::cache::cache_dir()`,
/// re-anchored off the binary directory so a read-only install (`.app`
/// bundle, `/Applications`) can still write logs. Falls back to a
/// CWD-relative `./data/logs` only when `data_dir()` is unavailable.
pub fn session_log_dir() -> PathBuf {
    match dirs::data_dir() {
        Some(dir) => dir.join("metabolopan").join("logs"),
        None => PathBuf::from("data/logs"),
    }
}

/// Returns `session_YYYYMMDD_HHMMSS_<pid>.log`. UTC timestamp.
pub fn session_log_path(started_at: DateTime<Utc>, pid: u32) -> PathBuf {
    let stem = started_at.format("session_%Y%m%d_%H%M%S").to_string();
    PathBuf::from(format!("{stem}_{pid}.log"))
}

/// Compact UTC time formatter: `HH:MM:SS.mmm UTC`. Used by the file-sink
/// fmt layer in place of the default ISO-8601 / system-time formats so
/// the on-disk log keeps a compact body (`HH:MM:SS.mmm`) parallel to the
/// in-window log pane. The suffix intentionally diverges by audience:
/// ` UTC` on disk (developer-grade, used in bug-report bundles read
/// across timezones) versus ` +HHMM` in the pane (live operator
/// wall-clock readability) — see the `fix-log-pane-local-time` change
/// (design D2) for the rationale.
///
/// Exposed publicly because `main.rs` assembles the fmt layer inline (the
/// generic-S `impl Layer<S>` return type would bake in S at the wrong
/// point in the layer chain — see D8 commentary on layer composition).
pub struct CompactUtcTime;

impl FormatTime for CompactUtcTime {
    fn format_time(&self, w: &mut Writer<'_>) -> std::fmt::Result {
        let now = Utc::now();
        write!(w, "{}", now.format("%H:%M:%S%.3f UTC"))
    }
}

/// Returns the fixed `EnvFilter` used by the file sink layer. Always
/// `FILE_SINK_DIRECTIVE`, regardless of `RUST_LOG`.
pub fn file_sink_env_filter() -> EnvFilter {
    EnvFilter::new(FILE_SINK_DIRECTIVE)
}

/// Opens (creates if absent) the session log file under `dir` and
/// returns its absolute path plus the open `File` handle. Caller wraps
/// the file in a `Mutex` and feeds it to `fmt::layer().with_writer(...)`.
///
/// Best-effort `create_dir_all(dir)` runs first. IO errors are
/// returned verbatim so the caller can WARN-and-continue per the
/// "missing file sink does not block startup" rule.
pub fn try_open_session_log(
    dir: &Path,
    started_at: DateTime<Utc>,
    pid: u32,
) -> io::Result<(PathBuf, File)> {
    fs::create_dir_all(dir)?;
    let path = dir.join(session_log_path(started_at, pid));
    let file = OpenOptions::new().create(true).append(true).open(&path)?;
    Ok((path, file))
}

/// Result of `clear_stale_session_logs`. `failures` carries per-file IO
/// errors so the caller can WARN each one without aborting startup.
#[derive(Debug, Default)]
pub struct CleanupReport {
    pub deleted: Vec<PathBuf>,
    pub retained: usize,
    pub failures: Vec<(PathBuf, io::Error)>,
    pub skipped_missing_dir: bool,
}

/// Removes `session_*.log` files older than `max_age_days` from `dir`.
/// All failure modes are surfaced through `CleanupReport` rather than
/// returned as `Result::Err` — startup must never abort on log cleanup.
pub fn clear_stale_session_logs(
    dir: &Path,
    now: DateTime<Utc>,
    max_age_days: u32,
) -> CleanupReport {
    let mut report = CleanupReport::default();
    if !dir.exists() {
        report.skipped_missing_dir = true;
        return report;
    }
    let read_dir = match fs::read_dir(dir) {
        Ok(rd) => rd,
        Err(e) => {
            report.failures.push((dir.to_path_buf(), e));
            return report;
        }
    };
    let threshold = now - Duration::days(i64::from(max_age_days));
    for entry in read_dir.flatten() {
        let path = entry.path();
        let name = match path.file_name().and_then(|s| s.to_str()) {
            Some(n) => n,
            None => continue,
        };
        if !is_session_log_name(name) {
            continue;
        }
        let mtime: SystemTime = match entry.metadata().and_then(|m| m.modified()) {
            Ok(t) => t,
            Err(e) => {
                report.failures.push((path, e));
                continue;
            }
        };
        let mtime_utc: DateTime<Utc> = mtime.into();
        if mtime_utc < threshold {
            match fs::remove_file(&path) {
                Ok(()) => report.deleted.push(path),
                Err(e) => report.failures.push((path, e)),
            }
        } else {
            report.retained += 1;
        }
    }
    report
}

/// Matches the `session_*.log` filename pattern emitted by
/// `session_log_path`. Tight enough to ignore `unrelated.txt`,
/// `session.txt`, `random.log`, etc.
fn is_session_log_name(name: &str) -> bool {
    name.starts_with("session_") && name.ends_with(".log")
}

#[cfg(test)]
mod tests {
    use super::*;
    use chrono::TimeZone;
    use filetime::{FileTime, set_file_mtime};
    use tempfile::tempdir;

    #[test]
    fn session_log_dir_is_under_data_dir() {
        // Re-anchored off the binary directory onto `dirs::data_dir()` so a
        // read-only install / `.app` bundle can still write logs; CWD-relative
        // fallback only when `data_dir()` is unavailable.
        match dirs::data_dir() {
            Some(base) => assert_eq!(session_log_dir(), base.join("metabolopan").join("logs")),
            None => assert_eq!(session_log_dir(), PathBuf::from("data/logs")),
        }
    }

    fn touch(path: &Path) {
        File::create(path).expect("create log fixture");
    }

    fn set_age(path: &Path, days_ago: i64) {
        let now = SystemTime::now();
        let target = now - std::time::Duration::from_secs((days_ago * 86_400) as u64);
        set_file_mtime(path, FileTime::from_system_time(target)).expect("set mtime");
    }

    #[test]
    fn session_log_path_format_is_stable() {
        let ts = Utc.with_ymd_and_hms(2026, 5, 24, 15, 30, 12).unwrap();
        let path = session_log_path(ts, 4242);
        assert_eq!(path.to_string_lossy(), "session_20260524_153012_4242.log");
    }

    #[test]
    fn clear_stale_session_logs_age_threshold() {
        let dir = tempdir().expect("tempdir");
        let stale = dir.path().join("session_20240101_000000_1.log");
        let fresh_edge = dir.path().join("session_20260517_000000_2.log");
        let fresh = dir.path().join("session_20260523_000000_3.log");
        touch(&stale);
        touch(&fresh_edge);
        touch(&fresh);
        set_age(&stale, 14);
        set_age(&fresh_edge, 6);
        set_age(&fresh, 1);

        let report = clear_stale_session_logs(dir.path(), Utc::now(), 7);

        assert_eq!(report.deleted, vec![stale.clone()]);
        assert_eq!(report.retained, 2);
        assert!(
            report.failures.is_empty(),
            "failures: {:?}",
            report.failures
        );
        assert!(!report.skipped_missing_dir);
        assert!(!stale.exists(), "stale file should be removed");
        assert!(fresh_edge.exists(), "edge-fresh file should be retained");
        assert!(fresh.exists(), "fresh file should be retained");
    }

    #[test]
    fn clear_stale_session_logs_missing_dir_no_error() {
        let dir = tempdir().expect("tempdir");
        let missing = dir.path().join("does_not_exist");
        let report = clear_stale_session_logs(&missing, Utc::now(), 7);
        assert!(report.skipped_missing_dir);
        assert!(report.failures.is_empty());
        assert!(report.deleted.is_empty());
        assert_eq!(report.retained, 0);
    }

    #[test]
    fn clear_stale_session_logs_per_file_failure_does_not_abort() {
        // Cross-platform "make one entry undeletable while leaving others
        // deletable": create a stale DIRECTORY whose name passes the
        // `session_*.log` filter. `fs::remove_file` returns an error when
        // asked to delete a directory (EISDIR on Unix, ERROR_ACCESS_DENIED
        // on Windows) — same effect as a permission-denied file deletion
        // without needing chmod / chflags / Administrator. Covers the spec
        // scenario "Per-file deletion failure does not abort cleanup".
        let dir = tempdir().expect("tempdir");
        let good = dir.path().join("session_good_1.log");
        let bad = dir.path().join("session_dir_should_fail_2.log");
        let fresh = dir.path().join("session_fresh_3.log");
        touch(&good);
        std::fs::create_dir(&bad).expect("create dir-named-as-log");
        touch(&fresh);
        set_age(&good, 30);
        // Setting mtime on a directory works on both Unix and Windows.
        set_age(&bad, 30);
        set_age(&fresh, 1);

        let report = clear_stale_session_logs(dir.path(), Utc::now(), 7);

        assert_eq!(
            report.deleted,
            vec![good.clone()],
            "only the regular stale file should be deleted; deleted={:?}",
            report.deleted
        );
        assert!(!good.exists(), "stale regular file should be gone");
        assert!(
            bad.exists(),
            "stale directory should still exist because remove_file failed on it"
        );
        assert!(fresh.exists(), "fresh file should be retained");
        assert_eq!(
            report.failures.len(),
            1,
            "exactly one failure expected; got {:?}",
            report.failures
        );
        assert_eq!(
            report.failures[0].0, bad,
            "failure path should be the directory"
        );
        assert_eq!(
            report.retained, 1,
            "the fresh file is retained (the failed entry is counted as a failure, not retained)"
        );
        assert!(!report.skipped_missing_dir);
    }

    #[test]
    fn clear_stale_session_logs_ignores_non_session_files() {
        let dir = tempdir().expect("tempdir");
        let session = dir.path().join("session_old_999.log");
        let unrelated = dir.path().join("unrelated.txt");
        let random_log = dir.path().join("random.log");
        touch(&session);
        touch(&unrelated);
        touch(&random_log);
        set_age(&session, 30);
        set_age(&unrelated, 30);
        set_age(&random_log, 30);

        let report = clear_stale_session_logs(dir.path(), Utc::now(), 7);

        assert_eq!(report.deleted, vec![session.clone()]);
        assert!(!session.exists());
        assert!(unrelated.exists(), "unrelated.txt must not be touched");
        assert!(random_log.exists(), "random.log must not be touched");
    }
}
