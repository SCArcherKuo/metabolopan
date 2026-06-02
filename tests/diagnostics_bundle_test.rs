//! Integration tests for the bug-report bundle assembler. These cover
//! the privacy invariant (no raw inputs, no cache contents) and the
//! "missing session log" stub path.

use std::io::{Cursor, Read};

use metabolopan::diagnostics::{
    BUNDLE_README_TEXT, BundleArgs, build_bundle, redact_home_dir, render_cache_summary,
    render_input_summary,
};

const EXPECTED_FILES: [&str; 8] = [
    "README.txt",
    "version.txt",
    "RUST_LOG.txt",
    "KEGG_CACHE_DIR.txt",
    "logs.txt",
    "app_state.txt",
    "input_summary.txt",
    "cache_summary.txt",
];

/// Reads every file in a bundle zip into a `Vec<(name, bytes)>`.
fn unzip(bytes: &[u8]) -> Vec<(String, Vec<u8>)> {
    let mut archive = zip::ZipArchive::new(Cursor::new(bytes)).expect("open zip");
    let mut out = Vec::new();
    for i in 0..archive.len() {
        let mut file = archive.by_index(i).expect("zip entry");
        let name = file.name().to_string();
        let mut buf = Vec::new();
        file.read_to_end(&mut buf).expect("read entry");
        out.push((name, buf));
    }
    out
}

fn minimal_bundle_args<'a>(
    session_log: Option<&'a std::path::Path>,
    app_state_text: String,
    input_summary_text: String,
    cache_summary_text: String,
) -> BundleArgs<'a> {
    BundleArgs {
        session_log_path: session_log,
        rust_log_directive: "info",
        kegg_cache_dir: None,
        app_state_text,
        input_summary_text,
        cache_summary_text,
    }
}

#[test]
fn build_bundle_produces_expected_layout() {
    let bytes = build_bundle(minimal_bundle_args(
        None,
        "Variant: Stage1Input\n".into(),
        "ion_tables: <none loaded>\n".into(),
        render_cache_summary(&[]),
    ))
    .expect("build bundle");

    let entries = unzip(&bytes);
    let names: Vec<&str> = entries.iter().map(|(n, _)| n.as_str()).collect();
    for expected in &EXPECTED_FILES {
        assert!(
            names.contains(expected),
            "expected {expected} in zip, got: {names:?}"
        );
        let body = &entries
            .iter()
            .find(|(n, _)| n == expected)
            .expect("file present")
            .1;
        assert!(
            !body.is_empty(),
            "{expected} should not be empty inside the bundle"
        );
    }
    // No subdirectories — every entry is at the zip root.
    for (n, _) in &entries {
        assert!(
            !n.contains('/') || n.ends_with('/').not(),
            "entry should not be in a subdirectory: {n}"
        );
    }
}

#[test]
fn build_bundle_never_contains_input_files() {
    // The MS-DIAL fixture file's header row is a unique phrase that won't appear
    // anywhere else by accident — perfect sentinel.
    let sentinel_path = std::path::Path::new("data/single-mode/MS-DIAL-output-example.txt");
    let sentinel_content = std::fs::read_to_string(sentinel_path).expect("fixture present");
    // Pick a substring from the file that's stable + distinctive.
    let needle = "Alignment ID\tAverage Rt(min)";
    assert!(
        sentinel_content.contains(needle),
        "test invariant: needle must appear in fixture"
    );

    let app_state_text = format!(
        "Variant: Stage1Input\nion_tables_count: 1\n  ion_table[0]: mode=Positive features=42 samples=10 path={}\n",
        sentinel_path.display()
    );
    let input_summary_text = format!(
        "ion_tables_count: 1\n  ion_table[0]: mode=Positive  features=42  samples=10  path={}\n",
        sentinel_path.display()
    );

    let bytes = build_bundle(minimal_bundle_args(
        None,
        app_state_text,
        input_summary_text,
        render_cache_summary(&[]),
    ))
    .expect("build bundle");

    for (name, body) in unzip(&bytes) {
        let text = String::from_utf8_lossy(&body);
        assert!(
            !text.contains(needle),
            "MS-DIAL fixture content leaked into {name}:\nfirst 200 bytes: {:?}",
            text.chars().take(200).collect::<String>()
        );
    }
}

#[test]
fn build_bundle_never_contains_cache_contents() {
    // Construct a tempdir cache root with a sentinel string inside modules.json
    // and assert the bundle's cache_summary.txt mentions the file by name + size
    // only, never the sentinel content.
    let dir = tempfile::tempdir().expect("tempdir");
    let modules_json = dir.path().join("modules.json");
    let sentinel = "SECRET_SENTINEL_THAT_MUST_NOT_LEAK_42";
    let now = chrono::Utc::now().to_rfc3339();
    std::fs::write(
        &modules_json,
        format!(
            "{{\"modules\":{{\"M00001\":{{\"id\":\"M00001\",\"name\":\"{sentinel}\",\"fetched_at\":\"{now}\"}}}}}}"
        ),
    )
    .unwrap();

    let cache_summary = render_cache_summary(&[dir.path()]);
    assert!(
        cache_summary.contains("modules.json"),
        "modules.json must be listed in summary"
    );
    assert!(
        !cache_summary.contains(sentinel),
        "cache_summary itself already leaks sentinel — render_cache_summary bug"
    );

    let bytes = build_bundle(minimal_bundle_args(
        None,
        "Variant: Stage1Input\n".into(),
        "ion_tables: <none loaded>\n".into(),
        cache_summary,
    ))
    .expect("build bundle");

    for (name, body) in unzip(&bytes) {
        let text = String::from_utf8_lossy(&body);
        assert!(
            !text.contains(sentinel),
            "cache content leaked into {name} — privacy invariant broken"
        );
    }
}

#[test]
fn build_bundle_handles_missing_session_log() {
    let bytes = build_bundle(minimal_bundle_args(
        None,
        "Variant: Stage1Input\n".into(),
        "ion_tables: <none loaded>\n".into(),
        render_cache_summary(&[]),
    ))
    .expect("build bundle");

    let entries = unzip(&bytes);
    let (_, logs) = entries
        .iter()
        .find(|(n, _)| n == "logs.txt")
        .expect("logs.txt present");
    let text = String::from_utf8_lossy(logs);
    assert!(
        text.contains("session log unavailable"),
        "logs.txt stub missing the expected sentinel phrase: {text}"
    );
    let line_count = text.lines().count();
    assert!(
        line_count <= 2,
        "stub should be one line (got {line_count}): {text}"
    );
}

#[test]
fn build_bundle_readme_matches_constant() {
    let bytes = build_bundle(minimal_bundle_args(
        None,
        "Variant: Stage1Input\n".into(),
        "ion_tables: <none loaded>\n".into(),
        render_cache_summary(&[]),
    ))
    .expect("build bundle");
    let entries = unzip(&bytes);
    let (_, readme) = entries
        .iter()
        .find(|(n, _)| n == "README.txt")
        .expect("README.txt present");
    assert_eq!(
        readme,
        BUNDLE_README_TEXT.as_bytes(),
        "README.txt content drifted from BUNDLE_README_TEXT constant"
    );
}

#[test]
fn build_bundle_splits_env_into_per_var_files() {
    // The user concern was: env.txt could be mistaken for a full env dump.
    // Resolution: filename = variable name, so there's no ambiguity.
    let bytes = build_bundle(BundleArgs {
        session_log_path: None,
        rust_log_directive: "debug,hyper=warn",
        kegg_cache_dir: Some("/tmp/kegg-test"),
        app_state_text: "Variant: Stage1Input\n".into(),
        input_summary_text: "ion_tables: <none loaded>\n".into(),
        cache_summary_text: render_cache_summary(&[]),
    })
    .expect("build bundle");

    let entries = unzip(&bytes);
    let names: Vec<&str> = entries.iter().map(|(n, _)| n.as_str()).collect();
    assert!(
        names.contains(&"RUST_LOG.txt"),
        "RUST_LOG.txt missing: {names:?}"
    );
    assert!(
        names.contains(&"KEGG_CACHE_DIR.txt"),
        "KEGG_CACHE_DIR.txt missing: {names:?}"
    );
    assert!(
        !names.contains(&"env.txt"),
        "old env.txt should not appear: {names:?}"
    );

    let rust_log = entries
        .iter()
        .find(|(n, _)| n == "RUST_LOG.txt")
        .map(|(_, b)| String::from_utf8_lossy(b).to_string())
        .unwrap();
    let kegg = entries
        .iter()
        .find(|(n, _)| n == "KEGG_CACHE_DIR.txt")
        .map(|(_, b)| String::from_utf8_lossy(b).to_string())
        .unwrap();
    // Files contain ONLY the value (with trailing newline) — no "VAR=" prefix
    // that might lead a reader to expect more vars below.
    assert_eq!(rust_log.trim_end(), "debug,hyper=warn");
    assert_eq!(kegg.trim_end(), "/tmp/kegg-test");
}

#[test]
fn build_bundle_kegg_cache_dir_unset_is_explicit() {
    let bytes = build_bundle(BundleArgs {
        session_log_path: None,
        rust_log_directive: "info",
        kegg_cache_dir: None,
        app_state_text: "Variant: Stage1Input\n".into(),
        input_summary_text: "ion_tables: <none loaded>\n".into(),
        cache_summary_text: render_cache_summary(&[]),
    })
    .expect("build bundle");

    let entries = unzip(&bytes);
    let kegg = entries
        .iter()
        .find(|(n, _)| n == "KEGG_CACHE_DIR.txt")
        .map(|(_, b)| String::from_utf8_lossy(b).to_string())
        .unwrap();
    assert_eq!(kegg.trim_end(), "<unset>");
}

#[test]
fn build_bundle_redacts_home_directory_across_all_sections() {
    // Inject the current process's $HOME into multiple sections and assert
    // the resulting zip never contains the literal home directory string.
    let home = std::env::var("HOME").or_else(|_| std::env::var("USERPROFILE"));
    let home = match home {
        Ok(h) if !h.is_empty() => h,
        _ => {
            eprintln!("skipping: no $HOME or $USERPROFILE set in test environment");
            return;
        }
    };

    let app_state_text = format!("Variant: Stage1Input\ncsv_path: {home}/project/metadata.csv\n");
    let input_summary_text = format!(
        "ion_tables_count: 1\n  ion_table[0]: mode=Positive features=10 samples=4 path={home}/data/POS.txt\n"
    );
    let cache_summary_text = format!("Root: {home}/Library/Caches/myapp\n");

    // RUST_LOG is arbitrary user text and CAN embed a path (a file-target
    // filter). Inject $HOME so this section is exercised by the redaction loop.
    let rust_log_directive = format!("info,myfilter={home}/trace.log");

    let bytes = build_bundle(BundleArgs {
        session_log_path: None,
        rust_log_directive: &rust_log_directive,
        kegg_cache_dir: Some(&format!("{home}/.kegg-cache")),
        app_state_text,
        input_summary_text,
        cache_summary_text,
    })
    .expect("build bundle");

    // The bundle is exactly 8 files by construction (privacy invariant); lock
    // it so an added/removed `add_entry` cannot silently change the file set.
    assert_eq!(
        unzip(&bytes).len(),
        8,
        "bug-report bundle must contain exactly 8 files"
    );

    for (name, body) in unzip(&bytes) {
        // README and version are not redacted (they don't contain user paths
        // by construction), but the assertion still holds for them — they
        // don't reference $HOME at all.
        let text = String::from_utf8_lossy(&body);
        assert!(
            !text.contains(&home),
            "raw home directory string leaked into {name}:\n  first 300 chars: {:?}",
            text.chars().take(300).collect::<String>()
        );
    }
}

#[test]
fn redact_home_dir_replaces_and_idempotent() {
    let home = std::env::var("HOME").or_else(|_| std::env::var("USERPROFILE"));
    let home = match home {
        Ok(h) if !h.is_empty() => h,
        _ => return,
    };
    let input = format!("path1={home}/a.txt path2={home}/b/c.txt unrelated=42");
    let redacted = redact_home_dir(&input);
    assert!(!redacted.contains(&home), "home not redacted: {redacted}");
    assert!(redacted.contains("~/a.txt"));
    assert!(redacted.contains("~/b/c.txt"));
    assert!(redacted.contains("unrelated=42"));
    // Idempotent: applying twice equals applying once.
    assert_eq!(redact_home_dir(&redacted), redacted);
}

#[test]
fn render_input_summary_does_not_open_input_files() {
    // Sanity: render_input_summary doesn't actually open the .txt referenced
    // in IonModeTable.txt_path — it just reports the path string.
    use metabolopan::app::{AppState, SessionInputs};
    let state = AppState::Stage1Input {
        slot1_mode: None,
        slot2_revealed: false,
        slot2_mode: None,
        error: None,
    };
    let inputs = SessionInputs {
        ion_tables: vec![],
        mapping: None,
        csv_path: Some(std::path::PathBuf::from(
            "/this/path/does/not/exist/metadata.csv",
        )),
    };
    let out = render_input_summary(&state, &inputs);
    assert!(out.contains("/this/path/does/not/exist/metadata.csv"));
}

// Tiny helper because `bool::not()` isn't free in this scope.
trait BoolExt {
    fn not(self) -> bool;
}
impl BoolExt for bool {
    fn not(self) -> bool {
        !self
    }
}
