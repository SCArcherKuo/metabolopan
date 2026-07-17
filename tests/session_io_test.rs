//! Integration tests for the `session-settings-io` capability.
//!
//! These tests are intentionally external to the `session_io` module
//! so they exercise the same `metabolopan::session_io::*` public API
//! that downstream callers (and the GUI) use.

use std::path::PathBuf;

use metabolopan::app::SessionSettings;
use metabolopan::session_io::{SCHEMA_VERSION, load_from_path};

/// Path to the golden fixture under `tests/fixtures/`.
/// The fixture represents `SessionSettings::default()` + the current
/// schema envelope; any drift between the in-code default and the
/// on-disk shape is caught here.
fn fixture_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests");
    p.push("fixtures");
    p.push("settings_default_v2.json");
    p
}

/// Path to the demoted v1 fixture (former `settings_default_v1.json`),
/// preserved verbatim under a different name so the schema-rejection path
/// has a real-world prior-release sample to load. Repurposed by
/// `add-rt-aware-dedup` at the v1 → v2 schema bump (same pattern as the
/// earlier v2 demotion; the former `settings_v2_rejected.json` was removed
/// because schema 2 is now the accepted current version).
fn v1_rejected_fixture_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests");
    p.push("fixtures");
    p.push("settings_v1_rejected.json");
    p
}

/// Triple-lock #1 of the schema-drift detection (per design D12): the
/// `SCHEMA_VERSION` constant. Any PR that bumps this MUST update this
/// test, which surfaces in review and forces a commit-message
/// explanation of why the bump happened.
#[test]
fn schema_version_is_locked_at_2() {
    assert_eq!(SCHEMA_VERSION, 2);
}

/// Triple-lock #2 of the schema-drift detection (per design D12): the
/// golden fixture roundtrip. Loading `settings_default_v2.json` and
/// comparing against `SessionSettings::default()` catches three classes
/// of drift in one assertion:
///
/// - Renaming a `SessionSettings` field without updating the fixture
///   → serde deny_unknown_fields rejects the load.
/// - Changing a default value in `src/app.rs` without updating the
///   fixture → `assert_eq!` reports the differing field.
/// - Adding a new field to `SessionSettings` → either fixture missing
///   it (serde error) or default-value drift on read.
#[test]
fn golden_default_snapshot_v2_roundtrips() {
    let snap = load_from_path(&fixture_path()).expect("golden fixture must load");
    assert_eq!(snap.schema_version, SCHEMA_VERSION);
    assert_eq!(snap.app_version, "0.0.0-test");
    assert_eq!(snap.saved_at, "2026-05-27T00:00:00Z");
    assert_eq!(snap.user_note, "");
    assert!(snap.input_files.is_empty());
    assert_eq!(snap.settings, SessionSettings::default());
}

/// Regression: the demoted v1 snapshot (former golden fixture) must be
/// rejected with `UnsupportedSchemaVersion { found: 1, expected: 2 }`.
/// Locks the v1-rejection contract in `session-settings-io`'s
/// "v1 snapshots are rejected (historical sample)" scenario.
#[test]
fn v1_snapshot_is_rejected_with_unsupported_schema_version() {
    use metabolopan::session_io::SnapshotError;
    let result = load_from_path(&v1_rejected_fixture_path());
    match result {
        Err(SnapshotError::UnsupportedSchemaVersion { found, expected }) => {
            assert_eq!(found, 1);
            assert_eq!(expected, SCHEMA_VERSION);
            assert_eq!(expected, 2);
        }
        other => {
            panic!("expected UnsupportedSchemaVersion {{ found: 1, expected: 2 }}, got {other:?}")
        }
    }
}

/// Triple-lock #3 of the schema-drift detection (per design D12):
/// the Stage 2 gate side-fix integration check. This test does not
/// load any fixture — it just drives the extracted
/// `check_group_membership` helper from the public API. Verifies the
/// "Treated" vs ["A", "B"] case from the proposal's example.
#[test]
fn stage2_gate_rejects_stale_numerator() {
    let groups = vec!["A".to_string(), "B".to_string()];
    let (num_ok, _den_ok) =
        metabolopan::ui::stage2_setup::check_group_membership(Some("Treated"), Some("A"), &groups);
    assert!(
        !num_ok,
        "Stage 2 gate must reject numerator that is not in current mapping.groups()"
    );
}
