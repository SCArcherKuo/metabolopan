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
    p.push("settings_default_v3.json");
    p
}

/// Path to the original v1 fixture, preserved verbatim under a different
/// name so the v1-rejection path has a real-world v1 sample to load (rather
/// than an inline string constant that could drift from the actual
/// prior-released file format). Added by `add-log-transform-and-scaling` at
/// the v1 → v2 schema bump.
fn v1_rejected_fixture_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests");
    p.push("fixtures");
    p.push("settings_v1_rejected.json");
    p
}

/// Path to the demoted v2 fixture (former `settings_default_v2.json`),
/// preserved verbatim under a different name so the v2-rejection path has
/// a real-world v2 sample to load. Added by `add-min-entry-size-filter` at
/// the v2 → v3 schema bump (same pattern as the v1 demotion).
fn v2_rejected_fixture_path() -> PathBuf {
    let mut p = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    p.push("tests");
    p.push("fixtures");
    p.push("settings_v2_rejected.json");
    p
}

/// Triple-lock #1 of the schema-drift detection (per design D12): the
/// `SCHEMA_VERSION` constant. Any PR that bumps this MUST update this
/// test, which surfaces in review and forces a commit-message
/// explanation of why the bump happened.
#[test]
fn schema_version_is_locked_at_3() {
    assert_eq!(SCHEMA_VERSION, 3);
}

/// Triple-lock #2 of the schema-drift detection (per design D12): the
/// golden fixture roundtrip. Loading `settings_default_v3.json` and
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
fn golden_default_snapshot_v3_roundtrips() {
    let snap = load_from_path(&fixture_path()).expect("golden fixture must load");
    assert_eq!(snap.schema_version, SCHEMA_VERSION);
    assert_eq!(snap.app_version, "0.0.0-test");
    assert_eq!(snap.saved_at, "2026-05-27T00:00:00Z");
    assert_eq!(snap.user_note, "");
    assert!(snap.input_files.is_empty());
    assert_eq!(snap.settings, SessionSettings::default());
}

/// Regression: an actual v1 snapshot (preserved from the pre-2026-05-27
/// schema) must be rejected with `UnsupportedSchemaVersion { found: 1,
/// expected: SCHEMA_VERSION }`. The `expected` value rolls forward with
/// each schema bump; locking it via `SCHEMA_VERSION` rather than a literal
/// keeps this test correct across future bumps.
#[test]
fn v1_snapshot_is_rejected_with_unsupported_schema_version() {
    use metabolopan::session_io::SnapshotError;
    let result = load_from_path(&v1_rejected_fixture_path());
    match result {
        Err(SnapshotError::UnsupportedSchemaVersion { found, expected }) => {
            assert_eq!(found, 1);
            assert_eq!(expected, SCHEMA_VERSION);
            assert_eq!(expected, 3);
        }
        other => {
            panic!("expected UnsupportedSchemaVersion {{ found: 1, expected: 3 }}, got {other:?}")
        }
    }
}

/// Regression: the demoted v2 snapshot (former golden fixture) must be
/// rejected with `UnsupportedSchemaVersion { found: 2, expected: 3 }`.
/// Locks the v2-rejection contract in `session-settings-io`'s
/// "v2 snapshots are rejected after the v2→v3 bump" scenario.
#[test]
fn v2_snapshot_is_rejected_with_unsupported_schema_version() {
    use metabolopan::session_io::SnapshotError;
    let result = load_from_path(&v2_rejected_fixture_path());
    match result {
        Err(SnapshotError::UnsupportedSchemaVersion { found, expected }) => {
            assert_eq!(found, 2);
            assert_eq!(expected, SCHEMA_VERSION);
            assert_eq!(expected, 3);
        }
        other => {
            panic!("expected UnsupportedSchemaVersion {{ found: 2, expected: 3 }}, got {other:?}")
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
