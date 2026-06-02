//! Save / load of every Stage 1–3 user-tunable parameter as a JSON file,
//! plus SHA-256 hashes of the currently-loaded MS-DIAL `.txt` +
//! metadata `.csv` inputs for drift detection.
//!
//! See `openspec/specs/session-settings-io/spec.md` for the normative
//! contract. The on-disk schema is owned by `schema::Snapshot`; the
//! `SCHEMA_VERSION` constant + golden fixture + version-rock test form
//! the triple-lock against schema drift.

pub mod error;
pub mod hash;
pub mod schema;
pub mod validate;

pub use error::SnapshotError;
pub use schema::{InputFileEntry, InputFileRole, SCHEMA_VERSION, Snapshot};
pub use validate::{ValidationResets, validate_against_inputs};

use std::fs;
use std::path::Path;

use crate::app::{SessionInputs, SessionSettings};
use crate::data::IonMode;

/// A per-role hash mismatch surfaced by `Snapshot::diff_input_hashes`.
/// `current = None` means the snapshot expected this role but the user
/// has no file loaded in it (e.g., snapshot saved dual-mode, user
/// currently single-mode).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct InputHashMismatch {
    pub role: InputFileRole,
    pub saved_name: String,
    pub saved_sha256: String,
    /// `Some((current_name, current_sha256))` when the current session
    /// has a file in the same role; `None` when the user does not.
    pub current: Option<(String, String)>,
}

/// Build a `Snapshot` from current session state.
///
/// Computes per-role hashes via `hash::sha256_file`. Roles whose source
/// is `None` at save time (in-memory / fixture-built ion tables with
/// `txt_path == None`; `inputs.csv_path == None`) are OMITTED from
/// `input_files` — see the "Snapshot save SHALL be available outside
/// `Initializing`" requirement in the spec.
pub fn from_session(
    settings: &SessionSettings,
    inputs: &SessionInputs,
    user_note: &str,
) -> Result<Snapshot, SnapshotError> {
    let mut input_files: Vec<InputFileEntry> = Vec::new();

    for table in &inputs.ion_tables {
        let Some(path) = table.txt_path.as_ref() else {
            continue;
        };
        let sha256 = hash::sha256_file(path).map_err(|source| SnapshotError::HashIo {
            path: path.clone(),
            source,
        })?;
        let name = path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();
        let role = match table.mode {
            IonMode::Positive => InputFileRole::Positive,
            IonMode::Negative => InputFileRole::Negative,
        };
        input_files.push(InputFileEntry { role, name, sha256 });
    }

    if let Some(path) = inputs.csv_path.as_ref() {
        let sha256 = hash::sha256_file(path).map_err(|source| SnapshotError::HashIo {
            path: path.clone(),
            source,
        })?;
        let name = path
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("")
            .to_string();
        input_files.push(InputFileEntry {
            role: InputFileRole::Metadata,
            name,
            sha256,
        });
    }

    Ok(Snapshot {
        schema_version: SCHEMA_VERSION,
        app_version: env!("CARGO_PKG_VERSION").to_string(),
        saved_at: chrono::Utc::now().to_rfc3339_opts(chrono::SecondsFormat::Secs, true),
        user_note: user_note.to_string(),
        input_files,
        settings: settings.clone(),
    })
}

/// Build a `Snapshot` and write it to `path` as pretty-printed UTF-8
/// JSON. Hash computation, serialise, and write all happen on the
/// calling thread (small budget: ≤ a few MB inputs, ≤ ~10 KB JSON).
///
/// On failure, the function MUST NOT leave a partial file at `path`.
/// `fs::write` is atomic on most platforms when the underlying write
/// fits in one syscall; for safety we never write before all hashes
/// succeed.
pub fn save_to_path(
    path: &Path,
    settings: &SessionSettings,
    inputs: &SessionInputs,
    user_note: &str,
) -> Result<(), SnapshotError> {
    let snapshot = from_session(settings, inputs, user_note)?;
    let json = serde_json::to_string_pretty(&snapshot).map_err(|e| SnapshotError::WriteIo {
        path: path.to_path_buf(),
        source: std::io::Error::other(e.to_string()),
    })?;
    fs::write(path, json).map_err(|source| SnapshotError::WriteIo {
        path: path.to_path_buf(),
        source,
    })?;
    Ok(())
}

/// Parse a settings JSON file and return the `Snapshot`. Validates
/// `schema_version` BEFORE deserialising the `settings` body so a
/// future-version payload cannot accidentally parse under the current
/// contract.
///
/// Error precedence: `Io` (read failure) → `JsonParse` (parse failure)
/// → `UnsupportedSchemaVersion`. The function does NOT mutate any
/// caller state — the caller decides what to do with the result.
pub fn load_from_path(path: &Path) -> Result<Snapshot, SnapshotError> {
    let bytes = fs::read(path).map_err(|source| SnapshotError::Io {
        path: path.to_path_buf(),
        source,
    })?;

    // First pass: parse to a generic Value so we can check
    // `schema_version` before deserialising the (potentially
    // future-shaped) `settings` block.
    let value: serde_json::Value = match serde_json::from_slice(&bytes) {
        Ok(v) => v,
        Err(e) => {
            return Err(SnapshotError::JsonParse {
                path: path.to_path_buf(),
                line: e.line(),
                column: e.column(),
                message: e.to_string(),
            });
        }
    };

    let found = value
        .get("schema_version")
        .and_then(|v| v.as_u64())
        .ok_or_else(|| SnapshotError::JsonParse {
            path: path.to_path_buf(),
            line: 0,
            column: 0,
            message: "missing or non-integer top-level field `schema_version`".to_string(),
        })?;
    if found != SCHEMA_VERSION as u64 {
        return Err(SnapshotError::UnsupportedSchemaVersion {
            found: found as u32,
            expected: SCHEMA_VERSION,
        });
    }

    // Second pass: full deserialise into `Snapshot`. `deny_unknown_fields`
    // catches typos / forward-incompat additions.
    let snapshot: Snapshot =
        serde_json::from_value(value).map_err(|e| SnapshotError::JsonParse {
            path: path.to_path_buf(),
            line: e.line(),
            column: e.column(),
            message: e.to_string(),
        })?;

    Ok(snapshot)
}

impl Snapshot {
    /// Compare this snapshot's `input_files` against the file hashes
    /// that `inputs` would produce right now.
    ///
    /// Semantics (per the spec):
    /// - Empty `input_files` (snapshot saved before any input loaded)
    ///   → empty result.
    /// - For each `entry` in `self.input_files`, locate the
    ///   corresponding current file by `role`. If found and hashes
    ///   differ → `Mismatch` with `current = Some(...)`. If the user
    ///   does not have a file in that role → `Mismatch` with
    ///   `current = None`.
    /// - Roles present in `inputs` but absent from `self.input_files`
    ///   are NOT reported (additive direction does not violate the
    ///   snapshot's reproducibility claim).
    pub fn diff_input_hashes(
        &self,
        inputs: &SessionInputs,
    ) -> Result<Vec<InputHashMismatch>, SnapshotError> {
        let mut out = Vec::new();

        for entry in &self.input_files {
            let current: Option<(String, String)> = match entry.role {
                InputFileRole::Positive | InputFileRole::Negative => {
                    let target_mode = if entry.role == InputFileRole::Positive {
                        IonMode::Positive
                    } else {
                        IonMode::Negative
                    };
                    let lookup = inputs
                        .ion_tables
                        .iter()
                        .find(|t| t.mode == target_mode)
                        .and_then(|t| t.txt_path.as_ref());
                    match lookup {
                        Some(p) => {
                            let sha256 =
                                hash::sha256_file(p).map_err(|source| SnapshotError::HashIo {
                                    path: p.clone(),
                                    source,
                                })?;
                            let name = p
                                .file_name()
                                .and_then(|s| s.to_str())
                                .unwrap_or("")
                                .to_string();
                            Some((name, sha256))
                        }
                        None => None,
                    }
                }
                InputFileRole::Metadata => match inputs.csv_path.as_ref() {
                    Some(p) => {
                        let sha256 =
                            hash::sha256_file(p).map_err(|source| SnapshotError::HashIo {
                                path: p.clone(),
                                source,
                            })?;
                        let name = p
                            .file_name()
                            .and_then(|s| s.to_str())
                            .unwrap_or("")
                            .to_string();
                        Some((name, sha256))
                    }
                    None => None,
                },
            };

            // Report a mismatch when either:
            // - current is None (snapshot expected, user doesn't have)
            // - current is Some(_) and sha256 differs from the saved value
            let mismatch = match &current {
                None => true,
                Some((_, current_sha)) => current_sha != &entry.sha256,
            };
            if mismatch {
                out.push(InputHashMismatch {
                    role: entry.role,
                    saved_name: entry.name.clone(),
                    saved_sha256: entry.sha256.clone(),
                    current,
                });
            }
        }

        Ok(out)
    }
}

impl SessionSettings {
    /// Apply a snapshot's settings to `*self`, honoring per-field
    /// resets recorded by `validate_against_inputs`.
    ///
    /// For each field listed in `resets` whose carrier is `Some(_)`,
    /// the corresponding field on `incoming` is set to `None` BEFORE
    /// the overwrite. Every other field on `incoming` overwrites
    /// `*self` verbatim.
    pub fn apply_snapshot(&mut self, mut incoming: SessionSettings, resets: &ValidationResets) {
        if resets.numerator.is_some() {
            incoming.numerator = None;
        }
        if resets.denominator.is_some() {
            incoming.denominator = None;
        }
        if resets.metadata_column.is_some() {
            incoming.metadata_column = None;
        }
        if resets.pqn_reference_group.is_some() {
            incoming.pqn_reference_group = None;
        }
        // Defensive coercion: Stage 2 DAM has no `NoCorrection` radio in
        // its UI (raw p-values across ~13 k features would flood the
        // result set with false positives). A hand-crafted or
        // future-version snapshot carrying `NoCorrection` here is
        // silently coerced back to BH; the user gets a WARN in logs.
        if matches!(
            incoming.dam_fdr_method,
            crate::dam::fdr::FdrMethod::NoCorrection
        ) {
            tracing::warn!(
                "snapshot dam_fdr_method=NoCorrection coerced to BenjaminiHochberg (Stage 2 UI never exposes None)"
            );
            incoming.dam_fdr_method = crate::dam::fdr::FdrMethod::BenjaminiHochberg;
        }
        *self = incoming;
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::{IonMode, IonModeTable, MetabolomicsTable};
    use ndarray::Array2;
    use std::io::Write;
    use std::path::PathBuf;
    use tempfile::{NamedTempFile, TempDir};

    fn empty_table() -> MetabolomicsTable {
        let intensity = Array2::<f64>::zeros((0, 0));
        MetabolomicsTable {
            annotated_count: 0,
            features: vec![],
            sample_cols: vec![],
            intensity_raw: intensity.clone(),
            intensity,
            excluded_cols: vec![],
        }
    }

    fn write_temp_with(name: &str, contents: &[u8]) -> (TempDir, PathBuf) {
        let dir = tempfile::tempdir().expect("tempdir");
        let path = dir.path().join(name);
        std::fs::write(&path, contents).expect("write temp file");
        (dir, path)
    }

    fn inputs_with(
        pos_path: Option<PathBuf>,
        neg_path: Option<PathBuf>,
        csv_path: Option<PathBuf>,
    ) -> SessionInputs {
        let mut ion_tables = Vec::new();
        if let Some(p) = pos_path {
            ion_tables.push(IonModeTable {
                mode: IonMode::Positive,
                table: empty_table(),
                txt_path: Some(p),
            });
        }
        if let Some(p) = neg_path {
            ion_tables.push(IonModeTable {
                mode: IonMode::Negative,
                table: empty_table(),
                txt_path: Some(p),
            });
        }
        SessionInputs {
            ion_tables,
            mapping: None,
            csv_path,
        }
    }

    fn build_snapshot_with_input_files(input_files: Vec<InputFileEntry>) -> Snapshot {
        Snapshot {
            schema_version: SCHEMA_VERSION,
            app_version: "0.0.0-test".to_string(),
            saved_at: "2026-05-26T00:00:00Z".to_string(),
            user_note: String::new(),
            input_files,
            settings: SessionSettings::default(),
        }
    }

    #[test]
    fn save_then_load_roundtrips() {
        let dir = tempfile::tempdir().expect("tempdir");
        let target = dir.path().join("snap.json");
        let settings = SessionSettings::default();
        let inputs = SessionInputs::default();

        save_to_path(&target, &settings, &inputs, "").expect("save");
        let snap = load_from_path(&target).expect("load");
        assert_eq!(snap.schema_version, SCHEMA_VERSION);
        assert_eq!(snap.user_note, "");
        assert!(snap.input_files.is_empty());
        assert_eq!(snap.settings, settings);
    }

    #[test]
    fn save_writes_pretty_printed_json() {
        let dir = tempfile::tempdir().expect("tempdir");
        let target = dir.path().join("snap.json");
        save_to_path(
            &target,
            &SessionSettings::default(),
            &SessionInputs::default(),
            "",
        )
        .expect("save");
        let text = std::fs::read_to_string(&target).expect("read back");
        // Pretty-printed JSON contains newlines + 2-space indent.
        assert!(text.contains('\n'), "expected newlines in pretty JSON");
        assert!(
            text.contains("\n  \""),
            "expected 2-space indent on top-level fields"
        );
    }

    #[test]
    fn load_rejects_schema_version_4() {
        let mut f = NamedTempFile::new().expect("tempfile");
        let body = r#"{
            "schema_version": 4,
            "app_version": "0.0.0",
            "saved_at": "2026-05-27T00:00:00Z",
            "user_note": "",
            "input_files": [],
            "settings": {}
        }"#;
        f.write_all(body.as_bytes()).expect("write");
        let err = load_from_path(f.path()).expect_err("should reject");
        match err {
            SnapshotError::UnsupportedSchemaVersion { found, expected } => {
                assert_eq!(found, 4);
                assert_eq!(expected, 3);
            }
            other => panic!("expected UnsupportedSchemaVersion, got {other:?}"),
        }
    }

    #[test]
    fn load_rejects_malformed_json() {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(b"{not valid json").expect("write");
        let err = load_from_path(f.path()).expect_err("should reject");
        assert!(
            matches!(err, SnapshotError::JsonParse { .. }),
            "expected JsonParse, got {err:?}"
        );
    }

    #[test]
    fn load_rejects_missing_schema_version() {
        let mut f = NamedTempFile::new().expect("tempfile");
        f.write_all(b"{\"foo\": 1}").expect("write");
        let err = load_from_path(f.path()).expect_err("should reject");
        match err {
            SnapshotError::JsonParse { message, .. } => {
                assert!(
                    message.contains("schema_version"),
                    "expected message about schema_version, got: {message}"
                );
            }
            other => panic!("expected JsonParse, got {other:?}"),
        }
    }

    #[test]
    fn diff_empty_input_files_returns_empty() {
        let snap = build_snapshot_with_input_files(vec![]);
        let inputs = SessionInputs::default();
        let diffs = snap.diff_input_hashes(&inputs).expect("diff");
        assert!(diffs.is_empty());
    }

    #[test]
    fn diff_all_match_returns_empty() {
        // Build a temp metadata file with known contents, hash it, then
        // construct a snapshot referencing the same hash.
        let (_dir, csv_path) = write_temp_with("metadata.csv", b"sample,group\nS01,A\n");
        let expected_sha = hash::sha256_file(&csv_path).expect("hash");
        let snap = build_snapshot_with_input_files(vec![InputFileEntry {
            role: InputFileRole::Metadata,
            name: "metadata.csv".to_string(),
            sha256: expected_sha,
        }]);
        let inputs = inputs_with(None, None, Some(csv_path));
        let diffs = snap.diff_input_hashes(&inputs).expect("diff");
        assert!(diffs.is_empty());
    }

    #[test]
    fn diff_one_mismatch_returns_one_record() {
        let (_dir, csv_path) = write_temp_with("metadata.csv", b"sample,group\nS01,A\n");
        let snap = build_snapshot_with_input_files(vec![InputFileEntry {
            role: InputFileRole::Metadata,
            name: "metadata.csv".to_string(),
            sha256: "0000000000000000000000000000000000000000000000000000000000000000".to_string(),
        }]);
        let inputs = inputs_with(None, None, Some(csv_path));
        let diffs = snap.diff_input_hashes(&inputs).expect("diff");
        assert_eq!(diffs.len(), 1);
        let m = &diffs[0];
        assert_eq!(m.role, InputFileRole::Metadata);
        assert_eq!(m.saved_name, "metadata.csv");
        assert!(m.current.is_some());
    }

    #[test]
    fn diff_snapshot_has_but_current_missing() {
        // Snapshot expects NEG + metadata; current has only metadata.
        let (_dir, csv_path) = write_temp_with("metadata.csv", b"sample,group\nS01,A\n");
        let csv_sha = hash::sha256_file(&csv_path).expect("hash");
        let snap = build_snapshot_with_input_files(vec![
            InputFileEntry {
                role: InputFileRole::Negative,
                name: "NEG.txt".to_string(),
                sha256: "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa"
                    .to_string(),
            },
            InputFileEntry {
                role: InputFileRole::Metadata,
                name: "metadata.csv".to_string(),
                sha256: csv_sha,
            },
        ]);
        let inputs = inputs_with(None, None, Some(csv_path));
        let diffs = snap.diff_input_hashes(&inputs).expect("diff");
        // Only NEG should be reported as mismatched; metadata matches.
        assert_eq!(diffs.len(), 1);
        assert_eq!(diffs[0].role, InputFileRole::Negative);
        assert!(diffs[0].current.is_none(), "current side should be None");
    }

    #[test]
    fn diff_current_has_but_snapshot_missing_not_reported() {
        // Snapshot has only metadata; current has POS + metadata.
        let (_dir1, csv_path) = write_temp_with("metadata.csv", b"sample,group\nS01,A\n");
        let (_dir2, pos_path) = write_temp_with("POS.txt", b"any positive payload");
        let csv_sha = hash::sha256_file(&csv_path).expect("hash");
        let snap = build_snapshot_with_input_files(vec![InputFileEntry {
            role: InputFileRole::Metadata,
            name: "metadata.csv".to_string(),
            sha256: csv_sha,
        }]);
        let inputs = inputs_with(Some(pos_path), None, Some(csv_path));
        let diffs = snap.diff_input_hashes(&inputs).expect("diff");
        // POS is current-only — must NOT appear.
        assert!(diffs.is_empty());
    }

    #[test]
    fn apply_snapshot_resets_invalid_then_overwrites() {
        let mut current = SessionSettings {
            fc_threshold: 1.1,
            numerator: Some("OldA".to_string()),
            ..SessionSettings::default()
        };
        let incoming = SessionSettings {
            fc_threshold: 9.9,
            numerator: Some("Treated".to_string()),
            ..SessionSettings::default()
        };
        let resets = ValidationResets {
            numerator: Some("Treated".to_string()),
            ..ValidationResets::default()
        };
        current.apply_snapshot(incoming, &resets);
        assert_eq!(current.numerator, None);
        assert_eq!(current.fc_threshold, 9.9);
    }

    #[test]
    fn apply_snapshot_with_no_resets_overwrites_verbatim() {
        let mut current = SessionSettings {
            fc_threshold: 1.1,
            ..SessionSettings::default()
        };
        let incoming = SessionSettings {
            fc_threshold: 9.9,
            numerator: Some("X".to_string()),
            ..SessionSettings::default()
        };
        current.apply_snapshot(incoming.clone(), &ValidationResets::default());
        assert_eq!(current.fc_threshold, 9.9);
        assert_eq!(current.numerator, Some("X".to_string()));
        assert_eq!(current, incoming);
    }

    #[test]
    fn apply_snapshot_coerces_dam_fdr_no_correction_to_bh() {
        // Adversarial / future-version snapshot carries NoCorrection in
        // the Stage 2 FDR field. The Stage 2 UI never offers None so
        // this would silently slip raw p-values into a 13-k-feature DAM
        // run; apply_snapshot must coerce back to BH.
        let mut current = SessionSettings::default();
        let incoming = SessionSettings {
            dam_fdr_method: crate::dam::fdr::FdrMethod::NoCorrection,
            // enrichment_fdr_method legitimately supports None — must be preserved.
            enrichment_fdr_method: crate::dam::fdr::FdrMethod::NoCorrection,
            ..SessionSettings::default()
        };
        current.apply_snapshot(incoming, &ValidationResets::default());
        assert_eq!(
            current.dam_fdr_method,
            crate::dam::fdr::FdrMethod::BenjaminiHochberg
        );
        assert_eq!(
            current.enrichment_fdr_method,
            crate::dam::fdr::FdrMethod::NoCorrection
        );
    }
}
