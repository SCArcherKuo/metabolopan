//! JSON schema types for the `session-settings-io` capability.
//!
//! `Snapshot` is the on-disk representation of a saved session. All field
//! names + variant names are part of the schema contract — renaming or
//! removing any of them requires a `SCHEMA_VERSION` bump (see the
//! "Snapshot schema version SHALL be a tracked invariant" requirement in
//! `openspec/specs/session-settings-io/spec.md`).

use serde::{Deserialize, Serialize};

use crate::app::SessionSettings;

/// Current schema version. Bumped exactly one when any field of
/// `Snapshot`, `InputFileEntry`, `SessionSettings`, or any enum
/// transitively referenced by `SessionSettings` changes name, is removed,
/// or is added. Triple-locked against drift by the
/// `tests/fixtures/settings_default_v3.json` golden fixture + the
/// version-rock test in `tests/session_io_test.rs`.
pub const SCHEMA_VERSION: u32 = 3;

/// Top-level shape of a saved settings JSON file.
#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Snapshot {
    pub schema_version: u32,
    pub app_version: String,
    /// UTC RFC3339 timestamp (e.g., `"2026-05-26T14:32:11Z"`); produced
    /// by `chrono::Utc::now().to_rfc3339_opts(SecondsFormat::Secs, true)`
    /// on save.
    pub saved_at: String,
    /// Always written as `""` on save. Users hand-edit the JSON if they
    /// want to attach a comment.
    pub user_note: String,
    /// Per-role hashes of inputs loaded at save time. May be empty when
    /// the snapshot was saved before any inputs were loaded; may be
    /// shorter than `inputs.ion_tables.len() + 1` when some `txt_path` /
    /// `csv_path` was `None` at save time.
    pub input_files: Vec<InputFileEntry>,
    pub settings: SessionSettings,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct InputFileEntry {
    pub role: InputFileRole,
    /// Basename of the source file (no parent directory).
    pub name: String,
    /// SHA-256 of the file bytes as a 64-character lowercase hex string.
    pub sha256: String,
}

/// Role of an input file in the saved session. Serde uses lowercase
/// strings in the JSON form (`"positive"`, `"negative"`, `"metadata"`);
/// in-code we work with the enum.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum InputFileRole {
    Positive,
    Negative,
    Metadata,
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dam::DamMethod;
    use crate::dam::fdr::FdrMethod;
    use crate::enrichment::EnrichmentDirection;
    use crate::normalize::{NormalizationMethod, PqnReference};

    /// Build a `SessionSettings` value with every field set to a
    /// non-default value, so a roundtrip test exercises every field's
    /// serde derive. Independent of `app.rs`'s in-test
    /// `non_default_settings` helper (we cannot import a `#[cfg(test)]`
    /// helper from another file).
    fn settings_with_every_field_non_default() -> SessionSettings {
        use crate::app::AnalysisMode;
        SessionSettings {
            analysis_mode: AnalysisMode::Module,
            kegg_species: Some("hsa".to_string()),
            organism_group_level: Some(2),
            organism_group: Some("Mammals".to_string()),
            min_group_overlap: 5,
            numerator: Some("treatment".to_string()),
            denominator: Some("control".to_string()),
            dam_method: DamMethod::Welch,
            drop_unknown: false,
            dedup_enabled: false,
            normalization: NormalizationMethod::Pqn {
                reference: PqnReference::Group("control".to_string()),
            },
            metadata_column: Some("dry_weight".to_string()),
            pqn_reference: PqnReference::Group("control".to_string()),
            pqn_reference_group: Some("control".to_string()),
            log_transform: false,
            dam_fdr_method: FdrMethod::BenjaminiYekutieli,
            fc_threshold: 4.0,
            fdr_threshold: 0.01,
            delta_threshold: 0.5,
            stage2_export_width_in: 6.0,
            stage2_export_height_in: 4.0,
            stage2_export_dpi: 600,
            direction: EnrichmentDirection::Up,
            top_n: 50,
            enrichment_fdr_threshold: 0.1,
            min_hit_count: 3,
            min_entry_size: 5,
            enrichment_fdr_method: FdrMethod::BenjaminiYekutieli,
            stage3_export_width_in: 5.0,
            stage3_export_height_in: 10.0,
            stage3_export_dpi: 600,
        }
    }

    fn sample_snapshot() -> Snapshot {
        Snapshot {
            schema_version: SCHEMA_VERSION,
            app_version: "0.0.0-test".to_string(),
            saved_at: "2026-05-26T00:00:00Z".to_string(),
            user_note: String::new(),
            input_files: vec![InputFileEntry {
                role: InputFileRole::Metadata,
                name: "metadata.csv".to_string(),
                sha256: "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
                    .to_string(),
            }],
            settings: settings_with_every_field_non_default(),
        }
    }

    #[test]
    fn snapshot_roundtrips_bit_equal() {
        let snap = sample_snapshot();
        let json = serde_json::to_string_pretty(&snap).expect("serialise");
        let parsed: Snapshot = serde_json::from_str(&json).expect("deserialise");
        assert_eq!(parsed.schema_version, snap.schema_version);
        assert_eq!(parsed.app_version, snap.app_version);
        assert_eq!(parsed.saved_at, snap.saved_at);
        assert_eq!(parsed.user_note, snap.user_note);
        assert_eq!(parsed.input_files, snap.input_files);
        assert_eq!(parsed.settings, snap.settings);
    }

    #[test]
    fn unit_variant_normalization_serialises_as_bare_string() {
        let mut s = settings_with_every_field_non_default();
        s.normalization = NormalizationMethod::Sum;
        let json = serde_json::to_string(&s).expect("serialise");
        assert!(
            json.contains("\"normalization\":\"Sum\""),
            "expected normalization to serialise as bare string, got: {json}"
        );
    }

    #[test]
    fn struct_variant_normalization_serialises_as_externally_tagged() {
        let mut s = settings_with_every_field_non_default();
        s.normalization = NormalizationMethod::Pqn {
            reference: PqnReference::AllSamples,
        };
        let json = serde_json::to_string(&s).expect("serialise");
        assert!(
            json.contains("\"normalization\":{\"Pqn\":{\"reference\":\"AllSamples\"}}"),
            "expected Pqn struct variant in externally-tagged form, got: {json}"
        );
    }

    #[test]
    fn deny_unknown_fields_rejects_extra_top_level_key() {
        // Build a valid snapshot, then inject an unknown top-level field.
        let snap = sample_snapshot();
        let mut value = serde_json::to_value(&snap).expect("to_value");
        if let serde_json::Value::Object(ref mut map) = value {
            map.insert(
                "future_field".to_string(),
                serde_json::Value::String("unexpected".to_string()),
            );
        }
        let json = serde_json::to_string(&value).expect("re-serialise");
        let parsed: Result<Snapshot, _> = serde_json::from_str(&json);
        assert!(
            parsed.is_err(),
            "deny_unknown_fields should reject extra top-level key, but parse succeeded"
        );
    }

    #[test]
    fn deny_unknown_fields_rejects_extra_input_file_key() {
        let mut snap = sample_snapshot();
        snap.input_files = vec![]; // we'll inject by hand below
        let mut value = serde_json::to_value(&snap).expect("to_value");
        // Inject an InputFileEntry-like object with an extra key.
        let entry = serde_json::json!({
            "role": "metadata",
            "name": "metadata.csv",
            "sha256": "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad",
            "future_field": "unexpected"
        });
        if let serde_json::Value::Object(ref mut map) = value {
            map.insert(
                "input_files".to_string(),
                serde_json::Value::Array(vec![entry]),
            );
        }
        let json = serde_json::to_string(&value).expect("re-serialise");
        let parsed: Result<Snapshot, _> = serde_json::from_str(&json);
        assert!(
            parsed.is_err(),
            "deny_unknown_fields should reject extra InputFileEntry key, but parse succeeded"
        );
    }

    #[test]
    fn input_file_role_serialises_as_lowercase() {
        let pos = serde_json::to_string(&InputFileRole::Positive).unwrap();
        let neg = serde_json::to_string(&InputFileRole::Negative).unwrap();
        let met = serde_json::to_string(&InputFileRole::Metadata).unwrap();
        assert_eq!(pos, "\"positive\"");
        assert_eq!(neg, "\"negative\"");
        assert_eq!(met, "\"metadata\"");
    }
}
