//! Validate snapshot settings against the currently-loaded session.
//!
//! See `openspec/specs/session-settings-io/spec.md` (Requirement:
//! `validate_against_inputs` SHALL identify input-dependent field
//! resets, and `apply_snapshot` SHALL overwrite settings with per-field
//! resets honored).

use crate::app::SessionSettings;
use crate::data::GroupMapping;
use crate::data::groups::UNASSIGNED;

/// The four input-dependent fields on `SessionSettings` that may need to
/// be reset to `None` at Load-apply time because the snapshot's value
/// does not exist in the current `GroupMapping`.
///
/// `Some(value)` means "the snapshot carried `Some(value)` and the
/// validator decided this field must be reset to `None` before
/// overwrite". The carried value is the SAVED value (used by the
/// confirm-modal renderer to display "Numerator group 'Treated' →
/// re-select required"). `None` means "no reset required for this
/// field" — either the snapshot had `None`, or its value is valid in
/// the current mapping.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct ValidationResets {
    pub numerator: Option<String>,
    pub denominator: Option<String>,
    pub metadata_column: Option<String>,
    pub pqn_reference_group: Option<String>,
}

impl ValidationResets {
    pub fn is_empty(&self) -> bool {
        self.numerator.is_none()
            && self.denominator.is_none()
            && self.metadata_column.is_none()
            && self.pqn_reference_group.is_none()
    }
}

/// Identify which of the four input-dependent fields on `settings`
/// reference values not present in `mapping`. Returns
/// `ValidationResets::default()` when `mapping == None` (no metadata
/// loaded yet — we have no reference data, so the snapshot's values are
/// accepted as-is; downstream gates catch problems later).
///
/// The group set used for `numerator`, `denominator`, and
/// `pqn_reference_group` excludes `UNASSIGNED` (matching how
/// `stage2_setup.rs` filters groups for the user-facing ComboBoxes).
/// The metadata column set is exactly `mapping.metadata_column_names()`
/// — already restricted to numeric columns by the loader.
pub fn validate_against_inputs(
    settings: &SessionSettings,
    mapping: Option<&GroupMapping>,
) -> ValidationResets {
    let Some(m) = mapping else {
        return ValidationResets::default();
    };

    let groups: Vec<String> = m.groups().into_iter().filter(|g| g != UNASSIGNED).collect();
    let in_groups = |name: &str| groups.iter().any(|g| g == name);

    let metadata_cols = m.metadata_column_names();
    let in_metadata_cols = |col: &str| metadata_cols.iter().any(|c| c == col);

    ValidationResets {
        numerator: settings
            .numerator
            .as_ref()
            .filter(|name| !in_groups(name))
            .cloned(),
        denominator: settings
            .denominator
            .as_ref()
            .filter(|name| !in_groups(name))
            .cloned(),
        metadata_column: settings
            .metadata_column
            .as_ref()
            .filter(|col| !in_metadata_cols(col))
            .cloned(),
        pqn_reference_group: settings
            .pqn_reference_group
            .as_ref()
            .filter(|g| !in_groups(g))
            .cloned(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::GroupMapping;
    use crate::data::groups::load_group_mapping;
    use std::io::Write;
    use tempfile::NamedTempFile;

    /// Build a minimal `GroupMapping` by writing a temp CSV and loading
    /// it via the public `load_group_mapping` — matches the unit-test
    /// pattern in `src/data/groups.rs`.
    fn mapping_with(groups: &[(&str, &str)], metadata_cols: &[(&str, Vec<f64>)]) -> GroupMapping {
        let mut header = String::from("sample,group");
        for (col, _) in metadata_cols {
            header.push(',');
            header.push_str(col);
        }
        header.push('\n');

        let mut body = String::new();
        for (i, (sample, group)) in groups.iter().enumerate() {
            body.push_str(sample);
            body.push(',');
            body.push_str(group);
            for (_, values) in metadata_cols {
                body.push(',');
                if let Some(v) = values.get(i) {
                    body.push_str(&v.to_string());
                }
            }
            body.push('\n');
        }

        let mut f = NamedTempFile::new().expect("create tempfile");
        f.write_all(header.as_bytes()).expect("write header");
        f.write_all(body.as_bytes()).expect("write body");

        let cols: Vec<String> = groups.iter().map(|(s, _)| s.to_string()).collect();
        load_group_mapping(f.path(), &cols).expect("load mapping")
    }

    fn baseline_settings() -> SessionSettings {
        SessionSettings::default()
    }

    #[test]
    fn no_mapping_yields_default_resets() {
        let mut s = baseline_settings();
        s.numerator = Some("Treated".to_string());
        s.metadata_column = Some("dry_weight".to_string());
        let r = validate_against_inputs(&s, None);
        assert_eq!(r, ValidationResets::default());
        assert!(r.is_empty());
    }

    #[test]
    fn stale_numerator_is_recorded() {
        let m = mapping_with(&[("s1", "A"), ("s2", "B")], &[]);
        let mut s = baseline_settings();
        s.numerator = Some("Treated".to_string());
        let r = validate_against_inputs(&s, Some(&m));
        assert_eq!(r.numerator, Some("Treated".to_string()));
        assert_eq!(r.denominator, None);
        assert_eq!(r.metadata_column, None);
        assert_eq!(r.pqn_reference_group, None);
        assert!(!r.is_empty());
    }

    #[test]
    fn valid_numerator_is_not_recorded() {
        let m = mapping_with(&[("s1", "A"), ("s2", "B")], &[]);
        let mut s = baseline_settings();
        s.numerator = Some("A".to_string());
        let r = validate_against_inputs(&s, Some(&m));
        assert_eq!(r.numerator, None);
        assert!(r.is_empty());
    }

    #[test]
    fn unassigned_is_excluded_from_groups_set() {
        // A mapping where "Unassigned" appears as a group label: the
        // validator must treat numerator="Unassigned" as invalid
        // because UNASSIGNED is filtered out of the group set.
        let m = mapping_with(&[("s1", "A"), ("s2", UNASSIGNED)], &[]);
        let mut s = baseline_settings();
        s.numerator = Some(UNASSIGNED.to_string());
        let r = validate_against_inputs(&s, Some(&m));
        assert_eq!(r.numerator, Some(UNASSIGNED.to_string()));
    }

    #[test]
    fn stale_metadata_column_is_recorded() {
        let m = mapping_with(
            &[("s1", "A"), ("s2", "B")],
            &[("dry_weight", vec![1.0, 2.0])],
        );
        let mut s = baseline_settings();
        s.metadata_column = Some("fresh_weight".to_string());
        let r = validate_against_inputs(&s, Some(&m));
        assert_eq!(r.metadata_column, Some("fresh_weight".to_string()));
    }

    #[test]
    fn valid_metadata_column_is_not_recorded() {
        let m = mapping_with(
            &[("s1", "A"), ("s2", "B")],
            &[("dry_weight", vec![1.0, 2.0])],
        );
        let mut s = baseline_settings();
        s.metadata_column = Some("dry_weight".to_string());
        let r = validate_against_inputs(&s, Some(&m));
        assert_eq!(r.metadata_column, None);
    }

    #[test]
    fn stale_pqn_reference_group_is_recorded() {
        let m = mapping_with(&[("s1", "A"), ("s2", "B")], &[]);
        let mut s = baseline_settings();
        s.pqn_reference_group = Some("MissingGroup".to_string());
        let r = validate_against_inputs(&s, Some(&m));
        assert_eq!(r.pqn_reference_group, Some("MissingGroup".to_string()));
    }

    #[test]
    fn all_four_resets_can_coexist() {
        let m = mapping_with(&[("s1", "A"), ("s2", "B")], &[]);
        let mut s = baseline_settings();
        s.numerator = Some("X".to_string());
        s.denominator = Some("Y".to_string());
        s.metadata_column = Some("Z".to_string());
        s.pqn_reference_group = Some("W".to_string());
        let r = validate_against_inputs(&s, Some(&m));
        assert_eq!(r.numerator, Some("X".to_string()));
        assert_eq!(r.denominator, Some("Y".to_string()));
        assert_eq!(r.metadata_column, Some("Z".to_string()));
        assert_eq!(r.pqn_reference_group, Some("W".to_string()));
        assert!(!r.is_empty());
    }
}
