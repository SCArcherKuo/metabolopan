//! Validate snapshot settings against the currently-loaded session.
//!
//! See the `session-settings-io` capability spec (Requirement:
//! `validate_against_inputs` SHALL identify input-dependent field
//! resets, and `apply_snapshot` SHALL overwrite settings with per-field
//! resets honored).

use crate::app::{AnalysisRoute, SessionSettings};
use crate::data::GroupMapping;
use crate::data::groups::UNASSIGNED;

/// The input-dependent fields on `SessionSettings` that may need to be reset at
/// Load-apply time because the snapshot's value is not valid in the current
/// session.
///
/// `Some(value)` means "the snapshot carried `Some(value)` and the
/// validator decided this field must be reset to `None` before
/// overwrite". The carried value is the SAVED value (used by the
/// confirm-modal renderer to display "Numerator group 'Treated' →
/// re-select required"). `None` means "no reset required for this
/// field" — either the snapshot had `None`, or its value is valid in
/// the current mapping.
///
/// [`ValidationResets::analysis_route`] is the exception to all of that: see
/// its own doc comment.
#[derive(Debug, Default, Clone, PartialEq, Eq)]
pub struct ValidationResets {
    pub numerator: Option<String>,
    pub denominator: Option<String>,
    pub metadata_column: Option<String>,
    pub pqn_reference_group: Option<String>,
    /// The saved coverage group selection, when ANY entry of it names a group
    /// absent from the current mapping. The whole list is recorded, not just
    /// the stale entries, so the confirm modal can show what the user chose.
    pub coverage_selected_groups: Option<Vec<String>>,
    /// The INCOMING route, when it differs from the session's current one.
    ///
    /// Not a field reset — `apply_snapshot` applies the incoming route
    /// verbatim, because a snapshot that says "coverage route" means it and
    /// discarding it would make a saved coverage session unrestorable. The
    /// entry exists so the confirm modal can warn that loading this file
    /// switches routes, and so the apply site knows it must also repair
    /// navigation (it will be sitting on a screen belonging to the route it
    /// just left).
    pub analysis_route: Option<AnalysisRoute>,
}

impl ValidationResets {
    pub fn is_empty(&self) -> bool {
        self.numerator.is_none()
            && self.denominator.is_none()
            && self.metadata_column.is_none()
            && self.pqn_reference_group.is_none()
            && self.coverage_selected_groups.is_none()
            && self.analysis_route.is_none()
    }
}

/// Identify which input-dependent fields on `settings` are not valid in the
/// current session.
///
/// The group set used for `numerator`, `denominator`, `pqn_reference_group`,
/// and `coverage_selected_groups` excludes `UNASSIGNED` (matching how
/// `stage2_setup.rs` filters groups for the user-facing ComboBoxes).
/// The metadata column set is exactly `mapping.metadata_column_names()`
/// — already restricted to numeric columns by the loader.
///
/// `mapping == None` (no metadata loaded yet) suppresses every mapping-derived
/// check — there is no reference data, so the snapshot's values are accepted
/// as-is and downstream gates handle any later mismatch. The `analysis_route`
/// check still runs: it validates against the session's own route, not against
/// the mapping.
///
/// Pure and side-effect-free.
pub fn validate_against_inputs(
    settings: &SessionSettings,
    mapping: Option<&GroupMapping>,
    current_route: AnalysisRoute,
) -> ValidationResets {
    // Runs whether or not a mapping exists.
    let analysis_route =
        (settings.analysis_route != current_route).then_some(settings.analysis_route);

    let Some(m) = mapping else {
        return ValidationResets {
            analysis_route,
            ..Default::default()
        };
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
        // ANY stale entry invalidates the WHOLE selection rather than silently
        // keeping the valid subset: a partly-recognised selection is a
        // selection the user never made, so make them re-confirm it. An empty
        // list is never invalid — it names no group.
        coverage_selected_groups: settings
            .coverage_selected_groups
            .as_ref()
            .filter(|list| list.iter().any(|g| !in_groups(g)))
            .cloned(),
        analysis_route,
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
        let r = validate_against_inputs(&s, None, AnalysisRoute::DamEnrichment);
        assert_eq!(r, ValidationResets::default());
        assert!(r.is_empty());
    }

    #[test]
    fn stale_numerator_is_recorded() {
        let m = mapping_with(&[("s1", "A"), ("s2", "B")], &[]);
        let mut s = baseline_settings();
        s.numerator = Some("Treated".to_string());
        let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
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
        let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
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
        let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
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
        let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
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
        let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
        assert_eq!(r.metadata_column, None);
    }

    #[test]
    fn stale_pqn_reference_group_is_recorded() {
        let m = mapping_with(&[("s1", "A"), ("s2", "B")], &[]);
        let mut s = baseline_settings();
        s.pqn_reference_group = Some("MissingGroup".to_string());
        let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
        assert_eq!(r.pqn_reference_group, Some("MissingGroup".to_string()));
    }

    /// A stale group in the coverage selection invalidates the WHOLE list, and
    /// the recorded value is the saved list so the modal can name it.
    #[test]
    fn stale_coverage_selected_groups_are_recorded_whole() {
        let m = mapping_with(&[("s1", "A"), ("s2", "B")], &[]);
        let mut s = baseline_settings();
        s.coverage_selected_groups = Some(vec!["Control".to_string(), "Treated".to_string()]);
        let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
        assert_eq!(
            r.coverage_selected_groups,
            Some(vec!["Control".to_string(), "Treated".to_string()])
        );

        // One valid + one stale is still wholly invalid: a partly-recognised
        // selection is a selection the user never made.
        s.coverage_selected_groups = Some(vec!["A".to_string(), "Treated".to_string()]);
        let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
        assert_eq!(
            r.coverage_selected_groups,
            Some(vec!["A".to_string(), "Treated".to_string()])
        );
    }

    /// A fully-valid selection, `None`, and `Some(vec![])` all pass: neither
    /// sentinel names a group, so neither can name a stale one.
    #[test]
    fn valid_or_empty_coverage_selection_is_not_recorded() {
        let m = mapping_with(&[("s1", "A"), ("s2", "B")], &[]);
        let mut s = baseline_settings();
        for sel in [
            Some(vec!["A".to_string(), "B".to_string()]),
            Some(vec![]),
            None,
        ] {
            s.coverage_selected_groups = sel.clone();
            let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
            assert_eq!(r.coverage_selected_groups, None, "for selection {sel:?}");
        }
    }

    /// `Unassigned` is never offerable, so a selection naming it is stale even
    /// though the mapping "has" that label.
    #[test]
    fn unassigned_in_the_coverage_selection_is_stale() {
        let m = mapping_with(&[("s1", "A"), ("s2", UNASSIGNED)], &[]);
        let mut s = baseline_settings();
        s.coverage_selected_groups = Some(vec![UNASSIGNED.to_string()]);
        let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
        assert_eq!(
            r.coverage_selected_groups,
            Some(vec![UNASSIGNED.to_string()])
        );
    }

    /// The route check does not depend on the mapping, so it still fires when
    /// `mapping == None` suppresses every other check. The recorded value is
    /// the INCOMING route.
    #[test]
    fn a_route_change_is_recorded_even_with_no_mapping() {
        let mut s = baseline_settings();
        s.analysis_route = AnalysisRoute::KeggCoverage;
        s.numerator = Some("Treated".to_string());

        let r = validate_against_inputs(&s, None, AnalysisRoute::DamEnrichment);
        assert_eq!(r.analysis_route, Some(AnalysisRoute::KeggCoverage));
        assert!(!r.is_empty());
        // Every mapping-derived check stayed suppressed.
        assert_eq!(r.numerator, None);
        assert_eq!(r.coverage_selected_groups, None);
    }

    /// Same route in and out: nothing to warn about, nothing to repair.
    #[test]
    fn a_matching_route_records_nothing() {
        let mut s = baseline_settings();
        s.analysis_route = AnalysisRoute::KeggCoverage;
        let r = validate_against_inputs(&s, None, AnalysisRoute::KeggCoverage);
        assert_eq!(r.analysis_route, None);
        assert!(r.is_empty());
    }

    #[test]
    fn all_four_resets_can_coexist() {
        let m = mapping_with(&[("s1", "A"), ("s2", "B")], &[]);
        let mut s = baseline_settings();
        s.numerator = Some("X".to_string());
        s.denominator = Some("Y".to_string());
        s.metadata_column = Some("Z".to_string());
        s.pqn_reference_group = Some("W".to_string());
        let r = validate_against_inputs(&s, Some(&m), AnalysisRoute::DamEnrichment);
        assert_eq!(r.numerator, Some("X".to_string()));
        assert_eq!(r.denominator, Some("Y".to_string()));
        assert_eq!(r.metadata_column, Some("Z".to_string()));
        assert_eq!(r.pqn_reference_group, Some("W".to_string()));
        assert!(!r.is_empty());
    }
}
