//! Building the coverage route's detected feature set: the group-presence
//! filter, then deduplication, per ion-mode table.
//!
//! Pure — no I/O, no network, no `tracing` from these bodies (the orchestrator
//! logs around the calls). Owner: the `kegg-coverage` capability.

use std::collections::HashSet;

use crate::data::GroupMapping;
use crate::data::groups::UNASSIGNED;
use crate::data::types::MetabolomicsTable;
use crate::dedup::{DedupReport, run_dedup};

/// Resolve `settings.coverage_selected_groups` against a mapping.
///
/// **`None` MUST NOT be resolved with `Option::unwrap_or_default()`.** That maps
/// `None` to the EMPTY list, which is the *opposite* of what it means: the OR
/// across selected groups would range over nothing, no feature would survive,
/// and `D` would come back empty. `None` means "not yet chosen ⇒ use every
/// group"; only `Some(vec![])` means "the user deliberately chose none".
///
/// The setup screen normally replaces `None` with `Some(<all groups>)` before a
/// run can start, so `None` should not reach here on the ordinary path — but
/// `apply_snapshot` deliberately writes `None` on a stale-selection reset, so
/// this fallback has to be correct on its own rather than relying on a render
/// having happened first.
///
/// `Unassigned` is never offerable as a selectable group, so it is filtered out
/// of the "use everything" fallback. A caller that explicitly passes it in
/// `Some(list)` gets it back — but no feature can be present in it, because
/// `group_of` returning `UNASSIGNED` is what "this sample has no group" means
/// and the setup screen never offers it.
pub fn selected_groups(selection: Option<&[String]>, mapping: &GroupMapping) -> Vec<String> {
    match selection {
        Some(list) => list.to_vec(),
        None => mapping
            .groups()
            .into_iter()
            .filter(|g| g != UNASSIGNED)
            .collect(),
    }
}

/// Is feature `f` present in the group whose columns are `cols` (indices into
/// `table.sample_cols`)?
///
/// "Has signal" is `is_finite() && > 0.0`, **not** merely "not NaN". The MS-DIAL
/// parser maps empty / `NA` / `null` cells to `NaN` but parses a literal `0` as
/// `0.0`, and MS-DIAL writes literal zeros for not-detected. Across all three
/// bundled fixtures not one sample cell is empty and about 8 % are literal `0`,
/// so a NaN-only test would be a no-op on real data.
///
/// The `count >= 1` clause is an unconditional floor, so a `threshold` of `0.0`
/// degrades to "present in at least one sample of the group" rather than
/// "vacuously present everywhere".
///
/// Reads `intensity_raw`, the immutable as-loaded matrix — never `intensity`.
/// No sample normalization is offered on the coverage route, so the working
/// matrix is never populated here; reading it would read whatever a previous
/// DAM run in the same session happened to leave behind.
fn present_in_group(
    table: &MetabolomicsTable,
    feature: usize,
    cols: &[usize],
    threshold: f64,
) -> bool {
    if cols.is_empty() {
        return false;
    }
    let count = cols
        .iter()
        .filter(|&&c| {
            let v = table.intensity_raw[[feature, c]];
            v.is_finite() && v > 0.0
        })
        .count();
    count >= 1 && (count as f64 / cols.len() as f64) >= threshold
}

/// Column indices of each selected group within THIS table's `sample_cols`.
///
/// A group with no columns in this table yields an empty vec, which
/// [`present_in_group`] rejects — so in dual mode a group that exists only in
/// POS simply never votes for a NEG feature, and NEG features can still survive
/// through any other selected group.
fn group_columns(
    table: &MetabolomicsTable,
    mapping: &GroupMapping,
    groups: &[String],
) -> Vec<Vec<usize>> {
    groups
        .iter()
        .map(|g| {
            table
                .sample_cols
                .iter()
                .enumerate()
                .filter(|(_, name)| mapping.group_of(name) == g)
                .map(|(idx, _)| idx)
                .collect()
        })
        .collect()
}

/// Indices of the features present in AT LEAST ONE selected group.
///
/// A plain OR across the selection: a compound seen in any one selected
/// condition is a compound the sample contains.
pub fn group_presence_survivors(
    table: &MetabolomicsTable,
    mapping: &GroupMapping,
    groups: &[String],
    threshold: f64,
) -> Vec<usize> {
    let cols_per_group = group_columns(table, mapping, groups);
    (0..table.features.len())
        .filter(|&f| {
            cols_per_group
                .iter()
                .any(|cols| present_in_group(table, f, cols, threshold))
        })
        .collect()
}

/// One ion-mode table's surviving features plus the funnel counts the two
/// filters produce.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct DetectedFeatures {
    /// Indices into `table.features`, ascending.
    pub kept: Vec<usize>,
    /// `table.features.len()`.
    pub raw_features: usize,
    /// Survivors of the group-presence filter. `None` when no metadata `.csv`
    /// was supplied — the stage did not run, so the renderer omits the term
    /// rather than printing a tautology.
    pub in_selected_groups: Option<usize>,
    /// `kept.len()`. Equals `in_selected_groups` (or `raw_features`) when
    /// deduplication is off.
    pub after_dedup: usize,
}

/// Build one ion-mode table's detected feature set: group-presence filter, then
/// deduplication.
///
/// **The group filter MUST precede deduplication.** The cascade picks one
/// survivor per InChIKey retention-time cluster on annotation quality alone; if
/// it ran first, an InChIKey whose highest-quality feature happens to appear
/// only in an unselected group (a QC pool, a solvent blank) would elect that
/// feature as the cluster's representative and then lose the whole compound at
/// the group stage — even though a perfectly good feature for it exists in the
/// selected groups. Filtering first lets the cascade choose a representative
/// from among the features that actually exist in the samples the user cares
/// about.
///
/// `mapping == None` (no metadata `.csv`) skips the group stage entirely: every
/// feature passes and no intensity value is read.
///
/// Features with no InChIKey are NOT removed here — they are excluded
/// structurally downstream, since `D` is the KEGG image of the surviving
/// InChIKeys and a feature without one contributes nothing to it.
/// `settings.drop_unknown` is never read on this route: offering a checkbox for
/// something that is not a choice would misrepresent it as one.
///
/// Returns the `DedupReport` alongside, for the Data tab's `Dedupe:` line and
/// audit download — on this route that report is the only surface on which
/// deduplication's effect is observable at all, since it provably cannot move
/// any reported coverage number.
pub fn detect_features(
    table: &MetabolomicsTable,
    mapping: Option<&GroupMapping>,
    groups: &[String],
    threshold: f64,
    dedup_enabled: bool,
    dedup_rt_tolerance_min: f64,
) -> (DetectedFeatures, Option<DedupReport>) {
    let raw_features = table.features.len();

    // ── 1. Group-presence filter ──
    let (surviving, in_selected_groups) = match mapping {
        Some(m) => {
            let s = group_presence_survivors(table, m, groups, threshold);
            let n = s.len();
            (s, Some(n))
        }
        None => ((0..raw_features).collect::<Vec<_>>(), None),
    };

    // ── 2. Deduplication ──
    let (kept, report) = if dedup_enabled {
        // `run_dedup` needs a contiguous slice, and it must see ONLY the group
        // survivors — that ordering is the whole point of this function. The
        // clone is a one-off over an already-filtered set inside a
        // multi-minute network-bound run.
        let subset: Vec<_> = surviving
            .iter()
            .map(|&i| table.features[i].clone())
            .collect();
        let (keep_local, report) = run_dedup(&subset, dedup_rt_tolerance_min);
        // `keep_local` indexes `subset`; map back to table indices, ascending.
        let mut kept: Vec<usize> = keep_local.iter().map(|&l| surviving[l]).collect();
        kept.sort_unstable();
        (kept, Some(report))
    } else {
        (surviving, None)
    };

    let after_dedup = kept.len();
    (
        DetectedFeatures {
            kept,
            raw_features,
            in_selected_groups,
            after_dedup,
        },
        report,
    )
}

/// The distinct InChIKeys of the kept features, in first-appearance order.
///
/// First-appearance rather than sorted so the resolver's request order tracks
/// the table's own order, which keeps the progress strip's `current` value
/// moving through the file the way the user would read it.
pub fn inchikeys_of(table: &MetabolomicsTable, kept: &[usize]) -> Vec<String> {
    let mut seen: HashSet<&str> = HashSet::new();
    let mut out = Vec::new();
    for &i in kept {
        if let Some(k) = table.features[i].inchikey.as_deref()
            && seen.insert(k)
        {
            out.push(k.to_string());
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::types::FeatureMeta;
    use ndarray::Array2;

    fn feature(id: &str, inchikey: Option<&str>) -> FeatureMeta {
        FeatureMeta {
            alignment_id: id.to_string(),
            metabolite_name: format!("M{id}"),
            inchikey: inchikey.map(str::to_string),
            adduct_type: None,
            average_rt_min: Some(1.0),
            average_mz: None,
            formula: None,
            smiles: None,
            fill_percent: None,
            ms_ms_matched: None,
            isotope_tracking_weight_number: None,
            total_score: None,
            sn_average: None,
        }
    }

    /// A table with the given per-feature intensity rows over `sample_cols`.
    fn table(
        sample_cols: &[&str],
        features: Vec<FeatureMeta>,
        rows: &[&[f64]],
    ) -> MetabolomicsTable {
        let n_f = features.len();
        let n_s = sample_cols.len();
        let flat: Vec<f64> = rows.iter().flat_map(|r| r.iter().copied()).collect();
        let raw = Array2::from_shape_vec((n_f, n_s), flat).expect("shape matches");
        let annotated_count = features.iter().filter(|f| f.inchikey.is_some()).count();
        MetabolomicsTable {
            features,
            sample_cols: sample_cols.iter().map(|s| s.to_string()).collect(),
            // `intensity` is deliberately DIFFERENT from `intensity_raw` so a
            // test reading the wrong matrix fails loudly.
            intensity: Array2::from_elem((n_f, n_s), f64::NAN),
            intensity_raw: raw,
            excluded_cols: vec![],
            annotated_count,
        }
    }

    fn mapping_from(rows: &[(&str, &str)], sample_cols: &[&str]) -> GroupMapping {
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().expect("tempfile");
        writeln!(f, "sample,group").unwrap();
        for (s, g) in rows {
            writeln!(f, "{s},{g}").unwrap();
        }
        let cols: Vec<String> = sample_cols.iter().map(|s| s.to_string()).collect();
        crate::data::groups::load_group_mapping(f.path(), &cols).expect("mapping loads")
    }

    const CTRL6: [&str; 6] = ["c1", "c2", "c3", "c4", "c5", "c6"];

    fn ctrl_mapping() -> GroupMapping {
        mapping_from(
            &[
                ("c1", "Control"),
                ("c2", "Control"),
                ("c3", "Control"),
                ("c4", "Control"),
                ("c5", "Control"),
                ("c6", "Control"),
            ],
            &CTRL6,
        )
    }

    fn survives(row: &[f64], threshold: f64) -> bool {
        let t = table(&CTRL6, vec![feature("1", Some("K1"))], &[row]);
        let m = ctrl_mapping();
        !group_presence_survivors(&t, &m, &["Control".to_string()], threshold).is_empty()
    }

    /// 3 of 6 at a 0.5 threshold is exactly at the boundary — and `>=` means it
    /// is present.
    #[test]
    fn exactly_at_the_threshold_is_present() {
        assert!(survives(&[0.0, 0.0, 15200.0, 18900.0, 0.0, 21000.0], 0.5));
    }

    /// 2 of 6 ≈ 0.333 is below 0.5.
    #[test]
    fn just_below_the_threshold_is_absent() {
        assert!(!survives(&[0.0, 0.0, 0.0, 18900.0, 0.0, 21000.0], 0.5));
    }

    /// Literal zeros are absences, not values. A NaN-only missing-value test
    /// would wrongly have admitted this feature — and ~8 % of real cells are
    /// literal zeros, so that mistake would be invisible in review and enormous
    /// in effect.
    #[test]
    fn literal_zeros_are_absences() {
        assert!(!survives(&[0.0; 6], 0.5));
        assert!(!survives(&[0.0; 6], 0.0));
    }

    /// NaN cells count against presence exactly as zeros do.
    #[test]
    fn nan_cells_are_absences_too() {
        let nan = f64::NAN;
        // Three real values out of six: present.
        assert!(survives(&[nan, nan, 15200.0, 18900.0, 0.0, 21000.0], 0.5));
        // Only two: absent.
        assert!(!survives(&[nan, nan, nan, nan, 15200.0, 21000.0], 0.5));
    }

    /// A threshold of 0.0 degrades to "at least one sample", NOT to "everything
    /// passes" — the `count >= 1` floor is unconditional.
    #[test]
    fn a_zero_threshold_still_requires_one_sample() {
        assert!(survives(&[0.0, 0.0, 0.0, 0.0, 0.0, 42.0], 0.0));
        assert!(!survives(&[0.0, 0.0, 0.0, 0.0, 0.0, 0.0], 0.0));
        assert!(!survives(&[f64::NAN; 6], 0.0));
    }

    /// Presence in ANY one selected group is enough.
    #[test]
    fn presence_in_any_one_selected_group_is_enough() {
        let cols = ["c1", "c2", "t1", "t2"];
        let m = mapping_from(
            &[
                ("c1", "Control"),
                ("c2", "Control"),
                ("t1", "Treated"),
                ("t2", "Treated"),
            ],
            &cols,
        );
        // Absent from Control, present in Treated.
        let t = table(
            &cols,
            vec![feature("1", Some("K1"))],
            &[&[0.0, 0.0, 500.0, 900.0]],
        );
        let groups = vec!["Control".to_string(), "Treated".to_string()];
        assert_eq!(group_presence_survivors(&t, &m, &groups, 0.5), vec![0]);
        // With only Control selected it is gone.
        assert!(group_presence_survivors(&t, &m, &["Control".to_string()], 0.5).is_empty());
    }

    /// A group with no columns in THIS table cannot veto it: `cols` is empty, so
    /// `present` is false for that group, and the feature survives through any
    /// other selected group that does have columns here.
    #[test]
    fn a_group_absent_from_this_table_does_not_veto_it() {
        let cols = ["c1", "c2"];
        // `Treated` exists in the mapping but has no column in this table.
        let m = mapping_from(
            &[("c1", "Control"), ("c2", "Control"), ("t1", "Treated")],
            &["c1", "c2", "t1"],
        );
        let t = table(&cols, vec![feature("1", Some("K1"))], &[&[500.0, 900.0]]);
        let groups = vec!["Control".to_string(), "Treated".to_string()];
        assert_eq!(group_presence_survivors(&t, &m, &groups, 0.5), vec![0]);
        // And a table where ONLY the absent group is selected keeps nothing.
        assert!(group_presence_survivors(&t, &m, &["Treated".to_string()], 0.5).is_empty());
    }

    /// The filter reads `intensity_raw`, so an earlier DAM run's normalization
    /// left in `intensity` cannot affect it. The fixture builder fills
    /// `intensity` with NaN, which would make every feature absent if read.
    #[test]
    fn the_raw_matrix_is_the_one_read() {
        let t = table(
            &CTRL6,
            vec![feature("1", Some("K1"))],
            &[&[100.0, 200.0, 300.0, 400.0, 500.0, 600.0]],
        );
        assert!(t.intensity.iter().all(|v| v.is_nan()));
        let m = ctrl_mapping();
        assert_eq!(
            group_presence_survivors(&t, &m, &["Control".to_string()], 0.5),
            vec![0]
        );
    }

    /// `None` means "use every non-`Unassigned` group", NOT the empty list.
    /// `unwrap_or_default()` here would silently empty `D`.
    #[test]
    fn a_none_selection_resolves_to_every_group() {
        let m = mapping_from(
            &[("c1", "Control"), ("t1", "Treated"), ("x1", UNASSIGNED)],
            &["c1", "t1", "x1"],
        );
        let mut got = selected_groups(None, &m);
        got.sort();
        assert_eq!(got, vec!["Control".to_string(), "Treated".to_string()]);
        assert!(!got.iter().any(|g| g == UNASSIGNED));
    }

    /// `Some(vec![])` is honoured verbatim — it means "deliberately none", which
    /// is a different state from "not yet chosen" and must stay distinguishable.
    #[test]
    fn an_empty_selection_is_honoured_verbatim() {
        let m = mapping_from(&[("c1", "Control")], &["c1"]);
        assert!(selected_groups(Some(&[]), &m).is_empty());

        let t = table(&["c1"], vec![feature("1", Some("K1"))], &[&[500.0]]);
        assert!(group_presence_survivors(&t, &m, &[], 0.5).is_empty());
    }

    /// **The case the filter-before-dedup order exists for.**
    ///
    /// One InChIKey, two features in the same RT cluster. Feature A has the
    /// higher `total_score` but sits only in the unselected `QC` group; feature
    /// B sits in the selected `Control` group. Filtering first lets B become the
    /// cluster's survivor, so the compound reaches `D`. Had the cascade run
    /// first, A would have won the cluster and then been dropped by the group
    /// filter — losing the compound entirely.
    #[test]
    fn the_group_filter_runs_before_deduplication() {
        let cols = ["c1", "c2", "q1", "q2"];
        let m = mapping_from(
            &[
                ("c1", "Control"),
                ("c2", "Control"),
                ("q1", "QC"),
                ("q2", "QC"),
            ],
            &cols,
        );
        let mut a = feature("A", Some("SAMEKEY"));
        a.total_score = Some(99.0);
        let mut b = feature("B", Some("SAMEKEY"));
        b.total_score = Some(10.0);
        // A: QC only. B: Control only. Same RT, so one cluster.
        let t = table(
            &cols,
            vec![a, b],
            &[&[0.0, 0.0, 800.0, 900.0], &[500.0, 600.0, 0.0, 0.0]],
        );

        let selected = vec!["Control".to_string()];
        let (detected, report) = detect_features(&t, Some(&m), &selected, 0.5, true, 0.1);

        assert_eq!(detected.in_selected_groups, Some(1), "only B passes step 1");
        assert_eq!(
            detected.kept,
            vec![1],
            "B survives as the cluster's champion"
        );
        assert_eq!(
            inchikeys_of(&t, &detected.kept),
            vec!["SAMEKEY".to_string()],
            "the compound reaches D through B"
        );
        assert!(report.is_some());

        // The counterfactual: cascade first over the WHOLE table elects A…
        let (keep_all, _) = run_dedup(&t.features, 0.1);
        assert_eq!(
            keep_all.iter().copied().collect::<Vec<_>>(),
            vec![0],
            "A wins the cluster on total_score when the group filter has not run"
        );
        // …and A is then removed by the group filter, losing the compound.
        assert!(
            !group_presence_survivors(&t, &m, &selected, 0.5).contains(&0),
            "A does not survive the group filter, so dedup-first loses SAMEKEY"
        );
    }

    /// With no `.csv`, step 1 is inert: every feature passes, the funnel term is
    /// `None`, and no intensity is read.
    #[test]
    fn no_mapping_skips_the_group_stage() {
        let t = table(
            &CTRL6,
            vec![feature("1", Some("K1")), feature("2", None)],
            &[&[0.0; 6], &[0.0; 6]],
        );
        let (detected, report) = detect_features(&t, None, &[], 0.5, false, 0.1);
        assert_eq!(detected.raw_features, 2);
        assert_eq!(detected.in_selected_groups, None);
        assert_eq!(detected.kept, vec![0, 1], "every feature passes step 1");
        assert!(report.is_none());
    }

    /// Deduplication off keeps every group survivor, and no report is produced.
    #[test]
    fn dedup_off_keeps_every_group_survivor() {
        let cols = ["c1", "c2"];
        let m = mapping_from(&[("c1", "Control"), ("c2", "Control")], &cols);
        let t = table(
            &cols,
            vec![feature("A", Some("SAMEKEY")), feature("B", Some("SAMEKEY"))],
            &[&[500.0, 600.0], &[700.0, 800.0]],
        );
        let selected = vec!["Control".to_string()];

        let (off, report_off) = detect_features(&t, Some(&m), &selected, 0.5, false, 0.1);
        assert_eq!(off.kept, vec![0, 1]);
        assert!(report_off.is_none());

        let (on, report_on) = detect_features(&t, Some(&m), &selected, 0.5, true, 0.1);
        assert_eq!(on.kept.len(), 1, "one champion per RT cluster");
        assert!(report_on.is_some());

        // The funnel MOVES…
        assert_ne!(off.after_dedup, on.after_dedup);
        // …but the InChIKey set — and therefore D — does not.
        assert_eq!(inchikeys_of(&t, &off.kept), inchikeys_of(&t, &on.kept));
    }

    /// Features with no InChIKey are not filtered out here; they simply
    /// contribute nothing to `D`. `drop_unknown` is never consulted.
    #[test]
    fn unannotated_features_pass_through_and_contribute_nothing() {
        let cols = ["c1", "c2"];
        let m = mapping_from(&[("c1", "Control"), ("c2", "Control")], &cols);
        let t = table(
            &cols,
            vec![feature("A", Some("K1")), feature("B", None)],
            &[&[500.0, 600.0], &[700.0, 800.0]],
        );
        let (detected, _) = detect_features(&t, Some(&m), &["Control".to_string()], 0.5, true, 0.1);
        assert_eq!(detected.kept, vec![0, 1], "both survive the two filters");
        assert_eq!(inchikeys_of(&t, &detected.kept), vec!["K1".to_string()]);
    }

    /// **T6-C.** With no `.csv` the group filter is inert, so the coverage
    /// route's surviving feature set must equal exactly the set `run_dam`
    /// iterates over for the same dedup settings.
    ///
    /// `run_dam`'s loop skips on three conditions: the dedup mask, then
    /// `drop_unknown`, then the pre-filter. The coverage route never reads
    /// `drop_unknown` and runs no pre-filter (both are two-group-comparison
    /// concerns), so with those two neutral the two sets must coincide — and
    /// that is exactly what "the coverage route reuses the DAM route's
    /// deduplication" has to mean to be worth saying.
    #[test]
    fn with_no_csv_the_kept_set_equals_run_dams_dedup_mask() {
        let cols = ["c1", "c2"];
        let mut a = feature("A", Some("K1"));
        a.total_score = Some(90.0);
        let mut b = feature("B", Some("K1"));
        b.total_score = Some(10.0);
        let mut c = feature("C", Some("K2"));
        c.average_rt_min = Some(9.0);
        let t = table(
            &cols,
            vec![a, b, c, feature("D", None)],
            &[&[1.0, 2.0], &[3.0, 4.0], &[5.0, 6.0], &[7.0, 8.0]],
        );

        for tol in [0.1, 1.0, 5.0] {
            let (detected, _) = detect_features(&t, None, &[], 0.5, true, tol);
            let (dam_mask, _) = run_dedup(&t.features, tol);
            let mut expected: Vec<usize> = dam_mask.into_iter().collect();
            expected.sort_unstable();
            assert_eq!(
                detected.kept, expected,
                "coverage and run_dam must iterate the same features at tol={tol}"
            );
        }
    }

    /// Changing the RT tolerance moves which features represent an InChIKey,
    /// never which InChIKeys exist — so every coverage number is untouched.
    #[test]
    fn the_rt_tolerance_never_changes_the_inchikey_set() {
        let cols = ["c1"];
        let mut a = feature("A", Some("K1"));
        a.average_rt_min = Some(1.0);
        a.total_score = Some(90.0);
        let mut b = feature("B", Some("K1"));
        b.average_rt_min = Some(4.0);
        b.total_score = Some(10.0);
        let t = table(&cols, vec![a, b], &[&[1.0], &[2.0]]);

        let (tight, _) = detect_features(&t, None, &[], 0.5, true, 0.1);
        let (loose, _) = detect_features(&t, None, &[], 0.5, true, 5.0);
        assert_ne!(tight.after_dedup, loose.after_dedup, "the funnel moves");
        assert_eq!(
            inchikeys_of(&t, &tight.kept),
            inchikeys_of(&t, &loose.kept),
            "the InChIKey set — and therefore D — does not"
        );
    }

    /// `inchikeys_of` deduplicates in first-appearance order.
    #[test]
    fn inchikeys_are_distinct_in_first_appearance_order() {
        let cols = ["c1"];
        let t = table(
            &cols,
            vec![
                feature("A", Some("ZZZ")),
                feature("B", None),
                feature("C", Some("AAA")),
                feature("D", Some("ZZZ")),
            ],
            &[&[1.0], &[1.0], &[1.0], &[1.0]],
        );
        assert_eq!(
            inchikeys_of(&t, &[0, 1, 2, 3]),
            vec!["ZZZ".to_string(), "AAA".to_string()]
        );
    }
}
