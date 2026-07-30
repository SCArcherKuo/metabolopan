//! Regression test for the coverage route's load-bearing claim: **deduplication
//! cannot change any reported coverage number.**
//!
//! `add-kegg-coverage-route` design decision D16 keeps the deduplication control
//! on the coverage setup screen *only* because it moves the funnel and the
//! exported metabolite names. Everything else — `D`, every row's `entry_size` /
//! `hits` / `coverage` / `share` / `hit_compounds`, `detected_total`,
//! `detected_in_entries`, `entries_without_compounds` — is invariant under it,
//! and the setup screen's grey sub-hint tells the user so in as many words.
//!
//! If that claim ever stops holding, the UI text becomes a false statement about
//! the user's own data. This test is what stops that silently.
//!
//! The network half is deliberately excluded: `D` is the KEGG image of the
//! surviving InChIKey SET, so if the SET is invariant then `D` is too, whatever
//! the resolver returns. The test therefore substitutes a deterministic
//! InChIKey → cpd mapping and asserts on `coverage::compute`'s full output.

use std::collections::{HashMap, HashSet};

use metabolopan::coverage::detect::{detect_features, inchikeys_of};
use metabolopan::coverage::{CoverageResult, compute};
use metabolopan::data::types::MetabolomicsTable;
use metabolopan::data::{GroupMapping, load_group_mapping, parse_msdial_txt};
use metabolopan::kegg::types::KeggCompoundSet;

/// A deterministic stand-in for the PubChem → KEGG resolver: every InChIKey maps
/// to one synthetic cpd ID derived from the key itself.
///
/// Determinism is what makes the comparison meaningful — a real resolver's
/// network variability would mask the very drift this test looks for.
fn cpd_of(inchikey: &str) -> String {
    // The 14-character skeleton block is what actually identifies the compound;
    // using it means two InChIKeys differing only in protonation collapse the
    // same way a real resolver would.
    let skeleton: String = inchikey.chars().take(14).collect();
    format!("C{skeleton}")
}

fn detected_set(table: &MetabolomicsTable, kept: &[usize]) -> HashSet<String> {
    inchikeys_of(table, kept)
        .iter()
        .map(|k| cpd_of(k))
        .collect()
}

/// A synthetic catalogue built from the detected compounds themselves, plus
/// entries the data cannot hit and one empty entry.
///
/// Built from the DEDUP-ON set so the two runs are scored against an identical
/// catalogue — deriving it separately per run would let the catalogue absorb a
/// difference this test exists to detect.
fn catalogue(detected: &HashSet<String>) -> Vec<KeggCompoundSet> {
    let mut ids: Vec<&String> = detected.iter().collect();
    ids.sort();
    let mut entries = Vec::new();
    // Overlapping entries of assorted sizes, so `share` values genuinely exceed
    // 1.0 in sum and `coverage` spans the range.
    for (i, chunk) in ids.chunks(7).enumerate().take(30) {
        entries.push(KeggCompoundSet {
            id: format!("map{i:05}"),
            name: format!("Synthetic entry {i}"),
            compounds: chunk
                .iter()
                .map(|s| (*s).clone())
                // One unreachable compound per entry, so no entry is trivially
                // 100 % covered.
                .chain(std::iter::once(format!("CUNREACHABLE{i}")))
                .collect(),
        });
    }
    // A zero-compound entry — every KEGG global/overview map is one.
    entries.push(KeggCompoundSet {
        id: "map01100".to_string(),
        name: "Metabolic pathways (synthetic global map)".to_string(),
        compounds: vec![],
    });
    entries
}

/// Assert the two results agree on EVERY field, naming the first divergence.
fn assert_results_identical(on: &CoverageResult, off: &CoverageResult) {
    assert_eq!(
        on.detected_total, off.detected_total,
        "|D| moved under deduplication"
    );
    assert_eq!(on.entries_total, off.entries_total);
    assert_eq!(on.entries_without_compounds, off.entries_without_compounds);
    assert_eq!(
        on.detected_in_entries, off.detected_in_entries,
        "detected_in_entries moved under deduplication"
    );
    assert_eq!(on.rows.len(), off.rows.len());
    for (a, b) in on.rows.iter().zip(&off.rows) {
        assert_eq!(a.entry_id, b.entry_id, "row order diverged");
        assert_eq!(a.entry_size, b.entry_size, "{} entry_size", a.entry_id);
        assert_eq!(a.hits, b.hits, "{} hits", a.entry_id);
        assert_eq!(a.coverage, b.coverage, "{} coverage", a.entry_id);
        assert_eq!(a.share, b.share, "{} share", a.entry_id);
        assert_eq!(
            a.hit_compounds, b.hit_compounds,
            "{} hit_compounds",
            a.entry_id
        );
    }
    // `PartialEq` over the whole struct, as a backstop against a field added
    // later that the per-field assertions above forget.
    assert_eq!(
        on, off,
        "CoverageResult diverged in a field not named above"
    );
}

/// Run the full detect → compute path twice on one table and compare.
fn assert_invariant_on(
    table: &MetabolomicsTable,
    mapping: Option<&GroupMapping>,
    groups: &[String],
) {
    let (on, report_on) = detect_features(table, mapping, groups, 0.5, true, 0.1);
    let (off, report_off) = detect_features(table, mapping, groups, 0.5, false, 0.1);

    // Preconditions: the fixture must actually contain duplicates, or the test
    // would pass vacuously.
    assert!(
        report_on.is_some() && report_off.is_none(),
        "a report is produced iff dedup ran"
    );
    let dropped = report_on.as_ref().map(|r| r.dropped.len()).unwrap_or(0);
    assert!(
        dropped > 0,
        "fixture has no same-InChIKey duplicates — the test would pass vacuously"
    );
    assert!(
        on.after_dedup < off.after_dedup,
        "the funnel must MOVE: {} vs {}",
        on.after_dedup,
        off.after_dedup
    );
    // The group stage is unaffected by dedup, which runs after it.
    assert_eq!(on.in_selected_groups, off.in_selected_groups);

    // The claim itself.
    let d_on = detected_set(table, &on.kept);
    let d_off = detected_set(table, &off.kept);
    assert_eq!(d_on, d_off, "D moved under deduplication");

    let entries = catalogue(&d_on);
    let result_on = compute(&d_on, &entries);
    let result_off = compute(&d_off, &entries);
    assert_results_identical(&result_on, &result_off);

    // And the catalogue is non-trivial, so the comparison had something to say.
    assert!(
        result_on.rows.iter().any(|r| r.hits > 0),
        "no entry was hit — the catalogue is not exercising the comparison"
    );
    assert!(result_on.entries_without_compounds >= 1);
}

/// Single-mode, no metadata `.csv` — the group stage is inert, so this isolates
/// deduplication as the only filter that ran.
#[test]
fn dedup_never_moves_a_coverage_number_single_mode() {
    let txt = std::path::Path::new("data/single-mode/MS-DIAL-output-example.txt");
    if !txt.exists() {
        eprintln!("skipping: single-mode fixture absent (Git LFS not materialised)");
        return;
    }
    let table = parse_msdial_txt(txt).expect("parse single-mode fixture");
    assert_invariant_on(&table, None, &[]);
}

/// With a metadata `.csv` and a real group selection, so the group filter runs
/// FIRST and deduplication then operates on its survivors — the ordering the
/// route specifies. The invariance must survive that composition too.
#[test]
fn dedup_never_moves_a_coverage_number_with_a_group_filter() {
    let txt = std::path::Path::new("data/single-mode/MS-DIAL-output-example.txt");
    let csv = std::path::Path::new("data/single-mode/metadata-example.csv");
    if !txt.exists() || !csv.exists() {
        eprintln!("skipping: single-mode fixtures absent");
        return;
    }
    let table = parse_msdial_txt(txt).expect("parse single-mode fixture");
    let mapping = load_group_mapping(csv, &table.sample_cols).expect("load mapping");

    let groups: Vec<String> = mapping
        .groups()
        .into_iter()
        .filter(|g| g != metabolopan::data::UNASSIGNED)
        .collect();
    assert!(
        !groups.is_empty(),
        "fixture mapping has no assignable groups"
    );

    assert_invariant_on(&table, Some(&mapping), &groups);
}

/// The RT tolerance decides how many features represent an InChIKey, never which
/// InChIKeys exist — so every reported number is untouched by it as well.
#[test]
fn the_rt_tolerance_never_moves_a_coverage_number() {
    let txt = std::path::Path::new("data/single-mode/MS-DIAL-output-example.txt");
    if !txt.exists() {
        eprintln!("skipping: single-mode fixture absent");
        return;
    }
    let table = parse_msdial_txt(txt).expect("parse single-mode fixture");

    let (tight, _) = detect_features(&table, None, &[], 0.5, true, 0.1);
    let (loose, _) = detect_features(&table, None, &[], 0.5, true, 30.0);
    assert!(
        tight.after_dedup > loose.after_dedup,
        "a 30-minute tolerance must collapse more clusters than a 0.1-minute one \
         ({} vs {})",
        tight.after_dedup,
        loose.after_dedup
    );

    let d_tight = detected_set(&table, &tight.kept);
    let d_loose = detected_set(&table, &loose.kept);
    assert_eq!(d_tight, d_loose, "D moved with the RT tolerance");

    let entries = catalogue(&d_tight);
    assert_results_identical(&compute(&d_tight, &entries), &compute(&d_loose, &entries));
}

/// The number of resolver inputs is identical too, so the two runs cost the same
/// in PubChem and KEGG requests. `run_coverage` hands the resolver exactly this
/// list.
#[test]
fn dedup_does_not_change_the_resolver_workload() {
    let txt = std::path::Path::new("data/single-mode/MS-DIAL-output-example.txt");
    if !txt.exists() {
        eprintln!("skipping: single-mode fixture absent");
        return;
    }
    let table = parse_msdial_txt(txt).expect("parse single-mode fixture");
    let (on, _) = detect_features(&table, None, &[], 0.5, true, 0.1);
    let (off, _) = detect_features(&table, None, &[], 0.5, false, 0.1);

    let keys_on: HashSet<String> = inchikeys_of(&table, &on.kept).into_iter().collect();
    let keys_off: HashSet<String> = inchikeys_of(&table, &off.kept).into_iter().collect();
    assert_eq!(
        keys_on, keys_off,
        "the InChIKey list handed to the resolver moved"
    );
}

/// The other half of D16: deduplication DOES change which MS-DIAL metabolite
/// names the CSV attaches to each compound. Without this the control would be
/// entirely inert and there would be no reason to keep it.
#[test]
fn dedup_does_change_the_exported_metabolite_names() {
    let txt = std::path::Path::new("data/single-mode/MS-DIAL-output-example.txt");
    if !txt.exists() {
        eprintln!("skipping: single-mode fixture absent");
        return;
    }
    let table = parse_msdial_txt(txt).expect("parse single-mode fixture");
    let (on, _) = detect_features(&table, None, &[], 0.5, true, 0.1);
    let (off, _) = detect_features(&table, None, &[], 0.5, false, 0.1);

    let names_of = |kept: &[usize]| -> HashMap<String, HashSet<String>> {
        let mut m: HashMap<String, HashSet<String>> = HashMap::new();
        for &i in kept {
            let f = &table.features[i];
            if let Some(k) = f.inchikey.as_deref() {
                m.entry(cpd_of(k))
                    .or_default()
                    .insert(f.metabolite_name.clone());
            }
        }
        m
    };
    let names_on = names_of(&on.kept);
    let names_off = names_of(&off.kept);

    assert_eq!(
        names_on.keys().collect::<HashSet<_>>(),
        names_off.keys().collect::<HashSet<_>>(),
        "the compound KEYS must match — only the name lists may differ"
    );
    assert_ne!(
        names_on, names_off,
        "deduplication must change at least one compound's name list, or the \
         control has no observable effect at all and D16's justification fails"
    );
}
