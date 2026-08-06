//! KEGG coverage survey: descriptive hit counts and coverage percentages over
//! a KEGG entry catalogue, with **no statistical test**.
//!
//! This is the coverage route's counterpart to `crate::enrichment`, and the
//! contrast is the point. `run_ora` asks "is this entry enriched among the
//! metabolites that changed", which needs a foreground `K` drawn from a
//! background `N` — and without a differential comparison there is no
//! defensible `K`. [`compute`] asks only "how much of this entry did I detect",
//! which is a count and a ratio.
//!
//! Owner: the `kegg-coverage` capability.

pub mod detect;
pub mod export;
pub mod types;

use std::collections::HashSet;

use crate::kegg::types::KeggCompoundSet;

pub use types::{CoverageResult, CoverageRow, CoverageSortKey};

/// Compute the coverage of every catalogue entry by the detected compound set.
///
/// `detected` is `D` — the detected, KEGG-mapped cpd IDs. `entries` is the
/// assembled Pathway or Module catalogue.
///
/// **Pure**: no I/O, no global state, no `tracing` events, no panics on any
/// well-typed input. Identical inputs return an identical `CoverageResult`,
/// row order included.
///
/// Returns one row for **every** entry, in the input order — including
/// `hits == 0` entries and `entry_size == 0` entries. It applies no filtering
/// of its own: `coverage_min_entry_size`, `min_hit_count`, and `top_n` are
/// result-screen filters re-applied live over these rows by [`displayed_rows`],
/// and there is no entry-identity-based exclusion at any layer. Fixing the
/// returned order here is what makes the determinism guarantee testable rather
/// than vacuous.
///
/// Takes NO name map: the `cpd -> metabolite name` mapping is a presentation
/// concern that travels on the CSV exporter's context, so `hit_compounds` holds
/// the bare cpd IDs the on-screen table renders.
///
/// KEGG's global/overview maps (`hsa01100` and the 12 other `br08901`
/// "Global and overview maps" entries) carry no `COMPOUND` section, so they
/// arrive as `compounds: []` and yield `entry_size == 0` rows. They are removed
/// by the minimum-entry-size floor like any other empty entry and are NOT
/// special-cased by ID here or anywhere else.
pub fn compute(detected: &HashSet<String>, entries: &[KeggCompoundSet]) -> CoverageResult {
    let detected_total = detected.len();
    // `|D|` is fixed across the whole result, so hoist the guarded reciprocal
    // rather than branching per row.
    let share_denom = detected_total as f64;

    let mut entries_without_compounds = 0usize;
    // Every detected compound reachable from ANY entry — the provenance count.
    // Accumulated over all entries before any display filter, so the result
    // screen's filters cannot move it.
    let mut reached: HashSet<&str> = HashSet::new();

    let mut rows = Vec::with_capacity(entries.len());
    for entry in entries {
        // Set semantics: `entry_size` is `|C|`, so a compound listed twice in a
        // KEGG record counts once — and counts once in `hits` too. Using the
        // raw slice length for one and set intersection for the other would let
        // `coverage` exceed 1.0.
        let compound_set: HashSet<&str> = entry.compounds.iter().map(String::as_str).collect();
        let entry_size = compound_set.len();
        if entry_size == 0 {
            entries_without_compounds += 1;
        }

        // Collected as borrows of `entry.compounds` (which outlives the loop) so
        // the same slice can feed both `reached` and the owned row field.
        let mut hit_refs: Vec<&str> = compound_set
            .iter()
            .copied()
            .filter(|c| detected.contains(*c))
            .collect();
        // HashSet iteration order is not stable across runs; sorting is what
        // makes the determinism guarantee hold for `hit_compounds` too.
        hit_refs.sort_unstable();
        let hits = hit_refs.len();
        reached.extend(hit_refs.iter().copied());
        let hit_compounds: Vec<String> = hit_refs.iter().map(|c| (*c).to_string()).collect();

        // Both denominators are guarded: an empty entry and an empty `D` are
        // ordinary inputs here, not error cases, and must yield 0.0 rather than
        // NaN — a NaN would propagate into the sort, the plot, and the CSV.
        let coverage = if entry_size == 0 {
            0.0
        } else {
            hits as f64 / entry_size as f64
        };
        let share = if detected_total == 0 {
            0.0
        } else {
            hits as f64 / share_denom
        };

        rows.push(CoverageRow {
            entry_id: entry.id.clone(),
            entry_name: entry.name.clone(),
            entry_size,
            hits,
            coverage,
            share,
            hit_compounds,
        });
    }

    CoverageResult {
        entries_total: rows.len(),
        rows,
        detected_total,
        entries_without_compounds,
        detected_in_entries: reached.len(),
    }
}

/// The complete display-filter chain, in one struct so the table, the dot plot,
/// and the CSV cannot drift apart.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DisplayFilters {
    /// Clamped to `>= 1` by the caller (see `clamp_coverage_min_entry_size`),
    /// so an `entry_size == 0` row can never be displayed, plotted, or exported.
    pub min_entry_size: usize,
    pub min_hit_count: usize,
    pub sort_key: CoverageSortKey,
    pub top_n: usize,
}

/// Apply the entire display pipeline: drop `entry_size < min_entry_size`, drop
/// `hits < min_hit_count`, sort by `sort_key`, truncate to `top_n`.
///
/// The results table, the coverage dot plot, and the CSV exporter MUST all
/// obtain their rows from here, and none may re-implement any part of the
/// chain. That is what makes "the rows drawn always equal the rows in the
/// table" true by construction rather than by three implementations agreeing —
/// the enrichment route specifies the equivalent chain in two places and relies
/// on them staying in step.
///
/// All keys sort descending except `EntryId` (ascending, lexicographic). Ties
/// break by descending `hits`, then ascending `entry_id`. `f64` comparisons go
/// through `total_cmp`, so the ordering is total and cannot panic on a NaN that
/// some future arithmetic change lets through.
pub fn displayed_rows(result: &CoverageResult, filters: DisplayFilters) -> Vec<&CoverageRow> {
    let mut rows: Vec<&CoverageRow> = result
        .rows
        .iter()
        .filter(|r| r.entry_size >= filters.min_entry_size && r.hits >= filters.min_hit_count)
        .collect();

    rows.sort_by(|a, b| {
        let primary = match filters.sort_key {
            CoverageSortKey::Coverage => b.coverage.total_cmp(&a.coverage),
            CoverageSortKey::Hits => b.hits.cmp(&a.hits),
            CoverageSortKey::EntrySize => b.entry_size.cmp(&a.entry_size),
            CoverageSortKey::EntryId => a.entry_id.cmp(&b.entry_id),
        };
        primary
            .then_with(|| b.hits.cmp(&a.hits))
            .then_with(|| a.entry_id.cmp(&b.entry_id))
    });

    rows.truncate(filters.top_n);
    rows
}

#[cfg(test)]
mod tests {
    use super::*;

    fn entry(id: &str, name: &str, compounds: &[&str]) -> KeggCompoundSet {
        KeggCompoundSet {
            id: id.to_string(),
            name: name.to_string(),
            compounds: compounds.iter().map(|s| s.to_string()).collect(),
        }
    }

    fn detected(ids: &[&str]) -> HashSet<String> {
        ids.iter().map(|s| s.to_string()).collect()
    }

    fn row<'a>(result: &'a CoverageResult, id: &str) -> &'a CoverageRow {
        result
            .rows
            .iter()
            .find(|r| r.entry_id == id)
            .unwrap_or_else(|| panic!("row {id} should be present"))
    }

    /// Partial coverage arithmetic: 18 of 42 compounds detected.
    #[test]
    fn partial_coverage_arithmetic() {
        let compounds: Vec<String> = (0..42).map(|i| format!("C{i:05}")).collect();
        let refs: Vec<&str> = compounds.iter().map(String::as_str).collect();
        let d = detected(&refs[..18]);
        let result = compute(&d, &[entry("p1", "Partial", &refs)]);

        let r = row(&result, "p1");
        assert_eq!(r.entry_size, 42);
        assert_eq!(r.hits, 18);
        assert!((r.coverage - 18.0 / 42.0).abs() < 1e-12);
        assert!((r.share - 1.0).abs() < 1e-12, "all 18 of D hit this entry");
        assert_eq!(r.hit_compounds.len(), 18);
    }

    /// An entry disjoint from `D` still gets a row — it is not silently dropped.
    #[test]
    fn zero_hit_entries_are_retained() {
        let d = detected(&["C00001", "C00002"]);
        let result = compute(&d, &[entry("p1", "Disjoint", &["C00031", "C00033"])]);

        assert_eq!(result.rows.len(), 1);
        let r = row(&result, "p1");
        assert_eq!(r.hits, 0);
        assert_eq!(r.coverage, 0.0);
        assert!(r.hit_compounds.is_empty());
    }

    /// A zero-compound entry (every KEGG global/overview map) divides by zero
    /// twice over and must produce 0.0 both times, plus be counted.
    #[test]
    fn zero_compound_entries_produce_a_row_without_nan() {
        let d = detected(&["C00001"]);
        let result = compute(
            &d,
            &[
                entry("hsa01100", "Metabolic pathways", &[]),
                entry("hsa00010", "Glycolysis", &["C00001", "C00002"]),
                entry("hsa01110", "Biosynthesis of secondary metabolites", &[]),
            ],
        );

        assert_eq!(result.entries_total, 3);
        assert_eq!(result.entries_without_compounds, 2);
        for id in ["hsa01100", "hsa01110"] {
            let r = row(&result, id);
            assert_eq!(r.entry_size, 0);
            assert_eq!(r.hits, 0);
            assert_eq!(r.coverage, 0.0);
            assert_eq!(r.share, 0.0);
            assert!(r.coverage.is_finite() && r.share.is_finite());
        }
    }

    /// An entry whose ID looks like a global map but which HAS compounds is
    /// treated exactly like any other entry — nothing keys on the identifier.
    #[test]
    fn no_entry_is_excluded_by_identity() {
        let d = detected(&["C00031", "C00092"]);
        let result = compute(&d, &[entry("hsa01100", "Metabolic pathways", &["C00031"])]);

        let r = row(&result, "hsa01100");
        assert_eq!(r.entry_size, 1);
        assert_eq!(r.hits, 1);
        assert_eq!(r.coverage, 1.0);
        assert_eq!(result.entries_without_compounds, 0);

        let shown = displayed_rows(
            &result,
            DisplayFilters {
                min_entry_size: 1,
                min_hit_count: 0,
                sort_key: CoverageSortKey::Coverage,
                top_n: 100,
            },
        );
        assert_eq!(shown.len(), 1);
        assert_eq!(shown[0].entry_id, "hsa01100");
    }

    /// An empty `D` guards the `share` denominator.
    #[test]
    fn empty_detected_set_does_not_panic() {
        let result = compute(
            &HashSet::new(),
            &[entry("p1", "A", &["C00001"]), entry("p2", "B", &[])],
        );

        assert_eq!(result.detected_total, 0);
        assert_eq!(result.detected_in_entries, 0);
        for r in &result.rows {
            assert_eq!(r.hits, 0);
            assert_eq!(r.coverage, 0.0);
            assert_eq!(r.share, 0.0);
        }
    }

    /// An empty entry slice is a valid input, not an error.
    #[test]
    fn empty_entry_slice_yields_an_empty_result() {
        let result = compute(&detected(&["C00001"]), &[]);
        assert!(result.rows.is_empty());
        assert_eq!(result.entries_total, 0);
        assert_eq!(result.entries_without_compounds, 0);
        assert_eq!(result.detected_total, 1);
        assert_eq!(result.detected_in_entries, 0);
    }

    /// Identical inputs give an identical result, row order and each row's
    /// `hit_compounds` order included. The `HashSet` intersection makes this
    /// non-trivial: without the explicit sort, `hit_compounds` would vary.
    #[test]
    fn determinism_on_identical_inputs() {
        let d = detected(&["C00031", "C00092", "C00001", "C00267"]);
        let entries = [
            entry(
                "p1",
                "A",
                &["C00031", "C00092", "C00001", "C00267", "C99999"],
            ),
            entry("p2", "B", &["C00267", "C00031"]),
            entry("p3", "C", &[]),
        ];
        assert_eq!(compute(&d, &entries), compute(&d, &entries));
        assert_eq!(
            compute(&d, &entries).rows[0].hit_compounds,
            vec!["C00001", "C00031", "C00092", "C00267"]
        );
    }

    /// Rows come back in the INPUT order, so sorting stays a display concern.
    #[test]
    fn rows_preserve_input_order() {
        let d = detected(&["C00001"]);
        let result = compute(
            &d,
            &[
                entry("zzz", "Last alphabetically, full coverage", &["C00001"]),
                entry("aaa", "First alphabetically, no hits", &["C99999"]),
            ],
        );
        assert_eq!(
            result.rows.iter().map(|r| &r.entry_id).collect::<Vec<_>>(),
            vec!["zzz", "aaa"]
        );
    }

    /// `detected_in_entries` counts `D` members reaching ANY entry — including
    /// entries a display filter would later remove.
    #[test]
    fn detected_in_entries_counts_compounds_reaching_any_entry() {
        let d = detected(&["C00001", "C00002", "C00003"]);
        let result = compute(
            &d,
            &[
                // C00003 is reachable ONLY through this 1-compound entry, which
                // the default min_entry_size of 3 would filter out of the table.
                entry("small", "Tiny", &["C00003"]),
                entry("big", "Large", &["C00001", "C00002", "C99998", "C99999"]),
            ],
        );

        assert_eq!(result.detected_total, 3);
        assert_eq!(result.detected_in_entries, 3);

        // Raising the floor changes which rows show, never the provenance count.
        let shown = displayed_rows(
            &result,
            DisplayFilters {
                min_entry_size: 3,
                min_hit_count: 0,
                sort_key: CoverageSortKey::Coverage,
                top_n: 100,
            },
        );
        assert_eq!(shown.len(), 1);
        assert_eq!(result.detected_in_entries, 3);
    }

    /// A compound listed twice in one KEGG record counts once in both `|C|` and
    /// `hits`, so `coverage` can never exceed 1.0.
    #[test]
    fn duplicate_compounds_in_an_entry_count_once() {
        let d = detected(&["C00031"]);
        let result = compute(&d, &[entry("p1", "Dup", &["C00031", "C00031", "C00092"])]);

        let r = row(&result, "p1");
        assert_eq!(r.entry_size, 2);
        assert_eq!(r.hits, 1);
        assert!(r.coverage <= 1.0);
        assert_eq!(r.hit_compounds, vec!["C00031"]);
    }

    // ---- displayed_rows ----

    fn mk_row(id: &str, entry_size: usize, hits: usize) -> CoverageRow {
        CoverageRow {
            entry_id: id.to_string(),
            entry_name: format!("Entry {id}"),
            entry_size,
            hits,
            coverage: if entry_size == 0 {
                0.0
            } else {
                hits as f64 / entry_size as f64
            },
            share: 0.0,
            hit_compounds: vec![],
        }
    }

    fn result_of(rows: Vec<CoverageRow>) -> CoverageResult {
        CoverageResult {
            entries_total: rows.len(),
            entries_without_compounds: rows.iter().filter(|r| r.entry_size == 0).count(),
            rows,
            detected_total: 0,
            detected_in_entries: 0,
        }
    }

    /// The default filter values applied to a mixed catalogue: the floor of 3
    /// removes the 0–2 band and coverage sorts descending.
    #[test]
    fn default_floor_removes_the_zero_to_two_band() {
        let result = result_of(vec![
            mk_row("e0", 0, 0),
            mk_row("e1", 1, 1),
            mk_row("e2", 2, 2),
            mk_row("e3", 3, 3),
            mk_row("e4", 4, 2),
        ]);
        let shown = displayed_rows(
            &result,
            DisplayFilters {
                min_entry_size: 3,
                min_hit_count: 0,
                sort_key: CoverageSortKey::Coverage,
                top_n: 100,
            },
        );
        assert_eq!(
            shown.iter().map(|r| &r.entry_id).collect::<Vec<_>>(),
            vec!["e3", "e4"]
        );
    }

    /// Equal coverage breaks by descending hits, then ascending entry_id.
    #[test]
    fn equal_coverage_breaks_by_hits_then_entry_id() {
        let result = result_of(vec![
            mk_row("b", 6, 3),   // coverage 0.5, hits 3
            mk_row("a", 42, 21), // coverage 0.5, hits 21
            mk_row("z", 6, 3),   // coverage 0.5, hits 3
            mk_row("c", 6, 3),   // coverage 0.5, hits 3
        ]);
        let shown = displayed_rows(
            &result,
            DisplayFilters {
                min_entry_size: 1,
                min_hit_count: 0,
                sort_key: CoverageSortKey::Coverage,
                top_n: 100,
            },
        );
        assert_eq!(
            shown.iter().map(|r| &r.entry_id).collect::<Vec<_>>(),
            vec!["a", "b", "c", "z"]
        );
    }

    /// Every sort key, including `EntryId`'s ascending exception.
    #[test]
    fn each_sort_key_orders_as_specified() {
        let result = result_of(vec![
            mk_row("m", 10, 1), // coverage 0.10
            mk_row("a", 4, 2),  // coverage 0.50
            mk_row("z", 20, 3), // coverage 0.15
        ]);
        let shown = |key| {
            displayed_rows(
                &result,
                DisplayFilters {
                    min_entry_size: 1,
                    min_hit_count: 0,
                    sort_key: key,
                    top_n: 100,
                },
            )
            .iter()
            .map(|r| r.entry_id.clone())
            .collect::<Vec<_>>()
        };
        assert_eq!(shown(CoverageSortKey::Coverage), vec!["a", "z", "m"]);
        assert_eq!(shown(CoverageSortKey::Hits), vec!["z", "a", "m"]);
        assert_eq!(shown(CoverageSortKey::EntrySize), vec!["z", "m", "a"]);
        assert_eq!(shown(CoverageSortKey::EntryId), vec!["a", "m", "z"]);
    }

    /// `min_hit_count` and `top_n` apply in the specified order: filter, sort,
    /// then truncate — so `top_n` takes the best rows, not the first ones.
    #[test]
    fn min_hit_count_then_sort_then_truncate() {
        let result = result_of(vec![
            mk_row("low", 10, 1), // coverage 0.10, hits 1
            mk_row("best", 4, 4), // coverage 1.00, hits 4
            mk_row("mid", 10, 5), // coverage 0.50, hits 5
        ]);
        let shown = displayed_rows(
            &result,
            DisplayFilters {
                min_entry_size: 1,
                min_hit_count: 2,
                sort_key: CoverageSortKey::Coverage,
                top_n: 1,
            },
        );
        assert_eq!(shown.len(), 1);
        assert_eq!(shown[0].entry_id, "best");
    }

    /// Sorting the same result twice gives the same order.
    #[test]
    fn sorting_is_stable_across_repeated_renders() {
        let result = result_of(vec![
            mk_row("a", 6, 3),
            mk_row("b", 6, 3),
            mk_row("c", 4, 2),
        ]);
        let filters = DisplayFilters {
            min_entry_size: 1,
            min_hit_count: 0,
            sort_key: CoverageSortKey::Coverage,
            top_n: 100,
        };
        assert_eq!(
            displayed_rows(&result, filters)
                .iter()
                .map(|r| &r.entry_id)
                .collect::<Vec<_>>(),
            displayed_rows(&result, filters)
                .iter()
                .map(|r| &r.entry_id)
                .collect::<Vec<_>>()
        );
    }

    /// A `min_entry_size` of 1 keeps every zero-compound entry out, which is
    /// what the hard minimum exists to guarantee.
    #[test]
    fn the_floor_at_its_minimum_still_excludes_empty_entries() {
        let result = result_of(vec![mk_row("empty", 0, 0), mk_row("real", 1, 0)]);
        let shown = displayed_rows(
            &result,
            DisplayFilters {
                min_entry_size: 1,
                min_hit_count: 0,
                sort_key: CoverageSortKey::Coverage,
                top_n: 100,
            },
        );
        assert_eq!(shown.len(), 1);
        assert_eq!(shown[0].entry_id, "real");
    }
}
