//! Over-representation analysis (ORA): per-entry hypergeometric test +
//! user-selected FDR correction (BY default for Stage 3 ORA; BH or no correction
//! opt-in) over all entries.
//! An "entry" is a `KeggCompoundSet` — either a KEGG pathway (pathway
//! mode) or a KEGG module (module mode); the math is identical for both.

use statrs::distribution::{DiscreteCDF, Hypergeometric};
use std::collections::HashSet;
use tracing::error;

use crate::dam::fdr::{FdrMethod, adjust_pvalues};
use crate::enrichment::types::{EnrichmentDirection, EnrichmentResult, EnrichmentRow};
use crate::kegg::KeggCompoundSet;

/// Compute ORA for every entry. Rows are sorted by ascending FDR; ties
/// broken by ascending `entry_id`. The function is pure (no I/O).
///
/// Universe `N = |universe|`, foreground `K = |dam_cpd|`, per-entry:
///
/// - `M_p = |entry.compounds ∩ universe|` (entry size restricted to
///   what we could have measured).
/// - `k_p = |dam_cpd ∩ entry.compounds|`.
/// - `p_value`: `1.0` if any of `k_p == 0 || M_p == 0 || K == 0 || N == 0`
///   (short-circuit to avoid undefined CDF arguments). Otherwise
///   `1 - HypergeometricCDF(k_p - 1; N, M_p, K)`.
///
/// `fdr_method` selects the FDR correction (BY default for Stage 3 ORA; BH or no correction opt-in). The
/// chosen method is corrected over ALL entries' p-values (including those
/// that short-circuited to 1.0) and recorded on `EnrichmentResult.fdr_method`.
/// BY is the more conservative choice when entries share compounds, which is
/// the typical situation for pathway/module ORA.
pub fn run_ora(
    universe: &HashSet<String>,
    dam_cpd: &HashSet<String>,
    entries: &[KeggCompoundSet],
    min_hit_count: usize,
    direction: EnrichmentDirection,
    fdr_method: FdrMethod,
    min_entry_size: usize,
) -> EnrichmentResult {
    let n_universe = universe.len();
    let k_size = dam_cpd.len();
    let empty_compound_count = entries.iter().filter(|e| e.compounds.is_empty()).count();

    // Empty universe: return an empty result. universe_size == 0 is the
    // sentinel callers check for the empty-state UI message. The pre-FDR
    // filter does not run (nothing to filter against), so the dropped
    // count is 0.
    if n_universe == 0 {
        return EnrichmentResult {
            universe_size: 0,
            dam_cpd_size: k_size,
            direction,
            min_hit_count,
            min_entry_size,
            entries_dropped_by_min_entry_size: 0,
            empty_compound_count,
            rows: vec![],
            fdr_method,
        };
    }

    // K ⊆ N invariant. K is built upstream from the same `feature_to_cpds`
    // map as N, so violations should be impossible by construction — but a
    // future refactor could break the invariant silently, producing
    // `Hypergeometric::new` errors that the per-entry loop maps to p=1.0
    // (all-non-significant) with no signal upstream. `debug_assert!` makes
    // any such regression fail loudly in dev / test; the per-run summary
    // `error!` log below provides the same signal in release builds. Placed
    // AFTER the empty-N short-circuit because callers may legitimately pass
    // a non-empty K with N=∅ as a degenerate "nothing to enrich against"
    // probe — that path returns empty rows without running Hypergeometric.
    debug_assert!(
        dam_cpd.is_subset(universe),
        "ORA invariant: K must be ⊆ N; got |K|={}, |N|={}",
        dam_cpd.len(),
        universe.len()
    );

    // Per-entry stats (before FDR).
    struct PreFdrRow {
        entry_id: String,
        entry_name: String,
        hits: usize,
        total: usize,
        expected: f64,
        enrichment_ratio: f64,
        p_value: f64,
        hit_kegg_ids: Vec<String>,
    }

    // Per-entry stats accumulated via for-loop (instead of iterator map) so we
    // can tally Hypergeometric domain errors across entries. domain_errors > 0
    // signals an upstream K ⊄ N violation (caught in debug by the assert
    // above; observable in release via the error! summary below).
    let mut pre: Vec<PreFdrRow> = Vec::with_capacity(entries.len());
    let mut domain_errors: usize = 0;
    let mut entries_dropped_by_min_entry_size: usize = 0;
    for p in entries {
        // entry_set is HashSet-deduped; reuse it for m_p so dup-cpd KEGG
        // entries don't inflate m_p relative to k_p (which is naturally
        // set-cardinality). Pre-2026-05-26 m_p iterated raw &p.compounds
        // (Vec, no dedup) — every duplicate occurrence in the COMPOUND block
        // inflated m_p, biased expected upward and enrichment_ratio
        // downward, and could in pathological cases push m_p > n_universe
        // and trigger the Hypergeometric Err arm.
        let entry_set: HashSet<&String> = p.compounds.iter().collect();
        let m_p: usize = entry_set
            .iter()
            .filter(|c| universe.contains(c.as_str()))
            .count();

        // Pre-FDR `min_entry_size` filter. Entries with too few testable
        // compounds in the universe never enter the FDR family — they
        // can't produce useful signal and would only dilute `m`. The
        // dropped count is surfaced via `entries_dropped_by_min_entry_size`
        // on the result so the Stage 3 UI can show the retention chain.
        if m_p < min_entry_size {
            entries_dropped_by_min_entry_size += 1;
            continue;
        }

        // k_p: entry ∩ dam_cpd. Also collect the actual hit ids.
        let mut hit_ids: Vec<String> = Vec::new();
        for c in dam_cpd {
            if entry_set.contains(c) {
                hit_ids.push(c.clone());
            }
        }
        hit_ids.sort();
        let k_p = hit_ids.len();

        let expected = if n_universe > 0 {
            (k_size as f64) * (m_p as f64) / (n_universe as f64)
        } else {
            0.0
        };
        let enrichment_ratio = if expected > 0.0 {
            (k_p as f64) / expected
        } else {
            f64::NAN
        };

        let (p_value, domain_err) = entry_pvalue(n_universe, m_p, k_size, k_p);
        if domain_err {
            domain_errors += 1;
        }

        pre.push(PreFdrRow {
            entry_id: p.id.clone(),
            entry_name: p.name.clone(),
            hits: k_p,
            total: m_p,
            expected,
            enrichment_ratio,
            p_value,
            hit_kegg_ids: hit_ids,
        });
    }

    if domain_errors > 0 {
        // This should never fire in normal operation (the debug_assert above
        // would have already caught K ⊄ N in dev/test). Surface a per-run
        // ERROR log if it does — signals an upstream invariant break that
        // would otherwise present as "all entries non-significant" with no
        // diagnostic.
        error!(
            domain_errors,
            total_entries = entries.len(),
            n_universe,
            k_size,
            "ORA hypergeometric domain errors: this indicates K ⊄ N or m_p > N (invariant violation upstream); affected entries assigned p=1.0"
        );
    }

    // FDR over all p-values (including k=0 rows with p=1.0). User-selected
    // BH (literature default) or BY (more conservative for shared-compound entries).
    let p_vec: Vec<f64> = pre.iter().map(|r| r.p_value).collect();
    let fdr_vec = adjust_pvalues(&p_vec, fdr_method);

    // Materialise rows with FDR, then sort.
    let mut rows: Vec<EnrichmentRow> = pre
        .drain(..)
        .zip(fdr_vec)
        .map(|(r, fdr)| EnrichmentRow {
            entry_id: r.entry_id,
            entry_name: r.entry_name,
            hits: r.hits,
            total: r.total,
            expected: r.expected,
            enrichment_ratio: r.enrichment_ratio,
            p_value: r.p_value,
            fdr,
            hit_kegg_ids: r.hit_kegg_ids,
        })
        .collect();

    rows.sort_by(|a, b| {
        a.fdr
            .partial_cmp(&b.fdr)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.entry_id.cmp(&b.entry_id))
    });

    EnrichmentResult {
        universe_size: n_universe,
        dam_cpd_size: k_size,
        direction,
        min_hit_count,
        min_entry_size,
        entries_dropped_by_min_entry_size,
        empty_compound_count,
        rows,
        fdr_method,
    }
}

/// Pure per-entry p-value decision, factored out of `run_ora`'s loop so the
/// domain-error path is unit-testable in debug without tripping `run_ora`'s
/// K ⊆ N `debug_assert!`. Owns the WHOLE decision and returns `(p_value, domain_err)`:
///
/// - zero-input short-circuit (`k_p == 0 || m_p == 0 || k_size == 0 || n_universe == 0`)
///   → `(1.0, false)` — a legitimate zero-hit (`k_p == 0`) row is NOT a domain error;
/// - `Hypergeometric::new` Err arm → `(1.0, true)` — reachable only when `k_size > n_universe`
///   (i.e. K ⊄ N), since `m_p ≤ n_universe` always (it is an intersection with the universe);
/// - otherwise → `(1 - CDF(k_p - 1; n_universe, m_p, k_size), false)`.
///
/// Carries no `debug_assert!`, so feeding `m_p > n_universe` / `k_size > n_universe`
/// exercises the `domain_err = true` arm directly in a debug `cargo test`.
fn entry_pvalue(n_universe: usize, m_p: usize, k_size: usize, k_p: usize) -> (f64, bool) {
    if k_p == 0 || m_p == 0 || k_size == 0 || n_universe == 0 {
        return (1.0, false);
    }
    // Hypergeometric in statrs: parameters are
    // (population N, successes-in-population K_pop, draws n).
    // We test "probability of seeing AT LEAST k_p" = 1 - P(X <= k_p - 1).
    match Hypergeometric::new(n_universe as u64, m_p as u64, k_size as u64) {
        Ok(dist) => (1.0 - dist.cdf((k_p as u64).saturating_sub(1)), false),
        Err(_) => (1.0, true),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::collections::HashSet;

    fn entry(id: &str, name: &str, compounds: &[&str]) -> KeggCompoundSet {
        KeggCompoundSet {
            id: id.to_string(),
            name: name.to_string(),
            compounds: compounds.iter().map(|s| s.to_string()).collect(),
        }
    }

    fn make_set<I: IntoIterator<Item = &'static str>>(items: I) -> HashSet<String> {
        items.into_iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn entry_pvalue_owns_short_circuit_and_flags_only_out_of_range() {
        // Out-of-range inputs are the ONLY ones that flag a domain error — and
        // they are reachable in a real run only when K ⊄ N. Tested here directly
        // because run_ora's debug_assert! would panic on such inputs first.
        let (p, err) = entry_pvalue(5, 6, 3, 2); // m_p > n_universe
        assert_eq!(p, 1.0);
        assert!(err, "m_p > n_universe must flag domain_err");
        let (p, err) = entry_pvalue(5, 3, 6, 2); // k_size > n_universe
        assert_eq!(p, 1.0);
        assert!(err, "k_size > n_universe must flag domain_err");

        // Every zero-input short-circuit yields (1.0, false): a legitimate
        // zero-hit (k_p == 0) row is NOT a domain error.
        for (n, m, k, kp) in [(5, 3, 2, 0), (5, 0, 2, 1), (5, 3, 0, 1), (0, 3, 2, 1)] {
            let (p, err) = entry_pvalue(n, m, k, kp);
            assert_eq!(
                p, 1.0,
                "short-circuit must yield p=1.0 for ({n},{m},{k},{kp})"
            );
            assert!(
                !err,
                "short-circuit must NOT flag domain_err for ({n},{m},{k},{kp})"
            );
        }

        // Valid in-range input: finite upper-tail p in (0, 1], no domain error.
        let (p, err) = entry_pvalue(100, 10, 20, 5);
        assert!(!err, "valid in-range input must not flag domain_err");
        assert!(p > 0.0 && p <= 1.0, "valid p must be in (0, 1], got {p}");
    }

    #[test]
    fn k_zero_short_circuits_to_p_one() {
        // K is empty: every entry short-circuits to p=1.0.
        let universe = make_set(["A", "B", "C", "D", "E"]);
        let dam_cpd: HashSet<String> = HashSet::new();
        let entries = vec![entry("p1", "Path 1", &["A", "B"])];
        let result = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiYekutieli,
            1,
        );
        assert_eq!(result.rows.len(), 1);
        assert_eq!(result.rows[0].p_value, 1.0);
        assert_eq!(result.rows[0].hits, 0);
        assert!(result.rows[0].expected == 0.0);
        assert!(result.rows[0].enrichment_ratio.is_nan());
    }

    #[test]
    fn m_zero_drops_entry_under_default_min_entry_size() {
        // Entry's compounds none in universe: M_p = 0. Under min_entry_size = 1
        // (effectively-no-filter default), m_p=0 < 1 still drops the entry.
        let universe = make_set(["A", "B"]);
        let dam_cpd = make_set(["A"]);
        let entries = vec![entry("p1", "Path 1", &["X", "Y"])];
        let result = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiYekutieli,
            1,
        );
        assert!(result.rows.is_empty(), "m_p=0 entry must be filtered out");
        assert_eq!(result.entries_dropped_by_min_entry_size, 1);
    }

    #[test]
    fn empty_universe_returns_no_rows() {
        let universe: HashSet<String> = HashSet::new();
        let dam_cpd = make_set(["A"]);
        let entries = vec![entry("p1", "Path 1", &["A"])];
        let result = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiYekutieli,
            1,
        );
        assert_eq!(result.universe_size, 0);
        assert!(result.rows.is_empty());
    }

    #[test]
    fn basic_enrichment_produces_finite_p() {
        // Universe N=10, K=3, entry has M=4 (A,B,E,F) of which 2 are in K (A,B).
        let universe = make_set(["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]);
        let dam_cpd = make_set(["A", "B", "C"]);
        let entries = vec![entry("p1", "Path 1", &["A", "B", "E", "F"])];
        let result = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiYekutieli,
            1,
        );
        let row = &result.rows[0];
        assert_eq!(row.hits, 2);
        assert_eq!(row.total, 4);
        assert!((row.expected - 1.2).abs() < 1e-9);
        assert!((row.enrichment_ratio - (2.0 / 1.2)).abs() < 1e-9);
        assert!(row.p_value.is_finite());
        assert!(row.p_value < 1.0);
        // FDR with only one entry equals p_value × c(1) = p_value × 1.
        assert!((row.fdr - row.p_value).abs() < 1e-9);
        // hit_kegg_ids sorted.
        assert_eq!(row.hit_kegg_ids, vec!["A".to_string(), "B".to_string()]);
    }

    #[test]
    fn rows_sorted_by_fdr_then_entry_id() {
        let universe = make_set(["A", "B", "C", "D"]);
        let dam_cpd = make_set(["A", "B"]);
        let entries = vec![
            entry("z1", "Z", &["C", "D"]),  // no hits
            entry("a1", "A", &["A", "B"]),  // 2 hits
            entry("a2", "A2", &["A", "B"]), // 2 hits — same p as a1
        ];
        let result = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiYekutieli,
            1,
        );
        // First two should be a1, a2 (tied by FDR, broken by id).
        // z1 (no hits, p=1) should be last.
        assert_eq!(result.rows[0].entry_id, "a1");
        assert_eq!(result.rows[1].entry_id, "a2");
        assert_eq!(result.rows[2].entry_id, "z1");
    }

    /// `min_hit_count` is inert inside `run_ora`: passing different values
    /// returns identical rows, p-values and adjusted values.
    ///
    /// It used to set a per-row `displayed` flag here, which froze a display
    /// filter to the moment of the run. The filter now lives with the
    /// consumers that draw and count; `run_ora` only records the value on the
    /// result as provenance. This test pins the inertness by CALLING TWICE —
    /// an earlier version asserted properties of the fixture and would have
    /// passed for any argument.
    #[test]
    fn min_hit_count_does_not_affect_run_ora_output() {
        let universe = make_set(["A", "B", "C", "D"]);
        let dam_cpd = make_set(["A", "B"]);
        let entries = vec![
            entry("p1", "Path 1", &["A", "B", "C"]), // 2 hits
            entry("p2", "Path 2", &["D"]),           // 0 hits
        ];
        let loose = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiYekutieli,
            1,
        );
        let strict = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            99,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiYekutieli,
            1,
        );
        assert_eq!(loose.rows.len(), strict.rows.len());
        for (a, b) in loose.rows.iter().zip(&strict.rows) {
            assert_eq!(a.entry_id, b.entry_id);
            assert_eq!(a.hits, b.hits);
            assert_eq!(a.p_value, b.p_value);
            assert_eq!(a.fdr, b.fdr);
        }
        // Only the provenance record differs.
        assert_eq!(loose.min_hit_count, 1);
        assert_eq!(strict.min_hit_count, 99);
    }

    #[test]
    fn fdr_m_equals_all_entries_including_zero_hit() {
        // 5 entries, only 1 has any hit. m for BY FDR MUST be 5, not 1.
        let universe = make_set(["A", "B", "C", "D", "E"]);
        let dam_cpd = make_set(["A"]);
        let entries = vec![
            entry("p1", "P1", &["A"]),
            entry("p2", "P2", &["B"]),
            entry("p3", "P3", &["C"]),
            entry("p4", "P4", &["D"]),
            entry("p5", "P5", &["E"]),
        ];
        let result = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiYekutieli,
            1,
        );
        // All 5 rows present.
        assert_eq!(result.rows.len(), 5);
        // p1's raw p was K=1, M=1, N=5, k=1 → 1 - hypercdf(0; 5, 1, 1) = 1 - 4/5 = 0.2.
        // After BY (c(5) = 1 + 1/2 + 1/3 + 1/4 + 1/5 ≈ 2.283), p1's
        // adjusted = 0.2 * c(5) * 5 / 1 = 2.283. Clamped to 1.0.
        let p1 = result.rows.iter().find(|r| r.entry_id == "p1").unwrap();
        assert!((p1.p_value - 0.2).abs() < 1e-9);
        assert!(p1.fdr >= 0.9, "BY-adjusted p must be much larger than raw");
    }

    #[test]
    fn empty_compound_count_surfaces_in_result() {
        // Module-mode case: some entries have no compound block. Under
        // min_entry_size = 1 those empty-compound entries (m_p = 0) are
        // dropped from `rows`, but `empty_compound_count` still counts them
        // as a diagnostic input statistic.
        let universe = make_set(["A", "B", "C"]);
        let dam_cpd = make_set(["A"]);
        let entries = vec![
            entry("p1", "Has compounds", &["A", "B"]),
            entry("p2", "No compounds 1", &[]),
            entry("p3", "No compounds 2", &[]),
        ];
        let result = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiYekutieli,
            1,
        );
        assert_eq!(result.empty_compound_count, 2);
        // p2 / p3 have m_p = 0 < min_entry_size = 1 → dropped from rows.
        assert!(result.rows.iter().all(|r| r.entry_id != "p2"));
        assert!(result.rows.iter().all(|r| r.entry_id != "p3"));
        // Only p1 (m_p = 2) survives.
        assert_eq!(result.rows.len(), 1);
        assert_eq!(result.rows[0].entry_id, "p1");
        // Both empty-compound entries appear in the dropped count along with
        // anyone else below the threshold.
        assert!(result.entries_dropped_by_min_entry_size >= 2);
    }

    /// An entry with `M_p = 2` is dropped when `min_entry_size = 3`. The
    /// dropped entry's would-be p-value does NOT enter the FDR family — `m`
    /// equals the surviving count.
    #[test]
    fn min_entry_size_drops_below_threshold() {
        let universe = make_set(["A", "B", "C", "D", "E"]);
        let dam_cpd = make_set(["A", "B"]);
        let entries = vec![
            entry("p1", "Path 1", &["A", "B", "C"]), // m_p = 3 → kept
            entry("p2", "Path 2", &["A", "D"]),      // m_p = 2 → dropped
        ];
        let result = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiHochberg,
            3,
        );
        assert_eq!(result.entries_dropped_by_min_entry_size, 1);
        assert_eq!(result.rows.len(), 1);
        assert_eq!(result.rows[0].entry_id, "p1");
        // p2's p-value MUST NOT enter the FDR family — the only surviving
        // entry's BH adjusted q is `p * m / rank` with m=1, so q equals raw p.
        assert!((result.rows[0].fdr - result.rows[0].p_value).abs() < 1e-9);
    }

    /// `min_entry_size = 1` is the "effectively no filter" choice. Only
    /// `M_p = 0` entries are dropped (mathematically un-testable). Entries
    /// with `M_p >= 1` are all retained, including zero-hit ones.
    #[test]
    fn min_entry_size_one_drops_only_zero_m_p_entries() {
        let universe = make_set(["A", "B", "C"]);
        let dam_cpd = make_set(["A"]);
        let entries = vec![
            entry("p1", "Has 1 in N", &["A"]),        // m_p = 1 → kept
            entry("p2", "Has 0 in N", &["X"]),        // m_p = 0 → dropped
            entry("p3", "Has 1 in N no hit", &["B"]), // m_p = 1, k_p = 0 → kept
        ];
        let result = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiHochberg,
            1,
        );
        assert_eq!(result.entries_dropped_by_min_entry_size, 1);
        assert_eq!(result.rows.len(), 2);
        let kept_ids: std::collections::HashSet<&str> =
            result.rows.iter().map(|r| r.entry_id.as_str()).collect();
        assert!(kept_ids.contains("p1"));
        assert!(kept_ids.contains("p3"));
        assert!(!kept_ids.contains("p2"));
    }

    /// Deliberately mixed fixture (M_p = 0/1/2/3/4) — with `min_entry_size = 3`
    /// exactly 3 entries are dropped and exactly 2 survive.
    #[test]
    fn entries_dropped_count_matches_expected() {
        let universe = make_set(["A", "B", "C", "D"]);
        let dam_cpd = make_set(["A", "B"]);
        let entries = vec![
            entry("p0", "m_p=0", &["X", "Y"]),
            entry("p1", "m_p=1", &["A", "X"]),
            entry("p2", "m_p=2", &["A", "B"]),
            entry("p3", "m_p=3", &["A", "B", "C"]),
            entry("p4", "m_p=4", &["A", "B", "C", "D"]),
        ];
        let result = run_ora(
            &universe,
            &dam_cpd,
            &entries,
            1,
            EnrichmentDirection::Both,
            FdrMethod::BenjaminiHochberg,
            3,
        );
        assert_eq!(result.entries_dropped_by_min_entry_size, 3);
        assert_eq!(result.rows.len(), 2);
    }
}
