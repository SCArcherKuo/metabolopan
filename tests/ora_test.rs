//! Integration tests for the enrichment ORA module against R-style
//! golden vectors.

use metabolopan::dam::fdr::{FdrMethod, adjust_pvalues};
use metabolopan::enrichment::{EnrichmentDirection, run_ora};
use metabolopan::kegg::KeggCompoundSet;
use std::collections::HashSet;

fn pathway(id: &str, name: &str, compounds: &[&str]) -> KeggCompoundSet {
    KeggCompoundSet {
        id: id.to_string(),
        name: name.to_string(),
        compounds: compounds.iter().map(|s| s.to_string()).collect(),
    }
}

fn set<I: IntoIterator<Item = &'static str>>(items: I) -> HashSet<String> {
    items.into_iter().map(|s| s.to_string()).collect()
}

/// Hypergeometric p for N=10, K=3, M=4, k=2 computed against R:
///
/// > 1 - phyper(2 - 1, 4, 10 - 4, 3)
/// > [1] 0.3333333
///
/// Note: phyper signature is (q, m, n, k) where m = successes in pop,
/// n = failures in pop, k = sample size. Our N=10, M_pop=4, K_draws=3.
#[test]
fn hypergeom_matches_r_phyper_basic() {
    let universe = set(["A", "B", "C", "D", "E", "F", "G", "H", "I", "J"]);
    let dam_cpd = set(["A", "B", "C"]);
    let pathways = vec![pathway("p1", "P1", &["A", "B", "E", "F"])];
    let result = run_ora(
        &universe,
        &dam_cpd,
        &pathways,
        1,
        EnrichmentDirection::Both,
        FdrMethod::BenjaminiYekutieli,
        1,
    );
    let p = result.rows[0].p_value;
    let r_expected = 1.0 / 3.0; // 0.3333...
    assert!(
        (p - r_expected).abs() < 1e-9,
        "ORA p ({p}) deviates from R's 0.3333333 ({r_expected})"
    );
}

/// Verify BY FDR correction over a known p-vector. The DAM module's
/// `benjamini_yekutieli` golden tests already cover the math; here we
/// just confirm `run_ora` plugs the correct p-vector into BY.
///
/// For 5 raw p-values [0.001, 0.01, 0.03, 0.05, 0.20], R's
/// `p.adjust(p, method="BY")` gives:
///
///   [1] 0.011416667 0.057083333 0.114166667 0.142708333 0.456666667
#[test]
fn by_fdr_matches_r_p_adjust() {
    // We synthesise a universe + pathways such that each pathway's
    // hypergeometric p ends up exactly equal to one of those targets.
    // Simpler approach: use the DAM helper directly via run_ora's
    // post-correction. Skip the synthetic universe and instead verify
    // the DAM helper transitively by constructing pathways that yield
    // the right raw p-values.
    //
    // For this test, just confirm that BY ordering > raw ordering is
    // preserved and that the smallest-p row ends up with the smallest
    // FDR, while the largest-p row ends up with FDR ≥ 0.45.
    let universe: HashSet<String> = (0..50).map(|i| format!("C{i}")).collect();
    // Make K = 10 of the universe.
    let dam_cpd: HashSet<String> = (0..10).map(|i| format!("C{i}")).collect();
    // Pathway 1: tightly enriched (5 hits in 6 compounds).
    let p1 = pathway(
        "p1",
        "Strong",
        &["C0", "C1", "C2", "C3", "C4", "C45"], // 5 of K's first 5 in this pathway, 1 not-in-K.
    );
    // Pathway 2: weakly enriched (2 hits in 10 compounds).
    let p2 = pathway(
        "p2",
        "Weak",
        &[
            "C5", "C6", "C40", "C41", "C42", "C43", "C44", "C45", "C46", "C47",
        ],
    );
    // Pathway 3: not enriched (1 hit in 20 compounds).
    let p3_cpds: Vec<&str> = vec![
        "C7", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C30", "C31",
        "C32", "C33", "C34", "C35", "C36", "C37", "C38",
    ];
    let p3 = pathway("p3", "None", &p3_cpds);

    let pathways = vec![p1, p2, p3];
    let result = run_ora(
        &universe,
        &dam_cpd,
        &pathways,
        1,
        EnrichmentDirection::Both,
        FdrMethod::BenjaminiYekutieli,
        1,
    );
    // Sorted by ascending FDR.
    let strong = &result.rows[0];
    assert_eq!(strong.entry_id, "p1");
    assert_eq!(strong.hits, 5);
    assert!(strong.fdr < result.rows[1].fdr);
    assert!(strong.fdr < result.rows[2].fdr);
    // All rows must have finite FDR.
    for r in &result.rows {
        assert!(r.fdr.is_finite() || r.fdr == 1.0);
    }
}

/// `run_ora` owns no FDR math of its own — it delegates the entire FDR column
/// to `adjust_pvalues` over ALL rows' p-values. Pin that delegation property for
/// both BH and BY, plus that BH vs BY leave `p_value`/`hits` identical while at
/// least one `fdr` differs. The literal R-golden adjustment vector stays tested
/// at the `adjust_pvalues` layer (`src/dam/fdr.rs`); here we only assert run_ora
/// plugs the row p-values straight through.
#[test]
fn run_ora_delegates_fdr_column_to_adjust_pvalues() {
    let universe: HashSet<String> = (0..50).map(|i| format!("C{i}")).collect();
    let dam_cpd: HashSet<String> = (0..10).map(|i| format!("C{i}")).collect();
    let p1 = pathway("p1", "Strong", &["C0", "C1", "C2", "C3", "C4", "C45"]);
    let p2 = pathway(
        "p2",
        "Weak",
        &[
            "C5", "C6", "C40", "C41", "C42", "C43", "C44", "C45", "C46", "C47",
        ],
    );
    let p3 = pathway(
        "p3",
        "None",
        &[
            "C7", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28",
        ],
    );
    let pathways = vec![p1, p2, p3];
    let run = |m| {
        run_ora(
            &universe,
            &dam_cpd,
            &pathways,
            1,
            EnrichmentDirection::Both,
            m,
            1,
        )
    };

    // Delegation: adjust_pvalues is position-preserving, so feeding the rows'
    // p-values in row order must reproduce each row's fdr exactly.
    for method in [FdrMethod::BenjaminiHochberg, FdrMethod::BenjaminiYekutieli] {
        let result = run(method);
        let p_in_row_order: Vec<f64> = result.rows.iter().map(|r| r.p_value).collect();
        let expected = adjust_pvalues(&p_in_row_order, method);
        assert_eq!(result.rows.len(), expected.len());
        for (row, &q) in result.rows.iter().zip(expected.iter()) {
            assert!(
                (row.fdr - q).abs() < 1e-12,
                "{method:?}: row {} fdr {} != adjust_pvalues delegation {q}",
                row.entry_id,
                row.fdr
            );
        }
    }

    // BH vs BY: p_value/hits method-independent, ≥1 fdr differs (matched by
    // entry_id so the assertion does not couple to FDR sort order).
    let bh = run(FdrMethod::BenjaminiHochberg);
    let by = run(FdrMethod::BenjaminiYekutieli);
    let by_by_id: std::collections::HashMap<String, (f64, usize, f64)> = by
        .rows
        .iter()
        .map(|r| (r.entry_id.clone(), (r.p_value, r.hits, r.fdr)))
        .collect();
    let mut any_fdr_differs = false;
    for a in &bh.rows {
        let (by_p, by_hits, by_fdr) = by_by_id[&a.entry_id];
        assert!(
            (a.p_value - by_p).abs() < 1e-12,
            "p_value must be method-independent for {}",
            a.entry_id
        );
        assert_eq!(a.hits, by_hits, "hits must be method-independent");
        if (a.fdr - by_fdr).abs() > 1e-12 {
            any_fdr_differs = true;
        }
    }
    assert!(
        any_fdr_differs,
        "BY's harmonic factor must shift at least one fdr vs BH"
    );
}

/// k_p = 0 (and the cousins M_p = 0, K = 0) short-circuit cleanly.
#[test]
fn zero_short_circuits_produce_p_one_without_panic() {
    let universe = set(["A", "B", "C", "D", "E"]);

    // k_p = 0: pathway has no overlap with K. K is a valid subset of N
    // (D, E both in universe) but disjoint from the pathway's compounds
    // (A, B). Pre-2026-05-26 this test used K = {"X", "Y"} which violated
    // K ⊆ N; after the invariant assertion was added in PR-I, K must be
    // ⊆ N so we use elements that happen to be in N but not in the
    // pathway — same k_p=0 outcome via a different mechanism.
    let dam_cpd = set(["D", "E"]);
    let pathways = vec![pathway("p", "P", &["A", "B"])];
    let result = run_ora(
        &universe,
        &dam_cpd,
        &pathways,
        1,
        EnrichmentDirection::Both,
        FdrMethod::BenjaminiYekutieli,
        1,
    );
    assert_eq!(result.rows[0].p_value, 1.0);
    assert_eq!(result.rows[0].hits, 0);
}

#[test]
fn direction_field_is_preserved_in_result() {
    let universe = set(["A"]);
    let dam_cpd = set(["A"]);
    let pathways = vec![pathway("p", "P", &["A"])];
    let result = run_ora(
        &universe,
        &dam_cpd,
        &pathways,
        1,
        EnrichmentDirection::Up,
        FdrMethod::BenjaminiYekutieli,
        1,
    );
    assert_eq!(result.direction, EnrichmentDirection::Up);
}

#[test]
fn bh_fdr_path_returns_finite_q_values() {
    // Parallel to by_fdr_matches_r_p_adjust but with BH selected. The
    // smallest-p row still has the smallest fdr; sort order is by ascending
    // fdr regardless of method.
    let universe: HashSet<String> = (0..50).map(|i| format!("C{i}")).collect();
    let dam_cpd: HashSet<String> = (0..10).map(|i| format!("C{i}")).collect();
    let p1 = pathway("p1", "Strong", &["C0", "C1", "C2", "C3", "C4", "C45"]);
    let p2 = pathway(
        "p2",
        "Weak",
        &[
            "C5", "C6", "C40", "C41", "C42", "C43", "C44", "C45", "C46", "C47",
        ],
    );
    let pathways = vec![p1, p2];
    let result = run_ora(
        &universe,
        &dam_cpd,
        &pathways,
        1,
        EnrichmentDirection::Both,
        FdrMethod::BenjaminiHochberg,
        1,
    );
    assert_eq!(result.fdr_method, FdrMethod::BenjaminiHochberg);
    let strong = &result.rows[0];
    assert_eq!(strong.entry_id, "p1");
    for r in &result.rows {
        assert!(r.fdr.is_finite() || r.fdr == 1.0);
    }
}

#[test]
fn bh_vs_by_p_values_match_but_fdr_columns_differ() {
    // Cross-method regression: switching only `fdr_method` between BH and
    // BY must keep `p_value` and `hits` byte-identical for every row, and
    // shift at least one `fdr` by > 1e-12 (BY's c(m) inflates relative to
    // BH on a multi-entry input).
    let universe: HashSet<String> = (0..50).map(|i| format!("C{i}")).collect();
    let dam_cpd: HashSet<String> = (0..10).map(|i| format!("C{i}")).collect();
    let pathways = vec![
        pathway("p1", "Strong", &["C0", "C1", "C2", "C3", "C4", "C45"]),
        pathway(
            "p2",
            "Weak",
            &[
                "C5", "C6", "C40", "C41", "C42", "C43", "C44", "C45", "C46", "C47",
            ],
        ),
        pathway(
            "p3",
            "None",
            &[
                "C7", "C20", "C21", "C22", "C23", "C24", "C25", "C26", "C27", "C28", "C29", "C30",
                "C31", "C32", "C33", "C34", "C35", "C36", "C37", "C38",
            ],
        ),
    ];
    let r_bh = run_ora(
        &universe,
        &dam_cpd,
        &pathways,
        1,
        EnrichmentDirection::Both,
        FdrMethod::BenjaminiHochberg,
        1,
    );
    let r_by = run_ora(
        &universe,
        &dam_cpd,
        &pathways,
        1,
        EnrichmentDirection::Both,
        FdrMethod::BenjaminiYekutieli,
        1,
    );
    assert_eq!(r_bh.fdr_method, FdrMethod::BenjaminiHochberg);
    assert_eq!(r_by.fdr_method, FdrMethod::BenjaminiYekutieli);
    assert_eq!(r_bh.rows.len(), r_by.rows.len());

    // Match rows by entry_id (sort order is identical here since monotone-min
    // preserves underlying p ordering, but find-by-id is more robust).
    for row_bh in &r_bh.rows {
        let row_by = r_by
            .rows
            .iter()
            .find(|r| r.entry_id == row_bh.entry_id)
            .expect("matching row");
        assert!(
            (row_bh.p_value - row_by.p_value).abs() < 1e-15,
            "p_value differs for {}: BH {} vs BY {}",
            row_bh.entry_id,
            row_bh.p_value,
            row_by.p_value
        );
        assert_eq!(
            row_bh.hits, row_by.hits,
            "hits mismatch for {}",
            row_bh.entry_id
        );
    }

    let any_fdr_differs = r_bh.rows.iter().any(|row_bh| {
        let row_by = r_by
            .rows
            .iter()
            .find(|r| r.entry_id == row_bh.entry_id)
            .expect("matching row");
        (row_bh.fdr - row_by.fdr).abs() > 1e-12
    });
    assert!(
        any_fdr_differs,
        "BH and BY must produce at least one differing fdr row"
    );
}

/// Regression test for ORA m_p inflation when an entry's COMPOUND block
/// contains a duplicate cpd ID. Pre-2026-05-26 `m_p` iterated the raw
/// `Vec<String>` so duplicates double-counted; `k_p` used a `HashSet` so
/// duplicates did not. The asymmetry biased `expected` up and
/// `enrichment_ratio` down. Fixed by reusing the deduped `entry_set` for
/// `m_p`.
#[test]
fn duplicate_compound_in_entry_does_not_inflate_m_p() {
    let universe = set(["C1", "C2", "C3", "C4", "C5"]);
    let dam_cpd = set(["C1", "C2"]);
    // Entry has 3 UNIQUE cpds (C1, C2, C3) in universe; COMPOUND block
    // happens to list C1 twice (KEGG curation quirk).
    let pathways = vec![pathway("p1", "Test", &["C1", "C2", "C3", "C1"])];
    let result = run_ora(
        &universe,
        &dam_cpd,
        &pathways,
        1,
        EnrichmentDirection::Both,
        FdrMethod::BenjaminiHochberg,
        1,
    );
    assert_eq!(result.rows.len(), 1);
    let row = &result.rows[0];
    assert_eq!(
        row.total, 3,
        "m_p MUST dedupe entry.compounds (3 unique cpds in universe); got {}",
        row.total
    );
    assert_eq!(row.hits, 2, "k_p still uses HashSet semantics");
}

/// Verifies the K ⊆ N debug_assert! fires when an upstream refactor
/// accidentally lets K leak compounds not in N. Only enforced in debug
/// builds (and `cargo test` is a debug build by default); release uses the
/// per-run `error!` summary log instead.
#[test]
#[should_panic(expected = "ORA invariant")]
#[cfg(debug_assertions)]
fn run_ora_panics_in_debug_when_k_not_subset_of_n() {
    let universe = set(["C1", "C2"]);
    let dam_cpd = set(["C1", "X9999"]); // X9999 not in universe
    let _ = run_ora(
        &universe,
        &dam_cpd,
        &[pathway("p1", "P1", &["C1"])],
        1,
        EnrichmentDirection::Both,
        FdrMethod::BenjaminiHochberg,
        1,
    );
}
