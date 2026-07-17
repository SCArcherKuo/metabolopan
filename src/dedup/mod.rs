//! InChIKey + retention-time deduplication for MS-DIAL features.
//!
//! Pure-function entry point `run_dedup` that takes `&[FeatureMeta]` plus an
//! `rt_tolerance_min` window and returns a kept-index set plus a `DedupReport`
//! enumerating every dup-loser. Features are grouped by InChIKey, each group is
//! partitioned into retention-time clusters (complete-linkage: each cluster's RT
//! span stays within the tolerance, so same-InChIKey features more than the
//! tolerance apart are kept as separate peaks), and the cascade picks one
//! survivor per cluster. Retention time
//! decides which features compete, never which one wins. The cascade is
//! documented in `cascade.rs`; adduct classification in `adduct.rs`. See the
//! `msdial-deduplication` capability spec for the formal contract.
//!
//! No I/O, no tracing, no global state. `run_dam` calls this function
//! pre-loop and consumes the kept set as an iteration-time skip; the
//! `MetabolomicsTable` matrices are never resized.

pub mod adduct;
pub mod cascade;
pub mod types;

pub use adduct::AdductClass;
pub use types::{CascadeStep, CascadeValue, DedupReport, DroppedFeature};

use std::collections::{HashMap, HashSet};

use crate::data::types::FeatureMeta;
use crate::dedup::cascade::{DecisionWinner, cascade_compare};

/// Group features by `Some(inchikey)`, partition each group into
/// retention-time clusters (`rt_clusters`), collapse each cluster to a single
/// cascade winner, and report every dup-loser with full provenance. A group
/// that spans more than one retention-time cluster therefore keeps more than
/// one survivor. Features with `inchikey.is_none()` bypass the clustering +
/// cascade entirely and flow into the kept set as `null_inchikey_passthrough`.
///
/// Returns `(kept_indices, report)` where
/// `kept_indices.len() + report.dropped.len() == features.len()` and
/// `report.kept_count == kept_indices.len()`.
///
/// Pure function: deterministic on identical inputs, no I/O, no tracing
/// events. Callers MAY emit tracing around the call (e.g. `run_dam`
/// logs the dropped count post-hoc) but this function does not.
pub fn run_dedup(features: &[FeatureMeta], rt_tolerance_min: f64) -> (HashSet<usize>, DedupReport) {
    let mut kept: HashSet<usize> = HashSet::new();
    let mut dropped: Vec<DroppedFeature> = Vec::new();
    let mut null_inchikey_passthrough: usize = 0;

    // Group non-null InChIKey indices in stable order by first-appearance,
    // and parallel record of group keys for deterministic iteration.
    let mut groups: HashMap<&str, Vec<usize>> = HashMap::new();
    let mut group_order: Vec<&str> = Vec::new();
    for (i, f) in features.iter().enumerate() {
        match f.inchikey.as_deref() {
            Some(k) => {
                let v = groups.entry(k).or_insert_with(|| {
                    group_order.push(k);
                    Vec::new()
                });
                v.push(i);
            }
            None => {
                kept.insert(i);
                null_inchikey_passthrough += 1;
            }
        }
    }

    for key in &group_order {
        let indices = &groups[key];
        // Fast path: a singleton InChIKey group (the common case on real data)
        // is always kept as-is, without the clustering allocation.
        if indices.len() == 1 {
            kept.insert(indices[0]);
            continue;
        }
        // Partition the InChIKey group into retention-time clusters; the cascade
        // runs INDEPENDENTLY per cluster. Complete-linkage bounds each cluster's
        // RT span to `rt_tolerance_min`, so same-InChIKey features more than the
        // tolerance apart land in different clusters and are each kept. See the
        // `msdial-deduplication` capability spec.
        for cluster in rt_clusters(indices, features, rt_tolerance_min) {
            if cluster.len() == 1 {
                kept.insert(cluster[0]);
                continue;
            }
            // Pairwise tournament over the cluster in ascending original-index
            // order: first index is the initial champion. For each subsequent
            // index, compare against the current champion; the loser is
            // recorded with full provenance and the winner becomes the new
            // champion.
            let mut champion = cluster[0];
            for &challenger in &cluster[1..] {
                let (winner, step, val_a, val_b) =
                    cascade_compare(&features[champion], &features[challenger]);
                let (kept_idx, lost_idx, lost_val, winner_val) = match winner {
                    DecisionWinner::A => (champion, challenger, val_b, val_a),
                    DecisionWinner::B => (challenger, champion, val_a, val_b),
                };
                dropped.push(DroppedFeature {
                    alignment_id: features[lost_idx].alignment_id.clone(),
                    inchikey: key.to_string(),
                    winner_alignment_id: features[kept_idx].alignment_id.clone(),
                    decided_at: step,
                    loser_value: lost_val,
                    winner_value: winner_val,
                });
                champion = kept_idx;
            }
            kept.insert(champion);
        }
    }

    let report = DedupReport {
        dropped,
        kept_count: kept.len(),
        null_inchikey_passthrough,
    };
    (kept, report)
}

/// Partition one InChIKey group's feature indices into retention-time clusters.
///
/// Features with `average_rt_min == Some(rt)` are sorted by RT and split by
/// **complete-linkage**, implemented as greedy left-anchored interval covering:
/// each cluster is anchored at its first (minimum) RT and a member joins while
/// `rt - cluster_start <= rt_tolerance_min` (inclusive boundary — a span exactly
/// equal to the tolerance keeps them together). In one dimension a set's span
/// (`max − min`) equals its maximum pairwise distance, so bounding each cluster's
/// span to the tolerance guarantees every pair inside a cluster is within the
/// tolerance — no chain can exceed it. When a boundary-straddling middle feature
/// could join either the lower or the upper cluster, the greedy rule keeps it in
/// the LOWER cluster (a deterministic left bias). Features with
/// `average_rt_min == None` form ONE separate "no-RT" cluster that never merges
/// with any RT-known cluster (unknown RT cannot be asserted to co-elute with a
/// known RT).
///
/// The returned clusters are ordered RT-known-first by ascending
/// cluster-minimum RT (the no-RT cluster, if any, is last), and each cluster's
/// members are in ascending original-feature-index order, so the caller's
/// tournament and `report.dropped` ordering are deterministic.
/// `rt_tolerance_min` is a UI-enforced precondition (finite, `> 0`); a
/// non-finite or `<= 0` value does not panic here (`total_cmp` handles `NaN`
/// and every gap comparison is well-defined) but yields an unspecified
/// partition.
fn rt_clusters(
    indices: &[usize],
    features: &[FeatureMeta],
    rt_tolerance_min: f64,
) -> Vec<Vec<usize>> {
    // Read each RT once and carry it alongside the index for the sort + walk.
    let mut rt_known: Vec<(usize, f64)> = Vec::new();
    let mut no_rt: Vec<usize> = Vec::new();
    for &i in indices {
        match features[i].average_rt_min {
            Some(rt) => rt_known.push((i, rt)),
            None => no_rt.push(i),
        }
    }

    // Ascending RT. The stable sort keeps equal-RT members in input
    // (ascending-index) order; `total_cmp` orders `NaN` without panicking. RT
    // alone is the sort key: cluster membership never depends on any secondary
    // tiebreak (equal RTs share a cluster since their span is 0, and each cluster
    // is re-sorted by original index below).
    rt_known.sort_by(|a, b| a.1.total_cmp(&b.1));

    // Complete-linkage via greedy left-anchoring: `cluster_start` is the RT of
    // the current cluster's first (minimum) member and is updated ONLY when a
    // new cluster opens. A member joins while its span from that anchor is within
    // the tolerance, bounding every cluster's span — and thus every pairwise
    // distance — to `rt_tolerance_min`. Keeping the anchor fixed while members
    // join is what packs a boundary-straddling feature into the LOWER cluster
    // (the deterministic left bias). A NaN or negative tolerance makes the join
    // test `false`, so each member re-anchors as its own cluster (panic-free,
    // unspecified).
    let mut clusters: Vec<Vec<usize>> = Vec::new();
    let mut cluster_start: Option<f64> = None;
    for &(i, rt) in &rt_known {
        match cluster_start {
            Some(start) if rt - start <= rt_tolerance_min => {
                clusters.last_mut().unwrap().push(i);
            }
            _ => {
                clusters.push(vec![i]);
                cluster_start = Some(rt);
            }
        }
    }

    // The tournament runs in ascending original-index order within each cluster;
    // the clusters themselves are already in ascending cluster-minimum RT order
    // (they were formed from the RT-sorted walk).
    for cluster in &mut clusters {
        cluster.sort_unstable();
    }
    // The no-RT members form one separate cluster, processed last.
    if !no_rt.is_empty() {
        no_rt.sort_unstable();
        clusters.push(no_rt);
    }
    clusters
}

#[cfg(test)]
mod tests {
    use super::*;

    fn feat(alignment_id: &str, inchikey: Option<&str>) -> FeatureMeta {
        FeatureMeta {
            alignment_id: alignment_id.to_string(),
            metabolite_name: String::new(),
            inchikey: inchikey.map(|s| s.to_string()),
            adduct_type: None,
            average_rt_min: None,
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

    fn feat_full(
        alignment_id: &str,
        inchikey: Option<&str>,
        ms_ms: Option<bool>,
        total: Option<f64>,
        adduct: Option<&str>,
        fill: Option<f64>,
    ) -> FeatureMeta {
        FeatureMeta {
            alignment_id: alignment_id.to_string(),
            metabolite_name: String::new(),
            inchikey: inchikey.map(|s| s.to_string()),
            adduct_type: adduct.map(|s| s.to_string()),
            average_rt_min: None,
            average_mz: None,
            formula: None,
            smiles: None,
            fill_percent: fill,
            ms_ms_matched: ms_ms,
            isotope_tracking_weight_number: None,
            total_score: total,
            sn_average: None,
        }
    }

    #[test]
    fn empty_input_empty_output() {
        let (kept, report) = run_dedup(&[], 0.1);
        assert!(kept.is_empty());
        assert!(report.dropped.is_empty());
        assert_eq!(report.kept_count, 0);
        assert_eq!(report.null_inchikey_passthrough, 0);
    }

    #[test]
    fn all_null_inchikey_passthrough() {
        let features: Vec<FeatureMeta> = (0..5).map(|i| feat(&format!("P{i}"), None)).collect();
        let (kept, report) = run_dedup(&features, 0.1);
        assert_eq!(kept.len(), 5);
        assert!(report.dropped.is_empty());
        assert_eq!(report.null_inchikey_passthrough, 5);
        assert_eq!(report.kept_count, 5);
    }

    #[test]
    fn singleton_kept_regardless_of_quality() {
        // A single feature with bad-quality annotation is still kept —
        // cascade-only semantics. The msms=false / isotope-pattern adduct
        // would lose if a competitor existed, but as a singleton it
        // passes through untouched.
        let features = vec![feat_full(
            "P0",
            Some("ABCD"),
            Some(false),
            Some(100.0),
            Some("[M+1]+"),
            Some(10.0),
        )];
        let (kept, report) = run_dedup(&features, 0.1);
        assert_eq!(kept.len(), 1);
        assert!(kept.contains(&0));
        assert!(report.dropped.is_empty());
        assert_eq!(report.kept_count, 1);
    }

    #[test]
    fn three_way_group_one_survivor() {
        // 3 features with InChIKey "X":
        //   A: msms=true, total=800
        //   B: msms=true, total=900   <- best on Total score
        //   C: msms=false             <- worst on MS/MS
        // Tournament: champion=A, challenge B -> B wins (TotalScore, A loses).
        //             champion=B, challenge C -> B wins (MsmsMatched, C loses).
        let features = vec![
            feat_full("A", Some("X"), Some(true), Some(800.0), None, None),
            feat_full("B", Some("X"), Some(true), Some(900.0), None, None),
            feat_full("C", Some("X"), Some(false), Some(950.0), None, None),
        ];
        let (kept, report) = run_dedup(&features, 0.1);
        assert_eq!(kept.len(), 1);
        assert!(kept.contains(&1));
        assert_eq!(report.dropped.len(), 2);
        // First dropped is A (lost to B on TotalScore).
        assert_eq!(report.dropped[0].alignment_id, "A");
        assert_eq!(report.dropped[0].winner_alignment_id, "B");
        assert_eq!(report.dropped[0].decided_at, CascadeStep::TotalScore);
        // Second dropped is C (lost to B on MsmsMatched).
        assert_eq!(report.dropped[1].alignment_id, "C");
        assert_eq!(report.dropped[1].winner_alignment_id, "B");
        assert_eq!(report.dropped[1].decided_at, CascadeStep::MsmsMatched);
        assert_eq!(report.kept_count, 1);
    }

    #[test]
    fn determinism_on_repeated_calls() {
        let features = vec![
            feat_full(
                "P0",
                Some("X"),
                Some(true),
                Some(800.0),
                Some("[M+H]+"),
                Some(95.0),
            ),
            feat_full(
                "P1",
                Some("X"),
                Some(true),
                Some(700.0),
                Some("[M+Na]+"),
                Some(70.0),
            ),
            feat_full(
                "P2",
                Some("Y"),
                Some(false),
                Some(500.0),
                Some("[M+H]+"),
                Some(50.0),
            ),
            feat_full(
                "P3",
                Some("Y"),
                Some(false),
                Some(600.0),
                Some("[M+H]+"),
                Some(60.0),
            ),
            feat_full(
                "P4", None, None, None, None, None, // null InChIKey, passthrough
            ),
        ];
        let (kept1, report1) = run_dedup(&features, 0.1);
        let (kept2, report2) = run_dedup(&features, 0.1);
        // HashSet equality is content-based; reports are Vec, so we need
        // their `dropped` to match in order.
        assert_eq!(kept1, kept2);
        assert_eq!(report1.dropped.len(), report2.dropped.len());
        for (a, b) in report1.dropped.iter().zip(report2.dropped.iter()) {
            assert_eq!(a.alignment_id, b.alignment_id);
            assert_eq!(a.winner_alignment_id, b.winner_alignment_id);
            assert_eq!(a.decided_at, b.decided_at);
        }
        assert_eq!(report1.kept_count, report2.kept_count);
        assert_eq!(
            report1.null_inchikey_passthrough,
            report2.null_inchikey_passthrough
        );
    }

    #[test]
    fn counts_add_up() {
        // 10 features: 3 groups of 2 (so 3 dropped) + 4 singletons.
        let features = vec![
            feat_full("A1", Some("X"), Some(true), Some(900.0), None, None),
            feat_full("A2", Some("X"), Some(true), Some(800.0), None, None),
            feat_full("B1", Some("Y"), Some(true), Some(900.0), None, None),
            feat_full("B2", Some("Y"), Some(true), Some(800.0), None, None),
            feat_full("C1", Some("Z"), Some(true), Some(900.0), None, None),
            feat_full("C2", Some("Z"), Some(true), Some(800.0), None, None),
            feat_full("S1", Some("S1k"), None, None, None, None),
            feat_full("S2", Some("S2k"), None, None, None, None),
            feat_full("S3", Some("S3k"), None, None, None, None),
            feat_full("N1", None, None, None, None, None),
        ];
        let (kept, report) = run_dedup(&features, 0.1);
        // 6 dup features -> 3 winners + 3 losers; 3 singletons + 1 null = 7 kept
        assert_eq!(kept.len(), 7);
        assert_eq!(report.dropped.len(), 3);
        assert_eq!(report.kept_count, 7);
        assert_eq!(report.null_inchikey_passthrough, 1);
        // Invariant
        assert_eq!(report.kept_count + report.dropped.len(), features.len());
    }

    #[test]
    fn null_inchikey_does_not_dedup_with_other_nulls() {
        // Two features with InChIKey == None should BOTH be kept, not
        // grouped as duplicates of each other.
        let features = vec![
            feat_full("N1", None, Some(true), Some(900.0), None, None),
            feat_full("N2", None, Some(false), Some(800.0), None, None),
        ];
        let (kept, report) = run_dedup(&features, 0.1);
        assert_eq!(kept.len(), 2);
        assert!(kept.contains(&0) && kept.contains(&1));
        assert!(report.dropped.is_empty());
        assert_eq!(report.null_inchikey_passthrough, 2);
    }

    // --- Retention-time clustering (add-rt-aware-dedup) ---

    /// A feature carrying an explicit `average_rt_min`; every cascade field is
    /// `None`, so a same-cluster tie resolves at level 4 (`alignment_id`).
    fn feat_rt(alignment_id: &str, inchikey: Option<&str>, rt: Option<f64>) -> FeatureMeta {
        let mut f = feat(alignment_id, inchikey);
        f.average_rt_min = rt;
        f
    }

    #[test]
    fn rt_within_tolerance_collapses_to_one_winner() {
        // Same InChIKey, RT 2.10 vs 2.13 (gap 0.03 <= 0.1) -> one cluster.
        let features = vec![
            feat_rt("A", Some("AAA"), Some(2.10)),
            feat_rt("B", Some("AAA"), Some(2.13)),
        ];
        let (kept, report) = run_dedup(&features, 0.1);
        assert_eq!(kept.len(), 1);
        assert_eq!(report.dropped.len(), 1);
        // Cascade fields all None -> alignment_id tiebreak; "A" < "B" keeps A.
        assert!(kept.contains(&0));
    }

    #[test]
    fn rt_beyond_tolerance_keeps_both() {
        // Same InChIKey, RT 2.10 vs 8.55 (gap 6.45 > 0.1) -> two clusters.
        let features = vec![
            feat_rt("A", Some("AAA"), Some(2.10)),
            feat_rt("B", Some("AAA"), Some(8.55)),
        ];
        let (kept, report) = run_dedup(&features, 0.1);
        assert_eq!(kept.len(), 2);
        assert!(kept.contains(&0) && kept.contains(&1));
        assert!(report.dropped.is_empty());
    }

    #[test]
    fn rt_boundary_gap_exactly_at_tolerance_stays_together() {
        // RT 2.0 vs 2.5 with tolerance 0.5 -> span == tolerance (exact in f64,
        // all dyadic) -> inclusive boundary keeps them in one cluster.
        let features = vec![
            feat_rt("A", Some("AAA"), Some(2.0)),
            feat_rt("B", Some("AAA"), Some(2.5)),
        ];
        let (kept, report) = run_dedup(&features, 0.5);
        assert_eq!(kept.len(), 1);
        assert_eq!(report.dropped.len(), 1);
    }

    #[test]
    fn rt_chain_wider_than_tolerance_is_split() {
        // Complete-linkage: a chain 0.00 / 0.08 / 0.16 has consecutive gaps
        // 0.08 ≤ 0.1 but a total span 0.16 > 0.1, so it does NOT merge into one
        // cluster (single-linkage would). Greedy left-anchoring yields
        // {A(0.00), B(0.08)} and {C(0.16)}, packing the middle 0.08 into the
        // LOWER cluster (left bias) rather than pairing it with 0.16.
        let features = vec![
            feat_rt("A", Some("AAA"), Some(0.00)),
            feat_rt("B", Some("AAA"), Some(0.08)),
            feat_rt("C", Some("AAA"), Some(0.16)),
        ];
        let (kept, report) = run_dedup(&features, 0.1);
        // {A,B} -> one survivor + one dropped; {C} -> kept as its own cluster.
        assert_eq!(report.dropped.len(), 1);
        assert_eq!(kept.len(), 2);
        assert!(kept.contains(&2)); // C (0.16) kept separately
        // Left bias: the single drop is inside the lower {A,B} cluster; C never
        // enters a two-member cluster, so it is neither the loser nor the winner.
        assert_ne!(report.dropped[0].alignment_id, "C");
        assert_ne!(report.dropped[0].winner_alignment_id, "C");
    }

    #[test]
    fn rt_none_features_form_one_separate_cluster() {
        // Two no-RT + one RT-known, same InChIKey: the two no-RT compete in one
        // cluster (one dropped); the RT-known feature is its own singleton (kept).
        let features = vec![
            feat_rt("A", Some("AAA"), None),
            feat_rt("B", Some("AAA"), None),
            feat_rt("C", Some("AAA"), Some(3.0)),
        ];
        let (kept, report) = run_dedup(&features, 0.1);
        assert_eq!(report.dropped.len(), 1);
        assert_eq!(kept.len(), 2);
        assert!(kept.contains(&0)); // A kept (alignment_id "A" < "B")
        assert!(kept.contains(&2)); // C (RT-known) always kept as its own cluster
        assert!(!kept.contains(&1)); // B dropped
    }

    #[test]
    fn all_none_rt_is_tolerance_invariant() {
        // Every member has average_rt_min == None -> one no-RT cluster regardless
        // of tolerance; the result is identical across tolerances (legacy
        // InChIKey-only behavior, byte-for-byte).
        let features = vec![
            feat_full("A", Some("X"), Some(true), Some(800.0), None, None),
            feat_full("B", Some("X"), Some(true), Some(900.0), None, None),
            feat_full("C", Some("X"), Some(false), Some(950.0), None, None),
        ];
        let (kept_small, rep_small) = run_dedup(&features, 0.001);
        let (kept_large, rep_large) = run_dedup(&features, f64::MAX);
        assert_eq!(kept_small, kept_large);
        assert_eq!(rep_small.dropped.len(), rep_large.dropped.len());
        assert_eq!(kept_small.len(), 1);
        assert!(kept_small.contains(&1)); // B wins on Total score among ms_ms=true
    }

    #[test]
    fn large_tolerance_merges_rt_known_but_keeps_no_rt_separate() {
        // Mixed group: two RT-known (far apart) + one no-RT. f64::MAX merges the
        // two RT-known into one cluster (one dropped) but the no-RT stays its own
        // cluster (kept) -- so a mixed group does NOT collapse to a single group.
        let features = vec![
            feat_rt("A", Some("AAA"), Some(1.0)),
            feat_rt("B", Some("AAA"), Some(9.0)),
            feat_rt("C", Some("AAA"), None),
        ];
        let (kept, report) = run_dedup(&features, f64::MAX);
        assert_eq!(report.dropped.len(), 1);
        assert_eq!(kept.len(), 2);
        assert!(kept.contains(&2)); // C (no-RT) always kept separately
    }

    #[test]
    fn equal_rt_determinism_by_alignment_id() {
        // Identical RT within one cluster; the tie resolves at alignment_id and
        // repeated calls agree.
        let features = vec![
            feat_rt("PEAK_0002", Some("AAA"), Some(2.0)),
            feat_rt("PEAK_0001", Some("AAA"), Some(2.0)),
        ];
        let (kept1, rep1) = run_dedup(&features, 0.1);
        let (kept2, rep2) = run_dedup(&features, 0.1);
        assert_eq!(kept1, kept2);
        assert_eq!(rep1.dropped.len(), rep2.dropped.len());
        assert_eq!(kept1.len(), 1);
        assert!(kept1.contains(&1)); // "PEAK_0001" (index 1) is lexicographically smaller
    }

    #[test]
    fn non_finite_or_nonpositive_tolerance_is_panic_free() {
        // run_dedup is a public pure function; a non-UI caller could pass a
        // non-finite or <= 0 tolerance. It must not panic; the count invariant
        // holds regardless (partition is otherwise unspecified).
        let features = vec![
            feat_rt("A", Some("AAA"), Some(2.0)),
            feat_rt("B", Some("AAA"), Some(2.05)),
            feat_rt("C", Some("AAA"), None),
        ];
        for tol in [f64::NAN, -1.0, 0.0, f64::INFINITY] {
            let (kept, report) = run_dedup(&features, tol);
            assert_eq!(kept.len() + report.dropped.len(), features.len());
        }
    }
}
