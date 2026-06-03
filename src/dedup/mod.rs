//! InChIKey deduplication for MS-DIAL features.
//!
//! Pure-function entry point `run_dedup` that takes `&[FeatureMeta]`
//! and returns a kept-index set plus a `DedupReport` enumerating every
//! dup-loser. The cascade is documented in `cascade.rs`; adduct
//! classification in `adduct.rs`. See the `msdial-deduplication` capability
//! spec for the formal contract.
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

/// Group features by `Some(inchikey)`, collapse each group to a single
/// cascade winner, and report every dup-loser with full provenance.
/// Features with `inchikey.is_none()` bypass the cascade entirely and
/// flow into the kept set as `null_inchikey_passthrough`.
///
/// Returns `(kept_indices, report)` where
/// `kept_indices.len() + report.dropped.len() == features.len()` and
/// `report.kept_count == kept_indices.len()`.
///
/// Pure function: deterministic on identical inputs, no I/O, no tracing
/// events. Callers MAY emit tracing around the call (e.g. `run_dam`
/// logs the dropped count post-hoc) but this function does not.
pub fn run_dedup(features: &[FeatureMeta]) -> (HashSet<usize>, DedupReport) {
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
        if indices.len() == 1 {
            kept.insert(indices[0]);
            continue;
        }
        // Pairwise tournament: first index is the initial champion. For
        // each subsequent index, compare against the current champion;
        // the loser is recorded with full provenance and the winner
        // becomes the new champion.
        let mut champion = indices[0];
        for &challenger in &indices[1..] {
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

    let report = DedupReport {
        dropped,
        kept_count: kept.len(),
        null_inchikey_passthrough,
    };
    (kept, report)
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
        let (kept, report) = run_dedup(&[]);
        assert!(kept.is_empty());
        assert!(report.dropped.is_empty());
        assert_eq!(report.kept_count, 0);
        assert_eq!(report.null_inchikey_passthrough, 0);
    }

    #[test]
    fn all_null_inchikey_passthrough() {
        let features: Vec<FeatureMeta> = (0..5).map(|i| feat(&format!("P{i}"), None)).collect();
        let (kept, report) = run_dedup(&features);
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
        let (kept, report) = run_dedup(&features);
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
        let (kept, report) = run_dedup(&features);
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
        let (kept1, report1) = run_dedup(&features);
        let (kept2, report2) = run_dedup(&features);
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
        let (kept, report) = run_dedup(&features);
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
        let (kept, report) = run_dedup(&features);
        assert_eq!(kept.len(), 2);
        assert!(kept.contains(&0) && kept.contains(&1));
        assert!(report.dropped.is_empty());
        assert_eq!(report.null_inchikey_passthrough, 2);
    }
}
