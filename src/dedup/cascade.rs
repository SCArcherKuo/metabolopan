//! Pairwise cascade comparator for the deduplication module.
//!
//! Implements the 6-level decision table specified by the
//! `msdial-deduplication` capability:
//!
//! | Level | Field                              | Rule                                                          |
//! |-------|------------------------------------|---------------------------------------------------------------|
//! | 1a    | `ms_ms_matched`                    | `Some(true)` > `Some(false)` > `None`                         |
//! | 1b    | `total_score`                      | larger wins (`f64::total_cmp`); `Some(_)` beats `None`        |
//! | 2     | adduct class                       | `Primary` < `NonPrimary` < `Dimer` < `Isotope`; Primary sub-rank |
//! | 3a    | `fill_percent`                     | same rule as `total_score`                                   |
//! | 3b    | `sn_average`                       | same rule as `total_score`                                   |
//! | 4     | `alignment_id`                     | lexicographic smaller wins (deterministic terminator)         |
//!
//! The MS-DIAL `Dot product` column was removed from the cascade: `Total score`
//! is the vendor-computed weighted composite of every spectral-similarity
//! metric (including dot products), so ranking on the raw dot product as well
//! double-counted the same signal.
//!
//! Pure function. No I/O, no tracing.

use std::cmp::Ordering;

use crate::data::types::FeatureMeta;
use crate::dedup::adduct::{AdductClass, classify, primary_subrank};
use crate::dedup::types::{CascadeStep, CascadeValue};

/// Which of the two input features won the pairwise cascade.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum DecisionWinner {
    A,
    B,
}

/// Compare two same-InChIKey features through the 6-level cascade.
/// Returns the winner, the step that decided, and the recorded value
/// on each side (for the audit report).
///
/// The cascade always terminates because level 4 (`alignment_id`) is
/// guaranteed to distinguish any two features — if even alignment IDs
/// match the inputs are bit-identical and the caller bears the choice
/// of which one to keep (we pick `A` and report `Tiebreak`).
pub fn cascade_compare(
    a: &FeatureMeta,
    b: &FeatureMeta,
) -> (
    DecisionWinner,
    CascadeStep,
    Option<CascadeValue>,
    Option<CascadeValue>,
) {
    // Level 1a — MS/MS matched. Three-way ordering: Some(true) > Some(false) > None.
    // Mirror the numeric levels' `.map(CascadeValue::Num)` convention: when a
    // side's `ms_ms_matched` was `None`, that side's value is the OUTER
    // `Option::None`, not `Some(Msms(None))`. The spec's "loser had None →
    // loser_value is None" rule applies uniformly across all cascade levels.
    if let Some(decision) = compare_msms(a.ms_ms_matched, b.ms_ms_matched) {
        return (
            decision,
            CascadeStep::MsmsMatched,
            a.ms_ms_matched.map(CascadeValue::Msms),
            b.ms_ms_matched.map(CascadeValue::Msms),
        );
    }
    // Level 1b — Total score (vendor-computed weighted composite of every
    // spectral-similarity metric, including dot products). Larger wins.
    if let Some(decision) = compare_optional_f64(a.total_score, b.total_score) {
        return (
            decision,
            CascadeStep::TotalScore,
            a.total_score.map(CascadeValue::Num),
            b.total_score.map(CascadeValue::Num),
        );
    }
    // Level 2 — Adduct class (with Primary sub-rank).
    let (class_a, sub_a) = adduct_rank(a);
    let (class_b, sub_b) = adduct_rank(b);
    if let Some(decision) = compare_adduct_rank((class_a, sub_a), (class_b, sub_b)) {
        return (
            decision,
            CascadeStep::AdductClass,
            Some(CascadeValue::Adduct {
                class: class_a,
                sub: sub_a,
            }),
            Some(CascadeValue::Adduct {
                class: class_b,
                sub: sub_b,
            }),
        );
    }
    // Level 3a — Fill %.
    if let Some(decision) = compare_optional_f64(a.fill_percent, b.fill_percent) {
        return (
            decision,
            CascadeStep::FillPercent,
            a.fill_percent.map(CascadeValue::Num),
            b.fill_percent.map(CascadeValue::Num),
        );
    }
    // Level 3b — S/N average.
    if let Some(decision) = compare_optional_f64(a.sn_average, b.sn_average) {
        return (
            decision,
            CascadeStep::SnAverage,
            a.sn_average.map(CascadeValue::Num),
            b.sn_average.map(CascadeValue::Num),
        );
    }
    // Level 4 — alignment_id (deterministic terminator).
    let winner = match a.alignment_id.cmp(&b.alignment_id) {
        Ordering::Less => DecisionWinner::A,
        // Equal includes the truly-identical case; we deterministically
        // keep A so the cascade always returns a decision.
        Ordering::Equal => DecisionWinner::A,
        Ordering::Greater => DecisionWinner::B,
    };
    (
        winner,
        CascadeStep::Tiebreak,
        Some(CascadeValue::AlignmentId(a.alignment_id.clone())),
        Some(CascadeValue::AlignmentId(b.alignment_id.clone())),
    )
}

fn compare_msms(a: Option<bool>, b: Option<bool>) -> Option<DecisionWinner> {
    // Encode as a u8 so we can use direct comparison: Some(true)=2,
    // Some(false)=1, None=0. Larger wins.
    let rank = |v: Option<bool>| match v {
        Some(true) => 2u8,
        Some(false) => 1u8,
        None => 0u8,
    };
    let ra = rank(a);
    let rb = rank(b);
    match ra.cmp(&rb) {
        Ordering::Greater => Some(DecisionWinner::A),
        Ordering::Less => Some(DecisionWinner::B),
        Ordering::Equal => None,
    }
}

fn compare_optional_f64(a: Option<f64>, b: Option<f64>) -> Option<DecisionWinner> {
    match (a, b) {
        (Some(x), Some(y)) => match x.total_cmp(&y) {
            Ordering::Greater => Some(DecisionWinner::A),
            Ordering::Less => Some(DecisionWinner::B),
            Ordering::Equal => None,
        },
        (Some(_), None) => Some(DecisionWinner::A),
        (None, Some(_)) => Some(DecisionWinner::B),
        (None, None) => None,
    }
}

fn adduct_rank(f: &FeatureMeta) -> (AdductClass, u8) {
    let class = classify(f.adduct_type.as_deref(), f.isotope_tracking_weight_number);
    let sub = match class {
        AdductClass::Primary => f
            .adduct_type
            .as_deref()
            .map(primary_subrank)
            // shouldn't happen — classify returns Primary only when adduct
            // string is Some and in the allowlist, but be safe.
            .unwrap_or(1),
        // For non-Primary classes the sub-rank has no defined meaning;
        // store 0 for stability.
        _ => 0,
    };
    (class, sub)
}

fn compare_adduct_rank(
    (class_a, sub_a): (AdductClass, u8),
    (class_b, sub_b): (AdductClass, u8),
) -> Option<DecisionWinner> {
    // Lower class rank wins. Class order: Primary < NonPrimary < Dimer < Isotope.
    let class_rank = |c: AdductClass| match c {
        AdductClass::Primary => 0u8,
        AdductClass::NonPrimary => 1u8,
        AdductClass::Dimer => 2u8,
        AdductClass::Isotope => 3u8,
    };
    let ra = class_rank(class_a);
    let rb = class_rank(class_b);
    match ra.cmp(&rb) {
        Ordering::Less => Some(DecisionWinner::A),
        Ordering::Greater => Some(DecisionWinner::B),
        Ordering::Equal => {
            // Same class — only Primary has a meaningful sub-rank.
            if class_a == AdductClass::Primary {
                match sub_a.cmp(&sub_b) {
                    Ordering::Less => Some(DecisionWinner::A),
                    Ordering::Greater => Some(DecisionWinner::B),
                    Ordering::Equal => None,
                }
            } else {
                None
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Build a feature with explicit cascade fields and everything else default.
    #[allow(clippy::too_many_arguments)]
    fn feat(
        alignment_id: &str,
        ms_ms: Option<bool>,
        total: Option<f64>,
        adduct: Option<&str>,
        iso_w: Option<i32>,
        fill: Option<f64>,
        sn: Option<f64>,
    ) -> FeatureMeta {
        FeatureMeta {
            alignment_id: alignment_id.to_string(),
            metabolite_name: String::new(),
            inchikey: None,
            adduct_type: adduct.map(|s| s.to_string()),
            average_rt_min: None,
            average_mz: None,
            formula: None,
            smiles: None,
            fill_percent: fill,
            ms_ms_matched: ms_ms,
            isotope_tracking_weight_number: iso_w,
            total_score: total,
            sn_average: sn,
        }
    }

    #[test]
    fn level_1a_ms_ms_true_beats_false() {
        let a = feat("A", Some(true), None, None, None, None, None);
        let b = feat("B", Some(false), None, None, None, None, None);
        let (winner, step, _, _) = cascade_compare(&a, &b);
        assert_eq!(winner, DecisionWinner::A);
        assert_eq!(step, CascadeStep::MsmsMatched);
    }

    #[test]
    fn level_1b_total_score_wins_on_ms_ms_tie() {
        // ms_ms_matched ties; total_score decides. total_score is the level-1b
        // numeric tiebreaker after the Dot product column was removed.
        let a = feat("A", Some(true), Some(0.85), None, None, None, None);
        let b = feat("B", Some(true), Some(0.72), None, None, None, None);
        let (winner, step, _, _) = cascade_compare(&a, &b);
        assert_eq!(winner, DecisionWinner::A);
        assert_eq!(step, CascadeStep::TotalScore);
    }

    #[test]
    fn level_2_adduct_class_beats_fill_percent() {
        // Both have all-equal upper fields; one is Primary [M+H]+, the
        // other is Dimer [2M+H]+. The Primary one wins on class even
        // though the Dimer one has higher Fill %.
        let a = feat(
            "A",
            Some(true),
            Some(0.85),
            Some("[M+H]+"),
            None,
            Some(20.0),
            None,
        );
        let b = feat(
            "B",
            Some(true),
            Some(0.85),
            Some("[2M+H]+"),
            None,
            Some(99.0),
            None,
        );
        let (winner, step, _, _) = cascade_compare(&a, &b);
        assert_eq!(winner, DecisionWinner::A);
        assert_eq!(step, CascadeStep::AdductClass);
    }

    #[test]
    fn level_3a_fill_percent_wins_on_adduct_tie() {
        let a = feat(
            "A",
            Some(true),
            Some(0.85),
            Some("[M+H]+"),
            None,
            Some(95.0),
            None,
        );
        let b = feat(
            "B",
            Some(true),
            Some(0.85),
            Some("[M+H]+"),
            None,
            Some(40.0),
            None,
        );
        let (winner, step, _, _) = cascade_compare(&a, &b);
        assert_eq!(winner, DecisionWinner::A);
        assert_eq!(step, CascadeStep::FillPercent);
    }

    #[test]
    fn level_3b_sn_breaks_fill_percent_tie() {
        let a = feat(
            "A",
            Some(true),
            Some(0.85),
            Some("[M+H]+"),
            None,
            Some(90.0),
            Some(120.0),
        );
        let b = feat(
            "B",
            Some(true),
            Some(0.85),
            Some("[M+H]+"),
            None,
            Some(90.0),
            Some(45.0),
        );
        let (winner, step, _, _) = cascade_compare(&a, &b);
        assert_eq!(winner, DecisionWinner::A);
        assert_eq!(step, CascadeStep::SnAverage);
    }

    #[test]
    fn level_4_alignment_id_lexicographic_terminator() {
        let a = feat("PEAK_0098", None, None, Some("[M+Na]+"), None, None, None);
        let b = feat("PEAK_0123", None, None, Some("[M+Na]+"), None, None, None);
        let (winner, step, _, _) = cascade_compare(&a, &b);
        assert_eq!(winner, DecisionWinner::A);
        assert_eq!(step, CascadeStep::Tiebreak);
    }

    #[test]
    fn some_beats_none_at_each_nullable_level() {
        // A has Some on every cascade field, B has None on all four.
        // Decision should short-circuit at level 1a (ms_ms_matched).
        let a = feat(
            "A",
            Some(true),
            Some(0.85),
            Some("[M+H]+"),
            Some(0),
            Some(90.0),
            Some(120.0),
        );
        let b = feat("B", None, None, None, None, None, None);
        let (winner, step, _, _) = cascade_compare(&a, &b);
        assert_eq!(winner, DecisionWinner::A);
        assert_eq!(step, CascadeStep::MsmsMatched);
    }

    #[test]
    fn total_score_none_loser_produces_none_loser_value() {
        // Spec scenario: "loser_value is None when loser was None on the
        // deciding field" — winner has Some(total_score), loser has None.
        // cascade.rs uses `.map(CascadeValue::Num)` so this yields outer None.
        let winner = feat("W", Some(true), Some(0.80), None, None, None, None);
        let loser = feat("L", Some(true), None, None, None, None, None);
        let (w, step, val_w, val_l) = cascade_compare(&winner, &loser);
        assert_eq!(w, DecisionWinner::A);
        assert_eq!(step, CascadeStep::TotalScore);
        assert!(
            val_l.is_none(),
            "loser_value must be None when loser had no total_score; got {val_l:?}"
        );
        assert!(
            matches!(val_w, Some(CascadeValue::Num(v)) if (v - 0.80).abs() < 1e-9),
            "winner_value must be Some(Num(0.80)); got {val_w:?}"
        );
    }

    #[test]
    fn primary_sub_rank_breaks_class_tie() {
        // Both Primary, but A is [M+H]+ (sub 0) and B is [M+Na]+ (sub 1).
        let a = feat(
            "A",
            Some(true),
            Some(0.85),
            Some("[M+H]+"),
            None,
            Some(90.0),
            Some(120.0),
        );
        let b = feat(
            "B",
            Some(true),
            Some(0.85),
            Some("[M+Na]+"),
            None,
            Some(90.0),
            Some(120.0),
        );
        let (winner, step, _, _) = cascade_compare(&a, &b);
        assert_eq!(winner, DecisionWinner::A);
        assert_eq!(step, CascadeStep::AdductClass);
    }
}
