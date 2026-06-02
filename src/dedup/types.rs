//! Public types for the deduplication audit report.

use crate::dedup::adduct::AdductClass;

/// Summary of one `run_dedup` invocation. The kept-index set returned
/// separately and `kept_count + dropped.len() == features.len()` is an
/// invariant.
#[derive(Clone, Debug)]
pub struct DedupReport {
    /// Every feature that lost a cascade comparison, in the order the
    /// pairwise tournament discovered them. Determinism is guaranteed:
    /// `run_dedup` called twice on the same slice produces a
    /// byte-identical `dropped` vector.
    pub dropped: Vec<DroppedFeature>,
    /// Number of features in the kept set returned alongside this
    /// report. Includes both dedup winners (per non-null InChIKey
    /// group) and null-InChIKey passthrough features.
    pub kept_count: usize,
    /// Subset of `kept_count` that consists of features with
    /// `inchikey.is_none()`. These do not participate in any cascade
    /// comparison.
    pub null_inchikey_passthrough: usize,
}

/// One dup-loser, with full provenance: which winner displaced it, at
/// which cascade level, and the values on both sides of that decision.
#[derive(Clone, Debug)]
pub struct DroppedFeature {
    pub alignment_id: String,
    pub inchikey: String,
    pub winner_alignment_id: String,
    pub decided_at: CascadeStep,
    pub loser_value: Option<CascadeValue>,
    pub winner_value: Option<CascadeValue>,
}

/// Which cascade level decided a given pairwise comparison. The
/// `Tiebreak` variant is the deterministic fallback (lexicographic
/// `alignment_id` compare) that guarantees the cascade always
/// terminates.
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CascadeStep {
    MsmsMatched,
    TotalScore,
    AdductClass,
    FillPercent,
    SnAverage,
    Tiebreak,
}

impl CascadeStep {
    /// Variant-name label for the dedup audit CSV's `decided_at` column.
    /// Exhaustive (no wildcard) so adding a `CascadeStep` variant forces a
    /// compile error HERE — at the enum's home — rather than silently
    /// mislabelling the audit CSV from a distant exporter (`move-labels-onto-types`).
    pub fn label(&self) -> &'static str {
        match self {
            CascadeStep::MsmsMatched => "MsmsMatched",
            CascadeStep::TotalScore => "TotalScore",
            CascadeStep::AdductClass => "AdductClass",
            CascadeStep::FillPercent => "FillPercent",
            CascadeStep::SnAverage => "SnAverage",
            CascadeStep::Tiebreak => "Tiebreak",
        }
    }
}

/// The deciding value captured for the audit. The discriminator MUST
/// match the `CascadeStep` per the `msdial-deduplication` spec:
/// numeric steps use `Num(_)`, `MsmsMatched` uses `Msms(_)`,
/// `AdductClass` uses `Adduct { class, sub }`, `Tiebreak` uses
/// `AlignmentId(_)`.
#[derive(Clone, Debug)]
pub enum CascadeValue {
    Num(f64),
    Adduct {
        class: AdductClass,
        sub: u8,
    },
    /// MS/MS-matched boolean as carried by the side whose value was
    /// `Some(_)`. When a side's `ms_ms_matched` was `None`, that side's
    /// `loser_value` / `winner_value` is the outer `Option::None` —
    /// mirroring the numeric levels' `Option<f64>.map(CascadeValue::Num)`
    /// convention. The inner `bool` is therefore always a known value.
    Msms(bool),
    AlignmentId(String),
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cascade_step_label_matches_variant_names() {
        assert_eq!(CascadeStep::MsmsMatched.label(), "MsmsMatched");
        assert_eq!(CascadeStep::TotalScore.label(), "TotalScore");
        assert_eq!(CascadeStep::AdductClass.label(), "AdductClass");
        assert_eq!(CascadeStep::FillPercent.label(), "FillPercent");
        assert_eq!(CascadeStep::SnAverage.label(), "SnAverage");
        assert_eq!(CascadeStep::Tiebreak.label(), "Tiebreak");
    }
}
