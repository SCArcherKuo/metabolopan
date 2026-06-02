//! Public types for the `normalize` module.

use std::fmt;

use serde::{Deserialize, Serialize};

/// Sample-axis normalization methods. `None` is the bit-equal passthrough used
/// by default; the other five methods correspond to the user-selectable
/// options in the Stage 2 setup screen.
///
/// Serde derives use externally-tagged representation (unit variants as bare
/// strings, struct variants as `{"Variant": {...}}`). This shape is part of
/// the `session-settings-io` capability's JSON schema contract — renaming a
/// variant requires a `SCHEMA_VERSION` bump there.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum NormalizationMethod {
    None,
    Sum,
    Median,
    Metadata { column: String },
    Quantile,
    Pqn { reference: PqnReference },
}

/// Where the PQN reference spectrum comes from. `AllSamples` mirrors the
/// classical Dieterle 2006 formulation but excludes `Unassigned` samples
/// (so the reference cohort stays consistent with how the rest of the app
/// treats Unassigned). `Group(name)` restricts the reference to one group.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum PqnReference {
    AllSamples,
    Group(String),
}

/// Configuration passed to [`crate::normalize::apply`]. A thin wrapper so
/// `&NormalizationConfig` is the call-site type rather than the bare enum;
/// future per-method tuning knobs (e.g. magnitude-preservation off) can land
/// here without breaking callers.
#[derive(Debug, Clone)]
pub struct NormalizationConfig {
    pub method: NormalizationMethod,
}

impl Default for NormalizationConfig {
    fn default() -> Self {
        Self {
            method: NormalizationMethod::None,
        }
    }
}

/// Errors that abort a normalization before the working matrix is returned.
/// Display messages are written to be rendered directly in the Stage 2 error
/// banner, so each variant names the offending sample / column / group /
/// value where applicable.
#[derive(Debug, Clone, PartialEq)]
pub enum NormalizationError {
    EmptyMatrix,
    ZeroFactor {
        sample: String,
        method: &'static str,
    },
    NanFactor {
        sample: String,
        method: &'static str,
    },
    MetadataColumnMissing {
        column: String,
    },
    MetadataValueMissing {
        sample: String,
        column: String,
    },
    MetadataValueNonPositive {
        sample: String,
        column: String,
        value: f64,
    },
    /// Surfaced by `apply_metadata` when a sample named in the per-mode
    /// `sample_cols` is absent from `mapping.sample_names` (e.g. a stale
    /// (mapping, tables) pair after a Stage 1 file re-pick). Returned instead
    /// of panicking so the skew surfaces as a recoverable Stage 2 banner error
    /// rather than aborting the GUI (`convert-defensive-panics-to-errors`).
    SampleNotInMapping {
        sample: String,
        column: String,
    },
    EmptyReferenceGroup {
        group: String,
    },
    ReferenceAllNan {
        method: &'static str,
    },
    /// Surfaced by the Metadata-column preflight when dropping samples
    /// without a metadata value would leave the chosen numerator or
    /// denominator group with fewer than 2 samples (DAM's per-group
    /// minimum). Names the offending group, how many samples remain, the
    /// minimum required, and the column the user picked.
    InsufficientSamplesAfterDrop {
        group: String,
        remaining: usize,
        required: usize,
        column: String,
    },
    /// Surfaced by `apply_pqn` when one or more samples produce a
    /// degenerate per-sample quotient median (NaN — no usable features
    /// against the reference spectrum, OR 0 — half-or-more of usable
    /// quotients are exactly 0, indicating a sparse / blank-like sample).
    /// Pre-2026-05-26 the implementation silently fell back to
    /// `factor = 1.0` for these samples (leaving them at sum-normalized
    /// scale while peers were PQN-scaled — producing artefactual
    /// differential abundance from scale mismatch). Now this is a hard
    /// error so the user decides whether to drop the offending samples or
    /// switch normalization method.
    PqnDegenerateSamples {
        samples: Vec<String>,
    },
}

impl fmt::Display for NormalizationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptyMatrix => write!(f, "Cannot normalize: intensity matrix is empty."),
            Self::ZeroFactor { sample, method } => write!(
                f,
                "Sample '{sample}' has a {method} of 0 — cannot apply {method} normalization."
            ),
            Self::NanFactor { sample, method } => write!(
                f,
                "Sample '{sample}' has no non-missing values — cannot apply {method} normalization."
            ),
            Self::MetadataColumnMissing { column } => write!(
                f,
                "Metadata column '{column}' is not present in the loaded metadata CSV."
            ),
            Self::MetadataValueMissing { sample, column } => write!(
                f,
                "Sample '{sample}' is missing a value in metadata column '{column}'."
            ),
            Self::MetadataValueNonPositive {
                sample,
                column,
                value,
            } => write!(
                f,
                "Sample '{sample}' has a non-positive value ({value}) in '{column}' — values must be > 0."
            ),
            Self::SampleNotInMapping { sample, column } => write!(
                f,
                "Sample '{sample}' is not present in the loaded metadata — cannot apply metadata \
                 normalization with column '{column}'. Re-load the metadata CSV so it covers every \
                 sample in the MS-DIAL .txt, or pick a different normalization method."
            ),
            Self::EmptyReferenceGroup { group } => write!(
                f,
                "Reference group '{group}' has no samples — cannot compute PQN reference spectrum."
            ),
            Self::ReferenceAllNan { method } => write!(
                f,
                "Reference distribution is all-NaN ({method}) — check that samples have measured values."
            ),
            Self::InsufficientSamplesAfterDrop {
                group,
                remaining,
                required,
                column,
            } => write!(
                f,
                "Dropping samples without a value in metadata column '{column}' leaves group '{group}' with only {remaining} sample(s); DAM requires at least {required}. Add values to the metadata CSV or pick a different normalization method."
            ),
            Self::PqnDegenerateSamples { samples } => {
                // Truncate the displayed list at 5 names for readability;
                // the full list is still in the struct field.
                let preview: String = if samples.len() > 5 {
                    let head: Vec<String> = samples.iter().take(5).cloned().collect();
                    format!("{} (and {} more)", head.join(", "), samples.len() - 5)
                } else {
                    samples.join(", ")
                };
                // After `filter-unassigned-samples-from-stage2`, Stage 2
                // callers narrow inputs to assigned samples only — so any
                // sample reaching this error path IS in the metadata CSV.
                // The pre-2026-05 hint "Drop them from the metadata CSV"
                // is therefore misleading (dropping a CSV row maps the
                // sample to UNASSIGNED, which is auto-filtered, so the
                // user would change the metadata, re-run, and hit the
                // same error). New guidance points at the two
                // actionable next steps.
                write!(
                    f,
                    "PQN: {} sample(s) have a degenerate quotient median (NaN or 0): {preview}. These samples have zero or NaN per-feature quotients against the reference. Pick a different normalization method, or check the affected sample columns in the MS-DIAL .txt.",
                    samples.len()
                )
            }
        }
    }
}

impl std::error::Error for NormalizationError {}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sample_not_in_mapping_display_names_sample_and_column() {
        let err = NormalizationError::SampleNotInMapping {
            sample: "Z9".into(),
            column: "dry_weight".into(),
        };
        let msg = format!("{err}");
        assert!(msg.contains("Z9"), "expected sample name 'Z9' in: {msg}");
        assert!(
            msg.contains("dry_weight"),
            "expected column name 'dry_weight' in: {msg}"
        );
        assert!(
            msg.contains("different normalization method"),
            "expected actionable hint in: {msg}"
        );
    }

    #[test]
    fn pqn_degenerate_display_carries_actionable_guidance() {
        let err = NormalizationError::PqnDegenerateSamples {
            samples: vec!["A1".into(), "A2".into()],
        };
        let msg = format!("{err}");
        assert!(msg.contains("A1"), "expected sample name 'A1' in: {msg}");
        assert!(msg.contains("A2"), "expected sample name 'A2' in: {msg}");
        assert!(
            msg.contains("different normalization method"),
            "expected actionable hint 'different normalization method' in: {msg}"
        );
        assert!(
            msg.contains("MS-DIAL .txt"),
            "expected actionable hint 'MS-DIAL .txt' in: {msg}"
        );
        // The old, misleading hint "Drop them from the metadata CSV" must
        // be gone (that path doesn't work post-Stage-2-boundary filter).
        assert!(
            !msg.contains("metadata CSV"),
            "old misleading hint 'metadata CSV' must be removed; got: {msg}"
        );
    }
}
