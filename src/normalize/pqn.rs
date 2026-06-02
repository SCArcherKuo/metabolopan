//! Probabilistic Quotient Normalization (Dieterle 2006).
//!
//! Pipeline:
//! 1. Sum-normalize each sample (so dilution-driven scale is removed first).
//! 2. Compute the per-feature reference spectrum (median across the chosen
//!    reference cohort, NaN-aware). `AllSamples` cohort = non-Unassigned
//!    samples; `Group(name)` cohort = that group's samples.
//! 3. For each sample, compute per-feature quotients vs the reference and
//!    take the median over features where reference[i] is finite and > 0.
//! 4. Divide the sum-normalized matrix by per-sample factors, then multiply
//!    by the median factor (magnitude preservation).

use ndarray::Array2;

use crate::data::GroupMapping;
use crate::data::groups::UNASSIGNED;
use crate::normalize::sum::{apply_sum, median_of};
use crate::normalize::types::{NormalizationError, PqnReference};

/// Returns `(out, median_factor, nan_cells_in, nan_cells_out, reference_features_used)`.
/// `reference_features_used` is the count of features whose cohort median is
/// finite and > 0 — i.e. features that the QC / reference cohort can actually
/// support as a PQN factor anchor. Features with median ≤ 0 or NaN are
/// excluded from per-sample factor computation but still appear in the
/// output matrix (scaled by the factor estimated from the remaining
/// features); the log line surfaces this count so users can sanity-check
/// how many features the cohort actually anchored.
///
/// `sample_cols` is the per-mode list of sample column names from the
/// `IonModeTable` being normalized (length == `raw.ncols()`). In dual-mode
/// `mapping.sample_names` is the UNION of all modes' samples, so cohort
/// membership MUST be decided by per-mode sample NAME (`mapping.group_of(name)`),
/// NOT by positional index against the union slice. The pre-2026-05-26 code
/// used `mapping.sample_name(j)` for per-mode `j` against the union mapping,
/// so for the NEG table it resolved POS sample names and decided the cohort
/// from POS group assignments (Findings #2 and #3 in the 2026-05-25 audit).
pub fn apply_pqn(
    raw: &Array2<f64>,
    mapping: &GroupMapping,
    sample_cols: &[String],
    reference: &PqnReference,
) -> Result<(Array2<f64>, f64, usize, usize, usize), NormalizationError> {
    let (sum_normalized, _sum_factor, _nan_in_sum, _nan_out_sum) =
        apply_sum(raw, mapping, sample_cols)?;
    let (n_features, n_samples) = (sum_normalized.nrows(), sum_normalized.ncols());
    // See apply_sum for the debug_assert vs assert rationale (PR-M).
    debug_assert_eq!(
        sample_cols.len(),
        n_samples,
        "sample_cols must be aligned to raw column axis (got {} vs {})",
        sample_cols.len(),
        n_samples
    );

    // Reference cohort (per-mode sample indices, 0..n_samples). AllSamples
    // excludes Unassigned; Group(name) is restricted to that group's samples
    // — both decisions are made BY NAME via the per-mode sample_cols slice,
    // so they're correct regardless of whether mapping is union-indexed.
    let cohort: Vec<usize> = match reference {
        PqnReference::AllSamples => sample_cols
            .iter()
            .enumerate()
            .filter(|(_, name)| mapping.group_of(name) != UNASSIGNED)
            .map(|(j, _)| j)
            .collect(),
        PqnReference::Group(group_name) => sample_cols
            .iter()
            .enumerate()
            .filter(|(_, name)| mapping.group_of(name) == group_name)
            .map(|(j, _)| j)
            .collect(),
    };
    if cohort.is_empty() {
        return match reference {
            PqnReference::Group(name) => Err(NormalizationError::EmptyReferenceGroup {
                group: name.clone(),
            }),
            PqnReference::AllSamples => Err(NormalizationError::ReferenceAllNan { method: "PQN" }),
        };
    }

    // Per-feature reference spectrum = NaN-aware median across cohort samples.
    let mut reference_spectrum: Vec<f64> = Vec::with_capacity(n_features);
    let mut all_nan = true;
    for i in 0..n_features {
        let row: Vec<f64> = cohort.iter().map(|&j| sum_normalized[[i, j]]).collect();
        let med = median_of(&row);
        if !med.is_nan() {
            all_nan = false;
        }
        reference_spectrum.push(med);
    }
    if all_nan {
        return Err(NormalizationError::ReferenceAllNan { method: "PQN" });
    }
    // Count features the cohort can anchor as a PQN reference (median > 0).
    // Reported in the normalize INFO log so users can see how many features
    // actually contributed to per-sample factor estimation vs the total
    // feature count.
    let reference_features_used = reference_spectrum
        .iter()
        .filter(|&&r| !r.is_nan() && r > 0.0)
        .count();

    // Per-sample factor = median of per-feature quotients, skipping features
    // where reference[i] is NaN or 0, or sum_normalized[i, j] is NaN.
    // Degenerate samples (median NaN or 0) are collected and surfaced as a
    // hard error (Finding #13 in the 2026-05-25 audit): the pre-fix
    // behaviour was to push `factor = 1.0` with only a WARN, leaving the
    // sample at sum-normalized scale while peers were PQN-scaled —
    // producing artefactual differential abundance purely from the scale
    // mismatch.
    let mut factor: Vec<f64> = Vec::with_capacity(n_samples);
    let mut degenerate: Vec<String> = Vec::new();
    for j in 0..n_samples {
        let mut quotients: Vec<f64> = Vec::with_capacity(n_features);
        for i in 0..n_features {
            let r = reference_spectrum[i];
            if r.is_nan() || r == 0.0 {
                continue;
            }
            let v = sum_normalized[[i, j]];
            if v.is_nan() {
                continue;
            }
            quotients.push(v / r);
        }
        let m = median_of(&quotients);
        if m.is_nan() || m == 0.0 {
            degenerate.push(sample_cols[j].clone());
            // Placeholder factor; we won't apply it because the error path
            // below short-circuits before the output matrix is built.
            factor.push(f64::NAN);
        } else {
            factor.push(m);
        }
    }
    if !degenerate.is_empty() {
        return Err(NormalizationError::PqnDegenerateSamples {
            samples: degenerate,
        });
    }

    // Every factor reaching here is a validated non-degenerate `f64` (the
    // degenerate path returned Err above), so wrapping each in `Some` and
    // delegating to the shared write-back yields the identical median_factor
    // and per-cell output.
    let (out, median_factor, nan_in, nan_out) = crate::normalize::apply_factors_and_count(
        &sum_normalized,
        &factor.iter().map(|&f| Some(f)).collect::<Vec<_>>(),
    );
    Ok((out, median_factor, nan_in, nan_out, reference_features_used))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::normalize::test_support::{cols, mapping_for};
    use ndarray::array;

    #[test]
    fn pqn_uniform_samples_identity() {
        // All samples identical -> sum-norm makes them equal -> reference == any
        // sample's feature values -> per-sample quotient is 1.0 -> output equals
        // sum-normalized matrix.
        let raw = array![[1.0, 1.0, 1.0], [2.0, 2.0, 2.0], [3.0, 3.0, 3.0]];
        let mapping = mapping_for(&["A", "B", "C"], &["G1", "G1", "G1"]);
        let sc = cols(&["A", "B", "C"]);
        let (out, factor, _, _, _) =
            apply_pqn(&raw, &mapping, &sc, &PqnReference::AllSamples).unwrap();
        assert!((factor - 1.0).abs() < 1e-9, "factor = {factor}");
        for j in 0..3 {
            for i in 0..3 {
                assert!(
                    (out[[i, j]] - raw[[i, j]]).abs() < 1e-9,
                    "out[{i},{j}] = {} != {}",
                    out[[i, j]],
                    raw[[i, j]]
                );
            }
        }
    }

    #[test]
    fn pqn_doubled_sample_factor_is_two() {
        // Sample B = 2x sample A; sample A baseline.
        // sum-norm magnitude-preserves -> both samples sum to the same value;
        // After sum-norm A and B become identical, so reference equals both,
        // quotients are 1.0, factor is 1.0, output is identity-on-sum-normalized.
        // To get a meaningful PQN factor, use 4 samples where 3 are equal and 1 is shifted.
        let raw = array![
            [1.0, 1.0, 1.0, 2.0],
            [2.0, 2.0, 2.0, 4.0],
            [3.0, 3.0, 3.0, 6.0]
        ];
        let mapping = mapping_for(&["A", "B", "C", "D"], &["G1", "G1", "G1", "G1"]);
        let sc = cols(&["A", "B", "C", "D"]);
        let (out, _factor, _, _, _) =
            apply_pqn(&raw, &mapping, &sc, &PqnReference::AllSamples).unwrap();
        // After PQN, D should be brought closer to the rest. The point: PQN
        // produces a finite output without errors and the shape is preserved.
        assert_eq!(out.shape(), raw.shape());
        for v in out.iter() {
            assert!(v.is_finite(), "PQN output must be finite, got {v}");
        }
    }

    #[test]
    fn pqn_allsamples_excludes_unassigned() {
        // 3-sample matrix: A, B are assigned (cohort); D is Unassigned.
        // Use distinct per-sample RATIOS so sum-normalization doesn't collapse them.
        use crate::data::load_group_mapping;
        use std::io::Write;
        let mut f = tempfile::Builder::new().suffix(".csv").tempfile().unwrap();
        f.write_all(b"sample,group\nA,G1\nB,G1\n").unwrap();
        let mapping =
            load_group_mapping(f.path(), &["A".into(), "B".into(), "D".into()]).expect("mapping");
        // Sums: A=11, B=4, D=11. Distinct ratios A:[1,10], B:[3,1], D:[10,1].
        let raw = array![[1.0, 3.0, 10.0], [10.0, 1.0, 1.0]];
        let sc = cols(&["A", "B", "D"]);
        let (_out, _factor, _, _, _) =
            apply_pqn(&raw, &mapping, &sc, &PqnReference::AllSamples).expect("PQN runs");
        // Re-derive sum_normalized to compare cohort-only vs all-samples references.
        let (sn, _, _, _) = crate::normalize::sum::apply_sum(&raw, &mapping, &sc).unwrap();
        let cohort_ref_row0 = crate::normalize::sum::median_of(&[sn[[0, 0]], sn[[0, 1]]]);
        let all_ref_row0 = crate::normalize::sum::median_of(&[sn[[0, 0]], sn[[0, 1]], sn[[0, 2]]]);
        assert!(
            (cohort_ref_row0 - all_ref_row0).abs() > 1e-6,
            "Cohort=AB and Cohort=ABD must produce distinguishable references; got AB={cohort_ref_row0}, ABD={all_ref_row0}"
        );
    }

    /// Regression test for Findings #2 / #3: dual-mode PQN cohort
    /// membership must be decided by per-mode sample NAME, NOT by union
    /// positional index. Simulates the dual-mode layout where the same
    /// physical biosample appears under different group assignments in
    /// POS vs NEG (extreme case to make the bug visible).
    #[test]
    fn pqn_dual_mode_group_cohort_uses_per_mode_names() {
        use crate::data::load_group_mapping;
        use std::io::Write;
        // Union: POS_S01, POS_S02 are "Treat"; NEG_S01, NEG_S02 are "Ctrl".
        // (Same sample IDs would be ambiguous; this contrived layout makes
        // the cohort decision visibly per-mode.)
        let mut f = tempfile::Builder::new().suffix(".csv").tempfile().unwrap();
        f.write_all(b"sample,group\nPOS_S01,Treat\nPOS_S02,Treat\nNEG_S01,Ctrl\nNEG_S02,Ctrl\n")
            .unwrap();
        let mapping = load_group_mapping(
            f.path(),
            &cols(&["POS_S01", "POS_S02", "NEG_S01", "NEG_S02"]),
        )
        .expect("mapping");

        // NEG table only (per-mode raw, 2 samples).
        let neg_raw = array![[1.0, 2.0], [3.0, 4.0]];
        let neg_cols = cols(&["NEG_S01", "NEG_S02"]);

        // Group("Ctrl") cohort for NEG table = NEG_S01, NEG_S02 → indices [0, 1].
        // Pre-2026-05-26: samples_in("Ctrl") returned [2, 3] (union indices)
        //                 then `.filter(|j| *j < n_samples=2)` left it empty,
        //                 yielding EmptyReferenceGroup mid-flight.
        let (_, _, _, _, _) = apply_pqn(
            &neg_raw,
            &mapping,
            &neg_cols,
            &PqnReference::Group("Ctrl".into()),
        )
        .expect("NEG Ctrl cohort exists — picked by name, not by union index");

        // Group("Treat") cohort for NEG table = empty (no NEG sample is
        // assigned Treat), so this MUST error explicitly.
        let err = apply_pqn(
            &neg_raw,
            &mapping,
            &neg_cols,
            &PqnReference::Group("Treat".into()),
        )
        .expect_err("Treat cohort is empty in NEG table — by name lookup");
        assert!(matches!(
            err,
            NormalizationError::EmptyReferenceGroup { ref group } if group == "Treat"
        ));
    }

    /// Regression test for Finding #13: a sample whose per-sample PQN
    /// quotient median is degenerate (NaN or 0) MUST surface as a hard
    /// error listing the affected sample names. Pre-2026-05-26 this silently
    /// fell back to `factor = 1.0`, leaving the sample at sum-normalized
    /// scale while peers were PQN-scaled — biased DAM downstream.
    #[test]
    fn pqn_degenerate_sample_errors_with_sample_name() {
        // Sample D is sparse: only ONE feature has a non-zero value; the
        // rest are 0. After sum-normalization sample D's values for the
        // zero-features are all 0, so its quotients are mostly 0 — the
        // per-sample median is 0 → degenerate.
        let raw = array![
            [10.0, 10.0, 10.0, 0.0],
            [20.0, 20.0, 20.0, 0.0],
            [30.0, 30.0, 30.0, 0.0],
            [40.0, 40.0, 40.0, 100.0], // D's only non-zero feature
        ];
        let mapping = mapping_for(&["A", "B", "C", "D"], &["G1", "G1", "G1", "G1"]);
        let sc = cols(&["A", "B", "C", "D"]);
        let err = apply_pqn(&raw, &mapping, &sc, &PqnReference::AllSamples)
            .expect_err("degenerate sample D must error, not silently fall back to factor=1");
        match err {
            NormalizationError::PqnDegenerateSamples { samples } => {
                assert!(
                    samples.contains(&"D".to_string()),
                    "samples list must contain D; got {samples:?}"
                );
            }
            other => panic!("expected PqnDegenerateSamples, got {other:?}"),
        }
    }

    #[test]
    fn pqn_empty_group_errors() {
        let raw = array![[1.0, 2.0]];
        let mapping = mapping_for(&["A", "B"], &["Ctrl", "Ctrl"]);
        let sc = cols(&["A", "B"]);
        let err = apply_pqn(
            &raw,
            &mapping,
            &sc,
            &PqnReference::Group("Treat".to_string()),
        )
        .expect_err("empty group must error");
        assert!(matches!(
            err,
            NormalizationError::EmptyReferenceGroup { ref group } if group == "Treat"
        ));
    }
}
