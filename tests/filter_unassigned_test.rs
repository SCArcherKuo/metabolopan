//! Stage 1 → Stage 2 boundary regression tests for
//! `filter-unassigned-samples-from-stage2`.
//!
//! Each test exercises the public `run_dam` entry point with inputs that
//! have been narrowed via the `without_unassigned_samples` helpers — the
//! same code path `start_dam` follows in the UI. We do not invoke
//! `start_dam` directly (it's private and bound to `App` state); the
//! helper-level tests live in `src/ui/stage2_setup.rs`.

use ndarray::Array2;
use std::io::Write;

use metabolopan::dam::{DamConfig, DamMethod, fdr::FdrMethod, run_dam};
use metabolopan::data::{
    FeatureMeta, GroupMapping, IonMode, IonModeTable, MetabolomicsTable, load_group_mapping,
    parse_msdial_txt,
};
use metabolopan::normalize::{NormalizationConfig, NormalizationMethod, PqnReference};

fn write_csv(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new()
        .suffix(".csv")
        .tempfile()
        .expect("tempfile");
    f.write_all(content.as_bytes()).expect("write fixture");
    f
}

fn cols(names: &[&str]) -> Vec<String> {
    names.iter().map(|s| s.to_string()).collect()
}

/// Base `DamConfig` for these boundary tests: Welch, no normalization,
/// drop-unknown on, dedup off, log-transform on, BH FDR. Each call overrides
/// `normalization` via struct-update (the field these tests actually vary).
fn base_cfg() -> DamConfig {
    DamConfig {
        method: DamMethod::Welch,
        normalization: NormalizationConfig::default(),
        drop_unknown: true,
        dedup_enabled: false,
        log_transform: true,
        fdr_method: FdrMethod::BenjaminiHochberg,
    }
}

fn empty_feature(id: &str) -> FeatureMeta {
    FeatureMeta {
        alignment_id: id.into(),
        metabolite_name: format!("met-{id}"),
        inchikey: Some(format!("KEY-{id}")),
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

/// Make a small synthetic table with `n_features` rows × the given sample
/// columns. Values are deterministic per (feature, sample) so the test
/// can compare two runs cell-by-cell.
fn make_synthetic_table(sample_cols_v: &[&str], n_features: usize) -> MetabolomicsTable {
    let n_samples = sample_cols_v.len();
    let mut intensity = Array2::<f64>::zeros((n_features, n_samples));
    // Distinct, finite, non-zero values: per-sample multiplier × per-feature
    // base so a t-test sees real signal.
    for i in 0..n_features {
        for j in 0..n_samples {
            // Multiplier 1.0 for "A*" samples, 3.0 for "B*" samples, 99.0
            // for "Blank*". Feature base = (i + 1) * 10.0.
            let mult = if sample_cols_v[j].starts_with('A') {
                1.0
            } else if sample_cols_v[j].starts_with('B') {
                3.0
            } else {
                99.0
            };
            intensity[[i, j]] = mult * ((i + 1) as f64) * 10.0;
        }
    }
    let features: Vec<FeatureMeta> = (0..n_features)
        .map(|i| empty_feature(&format!("F{i}")))
        .collect();
    MetabolomicsTable {
        annotated_count: features.iter().filter(|f| f.inchikey.is_some()).count(),
        features,
        sample_cols: cols(sample_cols_v),
        intensity_raw: intensity.clone(),
        intensity,
        excluded_cols: vec![],
    }
}

// =====================================================================
// T16 — Bug-report regression
// =====================================================================
//
// Reproduces `data/bug-report-2026-05-27_132310/`: the NEG fixture has 21
// sample columns, 2 of which (SolventBlank_NEG_01, SolventBlank_NEG_02)
// are NOT in `metadata.csv` and therefore land in the Unassigned bucket.
// Pre-fix: PQN reference = Group("QC") aborts with PqnDegenerateSamples
// because the per-sample factor loop iterates the blank columns.
// Post-fix: `without_unassigned_samples` drops them before `run_dam`
// sees the matrix — the run completes cleanly.

#[tokio::test]
async fn bug_report_regression_pqn_qc_reference_runs_without_degenerate_error() {
    let neg_path = std::path::Path::new("data/double-mode/data-negative.txt");
    if !neg_path.exists() {
        eprintln!("skipping: {} not present", neg_path.display());
        return;
    }
    let table = parse_msdial_txt(neg_path).expect("parse NEG fixture");
    let mapping_full = load_group_mapping(
        std::path::Path::new("data/double-mode/metadata.csv"),
        &table.sample_cols,
    )
    .expect("load mapping");

    // Pre-condition: this fixture must have at least one Unassigned sample
    // for the regression to be meaningful. (SolventBlank_NEG_01/_02.)
    assert!(
        !mapping_full
            .samples_in(metabolopan::data::UNASSIGNED)
            .is_empty(),
        "fixture invariant: NEG must carry at least one Unassigned sample"
    );

    // Apply the Stage 1 → Stage 2 boundary helpers (what start_dam does).
    let mapping = mapping_full.without_unassigned_samples();
    let mut table = table.without_unassigned_samples(&mapping_full);
    assert!(
        !table
            .sample_cols
            .contains(&"SolventBlank_NEG_01".to_string()),
        "boundary helper should have dropped SolventBlank_NEG_01"
    );
    assert!(
        !table
            .sample_cols
            .contains(&"SolventBlank_NEG_02".to_string()),
        "boundary helper should have dropped SolventBlank_NEG_02"
    );

    let result = run_dam(
        &mut table,
        &mapping,
        "Treatment",
        "Control",
        &DamConfig {
            normalization: NormalizationConfig {
                method: NormalizationMethod::Pqn {
                    reference: PqnReference::Group("QC".to_string()),
                },
            },
            ..base_cfg()
        },
        None,
    )
    .await;

    // The whole point: must NOT surface PqnDegenerateSamples.
    let res = result.expect("PQN with QC reference must succeed after boundary filter");
    assert!(
        !res.features.is_empty(),
        "expected at least one DAM feature after PQN ref=QC"
    );
}

// =====================================================================
// T17 — Bit-equal: auto-filter vs hand-pruned produce identical results
// =====================================================================
//
// Construct a small fixture in two parallel forms:
//   (a) auto-filter: full sample axis [A1, A2, Blank, B1, B2] + CSV that
//                    omits Blank → boundary helper drops Blank
//   (b) hand-pruned: sample axis [A1, A2, B1, B2] only + CSV for same
// Assert run_dam over (a) and (b) produce bit-equal DamFeature vectors.

#[tokio::test]
async fn bit_equal_auto_filter_matches_hand_pruned_input() {
    // Common metadata CSV: 4 assigned samples in two groups.
    let csv = "sample,group\nA1,A\nA2,A\nB1,B\nB2,B\n";

    // (a) Auto-filter path: full axis includes Blank, CSV omits it.
    let f_a = write_csv(csv);
    let mapping_a_full = load_group_mapping(f_a.path(), &cols(&["A1", "A2", "Blank", "B1", "B2"]))
        .expect("mapping a");
    let table_a = make_synthetic_table(&["A1", "A2", "Blank", "B1", "B2"], 8);
    let mapping_a = mapping_a_full.without_unassigned_samples();
    let mut table_a = table_a.without_unassigned_samples(&mapping_a_full);

    // (b) Hand-pruned path: axis already excludes Blank.
    let f_b = write_csv(csv);
    let mapping_b =
        load_group_mapping(f_b.path(), &cols(&["A1", "A2", "B1", "B2"])).expect("mapping b");
    let mut table_b = make_synthetic_table(&["A1", "A2", "B1", "B2"], 8);

    // Sanity: the two paths should now be observationally equivalent.
    assert_eq!(table_a.sample_cols, table_b.sample_cols);
    assert_eq!(table_a.intensity_raw.shape(), table_b.intensity_raw.shape());

    let norm = NormalizationConfig {
        method: NormalizationMethod::Sum,
    };
    let res_a = run_dam(
        &mut table_a,
        &mapping_a,
        "A",
        "B",
        &DamConfig {
            normalization: norm.clone(),
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("run_dam (auto-filter path) must succeed");
    let res_b = run_dam(
        &mut table_b,
        &mapping_b,
        "A",
        "B",
        &DamConfig {
            normalization: norm.clone(),
            ..base_cfg()
        },
        None,
    )
    .await
    .expect("run_dam (hand-pruned path) must succeed");

    assert_eq!(
        res_a.features.len(),
        res_b.features.len(),
        "feature counts must match"
    );
    for (i, (fa, fb)) in res_a.features.iter().zip(res_b.features.iter()).enumerate() {
        assert_eq!(fa.alignment_id, fb.alignment_id, "feature {i} id");
        // Means / medians / fold change should match exactly.
        assert_eq!(fa.numerator_mean, fb.numerator_mean, "feature {i} num_mean");
        assert_eq!(
            fa.denominator_mean, fb.denominator_mean,
            "feature {i} den_mean"
        );
        assert_eq!(
            fa.log2_fold_change, fb.log2_fold_change,
            "feature {i} log2FC"
        );
        // p-values: bit-equal via floating-point reordering equivalence.
        // assert_eq! works for finite values; NaN compares unequal so use a tolerant check.
        if fa.p_value.is_nan() {
            assert!(fb.p_value.is_nan(), "feature {i} NaN p_value mismatch");
        } else {
            assert_eq!(fa.p_value, fb.p_value, "feature {i} p_value");
        }
        if fa.p_adjusted.is_nan() {
            assert!(
                fb.p_adjusted.is_nan(),
                "feature {i} NaN p_adjusted mismatch"
            );
        } else {
            assert_eq!(fa.p_adjusted, fb.p_adjusted, "feature {i} p_adjusted");
        }
    }
}

// =====================================================================
// T18 — Dual-mode asymmetry
// =====================================================================
//
// POS and NEG carry DIFFERENT Unassigned sample-name sets. Each per-mode
// table should drop only its own per-mode Unassigned columns; the shared
// `assigned_mapping` should be identical across both calls.

#[test]
fn dual_mode_asymmetric_unassigned_filters_per_mode_independently() {
    // CSV: 4 assigned samples (2 POS + 2 NEG), no SolventBlank rows.
    let csv = "sample,biosample,group\nA_POS_01,bio-A,g1\nA_POS_02,bio-A,g1\nA_NEG_01,bio-A,g2\nA_NEG_02,bio-A,g2\n";
    let f = write_csv(csv);
    // Union sample axis: 2 POS assigned + 1 POS unassigned + 2 NEG assigned + 2 NEG unassigned.
    let union_cols = cols(&[
        "A_POS_01",
        "A_POS_02",
        "Blank_POS_01",
        "A_NEG_01",
        "A_NEG_02",
        "Blank_NEG_01",
        "Blank_NEG_02",
    ]);
    let mapping_full = load_group_mapping(f.path(), &union_cols).expect("mapping");

    // Per-mode tables — each carries only its own polarity's samples.
    let pos_table = make_synthetic_table(&["A_POS_01", "A_POS_02", "Blank_POS_01"], 3);
    let neg_table =
        make_synthetic_table(&["A_NEG_01", "A_NEG_02", "Blank_NEG_01", "Blank_NEG_02"], 3);
    let pos_it = IonModeTable {
        mode: IonMode::Positive,
        table: pos_table,
        txt_path: None,
    };
    let neg_it = IonModeTable {
        mode: IonMode::Negative,
        table: neg_table,
        txt_path: None,
    };

    // Boundary helpers (same path start_dam follows).
    let mapping_assigned: GroupMapping = mapping_full.without_unassigned_samples();
    let pos_assigned = pos_it.without_unassigned_samples(&mapping_full);
    let neg_assigned = neg_it.without_unassigned_samples(&mapping_full);

    // POS dropped Blank_POS_01 only.
    assert_eq!(
        pos_assigned.table.sample_cols,
        cols(&["A_POS_01", "A_POS_02"])
    );
    // NEG dropped Blank_NEG_01 + Blank_NEG_02 only.
    assert_eq!(
        neg_assigned.table.sample_cols,
        cols(&["A_NEG_01", "A_NEG_02"])
    );
    // The shared mapping carries both modes' assigned samples (4 total).
    assert_eq!(mapping_assigned.assigned_count(), 4);
    assert!(
        !mapping_assigned
            .groups()
            .contains(&"Unassigned".to_string())
    );
    // Per-mode sample axes are independent — each filter only consulted
    // its own per-mode `sample_cols`, never the union positional axis.
    assert!(
        !pos_assigned
            .table
            .sample_cols
            .contains(&"Blank_NEG_01".to_string())
    );
    assert!(
        !neg_assigned
            .table
            .sample_cols
            .contains(&"Blank_POS_01".to_string())
    );
}
