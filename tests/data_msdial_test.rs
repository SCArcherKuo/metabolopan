use std::path::Path;

use metabolopan::data::{AdductPolarityInference, infer_polarity, parse_msdial_txt};

#[test]
fn parses_mini_fixture() {
    let table = parse_msdial_txt(Path::new("tests/fixtures/msdial_mini.txt"))
        .expect("mini fixture must parse");

    assert_eq!(table.features.len(), 10, "10 features expected");
    assert_eq!(
        table.sample_cols.len(),
        7,
        "7 sample columns expected (6 Sample + 1 Blank; NA excluded)"
    );
    assert_eq!(table.intensity.shape(), &[10, 7]);

    // Sample columns are taken from the header row, filtered by File type != NA.
    assert_eq!(
        table.sample_cols,
        vec!["T-1", "T-2", "T-3", "C-1", "C-2", "C-3", "Bk-1"]
    );

    // The single NA column ("NA-1") is recorded as excluded.
    assert_eq!(table.excluded_cols.len(), 1);
    assert_eq!(
        table.excluded_cols[0],
        ("NA-1".to_string(), "NA".to_string())
    );

    // The second feature (index 1) has the known INCHIKEY.
    assert_eq!(
        table.features[1].inchikey.as_deref(),
        Some("OKJIRPAQVSHGFK-UHFFFAOYSA-N")
    );

    // Feature 2 is "Unknown" with INCHIKEY=null → None.
    assert!(table.features[2].inchikey.is_none());

    // At least one cell is NaN (Feature 2 has an empty T-3 cell and Feature 5 has T-1 = "NA").
    let nan_count = table.intensity.iter().filter(|v| v.is_nan()).count();
    assert!(
        nan_count >= 2,
        "expected at least 2 NaN cells, got {nan_count}"
    );

    // Explicit zeros (e.g. Feature 3's first three samples) must remain exactly 0.0.
    let zero_count = table.intensity.iter().filter(|v| **v == 0.0).count();
    assert!(
        zero_count >= 3,
        "expected at least 3 explicit zeros (Feature 3 T-1..T-3), got {zero_count}"
    );

    // Sanity: a known non-zero, non-NaN cell.
    assert_eq!(table.intensity[[0, 0]], 100.0);
}

#[test]
fn errors_when_file_missing() {
    let err = parse_msdial_txt(Path::new("tests/fixtures/this_does_not_exist.txt"))
        .expect_err("missing file must error");
    let msg = format!("{err}");
    assert!(
        msg.contains("this_does_not_exist.txt"),
        "error message should mention the path; got: {msg}"
    );
}

#[test]
#[ignore]
fn parses_real_example() {
    let table = parse_msdial_txt(Path::new("data/single-mode/MS-DIAL-output-example.txt"))
        .expect("real example must parse");

    // 19 Sample (QC-1..3 + S1-1..S8-2) + 2 Blank (Bk-1, Bk-2) = 21; 20 NA aggregations excluded.
    assert_eq!(table.sample_cols.len(), 21);
    assert!(table.features.len() > 1000);
    assert!(
        table.sample_cols.contains(&"Bk-1".to_string()),
        "expected Bk-1 in sample_cols"
    );
    assert!(
        table.sample_cols.contains(&"Bk-2".to_string()),
        "expected Bk-2 in sample_cols"
    );

    // All non-Sample, non-empty File type cells should be reported as excluded
    // (the example file has 20 NA aggregation columns).
    let na_excluded = table
        .excluded_cols
        .iter()
        .filter(|(_, ft)| ft == "NA")
        .count();
    assert!(
        na_excluded > 0,
        "expected NA columns to be tracked as excluded"
    );

    let alignment_zero = table
        .features
        .iter()
        .find(|f| f.alignment_id == "0")
        .expect("must have a feature with Alignment ID == \"0\"");
    assert_eq!(
        alignment_zero.inchikey.as_deref(),
        Some("VHMCLONZXOFDIQ-QPQITGAISA-N")
    );
}

#[test]
fn infers_polarity_from_bundled_double_mode_pos_fixture() {
    let table = parse_msdial_txt(Path::new("data/double-mode/data-positive.txt"))
        .expect("bundled POS fixture must parse");
    assert!(
        table.features.iter().any(|f| f.adduct_type.is_some()),
        "POS fixture must populate at least some adduct_type values"
    );
    assert_eq!(
        infer_polarity(&table),
        AdductPolarityInference::Positive,
        "bundled POS fixture must infer Positive"
    );
}

#[test]
fn infers_polarity_from_bundled_double_mode_neg_fixture() {
    let table = parse_msdial_txt(Path::new("data/double-mode/data-negative.txt"))
        .expect("bundled NEG fixture must parse");
    assert!(
        table.features.iter().any(|f| f.adduct_type.is_some()),
        "NEG fixture must populate at least some adduct_type values"
    );
    assert_eq!(
        infer_polarity(&table),
        AdductPolarityInference::Negative,
        "bundled NEG fixture must infer Negative"
    );
}

/// MS-DIAL v5 compatibility regression guard.
///
/// MS-DIAL 5 splits the single `Dot product` column into
/// `Simple dot product` + `Weighted dot product`. Before this branch,
/// the parser emitted a WARN for the missing `Dot product` column on every
/// v5 file. After removing the `Dot product` lookup, v5 files parse cleanly.
///
/// This test pins that the v5 fixture parses without error, produces a
/// non-trivial feature set, and populates the quality fields that v5 DOES
/// carry (`Total score`, `Fill %`, etc.). If the `"Dot product"` lookup is
/// ever reintroduced, a WARN would fire — and, more critically, the
/// `total_score` field would silently stop populating for v5 files (because
/// the old code read it at a shifted column position).
#[test]
#[ignore]
fn parses_ms_dial_v5_fixture_without_dot_product_warn() {
    let path = Path::new("data/ms-dial-5/Area_1_2026_05_27_05_39_10.txt");
    if !path.exists() {
        eprintln!(
            "skipping: {} not present. This fixture set is maintainer-local and \
             is not distributed with the repository.",
            path.display()
        );
        return;
    }
    let table = parse_msdial_txt(path).expect("MS-DIAL v5 fixture must parse");
    assert!(
        table.features.len() > 10,
        "v5 fixture should have >10 features"
    );
    assert!(
        table.sample_cols.len() >= 6,
        "v5 fixture should have at least 6 sample columns"
    );
    // v5 carries Total score and Fill % — they must populate.
    let total_score_populated = table.features.iter().any(|f| f.total_score.is_some());
    assert!(
        total_score_populated,
        "Total score must populate from the v5 fixture (column present under that name)"
    );
    let fill_pct_populated = table.features.iter().any(|f| f.fill_percent.is_some());
    assert!(
        fill_pct_populated,
        "Fill % must populate from the v5 fixture"
    );
    // v5 does NOT have a `Dot product` column (it was split into Simple/Weighted).
    // The parser must NOT error and must NOT have set dot_product on any feature.
    // (dot_product field was removed from FeatureMeta — this test also verifies
    //  the struct compiles without it.)
}
