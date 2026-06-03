use std::path::PathBuf;

use metabolopan::data::load_group_mapping;

fn fixture(name: &str) -> PathBuf {
    PathBuf::from("tests/fixtures").join(name)
}

fn write_tmp_csv(content: &str) -> tempfile::NamedTempFile {
    use std::io::Write;
    let mut f = tempfile::Builder::new()
        .suffix(".csv")
        .tempfile()
        .expect("tempfile");
    f.write_all(content.as_bytes()).expect("write tempfile");
    f
}

#[test]
fn parses_fixture_with_extra_row_ignored_and_unassigned_detected() {
    let sample_cols: Vec<String> = ["S01", "S02", "S03", "S04", "S99_missing"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let mapping =
        load_group_mapping(&fixture("groups_example.csv"), &sample_cols).expect("must parse");

    let groups = mapping.groups();
    assert!(groups.contains(&"ASAP".to_string()));
    assert!(groups.contains(&"CK".to_string()));
    assert!(groups.contains(&"Unassigned".to_string()));
    assert!(
        !groups.contains(&"Other".to_string()),
        "extra CSV row must be ignored"
    );

    assert_eq!(mapping.samples_in("ASAP"), vec![0, 1]);
    assert_eq!(mapping.samples_in("CK"), vec![2, 3]);
    assert_eq!(mapping.group_of("S99_missing"), "Unassigned");
    assert_eq!(mapping.assigned_count(), 4);
}

#[test]
fn rejects_empty_group() {
    let f = write_tmp_csv("sample,group\nS01,ASAP\nS05,\n");
    let err = load_group_mapping(f.path(), &["S01".to_string(), "S05".to_string()])
        .expect_err("must error on empty group");
    let msg = format!("{err}");
    assert!(msg.contains("S05"), "error must mention S05: {msg}");
    assert!(
        msg.to_lowercase().contains("empty"),
        "error must mention 'empty': {msg}"
    );
}

#[test]
fn rejects_duplicate_sample() {
    let f = write_tmp_csv("sample,group\nS01,ASAP\nS01,CK\n");
    let err = load_group_mapping(f.path(), &["S01".to_string()])
        .expect_err("must error on duplicate sample");
    let msg = format!("{err}");
    assert!(msg.contains("S01"), "error must mention S01: {msg}");
    assert!(
        msg.to_lowercase().contains("duplicate"),
        "error must mention 'duplicate': {msg}"
    );
}

#[test]
fn rejects_header_missing_sample_column() {
    // `sample_name` is not `sample`: name-based detection finds no `sample`
    // column, so this fails as a missing-required-column error (the prior
    // positional parser rejected it as an unrecognized header shape).
    let f = write_tmp_csv("sample_name,group\nS01,ASAP\n");
    let err = load_group_mapping(f.path(), &["S01".to_string()])
        .expect_err("must error on missing sample column");
    let msg = format!("{err}");
    assert!(msg.contains("sample"), "error must name the `sample` column: {msg}");
    assert!(
        msg.to_lowercase().contains("missing"),
        "error must say the column is missing: {msg}"
    );
    assert!(
        msg.contains("sample_name,group"),
        "error should echo the actual header: {msg}"
    );
}

#[test]
fn no_overlap_returns_all_unassigned() {
    let f = write_tmp_csv("sample,group\nX01,A\nX02,B\n");
    let mapping = load_group_mapping(f.path(), &["S01".to_string(), "S02".to_string()])
        .expect("must succeed");

    assert_eq!(mapping.groups(), vec!["Unassigned"]);
    assert_eq!(mapping.assigned_count(), 0);
}

#[test]
fn two_column_csv_has_no_metadata_columns() {
    let f = write_tmp_csv("sample,group\nS01,A\nS02,B\n");
    let mapping =
        load_group_mapping(f.path(), &["S01".to_string(), "S02".to_string()]).expect("parses");
    assert!(mapping.metadata_column_names().is_empty());
    assert!(mapping.metadata_values("dry_weight").is_none());
}

#[test]
fn extra_numeric_columns_parse_and_align() {
    let f = write_tmp_csv(
        "sample,group,dry_weight,dilution\nS01,A,12.4,1.0\nS02,A,11.8,1.0\nS03,B,12.1,2.0\nS04,B,11.5,2.0\n",
    );
    let mapping = load_group_mapping(
        f.path(),
        &[
            "S01".to_string(),
            "S02".to_string(),
            "S03".to_string(),
            "S04".to_string(),
        ],
    )
    .expect("parses");
    assert_eq!(
        mapping.metadata_column_names(),
        vec!["dry_weight".to_string(), "dilution".to_string()],
        "metadata_column_names preserves CSV header order"
    );
    assert_eq!(
        mapping.metadata_values("dry_weight").expect("present"),
        &[Some(12.4), Some(11.8), Some(12.1), Some(11.5)]
    );
    assert_eq!(
        mapping.metadata_values("dilution").expect("present"),
        &[Some(1.0), Some(1.0), Some(2.0), Some(2.0)]
    );
}

#[test]
fn empty_metadata_cell_becomes_none() {
    let f = write_tmp_csv("sample,group,dry_weight\nS01,A,\nS02,A,11.8\n");
    let mapping = load_group_mapping(f.path(), &["S01".to_string(), "S02".to_string()])
        .expect("parses (empty cells become None)");
    assert_eq!(
        mapping.metadata_values("dry_weight").expect("present"),
        &[None, Some(11.8)]
    );
}

#[test]
fn partially_numeric_metadata_column_is_dropped_silently() {
    // Per D1 of relax-metadata-csv-validation: a single non-empty non-numeric
    // cell disqualifies the whole column from the metadata API.
    let f = write_tmp_csv("sample,group,dry_weight\nS01,A,12.4\nS02,A,N/A\n");
    let mapping = load_group_mapping(f.path(), &["S01".to_string(), "S02".to_string()])
        .expect("non-numeric metadata cell must NOT error");
    assert!(
        mapping.metadata_column_names().is_empty(),
        "dry_weight (with N/A) must be excluded from the numeric API surface"
    );
    assert!(
        mapping.metadata_values("dry_weight").is_none(),
        "dropped column must return None from metadata_values"
    );
}

#[test]
fn label_column_alongside_numeric_yields_only_numeric_in_dropdown() {
    // Reproduces the user-reported bug: a `biosample` label column alongside
    // a real numeric `dry_weight` column must not block the CSV load; only
    // the numeric column appears in the API.
    let f = write_tmp_csv(
        "sample,group,biosample,dry_weight\n\
         S01,A,CTR-01,12.4\n\
         S02,A,CTR-02,11.8\n\
         S03,B,TRT-01,12.1\n\
         S04,B,TRT-02,11.5\n",
    );
    let mapping = load_group_mapping(
        f.path(),
        &[
            "S01".to_string(),
            "S02".to_string(),
            "S03".to_string(),
            "S04".to_string(),
        ],
    )
    .expect("CSV with mixed label + numeric metadata must parse");
    assert_eq!(
        mapping.metadata_column_names(),
        vec!["dry_weight".to_string()],
        "only the numeric column appears in the dropdown source"
    );
    assert!(mapping.metadata_values("biosample").is_none());
    assert_eq!(
        mapping.metadata_values("dry_weight").expect("present"),
        &[Some(12.4), Some(11.8), Some(12.1), Some(11.5)]
    );
}

#[test]
fn mixed_numeric_columns_preserve_header_order() {
    // Two numeric columns flanking a non-numeric one — the surviving numeric
    // columns appear in CSV header order, with the label column gone.
    let f = write_tmp_csv(
        "sample,group,dilution,batch,dry_weight\n\
         S01,A,1.0,B1,12.4\n\
         S02,A,1.0,B2,11.8\n",
    );
    let mapping =
        load_group_mapping(f.path(), &["S01".to_string(), "S02".to_string()]).expect("parses");
    assert_eq!(
        mapping.metadata_column_names(),
        vec!["dilution".to_string(), "dry_weight".to_string()],
        "non-numeric `batch` removed; surviving columns keep header order"
    );
}

#[test]
fn sample_missing_from_csv_gets_none_in_every_metadata_column() {
    // CSV has S01 and S02; sample_cols includes S99 which has no CSV row.
    let f = write_tmp_csv("sample,group,dry_weight,dilution\nS01,A,12.4,1.0\nS02,A,11.8,2.0\n");
    let mapping = load_group_mapping(
        f.path(),
        &["S01".to_string(), "S02".to_string(), "S99".to_string()],
    )
    .expect("parses");
    let dry = mapping.metadata_values("dry_weight").expect("present");
    let dil = mapping.metadata_values("dilution").expect("present");
    assert_eq!(dry.len(), 3);
    assert_eq!(dry[2], None, "S99 missing from CSV -> None in dry_weight");
    assert_eq!(dil[2], None, "S99 missing from CSV -> None in dilution");
    assert_eq!(mapping.group_of("S99"), "Unassigned");
}

#[test]
fn metadata_column_order_is_csv_header_order_not_alphabetical() {
    let f = write_tmp_csv("sample,group,total_protein,dry_weight,dilution\nS01,A,1.0,2.0,3.0\n");
    let mapping = load_group_mapping(f.path(), &["S01".to_string()]).expect("parses");
    assert_eq!(
        mapping.metadata_column_names(),
        vec![
            "total_protein".to_string(),
            "dry_weight".to_string(),
            "dilution".to_string(),
        ],
    );
}
