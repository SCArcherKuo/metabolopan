use std::io::Write;

use metabolopan::data::{IonMode, load_group_mapping};
use metabolopan::ui::stage1_input::{Stage1ValidationInput, validate_for_dam};

fn write_tmp_csv(content: &str) -> tempfile::NamedTempFile {
    let mut f = tempfile::Builder::new()
        .suffix(".csv")
        .tempfile()
        .expect("tempfile");
    f.write_all(content.as_bytes()).expect("write tempfile");
    f
}

fn sample_cols(names: &[&str]) -> Vec<String> {
    names.iter().map(|s| s.to_string()).collect()
}

fn load(csv: &str, samples: &[&str]) -> metabolopan::data::GroupMapping {
    let f = write_tmp_csv(csv);
    load_group_mapping(f.path(), &sample_cols(samples)).expect("load")
}

// After `reorder-gui-and-move-mode-to-stage3`, Stage 1 is mode-agnostic — the
// KEGG species / Group / fetch-loaded gates live on Stage 3 setup. This helper
// builds the slimmed validation input used by the universal + dual-mode
// integrity tests below.
fn input<'a>(
    mapping: Option<&'a metabolopan::data::GroupMapping>,
    slot1_sample_cols: &'a [String],
) -> Stage1ValidationInput<'a> {
    Stage1ValidationInput {
        route: metabolopan::app::AnalysisRoute::DamEnrichment,
        table_loaded: true,
        slot1_sample_cols,
        slot2_sample_cols: None,
        mapping,
        slot1_mode: Some(IonMode::Positive),
        slot2_revealed: false,
        slot2_mode: None,
    }
}

#[test]
fn passes_with_two_groups_two_samples_each() {
    let s1 = sample_cols(&["S1", "S2", "S3", "S4"]);
    let mapping = load(
        "sample,group\nS1,A\nS2,A\nS3,B\nS4,B\n",
        &["S1", "S2", "S3", "S4"],
    );
    assert!(validate_for_dam(input(Some(&mapping), &s1)).is_ok());
}

#[test]
fn fails_all_unassigned() {
    let s1 = sample_cols(&["S1", "S2"]);
    let mapping = load("sample,group\nX1,A\nX2,B\n", &["S1", "S2"]);
    let issues = validate_for_dam(input(Some(&mapping), &s1)).expect_err("must fail");
    assert!(
        issues
            .iter()
            .any(|s| s.contains("No samples in the metadata")),
        "expected 'No samples in the metadata' issue; got: {issues:?}"
    );
}

#[test]
fn fails_single_group() {
    let s1 = sample_cols(&["S1", "S2", "S3"]);
    let mapping = load("sample,group\nS1,A\nS2,A\nS3,A\n", &["S1", "S2", "S3"]);
    let issues = validate_for_dam(input(Some(&mapping), &s1)).expect_err("must fail");
    assert!(
        issues.iter().any(|s| s.contains("At least 2 groups")),
        "expected '< 2 groups' issue; got: {issues:?}"
    );
}

#[test]
fn fails_single_sample_in_group() {
    let s1 = sample_cols(&["S1", "S2", "S3", "S4"]);
    let mapping = load(
        "sample,group\nS1,A\nS2,B\nS3,B\nS4,B\n",
        &["S1", "S2", "S3", "S4"],
    );
    let issues = validate_for_dam(input(Some(&mapping), &s1)).expect_err("must fail");
    assert!(
        issues
            .iter()
            .any(|s| s.contains("`A`") && s.contains("1 sample")),
        "expected single-sample issue mentioning group A; got: {issues:?}"
    );
}

#[test]
fn lists_multiple_issues_together() {
    let s1 = sample_cols(&["S1", "S2"]);
    let mapping = load("sample,group\nS1,A\n", &["S1", "S2"]);
    let issues = validate_for_dam(input(Some(&mapping), &s1)).expect_err("must fail");
    assert!(
        issues.iter().any(|s| s.contains("At least 2 groups")),
        "expected '< 2 groups' issue; got: {issues:?}"
    );
    assert!(
        issues
            .iter()
            .any(|s| s.contains("`A`") && s.contains("1 sample")),
        "expected '< 2 samples' issue mentioning A; got: {issues:?}"
    );
}
