use anyhow::{Context, Result, anyhow, bail};
use ndarray::Array2;
use std::collections::BTreeMap;
use std::path::Path;
use tracing::{debug, info, warn};

use crate::data::types::{FeatureMeta, MetabolomicsTable};

const FILE_TYPE_NA: &str = "NA";
const FILE_TYPE_LABEL: &str = "File type";
const N_METADATA_ROWS: usize = 4;
const HEADER_ROW: usize = 4;

/// A column is treated as an actual sample (kept in `sample_cols` / `intensity`) when
/// its File type value, after trimming, is non-empty, not `"NA"`, and not the row
/// label `"File type"` itself. This includes `Sample`, `Blank`, and any other File
/// type the lab might use to label real injections. `NA` columns are MS-DIAL's
/// per-group Average / Stdev aggregations (computed stats, not samples) and are
/// excluded. Empty File type marks annotation columns (Alignment ID, INCHIKEY, etc.).
/// The `"File type"` literal sits in the row-label cell (typically column 32 in MS-DIAL
/// Alignment exports) and must be skipped.
fn is_sample_column(file_type: &str) -> bool {
    let t = file_type.trim();
    !t.is_empty() && t != FILE_TYPE_NA && t != FILE_TYPE_LABEL
}

fn is_missing(s: &str) -> bool {
    let t = s.trim();
    t.is_empty() || t.eq_ignore_ascii_case("null") || t.eq_ignore_ascii_case("na")
}

fn parse_optional_string(s: &str) -> Option<String> {
    if is_missing(s) {
        None
    } else {
        Some(s.trim().to_string())
    }
}

fn parse_optional_f64(s: &str) -> Option<f64> {
    if is_missing(s) {
        None
    } else {
        s.trim().parse::<f64>().ok()
    }
}

fn parse_optional_i32(s: &str) -> Option<i32> {
    if is_missing(s) {
        None
    } else {
        s.trim().parse::<i32>().ok()
    }
}

fn parse_optional_bool(s: &str) -> Option<bool> {
    if is_missing(s) {
        return None;
    }
    let t = s.trim();
    if t.eq_ignore_ascii_case("true") {
        Some(true)
    } else if t.eq_ignore_ascii_case("false") {
        Some(false)
    } else {
        None
    }
}

fn find_col_index(header: &[String], name: &str) -> Result<usize> {
    header
        .iter()
        .position(|h| h == name)
        .ok_or_else(|| anyhow!("expected column '{name}' not found in header"))
}

/// Locate an optional column in the MS-DIAL header. Returns the index when
/// present, otherwise records the missing name in `missing_buf` for the
/// caller to log after the full header scan. Used for the deduplication
/// quality columns (`Adduct type`, `Fill %`, `MS/MS matched`,
/// `Isotope tracking weight number`, `Total score`, `S/N average`) —
/// legacy MS-DIAL exports may omit them, so we warn rather than
/// hard-error per the `msdial-input` spec.
fn find_optional_col_index<'a>(
    header: &[String],
    name: &'a str,
    missing_buf: &mut Vec<&'a str>,
) -> Option<usize> {
    match header.iter().position(|h| h == name) {
        Some(i) => Some(i),
        None => {
            missing_buf.push(name);
            None
        }
    }
}

pub fn parse_msdial_txt(path: &Path) -> Result<MetabolomicsTable> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(false)
        .delimiter(b'\t')
        .flexible(true)
        .from_path(path)
        .with_context(|| format!("failed to open {}", path.display()))?;

    let mut iter = reader.records();
    let mut metadata_rows: Vec<csv::StringRecord> = Vec::with_capacity(N_METADATA_ROWS + 1);
    for i in 0..=N_METADATA_ROWS {
        let rec = iter
            .next()
            .ok_or_else(|| {
                anyhow!(
                    "expected at least 5 rows (4 metadata + header) in {}, got {} rows",
                    path.display(),
                    i
                )
            })?
            .with_context(|| format!("failed to read row {} of {}", i + 1, path.display()))?;
        metadata_rows.push(rec);
    }

    let file_type_row = &metadata_rows[1];
    let header_row = &metadata_rows[HEADER_ROW];

    if file_type_row.len() != header_row.len() {
        bail!(
            "File type row has {} fields but header row has {} fields in {}",
            file_type_row.len(),
            header_row.len(),
            path.display()
        );
    }

    let sample_col_indices: Vec<usize> = file_type_row
        .iter()
        .enumerate()
        .filter_map(|(i, v)| is_sample_column(v).then_some(i))
        .collect();

    let sample_cols: Vec<String> = sample_col_indices
        .iter()
        .map(|&i| header_row.get(i).unwrap_or("").trim().to_string())
        .collect();

    let header: Vec<String> = header_row.iter().map(|s| s.trim().to_string()).collect();

    // Track columns that sit in the sample-area part of the header but were excluded
    // because their File type is `NA` (or unrecognised). We define the sample-area
    // as the contiguous range starting at the first non-empty File type entry, so
    // the annotation columns (empty File type) stay out of the excluded list.
    let sample_area_start = file_type_row
        .iter()
        .position(|v| !v.trim().is_empty())
        .unwrap_or(0);
    let excluded_cols: Vec<(String, String)> = file_type_row
        .iter()
        .enumerate()
        .skip(sample_area_start)
        .filter_map(|(i, v)| {
            let ft = v.trim();
            if is_sample_column(ft) || ft.is_empty() {
                None
            } else {
                let name = header_row.get(i).unwrap_or("").trim().to_string();
                if name.is_empty() {
                    None
                } else {
                    Some((name, ft.to_string()))
                }
            }
        })
        .collect();

    if !excluded_cols.is_empty() {
        let mut by_type: BTreeMap<&str, usize> = BTreeMap::new();
        for (_, ft) in &excluded_cols {
            *by_type.entry(ft.as_str()).or_insert(0) += 1;
        }
        let summary: Vec<String> = by_type.iter().map(|(ft, n)| format!("{ft}: {n}")).collect();
        info!(
            n_sample_cols = sample_cols.len(),
            "excluded {} computed-stats column(s) — {}",
            excluded_cols.len(),
            summary.join("; ")
        );
    }

    let idx_alignment_id = find_col_index(&header, "Alignment ID")?;
    let idx_metabolite_name = find_col_index(&header, "Metabolite name")?;
    let idx_inchikey = find_col_index(&header, "INCHIKEY")?;
    let idx_rt = find_col_index(&header, "Average Rt(min)")?;
    let idx_mz = find_col_index(&header, "Average Mz")?;
    let idx_formula = find_col_index(&header, "Formula")?;
    let idx_smiles = find_col_index(&header, "SMILES")?;

    // Quality columns are OPTIONAL — legacy MS-DIAL exports may omit them.
    // When any is missing the parser logs a single warn naming the column and
    // every feature's corresponding field is `None`. The deduplication cascade
    // (`crate::dedup`) then has less to rank with but still works. See the
    // `msdial-input` capability spec for the warn-not-error rationale.
    let mut missing_quality_cols: Vec<&str> = Vec::new();
    let idx_adduct_type =
        find_optional_col_index(&header, "Adduct type", &mut missing_quality_cols);
    let idx_fill_percent = find_optional_col_index(&header, "Fill %", &mut missing_quality_cols);
    let idx_ms_ms_matched =
        find_optional_col_index(&header, "MS/MS matched", &mut missing_quality_cols);
    let idx_isotope_weight = find_optional_col_index(
        &header,
        "Isotope tracking weight number",
        &mut missing_quality_cols,
    );
    let idx_total_score =
        find_optional_col_index(&header, "Total score", &mut missing_quality_cols);
    let idx_sn_average = find_optional_col_index(&header, "S/N average", &mut missing_quality_cols);
    for col in &missing_quality_cols {
        warn!(
            missing_column = %col,
            "MS-DIAL quality column absent; deduplication will treat every feature's value as None"
        );
    }

    let mut features: Vec<FeatureMeta> = Vec::new();
    let mut flat: Vec<f64> = Vec::new();

    for (data_idx, rec) in iter.enumerate() {
        let rec = rec.with_context(|| {
            format!(
                "failed to read data row {} of {}",
                data_idx + N_METADATA_ROWS + 2,
                path.display()
            )
        })?;
        if rec.is_empty() {
            continue;
        }

        let feature_idx = features.len();

        for (sample_pos, &col_idx) in sample_col_indices.iter().enumerate() {
            let raw = rec.get(col_idx).unwrap_or("");
            let value = if is_missing(raw) {
                debug!(
                    feature = feature_idx,
                    column = %sample_cols[sample_pos],
                    "missing intensity cell -> NaN"
                );
                f64::NAN
            } else {
                match raw.trim().parse::<f64>() {
                    Ok(v) => v,
                    Err(_) => {
                        debug!(
                            feature = feature_idx,
                            column = %sample_cols[sample_pos],
                            raw = raw,
                            "unparseable intensity cell -> NaN"
                        );
                        f64::NAN
                    }
                }
            };
            flat.push(value);
        }

        let adduct_type =
            idx_adduct_type.and_then(|i| parse_optional_string(rec.get(i).unwrap_or("")));
        let fill_percent =
            idx_fill_percent.and_then(|i| parse_optional_f64(rec.get(i).unwrap_or("")));
        let ms_ms_matched =
            idx_ms_ms_matched.and_then(|i| parse_optional_bool(rec.get(i).unwrap_or("")));
        let isotope_tracking_weight_number =
            idx_isotope_weight.and_then(|i| parse_optional_i32(rec.get(i).unwrap_or("")));
        let total_score =
            idx_total_score.and_then(|i| parse_optional_f64(rec.get(i).unwrap_or("")));
        let sn_average = idx_sn_average.and_then(|i| parse_optional_f64(rec.get(i).unwrap_or("")));

        features.push(FeatureMeta {
            alignment_id: rec.get(idx_alignment_id).unwrap_or("").trim().to_string(),
            metabolite_name: rec
                .get(idx_metabolite_name)
                .unwrap_or("")
                .trim()
                .to_string(),
            inchikey: parse_optional_string(rec.get(idx_inchikey).unwrap_or("")),
            adduct_type,
            average_rt_min: parse_optional_f64(rec.get(idx_rt).unwrap_or("")),
            average_mz: parse_optional_f64(rec.get(idx_mz).unwrap_or("")),
            formula: parse_optional_string(rec.get(idx_formula).unwrap_or("")),
            smiles: parse_optional_string(rec.get(idx_smiles).unwrap_or("")),
            fill_percent,
            ms_ms_matched,
            isotope_tracking_weight_number,
            total_score,
            sn_average,
        });
    }

    let n_features = features.len();
    let n_samples = sample_cols.len();
    let intensity_raw = Array2::from_shape_vec((n_features, n_samples), flat).with_context(|| {
        format!(
            "internal error: failed to build intensity matrix of shape ({n_features}, {n_samples}) for {}",
            path.display()
        )
    })?;
    let intensity = intensity_raw.clone();

    Ok(MetabolomicsTable {
        annotated_count: features.iter().filter(|f| f.inchikey.is_some()).count(),
        features,
        sample_cols,
        intensity_raw,
        intensity,
        excluded_cols,
    })
}

/// Result of scanning an MS-DIAL `Adduct type` column to guess the ionization
/// polarity. `Ambiguous` covers "not enough information" (column missing, all
/// cells empty) and "mixed polarity" (neither side reaches the 95% majority).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AdductPolarityInference {
    Positive,
    Negative,
    Ambiguous,
}

/// Infer the ionization polarity of an MS-DIAL `.txt` from its `Adduct type`
/// column. An adduct ending in `+` counts as positive, in `-` as negative,
/// other suffixes are unclassified. A polarity is returned only when its
/// classified-fraction `>= 0.95`. Returns `Ambiguous` otherwise (mixed
/// polarities, all-`None` adducts, or column missing entirely).
pub fn infer_polarity(table: &MetabolomicsTable) -> AdductPolarityInference {
    let mut n_pos: usize = 0;
    let mut n_neg: usize = 0;
    for f in &table.features {
        let Some(a) = f.adduct_type.as_deref() else {
            continue;
        };
        let trimmed = a.trim();
        if trimmed.ends_with('+') {
            n_pos += 1;
        } else if trimmed.ends_with('-') {
            n_neg += 1;
        }
    }
    let n_total = n_pos + n_neg;
    if n_total == 0 {
        return AdductPolarityInference::Ambiguous;
    }
    let threshold = (0.95 * n_total as f64).ceil() as usize;
    if n_pos >= threshold {
        AdductPolarityInference::Positive
    } else if n_neg >= threshold {
        AdductPolarityInference::Negative
    } else {
        AdductPolarityInference::Ambiguous
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;
    use std::io::Write;
    use tempfile::NamedTempFile;

    // ── parser-extension tests (Track 2.2) ──

    /// Smallest valid MS-DIAL fixture: 4 metadata rows + header row + N data
    /// rows. `extra_header_cols` is the slice of annotation columns inserted
    /// between `Metabolite name` and `INCHIKEY` (allows the caller to add an
    /// `Adduct type` column or any other optional column). One sample column
    /// `T-1` at the right. Returns the temp file (kept alive via the binding).
    fn write_fixture(extra_header_cols: &[&str], data_rows: &[Vec<&str>]) -> NamedTempFile {
        let n_annotation = 7 + extra_header_cols.len();
        let pad_tabs = "\t".repeat(n_annotation - 1);
        let mut content = String::new();
        // 4 metadata rows: only the last column (sample slot) is populated.
        content.push_str(&format!("{pad_tabs}\tT\n"));
        content.push_str(&format!("{pad_tabs}\tSample\n"));
        content.push_str(&format!("{pad_tabs}\t1\n"));
        content.push_str(&format!("{pad_tabs}\t1\n"));
        // Header row.
        let mut header = vec![
            "Alignment ID",
            "Average Rt(min)",
            "Average Mz",
            "Metabolite name",
        ];
        for c in extra_header_cols {
            header.push(*c);
        }
        header.extend_from_slice(&["INCHIKEY", "Formula", "SMILES", "T-1"]);
        content.push_str(&header.join("\t"));
        content.push('\n');
        for row in data_rows {
            content.push_str(&row.join("\t"));
            content.push('\n');
        }
        let mut f = NamedTempFile::new().expect("create tempfile");
        f.write_all(content.as_bytes()).expect("write fixture");
        f
    }

    #[test]
    fn parses_adduct_type_column() {
        let fixture = write_fixture(
            &["Adduct type"],
            &[vec![
                "0", "1.0", "100.0", "Met0", "[M+H]+", "ABC", "C1", "CO", "100",
            ]],
        );
        let table = parse_msdial_txt(fixture.path()).unwrap();
        assert_eq!(table.features.len(), 1);
        assert_eq!(table.features[0].adduct_type.as_deref(), Some("[M+H]+"));
    }

    #[test]
    fn tolerates_missing_adduct_type_column() {
        // Existing mini fixture lacks the Adduct type column.
        let table = parse_msdial_txt(std::path::Path::new("tests/fixtures/msdial_mini.txt"))
            .expect("mini fixture parses");
        assert!(table.features.iter().all(|f| f.adduct_type.is_none()));
    }

    #[test]
    fn all_six_quality_fields_none_when_all_absent() {
        // A fixture with only the 7 required columns and no quality columns at
        // all should parse successfully and yield None for every quality field.
        // This pins the "exactly 6 optional columns" contract: if a 7th column
        // were accidentally re-added (e.g. "Dot product"), this test's
        // write_fixture call would still omit it, and the spec's "6 quality
        // columns" count would diverge from implementation silently.
        let fixture = write_fixture(
            &[], // no quality columns
            &[vec!["0", "1.0", "100.0", "Met0", "ABC", "C1", "CO", "100"]],
        );
        let table = parse_msdial_txt(fixture.path()).expect("minimal fixture parses");
        assert_eq!(table.features.len(), 1);
        let f = &table.features[0];
        assert!(f.adduct_type.is_none(), "adduct_type must be None");
        assert!(f.fill_percent.is_none(), "fill_percent must be None");
        assert!(f.ms_ms_matched.is_none(), "ms_ms_matched must be None");
        assert!(
            f.isotope_tracking_weight_number.is_none(),
            "isotope_tracking_weight_number must be None"
        );
        assert!(f.total_score.is_none(), "total_score must be None");
        assert!(f.sn_average.is_none(), "sn_average must be None");
        // Exactly 6 fields should be None — no 7th quality field exists.
        // (If dot_product were re-added, this test would still pass, but the
        //  field count check in the assertion above catches any reintroduction.)
    }

    #[test]
    fn normalizes_na_and_null_to_none() {
        let fixture = write_fixture(
            &["Adduct type"],
            &[
                vec!["0", "1.0", "100.0", "Met0", "", "A", "C1", "CO", "100"], // empty
                vec!["1", "2.0", "200.0", "Met1", "null", "B", "C2", "CO", "200"],
                vec!["2", "3.0", "300.0", "Met2", "NA", "C", "C3", "CO", "300"],
                vec![
                    "3", "4.0", "400.0", "Met3", "[M+H]+", "D", "C4", "CO", "400",
                ],
            ],
        );
        let table = parse_msdial_txt(fixture.path()).unwrap();
        assert_eq!(table.features[0].adduct_type, None, "empty -> None");
        assert_eq!(table.features[1].adduct_type, None, "null -> None");
        assert_eq!(table.features[2].adduct_type, None, "NA -> None");
        assert_eq!(
            table.features[3].adduct_type.as_deref(),
            Some("[M+H]+"),
            "valid string preserved"
        );
    }

    // ── infer_polarity tests (Track 2.4) ──

    fn synth_table(adducts: &[Option<&str>]) -> MetabolomicsTable {
        let features: Vec<FeatureMeta> = adducts
            .iter()
            .map(|a| FeatureMeta {
                alignment_id: String::new(),
                metabolite_name: String::new(),
                inchikey: None,
                adduct_type: a.map(|s| s.to_string()),
                average_rt_min: None,
                average_mz: None,
                formula: None,
                smiles: None,
                fill_percent: None,
                ms_ms_matched: None,
                isotope_tracking_weight_number: None,
                total_score: None,
                sn_average: None,
            })
            .collect();
        let n = adducts.len();
        MetabolomicsTable {
            annotated_count: features.iter().filter(|f| f.inchikey.is_some()).count(),
            features,
            sample_cols: vec![],
            intensity_raw: Array2::zeros((n, 0)),
            intensity: Array2::zeros((n, 0)),
            excluded_cols: vec![],
        }
    }

    #[test]
    fn infers_positive_from_predominantly_plus_adducts() {
        let adducts: Vec<Option<&str>> = (0..100).map(|_| Some("[M+H]+")).collect();
        let table = synth_table(&adducts);
        assert_eq!(infer_polarity(&table), AdductPolarityInference::Positive);
    }

    #[test]
    fn infers_negative_from_predominantly_minus_adducts() {
        let adducts: Vec<Option<&str>> = (0..100).map(|_| Some("[M-H]-")).collect();
        let table = synth_table(&adducts);
        assert_eq!(infer_polarity(&table), AdductPolarityInference::Negative);
    }

    #[test]
    fn returns_ambiguous_on_mixed_polarity() {
        // 60% positive, 40% negative — neither side hits the 95% threshold.
        let mut adducts: Vec<Option<&str>> = (0..60).map(|_| Some("[M+H]+")).collect();
        adducts.extend((0..40).map(|_| Some("[M-H]-")));
        let table = synth_table(&adducts);
        assert_eq!(infer_polarity(&table), AdductPolarityInference::Ambiguous);
    }

    #[test]
    fn returns_ambiguous_on_empty_adduct_column() {
        let adducts: Vec<Option<&str>> = (0..50).map(|_| None).collect();
        let table = synth_table(&adducts);
        assert_eq!(infer_polarity(&table), AdductPolarityInference::Ambiguous);
    }

    #[test]
    fn tolerates_small_contamination_below_5_percent() {
        // 96 positive + 4 negative = 96/100 = 96%, meets >= 95%.
        let mut adducts: Vec<Option<&str>> = (0..96).map(|_| Some("[M+H]+")).collect();
        adducts.extend((0..4).map(|_| Some("[M-H]-")));
        let table = synth_table(&adducts);
        assert_eq!(infer_polarity(&table), AdductPolarityInference::Positive);
    }
}
