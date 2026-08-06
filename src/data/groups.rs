use anyhow::{Context, Result, bail};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use tracing::warn;

pub const UNASSIGNED: &str = "Unassigned";

#[derive(Debug, Clone)]
pub struct GroupMapping {
    sample_names: Vec<String>,
    sample_col_to_group: Vec<String>,
    group_to_indices: HashMap<String, Vec<usize>>,
    /// Optional numeric metadata columns parsed from CSV columns other than
    /// `sample`, `group`, and `biosample`. Key = column name; value =
    /// per-sample `Option<f64>` aligned to `sample_names`.
    metadata_columns: HashMap<String, Vec<Option<f64>>>,
    /// CSV-header order of metadata column names (not alphabetical).
    metadata_order: Vec<String>,
    /// Optional biosample column, present when the CSV has a column named
    /// `biosample` (recognized by name in any position). Per-sample
    /// values aligned to `sample_names` (`None` when the CSV row lacks a
    /// value, or when the column is absent entirely). The outer `Option`
    /// distinguishes "column not present" from "column present but values
    /// missing for some samples".
    biosample: Option<Vec<Option<String>>>,
    /// CSV-only samples: the `sample` value of every CSV row that matched no
    /// entry in `sample_cols` (the rows ignored-with-a-WARN at load). Retained
    /// (not just logged) so Stage 1 can surface them in a banner. Sorted
    /// ascending, deduplicated (CSV `sample` values are unique — duplicate rows
    /// are rejected at parse), and disjoint from the `Unassigned` set (which is
    /// the opposite direction: `sample_cols` entries with no CSV row). Owner: the `group-mapping` capability.
    unmatched_csv_samples: Vec<String>,
}

impl GroupMapping {
    pub fn groups(&self) -> Vec<String> {
        let mut out: Vec<String> = self.group_to_indices.keys().cloned().collect();
        out.sort();
        out
    }

    pub fn samples_in(&self, group: &str) -> Vec<usize> {
        self.group_to_indices
            .get(group)
            .cloned()
            .unwrap_or_default()
    }

    pub fn group_of(&self, sample: &str) -> &str {
        for (i, name) in self.sample_names.iter().enumerate() {
            if name == sample {
                return self.sample_col_to_group[i].as_str();
            }
        }
        UNASSIGNED
    }

    pub fn groups_in_order(&self) -> Vec<(usize, String)> {
        self.sample_col_to_group
            .iter()
            .enumerate()
            .map(|(i, g)| (i, g.clone()))
            .collect()
    }

    pub fn assigned_count(&self) -> usize {
        self.sample_col_to_group
            .iter()
            .filter(|g| g.as_str() != UNASSIGNED)
            .count()
    }

    /// CSV-only samples: the `sample` value of every CSV row that named a
    /// sample absent from `sample_cols` (ignored-with-a-WARN at load). Sorted
    /// ascending, no repeats, and disjoint from the `Unassigned` set. Empty
    /// when every CSV row matched a sample column. Stage 1 lists these in a
    /// `theme::ERROR` banner (owner: the `stage1-ui` capability).
    pub fn unmatched_csv_samples(&self) -> &[String] {
        &self.unmatched_csv_samples
    }

    /// Names of the numeric metadata columns parsed from the CSV (every column
    /// other than `sample`, `group`, and `biosample`), in **CSV header order**
    /// (NOT alphabetically sorted). Empty when the CSV had only `sample,group`.
    pub fn metadata_column_names(&self) -> Vec<String> {
        self.metadata_order.clone()
    }

    /// Per-sample values for the named metadata column. The returned slice is
    /// aligned to the `sample_cols` passed to `load_group_mapping` (one entry
    /// per sample; `None` for samples missing from the CSV). Returns `None` when
    /// the column does not exist.
    pub fn metadata_values(&self, col: &str) -> Option<&[Option<f64>]> {
        self.metadata_columns.get(col).map(|v| v.as_slice())
    }

    /// Per-sample metadata value looked up by sample **NAME**, NOT positional
    /// index. Returns `Some(Some(v))` for "value present", `Some(None)` for
    /// "row present but value empty in CSV", `None` for "column missing OR
    /// sample not in the mapping". Designed for dual-mode where positional
    /// indices on a per-mode raw matrix do NOT correspond to the union-
    /// indexed `GroupMapping` — `apply_metadata` MUST use this accessor to
    /// avoid mixing per-mode and union indexing (Finding #1).
    pub fn metadata_value_of(&self, sample: &str, col: &str) -> Option<Option<f64>> {
        let values = self.metadata_columns.get(col)?;
        let idx = self.sample_names.iter().position(|s| s == sample)?;
        values.get(idx).copied()
    }

    /// Sample name at column index `j` in the MS-DIAL sample axis. Returns
    /// `None` when `j` is out of range. Used by normalization error messages
    /// to surface a meaningful sample identifier (e.g. "S03") rather than a
    /// bare column index.
    pub fn sample_name(&self, j: usize) -> Option<&str> {
        self.sample_names.get(j).map(|s| s.as_str())
    }

    /// Total number of samples on the MS-DIAL axis (== `sample_cols.len()`
    /// passed to `load_group_mapping`).
    pub fn sample_count(&self) -> usize {
        self.sample_names.len()
    }

    /// Biosample label for `sample` when the CSV provided a `biosample`
    /// column AND the row for `sample` had a non-empty value. Returns
    /// `None` when the CSV lacked the column, the cell was empty, or
    /// `sample` is absent from the CSV. Used by Stage 1's dual-mode
    /// cross-mode consistency check.
    /// Whether the CSV carried a `biosample` column. Used by Stage 1 dual-mode
    /// validation: a dual-mode run requires this to be `true` (per D4).
    pub fn has_biosample(&self) -> bool {
        self.biosample.is_some()
    }

    pub fn biosample_of(&self, sample: &str) -> Option<&str> {
        let bio = self.biosample.as_ref()?;
        let idx = self.sample_names.iter().position(|s| s == sample)?;
        bio.get(idx).and_then(|v| v.as_deref())
    }

    /// Count of DISTINCT non-null biosample labels among `group`'s samples, or
    /// `None` when the CSV carried no `biosample` column. Used by the Data tab's
    /// per-group `<group> (<n> samples, <b> biosamples)` line
    /// (`data-summary-panel`). Resolves the group's sample indices via
    /// `samples_in` → `sample_name` → `biosample_of` (the group is keyed by
    /// index, so a name-keyed accessor does not exist).
    pub fn biosample_count(&self, group: &str) -> Option<usize> {
        if !self.has_biosample() {
            return None;
        }
        let mut seen: std::collections::HashSet<&str> = std::collections::HashSet::new();
        for j in self.samples_in(group) {
            if let Some(name) = self.sample_name(j)
                && let Some(bio) = self.biosample_of(name)
            {
                seen.insert(bio);
            }
        }
        Some(seen.len())
    }

    /// Stage-1 → Stage-2 boundary helper. Returns a new owned mapping with
    /// every Unassigned sample removed from `sample_names`,
    /// `sample_col_to_group`, `group_to_indices` (the `UNASSIGNED` key
    /// disappears), every per-numeric-column `metadata_columns` slot, and
    /// `biosample`. Surviving groups are reindexed against the narrowed
    /// `sample_names` axis. The original mapping is not mutated. Owner: the `group-mapping` capability.
    pub fn without_unassigned_samples(&self) -> Self {
        let kept: Vec<usize> = self
            .sample_col_to_group
            .iter()
            .enumerate()
            .filter(|(_, g)| g.as_str() != UNASSIGNED)
            .map(|(i, _)| i)
            .collect();

        let sample_names: Vec<String> =
            kept.iter().map(|&i| self.sample_names[i].clone()).collect();
        let sample_col_to_group: Vec<String> = kept
            .iter()
            .map(|&i| self.sample_col_to_group[i].clone())
            .collect();

        let mut group_to_indices: HashMap<String, Vec<usize>> = HashMap::new();
        for (new_idx, g) in sample_col_to_group.iter().enumerate() {
            group_to_indices.entry(g.clone()).or_default().push(new_idx);
        }

        let metadata_columns: HashMap<String, Vec<Option<f64>>> = self
            .metadata_columns
            .iter()
            .map(|(name, values)| {
                // Bounds-checked: under the loader invariant every `i` in
                // `kept` is in range (so this equals `values[i]`); an
                // out-of-range slot degrades to `None` rather than panicking
                // (`convert-defensive-panics-to-errors`).
                let narrowed: Vec<Option<f64>> = kept
                    .iter()
                    .map(|&i| values.get(i).copied().flatten())
                    .collect();
                (name.clone(), narrowed)
            })
            .collect();

        let biosample: Option<Vec<Option<String>>> = self.biosample.as_ref().map(|values| {
            kept.iter()
                .map(|&i| values.get(i).cloned().flatten())
                .collect::<Vec<_>>()
        });

        Self {
            sample_names,
            sample_col_to_group,
            group_to_indices,
            metadata_columns,
            metadata_order: self.metadata_order.clone(),
            biosample,
            // Carried through unchanged: CSV-only samples are about the original
            // CSV load, not the (narrowed) sample axis. Stage 2+ never reads it.
            unmatched_csv_samples: self.unmatched_csv_samples.clone(),
        }
    }
}

pub fn load_group_mapping(csv_path: &Path, sample_cols: &[String]) -> Result<GroupMapping> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_path(csv_path)
        .with_context(|| format!("failed to open {}", csv_path.display()))?;

    let headers = reader
        .headers()
        .with_context(|| format!("failed to read header of {}", csv_path.display()))?;
    let header_vec: Vec<String> = headers.iter().map(|s| s.to_string()).collect();

    // Name-based column discovery (positions/order are unconstrained): the CSV
    // MUST contain columns named `sample` and `group` (exact, case-sensitive).
    // An optional `biosample` column is recognized by name in any position and
    // routed to a dedicated slot (never a metadata extra). Every other column
    // is an optional metadata extra. Owner: the `group-mapping` capability.
    let positions = |name: &str| -> Vec<usize> {
        header_vec
            .iter()
            .enumerate()
            .filter_map(|(i, h)| (h.as_str() == name).then_some(i))
            .collect::<Vec<_>>()
    };
    let header_str = header_vec.join(",");
    let sample_col_idx = match positions("sample").as_slice() {
        [idx] => *idx,
        [] => bail!(
            "metadata CSV is missing a required `sample` column; header was '{header_str}' in {}",
            csv_path.display()
        ),
        hits => bail!(
            "metadata CSV has {} `sample` columns; exactly one is required (header '{header_str}') in {}",
            hits.len(),
            csv_path.display()
        ),
    };
    let group_col_idx = match positions("group").as_slice() {
        [idx] => *idx,
        [] => bail!(
            "metadata CSV is missing a required `group` column; header was '{header_str}' in {}",
            csv_path.display()
        ),
        hits => bail!(
            "metadata CSV has {} `group` columns; exactly one is required (header '{header_str}') in {}",
            hits.len(),
            csv_path.display()
        ),
    };
    let biosample_col_idx: Option<usize> = match positions("biosample").as_slice() {
        [] => None,
        [idx] => Some(*idx),
        hits => bail!(
            "metadata CSV has {} `biosample` columns; at most one is allowed (header '{header_str}') in {}",
            hits.len(),
            csv_path.display()
        ),
    };

    // Metadata columns: every header position EXCEPT sample (0), group, and
    // biosample (when present). Preserve CSV header order via the index walk.
    let metadata_positions: Vec<(usize, String)> = header_vec
        .iter()
        .enumerate()
        .filter(|(i, _)| {
            *i != sample_col_idx && *i != group_col_idx && Some(*i) != biosample_col_idx
        })
        .map(|(i, h)| (i, h.clone()))
        .collect();
    let metadata_count = metadata_positions.len();

    // Accumulate per-sample state during the row scan. Per-metadata-column
    // numeric/non-numeric classification follows the existing "drop the
    // column if any non-empty cell fails numeric parse" rule.
    let mut from_csv: HashMap<String, String> = HashMap::new();
    let mut metadata_by_sample: HashMap<String, Vec<Option<f64>>> = HashMap::new();
    let mut biosample_by_sample: HashMap<String, Option<String>> = HashMap::new();
    let mut non_empty_counts: Vec<usize> = vec![0; metadata_count];
    let mut non_numeric_counts: Vec<usize> = vec![0; metadata_count];

    for (row_idx, rec) in reader.records().enumerate() {
        let row_no = row_idx + 2;
        let rec = rec
            .with_context(|| format!("failed to read row {} of {}", row_no, csv_path.display()))?;
        let min_fields = sample_col_idx.max(group_col_idx);
        if rec.len() <= min_fields {
            bail!(
                "row {} in {} has fewer than {} fields (need at least sample + group columns)",
                row_no,
                csv_path.display(),
                min_fields + 1
            );
        }
        let sample = rec.get(sample_col_idx).unwrap_or("").trim().to_string();
        let group = rec.get(group_col_idx).unwrap_or("").trim().to_string();
        if sample.is_empty() {
            bail!(
                "row {} in {} has empty sample name",
                row_no,
                csv_path.display()
            );
        }
        if group.is_empty() {
            bail!(
                "sample '{sample}' (row {}) has empty group in {}",
                row_no,
                csv_path.display()
            );
        }
        if from_csv.contains_key(&sample) {
            bail!(
                "duplicate sample '{sample}' (row {}) in {}",
                row_no,
                csv_path.display()
            );
        }

        let biosample_value = biosample_col_idx.and_then(|i| {
            let trimmed = rec.get(i).unwrap_or("").trim();
            if trimmed.is_empty() {
                None
            } else {
                Some(trimmed.to_string())
            }
        });

        let mut row_meta: Vec<Option<f64>> = Vec::with_capacity(metadata_count);
        for (col_offset, (header_idx, _)) in metadata_positions.iter().enumerate() {
            let raw = rec.get(*header_idx).unwrap_or("");
            let trimmed = raw.trim();
            if trimmed.is_empty() {
                row_meta.push(None);
            } else {
                non_empty_counts[col_offset] += 1;
                match trimmed.parse::<f64>() {
                    Ok(v) => row_meta.push(Some(v)),
                    Err(_) => {
                        // Column will be reclassified non-numeric below; per-row
                        // slot stays aligned via None so sibling numeric columns
                        // remain readable.
                        non_numeric_counts[col_offset] += 1;
                        row_meta.push(None);
                    }
                }
            }
        }

        from_csv.insert(sample.clone(), group);
        biosample_by_sample.insert(sample.clone(), biosample_value);
        metadata_by_sample.insert(sample, row_meta);
    }

    // Drop non-numeric metadata columns from the API surface with a WARN.
    let mut metadata_order: Vec<String> = Vec::with_capacity(metadata_count);
    for (col_offset, (_, col_name)) in metadata_positions.iter().enumerate() {
        if non_numeric_counts[col_offset] == 0 {
            metadata_order.push(col_name.clone());
        } else {
            warn!(
                "column '{col_name}' dropped from metadata normalization options: {failed} of {non_empty} cells failed numeric parse",
                failed = non_numeric_counts[col_offset],
                non_empty = non_empty_counts[col_offset],
            );
        }
    }

    let mut sample_col_to_group: Vec<String> = Vec::with_capacity(sample_cols.len());
    let mut group_to_indices: HashMap<String, Vec<usize>> = HashMap::new();
    let mut matched_csv_samples: HashSet<&str> = HashSet::new();

    for (i, name) in sample_cols.iter().enumerate() {
        let group = if let Some(g) = from_csv.get(name) {
            matched_csv_samples.insert(name.as_str());
            g.clone()
        } else {
            UNASSIGNED.to_string()
        };
        group_to_indices.entry(group.clone()).or_default().push(i);
        sample_col_to_group.push(group);
    }

    // Build per-numeric-column `Vec<Option<f64>>` aligned to `sample_cols`,
    // reading from `metadata_by_sample` via each surviving column's offset
    // within `metadata_positions`.
    let metadata_offset_by_name: HashMap<&str, usize> = metadata_positions
        .iter()
        .enumerate()
        .map(|(off, (_, name))| (name.as_str(), off))
        .collect();
    let mut metadata_columns: HashMap<String, Vec<Option<f64>>> = HashMap::new();
    for col_name in &metadata_order {
        let offset = metadata_offset_by_name[col_name.as_str()];
        let mut aligned: Vec<Option<f64>> = Vec::with_capacity(sample_cols.len());
        for name in sample_cols {
            let value = metadata_by_sample
                .get(name)
                .and_then(|row| row.get(offset).copied().flatten());
            aligned.push(value);
        }
        metadata_columns.insert(col_name.clone(), aligned);
    }

    // Build biosample aligned to `sample_cols` when the CSV had the column.
    let biosample = biosample_col_idx.map(|_| {
        sample_cols
            .iter()
            .map(|name| biosample_by_sample.get(name).cloned().unwrap_or(None))
            .collect::<Vec<Option<String>>>()
    });

    // CSV-only samples: rows that named a sample absent from `sample_cols`.
    // Ignored-with-a-WARN as before, but now ALSO retained on the mapping so
    // Stage 1 can list them in a banner. `from_csv` is keyed by unique sample
    // values (duplicate rows already rejected), so no repeats; sort for a
    // deterministic banner + reproducible tests.
    let mut unmatched_csv_samples: Vec<String> = Vec::new();
    for csv_sample in from_csv.keys() {
        if !matched_csv_samples.contains(csv_sample.as_str()) {
            warn!(
                "metadata sample '{csv_sample}' is not present in MS-DIAL sample columns; ignoring"
            );
            unmatched_csv_samples.push(csv_sample.clone());
        }
    }
    unmatched_csv_samples.sort();

    let mapping = GroupMapping {
        sample_names: sample_cols.to_vec(),
        sample_col_to_group,
        group_to_indices,
        metadata_columns,
        metadata_order,
        biosample,
        unmatched_csv_samples,
    };

    if mapping.assigned_count() == 0 {
        warn!(
            "no samples in {} match any MS-DIAL sample column; all samples will be Unassigned",
            csv_path.display()
        );
    }

    Ok(mapping)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_csv(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().expect("create tempfile");
        f.write_all(content.as_bytes()).expect("write fixture");
        f
    }

    fn sample_cols(names: &[&str]) -> Vec<String> {
        names.iter().map(|s| s.to_string()).collect()
    }

    #[test]
    fn parses_3_column_csv_with_biosample_at_position_1() {
        let f = write_csv(
            "sample,biosample,group\n\
             CTR_positive_01,CTR-01,control\n\
             CTR_positive_02,CTR-02,control\n\
             TM_positive_01,TM-01,treatment\n\
             TM_positive_02,TM-02,treatment\n",
        );
        let cols = sample_cols(&[
            "CTR_positive_01",
            "CTR_positive_02",
            "TM_positive_01",
            "TM_positive_02",
        ]);
        let m = load_group_mapping(f.path(), &cols).unwrap();
        assert_eq!(m.biosample_of("CTR_positive_01"), Some("CTR-01"));
        assert_eq!(m.biosample_of("TM_positive_02"), Some("TM-02"));
        assert_eq!(m.group_of("CTR_positive_01"), "control");
        assert_eq!(m.group_of("TM_positive_02"), "treatment");
        // metadata_column_names() is empty (no numeric extras beyond the
        // recognized sample/biosample/group trio).
        assert!(m.metadata_column_names().is_empty());
    }

    #[test]
    fn parses_2_column_csv_with_biosample_at_position_2() {
        let f = write_csv(
            "sample,group,biosample\n\
             S01,control,CTR-01\n\
             S02,treatment,TM-01\n",
        );
        let cols = sample_cols(&["S01", "S02"]);
        let m = load_group_mapping(f.path(), &cols).unwrap();
        assert_eq!(m.biosample_of("S01"), Some("CTR-01"));
        assert_eq!(m.biosample_of("S02"), Some("TM-01"));
        // biosample must NOT appear in metadata_column_names()
        assert!(m.metadata_column_names().is_empty());
        assert!(m.metadata_values("biosample").is_none());
    }

    #[test]
    fn biosample_of_returns_none_when_column_absent() {
        let f = write_csv("sample,group\nS01,control\nS02,treatment\n");
        let cols = sample_cols(&["S01", "S02"]);
        let m = load_group_mapping(f.path(), &cols).unwrap();
        assert_eq!(m.biosample_of("S01"), None);
        assert_eq!(m.biosample_of("S02"), None);
    }

    #[test]
    fn biosample_of_returns_none_for_empty_cell() {
        let f = write_csv(
            "sample,biosample,group\n\
             S01,,control\n\
             S02,TM-01,treatment\n",
        );
        let cols = sample_cols(&["S01", "S02"]);
        let m = load_group_mapping(f.path(), &cols).unwrap();
        assert_eq!(m.biosample_of("S01"), None);
        assert_eq!(m.biosample_of("S02"), Some("TM-01"));
    }

    #[test]
    fn biosample_of_returns_none_for_unknown_sample() {
        let f = write_csv("sample,biosample,group\nS01,A,control\n");
        let cols = sample_cols(&["S01", "S99"]);
        let m = load_group_mapping(f.path(), &cols).unwrap();
        assert_eq!(m.biosample_of("S99"), None);
        // S99 is also Unassigned (no CSV row); biosample slot is None at its index.
        assert_eq!(m.group_of("S99"), UNASSIGNED);
    }

    #[test]
    fn metadata_column_names_excludes_biosample_at_either_position() {
        // Position 1 (3-column form) plus numeric extra.
        let f = write_csv(
            "sample,biosample,group,dry_weight\n\
             S01,A,control,1.0\n\
             S02,B,control,2.0\n",
        );
        let m = load_group_mapping(f.path(), &sample_cols(&["S01", "S02"])).unwrap();
        assert_eq!(m.metadata_column_names(), vec!["dry_weight".to_string()]);
        assert!(m.metadata_values("biosample").is_none());

        // Position 2 (2-column form with biosample as a metadata-style extra).
        let f2 = write_csv(
            "sample,group,biosample,dry_weight\n\
             S01,control,A,1.0\n\
             S02,control,B,2.0\n",
        );
        let m2 = load_group_mapping(f2.path(), &sample_cols(&["S01", "S02"])).unwrap();
        assert_eq!(m2.metadata_column_names(), vec!["dry_weight".to_string()]);
        assert!(m2.metadata_values("biosample").is_none());
        assert_eq!(m2.biosample_of("S01"), Some("A"));
    }

    #[test]
    fn arbitrary_column_order_parses_by_name() {
        // group before sample (the prior positional rule rejected this).
        let f = write_csv("group,sample\nA,S01\nB,S02\n");
        let m = load_group_mapping(f.path(), &sample_cols(&["S01", "S02"])).unwrap();
        assert_eq!(m.group_of("S01"), "A");
        assert_eq!(m.group_of("S02"), "B");

        // A metadata column before the key columns.
        let f2 = write_csv("dry_weight,sample,group\n1.5,S01,A\n2.5,S02,A\n");
        let m2 = load_group_mapping(f2.path(), &sample_cols(&["S01", "S02"])).unwrap();
        assert_eq!(m2.group_of("S01"), "A");
        assert_eq!(m2.metadata_column_names(), vec!["dry_weight".to_string()]);
        assert_eq!(
            m2.metadata_values("dry_weight"),
            Some([Some(1.5), Some(2.5)].as_ref())
        );

        // group, biosample, sample — biosample recognized by name mid-header.
        let f3 = write_csv("group,biosample,sample\nA,BIO-1,S01\nA,BIO-2,S02\n");
        let m3 = load_group_mapping(f3.path(), &sample_cols(&["S01", "S02"])).unwrap();
        assert_eq!(m3.group_of("S01"), "A");
        assert_eq!(m3.biosample_of("S01"), Some("BIO-1"));
        assert_eq!(m3.biosample_of("S02"), Some("BIO-2"));
        assert!(m3.metadata_column_names().is_empty());
    }

    #[test]
    fn missing_sample_column_errors_naming_sample() {
        let f = write_csv("id,group\nS01,A\n");
        let err = load_group_mapping(f.path(), &sample_cols(&["S01"])).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("sample"), "msg: {msg}");
        assert!(
            msg.contains("id,group"),
            "header should be echoed; msg: {msg}"
        );
    }

    #[test]
    fn missing_group_column_errors_naming_group() {
        let f = write_csv("sample,treatment\nS01,A\n");
        let err = load_group_mapping(f.path(), &sample_cols(&["S01"])).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("group"), "msg: {msg}");
    }

    #[test]
    fn duplicate_key_column_is_rejected() {
        // duplicate `sample`
        let f = write_csv("sample,group,sample\nS01,A,X\n");
        let err = load_group_mapping(f.path(), &sample_cols(&["S01"])).unwrap_err();
        assert!(err.to_string().contains("sample"), "msg: {err}");

        // duplicate `group`
        let f2 = write_csv("sample,group,group\nS01,A,B\n");
        let err2 = load_group_mapping(f2.path(), &sample_cols(&["S01"])).unwrap_err();
        assert!(err2.to_string().contains("group"), "msg: {err2}");

        // duplicate `biosample`
        let f3 = write_csv("sample,biosample,group,biosample\nS01,X,A,Y\n");
        let err3 = load_group_mapping(f3.path(), &sample_cols(&["S01"])).unwrap_err();
        assert!(err3.to_string().contains("biosample"), "msg: {err3}");
    }

    #[test]
    fn without_unassigned_round_trip_on_zero_unassigned_mapping_is_clone() {
        // mapping with every sample assigned (no Unassigned) → narrowed == clone
        let f = write_csv(
            "sample,group,dry_weight\n\
             S01,A,1.0\n\
             S02,B,2.0\n",
        );
        let m_in = load_group_mapping(f.path(), &sample_cols(&["S01", "S02"])).unwrap();
        assert_eq!(m_in.assigned_count(), 2);

        let m_out = m_in.without_unassigned_samples();
        assert_eq!(m_out.groups(), m_in.groups());
        assert_eq!(m_out.assigned_count(), m_in.assigned_count());
        assert_eq!(m_out.sample_count(), m_in.sample_count());
        assert_eq!(m_out.samples_in("A"), m_in.samples_in("A"));
        assert_eq!(m_out.samples_in("B"), m_in.samples_in("B"));
        assert_eq!(m_out.group_of("S01"), m_in.group_of("S01"));
        assert_eq!(m_out.group_of("S02"), m_in.group_of("S02"));
        assert_eq!(
            m_out.metadata_values("dry_weight"),
            m_in.metadata_values("dry_weight")
        );
    }

    #[test]
    fn without_unassigned_drops_unassigned_samples_and_reindexes() {
        // 5 samples; C and E are Unassigned (missing from CSV).
        let f = write_csv(
            "sample,group\n\
             A,g1\n\
             B,g2\n\
             D,g1\n",
        );
        let m_in = load_group_mapping(f.path(), &sample_cols(&["A", "B", "C", "D", "E"])).unwrap();
        assert_eq!(m_in.assigned_count(), 3);
        assert_eq!(m_in.sample_count(), 5);
        assert_eq!(m_in.group_of("C"), UNASSIGNED);
        assert_eq!(m_in.group_of("E"), UNASSIGNED);

        let m_out = m_in.without_unassigned_samples();
        assert_eq!(m_out.assigned_count(), 3);
        assert_eq!(m_out.sample_count(), 3);
        // groups() returns alphabetical via HashMap-key sort; Unassigned gone.
        assert_eq!(m_out.groups(), vec!["g1".to_string(), "g2".to_string()]);
        // A → new idx 0, D → new idx 2 (B took 1); B → new idx 1.
        assert_eq!(m_out.samples_in("g1"), vec![0, 2]);
        assert_eq!(m_out.samples_in("g2"), vec![1]);
        // Group lookups for kept samples are stable.
        assert_eq!(m_out.group_of("A"), "g1");
        assert_eq!(m_out.group_of("B"), "g2");
        assert_eq!(m_out.group_of("D"), "g1");
        // Defensive sentinel for un-tracked samples.
        assert_eq!(m_out.group_of("C"), UNASSIGNED);
        assert_eq!(m_out.group_of("E"), UNASSIGNED);
    }

    #[test]
    fn without_unassigned_narrows_metadata_and_biosample_per_sample() {
        // 3 samples; C is Unassigned. dry_weight + biosample columns present.
        let f = write_csv(
            "sample,biosample,group,dry_weight\n\
             A,BIO-A,g1,10.0\n\
             B,BIO-B,g1,12.0\n",
        );
        let m_in = load_group_mapping(f.path(), &sample_cols(&["A", "B", "C"])).unwrap();
        assert_eq!(m_in.assigned_count(), 2);
        // Pre-filter slices include C's slot (Some(None)/None values).
        assert_eq!(
            m_in.metadata_values("dry_weight"),
            Some([Some(10.0), Some(12.0), None].as_ref())
        );

        let m_out = m_in.without_unassigned_samples();
        assert_eq!(
            m_out.metadata_column_names(),
            vec!["dry_weight".to_string()]
        );
        assert_eq!(
            m_out.metadata_values("dry_weight"),
            Some([Some(10.0), Some(12.0)].as_ref()) // C's None is dropped
        );
        assert_eq!(m_out.biosample_of("A"), Some("BIO-A"));
        assert_eq!(m_out.biosample_of("B"), Some("BIO-B"));
        // C is no longer tracked → biosample_of returns None (same semantics
        // as querying an unknown sample on m_in).
        assert_eq!(m_out.biosample_of("C"), None);
    }

    #[test]
    fn without_unassigned_does_not_mutate_source() {
        let f = write_csv("sample,group\nA,g1\nB,g1\n");
        let m_in = load_group_mapping(f.path(), &sample_cols(&["A", "B", "C"])).unwrap();
        let pre_groups = m_in.groups();
        let pre_assigned = m_in.assigned_count();
        let pre_sample_count = m_in.sample_count();
        let pre_g1 = m_in.samples_in("g1");

        let _m_out = m_in.without_unassigned_samples();

        assert_eq!(m_in.groups(), pre_groups);
        assert_eq!(m_in.assigned_count(), pre_assigned);
        assert_eq!(m_in.sample_count(), pre_sample_count);
        assert_eq!(m_in.samples_in("g1"), pre_g1);
    }

    #[test]
    fn existing_2_column_form_still_parses() {
        let f = write_csv(
            "sample,group,dry_weight\n\
             S01,A,1.0\n\
             S02,B,2.0\n",
        );
        let m = load_group_mapping(f.path(), &sample_cols(&["S01", "S02"])).unwrap();
        assert_eq!(m.groups(), vec!["A".to_string(), "B".to_string()]);
        assert_eq!(m.metadata_column_names(), vec!["dry_weight".to_string()]);
        assert_eq!(
            m.metadata_values("dry_weight"),
            Some([Some(1.0), Some(2.0)].as_ref())
        );
        // No biosample column → biosample_of returns None for everyone.
        assert_eq!(m.biosample_of("S01"), None);
    }

    #[test]
    fn biosample_count_distinct_within_group() {
        // control: S01+S02 share biosample BIO-A (technical replicates), S03 is
        // BIO-B → 2 distinct biosamples across 3 samples. treatment: 1 sample.
        let f = write_csv(
            "sample,biosample,group\n\
             S01,BIO-A,control\n\
             S02,BIO-A,control\n\
             S03,BIO-B,control\n\
             S04,BIO-C,treatment\n",
        );
        let m = load_group_mapping(f.path(), &sample_cols(&["S01", "S02", "S03", "S04"])).unwrap();
        assert_eq!(m.samples_in("control").len(), 3);
        assert_eq!(m.biosample_count("control"), Some(2));
        assert_eq!(m.biosample_count("treatment"), Some(1));
    }

    #[test]
    fn biosample_count_none_without_biosample_column() {
        let f = write_csv("sample,group\nS01,control\nS02,control\n");
        let m = load_group_mapping(f.path(), &sample_cols(&["S01", "S02"])).unwrap();
        assert_eq!(m.biosample_count("control"), None);
    }

    #[test]
    fn unmatched_csv_samples_retained_and_listed() {
        // S99 names a sample absent from sample_cols → CSV-only.
        let f = write_csv("sample,group\nS01,ASAP\nS99,Other\n");
        let m = load_group_mapping(f.path(), &sample_cols(&["S01"])).unwrap();
        // The row is still ignored from the mapping...
        assert!(!m.groups().contains(&"Other".to_string()));
        assert_eq!(m.group_of("S99"), UNASSIGNED);
        // ...but its name is now retained.
        assert_eq!(m.unmatched_csv_samples(), &["S99".to_string()]);
    }

    #[test]
    fn unmatched_csv_samples_empty_when_fully_matched() {
        let f = write_csv("sample,group\nS01,ASAP\nS02,CK\n");
        let m = load_group_mapping(f.path(), &sample_cols(&["S01", "S02"])).unwrap();
        assert!(m.unmatched_csv_samples().is_empty());
    }

    #[test]
    fn unmatched_csv_samples_lists_all_on_no_overlap() {
        // Every CSV sample is CSV-only (no overlap with sample_cols).
        let f = write_csv("sample,group\nX01,A\nX02,B\n");
        let m = load_group_mapping(f.path(), &sample_cols(&["S01", "S02"])).unwrap();
        assert_eq!(m.assigned_count(), 0);
        assert_eq!(
            m.unmatched_csv_samples(),
            &["X01".to_string(), "X02".to_string()]
        );
    }

    #[test]
    fn unmatched_csv_samples_sorted_deterministically() {
        // CSV row order Z9, A1, M5 — none in sample_cols → all CSV-only.
        // Byte-lexicographic sort regardless of row / HashMap iteration order.
        let f = write_csv("sample,group\nZ9,A\nA1,B\nM5,C\n");
        let m = load_group_mapping(f.path(), &sample_cols(&["S01"])).unwrap();
        assert_eq!(
            m.unmatched_csv_samples(),
            &["A1".to_string(), "M5".to_string(), "Z9".to_string()]
        );
    }

    #[test]
    fn unmatched_csv_samples_carried_through_without_unassigned() {
        // A, B assigned to g1; C is Unassigned (.txt-only); S99 is CSV-only.
        let f = write_csv("sample,group\nA,g1\nB,g1\nS99,Other\n");
        let m = load_group_mapping(f.path(), &sample_cols(&["A", "B", "C"])).unwrap();
        assert_eq!(m.unmatched_csv_samples(), &["S99".to_string()]);
        // The Stage 1 → Stage 2 narrowing preserves the CSV-only list verbatim.
        let narrowed = m.without_unassigned_samples();
        assert_eq!(narrowed.unmatched_csv_samples(), m.unmatched_csv_samples());
    }
}
