use ndarray::{Array2, Axis};

use crate::data::groups::{GroupMapping, UNASSIGNED};

#[derive(Debug, Clone)]
pub struct FeatureMeta {
    pub alignment_id: String,
    pub metabolite_name: String,
    pub inchikey: Option<String>,
    /// MS-DIAL `Adduct type` cell, e.g. `[M+H]+` / `[M-H]-`. `None` when the
    /// source file lacks the column or the cell is empty. Consumed by
    /// `infer_polarity` for the dual-mode adduct-column sanity check, and by
    /// `crate::dedup` to rank adduct quality during InChIKey deduplication.
    pub adduct_type: Option<String>,
    pub average_rt_min: Option<f64>,
    pub average_mz: Option<f64>,
    pub formula: Option<String>,
    pub smiles: Option<String>,
    /// MS-DIAL `Fill %`. Range 0..=100. Used by `crate::dedup` cascade level 3a
    /// (data-quality tiebreak — higher is better).
    pub fill_percent: Option<f64>,
    /// MS-DIAL `MS/MS matched`. `Some(true)` / `Some(false)` / `None`. Used by
    /// `crate::dedup` cascade level 1a (annotation confidence — True > False > None).
    pub ms_ms_matched: Option<bool>,
    /// MS-DIAL `Isotope tracking weight number`. `0` is the M0 monoisotopic peak;
    /// any positive value flags a natural-abundance isotope peak. Used by
    /// `crate::dedup` adduct classification to detect isotope rows.
    pub isotope_tracking_weight_number: Option<i32>,
    /// MS-DIAL `Total score`. Combined annotation-confidence score. Used by
    /// `crate::dedup` cascade level 1b (the vendor-computed weighted composite of
    /// every spectral-similarity metric, including dot products).
    pub total_score: Option<f64>,
    /// MS-DIAL `S/N average`. Used by `crate::dedup` cascade level 3b
    /// (data-quality tiebreak — higher is better).
    pub sn_average: Option<f64>,
}

#[derive(Debug)]
pub struct MetabolomicsTable {
    pub features: Vec<FeatureMeta>,
    pub sample_cols: Vec<String>,
    /// Raw intensity matrix as parsed from the MS-DIAL file. Never mutated
    /// after Stage 1 load. Stage 2's DAM runner uses this as the source for
    /// every per-run normalization.
    pub intensity_raw: Array2<f64>,
    /// Working intensity matrix. Equal to `intensity_raw` until Stage 2's
    /// `run_dam` writes a normalized matrix here. All downstream consumers
    /// (DAM stats, Stage 3 enrichment) read from this field.
    pub intensity: Array2<f64>,
    /// Columns that were intentionally excluded from `sample_cols` because their
    /// `File type` value is not `Sample`. Each entry is `(column_name, file_type)`.
    /// Used by the UI to surface what the parser dropped (e.g. `Blank` process blanks)
    /// so users don't think their data is missing.
    pub excluded_cols: Vec<(String, String)>,
    /// Count of features carrying a non-null `inchikey` ("annotated"), computed
    /// ONCE at construction (the Unknown count is `features.len() - annotated_count`).
    /// `without_unassigned_samples` preserves it verbatim (the feature axis is
    /// untouched). Replaces the old per-call `annotated_count()` scan. Owner: the `msdial-input` capability.
    pub annotated_count: usize,
}

impl MetabolomicsTable {
    /// Stage-1 → Stage-2 boundary helper. Returns a new owned table whose
    /// sample axis is narrowed to samples whose
    /// `mapping.group_of(name) != UNASSIGNED`. `features` and `excluded_cols`
    /// are cloned verbatim (the feature axis is untouched; `excluded_cols`
    /// is name-keyed, not index-keyed, so verbatim cloning is safe). Owner: the `msdial-input` capability.
    pub fn without_unassigned_samples(&self, mapping: &GroupMapping) -> Self {
        let kept: Vec<usize> = self
            .sample_cols
            .iter()
            .enumerate()
            .filter(|(_, name)| mapping.group_of(name) != UNASSIGNED)
            .map(|(i, _)| i)
            .collect();

        let sample_cols: Vec<String> = kept.iter().map(|&i| self.sample_cols[i].clone()).collect();
        let intensity_raw = self.intensity_raw.select(Axis(1), &kept);
        let intensity = self.intensity.select(Axis(1), &kept);

        Self {
            annotated_count: self.annotated_count,
            features: self.features.clone(),
            sample_cols,
            intensity_raw,
            intensity,
            excluded_cols: self.excluded_cols.clone(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::data::load_group_mapping;
    use ndarray::array;
    use std::io::Write;

    fn write_csv(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().expect("create tempfile");
        f.write_all(content.as_bytes()).expect("write fixture");
        f
    }

    fn cols(names: &[&str]) -> Vec<String> {
        names.iter().map(|s| s.to_string()).collect()
    }

    fn feature_meta(id: &str) -> FeatureMeta {
        FeatureMeta {
            alignment_id: id.into(),
            metabolite_name: format!("met-{id}"),
            inchikey: None,
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

    #[test]
    fn without_unassigned_round_trip_on_zero_unassigned_is_clone_equivalent() {
        let f = write_csv("sample,group\nA,g1\nB,g1\n");
        let mapping = load_group_mapping(f.path(), &cols(&["A", "B"])).unwrap();
        let intensity_raw = array![[1.0, 2.0], [3.0, 4.0]];
        let intensity = intensity_raw.clone();
        let t_in = MetabolomicsTable {
            annotated_count: 0,
            features: vec![feature_meta("0"), feature_meta("1")],
            sample_cols: cols(&["A", "B"]),
            intensity_raw: intensity_raw.clone(),
            intensity,
            excluded_cols: vec![],
        };

        let t_out = t_in.without_unassigned_samples(&mapping);
        assert_eq!(t_out.sample_cols, t_in.sample_cols);
        assert_eq!(t_out.intensity_raw.shape(), t_in.intensity_raw.shape());
        for ((a, b), c) in t_out
            .intensity_raw
            .iter()
            .zip(t_in.intensity_raw.iter())
            .zip(t_out.intensity.iter())
        {
            assert_eq!(a, b);
            assert_eq!(a, c);
        }
        assert_eq!(t_out.features.len(), t_in.features.len());
        assert_eq!(t_out.excluded_cols, t_in.excluded_cols);
    }

    #[test]
    fn without_unassigned_drops_unassigned_columns_and_preserves_features() {
        let f = write_csv("sample,group\nA,g1\nB,g1\nC,g2\n");
        // sample_cols include 2 extras (Blank01, Blank02) NOT in the CSV → Unassigned.
        let sample_cols_v = cols(&["A", "Blank01", "B", "Blank02", "C"]);
        let mapping = load_group_mapping(f.path(), &sample_cols_v).unwrap();
        let intensity_raw = array![[1.0, 10.0, 2.0, 20.0, 3.0], [4.0, 40.0, 5.0, 50.0, 6.0]];
        let intensity = intensity_raw.clone();
        let t_in = MetabolomicsTable {
            annotated_count: 0,
            features: vec![feature_meta("0"), feature_meta("1")],
            sample_cols: sample_cols_v,
            intensity_raw: intensity_raw.clone(),
            intensity,
            excluded_cols: vec![("NA-1".into(), "NA".into())],
        };

        let t_out = t_in.without_unassigned_samples(&mapping);
        assert_eq!(t_out.sample_cols, cols(&["A", "B", "C"]));
        assert_eq!(t_out.intensity_raw.shape(), &[2, 3]);
        // Verify the kept columns are 0, 2, 4 of the input.
        assert_eq!(t_out.intensity_raw[[0, 0]], 1.0);
        assert_eq!(t_out.intensity_raw[[0, 1]], 2.0);
        assert_eq!(t_out.intensity_raw[[0, 2]], 3.0);
        assert_eq!(t_out.intensity_raw[[1, 0]], 4.0);
        assert_eq!(t_out.intensity_raw[[1, 1]], 5.0);
        assert_eq!(t_out.intensity_raw[[1, 2]], 6.0);
        // intensity narrowed identically.
        for ((i, j), v) in t_out.intensity.indexed_iter() {
            assert_eq!(*v, t_out.intensity_raw[[i, j]]);
        }
        // Features and excluded_cols cloned verbatim.
        assert_eq!(t_out.features.len(), 2);
        assert_eq!(t_out.excluded_cols, vec![("NA-1".into(), "NA".into())]);
    }

    #[test]
    fn without_unassigned_preserves_nan_cells_at_new_positions() {
        let f = write_csv("sample,group\nA,g1\nC,g2\n");
        let sample_cols_v = cols(&["A", "Blank", "C"]);
        let mapping = load_group_mapping(f.path(), &sample_cols_v).unwrap();
        let intensity_raw = array![
            [1.0, 99.0, f64::NAN],
            [f64::NAN, 99.0, 6.0],
            [3.0, 99.0, 7.0]
        ];
        let intensity = intensity_raw.clone();
        let t_in = MetabolomicsTable {
            annotated_count: 0,
            features: vec![feature_meta("0"), feature_meta("1"), feature_meta("2")],
            sample_cols: sample_cols_v,
            intensity_raw,
            intensity,
            excluded_cols: vec![],
        };

        let t_out = t_in.without_unassigned_samples(&mapping);
        assert_eq!(t_out.sample_cols, cols(&["A", "C"]));
        // NaN at (0, 2) of input → (0, 1) of output (column 1 dropped).
        assert!(t_out.intensity_raw[[0, 1]].is_nan());
        // NaN at (1, 0) of input → (1, 0) of output.
        assert!(t_out.intensity_raw[[1, 0]].is_nan());
        // Real values preserved.
        assert_eq!(t_out.intensity_raw[[0, 0]], 1.0);
        assert_eq!(t_out.intensity_raw[[1, 1]], 6.0);
    }

    #[test]
    fn without_unassigned_does_not_mutate_source_table() {
        let f = write_csv("sample,group\nA,g1\n");
        let sample_cols_v = cols(&["A", "Blank"]);
        let mapping = load_group_mapping(f.path(), &sample_cols_v).unwrap();
        let intensity_raw = array![[1.0, 99.0]];
        let intensity = intensity_raw.clone();
        let t_in = MetabolomicsTable {
            annotated_count: 0,
            features: vec![feature_meta("0")],
            sample_cols: sample_cols_v.clone(),
            intensity_raw: intensity_raw.clone(),
            intensity,
            excluded_cols: vec![],
        };

        let _t_out = t_in.without_unassigned_samples(&mapping);
        assert_eq!(t_in.sample_cols, sample_cols_v);
        assert_eq!(t_in.intensity_raw.shape(), &[1, 2]);
        assert_eq!(t_in.intensity.shape(), &[1, 2]);
        assert_eq!(t_in.intensity_raw[[0, 0]], 1.0);
        assert_eq!(t_in.intensity_raw[[0, 1]], 99.0);
    }

    fn feature_with_key(id: &str, key: Option<&str>) -> FeatureMeta {
        FeatureMeta {
            inchikey: key.map(|k| k.to_string()),
            ..feature_meta(id)
        }
    }

    fn table_with_features(features: Vec<FeatureMeta>) -> MetabolomicsTable {
        let n = features.len();
        MetabolomicsTable {
            annotated_count: features.iter().filter(|f| f.inchikey.is_some()).count(),
            features,
            sample_cols: cols(&["A"]),
            intensity_raw: Array2::zeros((n, 1)),
            intensity: Array2::zeros((n, 1)),
            excluded_cols: vec![],
        }
    }

    #[test]
    fn annotated_count_mixed() {
        let t = table_with_features(vec![
            feature_with_key("0", Some("AAA")),
            feature_with_key("1", None),
            feature_with_key("2", Some("BBB")),
            feature_with_key("3", None),
        ]);
        assert_eq!(t.annotated_count, 2);
        assert_eq!(t.features.len() - t.annotated_count, 2); // Unknown
    }

    #[test]
    fn annotated_count_all_unknown_and_all_annotated() {
        let all_unknown = table_with_features(vec![
            feature_with_key("0", None),
            feature_with_key("1", None),
        ]);
        assert_eq!(all_unknown.annotated_count, 0);

        let all_annotated = table_with_features(vec![
            feature_with_key("0", Some("AAA")),
            feature_with_key("1", Some("BBB")),
        ]);
        assert_eq!(all_annotated.annotated_count, 2);
    }
}
