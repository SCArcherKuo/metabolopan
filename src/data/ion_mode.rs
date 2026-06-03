use std::fmt;
use std::path::PathBuf;

use anyhow::{Result, bail};
use serde::{Deserialize, Serialize};

use crate::data::groups::GroupMapping;
use crate::data::types::MetabolomicsTable;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum IonMode {
    Positive,
    Negative,
}

impl IonMode {
    /// The opposite ionization polarity. A property of the polarity type
    /// (used by the Stage 1 slot-2 auto-fill), not of the Stage 1 screen
    /// (`move-labels-onto-types` D3).
    pub fn opposite(self) -> IonMode {
        match self {
            IonMode::Positive => IonMode::Negative,
            IonMode::Negative => IonMode::Positive,
        }
    }
}

impl fmt::Display for IonMode {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            IonMode::Positive => "POS",
            IonMode::Negative => "NEG",
        })
    }
}

#[derive(Debug)]
pub struct IonModeTable {
    pub mode: IonMode,
    pub table: MetabolomicsTable,
    pub txt_path: Option<PathBuf>,
}

impl IonModeTable {
    /// Stage-1 → Stage-2 boundary helper. Returns a new owned `IonModeTable`
    /// whose inner `MetabolomicsTable` has been narrowed via
    /// `MetabolomicsTable::without_unassigned_samples`. `mode` and `txt_path`
    /// are cloned verbatim. See the `dual-mode-input` capability spec for
    /// the full contract.
    pub fn without_unassigned_samples(&self, mapping: &GroupMapping) -> Self {
        Self {
            mode: self.mode,
            table: self.table.without_unassigned_samples(mapping),
            txt_path: self.txt_path.clone(),
        }
    }
}

#[derive(Debug)]
pub struct IonModeTables(Vec<IonModeTable>);

impl IonModeTables {
    pub fn try_new(mut items: Vec<IonModeTable>) -> Result<Self> {
        match items.len() {
            0 => bail!("ion-mode bundle requires at least one IonModeTable, got 0"),
            1 => Ok(Self(items)),
            2 => {
                if items[0].mode == items[1].mode {
                    bail!(
                        "ion-mode bundle requires distinct modes; both slots are {}",
                        items[0].mode
                    );
                }
                if items[0].mode == IonMode::Negative {
                    items.swap(0, 1);
                }
                Ok(Self(items))
            }
            n => bail!("ion-mode bundle accepts 1 or 2 IonModeTables, got {n}"),
        }
    }

    pub fn into_inner(self) -> Vec<IonModeTable> {
        self.0
    }

    pub fn as_slice(&self) -> &[IonModeTable] {
        &self.0
    }

    pub fn len(&self) -> usize {
        self.0.len()
    }

    pub fn is_empty(&self) -> bool {
        self.0.is_empty()
    }

    pub fn iter(&self) -> std::slice::Iter<'_, IonModeTable> {
        self.0.iter()
    }
}

impl std::ops::Index<usize> for IonModeTables {
    type Output = IonModeTable;
    fn index(&self, i: usize) -> &IonModeTable {
        &self.0[i]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ndarray::Array2;

    fn make_table() -> MetabolomicsTable {
        let intensity = Array2::<f64>::zeros((0, 0));
        MetabolomicsTable {
            annotated_count: 0,
            features: vec![],
            sample_cols: vec![],
            intensity_raw: intensity.clone(),
            intensity,
            excluded_cols: vec![],
        }
    }

    fn entry(mode: IonMode) -> IonModeTable {
        IonModeTable {
            mode,
            table: make_table(),
            txt_path: None,
        }
    }

    #[test]
    fn display_uses_short_labels() {
        assert_eq!(IonMode::Positive.to_string(), "POS");
        assert_eq!(IonMode::Negative.to_string(), "NEG");
    }

    #[test]
    fn opposite_flips_polarity() {
        assert_eq!(IonMode::Positive.opposite(), IonMode::Negative);
        assert_eq!(IonMode::Negative.opposite(), IonMode::Positive);
    }

    #[test]
    fn try_new_accepts_single_entry() {
        let b = IonModeTables::try_new(vec![entry(IonMode::Positive)]).unwrap();
        assert_eq!(b.len(), 1);
        assert_eq!(b[0].mode, IonMode::Positive);
    }

    #[test]
    fn try_new_canonical_orders_positive_first() {
        let b = IonModeTables::try_new(vec![entry(IonMode::Negative), entry(IonMode::Positive)])
            .unwrap();
        assert_eq!(b.len(), 2);
        assert_eq!(b[0].mode, IonMode::Positive);
        assert_eq!(b[1].mode, IonMode::Negative);
    }

    #[test]
    fn try_new_already_canonical_unchanged() {
        let b = IonModeTables::try_new(vec![entry(IonMode::Positive), entry(IonMode::Negative)])
            .unwrap();
        assert_eq!(b[0].mode, IonMode::Positive);
        assert_eq!(b[1].mode, IonMode::Negative);
    }

    #[test]
    fn try_new_rejects_duplicate_modes() {
        let err = IonModeTables::try_new(vec![entry(IonMode::Positive), entry(IonMode::Positive)])
            .unwrap_err();
        assert!(err.to_string().contains("distinct modes"));
    }

    #[test]
    fn try_new_rejects_zero_and_three() {
        assert!(IonModeTables::try_new(vec![]).is_err());
        let three = vec![
            entry(IonMode::Positive),
            entry(IonMode::Negative),
            entry(IonMode::Positive),
        ];
        assert!(IonModeTables::try_new(three).is_err());
    }

    fn make_table_with_cols(names: &[&str]) -> MetabolomicsTable {
        let n = names.len();
        let intensity = Array2::<f64>::zeros((1, n));
        MetabolomicsTable {
            annotated_count: 0,
            features: vec![],
            sample_cols: names.iter().map(|s| s.to_string()).collect(),
            intensity_raw: intensity.clone(),
            intensity,
            excluded_cols: vec![],
        }
    }

    fn load_mapping_for(csv: &str, sample_cols: &[&str]) -> GroupMapping {
        use crate::data::load_group_mapping;
        use std::io::Write;
        let mut f = tempfile::NamedTempFile::new().expect("tempfile");
        f.write_all(csv.as_bytes()).expect("write fixture");
        let owned: Vec<String> = sample_cols.iter().map(|s| s.to_string()).collect();
        load_group_mapping(f.path(), &owned).expect("mapping")
    }

    #[test]
    fn without_unassigned_wrapper_preserves_mode_and_txt_path() {
        let mapping = load_mapping_for("sample,group\nS1,g1\nS2,g1\n", &["S1", "S2"]);
        let txt_path = Some(PathBuf::from("/tmp/data/neg.txt"));
        let it_in = IonModeTable {
            mode: IonMode::Negative,
            table: make_table_with_cols(&["S1", "S2"]),
            txt_path: txt_path.clone(),
        };

        let it_out = it_in.without_unassigned_samples(&mapping);
        assert_eq!(it_out.mode, IonMode::Negative);
        assert_eq!(it_out.txt_path, txt_path);
    }

    #[test]
    fn without_unassigned_wrapper_delegates_table_filtering() {
        let mapping = load_mapping_for("sample,group\nS1,g1\nS2,g1\n", &["S1", "Blank", "S2"]);
        let it_in = IonModeTable {
            mode: IonMode::Positive,
            table: make_table_with_cols(&["S1", "Blank", "S2"]),
            txt_path: None,
        };

        let it_out = it_in.without_unassigned_samples(&mapping);
        assert_eq!(
            it_out.table.sample_cols,
            vec!["S1".to_string(), "S2".into()]
        );
        assert_eq!(it_out.table.intensity_raw.ncols(), 2);
    }

    #[test]
    fn without_unassigned_wrapper_does_not_mutate_source() {
        let mapping = load_mapping_for("sample,group\nS1,g1\n", &["S1", "Blank"]);
        let pre_path = Some(PathBuf::from("/tmp/data/pos.txt"));
        let it_in = IonModeTable {
            mode: IonMode::Positive,
            table: make_table_with_cols(&["S1", "Blank"]),
            txt_path: pre_path.clone(),
        };

        let _it_out = it_in.without_unassigned_samples(&mapping);
        assert_eq!(it_in.mode, IonMode::Positive);
        assert_eq!(it_in.txt_path, pre_path);
        assert_eq!(
            it_in.table.sample_cols,
            vec!["S1".to_string(), "Blank".into()]
        );
        assert_eq!(it_in.table.intensity_raw.ncols(), 2);
    }
}
