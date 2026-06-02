//! Data types for Stage 3 enrichment analysis output.

use serde::{Deserialize, Serialize};

use crate::dam::fdr::FdrMethod;

/// Which DAM features feed the K set (the "foreground" of the
/// hypergeometric test).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum EnrichmentDirection {
    Up,
    Down,
    Both,
}

impl EnrichmentDirection {
    pub fn display_name(&self) -> &'static str {
        match self {
            Self::Up => "Up only",
            Self::Down => "Down only",
            Self::Both => "Both (Up + Down)",
        }
    }

    pub fn short_label(&self) -> &'static str {
        match self {
            Self::Up => "Up",
            Self::Down => "Down",
            Self::Both => "Both",
        }
    }
}

/// One row of the ORA result. Rows are returned sorted by ascending FDR;
/// ties broken by ascending `entry_id`.
#[derive(Debug, Clone, PartialEq)]
pub struct EnrichmentRow {
    /// Identifier of the compound-set entry (pathway ID like `gmx00010`
    /// in pathway mode, or module ID like `M00001` in module mode).
    pub entry_id: String,
    pub entry_name: String,
    /// Number of compounds in `dam_cpd` that hit this entry.
    pub hits: usize,
    /// Number of compounds in this entry that are also in the universe.
    pub total: usize,
    /// Expected hits under the null: `K * total / N` (or 0 if N == 0).
    pub expected: f64,
    /// `hits / expected` (or NaN if expected == 0).
    pub enrichment_ratio: f64,
    pub p_value: f64,
    pub fdr: f64,
    /// Compound IDs that hit this entry, sorted alphabetically.
    pub hit_kegg_ids: Vec<String>,
    /// Whether this row passes the post-FDR `min_hit_count` display
    /// filter (i.e. `hits >= min_hit_count`). FDR was computed over all
    /// rows regardless of this flag.
    pub displayed: bool,
}

/// Entry-level summary of the enrichment Run.
#[derive(Debug, Clone, PartialEq)]
pub struct EnrichmentResult {
    pub universe_size: usize,
    pub dam_cpd_size: usize,
    pub direction: EnrichmentDirection,
    pub min_hit_count: usize,
    /// Pre-FDR entry-size threshold applied in this run. Entries with
    /// `m_p < min_entry_size` were dropped before FDR — they do NOT
    /// appear in `rows` and do NOT contribute to the FDR family `m`.
    /// Added by `add-min-entry-size-filter`.
    pub min_entry_size: usize,
    /// Count of input entries dropped by the pre-FDR `min_entry_size`
    /// filter (those with `m_p < min_entry_size`). Stage 3 result panel
    /// surfaces this via `<retained> / <total>` arithmetic where
    /// `total = rows.len() + entries_dropped_by_min_entry_size`.
    pub entries_dropped_by_min_entry_size: usize,
    /// Count of input entries that had `compounds.is_empty()`. Most
    /// relevant for module mode, where some KEGG modules have no
    /// COMPOUND block. Pathway mode typically reports 0 here.
    pub empty_compound_count: usize,
    pub rows: Vec<EnrichmentRow>,
    /// FDR correction applied to the per-entry p values. Carried so the
    /// CSV exporter and downstream renderers can surface the choice.
    pub fdr_method: FdrMethod,
}
