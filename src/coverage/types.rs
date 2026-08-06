//! Data types for the KEGG coverage survey — the descriptive counterpart to
//! `crate::enrichment::types`.
//!
//! **The defining property of this module is negative.** Neither [`CoverageRow`]
//! nor [`CoverageResult`] carries a p-value, q-value, FDR method, expected
//! count, enrichment ratio, or direction, and neither reuses `EnrichmentRow` /
//! `EnrichmentResult`. That is a structural guarantee, not a convention: the
//! coverage route performs no statistical test, so no coverage surface — table,
//! dot plot, CSV, Data tab, log — has a field through which it could present
//! detectability as significance. Owner: the `kegg-coverage` capability.

use serde::{Deserialize, Serialize};

/// One row of the coverage survey: a single KEGG entry (pathway or module) and
/// how much of it the detected metabolome covers.
///
/// `coverage` and `share` answer different questions and are deliberately not
/// comparable: `coverage` is "what fraction of THIS ENTRY did I detect",
/// `share` is "what fraction of MY METABOLOME is in this entry". Only
/// `coverage` is displayed; `share` is exported only (owner: `kegg-coverage`).
#[derive(Debug, Clone, PartialEq)]
pub struct CoverageRow {
    /// Pathway ID (`hsa00010`) in Pathway mode, module ID (`M00001`) in Module
    /// mode.
    pub entry_id: String,
    pub entry_name: String,
    /// `|C|` — the number of distinct KEGG compounds this entry contains.
    /// Zero for KEGG's global/overview maps, which carry no `COMPOUND` section.
    pub entry_size: usize,
    /// `|D ∩ C|`.
    pub hits: usize,
    /// `hits / entry_size`, or `0.0` when `entry_size == 0`.
    pub coverage: f64,
    /// `hits / |D|`, or `0.0` when `|D| == 0`. Shares do NOT sum to 1.0 across
    /// rows — entries share compounds.
    pub share: f64,
    /// The sorted `D ∩ C` cpd IDs. Bare IDs: attaching the user's MS-DIAL
    /// metabolite names is a presentation concern handled by the CSV exporter.
    pub hit_compounds: Vec<String>,
}

/// The complete coverage survey: one row per catalogue entry plus the
/// provenance counts the result screen reports.
#[derive(Debug, Clone, PartialEq)]
pub struct CoverageResult {
    /// One row per entry given to [`compute`](super::compute), in the input
    /// order. Includes zero-hit and zero-compound entries — every filter is a
    /// display concern applied later by [`displayed_rows`](super::displayed_rows).
    pub rows: Vec<CoverageRow>,
    /// `|D|` — the detected, KEGG-mapped compound set.
    pub detected_total: usize,
    /// Every entry in the catalogue; equals `rows.len()`.
    pub entries_total: usize,
    /// Entries with `entry_size == 0`, reported so the result screen can say
    /// how much of the catalogue is unusable rather than swallowing it.
    pub entries_without_compounds: usize,
    /// `|D ∩ (⋃ over ALL entries of entry.compounds)|`.
    ///
    /// A **provenance** count over the whole catalogue, computed before any
    /// display filter — it answers "how much of my mapped metabolome does this
    /// catalogue know about at all", so raising `coverage_min_entry_size` or
    /// `min_hit_count` on the result screen MUST NOT move it.
    ///
    /// The name collides with `Stage3Funnel::detected_in_entries`, which is
    /// FOREGROUND-scoped on the enrichment route. Different quantities on
    /// different types; the coverage route never reads `Stage3Funnel`.
    pub detected_in_entries: usize,
}

/// Sort key for the coverage results table. Owner: the `kegg-coverage` capability.
///
/// All keys sort descending except [`CoverageSortKey::EntryId`], which sorts
/// ascending; ties break by descending hits then ascending entry id.
///
/// Reachable from `SessionSettings`, so the serialised variant names are part
/// of the on-disk snapshot contract.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq, Serialize, Deserialize)]
pub enum CoverageSortKey {
    /// `hits / entry_size` — the default, and the reason the minimum-entry-size
    /// floor exists (a 2-compound entry with 2 hits would otherwise top the
    /// table at 100%).
    #[default]
    Coverage,
    Hits,
    EntrySize,
    EntryId,
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The negative guarantee, asserted by construction: every field of
    /// `CoverageRow` and `CoverageResult` is named here, so adding a
    /// `p_value` / `fdr` / `expected` / `enrichment_ratio` / `direction`
    /// field fails to compile against these literals rather than silently
    /// giving a coverage surface a statistic to render.
    #[test]
    fn result_types_carry_no_statistical_field() {
        let CoverageRow {
            entry_id: _,
            entry_name: _,
            entry_size: _,
            hits: _,
            coverage: _,
            share: _,
            hit_compounds: _,
        } = CoverageRow {
            entry_id: "hsa00010".into(),
            entry_name: "Glycolysis".into(),
            entry_size: 3,
            hits: 1,
            coverage: 1.0 / 3.0,
            share: 0.5,
            hit_compounds: vec!["C00031".into()],
        };

        let CoverageResult {
            rows: _,
            detected_total: _,
            entries_total: _,
            entries_without_compounds: _,
            detected_in_entries: _,
        } = CoverageResult {
            rows: vec![],
            detected_total: 0,
            entries_total: 0,
            entries_without_compounds: 0,
            detected_in_entries: 0,
        };
    }

    /// `CoverageSortKey` derives `Default` (required by `#[default]`) and lands
    /// on `Coverage`.
    #[test]
    fn sort_key_defaults_to_coverage() {
        assert_eq!(CoverageSortKey::default(), CoverageSortKey::Coverage);
    }
}
