//! Build the `OrganismGroupIndex` from a loaded organism list.
//!
//! KEGG `/list/organism` returns one row per organism with a 4-level
//! semicolon-delimited taxonomy in column 4 (e.g.
//! `Eukaryotes;Animals;Mammals;Primates`). Module mode lets the user pick
//! one of the first three levels and a Group within that level; this
//! function precomputes the (level, group) → set-of-codes map so the UI
//! can populate the Group dropdown without re-walking the organism list.
//!
//! Level 1 has 2 candidates (Prokaryotes / Eukaryotes) in the current
//! KEGG dataset. Level 2 has 6 (Bacteria, Animals, Archaea, Fungi, Plants,
//! Protists). Level 3 has tens. Level 4 (e.g. Primates) is not exposed.

use chrono::{DateTime, Utc};
use std::collections::{HashMap, HashSet};
use tracing::warn;

use crate::kegg::types::{KeggOrganism, OrganismGroupIndex};

/// Build a fresh `OrganismGroupIndex` from the loaded organism list.
///
/// `fetched_at` is carried through from the caller (typically the
/// organisms cache's own timestamp) so the two on-disk caches share a
/// coherent reference time.
///
/// Organisms whose lineage has fewer than 3 semicolon-separated levels
/// (should not happen with current KEGG data, but defensive) are
/// inserted into the available levels and skipped for missing ones,
/// with a WARN log noting the malformed lineage.
pub fn build_organism_group_index(
    organisms: &[KeggOrganism],
    fetched_at: DateTime<Utc>,
) -> OrganismGroupIndex {
    let mut by_level: [HashMap<String, HashSet<String>>; 3] =
        [HashMap::new(), HashMap::new(), HashMap::new()];

    for org in organisms {
        let levels: Vec<&str> = org.lineage.split(';').map(str::trim).collect();
        if levels.len() < 3 {
            warn!(
                code = %org.code,
                lineage = %org.lineage,
                "organism lineage has fewer than 3 levels; indexing what's available"
            );
        }
        for (idx, level_name) in levels.iter().enumerate().take(3) {
            if level_name.is_empty() {
                continue;
            }
            by_level[idx]
                .entry(level_name.to_string())
                .or_default()
                .insert(org.code.clone());
        }
    }

    OrganismGroupIndex {
        fetched_at,
        by_level,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn org(code: &str, lineage: &str) -> KeggOrganism {
        KeggOrganism {
            code: code.to_string(),
            t_number: "T00000".to_string(),
            name: format!("Organism {code}"),
            lineage: lineage.to_string(),
        }
    }

    #[test]
    fn three_level_decomposition_inserts_into_each_level() {
        let orgs = vec![org("hsa", "Eukaryotes;Animals;Mammals;Primates")];
        let now = Utc::now();
        let index = build_organism_group_index(&orgs, now);
        assert!(index.by_level[0].contains_key("Eukaryotes"));
        assert!(index.by_level[1].contains_key("Animals"));
        assert!(index.by_level[2].contains_key("Mammals"));
        // Level 4 (Primates) is NOT indexed.
        assert!(!index.by_level[2].contains_key("Primates"));
        assert_eq!(index.fetched_at, now);
    }

    #[test]
    fn multiple_organisms_share_groups() {
        let orgs = vec![
            org("hsa", "Eukaryotes;Animals;Mammals;Primates"),
            org("mmu", "Eukaryotes;Animals;Mammals;Rodents"),
            org("ath", "Eukaryotes;Plants;Eudicots;Brassicales"),
        ];
        let index = build_organism_group_index(&orgs, Utc::now());
        // Level 1: all 3 under Eukaryotes.
        assert_eq!(index.by_level[0]["Eukaryotes"].len(), 3);
        // Level 2: 2 Animals, 1 Plants.
        assert_eq!(index.by_level[1]["Animals"].len(), 2);
        assert_eq!(index.by_level[1]["Plants"].len(), 1);
        // Level 3: 2 Mammals, 1 Eudicots.
        assert_eq!(index.by_level[2]["Mammals"].len(), 2);
        assert_eq!(index.by_level[2]["Eudicots"].len(), 1);
    }

    #[test]
    fn malformed_lineage_with_fewer_levels_still_indexes_available() {
        let orgs = vec![
            org("good", "Eukaryotes;Animals;Mammals;Primates"),
            org("bad", "Eukaryotes;Animals"), // only 2 levels
        ];
        let index = build_organism_group_index(&orgs, Utc::now());
        // Both inserted into level 1 + 2.
        assert_eq!(index.by_level[0]["Eukaryotes"].len(), 2);
        assert_eq!(index.by_level[1]["Animals"].len(), 2);
        // Only `good` in level 3.
        assert_eq!(index.by_level[2]["Mammals"].len(), 1);
    }

    #[test]
    fn empty_organisms_produces_empty_index() {
        let index = build_organism_group_index(&[], Utc::now());
        assert!(index.by_level[0].is_empty());
        assert!(index.by_level[1].is_empty());
        assert!(index.by_level[2].is_empty());
    }
}
