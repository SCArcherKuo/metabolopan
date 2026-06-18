use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct KeggOrganism {
    pub code: String,
    pub t_number: String,
    pub name: String,
    pub lineage: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct KeggCompoundSet {
    pub id: String,
    pub name: String,
    pub compounds: Vec<String>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SpeciesKegg {
    pub code: String,
    pub fetched_at: DateTime<Utc>,
    pub pathways: Vec<KeggCompoundSet>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct OrganismsCache {
    pub fetched_at: DateTime<Utc>,
    pub organisms: Vec<KeggOrganism>,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct KeggProgress {
    pub completed: usize,
    pub total: usize,
    pub current_pathway: String,
}

#[derive(Debug)]
pub enum KeggEvent {
    Progress(KeggProgress),
    Done(SpeciesKegg),
    Failed(String),
}

#[derive(Debug, Clone)]
pub enum KeggCacheScope {
    Organisms,
    Species(String),
    /// Invalidate the global KEGG modules cache (`modules.json`).
    Modules,
    /// Invalidate the derived organism-group precompute cache (`organism_groups.json`).
    OrganismGroups,
}

/// One entry in the global KEGG modules cache. `compounds` may be empty
/// for modules whose `/get` response lacks a COMPOUND block (e.g. some
/// signature/reaction-only modules); `complete_orgs` may be empty for
/// modules with no COMPLETE block (very new uncurated modules). Both
/// cases are valid KEGG responses and MUST be cached normally — only
/// transport failures skip the cache write.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct KeggModuleEntry {
    /// Human-readable name from `/list/module` (e.g. "Glycolysis
    /// (Embden-Meyerhof pathway), glucose => pyruvate").
    pub name: String,
    /// Compound IDs (`Cxxxxx`) in source order from the COMPOUND block.
    pub compounds: Vec<String>,
    /// Organism codes (3-6 char strings like `hsa`, `ath`) where this
    /// module is fully complete, per the COMPLETE block. Scientific and
    /// common names from KEGG are deliberately dropped — codes alone
    /// join back to `organisms.json` when display names are needed.
    pub complete_orgs: HashSet<String>,
    pub fetched_at: DateTime<Utc>,
}

/// Global cache of all KEGG modules. Stored at `<cache_dir>/modules.json`.
/// Per-entry `fetched_at` (NOT a single file-level timestamp) supports
/// incremental fills across many sessions, like the CID→cpd cache.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct KeggModulesCache {
    pub modules: HashMap<String, KeggModuleEntry>,
}

/// Progress event emitted by `fetch_modules_incremental` per module.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ModuleFetchProgress {
    pub completed: usize,
    pub total: usize,
    pub current_id: String,
    /// Estimated time-to-completion in seconds, derived from a rolling
    /// average of per-completion wall-clock duration. `None` until the
    /// rolling-average buffer has at least 5 samples (warmup window).
    pub eta_secs: Option<u64>,
}

/// Precomputed index from KEGG taxonomy lineage (the semicolon-delimited
/// `KeggOrganism.lineage`, reconstructed from the BRITE `br08601` hierarchy)
/// to organism codes. Three levels
/// are exposed in the UI: Level 1 (e.g. Eukaryotes / Prokaryotes), Level 2
/// (e.g. Animals, Bacteria), Level 3 (e.g. Mammals, Insects).
///
/// `by_level[N - 1]` is `{ group_name -> set of organism codes }`. Stored
/// at `<cache_dir>/organism_groups.json` and rewritten every time the
/// organism list is loaded (cache hit OR fresh fetch); its `fetched_at`
/// MUST match the organisms cache's `fetched_at` so the two stay coherent.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct OrganismGroupIndex {
    pub fetched_at: DateTime<Utc>,
    pub by_level: [HashMap<String, HashSet<String>>; 3],
}

/// One entry in the on-disk CID → KEGG compound cache. `cpd` is `None`
/// when the KEGG `/conv` endpoint confirmed there is no compound mapping
/// for this CID; missing keys in the outer cache map mean "not yet
/// queried".
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct CidCpdEntry {
    pub cpd: Option<String>,
    pub fetched_at: DateTime<Utc>,
}

/// Progress event emitted by `kegg::conv::resolve_cids_to_cpds` per chunk.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConvProgress {
    pub completed_batches: usize,
    pub total_batches: usize,
    pub from_cache: usize,
    pub fetched: usize,
    /// Total CID count for this resolver call; the UI can render
    /// progress as `(from_cache + fetched) / total_inputs`.
    pub total_inputs: usize,
}
