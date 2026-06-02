//! Data types for the PubChem PUG REST client and resolver.

use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};

/// One entry in the on-disk InChIKey → CIDs cache. `cids` is empty (not
/// missing) when PubChem confirmed there is no CID for this InChIKey;
/// that distinguishes "no match" from "not yet queried" (the latter
/// shows up as a missing key in the outer cache map).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct InchikeyCidsEntry {
    pub cids: Vec<String>,
    pub fetched_at: DateTime<Utc>,
}

/// Progress event emitted by `resolve_inchikeys_to_cids` per batch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct PubchemProgress {
    /// Number of network batches completed so far.
    pub completed_batches: usize,
    /// Total network batches that will be issued (excludes cache hits).
    pub total_batches: usize,
    /// Input count served entirely from cache (no network).
    pub from_cache: usize,
    /// Input count fetched fresh from PubChem in this resolver call.
    pub fetched: usize,
    /// Total InChIKey count for this resolver call; the UI can render
    /// progress as `(from_cache + fetched) / total_inputs`.
    pub total_inputs: usize,
}
