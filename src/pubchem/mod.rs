//! PubChem PUG REST client and resolver for InChIKey → CID mapping.
//!
//! Owner: the `pubchem-fetching` capability.

pub mod cache;
pub mod client;
pub mod resolve;
pub mod types;

pub use cache::clear_stale_locks;
pub use client::PubchemClient;
pub use resolve::resolve_inchikeys_to_cids;
pub use types::{InchikeyCidsEntry, PubchemProgress};
