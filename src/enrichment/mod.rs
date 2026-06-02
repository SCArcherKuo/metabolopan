//! Stage 3 enrichment analysis: hypergeometric ORA + BY FDR on KEGG
//! entries (pathways in pathway mode, modules in module mode), plus
//! CSV export.
//!
//! See the `enrichment-ora` capability in `openspec/specs/` for the
//! requirement contract.

pub mod export;
pub mod ora;
pub mod types;

pub use crate::dam::fdr::FdrMethod;
pub use export::{export_csv, export_csv_with_mode};
pub use ora::run_ora;
pub use types::{EnrichmentDirection, EnrichmentResult, EnrichmentRow};
