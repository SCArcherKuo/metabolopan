pub mod brunner_munzel;
pub mod export;
pub mod fdr;
pub mod filter;
pub mod run;
pub mod student;
pub mod transforms;
pub mod types;
pub mod welch;

pub use run::{DamConfig, classify_trend, run_dam};
pub use types::{DamFeature, DamMethod, DamProgress, DamResult, FcBasis, Trend};
