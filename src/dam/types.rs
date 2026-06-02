use serde::{Deserialize, Serialize};

pub use crate::dam::fdr::FdrMethod;
use crate::dedup::types::DedupReport;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum DamMethod {
    Welch,
    Student,
    BrunnerMunzel,
}

impl DamMethod {
    pub fn display_name(self) -> &'static str {
        match self {
            DamMethod::Welch => "Welch's t-test",
            DamMethod::Student => "Student's t-test",
            DamMethod::BrunnerMunzel => "Brunner-Munzel",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FcBasis {
    /// log2(mean(num) / mean(den)) on the working (post-normalization, pre-test-transform)
    /// matrix. Parametric tests (Welch / Student) without arcsinh use this.
    Mean,
    /// log2(median(num) / median(den)) on the working matrix. BM uses this.
    Median,
    /// (mean(arcsinh(num)) − mean(arcsinh(den))) / ln(2). Parametric tests with the
    /// `Log transformation` (arcsinh) variance-stabilisation step use this so the
    /// reported FC sign always agrees with the t-statistic sign — even on heavy-tailed
    /// data where Jensen's inequality would make the raw mean ratio diverge from the
    /// arcsinh-mean difference.
    ArcsinhMean,
}

impl FcBasis {
    pub fn label(self) -> &'static str {
        match self {
            FcBasis::Mean => "mean",
            FcBasis::Median => "median",
            FcBasis::ArcsinhMean => "arcsinh-mean",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Trend {
    Up,
    Down,
    NotSignificant,
}

impl Trend {
    pub fn label(self) -> &'static str {
        match self {
            Trend::Up => "up",
            Trend::Down => "down",
            Trend::NotSignificant => "ns",
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DamFeature {
    pub alignment_id: String,
    pub metabolite_name: String,
    pub inchikey: Option<String>,
    pub average_rt_min: Option<f64>,
    pub average_mz: Option<f64>,
    pub formula: Option<String>,
    pub smiles: Option<String>,
    pub numerator_mean: f64,
    pub denominator_mean: f64,
    pub numerator_median: f64,
    pub denominator_median: f64,
    pub fold_change: f64,
    pub log2_fold_change: f64,
    pub fc_basis: FcBasis,
    pub p_value: f64,
    pub p_adjusted: f64,
    pub neg_log10_p_adjusted: f64,
    pub effect_size: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DamProgress {
    pub completed: usize,
    pub total: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DamResult {
    pub method: DamMethod,
    pub numerator: String,
    pub denominator: String,
    pub features: Vec<DamFeature>,
    pub skipped: usize,
    /// FDR correction applied to the per-feature p values. Carried so the
    /// CSV exporter and downstream renderers can surface the choice.
    #[serde(default)]
    pub fdr_method: FdrMethod,
    /// Optional deduplication audit report. `Some(_)` when `run_dam` was
    /// invoked with `dedup_enabled = true`; `None` when dedup was opted
    /// out OR when this `DamResult` was deserialized (the field is not
    /// persisted across the JSON boundary — see `#[serde(skip)]`).
    /// Drives the Stage 2 threshold screen's "Download dedup audit
    /// (CSV)" button visibility.
    #[serde(skip)]
    pub dedup_report: Option<DedupReport>,
}
