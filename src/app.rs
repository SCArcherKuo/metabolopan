use std::path::PathBuf;
use std::sync::mpsc;
use std::time::Duration;

use chrono::{DateTime, Utc};
use eframe::{App as EframeApp, Frame};
use egui::{CentralPanel, Context, TopBottomPanel};
use tokio::runtime::Runtime;
use tracing::{error, info};

use crate::dam::fdr::FdrMethod;
use crate::dam::{DamMethod, DamResult};
use crate::data::{GroupMapping, IonMode, IonModeTable};
use crate::enrichment::{EnrichmentDirection, EnrichmentResult};
use crate::kegg::{
    ConvProgress, KeggClient, KeggEvent, KeggModulesCache, KeggOrganism, ModuleFetchProgress,
    OrganismsCache, SpeciesKegg,
};
use crate::logging::LogStore;
use crate::normalize::{NormalizationMethod, PqnReference};
use crate::pubchem::PubchemProgress;
use crate::theme;
use crate::ui::organism_group_selector::OrganismGroupSelectorState;
use crate::ui::species_selector::SpeciesSelectorState;
use crate::ui::{
    bottom_panel, initializing, log_pane, settings_modals, stage1_input, stage2_running,
    stage2_setup, stage2_threshold, stage3_result, stage3_running, stage3_setup, stepper,
};
use serde::{Deserialize, Serialize};
use std::collections::HashSet;

// Note: `AppState` cannot derive `Debug` because `egui::TextureHandle` (held by
// `Stage2DamThreshold`) does not implement `Debug`. Logging that needs to inspect
// state should do so per-variant.

/// Volcano render result: RGBA buffer plus the canvas dimensions it was drawn at.
/// Carrying the dims with the buffer lets the UI thread upload a texture of the
/// correct size even if the user changed the export inputs mid-flight.
pub type VolcanoRender = (Vec<u8>, u32, u32);

/// Top-level analysis mode picked at Stage 1. Pathway is the historical
/// default; Module enables the new KEGG-module ORA flow.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum AnalysisMode {
    Pathway,
    Module,
}

impl AnalysisMode {
    /// Singular lowercase noun for the entity this mode enriches over.
    /// Threaded into the dot-plot renderer (`DotplotOpts.entry_label`) so
    /// the empty-state placeholder and Y-axis title read correctly in
    /// Module mode — Pathway mode renderers used to hardcode "pathways".
    pub fn entry_label_singular(self) -> &'static str {
        match self {
            Self::Pathway => "pathway",
            Self::Module => "module",
        }
    }
}

/// Mode-specific analysis payload threaded through Stage 2 → Stage 3.
/// Replaces the prior `species_kegg: SpeciesKegg` field in downstream
/// states so module-mode runs can carry their own data without breaking
/// the pathway-mode contract.
#[allow(clippy::large_enum_variant)]
#[derive(Clone)]
pub enum AnalysisPayload {
    Pathway {
        species_kegg: SpeciesKegg,
    },
    Module {
        modules_pack: KeggModulesCache,
        group_level: u8,
        group_name: String,
        group_org_codes: HashSet<String>,
        min_group_overlap: usize,
    },
}

impl AnalysisPayload {
    pub fn mode(&self) -> AnalysisMode {
        match self {
            Self::Pathway { .. } => AnalysisMode::Pathway,
            Self::Module { .. } => AnalysisMode::Module,
        }
    }

    /// Convenience accessor for code paths that ONLY apply in pathway
    /// mode (e.g. species code shown in result panels). Returns `None`
    /// in module mode.
    pub fn pathway_species(&self) -> Option<&SpeciesKegg> {
        match self {
            Self::Pathway { species_kegg } => Some(species_kegg),
            Self::Module { .. } => None,
        }
    }
}

/// Central store for every user-tunable parameter across all stages.
/// Introduced by the `refactor-session-settings` change so that
/// settings survive transitions without each transition function
/// hand-copying field subsets. See the `app-shell` capability spec
/// for the normative field list and default values, and the
/// `refactor-session-settings` change's design notes (decision D11)
/// for the table of reset/clear API surfaces.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SessionSettings {
    // ── Stage 1 / mode-related ──
    pub analysis_mode: AnalysisMode,
    pub kegg_species: Option<String>,
    pub organism_group_level: Option<u8>,
    pub organism_group: Option<String>,
    pub min_group_overlap: usize,

    // ── Stage 2 DAM setup ──
    pub numerator: Option<String>,
    pub denominator: Option<String>,
    pub dam_method: DamMethod,
    pub drop_unknown: bool,
    pub dedup_enabled: bool,
    pub normalization: NormalizationMethod,
    pub metadata_column: Option<String>,
    pub pqn_reference: PqnReference,
    pub pqn_reference_group: Option<String>,
    /// Welch / Student pre-test arcsinh toggle. `true` = apply the project's
    /// "generalised log" (asinh) before the t-test; `false` = pass the working
    /// matrix to the t-test directly. BM ignores this field (rank-invariant
    /// under monotone arcsinh). Added in SCHEMA_VERSION 2 by
    /// add-log-transform-and-scaling.
    #[serde(default = "default_log_transform")]
    pub log_transform: bool,
    pub dam_fdr_method: FdrMethod,

    // ── Stage 2 DAM result thresholds + export size ──
    pub fc_threshold: f64,
    pub fdr_threshold: f64,
    pub delta_threshold: f64,
    pub stage2_export_width_in: f64,
    pub stage2_export_height_in: f64,
    pub stage2_export_dpi: u32,

    // ── Stage 3 setup ──
    pub direction: EnrichmentDirection,
    pub top_n: usize,
    pub enrichment_fdr_threshold: f64,
    pub min_hit_count: usize,
    /// Pre-FDR entry-size filter for Stage 3 ORA. Entries with `m_p`
    /// (universe-restricted compound count) below this threshold are
    /// dropped before FDR — they don't enter the FDR family `m`.
    /// Added in SCHEMA_VERSION 3 by add-min-entry-size-filter.
    #[serde(default = "default_min_entry_size")]
    pub min_entry_size: usize,
    pub enrichment_fdr_method: FdrMethod,

    // ── Stage 3 result export size ──
    pub stage3_export_width_in: f64,
    pub stage3_export_height_in: f64,
    pub stage3_export_dpi: u32,
}

/// Serde-default helper for the `log_transform` field. Returning `true` so
/// snapshots saved before SCHEMA_VERSION 2 — if any bypassed the strict version
/// gate — would fall back to the pre-change pipeline behaviour. The gate
/// normally rejects v1 outright; this is defensive.
fn default_log_transform() -> bool {
    true
}

/// Serde-default helper for the `min_entry_size` field. Returning `1` to
/// match `SessionSettings::default()`. Defensive against snapshots saved
/// before SCHEMA_VERSION 3 that bypass the strict version gate.
fn default_min_entry_size() -> usize {
    1
}

/// Initial Stage 3 dot-plot export height (inches), sized to the rows the
/// plot will actually show: `min(top_n, displayed_rows)` (at least 1) at the
/// per-row `0.3 in` rhythm plus a `1.0 in` base, clamped to `[2.0, 40.0]`.
/// Sizing to the displayed-row count (rather than raw `top_n`) stops a sparse
/// result — far fewer significant entries than `top_n` — from rendering in a
/// tall band of whitespace. Only sets the initial value at result-entry; the
/// result-screen Height field remains a user override.
///
/// `pub(crate)` so the Stage 3 result UI can re-fit the height on each
/// "Re-draw dot plot" against the live displayed-row count (the display
/// filters — Enrichment FDR threshold / Min hit count / Top N — change the
/// row count without a re-run, so the run-entry autosize alone goes stale).
pub(crate) fn stage3_autosize_height_in(top_n: usize, displayed_rows: usize) -> f64 {
    let effective = top_n.min(displayed_rows).max(1);
    // Lower clamp 2.0 (was 1.5): the dot plot now scales fonts/elements by the
    // fixed width (`common::design_scale_by_width`), so the legend renders at
    // full size even on a short canvas. 2.0 in (600 px @ 300 dpi) guarantees the
    // full-height colorbar + Hits block clears the canvas on 1–3-row results
    // that would otherwise size to 1.3–1.9 in and clip the last reference dot.
    ((effective as f64) * 0.3 + 1.0).clamp(2.0, 40.0)
}

impl Default for SessionSettings {
    fn default() -> Self {
        Self {
            // Stage 1 / mode
            analysis_mode: AnalysisMode::Pathway,
            kegg_species: None,
            organism_group_level: None,
            organism_group: None,
            min_group_overlap: 1,

            // Stage 2 setup
            numerator: None,
            denominator: None,
            dam_method: DamMethod::Student,
            drop_unknown: true,
            dedup_enabled: true,
            normalization: NormalizationMethod::None,
            metadata_column: None,
            pqn_reference: PqnReference::AllSamples,
            pqn_reference_group: None,
            log_transform: default_log_transform(),
            dam_fdr_method: FdrMethod::BenjaminiHochberg,

            // Stage 2 result
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            stage2_export_width_in: 3.5,
            stage2_export_height_in: 2.2,
            stage2_export_dpi: 300,

            // Stage 3 setup
            direction: EnrichmentDirection::Both,
            top_n: 20,
            enrichment_fdr_threshold: 0.05,
            min_hit_count: 1,
            min_entry_size: default_min_entry_size(),
            enrichment_fdr_method: FdrMethod::BenjaminiHochberg,

            // Stage 3 result (default height = top_n * 0.3 + 1.0 = 7.0)
            stage3_export_width_in: 3.5,
            stage3_export_height_in: 7.0,
            stage3_export_dpi: 300,
        }
    }
}

impl SessionSettings {
    /// Reset surface for "Back to DAM Setup" on Stage 2 DAM result. After
    /// `reorder-gui-and-move-mode-to-stage3` (Phase 2), this body is a
    /// NO-OP: every Stage 2 settings field is preserved across the
    /// transition. The method name remains so transition call sites stay
    /// spec-anchored; a future change can adjust the body without
    /// touching the call surface.
    pub fn reset_stage2_choices_on_change_comparison(&mut self) {}

    /// Reset surface for "< Back to DAM Result" on Stage 3 setup. After
    /// the smoke-test feedback during `reorder-gui-and-move-mode-to-stage3`
    /// (post-Phase 5), this body is a NO-OP: every settings field is
    /// preserved when the user navigates back to Stage 2 result, so they
    /// can revisit thresholds + export size without losing their picks.
    /// The method name remains so transition call sites stay
    /// spec-anchored.
    pub fn reset_for_back_to_stage2_threshold(&mut self) {}

    /// Reset surface for "Continue to Enrichment" on Stage 2 DAM result.
    /// After the smoke-test feedback (post-Phase 5), this body is a
    /// NO-OP: every Stage 3 settings field is preserved across Continue
    /// transitions, so the user's prior direction / top_n / FDR / min hit
    /// count carry over into the next Stage 3 setup session.
    pub fn reset_stage3_on_continue_to_enrichment(&mut self) {}

    /// Reset surface for "< Back to Input" on Stage 2 DAM setup. After
    /// the smoke-test feedback (post-Phase 5), this body is a NO-OP:
    /// every Stage 2 and Stage 3 settings field is preserved across the
    /// Back-to-Stage-1 transition. If the user later re-picks files at
    /// Stage 1 such that the preserved numerator / denominator groups
    /// no longer exist, the Stage 2 setup dropdown gates pick that up
    /// (it shows "— pick one —" until the user re-selects a valid
    /// group).
    pub fn reset_for_back_to_stage1(&mut self) {}

    /// Reset surface for the Pathway ↔ Module radio toggle on Stage 3
    /// setup. After `reorder-gui-and-move-mode-to-stage3` (Phase 2), the
    /// body is reduced to `self.analysis_mode = new_mode;` only — the
    /// inactive mode's selection (`kegg_species` for Pathway,
    /// `organism_group_level` / `organism_group` for Module) is
    /// PRESERVED so the user can toggle between modes without losing
    /// either side's prior selection. The paired `cache.clear_for_mode_switch`
    /// is now also a no-op (see `SessionCache` below); the two methods
    /// stay paired at call sites for the same spec-anchoring reason.
    pub fn reset_kegg_selection_for_mode_switch(&mut self, new_mode: AnalysisMode) {
        self.analysis_mode = new_mode;
    }
}

/// Loaded raw inputs that downstream stages read but do not modify.
/// `IonModeTable` slot count is 1 (single-mode) or 2 (dual-mode); the
/// `IonModeTables` newtype invariant (length 1–2, no duplicate mode,
/// Positive-first ordering) is enforced transiently inside
/// `promote_to_stage2`, NOT stored here — Stage 1 mid-pick states would
/// not satisfy that invariant yet.
#[derive(Debug, Default)]
pub struct SessionInputs {
    pub ion_tables: Vec<IonModeTable>,
    pub mapping: Option<GroupMapping>,
    pub csv_path: Option<PathBuf>,
}

/// Per-analysis fetched KEGG data, held in memory after the fetch
/// completes. Both pathway-mode (`species_kegg`) and module-mode
/// (`modules_pack` + `group_org_codes`) caches can coexist; today's
/// behavior clears the other side when the user toggles
/// `AnalysisMode` (via `clear_for_mode_switch`). Change #2 will stop
/// calling that method, allowing parallel retention.
///
/// `App::organisms` (the eager-loaded KEGG organism roster) is NOT in
/// this struct — it is session-immutable lookup data, loaded once at
/// startup, and lives directly on `App`. See the `app-shell` capability
/// spec (decision D10 / preamble) for the rule.
#[derive(Debug, Default, Clone)]
pub struct SessionCache {
    pub species_kegg: Option<SpeciesKegg>,
    pub modules_pack: Option<KeggModulesCache>,
    pub group_org_codes: Option<HashSet<String>>,
}

impl SessionCache {
    /// After `reorder-gui-and-move-mode-to-stage3` (Phase 2), this body
    /// is a NO-OP: Pathway and Module caches coexist for the lifetime of
    /// the session so the user can toggle modes without losing either
    /// mode's fetched data. The method name remains so transition call
    /// sites stay spec-anchored.
    pub fn clear_for_mode_switch(&mut self, new_mode: AnalysisMode) {
        let _ = new_mode;
    }
}

// `AppState` variants carry large analysis payloads (`MetabolomicsTable`,
// `DamResult`, `EnrichmentResult`, hash sets / maps, texture handles).
// Variant size disparity is inherent and irrelevant — there is exactly
// one `AppState` value alive in the process at a time, and transitions
// move ownership rather than allocating new memory for each variant.
#[allow(clippy::large_enum_variant)]
pub enum AppState {
    /// Startup splash. Eagerly loading the KEGG organism list (cache or
    /// fresh REST) before the rest of the UI is reachable. On success,
    /// transitions to `Stage1Input` with `app.organisms` populated. On
    /// failure, surfaces a Retry button (and "Use cached organisms"
    /// when a stale cache exists in `fallback_cache`).
    ///
    /// **Refactor note (`refactor-session-settings`)**: as of Phase 3,
    /// every `AppState` variant carries ONLY runtime artifacts — mpsc
    /// receivers, progress accumulators, computed outputs, textures, and
    /// screen-local error strings. Every user-tunable parameter lives on
    /// `App::settings`; every loaded input lives on `App::inputs`; every
    /// per-analysis fetched cache lives on `App::cache`. See the
    /// `app-shell` capability spec for the normative contract.
    Initializing {
        load_rx: mpsc::Receiver<Result<Vec<KeggOrganism>, String>>,
        /// Pre-loaded stale cache, surfaced as a fallback on failure.
        fallback_cache: Option<OrganismsCache>,
        /// Latest error from a previous load attempt, if any.
        last_error: Option<String>,
    },
    /// Stage 1 — file pickers + mode toggle + species/group selector.
    /// The variant carries only Stage-1-screen-local UI state for the
    /// file-pick radios; everything else (analysis mode, selected
    /// species/group, loaded `.txt`/`.csv`, KEGG cache) lives on
    /// `App::settings` / `App::inputs` / `App::cache`.
    Stage1Input {
        /// Slot #1 ionization-mode radio. `None` on fresh entry; user
        /// must pick before Start. Mirrors the file picker's per-slot UI.
        slot1_mode: Option<IonMode>,
        /// Whether the user has clicked "+ Add second ionization mode" to
        /// reveal slot #2. Stays `false` until clicked; an `×` removal
        /// sets it back to `false`.
        slot2_revealed: bool,
        /// Slot #2 ionization-mode radio. `None` until set by the user.
        slot2_mode: Option<IonMode>,
        error: Option<String>,
    },
    /// Stage 2 — DAM setup screen. Every user-tunable parameter
    /// (numerator/denominator/method/dedup/drop_unknown/normalization
    /// fields, FDR method) lives on `App::settings`; only the
    /// screen-local error string lives on the variant.
    Stage2DamSetup { error: Option<String> },
    /// Transient: DAM tasks running. Runtime artifacts only — every
    /// parameter read by the worker comes from `App::settings`,
    /// `App::inputs`, `App::cache` at spawn time.
    Stage2DamRunning {
        /// Per-mode result channel. Each spawned worker emits one
        /// message `(mode_idx, Result<DamResult, String>)`.
        result_rx: mpsc::Receiver<(usize, Result<DamResult, String>)>,
        /// One progress channel per mode (index matches
        /// `App::inputs.ion_tables`).
        progress_rxs: Vec<mpsc::Receiver<crate::dam::DamProgress>>,
        mode_completed: Vec<usize>,
        mode_total: Vec<usize>,
        /// Accumulator for per-mode terminal results.
        dam_results_accum: Vec<Option<Result<DamResult, String>>>,
    },
    /// Stage 2 — DAM result (threshold) screen. Carries the computed
    /// `dam_results` plus volcano-rendering runtime. Thresholds and
    /// export size are on `App::settings`.
    Stage2DamThreshold {
        dam_results: Vec<DamResult>,
        /// Which ion-mode tab is currently selected in the volcano
        /// area. Default = first loaded mode (`inputs.ion_tables[0].mode`).
        active_volcano_tab: IonMode,
        /// Per-mode volcano texture cache. Index matches
        /// `App::inputs.ion_tables`. Slider/threshold changes
        /// invalidate ALL entries.
        volcano_textures: Vec<Option<egui::TextureHandle>>,
        rendering: bool,
        render_rx: Option<mpsc::Receiver<(usize, Result<VolcanoRender, String>)>>,
    },
    /// Stage 3 — enrichment setup screen. Carries `dam_results` (to
    /// thread DAM output forward into the enrichment Run), the
    /// screen-local error string, and the OPTIONAL inline KEGG-fetch
    /// progress state for whichever mode-specific fetch is currently in
    /// flight. Direction / Top N / Min hit count / FDR method etc. live
    /// on `App::settings`.
    ///
    /// Why the fetch state lives here: per `reorder-gui-and-move-mode-to-stage3`
    /// D1+D2, the Mode toggle + species/group selector live on this
    /// screen, and a fetch triggered by the user picking a
    /// species/Group renders progress INLINE on this screen rather than
    /// transitioning to a separate `KeggFetching` / `ModulesFetching`
    /// variant. The two flat optionals (instead of one enum) allow
    /// either mode's fetch to coexist in flight; toggling mode while a
    /// fetch is in flight does not cancel it (D6).
    Stage3EnrichSetup {
        dam_results: Vec<DamResult>,
        error: Option<String>,
        /// `Some(_)` while a per-species pathway fetch is streaming. The
        /// species code being fetched is read from
        /// `App::settings.kegg_species`. Cleared back to `None` by the
        /// terminal-event handler on Done / Failed.
        kegg_fetch: Option<KeggFetchInFlight>,
        /// `Some(_)` while a module-mode bulk fetch is streaming. The
        /// Group/level being fetched is on `App::settings`. Cleared back
        /// to `None` on terminal event.
        modules_fetch: Option<ModulesFetchInFlight>,
    },
    /// Transient: Stage 3 orchestrator running (PubChem → KEGG conv →
    /// ORA). Carries `dam_results` (still needed for downstream display
    /// and for back-to-threshold) plus the 3-phase progress runtime.
    Stage3EnrichRunning {
        dam_results: Vec<DamResult>,
        phase: Stage3Phase,
        pubchem_progress_rx: mpsc::Receiver<PubchemProgress>,
        kegg_conv_progress_rx: mpsc::Receiver<ConvProgress>,
        result_rx: mpsc::Receiver<Result<Stage3RunOutput, String>>,
        pubchem_completed: usize,
        pubchem_total: usize,
        kegg_conv_completed: usize,
        kegg_conv_total: usize,
    },
    /// Stage 3 — enrichment result screen. Carries the computed output
    /// (`enrichment_result`, `feature_to_cpds`, `mapped_universe`,
    /// time-span tuples, `module_retention`, `dual_mode_breakdown`) and
    /// dot-plot rendering runtime. Export size is on `App::settings`.
    Stage3EnrichResult {
        dam_results: Vec<DamResult>,
        /// Track E: populated in module mode from
        /// `Stage3RunOutput.module_retention`. `None` in pathway mode
        /// (the result panel reads the species' own `code` / `fetched_at`
        /// from `App::cache.species_kegg` instead).
        module_retention: Option<crate::stage3::ModuleRetention>,
        enrichment_result: EnrichmentResult,
        mapped_universe: std::collections::HashSet<String>,
        feature_to_cpds: std::collections::HashMap<String, std::collections::HashSet<String>>,
        pubchem_time_span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>,
        kegg_conv_time_span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>,
        /// Dual-mode partition counts. `None` in single-mode runs.
        dual_mode_breakdown: Option<crate::stage3::DualModeBreakdown>,
        /// Provenance funnel counts (InChIKey → CID → KEGG-cpd, detected +
        /// foreground) surfaced by the Data tab (`add-bottom-panel-data-tab`).
        funnel: Stage3Funnel,
        dotplot_tex: Option<egui::TextureHandle>,
        rendering: bool,
        render_rx: Option<mpsc::Receiver<Result<DotplotRender, String>>>,
        refresh_state: RefreshState,
        /// Variant-internal "Start a new analysis?" confirmation flag.
        /// `true` while the loss-warning modal is open; set by the
        /// `Start a new analysis` button, cleared on Cancel, and consumed
        /// by `App::start_new_round` on Confirm. NOT part of the App-level
        /// modal mutual-exclusion family (it mirrors `refresh_state`, a
        /// central-panel confirm). Ephemeral UI state — deliberately NOT
        /// surfaced in the bug-report `app_state.txt` snapshot.
        confirming_new_round: bool,
        /// `true` once the user has hand-edited the `Height (in)` export
        /// field on this result screen. While `false`, every "Re-draw dot
        /// plot" re-fits the height to the live displayed-row count (so
        /// loosening the FDR threshold and redrawing grows the canvas
        /// instead of cramming rows into a stale autosize); once `true`,
        /// the user's height is honored verbatim. Reset to `false` on every
        /// fresh Run/Re-run (a new result state is built per run). Ephemeral
        /// UI state — not persisted, not in the bug-report snapshot.
        height_user_overridden: bool,
    },
}

/// Drain available KEGG progress events into the in-flight struct;
/// return the terminal event (Done / Failed) if one arrived this frame.
/// Used by `update()` while on `Stage3EnrichSetup`.
fn drain_kegg_progress(fetch: &mut KeggFetchInFlight) -> Option<KeggEvent> {
    loop {
        match fetch.progress_rx.try_recv() {
            Ok(KeggEvent::Progress(p)) => {
                fetch.completed = p.completed;
                fetch.total = p.total;
                fetch.current_pathway = p.current_pathway.clone();
            }
            Ok(other) => return Some(other),
            Err(mpsc::TryRecvError::Empty) | Err(mpsc::TryRecvError::Disconnected) => return None,
        }
    }
}

/// Drain available module-fetch progress events into the in-flight
/// struct; return the terminal event (Done / Failed) if one arrived this
/// frame.
fn drain_modules_progress(fetch: &mut ModulesFetchInFlight) -> Option<ModulesFetchEvent> {
    loop {
        match fetch.progress_rx.try_recv() {
            Ok(ModulesFetchEvent::Progress(p)) => {
                fetch.completed = p.completed;
                fetch.total = p.total;
                fetch.current_id = p.current_id.clone();
                fetch.eta_secs = p.eta_secs;
            }
            Ok(other) => return Some(other),
            Err(mpsc::TryRecvError::Empty) | Err(mpsc::TryRecvError::Disconnected) => return None,
        }
    }
}

/// Inline KEGG species-pathway fetch state held by `Stage3EnrichSetup`
/// while a fetch is streaming. Replaces the deleted `AppState::KeggFetching`
/// variant — instead of leaving the setup screen for a dedicated fetching
/// screen, the user stays on Stage 3 setup and a progress strip renders
/// inline. The terminal-event handler writes the fetched `SpeciesKegg`
/// into `App::cache.species_kegg` and clears this field to `None`.
pub struct KeggFetchInFlight {
    pub progress_rx: mpsc::Receiver<KeggEvent>,
    pub completed: usize,
    pub total: usize,
    pub current_pathway: String,
}

/// Inline module-mode bulk fetch state held by `Stage3EnrichSetup` while a
/// fetch is streaming. Module parallel of `KeggFetchInFlight`. On terminal
/// event the handler writes the fetched `KeggModulesCache` into
/// `App::cache.modules_pack` and clears this field to `None`.
pub struct ModulesFetchInFlight {
    pub progress_rx: mpsc::Receiver<ModulesFetchEvent>,
    pub completed: usize,
    pub total: usize,
    pub current_id: String,
    pub eta_secs: Option<u64>,
}

/// Event surfaced by the modules-fetch orchestrator. Mirrors the
/// `KeggEvent` shape for the pathway-mode flow.
#[derive(Debug)]
pub enum ModulesFetchEvent {
    Progress(ModuleFetchProgress),
    Done(KeggModulesCache),
    Failed(String),
}

/// Which phase of the Stage 3 Run is currently in flight.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stage3Phase {
    PubChem,
    KeggConv,
    Ora,
}

/// State of the refresh subsystem on the Stage 3 result screen.
#[derive(Debug, Default)]
pub enum RefreshState {
    #[default]
    Idle,
    ConfirmingPubchem,
    RefreshingPubchem {
        progress_rx: mpsc::Receiver<PubchemProgress>,
        result_rx: mpsc::Receiver<Result<Stage3RunOutput, String>>,
        completed: usize,
        total: usize,
    },
    ConfirmingKegg,
    RefreshingKegg {
        progress_rx: mpsc::Receiver<ConvProgress>,
        result_rx: mpsc::Receiver<Result<Stage3RunOutput, String>>,
        completed: usize,
        total: usize,
    },
}

/// Payload returned by the Stage 3 orchestrator task and consumed by
/// `handle_stage3_terminal` to populate `Stage3EnrichResult`.
pub struct Stage3RunOutput {
    pub enrichment_result: EnrichmentResult,
    pub mapped_universe: std::collections::HashSet<String>,
    pub feature_to_cpds: std::collections::HashMap<String, std::collections::HashSet<String>>,
    pub pubchem_time_span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>,
    pub kegg_conv_time_span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>,
    /// Populated in module mode only; `None` in pathway mode. Drives the
    /// Stage 3 result panel's module-mode "Data sources" copy (Group,
    /// Level, min_group_overlap, retention summary, time span).
    pub module_retention: Option<crate::stage3::ModuleRetention>,
    /// Populated in dual-mode (`dam_results.len() == 2`) runs; `None` for
    /// single-mode. Drives the dual-mode breakdown block in the Stage 3
    /// result panel surfacing per-mode partition of `|N|` and `|K|` plus
    /// the conflict-excluded count.
    pub dual_mode_breakdown: Option<crate::stage3::DualModeBreakdown>,
    /// Provenance funnel counts for the Data tab (`add-bottom-panel-data-tab`),
    /// computed as by-products of the universe/foreground construction — no new
    /// network calls. `detected_*` describe the full measurable-metabolome path
    /// (all DAM features); `foreground_*` describe the significant subset (K).
    /// `detected_in_entries` is the count of foreground cpds appearing in ≥ 1
    /// tested entry. Invariants (asserted in the orchestrator):
    /// `detected_inchikeys >= foreground_inchikeys`,
    /// `detected_cids >= foreground_cids`, `universe_size >= K`.
    pub funnel: Stage3Funnel,
}

impl Stage3RunOutput {
    /// Consume the orchestrator output into a fully-populated
    /// `AppState::Stage3EnrichResult` for the given `dam_results`. The five
    /// result-screen runtime fields start fresh (`dotplot_tex: None`,
    /// `rendering: false`, `render_rx: None`, `refresh_state: Idle`,
    /// `confirming_new_round: false`). Replaces the manual field-by-field splat
    /// at the `handle_stage3_terminal` success arm.
    fn into_result_state(self, dam_results: Vec<DamResult>) -> AppState {
        AppState::Stage3EnrichResult {
            dam_results,
            module_retention: self.module_retention,
            enrichment_result: self.enrichment_result,
            mapped_universe: self.mapped_universe,
            feature_to_cpds: self.feature_to_cpds,
            pubchem_time_span: self.pubchem_time_span,
            kegg_conv_time_span: self.kegg_conv_time_span,
            dual_mode_breakdown: self.dual_mode_breakdown,
            funnel: self.funnel,
            dotplot_tex: None,
            rendering: false,
            render_rx: None,
            refresh_state: RefreshState::Idle,
            confirming_new_round: false,
            // Fresh per-run state: the run-entry autosize (set in
            // `handle_stage3_terminal`) is authoritative until the user
            // hand-edits the Height field on this screen.
            height_user_overridden: false,
        }
    }
}

/// Read-only provenance funnel counts surfaced on `Stage3RunOutput.funnel`.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct Stage3Funnel {
    /// Distinct InChIKeys resolved in Phase 1 (deduped union of all DAM
    /// features' InChIKeys).
    pub detected_inchikeys: usize,
    /// Distinct PubChem CIDs collected across all resolved InChIKeys.
    pub detected_cids: usize,
    /// Distinct InChIKeys among significant DAM features (active direction).
    pub foreground_inchikeys: usize,
    /// Distinct PubChem CIDs reachable from the foreground InChIKeys.
    pub foreground_cids: usize,
    /// Foreground (`K`) cpds appearing in ≥ 1 tested entry's compound set.
    pub detected_in_entries: usize,
}

/// Dot plot render channel payload: RGBA buffer plus its dimensions.
pub type DotplotRender = (Vec<u8>, u32, u32);

/// Helper for transitions returning to `Stage1Input` from a downstream
/// state: derive the slot-radio fields from a loaded `ion_tables` so the
/// user's mode choices are preserved across navigation.
pub fn slot_fields_from(ion_tables: &[IonModeTable]) -> (Option<IonMode>, bool, Option<IonMode>) {
    (
        ion_tables.first().map(|it| it.mode),
        ion_tables.len() >= 2,
        ion_tables.get(1).map(|it| it.mode),
    )
}

impl Default for AppState {
    fn default() -> Self {
        AppState::Stage1Input {
            slot1_mode: None,
            slot2_revealed: false,
            slot2_mode: None,
            error: None,
        }
    }
}

/// Holds the organism list once loaded (lazy).
#[derive(Debug, Default)]
pub struct OrganismsLoad {
    pub state: OrganismsLoadState,
}

#[derive(Debug, Default)]
pub enum OrganismsLoadState {
    #[default]
    Idle,
    Loading {
        rx: mpsc::Receiver<OrganismsLoadResult>,
    },
    Loaded(Vec<KeggOrganism>),
    Failed(String),
}

#[derive(Debug)]
pub enum OrganismsLoadResult {
    Ok(Vec<KeggOrganism>),
    Err(String),
}

/// State of the `[Download bug report…]` modal flow. Closed until the
/// user clicks the log-pane button; Confirming while the privacy-list
/// modal is up; Saving while the background thread is assembling the
/// zip + writing it. Once Saving completes (success or error), returns
/// to Closed.
pub enum BundleModalState {
    Closed,
    Confirming,
    Saving {
        /// On success the worker sends `Ok((path, zip_size_bytes))` so the
        /// completion INFO event can report the byte size — spec scenario
        /// "User completes the export" mandates `(target path, zip byte size)`.
        result_rx: mpsc::Receiver<Result<(PathBuf, u64), String>>,
    },
}

/// State of the `[Save settings…]` modal flow. Closed until the user
/// clicks the log-pane button; Confirming while the pre-save summary
/// modal is up. The OS save-file dialog blocks the UI thread so no
/// dedicated `Saving` variant is needed (synchronous write).
#[derive(Debug, Default)]
pub enum SettingsSaveModalState {
    #[default]
    Closed,
    Confirming,
}

/// State of the `[Load settings…]` modal flow. The Confirming variant
/// is populated AFTER `session_io::load_from_path` has parsed the file
/// — the modal body needs the snapshot contents (saved_at, hash
/// mismatches, validation resets) to render.
///
/// Same rationale as `AppState` for the `large_enum_variant` allow:
/// only one value is alive at a time and transitions move ownership.
#[allow(clippy::large_enum_variant)]
#[derive(Debug, Default)]
pub enum SettingsLoadModalState {
    #[default]
    Closed,
    Confirming {
        snapshot: crate::session_io::Snapshot,
        mismatches: Vec<crate::session_io::InputHashMismatch>,
        resets: crate::session_io::ValidationResets,
        path: PathBuf,
    },
}

/// Shared error-toast modal surfaced by both save and load failures
/// (and any future capability that wants a one-line user-readable
/// error). CENTER_CENTER egui Window with a red title and a single
/// `OK` button.
#[derive(Debug, Default, Clone)]
pub enum ErrorModalState {
    #[default]
    Closed,
    Open {
        title: String,
        message: String,
    },
}

/// Returns the default save-dialog filename `bug-report-YYYY-MM-DD_HHMMSS.zip`
/// using local wall-clock time. Local time (not UTC) matches what the user
/// sees on their desktop when they look at the saved file.
pub fn default_bundle_filename() -> String {
    chrono::Local::now()
        .format("bug-report-%Y-%m-%d_%H%M%S.zip")
        .to_string()
}

pub struct App {
    // ── The four sibling fields introduced by the `refactor-session-settings`
    // change. See the `app-shell` capability spec preamble for the
    // normative contract. As of the Phase 2 commit of that refactor, these
    // fields are populated alongside the existing per-variant fields on
    // `AppState`; Phase 3 will slim `AppState` and migrate every callsite to
    // read from these.
    pub state: AppState,
    /// Every user-tunable parameter across all stages. Single source of
    /// truth for defaults (`SessionSettings::default()`); transitions
    /// mutate it via named reset APIs (`reset_*`), never inline bundles.
    pub settings: SessionSettings,
    /// Loaded raw inputs (`IonModeTable`s, `GroupMapping`, `csv_path`).
    /// Survives every state transition; never re-built from
    /// `AppState` variants.
    pub inputs: SessionInputs,
    /// Per-analysis fetched KEGG data (pathway-mode + module-mode
    /// caches can coexist; cleared on user-driven mode toggle via
    /// `SessionCache::clear_for_mode_switch`). Does NOT hold
    /// `App::organisms` — that is session-immutable lookup data.
    pub cache: SessionCache,

    // ── Session-immutable infrastructure fields (per D10 / spec preamble).
    // Neither settings nor input data nor analysis cache — they are
    // session-immutable lookup data and UI control state that must NOT be
    // folded into the four sibling fields above.
    pub log: LogStore,
    pub log_ui: log_pane::LogPaneState,
    pub rt: Runtime,
    pub kegg: KeggClient,
    /// Eager-loaded KEGG organism roster. Session-immutable; loaded once
    /// at startup from `organisms.json` cache or via `/list/organism`.
    /// Does NOT belong in `SessionCache` (does not vary per-analysis).
    pub organisms: OrganismsLoad,
    /// Pathway-mode species selector UI state (filter text, scroll,
    /// open state). UI ephemera, neither settings nor cache.
    pub species_selector: SpeciesSelectorState,
    /// Module-mode organism-group selector UI state. UI ephemera.
    pub organism_group_selector: OrganismGroupSelectorState,
    /// Absolute path to this process's session log file under
    /// `<bin>/data/logs/`. `None` when the file sink failed to
    /// initialise at startup; the bug-report bundle exporter renders
    /// `logs.txt` as a one-line stub in that case.
    pub session_log_path: Option<PathBuf>,
    /// `RUST_LOG` directive captured at startup. Surfaced in
    /// `env.txt` inside the bug-report bundle.
    pub rust_log_directive: String,
    /// Value of `KEGG_CACHE_DIR` at startup, if it was set.
    /// Surfaced in `env.txt` inside the bug-report bundle.
    pub kegg_cache_dir_env: Option<String>,
    /// Bug-report download modal state machine.
    pub bundle_modal: BundleModalState,
    /// Save-settings modal state machine.
    pub settings_save_modal: SettingsSaveModalState,
    /// Load-settings modal state machine.
    pub settings_load_modal: SettingsLoadModalState,
    /// Shared error toast for save / load failures.
    pub error_modal: ErrorModalState,
    /// Stepper step icons, decoded once and uploaded to GPU textures on the
    /// first stepper render (texture upload needs a live `egui::Context`,
    /// unavailable at `App::new`). UI plumbing analogous to the volcano /
    /// dot-plot `TextureHandle`s, not a four-sibling contract field (see
    /// `add-stepper-step-icons` design D4).
    pub stepper_icons: Option<crate::ui::stepper::StepperIcons>,
    /// `rat_face.png` texture, lazily uploaded on the first render of the
    /// bug-report confirm modal title bar.
    pub rat_face_tex: Option<egui::TextureHandle>,
    /// Whether the easter egg image popup window is currently open.
    pub show_rat_easter_egg: bool,
}

impl App {
    pub fn new(
        log: LogStore,
        filter_directive: String,
        rt: Runtime,
        session_log_path: Option<PathBuf>,
    ) -> Self {
        // Clear stale cache locks from any crashed prior process. These
        // are advisory locks for PubChem and KEGG conv cache writes; an
        // orphaned lock would otherwise block this session's writes for
        // up to 30 s per write. Errors are logged but non-fatal — a
        // permission-denied cache dir is recoverable; refusal to start
        // is not.
        if let Err(e) = crate::pubchem::clear_stale_locks() {
            tracing::warn!(error = %e, "failed to clear stale PubChem cache lock at startup");
        }
        if let Err(e) = crate::kegg::clear_stale_locks() {
            tracing::warn!(error = %e, "failed to clear stale KEGG cid_to_cpd cache lock at startup");
        }

        let kegg = KeggClient::new();

        // Eagerly load the organism list before the user interacts with
        // any UI (Track C / kegg-fetching spec). Initial state is
        // `Initializing` showing a splash; a cache hit is sub-frame,
        // a cold network fetch shows the spinner during `/list/organism`.
        // Read any existing cache (even stale) to provide a "Use cached"
        // fallback on failure.
        let fallback_cache = match crate::kegg::cache::read_organisms() {
            Ok(c) => c,
            Err(e) => {
                tracing::warn!(error = %e, "failed to read organisms cache for fallback");
                None
            }
        };

        let load_rx = spawn_eager_organism_load(&kegg, &rt);

        let kegg_cache_dir_env = std::env::var("KEGG_CACHE_DIR")
            .ok()
            .filter(|s| !s.is_empty());

        Self {
            state: AppState::Initializing {
                load_rx,
                fallback_cache,
                last_error: None,
            },
            // Three sibling structs introduced by `refactor-session-settings`.
            // Populated alongside `state` from this Phase 2 commit; Phase 3
            // will migrate every read/write of user-tunable parameters to
            // these and slim the `AppState` variants accordingly.
            settings: SessionSettings::default(),
            inputs: SessionInputs::default(),
            cache: SessionCache::default(),
            log,
            log_ui: log_pane::LogPaneState {
                active_tab: log_pane::BottomTab::default(),
                filter_directive: filter_directive.clone(),
                bundle_export_requested: false,
                settings_save_requested: false,
                settings_load_requested: false,
                refresh_pubchem_requested: false,
                refresh_kegg_conv_requested: false,
                rerun_enrichment_requested: false,
                refresh_catalogue_requested: false,
            },
            rt,
            kegg,
            organisms: OrganismsLoad::default(),
            species_selector: SpeciesSelectorState::default(),
            organism_group_selector: OrganismGroupSelectorState::default(),
            session_log_path,
            rust_log_directive: filter_directive,
            kegg_cache_dir_env,
            bundle_modal: BundleModalState::Closed,
            settings_save_modal: SettingsSaveModalState::Closed,
            settings_load_modal: SettingsLoadModalState::Closed,
            error_modal: ErrorModalState::Closed,
            stepper_icons: None,
            rat_face_tex: None,
            show_rat_easter_egg: false,
        }
    }

    /// Discard the current analysis and start a fresh session at Stage 1.
    /// Invoked by the `Start a new analysis` confirmation on the Stage 3
    /// result screen (see `stage3-ui` / `app-shell` specs).
    ///
    /// This is the only explicit, user-initiated FULL reset — distinct from
    /// navigation transitions, which preserve every setting/input/cache for
    /// the session lifetime via the six no-op reset APIs. It deliberately
    /// does NOT call `reset_for_back_to_stage1` (that is the preserve path);
    /// it performs a wholesale `= ::default()` on the three siblings and
    /// resets the selector UI control state too. `organisms` / `rt` / `kegg`
    /// are session-immutable infrastructure and are left untouched, so the
    /// `Initializing` organism-load splash does NOT re-run.
    pub fn start_new_round(&mut self) {
        self.settings = SessionSettings::default();
        self.inputs = SessionInputs::default();
        self.cache = SessionCache::default();
        self.species_selector = SpeciesSelectorState::default();
        self.organism_group_selector = OrganismGroupSelectorState::default();
        self.state = AppState::Stage1Input {
            slot1_mode: None,
            slot2_revealed: false,
            slot2_mode: None,
            error: None,
        };
        // Count-only log line (no input/sample names) so bug-report bundles
        // record the reset while staying privacy-safe.
        tracing::info!("new analysis round started — session reset");
    }

    /// Renders the bug-report confirm modal and the "Writing bundle…"
    /// progress modal. Polls the result channel each frame and
    /// transitions back to `Closed` on completion.
    pub fn show_bundle_modal(&mut self, ctx: &Context) {
        match &self.bundle_modal {
            BundleModalState::Closed => {}
            BundleModalState::Confirming => self.show_bundle_confirm_modal(ctx),
            BundleModalState::Saving { .. } => self.show_bundle_saving_modal(ctx),
        }
    }

    fn show_bundle_confirm_modal(&mut self, ctx: &Context) {
        use egui::{Align2, Window};

        let mut want_save = false;
        let mut want_cancel = false;
        let mut open_easter_egg = false;

        // Lazily upload rat_face texture; extract Copy IDs before the closure
        // to avoid holding a borrow on `self` across the `show` closure.
        let rat = {
            let tex = self.rat_face_tex.get_or_insert_with(|| {
                let bytes = include_bytes!("../assets/rat_face.png");
                let img = image::load_from_memory(bytes)
                    .expect("rat_face.png should decode")
                    .to_rgba8();
                let (w, h) = (img.width() as usize, img.height() as usize);
                let color = egui::ColorImage::from_rgba_unmultiplied([w, h], img.as_raw());
                ctx.load_texture("rat_face", color, egui::TextureOptions::LINEAR)
            });
            (tex.id(), tex.size_vec2())
        };

        Window::new("Download bug report")
            .title_bar(false)
            .collapsible(false)
            .resizable(false)
            .anchor(Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                // Custom title row: clickable rat icon + heading text.
                ui.add_space(4.0);
                ui.horizontal(|ui| {
                    ui.add_space(4.0);
                    let resp = ui
                        .add(
                            egui::Image::from_texture(egui::load::SizedTexture::new(rat.0, rat.1))
                                .fit_to_exact_size(egui::vec2(22.0, 22.0))
                                .sense(egui::Sense::click()),
                        )
                        .on_hover_cursor(egui::CursorIcon::PointingHand);
                    if resp.clicked() {
                        open_easter_egg = true;
                    }
                    ui.heading("Download bug report");
                });
                ui.add_space(4.0);
                ui.separator();
                ui.add_space(4.0);

                ui.label("This zip will contain the following files:");
                ui.add_space(4.0);
                for name in &[
                    "README.txt",
                    "version.txt",
                    "RUST_LOG.txt",
                    "KEGG_CACHE_DIR.txt",
                    "logs.txt",
                    "app_state.txt",
                    "input_summary.txt",
                    "cache_summary.txt",
                ] {
                    ui.label(format!("  • {name}"));
                }
                ui.add_space(8.0);
                ui.label(crate::diagnostics::BUNDLE_PRIVACY_LINE);
                ui.add_space(8.0);
                ui.horizontal(|ui| {
                    if ui.button("Cancel").clicked() {
                        want_cancel = true;
                    }
                    if ui.button("Save…").clicked() {
                        want_save = true;
                    }
                });
                ui.add_space(4.0);
            });

        if open_easter_egg {
            self.show_rat_easter_egg = true;
        }
        if want_cancel {
            self.bundle_modal = BundleModalState::Closed;
            return;
        }
        if want_save {
            let Some(path) = rfd::FileDialog::new()
                .add_filter("zip", &["zip"])
                .set_file_name(default_bundle_filename())
                .save_file()
            else {
                // User cancelled the OS save dialog — close without writing.
                self.bundle_modal = BundleModalState::Closed;
                return;
            };
            let result_rx = self.dispatch_bundle_build(path);
            self.bundle_modal = BundleModalState::Saving { result_rx };
        }
    }

    fn show_rat_easter_egg_window(&mut self, ctx: &Context) {
        if !self.show_rat_easter_egg {
            return;
        }
        let mut open = true;
        egui::Window::new("Rat Gallery")
            .open(&mut open)
            .default_width(1280.0)
            .resizable(true)
            .show(ctx, |ui| {
                ui.add(
                    egui::Image::new(
                        "https://drive.google.com/thumbnail?id=1EuOcfu4gizqi6DORjKiM5a9MYmBpfbmb&sz=s1280",
                    )
                    .shrink_to_fit(),
                );
            });
        if !open {
            self.show_rat_easter_egg = false;
        }
    }

    fn show_bundle_saving_modal(&mut self, ctx: &Context) {
        use egui::{Align2, Window};

        // Poll the result channel; if a result has arrived, log + close.
        let mut finished = None;
        if let BundleModalState::Saving { result_rx } = &self.bundle_modal
            && let Ok(msg) = result_rx.try_recv()
        {
            finished = Some(msg);
        }
        match finished {
            Some(Ok((path, size))) => {
                info!(
                    path = %path.display(),
                    size = size,
                    "bug report bundle written"
                );
                self.bundle_modal = BundleModalState::Closed;
                return;
            }
            Some(Err(e)) => {
                error!(error = %e, "bug report bundle write failed");
                self.bundle_modal = BundleModalState::Closed;
                return;
            }
            None => {}
        }

        Window::new("Writing bundle…")
            .collapsible(false)
            .resizable(false)
            .anchor(Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                ui.label("Assembling bug report bundle, please wait…");
            });
    }

    fn dispatch_bundle_build(
        &self,
        path: PathBuf,
    ) -> mpsc::Receiver<Result<(PathBuf, u64), String>> {
        use crate::diagnostics::{
            BundleArgs, build_bundle, render_app_state, render_cache_summary, render_input_summary,
        };

        let (tx, rx) = mpsc::channel();
        let app_state_text =
            render_app_state(&self.state, &self.settings, &self.inputs, &self.cache);
        let input_summary_text = render_input_summary(&self.state, &self.inputs);
        let kegg_cache_root = crate::kegg::cache::cache_dir();
        let pubchem_cache_root = crate::pubchem::cache::cache_dir();
        let cache_summary_text =
            render_cache_summary(&[kegg_cache_root.as_path(), pubchem_cache_root.as_path()]);
        let session_log_path = self.session_log_path.clone();
        let directive = self.rust_log_directive.clone();
        let kegg_env = self.kegg_cache_dir_env.clone();

        std::thread::spawn(move || {
            let result: Result<(PathBuf, u64), String> = (|| {
                let args = BundleArgs {
                    session_log_path: session_log_path.as_deref(),
                    rust_log_directive: &directive,
                    kegg_cache_dir: kegg_env.as_deref(),
                    app_state_text,
                    input_summary_text,
                    cache_summary_text,
                };
                let bytes = build_bundle(args).map_err(|e| e.to_string())?;
                let size = bytes.len() as u64;
                std::fs::write(&path, bytes).map_err(|e| e.to_string())?;
                Ok((path, size))
            })();
            let _ = tx.send(result);
        });

        rx
    }

    /// Render the shared error-toast modal (if Open).
    fn show_error_modal(&mut self, ctx: &Context) {
        use egui::{Align2, RichText, Window};

        let ErrorModalState::Open { title, message } = self.error_modal.clone() else {
            return;
        };

        let mut want_close = false;
        Window::new(RichText::new(&title).color(theme::ERROR))
            .collapsible(false)
            .resizable(false)
            .anchor(Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ctx, |ui| {
                ui.label(&message);
                ui.add_space(8.0);
                if ui.button("OK").clicked() {
                    want_close = true;
                }
            });
        if want_close {
            self.error_modal = ErrorModalState::Closed;
        }
    }

    pub(crate) fn open_error_toast_for_save(&mut self, e: &crate::session_io::SnapshotError) {
        self.error_modal = ErrorModalState::Open {
            title: "Cannot save settings".to_string(),
            message: e.to_string(),
        };
    }

    pub(crate) fn open_error_toast_for_load(&mut self, e: &crate::session_io::SnapshotError) {
        self.error_modal = ErrorModalState::Open {
            title: "Cannot load settings".to_string(),
            message: e.to_string(),
        };
    }

    /// Re-issue the eager organism load (user clicked Retry on the
    /// Initializing splash).
    pub fn retry_organism_load(&mut self) {
        if !matches!(self.state, AppState::Initializing { .. }) {
            tracing::warn!("retry_organism_load called outside Initializing; ignoring");
            return;
        }
        let load_rx = spawn_eager_organism_load(&self.kegg, &self.rt);
        if let AppState::Initializing {
            load_rx: rx,
            last_error,
            ..
        } = &mut self.state
        {
            *rx = load_rx;
            *last_error = None;
        }
    }

    /// Accept a stale cache (user clicked "Use cached organisms") and
    /// transition into Stage1Input.
    pub fn accept_fallback_cache(&mut self) {
        let prev = std::mem::take(&mut self.state);
        if let AppState::Initializing {
            fallback_cache: Some(cache),
            ..
        } = prev
        {
            info!(
                count = cache.organisms.len(),
                fetched_at = %cache.fetched_at,
                "accepted stale organisms cache as fallback"
            );
            self.organisms.state = OrganismsLoadState::Loaded(cache.organisms);
            self.state = AppState::default();
        } else {
            // Restore prior state if there was no fallback to accept.
            self.state = prev;
        }
    }

    /// Trigger an async load of the organism list. Idempotent (no-op while
    /// already loading or already loaded successfully).
    pub fn ensure_organisms_loading(&mut self) {
        match self.organisms.state {
            OrganismsLoadState::Loaded(_) | OrganismsLoadState::Loading { .. } => return,
            _ => {}
        }
        let (tx, rx) = mpsc::channel::<OrganismsLoadResult>();
        let client = self.kegg.clone();
        self.rt.spawn(async move {
            let result = match crate::kegg::list_organisms(&client).await {
                Ok(v) => OrganismsLoadResult::Ok(v),
                Err(e) => {
                    error!(error = %e, "list_organisms failed");
                    OrganismsLoadResult::Err(e.to_string())
                }
            };
            // Best-effort send; UI may have closed the channel if the user moved on.
            let _ = tx.send(result);
        });
        self.organisms.state = OrganismsLoadState::Loading { rx };
        info!("organisms load started");
    }

    /// Drain the eager-startup organism-load channel while in
    /// `AppState::Initializing`. On Ok: store the list on `app.organisms`
    /// and transition to `Stage1Input`. On Err: set `last_error` and stay
    /// in `Initializing` for the user to Retry.
    fn drain_initializing(&mut self) {
        let received = if let AppState::Initializing { load_rx, .. } = &self.state {
            load_rx.try_recv().ok()
        } else {
            None
        };
        let Some(result) = received else {
            return;
        };
        match result {
            Ok(organisms) => {
                info!(count = organisms.len(), "eager organisms load completed");
                self.organisms.state = OrganismsLoadState::Loaded(organisms);
                self.state = AppState::default();
            }
            Err(msg) => {
                error!(error = %msg, "eager organisms load failed");
                if let AppState::Initializing { last_error, .. } = &mut self.state {
                    *last_error = Some(msg);
                }
            }
        }
    }

    /// Drain the organism load channel if present.
    fn drain_organisms_load(&mut self) {
        let received = if let OrganismsLoadState::Loading { rx } = &self.organisms.state {
            rx.try_recv().ok()
        } else {
            None
        };
        if let Some(result) = received {
            match result {
                OrganismsLoadResult::Ok(v) => {
                    info!(count = v.len(), "organisms load completed");
                    self.organisms.state = OrganismsLoadState::Loaded(v);
                }
                OrganismsLoadResult::Err(msg) => {
                    self.organisms.state = OrganismsLoadState::Failed(msg);
                }
            }
        }
    }
}

impl EframeApp for App {
    fn update(&mut self, ctx: &Context, _frame: &mut Frame) {
        self.drain_organisms_load();
        self.drain_initializing();

        self.drain_stage3_setup_fetch();
        self.drain_stage2_running();
        self.drain_stage3_running();

        // Bottom panel — a two-tab container (Data | Log) per
        // `add-bottom-panel-data-tab`. The tab strip + per-tab visibility rules
        // live in `bottom_panel::show`; the Log tab preserves the prior pane
        // behaviour verbatim. Panel id kept as "log_pane" so the resizable
        // height the user dragged persists across the rename.
        TopBottomPanel::bottom("log_pane")
            .resizable(true)
            // ~15 em at egui's 14 px body font — a roomier default so the Data
            // tab's stage summary is visible without dragging.
            .default_height(210.0)
            .min_height(60.0)
            // Page-chrome secondary background (`#CED8D9`) so the panel reads as
            // a distinct surface over the central `BACKGROUND`. Override fill
            // only; keep the default panel insets/stroke.
            .frame(egui::Frame::side_top_panel(&ctx.style()).fill(theme::BACKGROUND_SECONDARY))
            .show(ctx, |ui| {
                bottom_panel::show(ui, self);
            });

        self.drain_modal_requests();
        self.render_modals(ctx);

        // Global stage stepper / breadcrumb — clickable back-navigation, above
        // the screen body. Rendered after the bottom panel so it claims the top
        // strip; skips itself during `Initializing` (no navigation possible).
        TopBottomPanel::top("stepper")
            // Same page-chrome secondary background as the bottom panel.
            .frame(egui::Frame::side_top_panel(&ctx.style()).fill(theme::BACKGROUND_SECONDARY))
            .show(ctx, |ui| {
                stepper::show(ui, self);
            });

        CentralPanel::default().show(ctx, |ui| {
            // Dispatch by current variant. Stage1Input and Stage3EnrichSetup
            // need app-level state (organism list, runtime) so we route
            // through `App` methods that own those references.
            match &self.state {
                AppState::Initializing { .. } => {
                    initializing::show(ui, self);
                }
                AppState::Stage1Input { .. } => {
                    stage1_input::show(ui, self);
                }
                AppState::Stage2DamSetup { .. } => {
                    stage2_setup::show(ui, self);
                }
                AppState::Stage2DamRunning { .. } => {
                    stage2_running::show(ui, self);
                }
                AppState::Stage2DamThreshold { .. } => {
                    stage2_threshold::show(ui, self);
                }
                AppState::Stage3EnrichSetup { .. } => {
                    stage3_setup::show(ui, self);
                }
                AppState::Stage3EnrichRunning { .. } => {
                    stage3_running::show(ui, self);
                }
                AppState::Stage3EnrichResult { .. } => {
                    stage3_result::show(ui, self);
                }
            }
        });

        // Ensure new log lines + KEGG progress surface promptly even while the
        // user is idle.
        ctx.request_repaint_after(Duration::from_millis(250));
    }
}

impl App {
    /// Drain any in-flight KEGG / module fetch progress while on Stage 3
    /// setup. Both `kegg_fetch` and `modules_fetch` live as optional
    /// fields on `Stage3EnrichSetup`; each is drained independently so
    /// that toggling mode mid-fetch does not interrupt either stream
    /// (per `reorder-gui-and-move-mode-to-stage3` D6).
    fn drain_stage3_setup_fetch(&mut self) {
        if let AppState::Stage3EnrichSetup {
            kegg_fetch,
            modules_fetch,
            ..
        } = &mut self.state
        {
            let kegg_terminal = kegg_fetch.as_mut().and_then(drain_kegg_progress);
            let modules_terminal = modules_fetch.as_mut().and_then(drain_modules_progress);
            if let Some(event) = kegg_terminal {
                self.handle_kegg_terminal_event(event);
            }
            if let Some(event) = modules_terminal {
                self.handle_modules_fetch_terminal_event(event);
            }
        }
    }

    /// Drain DAM progress + check for per-mode terminal results. We only
    /// call `handle_dam_terminal` once every spawned worker has reported
    /// (i.e. every slot in `dam_results_accum` is `Some`).
    fn drain_stage2_running(&mut self) {
        let mut all_done = false;
        if let AppState::Stage2DamRunning {
            progress_rxs,
            mode_completed,
            mode_total,
            result_rx,
            dam_results_accum,
            ..
        } = &mut self.state
        {
            for (idx, rx) in progress_rxs.iter().enumerate() {
                while let Ok(p) = rx.try_recv() {
                    if let Some(c) = mode_completed.get_mut(idx) {
                        *c = p.completed;
                    }
                    if let Some(t) = mode_total.get_mut(idx) {
                        *t = p.total;
                    }
                }
            }
            while let Ok((idx, res)) = result_rx.try_recv() {
                if let Some(slot) = dam_results_accum.get_mut(idx) {
                    *slot = Some(res);
                }
            }
            all_done =
                !dam_results_accum.is_empty() && dam_results_accum.iter().all(|x| x.is_some());
        }
        if all_done {
            self.handle_dam_terminal();
        }
    }

    /// Drain Stage 3 progress + check for terminal result.
    fn drain_stage3_running(&mut self) {
        if let AppState::Stage3EnrichRunning {
            phase,
            pubchem_progress_rx,
            kegg_conv_progress_rx,
            pubchem_completed,
            pubchem_total,
            kegg_conv_completed,
            kegg_conv_total,
            result_rx,
            ..
        } = &mut self.state
        {
            while let Ok(p) = pubchem_progress_rx.try_recv() {
                *pubchem_completed = p.from_cache + p.fetched;
                *pubchem_total = p.total_inputs.max(*pubchem_total);
                // Advance phase from PubChem to KeggConv when the first
                // PubChem batch completes — heuristic, but the running
                // screen's phase label is purely informational.
                if matches!(phase, Stage3Phase::PubChem)
                    && *pubchem_completed >= *pubchem_total
                    && *pubchem_total > 0
                {
                    *phase = Stage3Phase::KeggConv;
                }
            }
            while let Ok(p) = kegg_conv_progress_rx.try_recv() {
                *kegg_conv_completed = p.from_cache + p.fetched;
                *kegg_conv_total = p.total_inputs.max(*kegg_conv_total);
                if matches!(phase, Stage3Phase::KeggConv | Stage3Phase::PubChem) {
                    *phase = Stage3Phase::KeggConv;
                }
            }
            let terminal = result_rx.try_recv().ok();
            if let Some(msg) = terminal {
                self.handle_stage3_terminal(msg);
            }
        }
    }

    /// Drain the log-pane modal-request flags. Each flag is reset every
    /// frame regardless of whether the modal opens — see the
    /// mutual-exclusion invariant in the app-shell spec. If any modal is
    /// already non-Closed, drop the new request with a `warn!` line (the
    /// user dismisses the current modal first).
    fn drain_modal_requests(&mut self) {
        let any_modal_open = !matches!(self.bundle_modal, BundleModalState::Closed)
            || !matches!(self.settings_save_modal, SettingsSaveModalState::Closed)
            || !matches!(self.settings_load_modal, SettingsLoadModalState::Closed)
            || !matches!(self.error_modal, ErrorModalState::Closed);

        if self.log_ui.bundle_export_requested {
            self.log_ui.bundle_export_requested = false;
            if !any_modal_open {
                self.bundle_modal = BundleModalState::Confirming;
            } else {
                tracing::warn!("bundle export request dropped: another modal is open");
            }
        }

        if self.log_ui.settings_save_requested {
            self.log_ui.settings_save_requested = false;
            if !any_modal_open {
                self.settings_save_modal = SettingsSaveModalState::Confirming;
            } else {
                tracing::warn!("settings save request dropped: another modal is open");
            }
        }

        if self.log_ui.settings_load_requested {
            self.log_ui.settings_load_requested = false;
            if !any_modal_open {
                settings_modals::open_load(self);
            } else {
                tracing::warn!("settings load request dropped: another modal is open");
            }
        }
    }

    /// Render the four App-level modals (bundle, settings save, settings
    /// load, error) plus the rat easter-egg window, in the established
    /// order. Runs after `drain_modal_requests` so a flag set this frame
    /// opens its modal this frame.
    fn render_modals(&mut self, ctx: &Context) {
        self.show_bundle_modal(ctx);
        settings_modals::show_save(self, ctx);
        settings_modals::show_load(self, ctx);
        self.show_error_modal(ctx);
        self.show_rat_easter_egg_window(ctx);
    }

    /// Handle the terminal `Done` / `Failed` event of an in-flight KEGG fetch.
    /// Transitions out of `KeggFetching` and back to `Stage1Input` with the
    /// appropriate result populated. Caller MUST have already exited any
    /// borrow on the `progress_rx` field.
    /// Handle the terminal `Ok` / `Err` from a DAM background task. Transitions out
    /// of `Stage2DamRunning` to either `Stage2DamThreshold` (success, with default
    /// thresholds) or back to `Stage2DamSetup` (failure, with the error message and
    /// the previous group / method selections preserved).
    /// Called from `update()` once every spawned DAM worker has reported a
    /// terminal `Result` into `dam_results_accum`. Folds the per-mode results
    /// into either a `Stage2DamThreshold` transition (all modes succeeded) or
    /// a `Stage2DamSetup` transition with an error message naming the failing
    /// mode(s).
    fn handle_dam_terminal(&mut self) {
        let prev = std::mem::take(&mut self.state);
        let AppState::Stage2DamRunning {
            dam_results_accum, ..
        } = prev
        else {
            return;
        };

        // Partition: all-Ok → transition forward; any-Err → return to Setup
        // with an error string naming the failing mode(s).
        let mut ok_results: Vec<DamResult> = Vec::with_capacity(dam_results_accum.len());
        let mut errors: Vec<String> = Vec::new();
        for (idx, slot) in dam_results_accum.into_iter().enumerate() {
            match slot {
                Some(Ok(r)) => ok_results.push(r),
                Some(Err(e)) => {
                    let label = self
                        .inputs
                        .ion_tables
                        .get(idx)
                        .map(|it| it.mode.to_string())
                        .unwrap_or_else(|| format!("mode #{}", idx + 1));
                    errors.push(format!("{label}: {e}"));
                }
                None => {
                    errors.push(format!("mode #{} did not report a result", idx + 1));
                }
            }
        }

        if errors.is_empty() {
            info!(
                modes = ok_results.len(),
                fdr_method = ?self.settings.dam_fdr_method,
                "DAM complete"
            );
            let active_volcano_tab = self
                .inputs
                .ion_tables
                .first()
                .map(|it| it.mode)
                .unwrap_or(IonMode::Positive);
            let volcano_textures = vec![None; self.inputs.ion_tables.len()];
            self.state = AppState::Stage2DamThreshold {
                dam_results: ok_results,
                active_volcano_tab,
                volcano_textures,
                rendering: false,
                render_rx: None,
            };
        } else {
            let msg = errors.join("; ");
            error!(error = %msg, "DAM failed");
            self.state = AppState::Stage2DamSetup { error: Some(msg) };
        }
    }

    /// Handle the terminal `Done` / `Failed` from an in-flight module
    /// fetch. The current variant MUST be `Stage3EnrichSetup` with
    /// `modules_fetch.is_some()` — on Done writes `cache.modules_pack`
    /// and clears the in-flight slot to `None` (no AppState transition);
    /// on Failed writes the error string into `Stage3EnrichSetup.error`
    /// and clears the slot.
    fn handle_modules_fetch_terminal_event(&mut self, event: ModulesFetchEvent) {
        let AppState::Stage3EnrichSetup {
            modules_fetch,
            error,
            ..
        } = &mut self.state
        else {
            error!(
                "unexpected terminal modules event outside Stage3EnrichSetup; current state ignored"
            );
            return;
        };
        match event {
            ModulesFetchEvent::Done(cache) => {
                info!(
                    cached = cache.modules.len(),
                    group = ?self.settings.organism_group,
                    "modules fetch complete"
                );
                self.cache.modules_pack = Some(cache);
                *modules_fetch = None;
                *error = None;
            }
            ModulesFetchEvent::Failed(msg) => {
                error!(error = %msg, "modules fetch failed");
                self.cache.modules_pack = None;
                *modules_fetch = None;
                *error = Some(format!("KEGG modules fetch failed: {msg}"));
            }
            ModulesFetchEvent::Progress(_) => {
                error!("unexpected Progress event in modules terminal handler");
                *modules_fetch = None;
            }
        }
    }

    /// Handle the terminal `Done` / `Failed` from an in-flight species
    /// fetch. The current variant MUST be `Stage3EnrichSetup` with
    /// `kegg_fetch.is_some()` — on Done writes `cache.species_kegg` and
    /// clears the in-flight slot to `None` (no AppState transition); on
    /// Failed writes the error string into `Stage3EnrichSetup.error` and
    /// clears the slot.
    fn handle_kegg_terminal_event(&mut self, event: KeggEvent) {
        let species = self.settings.kegg_species.clone().unwrap_or_default();
        let AppState::Stage3EnrichSetup {
            kegg_fetch, error, ..
        } = &mut self.state
        else {
            error!(
                ?event,
                "unexpected terminal KEGG event outside Stage3EnrichSetup; current state ignored"
            );
            return;
        };
        match event {
            KeggEvent::Done(species_kegg) => {
                info!(code = %species, pathways = species_kegg.pathways.len(), "KEGG fetch complete");
                self.cache.species_kegg = Some(species_kegg);
                *kegg_fetch = None;
                *error = None;
            }
            KeggEvent::Failed(msg) => {
                error!(code = %species, error = %msg, "KEGG fetch failed");
                self.cache.species_kegg = None;
                *kegg_fetch = None;
                *error = Some(format!("KEGG fetch failed: {msg}"));
            }
            KeggEvent::Progress(_) => {
                error!("unexpected Progress event in KEGG terminal handler");
                *kegg_fetch = None;
            }
        }
    }

    /// Spawn the Stage 3 orchestrator and transition into
    /// `AppState::Stage3EnrichRunning`. Shared by `start_run` and `rerun` (both
    /// perform the identical channel-trio wiring + `run_stage3` spawn +
    /// `Stage3EnrichRunning` construction). `start_refresh` does NOT use this —
    /// it keeps its own in-place `RefreshState` transition + progress bridges.
    ///
    /// Log-silent: the only diagnostic that differed across callers was
    /// `start_run`'s `n_modes` `info!`, which stays at the `start_run` call site
    /// (emitted before this call). The orchestrator's terminal error is logged
    /// once by `handle_stage3_terminal`.
    pub(crate) fn spawn_stage3_run(
        &mut self,
        dam_results: Vec<DamResult>,
        target: AnalysisPayload,
        params: crate::stage3::Stage3Params,
        pubchem_total: usize,
    ) {
        let (pub_tx, pub_rx) = mpsc::channel();
        let (kegg_tx, kegg_rx) = mpsc::channel();
        let (result_tx, result_rx) = mpsc::channel::<Result<Stage3RunOutput, String>>();
        let kegg_client = self.kegg.clone();
        let dam_results_clone = dam_results.clone();
        self.rt.spawn(async move {
            let pubchem = crate::pubchem::PubchemClient::new();
            let r = crate::stage3::run_stage3(
                &pubchem,
                &kegg_client,
                &dam_results_clone,
                &target,
                params,
                pub_tx,
                kegg_tx,
            )
            .await;
            let _ = result_tx.send(r.map_err(|e| e.to_string()));
        });
        self.state = AppState::Stage3EnrichRunning {
            dam_results,
            phase: Stage3Phase::PubChem,
            pubchem_progress_rx: pub_rx,
            kegg_conv_progress_rx: kegg_rx,
            result_rx,
            pubchem_completed: 0,
            pubchem_total,
            kegg_conv_completed: 0,
            kegg_conv_total: 0,
        };
    }

    /// Handle the orchestrator's terminal result for Stage 3. On Ok,
    /// move into `Stage3EnrichResult` with the payload populated; on
    /// Err, return to `Stage3EnrichSetup` with the error preserved.
    fn handle_stage3_terminal(&mut self, msg: Result<Stage3RunOutput, String>) {
        let prev = std::mem::take(&mut self.state);
        let AppState::Stage3EnrichRunning { dam_results, .. } = prev else {
            return;
        };
        match msg {
            Ok(out) => {
                // Rows the dot plot will actually show — reused for the
                // run-complete log and the initial export-height auto-size.
                let displayed = out
                    .enrichment_result
                    .rows
                    .iter()
                    .filter(|r| r.fdr < self.settings.enrichment_fdr_threshold && r.displayed)
                    .count();
                info!(
                    universe = out.mapped_universe.len(),
                    mapped_features = out.feature_to_cpds.len(),
                    fdr_method = ?self.settings.enrichment_fdr_method,
                    sig_pathways = displayed,
                    "Stage 3 Run complete"
                );
                // Auto-size the initial export height to min(top_n, displayed)
                // so a sparse result doesn't render in a tall band of
                // whitespace. Result-entry only; the result-screen Height
                // field remains a user override.
                self.settings.stage3_export_height_in =
                    stage3_autosize_height_in(self.settings.top_n, displayed);
                self.state = out.into_result_state(dam_results);
            }
            Err(e) => {
                error!(error = %e, "Stage 3 Run failed");
                self.state = AppState::Stage3EnrichSetup {
                    dam_results,
                    error: Some(e),
                    kegg_fetch: None,
                    modules_fetch: None,
                };
            }
        }
    }
}

/// Spawn the eager organism-list load on the tokio runtime and return a
/// receiver the UI thread will drain each frame while in
/// `AppState::Initializing`. Reads from cache when present (sub-frame);
/// otherwise issues `GET /list/organism`. Either path also rewrites the
/// derived `organism_groups.json` precompute (see `kegg::list_organisms`).
fn spawn_eager_organism_load(
    client: &KeggClient,
    rt: &Runtime,
) -> mpsc::Receiver<Result<Vec<KeggOrganism>, String>> {
    let (tx, rx) = mpsc::channel::<Result<Vec<KeggOrganism>, String>>();
    let client = client.clone();
    rt.spawn(async move {
        let result = match crate::kegg::list_organisms(&client).await {
            Ok(v) => Ok(v),
            Err(e) => {
                error!(error = %e, "eager list_organisms failed");
                Err(e.to_string())
            }
        };
        let _ = tx.send(result);
    });
    info!("eager organisms load started");
    rx
}

#[cfg(test)]
mod tests {
    //! Unit tests for `SessionSettings` / `SessionCache` defaults, named
    //! reset APIs, and cache clearing behavior. The reset surfaces here
    //! lock in today's pre-refactor behavior bit-equally — see the
    //! `refactor-session-settings` change's design notes (D11) for
    //! the full inventory.
    use super::*;

    #[test]
    fn stage3_autosize_height_sparse_result_is_compact() {
        // top_n=20 but only 8 rows displayed → sized to 8, not 20.
        assert!((stage3_autosize_height_in(20, 8) - 3.4).abs() < 1e-9);
    }

    #[test]
    fn stage3_autosize_height_dense_result_capped_by_top_n() {
        // 50 pass but top_n=20 caps the plot → identical to the old formula.
        assert!((stage3_autosize_height_in(20, 50) - 7.0).abs() < 1e-9);
    }

    #[test]
    fn stage3_autosize_height_empty_result_clamps_to_minimum() {
        // 0 displayed → max(1) → 1.3 → clamped up to the 2.0 floor (raised from
        // 1.5 so the now-full-size legend fits a sparse-result canvas).
        assert!((stage3_autosize_height_in(20, 0) - 2.0).abs() < 1e-9);
        assert!((stage3_autosize_height_in(1, 1) - 2.0).abs() < 1e-9);
    }

    #[test]
    fn stage3_autosize_height_clamps_to_upper_bound() {
        // min(200, 200) = 200 → 61.0 → clamped to the 40.0 ceiling.
        assert!((stage3_autosize_height_in(200, 200) - 40.0).abs() < 1e-9);
    }

    /// Construct a `SessionSettings` with every field deliberately set
    /// to a non-default value. Used by the `reset_*` tests to verify
    /// that a method clears only its documented field surface and
    /// leaves the rest verbatim. Field values are chosen to be valid
    /// (e.g. `min_group_overlap = 5`, `top_n = 50`, FDR method = BY)
    /// so the test compiles even if a future change tightens
    /// validation.
    fn non_default_settings() -> SessionSettings {
        SessionSettings {
            // Stage 1 / mode
            analysis_mode: AnalysisMode::Module,
            kegg_species: Some("hsa".to_string()),
            organism_group_level: Some(2),
            organism_group: Some("Mammals".to_string()),
            min_group_overlap: 5,

            // Stage 2 setup
            numerator: Some("treatment".to_string()),
            denominator: Some("control".to_string()),
            dam_method: DamMethod::Welch,
            drop_unknown: false,
            dedup_enabled: false,
            normalization: NormalizationMethod::Sum,
            metadata_column: Some("dry_weight".to_string()),
            pqn_reference: PqnReference::Group("control".to_string()),
            pqn_reference_group: Some("control".to_string()),
            log_transform: false,
            dam_fdr_method: FdrMethod::BenjaminiYekutieli,

            // Stage 2 result
            fc_threshold: 4.0,
            fdr_threshold: 0.01,
            delta_threshold: 0.5,
            stage2_export_width_in: 6.0,
            stage2_export_height_in: 4.0,
            stage2_export_dpi: 600,

            // Stage 3 setup
            direction: EnrichmentDirection::Up,
            top_n: 50,
            enrichment_fdr_threshold: 0.1,
            min_hit_count: 3,
            min_entry_size: 5,
            enrichment_fdr_method: FdrMethod::BenjaminiYekutieli,

            // Stage 3 result
            stage3_export_width_in: 5.0,
            stage3_export_height_in: 10.0,
            stage3_export_dpi: 600,
        }
    }

    // ── Task 1.11(a): Default values match documented defaults ──

    #[test]
    fn default_settings_match_documented_values() {
        let d = SessionSettings::default();
        // Stage 1 / mode
        assert_eq!(d.analysis_mode, AnalysisMode::Pathway);
        assert_eq!(d.kegg_species, None);
        assert_eq!(d.organism_group_level, None);
        assert_eq!(d.organism_group, None);
        assert_eq!(d.min_group_overlap, 1);
        // Stage 2 setup
        assert_eq!(d.numerator, None);
        assert_eq!(d.denominator, None);
        assert_eq!(d.dam_method, DamMethod::Student);
        assert!(d.drop_unknown);
        assert!(d.dedup_enabled);
        assert_eq!(d.normalization, NormalizationMethod::None);
        assert_eq!(d.metadata_column, None);
        assert_eq!(d.pqn_reference, PqnReference::AllSamples);
        assert_eq!(d.pqn_reference_group, None);
        assert_eq!(d.dam_fdr_method, FdrMethod::BenjaminiHochberg);
        // Stage 2 result
        assert_eq!(d.fc_threshold, 2.0);
        assert_eq!(d.fdr_threshold, 0.05);
        assert_eq!(d.delta_threshold, 0.33);
        assert_eq!(d.stage2_export_width_in, 3.5);
        assert_eq!(d.stage2_export_height_in, 2.2);
        assert_eq!(d.stage2_export_dpi, 300);
        // Stage 3 setup
        assert_eq!(d.direction, EnrichmentDirection::Both);
        assert_eq!(d.top_n, 20);
        assert_eq!(d.enrichment_fdr_threshold, 0.05);
        assert_eq!(d.min_hit_count, 1);
        assert_eq!(d.enrichment_fdr_method, FdrMethod::BenjaminiHochberg);
        // Stage 3 result (default height = top_n * 0.3 + 1.0 = 7.0 for top_n = 20)
        assert_eq!(d.stage3_export_width_in, 3.5);
        assert_eq!(d.stage3_export_height_in, 7.0);
        assert_eq!(d.stage3_export_dpi, 300);
    }

    // ── Task 1.11(b): Each reset_* / clear_* touches only its surface ──

    #[test]
    fn reset_stage2_choices_on_change_comparison_preserves_all_fields() {
        // After `reorder-gui-and-move-mode-to-stage3` Phase 2, this method
        // is a no-op: pressing "Back to DAM Setup" preserves every Stage 2
        // settings field so the user can adjust one choice and re-run
        // without losing the others.
        let mut s = non_default_settings();
        let baseline = s.clone();
        s.reset_stage2_choices_on_change_comparison();
        assert_eq!(
            s, baseline,
            "no settings field MUST change on Back to DAM Setup"
        );
    }

    #[test]
    fn reset_for_back_to_stage2_threshold_preserves_all_fields() {
        // Post-smoke-test feedback: settings persist across the Stage 3
        // setup → Stage 2 result Back transition. The method body is a
        // no-op; this test locks that in.
        let mut s = non_default_settings();
        let baseline = s.clone();
        s.reset_for_back_to_stage2_threshold();
        assert_eq!(
            s, baseline,
            "Back to DAM Result MUST preserve every settings field"
        );
    }

    #[test]
    fn reset_stage3_on_continue_to_enrichment_preserves_all_fields() {
        // Post-smoke-test feedback: Stage 3 settings persist across
        // Continue-to-Enrichment transitions so the user's prior
        // direction / top_n / FDR / min hit count carry forward.
        let mut s = non_default_settings();
        let baseline = s.clone();
        s.reset_stage3_on_continue_to_enrichment();
        assert_eq!(
            s, baseline,
            "Continue to Enrichment MUST preserve every settings field"
        );
    }

    #[test]
    fn reset_for_back_to_stage1_preserves_all_fields() {
        // Post-smoke-test feedback: settings persist across the Stage 2
        // setup → Stage 1 Back transition. If the user re-picks files
        // such that the preserved numerator / denominator groups no
        // longer exist, the Stage 2 setup gate refuses to start DAM
        // until the user re-selects valid groups; this is enforced at
        // the UI gate level, not by clearing settings here.
        let mut s = non_default_settings();
        let baseline = s.clone();
        s.reset_for_back_to_stage1();
        assert_eq!(
            s, baseline,
            "Back to Input MUST preserve every settings field"
        );
    }

    // ── Task 1.11(c): Mode-switch cleanup in both directions ──

    #[test]
    fn reset_kegg_selection_pathway_to_module_only_changes_analysis_mode() {
        // After Phase 2: only `analysis_mode` flips; both modes' selection
        // fields are preserved so the user can toggle freely.
        let mut s = non_default_settings();
        s.analysis_mode = AnalysisMode::Pathway;
        s.kegg_species = Some("hsa".to_string());
        s.organism_group_level = Some(2);
        s.organism_group = Some("Mammals".to_string());
        s.min_group_overlap = 5;

        s.reset_kegg_selection_for_mode_switch(AnalysisMode::Module);

        assert_eq!(s.analysis_mode, AnalysisMode::Module);
        assert_eq!(s.kegg_species, Some("hsa".to_string())); // preserved
        assert_eq!(s.organism_group_level, Some(2));
        assert_eq!(s.organism_group, Some("Mammals".to_string()));
        assert_eq!(s.min_group_overlap, 5);
    }

    #[test]
    fn reset_kegg_selection_module_to_pathway_only_changes_analysis_mode() {
        let mut s = non_default_settings();
        s.analysis_mode = AnalysisMode::Module;
        s.kegg_species = Some("hsa".to_string());
        s.organism_group_level = Some(2);
        s.organism_group = Some("Mammals".to_string());
        s.min_group_overlap = 5;

        s.reset_kegg_selection_for_mode_switch(AnalysisMode::Pathway);

        assert_eq!(s.analysis_mode, AnalysisMode::Pathway);
        // Both selections preserved across the toggle.
        assert_eq!(s.kegg_species, Some("hsa".to_string()));
        assert_eq!(s.organism_group_level, Some(2));
        assert_eq!(s.organism_group, Some("Mammals".to_string()));
        assert_eq!(s.min_group_overlap, 5);
    }

    #[test]
    fn session_cache_clear_for_mode_switch_is_a_no_op() {
        // After Phase 2: Pathway and Module caches coexist for the
        // lifetime of the session — `clear_for_mode_switch` no longer
        // touches any field.
        let mut cache = SessionCache {
            species_kegg: None,
            modules_pack: None,
            group_org_codes: Some(std::collections::HashSet::from(["hsa".to_string()])),
        };
        let has_codes_before = cache.group_org_codes.is_some();
        cache.clear_for_mode_switch(AnalysisMode::Pathway);
        assert_eq!(cache.group_org_codes.is_some(), has_codes_before);
        cache.clear_for_mode_switch(AnalysisMode::Module);
        assert_eq!(cache.group_org_codes.is_some(), has_codes_before);
    }

    // ── add-new-analysis-reset: App::start_new_round ──────────────────

    /// Build an `App` for unit tests WITHOUT firing the eager organism
    /// load over the network. `App::new` spawns `list_organisms` onto the
    /// runtime, but a current-thread runtime never drives spawned tasks
    /// absent a `block_on`, so the task is queued and dropped untouched
    /// when the runtime is dropped at end of test — no request fires.
    fn test_app() -> App {
        let rt = tokio::runtime::Builder::new_current_thread()
            .build()
            .expect("build current-thread test runtime");
        App::new(LogStore::new(16), "info".to_string(), rt, None)
    }

    /// A fully-populated `Stage3EnrichResult` state — the variant
    /// `start_new_round` is reachable from in production.
    fn stage3_result_state() -> AppState {
        AppState::Stage3EnrichResult {
            dam_results: vec![],
            module_retention: None,
            enrichment_result: EnrichmentResult {
                universe_size: 0,
                dam_cpd_size: 0,
                direction: EnrichmentDirection::Both,
                min_hit_count: 1,
                min_entry_size: 1,
                entries_dropped_by_min_entry_size: 0,
                empty_compound_count: 0,
                rows: vec![],
                fdr_method: FdrMethod::BenjaminiHochberg,
            },
            mapped_universe: std::collections::HashSet::new(),
            feature_to_cpds: std::collections::HashMap::new(),
            pubchem_time_span: None,
            kegg_conv_time_span: None,
            dual_mode_breakdown: None,
            funnel: Stage3Funnel::default(),
            dotplot_tex: None,
            rendering: false,
            render_rx: None,
            refresh_state: RefreshState::Idle,
            confirming_new_round: false,
            height_user_overridden: false,
        }
    }

    #[test]
    fn start_new_round_resets_session_and_returns_to_stage1() {
        let mut app = test_app();
        // Dirty every resettable surface, mid-analysis.
        app.settings = non_default_settings();
        app.inputs = SessionInputs {
            ion_tables: vec![],
            mapping: None,
            csv_path: Some(std::path::PathBuf::from("/tmp/m.csv")),
        };
        app.cache = SessionCache {
            species_kegg: None,
            modules_pack: None,
            group_org_codes: Some(std::collections::HashSet::from(["hsa".to_string()])),
        };
        app.species_selector.filter = "soy".to_string();
        app.species_selector.picker_open = true;
        app.organism_group_selector.level = 1;
        app.state = stage3_result_state();

        app.start_new_round();

        // Settings fully back to documented defaults.
        assert_eq!(app.settings, SessionSettings::default());
        // Inputs cleared (no PartialEq on SessionInputs → field checks).
        assert!(app.inputs.ion_tables.is_empty());
        assert!(app.inputs.mapping.is_none());
        assert!(app.inputs.csv_path.is_none());
        // Cache cleared.
        assert!(app.cache.species_kegg.is_none());
        assert!(app.cache.modules_pack.is_none());
        assert!(app.cache.group_org_codes.is_none());
        // Selector UI control state back to fresh defaults.
        assert!(app.species_selector.filter.is_empty());
        assert!(!app.species_selector.picker_open);
        assert_eq!(app.organism_group_selector.level, 2);
        // Landed on a fresh Stage 1 input screen.
        assert!(matches!(
            app.state,
            AppState::Stage1Input {
                slot1_mode: None,
                slot2_revealed: false,
                slot2_mode: None,
                error: None,
            }
        ));
    }

    #[test]
    fn start_new_round_leaves_organism_roster_and_skips_initializing() {
        let mut app = test_app();
        app.organisms.state = OrganismsLoadState::Loaded(vec![]);
        app.state = stage3_result_state();

        app.start_new_round();

        // Session-immutable organism roster untouched (still Loaded, not
        // reset to Idle / re-fetched); and we do NOT bounce back through
        // the Initializing splash.
        assert!(matches!(app.organisms.state, OrganismsLoadState::Loaded(_)));
        assert!(!matches!(app.state, AppState::Initializing { .. }));
        assert!(matches!(app.state, AppState::Stage1Input { .. }));
    }

    #[test]
    fn confirming_new_round_defaults_false_and_toggle_is_sibling_safe() {
        // The flag defaults `false` at construction, and opening then
        // cancelling the confirm modal (flipping the flag) never touches
        // the session siblings — only Confirm (`start_new_round`) resets.
        let mut app = test_app();
        app.settings = non_default_settings();
        let settings_snapshot = app.settings.clone();
        app.state = stage3_result_state();

        // Default is false at construction.
        let AppState::Stage3EnrichResult {
            confirming_new_round,
            ..
        } = &app.state
        else {
            panic!("expected Stage3EnrichResult");
        };
        assert!(!confirming_new_round);

        // Open (button) then cancel (Cancel) — flip the flag both ways.
        if let AppState::Stage3EnrichResult {
            confirming_new_round,
            ..
        } = &mut app.state
        {
            *confirming_new_round = true;
            *confirming_new_round = false;
        }

        // The open/cancel toggle left every session sibling untouched.
        assert_eq!(app.settings, settings_snapshot);
        assert!(matches!(app.state, AppState::Stage3EnrichResult { .. }));
    }
}
