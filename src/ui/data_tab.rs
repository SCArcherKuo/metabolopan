//! The **Data** tab of the bottom panel — a read-only, stage-aware summary of
//! every per-stage data count (raw input, annotated/Unknown split, biosample
//! counts, dedup audit, pre-filter, DAM tallies, and the Stage 3 enrichment
//! provenance funnel). Defined by the `data-summary-panel` capability.
//!
//! The renderer branches on the current `AppState` discriminant and derives
//! every count on the frame from `App::{inputs, cache, settings, state}` — it
//! stores no counts of its own (the sole exception is the Stage 3 funnel
//! counts the orchestrator returns on `Stage3RunOutput`).

use std::collections::HashSet;

use chrono::{DateTime, Utc};
use egui::{RichText, ScrollArea};

use crate::app::{
    AnalysisMode, App, AppState, RefreshState, RunningPayload, SessionCache, SessionSettings,
    Stage3Funnel,
};
use crate::dam::{DamMethod, DamResult, Trend, classify_trend};
use crate::data::{GroupMapping, IonMode, IonModeTable, UNASSIGNED};
use crate::enrichment::types::EnrichmentResult;
use crate::kegg::{KeggModulesCache, SpeciesKegg};
use crate::stage3::{DualModeBreakdown, ModuleRetention};
use crate::theme;
use crate::ui::widgets::{kv_line, kv_line_colored, section_header};

/// Full ion-mode name for the slot header (`MS-DIAL data slot 1 (Positive)`).
fn mode_full(m: IonMode) -> &'static str {
    match m {
        IonMode::Positive => "Positive",
        IonMode::Negative => "Negative",
    }
}

/// Short DAM-method label for the `DAM data` block. Distinct from
/// `DamMethod::display_name()` (which returns the ASCII-hyphen `"Brunner-Munzel"`
/// with no " test" suffix and feeds the plot strips): here we want the en-dash
/// + " test" form matching the Stage 2 setup method radio.
fn method_label(m: DamMethod) -> &'static str {
    match m {
        DamMethod::Welch => "Welch's t-test",
        DamMethod::Student => "Student's t-test",
        DamMethod::BrunnerMunzel => "Brunner–Munzel test",
    }
}

/// `DAM data` block (`Comparison:` + `Statistical method:`), rendered on every
/// post-DAM stage (DAM Result + both Enrichment pages) from `dam_results[0]`.
fn render_dam_data_block(ui: &mut egui::Ui, dam_results: &[DamResult]) {
    let Some(dam) = dam_results.first() else {
        return;
    };
    section_header(ui, "DAM data");
    kv_line(
        ui,
        &format!("Comparison: {} / {}", dam.numerator, dam.denominator),
    );
    kv_line(
        ui,
        &format!("Statistical method: {}", method_label(dam.method)),
    );
    ui.add_space(6.0);
}

/// Render the Data-tab body. Takes `&mut App` because the relocated
/// `Download dedup audit (CSV)` action (Stage 2 result) performs a file-dialog
/// write inline, mirroring `log_pane::show` mutating `LogPaneState` flags.
/// What the user clicked in the Stage 3 result Cache-data block. Collected
/// during the immutable-borrow render pass, then translated to `LogPaneState`
/// request flags after the scroll closure (so `stage3_result::show` drains them
/// into the existing confirm-modal / re-run flow the same frame).
enum ResultCacheAction {
    None,
    RefreshCatalogue,
    RefreshPubchem,
    RefreshKegg,
    Rerun,
}

pub fn show(ui: &mut egui::Ui, app: &mut App) {
    let mut dedup_download = false;
    // Cache-refresh intents collected inside the scroll closure (which holds an
    // immutable borrow of `app`), acted on with `&mut app` afterwards — same
    // pattern as `dedup_download`.
    let mut setup_module_refresh = false;
    let mut setup_pathway_refresh = false;
    let mut organisms_refresh = false;
    let mut result_cache_action = ResultCacheAction::None;

    // Organism roster fetched-date + loading state for the Cache-data block
    // (mode-independent, Stage 3 only). While a refresh is in flight the state is
    // `Loading`; show the stashed prior timestamp so the date line persists.
    let (organisms_fetched_at, organisms_loading) = match &app.organisms.state {
        crate::app::OrganismsLoadState::Loaded { fetched_at, .. } => (Some(*fetched_at), false),
        crate::app::OrganismsLoadState::Loading { .. } => (
            app.organisms.refresh_stash.as_ref().map(|c| c.fetched_at),
            true,
        ),
        _ => (None, false),
    };

    // Session settings save / load toolbar — relocated from the Log pane
    // toolbar by `move-settings-buttons-to-data-tab`. Rendered above the
    // scroll area so it stays pinned while the stage summary scrolls. Save is
    // hidden during `Initializing`; Load is enabled only on Stage 1 (greyed
    // elsewhere with a hover hint). Clicks set the same `LogPaneState` flags
    // `App::update()` drains — see the `app-shell` spec.
    let enable_settings_save = !matches!(app.state, AppState::Initializing { .. });
    let enable_settings_load = matches!(app.state, AppState::Stage1Input { .. });
    if enable_settings_save {
        ui.horizontal(|ui| {
            if ui.button("Save settings…").clicked() {
                app.log_ui.settings_save_requested = true;
            }
            let resp = ui.add_enabled(enable_settings_load, egui::Button::new("Load settings…"));
            let resp = if !enable_settings_load {
                resp.on_disabled_hover_text(
                    "Loading settings is only available on the Stage 1 input screen.",
                )
            } else {
                resp
            };
            if resp.clicked() && enable_settings_load {
                app.log_ui.settings_load_requested = true;
            }
        });
        ui.separator();
    }

    ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            let inputs = &app.inputs;
            let settings = &app.settings;
            let cache = &app.cache;
            let mapping = inputs.mapping.as_ref();
            let (fc, fdr, delta) = (
                settings.fc_threshold,
                settings.fdr_threshold,
                settings.delta_threshold,
            );
            // Only a setup screen carries `modules_fetch`; `Running` does not.
            let module_fetch_in_flight = crate::app::setup_fetch_slots(&app.state)
                .is_some_and(|(_, modules_fetch)| modules_fetch.is_some());
            // The Cache-data block renders on the RUNNING states too, where a
            // refresh click would spawn a fetch no screen can receive.
            let busy = crate::app::is_busy(&app.state);
            match &app.state {
                AppState::Initializing { .. } => {
                    ui.label("Loading…");
                }
                // Nothing is loaded on the route chooser and nothing is loading
                // either, so the panel renders EMPTY. It previously borrowed
                // `Initializing`'s "Loading…" stub, which said something untrue:
                // the roster load has already finished by the time this screen
                // appears, and a reader waiting for it would wait forever.
                AppState::Stage0ChooseAnalysis => {}
                AppState::Stage1Input { .. } => {
                    render_stage1(ui, &inputs.ion_tables, mapping);
                }
                AppState::Stage2DamSetup { .. } | AppState::Stage2DamRunning { .. } => {
                    render_stage2_setup(ui, &inputs.ion_tables, mapping);
                }
                AppState::Stage2DamThreshold { dam_results, .. } => {
                    render_dam_data_block(ui, dam_results);
                    render_dam_slots(
                        ui,
                        &inputs.ion_tables,
                        fc,
                        fdr,
                        delta,
                        dam_results,
                        &mut dedup_download,
                        true,
                    );
                    render_metadata(ui, &inputs.ion_tables, mapping, false);
                }
                AppState::Stage3EnrichSetup { dam_results, .. }
                | AppState::Stage3EnrichRunning {
                    payload: RunningPayload::Enrichment(dam_results),
                    ..
                } => {
                    // Cache data → Enrichment data → DAM data → slots → metadata.
                    render_cache_block_setup(
                        ui,
                        settings,
                        cache,
                        module_fetch_in_flight,
                        organisms_fetched_at,
                        organisms_loading,
                        busy,
                        &mut setup_module_refresh,
                        &mut setup_pathway_refresh,
                        &mut organisms_refresh,
                    );
                    render_enrichment_setup_block(ui, settings, cache);
                    ui.add_space(8.0);
                    render_dam_data_block(ui, dam_results);
                    render_dam_slots(
                        ui,
                        &inputs.ion_tables,
                        fc,
                        fdr,
                        delta,
                        dam_results,
                        &mut dedup_download,
                        false,
                    );
                    render_metadata(ui, &inputs.ion_tables, mapping, false);
                }
                AppState::Stage3EnrichResult {
                    dam_results,
                    module_retention,
                    enrichment_result,
                    dual_mode_breakdown,
                    funnel,
                    pubchem_time_span,
                    kegg_conv_time_span,
                    refresh_state,
                    ..
                } => {
                    let busy = crate::app::is_busy(&app.state);
                    render_cache_block_result(
                        ui,
                        settings.analysis_mode,
                        cache.species_kegg.as_ref(),
                        module_retention.as_ref(),
                        *pubchem_time_span,
                        *kegg_conv_time_span,
                        refresh_state,
                        busy,
                        organisms_fetched_at,
                        organisms_loading,
                        &mut result_cache_action,
                        &mut organisms_refresh,
                    );
                    render_enrichment_result_block(
                        ui,
                        settings.analysis_mode,
                        settings.enrichment_fdr_threshold,
                        settings.kegg_species.as_deref(),
                        module_retention.as_ref(),
                        enrichment_result,
                        dual_mode_breakdown.as_ref(),
                        funnel,
                    );
                    ui.add_space(8.0);
                    render_dam_data_block(ui, dam_results);
                    render_dam_slots(
                        ui,
                        &inputs.ion_tables,
                        fc,
                        fdr,
                        delta,
                        dam_results,
                        &mut dedup_download,
                        false,
                    );
                    render_metadata(ui, &inputs.ion_tables, mapping, false);
                }
                // ── Coverage route ──
                //
                // No `DAM data` block on any of these: that block reports a
                // numerator/denominator comparison, a statistical method, and an
                // up/down/ns tally, none of which exist on a route that runs no
                // differential analysis. Rendering it with zeros would
                // misrepresent the run.
                AppState::Stage2CoverageSetup { .. }
                | AppState::Stage3EnrichRunning {
                    payload: RunningPayload::Coverage,
                    ..
                } => {
                    render_cache_block_setup(
                        ui,
                        settings,
                        cache,
                        module_fetch_in_flight,
                        organisms_fetched_at,
                        organisms_loading,
                        busy,
                        &mut setup_module_refresh,
                        &mut setup_pathway_refresh,
                        &mut organisms_refresh,
                    );
                    render_coverage_target_block(ui, settings, cache);
                    ui.add_space(8.0);
                    render_stage2_setup(ui, &inputs.ion_tables, mapping);
                }
                AppState::Stage3CoverageResult {
                    coverage_result,
                    funnel,
                    module_retention,
                    mode_partition,
                    dedup_reports,
                    pubchem_time_span,
                    kegg_conv_time_span,
                    ..
                } => {
                    render_coverage_cache_block(
                        ui,
                        settings.analysis_mode,
                        cache.species_kegg.as_ref(),
                        module_retention.as_ref(),
                        *pubchem_time_span,
                        *kegg_conv_time_span,
                        organisms_fetched_at,
                        organisms_loading,
                        &mut organisms_refresh,
                    );
                    render_coverage_result_block(
                        ui,
                        settings,
                        cache,
                        coverage_result,
                        funnel,
                        mode_partition.as_ref(),
                    );
                    ui.add_space(8.0);
                    render_coverage_slots(
                        ui,
                        &inputs.ion_tables,
                        dedup_reports,
                        &mut dedup_download,
                    );
                    render_metadata(ui, &inputs.ion_tables, mapping, false);
                }
            }
        });

    // Borrows on `app` are released; perform the `&mut app` actions outside the
    // scroll closure so they do not overlap the immutable reads above.
    if setup_pathway_refresh {
        crate::ui::stage3_setup::handle_species_refresh(app);
    }
    if setup_module_refresh
        && let (Some(level), Some(group), Some(org_codes)) = (
            app.settings.organism_group_level,
            app.settings.organism_group.clone(),
            app.cache.group_org_codes.clone(),
        )
    {
        crate::ui::stage3_setup::spawn_modules_fetch(app, level, group, org_codes, true);
    }
    match result_cache_action {
        // The Stage 3 result refresh / re-run buttons set request flags that
        // `stage3_result::show` drains into the existing confirm-modal flow
        // (the Data tab renders before the central panel each frame).
        ResultCacheAction::RefreshCatalogue => app.log_ui.refresh_catalogue_requested = true,
        ResultCacheAction::RefreshPubchem => app.log_ui.refresh_pubchem_requested = true,
        ResultCacheAction::RefreshKegg => app.log_ui.refresh_kegg_conv_requested = true,
        ResultCacheAction::Rerun => app.log_ui.rerun_enrichment_requested = true,
        ResultCacheAction::None => {}
    }
    // The organism-roster refresh button (setup OR result) sets a flag drained by
    // `stage3_setup::show` / `stage3_result::show`, which open the Stage-3-local
    // refresh confirm (NOT an App-level modal — see the `app-shell` spec).
    if organisms_refresh {
        app.log_ui.organisms_refresh_requested = true;
    }
    if dedup_download {
        // Same writer on both routes; only the source of the reports differs.
        match &app.state {
            AppState::Stage3CoverageResult { dedup_reports, .. } => {
                let reports: Vec<&crate::dedup::DedupReport> = dedup_reports.iter().collect();
                crate::ui::stage2_threshold::write_dedup_audit(
                    &reports,
                    app.inputs.ion_tables.as_slice(),
                );
            }
            _ => crate::ui::stage2_threshold::download_dedup_audit_csv(app),
        }
    }
}

// ── Stage 1 ──────────────────────────────────────────────────────────────

fn render_stage1(ui: &mut egui::Ui, ion_tables: &[IonModeTable], mapping: Option<&GroupMapping>) {
    for (i, it) in ion_tables.iter().enumerate() {
        slot_header(ui, i, it.mode);
        kv_line(
            ui,
            &format!(
                "Raw input: {} features, {} sample columns",
                it.table.features.len(),
                it.table.sample_cols.len()
            ),
        );
        ui.add_space(6.0);
    }
    render_metadata(ui, ion_tables, mapping, true);
}

// ── Stage 2 setup ────────────────────────────────────────────────────────

fn render_stage2_setup(
    ui: &mut egui::Ui,
    ion_tables: &[IonModeTable],
    mapping: Option<&GroupMapping>,
) {
    for (i, it) in ion_tables.iter().enumerate() {
        slot_header(ui, i, it.mode);
        raw_input_with_split(ui, &it.table);
        ui.add_space(6.0);
    }
    render_metadata(ui, ion_tables, mapping, false);
}

// ── per-slot DAM blocks (Stage 2 result + Stage 3) ─────────────────────────

#[allow(clippy::too_many_arguments)]
fn render_dam_slots(
    ui: &mut egui::Ui,
    ion_tables: &[IonModeTable],
    fc: f64,
    fdr: f64,
    delta: f64,
    dam_results: &[DamResult],
    dedup_download: &mut bool,
    show_dedup_button: bool,
) {
    for (i, it) in ion_tables.iter().enumerate() {
        slot_header(ui, i, it.mode);
        raw_input_with_split(ui, &it.table);
        if let Some(dam) = dam_results.get(i) {
            if let Some(report) = dam.dedup_report.as_ref() {
                kv_line(
                    ui,
                    &format!(
                        "Dedupe: {} dup-losers dropped, {} null-InChIKey passed through",
                        report.dropped.len(),
                        report.null_inchikey_passthrough
                    ),
                );
            }
            kv_line(ui, &format!("Pre-filter: {} features dropped", dam.skipped));
            kv_line(ui, &format!("DAM input: {} features", dam.features.len()));
            let (up, down) = count_trends(dam, fc, fdr, delta);
            let ns = dam.features.len().saturating_sub(up + down);
            kv_line(ui, &dam_line_str(up, down, ns));
        }
        ui.add_space(6.0);
    }

    // Relocated dedup-audit download (present only when at least one mode
    // produced a dedup report). The click sets a flag drained by `show`.
    if show_dedup_button
        && dam_results.iter().any(|r| r.dedup_report.is_some())
        && ui.button("Download dedup audit (CSV)").clicked()
    {
        *dedup_download = true;
    }
    ui.add_space(6.0);
}

// ── Coverage-route blocks ──────────────────────────────────────────────────

/// `Coverage data` on the coverage SETUP screens: the analysis target, and the
/// group selection when a `.csv` was supplied.
///
/// Named `Coverage data`, not `Enrichment data`: this route computes no
/// enrichment, and a header naming one would be the first thing a reader of a
/// bug report saw.
fn render_coverage_target_block(
    ui: &mut egui::Ui,
    settings: &SessionSettings,
    cache: &SessionCache,
) {
    section_header(ui, "Coverage data");
    kv_line(
        ui,
        &format!("Analysis mode: {}", mode_name(settings.analysis_mode)),
    );
    match settings.analysis_mode {
        AnalysisMode::Pathway => {
            kv_line(
                ui,
                &format!(
                    "Species: {}",
                    settings.kegg_species.as_deref().unwrap_or("—")
                ),
            );
            if let Some(sk) = &cache.species_kegg {
                kv_line(ui, &format!("Pathways in catalogue: {}", sk.pathways.len()));
            }
        }
        AnalysisMode::Module => {
            kv_line(
                ui,
                &format!(
                    "Group: {} (level {})",
                    settings.organism_group.as_deref().unwrap_or("—"),
                    settings
                        .organism_group_level
                        .map(|l| l.to_string())
                        .unwrap_or_else(|| "—".into())
                ),
            );
        }
    }
    // Group selection: the one input on this route that can silently REMOVE
    // features, so a bug report has to record it. Group names are the user's
    // own and already on screen, so they are safe in the bundle.
    match settings.coverage_selected_groups.as_ref() {
        Some(groups) if groups.is_empty() => {
            kv_line_colored(ui, "Sample groups: none selected", theme::WARNING);
        }
        Some(groups) => {
            kv_line(ui, &format!("Sample groups: {}", groups.join(", ")));
            kv_line(
                ui,
                &format!(
                    "Detected in >= {:.0}% of a group's samples",
                    settings.coverage_presence_threshold * 100.0
                ),
            );
        }
        None => {}
    }
    ui.add_space(6.0);
}

/// `Coverage data` on the RESULT screen: the target, the provenance funnel, and
/// — in dual mode — the per-mode partition.
///
/// The funnel has **no foreground branch**: there is no foreground on this
/// route, so no `foreground_*` value is rendered and no label uses the words
/// "foreground", "significant", or "universe".
fn render_coverage_result_block(
    ui: &mut egui::Ui,
    settings: &SessionSettings,
    cache: &SessionCache,
    result: &crate::coverage::CoverageResult,
    funnel: &crate::app::CoverageFunnel,
    partition: Option<&crate::stage3::CoverageModePartition>,
) {
    render_coverage_target_block(ui, settings, cache);

    section_header(ui, "Coverage funnel");
    kv_line(ui, &format!("Raw features: {}", funnel.raw_features));
    if let Some(n) = funnel.in_selected_groups {
        kv_line(ui, &format!("In selected groups: {n}"));
    }
    kv_line(ui, &format!("After deduplication: {}", funnel.after_dedup));
    kv_line(
        ui,
        &format!("Distinct InChIKeys: {}", funnel.detected_inchikeys),
    );
    kv_line(ui, &format!("Distinct CIDs: {}", funnel.detected_cids));
    kv_line(ui, &format!("KEGG compounds: {}", result.detected_total));
    kv_line(
        ui,
        &format!("In at least one entry: {}", result.detected_in_entries),
    );
    kv_line(
        ui,
        &format!(
            "Entries: {} ({} with no KEGG compounds)",
            result.entries_total, result.entries_without_compounds
        ),
    );

    // The Data tab is the SOLE surface for this partition — the results table
    // deliberately renders no per-mode columns.
    if let Some(p) = partition {
        section_header(ui, "Ionization modes");
        kv_line(ui, &format!("POS only: {}", p.pos_only));
        kv_line(ui, &format!("NEG only: {}", p.neg_only));
        kv_line(ui, &format!("In both: {}", p.in_both));
    }
    ui.add_space(6.0);
}

/// Per-slot MS-DIAL blocks for the coverage result, with the `Dedupe:` line and
/// the relocated audit download.
///
/// On this route that line carries extra weight: deduplication provably cannot
/// change any reported coverage number, so it and the audit are the ONLY
/// surfaces on which its effect is observable at all. Without them,
/// "the deduplication controls are inspectable" would mean "save a file and
/// open it in a spreadsheet".
fn render_coverage_slots(
    ui: &mut egui::Ui,
    ion_tables: &[IonModeTable],
    dedup_reports: &[crate::dedup::DedupReport],
    dedup_download: &mut bool,
) {
    for (i, it) in ion_tables.iter().enumerate() {
        slot_header(ui, i, it.mode);
        raw_input_with_split(ui, &it.table);
        if let Some(report) = dedup_reports.get(i) {
            kv_line(
                ui,
                &format!(
                    "Dedupe: {} dup-losers dropped, {} null-InChIKey passed through",
                    report.dropped.len(),
                    report.null_inchikey_passthrough
                ),
            );
        }
        ui.add_space(6.0);
    }
    // Absent when the run was performed with `dedup_enabled = false` — there is
    // no report to export.
    if !dedup_reports.is_empty() && ui.button("Download dedup audit (CSV)").clicked() {
        *dedup_download = true;
    }
    ui.add_space(6.0);
}

/// `Cache data` on the coverage result screen.
///
/// Same fetched-date lines as the enrichment result's block, WITHOUT any
/// refresh or re-run button: this route offers no PubChem/KEGG refresh action,
/// which is also why its state carries no `refresh_state`. Rendering the
/// buttons anyway would advertise an action the route cannot perform.
#[allow(clippy::too_many_arguments)]
fn render_coverage_cache_block(
    ui: &mut egui::Ui,
    mode: AnalysisMode,
    species: Option<&SpeciesKegg>,
    module_retention: Option<&ModuleRetention>,
    pubchem_time_span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>,
    kegg_conv_time_span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>,
    organisms_fetched_at: Option<DateTime<Utc>>,
    organisms_loading: bool,
    out_organisms_refresh: &mut bool,
) {
    section_header(ui, "Cache data");
    render_organism_cache_row(
        ui,
        organisms_fetched_at,
        organisms_loading,
        out_organisms_refresh,
    );
    match mode {
        AnalysisMode::Pathway => {
            if let Some(sk) = species {
                kv_line_colored(
                    ui,
                    &format!(
                        "KEGG pathways ({}): {}",
                        sk.code,
                        sk.fetched_at.format("%Y-%m-%d %H:%M UTC")
                    ),
                    theme::TEXT_SECONDARY,
                );
            }
        }
        AnalysisMode::Module => {
            if let Some(r) = module_retention {
                kv_line_colored(
                    ui,
                    &format!(
                        "KEGG modules fetched date: {} -> {}",
                        r.oldest_fetched_at.format("%Y-%m-%d"),
                        r.newest_fetched_at.format("%Y-%m-%d"),
                    ),
                    theme::TEXT_SECONDARY,
                );
            }
        }
    }
    kv_line_colored(
        ui,
        &pubchem_span_str(pubchem_time_span),
        theme::TEXT_SECONDARY,
    );
    kv_line_colored(
        ui,
        &kegg_conv_span_str(kegg_conv_time_span),
        theme::TEXT_SECONDARY,
    );
    ui.add_space(8.0);
}

// ── Stage 3 enrichment-data block (setup) ──────────────────────────────────

/// Title-Case mode name for the `Analysis mode:` line.
fn mode_name(m: AnalysisMode) -> &'static str {
    match m {
        AnalysisMode::Pathway => "Pathway",
        AnalysisMode::Module => "Module",
    }
}

fn render_enrichment_setup_block(
    ui: &mut egui::Ui,
    settings: &SessionSettings,
    cache: &SessionCache,
) {
    section_header(ui, "Enrichment data");
    kv_line(
        ui,
        &format!("Analysis mode: {}", mode_name(settings.analysis_mode)),
    );
    match settings.analysis_mode {
        AnalysisMode::Pathway => {
            if let Some(code) = settings.kegg_species.as_deref() {
                kv_line(ui, &format!("Species: {code}"));
            }
        }
        AnalysisMode::Module => {
            if let (Some(group), Some(level), Some(orgs)) = (
                settings.organism_group.as_deref(),
                settings.organism_group_level,
                cache.group_org_codes.as_ref(),
            ) {
                kv_line(
                    ui,
                    &format!("Group: {group} (Level {level}, {} organisms)", orgs.len()),
                );
            }
        }
    }
    match settings.analysis_mode {
        AnalysisMode::Pathway => match &cache.species_kegg {
            Some(sk) => {
                kv_line(ui, &format!("Pathways fetched: {}", sk.pathways.len()));
            }
            None => {
                kv_line(ui, "Pathways: not fetched yet");
            }
        },
        AnalysisMode::Module => match &cache.modules_pack {
            Some(pack) => {
                kv_line(ui, &format!("Modules fetched: {}", pack.modules.len()));
                if let Some(orgs) = &cache.group_org_codes {
                    // Memoize the ~573-module Group-overlap scan across frames via
                    // egui's cross-frame `Context::data`, keyed on every input that
                    // can change it; recompute only on a key miss.
                    let key = module_group_memo_key(
                        settings.min_group_overlap,
                        settings.organism_group.as_deref(),
                        settings.organism_group_level,
                        pack.modules.len(),
                    );
                    let id = egui::Id::new("data_tab.module_group_counts");
                    let cached = ui
                        .ctx()
                        .data(|d| d.get_temp::<(ModuleGroupKey, (usize, usize))>(id));
                    let (in_group, with_compounds) = match cached {
                        Some((k, v)) if k == key => v,
                        _ => {
                            let v = module_group_counts(pack, orgs, settings.min_group_overlap);
                            ui.ctx().data_mut(|d| d.insert_temp(id, (key, v)));
                            v
                        }
                    };
                    kv_line(ui, &format!("In selected Group: {in_group}"));
                    kv_line(ui, &with_compound_line(in_group, in_group - with_compounds));
                }
            }
            None => {
                kv_line(ui, "Modules: not fetched yet");
            }
        },
    }
}

// ── Stage 3 enrichment-data block (result) ─────────────────────────────────

#[allow(clippy::too_many_arguments)]
fn render_enrichment_result_block(
    ui: &mut egui::Ui,
    mode: AnalysisMode,
    fdr_threshold: f64,
    species: Option<&str>,
    module_retention: Option<&ModuleRetention>,
    result: &EnrichmentResult,
    breakdown: Option<&DualModeBreakdown>,
    funnel: &Stage3Funnel,
) {
    let entry_noun = format!("{}s", mode.entry_label_singular());
    let tested = result.rows.len();
    let min_entry = result.min_entry_size;

    section_header(ui, "Enrichment data");
    kv_line(ui, &format!("Analysis mode: {}", mode_name(mode)));
    match mode {
        AnalysisMode::Module => {
            if let Some(r) = module_retention {
                kv_line(
                    ui,
                    &format!(
                        "Group: {} (Level {}, {} organisms)",
                        r.group_name, r.group_level, r.group_org_count
                    ),
                );
            }
        }
        AnalysisMode::Pathway => {
            if let Some(code) = species {
                kv_line(ui, &format!("Species: {code}"));
            }
        }
    }
    match mode {
        AnalysisMode::Module => {
            if let Some(r) = module_retention {
                let empty = result.empty_compound_count;
                kv_line(ui, &format!("Modules fetched: {}", r.total_modules));
                kv_line(ui, &format!("In selected Group: {}", r.retained_modules));
                kv_line(ui, &with_compound_line(r.retained_modules, empty));
                kv_line(
                    ui,
                    &format!("Tested: {tested} (≥ {min_entry} compounds in universe)"),
                );
            }
        }
        AnalysisMode::Pathway => {
            let fetched = tested + result.entries_dropped_by_min_entry_size;
            kv_line(ui, &format!("Pathways fetched: {fetched}"));
            kv_line(
                ui,
                &format!("Tested: {tested} / {fetched} (≥ {min_entry} compounds in universe)"),
            );
        }
    }

    // Universe funnel (all tested features mapped to KEGG — the measurable metabolome).
    ui.add_space(6.0);
    ui.label(
        RichText::new("Universe — all tested features (measurable metabolome)")
            .strong()
            .color(theme::HEADING),
    );
    ui.label(format!("{} InChIKeys", funnel.detected_inchikeys));
    ui.label(format!("->{} PubChem CIDs", funnel.detected_cids));
    let universe_partition =
        breakdown.map(|b| (b.universe_pos_only, b.universe_neg_only, b.universe_in_both));
    ui.label(universe_kegg_line(result.universe_size, universe_partition));

    // Foreground funnel (significant subset, framed as "of which").
    ui.add_space(6.0);
    ui.label(
        RichText::new("Foreground — significant features (active direction)")
            .strong()
            .color(theme::HEADING),
    );
    ui.label(format!(
        "of which {} InChIKeys",
        funnel.foreground_inchikeys
    ));
    ui.label(format!("->{} PubChem CIDs", funnel.foreground_cids));
    let fg_partition = breakdown.map(|b| {
        (
            b.foreground_pos_only,
            b.foreground_neg_only,
            b.foreground_agree_both,
            b.foreground_excluded_conflict,
        )
    });
    ui.label(foreground_kegg_line(result.dam_cpd_size, fg_partition));

    // Coverage + significance.
    ui.add_space(6.0);
    ui.label(format!(
        "Detected in tested {entry_noun}: {} KEGG cpds",
        funnel.detected_in_entries
    ));
    let sig = result
        .rows
        .iter()
        .filter(|r| r.fdr < fdr_threshold && r.displayed)
        .count();
    ui.label(format!(
        "Significant {entry_noun}: {sig} (FDR < {fdr_threshold})"
    ));

    // Dual-mode K-source hint.
    if let Some(b) = breakdown
        && let Some(hint) = k_source_hint(
            b.foreground_pos_only,
            b.foreground_neg_only,
            b.foreground_agree_both,
        )
    {
        ui.colored_label(theme::WARNING, hint);
    }
}

// ── Stage 3 cache-data block (relocated from the screen bodies) ────────────

/// `Cache data` block for Stage 3 **setup**: the cache fetched-date line + the
/// `Refresh` button, relocated from the setup body. Records the click intent in
/// the `out_*` flags (the caller acts on `&mut app` after the scroll closure).
/// The mode-independent `KEGG organism list: <fetched_at>` line + `Refresh KEGG
/// organism list` button, shared by the setup and result Cache-data blocks. The
/// roster underlies both selectors, so it renders in both Pathway and Module
/// mode. `organisms_fetched_at` is the loaded roster's timestamp (or the stashed
/// prior timestamp while a refresh is in flight); when `None` (Idle / failed
/// eager load) the row is omitted. The button is disabled while a refresh is in
/// flight. Records the click intent in `out_organisms_refresh`.
fn render_organism_cache_row(
    ui: &mut egui::Ui,
    organisms_fetched_at: Option<DateTime<Utc>>,
    organisms_loading: bool,
    out_organisms_refresh: &mut bool,
) {
    let Some(ts) = organisms_fetched_at else {
        return;
    };
    kv_line_colored(
        ui,
        &format!("KEGG organism list: {}", ts.format("%Y-%m-%d %H:%M UTC")),
        theme::TEXT_SECONDARY,
    );
    let resp = ui
        .add_enabled(
            !organisms_loading,
            egui::Button::new("Refresh KEGG organism list"),
        )
        .on_disabled_hover_text("Refreshing organism list…");
    if resp.clicked() {
        *out_organisms_refresh = true;
    }
}

#[allow(clippy::too_many_arguments)]
fn render_cache_block_setup(
    ui: &mut egui::Ui,
    settings: &SessionSettings,
    cache: &SessionCache,
    module_fetch_in_flight: bool,
    organisms_fetched_at: Option<DateTime<Utc>>,
    organisms_loading: bool,
    // The shared busy predicate for the current state. This block renders on
    // the RUNNING states of both routes as well as the setup screens, and a
    // running state owns neither an in-flight fetch slot nor an `error` field
    // — so a refused refresh there could not even report itself. Disabling is
    // the only intelligible answer; see the `data-summary-panel` spec.
    busy: bool,
    out_module_refresh: &mut bool,
    out_pathway_refresh: &mut bool,
    out_organisms_refresh: &mut bool,
) {
    section_header(ui, "Cache data");
    // Organism roster (mode-independent) sits at the TOP of the Cache-data
    // block, above the mode-specific pathway/module entries.
    render_organism_cache_row(
        ui,
        organisms_fetched_at,
        organisms_loading,
        out_organisms_refresh,
    );
    match settings.analysis_mode {
        AnalysisMode::Pathway => {
            if let Some(sk) = &cache.species_kegg {
                kv_line_colored(
                    ui,
                    &format!(
                        "KEGG pathways ({}): {}",
                        sk.code,
                        sk.fetched_at.format("%Y-%m-%d %H:%M UTC")
                    ),
                    theme::TEXT_SECONDARY,
                );
                if ui
                    .add_enabled(!busy, egui::Button::new("Refresh KEGG pathway cache"))
                    .on_disabled_hover_text(
                        "Unavailable while an analysis or fetch is running — \
                         this screen cannot receive the result.",
                    )
                    .clicked()
                {
                    *out_pathway_refresh = true;
                }
            }
        }
        AnalysisMode::Module => {
            if let Some(pack) = &cache.modules_pack {
                if let (Some(oldest), Some(newest)) = (
                    pack.modules.values().map(|m| m.fetched_at).min(),
                    pack.modules.values().map(|m| m.fetched_at).max(),
                ) {
                    kv_line_colored(
                        ui,
                        &format!(
                            "KEGG modules fetched date: {} -> {}",
                            oldest.format("%Y-%m-%d"),
                            newest.format("%Y-%m-%d"),
                        ),
                        theme::TEXT_SECONDARY,
                    );
                }
                // Shown only when a Group is selected + cached and no fetch is in
                // flight (mirrors the prior setup-body `refresh_ready` gate).
                let refresh_ready = settings.organism_group_level.is_some()
                    && settings.organism_group.is_some()
                    && cache.group_org_codes.is_some()
                    && !module_fetch_in_flight;
                if refresh_ready
                    && ui
                        .add_enabled(!busy, egui::Button::new("Refresh KEGG module cache"))
                        .on_disabled_hover_text(
                            "Unavailable while an analysis or fetch is running — \
                             this screen cannot receive the result.",
                        )
                        .clicked()
                {
                    *out_module_refresh = true;
                }
            }
        }
    }
    ui.add_space(8.0);
}

/// `Cache data` block for Stage 3 **result**: the mode-aware cache date, the
/// PubChem + KEGG-conv fetched-date spans with their `Refresh` buttons (+ inline
/// progress), the `Re-run enrichment` button, and the re-run reminder note. All
/// relocated from the result body. Records intent in `out`.
#[allow(clippy::too_many_arguments)]
fn render_cache_block_result(
    ui: &mut egui::Ui,
    mode: AnalysisMode,
    species: Option<&SpeciesKegg>,
    module_retention: Option<&ModuleRetention>,
    pubchem_time_span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>,
    kegg_conv_time_span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>,
    refresh_state: &RefreshState,
    busy: bool,
    organisms_fetched_at: Option<DateTime<Utc>>,
    organisms_loading: bool,
    out: &mut ResultCacheAction,
    out_organisms_refresh: &mut bool,
) {
    section_header(ui, "Cache data");

    // Organism roster (mode-independent) sits at the TOP of the Cache-data
    // block, above the mode-specific catalogue / PubChem / KEGG-conv entries.
    render_organism_cache_row(
        ui,
        organisms_fetched_at,
        organisms_loading,
        out_organisms_refresh,
    );

    // Mode-aware fetched-date line.
    match mode {
        AnalysisMode::Pathway => {
            if let Some(sk) = species {
                kv_line_colored(
                    ui,
                    &format!(
                        "KEGG pathways ({}): {}",
                        sk.code,
                        sk.fetched_at.format("%Y-%m-%d %H:%M UTC")
                    ),
                    theme::TEXT_SECONDARY,
                );
            }
        }
        AnalysisMode::Module => {
            if let Some(r) = module_retention {
                kv_line_colored(
                    ui,
                    &format!(
                        "KEGG modules fetched date: {} -> {}",
                        r.oldest_fetched_at.format("%Y-%m-%d"),
                        r.newest_fetched_at.format("%Y-%m-%d"),
                    ),
                    theme::TEXT_SECONDARY,
                );
            }
        }
    }

    // Catalogue (module/pathway) refresh. The result state lacks the
    // catalogue-fetch progress infra, so clicking this navigates back to Stage 3
    // setup and re-fetches there (see `stage3_result::drain_cache_actions`).
    let catalogue_label = match mode {
        AnalysisMode::Pathway => "Refresh KEGG pathway cache",
        AnalysisMode::Module => "Refresh KEGG module cache",
    };
    if ui
        .add_enabled(!busy, egui::Button::new(catalogue_label))
        .clicked()
    {
        *out = ResultCacheAction::RefreshCatalogue;
    }

    // PubChem line + refresh + in-flight progress.
    kv_line_colored(
        ui,
        &pubchem_span_str(pubchem_time_span),
        theme::TEXT_SECONDARY,
    );
    if ui
        .add_enabled(!busy, egui::Button::new("Refresh PubChem cache"))
        .clicked()
    {
        *out = ResultCacheAction::RefreshPubchem;
    }
    if let RefreshState::RefreshingPubchem {
        completed, total, ..
    } = refresh_state
    {
        refresh_progress_row(ui, "Refreshing PubChem cache…", *completed, *total);
    }

    // KEGG conv line + refresh + in-flight progress.
    kv_line_colored(
        ui,
        &kegg_conv_span_str(kegg_conv_time_span),
        theme::TEXT_SECONDARY,
    );
    if ui
        .add_enabled(!busy, egui::Button::new("Refresh KEGG conv cache"))
        .clicked()
    {
        *out = ResultCacheAction::RefreshKegg;
    }
    if let RefreshState::RefreshingKegg {
        completed, total, ..
    } = refresh_state
    {
        refresh_progress_row(ui, "Refreshing KEGG conv cache…", *completed, *total);
    }

    // Re-run reminder note with the `Re-run enrichment` Secondary button
    // embedded inline (the default egui button style IS the §2 Secondary
    // component per `theme::install`; Primary is opt-in via `widgets`). The
    // bracketed mockup token `[Re-run enrichment]` is this button.
    ui.horizontal(|ui| {
        ui.label(RichText::new("Please").color(theme::TEXT_SECONDARY));
        if ui
            .add_enabled(!busy, egui::Button::new("Re-run enrichment"))
            .clicked()
        {
            *out = ResultCacheAction::Rerun;
        }
        ui.label(RichText::new("after refresh cache.").color(theme::TEXT_SECONDARY));
    });
    ui.add_space(8.0);
}

/// `PubChem CIDs fetched date: <min> -> <max> (<n> entries used)`, or the
/// `(no entries used)` form. Pure builder (unit-tested).
fn pubchem_span_str(span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>) -> String {
    match span {
        Some((min, max, n)) => format!(
            "PubChem CIDs fetched date: {} -> {} ({n} entries used)",
            min.format("%Y-%m-%d"),
            max.format("%Y-%m-%d"),
        ),
        None => "PubChem CIDs fetched date: (no entries used)".to_string(),
    }
}

/// `KEGG conv fetched date: <min> -> <max> (<n> entries used)`, or the
/// `(no entries used)` form. Pure builder (unit-tested).
fn kegg_conv_span_str(span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>) -> String {
    match span {
        Some((min, max, n)) => format!(
            "KEGG conv fetched date: {} -> {} ({n} entries used)",
            min.format("%Y-%m-%d"),
            max.format("%Y-%m-%d"),
        ),
        None => "KEGG conv fetched date: (no entries used)".to_string(),
    }
}

/// In-flight refresh progress row (label + bar), shown under the cache line the
/// refresh pertains to. Mirrors `stage3_result::refresh_progress_row`.
fn refresh_progress_row(ui: &mut egui::Ui, label: &str, completed: usize, total: usize) {
    let frac = if total == 0 {
        0.0
    } else {
        completed as f32 / total as f32
    };
    ui.label(label);
    crate::ui::widgets::progress_bar(
        ui,
        egui::ProgressBar::new(frac)
            .text(format!("{completed} / {total}"))
            .desired_width(420.0),
        theme::SURFACE,
    );
}

// ── shared helpers ──────────────────────────────────────────────────────────

fn slot_header(ui: &mut egui::Ui, idx: usize, mode: IonMode) {
    section_header(
        ui,
        format!("MS-DIAL data slot {} ({})", idx + 1, mode_full(mode)),
    );
}

fn raw_input_with_split(ui: &mut egui::Ui, table: &crate::data::MetabolomicsTable) {
    kv_line(
        ui,
        &raw_input_split_str(table.features.len(), table.annotated_count),
    );
}

/// Cache key for the Module-mode group-overlap memo (egui `Context::data`).
type ModuleGroupKey = (usize, Option<String>, Option<u8>, usize);

/// Build the [`ModuleGroupKey`] capturing every input `module_group_counts`
/// depends on: the overlap threshold, the selected Group identity (name + level
/// deterministically fix the org set via the session-immutable organism
/// roster), and the catalogue-snapshot size (`modules.len()` only changes on an
/// explicit fetch/refresh). Any change to a keyed input is a cache miss.
fn module_group_memo_key(
    min_overlap: usize,
    group: Option<&str>,
    level: Option<u8>,
    n_modules: usize,
) -> ModuleGroupKey {
    (min_overlap, group.map(str::to_owned), level, n_modules)
}

/// Count of in-Group modules and how many of those have a non-empty compound
/// list — the cheap (no-clone) counterpart of `assemble_module_entries` used
/// for the setup-screen funnel.
fn module_group_counts(
    pack: &KeggModulesCache,
    group_orgs: &HashSet<String>,
    min_overlap: usize,
) -> (usize, usize) {
    let mut in_group = 0;
    let mut with_compounds = 0;
    for entry in pack.modules.values() {
        if entry.complete_orgs.intersection(group_orgs).count() >= min_overlap {
            in_group += 1;
            if !entry.compounds.is_empty() {
                with_compounds += 1;
            }
        }
    }
    (in_group, with_compounds)
}

/// Up / down counts under the live thresholds, via the public `classify_trend`
/// seam (`ns` is derived at the call site as `kept − up − down`).
fn count_trends(dam: &DamResult, fc: f64, fdr: f64, delta: f64) -> (usize, usize) {
    let mut up = 0;
    let mut down = 0;
    for feat in &dam.features {
        match classify_trend(feat, fc, fdr, delta, dam.method) {
            Trend::Up => up += 1,
            Trend::Down => down += 1,
            Trend::NotSignificant => {}
        }
    }
    (up, down)
}

fn render_metadata(
    ui: &mut egui::Ui,
    ion_tables: &[IonModeTable],
    mapping: Option<&GroupMapping>,
    show_samples_and_unassigned: bool,
) {
    let Some(mapping) = mapping else {
        return;
    };
    section_header(ui, "Metadata");

    if show_samples_and_unassigned {
        render_coverage(ui, mapping, ion_tables);
    }

    kv_line(ui, "Groups:");
    for group in mapping.groups() {
        if group == UNASSIGNED {
            if show_samples_and_unassigned {
                let count = mapping.samples_in(&group).len();
                ui.colored_label(
                    theme::WARNING,
                    format!("  {group} ({count} sample{})", plural(count)),
                );
            }
            continue;
        }
        let bios = mapping.biosample_count(&group);
        let count = mapping.samples_in(&group).len();
        ui.label(group_row_str(
            &group,
            count,
            bios,
            show_samples_and_unassigned,
        ));
    }
}

fn render_coverage(ui: &mut egui::Ui, mapping: &GroupMapping, ion_tables: &[IonModeTable]) {
    let matched = mapping.assigned_count();
    let total = mapping.groups_in_order().len();
    let coverage = if ion_tables.len() == 2 {
        let count_assigned = |it: &IonModeTable| -> usize {
            it.table
                .sample_cols
                .iter()
                .filter(|s| mapping.group_of(s) != UNASSIGNED)
                .count()
        };
        let n1 = count_assigned(&ion_tables[0]);
        let n2 = count_assigned(&ion_tables[1]);
        let l1 = ion_tables[0].mode;
        let l2 = ion_tables[1].mode;
        format!("Matched {n1} {l1} + {n2} {l2} = {matched} / {total} sample columns from CSV")
    } else {
        format!("Matched {matched} / {total} sample columns from CSV")
    };
    if matched < total {
        ui.colored_label(theme::WARNING, coverage);
    } else {
        ui.label(coverage);
    }
}

// ── pure string builders (unit-tested) ─────────────────────────────────────

/// `Raw input: N features (A annotated · U unknown)` (Stage 2+).
fn raw_input_split_str(total: usize, annotated: usize) -> String {
    let unknown = total - annotated;
    format!("Raw input: {total} features ({annotated} annotated · {unknown} unknown)")
}

/// `DAM: U up, D down, N ns`.
fn dam_line_str(up: usize, down: usize, ns: usize) -> String {
    format!("DAM: {up} up, {down} down, {ns} ns")
}

/// `With compound list: <in_group − empty>  (−<empty> empty)`.
fn with_compound_line(in_group: usize, empty: usize) -> String {
    format!(
        "With compound list: {}  (−{empty} empty)",
        in_group.saturating_sub(empty)
    )
}

/// Universe KEGG-cpd line; dual-mode adds the `(POS-only …)` partition.
fn universe_kegg_line(n: usize, partition: Option<(usize, usize, usize)>) -> String {
    match partition {
        Some((a, b, c)) => {
            format!("->{n} KEGG cpds  (POS-only: {a}; NEG-only: {b}; in both: {c})")
        }
        None => format!("->{n} KEGG cpds"),
    }
}

/// Foreground KEGG-cpd line; dual-mode adds the per-mode + conflict partition.
fn foreground_kegg_line(k: usize, partition: Option<(usize, usize, usize, usize)>) -> String {
    match partition {
        Some((d, e, f, g)) => format!(
            "->{k} KEGG cpds  (sig POS-only: {d}; sig NEG-only: {e}; agree both: {f}; excluded by conflict: {g})"
        ),
        None => format!("->{k} KEGG cpds"),
    }
}

/// Dual-mode K-source hint derived from the foreground partition. A mode is
/// "silent" when it contributed nothing to K (its own-only count AND the
/// agree-both count are both zero).
fn k_source_hint(pos_only: usize, neg_only: usize, agree_both: usize) -> Option<String> {
    let pos_silent = pos_only == 0 && agree_both == 0;
    let neg_silent = neg_only == 0 && agree_both == 0;
    match (pos_silent, neg_silent) {
        (true, true) => Some("K is empty.".to_string()),
        (false, true) => {
            Some("K source: POS only (NEG had 0 sig features in the active direction).".to_string())
        }
        (true, false) => {
            Some("K source: NEG only (POS had 0 sig features in the active direction).".to_string())
        }
        (false, false) => None,
    }
}

/// A non-Unassigned group row. On Stage 1 (`show_samples = true`) the row shows
/// sample counts and, when a biosample column exists, biosample counts. On
/// Stage 2+ it shows biosample counts only (sample counts omitted).
fn group_row_str(
    group: &str,
    samples: usize,
    biosamples: Option<usize>,
    show_samples: bool,
) -> String {
    match (show_samples, biosamples) {
        (true, Some(b)) => format!("  {group} ({samples} samples, {b} biosamples)"),
        (true, None) => format!("  {group} ({samples} sample{})", plural(samples)),
        (false, Some(b)) => format!("  {group} ({b} biosamples)"),
        (false, None) => format!("  {group} ({samples} samples)"),
    }
}

fn plural(n: usize) -> &'static str {
    if n == 1 { "" } else { "s" }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn raw_input_split_str_formats_counts() {
        assert_eq!(
            raw_input_split_str(16956, 12805),
            "Raw input: 16956 features (12805 annotated · 4151 unknown)"
        );
        assert_eq!(
            raw_input_split_str(10, 10),
            "Raw input: 10 features (10 annotated · 0 unknown)"
        );
    }

    #[test]
    fn method_label_uses_short_en_dash_forms() {
        assert_eq!(method_label(DamMethod::Welch), "Welch's t-test");
        assert_eq!(method_label(DamMethod::Student), "Student's t-test");
        // En-dash + " test" (NOT DamMethod::display_name's ASCII "Brunner-Munzel").
        assert_eq!(
            method_label(DamMethod::BrunnerMunzel),
            "Brunner–Munzel test"
        );
    }

    #[test]
    fn cache_span_strs_format_both_branches() {
        use chrono::TimeZone;
        let a = Utc.with_ymd_and_hms(2026, 5, 22, 0, 0, 0).unwrap();
        let b = Utc.with_ymd_and_hms(2026, 5, 27, 0, 0, 0).unwrap();
        assert_eq!(
            pubchem_span_str(Some((a, b, 8812))),
            "PubChem CIDs fetched date: 2026-05-22 -> 2026-05-27 (8812 entries used)"
        );
        assert_eq!(
            pubchem_span_str(None),
            "PubChem CIDs fetched date: (no entries used)"
        );
        assert_eq!(
            kegg_conv_span_str(Some((a, b, 8515))),
            "KEGG conv fetched date: 2026-05-22 -> 2026-05-27 (8515 entries used)"
        );
        assert_eq!(
            kegg_conv_span_str(None),
            "KEGG conv fetched date: (no entries used)"
        );
    }

    #[test]
    fn dam_line_str_formats_tally() {
        assert_eq!(
            dam_line_str(557, 488, 5449),
            "DAM: 557 up, 488 down, 5449 ns"
        );
        assert_eq!(dam_line_str(0, 0, 0), "DAM: 0 up, 0 down, 0 ns");
    }

    #[test]
    fn group_row_str_stage1_with_biosamples() {
        assert_eq!(
            group_row_str("Control", 16, Some(8), true),
            "  Control (16 samples, 8 biosamples)"
        );
    }

    #[test]
    fn group_row_str_stage1_without_biosample_column_pluralizes() {
        assert_eq!(group_row_str("QC", 1, None, true), "  QC (1 sample)");
        assert_eq!(group_row_str("QC", 6, None, true), "  QC (6 samples)");
    }

    #[test]
    fn group_row_str_stage2_shows_biosamples_only() {
        assert_eq!(
            group_row_str("Treatment", 16, Some(8), false),
            "  Treatment (8 biosamples)"
        );
        assert_eq!(
            group_row_str("Treatment", 16, None, false),
            "  Treatment (16 samples)"
        );
    }

    #[test]
    fn mode_full_long_names() {
        assert_eq!(mode_full(IonMode::Positive), "Positive");
        assert_eq!(mode_full(IonMode::Negative), "Negative");
    }

    #[test]
    fn with_compound_line_subtracts_empties() {
        // 158 in group, 5 empty ->153 with compound list.
        assert_eq!(
            with_compound_line(158, 5),
            "With compound list: 153  (−5 empty)"
        );
        assert_eq!(
            with_compound_line(10, 0),
            "With compound list: 10  (−0 empty)"
        );
    }

    #[test]
    fn universe_kegg_line_dual_and_single() {
        assert_eq!(
            universe_kegg_line(453, Some((289, 111, 53))),
            "->453 KEGG cpds  (POS-only: 289; NEG-only: 111; in both: 53)"
        );
        assert_eq!(universe_kegg_line(453, None), "->453 KEGG cpds");
        // Partition sums to N.
        assert_eq!(289 + 111 + 53, 453);
    }

    #[test]
    fn foreground_kegg_line_dual_and_single() {
        assert_eq!(
            foreground_kegg_line(58, Some((40, 18, 0, 0))),
            "->58 KEGG cpds  (sig POS-only: 40; sig NEG-only: 18; agree both: 0; excluded by conflict: 0)"
        );
        assert_eq!(foreground_kegg_line(58, None), "->58 KEGG cpds");
    }

    #[test]
    fn k_source_hint_cases() {
        // Both contribute ->no hint.
        assert_eq!(k_source_hint(40, 18, 0), None);
        assert_eq!(k_source_hint(0, 0, 4), None); // agree-both keeps both live
        // NEG silent (no neg-only, no agree) ->POS only.
        assert_eq!(
            k_source_hint(40, 0, 0).as_deref(),
            Some("K source: POS only (NEG had 0 sig features in the active direction).")
        );
        // POS silent ->NEG only.
        assert_eq!(
            k_source_hint(0, 18, 0).as_deref(),
            Some("K source: NEG only (POS had 0 sig features in the active direction).")
        );
        // Both silent ->K is empty.
        assert_eq!(k_source_hint(0, 0, 0).as_deref(), Some("K is empty."));
    }

    #[test]
    fn module_group_memo_key_changes_with_each_input() {
        let base = module_group_memo_key(1, Some("Animals"), Some(2), 573);
        assert_eq!(
            base,
            module_group_memo_key(1, Some("Animals"), Some(2), 573)
        );
        assert_ne!(
            base,
            module_group_memo_key(2, Some("Animals"), Some(2), 573)
        ); // overlap
        assert_ne!(base, module_group_memo_key(1, Some("Plants"), Some(2), 573)); // group name
        assert_ne!(
            base,
            module_group_memo_key(1, Some("Animals"), Some(3), 573)
        ); // level
        assert_ne!(
            base,
            module_group_memo_key(1, Some("Animals"), Some(2), 574)
        ); // catalogue size
    }
}
