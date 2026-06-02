//! Stage 3 setup screen. Controls: Analysis Mode toggle, species/group
//! selector, inline KEGG fetch progress, direction selector, top N,
//! enrichment FDR threshold, FDR correction radio, minimum hit count,
//! Run Enrichment button.
//!
//! Per `reorder-gui-and-move-mode-to-stage3` the Mode toggle and the
//! mode-specific selector both live on this screen (no longer Stage 1),
//! and a fetch triggered by the user picking a species / Group writes
//! into `Stage3EnrichSetup.kegg_fetch` / `modules_fetch` rather than
//! transitioning to a dedicated fetching variant. Toggling mode while a
//! fetch is in flight does NOT cancel the in-flight fetch — both
//! `<mode>_fetch` slots can coexist (D2 + D6).

use egui::RichText;
use std::sync::mpsc;
use tokio::sync::mpsc as tokio_mpsc;
use tracing::{error, info, warn};

use crate::theme;

use crate::app::{
    AnalysisMode, AnalysisPayload, App, AppState, KeggFetchInFlight, ModulesFetchInFlight,
    OrganismsLoadState, SessionCache, SessionSettings,
};
use crate::dam::DamResult;
use crate::enrichment::EnrichmentDirection;
use crate::kegg::{self, KeggCacheScope, KeggEvent, KeggProgress};
use crate::stage3::Stage3Params;
use crate::ui::species_selector::{self, SpeciesSelectorEvent};

#[derive(Debug, Clone, Copy)]
enum Action {
    None,
    Run,
}

pub fn show(ui: &mut egui::Ui, app: &mut App) {
    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
        // Defensive state guard — the central-panel dispatcher only routes here
        // in this state. The DAM-result context line (`<num> vs <den>, <method>`
        // + the up/down/ns tally) moved to the bottom-panel Data tab's `DAM data`
        // + per-slot blocks (`data-summary-panel`), so the setup body opens
        // straight into the controls.
        if !matches!(&app.state, AppState::Stage3EnrichSetup { .. }) {
            return;
        }

        ui.heading(egui::RichText::new("Stage 3 — Enrichment Analysis").color(theme::HEADING));
        ui.add_space(6.0);
        // Back-navigation handled by the global stage stepper (`ui::stepper`).
        ui.add_space(8.0);

        // === Mode toggle ===
        render_mode_toggle(ui, app);
        ui.add_space(6.0);

        // === Mode-aware selector ===
        match app.settings.analysis_mode {
            AnalysisMode::Pathway => render_species_selector(ui, app),
            AnalysisMode::Module => render_organism_group_selector(ui, app),
        }
        ui.add_space(6.0);

        // === Inline KEGG fetch progress strip (active mode only) ===
        render_inline_fetch_progress(ui, app);
        ui.add_space(8.0);

        // === Direction / FDR / Min hit count ===
        // Top N moved to Stage 3 result per smoke-test feedback — the user
        // adjusts it there to iterate on the dot plot without navigating
        // back. See `src/ui/stage3_result.rs`.
        let settings = &mut app.settings;
        ui.label("Include DAM features with direction:");
        ui.radio_value(
            &mut settings.direction,
            EnrichmentDirection::Both,
            "Both (Up + Down)",
        );
        ui.radio_value(&mut settings.direction, EnrichmentDirection::Up, "Up only");
        ui.radio_value(
            &mut settings.direction,
            EnrichmentDirection::Down,
            "Down only",
        );
        ui.add_space(8.0);

        // Pre-FDR entry-size filter (position 6 per stage3-ui spec — adjacent to
        // Direction in the "what to test" cluster).
        ui.horizontal(|ui| {
            ui.label(format!(
                "Minimum number of compounds detected in a {}:",
                settings.analysis_mode.entry_label_singular()
            ));
            ui.add(
                egui::DragValue::new(&mut settings.min_entry_size)
                    .speed(1)
                    .range(1..=20),
            )
            .on_hover_text(
                "Drop pathways/modules with fewer than N measurable compounds before FDR. \
                 Reduces multiple-testing penalty.",
            );
        });
        ui.add_space(8.0);

        // `Enrichment FDR threshold` + `Minimum hit count` moved to the Stage 3
        // result screen (`add-bottom-panel-data-tab`) so they can be tuned and
        // re-applied without navigating back. The FDR correction *method* stays
        // here (it changes how q-values are computed, not just the display cut).

        // FDR correction radio. Independent of Stage 2's choice (per add-bh-fdr D3).
        // `None` is exposed here ONLY (Stage 2 setup hides it; an adversarial
        // snapshot is defensively coerced back to BH in `apply_snapshot`).
        // Stage 3 exposes the `None` (NoCorrection) variant → `include_none = true`.
        crate::ui::widgets::fdr_method_radios(ui, &mut settings.enrichment_fdr_method, true);
        ui.label(
            RichText::new(
                "None skips multiple-testing correction entirely — exploratory only, not for publication.",
            )
            .small()
            .color(theme::TEXT),
        );

        ui.add_space(12.0);
        if let AppState::Stage3EnrichSetup { error, .. } = &app.state
            && let Some(e) = error
        {
            ui.colored_label(theme::ERROR, e.clone());
        }

        // === Run button (Back moved to top in Phase 4) ===
        let run_enabled = run_button_enabled(app);
        let mut action = Action::None;
        // Disabled-state hint: explain *why* the button is disabled when a
        // fetch is in flight in the active mode (per stage3-ui spec scenario
        // "Inline progress strip during pathway fetch" and design D2).
        let disabled_hint = if !run_enabled {
            match app.settings.analysis_mode {
                AnalysisMode::Pathway
                    if matches!(
                        &app.state,
                        AppState::Stage3EnrichSetup {
                            kegg_fetch: Some(_),
                            ..
                        }
                    ) =>
                {
                    Some("Waiting for KEGG pathway fetch…")
                }
                AnalysisMode::Module
                    if matches!(
                        &app.state,
                        AppState::Stage3EnrichSetup {
                            modules_fetch: Some(_),
                            ..
                        }
                    ) =>
                {
                    Some("Waiting for KEGG modules fetch…")
                }
                _ => None,
            }
        } else {
            None
        };
        let resp = crate::ui::widgets::primary_button(ui, "Run Enrichment", run_enabled);
        let resp = if let Some(hint) = disabled_hint {
            resp.on_disabled_hover_text(hint)
        } else {
            resp
        };
        if resp.clicked() && run_enabled {
            action = Action::Run;
        }

        match action {
            Action::None => {}
            Action::Run => start_run(app),
        }
        });
}

/// Render the Pathway / Module radio. On change, update
/// `settings.analysis_mode` via the named API. Pathway and Module
/// selections coexist — the API is a near-no-op that only sets the mode
/// (per `reorder-gui-and-move-mode-to-stage3` D3).
fn render_mode_toggle(ui: &mut egui::Ui, app: &mut App) {
    let current_mode = app.settings.analysis_mode;
    let mut new_mode = current_mode;
    ui.horizontal(|ui| {
        ui.label(
            RichText::new("Analysis mode:")
                .strong()
                .color(theme::HEADING),
        );
        ui.radio_value(&mut new_mode, AnalysisMode::Pathway, "Pathway");
        ui.radio_value(&mut new_mode, AnalysisMode::Module, "Module");
    });

    if new_mode != current_mode {
        // The body of this API was reduced to `self.analysis_mode = new_mode;`
        // by Phase 2 of `reorder-gui-and-move-mode-to-stage3`, but the call
        // site stays so spec scenarios remain anchored to a named API.
        app.settings.reset_kegg_selection_for_mode_switch(new_mode);
        app.cache.clear_for_mode_switch(new_mode);
    }
}

fn render_species_selector(ui: &mut egui::Ui, app: &mut App) {
    let (organisms_view, loading, load_error) = match &app.organisms.state {
        OrganismsLoadState::Loaded(v) => (Some(v.as_slice()), false, None),
        OrganismsLoadState::Loading { .. } => (None, true, None),
        OrganismsLoadState::Failed(msg) => (None, false, Some(msg.as_str())),
        OrganismsLoadState::Idle => (None, false, None),
    };
    let current = app.settings.kegg_species.clone();
    let selector_enabled = matches!(&app.state, AppState::Stage3EnrichSetup { .. });

    let event = species_selector::show(
        ui,
        &mut app.species_selector,
        organisms_view,
        loading,
        load_error,
        current.as_deref(),
        selector_enabled,
    );

    // The `Cached <ts>` fetched-date line + the `Refresh KEGG pathway cache`
    // button moved to the bottom-panel Data tab's `Cache data` block
    // (`apply-ui-design-md-tweaks`); the body keeps only the selector.

    match event {
        SpeciesSelectorEvent::None => {}
        SpeciesSelectorEvent::OpenedAndNeedsLoad => {
            app.ensure_organisms_loading();
        }
        SpeciesSelectorEvent::Selected(code) => handle_species_selected(app, code),
    }
}

fn render_organism_group_selector(ui: &mut egui::Ui, app: &mut App) {
    let current_group = app.settings.organism_group.clone();
    let enabled = matches!(&app.state, AppState::Stage3EnrichSetup { .. });

    let index = match crate::kegg::cache::read_organism_group_index() {
        Ok(Some(idx)) => Some(idx),
        Ok(None) => None,
        Err(e) => {
            warn!(error = %e, "failed to read organism group index for Stage 3 setup selector");
            None
        }
    };

    let event = crate::ui::organism_group_selector::show(
        ui,
        &mut app.organism_group_selector,
        index.as_ref(),
        current_group.as_deref(),
        enabled,
    );

    // The `KEGG modules fetched date: …` span line + the `Refresh KEGG module
    // cache` button moved to the bottom-panel Data tab's `Cache data` block
    // (`apply-ui-design-md-tweaks`); the body keeps only the selector. The
    // initial fetch on Group selection still fires here.
    match event {
        crate::ui::organism_group_selector::OrganismGroupSelectorEvent::None => {}
        crate::ui::organism_group_selector::OrganismGroupSelectorEvent::LevelChanged(lvl) => {
            app.settings.organism_group_level = Some(lvl);
            app.settings.organism_group = None;
            app.cache.group_org_codes = None;
            // Switching Level invalidates the org-codes set for any prior
            // Group; the modules cache itself (global to all Groups) stays.
        }
        crate::ui::organism_group_selector::OrganismGroupSelectorEvent::GroupSelected {
            level,
            group,
            org_codes,
        } => {
            // Skip-if-already-cached: a redundant click on the already-selected
            // Group with a warm cache is inert (mirrors `handle_species_selected`'s
            // early return). When the cache is empty — e.g. after a settings load
            // — the guard does not apply and the click triggers the fetch.
            let already_cached = app.settings.organism_group.as_deref() == Some(group.as_str())
                && app.cache.modules_pack.is_some()
                && app.cache.group_org_codes.is_some();
            if !already_cached {
                spawn_modules_fetch(app, level, group, org_codes, false);
            }
        }
    }
}

/// Render the inline progress strip for whichever fetch the ACTIVE mode
/// has in flight. When the active mode's `<mode>_fetch == None`, nothing
/// renders.
fn render_inline_fetch_progress(ui: &mut egui::Ui, app: &mut App) {
    let AppState::Stage3EnrichSetup {
        kegg_fetch,
        modules_fetch,
        ..
    } = &app.state
    else {
        return;
    };
    match app.settings.analysis_mode {
        AnalysisMode::Pathway => {
            let Some(f) = kegg_fetch.as_ref() else {
                return;
            };
            let fraction = if f.total == 0 {
                0.0
            } else {
                f.completed as f32 / f.total as f32
            };
            ui.horizontal(|ui| {
                crate::ui::widgets::progress_bar(
                    ui,
                    egui::ProgressBar::new(fraction)
                        .text(format!("{} / {}", f.completed, f.total))
                        .desired_width(220.0),
                    crate::theme::SURFACE,
                );
                if !f.current_pathway.is_empty() {
                    ui.label(
                        RichText::new(format!("Fetching {}", f.current_pathway))
                            .small()
                            .color(theme::TEXT),
                    );
                }
            });
        }
        AnalysisMode::Module => {
            let Some(f) = modules_fetch.as_ref() else {
                return;
            };
            let fraction = if f.total == 0 {
                0.0
            } else {
                f.completed as f32 / f.total as f32
            };
            ui.horizontal(|ui| {
                crate::ui::widgets::progress_bar(
                    ui,
                    egui::ProgressBar::new(fraction)
                        .text(format!("{} / {}", f.completed, f.total))
                        .desired_width(220.0),
                    crate::theme::SURFACE,
                );
                let mut hint = if f.current_id.is_empty() {
                    String::new()
                } else {
                    format!("Fetching {}", f.current_id)
                };
                if let Some(eta) = f.eta_secs {
                    if !hint.is_empty() {
                        hint.push(' ');
                    }
                    hint.push_str(&format!("ETA {eta}s"));
                }
                if !hint.is_empty() {
                    ui.label(RichText::new(hint).small().color(theme::TEXT));
                }
            });
        }
    }
}

/// True iff the user can click `Run enrichment` right now: active mode's
/// selection is complete AND its cache is populated AND no fetch for that
/// mode is in flight.
fn run_button_enabled(app: &App) -> bool {
    let AppState::Stage3EnrichSetup {
        kegg_fetch,
        modules_fetch,
        ..
    } = &app.state
    else {
        return false;
    };
    match app.settings.analysis_mode {
        AnalysisMode::Pathway => {
            kegg_fetch.is_none()
                && app.settings.kegg_species.is_some()
                && app.cache.species_kegg.is_some()
        }
        AnalysisMode::Module => {
            modules_fetch.is_none()
                && app.settings.organism_group_level.is_some()
                && app.settings.organism_group.is_some()
                && app.cache.modules_pack.is_some()
                && app.cache.group_org_codes.is_some()
        }
    }
}

/// Look up a Group's organism-code set in the on-disk-derived index with NO
/// network call. Bounds-guarded: a `level` outside `1..=3` or an absent group
/// returns `None`. Used by the module rehydrate path to repopulate
/// `cache.group_org_codes` from a selection restored by a settings snapshot.
fn lookup_group_codes(
    index: &crate::kegg::OrganismGroupIndex,
    level: u8,
    group: &str,
) -> Option<std::collections::HashSet<String>> {
    if !(1..=3).contains(&level) {
        return None;
    }
    index.by_level[(level - 1) as usize].get(group).cloned()
}

/// What [`rehydrate_stage3_cache`] should do for a restored selection whose
/// cache is (partly) empty. Pure decision — no side effects — so it is
/// unit-testable without egui or network, mirroring the pure-core /
/// thin-wrapper split used by `build_stage3_run_inputs`.
#[derive(Debug, Clone, PartialEq, Eq)]
enum RehydrateAction {
    None,
    LoadPathway(String),
    LoadModule {
        level: u8,
        group: String,
        need_codes: bool,
        need_pack: bool,
    },
}

/// Decide how to reconcile a restored KEGG selection with the live cache on
/// entering Stage 3 setup. Returns a non-`None` action ONLY when the active
/// mode carries a selection whose cache is missing and no fetch for that mode is
/// in flight — so it is a no-op on the normal forward path (selection `None`)
/// and on back/forward navigation (cache already present).
fn rehydrate_action(
    settings: &SessionSettings,
    cache: &SessionCache,
    kegg_fetch_in_flight: bool,
    modules_fetch_in_flight: bool,
) -> RehydrateAction {
    match settings.analysis_mode {
        AnalysisMode::Pathway => match &settings.kegg_species {
            Some(code) if cache.species_kegg.is_none() && !kegg_fetch_in_flight => {
                RehydrateAction::LoadPathway(code.clone())
            }
            _ => RehydrateAction::None,
        },
        AnalysisMode::Module => match (settings.organism_group_level, &settings.organism_group) {
            (Some(level), Some(group)) => {
                let need_codes = cache.group_org_codes.is_none();
                let need_pack = cache.modules_pack.is_none() && !modules_fetch_in_flight;
                if need_codes || need_pack {
                    RehydrateAction::LoadModule {
                        level,
                        group: group.clone(),
                        need_codes,
                        need_pack,
                    }
                } else {
                    RehydrateAction::None
                }
            }
            _ => RehydrateAction::None,
        },
    }
}

/// Reconcile a restored KEGG selection with the live cache on entering Stage 3
/// setup. Called from the two transitions that construct `Stage3EnrichSetup`
/// (`stage2_threshold::continue_to_enrichment` and the `stepper` jump to step 3)
/// so that a snapshot-restored species/Group — whose fetched cache is NOT
/// restored by `apply_snapshot` — auto-loads and `Run Enrichment` enables
/// without manual re-selection. Transition-triggered (NOT per-frame) so a failed
/// fetch does not re-spawn in a loop; the in-flight slot set by the spawn helpers
/// blocks a double-spawn within one transition. A no-op unless the active mode
/// has a selection whose cache is empty (see [`rehydrate_action`]).
pub(crate) fn rehydrate_stage3_cache(app: &mut App) {
    let (kegg_in_flight, modules_in_flight) = match &app.state {
        AppState::Stage3EnrichSetup {
            kegg_fetch,
            modules_fetch,
            ..
        } => (kegg_fetch.is_some(), modules_fetch.is_some()),
        _ => return,
    };

    match rehydrate_action(&app.settings, &app.cache, kegg_in_flight, modules_in_flight) {
        RehydrateAction::None => {}
        // Reuses the species selector's load path: disk fast-path via
        // `kegg::cache::read_species`, else spawns the network fetch.
        RehydrateAction::LoadPathway(code) => handle_species_selected(app, code),
        RehydrateAction::LoadModule {
            level,
            group,
            need_codes,
            need_pack,
        } => {
            // Re-derive the Group's org-code set from the on-disk index — no
            // network. Skip silently if the index is unreadable or the stored
            // (level, group) is no longer present (the manual selector path
            // still works).
            let codes = match crate::kegg::cache::read_organism_group_index() {
                Ok(Some(index)) => lookup_group_codes(&index, level, &group),
                Ok(None) => None,
                Err(e) => {
                    warn!(error = %e, "failed to read organism group index for Stage 3 rehydrate");
                    None
                }
            };
            let Some(codes) = codes else {
                return;
            };
            if need_pack {
                // `spawn_modules_fetch` writes `group_org_codes` AND (on Done)
                // `modules_pack`, so the `need_codes` write is covered here too.
                spawn_modules_fetch(app, level, group, codes, false);
            } else if need_codes {
                // Defensive partial-state branch (pack present, codes missing) —
                // unreachable in the load flow where both are None together.
                app.cache.group_org_codes = Some(codes);
            }
        }
    }
}

fn handle_species_selected(app: &mut App, code: String) {
    let unchanged = app.settings.kegg_species.as_deref() == Some(code.as_str());
    if unchanged && app.cache.species_kegg.is_some() {
        return;
    }
    app.settings.kegg_species = Some(code.clone());
    app.cache.species_kegg = None;

    // Cache fast-path.
    match kegg::cache::read_species(&code) {
        Ok(Some(species_data)) => {
            info!(code = %code, pathways = species_data.pathways.len(), "loaded species from cache fast-path");
            app.cache.species_kegg = Some(species_data);
            return;
        }
        Ok(None) => {}
        Err(e) => {
            error!(code = %code, error = %e, "failed to read species cache; falling back to network");
        }
    }

    spawn_species_fetch(app, code);
}

// `pub(crate)` so the bottom-panel Data tab's relocated `Refresh KEGG pathway
// cache` button (Cache-data block) can trigger a force re-fetch directly.
pub(crate) fn handle_species_refresh(app: &mut App) {
    let code = match app.settings.kegg_species.as_ref() {
        Some(c) => c.clone(),
        None => return,
    };
    if let AppState::Stage3EnrichSetup { kegg_fetch, .. } = &app.state
        && kegg_fetch.is_some()
    {
        warn!("KEGG species fetch already in progress; ignoring re-fetch click");
        return;
    }
    if let Err(e) = kegg::invalidate_cache(KeggCacheScope::Species(code.clone())) {
        error!(code = %code, error = %e, "failed to invalidate species cache");
    }
    app.cache.species_kegg = None;
    spawn_species_fetch(app, code);
}

fn spawn_species_fetch(app: &mut App, code: String) {
    if app.inputs.ion_tables.is_empty() || app.inputs.mapping.is_none() {
        warn!(code = %code, "cannot start KEGG species fetch without table + mapping; aborting");
        return;
    }

    let (event_tx, event_rx) = mpsc::channel::<KeggEvent>();
    let (progress_tx, mut progress_rx) = tokio_mpsc::channel::<KeggProgress>(64);
    let event_tx_for_progress = event_tx.clone();
    let client = app.kegg.clone();
    let code_for_task = code.clone();

    app.rt.spawn(async move {
        while let Some(p) = progress_rx.recv().await {
            if event_tx_for_progress.send(KeggEvent::Progress(p)).is_err() {
                break;
            }
        }
    });

    app.rt.spawn(async move {
        info!(code = %code_for_task, "starting KEGG species fetch");
        let outcome = kegg::fetch_species_pathways(&client, &code_for_task, progress_tx).await;
        let event = match outcome {
            Ok(s) => KeggEvent::Done(s),
            Err(e) => {
                error!(code = %code_for_task, error = %e, "KEGG species fetch failed");
                KeggEvent::Failed(e.to_string())
            }
        };
        let _ = event_tx.send(event);
    });

    if let AppState::Stage3EnrichSetup { kegg_fetch, .. } = &mut app.state {
        *kegg_fetch = Some(KeggFetchInFlight {
            progress_rx: event_rx,
            completed: 0,
            total: 0,
            current_pathway: String::new(),
        });
    }
}

// `pub(crate)` so the bottom-panel Data tab's relocated `Refresh KEGG module
// cache` button (Cache-data block) can trigger a force re-fetch directly.
pub(crate) fn spawn_modules_fetch(
    app: &mut App,
    level: u8,
    group: String,
    org_codes: std::collections::HashSet<String>,
    force_refresh: bool,
) {
    if app.inputs.ion_tables.is_empty() || app.inputs.mapping.is_none() {
        warn!("spawn_modules_fetch called without complete inputs; aborting");
        return;
    }

    // Record the user's selection on settings/cache before the fetch
    // completes — Group selection persists across mode toggles per D3.
    app.settings.organism_group_level = Some(level);
    app.settings.organism_group = Some(group);
    app.cache.group_org_codes = Some(org_codes);

    let (event_tx, event_rx) = mpsc::channel::<crate::app::ModulesFetchEvent>();
    let (progress_tx, mut progress_rx) =
        tokio_mpsc::channel::<crate::kegg::ModuleFetchProgress>(64);
    let event_tx_for_progress = event_tx.clone();
    let kegg_client = app.kegg.clone();

    app.rt.spawn(async move {
        while let Some(p) = progress_rx.recv().await {
            if event_tx_for_progress
                .send(crate::app::ModulesFetchEvent::Progress(p))
                .is_err()
            {
                break;
            }
        }
    });

    app.rt.spawn(async move {
        let event = match crate::kegg::fetch_modules(&kegg_client, force_refresh, progress_tx).await
        {
            Ok(cache) => crate::app::ModulesFetchEvent::Done(cache),
            Err(e) => {
                error!(error = %e, "fetch_modules failed");
                crate::app::ModulesFetchEvent::Failed(e.to_string())
            }
        };
        let _ = event_tx.send(event);
    });

    if let AppState::Stage3EnrichSetup { modules_fetch, .. } = &mut app.state {
        *modules_fetch = Some(ModulesFetchInFlight {
            progress_rx: event_rx,
            completed: 0,
            total: 0,
            current_id: String::new(),
            eta_secs: None,
        });
    }
}

/// Build an `AnalysisPayload` from the current settings + cache. Returns
/// `None` if any required piece is missing (cache absent for the active
/// mode, or Module-mode group selection not made yet). Callers should
/// surface a user-facing error if they get `None` at a spawn site.
pub(crate) fn build_analysis_payload(
    settings: &SessionSettings,
    cache: &SessionCache,
) -> Option<AnalysisPayload> {
    match settings.analysis_mode {
        AnalysisMode::Pathway => cache
            .species_kegg
            .clone()
            .map(|species_kegg| AnalysisPayload::Pathway { species_kegg }),
        AnalysisMode::Module => {
            let modules_pack = cache.modules_pack.clone()?;
            let group_level = settings.organism_group_level?;
            let group_name = settings.organism_group.clone()?;
            let group_org_codes = cache.group_org_codes.clone()?;
            Some(AnalysisPayload::Module {
                modules_pack,
                group_level,
                group_name,
                group_org_codes,
                min_group_overlap: settings.min_group_overlap,
            })
        }
    }
}

/// Build the orchestrator inputs `(Stage3Params, target, pubchem_total)` from a
/// `dam_results` slice plus the current `settings`/`cache`. The SINGLE place that
/// constructs `Stage3Params` from `app.settings`, resolves the `AnalysisTarget`
/// via `build_analysis_payload`, and sums `pubchem_total` (InChIKey-bearing
/// features ACROSS all modes — additive, not deduplicated; design D3). Returns
/// `None` when `build_analysis_payload` is `None` (KEGG cache for the active mode
/// missing). `method` reads `dam_results[0].method` (all modes share method per
/// Stage 2 design). Both `force_refresh_*` default `false`; `start_refresh`
/// overrides them after the call. Shared by `build_stage3_spawn_inputs` (Run),
/// `rerun`, and `start_refresh`.
pub(crate) fn build_stage3_run_inputs(
    dam_results: &[DamResult],
    settings: &SessionSettings,
    cache: &SessionCache,
) -> Option<(Stage3Params, AnalysisPayload, usize)> {
    let target = build_analysis_payload(settings, cache)?;
    let params = Stage3Params {
        method: dam_results[0].method,
        fc_threshold: settings.fc_threshold,
        fdr_threshold: settings.fdr_threshold,
        delta_threshold: settings.delta_threshold,
        direction: settings.direction,
        min_hit_count: settings.min_hit_count,
        min_entry_size: settings.min_entry_size,
        fdr_method: settings.enrichment_fdr_method,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };
    let pubchem_total: usize = dam_results
        .iter()
        .map(|d| d.features.iter().filter(|f| f.inchikey.is_some()).count())
        .sum();
    Some((params, target, pubchem_total))
}

/// Setup-screen wrapper over [`build_stage3_run_inputs`]: destructures the
/// `AppState::Stage3EnrichSetup` variant and prepends the cloned `dam_results`.
/// Returns `None` if the state variant is wrong OR the KEGG cache for the active
/// mode is missing. Kept as `pub fn` with this exact signature so integration
/// tests (`tests/stage3_ui_spawn_test.rs`) can drive the UI plumbing without an
/// eframe app. The return tuple is `(dam_results_clone, params, target, pubchem_total)`;
/// `dam_results_clone` is NEVER collapsed to length 1 from a dual-mode source —
/// that bug is the reason this seam exists (`fix-stage3-ui-dual-mode-spawn`).
pub fn build_stage3_spawn_inputs(
    state: &AppState,
    settings: &SessionSettings,
    cache: &SessionCache,
) -> Option<(Vec<DamResult>, Stage3Params, AnalysisPayload, usize)> {
    let AppState::Stage3EnrichSetup { dam_results, .. } = state else {
        return None;
    };
    let (params, target, pubchem_total) = build_stage3_run_inputs(dam_results, settings, cache)?;
    Some((dam_results.clone(), params, target, pubchem_total))
}

fn start_run(app: &mut App) {
    // Compute spawn inputs from refs first (non-destructive). If the cache is
    // missing the helper returns None and we restore Stage3EnrichSetup with an
    // error message.
    let spawn_inputs = build_stage3_spawn_inputs(&app.state, &app.settings, &app.cache);

    let prev = std::mem::take(&mut app.state);
    let AppState::Stage3EnrichSetup { dam_results, .. } = prev else {
        return;
    };

    let Some((dam_results_clone, params, target_clone, pubchem_total)) = spawn_inputs else {
        app.state = AppState::Stage3EnrichSetup {
            dam_results,
            error: Some(
                "KEGG cache for the selected mode is missing; pick a species/Group again to re-fetch."
                    .to_string(),
            ),
            kegg_fetch: None,
            modules_fetch: None,
        };
        return;
    };

    // `start_run` emits the `n_modes` diagnostic BEFORE spawning (spec
    // requirement: bug-report bundles record `n_modes` even if the orchestrator
    // fails to start). `spawn_stage3_run` is log-silent, so `rerun` — which has
    // never emitted this line — stays byte-identical. (`dam_results` from `prev`
    // is consumed only by the cache-missing arm above; the running state's Vec
    // is `dam_results_clone`.)
    info!(
        direction = ?app.settings.direction,
        top_n = app.settings.top_n,
        enrichment_fdr_threshold = app.settings.enrichment_fdr_threshold,
        min_hit_count = app.settings.min_hit_count,
        n_modes = dam_results_clone.len(),
        pubchem_inputs = pubchem_total,
        "Stage 3 Run starting"
    );
    app.spawn_stage3_run(dam_results_clone, target_clone, params, pubchem_total);
}

#[cfg(test)]
mod build_run_inputs_tests {
    use super::*;
    use crate::dam::DamMethod;
    use crate::dam::fdr::FdrMethod;
    use crate::dam::types::{DamFeature, FcBasis};
    use crate::kegg::{KeggCompoundSet, SpeciesKegg};

    fn feat(inchikey: Option<&str>) -> DamFeature {
        DamFeature {
            alignment_id: "aid".into(),
            metabolite_name: "met".into(),
            inchikey: inchikey.map(|s| s.to_string()),
            average_rt_min: None,
            average_mz: None,
            formula: None,
            smiles: None,
            numerator_mean: 0.0,
            denominator_mean: 0.0,
            numerator_median: 0.0,
            denominator_median: 0.0,
            fold_change: 0.0,
            log2_fold_change: 0.0,
            fc_basis: FcBasis::Mean,
            p_value: 1.0,
            p_adjusted: 1.0,
            neg_log10_p_adjusted: 0.0,
            effect_size: None,
        }
    }

    fn dam(features: Vec<DamFeature>) -> DamResult {
        DamResult {
            method: DamMethod::Welch,
            numerator: "A".into(),
            denominator: "B".into(),
            features,
            skipped: 0,
            fdr_method: FdrMethod::BenjaminiHochberg,
            dedup_report: None,
        }
    }

    fn species_cache() -> SessionCache {
        SessionCache {
            species_kegg: Some(SpeciesKegg {
                code: "tst".into(),
                fetched_at: chrono::Utc::now(),
                pathways: vec![KeggCompoundSet {
                    id: "tst00001".into(),
                    name: "p1".into(),
                    compounds: vec!["C1".into()],
                }],
            }),
            ..SessionCache::default()
        }
    }

    #[test]
    fn maps_settings_and_sums_pubchem_total_across_modes() {
        // POS: 2 annotated + 1 unknown; NEG: 1 annotated → pubchem_total = 3.
        let dam_results = vec![
            dam(vec![feat(Some("K1")), feat(Some("K2")), feat(None)]),
            dam(vec![feat(Some("K3"))]),
        ];
        let settings = SessionSettings {
            min_entry_size: 7,
            min_hit_count: 4,
            ..SessionSettings::default()
        };
        let cache = species_cache();

        let (params, target, pubchem_total) =
            build_stage3_run_inputs(&dam_results, &settings, &cache).expect("cache present");

        // pubchem_total is the additive InChIKey-bearing count across BOTH modes.
        assert_eq!(pubchem_total, 3);
        // params maps from settings; `method` reads `dam_results[0]`.
        assert_eq!(params.method, DamMethod::Welch);
        assert_eq!(params.min_entry_size, 7);
        assert_eq!(params.min_hit_count, 4);
        assert_eq!(params.direction, settings.direction);
        assert_eq!(params.fdr_method, settings.enrichment_fdr_method);
        // Float fields are exact copies — compare bit patterns (clippy float_cmp).
        assert_eq!(
            params.fc_threshold.to_bits(),
            settings.fc_threshold.to_bits()
        );
        assert_eq!(
            params.fdr_threshold.to_bits(),
            settings.fdr_threshold.to_bits()
        );
        assert_eq!(
            params.delta_threshold.to_bits(),
            settings.delta_threshold.to_bits()
        );
        // Both force_refresh flags default false (start_refresh overrides after).
        assert!(!params.force_refresh_pubchem);
        assert!(!params.force_refresh_kegg_conv);
        assert!(matches!(target, AnalysisPayload::Pathway { .. }));
    }

    #[test]
    fn returns_none_when_species_kegg_absent() {
        let dam_results = vec![dam(vec![feat(Some("K1"))])];
        let settings = SessionSettings::default();
        let cache = SessionCache::default(); // species_kegg = None
        assert!(build_stage3_run_inputs(&dam_results, &settings, &cache).is_none());
    }
}

#[cfg(test)]
mod rehydrate_tests {
    use super::*;
    use crate::kegg::{KeggModulesCache, OrganismGroupIndex, SpeciesKegg};
    use std::collections::{HashMap, HashSet};

    fn index_with_animals() -> OrganismGroupIndex {
        let mut by_level: [HashMap<String, HashSet<String>>; 3] =
            [HashMap::new(), HashMap::new(), HashMap::new()];
        by_level[1].insert(
            "Animals".to_string(),
            HashSet::from(["hsa".to_string(), "mmu".to_string()]),
        );
        OrganismGroupIndex {
            fetched_at: chrono::Utc::now(),
            by_level,
        }
    }

    fn some_species_kegg() -> SpeciesKegg {
        SpeciesKegg {
            code: "hsa".into(),
            fetched_at: chrono::Utc::now(),
            pathways: vec![],
        }
    }

    fn some_modules_pack() -> KeggModulesCache {
        KeggModulesCache {
            modules: HashMap::new(),
        }
    }

    fn pathway_settings() -> SessionSettings {
        SessionSettings {
            analysis_mode: AnalysisMode::Pathway,
            kegg_species: Some("hsa".into()),
            ..SessionSettings::default()
        }
    }

    fn module_settings() -> SessionSettings {
        SessionSettings {
            analysis_mode: AnalysisMode::Module,
            organism_group_level: Some(2),
            organism_group: Some("Animals".into()),
            ..SessionSettings::default()
        }
    }

    // ── lookup_group_codes ──

    #[test]
    fn lookup_group_codes_hit_returns_the_set() {
        let index = index_with_animals();
        let codes = lookup_group_codes(&index, 2, "Animals").expect("Animals present at L2");
        assert_eq!(codes, HashSet::from(["hsa".to_string(), "mmu".to_string()]));
    }

    #[test]
    fn lookup_group_codes_absent_group_is_none() {
        let index = index_with_animals();
        assert!(lookup_group_codes(&index, 2, "Plants").is_none());
        // Right name but wrong level (Animals lives at L2, not L1/L3).
        assert!(lookup_group_codes(&index, 1, "Animals").is_none());
    }

    #[test]
    fn lookup_group_codes_out_of_range_level_is_none() {
        let index = index_with_animals();
        assert!(lookup_group_codes(&index, 0, "Animals").is_none());
        assert!(lookup_group_codes(&index, 4, "Animals").is_none());
    }

    // ── rehydrate_action: Pathway ──

    #[test]
    fn pathway_no_selection_is_noop() {
        // Default settings: Pathway mode, kegg_species = None (normal forward path).
        let settings = SessionSettings::default();
        let cache = SessionCache::default();
        assert_eq!(
            rehydrate_action(&settings, &cache, false, false),
            RehydrateAction::None
        );
    }

    #[test]
    fn pathway_warm_cache_is_noop() {
        // Selection restored AND species cache present (back/forward nav).
        let settings = pathway_settings();
        let cache = SessionCache {
            species_kegg: Some(some_species_kegg()),
            ..SessionCache::default()
        };
        assert_eq!(
            rehydrate_action(&settings, &cache, false, false),
            RehydrateAction::None
        );
    }

    #[test]
    fn pathway_selection_with_empty_cache_loads() {
        // The load-settings situation: species selected, cache empty, no fetch.
        let settings = pathway_settings();
        let cache = SessionCache::default();
        assert_eq!(
            rehydrate_action(&settings, &cache, false, false),
            RehydrateAction::LoadPathway("hsa".into())
        );
    }

    #[test]
    fn pathway_fetch_in_flight_is_noop() {
        let settings = pathway_settings();
        let cache = SessionCache::default();
        assert_eq!(
            rehydrate_action(&settings, &cache, true, false),
            RehydrateAction::None
        );
    }

    // ── rehydrate_action: Module — all four (need_codes, need_pack) combos ──

    #[test]
    fn module_no_selection_is_noop() {
        let settings = SessionSettings {
            analysis_mode: AnalysisMode::Module,
            ..SessionSettings::default()
        };
        let cache = SessionCache::default();
        assert_eq!(
            rehydrate_action(&settings, &cache, false, false),
            RehydrateAction::None
        );
    }

    #[test]
    fn module_both_empty_needs_codes_and_pack() {
        // (need_codes=true, need_pack=true): the load-settings situation.
        let settings = module_settings();
        let cache = SessionCache::default();
        assert_eq!(
            rehydrate_action(&settings, &cache, false, false),
            RehydrateAction::LoadModule {
                level: 2,
                group: "Animals".into(),
                need_codes: true,
                need_pack: true,
            }
        );
    }

    #[test]
    fn module_codes_present_pack_missing_needs_pack_only() {
        // (need_codes=false, need_pack=true).
        let settings = module_settings();
        let cache = SessionCache {
            group_org_codes: Some(HashSet::from(["hsa".to_string()])),
            ..SessionCache::default()
        };
        assert_eq!(
            rehydrate_action(&settings, &cache, false, false),
            RehydrateAction::LoadModule {
                level: 2,
                group: "Animals".into(),
                need_codes: false,
                need_pack: true,
            }
        );
    }

    #[test]
    fn module_pack_present_codes_missing_needs_codes_only() {
        // (need_codes=true, need_pack=false): the defensive partial-state arm —
        // unreachable in the load flow (both are None together) but covered here.
        let settings = module_settings();
        let cache = SessionCache {
            modules_pack: Some(some_modules_pack()),
            ..SessionCache::default()
        };
        assert_eq!(
            rehydrate_action(&settings, &cache, false, false),
            RehydrateAction::LoadModule {
                level: 2,
                group: "Animals".into(),
                need_codes: true,
                need_pack: false,
            }
        );
    }

    #[test]
    fn module_fetch_in_flight_clears_need_pack() {
        // (need_codes=true, need_pack=false) via an in-flight modules fetch.
        let settings = module_settings();
        let cache = SessionCache::default();
        assert_eq!(
            rehydrate_action(&settings, &cache, false, true),
            RehydrateAction::LoadModule {
                level: 2,
                group: "Animals".into(),
                need_codes: true,
                need_pack: false,
            }
        );
    }

    #[test]
    fn module_warm_cache_is_noop() {
        // (need_codes=false, need_pack=false): both caches present → None.
        let settings = module_settings();
        let cache = SessionCache {
            modules_pack: Some(some_modules_pack()),
            group_org_codes: Some(HashSet::from(["hsa".to_string()])),
            ..SessionCache::default()
        };
        assert_eq!(
            rehydrate_action(&settings, &cache, false, false),
            RehydrateAction::None
        );
    }
}
