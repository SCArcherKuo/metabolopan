//! Stage 3 setup screen. Controls: Analysis Mode toggle, species/group
//! selector, inline KEGG fetch progress, direction selector, minimum entry
//! size, FDR correction radio, Run Enrichment button. (Top N, the significance
//! threshold and minimum hit count live on the Stage 3 RESULT screen — they
//! were relocated by `add-bottom-panel-data-tab`.)
//!
//! Per `reorder-gui-and-move-mode-to-stage3` the Mode toggle and the
//! mode-specific selector both live on this screen (no longer Stage 1),
//! and a fetch triggered by the user picking a species / Group writes
//! into `Stage3EnrichSetup.kegg_fetch` / `modules_fetch` rather than
//! transitioning to a dedicated fetching variant. Toggling mode while a
//! fetch is IN FLIGHT cancels + clears that fetch (and blanks its
//! un-fetched selection) via `cancel_inflight_for_mode_switch`, so the two
//! modes never fetch at once and contend for the shared KEGG rate limit —
//! superseding the older "both slots coexist in flight (D2 + D6)" stance;
//! completed caches/selections still coexist across the toggle.

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

        // === Shared analysis-target block ===
        render_analysis_target(ui, app);

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
            // Method-dependent: under `NoCorrection` there is no FDR stage for
            // the drop to precede and no multiple-testing penalty to reduce, so
            // both sentences would be false. Read from SETTINGS here — this
            // screen describes the run the user is about to start, the
            // complement of the result screen's method-from-the-run rule.
            .on_hover_text(match settings.enrichment_fdr_method {
                crate::dam::fdr::FdrMethod::NoCorrection => {
                    "Drop pathways/modules with fewer than N measurable compounds \
                     before the significance test."
                }
                _ => {
                    "Drop pathways/modules with fewer than N measurable compounds before FDR. \
                     Reduces multiple-testing penalty."
                }
            });
        });
        ui.add_space(8.0);

        // `Enrichment FDR threshold` + `Minimum hit count` moved to the Stage 3
        // result screen (`add-bottom-panel-data-tab`) so they can be tuned and
        // re-applied without navigating back. The FDR correction *method* stays
        // here: it decides which significance quantity the run produces, not
        // just where the display cut falls.

        // FDR correction radio. Independent of Stage 2's choice (per add-bh-fdr D3).
        // `No correction` is exposed here ONLY (Stage 2 setup hides it; an
        // adversarial snapshot is defensively coerced back to BH in
        // `apply_snapshot`).
        // Stage 3 exposes the `NoCorrection` variant → `include_none = true`.
        crate::ui::widgets::fdr_method_radios(ui, &mut settings.enrichment_fdr_method, true);
        ui.label(
            RichText::new(
                "No correction skips multiple-testing correction entirely — exploratory only, not for publication.",
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
        let run_enabled = target_ready(app);
        let mut action = Action::None;
        // Disabled-state hint: explain *why* the button is disabled when a
        // fetch is in flight in the active mode (per stage3-ui spec scenario
        // "Inline progress strip during pathway fetch" and design D2).
        let disabled_hint = (!run_enabled).then(|| fetch_in_flight_hint(app)).flatten();
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

/// Why the active setup screen's Run button is disabled, when the reason is a
/// KEGG fetch in flight for the ACTIVE mode. `None` when no such fetch is
/// running (the caller has other reasons and its own texts for those).
///
/// Shared by `Run Enrichment` and `Run Coverage`: the two buttons gate on the
/// same fetch for the same reason, and the hint is the user's only explanation
/// of why a button they can see is inert. Reads the slots through
/// [`crate::app::setup_fetch_slots`], so it is variant-agnostic like the block
/// it belongs to.
pub(crate) fn fetch_in_flight_hint(app: &App) -> Option<&'static str> {
    let (kegg_fetch, modules_fetch) = crate::app::setup_fetch_slots(&app.state)?;
    match app.settings.analysis_mode {
        AnalysisMode::Pathway => kegg_fetch
            .is_some()
            .then_some("Waiting for KEGG pathway fetch…"),
        AnalysisMode::Module => modules_fetch
            .is_some()
            .then_some("Waiting for KEGG modules fetch…"),
    }
}

/// Render the Pathway / Module radio. On change, update
/// `settings.analysis_mode` via the named API. Pathway and Module
/// selections coexist — the API is a near-no-op that only sets the mode
/// (per `reorder-gui-and-move-mode-to-stage3` D3).
/// The shared **analysis-target block**: the Analysis Mode toggle, the
/// mode-aware target selector, and the inline KEGG fetch progress strip.
///
/// Rendered identically by every screen that selects an analysis target — the
/// Stage 3 enrichment setup screen and, on the KEGG coverage route, the
/// coverage setup screen. There MUST NOT be a second implementation: the two
/// screens differ only in the controls AROUND this block, and a divergence
/// between copies would surface as two subtly different species selectors
/// reached by two different routes (owner: the `stage3-ui` capability).
///
/// Variant-agnostic by construction — every slot access inside goes through
/// `crate::app::{is_target_setup, setup_fetch_slots, setup_fetch_slots_mut}`,
/// so a new setup variant is enabled by adding one arm there, not by editing
/// this function or any of the renderers it calls.
pub(crate) fn render_analysis_target(ui: &mut egui::Ui, app: &mut App) {
    render_mode_toggle(ui, app);
    ui.add_space(6.0);

    match app.settings.analysis_mode {
        AnalysisMode::Pathway => render_species_selector(ui, app),
        AnalysisMode::Module => {
            render_organism_group_selector(ui, app);
            // Module-mode-only Group-overlap threshold, directly under the
            // Level + Group picker (binds the existing min_group_overlap).
            render_min_group_overlap(ui, app);
        }
    }
    ui.add_space(6.0);

    // Inline KEGG fetch progress strip (active mode only).
    render_inline_fetch_progress(ui, app);
    ui.add_space(8.0);
}

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
        // Cancel the in-flight fetch on the mode being left (so it stops
        // hitting the shared KEGG client / rate limit) and blank an
        // incomplete selection. Completed caches + selections still persist.
        cancel_inflight_for_mode_switch(app, new_mode);
        // The body of this API was reduced to `self.analysis_mode = new_mode;`
        // by Phase 2 of `reorder-gui-and-move-mode-to-stage3`, but the call
        // site stays so spec scenarios remain anchored to a named API.
        app.settings.reset_kegg_selection_for_mode_switch(new_mode);
        app.cache.clear_for_mode_switch(new_mode);
    }
}

/// Cancel the in-flight fetch on the mode being left (so it stops contending
/// for the shared KEGG client) and, because an INCOMPLETE fetch leaves its
/// selection with no completed cache backing, reset that selection to blank.
/// A COMPLETED selection/cache (no fetch in flight) is left untouched, so the
/// `analysis-mode-toggle` coexistence contract still holds. `new_mode` is used
/// only for the diagnostic log.
pub(crate) fn cancel_inflight_for_mode_switch(app: &mut App, new_mode: AnalysisMode) {
    let (species_incomplete, group_incomplete) = crate::app::setup_fetch_slots(&app.state)
        .map(|(k, m)| (k.is_some(), m.is_some()))
        .unwrap_or((false, false));

    crate::app::abort_and_clear_setup_fetches(&mut app.state);

    if species_incomplete {
        info!(
            ?new_mode,
            "stopping in-flight KEGG species fetch: analysis mode was switched before it finished; resetting the incomplete species selection to blank"
        );
        app.settings.kegg_species = None;
        app.cache.species_kegg = None;
        app.species_selector = species_selector::SpeciesSelectorState::default();
    }
    if group_incomplete {
        info!(
            ?new_mode,
            "stopping in-flight KEGG module fetch: analysis mode was switched before it finished; resetting the incomplete group selection to blank"
        );
        app.settings.organism_group = None;
        app.settings.organism_group_level = None;
        app.cache.group_org_codes = None;
        app.organism_group_selector =
            crate::ui::organism_group_selector::OrganismGroupSelectorState::default();
    }
}

fn render_species_selector(ui: &mut egui::Ui, app: &mut App) {
    let (organisms_view, loading, load_error) = match &app.organisms.state {
        OrganismsLoadState::Loaded { organisms, .. } => (Some(organisms.as_slice()), false, None),
        OrganismsLoadState::Loading { .. } => (None, true, None),
        OrganismsLoadState::Failed(msg) => (None, false, Some(msg.as_str())),
        OrganismsLoadState::Idle => (None, false, None),
    };
    let current = app.settings.kegg_species.clone();
    let selector_enabled = crate::app::is_target_setup(&app.state);

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
    let enabled = crate::app::is_target_setup(&app.state);

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

/// Effective `DragValue` ceiling for `min_group_overlap`: the selected Group's
/// organism count (`group_len`) soft-capped at 20, but never below an
/// already-`stored` value so a hand-set value above the soft cap is preserved
/// rather than silently truncated (review M2). `group_len == usize::MAX` (no
/// Group selected/fetched) ⇒ only the soft cap applies. Callers then
/// `stored.clamp(1, eff_max)`, which only ever lowers a value to the Group's
/// actual organism count — the real domain limit.
fn min_group_overlap_eff_max(group_len: usize, stored: usize) -> usize {
    group_len.min(stored.max(20))
}

/// Module-mode-only control: the minimum number of the selected Group's
/// organisms a module must be complete in to enter the analysis. Binds the
/// existing `app.settings.min_group_overlap` (default 1; previously had no UI
/// and was pinned to its default). The drag ceiling is the selected Group's
/// organism count soft-capped at 20, but expands to preserve a value already
/// set above 20 (only reachable via a hand-edited session JSON) rather than
/// silently truncating it; the stored value is clamped to the Group's ACTUAL
/// organism count so an out-of-range leftover from a previously-larger Group
/// can't silently zero out module retention. Changing it re-filters the
/// already-fetched modules (no KEGG re-fetch) and the Data tab's
/// `In selected Group` count updates live (its memo keys on this value).
fn render_min_group_overlap(ui: &mut egui::Ui, app: &mut App) {
    // `usize::MAX` when no Group is selected/fetched ⇒ no group cap (soft cap 20
    // still applies via `.max(20)`). Read `app.cache` / write `app.settings`
    // directly here — this runs before `show`'s `let settings = &mut ...`
    // re-borrow, so reusing that alias would overlap borrows of `app`.
    let group_len = app
        .cache
        .group_org_codes
        .as_ref()
        .map(|c| c.len())
        .unwrap_or(usize::MAX);
    let eff_max = min_group_overlap_eff_max(group_len, app.settings.min_group_overlap);
    app.settings.min_group_overlap = app.settings.min_group_overlap.clamp(1, eff_max);
    ui.horizontal(|ui| {
        ui.label("Minimum number of group organisms a module must be complete in:");
        ui.add(
            egui::DragValue::new(&mut app.settings.min_group_overlap)
                .speed(1)
                .range(1..=eff_max),
        )
        .on_hover_text(
            "Keep only modules that are complete in at least N organisms of the selected group. \
             Higher values retain only modules conserved across the group.",
        );
    });
}

/// Render the inline progress strip for whichever fetch the ACTIVE mode
/// has in flight. When the active mode's `<mode>_fetch == None`, nothing
/// renders.
fn render_inline_fetch_progress(ui: &mut egui::Ui, app: &mut App) {
    let Some((kegg_fetch, modules_fetch)) = crate::app::setup_fetch_slots(&app.state) else {
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
/// Is the active mode's analysis target selected AND fetched, with no fetch in
/// flight?
///
/// Shared by `Run Enrichment` and `Run Coverage`: both gate on the same three
/// conditions for the same reason, and a second copy would let the two screens
/// disagree about when a target is usable.
pub(crate) fn target_ready(app: &App) -> bool {
    let Some((kegg_fetch, modules_fetch)) = crate::app::setup_fetch_slots(&app.state) else {
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
    // Through the shared accessor, so both setup screens rehydrate: a
    // coverage-route jump back to `Stage2CoverageSetup` restores a selection
    // whose cache is empty exactly as the enrichment setup screen does.
    let (kegg_in_flight, modules_in_flight) = match crate::app::setup_fetch_slots(&app.state) {
        Some((kegg_fetch, modules_fetch)) => (kegg_fetch.is_some(), modules_fetch.is_some()),
        None => return,
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

    // The selection changed: cancel any in-flight fetch for the PREVIOUS
    // species and clear its progress strip NOW, on every downstream path.
    // The cache fast-path below returns before `spawn_species_fetch` (whose
    // own abort-prior would otherwise never run), so a switch to an
    // already-cached species would leave the old fetch streaming and its
    // progress bar visible.
    if let Some((kegg_fetch, _)) = crate::app::setup_fetch_slots_mut(&mut app.state)
        && let Some(prev) = kegg_fetch.take()
    {
        info!(
            species = %code,
            "stopping the in-flight KEGG species fetch: a different species was selected"
        );
        prev.abort_tasks();
    }

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

/// Drain `LogPaneState.organisms_refresh_requested` and render the
/// organism-roster refresh confirm. `pub(crate)` for its ONE caller,
/// `App::drain_frame_dialogs` — the frame owns this, not a screen, because the
/// `Refresh KEGG organism list` button renders on five `AppState` variants and
/// a drain owned by two `show` functions left the request set on the other
/// three. Modelled on the `RefreshState` cache-refresh confirms: it is NOT an
/// App-level `*ModalState` and sits OUTSIDE the four-modal mutual-exclusion
/// family / `drain_modal_requests` (per the `app-shell` organism-roster refresh
/// requirement). On confirm it calls `App::handle_organisms_refresh`.
///
/// Rendered unconditionally, so an unanswered confirm follows the user across
/// navigation rather than vanishing and re-appearing; `App::start_new_round` is
/// the only thing that closes it for them.
pub(crate) fn drain_organisms_refresh_confirm(app: &mut App, ctx: &egui::Context) {
    if std::mem::take(&mut app.log_ui.organisms_refresh_requested)
        && matches!(app.organisms.state, OrganismsLoadState::Loaded { .. })
    {
        app.log_ui.organisms_refresh_confirm_open = true;
    }
    if !app.log_ui.organisms_refresh_confirm_open {
        return;
    }
    let mut want_refresh = false;
    let mut want_cancel = false;
    egui::Window::new("Refresh KEGG organism list")
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, egui::vec2(0.0, 0.0))
        .show(ctx, |ui| {
            ui.label(
                "Re-fetch the full KEGG organism roster from KEGG and rebuild the organism groups?",
            );
            // Disclosure at the point of decision: the user cannot otherwise
            // see that confirming may empty their current selection. It stops
            // there deliberately — a completed run's recorded target is taken
            // from the run itself, so this cannot relabel finished results.
            ui.label(
                "If KEGG has retired your selected species or Group, it is cleared from the \
                 current selection and must be re-picked before the next run.",
            );
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Cancel").clicked() {
                    want_cancel = true;
                }
                if ui.button("Refresh").clicked() {
                    want_refresh = true;
                }
            });
        });
    if want_cancel {
        app.log_ui.organisms_refresh_confirm_open = false;
    }
    if want_refresh {
        app.log_ui.organisms_refresh_confirm_open = false;
        app.handle_organisms_refresh();
    }
}

// `pub(crate)` so the bottom-panel Data tab's relocated `Refresh KEGG pathway
// cache` button (Cache-data block) can trigger a force re-fetch directly.
pub(crate) fn handle_species_refresh(app: &mut App) {
    let code = match app.settings.kegg_species.as_ref() {
        Some(c) => c.clone(),
        None => return,
    };
    // BEFORE `invalidate_cache`, not only inside `spawn_species_fetch`: this
    // function destroys the on-disk cache and clears `cache.species_kegg`
    // before delegating, so a precondition checked only at the spawn would
    // leave the user with neither cache and no fetch. Note the guard below
    // cannot cover this — `setup_fetch_slots` returns `None` off a setup
    // screen, so the `&&` short-circuits and never fires.
    if !crate::app::is_target_setup(&app.state) {
        warn!(
            "species cache refresh requested from a screen that owns no in-flight slot; refusing"
        );
        return;
    }
    if let Some((kegg_fetch, _)) = crate::app::setup_fetch_slots(&app.state)
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
    // Same precondition and the same reason as `spawn_modules_fetch`: a screen
    // must own the slot that will receive the result. Owner: the `kegg-fetching` capability.
    if !crate::app::is_target_setup(&app.state) {
        warn!(code = %code, "species fetch requested from a screen that owns no in-flight slot; refusing");
        if let Some(err) = crate::app::screen_error_mut(&mut app.state) {
            *err = Some("Cannot fetch KEGG pathways from this screen.".into());
        }
        return;
    }

    let (event_tx, event_rx) = mpsc::channel::<KeggEvent>();
    let (progress_tx, mut progress_rx) = tokio_mpsc::channel::<KeggProgress>(64);
    let event_tx_for_progress = event_tx.clone();
    let client = app.kegg.clone();
    let code_for_task = code.clone();

    let relay_handle = app
        .rt
        .spawn(async move {
            while let Some(p) = progress_rx.recv().await {
                if event_tx_for_progress.send(KeggEvent::Progress(p)).is_err() {
                    break;
                }
            }
        })
        .abort_handle();

    let fetch_handle = app
        .rt
        .spawn(async move {
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
        })
        .abort_handle();

    if let Some((kegg_fetch, _)) = crate::app::setup_fetch_slots_mut(&mut app.state) {
        // Cancel a prior in-flight fetch before replacing the slot so the
        // abandoned fetch stops issuing KEGG requests on the shared client
        // instead of running to completion and discarding its result.
        if let Some(prev) = kegg_fetch.take() {
            prev.abort_tasks();
        }
        *kegg_fetch = Some(KeggFetchInFlight {
            progress_rx: event_rx,
            completed: 0,
            total: 0,
            current_pathway: String::new(),
            fetch_handle,
            relay_handle,
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
    // The load-bearing precondition is that a screen can RECEIVE the result,
    // not what the user loaded: a KEGG catalogue fetch is keyed by Group and
    // reads nothing from the samples, and the coverage route's metadata `.csv`
    // is optional by design. Checked BEFORE spawning — the tasks below are
    // installed into a slot afterwards, so on a slotless state they would run
    // detached: undrained, their `AbortHandle`s never stored (so the in-flight
    // cancellation path can never reach them), holding `.modules.lock` for the
    // full 6-12 minute fetch. Owner: the `module-fetching` capability.
    if !crate::app::is_target_setup(&app.state) {
        warn!("module fetch requested from a screen that owns no in-flight slot; refusing");
        if let Some(err) = crate::app::screen_error_mut(&mut app.state) {
            *err = Some("Cannot refresh the KEGG module cache from this screen.".into());
        }
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

    let relay_handle = app
        .rt
        .spawn(async move {
            while let Some(p) = progress_rx.recv().await {
                if event_tx_for_progress
                    .send(crate::app::ModulesFetchEvent::Progress(p))
                    .is_err()
                {
                    break;
                }
            }
        })
        .abort_handle();

    let fetch_handle = app
        .rt
        .spawn(async move {
            let event =
                match crate::kegg::fetch_modules(&kegg_client, force_refresh, progress_tx).await {
                    Ok(cache) => crate::app::ModulesFetchEvent::Done(cache),
                    Err(e) => {
                        error!(error = %e, "fetch_modules failed");
                        crate::app::ModulesFetchEvent::Failed(e.to_string())
                    }
                };
            let _ = event_tx.send(event);
        })
        .abort_handle();

    if let Some((_, modules_fetch)) = crate::app::setup_fetch_slots_mut(&mut app.state) {
        // Cancel a prior in-flight fetch before replacing the slot. Aborting
        // the fetch task drops its future, running `ModulesFetchGuard::Drop`
        // and releasing `.modules.lock`.
        if let Some(prev) = modules_fetch.take() {
            prev.abort_tasks();
        }
        *modules_fetch = Some(ModulesFetchInFlight {
            progress_rx: event_rx,
            completed: 0,
            total: 0,
            current_id: String::new(),
            eta_secs: None,
            fetch_handle,
            relay_handle,
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

    // Cancel any in-flight setup-screen fetch before leaving the state, so a
    // still-streaming other-mode fetch (or a refresh) is not orphaned by the
    // transition. Covers the cache-missing arm below, which rebuilds setup
    // with both fetch slots `None`.
    if crate::app::is_busy(&app.state) {
        info!("stopping the in-flight setup fetch: starting the enrichment run");
    }
    crate::app::abort_in_flight(&app.state);
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
    app.spawn_stage3_run(crate::app::RunSpawn {
        payload: crate::app::RunPayloadSpec::Enrichment {
            dam_results: dam_results_clone,
            params,
        },
        target: target_clone,
        pubchem_total,
    });
}

#[cfg(test)]
mod build_run_inputs_tests {
    use super::*;
    use crate::dam::DamMethod;
    use crate::dam::fdr::FdrMethod;
    use crate::dam::types::{DamFeature, FcBasis};
    use crate::kegg::{KeggCompoundSet, SpeciesKegg};

    #[test]
    fn min_group_overlap_eff_max_soft_caps_and_preserves() {
        // Default / any ≤20 value under a large group → soft cap 20.
        assert_eq!(min_group_overlap_eff_max(157, 1), 20);
        // Small group → capped at the group's actual organism count (the
        // domain limit); a stale larger stored value then clamps down to it.
        assert_eq!(min_group_overlap_eff_max(4, 12), 4);
        assert_eq!(12usize.clamp(1, min_group_overlap_eff_max(4, 12)), 4);
        // A hand-set value above the soft cap under a big group is PRESERVED,
        // not silently truncated to 20 (review M2).
        assert_eq!(min_group_overlap_eff_max(100, 25), 25);
        assert_eq!(25usize.clamp(1, min_group_overlap_eff_max(100, 25)), 25);
        // No Group selected/fetched (usize::MAX) → only the soft cap applies.
        assert_eq!(min_group_overlap_eff_max(usize::MAX, 1), 20);
        assert_eq!(min_group_overlap_eff_max(usize::MAX, 25), 25);
    }

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

    #[test]
    fn reselect_species_aborts_and_clears_prior_inflight_fetch() {
        // Regression (found in smoke testing: select `csab` then `hsa`, csab's
        // progress bar kept running): switching species while a fetch is in
        // flight must cancel it and clear the progress strip on EVERY path —
        // including ones that never reach `spawn_species_fetch`'s own abort
        // (the cache fast-path, or — as exercised here — the early-return when
        // inputs are absent). The clear lives in `handle_species_selected`.
        let parked_rt = tokio::runtime::Builder::new_multi_thread()
            .worker_threads(1)
            .enable_all()
            .build()
            .expect("parked-task runtime");
        let fj = parked_rt.spawn(std::future::pending::<()>());
        let rj = parked_rt.spawn(std::future::pending::<()>());
        let fetch_handle = fj.abort_handle();
        let relay_handle = rj.abort_handle();
        let (_tx, progress_rx) = std::sync::mpsc::channel::<KeggEvent>();

        let app_rt = tokio::runtime::Builder::new_current_thread()
            .build()
            .expect("app runtime");
        let mut app = App::new(
            crate::logging::LogStore::new(16),
            "info".to_string(),
            app_rt,
            None,
        );
        app.settings.kegg_species = Some("csab".to_string());
        app.state = AppState::Stage3EnrichSetup {
            dam_results: vec![],
            error: None,
            kegg_fetch: Some(KeggFetchInFlight {
                progress_rx,
                completed: 0,
                total: 0,
                current_pathway: String::new(),
                fetch_handle,
                relay_handle,
            }),
            modules_fetch: None,
        };

        // Switch to a different, uncached species.
        //
        // This test used to lean on `spawn_species_fetch` early-returning
        // because `inputs` was empty, which isolated the clear by leaving the
        // slot at `None`. That precondition is gone: the guard is now "does a
        // screen own the slot" (`is_target_setup`), which this state satisfies,
        // so a fetch for the NEW species is installed. The load-bearing
        // assertion was always the abort, not the emptiness — the regression
        // this test was written for was csab's progress bar still running after
        // switching to hsa.
        handle_species_selected(&mut app, "zzqx_uncached_code".to_string());

        // Both of the PRIOR fetch's tasks were aborted, so nothing from csab
        // keeps streaming into the strip.
        assert!(parked_rt.block_on(fj).unwrap_err().is_cancelled());
        assert!(parked_rt.block_on(rj).unwrap_err().is_cancelled());
        // The screen still owns the slot; what it now holds is the newly
        // selected species' fetch, not the stale one whose handles just died.
        assert!(
            crate::app::setup_fetch_slots(&app.state).is_some(),
            "the setup screen must still own its fetch slots"
        );
        assert_eq!(
            app.settings.kegg_species.as_deref(),
            Some("zzqx_uncached_code"),
            "the new selection is recorded"
        );
    }

    /// A `Stage3EnrichRunning` state — owns no in-flight fetch slot and no
    /// `error` field, which is what makes it the one reachable place a spawn
    /// can be requested with nowhere to put the result.
    fn running_state(rt: &tokio::runtime::Runtime) -> AppState {
        fn dummy_rx<T>() -> std::sync::mpsc::Receiver<T> {
            let (_tx, rx) = std::sync::mpsc::channel();
            rx
        }
        AppState::Stage3EnrichRunning {
            payload: crate::app::RunningPayload::Coverage,
            phase: crate::app::Stage3Phase::PubChem,
            pubchem_progress_rx: dummy_rx(),
            kegg_conv_progress_rx: dummy_rx(),
            result_rx: dummy_rx(),
            pubchem_completed: 0,
            pubchem_total: 0,
            kegg_conv_completed: 0,
            kegg_conv_total: 0,
            run_handle: rt.spawn(std::future::pending::<()>()).abort_handle(),
        }
    }

    /// A Group pick with NO metadata `.csv` must still start the fetch.
    ///
    /// `inputs.ion_tables` is left EMPTY on purpose. Under the old guard
    /// (`ion_tables.is_empty() || mapping.is_none()`) that alone forced the
    /// early return, so this test also fails against a half-fix that removes
    /// only the `mapping` half — which is exactly the signal wanted. A KEGG
    /// catalogue fetch reads neither the samples nor the mapping.
    #[test]
    fn group_selection_without_a_mapping_installs_the_fetch_slot() {
        let rt = parked_app_rt();
        let mut app = App::new(
            crate::logging::LogStore::new(16),
            "info".to_string(),
            rt,
            None,
        );
        app.state = AppState::Stage2CoverageSetup {
            error: None,
            stale_groups_notice: None,
            kegg_fetch: None,
            modules_fetch: None,
        };
        assert!(app.inputs.mapping.is_none() && app.inputs.ion_tables.is_empty());

        spawn_modules_fetch(
            &mut app,
            2,
            "Plants".to_string(),
            std::collections::HashSet::new(),
            false,
        );

        let (_, modules_fetch) =
            crate::app::setup_fetch_slots(&app.state).expect("coverage setup owns the slots");
        assert!(
            modules_fetch.is_some(),
            "the fetch must start without a metadata CSV — the coverage route makes it optional"
        );
        assert_eq!(app.settings.organism_group.as_deref(), Some("Plants"));
    }

    /// A spawn requested from a state that owns no slot must start nothing.
    ///
    /// Without the precondition the tasks are spawned first and the slot
    /// installed afterwards inside an `if let` with no `else`, so they run
    /// detached: undrained, their `AbortHandle`s never stored, holding
    /// `.modules.lock` for the whole fetch.
    #[test]
    fn a_spawn_from_a_slotless_state_starts_nothing() {
        let rt = parked_app_rt();
        let mut app = App::new(
            crate::logging::LogStore::new(16),
            "info".to_string(),
            rt,
            None,
        );
        let handle_rt = parked_app_rt();
        app.state = running_state(&handle_rt);

        spawn_modules_fetch(
            &mut app,
            2,
            "Plants".to_string(),
            std::collections::HashSet::new(),
            false,
        );

        assert!(
            crate::app::setup_fetch_slots(&app.state).is_none(),
            "precondition of the test: this state owns no slot"
        );
        assert!(
            app.settings.organism_group.is_none(),
            "the guard must return BEFORE recording the selection, i.e. before spawning"
        );
    }

    /// A species refresh from a slotless state must not destroy the cache.
    ///
    /// `handle_species_refresh` invalidates the on-disk cache and clears
    /// `cache.species_kegg` before delegating to the spawn, and its own
    /// in-flight guard short-circuits to false off a setup screen. A
    /// precondition living only inside `spawn_species_fetch` would leave the
    /// user with neither cache and no fetch to refill them.
    #[test]
    fn a_species_refresh_from_a_slotless_state_keeps_the_cache() {
        let rt = parked_app_rt();
        let mut app = App::new(
            crate::logging::LogStore::new(16),
            "info".to_string(),
            rt,
            None,
        );
        let handle_rt = parked_app_rt();
        app.state = running_state(&handle_rt);
        app.settings.kegg_species = Some("ath".to_string());
        app.cache.species_kegg = Some(crate::kegg::SpeciesKegg {
            code: "ath".into(),
            fetched_at: chrono::Utc::now(),
            pathways: vec![],
        });

        handle_species_refresh(&mut app);

        assert!(
            app.cache.species_kegg.is_some(),
            "the in-memory catalogue must survive a refusal — it is not replaced by anything"
        );
    }

    fn parked_app_rt() -> tokio::runtime::Runtime {
        tokio::runtime::Builder::new_multi_thread()
            .worker_threads(1)
            .enable_all()
            .build()
            .expect("parked-task runtime")
    }

    fn empty_app() -> App {
        let app_rt = tokio::runtime::Builder::new_current_thread()
            .build()
            .expect("app runtime");
        App::new(
            crate::logging::LogStore::new(16),
            "info".to_string(),
            app_rt,
            None,
        )
    }

    #[test]
    fn mode_switch_blanks_incomplete_species_selection() {
        // Requirement: switching mode while the pathway species fetch is still
        // running must cancel it AND reset the (un-fetched) species selection
        // to blank — there is no completed cache to preserve.
        let rt = parked_app_rt();
        let fj = rt.spawn(std::future::pending::<()>());
        let rj = rt.spawn(std::future::pending::<()>());
        let (fetch_handle, relay_handle) = (fj.abort_handle(), rj.abort_handle());
        let (_tx, progress_rx) = std::sync::mpsc::channel::<KeggEvent>();

        let mut app = empty_app();
        app.settings.analysis_mode = AnalysisMode::Pathway;
        app.settings.kegg_species = Some("csab".to_string());
        app.state = AppState::Stage3EnrichSetup {
            dam_results: vec![],
            error: None,
            kegg_fetch: Some(KeggFetchInFlight {
                progress_rx,
                completed: 0,
                total: 0,
                current_pathway: String::new(),
                fetch_handle,
                relay_handle,
            }),
            modules_fetch: None,
        };

        cancel_inflight_for_mode_switch(&mut app, AnalysisMode::Module);

        assert_eq!(app.settings.kegg_species, None, "selection blanked");
        assert!(app.cache.species_kegg.is_none());
        assert!(matches!(
            &app.state,
            AppState::Stage3EnrichSetup {
                kegg_fetch: None,
                ..
            }
        ));
        assert!(rt.block_on(fj).unwrap_err().is_cancelled());
        assert!(rt.block_on(rj).unwrap_err().is_cancelled());
    }

    #[test]
    fn mode_switch_blanks_incomplete_group_selection() {
        let rt = parked_app_rt();
        let fj = rt.spawn(std::future::pending::<()>());
        let rj = rt.spawn(std::future::pending::<()>());
        let (fetch_handle, relay_handle) = (fj.abort_handle(), rj.abort_handle());
        let (_tx, progress_rx) = std::sync::mpsc::channel::<crate::app::ModulesFetchEvent>();

        let mut app = empty_app();
        app.settings.analysis_mode = AnalysisMode::Module;
        app.settings.organism_group_level = Some(3);
        app.settings.organism_group = Some("Mammals".to_string());
        app.organism_group_selector.level = 3;
        app.state = AppState::Stage3EnrichSetup {
            dam_results: vec![],
            error: None,
            kegg_fetch: None,
            modules_fetch: Some(ModulesFetchInFlight {
                progress_rx,
                completed: 0,
                total: 0,
                current_id: String::new(),
                eta_secs: None,
                fetch_handle,
                relay_handle,
            }),
        };

        cancel_inflight_for_mode_switch(&mut app, AnalysisMode::Pathway);

        assert_eq!(app.settings.organism_group, None, "group blanked");
        assert_eq!(app.settings.organism_group_level, None, "level blanked");
        assert!(app.cache.group_org_codes.is_none());
        assert_eq!(
            app.organism_group_selector.level, 2,
            "selector level reset to default"
        );
        assert!(matches!(
            &app.state,
            AppState::Stage3EnrichSetup {
                modules_fetch: None,
                ..
            }
        ));
        assert!(rt.block_on(fj).unwrap_err().is_cancelled());
        assert!(rt.block_on(rj).unwrap_err().is_cancelled());
    }

    #[test]
    fn mode_switch_preserves_completed_selection() {
        // Coexistence contract: with NO fetch in flight (completed cache), the
        // selection + cache MUST survive the toggle untouched.
        let mut app = empty_app();
        app.settings.kegg_species = Some("hsa".to_string());
        app.cache.species_kegg = Some(crate::kegg::SpeciesKegg {
            code: "hsa".into(),
            fetched_at: chrono::Utc::now(),
            pathways: vec![],
        });
        app.state = AppState::Stage3EnrichSetup {
            dam_results: vec![],
            error: None,
            kegg_fetch: None,
            modules_fetch: None,
        };

        cancel_inflight_for_mode_switch(&mut app, AnalysisMode::Module);

        assert_eq!(app.settings.kegg_species.as_deref(), Some("hsa"));
        assert!(app.cache.species_kegg.is_some());
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
