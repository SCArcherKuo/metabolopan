//! Global stage stepper / breadcrumb — the `Input › DAM Setup › DAM Result ›
//! Enrichment Analysis › Enrichment Result` row from every `UI-design.md`
//! mockup. Rendered once in a `TopBottomPanel::top` above the `CentralPanel`.
//!
//! The current stage is **bold**; already-reached stages are clickable buttons
//! that jump BACK to that stage; not-yet-reached stages are plain
//! (non-interactive) labels — NOT greyed/faded. You cannot jump forward (a
//! forward stage requires running the pipeline). Because the stepper carries
//! all back-navigation, the per-screen `< Back to …` buttons are removed.

use crate::app::{AnalysisRoute, App, AppState, RunningPayload, slot_fields_from};
use crate::dam::DamResult;
use crate::data::IonMode;

/// The `DamEnrichment` track — five stages, in order. Index = step number.
const STEP_LABELS_DAM: [&str; 5] = [
    "Input",
    "DAM Setup",
    "DAM Result",
    "Enrichment Analysis",
    "Enrichment Result",
];

/// The `KeggCoverage` track — three stages. The coverage route genuinely has
/// fewer stages (no threshold screen, no separate target-selection stage), so a
/// shorter track is the truthful rendering rather than padding it to five.
const STEP_LABELS_COVERAGE: [&str; 3] = ["Input", "Setup", "Coverage"];

/// The active route's step labels.
fn track_labels(route: AnalysisRoute) -> &'static [&'static str] {
    match route {
        AnalysisRoute::DamEnrichment => &STEP_LABELS_DAM,
        AnalysisRoute::KeggCoverage => &STEP_LABELS_COVERAGE,
    }
}

/// Map the current `AppState` to its step index within the ACTIVE route's
/// track. Transient running states map to their setup step. `None` means no
/// stepper is rendered.
///
/// `Stage1Input` and `Stage3EnrichRunning` are shared by both routes, so the
/// variant alone cannot decide — the index comes from `route` (equivalently,
/// for the running screen, from its `RunningPayload` discriminant).
///
/// A variant belonging to the OTHER route returns `None` rather than an index
/// into a track it is not on. That pairing is unreachable — a route is chosen
/// before Stage 1 and can only change by returning to the chooser, which
/// discards the state — so the arm is a backstop: hiding the stepper is a
/// visible symptom, whereas an index into the wrong track would render a
/// plausible-looking row whose clicks navigate somewhere else entirely.
fn current_step(state: &AppState, route: AnalysisRoute) -> Option<usize> {
    match state {
        // Pre-route screens: no navigation is possible before the app is ready
        // and before a route has been chosen.
        AppState::Initializing { .. } | AppState::Stage0ChooseAnalysis => None,
        AppState::Stage1Input { .. } => Some(0),
        AppState::Stage3EnrichRunning { payload, .. } => match payload {
            RunningPayload::Enrichment(_) => Some(3),
            RunningPayload::Coverage => Some(1),
        },
        AppState::Stage2DamSetup { .. } | AppState::Stage2DamRunning { .. } => {
            matches!(route, AnalysisRoute::DamEnrichment).then_some(1)
        }
        AppState::Stage2DamThreshold { .. } => {
            matches!(route, AnalysisRoute::DamEnrichment).then_some(2)
        }
        AppState::Stage3EnrichSetup { .. } => {
            matches!(route, AnalysisRoute::DamEnrichment).then_some(3)
        }
        AppState::Stage3EnrichResult { .. } => {
            matches!(route, AnalysisRoute::DamEnrichment).then_some(4)
        }
        AppState::Stage2CoverageSetup { .. } => {
            matches!(route, AnalysisRoute::KeggCoverage).then_some(1)
        }
        AppState::Stage3CoverageResult { .. } => {
            matches!(route, AnalysisRoute::KeggCoverage).then_some(2)
        }
    }
}

/// The three decoded stepper icons, uploaded once to GPU textures and cached on
/// `App` for the session (see [`StepperIcons::load`]). `input` is the step-0
/// (`Input`) icon; `intermediate` is shared by steps 1–3 (`DAM Setup` /
/// `DAM Result` / `Enrichment Analysis`); `eating` is the step-4
/// (`Enrichment Result`) icon. No `Debug` (egui `TextureHandle` lacks it, same
/// as `App`).
pub struct StepperIcons {
    input: egui::TextureHandle,
    intermediate: egui::TextureHandle,
    eating: egui::TextureHandle,
}

impl StepperIcons {
    /// Decode the three embedded PNGs and upload them as textures. Called once,
    /// lazily, on the first stepper render — texture upload needs a live
    /// `egui::Context`, unavailable at `App::new`.
    fn load(ctx: &egui::Context) -> Self {
        Self {
            input: load_icon_texture(
                ctx,
                "stepper_icon_input",
                include_bytes!("../../assets/icon.png"),
            ),
            intermediate: load_icon_texture(
                ctx,
                "stepper_icon_intermediate",
                include_bytes!("../../assets/rat_single.png"),
            ),
            eating: load_icon_texture(
                ctx,
                "stepper_icon_eating",
                include_bytes!("../../assets/rat_eating.png"),
            ),
        }
    }
}

/// Decode an embedded RGBA PNG and upload it as a linear-filtered texture — the
/// same `image` → `ColorImage` → `load_texture` path the volcano / dot-plot
/// renderers and the `main.rs` window icon use.
fn load_icon_texture(ctx: &egui::Context, name: &str, bytes: &[u8]) -> egui::TextureHandle {
    let image = image::load_from_memory(bytes)
        .expect("embedded stepper icon PNG should decode")
        .to_rgba8();
    let (w, h) = (image.width() as usize, image.height() as usize);
    let color = egui::ColorImage::from_rgba_unmultiplied([w, h], image.as_raw());
    ctx.load_texture(name, color, egui::TextureOptions::LINEAR)
}

/// Map a step index within a `track_len`-long track to its icon texture:
/// first step → `input`, last step → `eating`, everything between →
/// `intermediate`.
///
/// Stated positionally rather than as a fixed 0/1..=3/4 table because that ONE
/// rule reproduces both tracks: the five-step DAM track keeps its historical
/// `0 → input`, `1..=3 → intermediate`, `4 → eating` mapping, and the three-step
/// coverage track gets `0 → input`, `1 → intermediate`, `2 → eating` — the same
/// "input / working / done" progression, reusing the same three assets. No new
/// asset is introduced.
///
/// A `track_len` of 0 is not reachable (both tracks are const arrays), and
/// `step` is always a valid index into the active track.
fn icon_for_step(step: usize, track_len: usize, icons: &StepperIcons) -> egui::TextureId {
    if step == 0 {
        icons.input.id()
    } else if step + 1 >= track_len {
        icons.eating.id()
    } else {
        icons.intermediate.id()
    }
}

/// Render the stepper row and perform any back-navigation the user clicked.
/// No-op (renders nothing) during `Initializing` and on the route chooser.
pub fn show(ui: &mut egui::Ui, app: &mut App) {
    let route = app.settings.analysis_route;
    let Some(current) = current_step(&app.state, route) else {
        return;
    };
    let labels = track_labels(route);

    // Lazily upload the icons on the first render (texture upload needs a live
    // Context, unavailable at `App::new`). Done BEFORE painting any segment so
    // the current step's icon is present on the very first frame — no flash-in
    // (design D4 / H2). `current_icon` is a `Copy` `TextureId`, so the
    // `app.stepper_icons` borrow ends here and `app` stays free for navigation.
    let icons = app
        .stepper_icons
        .get_or_insert_with(|| StepperIcons::load(ui.ctx()));
    let current_icon = icon_for_step(current, labels.len(), icons);

    // Pin every segment to the icon-row height so the row stays vertically
    // aligned: the current (icon) segment is taller than the others, and egui's
    // single-pass horizontal layout would otherwise shift the segments rendered
    // BEFORE it upward (see `segmented_icon_row_height`).
    let row_h = crate::ui::widgets::segmented_icon_row_height(ui);

    let mut jump_to: Option<usize> = None;
    // §4 Primary "Segmented" control — the steps sit in a FILL_SECONDARY track;
    // the current step is the white "floating slider".
    crate::ui::widgets::segmented_track(ui, |ui| {
        for (i, label) in labels.iter().enumerate() {
            if i > 0 {
                // ASCII separator — the default egui font lacks the ›/→ glyphs
                // (same reason the back-arrow used `<`), so use ">".
                ui.label(">");
            }
            // Current step carries its stage icon; past/future are icon-less.
            // Past (i < current) is clickable (jump back); current is selected
            // but inert; future (i > current) is disabled. Every segment is
            // pinned to `row_h` so they align.
            let icon = (i == current).then_some(current_icon);
            let enabled = i <= current;
            let resp = crate::ui::widgets::segmented_tab_with_icon(
                ui,
                i == current,
                label,
                enabled,
                true,
                icon,
                row_h,
            );
            if i < current && resp.clicked() {
                jump_to = Some(i);
            }
        }
    });

    if let Some(target) = jump_to {
        if crate::app::needs_nav_confirm(&app.state) {
            // A long module fetch is in flight: defer the jump behind a
            // confirm modal instead of cancelling silently. Respect the
            // App-level modal mutual-exclusion family.
            if app.any_modal_open() {
                tracing::warn!("back-navigation request dropped: another modal is open");
            } else {
                app.pending_back_nav = Some(target);
            }
        } else {
            navigate_back_to(app, target);
        }
    }
}

/// Consume the current state and pull out its `dam_results` if the variant
/// carries them (DAM result + every Stage 3 state). `None` for states that
/// never hold a completed DAM result.
fn take_dam_results(prev: AppState) -> Option<Vec<DamResult>> {
    match prev {
        AppState::Stage2DamThreshold { dam_results, .. }
        | AppState::Stage3EnrichSetup { dam_results, .. }
        | AppState::Stage3EnrichResult { dam_results, .. } => Some(dam_results),
        // Shared running screen: only the enrichment payload carries results.
        AppState::Stage3EnrichRunning {
            payload: RunningPayload::Enrichment(dam_results),
            ..
        } => Some(dam_results),
        _ => None,
    }
}

/// Jump back to an earlier step, reconstructing the target `AppState` from the
/// persistent siblings (`settings` / `inputs` / `cache`) + the `dam_results`
/// carried on the current state. Forward runtime artifacts are dropped, exactly
/// like the (now-removed) per-screen Back buttons. `target` is always a step
/// strictly before the current one (the caller only enables past steps), and
/// any target that needs `dam_results` is only reachable from a state that
/// carries them, so the `expect`s below are unreachable by construction.
///
/// **A step index means nothing without the route.** Step 1 is `DAM Setup` on
/// one track and `Setup` (coverage) on the other, so the two tracks get separate
/// dispatch tables. A single route-blind table is exactly how a coverage
/// back-jump would land on the DAM screen: `1 => Stage2DamSetup` reads correct
/// until a second track exists.
pub(crate) fn navigate_back_to(app: &mut App, target: usize) {
    // Cancel any in-flight background task on the leaving state BEFORE it is
    // dropped, so a running enrichment/DAM run or a fetch/refresh is aborted
    // rather than orphaned. Silent for every state except a long module
    // fetch, which is gated by a confirm before this is reached (see `show`
    // and `App::show_back_nav_confirm_modal`).
    if crate::app::is_busy(&app.state) {
        tracing::info!(
            target_step = target,
            "stopping in-flight work: navigated back to an earlier step via the stepper"
        );
    }
    crate::app::abort_in_flight(&app.state);
    let prev = std::mem::take(&mut app.state);
    let dam_results = take_dam_results(prev);

    match app.settings.analysis_route {
        AnalysisRoute::DamEnrichment => back_to_dam_step(app, target, dam_results),
        AnalysisRoute::KeggCoverage => back_to_coverage_step(app, target),
    }
}

/// Back-navigation targets on the five-step `DamEnrichment` track.
fn back_to_dam_step(app: &mut App, target: usize, dam_results: Option<Vec<DamResult>>) {
    match target {
        // Input — settings reset is a documented no-op; slot radios are
        // re-derived from the loaded ion tables so mode choices persist.
        0 => back_to_input(app),
        // DAM Setup.
        1 => {
            app.settings.reset_stage2_choices_on_change_comparison();
            app.state = AppState::Stage2DamSetup { error: None };
        }
        // DAM Result — rebuild the threshold screen from the carried results.
        2 => {
            app.settings.reset_for_back_to_stage2_threshold();
            let dam_results = dam_results.expect("DAM Result reachable only with dam_results");
            let active_volcano_tab = app
                .inputs
                .ion_tables
                .first()
                .map(|it| it.mode)
                .unwrap_or(IonMode::Positive);
            let volcano_textures = vec![None; app.inputs.ion_tables.len()];
            app.state = AppState::Stage2DamThreshold {
                dam_results,
                active_volcano_tab,
                volcano_textures,
                rendering: false,
                render_rx: None,
            };
        }
        // Enrichment Analysis (setup).
        3 => {
            let dam_results =
                dam_results.expect("Enrichment Setup reachable only with dam_results");
            app.state = AppState::Stage3EnrichSetup {
                dam_results,
                error: None,
                kegg_fetch: None,
                modules_fetch: None,
            };
            // Re-derive the active mode's KEGG cache from a restored selection
            // whose cache is empty (mirrors `continue_to_enrichment`).
            crate::ui::stage3_setup::rehydrate_stage3_cache(app);
        }
        // Step 4 (Enrichment Result) is never a back-target (it is the last
        // stage); the state was already taken, so leave it at its `Default`.
        _ => {
            debug_assert!(
                false,
                "navigate_back_to called with non-back target {target}"
            );
        }
    }
}

/// Back-navigation targets on the three-step `KeggCoverage` track, exhaustively.
///
/// Only `reset_for_back_to_stage1` has a corresponding target here.
/// `reset_stage2_choices_on_change_comparison` and
/// `reset_for_back_to_stage2_threshold` name DAM-route screens this track does
/// not contain, so they are deliberately NOT invoked — their absence is not an
/// omission. All three bodies are no-ops anyway, so every settings field
/// survives every coverage-route jump either way.
fn back_to_coverage_step(app: &mut App, target: usize) {
    match target {
        0 => back_to_input(app),
        // Setup. Reachable from `Stage3CoverageResult` only: the running screen
        // maps to step 1 itself, so `Setup` is its CURRENT step and clicking it
        // is a no-op that never reaches here.
        1 => {
            app.state = AppState::Stage2CoverageSetup {
                error: None,
                stale_groups_notice: None,
                kegg_fetch: None,
                modules_fetch: None,
            };
            // Same rehydrate as the enrichment setup target — the screen hosts
            // the same selector, so a restored selection with an empty cache
            // must reload for `Run Coverage` to be immediately clickable.
            crate::ui::stage3_setup::rehydrate_stage3_cache(app);
        }
        // Step 2 (Coverage result) is the last stage, never a back-target.
        _ => {
            debug_assert!(
                false,
                "navigate_back_to called with non-back target {target}"
            );
        }
    }
}

/// Step 0 on both tracks: back to `Stage1Input` with the slot radios re-derived
/// from the loaded ion tables, so the user's ionization-mode choices survive.
fn back_to_input(app: &mut App) {
    app.settings.reset_for_back_to_stage1();
    let (slot1_mode, slot2_revealed, slot2_mode) = slot_fields_from(&app.inputs.ion_tables);
    app.state = AppState::Stage1Input {
        slot1_mode,
        slot2_revealed,
        slot2_mode,
        error: None,
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::app::{CoverageFunnel, RunningPayload};
    use crate::coverage::CoverageResult;

    const DAM: AnalysisRoute = AnalysisRoute::DamEnrichment;
    const COV: AnalysisRoute = AnalysisRoute::KeggCoverage;

    fn initializing() -> AppState {
        let (_tx, load_rx) = std::sync::mpsc::channel();
        AppState::Initializing {
            load_rx,
            fallback_cache: None,
            last_error: None,
        }
    }

    fn stage1() -> AppState {
        AppState::Stage1Input {
            slot1_mode: None,
            slot2_revealed: false,
            slot2_mode: None,
            error: None,
        }
    }

    fn dam_threshold() -> AppState {
        AppState::Stage2DamThreshold {
            dam_results: vec![],
            active_volcano_tab: IonMode::Positive,
            volcano_textures: vec![],
            rendering: false,
            render_rx: None,
        }
    }

    fn enrich_setup() -> AppState {
        AppState::Stage3EnrichSetup {
            dam_results: vec![],
            error: None,
            kegg_fetch: None,
            modules_fetch: None,
        }
    }

    fn coverage_setup() -> AppState {
        AppState::Stage2CoverageSetup {
            error: None,
            stale_groups_notice: None,
            kegg_fetch: None,
            modules_fetch: None,
        }
    }

    fn coverage_result() -> AppState {
        AppState::Stage3CoverageResult {
            coverage_result: CoverageResult {
                rows: vec![],
                detected_total: 0,
                entries_total: 0,
                entries_without_compounds: 0,
                detected_in_entries: 0,
            },
            funnel: CoverageFunnel::default(),
            cpd_to_names: std::collections::HashMap::new(),
            module_retention: None,
            mode_partition: None,
            dedup_reports: vec![],
            pubchem_time_span: None,
            kegg_conv_time_span: None,
            dotplot_tex: None,
            rendering: false,
            render_rx: None,
            confirming_new_round: false,
            height_user_overridden: false,
        }
    }

    /// A `Stage3EnrichRunning` carrying `payload`. Needs a tokio runtime for the
    /// `AbortHandle`, which the step-index tests do not otherwise require.
    fn running(payload: RunningPayload, rt: &tokio::runtime::Runtime) -> AppState {
        let (_tx, rx1) = std::sync::mpsc::channel();
        let (_tx2, rx2) = std::sync::mpsc::channel();
        let (_tx3, rx3) = std::sync::mpsc::channel();
        AppState::Stage3EnrichRunning {
            payload,
            phase: crate::app::Stage3Phase::PubChem,
            pubchem_progress_rx: rx1,
            kegg_conv_progress_rx: rx2,
            result_rx: rx3,
            pubchem_completed: 0,
            pubchem_total: 0,
            kegg_conv_completed: 0,
            kegg_conv_total: 0,
            run_handle: rt.spawn(std::future::pending::<()>()).abort_handle(),
        }
    }

    fn test_rt() -> tokio::runtime::Runtime {
        tokio::runtime::Builder::new_current_thread()
            .build()
            .expect("test runtime")
    }

    /// The five-step DAM track, unchanged from before the route split.
    #[test]
    fn current_step_maps_dam_route_states_to_indices() {
        let rt = test_rt();
        assert_eq!(current_step(&initializing(), DAM), None);
        assert_eq!(current_step(&AppState::Stage0ChooseAnalysis, DAM), None);
        assert_eq!(current_step(&stage1(), DAM), Some(0));
        assert_eq!(
            current_step(&AppState::Stage2DamSetup { error: None }, DAM),
            Some(1)
        );
        assert_eq!(current_step(&dam_threshold(), DAM), Some(2));
        assert_eq!(current_step(&enrich_setup(), DAM), Some(3));
        assert_eq!(
            current_step(&running(RunningPayload::Enrichment(vec![]), &rt), DAM),
            Some(3)
        );
    }

    /// The three-step coverage track. `Stage1Input` and the shared running
    /// screen resolve through the route, not the variant.
    #[test]
    fn current_step_maps_coverage_route_states_to_indices() {
        let rt = test_rt();
        assert_eq!(current_step(&initializing(), COV), None);
        assert_eq!(current_step(&AppState::Stage0ChooseAnalysis, COV), None);
        assert_eq!(current_step(&stage1(), COV), Some(0));
        assert_eq!(current_step(&coverage_setup(), COV), Some(1));
        assert_eq!(
            current_step(&running(RunningPayload::Coverage, &rt), COV),
            Some(1)
        );
        assert_eq!(current_step(&coverage_result(), COV), Some(2));
    }

    /// The shared running screen resolves through its payload, so the SAME
    /// variant lands on step 3 of the five-step track or step 1 of the
    /// three-step one.
    #[test]
    fn the_shared_running_screen_maps_to_its_own_routes_step() {
        let rt = test_rt();
        assert_eq!(
            current_step(&running(RunningPayload::Enrichment(vec![]), &rt), DAM),
            Some(3)
        );
        assert_eq!(
            current_step(&running(RunningPayload::Coverage, &rt), COV),
            Some(1)
        );
    }

    /// A variant belonging to the other route hides the stepper rather than
    /// indexing into a track it is not on. Unreachable in practice; asserted so
    /// the failure mode stays "no stepper" rather than "wrong stepper".
    #[test]
    fn a_foreign_variant_renders_no_stepper() {
        assert_eq!(current_step(&coverage_setup(), DAM), None);
        assert_eq!(current_step(&coverage_result(), DAM), None);
        assert_eq!(
            current_step(&AppState::Stage2DamSetup { error: None }, COV),
            None
        );
        assert_eq!(current_step(&dam_threshold(), COV), None);
        assert_eq!(current_step(&enrich_setup(), COV), None);
    }

    /// Each route gets its own track, and the coverage track offers no DAM step
    /// — which is what makes a route switch impossible via back-navigation.
    #[test]
    fn each_route_has_its_own_track() {
        assert_eq!(
            track_labels(DAM),
            [
                "Input",
                "DAM Setup",
                "DAM Result",
                "Enrichment Analysis",
                "Enrichment Result"
            ]
        );
        assert_eq!(track_labels(COV), ["Input", "Setup", "Coverage"]);
        for label in track_labels(COV) {
            assert!(
                !track_labels(DAM).contains(label) || *label == "Input",
                "coverage track must not offer a DAM-only step: {label}"
            );
        }
    }

    /// `icon_for_step` is positional, so one rule serves both tracks: first →
    /// input, last → eating, middle → intermediate. The three textures are
    /// distinct so a mis-map cannot accidentally pass.
    #[test]
    fn icon_for_step_maps_indices_to_assets() {
        let ctx = egui::Context::default();
        let mk = |name: &str| {
            ctx.load_texture(
                name,
                egui::ColorImage::from_rgba_unmultiplied([1, 1], &[255, 255, 255, 255]),
                egui::TextureOptions::LINEAR,
            )
        };
        let icons = StepperIcons {
            input: mk("input"),
            intermediate: mk("intermediate"),
            eating: mk("eating"),
        };
        // Five-step DAM track — the historical 0 / 1..=3 / 4 mapping.
        assert_eq!(icon_for_step(0, 5, &icons), icons.input.id());
        assert_eq!(icon_for_step(1, 5, &icons), icons.intermediate.id());
        assert_eq!(icon_for_step(2, 5, &icons), icons.intermediate.id());
        assert_eq!(icon_for_step(3, 5, &icons), icons.intermediate.id());
        assert_eq!(icon_for_step(4, 5, &icons), icons.eating.id());
        // Three-step coverage track — same progression, same three assets.
        assert_eq!(icon_for_step(0, 3, &icons), icons.input.id());
        assert_eq!(icon_for_step(1, 3, &icons), icons.intermediate.id());
        assert_eq!(icon_for_step(2, 3, &icons), icons.eating.id());
        // The three icons really are distinct textures.
        assert_ne!(icons.input.id(), icons.intermediate.id());
        assert_ne!(icons.intermediate.id(), icons.eating.id());
    }

    /// `take_dam_results` reaches through the running screen's payload, and
    /// finds nothing on a coverage run — an `Enrichment(vec![])` fallback would
    /// have been indistinguishable from a real empty result.
    #[test]
    fn take_dam_results_reads_through_the_running_payload() {
        let rt = test_rt();
        assert!(take_dam_results(running(RunningPayload::Enrichment(vec![]), &rt)).is_some());
        assert!(take_dam_results(running(RunningPayload::Coverage, &rt)).is_none());
        assert!(take_dam_results(coverage_setup()).is_none());
        assert!(take_dam_results(coverage_result()).is_none());
        assert!(take_dam_results(AppState::Stage0ChooseAnalysis).is_none());
    }
}
