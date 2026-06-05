//! Global stage stepper / breadcrumb — the `Input › DAM Setup › DAM Result ›
//! Enrichment Analysis › Enrichment Result` row from every `UI-design.md`
//! mockup. Rendered once in a `TopBottomPanel::top` above the `CentralPanel`.
//!
//! The current stage is **bold**; already-reached stages are clickable buttons
//! that jump BACK to that stage; not-yet-reached stages are plain
//! (non-interactive) labels — NOT greyed/faded. You cannot jump forward (a
//! forward stage requires running the pipeline). Because the stepper carries
//! all back-navigation, the per-screen `< Back to …` buttons are removed.

use crate::app::{App, AppState, slot_fields_from};
use crate::dam::DamResult;
use crate::data::IonMode;

/// The five stages, in order. Index = step number (0..=4).
const STEP_LABELS: [&str; 5] = [
    "Input",
    "DAM Setup",
    "DAM Result",
    "Enrichment Analysis",
    "Enrichment Result",
];

/// Map the current `AppState` to its step index (0..=4). Transient running
/// states map to their setup step. `None` during `Initializing` (no stepper).
fn current_step(state: &AppState) -> Option<usize> {
    match state {
        AppState::Initializing { .. } => None,
        AppState::Stage1Input { .. } => Some(0),
        AppState::Stage2DamSetup { .. } | AppState::Stage2DamRunning { .. } => Some(1),
        AppState::Stage2DamThreshold { .. } => Some(2),
        AppState::Stage3EnrichSetup { .. } | AppState::Stage3EnrichRunning { .. } => Some(3),
        AppState::Stage3EnrichResult { .. } => Some(4),
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

/// Map a stepper step index (0..=4) to its icon texture: `0 → input`,
/// `1..=3 → intermediate` (shared), `4 → eating`. The `_` arm folds onto the
/// step-4 asset; the only real inputs are 0..=4 (one per `STEP_LABELS` entry).
fn icon_for_step(step: usize, icons: &StepperIcons) -> egui::TextureId {
    match step {
        0 => icons.input.id(),
        1..=3 => icons.intermediate.id(),
        _ => icons.eating.id(),
    }
}

/// Render the stepper row and perform any back-navigation the user clicked.
/// No-op (renders nothing) during `Initializing`.
pub fn show(ui: &mut egui::Ui, app: &mut App) {
    let Some(current) = current_step(&app.state) else {
        return;
    };

    // Lazily upload the icons on the first render (texture upload needs a live
    // Context, unavailable at `App::new`). Done BEFORE painting any segment so
    // the current step's icon is present on the very first frame — no flash-in
    // (design D4 / H2). `current_icon` is a `Copy` `TextureId`, so the
    // `app.stepper_icons` borrow ends here and `app` stays free for navigation.
    let icons = app
        .stepper_icons
        .get_or_insert_with(|| StepperIcons::load(ui.ctx()));
    let current_icon = icon_for_step(current, icons);

    // Pin every segment to the icon-row height so the row stays vertically
    // aligned: the current (icon) segment is taller than the others, and egui's
    // single-pass horizontal layout would otherwise shift the segments rendered
    // BEFORE it upward (see `segmented_icon_row_height`).
    let row_h = crate::ui::widgets::segmented_icon_row_height(ui);

    let mut jump_to: Option<usize> = None;
    // §4 Primary "Segmented" control — the steps sit in a FILL_SECONDARY track;
    // the current step is the white "floating slider".
    crate::ui::widgets::segmented_track(ui, |ui| {
        for (i, label) in STEP_LABELS.iter().enumerate() {
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
        | AppState::Stage3EnrichRunning { dam_results, .. }
        | AppState::Stage3EnrichResult { dam_results, .. } => Some(dam_results),
        _ => None,
    }
}

/// Jump back to an earlier step, reconstructing the target `AppState` from the
/// persistent siblings (`settings` / `inputs` / `cache`) + the `dam_results`
/// carried on the current state. Forward runtime artifacts are dropped, exactly
/// like the (now-removed) per-screen Back buttons. `target` is always a step
/// strictly before the current one (the caller only enables past steps), and
/// any target that needs `dam_results` is only reachable from a state that
/// carries them, so the `expect` below is unreachable by construction.
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

    match target {
        // Input — settings reset is a documented no-op; slot radios are
        // re-derived from the loaded ion tables so mode choices persist.
        0 => {
            app.settings.reset_for_back_to_stage1();
            let (slot1_mode, slot2_revealed, slot2_mode) = slot_fields_from(&app.inputs.ion_tables);
            app.state = AppState::Stage1Input {
                slot1_mode,
                slot2_revealed,
                slot2_mode,
                error: None,
            };
        }
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
        // stage); restore the taken state defensively if somehow requested.
        _ => {
            debug_assert!(
                false,
                "navigate_back_to called with non-back target {target}"
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// `current_step` maps each `AppState` to its stepper index. Covers the
    /// cheaply-constructible canonical states (heavy fields use `vec![]` /
    /// `None`), exercising every distinct step value None / 0 / 1 / 2 / 3. The
    /// `*Running` aliases and `Stage3EnrichResult` (step 4) reuse the same
    /// trivial match arms but need a full egui/orchestrator harness to build;
    /// the match is exhaustive so a dropped arm fails to compile rather than
    /// silently mis-mapping.
    #[test]
    fn current_step_maps_states_to_indices() {
        let (_tx, load_rx) = std::sync::mpsc::channel();
        let initializing = AppState::Initializing {
            load_rx,
            fallback_cache: None,
            last_error: None,
        };
        assert_eq!(current_step(&initializing), None);

        let stage1 = AppState::Stage1Input {
            slot1_mode: None,
            slot2_revealed: false,
            slot2_mode: None,
            error: None,
        };
        assert_eq!(current_step(&stage1), Some(0));

        let stage2_setup = AppState::Stage2DamSetup { error: None };
        assert_eq!(current_step(&stage2_setup), Some(1));

        let stage2_threshold = AppState::Stage2DamThreshold {
            dam_results: vec![],
            active_volcano_tab: IonMode::Positive,
            volcano_textures: vec![],
            rendering: false,
            render_rx: None,
        };
        assert_eq!(current_step(&stage2_threshold), Some(2));

        let stage3_setup = AppState::Stage3EnrichSetup {
            dam_results: vec![],
            error: None,
            kegg_fetch: None,
            modules_fetch: None,
        };
        assert_eq!(current_step(&stage3_setup), Some(3));
    }

    /// `icon_for_step` maps each stepper index to its asset texture: 0 → input,
    /// 1–3 → the shared intermediate icon, 4 → eating. The three textures are
    /// distinct so a mis-map can't accidentally pass.
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
        assert_eq!(icon_for_step(0, &icons), icons.input.id());
        assert_eq!(icon_for_step(1, &icons), icons.intermediate.id());
        assert_eq!(icon_for_step(2, &icons), icons.intermediate.id());
        assert_eq!(icon_for_step(3, &icons), icons.intermediate.id());
        assert_eq!(icon_for_step(4, &icons), icons.eating.id());
        // The three icons really are distinct textures.
        assert_ne!(icons.input.id(), icons.intermediate.id());
        assert_ne!(icons.intermediate.id(), icons.eating.id());
    }
}
