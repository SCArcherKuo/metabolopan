//! Stage 3 result screen — dot plot preview, export buttons, time-span
//! display, and two refresh buttons with confirmation modals.

use egui::RichText;
use std::sync::mpsc;
use tracing::{error, info, warn};

use crate::app::{
    AnalysisMode, App, AppState, DrawnPlot, EnrichDisplayFilters, RefreshState, Stage3RunOutput,
};
use crate::enrichment::{RowSelection, export_csv_with_mode};
use crate::plot::{DotplotOpts, export_dotplot_png, render_dotplot};
use crate::pubchem::PubchemClient;
use crate::stage3::run_stage3;
use crate::theme;

#[derive(Debug, Clone, Copy)]
enum Action {
    None,
    Redraw,
    DownloadPng,
    /// `Download enrichment results (CSV)` — the rows the figure is drawn from.
    DownloadFigureCsv,
    /// `Download all results (CSV)` — every surviving row.
    DownloadAllCsv,
    ConfirmRefreshPubchem,
    ConfirmRefreshKegg,
    CancelRefresh,
    RequestNewRound,
    ConfirmNewRound,
    CancelNewRound,
}

pub fn show(ui: &mut egui::Ui, app: &mut App) {
    drain_render(app, ui.ctx());
    drain_refresh(app);
    drain_cache_actions(app);

    // Snapshot read-only fields we need to render. Built by a pure function so
    // the run-vs-settings sourcing rule below is assertable without egui — see
    // the `stage3-ui` result-screen requirement.
    let snap = result_snap(&app.state, &app.settings);
    let Some(snap) = snap else {
        return;
    };
    show_inner(ui, app, snap);
}

/// Build the render snapshot from a result state.
///
/// **Method from the run, threshold from the controls.** `fdr_method` is read
/// from `EnrichmentResult`, the method that produced the numbers on display;
/// `fdr_threshold` is read from settings, because it is a live display filter
/// the user tunes on this screen. Sourcing the method from settings would let a
/// label describe an export it does not match — the exporter reads the result.
/// See the `stage3-ui` capability spec.
fn result_snap(state: &AppState, settings: &crate::app::SessionSettings) -> Option<ResultSnap> {
    match state {
        AppState::Stage3EnrichResult {
            dam_results,
            enrichment_result,
            feature_to_cpds,
            refresh_state,
            rendering,
            dotplot,
            confirming_new_round,
            height_user_overridden,
            ..
        } => {
            let mode = settings.analysis_mode;
            let enrichment_fdr_threshold = settings.enrichment_fdr_threshold;
            // Identity, cache freshness, the Refresh / Re-run buttons, and the
            // provenance counts all moved to the bottom-panel Data tab
            // (`data-summary-panel`). The result body keeps only the
            // display-filter inputs, export size, draw, preview, and downloads.
            Some(ResultSnap {
                mode,
                // Scope counts for the refresh-confirmation modal copy.
                dam_features_total: dam_results.iter().map(|d| d.features.len()).sum(),
                mapped_features: feature_to_cpds.len(),
                export_w_in: settings.stage3_export_width_in,
                export_h_in: settings.stage3_export_height_in,
                export_dpi: settings.stage3_export_dpi,
                top_n: settings.top_n,
                fdr_threshold: enrichment_fdr_threshold,
                // Method from the RUN, threshold from the controls.
                fdr_method: enrichment_result.fdr_method,
                min_hit_count: settings.min_hit_count,
                refresh_state_kind: refresh_state_kind(refresh_state),
                confirming_new_round: *confirming_new_round,
                rendering: *rendering,
                drawn_filters: dotplot.as_ref().map(|d| d.filters),
                height_user_overridden: *height_user_overridden,
            })
        }
        _ => None,
    }
}

/// Label for the significance-threshold DragValue.
///
/// Names the quantity the threshold is compared against, so it follows the
/// method that produced that quantity. Under `NoCorrection` the value compared
/// is a raw p-value; calling it an FDR would assert a correction that was not
/// performed. See the `stage3-ui` capability spec.
fn threshold_label(method: crate::dam::fdr::FdrMethod) -> String {
    format!("Enrichment {} threshold:", method.metric_label())
}

/// Hover text for the `Minimum hit count` DragValue.
///
/// The `after FDR` clause is DROPPED under `NoCorrection`, not re-pointed at
/// another stage: there is no FDR stage for the filter to follow, and naming a
/// different one would assert a second thing that did not happen.
fn min_hit_tooltip(method: crate::dam::fdr::FdrMethod) -> &'static str {
    match method {
        crate::dam::fdr::FdrMethod::NoCorrection => {
            "Hide pathways/modules with fewer than N hits. \
             Display-only; does not change p-values."
        }
        _ => {
            "Hide pathways/modules with fewer than N hits after FDR. \
             Display-only; does not change p-values."
        }
    }
}

/// **The one predicate.** Whether the held texture still describes the live
/// filter values — `None` meaning no texture is held, which is never live.
///
/// Three things read it: the draw button's label, the decision to blit, and the
/// not-yet-drawn prompt (which is templated FROM the label). They must read it
/// from one place; deriving any of them from raw texture-PRESENCE instead puts
/// `Click "Re-draw dot plot" to render the plot.` beside an empty preview on the
/// discard frame, offering to re-draw a figure that is no longer there.
fn texture_is_live(drawn: Option<EnrichDisplayFilters>, live: EnrichDisplayFilters) -> bool {
    drawn == Some(live)
}

/// The draw button's label, from the one predicate above.
fn draw_button_label(texture_is_live: bool) -> &'static str {
    if texture_is_live {
        "Re-draw dot plot"
    } else {
        "Draw dot plot"
    }
}

/// The not-yet-drawn prompt, templated from the button's label so the two can
/// never name different buttons.
fn not_yet_drawn_prompt(button_label: &str) -> String {
    format!("Click \"{button_label}\" to render the plot.")
}

fn show_inner(ui: &mut egui::Ui, app: &mut App, snap: ResultSnap) {
    let mut new_w_in = snap.export_w_in;
    let mut new_h_in = snap.export_h_in;
    let mut new_dpi = snap.export_dpi;
    let mut new_top_n = snap.top_n;
    let mut new_fdr_threshold = snap.fdr_threshold;
    let mut new_min_hit_count = snap.min_hit_count;
    let mut action = Action::None;
    // Whether the held texture still describes the live filter values. Computed
    // inside the body, once, above the draw button; read again after the body to
    // perform the discard. Declared here because the body is a closure, exactly
    // as `action` is.
    let mut texture_is_live = false;

    let refresh_inflight = !matches!(snap.refresh_state_kind, RefreshKind::Idle);
    let busy = refresh_inflight || snap.rendering;

    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            // Heading: mode-specific suffix (Title Case per Phase 4).
            let mode_suffix = match snap.mode {
                AnalysisMode::Pathway => "Pathway Mode",
                AnalysisMode::Module => "Module Mode",
            };
            ui.heading(
                egui::RichText::new(format!("Stage 3 — Enrichment Result · {mode_suffix}"))
                    .color(theme::HEADING),
            );
            ui.add_space(6.0);
            // Back-navigation handled by the global stage stepper (`ui::stepper`).
            ui.add_space(6.0);

            // The "Data sources for this run" panel — identity / cache
            // fetched-date spans / the `Refresh PubChem cache` + `Refresh KEGG
            // conv cache` buttons (with inline progress) / the `Re-run
            // enrichment` button — moved to the bottom-panel Data tab's `Cache
            // data` block (`data-summary-panel`). The Data-tab buttons set
            // `LogPaneState` request flags that `drain_cache_actions` (top of
            // this fn) feeds into the existing confirm-modal / re-run flow, so
            // the modal + `RefreshState` machinery below is unchanged.

            // Significance threshold + Minimum hit count — relocated from
            // Stage 3 setup (`add-bottom-panel-data-tab`). Both are live display
            // filters: editing either DISCARDS the figure, and the next
            // `Draw dot plot` re-applies it without re-spawning the
            // orchestrator. Both reach one of the two CSV downloads — the
            // filtered one, whose definition is the figure's row set — and
            // neither reaches the other (see the `enrichment-ora` capability).
            // The threshold's label names the quantity it is compared against,
            // so it follows the RUN's method (`snap.fdr_method`), never
            // `app.settings.enrichment_fdr_method`.
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                ui.label(threshold_label(snap.fdr_method));
                ui.add(
                    egui::DragValue::new(&mut new_fdr_threshold)
                        .speed(0.001)
                        .range(0.0001..=1.0),
                );
            });
            ui.horizontal(|ui| {
                ui.label("Minimum hit count:");
                ui.add(
                    egui::DragValue::new(&mut new_min_hit_count)
                        .speed(1)
                        .range(1..=10),
                )
                .on_hover_text(min_hit_tooltip(snap.fdr_method));
            });

            // Top N input — moved from Stage 3 setup to Stage 3 result so
            // the user can iterate on the dot plot's display cap without
            // navigating back. Mode-aware label.
            ui.add_space(8.0);
            let top_n_label = match snap.mode {
                AnalysisMode::Pathway => "Top N pathways:",
                AnalysisMode::Module => "Top N modules:",
            };
            ui.horizontal(|ui| {
                ui.label(top_n_label);
                ui.add(egui::DragValue::new(&mut new_top_n).speed(1).range(1..=100));
            });

            ui.add_space(12.0);
            crate::ui::widgets::png_export_size_controls(
                ui,
                &mut new_w_in,
                &mut new_h_in,
                &mut new_dpi,
            );
            // The `-> W × H px` hint was removed (`apply-ui-design-md-tweaks`);
            // the pixel dimensions are implied by the Width / Height / DPI
            // inputs above. `export_pixels` is still used by the render + PNG
            // export sites.

            // ── The one comparison ──
            //
            // Every display filter this frame has been edited into the `new_*`
            // locals above; the write-back to `app.settings` happens BELOW the
            // preview, so `app.settings` still holds last frame's values here.
            // Comparing against the locals is what makes the discard land on
            // the SAME frame as the change, matching the coverage screen, whose
            // controls write settings directly above its preview.
            //
            // THREE things read the result — the button label, the decision to
            // blit, and the not-yet-drawn prompt, which is templated FROM the
            // label. They must read it from one place. Deriving the label from
            // raw texture-PRESENCE instead would put `Click "Re-draw dot plot"
            // to render the plot.` beside an empty preview on the discard
            // frame, offering to re-draw a figure that is no longer there: the
            // texture is still `Some` until the dispatch block below clears it.
            let live_filters = EnrichDisplayFilters {
                fdr_threshold: new_fdr_threshold,
                min_hit_count: new_min_hit_count,
                top_n: new_top_n,
            };
            texture_is_live = self::texture_is_live(snap.drawn_filters, live_filters);

            ui.add_space(8.0);
            let button_label = draw_button_label(texture_is_live);
            // Always drawable. The empty case is not a UI state: the renderer
            // draws its own "No <entry>s passed …" placeholder INTO the image, which
            // is what the PNG export carries too. Gating the button on a live
            // emptiness check made the screen react to one control change (the
            // one that reaches zero) and ignore every other, and trapped the user
            // in a state whose only advertised exit was a full re-run.
            let can_draw = !busy;
            ui.horizontal(|ui| {
                if crate::ui::widgets::primary_button(ui, button_label, can_draw).clicked() {
                    action = Action::Redraw;
                }
                if snap.rendering {
                    ui.spinner();
                    ui.label("Rendering…");
                }
            });

            ui.add_space(6.0);
            ui.label(RichText::new("Dot plot").strong().color(theme::HEADING));
            if texture_is_live
                && let AppState::Stage3EnrichResult {
                    dotplot: Some(plot),
                    ..
                } = &app.state
            {
                let size = plot.tex.size_vec2();
                ui.add(egui::Image::new(&plot.tex).fit_to_exact_size(size));
            } else if !snap.rendering {
                // The third state renders NEITHER: while a render is in flight
                // and no figure is shown, the area is empty and the progress
                // lives beside the draw button. A prompt here would invite a
                // click on a button that render has disabled.
                ui.label(
                    RichText::new(not_yet_drawn_prompt(button_label))
                        .small()
                        .color(theme::TEXT),
                );
            }

            ui.add_space(8.0);
            if ui
                .add_enabled(!busy, egui::Button::new("Download dot plot PNG"))
                .clicked()
            {
                action = Action::DownloadPng;
            }
            // Two downloads, because one cannot mean both. This screen had a
            // single button that meant "the filtered rows" and then meant
            // "everything", each time silently. Same pair, and the same names,
            // as Stage 2.
            if ui
                .add_enabled(
                    !busy,
                    egui::Button::new("Download enrichment results (CSV)"),
                )
                .clicked()
            {
                action = Action::DownloadFigureCsv;
            }
            if ui
                .add_enabled(!busy, egui::Button::new("Download all results (CSV)"))
                .clicked()
            {
                action = Action::DownloadAllCsv;
            }
            // Start a new analysis — discards this analysis and returns to
            // Stage 1 after a loss-warning confirmation. `!busy`-gated so it
            // can't fire mid refresh/render (keeps `rendering` unreachable
            // while the confirm modal is open).
            if ui
                .add_enabled(!busy, egui::Button::new("Start a new analysis"))
                .clicked()
            {
                action = Action::RequestNewRound;
            }
            // Refresh-in-flight progress now renders inline under each cache
            // line in the data-sources panel above (not here).
        });

    // Confirmation modal.
    if matches!(
        snap.refresh_state_kind,
        RefreshKind::ConfirmingPubchem | RefreshKind::ConfirmingKegg
    ) {
        let is_pubchem = matches!(snap.refresh_state_kind, RefreshKind::ConfirmingPubchem);
        let (title, body, n) = if is_pubchem {
            (
                "Refresh PubChem cache?",
                "This will re-fetch every InChIKey from PubChem (~3-5 min on 13k features).",
                snap.dam_features_total,
            )
        } else {
            (
                "Refresh KEGG conv cache?",
                "This will re-fetch every CID from KEGG /conv/compound/pubchem (~1 min on 10k CIDs).",
                snap.mapped_features,
            )
        };
        egui::Window::new(title)
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ui.ctx(), |ui| {
                ui.label(format!("{body} ({n} entries scoped for this run)"));
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    if ui.button("Continue").clicked() {
                        action = if is_pubchem {
                            Action::ConfirmRefreshPubchem
                        } else {
                            Action::ConfirmRefreshKegg
                        };
                    }
                    if ui.button("Cancel").clicked() {
                        action = Action::CancelRefresh;
                    }
                });
            });
    }

    // New-analysis confirmation modal. Variant-internal (mirrors the refresh
    // confirm modal above), gated on `refresh_state == Idle` so it can never
    // co-show with the refresh modal; the trigger button is `!busy`-gated, so
    // `rendering == true` is unreachable while this modal is open.
    if snap.confirming_new_round && matches!(snap.refresh_state_kind, RefreshKind::Idle) {
        egui::Window::new(RichText::new("Start a new analysis?").color(theme::HEADING))
            .collapsible(false)
            .resizable(false)
            .anchor(egui::Align2::CENTER_CENTER, [0.0, 0.0])
            .show(ui.ctx(), |ui| {
                ui.label(
                    RichText::new(
                        "This clears the loaded files and resets all parameters to defaults.",
                    )
                    .color(theme::TEXT),
                );
                ui.label(
                    RichText::new(
                        "The current DAM / enrichment results and any un-downloaded plots or CSV \
                         will be lost. This cannot be undone.",
                    )
                    .color(theme::WARNING),
                );
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    if ui.button("Start over").clicked() {
                        action = Action::ConfirmNewRound;
                    }
                    if ui.button("Cancel").clicked() {
                        action = Action::CancelNewRound;
                    }
                });
            });
    }

    // Write back export-size + top_n changes. Top N is a dot-plot
    // display cap (see plot/dotplot.rs); changing it on the result page
    // does not require a re-run of the orchestrator — only a fresh
    // `Draw dot plot`, which the user triggers explicitly.
    app.settings.stage3_export_width_in = new_w_in;
    app.settings.stage3_export_height_in = new_h_in;
    app.settings.stage3_export_dpi = new_dpi;
    app.settings.top_n = new_top_n;
    // The relocated significance threshold + Minimum hit count are written back
    // like Top N. All three are LIVE display filters on the FIGURE: the Data
    // tab's significant-count reads them immediately, moving any of them
    // discards the texture, and the next draw renders from them.
    // `Download enrichment results (CSV)` reads the first two;
    // `Download all results (CSV)` reads none of them, and `Top N` reaches
    // neither file. No re-run — none can move `m` or any p-value.
    app.settings.enrichment_fdr_threshold = new_fdr_threshold;
    app.settings.min_hit_count = new_min_hit_count;

    // Invalidation. The comparison was made above the draw button so the label,
    // the blit and the prompt could all agree on the same frame; the texture
    // itself is dropped here, once the borrow of `app.state` the preview held
    // has ended. `DrawnPlot` binds the texture to its filters, so this takes
    // both — there is no way to leave a filter record behind.
    if !texture_is_live
        && let AppState::Stage3EnrichResult { dotplot, .. } = &mut app.state
        && dotplot.take().is_some()
    {
        info!("display filter changed; dot plot texture discarded");
    }

    // A same-frame change to the Height DragValue is a manual override: from
    // here on, redraws honor the user's height instead of re-fitting it to the
    // displayed-row count. A redraw-driven auto-fit writes `settings` in a
    // PRIOR frame, so `snap.export_h_in` already equals `new_h_in` and this
    // stays false for it — only a genuine user drag this frame trips it.
    if !snap.height_user_overridden
        && new_h_in != snap.export_h_in
        && let AppState::Stage3EnrichResult {
            height_user_overridden,
            ..
        } = &mut app.state
    {
        *height_user_overridden = true;
    }

    // Dispatch. (Re-run + the two cache-refresh REQUESTS now arrive via
    // `drain_cache_actions` from the Data-tab Cache block, not the body — those
    // call `rerun` / `request_refresh` directly; the confirm-modal Confirm /
    // Cancel still flow through `action`.)
    match action {
        Action::None => {}
        Action::Redraw => spawn_render(app),
        Action::DownloadPng => download_png(app),
        Action::DownloadFigureCsv => download_csv(app, CsvButton::Figure),
        Action::DownloadAllCsv => download_csv(app, CsvButton::All),
        Action::ConfirmRefreshPubchem => start_refresh(app, true),
        Action::ConfirmRefreshKegg => start_refresh(app, false),
        Action::CancelRefresh => cancel_refresh(app),
        Action::RequestNewRound => set_confirming_new_round(app, true),
        Action::ConfirmNewRound => app.start_new_round(),
        Action::CancelNewRound => set_confirming_new_round(app, false),
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum RefreshKind {
    Idle,
    ConfirmingPubchem,
    ConfirmingKegg,
    RefreshingPubchem,
    RefreshingKegg,
}

fn refresh_state_kind(s: &RefreshState) -> RefreshKind {
    match s {
        RefreshState::Idle => RefreshKind::Idle,
        RefreshState::ConfirmingPubchem => RefreshKind::ConfirmingPubchem,
        RefreshState::ConfirmingKegg => RefreshKind::ConfirmingKegg,
        RefreshState::RefreshingPubchem { .. } => RefreshKind::RefreshingPubchem,
        RefreshState::RefreshingKegg { .. } => RefreshKind::RefreshingKegg,
    }
}

/// Helper: extract `module_retention` from the current Stage3EnrichResult
/// state without holding a borrow across the snapshot construction.
struct ResultSnap {
    mode: AnalysisMode,
    /// Total DAM features across modes — scope count in the PubChem refresh
    /// confirmation modal ("(N entries scoped for this run)").
    dam_features_total: usize,
    /// Features mapped to ≥1 cpd — scope count in the KEGG conv refresh modal.
    mapped_features: usize,
    export_w_in: f64,
    export_h_in: f64,
    export_dpi: u32,
    /// Top N input — changing this on the result page re-targets the next
    /// "Draw / Re-draw dot plot" render.
    top_n: usize,
    /// Significance threshold — relocated from Stage 3 setup
    /// (`add-bottom-panel-data-tab`). A live display cutoff read from SETTINGS,
    /// because the user tunes it on this screen (re-draw, not re-run).
    fdr_threshold: f64,
    /// Correction method of the run being displayed, read from
    /// `EnrichmentResult` — NOT from settings. Every significance noun on this
    /// screen derives from it, so a label can never describe an export the run
    /// did not produce. Pairs with `fdr_threshold` above: method from the run,
    /// threshold from the controls.
    fdr_method: crate::dam::fdr::FdrMethod,
    /// Minimum hit count — relocated from Stage 3 setup. A live display filter
    /// on hit count, re-applied on each draw. It reaches ONE of the two CSV
    /// downloads — the filtered one, whose definition is the figure's row set —
    /// and not the other; and it is consulted after the FDR correction, so it
    /// moves no statistic.
    min_hit_count: usize,
    refresh_state_kind: RefreshKind,
    /// Variant-internal "Start a new analysis?" flag, mirrored here so the
    /// modal-render block (which only has `snap`) can decide whether to show.
    confirming_new_round: bool,
    rendering: bool,
    /// The display filters the currently-held texture was rendered FROM, or
    /// `None` when no texture is held. Compared against the live values to
    /// decide whether the figure still describes what the controls say — see
    /// `show_inner`, which makes that comparison ONCE.
    drawn_filters: Option<EnrichDisplayFilters>,
    /// Whether the user has hand-edited the Height field this result-state.
    /// Mirrored here so the writeback block can detect a fresh height drag
    /// (`new_h_in != export_h_in`) and latch the variant flag.
    height_user_overridden: bool,
}

/// Count rows the dot plot will actually draw before the `top_n` truncation —
/// the SAME predicate as `plot::dotplot`
/// (`r.hits >= min_hit_count && r.fdr < threshold`) and the result-screen
/// significant-entry count. Both filters are live: they are read from settings
/// at draw time, so the canvas height tracks the rows actually on screen.
fn displayed_row_count(
    rows: &[crate::enrichment::types::EnrichmentRow],
    fdr_threshold: f64,
    min_hit_count: usize,
) -> usize {
    rows.iter()
        .filter(|r| r.hits >= min_hit_count && r.fdr < fdr_threshold)
        .count()
}

/// Export height (inches) for the next dot-plot render. While the user has not
/// hand-edited the Height field (`overridden == false`), re-fit it to the live
/// displayed-row count so changing any display filter (the significance
/// threshold, min hit count, or Top N) and re-drawing grows or
/// shrinks the canvas to match — each of those changes now also discarding the
/// previous figure, so the two can never be seen side by side — fixing the
/// stale-autosize squish where rows were crammed into a height computed for the
/// previous run's count. Once the user overrides Height, honor it verbatim.
fn effective_dotplot_height_in(
    overridden: bool,
    current_height_in: f64,
    top_n: usize,
    fdr_threshold: f64,
    min_hit_count: usize,
    rows: &[crate::enrichment::types::EnrichmentRow],
) -> f64 {
    if overridden {
        current_height_in
    } else {
        crate::app::stage3_autosize_height_in(
            top_n,
            displayed_row_count(rows, fdr_threshold, min_hit_count),
        )
    }
}

fn spawn_render(app: &mut App) {
    let export_width_in = app.settings.stage3_export_width_in;
    let export_dpi = app.settings.stage3_export_dpi;
    // The filters this render is being launched with. They travel with the
    // finished buffer so `drain_render` can drop a render that completes after
    // the user has moved past it — the filter controls stay enabled during a
    // render, so that is reachable, and when it happens there is no texture on
    // screen for the ordinary invalidation to clear.
    let filters = app.settings.enrichment_display_filters();
    let enrichment_fdr_threshold = filters.fdr_threshold;
    let min_hit_count = filters.min_hit_count;
    let top_n = filters.top_n;
    let AppState::Stage3EnrichResult {
        enrichment_result,
        rendering,
        render_rx,
        height_user_overridden,
        ..
    } = &mut app.state
    else {
        return;
    };
    // Method from the RUN, not from settings — the plot's chrome must name
    // what produced the values it draws (`enrichment-dot-plot` requires the
    // same `FdrMethod` that fed the `run_ora` call which produced this result).
    let fdr_method = enrichment_result.fdr_method;
    // Re-fit the export height to the rows this draw will actually show (unless
    // the user has hand-set Height), then persist it so the Height DragValue
    // reflects the fitted value next frame. This is what makes a redraw after
    // loosening the FDR threshold grow the canvas instead of cramming rows into
    // a height autosized for the previous run.
    let export_height_in = effective_dotplot_height_in(
        *height_user_overridden,
        app.settings.stage3_export_height_in,
        top_n,
        enrichment_fdr_threshold,
        min_hit_count,
        &enrichment_result.rows,
    );
    app.settings.stage3_export_height_in = export_height_in;
    let (w_px, h_px) =
        crate::ui::widgets::export_pixels(export_width_in, export_height_in, export_dpi);
    let opts = DotplotOpts {
        width_px: w_px,
        height_px: h_px,
        fdr_threshold: enrichment_fdr_threshold,
        min_hit_count,
        top_n,
        fdr_method,
        entry_label: app.settings.analysis_mode.entry_label_singular(),
    };
    let result_clone = enrichment_result.clone();
    let (tx, rx) =
        mpsc::channel::<Result<crate::app::DotplotRenderOf<EnrichDisplayFilters>, String>>();
    *render_rx = Some(rx);
    *rendering = true;
    info!(width_px = w_px, height_px = h_px, "rendering dot plot");
    app.rt.spawn(async move {
        let r = tokio::task::spawn_blocking(move || render_dotplot(&result_clone, &opts))
            .await
            .map_err(|e| e.to_string())
            .and_then(|res| res.map_err(|e| e.to_string()))
            .map(|buf| ((buf, w_px, h_px), filters));
        let _ = tx.send(r);
    });
}

fn drain_render(app: &mut App, ctx: &egui::Context) {
    let r = {
        let AppState::Stage3EnrichResult {
            rendering,
            render_rx,
            ..
        } = &mut app.state
        else {
            return;
        };
        if !*rendering {
            return;
        }
        let Some(rx) = render_rx else { return };
        let Ok(msg) = rx.try_recv() else { return };
        *rendering = false;
        *render_rx = None;
        match msg {
            Ok(triple) => Some(triple),
            Err(e) => {
                error!(error = %e, "dot plot render failed");
                None
            }
        }
    };
    let Some(((buf, w, h), filters)) = r else {
        return;
    };
    // A render that finished after a filter moved describes values the user has
    // already left. Drop it rather than install it — the same comparison the
    // frame body makes, at the one point where there is no texture on screen for
    // that comparison to reach.
    if filters != app.settings.enrichment_display_filters() {
        info!("display filter changed during render; finished dot plot discarded");
        return;
    }
    let img = egui::ColorImage::from_rgba_unmultiplied([w as usize, h as usize], &buf);
    let tex = ctx.load_texture("dotplot", img, egui::TextureOptions::LINEAR);
    if let AppState::Stage3EnrichResult { dotplot, .. } = &mut app.state {
        *dotplot = Some(DrawnPlot { tex, filters });
        info!(width_px = w, height_px = h, "dot plot texture uploaded");
    }
}

fn drain_refresh(app: &mut App) {
    // Drain progress events from the active refresh into completed/total.
    let mut terminal: Option<Result<Stage3RunOutput, String>> = None;

    if let AppState::Stage3EnrichResult { refresh_state, .. } = &mut app.state {
        match refresh_state {
            RefreshState::RefreshingPubchem {
                progress_rx,
                result_rx,
                completed,
                total,
                ..
            } => {
                while let Ok(p) = progress_rx.try_recv() {
                    *completed = p.from_cache + p.fetched;
                    *total = p.total_inputs.max(*total);
                }
                if let Ok(t) = result_rx.try_recv() {
                    terminal = Some(t);
                }
            }
            RefreshState::RefreshingKegg {
                progress_rx,
                result_rx,
                completed,
                total,
                ..
            } => {
                while let Ok(p) = progress_rx.try_recv() {
                    *completed = p.from_cache + p.fetched;
                    *total = p.total_inputs.max(*total);
                }
                if let Ok(t) = result_rx.try_recv() {
                    terminal = Some(t);
                }
            }
            _ => {}
        }
    }

    if let Some(t) = terminal
        && let AppState::Stage3EnrichResult {
            refresh_state,
            pubchem_time_span,
            kegg_conv_time_span,
            ..
        } = &mut app.state
    {
        match t {
            Ok(out) => {
                info!("refresh complete; updating time-span fields");
                *pubchem_time_span = out.pubchem_time_span;
                *kegg_conv_time_span = out.kegg_conv_time_span;
            }
            Err(e) => {
                error!(error = %e, "refresh failed");
            }
        }
        *refresh_state = RefreshState::Idle;
    }
}

fn download_png(app: &App) {
    let export_width_in = app.settings.stage3_export_width_in;
    let export_dpi = app.settings.stage3_export_dpi;
    let enrichment_fdr_threshold = app.settings.enrichment_fdr_threshold;
    let min_hit_count = app.settings.min_hit_count;
    let top_n = app.settings.top_n;
    let AppState::Stage3EnrichResult {
        enrichment_result,
        height_user_overridden,
        ..
    } = &app.state
    else {
        return;
    };
    // Method from the RUN, not from settings — the plot's chrome must name
    // what produced the values it draws (`enrichment-dot-plot` requires the
    // same `FdrMethod` that fed the `run_ora` call which produced this result).
    let fdr_method = enrichment_result.fdr_method;
    // Size the export the same way the preview does. A filter change discards
    // the preview, so this button is the one route to a figure at the current
    // values without pressing Draw first — and it must produce the figure the
    // preview WOULD show, not the one it last did. `&App` is immutable here, so
    // the fitted value is used locally and not persisted — the next draw writes
    // it back to settings.
    let export_height_in = effective_dotplot_height_in(
        *height_user_overridden,
        app.settings.stage3_export_height_in,
        top_n,
        enrichment_fdr_threshold,
        min_hit_count,
        &enrichment_result.rows,
    );
    let Some(path) = crate::ui::widgets::save_dialog("PNG", "png", "enrichment-dotplot.png") else {
        return;
    };
    let (w_px, h_px) =
        crate::ui::widgets::export_pixels(export_width_in, export_height_in, export_dpi);
    let opts = DotplotOpts {
        width_px: w_px,
        height_px: h_px,
        fdr_threshold: enrichment_fdr_threshold,
        min_hit_count,
        top_n,
        fdr_method,
        entry_label: app.settings.analysis_mode.entry_label_singular(),
    };
    if let Err(e) = export_dotplot_png(enrichment_result, &opts, export_dpi, &path) {
        error!(error = %e, "dot plot PNG export failed");
    } else {
        info!(
            path = %path.display(),
            width_px = w_px,
            height_px = h_px,
            dpi = export_dpi,
            "dot plot PNG exported"
        );
    }
}

/// Which of the screen's two CSV buttons was pressed.
///
/// Named rather than a bare `bool` so a call site cannot say "filtered" while
/// meaning something else. Each arm carries its own default filename: nothing
/// INSIDE either file records which button wrote it, so the filename is the only
/// distinction the user carries away.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CsvButton {
    /// `Download enrichment results (CSV)` — the rows the figure is drawn from.
    Figure,
    /// `Download all results (CSV)` — every surviving row.
    All,
}

impl CsvButton {
    /// Deliberately NOT a bare `all_results.csv` for the complete export: Stage 2
    /// already proposes that name for an unrelated file, and the two would
    /// overwrite each other in one download folder.
    fn default_filename(self) -> &'static str {
        match self {
            CsvButton::Figure => "enrichment.csv",
            CsvButton::All => "enrichment_all_results.csv",
        }
    }
}

fn download_csv(app: &App, button: CsvButton) {
    // `Top N` reaches NEITHER file. It is an ordering cap — "how many of the
    // ranked rows fit on the axis" — so a file bounded by it would have a row
    // count meaning "the twenty I happened to be looking at". The other two are
    // per-row tests, which is what makes them meaningful in a file.
    let selection = match button {
        CsvButton::All => RowSelection::All,
        CsvButton::Figure => RowSelection::Figure {
            fdr_threshold: app.settings.enrichment_fdr_threshold,
            min_hit_count: app.settings.min_hit_count,
        },
    };
    let AppState::Stage3EnrichResult {
        enrichment_result,
        dual_mode_breakdown,
        module_retention,
        ..
    } = &app.state
    else {
        return;
    };
    let Some(path) = crate::ui::widgets::save_dialog("CSV", "csv", button.default_filename())
    else {
        return;
    };
    let mut file = match std::fs::File::create(&path) {
        Ok(f) => f,
        Err(e) => {
            error!(path = %path.display(), error = %e, "could not create CSV file");
            return;
        }
    };
    let is_dual = dual_mode_breakdown.is_some();
    // Module-mode runs self-document the Group-overlap threshold used; Pathway
    // runs (`module_retention == None`) pass `None` so the line is omitted.
    // BOTH buttons supply it: the two files' comment blocks are required to be
    // identical, so a site that omitted it would break that rather than differ.
    let min_group_overlap = module_retention.as_ref().map(|r| r.min_group_overlap);
    // Which rows are written is a property of the BUTTON and of nothing else.
    if let Err(e) = export_csv_with_mode(
        &mut file,
        enrichment_result,
        is_dual,
        min_group_overlap,
        selection,
    ) {
        error!(error = %e, "enrichment CSV export failed");
    } else {
        info!(path = %path.display(), is_dual, ?button, "enrichment CSV exported");
    }
}

/// `pub(crate)` purely so `app.rs`'s test module can drive it — this file's own
/// test module holds pure helpers and no `App` fixture. Same move, and same
/// reason, as `build_stage3_spawn_inputs`.
pub(crate) fn rerun(app: &mut App) {
    // Establish that the run CAN be assembled before consuming the state that
    // holds the previous one. `Stage3EnrichResult` is the sole holder of the
    // enrichment output and there is no route back to a discarded run, so
    // transitioning first and discovering afterwards costs the user everything
    // and returns nothing. `start_refresh` below is the local precedent for the
    // borrow-first shape. See the `data-summary-panel` capability spec.
    {
        let AppState::Stage3EnrichResult { dam_results, .. } = &app.state else {
            return;
        };
        if crate::ui::stage3_setup::build_stage3_run_inputs(dam_results, &app.settings, &app.cache)
            .is_none()
        {
            // Deliberately does NOT say "cache missing": in the state that
            // actually produces this, the cache is present and the target
            // SELECTION is gone (a retired organism clears it).
            warn!(
                "re-run refused: the analysis target for the current mode is incomplete; state left intact"
            );
            return;
        }
    }

    let prev = std::mem::take(&mut app.state);
    let AppState::Stage3EnrichResult { dam_results, .. } = prev else {
        return;
    };
    let Some((params, target, pubchem_total)) =
        crate::ui::stage3_setup::build_stage3_run_inputs(&dam_results, &app.settings, &app.cache)
    else {
        // Unreachable: the guard above established this succeeds. Restoring the
        // result is impossible here (its other fields were consumed by the take),
        // so the honest fallback is the setup screen with an accurate message.
        error!("re-run inputs vanished between the guard and the spawn");
        app.state = AppState::Stage3EnrichSetup {
            dam_results,
            error: Some(
                "The analysis target became unavailable while starting the re-run.".to_string(),
            ),
            kegg_fetch: None,
            modules_fetch: None,
        };
        return;
    };
    app.spawn_stage3_run(crate::app::RunSpawn {
        payload: crate::app::RunPayloadSpec::Enrichment {
            dam_results,
            params,
        },
        target,
        pubchem_total,
    });
}

/// Drain the Data-tab Cache-block request flags (set by `data_tab::show`, which
/// renders in the bottom panel BEFORE this central-panel fn each frame) into the
/// existing confirm-modal / re-run flow. Only honoured while no refresh is
/// already in flight (the Data-tab buttons are disabled when busy; this is
/// defense-in-depth). Re-run transitions the state to `Stage3EnrichRunning`, so
/// the snapshot below returns `None` and the screen yields to the running view.
fn drain_cache_actions(app: &mut App) {
    let catalogue = std::mem::take(&mut app.log_ui.refresh_catalogue_requested);
    let pubchem = std::mem::take(&mut app.log_ui.refresh_pubchem_requested);
    let kegg = std::mem::take(&mut app.log_ui.refresh_kegg_conv_requested);
    let rerun_req = std::mem::take(&mut app.log_ui.rerun_enrichment_requested);
    let idle = matches!(
        &app.state,
        AppState::Stage3EnrichResult {
            refresh_state: RefreshState::Idle,
            ..
        }
    );
    if !idle {
        return;
    }
    if catalogue {
        refresh_catalogue_via_setup(app);
    } else if pubchem {
        request_refresh(app, RefreshKind::ConfirmingPubchem);
    } else if kegg {
        request_refresh(app, RefreshKind::ConfirmingKegg);
    } else if rerun_req {
        rerun(app);
    }
}

/// Result-page catalogue (module / pathway) refresh. The `Stage3EnrichResult`
/// state has no `kegg_fetch` / `modules_fetch` in-flight infrastructure (that
/// lives on `Stage3EnrichSetup`), so navigate back to Stage 3 setup — carrying
/// `dam_results`, preserving the Group/species selection — and trigger the
/// force re-fetch there, where the fetch's progress strip + `Run Enrichment`
/// live. Only reached from `drain_cache_actions` after the `idle` guard, so the
/// state is always `Stage3EnrichResult` here.
pub(crate) fn refresh_catalogue_via_setup(app: &mut App) {
    // Establish that the re-fetch CAN start before consuming the state that
    // holds the completed run. The target selection is what the fetch is keyed
    // by, and a retired organism clears it while leaving the fetched catalogue
    // intact — so this is NOT the same predicate `rerun` uses, which
    // additionally needs that catalogue. See the `data-summary-panel` spec.
    let target_ready = match app.settings.analysis_mode {
        AnalysisMode::Module => {
            app.settings.organism_group_level.is_some()
                && app.settings.organism_group.is_some()
                && app.cache.group_org_codes.is_some()
        }
        AnalysisMode::Pathway => app.settings.kegg_species.is_some(),
    };
    if !target_ready {
        warn!(
            mode = ?app.settings.analysis_mode,
            "catalogue refresh refused: no analysis target selected; state left intact"
        );
        return;
    }

    let prev = std::mem::take(&mut app.state);
    if let AppState::Stage3EnrichResult { dam_results, .. } = prev {
        app.state = AppState::Stage3EnrichSetup {
            dam_results,
            error: None,
            kegg_fetch: None,
            modules_fetch: None,
        };
        match app.settings.analysis_mode {
            AnalysisMode::Module => {
                if let (Some(level), Some(group), Some(org_codes)) = (
                    app.settings.organism_group_level,
                    app.settings.organism_group.clone(),
                    app.cache.group_org_codes.clone(),
                ) {
                    crate::ui::stage3_setup::spawn_modules_fetch(
                        app, level, group, org_codes, true,
                    );
                }
            }
            AnalysisMode::Pathway => crate::ui::stage3_setup::handle_species_refresh(app),
        }
    }
}

fn request_refresh(app: &mut App, kind: RefreshKind) {
    if let AppState::Stage3EnrichResult { refresh_state, .. } = &mut app.state {
        *refresh_state = match kind {
            RefreshKind::ConfirmingPubchem => RefreshState::ConfirmingPubchem,
            RefreshKind::ConfirmingKegg => RefreshState::ConfirmingKegg,
            _ => return,
        };
    }
}

fn cancel_refresh(app: &mut App) {
    if let AppState::Stage3EnrichResult { refresh_state, .. } = &mut app.state {
        *refresh_state = RefreshState::Idle;
    }
}

/// Set the variant-internal `confirming_new_round` flag on the current
/// Stage 3 result state — `true` opens the "Start a new analysis?" confirm
/// modal, `false` closes it (Cancel). No-op off the result state. The actual
/// reset on confirm is `App::start_new_round`, not this helper.
fn set_confirming_new_round(app: &mut App, value: bool) {
    if let AppState::Stage3EnrichResult {
        confirming_new_round,
        ..
    } = &mut app.state
    {
        *confirming_new_round = value;
    }
}

fn start_refresh(app: &mut App, is_pubchem: bool) {
    // Build the orchestrator inputs via the shared seam, then set the
    // force-refresh flags for whichever cache is being refreshed. `start_refresh`
    // keeps its own in-place `RefreshState` transition (it does NOT enter
    // `Stage3EnrichRunning`) + the progress-bridge threads, so it does NOT use
    // `App::spawn_stage3_run`.
    let (dam_results, params, target, pubchem_total) = {
        let AppState::Stage3EnrichResult { dam_results, .. } = &app.state else {
            return;
        };
        let Some((mut params, target, pubchem_total)) =
            crate::ui::stage3_setup::build_stage3_run_inputs(
                dam_results,
                &app.settings,
                &app.cache,
            )
        else {
            error!("start_refresh: AnalysisPayload could not be built; cache missing");
            return;
        };
        params.force_refresh_pubchem = is_pubchem;
        params.force_refresh_kegg_conv = !is_pubchem;
        (dam_results.clone(), params, target, pubchem_total)
    };

    let (pub_tx, pub_rx_for_state) = mpsc::channel();
    let (kegg_tx, kegg_rx_for_state) = mpsc::channel();
    let (pub_tx_run, pub_rx_run) = mpsc::channel();
    let (kegg_tx_run, kegg_rx_run) = mpsc::channel();
    let (result_tx, result_rx) = mpsc::channel::<Result<Stage3RunOutput, String>>();

    // Bridge: forward run-side progress to the state-side channels.
    let _ = std::thread::Builder::new()
        .name("pubchem-progress-bridge".into())
        .spawn(move || {
            while let Ok(ev) = pub_rx_run.recv() {
                let _ = pub_tx.send(ev);
            }
        });
    let _ = std::thread::Builder::new()
        .name("kegg-progress-bridge".into())
        .spawn(move || {
            while let Ok(ev) = kegg_rx_run.recv() {
                let _ = kegg_tx.send(ev);
            }
        });

    let kegg_client = app.kegg.clone();
    let run_handle = app
        .rt
        .spawn(async move {
            let pubchem = PubchemClient::new();
            let r = run_stage3(
                &pubchem,
                &kegg_client,
                &dam_results,
                &target,
                params,
                pub_tx_run,
                kegg_tx_run,
            )
            .await;
            let _ = result_tx.send(r.map_err(|e| e.to_string()));
        })
        .abort_handle();

    if let AppState::Stage3EnrichResult { refresh_state, .. } = &mut app.state {
        *refresh_state = if is_pubchem {
            RefreshState::RefreshingPubchem {
                progress_rx: pub_rx_for_state,
                result_rx,
                completed: 0,
                total: pubchem_total,
                run_handle,
            }
        } else {
            RefreshState::RefreshingKegg {
                progress_rx: kegg_rx_for_state,
                result_rx,
                completed: 0,
                total: 0,
                run_handle,
            }
        };
    }
}

#[cfg(test)]
mod tests {
    use crate::dam::fdr::FdrMethod;

    /// A minimal `Stage3EnrichResult` carrying a chosen correction method.
    ///
    /// Duplicated rather than shared with `app.rs`'s `stage3_result_state()`,
    /// which is private to that module's test scope.
    fn result_state_with(method: FdrMethod) -> AppState {
        AppState::Stage3EnrichResult {
            dam_results: vec![],
            module_retention: None,
            enrichment_result: crate::enrichment::types::EnrichmentResult {
                universe_size: 0,
                dam_cpd_size: 0,
                direction: crate::enrichment::types::EnrichmentDirection::Both,
                min_hit_count: 1,
                min_entry_size: 1,
                entries_dropped_by_min_entry_size: 0,
                empty_compound_count: 0,
                rows: vec![],
                fdr_method: method,
            },
            mapped_universe: std::collections::HashSet::new(),
            feature_to_cpds: std::collections::HashMap::new(),
            pubchem_time_span: None,
            kegg_conv_time_span: None,
            dual_mode_breakdown: None,
            funnel: crate::app::Stage3Funnel::default(),
            dotplot: None,
            rendering: false,
            render_rx: None,
            refresh_state: RefreshState::Idle,
            confirming_new_round: false,
            height_user_overridden: false,
        }
    }

    /// **Method from the run, threshold from the controls.**
    ///
    /// The snapshot must take its correction method from `EnrichmentResult` —
    /// the method that produced the numbers on display — and its threshold from
    /// settings, which is the live control the user tunes on this screen. This
    /// is the half of the sourcing rule a test can see; that each label then
    /// passes `snap.fdr_method` rather than the `app.settings` value also in
    /// scope inside the render closure is not mechanically checkable here.
    #[test]
    fn snapshot_takes_method_from_the_run_and_threshold_from_settings() {
        let state = result_state_with(FdrMethod::BenjaminiYekutieli);
        // Settings deliberately disagree with the run.
        let settings = crate::app::SessionSettings {
            enrichment_fdr_method: FdrMethod::NoCorrection,
            enrichment_fdr_threshold: 0.02,
            ..Default::default()
        };

        let snap = result_snap(&state, &settings).expect("result state yields a snapshot");

        assert_eq!(snap.fdr_method, FdrMethod::BenjaminiYekutieli);
        assert_eq!(snap.fdr_threshold, 0.02);
    }

    #[test]
    fn result_snap_is_none_off_the_result_screen() {
        let settings = crate::app::SessionSettings::default();
        assert!(result_snap(&AppState::default(), &settings).is_none());
    }

    /// The three label builders name the QUANTITY, so they follow the method.
    #[test]
    fn label_builders_follow_the_method() {
        assert_eq!(
            threshold_label(FdrMethod::BenjaminiYekutieli),
            "Enrichment FDR threshold:"
        );
        assert_eq!(
            threshold_label(FdrMethod::NoCorrection),
            "Enrichment p-value threshold:"
        );

        // The `after FDR` clause is DROPPED under NoCorrection, not re-pointed:
        // there is no FDR stage for the filter to follow.
        assert!(min_hit_tooltip(FdrMethod::BenjaminiHochberg).contains("after FDR"));
        assert!(!min_hit_tooltip(FdrMethod::NoCorrection).contains("FDR"));
        assert!(!min_hit_tooltip(FdrMethod::NoCorrection).contains("after"));
    }

    use super::*;
    use crate::enrichment::types::EnrichmentRow;

    /// `has_hit == true` gives the row one hit, so it clears `min_hit_count = 1`
    /// and is filtered out by `min_hit_count = 2`.
    fn row(fdr: f64, has_hit: bool) -> EnrichmentRow {
        EnrichmentRow {
            entry_id: "p".to_string(),
            entry_name: "n".to_string(),
            hits: usize::from(has_hit),
            total: 5,
            expected: 1.0,
            enrichment_ratio: 1.0,
            p_value: fdr,
            fdr,
            hit_kegg_ids: vec![],
        }
    }

    #[test]
    fn displayed_row_count_matches_plot_predicate() {
        let rows = vec![
            row(0.01, true),  // significant + has a hit → counts
            row(0.01, false), // no hits → below any min_hit_count → excluded
            row(0.5, true),   // not significant at 0.05 → excluded
        ];
        assert_eq!(displayed_row_count(&rows, 0.05, 1), 1);
        // Loosen the threshold: the 0.5 row now passes too.
        assert_eq!(displayed_row_count(&rows, 1.0, 1), 2);
        // Raise min_hit_count instead — LIVE, no re-run: both one-hit rows drop.
        assert_eq!(displayed_row_count(&rows, 1.0, 2), 0);
    }

    #[test]
    fn effective_height_refits_when_filter_change_reveals_rows() {
        // 20 rows significant only once the FDR threshold is loosened past 0.5
        // — the exact Plot 1 vs Plot 2 scenario (FDR 0.05 → 0 rows; loosened →
        // 20 rows). Auto-fit is on (`overridden == false`).
        let rows: Vec<EnrichmentRow> = (0..20).map(|_| row(0.5, true)).collect();

        // At FDR 0.05: 0 displayed → autosize floors to 2.0 in (Plot 1's squish).
        let h_strict = effective_dotplot_height_in(false, 7.0, 20, 0.05, 1, &rows);
        assert!((h_strict - 2.0).abs() < 1e-9, "got {h_strict}");

        // Loosened to 1.0: 20 displayed (capped by top_n=20) → grows to 7.0 in
        // (Plot 2's correctly-sized, label-wrapping canvas) WITHOUT a re-run.
        let h_loose = effective_dotplot_height_in(false, 2.0, 20, 1.0, 1, &rows);
        assert!((h_loose - 7.0).abs() < 1e-9, "got {h_loose}");
    }

    #[test]
    fn effective_height_honors_user_override_verbatim() {
        let rows: Vec<EnrichmentRow> = (0..20).map(|_| row(0.5, true)).collect();
        // overridden == true → the hand-set height is returned regardless of
        // how many rows the loosened threshold would display.
        let h = effective_dotplot_height_in(true, 12.5, 20, 1.0, 1, &rows);
        assert!((h - 12.5).abs() < 1e-9, "got {h}");
    }
}

#[cfg(test)]
mod invalidation_tests {
    use super::*;

    fn f(fdr_threshold: f64, min_hit_count: usize, top_n: usize) -> EnrichDisplayFilters {
        EnrichDisplayFilters {
            fdr_threshold,
            min_hit_count,
            top_n,
        }
    }

    #[test]
    fn a_texture_is_live_only_while_every_filter_still_matches() {
        let drawn = f(0.05, 1, 20);
        assert!(texture_is_live(Some(drawn), drawn));
        // Each of the three, alone, is enough to kill it.
        assert!(!texture_is_live(Some(drawn), f(0.01, 1, 20)));
        assert!(!texture_is_live(Some(drawn), f(0.05, 3, 20)));
        assert!(!texture_is_live(Some(drawn), f(0.05, 1, 30)));
        // No texture is never live.
        assert!(!texture_is_live(None, drawn));
    }

    #[test]
    fn changing_a_filter_back_still_costs_a_re_render() {
        // The user raises `min_hit_count` and lowers it again. The predicate is
        // true again on the values alone — but by then the FIRST change has
        // already taken the texture to `None` (`show_inner`'s discard block), so
        // what the screen actually asks is `texture_is_live(None, ...)`.
        //
        // This is D1's accepted cost, pinned rather than assumed: a label would
        // have cleared itself here, a discard cannot un-discard.
        let drawn = f(0.05, 1, 20);
        assert!(!texture_is_live(Some(drawn), f(0.05, 3, 20)));
        assert!(!texture_is_live(None, drawn));
    }

    #[test]
    fn export_size_changes_do_not_invalidate() {
        // Width / Height / DPI are not display filters: they change the canvas,
        // not which rows are on it. They are absent from the tuple by design,
        // and `Re-draw dot plot` is reachable only through them.
        let drawn = f(0.05, 1, 20);
        assert!(texture_is_live(Some(drawn), drawn));
        assert_eq!(draw_button_label(true), "Re-draw dot plot");
    }

    #[test]
    fn the_discard_frame_never_offers_to_re_draw_what_is_gone() {
        // The one-frame artifact the same-frame design creates: on the discard
        // frame the texture is still `Some`, so anything reading raw presence
        // would say "Re-draw". Label and prompt both read the predicate instead.
        let drawn = f(0.05, 1, 20);
        let live = f(0.05, 4, 20);
        let label = draw_button_label(texture_is_live(Some(drawn), live));
        assert_eq!(label, "Draw dot plot");
        assert_eq!(
            not_yet_drawn_prompt(label),
            "Click \"Draw dot plot\" to render the plot."
        );
    }

    #[test]
    fn a_late_render_is_rejected_by_the_same_comparison() {
        // `drain_render` compares the filters the finished render carries
        // against the live ones. Same operator, at the one point where there is
        // no texture on screen for the frame body's comparison to reach.
        let launched_with = f(0.05, 1, 20);
        let live_now = f(0.05, 4, 20);
        assert_ne!(launched_with, live_now, "the render is stale on arrival");
        assert_eq!(launched_with, launched_with, "an unchanged render installs");
    }

    #[test]
    fn the_two_csv_buttons_propose_different_default_filenames() {
        // Nothing INSIDE either file records which button wrote it — the comment
        // block is a closed ordered contract this change adds no line to — so
        // the filename is the only distinction that survives the download.
        assert_ne!(
            CsvButton::Figure.default_filename(),
            CsvButton::All.default_filename()
        );
        assert_eq!(CsvButton::Figure.default_filename(), "enrichment.csv");
        // NOT a bare `all_results.csv`: Stage 2 already proposes that name for
        // an unrelated file, and the two would collide in one download folder.
        assert_eq!(
            CsvButton::All.default_filename(),
            "enrichment_all_results.csv"
        );
    }
}
