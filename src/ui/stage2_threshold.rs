use egui::RichText;
use std::sync::mpsc;
use tracing::{error, info};

use crate::app::{App, AppState, VolcanoRender};
use crate::dam::DamMethod;
use crate::data::IonMode;
use crate::plot::{VolcanoOpts, export_volcano_png, render_volcano};
use crate::theme;

#[derive(Debug, Clone, Copy)]
enum Action {
    None,
    Redraw,
    DownloadPng,
    DownloadDamCsv,
    DownloadAllCsv,
    ContinueToEnrichment,
    SwitchTab(IonMode),
}

pub fn show(ui: &mut egui::Ui, app: &mut App) {
    // First, drain any in-flight volcano render.
    drain_render(app, ui.ctx());

    // Snapshot the threshold-state fields we need to render UI, then release the borrow.
    let snapshot = match &app.state {
        AppState::Stage2DamThreshold {
            dam_results,
            rendering,
            active_volcano_tab,
            volcano_textures,
            ..
        } => {
            let ion_tables = app.inputs.ion_tables.as_slice();
            let active_idx = active_idx_for(ion_tables, *active_volcano_tab);
            let dam = &dam_results[active_idx];
            // Per-mode DAM up/down/ns tallies and the dedup summary moved to the
            // bottom-panel Data tab (`data-summary-panel`); the threshold body no
            // longer renders them.
            Some(Snapshot {
                method: dam.method,
                numerator: dam.numerator.clone(),
                denominator: dam.denominator.clone(),
                fc_threshold: app.settings.fc_threshold,
                fdr_threshold: app.settings.fdr_threshold,
                delta_threshold: app.settings.delta_threshold,
                export_width_in: app.settings.stage2_export_width_in,
                export_height_in: app.settings.stage2_export_height_in,
                export_dpi: app.settings.stage2_export_dpi,
                rendering: *rendering,
                has_texture: volcano_textures
                    .get(active_idx)
                    .and_then(|t| t.as_ref())
                    .is_some(),
                tabs: ion_tables.iter().map(|it| it.mode).collect(),
                active_tab: *active_volcano_tab,
            })
        }
        _ => None,
    };
    let Some(snap) = snapshot else {
        return;
    };

    // Build UI; collect any pending action.
    let mut action = Action::None;
    let mut new_fc = snap.fc_threshold;
    let mut new_fdr = snap.fdr_threshold;
    let mut new_delta = snap.delta_threshold;
    let mut new_w_in = snap.export_width_in;
    let mut new_h_in = snap.export_height_in;
    let mut new_dpi = snap.export_dpi;

    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.heading(
                egui::RichText::new(format!(
                    "Stage 2 — DAM Result: {} vs {}",
                    snap.numerator, snap.denominator
                ))
                .color(theme::HEADING),
            );
            ui.add_space(6.0);
            // Back-navigation to DAM Setup is handled by the global stage
            // stepper (`ui::stepper`); no per-screen Back button.
            ui.add_space(8.0);
            // The `Method: <method>` line moved to the Data tab's `DAM data`
            // block (`apply-ui-design-md-tweaks`); the body opens into the
            // threshold controls.

            // Dedup summary + the `Download dedup audit (CSV)` button moved to
            // the bottom-panel Data tab (`data-summary-panel`).

            let bm = snap.method == DamMethod::BrunnerMunzel;
            ui.horizontal(|ui| {
                ui.label("Fold change threshold:");
                ui.add(
                    egui::DragValue::new(&mut new_fc)
                        .speed(0.1)
                        .range(1.0..=1024.0),
                );
                ui.label(
                    RichText::new(format!(
                        "(uses |log2(FC)| ≥ {:.3} internally)",
                        new_fc.log2()
                    ))
                    .small()
                    .color(theme::TEXT),
                );
            });
            ui.horizontal(|ui| {
                ui.label("FDR threshold:");
                ui.add(
                    egui::DragValue::new(&mut new_fdr)
                        .speed(0.001)
                        .range(0.0001..=1.0),
                );
            });
            if bm {
                ui.horizontal(|ui| {
                    ui.label("δ threshold (Cliff's delta):");
                    ui.add(
                        egui::DragValue::new(&mut new_delta)
                            .speed(0.01)
                            .range(0.0..=1.0),
                    );
                });
            }

            // Per-mode DAM up/down/ns tally moved to the Data tab.

            ui.add_space(12.0);
            crate::ui::widgets::png_export_size_controls(
                ui,
                &mut new_w_in,
                &mut new_h_in,
                &mut new_dpi,
            );
            // ── Volcano area ──
            ui.add_space(8.0);
            // Tab bar in dual-mode (hidden in single-mode).
            if snap.tabs.len() >= 2 {
                // §4 Primary "Segmented" control — POS/NEG is a major view-mode
                // switch (a true 2-way toggle, the segmented control's home turf).
                crate::ui::widgets::segmented_track(ui, |ui| {
                    for mode in &snap.tabs {
                        let selected = *mode == snap.active_tab;
                        if crate::ui::widgets::segmented_tab(
                            ui,
                            selected,
                            &format!("{mode}"),
                            true,
                            true,
                        )
                        .clicked()
                            && !selected
                        {
                            action = Action::SwitchTab(*mode);
                        }
                    }
                });
                ui.add_space(4.0);
            }

            let button_label = if snap.has_texture {
                "Re-draw volcano"
            } else {
                "Draw volcano"
            };
            ui.horizontal(|ui| {
                if crate::ui::widgets::primary_button(ui, button_label, !snap.rendering).clicked() {
                    action = Action::Redraw;
                }
                if snap.rendering {
                    ui.spinner();
                    ui.label("Rendering…");
                }
            });

            // Volcano preview.
            ui.add_space(4.0);
            ui.label(
                RichText::new("Preview figure")
                    .strong()
                    .color(theme::HEADING),
            );
            if snap.has_texture
                && let AppState::Stage2DamThreshold {
                    volcano_textures,
                    active_volcano_tab,
                    ..
                } = &app.state
                && let Some(tex) = volcano_textures
                    .get(active_idx_for(
                        app.inputs.ion_tables.as_slice(),
                        *active_volcano_tab,
                    ))
                    .and_then(|t| t.as_ref())
            {
                let size = tex.size_vec2();
                ui.add(egui::Image::new(tex).fit_to_exact_size(size));
            } else if !snap.rendering {
                ui.label(
                    RichText::new(format!("Click \"{button_label}\" to render the plot."))
                        .small()
                        .color(theme::TEXT),
                );
            }

            ui.add_space(8.0);
            if ui.button("Download volcano PNG").clicked() {
                action = Action::DownloadPng;
            }
            ui.horizontal(|ui| {
                if ui.button("Download DAM (CSV)").clicked() {
                    action = Action::DownloadDamCsv;
                }
                if ui.button("Download all results (CSV)").clicked() {
                    action = Action::DownloadAllCsv;
                }
            });
            if crate::ui::widgets::primary_button(ui, "Continue to Enrichment", true).clicked() {
                action = Action::ContinueToEnrichment;
            }
        });

    // Write back any threshold / export changes the user dragged. Threshold
    // changes invalidate the volcano texture cache.
    let thresholds_changed = (app.settings.fc_threshold != new_fc)
        || (app.settings.fdr_threshold != new_fdr)
        || (app.settings.delta_threshold != new_delta);
    app.settings.fc_threshold = new_fc;
    app.settings.fdr_threshold = new_fdr;
    app.settings.delta_threshold = new_delta;
    app.settings.stage2_export_width_in = new_w_in;
    app.settings.stage2_export_height_in = new_h_in;
    app.settings.stage2_export_dpi = new_dpi;
    if thresholds_changed
        && let AppState::Stage2DamThreshold {
            volcano_textures, ..
        } = &mut app.state
    {
        for slot in volcano_textures.iter_mut() {
            *slot = None;
        }
    }

    // Dispatch the action.
    match action {
        Action::None => {}
        Action::Redraw => spawn_render(app),
        Action::DownloadPng => download_png(app),
        Action::DownloadDamCsv => download_csv(app, true),
        Action::DownloadAllCsv => download_csv(app, false),
        Action::ContinueToEnrichment => continue_to_enrichment(app),
        Action::SwitchTab(mode) => {
            if let AppState::Stage2DamThreshold {
                active_volcano_tab, ..
            } = &mut app.state
            {
                *active_volcano_tab = mode;
            }
        }
    }
}

/// Resolve the active tab to an index into `app.inputs.ion_tables` /
/// `dam_results`. Falls back to 0 if the requested mode isn't present.
fn active_idx_for(ion_tables: &[crate::data::IonModeTable], mode: IonMode) -> usize {
    ion_tables
        .iter()
        .position(|it| it.mode == mode)
        .unwrap_or(0)
}

struct Snapshot {
    method: DamMethod,
    numerator: String,
    denominator: String,
    fc_threshold: f64,
    fdr_threshold: f64,
    delta_threshold: f64,
    export_width_in: f64,
    export_height_in: f64,
    export_dpi: u32,
    rendering: bool,
    has_texture: bool,
    tabs: Vec<IonMode>,
    active_tab: IonMode,
}

fn spawn_render(app: &mut App) {
    let fc_threshold = app.settings.fc_threshold;
    let fdr_threshold = app.settings.fdr_threshold;
    let delta_threshold = app.settings.delta_threshold;
    let export_width_in = app.settings.stage2_export_width_in;
    let export_height_in = app.settings.stage2_export_height_in;
    let export_dpi = app.settings.stage2_export_dpi;
    let ion_tables_for_idx: Vec<IonMode> = app.inputs.ion_tables.iter().map(|it| it.mode).collect();

    let AppState::Stage2DamThreshold {
        dam_results,
        rendering,
        render_rx,
        active_volcano_tab,
        ..
    } = &mut app.state
    else {
        return;
    };
    let active_idx = ion_tables_for_idx
        .iter()
        .position(|m| *m == *active_volcano_tab)
        .unwrap_or(0);
    // Spec scenario "DAM success → Stage 2 threshold": volcano annotation
    // strip reads `dam_results[i].fdr_method` (immutable snapshot from the
    // DAM run), NOT `app.settings.dam_fdr_method`. Prevents future user
    // changes to the Stage 2 setup radio from silently relabeling an
    // already-completed DAM run.
    let fdr_method = dam_results[active_idx].fdr_method;
    let (w_px, h_px) =
        crate::ui::widgets::export_pixels(export_width_in, export_height_in, export_dpi);
    let opts = VolcanoOpts {
        width_px: w_px,
        height_px: h_px,
        fc_threshold,
        fdr_threshold,
        delta_threshold,
        fdr_method,
    };
    let method = dam_results[active_idx].method;
    let result_clone = dam_results[active_idx].clone();
    let mode_for_log = *active_volcano_tab;
    let (tx, rx) = mpsc::channel::<(usize, Result<VolcanoRender, String>)>();
    *render_rx = Some(rx);
    *rendering = true;
    match method {
        DamMethod::BrunnerMunzel => info!(
            method = ?method,
            mode = %mode_for_log,
            fc = fc_threshold,
            fdr = fdr_threshold,
            delta = delta_threshold,
            width_px = w_px,
            height_px = h_px,
            "rendering volcano"
        ),
        DamMethod::Welch | DamMethod::Student => info!(
            method = ?method,
            mode = %mode_for_log,
            fc = fc_threshold,
            fdr = fdr_threshold,
            width_px = w_px,
            height_px = h_px,
            "rendering volcano"
        ),
    }
    app.rt.spawn(async move {
        let r = tokio::task::spawn_blocking(move || render_volcano(&result_clone, &opts))
            .await
            .map_err(|e| e.to_string())
            .and_then(|res| res.map_err(|e| e.to_string()))
            .map(|buf| (buf, w_px, h_px));
        let _ = tx.send((active_idx, r));
    });
}

fn drain_render(app: &mut App, ctx: &egui::Context) {
    let render = {
        let AppState::Stage2DamThreshold {
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
        let Some(rx) = render_rx else {
            return;
        };
        let Ok((idx, msg)) = rx.try_recv() else {
            return;
        };
        *rendering = false;
        *render_rx = None;
        match msg {
            Ok(triple) => Some((idx, triple)),
            Err(e) => {
                error!(error = %e, "volcano render failed");
                None
            }
        }
    };
    let Some((idx, (buffer, w_px, h_px))) = render else {
        return;
    };
    let img = egui::ColorImage::from_rgba_unmultiplied([w_px as usize, h_px as usize], &buffer);
    let handle = ctx.load_texture("volcano", img, egui::TextureOptions::LINEAR);
    if let AppState::Stage2DamThreshold {
        volcano_textures, ..
    } = &mut app.state
        && let Some(slot) = volcano_textures.get_mut(idx)
    {
        *slot = Some(handle);
        info!(
            mode_idx = idx,
            width_px = w_px,
            height_px = h_px,
            "volcano texture uploaded"
        );
    }
}

fn download_png(app: &App) {
    let fc_threshold = app.settings.fc_threshold;
    let fdr_threshold = app.settings.fdr_threshold;
    let delta_threshold = app.settings.delta_threshold;
    let export_width_in = app.settings.stage2_export_width_in;
    let export_height_in = app.settings.stage2_export_height_in;
    let export_dpi = app.settings.stage2_export_dpi;

    let AppState::Stage2DamThreshold {
        dam_results,
        active_volcano_tab,
        ..
    } = &app.state
    else {
        return;
    };
    let ion_tables = app.inputs.ion_tables.as_slice();
    let active_idx = active_idx_for(ion_tables, *active_volcano_tab);
    let dam_result = &dam_results[active_idx];
    // Spec scenario "DAM success → Stage 2 threshold": CSV-export `# FDR:`
    // comment line reads from `dam_results[i].fdr_method` (the value used
    // at run_dam time), NOT from `app.settings.dam_fdr_method`. Same
    // immutable-snapshot reasoning as `spawn_render`.
    let fdr_method = dam_result.fdr_method;
    // Per-tab default filename: volcano-pos.png / volcano-neg.png in dual-mode;
    // bare volcano.png in single-mode.
    let default_name = if dam_results.len() >= 2 {
        match *active_volcano_tab {
            IonMode::Positive => "volcano-pos.png",
            IonMode::Negative => "volcano-neg.png",
        }
    } else {
        "volcano.png"
    };
    let Some(path) = crate::ui::widgets::save_dialog("PNG", "png", default_name) else {
        return;
    };
    let (w_px, h_px) =
        crate::ui::widgets::export_pixels(export_width_in, export_height_in, export_dpi);
    let opts = VolcanoOpts {
        width_px: w_px,
        height_px: h_px,
        fc_threshold,
        fdr_threshold,
        delta_threshold,
        fdr_method,
    };
    if let Err(e) = export_volcano_png(dam_result, &opts, export_dpi, &path) {
        error!(error = %e, "PNG export failed");
    } else {
        info!(
            path = %path.display(),
            mode = %*active_volcano_tab,
            width_in = export_width_in,
            height_in = export_height_in,
            dpi = export_dpi,
            width_px = w_px,
            height_px = h_px,
            "PNG exported"
        );
    }
}

fn download_csv(app: &App, only_dam: bool) {
    let fc_threshold = app.settings.fc_threshold;
    let fdr_threshold = app.settings.fdr_threshold;
    let delta_threshold = app.settings.delta_threshold;
    let AppState::Stage2DamThreshold { dam_results, .. } = &app.state else {
        return;
    };
    let ion_tables = app.inputs.ion_tables.as_slice();
    let default_name = if only_dam {
        "dam.csv"
    } else {
        "all_results.csv"
    };
    let Some(path) = crate::ui::widgets::save_dialog("CSV", "csv", default_name) else {
        return;
    };
    let file = match std::fs::File::create(&path) {
        Ok(f) => f,
        Err(e) => {
            error!(path = %path.display(), error = %e, "could not create CSV file");
            return;
        }
    };
    let result = if only_dam {
        if dam_results.len() >= 2 {
            crate::dam::export::export_dam_csv_multi(
                file,
                ion_tables,
                dam_results,
                fc_threshold,
                fdr_threshold,
                delta_threshold,
            )
        } else {
            crate::dam::export::export_dam_csv(
                file,
                &dam_results[0],
                fc_threshold,
                fdr_threshold,
                delta_threshold,
            )
        }
    } else if dam_results.len() >= 2 {
        crate::dam::export::export_all_csv_multi(
            file,
            ion_tables,
            dam_results,
            fc_threshold,
            fdr_threshold,
            delta_threshold,
        )
    } else {
        crate::dam::export::export_all_csv(
            file,
            &dam_results[0],
            fc_threshold,
            fdr_threshold,
            delta_threshold,
        )
    };
    if let Err(e) = result {
        error!(error = %e, "CSV export failed");
    } else {
        info!(path = %path.display(), only_dam, "CSV exported");
    }
}

/// Write the deduplication audit CSV to a user-chosen path. In single-mode
/// (one `DamResult`) writes the single report. In dual-mode, concatenates
/// per-mode reports separated by `# Mode: <mode>` header lines.
/// Download the dedup-audit CSV for the current `Stage2DamThreshold` run.
/// `pub(crate)` so the Data tab (`ui::data_tab`) can trigger the same export —
/// the dedup summary + download button relocated there by
/// `add-bottom-panel-data-tab`.
pub(crate) fn download_dedup_audit_csv(app: &App) {
    let AppState::Stage2DamThreshold { dam_results, .. } = &app.state else {
        return;
    };
    let ion_tables = app.inputs.ion_tables.as_slice();
    let any_report = dam_results.iter().any(|r| r.dedup_report.is_some());
    if !any_report {
        return;
    }
    let Some(path) = crate::ui::widgets::save_dialog("CSV", "csv", "dedup_audit.csv") else {
        return;
    };
    let mut file = match std::fs::File::create(&path) {
        Ok(f) => f,
        Err(e) => {
            error!(path = %path.display(), error = %e, "could not create dedup audit CSV file");
            return;
        }
    };
    use std::io::Write;
    for (idx, dr) in dam_results.iter().enumerate() {
        let Some(report) = dr.dedup_report.as_ref() else {
            continue;
        };
        // In dual-mode, prefix each report with a `# Mode:` discriminator.
        if dam_results.len() >= 2 {
            let mode = ion_tables
                .get(idx)
                .map(|it| it.mode.to_string())
                .unwrap_or_else(|| format!("mode-{}", idx + 1));
            if idx > 0 {
                let _ = writeln!(file);
            }
            if let Err(e) = writeln!(&mut file, "# Mode: {mode}") {
                error!(error = %e, "could not write Mode header to dedup audit CSV");
                return;
            }
        }
        if let Err(e) = crate::dam::export::export_dedup_audit_csv(&mut file, report) {
            error!(error = %e, "dedup audit CSV export failed");
            return;
        }
    }
    info!(path = %path.display(), modes = dam_results.len(), "dedup audit CSV exported");
}

fn continue_to_enrichment(app: &mut App) {
    let prev = std::mem::take(&mut app.state);
    let AppState::Stage2DamThreshold { dam_results, .. } = prev else {
        return;
    };
    // Per D11 / spec: hard-resets every Stage 3 settings field back to
    // default EVERY Continue (matches today's `continue_to_enrichment`
    // hard-coding the Stage 3 defaults at the call site).
    app.settings.reset_stage3_on_continue_to_enrichment();
    app.state = AppState::Stage3EnrichSetup {
        dam_results,
        error: None,
        kegg_fetch: None,
        modules_fetch: None,
    };
    // Reconcile a snapshot-restored KEGG selection (species/Group) with the
    // empty cache so `Run Enrichment` auto-enables without manual re-selection.
    crate::ui::stage3_setup::rehydrate_stage3_cache(app);
}
