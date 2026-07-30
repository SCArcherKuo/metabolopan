//! Stage 3 (coverage route) — the coverage result screen.
//!
//! Provenance funnel, the interpretation statement, four live display filters,
//! the sortable results table, the dot plot with its PNG export, the CSV
//! download, and `Start a new analysis`.
//!
//! Every filter on this screen is a pure filter over the `CoverageResult.rows`
//! already in hand. None of them re-runs `coverage::compute`, issues a PubChem
//! request, or issues a KEGG request — there is nothing on this screen that
//! requires recomputation.
//!
//! See the `coverage-ui` and `kegg-coverage` capability specs.

use egui::RichText;
use std::sync::mpsc;
use tracing::{error, info};

use crate::app::{AnalysisMode, App, AppState, CoverageFunnel};
use crate::coverage::{CoverageResult, CoverageRow, CoverageSortKey, displayed_rows};
use crate::plot::{CoverageDotplotOpts, export_coverage_dotplot_png, render_coverage_dotplot};
use crate::theme;
use crate::ui::widgets::primary_button;

/// The inline interpretation statement. Pinned verbatim by the `kegg-coverage`
/// capability's no-statistical-test requirement.
pub(crate) const INTERPRETATION: &str = "Descriptive coverage only — no statistical test is \
    performed. A high coverage percentage means many of this entry's compounds were detected, \
    which reflects both biology and what your method can detect.";

/// How many cpd IDs the `Compounds` cell shows before the `…(+N)` marker.
const INLINE_COMPOUNDS: usize = 4;
/// How many the hover shows — more than the cell, per the spec; the full list
/// is in the CSV.
const HOVER_COMPOUNDS: usize = 40;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Action {
    None,
    Redraw,
    DownloadPng,
    DownloadCsv,
    RequestNewRound,
    ConfirmNewRound,
    CancelNewRound,
}

pub fn show(ui: &mut egui::Ui, app: &mut App) {
    drain_render(app, ui.ctx());

    // Snapshot everything the body reads, so the `&app.state` borrow ends
    // before the filter controls take `&mut app.settings`.
    let snap = match &app.state {
        AppState::Stage3CoverageResult {
            coverage_result,
            funnel,
            rendering,
            confirming_new_round,
            dotplot_tex,
            ..
        } => Some((
            coverage_result.clone(),
            *funnel,
            *rendering,
            *confirming_new_round,
            dotplot_tex.clone(),
        )),
        _ => None,
    };
    let Some((result, funnel, rendering, confirming_new_round, dotplot_tex)) = snap else {
        return;
    };

    let mut action = Action::None;
    let mut new_w_in = app.settings.stage3_export_width_in;
    let mut new_h_in = app.settings.stage3_export_height_in;
    let mut new_dpi = app.settings.stage3_export_dpi;

    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.heading(RichText::new("Stage 3 — Coverage").color(theme::HEADING));
            ui.add_space(6.0);

            ui.label(RichText::new(funnel_line(&funnel, &result)).color(theme::TEXT));
            ui.add_space(6.0);
            ui.label(RichText::new(INTERPRETATION).color(theme::TEXT_SECONDARY));
            ui.add_space(10.0);

            render_filters(ui, app, &result);
            ui.add_space(10.0);

            // Recomputed AFTER the filter controls so a change lands on the
            // same frame the user made it, not the next one.
            let filters = app.settings.coverage_display_filters();
            let rows = displayed_rows(&result, filters);
            render_table(ui, app, &rows);
            ui.add_space(10.0);

            // ── Dot plot + exports ──
            //
            // Every one of these is disabled while a render is in flight; the
            // filter controls above are deliberately NOT, because they mutate
            // settings only and the next render picks up the new values.
            crate::ui::widgets::png_export_size_controls(
                ui,
                &mut new_w_in,
                &mut new_h_in,
                &mut new_dpi,
            );
            // The label matches the enrichment result screen byte-for-byte: two
            // screens with the same function must not use different words.
            if primary_button(ui, "Draw dot plot", !rendering).clicked() {
                action = Action::Redraw;
            }
            if let Some(tex) = dotplot_tex.as_ref() {
                let size = tex.size_vec2();
                let avail = ui.available_width().max(1.0);
                let scale = (avail / size.x).min(1.0);
                ui.add(egui::Image::new(tex).fit_to_exact_size(size * scale));
            }
            ui.add_space(6.0);
            if ui
                .add_enabled(!rendering, egui::Button::new("Download dot plot PNG"))
                .clicked()
            {
                action = Action::DownloadPng;
            }
            if ui
                .add_enabled(!rendering, egui::Button::new("Download coverage CSV"))
                .clicked()
            {
                action = Action::DownloadCsv;
            }

            ui.add_space(10.0);
            // `!rendering`-gated so it cannot fire mid-render, which is what
            // keeps `rendering == true` unreachable while the confirm is open.
            if ui
                .add_enabled(!rendering, egui::Button::new("Start a new analysis"))
                .clicked()
            {
                action = Action::RequestNewRound;
            }
        });

    if confirming_new_round {
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
                        "The current coverage results and any un-downloaded plots or CSV will \
                         be lost. This cannot be undone.",
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

    // The Height field is a user override once hand-edited, exactly as on the
    // enrichment result screen: while it is untouched, each draw re-fits it to
    // the live displayed-row count.
    if new_h_in != app.settings.stage3_export_height_in
        && let AppState::Stage3CoverageResult {
            height_user_overridden,
            ..
        } = &mut app.state
    {
        *height_user_overridden = true;
    }
    app.settings.stage3_export_width_in = new_w_in;
    app.settings.stage3_export_height_in = new_h_in;
    app.settings.stage3_export_dpi = new_dpi;

    match action {
        Action::None => {}
        Action::Redraw => spawn_render(app),
        Action::DownloadPng => download_png(app),
        Action::DownloadCsv => download_csv(app),
        Action::RequestNewRound => set_confirming_new_round(app, true),
        Action::ConfirmNewRound => app.start_new_round(),
        Action::CancelNewRound => set_confirming_new_round(app, false),
    }
}

fn set_confirming_new_round(app: &mut App, value: bool) {
    if let AppState::Stage3CoverageResult {
        confirming_new_round,
        ..
    } = &mut app.state
    {
        *confirming_new_round = value;
    }
}

/// Assemble the renderer options from the live settings and state.
///
/// One builder for the preview and the export, so a Download taken right after
/// a filter change — with no intervening Draw — produces the same image the
/// plot would show.
fn build_opts(app: &App, w_px: u32, h_px: u32) -> CoverageDotplotOpts {
    CoverageDotplotOpts {
        width_px: w_px,
        height_px: h_px,
        filters: app.settings.coverage_display_filters(),
        mode_label: match app.settings.analysis_mode {
            AnalysisMode::Pathway => "Pathway".to_string(),
            AnalysisMode::Module => "Module".to_string(),
        },
        target_label: target_label(app),
        detected_total: match &app.state {
            AppState::Stage3CoverageResult {
                coverage_result, ..
            } => coverage_result.detected_total,
            _ => 0,
        },
        group_record: group_record_counts(app),
    }
}

/// Species code in Pathway mode; `"<Level> / <Group>"` in Module mode.
fn target_label(app: &App) -> String {
    match app.settings.analysis_mode {
        AnalysisMode::Pathway => app
            .settings
            .kegg_species
            .clone()
            .unwrap_or_else(|| "—".to_string()),
        AnalysisMode::Module => format!(
            "{} / {}",
            app.settings
                .organism_group_level
                .map(|l| l.to_string())
                .unwrap_or_else(|| "—".into()),
            app.settings.organism_group.as_deref().unwrap_or("—")
        ),
    }
}

/// The compact `(selected, total, threshold)` the annotation strip renders.
/// `None` when no metadata `.csv` was supplied.
fn group_record_counts(app: &App) -> Option<(usize, usize, f64)> {
    let mapping = app.inputs.mapping.as_ref()?;
    let total = mapping
        .groups()
        .into_iter()
        .filter(|g| g != crate::data::groups::UNASSIGNED)
        .count();
    let selected = app
        .settings
        .coverage_selected_groups
        .as_ref()
        .map_or(total, Vec::len);
    Some((selected, total, app.settings.coverage_presence_threshold))
}

/// The full `(selected, all, threshold)` the CSV records. Unlike the plot's
/// compact form this carries the NAMES: a reader comparing two exports needs to
/// see which groups were dropped, not just how many.
fn group_record_names(app: &App) -> Option<(Vec<String>, Vec<String>, f64)> {
    let mapping = app.inputs.mapping.as_ref()?;
    let all: Vec<String> = mapping
        .groups()
        .into_iter()
        .filter(|g| g != crate::data::groups::UNASSIGNED)
        .collect();
    let selected = app
        .settings
        .coverage_selected_groups
        .clone()
        .unwrap_or_else(|| all.clone());
    Some((selected, all, app.settings.coverage_presence_threshold))
}

/// Re-fit the export height to the rows this draw will actually show, unless
/// the user has hand-set it. Same rule and same helper as the enrichment result
/// screen: without it, loosening a filter and redrawing crams more rows into a
/// height sized for the previous view.
fn effective_height_in(app: &App) -> f64 {
    let AppState::Stage3CoverageResult {
        coverage_result,
        height_user_overridden,
        ..
    } = &app.state
    else {
        return app.settings.stage3_export_height_in;
    };
    if *height_user_overridden {
        return app.settings.stage3_export_height_in;
    }
    let displayed = displayed_rows(coverage_result, app.settings.coverage_display_filters()).len();
    crate::app::stage3_autosize_height_in(app.settings.top_n, displayed)
}

fn spawn_render(app: &mut App) {
    let h_in = effective_height_in(app);
    app.settings.stage3_export_height_in = h_in;
    let (w_px, h_px) = crate::ui::widgets::export_pixels(
        app.settings.stage3_export_width_in,
        h_in,
        app.settings.stage3_export_dpi,
    );
    let opts = build_opts(app, w_px, h_px);

    let AppState::Stage3CoverageResult {
        coverage_result,
        rendering,
        render_rx,
        ..
    } = &mut app.state
    else {
        return;
    };
    let result_clone = coverage_result.clone();
    let (tx, rx) = mpsc::channel::<Result<crate::app::DotplotRender, String>>();
    *render_rx = Some(rx);
    *rendering = true;
    info!(
        width_px = w_px,
        height_px = h_px,
        "rendering coverage dot plot"
    );
    app.rt.spawn(async move {
        let r = tokio::task::spawn_blocking(move || render_coverage_dotplot(&result_clone, &opts))
            .await
            .map_err(|e| e.to_string())
            .and_then(|res| res.map_err(|e| e.to_string()))
            .map(|buf| (buf, w_px, h_px));
        let _ = tx.send(r);
    });
}

fn drain_render(app: &mut App, ctx: &egui::Context) {
    let r = {
        let AppState::Stage3CoverageResult {
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
                error!(error = %e, "coverage dot plot render failed");
                None
            }
        }
    };
    let Some((buf, w, h)) = r else { return };
    let img = egui::ColorImage::from_rgba_unmultiplied([w as usize, h as usize], &buf);
    let handle = ctx.load_texture("coverage_dotplot", img, egui::TextureOptions::LINEAR);
    if let AppState::Stage3CoverageResult { dotplot_tex, .. } = &mut app.state {
        *dotplot_tex = Some(handle);
        info!(
            width_px = w,
            height_px = h,
            "coverage dot plot texture uploaded"
        );
    }
}

fn download_png(app: &App) {
    let dpi = app.settings.stage3_export_dpi;
    let (w_px, h_px) = crate::ui::widgets::export_pixels(
        app.settings.stage3_export_width_in,
        effective_height_in(app),
        dpi,
    );
    let opts = build_opts(app, w_px, h_px);
    let AppState::Stage3CoverageResult {
        coverage_result, ..
    } = &app.state
    else {
        return;
    };
    let Some(path) = crate::ui::widgets::save_dialog("PNG", "png", "coverage-dotplot.png") else {
        return;
    };
    if let Err(e) = export_coverage_dotplot_png(coverage_result, &opts, dpi, &path) {
        error!(error = %e, "coverage dot plot PNG export failed");
    } else {
        info!(
            path = %path.display(),
            width_px = w_px,
            height_px = h_px,
            dpi,
            "coverage dot plot PNG exported"
        );
    }
}

fn download_csv(app: &App) {
    let target_label = target_label(app);
    let group_record = group_record_names(app);
    let filters = app.settings.coverage_display_filters();
    let mode = app.settings.analysis_mode;
    let AppState::Stage3CoverageResult {
        coverage_result,
        cpd_to_names,
        ..
    } = &app.state
    else {
        return;
    };
    let Some(path) = crate::ui::widgets::save_dialog("CSV", "csv", "coverage.csv") else {
        return;
    };
    let mut file = match std::fs::File::create(&path) {
        Ok(f) => f,
        Err(e) => {
            error!(path = %path.display(), error = %e, "could not create coverage CSV file");
            return;
        }
    };
    let ctx = crate::coverage::export::CoverageExportContext {
        mode,
        target_label,
        group_record,
        cpd_to_names,
        filters,
    };
    if let Err(e) = crate::coverage::export::export_coverage_csv(&mut file, coverage_result, &ctx) {
        error!(error = %e, "coverage CSV export failed");
    } else {
        info!(path = %path.display(), "coverage CSV exported");
    }
}

/// The provenance funnel: the resolution chain from raw features to compounds
/// reaching an entry, one term per stage in the order the stages run.
///
/// The `in selected groups` term is OMITTED entirely when no metadata `.csv`
/// was supplied — the stage did not run, and repeating the raw count there
/// would print a tautology as if it were a measurement.
///
/// Every term is sourced from `CoverageFunnel` and `CoverageResult`, never from
/// `Stage3Funnel`, which carries a foreground branch this route has no meaning
/// for.
pub(crate) fn funnel_line(funnel: &CoverageFunnel, result: &CoverageResult) -> String {
    let mut parts = vec![format!("{} features", thousands(funnel.raw_features))];
    if let Some(n) = funnel.in_selected_groups {
        parts.push(format!("{} in selected groups", thousands(n)));
    }
    parts.push(format!(
        "{} after deduplication",
        thousands(funnel.after_dedup)
    ));
    parts.push(format!(
        "{} InChIKeys",
        thousands(funnel.detected_inchikeys)
    ));
    parts.push(format!("{} CIDs", thousands(funnel.detected_cids)));
    parts.push(format!(
        "{} KEGG compounds",
        thousands(result.detected_total)
    ));
    parts.push(format!(
        "{} in at least one entry",
        thousands(result.detected_in_entries)
    ));
    parts.join(" -> ")
}

/// Thousands separators, so a five-digit feature count is readable at a glance.
fn thousands(n: usize) -> String {
    let s = n.to_string();
    let bytes = s.as_bytes();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

/// The four live display filters. All write settings only; the next render
/// picks them up, so they stay enabled even while a plot render is in flight.
fn render_filters(ui: &mut egui::Ui, app: &mut App, result: &CoverageResult) {
    let settings = &mut app.settings;

    ui.horizontal(|ui| {
        ui.label("Sort by");
        // Writes the SAME field the clickable column headers write; there is no
        // second source of truth for the sort key, so the two can never
        // disagree.
        egui::ComboBox::from_id_salt("coverage_sort_by")
            .selected_text(sort_key_label(settings.coverage_sort_key))
            .show_ui(ui, |ui| {
                for key in OFFERED_SORT_KEYS {
                    ui.selectable_value(&mut settings.coverage_sort_key, key, sort_key_label(key));
                }
            });
    });

    ui.horizontal(|ui| {
        ui.label("Top N entries:");
        ui.add(
            egui::DragValue::new(&mut settings.top_n)
                .speed(1)
                .range(1..=500),
        );
    });

    ui.horizontal(|ui| {
        ui.label("Minimum hit count:");
        ui.add(
            egui::DragValue::new(&mut settings.min_hit_count)
                .speed(1)
                .range(0..=100),
        );
    });

    ui.horizontal(|ui| {
        ui.label("Minimum entry size:");
        // Hard minimum 1 at the input, matching the load-boundary clamp, so a
        // zero-compound entry can never be displayed, plotted, or exported.
        ui.add(
            egui::DragValue::new(&mut settings.coverage_min_entry_size)
                .speed(1)
                .range(crate::app::MIN_COVERAGE_ENTRY_SIZE..=200),
        );
    });

    // What the floor is hiding. Without this the ~20 % of a species catalogue
    // that KEGG annotates with no compounds — every global/overview map among
    // them — would vanish with no account of it anywhere.
    if result.entries_without_compounds > 0 {
        ui.label(
            RichText::new(format!(
                "{} of {} entries have no compounds in KEGG and are never shown.",
                result.entries_without_compounds, result.entries_total
            ))
            .small()
            .color(theme::TEXT_SECONDARY),
        );
    }
}

/// The sort keys the UI currently offers.
///
/// `CoverageSortKey` also defines `EntrySize` and `EntryId`, and
/// `displayed_rows` orders by them correctly — they are simply not offered yet.
/// The two that are offered both answer a question about the RESULT ("how much
/// did I see", "how much is there"); `EntryId` is an alphabetical listing and
/// `EntrySize` reorders by a property of the catalogue rather than of the data,
/// and neither has earned a slot in a four-control filter row. Re-offering
/// either one is adding it to this array — the enum, the chain, and the
/// serialised snapshot form already support them.
const OFFERED_SORT_KEYS: [CoverageSortKey; 2] = [CoverageSortKey::Coverage, CoverageSortKey::Hits];

fn sort_key_label(key: CoverageSortKey) -> &'static str {
    match key {
        CoverageSortKey::Coverage => "Coverage",
        CoverageSortKey::Hits => "Hits",
        CoverageSortKey::EntrySize => "Entry size",
        CoverageSortKey::EntryId => "Entry ID",
    }
}

/// The results table.
///
/// `Hits` always renders BOTH numbers, and `Coverage` is never rendered without
/// `Hits` beside it: a coverage percentage divorced from its denominator is the
/// single most misreadable number on this screen, so keeping the fraction
/// adjacent is a display rule rather than a styling choice.
fn render_table(ui: &mut egui::Ui, app: &mut App, rows: &[&CoverageRow]) {
    let active = app.settings.coverage_sort_key;
    let mut clicked: Option<CoverageSortKey> = None;

    egui::ScrollArea::both()
        .id_salt("coverage_table")
        .max_height(360.0)
        .auto_shrink([false, false])
        .show(ui, |ui| {
            egui::Grid::new("coverage_rows")
                .striped(true)
                .num_columns(5)
                .show(ui, |ui| {
                    // Only the OFFERED keys get a clickable header, so the two
                    // sort surfaces agree on what is sortable: a header that
                    // set a key the `Sort by` selector cannot show would be a
                    // state the user could enter and not leave.
                    ui.label(RichText::new("Entry").color(theme::HEADING));
                    ui.label(RichText::new("Name").color(theme::HEADING));
                    if sort_header(ui, "Hits", CoverageSortKey::Hits, active) {
                        clicked = Some(CoverageSortKey::Hits);
                    }
                    if sort_header(ui, "Coverage", CoverageSortKey::Coverage, active) {
                        clicked = Some(CoverageSortKey::Coverage);
                    }
                    ui.label(RichText::new("Compounds").color(theme::HEADING));
                    ui.end_row();

                    for row in rows {
                        ui.label(&row.entry_id);
                        ui.label(&row.entry_name).on_hover_text(&row.entry_name);
                        ui.label(format!("{} / {}", row.hits, row.entry_size));
                        ui.label(format!("{:.1} %", row.coverage * 100.0));
                        ui.label(compounds_cell(&row.hit_compounds))
                            .on_hover_text(compounds_hover(&row.hit_compounds));
                        ui.end_row();
                    }
                });

            if rows.is_empty() {
                ui.label(
                    RichText::new("No entries match the current filters.")
                        .color(theme::TEXT_SECONDARY),
                );
            }
        });

    if let Some(key) = clicked {
        app.settings.coverage_sort_key = key;
    }
}

/// A clickable column header, marked when it is the active sort column.
/// Returns `true` on click.
///
/// The active marker is the WORD `desc`, not a `▾` glyph: the default egui font
/// has no arrow glyph and rendered it as a tofu box on macOS — the same
/// constraint that forces the ASCII `>` in the stepper. Spelling the direction
/// out also says which way the column is sorted, which the arrow only implied.
/// Every offered key sorts descending (see `kegg-coverage`), so the word is
/// currently constant; it is built from the key rather than hard-coded so a
/// future ascending key cannot silently mislabel itself.
fn sort_header(
    ui: &mut egui::Ui,
    label: &str,
    key: CoverageSortKey,
    active: CoverageSortKey,
) -> bool {
    let text = if key == active {
        RichText::new(format!("{label} ({})", sort_direction_label(key))).color(theme::PRIMARY)
    } else {
        RichText::new(label).color(theme::HEADING)
    };
    ui.add(egui::Label::new(text).sense(egui::Sense::click()))
        .clicked()
}

/// Which way a key sorts. All keys sort descending except `EntryId`.
fn sort_direction_label(key: CoverageSortKey) -> &'static str {
    match key {
        CoverageSortKey::EntryId => "asc",
        _ => "desc",
    }
}

/// The inline `Compounds` cell: KEGG cpd IDs ONLY (metabolite names live in the
/// CSV), truncated with a marker naming how many were omitted.
pub(crate) fn compounds_cell(ids: &[String]) -> String {
    render_compounds(ids, INLINE_COMPOUNDS)
}

/// The hover rendering — more IDs than the cell, per the spec. The full list is
/// always available in the CSV.
pub(crate) fn compounds_hover(ids: &[String]) -> String {
    render_compounds(ids, HOVER_COMPOUNDS)
}

fn render_compounds(ids: &[String], limit: usize) -> String {
    if ids.is_empty() {
        return "—".to_string();
    }
    if ids.len() <= limit {
        return ids.join(", ");
    }
    format!("{}, …(+{})", ids[..limit].join(", "), ids.len() - limit)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The statement names what the numbers are NOT, which is the point of
    /// pinning it: an edit that drops "no statistical test" would leave a
    /// coverage percentage looking like a result.
    #[test]
    fn interpretation_states_no_statistical_test() {
        let norm = INTERPRETATION
            .split_whitespace()
            .collect::<Vec<_>>()
            .join(" ");
        assert_eq!(
            norm,
            "Descriptive coverage only — no statistical test is performed. A high coverage \
             percentage means many of this entry's compounds were detected, which reflects \
             both biology and what your method can detect."
        );
    }

    fn result_with(detected_total: usize, detected_in_entries: usize) -> CoverageResult {
        CoverageResult {
            rows: vec![],
            detected_total,
            entries_total: 0,
            entries_without_compounds: 0,
            detected_in_entries,
        }
    }

    /// The full chain, with thousands separators.
    #[test]
    fn funnel_line_reports_every_stage_in_order() {
        let funnel = CoverageFunnel {
            raw_features: 12431,
            in_selected_groups: Some(9004),
            after_dedup: 2317,
            detected_inchikeys: 1974,
            detected_cids: 4187,
        };
        assert_eq!(
            funnel_line(&funnel, &result_with(391, 264)),
            "12,431 features -> 9,004 in selected groups -> 2,317 after deduplication \
             -> 1,974 InChIKeys -> 4,187 CIDs -> 391 KEGG compounds -> 264 in at least one entry"
        );
    }

    /// With no `.csv` the group term is OMITTED — not rendered with the raw
    /// count repeated, which would read as a measurement that never happened.
    #[test]
    fn funnel_line_omits_the_group_term_when_the_stage_did_not_run() {
        let funnel = CoverageFunnel {
            raw_features: 12431,
            in_selected_groups: None,
            after_dedup: 2317,
            detected_inchikeys: 1974,
            detected_cids: 4187,
        };
        let line = funnel_line(&funnel, &result_with(391, 264));
        assert!(!line.contains("in selected groups"));
        assert!(line.starts_with("12,431 features -> 2,317 after deduplication"));
    }

    #[test]
    fn thousands_groups_digits() {
        assert_eq!(thousands(0), "0");
        assert_eq!(thousands(999), "999");
        assert_eq!(thousands(1000), "1,000");
        assert_eq!(thousands(12431), "12,431");
        assert_eq!(thousands(1234567), "1,234,567");
    }

    fn ids(n: usize) -> Vec<String> {
        (0..n).map(|i| format!("C{i:05}")).collect()
    }

    /// The cell truncates with a marker accounting for the rest; the hover
    /// shows more than the cell.
    #[test]
    fn compounds_cell_truncates_and_accounts_for_the_remainder() {
        let eighteen = ids(18);
        let cell = compounds_cell(&eighteen);
        assert!(cell.ends_with("…(+14)"), "got {cell}");
        assert_eq!(cell.matches("C0").count(), INLINE_COMPOUNDS);

        let hover = compounds_hover(&eighteen);
        assert!(!hover.contains("…"), "18 fits under the hover limit");
        assert_eq!(hover.matches("C0").count(), 18);
    }

    /// No metabolite name ever reaches the cell — those are a CSV concern.
    #[test]
    fn compounds_cell_renders_ids_only() {
        assert_eq!(compounds_cell(&ids(2)), "C00000, C00001");
        assert_eq!(compounds_cell(&[]), "—");
    }
}
