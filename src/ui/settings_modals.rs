//! Settings save / load modals, relocated out of `src/app.rs` so that file
//! stays the state/dispatch layer and `src/ui/` owns the render layer (the
//! same `src/ui/<screen>::show(ui, app)` split the rest of the app follows;
//! `decompose-god-functions` D2).
//!
//! The three entry points take `&mut App`:
//! - [`show_save`] — the Save-settings confirm modal (rendered each frame when
//!   `app.settings_save_modal == Confirming`), invoked from `App::render_modals`.
//! - [`show_load`] — the Load-settings confirm modal (rendered each frame when
//!   `app.settings_load_modal == Confirming`), invoked from `App::render_modals`.
//! - [`open_load`] — opens the OS file picker and transitions
//!   `settings_load_modal` into `Confirming`, invoked from
//!   `App::drain_modal_requests` when the load flag is set.
//!
//! Both confirm modals build their mode / DAM / normalization / enrichment
//! preview through the single [`settings_summary_lines`] formatter (D3) so the
//! two summaries can never drift.

use egui::Context;
use tracing::{error, info};

use crate::app::{
    AnalysisMode, App, SessionSettings, SettingsLoadModalState, SettingsSaveModalState,
};
use crate::theme;

/// The four shared preview lines both settings modals render. Each modal keeps
/// its own label prefix (`"Analysis mode: "` vs `"  Analysis mode:  "`); only
/// the *values* are shared, so the formatter owns the value formatting and the
/// call site owns the layout.
pub(crate) struct SettingsSummary {
    pub mode: String,
    pub dam: String,
    pub normalization: String,
    pub enrichment: String,
}

/// Format the mode / DAM / normalization / enrichment summary values for a
/// `SessionSettings`. Single source of truth shared by the save and load
/// modals; the strings are byte-identical to the pre-`decompose-god-functions`
/// inline copies (the two diverged only in label prefix, which stays at the
/// call site).
pub(crate) fn settings_summary_lines(s: &SessionSettings) -> SettingsSummary {
    let mode = match s.analysis_mode {
        AnalysisMode::Pathway => {
            format!("Pathway ({})", s.kegg_species.as_deref().unwrap_or("—"))
        }
        AnalysisMode::Module => format!(
            "Module ({} / {})",
            s.organism_group_level
                .map(|l| l.to_string())
                .unwrap_or_else(|| "—".to_string()),
            s.organism_group.as_deref().unwrap_or("—"),
        ),
    };
    let dam = format!(
        "{} · FDR({}) · FC≥{} · q<{}",
        s.dam_method.display_name(),
        s.dam_fdr_method.short_label(),
        s.fc_threshold,
        s.fdr_threshold,
    );
    let normalization = format!("{:?}", s.normalization);
    let enrichment = format!(
        "{} · FDR({}) · top {} · min hit {}",
        s.direction.short_label(),
        s.enrichment_fdr_method.short_label(),
        s.top_n,
        s.min_hit_count,
    );
    SettingsSummary {
        mode,
        dam,
        normalization,
        enrichment,
    }
}

pub(crate) fn show_save(app: &mut App, ctx: &Context) {
    use egui::{Align2, Window};

    if !matches!(app.settings_save_modal, SettingsSaveModalState::Confirming) {
        return;
    }

    let mut want_save = false;
    let mut want_cancel = false;

    let summary = settings_summary_lines(&app.settings);
    let input_basenames: Vec<String> = {
        let mut v = Vec::new();
        for t in &app.inputs.ion_tables {
            if let Some(p) = &t.txt_path
                && let Some(name) = p.file_name().and_then(|s| s.to_str())
            {
                v.push(format!("  • {} ({})", name, t.mode));
            }
        }
        if let Some(p) = &app.inputs.csv_path
            && let Some(name) = p.file_name().and_then(|s| s.to_str())
        {
            v.push(format!("  • {} (metadata)", name));
        }
        v
    };

    Window::new("Save current settings to file")
        .collapsible(false)
        .resizable(false)
        .anchor(Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            ui.label("This JSON file will contain the following:");
            ui.add_space(4.0);
            ui.label(format!("Analysis mode: {}", summary.mode));
            ui.label(format!("DAM: {}", summary.dam));
            ui.label(format!("Normalization: {}", summary.normalization));
            ui.label(format!("Enrichment: {}", summary.enrichment));
            if !input_basenames.is_empty() {
                ui.add_space(4.0);
                ui.label("Input files (SHA-256 only — file contents are NOT stored):");
                for b in &input_basenames {
                    ui.label(b);
                }
            }
            ui.add_space(8.0);
            ui.label(
                egui::RichText::new(
                    "`user_note` will be written as an empty string. \
                     You can open the JSON in any text editor afterwards to \
                     add a comment.",
                )
                .small()
                .color(theme::TEXT),
            );
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Cancel").clicked() {
                    want_cancel = true;
                }
                if ui.button("Save…").clicked() {
                    want_save = true;
                }
            });
        });

    if want_cancel {
        app.settings_save_modal = SettingsSaveModalState::Closed;
        return;
    }
    if want_save {
        let default_name = chrono::Local::now()
            .format("metabolopan-settings-%Y-%m-%d_%H%M%S.json")
            .to_string();
        let Some(path) = rfd::FileDialog::new()
            .add_filter("json", &["json"])
            .set_file_name(default_name)
            .save_file()
        else {
            app.settings_save_modal = SettingsSaveModalState::Closed;
            return;
        };
        match crate::session_io::save_to_path(&path, &app.settings, &app.inputs, "") {
            Ok(()) => {
                let bytes = std::fs::metadata(&path).map(|m| m.len()).unwrap_or(0);
                info!(
                    path = %path.display(),
                    bytes = bytes,
                    "settings snapshot written"
                );
            }
            Err(e) => {
                error!(
                    path = %path.display(),
                    error = %e,
                    "settings snapshot write failed"
                );
                app.open_error_toast_for_save(&e);
            }
        }
        app.settings_save_modal = SettingsSaveModalState::Closed;
    }
}

/// Open the OS file picker and, on success, transition to the Load
/// confirm modal carrying the parsed snapshot + hash mismatches +
/// validation resets.
pub(crate) fn open_load(app: &mut App) {
    let Some(path) = rfd::FileDialog::new()
        .add_filter("json", &["json"])
        .pick_file()
    else {
        return; // user cancelled the OS dialog
    };

    let snapshot = match crate::session_io::load_from_path(&path) {
        Ok(s) => s,
        Err(e) => {
            error!(
                path = %path.display(),
                error = %e,
                "settings snapshot load failed"
            );
            app.open_error_toast_for_load(&e);
            return;
        }
    };

    let mismatches = match snapshot.diff_input_hashes(&app.inputs) {
        Ok(m) => m,
        Err(e) => {
            error!(
                path = %path.display(),
                error = %e,
                "settings snapshot load failed"
            );
            app.open_error_toast_for_load(&e);
            return;
        }
    };

    let resets =
        crate::session_io::validate_against_inputs(&snapshot.settings, app.inputs.mapping.as_ref());

    app.settings_load_modal = SettingsLoadModalState::Confirming {
        snapshot,
        mismatches,
        resets,
        path,
    };
}

/// Render the load-settings confirm modal (if Confirming) per
/// `app-shell` spec Load UX.
pub(crate) fn show_load(app: &mut App, ctx: &Context) {
    use egui::{Align2, RichText, Window};

    let SettingsLoadModalState::Confirming {
        snapshot,
        mismatches,
        resets,
        path,
    } = &app.settings_load_modal
    else {
        return;
    };

    let mut want_apply = false;
    let mut want_cancel = false;

    let filename = path
        .file_name()
        .and_then(|s| s.to_str())
        .unwrap_or("<unknown>")
        .to_string();
    let saved_at_local = chrono::DateTime::parse_from_rfc3339(&snapshot.saved_at)
        .map(|dt| {
            dt.with_timezone(&chrono::Local)
                .format("%Y-%m-%d %H:%M:%S %z")
                .to_string()
        })
        .unwrap_or_else(|_| snapshot.saved_at.clone());

    let current_app_version = env!("CARGO_PKG_VERSION");
    let show_version_hint = snapshot.app_version != current_app_version;

    let summary = settings_summary_lines(&snapshot.settings);

    // Build reset rows (snapshot of borrowed data — drop borrow
    // before mutating `app.settings`).
    let reset_rows: Vec<(&'static str, String)> = {
        let mut rows = Vec::new();
        if let Some(v) = &resets.numerator {
            rows.push(("Numerator group", v.clone()));
        }
        if let Some(v) = &resets.denominator {
            rows.push(("Denominator group", v.clone()));
        }
        if let Some(v) = &resets.metadata_column {
            rows.push(("Metadata column", v.clone()));
        }
        if let Some(v) = &resets.pqn_reference_group {
            rows.push(("PQN reference group", v.clone()));
        }
        rows
    };
    let mismatch_rows: Vec<(String, String, Option<String>)> = mismatches
        .iter()
        .map(|m| {
            let saved_prefix = m.saved_sha256.get(..8).unwrap_or("").to_string();
            let current_label = m
                .current
                .as_ref()
                .map(|(_, sha)| sha.get(..8).unwrap_or("").to_string());
            (m.saved_name.clone(), saved_prefix, current_label)
        })
        .collect();
    let inputs_empty = app.inputs.ion_tables.is_empty() && app.inputs.csv_path.is_none();
    let user_note = snapshot.user_note.clone();
    let saved_app_version = snapshot.app_version.clone();

    Window::new("Load settings from file")
        .collapsible(false)
        .resizable(false)
        .anchor(Align2::CENTER_CENTER, [0.0, 0.0])
        .show(ctx, |ui| {
            ui.label(format!("File:     {}", filename));
            ui.label(format!("Saved at: {}", saved_at_local));
            if show_version_hint {
                ui.label(format!(
                    "App ver:  {}  (current: {})",
                    saved_app_version, current_app_version
                ));
            }
            if !user_note.is_empty() {
                ui.label(format!("Note:     {:?}", user_note));
            }

            ui.add_space(6.0);
            ui.label(
                RichText::new("Settings preview:")
                    .strong()
                    .color(theme::HEADING),
            );
            ui.label(format!("  Analysis mode:  {}", summary.mode));
            ui.label(format!("  DAM:            {}", summary.dam));
            ui.label(format!("  Normalization:  {}", summary.normalization));
            ui.label(format!("  Enrichment:     {}", summary.enrichment));

            if !reset_rows.is_empty() {
                ui.add_space(6.0);
                ui.colored_label(
                    theme::WARNING,
                    format!(
                        "⚠ {} field(s) will be reset because they don't match the loaded metadata:",
                        reset_rows.len()
                    ),
                );
                for (label, value) in &reset_rows {
                    ui.label(format!(
                        "    • {} \"{}\"  -> re-select required",
                        label, value
                    ));
                }
            }

            if !mismatch_rows.is_empty() {
                ui.add_space(6.0);
                ui.colored_label(
                    theme::WARNING,
                    format!(
                        "⚠ {} input file(s) differ from the saved snapshot:",
                        mismatch_rows.len()
                    ),
                );
                for (name, saved, current) in &mismatch_rows {
                    let current_text = match current {
                        Some(c) => format!("current={}…", c),
                        None => "current=<not loaded>".to_string(),
                    };
                    ui.label(format!(
                        "    • {}  saved={}…  {}",
                        name, saved, current_text
                    ));
                }
                ui.label(
                    RichText::new("The settings will still apply if you continue.")
                        .small()
                        .color(theme::TEXT),
                );
            } else if inputs_empty && !snapshot.input_files.is_empty() {
                // Reachable when snapshot expected files but user has none — diff already
                // surfaces those as mismatch_rows; this branch is for symmetric clarity
                // when input_files is non-empty but every entry happens to land in
                // mismatch_rows.
            } else if inputs_empty {
                ui.add_space(6.0);
                ui.colored_label(
                    theme::TEXT,
                    "ℹ No input files currently loaded; hash check skipped.",
                );
            }

            ui.add_space(8.0);
            ui.horizontal(|ui| {
                if ui.button("Cancel").clicked() {
                    want_cancel = true;
                }
                if ui.button("Apply settings").clicked() {
                    want_apply = true;
                }
            });
        });

    if want_cancel {
        app.settings_load_modal = SettingsLoadModalState::Closed;
        return;
    }
    if want_apply {
        // Take ownership so we can move snapshot.settings without
        // hitting the borrow checker.
        let taken = std::mem::take(&mut app.settings_load_modal);
        if let SettingsLoadModalState::Confirming {
            snapshot,
            mismatches,
            resets,
            path,
        } = taken
        {
            let m_count = mismatches.len();
            let r_count = [
                &resets.numerator,
                &resets.denominator,
                &resets.metadata_column,
                &resets.pqn_reference_group,
            ]
            .iter()
            .filter(|x| x.is_some())
            .count();
            app.settings.apply_snapshot(snapshot.settings, &resets);
            info!(
                path = %path.display(),
                mismatches = m_count,
                resets = r_count,
                "settings snapshot loaded"
            );
        }
        // settings_load_modal is now Closed (default) via mem::take.
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::app::SessionSettings;

    #[test]
    fn settings_summary_lines_formats_mode_and_normalization() {
        // Pathway mode renders the species code (or "—" when unset).
        let mut s = SessionSettings {
            analysis_mode: AnalysisMode::Pathway,
            kegg_species: Some("hsa".to_string()),
            ..SessionSettings::default()
        };
        let out = settings_summary_lines(&s);
        assert_eq!(out.mode, "Pathway (hsa)");
        // Normalization is the Debug rendering of the chosen method.
        assert_eq!(out.normalization, format!("{:?}", s.normalization));
        // DAM / enrichment carry the documented separator vocabulary so a
        // future label rename can't silently drop a field.
        assert!(out.dam.contains("FDR(") && out.dam.contains("FC≥") && out.dam.contains("q<"));
        assert!(
            out.enrichment.contains("FDR(")
                && out.enrichment.contains("top ")
                && out.enrichment.contains("min hit ")
        );

        // Pathway mode with no species selected falls back to the em-dash.
        s.kegg_species = None;
        assert_eq!(settings_summary_lines(&s).mode, "Pathway (—)");

        // Module mode renders the "(level / group)" pair, em-dash when unset.
        s.analysis_mode = AnalysisMode::Module;
        s.organism_group_level = None;
        s.organism_group = None;
        assert_eq!(settings_summary_lines(&s).mode, "Module (— / —)");
    }
}
