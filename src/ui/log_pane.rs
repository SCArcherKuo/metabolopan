use egui::{Color32, RichText, ScrollArea, TextStyle};
use std::time::SystemTime;
use tracing::Level;

use crate::logging::LogStore;
use crate::theme;

/// Which tab of the bottom panel is currently selected. The bottom panel is a
/// two-tab container (`add-bottom-panel-data-tab`): a stage-aware **Data**
/// summary (default) and the **Log** event stream. UI ephemera — held on
/// `LogPaneState`, never persisted (not part of the `session-settings-io`
/// schema).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum BottomTab {
    #[default]
    Data,
    Log,
}

#[derive(Debug, Default)]
pub struct LogPaneState {
    /// Active bottom-panel tab. Defaults to `Data` so the per-stage summary is
    /// visible from first paint (per the `data-summary-panel` capability).
    pub active_tab: BottomTab,
    pub filter_directive: String,
    /// Set to `true` when the user clicks `[Download bug report…]`.
    /// `App::update()` consumes the flag and opens the confirm modal.
    /// Reset to `false` after consumption so the next frame doesn't
    /// re-open. Field lives alongside `filter_directive` on the same
    /// struct, per design D6 / D8.
    pub bundle_export_requested: bool,
    /// Set to `true` when the user clicks `[Save settings…]` in the Data
    /// tab toolbar (relocated from the Log pane by
    /// `move-settings-buttons-to-data-tab`).
    /// `App::update()` drains the flag every frame regardless of modal
    /// outcome (mutual-exclusion rule may drop the request); see the
    /// `app-shell` spec.
    pub settings_save_requested: bool,
    /// Set to `true` when the user clicks `[Load settings…]` in the Data
    /// tab toolbar (relocated from the Log pane by
    /// `move-settings-buttons-to-data-tab`).
    /// `App::update()` drains the flag every frame; only honored when
    /// `AppState == Stage1Input` (Load button is greyed elsewhere).
    pub settings_load_requested: bool,
    /// Set to `true` when the user clicks `Refresh PubChem cache` in the Data
    /// tab's Stage 3 result Cache-data block (`apply-ui-design-md-tweaks`).
    /// `stage3_result::show` drains it into the existing confirm-modal refresh
    /// flow the same frame (Data tab renders before the central panel).
    pub refresh_pubchem_requested: bool,
    /// Set to `true` when the user clicks `Refresh KEGG conv cache` in the Data
    /// tab's Stage 3 result Cache-data block. Drained by `stage3_result::show`.
    pub refresh_kegg_conv_requested: bool,
    /// Set to `true` when the user clicks `Re-run enrichment` in the Data tab's
    /// Stage 3 result Cache-data block. Drained by `stage3_result::show`.
    pub rerun_enrichment_requested: bool,
    /// Set to `true` when the user clicks the mode-aware `Refresh KEGG
    /// module/pathway cache` button in the Data tab's Stage 3 result Cache-data
    /// block. The result state lacks the catalogue-fetch progress infra, so
    /// `stage3_result::show` drains it by navigating back to Stage 3 setup
    /// (Group/species preserved) and triggering the force re-fetch there.
    pub refresh_catalogue_requested: bool,
    /// Set to `true` when the user clicks `Refresh KEGG organism list` in the
    /// Data tab's Cache-data block (the roster is mode- and route-independent,
    /// so the button renders on five `AppState` variants). Drained by
    /// `App::drain_frame_dialogs` every frame, which opens the refresh confirm
    /// (NOT an App-level modal — see the `app-shell` organism-roster refresh
    /// requirement). Frame-owned rather than screen-owned precisely so the set
    /// of screens that can produce it equals the set that consumes it.
    pub organisms_refresh_requested: bool,
    /// `true` while the organism-roster refresh confirm dialog is open. Lives
    /// here (like `RefreshState`) rather than as an App-level `*ModalState`, so
    /// it stays outside the four-modal mutual-exclusion rule and
    /// `drain_modal_requests` — App-owned, but not one of the four.
    ///
    /// Rendered unconditionally by `App::drain_frame_dialogs`, so an unanswered
    /// confirm follows the user across navigation instead of ghosting;
    /// `App::start_new_round` is the only thing that clears it unanswered.
    pub organisms_refresh_confirm_open: bool,
}

fn format_timestamp(ts: SystemTime) -> String {
    let dt: chrono::DateTime<chrono::Local> = ts.into();
    dt.format("%H:%M:%S%.3f %z").to_string()
}

fn level_color(level: Level) -> Color32 {
    match level {
        Level::ERROR => theme::ERROR,
        Level::WARN => theme::WARNING,
        // INFO is the body-text colour on the new light BACKGROUND. The old
        // (220,220,220) was light-grey-on-dark and would be unreadable here.
        Level::INFO => theme::TEXT,
        Level::DEBUG => theme::TEXT_SECONDARY,
        Level::TRACE => theme::TEXT_DISABLED,
    }
}

pub fn show(
    ui: &mut egui::Ui,
    store: &LogStore,
    state: &mut LogPaneState,
    show_bundle_button: bool,
) {
    ui.horizontal(|ui| {
        if ui.button("Clear").clicked() {
            store.clear();
        }
        // Per app-shell spec R3: button is hidden when the AppState is
        // `Initializing` (no analysis state to dump, session_log_path
        // may not be bound yet).
        if show_bundle_button && ui.button("Download bug report…").clicked() {
            state.bundle_export_requested = true;
        }
        // `Save settings…` / `Load settings…` now live in the Data tab's
        // top toolbar row (`move-settings-buttons-to-data-tab`).
        ui.separator();
        ui.label(
            RichText::new(format!("filter: {}", state.filter_directive))
                .small()
                .color(theme::TEXT),
        );
    });
    ui.separator();

    let lines = store.snapshot();
    ScrollArea::vertical()
        .stick_to_bottom(true)
        .auto_shrink([false, false])
        .show(ui, |ui| {
            for line in &lines {
                let level_color = level_color(line.level);
                let ts = format_timestamp(line.timestamp);
                let header = format!(
                    "{ts} [{lvl}] {target}",
                    lvl = line.level,
                    target = line.target
                );
                ui.horizontal_wrapped(|ui| {
                    ui.label(
                        RichText::new(header)
                            .monospace()
                            .text_style(TextStyle::Monospace)
                            .color(level_color),
                    );
                    ui.label(
                        RichText::new(format!(" — {}", line.message))
                            .monospace()
                            .text_style(TextStyle::Monospace),
                    );
                });
            }
        });
}

#[cfg(test)]
mod tests {
    use super::format_timestamp;
    use std::time::SystemTime;

    #[test]
    fn format_timestamp_emits_local_time_with_offset() {
        let now = SystemTime::now();
        let rendered = format_timestamp(now);

        let re = regex::Regex::new(r"^\d{2}:\d{2}:\d{2}\.\d{3} [+-]\d{4}$").unwrap();
        assert!(
            re.is_match(&rendered),
            "format_timestamp output {rendered:?} does not match HH:MM:SS.mmm ±HHMM"
        );

        let expected = chrono::DateTime::<chrono::Local>::from(now)
            .format("%H:%M:%S%.3f %z")
            .to_string();
        assert_eq!(
            rendered, expected,
            "format_timestamp must agree with chrono::Local rendering for the same instant"
        );
    }
}
