//! The bottom panel — a two-tab container (`add-bottom-panel-data-tab`).
//!
//! Renders a `Data` | `Log` tab strip at the top, then dispatches to the
//! stage-aware Data summary (`data_tab::show`, default-selected) or the Log
//! event stream (`log_pane::show`). The Data tab hosts the `Save settings…` /
//! `Load settings…` toolbar (`move-settings-buttons-to-data-tab`); the Log tab
//! keeps `[Clear]` + `[Download bug report…]`. The active tab lives on
//! `LogPaneState.active_tab`.

use crate::app::{App, AppState};
use crate::ui::widgets::{segmented_tab, segmented_track};
use crate::ui::{data_tab, log_pane};

use super::log_pane::BottomTab;

/// Render the bottom panel: tab strip + the active tab's body.
pub fn show(ui: &mut egui::Ui, app: &mut App) {
    // §4 Secondary "Segmented" tabs — Data/Log is card-internal navigation
    // (the page-level stepper is the Primary segmented set).
    //
    // To sit the track flush against the separator, the gap between two stacked
    // widgets is governed by `item_spacing.y` as read when the UPPER widget is
    // laid out. So zero it BEFORE the track (removing its trailing gap above the
    // rule), then restore it BEFORE the separator so the separator → body gap
    // below the rule stays normal.
    let prev_spacing = ui.spacing().item_spacing.y;
    ui.spacing_mut().item_spacing.y = 0.0;
    segmented_track(ui, |ui| {
        let active = app.log_ui.active_tab;
        if segmented_tab(ui, active == BottomTab::Data, "Data", true, false).clicked() {
            app.log_ui.active_tab = BottomTab::Data;
        }
        if segmented_tab(ui, active == BottomTab::Log, "Log", true, false).clicked() {
            app.log_ui.active_tab = BottomTab::Log;
        }
    });
    ui.spacing_mut().item_spacing.y = prev_spacing;
    ui.add(egui::Separator::default().spacing(2.0));

    match app.log_ui.active_tab {
        BottomTab::Data => data_tab::show(ui, app),
        BottomTab::Log => {
            // The bug-report button is hidden only during `Initializing`.
            // The `Save settings…` / `Load settings…` buttons have moved to
            // the Data tab's top toolbar (`move-settings-buttons-to-data-tab`),
            // so the Log tab no longer needs their gating booleans.
            let show_bundle_button = !matches!(app.state, AppState::Initializing { .. });
            log_pane::show(ui, &app.log, &mut app.log_ui, show_bundle_button);
        }
    }
}
