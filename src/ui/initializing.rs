//! Startup splash renderer. Shown while the eager KEGG organism list
//! load is in flight (and during Retry attempts) per Track C of the
//! `add-kegg-module-ora` change.
//!
//! Behaviour:
//! - Loading (no `last_error`): centered spinner + "Initializing —
//!   fetching KEGG organism list…".
//! - Error (`last_error` set): centered error message + "Retry" button.
//!   If `fallback_cache` is also present, an additional "Use cached
//!   organisms (N days old)" button is shown.

use egui::{Align, Layout, RichText};

use crate::app::{App, AppState};
use crate::theme;

pub fn show(ui: &mut egui::Ui, app: &mut App) {
    // Pull out display state without holding a long borrow.
    let (loading, last_error, fallback_age_days) = match &app.state {
        AppState::Initializing {
            last_error,
            fallback_cache,
            ..
        } => {
            let age = fallback_cache.as_ref().map(|c| {
                let delta = chrono::Utc::now() - c.fetched_at;
                delta.num_days().max(0)
            });
            (last_error.is_none(), last_error.clone(), age)
        }
        _ => return,
    };

    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.with_layout(Layout::top_down(Align::Center), |ui| {
                ui.add_space(120.0);

                if loading {
                    ui.spinner();
                    ui.add_space(8.0);
                    ui.label(
                        RichText::new("Initializing — fetching KEGG organism list…")
                            .strong()
                            .size(16.0),
                    );
                    ui.add_space(4.0);
                    ui.label(
                        RichText::new("(first launch may take ~30 s; cached afterwards)")
                            .color(theme::TEXT)
                            .small(),
                    );
                } else if let Some(err) = last_error.as_ref() {
                    ui.label(
                        RichText::new("Couldn't load KEGG organism list")
                            .strong()
                            .size(16.0)
                            .color(theme::ERROR),
                    );
                    ui.add_space(6.0);
                    ui.label(RichText::new(err).color(theme::ERROR));
                    ui.add_space(12.0);

                    ui.horizontal(|ui| {
                        ui.add_space(80.0);
                        if ui.button("Retry").clicked() {
                            app.retry_organism_load();
                        }
                        if let Some(days) = fallback_age_days {
                            let label = if days == 0 {
                                "Use cached organisms (today)".to_string()
                            } else if days == 1 {
                                "Use cached organisms (1 day old)".to_string()
                            } else {
                                format!("Use cached organisms ({days} days old)")
                            };
                            if ui.button(label).clicked() {
                                app.accept_fallback_cache();
                            }
                        }
                    });
                }
            });
        });
}
