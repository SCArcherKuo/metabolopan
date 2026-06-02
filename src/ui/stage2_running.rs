use crate::app::{App, AppState};

pub fn show(ui: &mut egui::Ui, app: &App) {
    let AppState::Stage2DamRunning {
        mode_completed,
        mode_total,
        ..
    } = &app.state
    else {
        return;
    };
    let method = app.settings.dam_method;
    let ion_tables = app.inputs.ion_tables.as_slice();
    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.add_space(60.0);
            ui.vertical_centered(|ui| {
                ui.heading(egui::RichText::new("Running DAM").color(crate::theme::HEADING));
                ui.add_space(12.0);
                ui.label(format!("Running DAM ({})…", method.display_name()));
                ui.add_space(8.0);

                let n = ion_tables.len();
                for (idx, it) in ion_tables.iter().enumerate() {
                    let completed = mode_completed.get(idx).copied().unwrap_or(0);
                    let total = mode_total.get(idx).copied().unwrap_or(0);
                    let fraction = if total == 0 {
                        0.0
                    } else {
                        completed as f32 / total as f32
                    };
                    let label_text = if n >= 2 {
                        let mode = it.mode;
                        format!("{mode}: {completed} / {total} features")
                    } else {
                        format!("{completed} / {total} features")
                    };
                    crate::ui::widgets::progress_bar(
                        ui,
                        egui::ProgressBar::new(fraction)
                            .text(label_text)
                            .desired_width(380.0),
                        if fraction >= 1.0 {
                            crate::theme::SUCCESS
                        } else {
                            crate::theme::PRIMARY
                        },
                    );
                    if idx + 1 < n {
                        ui.add_space(4.0);
                    }
                }
            });
        });
}
