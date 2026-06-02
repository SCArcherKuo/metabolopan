//! Stage 3 running screen — 3-phase progress (PubChem → KEGG conv → ORA).

use crate::app::{App, AppState, Stage3Phase};
use crate::theme;

pub fn show(ui: &mut egui::Ui, app: &App) {
    let AppState::Stage3EnrichRunning {
        phase,
        pubchem_completed,
        pubchem_total,
        kegg_conv_completed,
        kegg_conv_total,
        ..
    } = &app.state
    else {
        return;
    };

    let mode_suffix = match app.settings.analysis_mode {
        crate::app::AnalysisMode::Pathway => "Pathway mode",
        crate::app::AnalysisMode::Module => "Module mode",
    };

    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.add_space(60.0);
            ui.vertical_centered(|ui| match phase {
                Stage3Phase::PubChem => {
                    ui.heading(
                        egui::RichText::new(format!(
                            "Phase 1/3: Resolving PubChem CIDs… · {mode_suffix}"
                        ))
                        .color(theme::HEADING),
                    );
                    ui.add_space(12.0);
                    let fraction = if *pubchem_total == 0 {
                        0.0
                    } else {
                        *pubchem_completed as f32 / *pubchem_total as f32
                    };
                    crate::ui::widgets::progress_bar(
                        ui,
                        egui::ProgressBar::new(fraction)
                            .text(format!("{pubchem_completed} / {pubchem_total} InChIKeys"))
                            .desired_width(420.0),
                        if fraction >= 1.0 {
                            crate::theme::SUCCESS
                        } else {
                            crate::theme::PRIMARY
                        },
                    );
                }
                Stage3Phase::KeggConv => {
                    ui.heading(
                        egui::RichText::new(format!(
                            "Phase 2/3: Resolving KEGG compound IDs… · {mode_suffix}"
                        ))
                        .color(theme::HEADING),
                    );
                    ui.add_space(12.0);
                    let fraction = if *kegg_conv_total == 0 {
                        0.0
                    } else {
                        *kegg_conv_completed as f32 / *kegg_conv_total as f32
                    };
                    crate::ui::widgets::progress_bar(
                        ui,
                        egui::ProgressBar::new(fraction)
                            .text(format!("{kegg_conv_completed} / {kegg_conv_total} CIDs"))
                            .desired_width(420.0),
                        if fraction >= 1.0 {
                            crate::theme::SUCCESS
                        } else {
                            crate::theme::PRIMARY
                        },
                    );
                }
                Stage3Phase::Ora => {
                    ui.heading(
                        egui::RichText::new(format!("Phase 3/3: Running ORA… · {mode_suffix}"))
                            .color(theme::HEADING),
                    );
                    ui.add_space(12.0);
                    ui.spinner();
                }
            });
        });
}
