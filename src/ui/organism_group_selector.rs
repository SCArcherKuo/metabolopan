//! Stage 1 widget for the module-mode "pick a taxonomy Level (1/2/3) +
//! pick a Group within that level" flow. Backed by the
//! `OrganismGroupIndex` precomputed at startup by `kegg::list_organisms`.
//!
//! Default Level = 2 (6 candidates: Animals, Archaea, Bacteria, Fungi,
//! Plants, Protists — the most useful biological partition). Level 1 has
//! 2 candidates (Prokaryotes/Eukaryotes); Level 3 has tens.

use egui::RichText;
use std::collections::HashSet;

use crate::theme;

use crate::kegg::OrganismGroupIndex;

/// Per-frame UI state for the organism-group selector.
#[derive(Debug)]
pub struct OrganismGroupSelectorState {
    /// Currently-active Level (1, 2, or 3). Defaults to 2.
    pub level: u8,
}

impl Default for OrganismGroupSelectorState {
    fn default() -> Self {
        Self { level: 2 }
    }
}

#[derive(Debug)]
pub enum OrganismGroupSelectorEvent {
    None,
    LevelChanged(u8),
    GroupSelected {
        level: u8,
        group: String,
        org_codes: HashSet<String>,
    },
}

/// Render the Level radio + Group dropdown. `index` is the precomputed
/// (level, group) → set-of-codes table built at startup; `current_group`
/// is the user's existing selection (highlighted in the dropdown). The
/// `enabled` flag locks the controls while a module fetch is in flight.
pub fn show(
    ui: &mut egui::Ui,
    state: &mut OrganismGroupSelectorState,
    index: Option<&OrganismGroupIndex>,
    current_group: Option<&str>,
    enabled: bool,
) -> OrganismGroupSelectorEvent {
    let mut event = OrganismGroupSelectorEvent::None;

    ui.label(
        RichText::new("Organism group")
            .strong()
            .color(theme::HEADING),
    );

    let Some(index) = index else {
        ui.label(
            RichText::new("(organism group index not loaded)")
                .color(theme::TEXT)
                .small(),
        );
        return event;
    };

    // ── Level radio ──
    ui.horizontal(|ui| {
        ui.label("Level:");
        for lvl in 1u8..=3 {
            let mut selected = state.level == lvl;
            if ui
                .add_enabled(enabled, egui::RadioButton::new(selected, format!("{lvl}")))
                .clicked()
                && !selected
            {
                selected = true;
                state.level = lvl;
                event = OrganismGroupSelectorEvent::LevelChanged(lvl);
            }
            let _ = selected;
        }
    });

    // ── Group dropdown ──
    let level_idx = (state.level.clamp(1, 3) - 1) as usize;
    let groups_map = &index.by_level[level_idx];
    let mut groups: Vec<(&String, &HashSet<String>)> = groups_map.iter().collect();
    groups.sort_by(|a, b| a.0.cmp(b.0));

    let button_label = match current_group {
        Some(g) => match groups_map.get(g) {
            Some(set) => format!("{} ({} organisms)", g, set.len()),
            None => g.to_string(),
        },
        None => format!(
            "Select a group (Level {} · {} candidates)",
            state.level,
            groups.len()
        ),
    };

    ui.horizontal(|ui| {
        // Primary dropdown — the Group selector is the core Module-mode choice.
        crate::ui::widgets::primary_dropdown(ui, |ui| {
            egui::ComboBox::from_id_salt("organism_group_combo")
                .selected_text(button_label)
                .show_ui(ui, |ui| {
                    if !enabled {
                        ui.disable();
                    }
                    for (group_name, codes) in &groups {
                        let label_text = format!("{} ({} organisms)", group_name, codes.len());
                        let is_current = current_group == Some(group_name.as_str());
                        // Selected item sits on the opaque PRIMARY fill — white text.
                        let label = if is_current {
                            egui::RichText::new(label_text).color(crate::theme::ON_PRIMARY)
                        } else {
                            egui::RichText::new(label_text)
                        };
                        // Re-emit `GroupSelected` even for the already-current
                        // Group (no `!is_current` suppression) so a restored
                        // Group whose cache is empty can be reloaded by clicking
                        // it — symmetric with the species selector. The Stage 3
                        // setup screen applies a skip-if-already-cached guard.
                        if ui
                            .add_enabled(enabled, egui::SelectableLabel::new(is_current, label))
                            .clicked()
                        {
                            event = OrganismGroupSelectorEvent::GroupSelected {
                                level: state.level,
                                group: (*group_name).clone(),
                                org_codes: (*codes).clone(),
                            };
                        }
                    }
                });
        });
    });

    event
}
