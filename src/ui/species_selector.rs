use egui::RichText;

use crate::kegg::KeggOrganism;
use crate::theme;

/// Per-frame UI state for the species selector (filter text, picker open flag).
#[derive(Debug, Default)]
pub struct SpeciesSelectorState {
    pub filter: String,
    pub picker_open: bool,
}

#[derive(Debug)]
pub enum SpeciesSelectorEvent {
    None,
    /// User selected this organism code.
    Selected(String),
    /// Picker was just opened and the organism list is not yet loaded — caller
    /// should kick off the lazy load.
    OpenedAndNeedsLoad,
}

/// Render the species selector.
///
/// - `organisms`: the loaded organism list, or `None` if not yet loaded.
/// - `loading`: true when an async load is in flight.
/// - `load_error`: human-readable error message from the last failed load.
/// - `current`: the currently selected organism code, if any.
/// - `enabled`: whether the picker button responds to clicks. While a KEGG
///   fetch is in flight we want the button locked, but the picker window (if
///   already open) should still be closeable.
pub fn show(
    ui: &mut egui::Ui,
    state: &mut SpeciesSelectorState,
    organisms: Option<&[KeggOrganism]>,
    loading: bool,
    load_error: Option<&str>,
    current: Option<&str>,
    enabled: bool,
) -> SpeciesSelectorEvent {
    let mut event = SpeciesSelectorEvent::None;

    let button_label = match current {
        Some(code) => match organisms.and_then(|os| os.iter().find(|o| o.code == code)) {
            Some(org) => format!("{} — {}", org.code, org.name),
            None => code.to_string(),
        },
        None => "Choose KEGG species".to_string(),
    };

    ui.horizontal(|ui| {
        // §2 Primary button — the KEGG-species selector is a core action.
        let button_response = crate::ui::widgets::primary_button_sized(
            ui,
            &button_label,
            enabled,
            egui::vec2(280.0, 24.0),
        );
        if button_response.clicked() {
            state.picker_open = !state.picker_open;
            if state.picker_open && organisms.is_none() && !loading {
                event = SpeciesSelectorEvent::OpenedAndNeedsLoad;
            }
        }
        // The pathway re-fetch action is the `[Refresh KEGG pathway cache]`
        // button rendered by the caller below the "Cached … ago" line, so it
        // sits adjacent to the cache-age it replaces (not an icon next to the
        // picker).
    });

    if state.picker_open {
        let ctx = ui.ctx().clone();
        let mut open = state.picker_open;
        egui::Window::new("Select KEGG species")
            .id(egui::Id::new("species_selector_window"))
            .open(&mut open)
            .resizable(true)
            .collapsible(false)
            .default_size([520.0, 520.0])
            .min_width(360.0)
            .min_height(220.0)
            .show(&ctx, |ui| {
                let body_event = picker_body(ui, state, organisms, loading, load_error, current);
                if !matches!(body_event, SpeciesSelectorEvent::None) {
                    event = body_event;
                }
                // Selecting an organism closes the window.
                if matches!(&event, SpeciesSelectorEvent::Selected(_)) {
                    // open will be flipped to false below.
                }
            });
        // Close the window when the user clicked an item.
        if matches!(&event, SpeciesSelectorEvent::Selected(_)) {
            open = false;
        }
        state.picker_open = open;
    }

    event
}

fn picker_body(
    ui: &mut egui::Ui,
    state: &mut SpeciesSelectorState,
    organisms: Option<&[KeggOrganism]>,
    loading: bool,
    load_error: Option<&str>,
    current: Option<&str>,
) -> SpeciesSelectorEvent {
    let mut event = SpeciesSelectorEvent::None;

    if loading {
        ui.horizontal(|ui| {
            ui.spinner();
            ui.label("Loading KEGG organism list…");
        });
        return event;
    }

    if let Some(err) = load_error {
        ui.colored_label(theme::ERROR, err);
        if ui.button("Retry").clicked() {
            event = SpeciesSelectorEvent::OpenedAndNeedsLoad;
        }
        return event;
    }

    let organisms = match organisms {
        Some(o) => o,
        None => {
            ui.label("(organism list not loaded)");
            return event;
        }
    };

    ui.horizontal(|ui| {
        ui.label("Search:");
        ui.add(egui::TextEdit::singleline(&mut state.filter).desired_width(f32::INFINITY));
    });

    let needle = state.filter.to_ascii_lowercase();
    let total = organisms.len();
    let filtered: Vec<&KeggOrganism> = if needle.is_empty() {
        organisms.iter().collect()
    } else {
        organisms
            .iter()
            .filter(|o| {
                o.code.to_ascii_lowercase().contains(&needle)
                    || o.name.to_ascii_lowercase().contains(&needle)
                    || o.lineage.to_ascii_lowercase().contains(&needle)
            })
            .collect()
    };

    ui.label(
        RichText::new(format!("{} / {} matches", filtered.len(), total))
            .small()
            .color(theme::TEXT),
    );

    ui.separator();

    // The ScrollArea fills the remaining vertical space inside the resizable
    // Window — drag the window's bottom-right corner to see more rows.
    egui::ScrollArea::vertical()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            // Selected row renders as a §3.3 Primary item: opaque PRIMARY fill
            // + white label (vs the translucent default text-selection halo).
            ui.style_mut().visuals.selection.bg_fill = theme::PRIMARY;
            const MAX_ROWS: usize = 1000;
            for org in filtered.iter().take(MAX_ROWS) {
                let is_current = Some(org.code.as_str()) == current;
                let label_text = format!("{} — {}", org.code, org.name);
                let label = if is_current {
                    RichText::new(label_text).color(theme::ON_PRIMARY)
                } else {
                    RichText::new(label_text)
                };
                if ui
                    .selectable_label(is_current, label)
                    .on_hover_text(&org.lineage)
                    .clicked()
                {
                    event = SpeciesSelectorEvent::Selected(org.code.clone());
                }
            }
            if filtered.len() > MAX_ROWS {
                ui.label(
                    RichText::new(format!(
                        "… {} more — refine your search",
                        filtered.len() - MAX_ROWS
                    ))
                    .small()
                    .color(theme::TEXT),
                );
            }
        });

    event
}
