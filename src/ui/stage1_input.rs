use egui::{Color32, RichText};
use std::path::PathBuf;
use tracing::{error, info, warn};

use crate::app::{App, AppState};
use crate::data::GroupMapping;
use crate::data::groups::{UNASSIGNED, load_group_mapping};
use crate::data::msdial::parse_msdial_txt;
use crate::data::{AdductPolarityInference, IonMode, IonModeTable, infer_polarity};
use crate::theme;
use crate::ui::widgets::{file_pick_button, primary_button};

/// Errors that block advancement from Stage 1 to Stage 2. After
/// `reorder-gui-and-move-mode-to-stage3`, Stage 1 is mode-agnostic:
/// Analysis Mode + species/group selector + their related gates live on
/// Stage 3 setup, so this struct only carries universal + dual-mode
/// integrity fields.
#[derive(Default)]
pub struct Stage1ValidationInput<'a> {
    pub table_loaded: bool,
    pub slot1_sample_cols: &'a [String],
    pub slot2_sample_cols: Option<&'a [String]>,
    pub mapping: Option<&'a GroupMapping>,
    pub slot1_mode: Option<IonMode>,
    pub slot2_revealed: bool,
    pub slot2_mode: Option<IonMode>,
}

pub fn validate_for_dam(input: Stage1ValidationInput<'_>) -> Result<(), Vec<String>> {
    let mut issues: Vec<String> = Vec::new();

    if !input.table_loaded {
        // No MS-DIAL file yet — keep "Continue to DAM" disabled, but render no
        // nag (the file pickers above are self-explanatory). An empty `Err`
        // shows nothing while still failing the `matches!(Ok)` gate.
        return Err(issues);
    }

    let mapping = match input.mapping {
        Some(m) => m,
        None => {
            // No metadata CSV picked yet — keep "Continue to DAM" disabled but
            // render NO nag (the `.csv` picker above is self-explanatory,
            // mirroring the no-`.txt` branch). Always fail the gate: returning
            // `Ok` here would wrongly enable Continue in single-mode. A `.csv`
            // that was picked but FAILED to parse surfaces its red error
            // elsewhere; this branch only covers "not provided yet".
            check_dual_mode_rules(&input, &mut issues);
            return Err(issues);
        }
    };

    if mapping.assigned_count() == 0 {
        issues.push(
            "No samples in the metadata match the MS-DIAL .txt. Check the `sample` column."
                .to_string(),
        );
    }

    let assignable: Vec<String> = mapping
        .groups()
        .into_iter()
        .filter(|g| g != UNASSIGNED)
        .collect();

    if assignable.len() < 2 {
        let found = if assignable.is_empty() {
            "none".to_string()
        } else {
            assignable.join(", ")
        };
        issues.push(format!(
            "At least 2 groups are needed for DAM. Found: {found}."
        ));
    }

    for group in &assignable {
        let n = mapping.samples_in(group).len();
        if n < 2 {
            issues.push(format!(
                "Group `{group}` has {n} sample(s). At least 2 samples per group are required."
            ));
        }
    }

    check_dual_mode_rules(&input, &mut issues);

    if issues.is_empty() {
        Ok(())
    } else {
        Err(issues)
    }
}

fn check_dual_mode_rules(input: &Stage1ValidationInput<'_>, issues: &mut Vec<String>) {
    if !input.table_loaded {
        return;
    }

    if input.slot1_mode.is_none() {
        issues.push("Pick the ionization mode for the first .txt.".to_string());
    }

    if !input.slot2_revealed {
        return;
    }

    let slot2_loaded = input.slot2_sample_cols.is_some();

    if !slot2_loaded {
        issues.push("Upload the second MS-DIAL .txt or remove slot #2.".to_string());
    }

    if slot2_loaded && input.slot2_mode.is_none() {
        issues.push("Pick the ionization mode for the second .txt.".to_string());
    }

    if let (Some(m1), Some(m2)) = (input.slot1_mode, input.slot2_mode)
        && m1 == m2
    {
        issues.push(format!(
            "POS and NEG must be different — both slots are set to {m1}."
        ));
    }

    let Some(mapping) = input.mapping else {
        return;
    };

    if slot2_loaded && !mapping.has_biosample() {
        issues.push(
            "Dual-mode requires a 'biosample' column in the metadata CSV. Add it or remove the second .txt file.".to_string(),
        );
    }

    let Some(slot2_cols) = input.slot2_sample_cols else {
        return;
    };
    let (Some(m1), Some(m2)) = (input.slot1_mode, input.slot2_mode) else {
        return;
    };
    if m1 == m2 {
        return;
    }

    let slot1_label = m1.to_string();
    let slot2_label = m2.to_string();

    let mut per_group: std::collections::BTreeMap<String, (usize, usize)> =
        std::collections::BTreeMap::new();
    for s in input.slot1_sample_cols {
        let g = mapping.group_of(s);
        if g != UNASSIGNED {
            per_group.entry(g.to_string()).or_insert((0, 0)).0 += 1;
        }
    }
    for s in slot2_cols {
        let g = mapping.group_of(s);
        if g != UNASSIGNED {
            per_group.entry(g.to_string()).or_insert((0, 0)).1 += 1;
        }
    }
    for (g, (n1, n2)) in &per_group {
        if *n1 == 0 && *n2 == 0 {
            continue;
        }
        if *n1 < 2 || *n2 < 2 {
            issues.push(format!(
                "Group '{g}' has {n1} sample(s) in {slot1_label} but {n2} in {slot2_label} — both modes need ≥ 2."
            ));
        }
    }

    let dup_check = |samples: &[String], label: &str, issues: &mut Vec<String>| {
        let mut counts: std::collections::BTreeMap<String, usize> =
            std::collections::BTreeMap::new();
        for s in samples {
            if let Some(bio) = mapping.biosample_of(s) {
                *counts.entry(bio.to_string()).or_insert(0) += 1;
            }
        }
        for (bio, n) in counts.iter().filter(|(_, n)| **n > 1) {
            issues.push(format!(
                "Biosample '{bio}' appears in {n} {label} rows — must be unique per mode."
            ));
        }
    };
    dup_check(input.slot1_sample_cols, &slot1_label, issues);
    dup_check(slot2_cols, &slot2_label, issues);

    let mut bio_groups: std::collections::BTreeMap<String, Vec<(&str, String)>> =
        std::collections::BTreeMap::new();
    for s in input.slot1_sample_cols {
        if let Some(bio) = mapping.biosample_of(s) {
            let g = mapping.group_of(s);
            if g != UNASSIGNED {
                bio_groups
                    .entry(bio.to_string())
                    .or_default()
                    .push((slot1_label.as_str(), g.to_string()));
            }
        }
    }
    for s in slot2_cols {
        if let Some(bio) = mapping.biosample_of(s) {
            let g = mapping.group_of(s);
            if g != UNASSIGNED {
                bio_groups
                    .entry(bio.to_string())
                    .or_default()
                    .push((slot2_label.as_str(), g.to_string()));
            }
        }
    }
    for (bio, entries) in &bio_groups {
        if entries.len() < 2 {
            continue;
        }
        let unique_groups: std::collections::BTreeSet<&str> =
            entries.iter().map(|(_, g)| g.as_str()).collect();
        if unique_groups.len() > 1 {
            let pos_g = entries
                .iter()
                .find(|(l, _)| *l == slot1_label)
                .map(|(_, g)| g.clone())
                .unwrap_or_default();
            let neg_g = entries
                .iter()
                .find(|(l, _)| *l == slot2_label)
                .map(|(_, g)| g.clone())
                .unwrap_or_default();
            issues.push(format!(
                "Biosample '{bio}' is in group '{pos_g}' in {slot1_label} but '{neg_g}' in {slot2_label}."
            ));
        }
    }
}

/// Slot #1 auto-fill helper (D1): maps the `infer_polarity` result to the value
/// to write into `slot1_mode` on a fresh `.txt` load. `Ambiguous` returns `None`
/// so the existing "Could not auto-detect…" grey hint stays.
fn decide_slot1_mode_on_file_load(inferred: AdductPolarityInference) -> Option<IonMode> {
    match inferred {
        AdductPolarityInference::Positive => Some(IonMode::Positive),
        AdductPolarityInference::Negative => Some(IonMode::Negative),
        AdductPolarityInference::Ambiguous => None,
    }
}

/// Slot #2 auto-fill helper (D3 trigger #3 + D6 "write-always" convention).
///
/// Returns the value the caller MUST unconditionally write into `slot2_mode`.
/// No-op cases (slot 2 not revealed, slot 1 is `None`, slot 2 already has the
/// correct opposite of slot 1) return `slot2_before` verbatim so the caller's
/// `*slot2_mode = helper(...)` is a self-assignment.
fn decide_slot2_mode_on_slot1_change(
    slot1_new: Option<IonMode>,
    slot2_revealed: bool,
    slot2_before: Option<IonMode>,
) -> Option<IonMode> {
    if !slot2_revealed {
        return slot2_before;
    }
    let Some(s1) = slot1_new else {
        return slot2_before;
    };
    let opposite = s1.opposite();
    match slot2_before {
        None => Some(opposite),
        Some(s2) if s2 == s1 => Some(opposite),
        Some(_) => slot2_before,
    }
}

pub fn show(ui: &mut egui::Ui, app: &mut App) {
    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            ui.heading(egui::RichText::new("Stage 1 — Input").color(theme::HEADING));
            ui.add_space(8.0);

            // === MS-DIAL .txt picker(s) ===
            render_txt_picker(ui, app);

            ui.add_space(8.0);

            // === group .csv picker ===
            render_csv_picker(ui, app);

            // The CSV→sample coverage line and the per-group sample / biosample
            // list moved to the bottom-panel Data tab (`data-summary-panel`).
            // The Stage 1 body keeps only the pickers + Continue gate.

            // === CSV-only samples banner (persistent, non-blocking) ===
            // Metadata rows whose `sample` matched no .txt column are ignored
            // (and WARN-logged) but retained on the mapping; surface them loudly
            // here so a casing/underscore typo isn't silently dropped (turning a
            // 6-vs-6 comparison into 6-vs-5). Distinct from Unassigned (.txt
            // columns with no CSV row, shown yellow in the Data tab). Does NOT
            // gate Continue — see the `stage1-ui` capability spec.
            if let Some(mapping) = app.inputs.mapping.as_ref() {
                let unmatched = mapping.unmatched_csv_samples();
                if !unmatched.is_empty() {
                    ui.add_space(8.0);
                    egui::Frame::group(ui.style())
                        .stroke(egui::Stroke::new(1.0, theme::ERROR))
                        .show(ui, |ui| {
                            ui.colored_label(
                                theme::ERROR,
                                RichText::new(format!(
                                    "{} metadata row(s) name samples not found in the MS-DIAL .txt — these rows are ignored:",
                                    unmatched.len()
                                ))
                                .strong(),
                            );
                            ui.colored_label(theme::ERROR, unmatched.join(", "));
                            ui.colored_label(
                                theme::ERROR,
                                "Fix the `sample` column in your .csv (check casing / underscores) or remove these rows.",
                            );
                        });
                }
            }

            ui.add_space(12.0);

            // === Error block ===
            if let AppState::Stage1Input { error: Some(e), .. } = &app.state {
                ui.colored_label(theme::ERROR, e.clone());
            }

            // === Validation + Continue button ===
            //
            // Stage 1 is mode-agnostic after `reorder-gui-and-move-mode-to-stage3`:
            // Analysis Mode + species/group selector live on Stage 3 setup. The
            // gate only enforces universal + dual-mode integrity checks.
            let validation_result = if let AppState::Stage1Input {
                slot1_mode,
                slot2_revealed,
                slot2_mode,
                ..
            } = &app.state
            {
                let ion_tables = app.inputs.ion_tables.as_slice();
                let slot1_sample_cols: &[String] = ion_tables
                    .first()
                    .map(|it| it.table.sample_cols.as_slice())
                    .unwrap_or(&[]);
                let slot2_sample_cols: Option<&[String]> =
                    ion_tables.get(1).map(|it| it.table.sample_cols.as_slice());
                Some(validate_for_dam(Stage1ValidationInput {
                    table_loaded: !ion_tables.is_empty(),
                    slot1_sample_cols,
                    slot2_sample_cols,
                    mapping: app.inputs.mapping.as_ref(),
                    slot1_mode: *slot1_mode,
                    slot2_revealed: *slot2_revealed,
                    slot2_mode: *slot2_mode,
                }))
            } else {
                None
            };

            if let Some(Err(ref issues)) = validation_result {
                for issue in issues {
                    ui.colored_label(theme::ERROR, issue);
                }
            }

            let can_start = matches!(validation_result, Some(Ok(())));
            let start_button = primary_button(ui, "Continue to DAM", can_start);
            if start_button.clicked() && can_start {
                promote_to_stage2(app);
            }
        });
}

fn render_txt_picker(ui: &mut egui::Ui, app: &mut App) {
    render_slot(ui, app, 0);

    let (slot1_loaded, slot2_revealed) = match &app.state {
        AppState::Stage1Input { slot2_revealed, .. } => {
            (!app.inputs.ion_tables.is_empty(), *slot2_revealed)
        }
        _ => (false, false),
    };

    if slot1_loaded && !slot2_revealed {
        ui.add_space(4.0);
        if ui.button("+ Add second ionization mode").clicked()
            && let AppState::Stage1Input {
                slot1_mode,
                slot2_revealed,
                slot2_mode,
                ..
            } = &mut app.state
        {
            *slot2_revealed = true;
            // Trigger #1 (D3): auto-fill slot 2 to opposite of slot 1 on reveal.
            // Intentionally inlined per D3 — the logic is a single conditional.
            if let Some(s1) = *slot1_mode {
                *slot2_mode = Some(s1.opposite());
            }
        }
    } else if slot2_revealed {
        ui.add_space(8.0);
        render_slot(ui, app, 1);
    }
}

fn render_slot(ui: &mut egui::Ui, app: &mut App, slot_idx: usize) {
    // ── Snapshot reads ──
    let it = app.inputs.ion_tables.get(slot_idx);
    let loaded_path = it.and_then(|x| x.txt_path.clone());
    let (current_mode_radio, disabled_mode) = match &app.state {
        AppState::Stage1Input {
            slot1_mode,
            slot2_mode,
            ..
        } => {
            let mode_radio = if slot_idx == 0 {
                *slot1_mode
            } else {
                *slot2_mode
            };
            let disabled = if slot_idx == 1 { *slot1_mode } else { None };
            (mode_radio, disabled)
        }
        _ => return,
    };
    let detected = it.map(|x| infer_polarity(&x.table));

    let mut picked: Option<PathBuf> = None;
    let mut remove_clicked = false;
    let mut chosen_mode = current_mode_radio;

    // ── Row 1: picker + filename + (slot #2) remove ──
    ui.horizontal(|ui| {
        let label = if slot_idx == 0 {
            "Choose MS-DIAL .txt"
        } else {
            "Choose MS-DIAL .txt (slot #2)"
        };
        picked = file_pick_button(ui, label, "MS-DIAL txt", "txt", true);
        if let Some(p) = &loaded_path {
            display_path_label(ui, p);
        }
        if slot_idx == 1
            && ui
                .button("×")
                .on_hover_text("Remove slot #2 and return to single-mode")
                .clicked()
        {
            remove_clicked = true;
        }
    });

    // ── Row 2: mode radio ──
    ui.horizontal(|ui| {
        ui.label("Mode:");
        let pos_enabled = disabled_mode != Some(IonMode::Positive);
        let neg_enabled = disabled_mode != Some(IonMode::Negative);
        let pos_resp = ui.add_enabled(
            pos_enabled,
            egui::RadioButton::new(chosen_mode == Some(IonMode::Positive), "Positive"),
        );
        if pos_resp.clicked() && pos_enabled {
            chosen_mode = Some(IonMode::Positive);
        }
        if !pos_enabled {
            pos_resp.on_disabled_hover_text("Mode already assigned to slot #1");
        }
        let neg_resp = ui.add_enabled(
            neg_enabled,
            egui::RadioButton::new(chosen_mode == Some(IonMode::Negative), "Negative"),
        );
        if neg_resp.clicked() && neg_enabled {
            chosen_mode = Some(IonMode::Negative);
        }
        if !neg_enabled {
            neg_resp.on_disabled_hover_text("Mode already assigned to slot #1");
        }
    });

    // Per-slot feature / sample-column counts and the excluded computed-stats
    // columns now live in the bottom-panel Data tab (`data-summary-panel`), not
    // inline here — the Stage 1 body shows only the picker + mode radio.

    if let Some(d) = detected {
        render_adduct_hint(ui, chosen_mode, d);
    }

    // ── Apply: mode change ──
    if chosen_mode != current_mode_radio {
        if let AppState::Stage1Input {
            slot1_mode,
            slot2_revealed,
            slot2_mode,
            ..
        } = &mut app.state
        {
            if slot_idx == 0 {
                *slot1_mode = chosen_mode;
                // D3 trigger #3: slot 1 mode changed → auto-fill / flip slot 2
                // if it's unset or conflicts. Helper returns slot2_mode verbatim
                // for no-op cases (write-always convention per D6).
                *slot2_mode =
                    decide_slot2_mode_on_slot1_change(*slot1_mode, *slot2_revealed, *slot2_mode);
            } else {
                *slot2_mode = chosen_mode;
            }
        }
        if let Some(m) = chosen_mode
            && let Some(it) = app.inputs.ion_tables.get_mut(slot_idx)
        {
            it.mode = m;
        }
        // D3 trigger #3 follow-up: when slot 0 mode changes, slot 2's
        // `IonModeTable.mode` (if loaded) must follow the new slot2_mode
        // so the data-model stays consistent with the radio.
        if slot_idx == 0
            && let AppState::Stage1Input { slot2_mode: s2, .. } = &app.state
            && let Some(new_s2_mode) = *s2
            && let Some(it1) = app.inputs.ion_tables.get_mut(1)
        {
            it1.mode = new_s2_mode;
        }
    }

    // ── Apply: × remove (slot #2 only) ──
    if remove_clicked {
        let mut shrunk = false;
        if app.inputs.ion_tables.len() > 1 {
            app.inputs.ion_tables.truncate(1);
            shrunk = true;
        }
        if let AppState::Stage1Input {
            slot2_revealed,
            slot2_mode,
            ..
        } = &mut app.state
        {
            *slot2_revealed = false;
            *slot2_mode = None;
        }
        if shrunk {
            // Mapping was loaded against the union of both modes' sample
            // columns; rebuild against the now-single-mode sample list.
            reload_mapping_after_ion_tables_change(app);
        }
    }

    // ── Apply: file pick ──
    if let Some(picked) = picked {
        let was_slot2_load = slot_idx == 1;
        apply_picked_file_to_slot(app, slot_idx, picked);
        if was_slot2_load {
            reload_mapping_after_ion_tables_change(app);
        }
    }
}

/// Parse the picked `.txt` and stash it into `app.inputs.ion_tables[slot_idx]`.
/// Re-picking slot #1 invalidates the CSV mapping and the species KEGG cache.
///
/// On parse success, runs `infer_polarity` and applies the auto-fill rules
/// (auto-infer-stage1-ion-mode §2.1 / §2.2):
/// - Slot 0: writes `slot1_mode` from the inferred polarity (D1); `IonModeTable.mode`
///   falls back to `Positive` when Ambiguous, preserving pre-change placeholder behavior.
/// - Slot 1: if slot 1 has a mode (auto-filled or manual), writes `slot2_mode` to
///   `Some(opposite(slot1))` (D3 trigger #2); `IonModeTable.mode` follows.
///   When slot 1 is unset, `slot2_mode` is left alone and `IonModeTable.mode`
///   defaults to `Negative` (literal pre-change fallback per design D6 H3 lock-in).
fn apply_picked_file_to_slot(app: &mut App, slot_idx: usize, picked: PathBuf) {
    match parse_msdial_txt(&picked) {
        Ok(t) => {
            info!(
                path = %picked.display(),
                slot = slot_idx,
                features = t.features.len(),
                samples = t.sample_cols.len(),
                "loaded MS-DIAL .txt"
            );
            // ── Auto-infer + write the radio state ──
            // For slot 0 we re-stamp slot1_mode from infer_polarity (overwriting any
            // prior value — re-picking implies the prior value referred to a now-gone file
            // per D2). For slot 1 we mirror the opposite of slot1_mode into slot2_mode
            // (D3 trigger #2), overwriting any prior slot 2 value.
            let inferred = infer_polarity(&t);
            let prior_slot1_mode = match &app.state {
                AppState::Stage1Input { slot1_mode, .. } => *slot1_mode,
                _ => None,
            };
            let new_slot1_mode = if slot_idx == 0 {
                decide_slot1_mode_on_file_load(inferred)
            } else {
                prior_slot1_mode
            };
            if slot_idx == 0
                && let AppState::Stage1Input {
                    slot1_mode,
                    slot2_revealed,
                    slot2_mode,
                    ..
                } = &mut app.state
            {
                *slot1_mode = new_slot1_mode;
                // D3 trigger #3 on the auto-fill path: slot 1's mode just changed
                // via re-pick / file load (per spec line 15: "via auto-fill on
                // file load OR via the user manually clicking the radio").
                // Without this call, a re-pick that flips slot 1's polarity
                // would leave slot 2 conflicting (e.g. POS slot1 → reveal slot 2
                // (auto NEG) → re-pick slot 1 with NEG file → both Negative).
                *slot2_mode =
                    decide_slot2_mode_on_slot1_change(*slot1_mode, *slot2_revealed, *slot2_mode);
            }
            if slot_idx == 1
                && let Some(s1) = new_slot1_mode
                && let AppState::Stage1Input { slot2_mode, .. } = &mut app.state
            {
                *slot2_mode = Some(s1.opposite());
            }

            // Placeholder mode for the non-Option `IonModeTable.mode` field.
            // Per §2.5: the block is retained (not deleted) because the field
            // requires *some* value; the fallback rules now derive from the
            // auto-fill outcomes above.
            //   - slot 0: unwrap_or(Positive) — fires only when Ambiguous.
            //   - slot 1: opposite-of-slot1 if slot1 is set; else Negative (pre-change literal).
            let placeholder_mode = if slot_idx == 0 {
                new_slot1_mode.unwrap_or(IonMode::Positive)
            } else {
                new_slot1_mode
                    .map(IonMode::opposite)
                    .unwrap_or(IonMode::Negative)
            };
            let new_table = IonModeTable {
                mode: placeholder_mode,
                table: t,
                txt_path: Some(picked),
            };
            if slot_idx >= app.inputs.ion_tables.len() {
                app.inputs.ion_tables.push(new_table);
            } else {
                app.inputs.ion_tables[slot_idx] = new_table;
            }
            // Slot 1 IonModeTable.mode sync after a slot 0 re-pick that may
            // have flipped slot2_mode via the D3 trigger #3 call above. Mirrors
            // the render_slot manual-click path so data-model stays consistent
            // with the radio for both auto-fill and manual transitions.
            if slot_idx == 0
                && let AppState::Stage1Input { slot2_mode: s2, .. } = &app.state
                && let Some(new_s2_mode) = *s2
                && let Some(it1) = app.inputs.ion_tables.get_mut(1)
            {
                it1.mode = new_s2_mode;
            }
            if slot_idx == 0 {
                app.inputs.mapping = None;
                app.inputs.csv_path = None;
                app.cache.species_kegg = None;
            }
            if let AppState::Stage1Input { error, .. } = &mut app.state {
                *error = None;
            }
        }
        Err(e) => {
            error!(path = %picked.display(), slot = slot_idx, error = %e, "failed to parse MS-DIAL .txt");
            let err_msg = format!("Failed to parse .txt ({}): {e}", picked.display());
            if slot_idx == 0 {
                app.inputs.ion_tables.clear();
                app.inputs.mapping = None;
                app.inputs.csv_path = None;
                app.cache.species_kegg = None;
                if let AppState::Stage1Input {
                    slot1_mode,
                    slot2_revealed,
                    slot2_mode,
                    error,
                    ..
                } = &mut app.state
                {
                    *slot1_mode = None;
                    *slot2_revealed = false;
                    *slot2_mode = None;
                    *error = Some(err_msg);
                }
            } else {
                app.inputs.ion_tables.truncate(1);
                if let AppState::Stage1Input { error, .. } = &mut app.state {
                    *error = Some(err_msg);
                }
            }
        }
    }
}

// The (radio, detect) → hint contract is a pure function of its two inputs
// regardless of how the radio reached its current value. After
// `auto-infer-stage1-ion-mode`, the (Unset, Positive) and (Unset, Negative)
// arms are unreachable on the happy path for slot #1 (auto-fill commits
// `slot1_mode = Some(X)` in the same frame as parse-success), so the only
// `Unset` arm that fires regularly is (Unset, Ambiguous). The other two stay
// as defensive backstops — don't delete them.
fn render_adduct_hint(
    ui: &mut egui::Ui,
    radio: Option<IonMode>,
    detected: AdductPolarityInference,
) {
    use AdductPolarityInference::*;
    let (text, color): (Option<String>, Color32) = match (radio, detected) {
        (None, Positive) => (
            Some("Detected: Positive (from Adduct type column). Pick a mode to continue.".into()),
            theme::TEXT,
        ),
        (None, Negative) => (
            Some("Detected: Negative (from Adduct type column). Pick a mode to continue.".into()),
            theme::TEXT,
        ),
        (None, Ambiguous) => (
            Some(
                "Could not auto-detect mode from Adduct type column. Pick a mode to continue."
                    .into(),
            ),
            theme::TEXT,
        ),
        (Some(IonMode::Positive), Negative) => (
            Some("Adduct column says Negative but you selected Positive. Please confirm.".into()),
            theme::WARNING,
        ),
        (Some(IonMode::Negative), Positive) => (
            Some("Adduct column says Positive but you selected Negative. Please confirm.".into()),
            theme::WARNING,
        ),
        _ => (None, theme::TEXT),
    };
    if let Some(t) = text {
        ui.label(RichText::new(t).small().color(color));
    }
}

/// Union of all loaded `IonModeTable.table.sample_cols`, deduplicated.
fn union_sample_cols(ion_tables: &[IonModeTable]) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    let mut seen: std::collections::HashSet<String> = std::collections::HashSet::new();
    for it in ion_tables {
        for s in &it.table.sample_cols {
            if seen.insert(s.clone()) {
                out.push(s.clone());
            }
        }
    }
    out
}

/// Re-load the mapping at `app.inputs.csv_path` against the current
/// `app.inputs.ion_tables` union.
fn reload_mapping_after_ion_tables_change(app: &mut App) {
    let Some(csv_path) = app.inputs.csv_path.clone() else {
        return;
    };
    if app.inputs.ion_tables.is_empty() {
        return;
    }
    let sample_cols = union_sample_cols(&app.inputs.ion_tables);
    match load_group_mapping(&csv_path, &sample_cols) {
        Ok(m) => {
            info!(
                path = %csv_path.display(),
                assigned = m.assigned_count(),
                groups = m.groups().len(),
                "reloaded group mapping after ion_tables change"
            );
            app.inputs.mapping = Some(m);
            if let AppState::Stage1Input { error, .. } = &mut app.state {
                *error = None;
            }
        }
        Err(e) => {
            error!(path = %csv_path.display(), error = %e, "failed to reload group .csv after ion_tables change");
            app.inputs.mapping = None;
            if let AppState::Stage1Input { error, .. } = &mut app.state {
                *error = Some(format!("Failed to parse .csv: {e}"));
            }
        }
    }
}

fn render_csv_picker(ui: &mut egui::Ui, app: &mut App) {
    let csv_enabled = !app.inputs.ion_tables.is_empty();
    let displayed_path = app.inputs.csv_path.clone();
    let mut picked: Option<PathBuf> = None;
    ui.horizontal(|ui| {
        picked = file_pick_button(ui, "Choose group .csv", "metadata csv", "csv", csv_enabled);
        if let Some(p) = &displayed_path {
            display_path_label(ui, p);
        }
    });

    if let Some(picked) = picked {
        let sample_cols = union_sample_cols(&app.inputs.ion_tables);
        match load_group_mapping(&picked, &sample_cols) {
            Ok(m) => {
                info!(
                    path = %picked.display(),
                    assigned = m.assigned_count(),
                    groups = m.groups().len(),
                    "loaded group mapping"
                );
                app.inputs.mapping = Some(m);
                app.inputs.csv_path = Some(picked);
                if let AppState::Stage1Input { error, .. } = &mut app.state {
                    *error = None;
                }
            }
            Err(e) => {
                error!(path = %picked.display(), error = %e, "failed to parse group .csv");
                app.inputs.mapping = None;
                app.inputs.csv_path = Some(picked);
                if let AppState::Stage1Input { error, .. } = &mut app.state {
                    *error = Some(format!("Failed to parse .csv: {e}"));
                }
            }
        }
    }
}

fn promote_to_stage2(app: &mut App) {
    let prev = std::mem::take(&mut app.state);
    if !matches!(prev, AppState::Stage1Input { .. }) {
        app.state = prev;
        return;
    }

    // Run the data-model backstop on inputs.ion_tables: enforce length 1–2,
    // no-duplicate-mode, and canonical Positive-first ordering. Validation
    // in `validate_for_dam` already enforces these, so any failure here is
    // unreachable — fall back to default state if it fires.
    let ion_tables = std::mem::take(&mut app.inputs.ion_tables);
    match crate::data::IonModeTables::try_new(ion_tables) {
        Ok(v) => {
            app.inputs.ion_tables = v.into_inner();
        }
        Err(e) => {
            warn!(
                error = %e,
                "promote_to_stage2: ion_tables failed IonModeTables::try_new (validation should have prevented this); restoring state"
            );
            app.state = AppState::default();
            return;
        }
    }

    info!(mode = ?app.settings.analysis_mode, "transitioning to Stage 2");
    app.state = AppState::Stage2DamSetup { error: None };
}

fn display_path_label(ui: &mut egui::Ui, p: &std::path::Path) {
    let basename = p
        .file_name()
        .map(|s| s.to_string_lossy().to_string())
        .unwrap_or_else(|| p.display().to_string());
    // Picked-file basename in success green (初葉綠) so a clean load reads at a
    // glance; full path still on hover. Shared by all three pickers.
    let response = ui.colored_label(theme::SUCCESS, basename);
    response.on_hover_text(p.display().to_string());
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    use tempfile::NamedTempFile;

    fn write_csv(content: &str) -> NamedTempFile {
        let mut f = NamedTempFile::new().expect("create tempfile");
        f.write_all(content.as_bytes()).expect("write fixture");
        f
    }

    fn cols(names: &[&str]) -> Vec<String> {
        names.iter().map(|s| s.to_string()).collect()
    }

    fn balanced_dual_mode_csv() -> NamedTempFile {
        write_csv(
            "sample,biosample,group\n\
             P1,BIO-1,ctrl\n\
             P2,BIO-2,ctrl\n\
             P3,BIO-3,treat\n\
             P4,BIO-4,treat\n\
             N1,BIO-1,ctrl\n\
             N2,BIO-2,ctrl\n\
             N3,BIO-3,treat\n\
             N4,BIO-4,treat\n",
        )
    }

    fn pos_cols() -> Vec<String> {
        cols(&["P1", "P2", "P3", "P4"])
    }
    fn neg_cols() -> Vec<String> {
        cols(&["N1", "N2", "N3", "N4"])
    }

    fn balanced_mapping() -> GroupMapping {
        let f = balanced_dual_mode_csv();
        let union: Vec<String> = pos_cols().into_iter().chain(neg_cols()).collect();
        load_group_mapping(f.path(), &union).expect("balanced fixture mapping loads")
    }

    fn dual_mode_input<'a>(
        mapping: &'a GroupMapping,
        slot1: &'a [String],
        slot2: &'a [String],
    ) -> Stage1ValidationInput<'a> {
        Stage1ValidationInput {
            table_loaded: true,
            slot1_sample_cols: slot1,
            slot2_sample_cols: Some(slot2),
            mapping: Some(mapping),
            slot1_mode: Some(IonMode::Positive),
            slot2_revealed: true,
            slot2_mode: Some(IonMode::Negative),
        }
    }

    #[test]
    fn dual_mode_happy_path_passes() {
        let m = balanced_mapping();
        let s1 = pos_cols();
        let s2 = neg_cols();
        let res = validate_for_dam(dual_mode_input(&m, &s1, &s2));
        assert!(res.is_ok(), "expected Ok, got {res:?}");
    }

    #[test]
    fn dual_mode_with_2_column_csv_blocks_with_specific_message() {
        let f = write_csv(
            "sample,group\n\
             P1,ctrl\nP2,ctrl\nP3,treat\nP4,treat\n\
             N1,ctrl\nN2,ctrl\nN3,treat\nN4,treat\n",
        );
        let union: Vec<String> = pos_cols().into_iter().chain(neg_cols()).collect();
        let m = load_group_mapping(f.path(), &union).unwrap();
        let s1 = pos_cols();
        let s2 = neg_cols();
        let issues = validate_for_dam(dual_mode_input(&m, &s1, &s2)).unwrap_err();
        assert!(
            issues
                .iter()
                .any(|s| s.contains("Dual-mode requires a 'biosample' column")),
            "issues: {issues:?}"
        );
    }

    #[test]
    fn dual_mode_same_mode_both_slots_blocks() {
        let m = balanced_mapping();
        let s1 = pos_cols();
        let s2 = neg_cols();
        let mut input = dual_mode_input(&m, &s1, &s2);
        input.slot1_mode = Some(IonMode::Positive);
        input.slot2_mode = Some(IonMode::Positive);
        let issues = validate_for_dam(input).unwrap_err();
        assert!(
            issues
                .iter()
                .any(|s| s.contains("POS and NEG must be different")),
            "issues: {issues:?}"
        );
    }

    #[test]
    fn dual_mode_group_count_mismatch_blocks_per_mode() {
        let m = balanced_mapping();
        let s1 = pos_cols();
        let s2 = cols(&["N1", "N2", "N3"]);
        let issues = validate_for_dam(dual_mode_input(&m, &s1, &s2)).unwrap_err();
        assert!(
            issues.iter().any(|s| s.contains("Group 'treat'")
                && s.contains("2 sample")
                && s.contains("POS")
                && s.contains("1 in NEG")),
            "issues: {issues:?}"
        );
    }

    #[test]
    fn dual_mode_biosample_appears_twice_in_same_mode_blocks() {
        let f = write_csv(
            "sample,biosample,group\n\
             P1,BIO-1,ctrl\nP2,BIO-1,ctrl\nP3,BIO-2,treat\nP4,BIO-3,treat\n\
             N1,BIO-1,ctrl\nN2,BIO-2,ctrl\nN3,BIO-3,treat\nN4,BIO-4,treat\n",
        );
        let union: Vec<String> = pos_cols().into_iter().chain(neg_cols()).collect();
        let m = load_group_mapping(f.path(), &union).unwrap();
        let s1 = pos_cols();
        let s2 = neg_cols();
        let issues = validate_for_dam(dual_mode_input(&m, &s1, &s2)).unwrap_err();
        assert!(
            issues.iter().any(|s| s.contains("Biosample 'BIO-1'")
                && s.contains("2 POS rows")
                && s.contains("must be unique")),
            "issues: {issues:?}"
        );
    }

    #[test]
    fn dual_mode_biosample_group_mismatch_blocks_with_names() {
        let f = write_csv(
            "sample,biosample,group\n\
             P1,BIO-1,ctrl\nP2,BIO-2,ctrl\nP3,BIO-3,treat\nP4,BIO-4,treat\n\
             N1,BIO-1,treat\nN2,BIO-2,ctrl\nN3,BIO-3,treat\nN4,BIO-4,treat\n",
        );
        let union: Vec<String> = pos_cols().into_iter().chain(neg_cols()).collect();
        let m = load_group_mapping(f.path(), &union).unwrap();
        let s1 = pos_cols();
        let s2 = neg_cols();
        let issues = validate_for_dam(dual_mode_input(&m, &s1, &s2)).unwrap_err();
        assert!(
            issues.iter().any(|s| s.contains("Biosample 'BIO-1'")
                && s.contains("group 'ctrl'")
                && s.contains("POS")
                && s.contains("'treat'")
                && s.contains("NEG")),
            "issues: {issues:?}"
        );
    }

    #[test]
    fn single_mode_with_3_column_csv_passes_biosample_ignored() {
        let m = balanced_mapping();
        let s1 = pos_cols();
        let input = Stage1ValidationInput {
            table_loaded: true,
            slot1_sample_cols: &s1,
            slot2_sample_cols: None,
            mapping: Some(&m),
            slot1_mode: Some(IonMode::Positive),
            slot2_revealed: false,
            slot2_mode: None,
        };
        assert!(validate_for_dam(input).is_ok());
    }

    #[test]
    fn missing_csv_blocks_the_gate_silently_without_a_nag() {
        // No metadata CSV picked yet: the gate MUST fail (so "Continue to DAM"
        // stays disabled) but emit NO "Upload …" / "Choose …" nag — mirroring
        // the no-`.txt` branch. Regression guard for the removed nag.
        let s1 = pos_cols();
        let input = Stage1ValidationInput {
            table_loaded: true,
            slot1_sample_cols: &s1,
            slot2_sample_cols: None,
            mapping: None,
            slot1_mode: Some(IonMode::Positive),
            slot2_revealed: false,
            slot2_mode: None,
        };
        let issues = validate_for_dam(input).expect_err("missing CSV must fail the gate");
        assert!(
            !issues
                .iter()
                .any(|s| s.contains("Upload the group mapping") || s.contains("Choose the group")),
            "no-CSV gate must render no nag, got: {issues:?}"
        );
    }

    #[test]
    fn slot1_mode_unset_blocks_in_single_mode() {
        let m = balanced_mapping();
        let s1 = pos_cols();
        let mut input = Stage1ValidationInput {
            table_loaded: true,
            slot1_sample_cols: &s1,
            slot2_sample_cols: None,
            mapping: Some(&m),
            slot1_mode: None,
            slot2_revealed: false,
            slot2_mode: None,
        };
        let issues = validate_for_dam(Stage1ValidationInput {
            slot1_mode: None,
            ..input
        })
        .unwrap_err();
        assert!(
            issues
                .iter()
                .any(|s| s.contains("Pick the ionization mode for the first .txt")),
            "issues: {issues:?}"
        );
        input.slot1_mode = Some(IonMode::Positive);
        assert!(validate_for_dam(input).is_ok());
    }

    #[test]
    fn slot2_revealed_without_file_blocks() {
        let m = balanced_mapping();
        let s1 = pos_cols();
        let input = Stage1ValidationInput {
            table_loaded: true,
            slot1_sample_cols: &s1,
            slot2_sample_cols: None,
            mapping: Some(&m),
            slot1_mode: Some(IonMode::Positive),
            slot2_revealed: true,
            slot2_mode: None,
        };
        let issues = validate_for_dam(input).unwrap_err();
        assert!(
            issues
                .iter()
                .any(|s| s.contains("Upload the second MS-DIAL .txt or remove slot #2")),
            "issues: {issues:?}"
        );
    }

    #[test]
    fn existing_single_mode_validation_passes_when_complete() {
        let f = write_csv("sample,group\nP1,ctrl\nP2,ctrl\nP3,treat\nP4,treat\n");
        let s1 = pos_cols();
        let m = load_group_mapping(f.path(), &s1).unwrap();
        let input = Stage1ValidationInput {
            table_loaded: true,
            slot1_sample_cols: &s1,
            slot2_sample_cols: None,
            mapping: Some(&m),
            slot1_mode: Some(IonMode::Positive),
            slot2_revealed: false,
            slot2_mode: None,
        };
        assert!(validate_for_dam(input).is_ok());
    }

    // ─── Auto-infer helpers (auto-infer-stage1-ion-mode §1.3) ───

    #[test]
    fn decide_slot1_mode_on_file_load_maps_each_variant() {
        assert_eq!(
            decide_slot1_mode_on_file_load(AdductPolarityInference::Positive),
            Some(IonMode::Positive)
        );
        assert_eq!(
            decide_slot1_mode_on_file_load(AdductPolarityInference::Negative),
            Some(IonMode::Negative)
        );
        assert_eq!(
            decide_slot1_mode_on_file_load(AdductPolarityInference::Ambiguous),
            None
        );
    }

    // ── 9-cell matrix for slot2_revealed = true ──
    // (slot1_new ∈ {None, Some(Pos), Some(Neg)}) × (slot2_before ∈ {None, Some(Pos), Some(Neg)})

    #[test]
    fn slot2_helper_revealed_slot1_none_returns_slot2_before_unchanged() {
        // slot1=None row collapses to "return slot2_before verbatim" for all three slot2_before
        // values (per tasks 1.3 — each cell asserted explicitly even though they all share the
        // same code path).
        assert_eq!(decide_slot2_mode_on_slot1_change(None, true, None), None);
        assert_eq!(
            decide_slot2_mode_on_slot1_change(None, true, Some(IonMode::Positive)),
            Some(IonMode::Positive)
        );
        assert_eq!(
            decide_slot2_mode_on_slot1_change(None, true, Some(IonMode::Negative)),
            Some(IonMode::Negative)
        );
    }

    #[test]
    fn slot2_helper_revealed_slot1_positive_fills_or_flips() {
        // slot2 unset → fill with opposite (Negative).
        assert_eq!(
            decide_slot2_mode_on_slot1_change(Some(IonMode::Positive), true, None),
            Some(IonMode::Negative)
        );
        // slot2 conflicts (also Positive) → flip to Negative (D3 trigger #3 conflict resolution).
        assert_eq!(
            decide_slot2_mode_on_slot1_change(
                Some(IonMode::Positive),
                true,
                Some(IonMode::Positive)
            ),
            Some(IonMode::Negative)
        );
        // slot2 already Negative (no-op) → return Negative unchanged.
        assert_eq!(
            decide_slot2_mode_on_slot1_change(
                Some(IonMode::Positive),
                true,
                Some(IonMode::Negative)
            ),
            Some(IonMode::Negative)
        );
    }

    #[test]
    fn slot2_helper_revealed_slot1_negative_fills_or_flips() {
        // slot2 unset → fill with opposite (Positive).
        assert_eq!(
            decide_slot2_mode_on_slot1_change(Some(IonMode::Negative), true, None),
            Some(IonMode::Positive)
        );
        // slot2 already Positive (no-op) → return Positive unchanged.
        assert_eq!(
            decide_slot2_mode_on_slot1_change(
                Some(IonMode::Negative),
                true,
                Some(IonMode::Positive)
            ),
            Some(IonMode::Positive)
        );
        // slot2 conflicts (also Negative) → flip to Positive.
        assert_eq!(
            decide_slot2_mode_on_slot1_change(
                Some(IonMode::Negative),
                true,
                Some(IonMode::Negative)
            ),
            Some(IonMode::Positive)
        );
    }

    // ── slot2_revealed = false: helper MUST return slot2_before verbatim for every (slot1_new, slot2_before) pair ──

    #[test]
    fn slot2_helper_not_revealed_returns_slot2_before_for_every_combo() {
        for slot1_new in [None, Some(IonMode::Positive), Some(IonMode::Negative)] {
            for slot2_before in [None, Some(IonMode::Positive), Some(IonMode::Negative)] {
                assert_eq!(
                    decide_slot2_mode_on_slot1_change(slot1_new, false, slot2_before),
                    slot2_before,
                    "slot2_revealed=false should always return slot2_before; \
                     got mismatch for slot1_new={slot1_new:?}, slot2_before={slot2_before:?}"
                );
            }
        }
    }
}
