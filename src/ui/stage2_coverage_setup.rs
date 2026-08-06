//! Stage 2 (coverage route) — the coverage setup screen.
//!
//! Renders, in the order the filters actually run: the sample-group
//! multi-select with its detection threshold, the deduplication block, the
//! SHARED analysis-target block (mode toggle + selector + inline fetch
//! progress), and `Run Coverage`.
//!
//! It deliberately renders no direction radio, no FDR control, no normalization
//! control, no numerator/denominator picker, no statistical-method radio, and
//! no `Drop Unknown features` checkbox — Owner: the `coverage-ui` capability.

use egui::RichText;

use crate::app::{App, AppState};
use crate::data::groups::UNASSIGNED;
use crate::theme;
use crate::ui::widgets::primary_button;

/// Hover text for the detection-threshold input.
const THRESHOLD_TOOLTIP: &str = "A feature counts as detected in a group when it has a real \
    measured intensity (present and above zero) in at least this fraction of that group's \
    samples. Features detected in no selected group are excluded.";

/// Grey hint under the group checkboxes.
///
/// It names QC pools and solvent blanks but the screen never GUESSES which
/// groups those are: group naming is the user's, and a wrong guess silently
/// discards real data.
const GROUP_HINT: &str =
    "Uncheck QC pools or solvent blanks so their compounds do not enter the results.";

/// Grey sub-hint under the deduplication controls.
///
/// Normative copy, not a suggestion. `kegg-coverage` proves deduplication
/// cannot move any reported number; a label that let a user believe otherwise
/// would be exactly the failure this change deleted the overview-map toggle to
/// avoid.
const DEDUP_HINT: &str = "Deduplication groups features by InChIKey, so it cannot change which \
    compounds are found or any coverage number. It changes only which metabolite names are \
    listed for each compound in the exported CSV.";

/// Disabled-hover text for the coverage-only Run gate.
const NO_GROUPS_HINT: &str = "Select at least one sample group.";

pub fn show(ui: &mut egui::Ui, app: &mut App) {
    let mut run_clicked = false;
    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            if !matches!(&app.state, AppState::Stage2CoverageSetup { .. }) {
                return;
            }
            ui.heading(RichText::new("Stage 2 — Setup").color(theme::HEADING));
            ui.add_space(8.0);

            render_group_block(ui, app);
            render_dedup_block(ui, app);

            // The T3-extracted shared surface: one call, no route-specific copy.
            crate::ui::stage3_setup::render_analysis_target(ui, app);

            run_clicked = render_run_button(ui, app);

            if let AppState::Stage2CoverageSetup { error, .. } = &app.state
                && let Some(msg) = error
            {
                ui.add_space(8.0);
                ui.colored_label(theme::ERROR, msg.clone());
            }
        });

    if run_clicked {
        start_coverage_run(app);
    }
}

/// The offered groups: every group in the mapping except `Unassigned`.
///
/// `Unassigned` is the ABSENCE of a group assignment, not a condition the user
/// chose to measure, so it is never offerable — which is also what lets the
/// presence filter treat it as unreachable.
fn offered_groups(app: &App) -> Vec<String> {
    app.inputs
        .mapping
        .as_ref()
        .map(|m| m.groups().into_iter().filter(|g| g != UNASSIGNED).collect())
        .unwrap_or_default()
}

/// Reconcile `coverage_selected_groups` with the groups currently on offer.
///
/// Two rules, and they must not be collapsed into one:
///
/// - **`None` → `Some(<every offered group>)`.** `None` means "not yet chosen",
///   so the default behaviour is "use everything I measured".
/// - **`Some(list)` with ANY entry absent from the offered set → reset to
///   `None`** (which the first rule then re-initialises), plus a persistent
///   notice naming the dropped groups.
///
/// The staleness test is **per-entry membership**, never a length or set
/// comparison. `list.len() != offered.len()` or `list != offered` would reset
/// the deliberately-none state every frame, re-check every box, and make the
/// Run gate unreachable.
///
/// `Some(vec![])` is never stale — it names no group, so "any entry is absent"
/// is vacuously false — and it is never re-initialised, because it is not
/// `None`. That is the whole reason the field is an `Option<Vec<_>>`.
///
/// Without the stale branch, swapping the metadata `.csv` at Stage 1 and
/// returning here is silently wrong: settings persist across every navigation
/// transition and `reset_for_back_to_stage1` is a no-op, so the field keeps the
/// OLD group names. Not `None`, so the initialise rule does not fire; not
/// `Some(vec![])`, so the Run gate does not fire; and every selected group then
/// matches zero sample columns — an EMPTY `D` after a multi-minute run.
fn reconcile_selection(app: &mut App, offered: &[String]) {
    if let Some(list) = app.settings.coverage_selected_groups.as_ref() {
        let dropped: Vec<String> = list
            .iter()
            .filter(|g| !offered.iter().any(|o| o == *g))
            .cloned()
            .collect();
        if !dropped.is_empty() {
            let notice = format!(
                "Group selection reset: {} {} not in the current metadata.",
                dropped
                    .iter()
                    .map(|g| format!("'{g}'"))
                    .collect::<Vec<_>>()
                    .join(", "),
                if dropped.len() == 1 { "is" } else { "are" }
            );
            tracing::warn!(
                dropped = dropped.len(),
                "coverage group selection reset: groups absent from the current metadata"
            );
            app.settings.coverage_selected_groups = None;
            // Held on the variant, not recomputed: the guard repairs its own
            // trigger condition on the frame it fires, so a condition-derived
            // warning would be visible for exactly one frame.
            if let AppState::Stage2CoverageSetup {
                stale_groups_notice,
                ..
            } = &mut app.state
            {
                *stale_groups_notice = Some(notice);
            }
        }
    }
    if app.settings.coverage_selected_groups.is_none() {
        app.settings.coverage_selected_groups = Some(offered.to_vec());
    }
}

/// The `Sample groups` block — rendered only when a metadata `.csv` was loaded.
fn render_group_block(ui: &mut egui::Ui, app: &mut App) {
    if app.inputs.mapping.is_none() {
        return;
    }
    let offered = offered_groups(app);
    reconcile_selection(app, &offered);

    let counts: Vec<usize> = {
        let mapping = app.inputs.mapping.as_ref().expect("checked above");
        offered
            .iter()
            .map(|g| mapping.samples_in(g).len())
            .collect()
    };

    crate::ui::widgets::section_header(ui, "Sample groups");

    let mut selection = app
        .settings
        .coverage_selected_groups
        .clone()
        .unwrap_or_default();
    let mut changed = false;

    for (g, n) in offered.iter().zip(&counts) {
        let mut checked = selection.iter().any(|s| s == g);
        if ui
            .checkbox(&mut checked, format!("{g} ({n} samples)"))
            .changed()
        {
            changed = true;
            if checked {
                selection.push(g.clone());
            } else {
                selection.retain(|s| s != g);
            }
        }
    }

    ui.horizontal(|ui| {
        if ui.button("Select all").clicked() {
            selection = offered.clone();
            changed = true;
        }
        if ui.button("Select none").clicked() {
            selection.clear();
            changed = true;
        }
    });

    if changed {
        // Re-derive from `offered` so the persisted list keeps the offered
        // order (alphabetical, per `GroupMapping::groups`) regardless of click
        // order — the saved snapshot is then click-order-independent.
        app.settings.coverage_selected_groups = Some(
            offered
                .iter()
                .filter(|g| selection.iter().any(|s| s == *g))
                .cloned()
                .collect(),
        );
        // The user has now made a selection of their own; the reset notice has
        // served its purpose.
        if let AppState::Stage2CoverageSetup {
            stale_groups_notice,
            ..
        } = &mut app.state
        {
            *stale_groups_notice = None;
        }
    }

    if let AppState::Stage2CoverageSetup {
        stale_groups_notice: Some(notice),
        ..
    } = &app.state
    {
        ui.colored_label(theme::WARNING, notice.clone());
    }

    ui.label(RichText::new(GROUP_HINT).small().color(theme::TEXT));
    ui.add_space(6.0);

    // Detection threshold — a percentage on screen, a fraction in settings.
    ui.horizontal(|ui| {
        ui.label("Detected in at least");
        let mut percent = app.settings.coverage_presence_threshold * 100.0;
        let resp = ui.add(
            egui::DragValue::new(&mut percent)
                .speed(1.0)
                .range(0.0..=100.0)
                .suffix(" %"),
        );
        if resp.changed() {
            app.settings.coverage_presence_threshold =
                crate::app::clamp_coverage_presence_threshold(percent / 100.0);
        }
        resp.on_hover_text(THRESHOLD_TOOLTIP);
        ui.label("of a group's samples");
    });
    ui.add_space(10.0);
}

/// The `Deduplication` block. Same labels and the same bound settings fields as
/// the Stage 2 DAM setup screen — it is the same operation on the same data.
///
/// Rendered BELOW the group block, matching the order the filters actually run
/// in (group presence → deduplication).
fn render_dedup_block(ui: &mut egui::Ui, app: &mut App) {
    crate::ui::widgets::section_header(ui, "Deduplication");
    let settings = &mut app.settings;
    ui.checkbox(
        &mut settings.dedup_enabled,
        "Deduplicate features by InChIKey (keep best per compound)",
    );
    if settings.dedup_enabled {
        ui.horizontal(|ui| {
            ui.label("RT tolerance (min)");
            ui.add(
                egui::DragValue::new(&mut settings.dedup_rt_tolerance_min)
                    .speed(0.01)
                    .range(crate::app::MIN_DEDUP_RT_TOLERANCE_MIN..=f64::MAX),
            )
            .on_hover_text(
                "Same-InChIKey features more than this far apart in retention \
                 time (± minutes) are kept as separate peaks instead of \
                 deduplicated.",
            );
        });
    }
    ui.label(RichText::new(DEDUP_HINT).small().color(theme::TEXT));
    ui.add_space(10.0);
}

/// `true` when the coverage-only gate fires: a `.csv` is loaded and the user has
/// deliberately unchecked every group.
///
/// `None` never reaches here — `reconcile_selection` has already replaced it
/// with the all-groups default on the same frame.
fn no_groups_selected(app: &App) -> bool {
    app.inputs.mapping.is_some()
        && app
            .settings
            .coverage_selected_groups
            .as_ref()
            .is_some_and(|s| s.is_empty())
}

/// `Run Coverage`, gated on the same conditions as `Run Enrichment` plus the
/// zero-selected-groups case. Returns `true` when clicked.
fn render_run_button(ui: &mut egui::Ui, app: &mut App) -> bool {
    let empty_selection = no_groups_selected(app);
    let enabled = crate::ui::stage3_setup::target_ready(app) && !empty_selection;
    // The empty-selection reason wins: it is the one the user can act on
    // without leaving the screen.
    let hint = if empty_selection {
        Some(NO_GROUPS_HINT)
    } else if enabled {
        None
    } else {
        crate::ui::stage3_setup::fetch_in_flight_hint(app)
    };

    let resp = primary_button(ui, "Run Coverage", enabled);
    let resp = match hint {
        Some(h) => resp.on_disabled_hover_text(h),
        None => resp,
    };
    resp.clicked() && enabled
}

/// Build the spawn request and hand off to the SHARED spawn helper.
///
/// A second spawn helper is deliberately not introduced: the channel-trio
/// wiring, the `AbortHandle` capture, and the `Stage3EnrichRunning` literal are
/// exactly the plumbing whose per-call-site duplication caused the
/// `fix-stage3-ui-dual-mode-spawn` regression.
fn start_coverage_run(app: &mut App) {
    // Matches the enrichment route's `start_run`: the transition below replaces
    // `app.state`, so a still-streaming other-mode fetch (or a refresh) would be
    // orphaned by it — its `AbortHandle` goes out of scope with the old state.
    // Unreachable today via the Run gate, exactly as on the enrichment route,
    // which carries this guard anyway.
    if crate::app::is_busy(&app.state) {
        tracing::info!("stopping the in-flight setup fetch: starting the coverage run");
    }
    crate::app::abort_in_flight(&app.state);

    let Some(target) = crate::ui::stage3_setup::build_analysis_payload(&app.settings, &app.cache)
    else {
        if let AppState::Stage2CoverageSetup { error, .. } = &mut app.state {
            *error = Some("KEGG cache for the selected mode is missing.".to_string());
        }
        return;
    };

    // Same definition as the enrichment route's: the additive
    // InChIKey-bearing-feature count across ALL modes, which is what the
    // resolver's undeduped caller-list length reports back, so the progress
    // bar's seed and its denominator agree.
    let pubchem_total: usize = app
        .inputs
        .ion_tables
        .iter()
        .map(|t| {
            t.table
                .features
                .iter()
                .filter(|f| f.inchikey.is_some())
                .count()
        })
        .sum();

    let params = crate::stage3::CoverageParams {
        selected_groups: app.settings.coverage_selected_groups.clone(),
        presence_threshold: app.settings.coverage_presence_threshold,
        dedup_enabled: app.settings.dedup_enabled,
        dedup_rt_tolerance_min: app.settings.dedup_rt_tolerance_min,
        force_refresh_pubchem: false,
        force_refresh_kegg_conv: false,
    };

    // The group filter and deduplication run HERE, synchronously, before the
    // spawn: they are pure and fast, and doing them now is what keeps the
    // intensity matrices out of the spawned future entirely.
    let prepared = crate::stage3::prepare_features(
        crate::stage3::CoverageInputs {
            ion_tables: &app.inputs.ion_tables,
            mapping: app.inputs.mapping.as_ref(),
        },
        &params,
    );

    tracing::info!(
        mode = ?app.settings.analysis_mode,
        n_modes = app.inputs.ion_tables.len(),
        pubchem_inputs = pubchem_total,
        inchikeys = prepared.inchikeys.len(),
        selected_groups = params.selected_groups.as_ref().map_or(0, Vec::len),
        "Coverage Run starting"
    );

    app.spawn_stage3_run(crate::app::RunSpawn {
        payload: crate::app::RunPayloadSpec::Coverage(Box::new(crate::app::CoverageSpawn {
            prepared,
            force_refresh_pubchem: params.force_refresh_pubchem,
            force_refresh_kegg_conv: params.force_refresh_kegg_conv,
        })),
        target,
        pubchem_total,
    });
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The pinned copy. These three strings are the screen's entire explanation
    /// of what it does and does not do, so an edit should be a deliberate act.
    #[test]
    fn pinned_copy_matches_the_spec() {
        assert_eq!(
            GROUP_HINT,
            "Uncheck QC pools or solvent blanks so their compounds do not enter the results."
        );
        assert_eq!(NO_GROUPS_HINT, "Select at least one sample group.");
        // Whitespace-normalised comparison: the source wraps these for
        // readability, and the wrapping is typography, not contract.
        let norm = |s: &str| s.split_whitespace().collect::<Vec<_>>().join(" ");
        assert_eq!(
            norm(DEDUP_HINT),
            "Deduplication groups features by InChIKey, so it cannot change which compounds \
             are found or any coverage number. It changes only which metabolite names are \
             listed for each compound in the exported CSV."
        );
        assert_eq!(
            norm(THRESHOLD_TOOLTIP),
            "A feature counts as detected in a group when it has a real measured intensity \
             (present and above zero) in at least this fraction of that group's samples. \
             Features detected in no selected group are excluded."
        );
    }
}
