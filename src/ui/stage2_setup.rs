use egui::RichText;
use std::sync::mpsc;
use tracing::{error, info, warn};

use crate::app::{App, AppState, SessionSettings};
use crate::dam;
use crate::dam::DamMethod;
use crate::data::groups::UNASSIGNED;
use crate::normalize::{NormalizationMethod, PqnReference};
use crate::theme;

/// Radio-tag enum: lets us bind a radio group to a discriminant without
/// caring about the inner payload (which is stored in `metadata_column` /
/// `pqn_reference_group` instead).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum NormKind {
    None,
    Sum,
    Median,
    Metadata,
    Quantile,
    Pqn,
}

impl NormKind {
    fn from_method(m: &NormalizationMethod) -> Self {
        match m {
            NormalizationMethod::None => Self::None,
            NormalizationMethod::Sum => Self::Sum,
            NormalizationMethod::Median => Self::Median,
            NormalizationMethod::Metadata { .. } => Self::Metadata,
            NormalizationMethod::Quantile => Self::Quantile,
            NormalizationMethod::Pqn { .. } => Self::Pqn,
        }
    }
}

pub fn show(ui: &mut egui::Ui, app: &mut App) {
    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
        // Stay on the Stage 2 setup screen.
        if !matches!(&app.state, AppState::Stage2DamSetup { .. }) {
            return;
        }
        let error_text: Option<String> = match &app.state {
            AppState::Stage2DamSetup { error } => error.clone(),
            _ => None,
        };

        // Inputs (immutable read).
        let mapping = match app.inputs.mapping.as_ref() {
            Some(m) => m,
            None => return,
        };
        let ion_tables = app.inputs.ion_tables.as_slice();
        if ion_tables.is_empty() {
            return;
        }

        // User-tunable settings (mutable, destructured so the borrow checker
        // allows independent &mut per field).
        let SessionSettings {
            numerator,
            denominator,
            dam_method,
            drop_unknown,
            dedup_enabled,
            dedup_rt_tolerance_min,
            normalization,
            metadata_column,
            pqn_reference,
            pqn_reference_group,
            log_transform,
            dam_fdr_method,
            ..
        } = &mut app.settings;

        ui.heading(egui::RichText::new("Stage 2 — DAM Setup").color(theme::HEADING));
        ui.add_space(6.0);
        // Back-navigation is handled by the global stage stepper
        // (`ui::stepper`); no per-screen Back button.
        ui.add_space(8.0);

        // The per-mode input summary (`Input: N features (A annotated · U Unknown)`)
        // moved to the bottom-panel Data tab (`data-summary-panel`); the setup
        // body now starts directly with the DAM controls.

        // Group dropdowns.
        let groups: Vec<String> = mapping
            .groups()
            .into_iter()
            .filter(|g| g != UNASSIGNED)
            .collect();

        ui.horizontal(|ui| {
            ui.label("Numerator group:");
            // Primary dropdown — Num/Den is the core comparison choice.
            crate::ui::widgets::primary_dropdown(ui, |ui| {
                egui::ComboBox::from_id_salt("num_group")
                    .selected_text(numerator.as_deref().unwrap_or("— pick one —"))
                    .show_ui(ui, |ui| {
                        // Both dropdowns list every group; a pick that collides
                        // with the sibling field swaps the two (see
                        // `apply_group_pick`).
                        for g in &groups {
                            let selected = numerator.as_deref() == Some(g.as_str());
                            let label = group_item_label(g, selected);
                            if ui.selectable_label(selected, label).clicked() {
                                apply_group_pick(numerator, denominator, g, &groups);
                            }
                        }
                    });
            });
        });

        ui.horizontal(|ui| {
            ui.label("Denominator group:");
            crate::ui::widgets::primary_dropdown(ui, |ui| {
                egui::ComboBox::from_id_salt("den_group")
                    .selected_text(denominator.as_deref().unwrap_or("— pick one —"))
                    .show_ui(ui, |ui| {
                        // Symmetric to the numerator dropdown: full list, swap on
                        // collision with the numerator.
                        for g in &groups {
                            let selected = denominator.as_deref() == Some(g.as_str());
                            let label = group_item_label(g, selected);
                            if ui.selectable_label(selected, label).clicked() {
                                apply_group_pick(denominator, numerator, g, &groups);
                            }
                        }
                    });
            });
        });

        ui.add_space(8.0);

        // ── Processing-order block ─────────────────────────────────────
        // Controls below are listed top-to-bottom in the order they fire
        // inside `run_dam` (see `src/dam/run.rs`):
        //   ① sample normalization — table-level, rewrites `intensity`
        //      from `intensity_raw` before any per-feature work.
        //   ② InChIKey dedup — computes the kept-index mask; the per-
        //      feature loop's step 3a then skips dedup-losers BEFORE
        //      the Drop-Unknown branch fires.
        //   ③ Drop Unknown — loop step 3b: skip features whose
        //      InChIKey is None.
        //   ④ Statistical method — loop step 3e: Welch / Student / BM.
        //   ⑤ FDR correction — applied once to the family of p-values
        //      AFTER the per-feature loop completes.
        // The order was set by `reorder-stage2-setup-by-processing-order`
        // (archived 2026-05-26 follow-up to `add-session-settings-io`).

        // ① Sample normalization radio + conditional sub-controls.
        ui.label("Sample normalization:");
        let mut kind = NormKind::from_method(normalization);
        let metadata_cols = mapping.metadata_column_names();
        let has_metadata = !metadata_cols.is_empty();
        let analyzed_groups: Vec<String> = mapping
            .groups()
            .into_iter()
            .filter(|g| g != UNASSIGNED)
            .collect();

        ui.radio_value(&mut kind, NormKind::None, "None");
        ui.radio_value(&mut kind, NormKind::Sum, "Sum");
        ui.radio_value(&mut kind, NormKind::Median, "Median");
        ui.horizontal(|ui| {
            ui.add_enabled_ui(has_metadata, |ui| {
                ui.radio_value(&mut kind, NormKind::Metadata, "Metadata column");
            });
            if !has_metadata {
                ui.label(
                    RichText::new("(add numeric columns to your metadata CSV, e.g. dry_weight)")
                        .small()
                        .color(theme::TEXT),
                );
            }
        });
        ui.radio_value(&mut kind, NormKind::Quantile, "Quantile Normalization");
        ui.radio_value(
            &mut kind,
            NormKind::Pqn,
            "Probabilistic Quotient Normalization (PQN)",
        );

        // Apply the radio selection back to `normalization`, picking up the
        // appropriate sub-state from the persistent fields.
        *normalization = match kind {
            NormKind::None => NormalizationMethod::None,
            NormKind::Sum => NormalizationMethod::Sum,
            NormKind::Median => NormalizationMethod::Median,
            NormKind::Metadata => NormalizationMethod::Metadata {
                column: metadata_column.clone().unwrap_or_default(),
            },
            NormKind::Quantile => NormalizationMethod::Quantile,
            NormKind::Pqn => NormalizationMethod::Pqn {
                reference: match pqn_reference {
                    PqnReference::AllSamples => PqnReference::AllSamples,
                    PqnReference::Group(_) => {
                        PqnReference::Group(pqn_reference_group.clone().unwrap_or_default())
                    }
                },
            },
        };

        // Conditional sub-controls for Metadata.
        if kind == NormKind::Metadata && has_metadata {
            ui.horizontal(|ui| {
                ui.add_space(20.0);
                ui.label("Column:");
                let selected_label = metadata_column.as_deref().unwrap_or("— pick one —");
                egui::ComboBox::from_id_salt("normalization_metadata_column")
                    .selected_text(selected_label)
                    .show_ui(ui, |ui| {
                        for col in &metadata_cols {
                            ui.selectable_value(metadata_column, Some(col.clone()), col);
                        }
                    });
            });
            // Pre-DAM warning: any sample without a value in this column will be
            // dropped from the analysis. Shown as soon as the user picks a column
            // with incomplete coverage so the dropped set is visible before they
            // click Start DAM.
            if let Some(col) = metadata_column.as_deref() {
                let dropped = crate::normalize::metadata::dropped_samples(mapping, col);
                if !dropped.is_empty() {
                    ui.horizontal(|ui| {
                        ui.add_space(20.0);
                        ui.colored_label(
                            theme::WARNING,
                            format!(
                                "⚠ {} sample(s) without a value in '{col}' will be dropped: {}",
                                dropped.len(),
                                dropped.join(", ")
                            ),
                        );
                    });
                }
            }
        }

        // Conditional sub-controls for PQN.
        if kind == NormKind::Pqn {
            ui.horizontal(|ui| {
                ui.add_space(20.0);
                ui.label("Reference:");
                let mut is_group = matches!(pqn_reference, PqnReference::Group(_));
                ui.radio_value(&mut is_group, false, "All samples");
                ui.radio_value(&mut is_group, true, "Group");
                // Apply selection back to the persistent state.
                *pqn_reference = if is_group {
                    PqnReference::Group(pqn_reference_group.clone().unwrap_or_default())
                } else {
                    PqnReference::AllSamples
                };
                if is_group {
                    let selected_label = pqn_reference_group.as_deref().unwrap_or("— pick one —");
                    egui::ComboBox::from_id_salt("normalization_pqn_group")
                        .selected_text(selected_label)
                        .show_ui(ui, |ui| {
                            for g in &analyzed_groups {
                                ui.selectable_value(pqn_reference_group, Some(g.clone()), g);
                            }
                        });
                }
            });
        }

        ui.add_space(8.0);

        // ② InChIKey deduplication checkbox.
        ui.checkbox(
            dedup_enabled,
            "Deduplicate features by InChIKey (keep best per compound)",
        );

        // ②a Retention-time tolerance — a sub-control of the dedup checkbox,
        // enabled only when dedup is on. Each retention-time cluster's RT span is
        // bounded by this (complete-linkage), so same-InChIKey features more than
        // this far apart (± minutes) are kept as separate peaks rather than
        // deduplicated. The value is kept strictly positive by the DragValue
        // floor plus the `SessionSettings::apply_snapshot` clamp (both at
        // `crate::app::MIN_DEDUP_RT_TOLERANCE_MIN`), so it is NOT re-clamped
        // here every frame. See the `msdial-deduplication` capability.
        ui.horizontal(|ui| {
            ui.label("RT tolerance (min)");
            let resp = ui.add_enabled(
                *dedup_enabled,
                egui::DragValue::new(dedup_rt_tolerance_min)
                    .speed(0.01)
                    .range(crate::app::MIN_DEDUP_RT_TOLERANCE_MIN..=f64::MAX),
            );
            if *dedup_enabled {
                resp.on_hover_text(
                    "Same-InChIKey features more than this far apart in retention \
                     time (± minutes) are kept as separate peaks instead of \
                     deduplicated.",
                );
            } else {
                resp.on_disabled_hover_text(
                    "Retention-time tolerance applies only when deduplication is enabled.",
                );
            }
        });

        // ③ Drop unknown checkbox.
        ui.checkbox(drop_unknown, "Drop unknown features (no InChIKey)");

        ui.add_space(8.0);

        // ③d Log transformation checkbox. Disabled when BM is selected. Per the
        // user's preference, the checkbox is also auto-UNCHECKED on every BM
        // render (write-through `log_transform = false`) so the disabled state
        // reads as "off and locked" rather than "on but locked" — the latter is
        // visually ambiguous about whether arcsinh would apply. When the user
        // switches back to Welch/Student the checkbox stays unchecked and they
        // must explicitly re-check it if they want arcsinh on (no prior-value
        // restoration). See dam-analysis BM-bypass requirement + stage2-ui
        // "BM disables and auto-clears..." scenario.
        let log_enabled = *dam_method != DamMethod::BrunnerMunzel;
        if !log_enabled {
            *log_transform = false;
        }
        let log_resp = ui.add_enabled(
            log_enabled,
            egui::Checkbox::new(log_transform, "Log transformation"),
        );
        if !log_enabled {
            log_resp.on_disabled_hover_text(
                "Brunner–Munzel is rank-based; arcsinh is monotone, so the toggle has no effect on the result.",
            );
        }
        ui.label(
            RichText::new("(arcsinh — generalised log)")
                .small()
                .color(theme::TEXT),
        );

        ui.add_space(8.0);

        // ④ Statistical method radio.
        ui.label("Statistical method:");
        ui.radio_value(
            dam_method,
            DamMethod::Student,
            "Student's t-test (equal variances)",
        );
        ui.radio_value(
            dam_method,
            DamMethod::Welch,
            "Welch's t-test (unequal variances)",
        );
        ui.radio_value(
            dam_method,
            DamMethod::BrunnerMunzel,
            "Brunner–Munzel test (non-parametric)",
        );

        ui.add_space(8.0);

        // ⑤ FDR correction radio. Independent of Stage 3's choice (per D3).
        // Stage 2 hides the `None` variant (Stage-3-only) → `include_none = false`.
        crate::ui::widgets::fdr_method_radios(ui, dam_fdr_method, false);
        ui.label(
            RichText::new(
                "BH is the literature default. \
                 BY is more conservative when features are correlated.",
            )
            .small()
            .color(theme::TEXT),
        );

        ui.add_space(12.0);

        if let Some(e) = &error_text {
            ui.colored_label(theme::ERROR, e.clone());
        }

        // Gate: groups picked AND (Metadata → column picked AND valid) AND
        // (PQN Group → reference group picked AND valid).
        let metadata_ok = !matches!(kind, NormKind::Metadata)
            || metadata_column
                .as_deref()
                .is_some_and(|c| metadata_cols.iter().any(|m| m == c));
        let pqn_group_ok = match (&kind, &*pqn_reference) {
            (NormKind::Pqn, PqnReference::Group(_)) => pqn_reference_group
                .as_deref()
                .is_some_and(|g| analyzed_groups.iter().any(|x| x == g)),
            _ => true,
        };
        // Side fix from `add-session-settings-io`: the gate trusted
        // `is_some()` without checking group membership. After
        // `reorder-gui-and-move-mode-to-stage3` made every reset API no-op,
        // a preserved `numerator` whose value disappeared from the new
        // metadata would sail past the gate and feed `run_dam` a 0-sample
        // group (NaN p-values throughout). Loading a settings snapshot from
        // an unrelated session amplifies this. Now: both groups must
        // additionally appear in `mapping.groups()` (excluding UNASSIGNED).
        let (num_in_groups, den_in_groups) =
            check_group_membership(numerator.as_deref(), denominator.as_deref(), &groups);
        let can_start =
            num_in_groups && den_in_groups && numerator != denominator && metadata_ok && pqn_group_ok;

        // Disabled-hover hint: specifically call out the group-membership
        // failure when both groups are Some but at least one is stale (the
        // ComboBox shows the persisted value via `selected_text` fallback,
        // but the dropdown list doesn't include it — without this hint the
        // user has no signal explaining why Start DAM is greyed).
        let groups_hint: Option<&'static str> = if !can_start
            && numerator.is_some()
            && denominator.is_some()
            && !(num_in_groups && den_in_groups)
        {
            Some("Numerator/denominator group not present in the loaded metadata.")
        } else {
            None
        };

        // Inline warning label above the Start DAM button — surfaces the
        // group-membership failure without requiring the user to hover.
        // Tooltip stays as the secondary affordance. Matches the existing
        // "dedup audit" / "metadata-column dropped samples" inline-warning
        // pattern in this file.
        if let Some(hint) = groups_hint {
            ui.colored_label(theme::WARNING, format!("⚠ {hint}"));
        }

        let mut start_clicked = false;
        let resp = crate::ui::widgets::primary_button(ui, "Start DAM", can_start);
        let resp = if let Some(hint) = groups_hint {
            resp.on_disabled_hover_text(hint)
        } else {
            resp
        };
        if resp.clicked() && can_start {
            start_clicked = true;
        }

        // Drop the &mut app.settings borrow before dispatching transitions.
        drop(metadata_cols);
        drop(analyzed_groups);
        // (settings destructure ends with the function scope; the bindings
        // above are mutable references into `app.settings`, which are
        // released as we exit the helper closures.)

        if start_clicked {
            start_dam(app);
        }
        });
}

fn start_dam(app: &mut App) {
    // Validate destination state shape and snapshot settings before mutating.
    if !matches!(&app.state, AppState::Stage2DamSetup { .. }) {
        return;
    }
    let original_ion_tables = app.inputs.ion_tables.as_slice();
    if original_ion_tables.is_empty() {
        warn!("start_dam called with empty ion_tables; aborting");
        return;
    }
    let original_mapping = match app.inputs.mapping.as_ref() {
        Some(m) => m.clone(),
        None => return,
    };

    // Stage 1 → Stage 2 boundary: narrow inputs to assigned samples only.
    // Owned by start_dam for the duration of this call; `app.inputs` is
    // never mutated.
    let (mapping, ion_tables) =
        match build_stage2_boundary_view(&original_mapping, original_ion_tables) {
            BoundaryResult::Ok {
                mapping,
                ion_tables,
            } => (mapping, ion_tables),
            BoundaryResult::Stage1GateRegression => {
                tracing::error!(
                    assigned_count = 0_usize,
                    "Stage 1 gate regression — no samples have group assignments at start_dam"
                );
                app.error_modal = crate::app::ErrorModalState::Open {
                    title: "Cannot start DAM".to_string(),
                    message: "No samples have group assignments. Re-pick the metadata CSV."
                        .to_string(),
                };
                return;
            }
        };

    let (numerator, denominator) = match (
        app.settings.numerator.clone(),
        app.settings.denominator.clone(),
    ) {
        (Some(n), Some(d)) => (n, d),
        _ => return,
    };
    let method = app.settings.dam_method;
    let drop_unknown = app.settings.drop_unknown;
    let dedup_enabled = app.settings.dedup_enabled;
    let dedup_rt_tolerance_min = app.settings.dedup_rt_tolerance_min;
    let log_transform = app.settings.log_transform;
    let fdr_method = app.settings.dam_fdr_method;
    let normalization = app.settings.normalization.clone();
    let metadata_column = app.settings.metadata_column.clone();
    let pqn_reference = app.settings.pqn_reference.clone();
    let pqn_reference_group = app.settings.pqn_reference_group.clone();

    let n_modes = ion_tables.len();
    info!(
        method = ?method,
        numerator = %numerator,
        denominator = %denominator,
        drop_unknown,
        n_modes,
        "starting DAM"
    );

    // Assemble the user's NormalizationConfig from the per-field state.
    let normalization_method = match &normalization {
        crate::normalize::NormalizationMethod::Metadata { column: _ } => {
            crate::normalize::NormalizationMethod::Metadata {
                column: metadata_column.clone().unwrap_or_default(),
            }
        }
        crate::normalize::NormalizationMethod::Pqn { reference: _ } => {
            let reference = match &pqn_reference {
                crate::normalize::PqnReference::Group(_) => crate::normalize::PqnReference::Group(
                    pqn_reference_group.clone().unwrap_or_default(),
                ),
                crate::normalize::PqnReference::AllSamples => {
                    crate::normalize::PqnReference::AllSamples
                }
            };
            crate::normalize::NormalizationMethod::Pqn { reference }
        }
        other => other.clone(),
    };
    let norm_config = crate::normalize::NormalizationConfig {
        method: normalization_method,
    };

    // Metadata-column preflight (per-mode): if dropping samples without a
    // value would leave either group below the 2-sample DAM minimum in ANY
    // mode, surface the error before doing any normalization work.
    // PR-L (Finding #4 from the PR-H/J review): collect failures across
    // ALL modes before bailing, so dual-mode users see every problem in
    // one cycle instead of having to fix POS, retry, then discover NEG.
    // PR-N (Finding #5): replace the hand-rolled
    // samples_in/sample_name/col_values.get nested lookup with
    // mapping.metadata_value_of(name, column) — the by-name accessor PR-H
    // introduced for exactly this dual-mode pattern.
    if let crate::normalize::NormalizationMethod::Metadata { column } = &norm_config.method {
        const MIN_PER_GROUP: usize = 2;
        if mapping.metadata_values(column).is_some() {
            let mut failures: Vec<String> = Vec::new();
            for it in &ion_tables {
                let count_remaining = |group_name: &str| -> usize {
                    it.table
                        .sample_cols
                        .iter()
                        .filter(|name| mapping.group_of(name) == group_name)
                        .filter_map(|name| mapping.metadata_value_of(name, column).flatten())
                        .count()
                };
                let num_remaining = count_remaining(&numerator);
                let den_remaining = count_remaining(&denominator);
                let failed = if num_remaining < MIN_PER_GROUP {
                    Some(("numerator", numerator.clone(), num_remaining))
                } else if den_remaining < MIN_PER_GROUP {
                    Some(("denominator", denominator.clone(), den_remaining))
                } else {
                    None
                };
                if let Some((label, group_name, remaining)) = failed {
                    let err = crate::normalize::NormalizationError::InsufficientSamplesAfterDrop {
                        group: format!("{label} '{group_name}'"),
                        remaining,
                        required: MIN_PER_GROUP,
                        column: column.clone(),
                    };
                    error!(error = %err, mode = %it.mode, "metadata preflight failed");
                    failures.push(if n_modes >= 2 {
                        format!("{}: {err}", it.mode)
                    } else {
                        format!("{err}")
                    });
                }
            }
            if !failures.is_empty() {
                let msg = failures.join("\n");
                app.state = AppState::Stage2DamSetup { error: Some(msg) };
                return;
            }
        }
    }

    // Synchronous per-mode normalization smoke-test: invoke `validate` on every
    // ion table BEFORE spawning any tokio worker; surface every per-mode
    // failure in one error message (PR-L, Finding #4). Pre-2026-05-26 this
    // hardcoded `ion_tables[0]` (POS only) — Finding #12 in the 2026-05-25
    // audit, fixed in PR-H. `validate` is the non-logging sibling of `apply`:
    // it runs the identical per-mode check but does NOT emit the `normalize:`
    // INFO line, so each mode's normalization is logged once (by the real DAM
    // worker that recomputes it inside its own async task) instead of twice.
    // The check's result is discarded either way. Cost is millisecond-scale
    // vs the value of catching NEG-only failures (e.g. PQN ReferenceAllNan
    // on NEG, Median ZeroFactor on a sparse NEG sample) at the Setup gate
    // instead of mid-flight in an async worker.
    let mut apply_failures: Vec<String> = Vec::new();
    for it in &ion_tables {
        if let Err(e) = crate::normalize::validate(
            &norm_config,
            &it.table.intensity_raw,
            &mapping,
            &it.table.sample_cols,
        ) {
            error!(error = %e, mode = %it.mode, "normalization preflight failed");
            apply_failures.push(if n_modes >= 2 {
                format!("{}: {e}", it.mode)
            } else {
                format!("{e}")
            });
        }
    }
    if !apply_failures.is_empty() {
        let msg = apply_failures.join("\n");
        app.state = AppState::Stage2DamSetup { error: Some(msg) };
        return;
    }

    // Per-mode fan-out: one tokio worker per IonModeTable.
    let (tx, rx) = mpsc::channel::<(usize, Result<dam::DamResult, String>)>();
    let mut mode_total: Vec<usize> = Vec::with_capacity(n_modes);
    let mut progress_rxs: Vec<mpsc::Receiver<dam::DamProgress>> = Vec::with_capacity(n_modes);
    let mut worker_handles: Vec<tokio::task::AbortHandle> = Vec::with_capacity(n_modes);

    // Bundle the DAM configuration values into one `DamConfig` before the
    // fan-out. Each mode's spawned worker gets a clone — the config is identical
    // across modes; only the per-call I/O (table / mapping / num / den) differs.
    // Replaces the per-task copies the prior positional signature needed
    // (`introduce-dam-config-struct` D3).
    let dam_config = dam::DamConfig {
        method,
        normalization: norm_config,
        drop_unknown,
        dedup_enabled,
        dedup_rt_tolerance_min,
        log_transform,
        fdr_method,
    };
    for (idx, it) in ion_tables.iter().enumerate() {
        mode_total.push(it.table.features.len());
        let table_clone = clone_table(&it.table);
        let mapping_clone = clone_mapping(&mapping);
        let num_clone = numerator.clone();
        let den_clone = denominator.clone();
        let config_for_task = dam_config.clone();
        let tx_for_task = tx.clone();

        let (per_mode_prog_tx, per_mode_prog_rx) = mpsc::channel::<dam::DamProgress>();
        progress_rxs.push(per_mode_prog_rx);

        let handle = app
            .rt
            .spawn(async move {
                let mut table_for_dam = table_clone;
                let result = match dam::run_dam(
                    &mut table_for_dam,
                    &mapping_clone,
                    &num_clone,
                    &den_clone,
                    &config_for_task,
                    Some(per_mode_prog_tx),
                )
                .await
                {
                    Ok(r) => Ok(r),
                    Err(e) => {
                        error!(mode_idx = idx, error = %e, "DAM failed");
                        Err(e.to_string())
                    }
                };
                let _ = tx_for_task.send((idx, result));
            })
            .abort_handle();
        worker_handles.push(handle);
    }
    drop(tx);

    let mode_completed = vec![0usize; n_modes];
    let dam_results_accum: Vec<Option<Result<dam::DamResult, String>>> =
        std::iter::repeat_with(|| None).take(n_modes).collect();

    app.state = AppState::Stage2DamRunning {
        result_rx: rx,
        progress_rxs,
        mode_completed,
        mode_total,
        dam_results_accum,
        worker_handles,
    };
}

/// Outcome of the Stage 1 → Stage 2 boundary narrowing. `Ok` carries the
/// assigned-only mapping + per-mode tables for the rest of `start_dam` to
/// consume. `Stage1GateRegression` is the defensive case (every sample is
/// Unassigned); reachable only if the Stage 1 gate that requires
/// `assigned_count > 0` has regressed. The caller (`start_dam`) is
/// responsible for the user-visible side effects (error toast then bail).
#[allow(clippy::large_enum_variant)]
pub(crate) enum BoundaryResult {
    Ok {
        mapping: crate::data::GroupMapping,
        ion_tables: Vec<crate::data::IonModeTable>,
    },
    Stage1GateRegression,
}

/// Build the Stage 1 → Stage 2 boundary view. Pure function: narrows
/// `original_mapping` and `original_ion_tables` via the
/// `without_unassigned_samples` helpers, emits a `tracing::info!` audit
/// log when at least one Unassigned sample was dropped (count only — no
/// names, per the `stage1-ui` "Bug-report bundle does not leak Unassigned
/// sample names" requirement), and fires a `debug_assert!` to catch
/// Stage 1 gate regressions loudly in debug / test builds. Returns
/// `BoundaryResult::Stage1GateRegression` in release builds when every
/// sample is Unassigned so the caller can degrade gracefully (error log
/// then user-visible toast). Mirrors ORA's K ⊆ N invariant pattern; see
/// the project's numerical / biological conventions.
pub(crate) fn build_stage2_boundary_view(
    original_mapping: &crate::data::GroupMapping,
    original_ion_tables: &[crate::data::IonModeTable],
) -> BoundaryResult {
    let unassigned_count = original_mapping.samples_in(UNASSIGNED).len();
    let mapping = original_mapping.without_unassigned_samples();

    debug_assert!(
        mapping.assigned_count() > 0,
        "Stage 1 gate regression: assigned_count == 0 reached start_dam"
    );
    if mapping.assigned_count() == 0 {
        return BoundaryResult::Stage1GateRegression;
    }

    if unassigned_count > 0 {
        info!(
            dropped_count = unassigned_count,
            "dropping {} Unassigned sample(s) before Stage 2 processing", unassigned_count
        );
    }

    let ion_tables: Vec<crate::data::IonModeTable> = original_ion_tables
        .iter()
        .map(|it| it.without_unassigned_samples(original_mapping))
        .collect();

    BoundaryResult::Ok {
        mapping,
        ion_tables,
    }
}

fn clone_table(t: &crate::data::MetabolomicsTable) -> crate::data::MetabolomicsTable {
    crate::data::MetabolomicsTable {
        annotated_count: t.annotated_count,
        features: t.features.clone(),
        sample_cols: t.sample_cols.clone(),
        intensity_raw: t.intensity_raw.clone(),
        intensity: t.intensity.clone(),
        excluded_cols: t.excluded_cols.clone(),
    }
}

fn clone_mapping(m: &crate::data::GroupMapping) -> crate::data::GroupMapping {
    m.clone()
}

/// Side-fix helper extracted from `Start DAM`'s gate so it is testable
/// without spinning up an egui context. Returns `(num_in_groups,
/// den_in_groups)`. `groups` MUST already have `UNASSIGNED` filtered out
/// (matches how the dropdowns build their option list).
///
/// Both predicates return `false` when the corresponding name is `None`
/// — the existing "groups not picked" gate covers that path; this
/// helper just adds the membership requirement on top.
pub fn check_group_membership(
    numerator: Option<&str>,
    denominator: Option<&str>,
    groups: &[String],
) -> (bool, bool) {
    let num_ok = numerator.is_some_and(|n| groups.iter().any(|g| g == n));
    let den_ok = denominator.is_some_and(|d| groups.iter().any(|g| g == d));
    (num_ok, den_ok)
}

/// Apply a group pick to `picked` (the field whose dropdown the user clicked
/// in), swapping with `sibling` on collision and guarding against relocating a
/// stale value.
///
/// Rule: set `*picked = Some(g)`. If `*sibling == Some(g)` (a colliding pick),
/// move `picked`'s PREVIOUS value into `*sibling` — but only when that previous
/// value is a current `groups` member; otherwise clear `*sibling` to `None`.
/// The single `prev.filter(is_member)` handles all three previous-value cases:
/// a valid group (full swap), a stale group (cleared — never relocated into the
/// sibling), and `None` (cleared). A non-colliding pick leaves `sibling`
/// untouched. This makes a two-group direction flip a single click while
/// preserving the `numerator != denominator` invariant.
///
/// `groups` MUST be the already-`UNASSIGNED`-filtered slice the dropdowns
/// render from (the caller filters it where the option list is built); the
/// helper trusts this and never re-checks for `UNASSIGNED`, so passing a raw
/// `mapping.groups()` would let `UNASSIGNED` become a relocatable sibling value.
/// Label for one item in a Primary group dropdown: white `ON_PRIMARY` text
/// when selected (it sits on the opaque `PRIMARY` selected-item fill), default
/// body text otherwise. Keeps the §3.3 Primary selected-item readable.
fn group_item_label(g: &str, selected: bool) -> RichText {
    let text = RichText::new(g);
    if selected {
        text.color(crate::theme::ON_PRIMARY)
    } else {
        text
    }
}

pub fn apply_group_pick(
    picked: &mut Option<String>,
    sibling: &mut Option<String>,
    g: &str,
    groups: &[String],
) {
    let collision = sibling.as_deref() == Some(g);
    let prev = picked.replace(g.to_string());
    if collision {
        *sibling = prev.filter(|v| groups.iter().any(|x| x == v));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // `clamp_rt_tolerance` + its floor live with the `SessionSettings` field in
    // `crate::app` (applied at the persistence boundary AND this UI).
    use crate::app::{MIN_DEDUP_RT_TOLERANCE_MIN, clamp_rt_tolerance};

    #[test]
    fn clamp_rt_tolerance_rejects_nonpositive_and_nonfinite() {
        assert_eq!(clamp_rt_tolerance(0.0), MIN_DEDUP_RT_TOLERANCE_MIN);
        assert_eq!(clamp_rt_tolerance(-1.0), MIN_DEDUP_RT_TOLERANCE_MIN);
        assert_eq!(clamp_rt_tolerance(f64::NAN), MIN_DEDUP_RT_TOLERANCE_MIN);
        assert_eq!(
            clamp_rt_tolerance(f64::INFINITY),
            MIN_DEDUP_RT_TOLERANCE_MIN
        );
        assert_eq!(clamp_rt_tolerance(0.0005), MIN_DEDUP_RT_TOLERANCE_MIN);
    }

    #[test]
    fn clamp_rt_tolerance_passes_valid_values() {
        assert_eq!(clamp_rt_tolerance(0.1), 0.1);
        assert_eq!(
            clamp_rt_tolerance(MIN_DEDUP_RT_TOLERANCE_MIN),
            MIN_DEDUP_RT_TOLERANCE_MIN
        );
        assert_eq!(clamp_rt_tolerance(5.0), 5.0);
    }

    #[test]
    fn default_dedup_rt_tolerance_is_point_one() {
        assert_eq!(SessionSettings::default().dedup_rt_tolerance_min, 0.1);
    }

    #[test]
    fn both_groups_present_pass() {
        let groups = vec!["A".to_string(), "B".to_string()];
        let (num, den) = check_group_membership(Some("A"), Some("B"), &groups);
        assert!(num);
        assert!(den);
    }

    #[test]
    fn stale_numerator_fails_only_numerator() {
        let groups = vec!["A".to_string(), "B".to_string()];
        let (num, den) = check_group_membership(Some("Treated"), Some("A"), &groups);
        assert!(!num);
        assert!(den);
    }

    #[test]
    fn stale_denominator_fails_only_denominator() {
        let groups = vec!["A".to_string(), "B".to_string()];
        let (num, den) = check_group_membership(Some("A"), Some("Control"), &groups);
        assert!(num);
        assert!(!den);
    }

    #[test]
    fn none_inputs_both_fail() {
        let groups = vec!["A".to_string(), "B".to_string()];
        let (num, den) = check_group_membership(None, None, &groups);
        assert!(!num);
        assert!(!den);
    }

    #[test]
    fn empty_groups_list_always_fails() {
        let groups: Vec<String> = vec![];
        let (num, den) = check_group_membership(Some("A"), Some("B"), &groups);
        assert!(!num);
        assert!(!den);
    }

    // --- apply_group_pick tests (swap-on-collision + stale guard) ---

    #[test]
    fn pick_swaps_on_collision_numerator_side() {
        // (a) two-group flip: num=A, den=B; pick B in numerator -> (B, A).
        let groups = cols_vec(&["A", "B"]);
        let mut num = Some("A".to_string());
        let mut den = Some("B".to_string());
        apply_group_pick(&mut num, &mut den, "B", &groups);
        assert_eq!(num.as_deref(), Some("B"));
        assert_eq!(den.as_deref(), Some("A"));
    }

    #[test]
    fn pick_swaps_on_collision_denominator_side() {
        // (b) symmetric flip via the denominator dropdown: same (A,B) start,
        // pick A in denominator (picked=den, sibling=num) -> (B, A).
        let groups = cols_vec(&["A", "B"]);
        let mut num = Some("A".to_string());
        let mut den = Some("B".to_string());
        apply_group_pick(&mut den, &mut num, "A", &groups);
        assert_eq!(num.as_deref(), Some("B"));
        assert_eq!(den.as_deref(), Some("A"));
    }

    #[test]
    fn non_colliding_pick_leaves_sibling_unchanged() {
        // (c) num=A, den=B; pick C in numerator -> (C, B), sibling untouched.
        let groups = cols_vec(&["A", "B", "C"]);
        let mut num = Some("A".to_string());
        let mut den = Some("B".to_string());
        apply_group_pick(&mut num, &mut den, "C", &groups);
        assert_eq!(num.as_deref(), Some("C"));
        assert_eq!(den.as_deref(), Some("B"));
    }

    #[test]
    fn stale_previous_value_not_relocated_into_sibling() {
        // (d) stale guard: num="Treated" (not in groups), den="A"; pick A in
        // numerator -> (A, None). The stale value is cleared, NOT moved to den.
        let groups = cols_vec(&["A", "B"]);
        let mut num = Some("Treated".to_string());
        let mut den = Some("A".to_string());
        apply_group_pick(&mut num, &mut den, "A", &groups);
        assert_eq!(num.as_deref(), Some("A"));
        assert_eq!(den, None);
    }

    #[test]
    fn empty_picked_field_collision_clears_sibling() {
        // (e) num=None, den="A"; pick A in numerator -> (A, None). No valid
        // previous value to swap in, so the sibling is cleared.
        let groups = cols_vec(&["A", "B"]);
        let mut num: Option<String> = None;
        let mut den = Some("A".to_string());
        apply_group_pick(&mut num, &mut den, "A", &groups);
        assert_eq!(num.as_deref(), Some("A"));
        assert_eq!(den, None);
    }

    #[test]
    fn single_group_ping_pong_never_sets_both() {
        // (f) groups=[A]: alternating picks can never leave both fields set.
        let groups = cols_vec(&["A"]);
        let mut num: Option<String> = None;
        let mut den: Option<String> = None;

        apply_group_pick(&mut num, &mut den, "A", &groups);
        assert_eq!(num.as_deref(), Some("A"));
        assert_eq!(den, None);
        assert!(!(num.is_some() && den.is_some()));

        // Pick A in the denominator dropdown: collides with num, swaps.
        apply_group_pick(&mut den, &mut num, "A", &groups);
        assert_eq!(den.as_deref(), Some("A"));
        assert_eq!(num, None);
        assert!(!(num.is_some() && den.is_some()));
    }

    // --- build_stage2_boundary_view tests ---

    use crate::data::{IonMode, IonModeTable, MetabolomicsTable, load_group_mapping};
    use ndarray::Array2;
    use std::io::Write;

    fn write_csv(content: &str) -> tempfile::NamedTempFile {
        let mut f = tempfile::NamedTempFile::new().expect("create tempfile");
        f.write_all(content.as_bytes()).expect("write fixture");
        f
    }

    fn cols_vec(names: &[&str]) -> Vec<String> {
        names.iter().map(|s| s.to_string()).collect()
    }

    fn make_table_with_cols(names: &[&str]) -> MetabolomicsTable {
        let n = names.len();
        let intensity = Array2::<f64>::zeros((1, n));
        MetabolomicsTable {
            annotated_count: 0,
            features: vec![],
            sample_cols: cols_vec(names),
            intensity_raw: intensity.clone(),
            intensity,
            excluded_cols: vec![],
        }
    }

    fn make_ion_table(mode: IonMode, names: &[&str]) -> IonModeTable {
        IonModeTable {
            mode,
            table: make_table_with_cols(names),
            txt_path: None,
        }
    }

    #[test]
    fn boundary_view_narrows_mapping_and_per_mode_tables() {
        let f = write_csv("sample,group\nA,g1\nB,g1\nC,g2\nD,g2\n");
        let mapping = load_group_mapping(
            f.path(),
            &cols_vec(&["A", "Blank1", "B", "C", "Blank2", "D"]),
        )
        .unwrap();
        let ion_tables = vec![make_ion_table(
            IonMode::Positive,
            &["A", "Blank1", "B", "C", "Blank2", "D"],
        )];

        let result = build_stage2_boundary_view(&mapping, &ion_tables);
        match result {
            BoundaryResult::Ok {
                mapping: m_out,
                ion_tables: tables_out,
            } => {
                assert_eq!(m_out.assigned_count(), 4);
                assert!(!m_out.groups().contains(&"Unassigned".to_string()));
                assert_eq!(tables_out.len(), 1);
                assert_eq!(
                    tables_out[0].table.sample_cols,
                    cols_vec(&["A", "B", "C", "D"])
                );
            }
            BoundaryResult::Stage1GateRegression => panic!("expected Ok, got regression"),
        }
    }

    #[test]
    fn boundary_view_emits_info_log_when_unassigned_present() {
        use crate::logging::{LogLayer, LogStore};
        use tracing::Level;
        use tracing_subscriber::Registry;
        use tracing_subscriber::layer::SubscriberExt;

        let store = LogStore::new(100);
        let subscriber = Registry::default().with(LogLayer::new(store.clone()));

        let f = write_csv("sample,group\nA,g1\nB,g1\nC,g2\nD,g2\n");
        let mapping = load_group_mapping(
            f.path(),
            &cols_vec(&["A", "Blank1", "B", "C", "Blank2", "D"]),
        )
        .unwrap();
        let ion_tables = vec![make_ion_table(
            IonMode::Positive,
            &["A", "Blank1", "B", "C", "Blank2", "D"],
        )];

        tracing::subscriber::with_default(subscriber, || {
            let _ = build_stage2_boundary_view(&mapping, &ion_tables);
        });

        let info_events: Vec<_> = store
            .snapshot()
            .into_iter()
            .filter(|l| l.level == Level::INFO && l.message.contains("Unassigned sample"))
            .collect();
        assert_eq!(
            info_events.len(),
            1,
            "expected exactly one INFO event mentioning Unassigned samples; got: {:?}",
            info_events.iter().map(|l| &l.message).collect::<Vec<_>>()
        );
        assert!(
            info_events[0].message.contains("dropping 2"),
            "expected 'dropping 2' in message body; got: {}",
            info_events[0].message
        );
        // Privacy: no sample NAMES in the log body.
        assert!(
            !info_events[0].message.contains("Blank1"),
            "log must not name samples; got: {}",
            info_events[0].message
        );
        assert!(
            !info_events[0].message.contains("Blank2"),
            "log must not name samples; got: {}",
            info_events[0].message
        );
    }

    #[test]
    fn boundary_view_does_not_emit_info_log_when_no_unassigned() {
        use crate::logging::{LogLayer, LogStore};
        use tracing::Level;
        use tracing_subscriber::Registry;
        use tracing_subscriber::layer::SubscriberExt;

        let store = LogStore::new(100);
        let subscriber = Registry::default().with(LogLayer::new(store.clone()));

        let f = write_csv("sample,group\nA,g1\nB,g1\n");
        let mapping = load_group_mapping(f.path(), &cols_vec(&["A", "B"])).unwrap();
        let ion_tables = vec![make_ion_table(IonMode::Positive, &["A", "B"])];

        tracing::subscriber::with_default(subscriber, || {
            let _ = build_stage2_boundary_view(&mapping, &ion_tables);
        });

        let info_events: Vec<_> = store
            .snapshot()
            .into_iter()
            .filter(|l| l.level == Level::INFO && l.message.contains("Unassigned sample"))
            .collect();
        assert!(
            info_events.is_empty(),
            "expected no INFO event for zero-unassigned case; got: {:?}",
            info_events.iter().map(|l| &l.message).collect::<Vec<_>>()
        );
    }

    #[test]
    #[should_panic(expected = "Stage 1 gate regression")]
    fn boundary_view_panics_in_debug_when_all_samples_unassigned() {
        // Mapping built from a CSV that has zero overlap with sample_cols
        // → every sample is Unassigned → assigned_count == 0.
        let f = write_csv("sample,group\nX99,g1\n");
        let mapping = load_group_mapping(f.path(), &cols_vec(&["A", "B"])).unwrap();
        assert_eq!(mapping.assigned_count(), 0);

        let ion_tables = vec![make_ion_table(IonMode::Positive, &["A", "B"])];
        // debug_assert! panics here in debug / test builds (cargo test runs debug).
        let _ = build_stage2_boundary_view(&mapping, &ion_tables);
    }
}
