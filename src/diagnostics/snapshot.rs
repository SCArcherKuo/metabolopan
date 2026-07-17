//! Human-readable dump of the current `AppState` + `SessionSettings` /
//! `SessionInputs` / `SessionCache` for the bug-report bundle's
//! `app_state.txt`. Match is exhaustive on `AppState` — adding a new
//! variant will be a compile error here, which is the intended safety
//! net (do NOT add `_ =>` arms).

use std::collections::HashMap;
use std::fmt::Write;

use crate::app::{AppState, SessionCache, SessionInputs, SessionSettings};
use crate::data::{GroupMapping, IonModeTable};

/// Returns a plain-text block describing the current `AppState` variant
/// and the user-facing parameters carried by the sibling
/// `SessionSettings` / `SessionInputs` / `SessionCache` structs.
/// Excludes raw intensity matrices, `TextureHandle`s, and
/// `mpsc::Receiver`s.
pub fn render_app_state(
    state: &AppState,
    settings: &SessionSettings,
    inputs: &SessionInputs,
    cache: &SessionCache,
) -> String {
    let mut s = String::new();

    // ── Settings (every user-tunable parameter, regardless of variant) ──
    writeln!(s, "[settings]").ok();
    writeln!(s, "analysis_mode: {:?}", settings.analysis_mode).ok();
    writeln!(
        s,
        "kegg_species: {}",
        opt_str(settings.kegg_species.as_deref())
    )
    .ok();
    writeln!(
        s,
        "organism_group_level: {}",
        opt_dbg(settings.organism_group_level.as_ref())
    )
    .ok();
    writeln!(
        s,
        "organism_group: {}",
        opt_str(settings.organism_group.as_deref())
    )
    .ok();
    writeln!(s, "min_group_overlap: {}", settings.min_group_overlap).ok();
    writeln!(s, "numerator: {}", opt_str(settings.numerator.as_deref())).ok();
    writeln!(
        s,
        "denominator: {}",
        opt_str(settings.denominator.as_deref())
    )
    .ok();
    writeln!(s, "dam_method: {:?}", settings.dam_method).ok();
    writeln!(s, "drop_unknown: {}", settings.drop_unknown).ok();
    writeln!(s, "dedup_enabled: {}", settings.dedup_enabled).ok();
    writeln!(
        s,
        "dedup_rt_tolerance_min: {}",
        settings.dedup_rt_tolerance_min
    )
    .ok();
    writeln!(s, "log_transform: {}", settings.log_transform).ok();
    writeln!(s, "normalization: {:?}", settings.normalization).ok();
    writeln!(
        s,
        "metadata_column: {}",
        opt_str(settings.metadata_column.as_deref())
    )
    .ok();
    writeln!(s, "pqn_reference: {:?}", settings.pqn_reference).ok();
    writeln!(
        s,
        "pqn_reference_group: {}",
        opt_str(settings.pqn_reference_group.as_deref())
    )
    .ok();
    writeln!(s, "dam_fdr_method: {:?}", settings.dam_fdr_method).ok();
    writeln!(s, "fc_threshold: {}", settings.fc_threshold).ok();
    writeln!(s, "fdr_threshold: {}", settings.fdr_threshold).ok();
    writeln!(s, "delta_threshold: {}", settings.delta_threshold).ok();
    writeln!(
        s,
        "stage2_export: {}in x {}in @ {}dpi",
        settings.stage2_export_width_in,
        settings.stage2_export_height_in,
        settings.stage2_export_dpi
    )
    .ok();
    writeln!(s, "direction: {:?}", settings.direction).ok();
    writeln!(s, "top_n: {}", settings.top_n).ok();
    writeln!(
        s,
        "enrichment_fdr_threshold: {}",
        settings.enrichment_fdr_threshold
    )
    .ok();
    writeln!(s, "min_hit_count: {}", settings.min_hit_count).ok();
    writeln!(s, "min_entry_size: {}", settings.min_entry_size).ok();
    writeln!(
        s,
        "enrichment_fdr_method: {:?}",
        settings.enrichment_fdr_method
    )
    .ok();
    writeln!(
        s,
        "stage3_export: {}in x {}in @ {}dpi",
        settings.stage3_export_width_in,
        settings.stage3_export_height_in,
        settings.stage3_export_dpi
    )
    .ok();

    // ── Inputs (loaded raw data summary, no contents) ──
    writeln!(s, "\n[inputs]").ok();
    write_ion_tables(&mut s, &inputs.ion_tables);
    writeln!(
        s,
        "csv_path: {}",
        inputs
            .csv_path
            .as_ref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "<none>".into())
    )
    .ok();
    write_mapping(&mut s, inputs.mapping.as_ref());

    // ── Cache (per-analysis fetched KEGG data) ──
    writeln!(s, "\n[cache]").ok();
    writeln!(
        s,
        "species_kegg_loaded: {} (code: {})",
        cache.species_kegg.is_some(),
        cache
            .species_kegg
            .as_ref()
            .map(|sk| sk.code.as_str())
            .unwrap_or("<n/a>")
    )
    .ok();
    writeln!(
        s,
        "modules_pack_loaded: {} (modules: {})",
        cache.modules_pack.is_some(),
        cache.modules_pack.as_ref().map_or(0, |p| p.modules.len())
    )
    .ok();
    writeln!(
        s,
        "group_org_codes_count: {}",
        cache.group_org_codes.as_ref().map_or(0, |c| c.len())
    )
    .ok();

    // ── Variant + variant-local runtime fields ──
    writeln!(s, "\n[state]").ok();
    match state {
        AppState::Initializing {
            fallback_cache,
            last_error,
            ..
        } => {
            writeln!(s, "Variant: Initializing").ok();
            writeln!(s, "fallback_cache_present: {}", fallback_cache.is_some()).ok();
            writeln!(s, "last_error: {}", opt_str(last_error.as_deref())).ok();
        }
        AppState::Stage1Input {
            slot1_mode,
            slot2_revealed,
            slot2_mode,
            error,
        } => {
            writeln!(s, "Variant: Stage1Input").ok();
            writeln!(s, "slot1_mode: {}", opt_dbg(slot1_mode.as_ref())).ok();
            writeln!(s, "slot2_revealed: {slot2_revealed}").ok();
            writeln!(s, "slot2_mode: {}", opt_dbg(slot2_mode.as_ref())).ok();
            writeln!(s, "error: {}", opt_str(error.as_deref())).ok();
        }
        AppState::Stage2DamSetup { error } => {
            writeln!(s, "Variant: Stage2DamSetup").ok();
            writeln!(s, "error: {}", opt_str(error.as_deref())).ok();
        }
        AppState::Stage2DamRunning {
            mode_completed,
            mode_total,
            ..
        } => {
            writeln!(s, "Variant: Stage2DamRunning").ok();
            for (idx, (done, total)) in mode_completed.iter().zip(mode_total.iter()).enumerate() {
                writeln!(s, "mode_{idx}_progress: {done}/{total}").ok();
            }
        }
        AppState::Stage2DamThreshold {
            dam_results,
            active_volcano_tab,
            rendering,
            ..
        } => {
            writeln!(s, "Variant: Stage2DamThreshold").ok();
            writeln!(s, "active_volcano_tab: {active_volcano_tab:?}").ok();
            writeln!(s, "rendering: {rendering}").ok();
            write_dam_results(&mut s, dam_results);
        }
        AppState::Stage3EnrichSetup {
            dam_results,
            error,
            kegg_fetch,
            modules_fetch,
        } => {
            writeln!(s, "Variant: Stage3EnrichSetup").ok();
            writeln!(s, "error: {}", opt_str(error.as_deref())).ok();
            match kegg_fetch.as_ref() {
                None => {
                    writeln!(s, "kegg_fetch: <none in flight>").ok();
                }
                Some(f) => {
                    writeln!(
                        s,
                        "kegg_fetch: in_flight completed={} total={} current_pathway={}",
                        f.completed, f.total, f.current_pathway
                    )
                    .ok();
                }
            }
            match modules_fetch.as_ref() {
                None => {
                    writeln!(s, "modules_fetch: <none in flight>").ok();
                }
                Some(f) => {
                    writeln!(
                        s,
                        "modules_fetch: in_flight completed={} total={} current_id={} eta_secs={}",
                        f.completed,
                        f.total,
                        f.current_id,
                        f.eta_secs
                            .map(|v| v.to_string())
                            .unwrap_or_else(|| "<n/a>".into())
                    )
                    .ok();
                }
            }
            write_dam_results(&mut s, dam_results);
        }
        AppState::Stage3EnrichRunning {
            dam_results,
            phase,
            pubchem_completed,
            pubchem_total,
            kegg_conv_completed,
            kegg_conv_total,
            ..
        } => {
            writeln!(s, "Variant: Stage3EnrichRunning").ok();
            writeln!(s, "phase: {phase:?}").ok();
            writeln!(s, "pubchem_progress: {pubchem_completed}/{pubchem_total}").ok();
            writeln!(
                s,
                "kegg_conv_progress: {kegg_conv_completed}/{kegg_conv_total}"
            )
            .ok();
            write_dam_results(&mut s, dam_results);
        }
        AppState::Stage3EnrichResult {
            dam_results,
            module_retention,
            enrichment_result,
            mapped_universe,
            feature_to_cpds,
            pubchem_time_span,
            kegg_conv_time_span,
            dual_mode_breakdown,
            rendering,
            refresh_state,
            ..
        } => {
            writeln!(s, "Variant: Stage3EnrichResult").ok();
            writeln!(s, "universe_size: {}", mapped_universe.len()).ok();
            writeln!(s, "entries_tested: {}", enrichment_result.rows.len()).ok();
            let entries_hit = enrichment_result.rows.iter().filter(|r| r.hits > 0).count();
            writeln!(s, "entries_hit: {entries_hit}").ok();
            writeln!(s, "feature_to_cpds_size: {}", feature_to_cpds.len()).ok();
            writeln!(
                s,
                "enrichment_result_fdr_method: {:?}",
                enrichment_result.fdr_method
            )
            .ok();
            writeln!(
                s,
                "pubchem_time_span: {}",
                fmt_time_span(pubchem_time_span.as_ref())
            )
            .ok();
            writeln!(
                s,
                "kegg_conv_time_span: {}",
                fmt_time_span(kegg_conv_time_span.as_ref())
            )
            .ok();
            writeln!(
                s,
                "module_retention_present: {}",
                module_retention.is_some()
            )
            .ok();
            writeln!(
                s,
                "dual_mode_breakdown_present: {}",
                dual_mode_breakdown.is_some()
            )
            .ok();
            writeln!(s, "rendering: {rendering}").ok();
            writeln!(s, "refresh_state: {refresh_state:?}").ok();
            write_dam_results(&mut s, dam_results);
        }
    }
    s
}

fn opt_str(o: Option<&str>) -> String {
    o.map(|s| s.to_string()).unwrap_or_else(|| "<none>".into())
}

fn opt_dbg<T: std::fmt::Debug>(o: Option<&T>) -> String {
    o.map(|v| format!("{v:?}"))
        .unwrap_or_else(|| "<none>".into())
}

fn write_ion_tables(s: &mut String, ion_tables: &[IonModeTable]) {
    if ion_tables.is_empty() {
        writeln!(s, "ion_tables: <no table loaded>").ok();
        return;
    }
    writeln!(s, "ion_tables_count: {}", ion_tables.len()).ok();
    for (i, t) in ion_tables.iter().enumerate() {
        let path = t
            .txt_path
            .as_ref()
            .map(|p| p.display().to_string())
            .unwrap_or_else(|| "<no path>".into());
        writeln!(
            s,
            "  ion_table[{i}]: mode={:?} features={} samples={} path={}",
            t.mode,
            t.table.features.len(),
            t.table.sample_cols.len(),
            path
        )
        .ok();
    }
}

fn write_mapping(s: &mut String, mapping: Option<&GroupMapping>) {
    match mapping {
        None => {
            writeln!(s, "mapping: <no mapping loaded>").ok();
        }
        Some(m) => {
            writeln!(s, "mapping_sample_count: {}", m.sample_count()).ok();
            writeln!(s, "mapping_has_biosample: {}", m.has_biosample()).ok();
            let groups = m.groups_in_order();
            let mut group_counts: HashMap<String, usize> = HashMap::new();
            for (_, g) in &groups {
                *group_counts.entry(g.clone()).or_insert(0) += 1;
            }
            for (g, c) in &group_counts {
                writeln!(s, "  group[{g}]: {c} samples").ok();
            }
            let cols = m.metadata_column_names();
            writeln!(s, "mapping_metadata_columns: {cols:?}").ok();
        }
    }
}

fn write_dam_results(s: &mut String, dam_results: &[crate::dam::DamResult]) {
    writeln!(s, "dam_results_count: {}", dam_results.len()).ok();
    for (i, dr) in dam_results.iter().enumerate() {
        writeln!(
            s,
            "  dam_result[{i}]: method={:?} num={} den={} features={} skipped={} fdr_method={:?}",
            dr.method,
            dr.numerator,
            dr.denominator,
            dr.features.len(),
            dr.skipped,
            dr.fdr_method
        )
        .ok();
    }
}

fn fmt_time_span(
    span: Option<&(
        chrono::DateTime<chrono::Utc>,
        chrono::DateTime<chrono::Utc>,
        usize,
    )>,
) -> String {
    match span {
        None => "<n/a>".into(),
        Some((min, max, n)) => format!("min={min} max={max} entries={n}"),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Drift-guard (`harden-diagnostics-redaction` D2): the `[settings]` block
    /// is a hand-written `writeln!` per `SessionSettings` field, with no
    /// compile-time completeness check (the struct is `Serialize`, not
    /// pattern-matched). Adding a field would silently drop it from every bug
    /// report. This test fails loudly if a serde field has no rendered label.
    #[test]
    fn snapshot_renders_every_session_settings_field() {
        let settings = SessionSettings::default();
        let json = serde_json::to_value(&settings).expect("SessionSettings serialises");
        let obj = json.as_object().expect("SessionSettings is a flat object");
        let inputs = SessionInputs::default();
        let cache = SessionCache::default();
        let rendered = render_app_state(&AppState::default(), &settings, &inputs, &cache);

        // The six export-size fields render under two composite labels rather
        // than a per-field `<name>:` line. Every other field MUST appear as a
        // literal `<name>:` label, so a genuinely new field fails loudly.
        const COMPOSITE: &[(&str, &str)] = &[
            ("stage2_export_width_in", "stage2_export:"),
            ("stage2_export_height_in", "stage2_export:"),
            ("stage2_export_dpi", "stage2_export:"),
            ("stage3_export_width_in", "stage3_export:"),
            ("stage3_export_height_in", "stage3_export:"),
            ("stage3_export_dpi", "stage3_export:"),
        ];
        for key in obj.keys() {
            let composite = COMPOSITE.iter().find(|(f, _)| f == key).map(|(_, a)| *a);
            let owned = format!("{key}:");
            let needle = composite.unwrap_or(owned.as_str());
            // Line-anchored (a trimmed line must START with the label) so e.g.
            // `fdr_threshold:` cannot spuriously match `enrichment_fdr_threshold:`.
            assert!(
                rendered.lines().any(|l| l.trim_start().starts_with(needle)),
                "SessionSettings field `{key}` is missing from app_state.txt [settings] block \
                 — add a writeln! for it in render_app_state (or, if it is a composite export \
                 field, extend the COMPOSITE allowlist)"
            );
        }
    }
}
