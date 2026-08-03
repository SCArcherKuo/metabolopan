//! Stage 3 orchestrator: drives PubChem → KEGG conv → ORA in sequence,
//! emitting progress events to the running screen and assembling the
//! `Stage3RunOutput` consumed by `Stage3EnrichResult`.
//!
//! See the `stage3-ui` and `enrichment-ora` capabilities for the
//! end-to-end contract.

use anyhow::Result;
use chrono::{DateTime, Utc};
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::sync::mpsc;
use tracing::{debug, info};

use crate::app::{AnalysisMode, AnalysisPayload, CoverageFunnel, Stage3Funnel, Stage3RunOutput};
use crate::coverage::CoverageResult;
use crate::dam::fdr::FdrMethod;
use crate::dam::run::classify_trend;
use crate::dam::{DamMethod, DamResult};
use crate::data::{GroupMapping, IonModeTable};
use crate::dedup::DedupReport;
use crate::enrichment::run_ora;
use crate::enrichment::types::{EnrichmentDirection, EnrichmentResult};
use crate::kegg::{
    ConvProgress, KeggClient, KeggCompoundSet, KeggModulesCache, resolve_cids_to_cpds,
};
use crate::kegg::{cache as kegg_cache, types::CidCpdEntry};
use crate::pubchem::cache as pubchem_cache;
use crate::pubchem::types::InchikeyCidsEntry;
use crate::pubchem::{PubchemClient, PubchemProgress, resolve_inchikeys_to_cids};

/// Helper: log a sample of compound IDs (first N) without dumping the
/// entire set into the log pane. Used by Stage 3 diagnostics.
fn fmt_cpd_sample(cpds: impl Iterator<Item = impl AsRef<str>>, max: usize) -> String {
    let mut v: Vec<String> = cpds.take(max).map(|s| s.as_ref().to_string()).collect();
    v.sort();
    v.join(", ")
}

/// `AnalysisTarget` is the orchestrator-facing alias for `AnalysisPayload`
/// — they carry identical mode-specific data; the alias documents intent
/// when threaded into `run_stage3`. Pathway-mode runs operate on a
/// `SpeciesKegg`'s pathways; module-mode runs filter a `KeggModulesCache`
/// by Group overlap, then map retained modules into `KeggCompoundSet`s.
pub type AnalysisTarget = AnalysisPayload;

/// Retention summary surfaced on `Stage3RunOutput.module_retention` in
/// module mode (`None` in pathway mode). Drives the Stage 3 result
/// panel's mode-specific "Data sources for this run" copy.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct ModuleRetention {
    /// Total modules in the cache (i.e. `/list/module` count at fetch time).
    pub total_modules: usize,
    /// Modules surviving the `|complete_orgs ∩ group_orgs| >= min_group_overlap` filter.
    pub retained_modules: usize,
    /// Active threshold used for the Group filter. Default 1 (permissive ∃).
    pub min_group_overlap: usize,
    /// Selected Level (1..=3). Validated at construction.
    pub group_level: u8,
    /// Selected Group name (e.g. "Animals").
    pub group_name: String,
    /// `|group_orgs|` — count of organism codes in the selected Group.
    pub group_org_count: usize,
    /// Min `fetched_at` across retained modules. Used by the result
    /// panel's "Modules cache time span" copy, parallel to the
    /// `pubchem_time_span` / `kegg_conv_time_span` tuples.
    pub oldest_fetched_at: DateTime<Utc>,
    pub newest_fetched_at: DateTime<Utc>,
}

/// Per-mode partition counts surfaced on `Stage3RunOutput.dual_mode_breakdown`
/// in dual-mode runs (`None` in single-mode). Drives the Stage 3 result
/// panel's dual-mode breakdown block. The arithmetic SHALL satisfy:
/// `|N| == universe_pos_only + universe_neg_only + universe_in_both` and
/// `|K| == foreground_pos_only + foreground_neg_only + foreground_agree_both`.
/// `foreground_excluded_conflict` is reported separately as a diagnostic
/// count (cpds dropped by the conflict-only-strict union rule are by
/// definition not in K).
#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct DualModeBreakdown {
    pub universe_pos_only: usize,
    pub universe_neg_only: usize,
    pub universe_in_both: usize,
    pub foreground_pos_only: usize,
    pub foreground_neg_only: usize,
    pub foreground_agree_both: usize,
    pub foreground_excluded_conflict: usize,
}

impl DualModeBreakdown {
    /// Tally a cpd into the universe (N) partition by which modes reach it.
    /// Moved verbatim from `build_dam_cpd_dual`'s former in-loop counter.
    fn tally_universe(&mut self, pos: ModeTrend, neg: ModeTrend) {
        match (pos != ModeTrend::Absent, neg != ModeTrend::Absent) {
            (true, false) => self.universe_pos_only += 1,
            (false, true) => self.universe_neg_only += 1,
            (true, true) => self.universe_in_both += 1,
            (false, false) => {}
        }
    }

    /// Tally an in-K cpd into the foreground partition by which modes flag it
    /// significant (Up/Down). Moved verbatim from the former in-loop counter.
    fn tally_in_k(&mut self, pos: ModeTrend, neg: ModeTrend) {
        let pos_sig = matches!(pos, ModeTrend::Up | ModeTrend::Down);
        let neg_sig = matches!(neg, ModeTrend::Up | ModeTrend::Down);
        match (pos_sig, neg_sig) {
            (true, false) => self.foreground_pos_only += 1,
            (false, true) => self.foreground_neg_only += 1,
            (true, true) => self.foreground_agree_both += 1,
            (false, false) => {}
        }
    }
}

/// Parameters carried into the orchestrator from `Stage3EnrichRunning`.
#[derive(Debug, Clone, Copy)]
pub struct Stage3Params {
    pub method: DamMethod,
    pub fc_threshold: f64,
    pub fdr_threshold: f64,
    pub delta_threshold: f64,
    pub direction: EnrichmentDirection,
    pub min_hit_count: usize,
    /// Pre-FDR entry-size threshold passed through to `run_ora`. Entries
    /// with `m_p < min_entry_size` are dropped before FDR.
    pub min_entry_size: usize,
    pub fdr_method: FdrMethod,
    pub force_refresh_pubchem: bool,
    pub force_refresh_kegg_conv: bool,
}

/// Output of [`resolve_detected_compounds`] — the InChIKey → CID → KEGG-cpd
/// chain both analysis routes share.
///
/// Carries `all_cids` rather than a bare count because callers need the list
/// itself: `compute_kegg_conv_time_span` reads the cache entry per CID to
/// derive the fetched-date span the Data tab renders.
///
/// The two provenance counts are DERIVED, not stored, so they cannot drift
/// from their sources:
/// - detected InChIKeys = `inchikey_to_cids.len()` — the resolved map's key
///   count, NOT the input slice's length. Single-mode callers pass an
///   undeduped per-feature list, so its length would overcount.
/// - detected CIDs = `all_cids.len()` — unique CIDs submitted to KEGG `/conv`.
pub struct ResolvedCompounds {
    pub inchikey_to_cids: HashMap<String, Vec<String>>,
    pub cid_to_cpd: HashMap<String, Option<String>>,
    /// Sorted, deduplicated CIDs across every resolved InChIKey — the exact
    /// slice handed to `resolve_cids_to_cpds`.
    pub all_cids: Vec<String>,
}

/// Resolve a list of InChIKeys to KEGG compound IDs: PubChem InChIKey → CID,
/// then KEGG `/conv` CID → cpd.
///
/// Extracted from `run_stage3` so the enrichment and coverage routes share one
/// implementation of the slow, network-bound part of the pipeline — batching,
/// retry, cache read/write, force-refresh, and progress emission must not drift
/// between them (see the `stage3-ui` capability spec).
///
/// **Callers do not deduplicate.** `resolve_inchikeys_to_cids` and
/// `resolve_cids_to_cpds` both dedupe internally via
/// `crate::seq::dedupe_preserve_order`, so the input contract is "a list which
/// may repeat" and is identical for both routes. The enrichment route therefore
/// keeps passing `collect_inchikeys`'s output unchanged, including its
/// deliberately-undeduped single-mode branch.
pub async fn resolve_detected_compounds(
    pubchem_client: &PubchemClient,
    kegg_client: &KeggClient,
    inchikeys: &[String],
    force_refresh_pubchem: bool,
    force_refresh_kegg_conv: bool,
    pubchem_progress_tx: mpsc::Sender<PubchemProgress>,
    kegg_conv_progress_tx: mpsc::Sender<ConvProgress>,
) -> Result<ResolvedCompounds> {
    let inchikey_to_cids = resolve_inchikeys_to_cids(
        pubchem_client,
        inchikeys,
        force_refresh_pubchem,
        Some(pubchem_progress_tx),
    )
    .await?;

    // ── Phase 2: collect unique CIDs and resolve → cpd IDs ──
    let mut all_cids: Vec<String> = inchikey_to_cids
        .values()
        .flat_map(|v| v.iter().cloned())
        .collect();
    all_cids.sort();
    all_cids.dedup();

    let inchikeys_with_cids = inchikey_to_cids.values().filter(|v| !v.is_empty()).count();
    info!(
        inchikeys_resolved_to_any_cid = inchikeys_with_cids,
        unique_cids_to_lookup = all_cids.len(),
        "Phase 1 PubChem resolution complete — Phase 2 KEGG /conv begins"
    );

    let cid_to_cpd = resolve_cids_to_cpds(
        kegg_client,
        &all_cids,
        force_refresh_kegg_conv,
        Some(kegg_conv_progress_tx),
    )
    .await?;

    Ok(ResolvedCompounds {
        inchikey_to_cids,
        cid_to_cpd,
        all_cids,
    })
}

/// Run the full Stage 3 pipeline. Pure async function — no UI awareness.
/// Errors abort the Run; previously-written cache entries persist.
///
/// `target` selects the entry catalogue ORA tests against:
/// - `AnalysisTarget::Pathway` (pathway mode): every pathway in the
///   selected `SpeciesKegg`.
/// - `AnalysisTarget::Module` (module mode): modules from the
///   `KeggModulesCache` whose `complete_orgs` overlap the Group's
///   `org_codes` set by at least `min_group_overlap`.
///
/// The PubChem and KEGG conv phases are mode-agnostic — they always
/// produce the measurable-metabolome universe `N`.
pub async fn run_stage3(
    pubchem_client: &PubchemClient,
    kegg_client: &KeggClient,
    dam_results: &[DamResult],
    target: &AnalysisTarget,
    params: Stage3Params,
    pubchem_progress_tx: mpsc::Sender<PubchemProgress>,
    kegg_conv_progress_tx: mpsc::Sender<ConvProgress>,
) -> Result<Stage3RunOutput> {
    if dam_results.is_empty() {
        anyhow::bail!("run_stage3 requires at least one DamResult");
    }
    // Lock the 1-or-2-mode contract one layer earlier than `build_dam_cpd_from_trends`'s
    // own debug_assert. The orchestrator's DualModeBreakdown partition is hardwired
    // 2-way POS/NEG; anything beyond len()=2 would silently undercount. Catches
    // future UI regressions (`fix-stage3-ui-dual-mode-spawn` Finding §5).
    debug_assert!(
        dam_results.len() <= 2,
        "run_stage3: dam_results.len() = {} exceeds the 1-or-2 mode contract",
        dam_results.len()
    );
    let n_modes = dam_results.len();
    let is_dual = n_modes >= 2;

    // ── Phase 1: resolve every annotated DAM feature's InChIKey → CIDs ──
    // Single-mode: preserve the prior ordering (no dedup) so behavior is
    // bit-equal to single-mode regression baselines. Dual-mode: union across
    // modes + sort + dedup so the single PubChem call covers both modes
    // without duplicate work (cache is keyed per InChIKey).
    let inchikeys: Vec<String> = collect_inchikeys(dam_results, is_dual);

    info!(
        mode = ?target.mode(),
        n_modes,
        dam_features_total = dam_results.iter().map(|d| d.features.len()).sum::<usize>(),
        inchikeys_to_resolve = inchikeys.len(),
        direction = ?params.direction,
        fc_threshold = params.fc_threshold,
        fdr_threshold = params.fdr_threshold,
        delta_threshold = params.delta_threshold,
        min_hit_count = params.min_hit_count,
        "Stage 3 Run starting — Phase 1 PubChem resolution begins"
    );

    // Phases 1 + 2 are identical on both analysis routes, so they live in the
    // shared resolver rather than here (see the `stage3-ui` capability spec).
    let ResolvedCompounds {
        inchikey_to_cids,
        cid_to_cpd,
        all_cids,
    } = resolve_detected_compounds(
        pubchem_client,
        kegg_client,
        &inchikeys,
        params.force_refresh_pubchem,
        params.force_refresh_kegg_conv,
        pubchem_progress_tx,
        kegg_conv_progress_tx,
    )
    .await?;

    // ── Build per-mode feature_to_cpds maps (multi-mapping rule D8) + the
    //    measurable-metabolome universe N (union of mapped cpd sets, D1).
    //    Single-mode is the degenerate length-1 case, bit-equal to pre-7.3. ──
    let (per_mode_feature_to_cpds, mapped_universe) =
        build_per_mode_maps(dam_results, &inchikey_to_cids, &cid_to_cpd);

    let cids_with_cpd = cid_to_cpd.values().filter(|v| v.is_some()).count();
    info!(
        cids_resolved_to_cpd = cids_with_cpd,
        cids_no_cpd_mapping = all_cids.len() - cids_with_cpd,
        features_mapped_to_any_cpd = per_mode_feature_to_cpds
            .iter()
            .map(|m| m.len())
            .sum::<usize>(),
        universe_size = mapped_universe.len(),
        "Phase 2 KEGG /conv complete — measurable-metabolome universe N built"
    );

    // ── Build K + dual-mode breakdown ──
    // K-build runs through the SAME conflict-only-strict helper for single- and
    // dual-mode: single-mode is the degenerate length-1 case, so it now ALSO
    // applies the conflict rule — a compound reached by both an Up and a Down
    // feature (intra-mode conflict) is excluded from K. This is a deliberate
    // consistency fix (refine-stage3-dual-mode-internals): pre-change single-mode
    // `build_dam_cpd` kept such compounds. The pos/neg `DualModeBreakdown` is a
    // dual-mode-only artifact, surfaced (`Some`) only when `is_dual`; single-mode
    // keeps `dual_mode_breakdown = None`.
    // Classify every feature ONCE (per mode, indexed by feature position) and
    // reuse it for both the K-build and the foreground funnel — removing the
    // redundant second `classify_trend` pass the funnel used to run
    // (refine-stage3-dual-mode-internals D2).
    let per_mode_trends: Vec<Vec<crate::dam::types::Trend>> = dam_results
        .iter()
        .map(|dam| {
            dam.features
                .iter()
                .map(|f| {
                    classify_trend(
                        f,
                        params.fc_threshold,
                        params.fdr_threshold,
                        params.delta_threshold,
                        params.method,
                    )
                })
                .collect()
        })
        .collect();

    let (dam_cpd, breakdown, conflict_sample) = build_dam_cpd_from_trends(
        dam_results,
        &per_mode_feature_to_cpds,
        &per_mode_trends,
        &mapped_universe,
        params.direction,
    );
    let n_excluded_conflict = breakdown.foreground_excluded_conflict;
    let dual_mode_breakdown = if is_dual { Some(breakdown) } else { None };

    // ── Assemble entries based on mode (Track E) ──
    let (entries, module_retention): (Vec<KeggCompoundSet>, Option<ModuleRetention>) = match target
    {
        AnalysisTarget::Pathway { species_kegg } => (species_kegg.pathways.clone(), None),
        AnalysisTarget::Module {
            modules_pack,
            group_level,
            group_name,
            group_org_codes,
            min_group_overlap,
        } => assemble_module_entries(
            modules_pack,
            *group_level,
            group_name,
            group_org_codes,
            *min_group_overlap,
        ),
    };

    // Emit the K / universe / conflict diagnostics together with the
    // entry-assembly summary, now that `entries` is built. The K logs precede
    // the entry logs (identical global order to before — entry assembly itself
    // emits no events).
    log_k_diagnostics(
        &dam_cpd,
        dual_mode_breakdown.as_ref(),
        n_excluded_conflict,
        &conflict_sample,
        &entries,
        &module_retention,
    );

    // ── Phase 3: ORA (synchronous; instant) ──
    let entries_total = entries.len();
    let enrichment_result = run_ora(
        &mapped_universe,
        &dam_cpd,
        &entries,
        params.min_hit_count,
        params.direction,
        params.fdr_method,
        params.min_entry_size,
    );
    info!(
        min_entry_size = params.min_entry_size,
        entries_total,
        entries_dropped_by_min_entry_size = enrichment_result.entries_dropped_by_min_entry_size,
        entries_tested = enrichment_result.rows.len(),
        "Pre-FDR min_entry_size filter applied"
    );

    // ── Diagnostic: which K cpds appear in any entry vs none? ──
    let mut k_covered: HashSet<&String> = HashSet::new();
    for row in &enrichment_result.rows {
        for hit in &row.hit_kegg_ids {
            k_covered.insert(hit);
        }
    }
    let k_uncovered: Vec<&String> = dam_cpd.iter().filter(|c| !k_covered.contains(c)).collect();
    log_ora_diagnostics(
        &enrichment_result,
        &dam_cpd,
        target,
        &k_covered,
        &k_uncovered,
    );

    // ── Compute time spans from caches (D16) ──
    let pubchem_time_span =
        compute_pubchem_time_span(&pubchem_cache::read_cache().unwrap_or_default(), &inchikeys);
    let kegg_conv_time_span = compute_kegg_conv_time_span(
        &kegg_cache::read_cid_to_cpd_cache().unwrap_or_default(),
        &all_cids,
    );

    // For Stage3RunOutput.feature_to_cpds (a single map): single-mode uses
    // per_mode_feature_to_cpds[0] verbatim (bit-equal). Dual-mode merges the
    // per-mode maps by inchikey, taking the cpd-set union.
    let feature_to_cpds: HashMap<String, HashSet<String>> = if is_dual {
        merge_feature_to_cpds(&per_mode_feature_to_cpds)
    } else {
        per_mode_feature_to_cpds
            .into_iter()
            .next()
            .unwrap_or_default()
    };

    // ── Provenance funnel counts for the Data tab (read-only by-products of
    //    the universe/foreground construction; no new network calls). ──
    let funnel = build_funnel(
        dam_results,
        &per_mode_trends,
        &inchikey_to_cids,
        &all_cids,
        &params,
        k_covered.len(),
    );
    // Funnel monotonicity (mirrors the K ⊆ N invariant pattern): loud in
    // dev/test, never panics in release (the counts are presentation-only).
    debug_assert!(
        funnel.detected_inchikeys >= funnel.foreground_inchikeys,
        "funnel: detected_inchikeys {} < foreground_inchikeys {}",
        funnel.detected_inchikeys,
        funnel.foreground_inchikeys
    );
    debug_assert!(
        funnel.detected_cids >= funnel.foreground_cids,
        "funnel: detected_cids {} < foreground_cids {}",
        funnel.detected_cids,
        funnel.foreground_cids
    );
    debug_assert!(
        mapped_universe.len() >= dam_cpd.len(),
        "funnel: universe {} < K {}",
        mapped_universe.len(),
        dam_cpd.len()
    );

    Ok(Stage3RunOutput {
        enrichment_result,
        mapped_universe,
        feature_to_cpds,
        pubchem_time_span,
        kegg_conv_time_span,
        module_retention,
        dual_mode_breakdown,
        funnel,
    })
}

/// The loaded session inputs a coverage run reads.
///
/// Bundled rather than passed as two parameters; they also travel together
/// everywhere — the mapping is only ever interpreted against these tables'
/// sample columns.
#[derive(Debug, Clone, Copy)]
pub struct CoverageInputs<'a> {
    pub ion_tables: &'a [IonModeTable],
    /// `None` when no metadata `.csv` was supplied — fully supported on this
    /// route, and the reason the group-presence stage is optional.
    pub mapping: Option<&'a GroupMapping>,
}

/// Per-run parameters for [`run_coverage`].
///
/// Deliberately NOT a reuse of [`Stage3Params`]: that struct is more than half
/// statistics (`method`, `fc_threshold`, `fdr_threshold`, `delta_threshold`,
/// `direction`, `fdr_method`), none of which exists on a route that performs no
/// test. Sharing it would put six permanently-unread fields on every coverage
/// run and invite exactly the confusion the route's `no statistical test`
/// guarantee exists to prevent.
///
/// The two display filters (`min_hit_count`, `top_n`) are absent for a
/// different reason: they are re-applied live on the result screen over the
/// rows already in hand, so the run never reads them.
#[derive(Debug, Clone, PartialEq)]
pub struct CoverageParams {
    /// `settings.coverage_selected_groups`, verbatim. `None` ("not yet chosen")
    /// and `Some(vec![])` ("deliberately none") are DIFFERENT and are resolved
    /// by `coverage::detect::selected_groups`, never by `unwrap_or_default`.
    pub selected_groups: Option<Vec<String>>,
    pub presence_threshold: f64,
    pub dedup_enabled: bool,
    pub dedup_rt_tolerance_min: f64,
    pub force_refresh_pubchem: bool,
    pub force_refresh_kegg_conv: bool,
}

/// Everything a coverage run needs from the loaded tables, extracted BEFORE the
/// orchestrator is spawned.
///
/// The split is here because every table-dependent step of a coverage run — the
/// group-presence filter, deduplication, and reading each surviving feature's
/// metabolite name — is pure, synchronous, and fast, while the only slow step
/// is the network resolver, which never touches a table. Handing the spawned
/// task this struct instead of the tables means the async future never holds
/// the intensity matrices: no multi-megabyte clone per run, and
/// `MetabolomicsTable` does not have to become `Clone` (which would make an
/// expensive copy easy to write by accident anywhere else in the app).
#[derive(Debug, Clone)]
pub struct PreparedFeatures {
    /// Distinct InChIKeys to resolve — the union across modes, in
    /// first-appearance order.
    pub inchikeys: Vec<String>,
    /// Per ion-mode table, `(inchikey, metabolite_name)` for each SURVIVING
    /// annotated feature.
    ///
    /// Per table rather than unioned because one InChIKey may carry different
    /// MS-DIAL names in POS and NEG, and the CSV lists every distinct one. Only
    /// the survivors appear, which is what makes deduplication observable in
    /// the exported names — its one effect on this route.
    pub per_mode_annotations: Vec<Vec<(String, String)>>,
    pub raw_features: usize,
    /// `None` when no metadata `.csv` was supplied — the stage did not run.
    pub in_selected_groups: Option<usize>,
    pub after_dedup: usize,
    /// One per ion-mode table, in order; empty when dedup was off.
    pub dedup_reports: Vec<DedupReport>,
}

/// Apply the group-presence filter and deduplication to every loaded table and
/// extract what the run needs from them.
///
/// Synchronous and cheap — call it on the UI thread immediately before
/// spawning. Emits one count-only `info!` per table for the group filter: no
/// sample names and no group names beyond those already on screen, so the event
/// is safe inside a bug-report bundle. Every other exclusion on this route is
/// surfaced in both the UI and the log; this one matches.
pub fn prepare_features(inputs: CoverageInputs<'_>, params: &CoverageParams) -> PreparedFeatures {
    let CoverageInputs {
        ion_tables,
        mapping,
    } = inputs;

    // Resolve the selection ONCE, against the mapping, so every table applies
    // the same group list. With no mapping there is nothing to resolve against
    // and the filter is inert anyway.
    let groups: Vec<String> = match mapping {
        Some(m) => crate::coverage::detect::selected_groups(params.selected_groups.as_deref(), m),
        None => Vec::new(),
    };

    let mut raw_features = 0usize;
    let mut in_selected_groups: Option<usize> = None;
    let mut after_dedup = 0usize;
    let mut dedup_reports: Vec<DedupReport> = Vec::new();
    let mut per_mode_annotations: Vec<Vec<(String, String)>> = Vec::new();

    for table in ion_tables {
        let (detected, report) = crate::coverage::detect::detect_features(
            &table.table,
            mapping,
            &groups,
            params.presence_threshold,
            params.dedup_enabled,
            params.dedup_rt_tolerance_min,
        );
        if let Some(surviving) = detected.in_selected_groups {
            info!(
                mode = ?table.mode,
                features_total = detected.raw_features,
                removed_by_group_filter = detected.raw_features - surviving,
                surviving,
                "coverage group-presence filter applied"
            );
            *in_selected_groups.get_or_insert(0) += surviving;
        }
        raw_features += detected.raw_features;
        after_dedup += detected.after_dedup;
        if let Some(r) = report {
            dedup_reports.push(r);
        }
        per_mode_annotations.push(
            detected
                .kept
                .iter()
                .filter_map(|&i| {
                    let f = &table.table.features[i];
                    f.inchikey
                        .as_ref()
                        .map(|k| (k.clone(), f.metabolite_name.clone()))
                })
                .collect(),
        );
    }

    let all_keys: Vec<String> = per_mode_annotations
        .iter()
        .flatten()
        .map(|(k, _)| k.clone())
        .collect();
    let inchikeys = crate::seq::dedupe_preserve_order(&all_keys);

    PreparedFeatures {
        inchikeys,
        per_mode_annotations,
        raw_features,
        in_selected_groups,
        after_dedup,
        dedup_reports,
    }
}

/// Which ionization modes a detected compound was reached through.
///
/// Computed for the Data tab only — the Data tab is its sole surface. It MUST
/// NOT affect membership in `D`, which is the plain union `D_pos ∪ D_neg`:
/// with no differential comparison there is no directional verdict that could
/// contradict another, so this route has no conflict rule to apply.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct CoverageModePartition {
    pub pos_only: usize,
    pub neg_only: usize,
    pub in_both: usize,
}

/// Output of [`run_coverage`] — the coverage route's counterpart to
/// [`Stage3RunOutput`].
pub struct CoverageRunOutput {
    pub coverage_result: CoverageResult,
    pub funnel: CoverageFunnel,
    /// cpd ID → the user's MS-DIAL metabolite names for the features that
    /// resolved to it, sorted and deduplicated. Consumed only by the CSV
    /// exporter (`CoverageExportContext`); the on-screen table renders bare
    /// cpd IDs, which is why `coverage::compute` never sees this map.
    pub cpd_to_names: HashMap<String, Vec<String>>,
    /// Module mode only, exactly as on the enrichment route.
    pub module_retention: Option<ModuleRetention>,
    /// Pathway mode only: the species code this run was performed against,
    /// captured from the run's own target. `None` in Module mode, where the
    /// target components already ride on `module_retention`.
    pub target_species: Option<String>,
    /// `None` in single-mode runs.
    pub mode_partition: Option<CoverageModePartition>,
    /// One `DedupReport` per ion-mode table, in `ion_tables` order. Empty when
    /// `dedup_enabled` is false.
    pub dedup_reports: Vec<DedupReport>,
    pub pubchem_time_span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>,
    pub kegg_conv_time_span: Option<(DateTime<Utc>, DateTime<Utc>, usize)>,
}

/// Run the coverage survey: resolve the prepared feature set to KEGG compounds
/// through the SHARED resolver, assemble the entry catalogue from the SHARED
/// `AnalysisTarget`, and compute descriptive coverage.
///
/// No statistical test is performed and none can be: `CoverageResult` has no
/// field to carry one.
///
/// Takes [`PreparedFeatures`] rather than the tables — see that type for why.
/// `force_refresh` is `(pubchem, kegg_conv)`, paired so the signature stays
/// inside clippy's argument limit; the two flags are always set together.
pub async fn run_coverage(
    pubchem_client: &PubchemClient,
    kegg_client: &KeggClient,
    prepared: PreparedFeatures,
    target: &AnalysisTarget,
    force_refresh: (bool, bool),
    pubchem_progress_tx: mpsc::Sender<PubchemProgress>,
    kegg_conv_progress_tx: mpsc::Sender<ConvProgress>,
) -> Result<CoverageRunOutput> {
    let (force_refresh_pubchem, force_refresh_kegg_conv) = force_refresh;
    let PreparedFeatures {
        inchikeys,
        per_mode_annotations,
        raw_features,
        in_selected_groups,
        after_dedup,
        dedup_reports,
    } = prepared;

    info!(
        mode = ?target.mode(),
        n_modes = per_mode_annotations.len(),
        raw_features,
        after_dedup,
        inchikeys_to_resolve = inchikeys.len(),
        "Coverage run starting"
    );

    // ── The shared resolver: InChIKey → CID → KEGG cpd ──
    let resolved = resolve_detected_compounds(
        pubchem_client,
        kegg_client,
        &inchikeys,
        force_refresh_pubchem,
        force_refresh_kegg_conv,
        pubchem_progress_tx,
        kegg_conv_progress_tx,
    )
    .await?;

    // ── D, and the cpd → metabolite-name map the CSV needs ──
    let cpds_of = |key: &str| -> Vec<String> {
        resolved
            .inchikey_to_cids
            .get(key)
            .into_iter()
            .flatten()
            .filter_map(|cid| resolved.cid_to_cpd.get(cid).and_then(|c| c.clone()))
            .collect()
    };
    let (detected_cpds, cpd_to_names) = build_detected_and_names(&per_mode_annotations, cpds_of);

    // Per-mode partition — Data tab only, never a membership rule.
    let per_mode_cpds: Vec<HashSet<String>> = per_mode_annotations
        .iter()
        .map(|anns| anns.iter().flat_map(|(k, _)| cpds_of(k)).collect())
        .collect();
    let mode_partition = partition_by_mode(&per_mode_cpds);
    debug_assert_eq!(
        per_mode_cpds.iter().flatten().collect::<HashSet<_>>().len(),
        detected_cpds.len(),
        "D must be the plain union of the per-mode compound sets"
    );

    // ── Entries: the SAME catalogue assembly the enrichment route uses ──
    let (entries, module_retention): (Vec<KeggCompoundSet>, Option<ModuleRetention>) = match target
    {
        AnalysisTarget::Pathway { species_kegg } => (species_kegg.pathways.clone(), None),
        AnalysisTarget::Module {
            modules_pack,
            group_level,
            group_name,
            group_org_codes,
            min_group_overlap,
        } => assemble_module_entries(
            modules_pack,
            *group_level,
            group_name,
            group_org_codes,
            *min_group_overlap,
        ),
    };

    let coverage_result = crate::coverage::compute(&detected_cpds, &entries);

    let funnel = CoverageFunnel {
        raw_features,
        in_selected_groups,
        after_dedup,
        // The resolved map's key count, not the input slice's length — the same
        // derivation the enrichment funnel uses.
        detected_inchikeys: resolved.inchikey_to_cids.len(),
        detected_cids: resolved.all_cids.len(),
    };

    let pubchem_cache_map = pubchem_cache::read_cache().unwrap_or_default();
    let kegg_conv_cache_map = kegg_cache::read_cid_to_cpd_cache().unwrap_or_default();
    let pubchem_time_span = compute_pubchem_time_span(&pubchem_cache_map, &inchikeys);
    let kegg_conv_time_span = compute_kegg_conv_time_span(&kegg_conv_cache_map, &resolved.all_cids);

    info!(
        entries_total = coverage_result.entries_total,
        entries_without_compounds = coverage_result.entries_without_compounds,
        detected_total = coverage_result.detected_total,
        detected_in_entries = coverage_result.detected_in_entries,
        "Coverage run complete"
    );

    Ok(CoverageRunOutput {
        coverage_result,
        funnel,
        cpd_to_names,
        module_retention,
        target_species: captured_species(target),
        mode_partition,
        dedup_reports,
        pubchem_time_span,
        kegg_conv_time_span,
    })
}

/// The species code a coverage run was performed against, taken from the run's
/// own `AnalysisTarget`. `None` in Module mode, where `ModuleRetention` already
/// carries the target's components.
///
/// A free function rather than an inline `match` so it is unit-testable without
/// a runtime: `run_coverage` is `async` and takes two live clients, and this is
/// the one assertion that distinguishes capturing from the run's target from
/// capturing from `app.settings` when the result state is built. Note this
/// function — like `run_coverage` itself — has no access to `SessionSettings`,
/// so the capture cannot silently fall back to them.
pub(crate) fn captured_species(target: &AnalysisTarget) -> Option<String> {
    target.pathway_species().map(|sk| sk.code.clone())
}

impl CoverageRunOutput {
    /// Consume the orchestrator output into a fully-populated
    /// `AppState::Stage3CoverageResult`. The result-screen runtime fields start
    /// fresh, mirroring `Stage3RunOutput::into_result_state`.
    pub(crate) fn into_result_state(self) -> crate::app::AppState {
        crate::app::AppState::Stage3CoverageResult {
            coverage_result: self.coverage_result,
            funnel: self.funnel,
            cpd_to_names: self.cpd_to_names,
            module_retention: self.module_retention,
            target_species: self.target_species,
            mode_partition: self.mode_partition,
            dedup_reports: self.dedup_reports,
            pubchem_time_span: self.pubchem_time_span,
            kegg_conv_time_span: self.kegg_conv_time_span,
            dotplot_tex: None,
            rendering: false,
            render_rx: None,
            confirming_new_round: false,
            // Fresh per-run state: the run-entry autosize is authoritative
            // until the user hand-edits the Height field on this screen.
            height_user_overridden: false,
        }
    }
}

/// Build `D` and the cpd → metabolite-name map from the SURVIVING features.
///
/// `per_mode_annotations` already contains only the features that survived both
/// filters (see [`PreparedFeatures`]), which is precisely how deduplication
/// earns its place on this route: it elects a different representative per
/// InChIKey, so `C00031 (D-Glucose / Glucose)` becomes `C00031 (D-Glucose)`.
/// Including the dup-losers would make the dedup control inert in its ONE
/// observable effect — `D` and every number derived from it are invariant under
/// it.
///
/// Per table rather than over the unioned InChIKey list, because one InChIKey
/// may carry different MS-DIAL names in POS and NEG and the CSV lists every
/// distinct one.
///
/// Names come back sorted and deduplicated, so the CSV cell is deterministic.
fn build_detected_and_names(
    per_mode_annotations: &[Vec<(String, String)>],
    cpds_of: impl Fn(&str) -> Vec<String>,
) -> (HashSet<String>, HashMap<String, Vec<String>>) {
    let mut detected: HashSet<String> = HashSet::new();
    let mut names: HashMap<String, HashSet<String>> = HashMap::new();
    for anns in per_mode_annotations {
        for (key, metabolite_name) in anns {
            for cpd in cpds_of(key) {
                detected.insert(cpd.clone());
                names
                    .entry(cpd)
                    .or_default()
                    .insert(metabolite_name.clone());
            }
        }
    }
    let names = names
        .into_iter()
        .map(|(cpd, set)| {
            let mut v: Vec<String> = set.into_iter().collect();
            v.sort();
            (cpd, v)
        })
        .collect();
    (detected, names)
}

/// Partition the detected compounds by which ionization mode reached them.
///
/// `None` for a single-mode run: there is nothing to partition, and reporting
/// "100 % POS-only" would be a tautology dressed as a finding.
///
/// **Descriptive only.** `D` is the plain union `D_pos ∪ D_neg`, and this
/// function is deliberately incapable of changing it — it takes the per-mode
/// sets and returns three counts. On the enrichment route the equivalent
/// dual-mode logic applies a conflict rule that EXCLUDES compounds; here there
/// is no differential comparison, so there is no directional verdict that could
/// contradict another and nothing to exclude.
fn partition_by_mode(per_mode_cpds: &[HashSet<String>]) -> Option<CoverageModePartition> {
    let [pos, neg] = per_mode_cpds else {
        return None;
    };
    Some(CoverageModePartition {
        pos_only: pos.difference(neg).count(),
        neg_only: neg.difference(pos).count(),
        in_both: pos.intersection(neg).count(),
    })
}

/// Phase 1 — collect the InChIKeys to resolve against PubChem. Single-mode
/// preserves the prior ordering (no dedup) so behaviour is bit-equal to
/// single-mode regression baselines; dual-mode unions across modes + sorts +
/// dedups so the single PubChem call covers both modes without duplicate work
/// (the cache is keyed per InChIKey).
fn collect_inchikeys(dam_results: &[DamResult], is_dual: bool) -> Vec<String> {
    if is_dual {
        let mut set: HashSet<String> = HashSet::new();
        for dam in dam_results {
            for f in &dam.features {
                if let Some(k) = &f.inchikey {
                    set.insert(k.clone());
                }
            }
        }
        let mut v: Vec<String> = set.into_iter().collect();
        v.sort();
        v
    } else {
        dam_results[0]
            .features
            .iter()
            .filter_map(|f| f.inchikey.clone())
            .collect()
    }
}

/// Build the per-mode `feature → cpds` maps (multi-mapping rule D8) and the
/// measurable-metabolome universe N (union of mapped cpd sets across modes,
/// D1). Single-mode is the degenerate length-1 case: the only map entry equals
/// the pre-Track-7.3 `feature_to_cpds` map and N equals its flattened cpd set
/// — bit-equal to the pre-7.3 behaviour.
fn build_per_mode_maps(
    dam_results: &[DamResult],
    inchikey_to_cids: &HashMap<String, Vec<String>>,
    cid_to_cpd: &HashMap<String, Option<String>>,
) -> (Vec<HashMap<String, HashSet<String>>>, HashSet<String>) {
    let per_mode_feature_to_cpds: Vec<HashMap<String, HashSet<String>>> = dam_results
        .iter()
        .map(|dam| build_feature_to_cpds_for_mode(dam, inchikey_to_cids, cid_to_cpd))
        .collect();
    let mapped_universe: HashSet<String> = per_mode_feature_to_cpds
        .iter()
        .flat_map(|m| m.values().flatten().cloned())
        .collect();
    (per_mode_feature_to_cpds, mapped_universe)
}

/// Emit the K / universe / conflict diagnostics followed by the entry-assembly
/// summary. Called once `entries` is built; the K logs precede the entry logs,
/// matching the pre-decomposition global order (entry assembly is silent).
/// Reproduces the exact `info!` events — same targets, fields, and order.
fn log_k_diagnostics(
    dam_cpd: &HashSet<String>,
    dual_mode_breakdown: Option<&DualModeBreakdown>,
    n_excluded_conflict: usize,
    conflict_sample: &[String],
    entries: &[KeggCompoundSet],
    module_retention: &Option<ModuleRetention>,
) {
    info!(
        foreground_k_size = dam_cpd.len(),
        k_sample = %fmt_cpd_sample(dam_cpd.iter().map(|s| s.as_str()), 20),
        "Foreground K (DAM-significant cpds matching direction) built"
    );
    if let Some(b) = dual_mode_breakdown {
        info!(
            universe_pos_only = b.universe_pos_only,
            universe_neg_only = b.universe_neg_only,
            universe_in_both = b.universe_in_both,
            foreground_pos_only = b.foreground_pos_only,
            foreground_neg_only = b.foreground_neg_only,
            foreground_agree_both = b.foreground_agree_both,
            foreground_excluded_conflict = b.foreground_excluded_conflict,
            "Dual-mode breakdown of N and K"
        );
    }
    // `conflict_sample` is truncated to 100 IDs in `build_dam_cpd_from_trends` for log
    // payload size; the true count lives on `n_excluded_conflict`. Logging
    // `conflict_sample.len()` would silently cap the reported number at 100 even
    // when hundreds of cpds were excluded — surface the count as `n_conflicts`
    // and keep `sample_size` as the (capped) ID-list length. Fires for BOTH
    // single- and dual-mode now that single-mode also applies the conflict rule.
    if n_excluded_conflict > 0 {
        info!(
            n_conflicts = n_excluded_conflict,
            sample_size = conflict_sample.len(),
            sample = %fmt_cpd_sample(conflict_sample.iter().map(|s| s.as_str()), 20),
            "cpds excluded from K by the conflict-only-strict rule"
        );
    }
    if dam_cpd.is_empty() {
        info!(
            "K = 0 — no DAM features pass the active trend filter under \
             these thresholds. Result will have all p_value = 1.0. Consider \
             loosening FC / FDR / δ thresholds or switching direction."
        );
    }

    // Log mode-specific entry-assembly results.
    match module_retention {
        Some(r) => {
            info!(
                total_modules_in_cache = r.total_modules,
                retained_after_group_filter = r.retained_modules,
                group_name = %r.group_name,
                group_level = r.group_level,
                group_org_count = r.group_org_count,
                min_group_overlap = r.min_group_overlap,
                "Module entries assembled (Group filter applied)"
            );
        }
        None => {
            info!(pathway_entries = entries.len(), "Pathway entries assembled");
        }
    }
    let entries_with_compounds = entries.iter().filter(|e| !e.compounds.is_empty()).count();
    info!(
        entries_total = entries.len(),
        entries_with_compounds,
        entries_empty_compound_list = entries.len() - entries_with_compounds,
        "Entry catalogue ready for ORA"
    );
}

/// Emit the post-ORA K-coverage diagnostics: the `info!` summary, the
/// uncovered-K sample, the no-hits mode hint, and the per-entry `debug!` hit
/// detail (capped at 50). `k_covered` / `k_uncovered` are derived in
/// `run_stage3` (k_covered is reused by the funnel). Reproduces the exact
/// events — same targets, fields, order, DEBUG-gating, and `.take(50)`.
fn log_ora_diagnostics(
    enrichment_result: &EnrichmentResult,
    dam_cpd: &HashSet<String>,
    target: &AnalysisTarget,
    k_covered: &HashSet<&String>,
    k_uncovered: &[&String],
) {
    let entries_with_hits = enrichment_result.rows.iter().filter(|r| r.hits > 0).count();
    info!(
        entries_tested = enrichment_result.rows.len(),
        entries_with_hits,
        k_compounds_covered = k_covered.len(),
        k_compounds_uncovered = k_uncovered.len(),
        "Phase 3 ORA complete"
    );
    if !k_uncovered.is_empty() {
        let sample = fmt_cpd_sample(k_uncovered.iter().map(|s| s.as_str()), 20);
        info!(
            uncovered_k_sample = %sample,
            "K cpds NOT covered by any entry's compound list — these cpds \
             contributed to K but didn't intersect with the active catalogue"
        );
    }
    if entries_with_hits == 0 && !dam_cpd.is_empty() {
        let mode_hint = match target.mode() {
            AnalysisMode::Module => {
                "Try pathway mode or a broader Group; module compound \
                 catalogues are narrow and may not overlap with your K"
            }
            AnalysisMode::Pathway => {
                "Try a more central species (broader pathway catalogue) or \
                 loosen DAM thresholds to expand K"
            }
        };
        info!(
            mode_hint = %mode_hint,
            "K is non-empty but NO entry contains any K cpd"
        );
    }
    // Per-entry hit detail at DEBUG (RUST_LOG=debug to surface).
    for row in enrichment_result
        .rows
        .iter()
        .filter(|r| r.hits > 0)
        .take(50)
    {
        debug!(
            entry_id = %row.entry_id,
            entry_name = %row.entry_name,
            hits = row.hits,
            total = row.total,
            p_value = row.p_value,
            fdr = row.fdr,
            hit_cpds = %row.hit_kegg_ids.join(","),
            "Entry with hits"
        );
    }
}

/// Build the Data-tab provenance funnel counts (read-only by-products of the
/// universe/foreground construction; no new network calls). `detected_inchikeys`
/// uses the resolved map's key count rather than the raw `inchikeys` Vec
/// (single-mode `inchikeys` is NOT deduped, so its length would overcount; the
/// map has one key per unique input). `detected_in_entries` is the number of K
/// cpds covered by some entry (`k_covered.len()`, passed by the caller).
fn build_funnel(
    dam_results: &[DamResult],
    per_mode_trends: &[Vec<crate::dam::types::Trend>],
    inchikey_to_cids: &HashMap<String, Vec<String>>,
    all_cids: &[String],
    params: &Stage3Params,
    k_covered_len: usize,
) -> Stage3Funnel {
    let fg_inchikeys =
        foreground_inchikey_set_from_trends(dam_results, per_mode_trends, params.direction);
    let mut fg_cids: HashSet<&str> = HashSet::new();
    for k in &fg_inchikeys {
        if let Some(cids) = inchikey_to_cids.get(k) {
            for c in cids {
                fg_cids.insert(c.as_str());
            }
        }
    }
    Stage3Funnel {
        detected_inchikeys: inchikey_to_cids.len(),
        detected_cids: all_cids.len(),
        foreground_inchikeys: fg_inchikeys.len(),
        foreground_cids: fg_cids.len(),
        detected_in_entries: k_covered_len,
    }
}

/// Distinct InChIKeys among DAM features classified significant in the active
/// direction (across all modes). Pre-conflict — the conflict-only-strict
/// exclusion happens at the cpd/K level and is reported separately via
/// `DualModeBreakdown.foreground_excluded_conflict`. Used for the Data-tab
/// provenance funnel's foreground stage (`add-bottom-panel-data-tab`).
fn foreground_inchikey_set_from_trends(
    dam_results: &[DamResult],
    per_mode_trends: &[Vec<crate::dam::types::Trend>],
    direction: EnrichmentDirection,
) -> HashSet<String> {
    let mut set: HashSet<String> = HashSet::new();
    for (mode_idx, dam) in dam_results.iter().enumerate() {
        for (feat_idx, f) in dam.features.iter().enumerate() {
            let trend = per_mode_trends[mode_idx][feat_idx];
            if matches_direction(trend, direction)
                && let Some(k) = &f.inchikey
            {
                set.insert(k.clone());
            }
        }
    }
    set
}

/// Test-only adapter preserving the pre-D2 `foreground_inchikey_set` call shape
/// (classify from `(method, fc, fdr, delta)`), delegating to the trend-threaded
/// production fn. Production builds `per_mode_trends` once in `run_stage3`.
#[cfg(test)]
fn foreground_inchikey_set(
    dam_results: &[DamResult],
    method: DamMethod,
    fc: f64,
    fdr: f64,
    delta: f64,
    direction: EnrichmentDirection,
) -> HashSet<String> {
    let per_mode_trends: Vec<Vec<crate::dam::types::Trend>> = dam_results
        .iter()
        .map(|dam| {
            dam.features
                .iter()
                .map(|f| classify_trend(f, fc, fdr, delta, method))
                .collect()
        })
        .collect();
    foreground_inchikey_set_from_trends(dam_results, &per_mode_trends, direction)
}

/// Per-mode counterpart of the original inline `feature_to_cpds` build.
/// Pure function — no I/O. Used by both single-mode (length-1 Vec) and
/// dual-mode (length-2 Vec) paths in `run_stage3`.
fn build_feature_to_cpds_for_mode(
    dam: &DamResult,
    inchikey_to_cids: &HashMap<String, Vec<String>>,
    cid_to_cpd: &HashMap<String, Option<String>>,
) -> HashMap<String, HashSet<String>> {
    let mut out: HashMap<String, HashSet<String>> = HashMap::new();
    for feat in &dam.features {
        let Some(inchikey) = &feat.inchikey else {
            continue;
        };
        let cids = match inchikey_to_cids.get(inchikey) {
            Some(v) => v,
            None => continue,
        };
        let mut cpds: HashSet<String> = HashSet::new();
        for cid in cids {
            if let Some(Some(cpd)) = cid_to_cpd.get(cid) {
                cpds.insert(cpd.clone());
            }
        }
        if !cpds.is_empty() {
            out.insert(inchikey.clone(), cpds);
        }
    }
    out
}

/// Merge per-mode `feature_to_cpds` maps into a single map. When the same
/// inchikey appears in multiple modes the cpd-sets are unioned (in practice
/// they are equal — PubChem + KEGG conv are deterministic per inchikey —
/// but the union is the safe definition).
fn merge_feature_to_cpds(
    per_mode: &[HashMap<String, HashSet<String>>],
) -> HashMap<String, HashSet<String>> {
    let mut out: HashMap<String, HashSet<String>> = HashMap::new();
    for m in per_mode {
        for (k, v) in m {
            out.entry(k.clone()).or_default().extend(v.iter().cloned());
        }
    }
    out
}

/// Per-mode aggregated trend (design D5). `Up`/`Down` mean every contributing
/// feature in that mode pointed the same way; `Conflict` means both Up and
/// Down features mapped to the same cpd in a single mode (same-InChIKey-
/// different-trends edge case); `NS` means no significant features; `Absent`
/// means the cpd is unreachable from this mode at all.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ModeTrend {
    Up,
    Down,
    Ns,
    Conflict,
    Absent,
}

fn aggregate_trends(trends: &[crate::dam::types::Trend]) -> ModeTrend {
    use crate::dam::types::Trend;
    if trends.is_empty() {
        return ModeTrend::Absent;
    }
    let has_up = trends.contains(&Trend::Up);
    let has_down = trends.contains(&Trend::Down);
    match (has_up, has_down) {
        (true, true) => ModeTrend::Conflict,
        (true, false) => ModeTrend::Up,
        (false, true) => ModeTrend::Down,
        (false, false) => ModeTrend::Ns,
    }
}

/// Conflict-only-strict membership rule (design D5). For a cpd whose per-mode
/// aggregated trends are `trends`, it enters K iff:
/// - direction = Up: any mode says Up AND no mode says Down AND no mode is Conflict
/// - direction = Down: symmetric
/// - direction = Both: at least one mode says Up or Down AND no Conflict AND
///   not both Up and Down across modes (inter-mode opposing → excluded)
#[cfg(test)]
fn is_in_k_dual(trends: &[ModeTrend], direction: EnrichmentDirection) -> bool {
    let has_up = trends.contains(&ModeTrend::Up);
    let has_down = trends.contains(&ModeTrend::Down);
    let has_conflict = trends.contains(&ModeTrend::Conflict);
    if has_conflict {
        return false;
    }
    match direction {
        EnrichmentDirection::Up => has_up && !has_down,
        EnrichmentDirection::Down => has_down && !has_up,
        EnrichmentDirection::Both => (has_up || has_down) && !(has_up && has_down),
    }
}

/// Whether a cpd was eliminated specifically by the conflict rule (i.e. has
/// signal in the active direction but was kept out of K because some mode
/// pointed the opposite way or had intra-mode Conflict). Used for the
/// `foreground_excluded_conflict` diagnostic count + conflict-sample log.
#[cfg(test)]
fn excluded_by_conflict(trends: &[ModeTrend], direction: EnrichmentDirection) -> bool {
    let has_up = trends.contains(&ModeTrend::Up);
    let has_down = trends.contains(&ModeTrend::Down);
    let has_conflict = trends.contains(&ModeTrend::Conflict);
    // A `Conflict` mode-trend means the SAME InChIKey produced both an Up and
    // a Down DAM feature in one ionization mode — so that mode unambiguously
    // had signal in BOTH directions. `aggregate_trends` collapses (Up, Down)
    // to `Conflict` (a separate variant), which means
    // `trends.contains(&ModeTrend::Up)` is FALSE for a (Conflict, Absent)
    // cpd. Without including `has_conflict` here, those cpds would be
    // dropped from K via `is_in_k_dual` but never counted in
    // `foreground_excluded_conflict` — under-reporting the diagnostic.
    let signal_in_dir = match direction {
        EnrichmentDirection::Up => has_up || has_conflict,
        EnrichmentDirection::Down => has_down || has_conflict,
        EnrichmentDirection::Both => has_up || has_down || has_conflict,
    };
    signal_in_dir && (has_conflict || (has_up && has_down))
}

/// The three-way verdict for a compound's per-mode aggregated trends under the
/// conflict-only-strict union rule. Exhaustive so the K-assembly loop must
/// handle all three; replaces the former `if is_in_k_dual { … } else if
/// excluded_by_conflict { … }` pair.
enum DualMembership {
    /// In K (≥ 1 mode signals the active direction, none opposes, none conflicts).
    InK,
    /// Had signal in the active direction but kept out of K by the conflict
    /// rule (opposing mode or intra-mode Conflict) — counted in
    /// `foreground_excluded_conflict` + the conflict-sample log.
    ExcludedByConflict,
    /// In N but neither in K nor conflict-excluded.
    Neither,
}

/// Single exhaustive classifier replacing `is_in_k_dual` + `excluded_by_conflict`
/// (kept `#[cfg(test)]` as reference oracles). Computes `has_up`/`has_down`/
/// `has_conflict` once and applies the EXACT prior predicate bodies, with
/// `InK` winning over `ExcludedByConflict` (mirroring the former `if/else if`
/// precedence) — so the partition is bit-identical to the prior pair.
fn classify_dual_membership(
    trends: &[ModeTrend],
    direction: EnrichmentDirection,
) -> DualMembership {
    let has_up = trends.contains(&ModeTrend::Up);
    let has_down = trends.contains(&ModeTrend::Down);
    let has_conflict = trends.contains(&ModeTrend::Conflict);

    // InK — the conflict-only-strict membership rule (former `is_in_k_dual`).
    let in_k = !has_conflict
        && match direction {
            EnrichmentDirection::Up => has_up && !has_down,
            EnrichmentDirection::Down => has_down && !has_up,
            EnrichmentDirection::Both => (has_up || has_down) && !(has_up && has_down),
        };
    if in_k {
        return DualMembership::InK;
    }

    // ExcludedByConflict — had signal in the active direction (the
    // `signal_in_dir = has_up || has_conflict` term carries intra-mode
    // `(Conflict, Absent)` cpds, former `excluded_by_conflict`).
    let signal_in_dir = match direction {
        EnrichmentDirection::Up => has_up || has_conflict,
        EnrichmentDirection::Down => has_down || has_conflict,
        EnrichmentDirection::Both => has_up || has_down || has_conflict,
    };
    if signal_in_dir && (has_conflict || (has_up && has_down)) {
        DualMembership::ExcludedByConflict
    } else {
        DualMembership::Neither
    }
}

/// Dual-mode K assembly. Walks each cpd in N, computes per-mode aggregated
/// trends, applies the conflict-only-strict rule, and accumulates
/// `DualModeBreakdown` partition counts. Returns (K, breakdown, conflict_sample).
/// The conflict_sample is sorted + capped at 100 ids for the log message —
/// callers should further trim before display.
fn build_dam_cpd_from_trends(
    dam_results: &[DamResult],
    per_mode_feature_to_cpds: &[HashMap<String, HashSet<String>>],
    per_mode_trends: &[Vec<crate::dam::types::Trend>],
    universe: &HashSet<String>,
    direction: EnrichmentDirection,
) -> (HashSet<String>, DualModeBreakdown, Vec<String>) {
    use crate::dam::types::Trend;
    let n_modes = dam_results.len();

    // DualModeBreakdown is a 2-way POS/NEG partition (indices 0 = POS, 1 = NEG
    // per `IonModeTables::try_new` canonical ordering). The spec contract on
    // the struct doc — `|K| == foreground_pos_only + foreground_neg_only +
    // foreground_agree_both` — holds only when n_modes ∈ {1, 2}. With n_modes
    // ≥ 3 a cpd with signal only in mode 2+ would enter K via `is_in_k_dual`
    // (which sees all modes via `trends.contains(...)`) but contribute nothing
    // to the partition counters (pos/neg below are hardcoded to indices 0/1),
    // silently breaking the invariant. The app guarantees n_modes ≤ 2 via
    // IonModeTables; this assert pins the contract so a future multi-mode
    // extension fails loudly in dev/test instead of leaking inconsistent
    // breakdowns.
    debug_assert!(
        n_modes <= 2,
        "build_dam_cpd_from_trends: DualModeBreakdown assumes n_modes ∈ {{1, 2}}; got {n_modes}",
    );

    // Per mode: cpd → list of per-feature trends that touched this cpd.
    let mut per_mode_cpd_trends: Vec<HashMap<String, Vec<Trend>>> = vec![HashMap::new(); n_modes];
    for (idx, dam) in dam_results.iter().enumerate() {
        let f2c = &per_mode_feature_to_cpds[idx];
        for (feat_idx, feat) in dam.features.iter().enumerate() {
            let Some(inchikey) = &feat.inchikey else {
                continue;
            };
            let Some(cpds) = f2c.get(inchikey) else {
                continue;
            };
            let trend = per_mode_trends[idx][feat_idx];
            for c in cpds {
                per_mode_cpd_trends[idx]
                    .entry(c.clone())
                    .or_default()
                    .push(trend);
            }
        }
    }

    // Walk N, aggregate per mode, apply the rule via the single exhaustive
    // classifier, and accumulate partitions through the breakdown's own tally.
    let mut k: HashSet<String> = HashSet::new();
    let mut breakdown = DualModeBreakdown::default();
    let mut conflicts: Vec<String> = Vec::new();

    // n_modes ≤ 2 (asserted at fn entry). Indices: 0 = POS, 1 = NEG.
    for c in universe {
        let trends: Vec<ModeTrend> = (0..n_modes)
            .map(|i| match per_mode_cpd_trends[i].get(c) {
                Some(ts) => aggregate_trends(ts),
                None => ModeTrend::Absent,
            })
            .collect();
        let pos = trends.first().copied().unwrap_or(ModeTrend::Absent);
        let neg = trends.get(1).copied().unwrap_or(ModeTrend::Absent);

        breakdown.tally_universe(pos, neg);

        match classify_dual_membership(&trends, direction) {
            DualMembership::InK => {
                k.insert(c.clone());
                breakdown.tally_in_k(pos, neg);
            }
            DualMembership::ExcludedByConflict => {
                breakdown.foreground_excluded_conflict += 1;
                conflicts.push(c.clone());
            }
            DualMembership::Neither => {}
        }
    }
    conflicts.sort();
    let conflict_sample: Vec<String> = conflicts.into_iter().take(100).collect();

    (k, breakdown, conflict_sample)
}

/// Test-only adapter preserving the pre-D2 `build_dam_cpd_dual` call shape:
/// classifies each feature from `(method, fc, fdr, delta)` and delegates to
/// `build_dam_cpd_from_trends`. Production threads precomputed `per_mode_trends`
/// instead (see `run_stage3`), so this re-classification lives only in tests.
#[cfg(test)]
#[allow(clippy::too_many_arguments)]
fn build_dam_cpd_dual(
    dam_results: &[DamResult],
    per_mode_feature_to_cpds: &[HashMap<String, HashSet<String>>],
    universe: &HashSet<String>,
    method: DamMethod,
    fc_threshold: f64,
    fdr_threshold: f64,
    delta_threshold: f64,
    direction: EnrichmentDirection,
) -> (HashSet<String>, DualModeBreakdown, Vec<String>) {
    let per_mode_trends: Vec<Vec<crate::dam::types::Trend>> = dam_results
        .iter()
        .map(|dam| {
            dam.features
                .iter()
                .map(|f| classify_trend(f, fc_threshold, fdr_threshold, delta_threshold, method))
                .collect()
        })
        .collect();
    build_dam_cpd_from_trends(
        dam_results,
        per_mode_feature_to_cpds,
        &per_mode_trends,
        universe,
        direction,
    )
}

/// Filter `modules_pack` by `|complete_orgs ∩ group_orgs| >= min_group_overlap`,
/// map retained entries to `KeggCompoundSet`, and compute a `ModuleRetention`
/// summary including the time span across retained modules' `fetched_at`.
fn assemble_module_entries(
    modules_pack: &KeggModulesCache,
    group_level: u8,
    group_name: &str,
    group_orgs: &HashSet<String>,
    min_group_overlap: usize,
) -> (Vec<KeggCompoundSet>, Option<ModuleRetention>) {
    let total_modules = modules_pack.modules.len();
    let mut entries: Vec<KeggCompoundSet> = Vec::new();
    let mut oldest: Option<DateTime<Utc>> = None;
    let mut newest: Option<DateTime<Utc>> = None;

    // Sort by module ID for deterministic output (HashMap iteration is
    // non-deterministic).
    let mut ids: Vec<&String> = modules_pack.modules.keys().collect();
    ids.sort();
    for id in ids {
        let entry = &modules_pack.modules[id];
        let overlap = entry.complete_orgs.intersection(group_orgs).count();
        if overlap < min_group_overlap {
            continue;
        }
        entries.push(KeggCompoundSet {
            id: id.clone(),
            name: entry.name.clone(),
            compounds: entry.compounds.clone(),
        });
        oldest = Some(oldest.map_or(entry.fetched_at, |o| o.min(entry.fetched_at)));
        newest = Some(newest.map_or(entry.fetched_at, |n| n.max(entry.fetched_at)));
    }

    let retention = ModuleRetention {
        total_modules,
        retained_modules: entries.len(),
        min_group_overlap,
        group_level,
        group_name: group_name.to_string(),
        group_org_count: group_orgs.len(),
        // When no modules survive the filter (empty Group, no overlap),
        // both timestamps default to Utc::now() — the UI surfaces an
        // empty-result message in this case, so the timestamps are
        // never user-visible.
        oldest_fetched_at: oldest.unwrap_or_else(Utc::now),
        newest_fetched_at: newest.unwrap_or_else(Utc::now),
    };
    (entries, Some(retention))
}

/// Compute the DAM cpd set K under the active direction filter, WITHOUT the
/// conflict-only-strict exclusion. This is the **pre-fix** single-mode K helper,
/// retained `#[cfg(test)]` as a reference oracle: production single-mode now
/// routes through `build_dam_cpd_dual` (degenerate length-1), which DOES apply
/// the conflict rule. Tests use this to document where the two agree (no-conflict
/// inputs) and diverge (intra-mode-conflict inputs).
#[cfg(test)]
fn build_dam_cpd(
    dam_result: &DamResult,
    feature_to_cpds: &HashMap<String, HashSet<String>>,
    method: DamMethod,
    fc_threshold: f64,
    fdr_threshold: f64,
    delta_threshold: f64,
    direction: EnrichmentDirection,
) -> HashSet<String> {
    let mut k = HashSet::new();
    for feat in &dam_result.features {
        let trend = classify_trend(feat, fc_threshold, fdr_threshold, delta_threshold, method);
        let include = matches_direction(trend, direction);
        if !include {
            continue;
        }
        let Some(inchikey) = &feat.inchikey else {
            continue;
        };
        if let Some(cpds) = feature_to_cpds.get(inchikey) {
            for c in cpds {
                k.insert(c.clone());
            }
        }
    }
    k
}

fn matches_direction(trend: crate::dam::types::Trend, direction: EnrichmentDirection) -> bool {
    use crate::dam::types::Trend;
    matches!(
        (direction, trend),
        (EnrichmentDirection::Up, Trend::Up)
            | (EnrichmentDirection::Down, Trend::Down)
            | (EnrichmentDirection::Both, Trend::Up | Trend::Down)
    )
}

fn compute_pubchem_time_span(
    cache: &HashMap<String, InchikeyCidsEntry>,
    inchikeys: &[String],
) -> Option<(DateTime<Utc>, DateTime<Utc>, usize)> {
    // Single-mode callers pass a per-feature InChIKey vec that may contain
    // duplicates (two features sharing the same annotation); dedup defensively
    // so `count` reports the number of distinct cached entries, not the number
    // of feature rows.
    let mut seen: HashSet<&String> = HashSet::new();
    let mut min_ts: Option<DateTime<Utc>> = None;
    let mut max_ts: Option<DateTime<Utc>> = None;
    let mut count: usize = 0;
    for k in inchikeys {
        if !seen.insert(k) {
            continue;
        }
        if let Some(entry) = cache.get(k) {
            count += 1;
            min_ts = Some(min_ts.map_or(entry.fetched_at, |m| m.min(entry.fetched_at)));
            max_ts = Some(max_ts.map_or(entry.fetched_at, |m| m.max(entry.fetched_at)));
        }
    }
    match (min_ts, max_ts) {
        (Some(a), Some(b)) => Some((a, b, count)),
        _ => None,
    }
}

fn compute_kegg_conv_time_span(
    cache: &HashMap<String, CidCpdEntry>,
    cids: &[String],
) -> Option<(DateTime<Utc>, DateTime<Utc>, usize)> {
    let mut min_ts: Option<DateTime<Utc>> = None;
    let mut max_ts: Option<DateTime<Utc>> = None;
    let mut count: usize = 0;
    for c in cids {
        if let Some(entry) = cache.get(c) {
            count += 1;
            min_ts = Some(min_ts.map_or(entry.fetched_at, |m| m.min(entry.fetched_at)));
            max_ts = Some(max_ts.map_or(entry.fetched_at, |m| m.max(entry.fetched_at)));
        }
    }
    match (min_ts, max_ts) {
        (Some(a), Some(b)) => Some((a, b, count)),
        _ => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dam::types::{DamFeature, FcBasis, Trend};
    use tracing_test::traced_test;

    /// `captured_species` is the whole point of the capture: it reads the run's
    /// own target and has no access to `SessionSettings`, so a completed run's
    /// identity cannot be re-derived from a selection the user may since have
    /// lost. Asserted here rather than through `run_coverage`, which is `async`
    /// and needs two live clients — this is the seam that makes the property
    /// testable at all.
    #[test]
    fn captured_species_reads_the_runs_own_target() {
        let target = AnalysisTarget::Pathway {
            species_kegg: crate::kegg::SpeciesKegg {
                code: "gmx".to_string(),
                fetched_at: Utc::now(),
                pathways: vec![],
            },
        };
        assert_eq!(captured_species(&target), Some("gmx".to_string()));
    }

    /// Module mode captures nothing new — `ModuleRetention` already carries the
    /// Group name and level, so a second copy would be two sources for one fact.
    #[test]
    fn captured_species_is_none_in_module_mode() {
        let target = AnalysisTarget::Module {
            modules_pack: crate::kegg::KeggModulesCache {
                modules: std::collections::HashMap::new(),
            },
            group_level: 2,
            group_name: "Plants".to_string(),
            group_org_codes: std::collections::HashSet::new(),
            min_group_overlap: 1,
        };
        assert_eq!(captured_species(&target), None);
    }

    // Decision A (logs are contract): the conflict-only-strict exclusion INFO log
    // is verified directly. `log_k_diagnostics` is private (only reachable from
    // inside this module — `run_stage3` would need a fully resolved &[DamResult] +
    // network-derived maps), so the test calls it directly with
    // `n_excluded_conflict > 0` and a `conflict_sample` naming a known cpd, and
    // asserts the event fires and its `sample` field names that cpd.
    #[traced_test]
    #[test]
    fn conflict_sample_info_log_fires_and_names_excluded_cpd() {
        let dam_cpd: std::collections::HashSet<String> =
            ["C00001".to_string()].into_iter().collect();
        let conflict_sample = vec!["C00002".to_string()];
        let entries: Vec<crate::kegg::KeggCompoundSet> = Vec::new();
        log_k_diagnostics(&dam_cpd, None, 1, &conflict_sample, &entries, &None);
        assert!(
            logs_contain("cpds excluded from K by the conflict-only-strict rule"),
            "conflict-exclusion INFO event must fire when n_excluded_conflict > 0"
        );
        assert!(
            logs_contain("C00002"),
            "the conflict event's `sample` field must name the excluded cpd"
        );
    }

    fn synth_feature(
        inchikey: Option<&str>,
        p_adjusted: f64,
        log2_fc: f64,
        effect_size: Option<f64>,
    ) -> DamFeature {
        DamFeature {
            alignment_id: "x".into(),
            metabolite_name: "x".into(),
            inchikey: inchikey.map(str::to_string),
            average_rt_min: None,
            average_mz: None,
            formula: None,
            smiles: None,
            numerator_mean: 0.0,
            denominator_mean: 0.0,
            numerator_median: 0.0,
            denominator_median: 0.0,
            fold_change: 0.0,
            log2_fold_change: log2_fc,
            fc_basis: FcBasis::Mean,
            p_value: 0.0,
            p_adjusted,
            neg_log10_p_adjusted: 0.0,
            effect_size,
        }
    }

    #[test]
    fn matches_direction_up_only() {
        assert!(matches_direction(Trend::Up, EnrichmentDirection::Up));
        assert!(!matches_direction(Trend::Down, EnrichmentDirection::Up));
        assert!(!matches_direction(
            Trend::NotSignificant,
            EnrichmentDirection::Up
        ));
    }

    #[test]
    fn matches_direction_down_only() {
        assert!(!matches_direction(Trend::Up, EnrichmentDirection::Down));
        assert!(matches_direction(Trend::Down, EnrichmentDirection::Down));
    }

    #[test]
    fn matches_direction_both() {
        assert!(matches_direction(Trend::Up, EnrichmentDirection::Both));
        assert!(matches_direction(Trend::Down, EnrichmentDirection::Both));
        assert!(!matches_direction(
            Trend::NotSignificant,
            EnrichmentDirection::Both
        ));
    }

    fn dam_result_with(features: Vec<DamFeature>) -> DamResult {
        DamResult {
            method: DamMethod::Welch,
            numerator: "T".into(),
            denominator: "C".into(),
            features,
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiHochberg,
            dedup_report: None,
        }
    }

    #[test]
    fn foreground_inchikey_set_respects_direction_and_dedups() {
        // KA Up, KB Down, KC ns, None=Up-but-no-key, duplicate KA Up.
        let dam = dam_result_with(vec![
            synth_feature(Some("KA"), 0.01, 2.0, None),
            synth_feature(Some("KB"), 0.01, -2.0, None),
            synth_feature(Some("KC"), 0.50, 2.0, None),
            synth_feature(None, 0.01, 2.0, None),
            synth_feature(Some("KA"), 0.01, 2.0, None), // dup InChIKey
        ]);
        let results = [dam];

        let both = foreground_inchikey_set(
            &results,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Both,
        );
        assert_eq!(both, ["KA", "KB"].iter().map(|s| s.to_string()).collect());

        let up = foreground_inchikey_set(
            &results,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert_eq!(up, ["KA"].iter().map(|s| s.to_string()).collect());

        let down = foreground_inchikey_set(
            &results,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Down,
        );
        assert_eq!(down, ["KB"].iter().map(|s| s.to_string()).collect());
    }

    #[test]
    fn foreground_inchikey_set_unions_across_modes() {
        let pos = dam_result_with(vec![synth_feature(Some("KA"), 0.01, 2.0, None)]); // Up
        let neg = dam_result_with(vec![
            synth_feature(Some("KA"), 0.01, 2.0, None), // Up (same key as POS)
            synth_feature(Some("KZ"), 0.01, 2.0, None), // Up (NEG-only)
        ]);
        let set = foreground_inchikey_set(
            &[pos, neg],
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Both,
        );
        // Union, deduped: KA (both modes) + KZ (NEG only).
        assert_eq!(set, ["KA", "KZ"].iter().map(|s| s.to_string()).collect());
    }

    #[test]
    fn build_dam_cpd_respects_direction_up() {
        let result = DamResult {
            method: DamMethod::Welch,
            numerator: "T".into(),
            denominator: "C".into(),
            features: vec![
                synth_feature(Some("KA"), 0.01, 2.0, None),  // Up
                synth_feature(Some("KB"), 0.01, -2.0, None), // Down
                synth_feature(Some("KC"), 0.50, 2.0, None),  // ns
                synth_feature(None, 0.01, 2.0, None),        // Up but no InChIKey
            ],
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            dedup_report: None,
        };
        let mut feature_to_cpds: HashMap<String, HashSet<String>> = HashMap::new();
        feature_to_cpds.insert(
            "KA".to_string(),
            ["C00001".to_string()].into_iter().collect(),
        );
        feature_to_cpds.insert(
            "KB".to_string(),
            ["C00002".to_string()].into_iter().collect(),
        );
        feature_to_cpds.insert(
            "KC".to_string(),
            ["C00003".to_string()].into_iter().collect(),
        );

        let k = build_dam_cpd(
            &result,
            &feature_to_cpds,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert_eq!(k, ["C00001".to_string()].into_iter().collect());

        let k_both = build_dam_cpd(
            &result,
            &feature_to_cpds,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Both,
        );
        assert_eq!(
            k_both,
            ["C00001".to_string(), "C00002".to_string()]
                .into_iter()
                .collect()
        );
    }

    #[test]
    fn compute_time_span_returns_none_when_no_overlap() {
        let cache: HashMap<String, InchikeyCidsEntry> = HashMap::new();
        let inchikeys = vec!["A".to_string()];
        assert!(compute_pubchem_time_span(&cache, &inchikeys).is_none());
    }

    #[test]
    fn compute_time_span_finds_min_max() {
        use chrono::TimeZone;
        let mut cache: HashMap<String, InchikeyCidsEntry> = HashMap::new();
        let t1 = Utc.with_ymd_and_hms(2026, 3, 1, 12, 0, 0).unwrap();
        let t2 = Utc.with_ymd_and_hms(2026, 5, 22, 12, 0, 0).unwrap();
        cache.insert(
            "A".to_string(),
            InchikeyCidsEntry {
                cids: vec![],
                fetched_at: t1,
            },
        );
        cache.insert(
            "B".to_string(),
            InchikeyCidsEntry {
                cids: vec![],
                fetched_at: t2,
            },
        );
        let inchikeys = vec!["A".to_string(), "B".to_string()];
        let span = compute_pubchem_time_span(&cache, &inchikeys).unwrap();
        assert_eq!(span.0, t1);
        assert_eq!(span.1, t2);
        assert_eq!(span.2, 2);
    }

    /// Single-mode callers pass a per-feature InChIKey vec that may carry
    /// duplicates (two DAM features sharing one annotation). The count returned
    /// must be the number of DISTINCT cached entries, not the number of input
    /// rows — otherwise the Stage 3 result panel over-reports "N entries" for
    /// the PubChem cache time span.
    #[test]
    fn compute_time_span_dedupes_duplicate_inputs() {
        use chrono::TimeZone;
        let mut cache: HashMap<String, InchikeyCidsEntry> = HashMap::new();
        let t1 = Utc.with_ymd_and_hms(2026, 3, 1, 12, 0, 0).unwrap();
        cache.insert(
            "A".to_string(),
            InchikeyCidsEntry {
                cids: vec![],
                fetched_at: t1,
            },
        );
        let inchikeys = vec!["A".to_string(), "A".to_string(), "A".to_string()];
        let span = compute_pubchem_time_span(&cache, &inchikeys).unwrap();
        assert_eq!(span.2, 1, "duplicate inputs must not inflate count");
    }

    // ── Track 7.8: dual-mode N/K math ──

    /// Build a `(DamResult, f2c)` pair from a list of
    /// `(inchikey, log2fc, p_adj, &[cpd...])`. Method is fixed at Welch.
    fn mk_mode(
        rows: &[(&str, f64, f64, &[&str])],
    ) -> (DamResult, HashMap<String, HashSet<String>>) {
        let mut features = Vec::new();
        let mut f2c: HashMap<String, HashSet<String>> = HashMap::new();
        for (ik, log2fc, p_adj, cpds) in rows {
            features.push(synth_feature(Some(ik), *p_adj, *log2fc, None));
            f2c.insert(
                (*ik).to_string(),
                cpds.iter().map(|c| (*c).to_string()).collect(),
            );
        }
        let dam = DamResult {
            method: DamMethod::Welch,
            numerator: "T".into(),
            denominator: "C".into(),
            features,
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiHochberg,
            dedup_report: None,
        };
        (dam, f2c)
    }

    fn n_union(maps: &[&HashMap<String, HashSet<String>>]) -> HashSet<String> {
        maps.iter()
            .flat_map(|m| m.values().flatten().cloned())
            .collect()
    }

    /// Single-mode surfaces no DualModeBreakdown (`dual_mode_breakdown == None`).
    /// Production single-mode now ROUTES THROUGH `build_dam_cpd_dual` (degenerate
    /// length-1) and computes a breakdown internally, but the orchestrator wraps
    /// it in `Some` only when `is_dual` — so single-mode still surfaces `None`.
    /// This test uses the `#[cfg(test)]` `build_dam_cpd` reference oracle to check
    /// the no-conflict K is unchanged.
    #[test]
    fn single_mode_dam_result_produces_no_breakdown() {
        let (dam, f2c) = mk_mode(&[("KA", 2.0, 0.01, &["C00001"])]);
        let n = n_union(&[&f2c]);
        // `build_dam_cpd` is the single-mode path. There's no breakdown to
        // return — only the K set. We verify K is sensible and rely on the
        // dispatch in `run_stage3` (which always uses None for len==1).
        let k = build_dam_cpd(
            &dam,
            &f2c,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert_eq!(k, ["C00001".to_string()].into_iter().collect());
        // Also: aggregator on a length-1 universe of length-1 features
        // produces Up, not Conflict.
        let trends: Vec<crate::dam::types::Trend> = dam
            .features
            .iter()
            .map(|f| classify_trend(f, 2.0, 0.05, 0.33, DamMethod::Welch))
            .collect();
        assert_eq!(aggregate_trends(&trends), ModeTrend::Up);
        let _ = n;
    }

    /// N in dual-mode covers POS-only, NEG-only, and both-mode cpds.
    #[test]
    fn dual_mode_union_n_covers_pos_only_neg_only_and_both() {
        let (pos, f2c_pos) = mk_mode(&[
            ("KA", 2.0, 0.01, &["C_BOTH"]),
            ("KB", 2.0, 0.01, &["C_POS"]),
        ]);
        let (neg, f2c_neg) = mk_mode(&[
            ("KA", 2.0, 0.01, &["C_BOTH"]),
            ("KC", 2.0, 0.01, &["C_NEG"]),
        ]);
        let universe = n_union(&[&f2c_pos, &f2c_neg]);
        assert_eq!(universe.len(), 3);
        let (_, breakdown, _) = build_dam_cpd_dual(
            &[pos, neg],
            &[f2c_pos, f2c_neg],
            &universe,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert_eq!(breakdown.universe_pos_only, 1);
        assert_eq!(breakdown.universe_neg_only, 1);
        assert_eq!(breakdown.universe_in_both, 1);
    }

    /// K excludes hard inter-mode direction conflict: same cpd Up in POS,
    /// Down in NEG → excluded.
    #[test]
    fn dual_mode_k_strict_excludes_hard_direction_conflict() {
        let (pos, f2c_pos) = mk_mode(&[("KA", 2.0, 0.01, &["C00001"])]);
        let (neg, f2c_neg) = mk_mode(&[("KA", -2.0, 0.01, &["C00001"])]);
        let universe = n_union(&[&f2c_pos, &f2c_neg]);
        let (k, breakdown, conflicts) = build_dam_cpd_dual(
            &[pos, neg],
            &[f2c_pos, f2c_neg],
            &universe,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert!(
            k.is_empty(),
            "C00001 should be excluded by conflict; K={k:?}"
        );
        assert_eq!(breakdown.foreground_excluded_conflict, 1);
        assert!(conflicts.contains(&"C00001".to_string()));
    }

    /// K keeps a cpd that is significant in one mode and NS / Absent in the
    /// other (no opposing-direction signal anywhere).
    #[test]
    fn dual_mode_k_keeps_one_mode_significant_other_ns_or_absent() {
        // POS: KA Up → C00001. NEG: KX NS (large p_adj) → also maps to C00001.
        let (pos, f2c_pos) = mk_mode(&[("KA", 2.0, 0.01, &["C00001"])]);
        let (neg, f2c_neg) = mk_mode(&[("KX", 2.0, 0.50, &["C00001"])]);
        let universe = n_union(&[&f2c_pos, &f2c_neg]);
        let (k, _, _) = build_dam_cpd_dual(
            &[pos.clone(), neg.clone()],
            &[f2c_pos.clone(), f2c_neg.clone()],
            &universe,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert!(
            k.contains("C00001"),
            "K should keep one-mode-signal cpd; K={k:?}"
        );

        // NEG absent for C00002: only POS has any feature mapping there.
        let (pos2, f2c_pos2) = mk_mode(&[("KA", 2.0, 0.01, &["C00002"])]);
        let (neg2, f2c_neg2) = mk_mode(&[("KZ", 0.1, 0.50, &[])]);
        let universe2 = n_union(&[&f2c_pos2, &f2c_neg2]);
        let (k2, _, _) = build_dam_cpd_dual(
            &[pos2, neg2],
            &[f2c_pos2, f2c_neg2],
            &universe2,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert!(k2.contains("C00002"));
    }

    /// K excludes a cpd flagged as Conflict within a single mode (same
    /// InChIKey having Up + Down features touching it — the same-InChIKey-
    /// different-trends edge case in D5).
    #[test]
    fn dual_mode_k_excludes_intra_mode_conflict() {
        // Same InChIKey duplicated in POS with opposite log2fcs. Both pass
        // significance; aggregate_trends → Conflict.
        let (pos, f2c_pos) = mk_mode(&[
            ("KA", 2.0, 0.01, &["C00001"]),
            ("KA", -2.0, 0.01, &["C00001"]),
        ]);
        let (neg, f2c_neg) = mk_mode(&[("KB", 2.0, 0.01, &["C00001"])]);
        let universe = n_union(&[&f2c_pos, &f2c_neg]);
        let (k, breakdown, _) = build_dam_cpd_dual(
            &[pos, neg],
            &[f2c_pos, f2c_neg],
            &universe,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert!(
            !k.contains("C00001"),
            "intra-mode Conflict → excluded; K={k:?}"
        );
        assert_eq!(breakdown.foreground_excluded_conflict, 1);
    }

    /// K excludes a cpd with intra-POS Conflict and NEG Absent. `has_up` is
    /// false (Conflict is its own variant), so pre-2026-05-26 the diagnostic
    /// counter incorrectly returned 0. The fix makes `signal_in_dir` include
    /// `has_conflict` (Conflict ⇒ Up signal was present in that mode).
    #[test]
    fn dual_mode_intra_pos_conflict_with_neg_absent_counts_as_excluded() {
        // POS: same InChIKey produces Up + Down features → ModeTrend::Conflict.
        // NEG: no features touch C00099 at all → ModeTrend::Absent.
        let (pos, f2c_pos) = mk_mode(&[
            ("KA", 2.0, 0.01, &["C00099"]),
            ("KA", -2.0, 0.01, &["C00099"]),
        ]);
        let (neg, f2c_neg) = mk_mode(&[("KZ", 2.0, 0.5, &[])]); // NEG has no feature for C00099
        let universe = n_union(&[&f2c_pos, &f2c_neg]);
        assert!(universe.contains("C00099"), "C00099 must be in universe");
        let (k, breakdown, conflicts) = build_dam_cpd_dual(
            &[pos, neg],
            &[f2c_pos, f2c_neg],
            &universe,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert!(
            !k.contains("C00099"),
            "intra-POS Conflict blocks K membership; K={k:?}"
        );
        assert_eq!(
            breakdown.foreground_excluded_conflict, 1,
            "POS=Conflict + NEG=Absent must increment the diagnostic counter"
        );
        assert!(
            conflicts.contains(&"C00099".to_string()),
            "C00099 must appear in the conflict_sample log list"
        );
    }

    /// Same Conflict-only-no-NEG-signal case but for direction=Down. Verifies
    /// the symmetric branch of `signal_in_dir`.
    #[test]
    fn dual_mode_intra_pos_conflict_with_neg_absent_counts_as_excluded_for_down() {
        let (pos, f2c_pos) = mk_mode(&[
            ("KA", 2.0, 0.01, &["C00099"]),
            ("KA", -2.0, 0.01, &["C00099"]),
        ]);
        let (neg, f2c_neg) = mk_mode(&[("KZ", 2.0, 0.5, &[])]);
        let universe = n_union(&[&f2c_pos, &f2c_neg]);
        let (_, breakdown, _) = build_dam_cpd_dual(
            &[pos, neg],
            &[f2c_pos, f2c_neg],
            &universe,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Down,
        );
        assert_eq!(breakdown.foreground_excluded_conflict, 1);
    }

    /// Inter-mode conflict via different InChIKeys mapping to the same cpd:
    /// POS has KA→C00001 Up, NEG has KB→C00001 Down. Aggregator says
    /// pos=Up, neg=Down → excluded.
    #[test]
    fn dual_mode_k_excludes_inter_mode_conflict_via_multiple_inchikeys() {
        let (pos, f2c_pos) = mk_mode(&[("KA", 2.0, 0.01, &["C00001"])]);
        let (neg, f2c_neg) = mk_mode(&[("KB", -2.0, 0.01, &["C00001"])]);
        let universe = n_union(&[&f2c_pos, &f2c_neg]);
        let (k, breakdown, _) = build_dam_cpd_dual(
            &[pos, neg],
            &[f2c_pos, f2c_neg],
            &universe,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Both,
        );
        assert!(!k.contains("C00001"));
        assert_eq!(breakdown.foreground_excluded_conflict, 1);
    }

    /// Universe + foreground partition counts arithmetically match the cpd
    /// totals across a deliberately mixed fixture.
    #[test]
    fn dual_mode_breakdown_counts_match_partitions() {
        // Build modes that produce ALL partition classes:
        // - C_AGREE  : POS Up, NEG Up    → foreground_agree_both
        // - C_POSONL : POS Up, NEG NS    → foreground_pos_only
        // - C_NEGONL : POS Absent, NEG Up → foreground_neg_only
        // - C_CONFL  : POS Up, NEG Down  → foreground_excluded_conflict
        // - C_NS_BTH : POS NS, NEG NS    → in universe_in_both, not in K
        let (pos, f2c_pos) = mk_mode(&[
            ("KA", 2.0, 0.01, &["C_AGREE"]),
            ("KB", 2.0, 0.01, &["C_POSONL"]),
            ("KC", 2.0, 0.01, &["C_CONFL"]),
            ("KD", 2.0, 0.50, &["C_NS_BTH"]),
        ]);
        let (neg, f2c_neg) = mk_mode(&[
            ("KA2", 2.0, 0.01, &["C_AGREE"]),
            ("KB2", 2.0, 0.50, &["C_POSONL"]),
            ("KE", 2.0, 0.01, &["C_NEGONL"]),
            ("KC2", -2.0, 0.01, &["C_CONFL"]),
            ("KF", 0.5, 0.50, &["C_NS_BTH"]),
        ]);
        let universe = n_union(&[&f2c_pos, &f2c_neg]);
        let (_k, breakdown, _) = build_dam_cpd_dual(
            &[pos, neg],
            &[f2c_pos, f2c_neg],
            &universe,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        // Universe partition: 0 POS-only (every POS cpd also touched in NEG),
        // 1 NEG-only (C_NEGONL), 4 in-both (C_AGREE/POSONL/CONFL/NS_BTH).
        assert_eq!(breakdown.universe_pos_only, 0);
        assert_eq!(breakdown.universe_neg_only, 1);
        assert_eq!(breakdown.universe_in_both, 4);
        let universe_total =
            breakdown.universe_pos_only + breakdown.universe_neg_only + breakdown.universe_in_both;
        assert_eq!(universe_total, universe.len());
        // Foreground partition: 1 agree_both, 1 pos_only, 1 neg_only, 1 excluded.
        assert_eq!(breakdown.foreground_agree_both, 1);
        assert_eq!(breakdown.foreground_pos_only, 1);
        assert_eq!(breakdown.foreground_neg_only, 1);
        assert_eq!(breakdown.foreground_excluded_conflict, 1);
    }

    /// In dual-mode with the NEG mode empty (no features mapping to cpds),
    /// the dual-mode K equals what the single-mode path would have produced
    /// for POS alone. This is the "regression bridge" test.
    #[test]
    fn dual_mode_ora_matches_single_mode_when_neg_is_empty() {
        let (pos, f2c_pos) = mk_mode(&[
            ("KA", 2.0, 0.01, &["C00001"]),
            ("KB", -2.0, 0.01, &["C00002"]),
            ("KC", 2.0, 0.50, &["C00003"]),
        ]);
        // NEG has no features at all → its f2c is empty.
        let neg = DamResult {
            method: DamMethod::Welch,
            numerator: "T".into(),
            denominator: "C".into(),
            features: vec![],
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiHochberg,
            dedup_report: None,
        };
        let f2c_neg: HashMap<String, HashSet<String>> = HashMap::new();
        let universe = n_union(&[&f2c_pos, &f2c_neg]);

        let single_k = build_dam_cpd(
            &pos,
            &f2c_pos,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Both,
        );
        let (dual_k, _, _) = build_dam_cpd_dual(
            &[pos, neg],
            &[f2c_pos, f2c_neg],
            &universe,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Both,
        );
        assert_eq!(single_k, dual_k);
    }

    /// Behavior fix (refine-stage3-dual-mode-internals): single-mode K now
    /// EXCLUDES a compound reached by both an Up and a Down feature (intra-mode
    /// conflict), because production single-mode routes through the degenerate
    /// length-1 `build_dam_cpd_dual` (which applies the conflict-only-strict
    /// rule). The pre-fix `build_dam_cpd` oracle KEPT such a compound — this
    /// test pins the divergence as the intended consistency fix.
    #[test]
    fn single_mode_excludes_intra_mode_conflict_cpd() {
        // KUP (Up) maps to C_CONFLICT + C_CLEAN; KDOWN (Down) maps to
        // C_CONFLICT. So C_CONFLICT is reached by opposite-trend features.
        let (dam, f2c) = mk_mode(&[
            ("KUP", 2.0, 0.01, &["C_CONFLICT", "C_CLEAN"]),
            ("KDOWN", -2.0, 0.01, &["C_CONFLICT"]),
        ]);
        let universe = n_union(&[&f2c]);

        // Production single-mode path = degenerate length-1 dual.
        let (dual_k, breakdown, conflicts) = build_dam_cpd_dual(
            std::slice::from_ref(&dam),
            std::slice::from_ref(&f2c),
            &universe,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert!(
            !dual_k.contains("C_CONFLICT"),
            "conflicting cpd must be EXCLUDED from single-mode K (the fix)"
        );
        assert!(dual_k.contains("C_CLEAN"), "clean Up cpd stays in K");
        assert_eq!(breakdown.foreground_excluded_conflict, 1);
        assert_eq!(conflicts, vec!["C_CONFLICT".to_string()]);

        // Pre-fix oracle KEPT the conflicting cpd — documents the divergence.
        let pre_fix_k = build_dam_cpd(
            &dam,
            &f2c,
            DamMethod::Welch,
            2.0,
            0.05,
            0.33,
            EnrichmentDirection::Up,
        );
        assert!(
            pre_fix_k.contains("C_CONFLICT"),
            "pre-fix single-mode KEPT the conflicting cpd (intended behavior change)"
        );
        assert!(pre_fix_k.contains("C_CLEAN"));
    }

    /// `classify_dual_membership` reproduces the prior `(is_in_k_dual,
    /// excluded_by_conflict)` verdict (InK-wins precedence) for every ModeTrend
    /// combination of length 1 and 2 × every direction.
    #[test]
    fn classify_dual_membership_matches_prior_predicates() {
        let variants = [
            ModeTrend::Up,
            ModeTrend::Down,
            ModeTrend::Ns,
            ModeTrend::Conflict,
            ModeTrend::Absent,
        ];
        let mut cases: Vec<Vec<ModeTrend>> = Vec::new();
        for &a in &variants {
            cases.push(vec![a]);
            for &b in &variants {
                cases.push(vec![a, b]);
            }
        }
        for dir in [
            EnrichmentDirection::Up,
            EnrichmentDirection::Down,
            EnrichmentDirection::Both,
        ] {
            for trends in &cases {
                let expected = if is_in_k_dual(trends, dir) {
                    "InK"
                } else if excluded_by_conflict(trends, dir) {
                    "ExcludedByConflict"
                } else {
                    "Neither"
                };
                let got = match classify_dual_membership(trends, dir) {
                    DualMembership::InK => "InK",
                    DualMembership::ExcludedByConflict => "ExcludedByConflict",
                    DualMembership::Neither => "Neither",
                };
                assert_eq!(
                    got, expected,
                    "trends={trends:?} dir={dir:?}: classifier diverged from oracle"
                );
            }
        }
        // The intra-mode (Conflict, Absent) under Up is ExcludedByConflict.
        assert!(matches!(
            classify_dual_membership(
                &[ModeTrend::Conflict, ModeTrend::Absent],
                EnrichmentDirection::Up
            ),
            DualMembership::ExcludedByConflict
        ));
    }

    /// `DualModeBreakdown::tally_universe` / `tally_in_k` partition arithmetic.
    #[test]
    fn dual_mode_breakdown_tally_arithmetic() {
        let mut u = DualModeBreakdown::default();
        u.tally_universe(ModeTrend::Up, ModeTrend::Absent); // pos only
        u.tally_universe(ModeTrend::Absent, ModeTrend::Down); // neg only
        u.tally_universe(ModeTrend::Up, ModeTrend::Ns); // both present
        u.tally_universe(ModeTrend::Absent, ModeTrend::Absent); // neither
        assert_eq!(u.universe_pos_only, 1);
        assert_eq!(u.universe_neg_only, 1);
        assert_eq!(u.universe_in_both, 1);

        let mut g = DualModeBreakdown::default();
        g.tally_in_k(ModeTrend::Up, ModeTrend::Ns); // pos sig only (Ns not sig)
        g.tally_in_k(ModeTrend::Ns, ModeTrend::Down); // neg sig only
        g.tally_in_k(ModeTrend::Up, ModeTrend::Down); // both sig
        g.tally_in_k(ModeTrend::Ns, ModeTrend::Absent); // neither sig
        assert_eq!(g.foreground_pos_only, 1);
        assert_eq!(g.foreground_neg_only, 1);
        assert_eq!(g.foreground_agree_both, 1);
    }

    // ── Coverage route: dual-mode union + partition (T6-D) ──

    fn cpd_set(ids: &[&str]) -> HashSet<String> {
        ids.iter().map(|s| s.to_string()).collect()
    }

    /// The `(inchikey, metabolite_name)` annotation list one ion-mode table
    /// contributes — exactly what `PreparedFeatures` carries.
    fn anns(pairs: &[(&str, &str)]) -> Vec<(String, String)> {
        pairs
            .iter()
            .map(|(k, n)| ((*k).to_string(), (*n).to_string()))
            .collect()
    }

    /// `D` is the plain set union across modes — no conflict rule, no
    /// strictness parameter, no direction-based exclusion. With no differential
    /// comparison there is no directional verdict that could contradict
    /// another, so there is nothing a conflict rule could act on.
    #[test]
    fn coverage_dual_mode_unions_without_a_conflict_rule() {
        let pos = cpd_set(&["C00001", "C00002"]);
        let neg = cpd_set(&["C00002", "C00007"]);
        let union: HashSet<String> = pos.union(&neg).cloned().collect();
        assert_eq!(union, cpd_set(&["C00001", "C00002", "C00007"]));

        // Nothing is excluded for appearing in only one mode or in both — the
        // union's size is exactly pos_only + neg_only + in_both.
        let p = partition_by_mode(&[pos, neg]).expect("dual mode partitions");
        assert_eq!(p.pos_only + p.neg_only + p.in_both, union.len());
    }

    /// The three partition buckets, on a set with one compound of each kind.
    #[test]
    fn coverage_mode_partition_counts_each_bucket() {
        let p = partition_by_mode(&[
            cpd_set(&["C00001", "C00002"]),
            cpd_set(&["C00002", "C00007"]),
        ])
        .expect("dual mode partitions");
        assert_eq!(p.pos_only, 1, "C00001");
        assert_eq!(p.neg_only, 1, "C00007");
        assert_eq!(p.in_both, 1, "C00002");
    }

    /// A single-mode run has nothing to partition. Reporting "100 % POS-only"
    /// would be a tautology dressed as a finding.
    #[test]
    fn coverage_mode_partition_is_none_in_single_mode() {
        assert_eq!(partition_by_mode(&[cpd_set(&["C00001"])]), None);
        assert_eq!(partition_by_mode(&[]), None);
    }

    /// Disjoint modes, and one mode empty — the degenerate ends of the range.
    #[test]
    fn coverage_mode_partition_handles_disjoint_and_empty_modes() {
        let disjoint = partition_by_mode(&[cpd_set(&["C00001"]), cpd_set(&["C00007"])])
            .expect("dual mode partitions");
        assert_eq!(
            (disjoint.pos_only, disjoint.neg_only, disjoint.in_both),
            (1, 1, 0)
        );

        let neg_empty = partition_by_mode(&[cpd_set(&["C00001", "C00002"]), cpd_set(&[])])
            .expect("dual mode partitions");
        assert_eq!(
            (neg_empty.pos_only, neg_empty.neg_only, neg_empty.in_both),
            (2, 0, 0)
        );
    }

    /// The one observable effect of deduplication on this route: which MS-DIAL
    /// metabolite name the CSV attaches to a compound.
    ///
    /// Two features share an InChIKey and therefore a cpd, but carry different
    /// names. With both kept, the CSV lists both; with the cascade having
    /// elected one, it lists one. `D` is identical in both cases — which is why
    /// this test, and the funnel, are the entire justification for keeping the
    /// dedup control on the coverage route at all (design D16).
    #[test]
    fn dedup_changes_the_exported_names_but_never_d() {
        let cpds_of = |key: &str| match key {
            "GLUCOSEKEY" => vec!["C00031".to_string()],
            "CITRATEKEY" => vec!["C00158".to_string()],
            _ => vec![],
        };

        // Dedup off: both same-InChIKey features survive, so both names travel.
        let off = anns(&[
            ("GLUCOSEKEY", "D-Glucose"),
            ("GLUCOSEKEY", "Glucose"),
            ("CITRATEKEY", "Citrate"),
        ]);
        // Dedup on: the cascade elected one, so only its name does.
        let on = anns(&[("GLUCOSEKEY", "D-Glucose"), ("CITRATEKEY", "Citrate")]);

        let (d_off, names_off) = build_detected_and_names(&[off], cpds_of);
        let (d_on, names_on) = build_detected_and_names(&[on], cpds_of);

        assert_eq!(d_off, d_on, "D is invariant under deduplication");
        assert_eq!(
            names_off.get("C00031"),
            Some(&vec!["D-Glucose".to_string(), "Glucose".to_string()]),
            "names sorted and both listed"
        );
        assert_eq!(
            names_on.get("C00031"),
            Some(&vec!["D-Glucose".to_string()]),
            "only the elected representative's name survives"
        );
    }

    /// One InChIKey named differently in POS and NEG contributes BOTH names —
    /// which is why the map is built per table rather than over the unioned key
    /// list.
    #[test]
    fn a_compound_named_differently_per_mode_lists_both_names() {
        let cpds_of = |key: &str| {
            if key == "GLUCOSEKEY" {
                vec!["C00031".to_string()]
            } else {
                vec![]
            }
        };
        let (d, names) = build_detected_and_names(
            &[
                anns(&[("GLUCOSEKEY", "Glucose (POS)")]),
                anns(&[("GLUCOSEKEY", "Glucose (NEG)")]),
            ],
            cpds_of,
        );
        assert_eq!(d.len(), 1);
        assert_eq!(
            names.get("C00031"),
            Some(&vec![
                "Glucose (NEG)".to_string(),
                "Glucose (POS)".to_string()
            ])
        );
    }

    /// An annotated feature whose InChIKey resolves to no cpd contributes
    /// nothing. (A feature with no InChIKey at all never reaches the annotation
    /// list in the first place — that exclusion is structural, one layer up.)
    #[test]
    fn unresolvable_features_contribute_nothing() {
        let cpds_of = |key: &str| {
            if key == "CITRATEKEY" {
                vec!["C00158".to_string()]
            } else {
                vec![]
            }
        };
        let (d, names) = build_detected_and_names(
            &[anns(&[("NOCPDKEY", "Unmapped"), ("CITRATEKEY", "Citrate")])],
            cpds_of,
        );
        assert_eq!(d, cpd_set(&["C00158"]));
        assert_eq!(names.len(), 1);
    }

    /// `CoverageParams` carries no statistical field. Asserted by construction,
    /// like the `CoverageResult` guarantee: a `direction` or `fdr_method` added
    /// here would fail to compile against this literal rather than quietly
    /// giving the route a knob it must not have.
    #[test]
    fn coverage_params_carry_no_statistical_field() {
        let CoverageParams {
            selected_groups: _,
            presence_threshold: _,
            dedup_enabled: _,
            dedup_rt_tolerance_min: _,
            force_refresh_pubchem: _,
            force_refresh_kegg_conv: _,
        } = CoverageParams {
            selected_groups: None,
            presence_threshold: 0.5,
            dedup_enabled: true,
            dedup_rt_tolerance_min: 0.1,
            force_refresh_pubchem: false,
            force_refresh_kegg_conv: false,
        };
    }
}
