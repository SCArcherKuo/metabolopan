pub mod cache;
pub mod client;
pub mod conv;
pub mod groups;
pub mod types;

use anyhow::Result;
use tokio::sync::mpsc;
use tracing::info;

pub use cache::{cache_dir, clear_stale_locks, invalidate_cache};
pub use client::KeggClient;
pub use conv::resolve_cids_to_cpds;
pub use groups::build_organism_group_index;
pub use types::{
    CidCpdEntry, ConvProgress, KeggCacheScope, KeggCompoundSet, KeggEvent, KeggModuleEntry,
    KeggModulesCache, KeggOrganism, KeggProgress, ModuleFetchProgress, OrganismGroupIndex,
    OrganismsCache, SpeciesKegg,
};

/// Public entry point used by the GUI. Reads the organism cache if present;
/// otherwise issues the HTTP request and writes the cache before returning.
///
/// Either path also rewrites the derived `organism_groups.json` precompute
/// cache so the two stay coherent (Track C / kegg-fetching spec). The
/// group-index `fetched_at` MUST equal the organisms cache's `fetched_at`.
pub async fn list_organisms(client: &KeggClient) -> Result<Vec<KeggOrganism>> {
    if let Some(cached) = cache::read_organisms()? {
        info!(
            count = cached.organisms.len(),
            fetched_at = %cached.fetched_at,
            "loaded KEGG organism list from cache"
        );
        let index = build_organism_group_index(&cached.organisms, cached.fetched_at);
        if let Err(e) = cache::write_organism_group_index(&index) {
            tracing::warn!(error = %e, "failed to persist organism group index");
        }
        return Ok(cached.organisms);
    }

    let organisms = client.list_organisms().await?;
    let fetched_at = chrono::Utc::now();
    let cache_entry = OrganismsCache {
        fetched_at,
        organisms: organisms.clone(),
    };
    cache::write_organisms(&cache_entry)?;
    let index = build_organism_group_index(&organisms, fetched_at);
    if let Err(e) = cache::write_organism_group_index(&index) {
        tracing::warn!(error = %e, "failed to persist organism group index");
    }
    info!(
        count = organisms.len(),
        "fetched KEGG organism list from REST and cached"
    );
    Ok(organisms)
}

/// Public entry point for fetching a species' pathway/compound data. Reads
/// from cache when available (no progress events emitted); otherwise hits
/// the REST API, streams progress through `progress_tx`, and persists the
/// result on success.
pub async fn fetch_species_pathways(
    client: &KeggClient,
    code: &str,
    progress_tx: mpsc::Sender<KeggProgress>,
) -> Result<SpeciesKegg> {
    if let Some(cached) = cache::read_species(code)? {
        info!(
            code = %code,
            pathways = cached.pathways.len(),
            fetched_at = %cached.fetched_at,
            "loaded KEGG species data from cache"
        );
        return Ok(cached);
    }

    let species = client.fetch_species_pathways(code, progress_tx).await?;
    cache::write_species(&species)?;
    info!(
        code = %code,
        pathways = species.pathways.len(),
        "fetched KEGG species data from REST and cached"
    );
    Ok(species)
}

/// Public entry point for fetching the global KEGG modules cache.
///
/// Flow:
/// 1. Acquire the long-running `.modules.lock` (waits up to 30 min for any
///    concurrent fetch; treats stale heartbeat >90 s as orphaned).
/// 2. Read existing `modules.json` (empty if missing).
/// 3. `GET /list/module` to enumerate all ~573 module IDs (the `M01063`
///    upper bound is the last ID, not the count — the range is sparse).
/// 4. Compute `missing_ids` = listed IDs absent from cache (OR all when
///    `force_refresh = true`).
/// 5. Fetch each missing module via `get_module_detail` with the 334 ms
///    throttle + 403/5xx retry policy; emit `ModuleFetchProgress` per
///    module.
/// 6. Merge fresh entries into cache; on `force_refresh`, also prune
///    cache entries no longer in `/list/module`.
/// 7. Persist atomically; release the lock (RAII guard's Drop).
pub async fn fetch_modules(
    client: &KeggClient,
    force_refresh: bool,
    progress_tx: mpsc::Sender<types::ModuleFetchProgress>,
) -> Result<KeggModulesCache> {
    let mut guard = cache::acquire_modules_fetch_lock()?;
    guard.heartbeat()?;

    let mut existing = cache::read_modules_cache()?;
    let listing = client.list_modules().await?;
    guard.heartbeat()?;

    let listed_ids: std::collections::HashSet<String> =
        listing.iter().map(|(id, _)| id.clone()).collect();
    let missing_ids: Vec<String> = if force_refresh {
        listing.iter().map(|(id, _)| id.clone()).collect()
    } else {
        listing
            .iter()
            .filter_map(|(id, _)| {
                if existing.modules.contains_key(id) {
                    None
                } else {
                    Some(id.clone())
                }
            })
            .collect()
    };
    info!(
        listed = listing.len(),
        missing = missing_ids.len(),
        force_refresh,
        "starting KEGG modules fetch"
    );

    // Emit an initial progress event so the ModulesFetching UI shows the
    // actual missing-count immediately (instead of a misleading "0 / 0"
    // while waiting for the first per-module fetch). For warm-cache
    // (`missing_ids.is_empty()`), this is the only event — the loop below
    // runs zero iterations and the orchestrator returns immediately.
    let initial_progress = types::ModuleFetchProgress {
        completed: 0,
        total: missing_ids.len(),
        current_id: if missing_ids.is_empty() {
            "(cache complete)".to_string()
        } else {
            "(starting fetch)".to_string()
        },
        eta_secs: None,
    };
    let _ = progress_tx.send(initial_progress).await;

    // Heartbeat-aware fetch loop. We pre-build a sub-channel and tee
    // progress through so the heartbeat thread can also kick.
    let entries = client
        .fetch_modules_incremental(&missing_ids, progress_tx)
        .await?;

    // Merge fresh entries into the cache.
    for (id, entry) in missing_ids.iter().zip(entries) {
        existing.modules.insert(id.clone(), entry);
    }
    // On force_refresh, also prune entries no longer in /list/module.
    if force_refresh {
        existing.modules.retain(|id, _| listed_ids.contains(id));
    }

    guard.heartbeat()?;
    cache::write_modules_cache(&existing)?;
    info!(
        cached = existing.modules.len(),
        "KEGG modules cache persisted"
    );
    Ok(existing)
}
