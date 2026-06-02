//! High-level cache-first resolver: turn a list of InChIKeys into
//! `InChIKey → Vec<CID>`, fetching missing entries from PubChem in
//! batches and writing positive + negative answers back to the persistent
//! cache.

use anyhow::Result;
use chrono::Utc;
use std::collections::HashMap;
use std::sync::mpsc;

use crate::pubchem::cache::{read_cache, write_cache};
use crate::pubchem::client::{PUBCHEM_MAX_BATCH, PubchemClient};
use crate::pubchem::types::{InchikeyCidsEntry, PubchemProgress};

/// Cache-first resolution of InChIKey → Vec<CID>.
///
/// - `force_refresh = true` bypasses cache reads (so refresh buttons can
///   re-fetch every entry), but still writes fresh entries.
/// - The cache is written incrementally per batch (atomic file write each
///   time) so a mid-run crash preserves progress made so far.
/// - `progress_tx`, when `Some`, receives one event per batch (including
///   a final event at completion). Failures to send are ignored so a
///   dropped UI receiver does not abort the resolver.
pub async fn resolve_inchikeys_to_cids(
    client: &PubchemClient,
    inchikeys: &[String],
    force_refresh: bool,
    progress_tx: Option<mpsc::Sender<PubchemProgress>>,
) -> Result<HashMap<String, Vec<String>>> {
    let total_inputs = inchikeys.len();

    // Dedupe inputs while preserving the first occurrence order. Callers
    // pass the per-feature InChIKey list which may have duplicates (e.g.
    // two features with the same compound annotation); we only need to
    // query each unique InChIKey once.
    let unique_inputs = crate::seq::dedupe_preserve_order(inchikeys);

    // Load the existing cache.
    let mut cache = read_cache()?;

    let mut result: HashMap<String, Vec<String>> = HashMap::new();
    let mut to_fetch: Vec<String> = Vec::new();

    if force_refresh {
        // Re-fetch everything; ignore cache reads. Cache writes still happen.
        to_fetch.extend(unique_inputs.iter().cloned());
    } else {
        for k in &unique_inputs {
            if let Some(entry) = cache.get(k) {
                result.insert(k.clone(), entry.cids.clone());
            } else {
                to_fetch.push(k.clone());
            }
        }
    }

    let from_cache = unique_inputs.len() - to_fetch.len();
    let total_batches = to_fetch.len().div_ceil(PUBCHEM_MAX_BATCH);
    let mut fetched = 0usize;

    // Emit an initial progress event so the UI can render "0 / N" before
    // the first batch lands.
    send_progress(
        &progress_tx,
        PubchemProgress {
            completed_batches: 0,
            total_batches,
            from_cache,
            fetched,
            total_inputs,
        },
    );

    for (batch_idx, batch) in to_fetch.chunks(PUBCHEM_MAX_BATCH).enumerate() {
        let pairs = client.post_inchikeys_to_cids(batch).await?;
        let now = Utc::now();

        // Write batch results to in-memory map AND persistent cache.
        for (inchikey, cids) in &pairs {
            cache.insert(
                inchikey.clone(),
                InchikeyCidsEntry {
                    cids: cids.clone(),
                    fetched_at: now,
                },
            );
            result.insert(inchikey.clone(), cids.clone());
        }
        write_cache(&cache)?;

        fetched += batch.len();
        send_progress(
            &progress_tx,
            PubchemProgress {
                completed_batches: batch_idx + 1,
                total_batches,
                from_cache,
                fetched,
                total_inputs,
            },
        );
    }

    // For callers passing duplicate InChIKeys, ensure the returned map's
    // values reflect the deduped results (already the case since we keyed
    // by InChIKey throughout).
    Ok(result)
}

fn send_progress(tx: &Option<mpsc::Sender<PubchemProgress>>, event: PubchemProgress) {
    if let Some(s) = tx
        && s.send(event).is_err()
    {
        // Receiver dropped; ignore.
    }
}

// (the `dedupe_preserve_order` test now lives in `src/seq.rs`, the helper's home)
