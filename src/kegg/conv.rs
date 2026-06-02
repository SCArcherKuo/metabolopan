//! High-level cache-first resolver: turn a list of PubChem CIDs into
//! `CID → Option<cpd_id>` via the KEGG `/conv/compound/pubchem` endpoint,
//! with a persistent per-entry-timestamped cache.

use anyhow::Result;
use chrono::Utc;
use std::collections::HashMap;
use std::sync::mpsc;

use crate::kegg::cache::{read_cid_to_cpd_cache, write_cid_to_cpd_cache};
use crate::kegg::client::{KEGG_CONV_MAX_BATCH, KeggClient};
use crate::kegg::types::{CidCpdEntry, ConvProgress};

/// Cache-first resolution of CID → Option<cpd_id>.
///
/// - `force_refresh = true` bypasses cache reads but still writes fresh
///   entries (so the refresh button can re-fetch every entry).
/// - Cache is written incrementally per batch.
/// - `progress_tx`, when `Some`, receives one event per batch (including
///   a final event at completion). Receiver-dropped sends are ignored.
pub async fn resolve_cids_to_cpds(
    client: &KeggClient,
    cids: &[String],
    force_refresh: bool,
    progress_tx: Option<mpsc::Sender<ConvProgress>>,
) -> Result<HashMap<String, Option<String>>> {
    let total_inputs = cids.len();
    let unique_inputs = crate::seq::dedupe_preserve_order(cids);

    let mut cache = read_cid_to_cpd_cache()?;
    let mut result: HashMap<String, Option<String>> = HashMap::new();
    let mut to_fetch: Vec<String> = Vec::new();

    if force_refresh {
        to_fetch.extend(unique_inputs.iter().cloned());
    } else {
        for c in &unique_inputs {
            if let Some(entry) = cache.get(c) {
                result.insert(c.clone(), entry.cpd.clone());
            } else {
                to_fetch.push(c.clone());
            }
        }
    }

    let from_cache = unique_inputs.len() - to_fetch.len();
    let total_batches = to_fetch.len().div_ceil(KEGG_CONV_MAX_BATCH);
    let mut fetched = 0usize;

    send_progress(
        &progress_tx,
        ConvProgress {
            completed_batches: 0,
            total_batches,
            from_cache,
            fetched,
            total_inputs,
        },
    );

    for (batch_idx, batch) in to_fetch.chunks(KEGG_CONV_MAX_BATCH).enumerate() {
        let pairs = client.conv_compound_pubchem(batch).await?;
        let now = Utc::now();

        for (cid, cpd) in &pairs {
            cache.insert(
                cid.clone(),
                CidCpdEntry {
                    cpd: cpd.clone(),
                    fetched_at: now,
                },
            );
            result.insert(cid.clone(), cpd.clone());
        }
        write_cid_to_cpd_cache(&cache)?;

        fetched += batch.len();
        send_progress(
            &progress_tx,
            ConvProgress {
                completed_batches: batch_idx + 1,
                total_batches,
                from_cache,
                fetched,
                total_inputs,
            },
        );
    }

    Ok(result)
}

fn send_progress(tx: &Option<mpsc::Sender<ConvProgress>>, event: ConvProgress) {
    if let Some(s) = tx
        && s.send(event).is_err()
    {
        // Receiver dropped; ignore.
    }
}

// (the `dedupe_preserve_order` test now lives in `src/seq.rs`, the helper's home)
