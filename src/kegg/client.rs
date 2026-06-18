use anyhow::{Context, Result};
use chrono::Utc;
use regex::Regex;
use reqwest::Url;
use std::sync::OnceLock;
use std::time::Duration;
use tokio::sync::mpsc;
use tokio::time::sleep;
use tracing::{debug, info, warn};

use crate::kegg::types::{
    KeggCompoundSet, KeggModuleEntry, KeggOrganism, KeggProgress, ModuleFetchProgress, SpeciesKegg,
};

const DEFAULT_BASE_URL: &str = "https://rest.kegg.jp";
/// Throttle between consecutive HTTP requests against rest.kegg.jp.
/// KEGG documents a soft cap of ~3 requests/second; 334 ms keeps us
/// under that. Earlier 50 ms was empirically too aggressive for /conv,
/// which started returning 403 (rate-limited) after a Stage 1 fetch.
const PER_REQUEST_DELAY: Duration = Duration::from_millis(334);

fn compound_id_re() -> &'static Regex {
    static RE: OnceLock<Regex> = OnceLock::new();
    RE.get_or_init(|| Regex::new(r"^C\d{5}$").expect("static regex compiles"))
}

#[derive(Debug, Clone)]
pub struct KeggClient {
    http: reqwest::Client,
    base_url: Url,
}

impl Default for KeggClient {
    fn default() -> Self {
        Self::new()
    }
}

impl KeggClient {
    pub fn new() -> Self {
        let base_url = Url::parse(DEFAULT_BASE_URL).expect("default base URL parses");
        Self::with_url_and_client(base_url, default_http_client())
    }

    pub fn with_base_url(base_url: Url) -> Self {
        Self::with_url_and_client(base_url, default_http_client())
    }

    fn with_url_and_client(base_url: Url, http: reqwest::Client) -> Self {
        Self { http, base_url }
    }

    pub fn base_url(&self) -> &Url {
        &self.base_url
    }

    fn endpoint(&self, path: &str) -> Result<Url> {
        // base_url may or may not end in '/'; join correctly.
        let mut joined = self.base_url.clone();
        let new_path = if joined.path().ends_with('/') {
            format!("{}{}", joined.path(), path.trim_start_matches('/'))
        } else {
            format!("{}/{}", joined.path(), path.trim_start_matches('/'))
        };
        joined.set_path(&new_path);
        Ok(joined)
    }

    /// One non-retrying KEGG GET: send → `error_for_status` → body text.
    /// Shared by the list/single-GET endpoints (`/list/*`, `/get/<pathway>`);
    /// these are single-attempt by design — the retrying loop is
    /// [`Self::get_with_retry`].
    async fn simple_get(&self, url: &Url) -> Result<String> {
        self.http
            .get(url.clone())
            .send()
            .await
            .with_context(|| format!("GET {url} failed"))?
            .error_for_status()
            .with_context(|| format!("GET {url} returned an error status"))?
            .text()
            .await
            .with_context(|| format!("failed to read body from {url}"))
    }

    pub async fn list_organisms(&self) -> Result<Vec<KeggOrganism>> {
        // The `/list/organism` special database was removed from the official
        // host (now HTTP 400); the BRITE "KEGG Organism" hierarchy `br08601`
        // carries the same roster with the taxonomy needed for the Group tree.
        let url = self.endpoint("get/br:br08601")?;
        let text = self.simple_get(&url).await?;
        Ok(parse_brite_organism_hierarchy(&text))
    }

    pub async fn list_pathways(&self, organism: &str) -> Result<Vec<(String, String)>> {
        let url = self.endpoint(&format!("list/pathway/{organism}"))?;
        let text = self.simple_get(&url).await?;
        Ok(parse_pathway_list(&text))
    }

    pub async fn get_pathway_detail(&self, pathway_id: &str) -> Result<String> {
        let url = self.endpoint(&format!("get/{pathway_id}"))?;
        self.simple_get(&url)
            .await
            .with_context(|| format!("failed to fetch KEGG pathway detail for {pathway_id}"))
    }

    pub async fn fetch_species_pathways(
        &self,
        organism: &str,
        progress_tx: mpsc::Sender<KeggProgress>,
    ) -> Result<SpeciesKegg> {
        info!(code = %organism, "fetching KEGG pathway list");
        let pathway_index = self.list_pathways(organism).await?;
        let total = pathway_index.len();
        info!(code = %organism, total, "starting per-pathway compound fetch");

        let mut pathways: Vec<KeggCompoundSet> = Vec::with_capacity(total);
        let mut channel_closed_logged = false;

        for (i, (id, name)) in pathway_index.iter().enumerate() {
            if i > 0 {
                sleep(PER_REQUEST_DELAY).await;
            }
            let detail = self
                .get_pathway_detail(id)
                .await
                .with_context(|| format!("failed to fetch detail for pathway {id}"))?;
            let compounds = parse_compound_ids(&detail);
            debug!(
                pathway = %id,
                compounds = compounds.len(),
                "parsed compound IDs for pathway"
            );
            pathways.push(KeggCompoundSet {
                id: id.clone(),
                name: name.clone(),
                compounds,
            });

            let progress = KeggProgress {
                completed: i + 1,
                total,
                current_pathway: id.clone(),
            };
            if let Err(e) = progress_tx.try_send(progress) {
                match e {
                    mpsc::error::TrySendError::Full(_) => {
                        // Drop progress event; UI will catch up on the next one.
                    }
                    mpsc::error::TrySendError::Closed(_) => {
                        if !channel_closed_logged {
                            warn!(code = %organism, "KEGG progress channel closed mid-fetch; continuing without progress updates");
                            channel_closed_logged = true;
                        }
                    }
                }
            }
        }

        Ok(SpeciesKegg {
            code: organism.to_string(),
            fetched_at: Utc::now(),
            pathways,
        })
    }

    /// GET `/list/module`. Returns `(module_id, name)` tuples. Module IDs
    /// in `/list/module` have NO prefix (bare `M00001`-`M01063`); the
    /// parser does NOT strip `md:`.
    pub async fn list_modules(&self) -> Result<Vec<(String, String)>> {
        let url = self.endpoint("list/module")?;
        let text = self.simple_get(&url).await?;
        Ok(parse_module_list(&text))
    }

    /// GET `/get/<module-id>`. Returns the raw flat-file body for
    /// downstream COMPOUND / COMPLETE parsing. Uses the same 403/5xx
    /// retry semantics as `/conv` (403 = rate-limit signal).
    pub async fn get_module_detail(&self, module_id: &str) -> Result<String> {
        let url = self.endpoint(&format!("get/{module_id}"))?;
        // 404 is fatal for a module GET; no throttle here — the caller
        // `fetch_modules_incremental` owns the inter-request sleep.
        match self
            .get_with_retry(&url, NotFound::Fatal, "KEGG /get/<module>")
            .await?
        {
            Some(body) => Ok(body),
            None => unreachable!("NotFound::Fatal never yields the all-no-match case"),
        }
    }

    /// One KEGG GET attempt, classified into a retry [`Attempt`]. Shared by the
    /// `/get/<module>` and `/conv` retry loops; the only per-call differences
    /// are the `label` (woven into messages) and the 404 policy (`on_404`).
    async fn send_once(&self, url: &Url, on_404: NotFound, label: &str) -> Attempt {
        let response = match self.http.get(url.clone()).send().await {
            Ok(r) => r,
            Err(e) if e.is_timeout() => {
                return Attempt::Retry(anyhow::anyhow!("{label} timed out: {e}"));
            }
            Err(e) => return Attempt::Fatal(anyhow::anyhow!("{label} request failed: {e}")),
        };
        let status = response.status();
        if status.as_u16() == 404 {
            return match on_404 {
                NotFound::AllNoMatch => Attempt::AllNoMatch,
                NotFound::Fatal => Attempt::Fatal(anyhow::anyhow!("{label} returned {status}")),
            };
        }
        // KEGG returns 403 as a rate-limit signal — must retry, not abort.
        if status.as_u16() == 403 {
            return Attempt::RateLimited(anyhow::anyhow!("{label} returned 403 (rate-limited)"));
        }
        if status.is_server_error() {
            return Attempt::Retry(anyhow::anyhow!("{label} returned {status}"));
        }
        if !status.is_success() {
            return Attempt::Fatal(anyhow::anyhow!("{label} returned {status}"));
        }
        match response.text().await {
            Ok(body) => Attempt::Ok(body),
            Err(e) => Attempt::Fatal(anyhow::anyhow!("{label} body read failed: {e}")),
        }
    }

    /// The shared 5-attempt KEGG retry loop: 403 → long (rate-limit) backoff,
    /// 5xx / timeout → short (network) backoff (both env-overridable). Returns
    /// `Ok(Some(body))` on success, `Ok(None)` for the 404→all-no-match case
    /// (only under `NotFound::AllNoMatch`), and `Err` on a fatal status or once
    /// the attempts are exhausted. The post-success throttle stays at the call
    /// site (only `/conv` sleeps), never inside this loop.
    async fn get_with_retry(
        &self,
        url: &Url,
        on_404: NotFound,
        label: &str,
    ) -> Result<Option<String>> {
        let mut last_err: Option<anyhow::Error> = None;
        for attempt in 1..=KEGG_CONV_MAX_ATTEMPTS {
            match self.send_once(url, on_404, label).await {
                Attempt::Ok(body) => return Ok(Some(body)),
                Attempt::AllNoMatch => return Ok(None),
                Attempt::RateLimited(e) => {
                    last_err = Some(e);
                    if attempt < KEGG_CONV_MAX_ATTEMPTS {
                        let backoff = kegg_conv_403_backoff();
                        warn!(
                            attempt,
                            max = KEGG_CONV_MAX_ATTEMPTS,
                            url = %url,
                            backoff_secs = backoff.as_secs(),
                            "{label} 403 (rate-limited); backing off"
                        );
                        sleep(backoff).await;
                    }
                }
                Attempt::Retry(e) => {
                    last_err = Some(e);
                    if attempt < KEGG_CONV_MAX_ATTEMPTS {
                        warn!(
                            attempt,
                            max = KEGG_CONV_MAX_ATTEMPTS,
                            url = %url,
                            "{label} transient error; retrying"
                        );
                        sleep(kegg_conv_network_backoff()).await;
                    }
                }
                Attempt::Fatal(e) => return Err(e),
            }
        }
        Err(last_err.unwrap_or_else(|| anyhow::anyhow!("{label} exhausted retries")))
    }

    /// Iterate `missing_ids`, fetch each via `get_module_detail`, parse
    /// both COMPOUND and COMPLETE blocks, and emit progress per module.
    /// Honours the 334 ms throttle between requests (KEGG etiquette).
    /// On any unrecoverable error from a single GET, returns early —
    /// modules fetched before the failure are NOT lost (the caller
    /// receives them in the returned Vec via `?` propagation only if
    /// the caller chooses to short-circuit; this fn itself bails on
    /// first hard failure).
    pub async fn fetch_modules_incremental(
        &self,
        missing_ids: &[String],
        progress_tx: tokio::sync::mpsc::Sender<ModuleFetchProgress>,
    ) -> Result<Vec<KeggModuleEntry>> {
        let total = missing_ids.len();
        let mut out: Vec<KeggModuleEntry> = Vec::with_capacity(total);
        let mut completed: usize = 0;
        // Rolling-average window for ETA. Capacity 10 per spec recommendation.
        let mut sample_buf: std::collections::VecDeque<f64> =
            std::collections::VecDeque::with_capacity(10);
        let mut channel_closed_logged = false;

        for (i, id) in missing_ids.iter().enumerate() {
            if i > 0 {
                sleep(PER_REQUEST_DELAY).await;
            }
            let started_at = std::time::Instant::now();
            let detail = self
                .get_module_detail(id)
                .await
                .with_context(|| format!("failed to fetch module {id}"))?;
            let elapsed = started_at.elapsed().as_secs_f64();
            if sample_buf.len() >= 10 {
                sample_buf.pop_front();
            }
            sample_buf.push_back(elapsed);

            let name = parse_module_name_from_detail(&detail);
            let compounds = parse_compound_ids(&detail);
            let complete_orgs = parse_complete_orgs(&detail);
            let entry = KeggModuleEntry {
                name,
                compounds,
                complete_orgs,
                fetched_at: Utc::now(),
            };
            out.push(entry);
            completed += 1;

            // ETA: only emit after ≥5 samples (warmup window). Include
            // retry-extended durations (Option A per design D16) — over-
            // estimation on a 12-min wait is friendlier than under-.
            let eta_secs = if sample_buf.len() >= 5 {
                let avg = sample_buf.iter().sum::<f64>() / sample_buf.len() as f64;
                let remaining = total.saturating_sub(completed) as f64;
                Some((remaining * avg).round() as u64)
            } else {
                None
            };

            let progress = ModuleFetchProgress {
                completed,
                total,
                current_id: id.clone(),
                eta_secs,
            };
            if let Err(e) = progress_tx.try_send(progress) {
                match e {
                    tokio::sync::mpsc::error::TrySendError::Full(_) => {
                        // Drop; UI will catch up on next tick.
                    }
                    tokio::sync::mpsc::error::TrySendError::Closed(_) => {
                        if !channel_closed_logged {
                            warn!("module fetch progress channel closed mid-fetch");
                            channel_closed_logged = true;
                        }
                    }
                }
            }
        }
        Ok(out)
    }

    /// GET `/conv/compound/pubchem:CID1+pubchem:CID2+...` for a batch of
    /// up to `KEGG_CONV_MAX_BATCH` PubChem CIDs. Returns one tuple per
    /// input CID: `(cid, Some(cpd))` when a `cpd:` mapping was returned,
    /// `(cid, None)` when the CID was queried but no `cpd:` line came
    /// back (filters out `glycan:`, `dr:`, etc.).
    ///
    /// KEGG returns HTTP 403 as a rate-limit signal (NOT a permanent
    /// forbidden). This method retries up to `KEGG_CONV_MAX_ATTEMPTS`
    /// times with `KEGG_CONV_403_BACKOFF` between attempts; 5xx / timeout
    /// retries use the same loop with a shorter `KEGG_CONV_NETWORK_BACKOFF`.
    /// 404 is treated as "all-no-match" (no retry).
    pub async fn conv_compound_pubchem(
        &self,
        cids: &[String],
    ) -> Result<Vec<(String, Option<String>)>> {
        if cids.is_empty() {
            return Ok(vec![]);
        }
        if cids.len() > KEGG_CONV_MAX_BATCH {
            anyhow::bail!(
                "batch size {} exceeds KEGG /conv maximum {}",
                cids.len(),
                KEGG_CONV_MAX_BATCH
            );
        }

        let joined = cids
            .iter()
            .map(|c| format!("pubchem:{c}"))
            .collect::<Vec<_>>()
            .join("+");
        let url = self.endpoint(&format!("conv/compound/{joined}"))?;

        // `/conv` treats 404 as "every input unmatched"; on success it applies
        // the 334 ms post-request throttle (KEGG etiquette) before returning.
        match self
            .get_with_retry(&url, NotFound::AllNoMatch, "KEGG /conv")
            .await?
        {
            Some(body) => {
                sleep(PER_REQUEST_DELAY).await;
                Ok(parse_conv_response(cids, &body))
            }
            None => Ok(no_match_for_all(cids)),
        }
    }
}

/// KEGG `/conv` accepts at most ~10 IDs per request.
pub const KEGG_CONV_MAX_BATCH: usize = 10;

/// Maximum attempts (initial + retries) for a single `/conv` batch.
const KEGG_CONV_MAX_ATTEMPTS: u32 = 5;
/// Default backoff between attempts when KEGG returns 403 (rate-limited).
/// Tests can override via the `KEGG_CONV_403_BACKOFF_MS` env var.
const KEGG_CONV_403_BACKOFF_DEFAULT: Duration = Duration::from_secs(5);
/// Default backoff between attempts on network/5xx errors. Tests can
/// override via `KEGG_CONV_NETWORK_BACKOFF_MS`.
const KEGG_CONV_NETWORK_BACKOFF_DEFAULT: Duration = Duration::from_secs(3);

fn kegg_conv_403_backoff() -> Duration {
    std::env::var("KEGG_CONV_403_BACKOFF_MS")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .map(Duration::from_millis)
        .unwrap_or(KEGG_CONV_403_BACKOFF_DEFAULT)
}

fn kegg_conv_network_backoff() -> Duration {
    std::env::var("KEGG_CONV_NETWORK_BACKOFF_MS")
        .ok()
        .and_then(|s| s.parse::<u64>().ok())
        .map(Duration::from_millis)
        .unwrap_or(KEGG_CONV_NETWORK_BACKOFF_DEFAULT)
}

/// 404 policy for a KEGG GET: `/get/<module>` treats 404 as fatal; `/conv`
/// treats it as "every input is unmatched" (no retry).
#[derive(Clone, Copy)]
enum NotFound {
    Fatal,
    AllNoMatch,
}

/// One classified KEGG GET attempt (shared by the `/get/<module>` and `/conv`
/// retry loops via `send_once` / `get_with_retry`).
enum Attempt {
    Ok(String),
    /// 404 under `NotFound::AllNoMatch` — every input is a confirmed no-match.
    AllNoMatch,
    /// HTTP 403 — KEGG's rate-limit signal. Retry with long backoff.
    RateLimited(anyhow::Error),
    /// 5xx / timeout. Retry with short backoff.
    Retry(anyhow::Error),
    /// Unrecoverable (4xx other than 403, 404-as-fatal, body read error).
    Fatal(anyhow::Error),
}

fn no_match_for_all(cids: &[String]) -> Vec<(String, Option<String>)> {
    cids.iter().map(|c| (c.clone(), None)).collect()
}

/// Parse a KEGG `/conv/compound/pubchem` response. Each line is
/// `pubchem:CID\tTARGET:ID`. We only keep `cpd:` targets; `glycan:`,
/// `dr:`, etc. are filtered out. Inputs missing from the response are
/// returned with `None`.
fn parse_conv_response(inputs: &[String], body: &str) -> Vec<(String, Option<String>)> {
    use std::collections::HashMap;
    use std::collections::hash_map::Entry;
    let mut matches: HashMap<String, String> = HashMap::new();
    for line in body.lines() {
        let parts: Vec<&str> = line.split('\t').collect();
        if parts.len() != 2 {
            continue;
        }
        let pubchem_part = parts[0].trim();
        let target_part = parts[1].trim();
        let cid = match pubchem_part.strip_prefix("pubchem:") {
            Some(s) => s.to_string(),
            None => continue,
        };
        let cpd = match target_part.strip_prefix("cpd:") {
            Some(s) => s.to_string(),
            None => continue, // filter out glycan:, dr:, etc.
        };
        // Deterministic first-write-wins. CidCpdEntry.cpd is Option<String>
        // (single value), so if KEGG ever returns multiple distinct cpd lines
        // for one CID we can't represent both — but at least surface the
        // conflict via WARN so we know when it happens in the wild.
        // Same-cpd duplicates stay silent.
        match matches.entry(cid.clone()) {
            Entry::Vacant(e) => {
                e.insert(cpd);
            }
            Entry::Occupied(e) => {
                let kept = e.get();
                if kept != &cpd {
                    warn!(
                        cid = %cid,
                        kept = %kept,
                        dropped = %cpd,
                        "KEGG /conv returned multiple cpd entries for one CID; keeping first"
                    );
                }
            }
        }
    }
    inputs
        .iter()
        .map(|cid| (cid.clone(), matches.remove(cid)))
        .collect()
}

fn default_http_client() -> reqwest::Client {
    // Same UA value as before (`<crate>/<version>`) → KEGG outbound bytes
    // unchanged; the shared builder simply centralises it.
    crate::cache_io::http_client(Duration::from_secs(30))
}

/// Parse the KEGG BRITE "KEGG Organism" hierarchy (`GET /get/br:br08601`) into
/// a flat organism list. The `/list/organism` special database was removed from
/// the official host (it now returns HTTP 400); `br08601` carries the same
/// roster as a taxonomy hierarchy whose leading character encodes the level:
/// `A` (Eukaryotes / Prokaryotes) → `B` (Animals / Bacteria / …) → `C`
/// (Mammals / …) → `D` (Primates / …) → `E` (leaf organism). `A`/`B`/`C`/`D`
/// lines name a group with a trailing ` (<count>)` occurrence count (stripped);
/// `E` leaves are `code  Scientific name (common name)`.
///
/// The reconstructed `lineage` is the accumulated `A`/`B`/`C`/`D` group names
/// joined by `;`, matching the semicolon-delimited form
/// [`build_organism_group_index`](crate::kegg::build_organism_group_index)
/// consumes. BRITE carries no KEGG T-numbers, so `t_number` is synthesized as
/// `T_{code}` (an inert placeholder; nothing downstream reads it as a real
/// identifier). Non-data lines (`+`, `!`, `#`, blank) are skipped.
pub fn parse_brite_organism_hierarchy(text: &str) -> Vec<KeggOrganism> {
    let mut out = Vec::new();
    // Current group name at each of the 4 taxonomy levels (A, B, C, D).
    let mut levels: [String; 4] = [String::new(), String::new(), String::new(), String::new()];
    for line in text.lines() {
        let mut chars = line.chars();
        let Some(prefix) = chars.next() else {
            continue; // blank line
        };
        let body = chars.as_str().trim();
        match prefix {
            'A' | 'B' | 'C' | 'D' => {
                let depth = (prefix as u8 - b'A') as usize;
                levels[depth] = strip_brite_count_suffix(body).to_string();
                // Clear deeper levels so a sibling branch can't inherit a stale
                // descendant name from the previous branch.
                for deeper in levels.iter_mut().skip(depth + 1) {
                    deeper.clear();
                }
            }
            'E' => {
                let mut parts = body.splitn(2, char::is_whitespace);
                let code = parts.next().unwrap_or("").trim().to_string();
                if code.is_empty() {
                    warn!(line = %line, "skipping BRITE organism leaf with empty code");
                    continue;
                }
                let name = parts.next().unwrap_or("").trim().to_string();
                let lineage = levels
                    .iter()
                    .filter(|l| !l.is_empty())
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(";");
                out.push(KeggOrganism {
                    t_number: format!("T_{code}"),
                    code,
                    name,
                    lineage,
                });
            }
            // '+' (column header), '!' (section marker), '#' (metadata), and any
            // other leading marker are non-data lines; skip.
            _ => continue,
        }
    }
    out
}

/// Strip a trailing ` (<digits>)` occurrence-count suffix that BRITE appends to
/// group names (e.g. `Animals (846)` → `Animals`). Only an all-digit
/// parenthesized suffix is removed, so a group name with other trailing parens
/// is left intact; leaf common-name parens (e.g. `(human)`) never reach here.
fn strip_brite_count_suffix(s: &str) -> &str {
    if let Some(open) = s.rfind(" (")
        && let Some(inner) = s[open + 2..].strip_suffix(')')
        && !inner.is_empty()
        && inner.bytes().all(|b| b.is_ascii_digit())
    {
        return s[..open].trim_end();
    }
    s
}

/// Parse `/list/module` response (tab-delimited 2-column: `M00001\tname`).
/// IDs in `/list/module` have NO `md:` prefix — DO NOT attempt to strip.
pub fn parse_module_list(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let mut parts = line.splitn(2, '\t');
        let id = parts.next().unwrap_or("").trim().to_string();
        let name = parts.next().unwrap_or("").trim().to_string();
        if id.is_empty() {
            warn!(line = %line, "skipping malformed module line (empty id)");
            continue;
        }
        out.push((id, name));
    }
    out
}

/// Extract a module's human-readable name from its `/get/<module>`
/// response by reading the `NAME` line. Falls back to an empty string
/// when no NAME line is present (defensive).
pub fn parse_module_name_from_detail(detail: &str) -> String {
    for line in detail.lines() {
        if let Some(rest) = line.strip_prefix("NAME") {
            return rest.trim().to_string();
        }
    }
    String::new()
}

pub fn parse_pathway_list(text: &str) -> Vec<(String, String)> {
    let mut out = Vec::new();
    for line in text.lines() {
        if line.trim().is_empty() {
            continue;
        }
        let mut parts = line.splitn(2, '\t');
        let raw_id = parts.next().unwrap_or("").trim();
        let name = parts.next().unwrap_or("").trim().to_string();
        let id = raw_id.strip_prefix("path:").unwrap_or(raw_id).to_string();
        if id.is_empty() {
            warn!(line = %line, "skipping malformed pathway line (empty id)");
            continue;
        }
        out.push((id, name));
    }
    out
}

/// Generic block-walker for KEGG `/get/{id}` flat-file responses. Scans
/// for `<keyword>` at column 0, collects the first whitespace-separated
/// token from that line and every indented continuation line, and
/// terminates at the next column-0 word-character keyword. Lines whose
/// column 0 is punctuation (`///`, etc.) do NOT terminate the block;
/// blank lines are treated as block-internal.
///
/// This generalizes the historical `parse_compound_ids` machine; the two
/// public wrappers (`parse_compound_ids` for `^C\d{5}$`-validated COMPOUND
/// IDs, `parse_complete_orgs` for unvalidated COMPLETE org codes) call it
/// with the appropriate keyword and post-process the tokens.
pub fn parse_block_first_tokens(detail: &str, keyword: &str) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    let mut inside = false;

    for line in detail.lines() {
        if inside {
            let mut chars = line.chars();
            match chars.next() {
                Some(c) if c.is_whitespace() => {
                    if let Some(token) = first_token(line) {
                        out.push(token.to_string());
                    }
                }
                Some(c) if is_word_char(c) => {
                    inside = false;
                }
                Some(_) => {
                    // Punctuation in column 0 — do NOT end the block.
                }
                None => {
                    // Empty line — block-internal, continue.
                }
            }
        }
        if !inside && line.starts_with(keyword) {
            // Require the keyword to be followed by whitespace or end-of-line
            // so e.g. "COMPLETEFOO" doesn't false-match "COMPLETE".
            let after = &line[keyword.len()..];
            if !after.is_empty() && !after.starts_with(char::is_whitespace) {
                continue;
            }
            inside = true;
            let rest = after.trim_start();
            if let Some(token) = first_token(rest) {
                out.push(token.to_string());
            }
        }
    }

    out
}

/// COMPOUND-block parser: returns compound IDs matching `^C\d{5}$`.
/// Non-matching tokens (e.g. `///`-line artifacts, malformed entries)
/// are dropped with a DEBUG log. Behaviour is identical to the
/// pre-refactor `parse_compound_ids` — implemented on top of
/// `parse_block_first_tokens` for sharing.
pub fn parse_compound_ids(detail: &str) -> Vec<String> {
    let re = compound_id_re();
    parse_block_first_tokens(detail, "COMPOUND")
        .into_iter()
        .filter(|token| {
            if re.is_match(token) {
                true
            } else {
                debug!(token = %token, "dropping non-`C#####` token from COMPOUND block");
                false
            }
        })
        .collect()
}

/// COMPLETE-block parser: returns organism codes (3-6 char strings like
/// `hsa`, `ath`). NO regex validation — KEGG org codes vary in length
/// and shape. Continuation lines look like
/// `<spaces><code>  <Scientific name> (<common name>)`; we keep the
/// leading code and drop everything else.
pub fn parse_complete_orgs(detail: &str) -> std::collections::HashSet<String> {
    parse_block_first_tokens(detail, "COMPLETE")
        .into_iter()
        .collect()
}

fn first_token(s: &str) -> Option<&str> {
    s.split_whitespace().next()
}

fn is_word_char(c: char) -> bool {
    c.is_ascii_alphanumeric() || c == '_'
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_compound_ids_single_line() {
        let detail = "ENTRY       gmx00010\nNAME        Glycolysis\nCOMPOUND    C00001 Water\nREFERENCE   ...\n";
        assert_eq!(parse_compound_ids(detail), vec!["C00001"]);
    }

    #[test]
    fn parse_compound_ids_multi_line() {
        let detail = "\
NAME        Foo
COMPOUND    C00001 Water
            C00002 ATP
            C00009 Orthophosphate
REFERENCE   bar
";
        assert_eq!(
            parse_compound_ids(detail),
            vec!["C00001", "C00002", "C00009"]
        );
    }

    #[test]
    fn parse_compound_ids_no_block() {
        let detail = "ENTRY       gmx99999\nNAME        Synthetic\nREFERENCE   ...\n";
        assert!(parse_compound_ids(detail).is_empty());
    }

    #[test]
    fn parse_compound_ids_filters_non_matching_token() {
        let detail = "COMPOUND    X12345 Not a compound\n            C00001 Water\nREFERENCE x\n";
        assert_eq!(parse_compound_ids(detail), vec!["C00001"]);
    }

    #[test]
    fn parse_compound_ids_punctuation_terminator_is_not_a_terminator() {
        // `///` begins with punctuation, not a word char — block does NOT end.
        let detail = "COMPOUND    C00001 Water\n///\n";
        assert_eq!(parse_compound_ids(detail), vec!["C00001"]);
    }

    #[test]
    fn parse_compound_ids_block_ends_on_word_keyword() {
        let detail = "\
COMPOUND    C00001 Water
            C00002 ATP
REFERENCE   x
            C99999 should NOT be picked up
";
        assert_eq!(parse_compound_ids(detail), vec!["C00001", "C00002"]);
    }

    #[test]
    fn parse_pathway_list_strips_path_prefix() {
        let text =
            "path:gmx00010\tGlycolysis - Glycine max (soybean)\npath:gmx00020\tCitrate cycle\n";
        let out = parse_pathway_list(text);
        assert_eq!(out.len(), 2);
        assert_eq!(out[0].0, "gmx00010");
        assert_eq!(out[1].0, "gmx00020");
    }

    /// A minimal `br08601` BRITE excerpt with the real wrapper lines, two
    /// sibling level-4 branches under one level-3, and a second level-2 branch
    /// to exercise the deeper-level clearing logic.
    const BRITE_SAMPLE: &str = "\
+E\tKEGG Organism
!
AEukaryotes (1307)
B  Animals (846)
C    Mammals (184)
D      Primates (33)
E        hsa  Homo sapiens (human)
E        ptr  Pan troglodytes (chimpanzee)
D      Rodents (29)
E        mmu  Mus musculus (house mouse)
B  Plants (200)
C    Eudicots (100)
D      Brassicales (10)
E        ath  Arabidopsis thaliana (thale cress)
#Last updated: June 18, 2026
";

    #[test]
    fn parse_brite_reconstructs_lineage_and_synthesizes_t_number() {
        let out = parse_brite_organism_hierarchy(BRITE_SAMPLE);
        assert_eq!(out.len(), 4);
        // First leaf: full 4-level lineage, count suffixes stripped, common-name
        // parens preserved, synthetic `T_{code}`.
        assert_eq!(out[0].code, "hsa");
        assert_eq!(out[0].t_number, "T_hsa");
        assert_eq!(out[0].name, "Homo sapiens (human)");
        assert_eq!(out[0].lineage, "Eukaryotes;Animals;Mammals;Primates");
    }

    #[test]
    fn parse_brite_clears_deeper_levels_across_sibling_branches() {
        let out = parse_brite_organism_hierarchy(BRITE_SAMPLE);
        let by_code = |c: &str| out.iter().find(|o| o.code == c).expect("present").clone();
        // Sibling D-branch under the same C: lineage swaps Primates → Rodents.
        assert_eq!(by_code("mmu").lineage, "Eukaryotes;Animals;Mammals;Rodents");
        // New B-branch (Plants) must not inherit Animals' C/D descendants.
        assert_eq!(
            by_code("ath").lineage,
            "Eukaryotes;Plants;Eudicots;Brassicales"
        );
    }

    #[test]
    fn parse_brite_skips_non_data_lines() {
        // The `+`, `!`, and `#` wrapper lines in BRITE_SAMPLE produce no records.
        let out = parse_brite_organism_hierarchy(BRITE_SAMPLE);
        assert!(out.iter().all(|o| !o.code.is_empty()));
        assert_eq!(out.len(), 4); // only the four E leaves
    }

    #[test]
    fn strip_brite_count_suffix_only_strips_digit_parens() {
        assert_eq!(strip_brite_count_suffix("Animals (846)"), "Animals");
        assert_eq!(strip_brite_count_suffix("Eukaryotes (1307)"), "Eukaryotes");
        // A non-digit parenthetical is left intact (defensive).
        assert_eq!(strip_brite_count_suffix("Group (TBD)"), "Group (TBD)");
        assert_eq!(strip_brite_count_suffix("NoParens"), "NoParens");
    }

    #[test]
    fn iso8601_round_trip_has_utc_suffix() {
        use chrono::Utc;
        let species = SpeciesKegg {
            code: "gmx".into(),
            fetched_at: Utc::now(),
            pathways: vec![],
        };
        let s = serde_json::to_string(&species).unwrap();
        // chrono's default serde format emits RFC-3339; assert UTC suffix is one of Z / +00:00.
        assert!(
            s.contains("Z\"") || s.contains("+00:00\""),
            "expected UTC suffix in {s}"
        );
        let back: SpeciesKegg = serde_json::from_str(&s).unwrap();
        assert_eq!(back, species);
    }

    #[test]
    fn parse_conv_standard_response() {
        let body = "pubchem:5793\tcpd:C00031\npubchem:99999\tcpd:C00074\n";
        let inputs = vec!["5793".to_string(), "12345".to_string(), "99999".to_string()];
        let parsed = parse_conv_response(&inputs, body);
        assert_eq!(parsed.len(), 3);
        assert_eq!(parsed[0], ("5793".to_string(), Some("C00031".to_string())));
        assert_eq!(parsed[1], ("12345".to_string(), None));
        assert_eq!(parsed[2], ("99999".to_string(), Some("C00074".to_string())));
    }

    #[test]
    fn parse_conv_filters_non_cpd_lines() {
        // Mix of cpd, glycan, dr — only cpd lines should survive.
        let body = "pubchem:5793\tcpd:C00031\n\
                    pubchem:11111\tglycan:G00100\n\
                    pubchem:22222\tdr:D00001\n";
        let inputs = vec!["5793".to_string(), "11111".to_string(), "22222".to_string()];
        let parsed = parse_conv_response(&inputs, body);
        assert_eq!(parsed[0].1.as_deref(), Some("C00031"));
        assert_eq!(parsed[1].1, None, "glycan line must be filtered");
        assert_eq!(parsed[2].1, None, "drug line must be filtered");
    }

    #[test]
    fn parse_conv_empty_body_all_no_match() {
        let inputs = vec!["A".to_string(), "B".to_string()];
        let parsed = parse_conv_response(&inputs, "");
        assert_eq!(parsed.len(), 2);
        assert!(parsed[0].1.is_none());
        assert!(parsed[1].1.is_none());
    }

    #[test]
    fn parse_conv_preserves_input_order() {
        let body = "pubchem:C\tcpd:C00003\npubchem:A\tcpd:C00001\n";
        let inputs = vec!["A".to_string(), "B".to_string(), "C".to_string()];
        let parsed = parse_conv_response(&inputs, body);
        assert_eq!(parsed[0], ("A".to_string(), Some("C00001".to_string())));
        assert_eq!(parsed[1], ("B".to_string(), None));
        assert_eq!(parsed[2], ("C".to_string(), Some("C00003".to_string())));
    }

    // ─── Track B: generic block parser + module parsers ──────────────────

    #[test]
    fn parse_block_first_tokens_multi_line() {
        let detail = "\
COMPLETE    hsa  Homo sapiens (human)
            ptr  Pan troglodytes (chimpanzee)
            pps  Pan paniscus (bonobo)
REACTION    R01786
";
        let tokens = parse_block_first_tokens(detail, "COMPLETE");
        assert_eq!(tokens, vec!["hsa", "ptr", "pps"]);
    }

    #[test]
    fn parse_block_first_tokens_punctuation_does_not_terminate() {
        let detail = "\
COMPOUND    C00001 Water
///
";
        let tokens = parse_block_first_tokens(detail, "COMPOUND");
        // `///` is punctuation in column 0, NOT a word char, so block
        // continues — but `///` is also not a continuation (no leading
        // whitespace), so no token is harvested from that line.
        assert_eq!(tokens, vec!["C00001"]);
    }

    #[test]
    fn parse_block_first_tokens_missing_keyword() {
        let detail = "ENTRY M00001\nNAME Glycolysis\nREFERENCE foo\n";
        let tokens = parse_block_first_tokens(detail, "COMPLETE");
        assert!(tokens.is_empty());
    }

    #[test]
    fn parse_block_first_tokens_keyword_prefix_match_is_strict() {
        // "COMPLETEFOO" must NOT false-match "COMPLETE".
        let detail = "COMPLETEFOO    not-a-match\nREFERENCE bar\n";
        let tokens = parse_block_first_tokens(detail, "COMPLETE");
        assert!(tokens.is_empty());
    }

    #[test]
    fn parse_complete_orgs_collects_org_codes_unvalidated() {
        let detail = "\
COMPLETE    hsa  Homo sapiens (human)
            cnb  Cryptococcus deneoformans
            longorgcode  Some Organism
REACTION    R01
";
        let set = parse_complete_orgs(detail);
        assert_eq!(set.len(), 3);
        assert!(set.contains("hsa"));
        assert!(set.contains("cnb"));
        assert!(set.contains("longorgcode"));
    }

    #[test]
    fn parse_complete_orgs_empty_when_block_absent() {
        let detail = "ENTRY M00099\nNAME Synthetic\nDEFINITION foo\n";
        let set = parse_complete_orgs(detail);
        assert!(set.is_empty());
    }

    #[test]
    fn parse_compound_ids_via_generalized_parser_still_validates() {
        // The wrapper still drops non-`C\d{5}$` tokens, exercising the
        // post-processing filter on top of parse_block_first_tokens.
        let detail = "\
COMPOUND    C00001 Water
            X12345 not-a-compound
            C00002 ATP
REFERENCE x
";
        let ids = parse_compound_ids(detail);
        assert_eq!(ids, vec!["C00001", "C00002"]);
    }

    #[test]
    fn parse_module_list_standard_response() {
        let text = "\
M00001\tGlycolysis (Embden-Meyerhof pathway), glucose => pyruvate
M00002\tGlycolysis, core module
M00010\tCitrate cycle, first carbon oxidation
";
        let out = parse_module_list(text);
        assert_eq!(out.len(), 3);
        assert_eq!(out[0].0, "M00001");
        assert!(out[0].1.starts_with("Glycolysis (Embden"));
        assert_eq!(out[2].0, "M00010");
    }

    #[test]
    fn parse_module_list_ids_have_no_md_prefix() {
        // `/list/module` returns bare IDs (NO `md:`). Confirm we don't
        // accidentally strip a prefix that isn't there or add one.
        let text = "M00001\tName\n";
        let out = parse_module_list(text);
        assert_eq!(out[0].0, "M00001");
    }

    #[test]
    fn parse_module_name_from_detail_finds_name_line() {
        let detail = "\
ENTRY       M00001            Pathway   Module
NAME        Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate
DEFINITION  (K00001) ...
";
        let name = parse_module_name_from_detail(detail);
        assert_eq!(
            name,
            "Glycolysis (Embden-Meyerhof pathway), glucose => pyruvate"
        );
    }

    #[test]
    fn parse_module_name_from_detail_missing_returns_empty() {
        let detail = "ENTRY M00099\nDEFINITION (K0)\n";
        let name = parse_module_name_from_detail(detail);
        assert_eq!(name, "");
    }

    #[test]
    fn parse_conv_response_first_write_wins_on_conflicting_cpd() {
        // KEGG returns two distinct cpd entries for the same CID. The cache
        // type (CidCpdEntry.cpd: Option<String>) holds only one value, so
        // we deterministically keep the first; the WARN log surfaces the
        // dropped value. This test only asserts behaviour, not the WARN
        // itself (would require tracing_test as a dev-dep).
        let body = "pubchem:5793\tcpd:C00031\npubchem:5793\tcpd:C00267\n";
        let inputs = vec!["5793".to_string()];
        let result = parse_conv_response(&inputs, body);
        assert_eq!(result.len(), 1);
        assert_eq!(result[0].0, "5793");
        assert_eq!(
            result[0].1,
            Some("C00031".to_string()),
            "first cpd must be kept (deterministic), not last"
        );
    }

    #[test]
    fn parse_conv_response_silent_on_identical_duplicate_cpd() {
        // Same cpd duplicated for the same CID → no conflict, no WARN.
        let body = "pubchem:5793\tcpd:C00031\npubchem:5793\tcpd:C00031\n";
        let inputs = vec!["5793".to_string()];
        let result = parse_conv_response(&inputs, body);
        assert_eq!(result[0].1, Some("C00031".to_string()));
    }

    #[test]
    fn parse_conv_response_handles_distinct_cids_normally() {
        // Two CIDs, each maps to one cpd. No conflict, no WARN.
        let body = "pubchem:5793\tcpd:C00031\npubchem:100\tcpd:C00100\n";
        let inputs = vec!["5793".to_string(), "100".to_string(), "999".to_string()];
        let result = parse_conv_response(&inputs, body);
        assert_eq!(result.len(), 3);
        assert_eq!(result[0].1, Some("C00031".to_string()));
        assert_eq!(result[1].1, Some("C00100".to_string()));
        assert_eq!(result[2].1, None, "unmatched CID must be None");
    }
}
