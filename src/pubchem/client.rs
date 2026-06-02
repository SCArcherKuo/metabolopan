//! PubChem PUG REST client. Batches InChIKeys against the
//! `compound/inchikey/property/InChIKey/CSV` POST endpoint.

use anyhow::{Context, Result, anyhow};
use std::collections::HashMap;
use std::time::Duration;

const DEFAULT_BASE_URL: &str = "https://pubchem.ncbi.nlm.nih.gov/rest/pug";
const ENDPOINT_PATH: &str = "/compound/inchikey/property/InChIKey/CSV";

/// Maximum number of InChIKeys to send in a single POST body. PubChem's
/// soft limit; 200 keeps the request payload < 8 KB on typical InChIKey
/// lengths.
pub const PUBCHEM_MAX_BATCH: usize = 200;

/// HTTP timeout per request (single attempt).
const REQUEST_TIMEOUT: Duration = Duration::from_secs(60);

/// Backoff between attempt 1 and attempt 2 on 5xx / timeout.
const RETRY_BACKOFF: Duration = Duration::from_millis(250);

pub struct PubchemClient {
    http: reqwest::Client,
    base_url: String,
}

impl PubchemClient {
    pub fn new() -> Self {
        Self {
            http: crate::cache_io::http_client(REQUEST_TIMEOUT),
            base_url: DEFAULT_BASE_URL.to_string(),
        }
    }

    /// Override the base URL (used by tests with wiremock).
    pub fn with_base_url(base_url: impl Into<String>) -> Self {
        Self {
            http: crate::cache_io::http_client(REQUEST_TIMEOUT),
            base_url: base_url.into(),
        }
    }

    /// POST a batch of up to `PUBCHEM_MAX_BATCH` InChIKeys to PubChem and
    /// return the parsed result. Every input InChIKey appears in the
    /// returned vector, including those PubChem found no CID for (with
    /// `Vec<CID> = []`).
    ///
    /// Retries once on HTTP 5xx or timeout with a `RETRY_BACKOFF` pause.
    /// HTTP 200 with empty body AND HTTP 404 are both treated as confirmed
    /// no-match for every input in the batch (no retry, no error) — the
    /// 404 path mirrors the KEGG /conv client's `ConvAttempt::AllNoMatch`.
    /// Other 4xx statuses are surfaced as `SendError::Fatal`.
    pub async fn post_inchikeys_to_cids(
        &self,
        inchikeys: &[String],
    ) -> Result<Vec<(String, Vec<String>)>> {
        if inchikeys.is_empty() {
            return Ok(vec![]);
        }
        if inchikeys.len() > PUBCHEM_MAX_BATCH {
            anyhow::bail!(
                "batch size {} exceeds PubChem maximum {}",
                inchikeys.len(),
                PUBCHEM_MAX_BATCH
            );
        }

        let url = format!("{}{}", self.base_url, ENDPOINT_PATH);
        let body = format!("inchikey={}", inchikeys.join(","));

        let body_text = match self.send_once(&url, &body).await {
            Ok(text) => text,
            Err(SendError::AllNoMatch) => {
                // Every input is a confirmed no-match; identical to the
                // empty-body 200 path inside parse_csv_response.
                return Ok(inchikeys.iter().map(|k| (k.clone(), vec![])).collect());
            }
            Err(SendError::Retryable(_)) => {
                tokio::time::sleep(RETRY_BACKOFF).await;
                match self.send_once(&url, &body).await {
                    Ok(text) => text,
                    Err(SendError::AllNoMatch) => {
                        return Ok(inchikeys.iter().map(|k| (k.clone(), vec![])).collect());
                    }
                    Err(SendError::Retryable(e)) => return Err(e),
                    Err(SendError::Fatal(e)) => return Err(e),
                }
            }
            Err(SendError::Fatal(e)) => return Err(e),
        };

        Ok(parse_csv_response(inchikeys, &body_text))
    }

    /// Issue one HTTP attempt and classify the result.
    async fn send_once(&self, url: &str, body: &str) -> std::result::Result<String, SendError> {
        let response = self
            .http
            .post(url)
            .header("Content-Type", "application/x-www-form-urlencoded")
            .body(body.to_string())
            .send()
            .await;

        let response = match response {
            Ok(r) => r,
            Err(e) if e.is_timeout() => {
                return Err(SendError::Retryable(anyhow!("PubChem POST timed out: {e}")));
            }
            Err(e) => return Err(SendError::Fatal(anyhow!("PubChem POST failed: {e}"))),
        };

        let status = response.status();
        if status.as_u16() == 404 {
            // PubChem PUG REST returns 404 (PUGREST.NotFound) when the entire
            // batch matches zero CIDs. Treat as confirmed no-match for every
            // input, identical to the empty-body 200 case. Mirrors KEGG /conv
            // (`src/kegg/client.rs:451`).
            return Err(SendError::AllNoMatch);
        }
        if status.is_server_error() {
            return Err(SendError::Retryable(anyhow!(
                "PubChem POST returned {}",
                status
            )));
        }
        if !status.is_success() {
            return Err(SendError::Fatal(anyhow!(
                "PubChem POST returned {}",
                status
            )));
        }
        response
            .text()
            .await
            .context("read PubChem response body")
            .map_err(SendError::Fatal)
    }
}

impl Default for PubchemClient {
    fn default() -> Self {
        Self::new()
    }
}

enum SendError {
    Retryable(anyhow::Error),
    Fatal(anyhow::Error),
    /// HTTP 404 = PubChem confirms zero matches for every InChIKey in the
    /// batch (documented PUGREST.NotFound response). Semantically identical
    /// to an HTTP 200 with empty body; caller should map each input to
    /// `(key, vec![])`. Mirrors `ConvAttempt::AllNoMatch` in the KEGG /conv
    /// client at `src/kegg/client.rs:451`.
    AllNoMatch,
}

/// Parse the CSV body returned by PubChem `property/InChIKey/CSV`.
///
/// Format (one header row, then one data row per match):
///
/// ```text
/// "InChIKey","CID"
/// "ZSLZBFCDCINBPY-ZSJDYOACSA-N","5793"
/// "FOOBARBAZQUX-XXXX-N","12345"
/// "FOOBARBAZQUX-XXXX-N","67890"
/// ```
///
/// Multiple CIDs for one InChIKey appear as separate rows; we aggregate.
/// Input InChIKeys with no rows in the response are returned with an
/// empty CID list (confirmed "no match").
fn parse_csv_response(inputs: &[String], body: &str) -> Vec<(String, Vec<String>)> {
    let mut matches: HashMap<String, Vec<String>> = HashMap::new();
    let mut lines = body.lines();

    // Parse header to figure out column order. PubChem returns
    // `"CID","InChIKey"` (the actual server response), but the official
    // docs and some mirrors return `"InChIKey","CID"`. We must NOT assume.
    let Some(header) = lines.next() else {
        // Empty body → no matches; caller will still emit all inputs as no-match.
        return inputs.iter().map(|k| (k.clone(), vec![])).collect();
    };
    let header_cols: Vec<String> = header
        .split(',')
        .map(|s| s.trim().trim_matches('"').to_ascii_lowercase())
        .collect();
    let inchikey_idx = header_cols.iter().position(|c| c == "inchikey");
    let cid_idx = header_cols.iter().position(|c| c == "cid");
    let (Some(inchikey_idx), Some(cid_idx)) = (inchikey_idx, cid_idx) else {
        // Header doesn't have the columns we expect — treat as all-no-match
        // rather than misparsing. Logging is the caller's job.
        return inputs.iter().map(|k| (k.clone(), vec![])).collect();
    };

    for line in lines {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let parts: Vec<&str> = line.split(',').collect();
        if parts.len() <= inchikey_idx.max(cid_idx) {
            continue;
        }
        let inchikey = parts[inchikey_idx].trim().trim_matches('"').to_string();
        let cid = parts[cid_idx].trim().trim_matches('"').to_string();
        if inchikey.is_empty() || cid.is_empty() {
            continue;
        }
        matches.entry(inchikey).or_default().push(cid);
    }

    inputs
        .iter()
        .map(|inchikey| {
            let cids = matches.remove(inchikey).unwrap_or_default();
            (inchikey.clone(), cids)
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_standard_response() {
        let body = "\"InChIKey\",\"CID\"\n\"K1\",\"5793\"\n\"K2\",\"100\"\n\"K2\",\"200\"\n";
        let inputs = vec!["K1".to_string(), "K2".to_string(), "K3".to_string()];
        let parsed = parse_csv_response(&inputs, body);

        assert_eq!(parsed.len(), 3);
        assert_eq!(parsed[0].0, "K1");
        assert_eq!(parsed[0].1, vec!["5793".to_string()]);
        assert_eq!(parsed[1].0, "K2");
        assert_eq!(parsed[1].1, vec!["100".to_string(), "200".to_string()]);
        assert_eq!(parsed[2].0, "K3");
        assert!(parsed[2].1.is_empty(), "K3 had no row, should be no-match");
    }

    #[test]
    fn parse_empty_body_treats_all_inputs_as_no_match() {
        let inputs = vec!["K1".to_string(), "K2".to_string()];
        let parsed = parse_csv_response(&inputs, "");
        assert_eq!(parsed.len(), 2);
        assert!(parsed[0].1.is_empty());
        assert!(parsed[1].1.is_empty());
    }

    #[test]
    fn parse_quoted_fields_stripped() {
        let body = "\"InChIKey\",\"CID\"\n\"K1\",\"5793\"\n";
        let inputs = vec!["K1".to_string()];
        let parsed = parse_csv_response(&inputs, body);
        assert_eq!(parsed[0].1, vec!["5793".to_string()]);
    }

    #[test]
    fn parse_preserves_input_order_and_count() {
        // Inputs intentionally not in alphabetical order; response in
        // arbitrary order. Output order MUST match input order.
        let body = "\"InChIKey\",\"CID\"\n\"K3\",\"3\"\n\"K1\",\"1\"\n";
        let inputs = vec!["K1".to_string(), "K2".to_string(), "K3".to_string()];
        let parsed = parse_csv_response(&inputs, body);
        assert_eq!(parsed[0].0, "K1");
        assert_eq!(parsed[0].1, vec!["1".to_string()]);
        assert_eq!(parsed[1].0, "K2");
        assert!(parsed[1].1.is_empty());
        assert_eq!(parsed[2].0, "K3");
        assert_eq!(parsed[2].1, vec!["3".to_string()]);
    }

    #[test]
    fn parse_actual_pubchem_header_order_cid_first() {
        // The REAL PubChem PUG REST endpoint returns "CID","InChIKey"
        // (NOT "InChIKey","CID"). Our parser must handle either order.
        let body = "\"CID\",\"InChIKey\"\n\
                    7577525,\"YJRSMJTVXWBFJJ-UHFFFAOYSA-N\"\n\
                    16406293,\"VHMCLONZXOFDIQ-QPQITGAISA-N\"\n";
        let inputs = vec![
            "YJRSMJTVXWBFJJ-UHFFFAOYSA-N".to_string(),
            "VHMCLONZXOFDIQ-QPQITGAISA-N".to_string(),
            "MISSING-NOMATCH-N".to_string(),
        ];
        let parsed = parse_csv_response(&inputs, body);
        assert_eq!(parsed[0].0, "YJRSMJTVXWBFJJ-UHFFFAOYSA-N");
        assert_eq!(parsed[0].1, vec!["7577525".to_string()]);
        assert_eq!(parsed[1].0, "VHMCLONZXOFDIQ-QPQITGAISA-N");
        assert_eq!(parsed[1].1, vec!["16406293".to_string()]);
        assert_eq!(parsed[2].0, "MISSING-NOMATCH-N");
        assert!(parsed[2].1.is_empty());
    }

    #[test]
    fn parse_unrecognised_header_treats_all_as_no_match() {
        // Defensive: if neither column header matches, don't misparse.
        let body = "\"Foo\",\"Bar\"\n\"VALUE-A\",\"VALUE-B\"\n";
        let inputs = vec!["K1".to_string()];
        let parsed = parse_csv_response(&inputs, body);
        assert_eq!(parsed.len(), 1);
        assert!(parsed[0].1.is_empty());
    }
}
