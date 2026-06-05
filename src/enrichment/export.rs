//! CSV exporter for Stage 3 enrichment results.

use anyhow::{Context, Result};
use std::io::Write;

use crate::csv_fmt::fmt_csv_f64;
use crate::enrichment::types::EnrichmentResult;

/// Column header line. Order MUST match the spec.
const HEADER: &str = "EntryID,EntryName,Hits,Total,Expected,EnrichmentRatio,PValue,FDR,HitKeggIDs";

/// Export the enrichment result as CSV.
///
/// - `only_displayed = true`: skip rows with `displayed == false` (i.e.
///   those filtered out by the user's `min_hit_count` setting).
/// - `only_displayed = false`: write every row in `result.rows`.
///
/// `HitKeggIDs` is semicolon-separated; `hit_kegg_ids` is already sorted
/// alphabetically by `run_ora`. NaN cells render as empty.
pub fn export_csv<W: Write>(
    writer: &mut W,
    result: &EnrichmentResult,
    only_displayed: bool,
) -> Result<()> {
    export_csv_with_mode(writer, result, only_displayed, false, None)
}

/// Variant of `export_csv` that prepends a `# Mode: dual (POS+NEG)` comment
/// line when `is_dual = true`, and appends a `# MinGroupOverlap: N` line
/// (after `# MinEntrySize:`) when `min_group_overlap = Some(n)` — i.e. for
/// Module-mode runs, where `n` is the run's
/// `module_retention.min_group_overlap`. Pathway runs pass `None`, so the line
/// is omitted and the output stays bit-equal to the pre-change exporter.
/// Single-mode/Pathway callers via `export_csv` get bit-equal output.
pub fn export_csv_with_mode<W: Write>(
    writer: &mut W,
    result: &EnrichmentResult,
    only_displayed: bool,
    is_dual: bool,
    min_group_overlap: Option<usize>,
) -> Result<()> {
    if is_dual {
        writeln!(writer, "# Mode: dual (POS+NEG)").context("write Mode tag line")?;
    }
    // Leading FDR-method tag line. Downstream parsers can strip via
    // comment='#' (pandas) or equivalent; the column header row is unchanged.
    writeln!(writer, "# FDR: {}", result.fdr_method.short_label()).context("write FDR tag line")?;
    // Pre-FDR `min_entry_size` filter tag line. Added in v3 by
    // add-min-entry-size-filter so the CSV self-documents which entries
    // were dropped before FDR.
    writeln!(writer, "# MinEntrySize: {}", result.min_entry_size)
        .context("write MinEntrySize tag line")?;
    // Module-mode Group-overlap filter tag line (Module mode only — Pathway
    // passes `None`). Self-documents the `min_group_overlap` the run used.
    if let Some(overlap) = min_group_overlap {
        writeln!(writer, "# MinGroupOverlap: {overlap}")
            .context("write MinGroupOverlap tag line")?;
    }
    // Rows go through `csv::Writer` for RFC-4180-correct quoting (the prior
    // hand-rolled `csv_escape` missed bare `\r` and surrounding whitespace).
    // The csv crate's default `\n` terminator matches the previous `writeln!`,
    // so the output is byte-identical for names without quote-forcing
    // characters. The comment lines above were written raw first so they are
    // never quoted.
    let mut w = csv::Writer::from_writer(writer);
    w.write_record(HEADER.split(','))
        .context("write CSV header")?;
    for row in &result.rows {
        if only_displayed && !row.displayed {
            continue;
        }
        w.write_record([
            row.entry_id.clone(),
            row.entry_name.clone(),
            row.hits.to_string(),
            row.total.to_string(),
            fmt_csv_f64(row.expected),
            fmt_csv_f64(row.enrichment_ratio),
            fmt_csv_f64(row.p_value),
            fmt_csv_f64(row.fdr),
            row.hit_kegg_ids.join(";"),
        ])
        .context("write CSV row")?;
    }
    w.flush().context("flush CSV writer")?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enrichment::types::{EnrichmentDirection, EnrichmentResult, EnrichmentRow};

    fn sample_result() -> EnrichmentResult {
        EnrichmentResult {
            universe_size: 100,
            dam_cpd_size: 10,
            direction: EnrichmentDirection::Both,
            min_hit_count: 1,
            min_entry_size: 1,
            entries_dropped_by_min_entry_size: 0,
            empty_compound_count: 0,
            rows: vec![
                EnrichmentRow {
                    entry_id: "p1".into(),
                    entry_name: "Glycolysis, simple".into(), // has comma
                    hits: 3,
                    total: 8,
                    expected: 0.8,
                    enrichment_ratio: 3.75,
                    p_value: 0.01,
                    fdr: 0.02,
                    hit_kegg_ids: vec!["C00031".into(), "C00074".into(), "C00103".into()],
                    displayed: true,
                },
                EnrichmentRow {
                    entry_id: "p2".into(),
                    entry_name: "Short".into(),
                    hits: 0,
                    total: 5,
                    expected: 0.0,
                    enrichment_ratio: f64::NAN,
                    p_value: 1.0,
                    fdr: 1.0,
                    hit_kegg_ids: vec![],
                    displayed: false,
                },
            ],
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        }
    }

    #[test]
    fn header_is_fixed() {
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result(), false).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // First three lines: `# FDR: …`, `# MinEntrySize: …`, then header.
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], "# FDR: BY");
        assert_eq!(lines[1], "# MinEntrySize: 1");
        assert_eq!(lines[2], HEADER);
    }

    #[test]
    fn min_entry_size_tag_line_reflects_value() {
        // Default fixture has min_entry_size = 1.
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result(), false).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(
            s.contains("\n# MinEntrySize: 1\n"),
            "expected MinEntrySize=1 tag; got: {s}"
        );

        // Bump the result's min_entry_size and re-export.
        let mut r = sample_result();
        r.min_entry_size = 5;
        let mut buf2 = Vec::new();
        export_csv(&mut buf2, &r, false).unwrap();
        let s2 = String::from_utf8(buf2).unwrap();
        assert!(
            s2.contains("\n# MinEntrySize: 5\n"),
            "expected MinEntrySize=5 tag; got: {s2}"
        );
    }

    #[test]
    fn min_group_overlap_tag_line_module_vs_pathway() {
        // Module mode: `Some(3)` emits `# MinGroupOverlap: 3` immediately
        // after `# MinEntrySize:`, before the header row.
        let mut buf = Vec::new();
        export_csv_with_mode(&mut buf, &sample_result(), false, false, Some(3)).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], "# FDR: BY");
        assert_eq!(lines[1], "# MinEntrySize: 1");
        assert_eq!(lines[2], "# MinGroupOverlap: 3");
        assert_eq!(lines[3], HEADER);

        // Pathway mode: `None` emits no MinGroupOverlap line and is
        // byte-identical to the 3-arg `export_csv` wrapper.
        let mut buf_none = Vec::new();
        export_csv_with_mode(&mut buf_none, &sample_result(), false, false, None).unwrap();
        let mut buf_wrapper = Vec::new();
        export_csv(&mut buf_wrapper, &sample_result(), false).unwrap();
        assert_eq!(
            buf_none, buf_wrapper,
            "None must match the 3-arg export_csv wrapper byte-for-byte"
        );
        assert!(
            !String::from_utf8(buf_none)
                .unwrap()
                .contains("# MinGroupOverlap"),
            "Pathway export must not contain a MinGroupOverlap line"
        );

        // Dual-mode Pathway (`is_dual = true`, `None`): unchanged order
        // `# Mode:` / `# FDR:` / `# MinEntrySize:` / header, no overlap line.
        let mut buf_dual = Vec::new();
        export_csv_with_mode(&mut buf_dual, &sample_result(), false, true, None).unwrap();
        let dual = String::from_utf8(buf_dual).unwrap();
        let dlines: Vec<&str> = dual.lines().collect();
        assert_eq!(dlines[0], "# Mode: dual (POS+NEG)");
        assert_eq!(dlines[1], "# FDR: BY");
        assert_eq!(dlines[2], "# MinEntrySize: 1");
        assert_eq!(dlines[3], HEADER);
        assert!(!dual.contains("# MinGroupOverlap"));
    }

    #[test]
    fn fdr_tag_line_reflects_method() {
        // Default fixture is BY → leading `# FDR: BY`.
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result(), false).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(
            s.starts_with("# FDR: BY\n"),
            "expected BY tag; got first 20: {:?}",
            &s[..s.len().min(20)]
        );

        // Switch to BH and re-export.
        let mut r = sample_result();
        r.fdr_method = crate::dam::fdr::FdrMethod::BenjaminiHochberg;
        let mut buf2 = Vec::new();
        export_csv(&mut buf2, &r, false).unwrap();
        let s2 = String::from_utf8(buf2).unwrap();
        assert!(
            s2.starts_with("# FDR: BH\n"),
            "expected BH tag; got first 20: {:?}",
            &s2[..s2.len().min(20)]
        );
    }

    #[test]
    fn comment_aware_parser_recovers_same_row_count_regardless_of_method() {
        // Simulate `pandas.read_csv(..., comment='#')` by stripping all lines
        // starting with `#` and counting data rows. Both methods must yield
        // the same data shape.
        fn data_lines(s: &str) -> Vec<&str> {
            s.lines().filter(|l| !l.starts_with('#')).collect()
        }
        let mut r = sample_result();
        let mut buf_by = Vec::new();
        export_csv(&mut buf_by, &r, false).unwrap();
        let s_by = String::from_utf8(buf_by).unwrap();

        r.fdr_method = crate::dam::fdr::FdrMethod::BenjaminiHochberg;
        let mut buf_bh = Vec::new();
        export_csv(&mut buf_bh, &r, false).unwrap();
        let s_bh = String::from_utf8(buf_bh).unwrap();

        let by_lines = data_lines(&s_by);
        let bh_lines = data_lines(&s_bh);
        assert_eq!(
            by_lines.len(),
            bh_lines.len(),
            "row counts must match across methods"
        );
        // First non-comment line is the header in both.
        assert_eq!(by_lines[0], HEADER);
        assert_eq!(bh_lines[0], HEADER);
    }

    #[test]
    fn semicolon_joined_hit_ids() {
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result(), false).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("C00031;C00074;C00103"));
    }

    #[test]
    fn nan_renders_as_empty_cell() {
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result(), false).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // p2 row has NaN ratio; expect ",," for that cell (empty between commas).
        let p2_line = s.lines().find(|l| l.starts_with("p2,")).unwrap();
        // p2 row breakdown: p2,Short,0,5,0,,1,1,
        // (expected=0, ratio NaN→empty, p=1, fdr=1, ids empty)
        assert!(
            p2_line.contains(",,1,1,"),
            "expected empty NaN cell in: {p2_line}"
        );
    }

    #[test]
    fn name_with_comma_is_quoted() {
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result(), false).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("\"Glycolysis, simple\""));
    }

    #[test]
    fn only_displayed_skips_filtered_rows() {
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result(), true).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // p1 is displayed=true, p2 is displayed=false → only p1 line.
        assert!(s.contains("p1,"));
        assert!(!s.contains("p2,"));
    }

    #[test]
    fn name_with_carriage_return_is_quoted() {
        // The intended escaping fix: a bare `\r` in an entry name forces
        // RFC-4180 quoting (the old hand-rolled escaper missed `\r`).
        let mut r = sample_result();
        r.rows[0].entry_name = "Glyco\rlysis".into();
        let mut buf = Vec::new();
        export_csv(&mut buf, &r, false).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(
            s.contains("\"Glyco\rlysis\""),
            "a name with a bare CR must be RFC-4180-quoted; got: {s:?}"
        );
    }

    #[test]
    fn comment_lines_are_not_quoted() {
        // Single-mode: first line is the FDR tag, written raw (unquoted).
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result(), false).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let first = s.lines().next().unwrap();
        assert!(first.starts_with("# FDR: "), "got: {first}");
        assert!(!first.starts_with('"'), "comment line must not be quoted");

        // Dual-mode: first line is the Mode tag, also raw.
        let mut buf2 = Vec::new();
        export_csv_with_mode(&mut buf2, &sample_result(), false, true, None).unwrap();
        let s2 = String::from_utf8(buf2).unwrap();
        let first2 = s2.lines().next().unwrap();
        assert!(
            first2.starts_with("# Mode: dual (POS+NEG)"),
            "got: {first2}"
        );
        assert!(!first2.starts_with('"'), "comment line must not be quoted");
    }
}
