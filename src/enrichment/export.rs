//! CSV exporter for Stage 3 enrichment results.

use anyhow::{Context, Result};
use std::io::Write;

use crate::csv_fmt::fmt_csv_f64;
use crate::enrichment::types::EnrichmentResult;

/// Column header line for a CORRECTED run. Order MUST match the spec.
const HEADER: &str = "EntryID,EntryName,Hits,Total,Expected,EnrichmentRatio,PValue,FDR,HitKeggIDs";

/// Header for an uncorrected run — the `FDR` column is omitted.
const HEADER_NO_CORRECTION: &str =
    "EntryID,EntryName,Hits,Total,Expected,EnrichmentRatio,PValue,HitKeggIDs";

/// The column header for `method`.
///
/// Exactly one thing varies this header: the correction method. Under
/// `NoCorrection`, `adjust_pvalues` returns its input unchanged, so an `FDR`
/// column would repeat `PValue` byte-for-byte under a name no correction
/// earned — and two identically-valued columns read as two independent
/// measurements. Nothing else varies it: not the analysis mode, not single vs
/// dual input, not any display filter. Owner: the `enrichment-ora` capability.
fn header_for(method: crate::dam::fdr::FdrMethod) -> &'static str {
    match method {
        crate::dam::fdr::FdrMethod::NoCorrection => HEADER_NO_CORRECTION,
        _ => HEADER,
    }
}

/// Which rows an export writes.
///
/// Named rather than a bare `bool` so a call site cannot say "filtered" while
/// meaning something else. The Stage 3 result screen has two download buttons
/// and this is the only thing that differs between the files they write: the
/// comment block, the header and the column set are identical.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum RowSelection {
    /// Every surviving row — what `Download all results (CSV)` writes.
    All,
    /// The rows the FIGURE is drawn from — what
    /// `Download enrichment results (CSV)` writes.
    ///
    /// **`top_n` is absent on purpose, and it is the only exclusion.** It is an
    /// ordering cap rather than a per-row test: it answers how many of the
    /// ranked rows fit on an axis, so a file bounded by it would have a row
    /// count meaning "the twenty I happened to be looking at". Both fields here
    /// are tests each row passes or fails on its own values, which is what makes
    /// them meaningful in a file. So this selection is the plot's row set
    /// BEFORE truncation — a superset of what is drawn, by exactly the rows
    /// `top_n` cut.
    Figure {
        fdr_threshold: f64,
        min_hit_count: usize,
    },
}

impl RowSelection {
    /// Whether `row` belongs in this selection.
    fn admits(self, row: &crate::enrichment::types::EnrichmentRow) -> bool {
        match self {
            RowSelection::All => true,
            RowSelection::Figure {
                fdr_threshold,
                min_hit_count,
            } => row.hits >= min_hit_count && row.fdr < fdr_threshold,
        }
    }
}

/// Export every surviving row of the enrichment result as CSV.
///
/// `min_entry_size` bounds this file as it bounds the filtered one: it drops
/// entries from the result itself, before any p-value is computed, so no
/// consumer ever sees them.
///
/// `HitKeggIDs` is semicolon-separated; `hit_kegg_ids` is already sorted
/// alphabetically by `run_ora`. NaN cells render as empty.
pub fn export_csv<W: Write>(writer: &mut W, result: &EnrichmentResult) -> Result<()> {
    export_csv_with_mode(writer, result, false, None, RowSelection::All)
}

/// Variant of `export_csv` that prepends a `# Mode: dual (POS+NEG)` comment
/// line when `is_dual = true`, and appends a `# MinGroupOverlap: N` line
/// (after `# MinEntrySize:`) when `min_group_overlap = Some(n)` — i.e. for
/// Module-mode runs, where `n` is the run's
/// `module_retention.min_group_overlap`. Pathway runs pass `None`, so the line
/// is omitted and the output stays bit-equal to the pre-change exporter.
/// Single-mode/Pathway callers via `export_csv` get bit-equal output.
/// `selection` names which rows are written. It changes nothing else: both
/// selections emit the identical comment block, header and column set, so two
/// files from one run differ only in how many data rows they carry.
pub fn export_csv_with_mode<W: Write>(
    writer: &mut W,
    result: &EnrichmentResult,
    is_dual: bool,
    min_group_overlap: Option<usize>,
    selection: RowSelection,
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
    let uncorrected = matches!(result.fdr_method, crate::dam::fdr::FdrMethod::NoCorrection);
    let mut w = csv::Writer::from_writer(writer);
    w.write_record(header_for(result.fdr_method).split(','))
        .context("write CSV header")?;
    for row in result.rows.iter().filter(|r| selection.admits(r)) {
        let mut fields = vec![
            row.entry_id.clone(),
            row.entry_name.clone(),
            row.hits.to_string(),
            row.total.to_string(),
            fmt_csv_f64(row.expected),
            fmt_csv_f64(row.enrichment_ratio),
            fmt_csv_f64(row.p_value),
        ];
        if !uncorrected {
            fields.push(fmt_csv_f64(row.fdr));
        }
        fields.push(row.hit_kegg_ids.join(";"));
        w.write_record(&fields).context("write CSV row")?;
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
                },
            ],
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        }
    }

    #[test]
    fn header_is_fixed_for_a_corrected_run() {
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result()).unwrap();
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
        export_csv(&mut buf, &sample_result()).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(
            s.contains("\n# MinEntrySize: 1\n"),
            "expected MinEntrySize=1 tag; got: {s}"
        );

        // Bump the result's min_entry_size and re-export.
        let mut r = sample_result();
        r.min_entry_size = 5;
        let mut buf2 = Vec::new();
        export_csv(&mut buf2, &r).unwrap();
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
        export_csv_with_mode(
            &mut buf,
            &sample_result(),
            false,
            Some(3),
            RowSelection::All,
        )
        .unwrap();
        let s = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], "# FDR: BY");
        assert_eq!(lines[1], "# MinEntrySize: 1");
        assert_eq!(lines[2], "# MinGroupOverlap: 3");
        assert_eq!(lines[3], HEADER);

        // Pathway mode: `None` emits no MinGroupOverlap line and is
        // byte-identical to the 3-arg `export_csv` wrapper.
        let mut buf_none = Vec::new();
        export_csv_with_mode(
            &mut buf_none,
            &sample_result(),
            false,
            None,
            RowSelection::All,
        )
        .unwrap();
        let mut buf_wrapper = Vec::new();
        export_csv(&mut buf_wrapper, &sample_result()).unwrap();
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
        export_csv_with_mode(
            &mut buf_dual,
            &sample_result(),
            true,
            None,
            RowSelection::All,
        )
        .unwrap();
        let dual = String::from_utf8(buf_dual).unwrap();
        let dlines: Vec<&str> = dual.lines().collect();
        assert_eq!(dlines[0], "# Mode: dual (POS+NEG)");
        assert_eq!(dlines[1], "# FDR: BY");
        assert_eq!(dlines[2], "# MinEntrySize: 1");
        assert_eq!(dlines[3], HEADER);
        assert!(!dual.contains("# MinGroupOverlap"));
    }

    /// `RowSelection::All` applies NO display filter — not `min_hit_count`, not
    /// the significance threshold, not `top_n`.
    ///
    /// This is what `Download all results (CSV)` writes, the same contract the
    /// Stage 2 button of that name has. The screen had ONE button once, and it
    /// meant this and then meant the filtered set, each time silently; the fix
    /// was a second button rather than a better single meaning.
    #[test]
    fn the_all_results_selection_applies_no_display_filter() {
        let r = sample_result();
        let mut buf = Vec::new();
        export_csv(&mut buf, &r).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // The fixture's two rows have 3 and 0 hits. Both are written, including
        // the zero-hit row that every positive `min_hit_count` would have hidden.
        assert_eq!(data_rows(&s), r.rows.len());
        assert!(s.contains("p2"), "the zero-hit row must still be exported");
    }

    /// Data rows only — the comment block and the header stripped.
    fn data_rows(csv: &str) -> usize {
        csv.lines()
            .filter(|l| !l.starts_with('#') && !l.starts_with("EntryID"))
            .count()
    }

    /// A result whose rows straddle both per-row filters AND exceed any
    /// plausible `top_n`, so one fixture can prove all three points at once.
    fn wide_result() -> EnrichmentResult {
        let mut r = sample_result();
        r.rows = (0..12)
            .map(|i| EnrichmentRow {
                entry_id: format!("keep{i}"),
                entry_name: format!("Keeper {i}"),
                hits: 5,
                total: 20,
                expected: 1.0,
                enrichment_ratio: 5.0,
                p_value: 0.001,
                fdr: 0.001,
                hit_kegg_ids: vec!["C00031".into()],
            })
            // Fails `min_hit_count` but passes the threshold.
            .chain((0..4).map(|i| EnrichmentRow {
                entry_id: format!("thin{i}"),
                entry_name: format!("Too few hits {i}"),
                hits: 1,
                total: 20,
                expected: 1.0,
                enrichment_ratio: 1.0,
                p_value: 0.001,
                fdr: 0.001,
                hit_kegg_ids: vec![],
            }))
            // Passes `min_hit_count` but fails the threshold.
            .chain((0..7).map(|i| EnrichmentRow {
                entry_id: format!("dull{i}"),
                entry_name: format!("Not significant {i}"),
                hits: 9,
                total: 20,
                expected: 1.0,
                enrichment_ratio: 1.0,
                p_value: 0.9,
                fdr: 0.9,
                hit_kegg_ids: vec![],
            }))
            .collect();
        r
    }

    /// The pair a user takes from ONE run: identical comment block, header and
    /// column set; different row sets, and nothing else.
    #[test]
    fn the_two_selections_differ_in_rows_and_in_nothing_else() {
        let r = wide_result();
        let figure = RowSelection::Figure {
            fdr_threshold: 0.05,
            min_hit_count: 3,
        };

        let mut all_buf = Vec::new();
        export_csv_with_mode(&mut all_buf, &r, false, Some(3), RowSelection::All).unwrap();
        let all = String::from_utf8(all_buf).unwrap();

        let mut fig_buf = Vec::new();
        export_csv_with_mode(&mut fig_buf, &r, false, Some(3), figure).unwrap();
        let fig = String::from_utf8(fig_buf).unwrap();

        // Everything above the first data row is byte-identical.
        let head = |s: &str| {
            s.lines()
                .take_while(|l| l.starts_with('#') || l.starts_with("EntryID"))
                .collect::<Vec<_>>()
                .join("\n")
        };
        assert_eq!(head(&all), head(&fig));

        assert_eq!(data_rows(&all), 23, "every surviving row");
        assert_eq!(data_rows(&fig), 12, "both per-row filters applied");
        assert!(!fig.contains("thin0"), "min_hit_count excludes it");
        assert!(!fig.contains("dull0"), "the threshold excludes it");
        assert!(fig.contains("keep11"));
    }

    /// `top_n` reaches NEITHER file, and it is the only exclusion. It caps how
    /// many ranked rows fit on an axis; a file bounded by it would have a row
    /// count meaning "the twenty I happened to be looking at".
    #[test]
    fn the_display_cap_reaches_neither_file() {
        let r = wide_result();
        let mut fig_buf = Vec::new();
        export_csv_with_mode(
            &mut fig_buf,
            &r,
            false,
            None,
            RowSelection::Figure {
                fdr_threshold: 0.05,
                min_hit_count: 3,
            },
        )
        .unwrap();
        let fig = String::from_utf8(fig_buf).unwrap();
        // With `top_n = 5` the plot would draw 5 of these 12. The file has all
        // 12 — the figure's row set BEFORE truncation.
        let top_n_the_plot_would_use = 5;
        assert_eq!(data_rows(&fig), 12);
        assert!(
            data_rows(&fig) > top_n_the_plot_would_use,
            "the file must carry more rows than the axis would show, or the \
             point is untested"
        );
    }

    /// Reachable with a strict threshold, and a natural first bug: the file is
    /// written, not suppressed and not empty.
    #[test]
    fn a_filtered_export_with_no_passing_rows_is_a_header_not_an_empty_file() {
        let r = wide_result();
        let mut buf = Vec::new();
        export_csv_with_mode(
            &mut buf,
            &r,
            false,
            None,
            RowSelection::Figure {
                fdr_threshold: 1e-12,
                min_hit_count: 1,
            },
        )
        .unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert_eq!(data_rows(&s), 0);
        assert!(s.starts_with("# FDR: BY\n# MinEntrySize: 1\n"));
        assert!(s.contains(HEADER), "the column header is still written");
    }

    /// An uncorrected export drops the `FDR` column entirely — it would have
    /// repeated `PValue` byte-for-byte under a name no correction earned, and
    /// two identically-valued columns read as two independent measurements.
    #[test]
    fn uncorrected_export_omits_the_fdr_column() {
        let mut r = sample_result();
        r.fdr_method = crate::dam::fdr::FdrMethod::NoCorrection;
        // `adjust_pvalues` passes p through for this method; mirror that here so
        // the duplication the omission prevents is actually present in the data.
        for row in &mut r.rows {
            row.fdr = row.p_value;
        }
        let mut buf = Vec::new();
        export_csv(&mut buf, &r).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = s.lines().collect();

        assert_eq!(lines[0], "# FDR: NoCorrection");
        let header = lines.iter().find(|l| l.starts_with("EntryID")).unwrap();
        assert_eq!(
            *header,
            "EntryID,EntryName,Hits,Total,Expected,EnrichmentRatio,PValue,HitKeggIDs"
        );
        assert!(!header.contains("FDR"));

        // Every data row is eight fields too, and no p-value appears twice.
        // Parsed with the csv reader, not `split(',')` — entry names may contain
        // a quoted comma, which is exactly why the exporter uses `csv::Writer`.
        let body: String = lines
            .iter()
            .filter(|l| !l.starts_with('#'))
            .map(|l| format!("{l}\n"))
            .collect();
        let mut rdr = csv::Reader::from_reader(body.as_bytes());
        assert_eq!(rdr.headers().unwrap().len(), 8);
        for rec in rdr.records() {
            let rec = rec.unwrap();
            assert_eq!(rec.len(), 8, "row: {rec:?}");
            let p = &rec[6];
            assert_eq!(
                rec.iter().filter(|f| *f == p).count(),
                1,
                "p-value repeated under a second name: {rec:?}"
            );
        }
    }

    #[test]
    fn fdr_tag_line_reflects_method() {
        // Default fixture is BY → leading `# FDR: BY`.
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result()).unwrap();
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
        export_csv(&mut buf2, &r).unwrap();
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
        export_csv(&mut buf_by, &r).unwrap();
        let s_by = String::from_utf8(buf_by).unwrap();

        r.fdr_method = crate::dam::fdr::FdrMethod::BenjaminiHochberg;
        let mut buf_bh = Vec::new();
        export_csv(&mut buf_bh, &r).unwrap();
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
        export_csv(&mut buf, &sample_result()).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("C00031;C00074;C00103"));
    }

    #[test]
    fn nan_renders_as_empty_cell() {
        let mut buf = Vec::new();
        export_csv(&mut buf, &sample_result()).unwrap();
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
        export_csv(&mut buf, &sample_result()).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains("\"Glycolysis, simple\""));
    }

    #[test]
    fn name_with_carriage_return_is_quoted() {
        // The intended escaping fix: a bare `\r` in an entry name forces
        // RFC-4180 quoting (the old hand-rolled escaper missed `\r`).
        let mut r = sample_result();
        r.rows[0].entry_name = "Glyco\rlysis".into();
        let mut buf = Vec::new();
        export_csv(&mut buf, &r).unwrap();
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
        export_csv(&mut buf, &sample_result()).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let first = s.lines().next().unwrap();
        assert!(first.starts_with("# FDR: "), "got: {first}");
        assert!(!first.starts_with('"'), "comment line must not be quoted");

        // Dual-mode: first line is the Mode tag, also raw.
        let mut buf2 = Vec::new();
        export_csv_with_mode(&mut buf2, &sample_result(), true, None, RowSelection::All).unwrap();
        let s2 = String::from_utf8(buf2).unwrap();
        let first2 = s2.lines().next().unwrap();
        assert!(
            first2.starts_with("# Mode: dual (POS+NEG)"),
            "got: {first2}"
        );
        assert!(!first2.starts_with('"'), "comment line must not be quoted");
    }
}
