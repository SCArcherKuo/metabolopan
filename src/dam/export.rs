//! CSV exporters for DAM results. Two flavours:
//!
//! - `export_dam_csv` — only features whose trend (under the active thresholds) is Up
//!   or Down.
//! - `export_all_csv` — every feature in `DamResult.features`, regardless of trend.
//!
//! ±∞ numeric values are written as `inf` / `-inf`. NaN values are written as empty
//! cells (more Excel- and R-friendly than the literal string `NaN`).

use anyhow::{Context, Result};
use std::io::Write;

use crate::dam::run::classify_trend;
use crate::dam::types::{DamMethod, DamResult, FcBasis, Trend};
use crate::data::IonModeTable;
use crate::dedup::AdductClass;
use crate::dedup::types::{CascadeValue, DedupReport};

pub fn export_dam_csv(
    writer: impl Write,
    result: &DamResult,
    fc_threshold: f64,
    fdr_threshold: f64,
    delta_threshold: f64,
) -> Result<()> {
    write_csv(
        writer,
        result,
        fc_threshold,
        fdr_threshold,
        delta_threshold,
        true,
    )
}

pub fn export_all_csv(
    writer: impl Write,
    result: &DamResult,
    fc_threshold: f64,
    fdr_threshold: f64,
    delta_threshold: f64,
) -> Result<()> {
    write_csv(
        writer,
        result,
        fc_threshold,
        fdr_threshold,
        delta_threshold,
        false,
    )
}

/// Dual-mode variant of `export_dam_csv`. Emits a `# Mode: dual (POS+NEG)`
/// comment line ahead of the `# FDR:` line, prepends a `Mode` column to the
/// header, and writes rows mode-by-mode in `ion_tables` order (slot 0 first,
/// slot 1 second). Trend filtering matches the single-mode `export_dam_csv`.
pub fn export_dam_csv_multi(
    writer: impl Write,
    ion_tables: &[IonModeTable],
    results: &[DamResult],
    fc_threshold: f64,
    fdr_threshold: f64,
    delta_threshold: f64,
) -> Result<()> {
    write_csv_multi(
        writer,
        ion_tables,
        results,
        fc_threshold,
        fdr_threshold,
        delta_threshold,
        true,
    )
}

/// Dual-mode variant of `export_all_csv`. Emits a `# Mode: dual (POS+NEG)`
/// header, the leading `Mode` column, and every feature row per mode in
/// `ion_tables` order regardless of trend.
pub fn export_all_csv_multi(
    writer: impl Write,
    ion_tables: &[IonModeTable],
    results: &[DamResult],
    fc_threshold: f64,
    fdr_threshold: f64,
    delta_threshold: f64,
) -> Result<()> {
    write_csv_multi(
        writer,
        ion_tables,
        results,
        fc_threshold,
        fdr_threshold,
        delta_threshold,
        false,
    )
}

fn write_csv(
    mut writer: impl Write,
    result: &DamResult,
    fc_threshold: f64,
    fdr_threshold: f64,
    delta_threshold: f64,
    only_dam: bool,
) -> Result<()> {
    // Leading FDR-method tag line. Downstream parsers can strip via comment='#'
    // (pandas) or equivalent; the column header row is unchanged.
    writeln!(writer, "# FDR: {}", result.fdr_method.short_label())
        .context("failed to write FDR tag line")?;
    let mut w = csv::Writer::from_writer(writer);
    w.write_record([
        "AlignmentID",
        "MetaboliteName",
        "INCHIKEY",
        "Formula",
        "SMILES",
        "RT(min)",
        "Mz",
        "NumeratorMean",
        "DenominatorMean",
        "NumeratorMedian",
        "DenominatorMedian",
        "FoldChange",
        "Log2FoldChange",
        "FCBasis",
        "PValue",
        "FDR",
        "NegLog10FDR",
        "EffectSize",
        "Trend",
    ])
    .context("failed to write CSV header")?;

    for feat in &result.features {
        let trend = classify_trend(
            feat,
            fc_threshold,
            fdr_threshold,
            delta_threshold,
            result.method,
        );
        if only_dam && trend == Trend::NotSignificant {
            continue;
        }
        w.write_record([
            feat.alignment_id.as_str(),
            feat.metabolite_name.as_str(),
            feat.inchikey.as_deref().unwrap_or(""),
            feat.formula.as_deref().unwrap_or(""),
            feat.smiles.as_deref().unwrap_or(""),
            &fmt_opt(feat.average_rt_min),
            &fmt_opt(feat.average_mz),
            &fmt_f(feat.numerator_mean),
            &fmt_f(feat.denominator_mean),
            &fmt_f(feat.numerator_median),
            &fmt_f(feat.denominator_median),
            &fmt_f(feat.fold_change),
            &fmt_f(feat.log2_fold_change),
            fc_basis_label(feat.fc_basis),
            &fmt_f(feat.p_value),
            &fmt_f(feat.p_adjusted),
            &fmt_f(feat.neg_log10_p_adjusted),
            &fmt_opt_inner(feat.effect_size, result.method),
            trend.label(),
        ])
        .context("failed to write CSV row")?;
    }
    w.flush()?;
    Ok(())
}

fn write_csv_multi(
    mut writer: impl Write,
    ion_tables: &[IonModeTable],
    results: &[DamResult],
    fc_threshold: f64,
    fdr_threshold: f64,
    delta_threshold: f64,
    only_dam: bool,
) -> Result<()> {
    // Leading dual-mode comment line, then the shared FDR tag. We assume
    // every DamResult in `results` was produced under the same FDR method
    // (the Stage 2 setup screen uses a single radio for both modes); we
    // emit the first result's method for the tag.
    writeln!(writer, "# Mode: dual (POS+NEG)").context("failed to write Mode tag line")?;
    let fdr_label = results
        .first()
        .map(|r| r.fdr_method.short_label())
        .unwrap_or("BH");
    writeln!(writer, "# FDR: {fdr_label}").context("failed to write FDR tag line")?;
    let mut w = csv::Writer::from_writer(writer);
    w.write_record([
        "Mode",
        "AlignmentID",
        "MetaboliteName",
        "INCHIKEY",
        "Formula",
        "SMILES",
        "RT(min)",
        "Mz",
        "NumeratorMean",
        "DenominatorMean",
        "NumeratorMedian",
        "DenominatorMedian",
        "FoldChange",
        "Log2FoldChange",
        "FCBasis",
        "PValue",
        "FDR",
        "NegLog10FDR",
        "EffectSize",
        "Trend",
    ])
    .context("failed to write CSV header")?;
    for (idx, result) in results.iter().enumerate() {
        let mode_label = ion_tables
            .get(idx)
            .map(|it| it.mode.to_string())
            .unwrap_or_else(|| format!("mode{idx}"));
        for feat in &result.features {
            let trend = classify_trend(
                feat,
                fc_threshold,
                fdr_threshold,
                delta_threshold,
                result.method,
            );
            if only_dam && trend == Trend::NotSignificant {
                continue;
            }
            w.write_record([
                mode_label.as_str(),
                feat.alignment_id.as_str(),
                feat.metabolite_name.as_str(),
                feat.inchikey.as_deref().unwrap_or(""),
                feat.formula.as_deref().unwrap_or(""),
                feat.smiles.as_deref().unwrap_or(""),
                &fmt_opt(feat.average_rt_min),
                &fmt_opt(feat.average_mz),
                &fmt_f(feat.numerator_mean),
                &fmt_f(feat.denominator_mean),
                &fmt_f(feat.numerator_median),
                &fmt_f(feat.denominator_median),
                &fmt_f(feat.fold_change),
                &fmt_f(feat.log2_fold_change),
                fc_basis_label(feat.fc_basis),
                &fmt_f(feat.p_value),
                &fmt_f(feat.p_adjusted),
                &fmt_f(feat.neg_log10_p_adjusted),
                &fmt_opt_inner(feat.effect_size, result.method),
                trend.label(),
            ])
            .context("failed to write CSV row")?;
        }
    }
    w.flush()?;
    Ok(())
}

fn fmt_f(v: f64) -> String {
    crate::csv_fmt::fmt_csv_f64(v)
}

fn fmt_opt(v: Option<f64>) -> String {
    match v {
        None => String::new(),
        Some(x) => fmt_f(x),
    }
}

fn fmt_opt_inner(v: Option<f64>, method: DamMethod) -> String {
    match (v, method) {
        // Parametric tests (Welch and Student) never produce an effect size; emit empty.
        (_, DamMethod::Welch | DamMethod::Student) => String::new(),
        (None, DamMethod::BrunnerMunzel) => String::new(),
        (Some(x), DamMethod::BrunnerMunzel) => fmt_f(x),
    }
}

fn fc_basis_label(b: FcBasis) -> &'static str {
    b.label()
}

/// Write the deduplication audit CSV. Two leading `#` comment lines
/// (banner + counts) sit above the column header; one data row per
/// dropped feature follows. Downstream parsers can strip the comment
/// lines via `comment='#'` (pandas) or equivalent.
///
/// Format defined by the `msdial-deduplication` capability spec.
pub fn export_dedup_audit_csv(writer: &mut impl Write, report: &DedupReport) -> Result<()> {
    writeln!(writer, "# Deduplication audit — generated by metabolopan")
        .context("failed to write dedup audit banner")?;
    writeln!(
        writer,
        "# Total dropped: {}; total kept: {}; null-InChIKey passthrough: {}",
        report.dropped.len(),
        report.kept_count,
        report.null_inchikey_passthrough
    )
    .context("failed to write dedup audit counts")?;
    let mut w = csv::Writer::from_writer(writer);
    w.write_record([
        "dropped_alignment_id",
        "inchikey",
        "winner_alignment_id",
        "decided_at",
        "loser_value",
        "winner_value",
    ])
    .context("failed to write dedup audit header")?;
    for d in &report.dropped {
        w.write_record([
            d.alignment_id.as_str(),
            d.inchikey.as_str(),
            d.winner_alignment_id.as_str(),
            d.decided_at.label(),
            &format_cascade_value(d.loser_value.as_ref()),
            &format_cascade_value(d.winner_value.as_ref()),
        ])
        .context("failed to write dedup audit row")?;
    }
    w.flush()?;
    Ok(())
}

fn format_cascade_value(v: Option<&CascadeValue>) -> String {
    match v {
        None => String::new(),
        Some(CascadeValue::Num(x)) => fmt_f(*x),
        Some(CascadeValue::Adduct { class, sub }) => match class {
            AdductClass::Primary => format!("Primary:{sub}"),
            other => other.label().to_string(),
        },
        Some(CascadeValue::Msms(true)) => "True".to_string(),
        Some(CascadeValue::Msms(false)) => "False".to_string(),
        Some(CascadeValue::AlignmentId(s)) => s.clone(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dam::types::DamFeature;
    use crate::dedup::types::CascadeStep;

    fn synth(name: &str, p_adj: f64, log2_fc: f64) -> DamFeature {
        DamFeature {
            alignment_id: name.into(),
            metabolite_name: name.into(),
            inchikey: Some("AAAA-BBBB".into()),
            average_rt_min: Some(1.0),
            average_mz: Some(100.0),
            formula: Some("C2H4O".into()),
            smiles: Some("CCO".into()),
            numerator_mean: 10.0,
            denominator_mean: 1.0,
            numerator_median: 10.0,
            denominator_median: 1.0,
            fold_change: 10.0,
            log2_fold_change: log2_fc,
            fc_basis: FcBasis::Mean,
            p_value: 0.01,
            p_adjusted: p_adj,
            neg_log10_p_adjusted: -p_adj.log10(),
            effect_size: None,
        }
    }

    #[test]
    fn dam_only_excludes_ns_rows() {
        let result = DamResult {
            method: DamMethod::Welch,
            numerator: "A".into(),
            denominator: "B".into(),
            features: vec![
                synth("up", 0.01, 2.0),
                synth("ns", 0.5, 2.0),   // not significant: FDR too large
                synth("ns2", 0.01, 0.1), // not significant: FC too small
            ],
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            dedup_report: None,
        };
        let mut buf = Vec::<u8>::new();
        export_dam_csv(&mut buf, &result, 2.0, 0.05, 0.33).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let lines: Vec<&str> = s.lines().collect();
        assert_eq!(lines[0], "# FDR: BY", "leading FDR tag line");
        assert_eq!(lines.len(), 3, "FDR tag + header + 1 row only; got {s}");
        assert!(lines[2].starts_with("up,"), "got {}", lines[2]);
    }

    #[test]
    fn fdr_tag_line_reflects_method() {
        // BH-tagged result emits `# FDR: BH` as the first line.
        let mut result = DamResult {
            method: DamMethod::Welch,
            numerator: "A".into(),
            denominator: "B".into(),
            features: vec![synth("up", 0.01, 2.0)],
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiHochberg,
            dedup_report: None,
        };
        let mut buf = Vec::<u8>::new();
        export_dam_csv(&mut buf, &result, 2.0, 0.05, 0.33).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(
            s.starts_with("# FDR: BH\n"),
            "BH tag; got first 20 chars: {:?}",
            &s[..s.len().min(20)]
        );

        // Flip to BY and re-export — first line changes to `# FDR: BY`.
        result.fdr_method = crate::dam::fdr::FdrMethod::BenjaminiYekutieli;
        let mut buf2 = Vec::<u8>::new();
        export_all_csv(&mut buf2, &result, 2.0, 0.05, 0.33).unwrap();
        let s2 = String::from_utf8(buf2).unwrap();
        assert!(
            s2.starts_with("# FDR: BY\n"),
            "BY tag; got first 20 chars: {:?}",
            &s2[..s2.len().min(20)]
        );
    }

    #[test]
    fn inf_serialises_as_inf() {
        let mut feat = synth("x", 0.01, f64::INFINITY);
        feat.fold_change = f64::INFINITY;
        let result = DamResult {
            method: DamMethod::Welch,
            numerator: "A".into(),
            denominator: "B".into(),
            features: vec![feat],
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            dedup_report: None,
        };
        let mut buf = Vec::<u8>::new();
        export_all_csv(&mut buf, &result, 2.0, 0.05, 0.33).unwrap();
        let s = String::from_utf8(buf).unwrap();
        assert!(s.contains(",inf,"), "expected inf token; got:\n{s}");
    }

    #[test]
    fn nan_serialises_as_empty() {
        let mut feat = synth("x", f64::NAN, 1.5);
        feat.p_value = f64::NAN;
        feat.neg_log10_p_adjusted = f64::NAN;
        let result = DamResult {
            method: DamMethod::Welch,
            numerator: "A".into(),
            denominator: "B".into(),
            features: vec![feat],
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            dedup_report: None,
        };
        let mut buf = Vec::<u8>::new();
        export_all_csv(&mut buf, &result, 2.0, 0.05, 0.33).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // Expect "" between commas for the NaN columns.
        assert!(s.contains(",,"), "expected empty cell for NaN: {s}");
    }

    // ── export_dedup_audit_csv tests ──

    use crate::dedup::types::DroppedFeature;

    fn drop(
        alignment_id: &str,
        inchikey: &str,
        winner: &str,
        decided_at: CascadeStep,
        loser: Option<CascadeValue>,
        winner_val: Option<CascadeValue>,
    ) -> DroppedFeature {
        DroppedFeature {
            alignment_id: alignment_id.into(),
            inchikey: inchikey.into(),
            winner_alignment_id: winner.into(),
            decided_at,
            loser_value: loser,
            winner_value: winner_val,
        }
    }

    #[test]
    fn dedup_audit_csv_header_and_shape() {
        let report = DedupReport {
            dropped: vec![
                drop(
                    "PEAK_0123",
                    "BQJC...",
                    "PEAK_0098",
                    CascadeStep::TotalScore,
                    Some(CascadeValue::Num(712.5)),
                    Some(CascadeValue::Num(891.2)),
                ),
                drop(
                    "PEAK_0456",
                    "XLY...",
                    "PEAK_0455",
                    CascadeStep::AdductClass,
                    Some(CascadeValue::Adduct {
                        class: AdductClass::NonPrimary,
                        sub: 0,
                    }),
                    Some(CascadeValue::Adduct {
                        class: AdductClass::Primary,
                        sub: 0,
                    }),
                ),
            ],
            kept_count: 50,
            null_inchikey_passthrough: 10,
        };
        let mut buf = Vec::<u8>::new();
        export_dedup_audit_csv(&mut buf, &report).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // 5 total lines: 2 comments + 1 header + 2 data
        assert_eq!(s.lines().count(), 5);
        let mut lines = s.lines();
        assert!(lines.next().unwrap().starts_with("# Deduplication audit"));
        assert_eq!(
            lines.next().unwrap(),
            "# Total dropped: 2; total kept: 50; null-InChIKey passthrough: 10"
        );
        assert_eq!(
            lines.next().unwrap(),
            "dropped_alignment_id,inchikey,winner_alignment_id,decided_at,loser_value,winner_value"
        );
    }

    #[test]
    fn dedup_audit_csv_numeric_rendering() {
        let report = DedupReport {
            dropped: vec![drop(
                "L",
                "K",
                "W",
                CascadeStep::TotalScore,
                Some(CascadeValue::Num(712.5)),
                Some(CascadeValue::Num(891.2)),
            )],
            kept_count: 1,
            null_inchikey_passthrough: 0,
        };
        let mut buf = Vec::<u8>::new();
        export_dedup_audit_csv(&mut buf, &report).unwrap();
        let s = String::from_utf8(buf).unwrap();
        // 4th line is the only data row: cells split by commas.
        let row = s.lines().nth(3).unwrap();
        let cells: Vec<&str> = row.split(',').collect();
        assert_eq!(cells[3], "TotalScore");
        assert_eq!(cells[4], "712.5");
        assert_eq!(cells[5], "891.2");
    }

    #[test]
    fn dedup_audit_csv_adduct_class_with_sub_rank() {
        let report = DedupReport {
            dropped: vec![drop(
                "L",
                "K",
                "W",
                CascadeStep::AdductClass,
                Some(CascadeValue::Adduct {
                    class: AdductClass::NonPrimary,
                    sub: 0,
                }),
                Some(CascadeValue::Adduct {
                    class: AdductClass::Primary,
                    sub: 0,
                }),
            )],
            kept_count: 1,
            null_inchikey_passthrough: 0,
        };
        let mut buf = Vec::<u8>::new();
        export_dedup_audit_csv(&mut buf, &report).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let row = s.lines().nth(3).unwrap();
        let cells: Vec<&str> = row.split(',').collect();
        assert_eq!(cells[3], "AdductClass");
        // NonPrimary has no sub-rank suffix; Primary always carries one.
        assert_eq!(cells[4], "NonPrimary");
        assert_eq!(cells[5], "Primary:0");
    }

    #[test]
    fn dedup_audit_csv_nan_renders_empty() {
        let report = DedupReport {
            dropped: vec![drop(
                "L",
                "K",
                "W",
                CascadeStep::TotalScore,
                Some(CascadeValue::Num(f64::NAN)),
                Some(CascadeValue::Num(f64::NAN)),
            )],
            kept_count: 1,
            null_inchikey_passthrough: 0,
        };
        let mut buf = Vec::<u8>::new();
        export_dedup_audit_csv(&mut buf, &report).unwrap();
        let s = String::from_utf8(buf).unwrap();
        let row = s.lines().nth(3).unwrap();
        // Cells 4 and 5 should be empty (no "NaN" literal anywhere)
        assert!(
            !row.contains("NaN"),
            "expected no NaN literal in audit CSV: {row}"
        );
        assert!(row.ends_with(",,"), "expected trailing empty cells: {row}");
    }
}
