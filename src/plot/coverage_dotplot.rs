//! Coverage dot plot renderer.
//!
//! Same chrome as the volcano and enrichment plots — white background, black
//! axes, no minor gridlines, scale-everything baseline, pHYs DPI metadata — on
//! the shared `crate::plot::common` scaffold.
//!
//! **The marker encoding deliberately inverts the enrichment plot's.** There it
//! is size = hits, colour = FDR; here it is size = `entry_size`, colour =
//! `hits`, with `coverage` on X. And there is no reference line at any X value:
//! the enrichment plot marks `enrichment_ratio = 1.0` because that is the null
//! expectation, and this route has no null to mark.
//!
//! theme: scope-excluded. YlOrRd = biological convention, not brand palette.
//!
//! Owner: the `coverage-dot-plot` capability.

use anyhow::{Result, anyhow};
use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};
use std::path::Path;

use crate::coverage::{
    CoverageResult, CoverageRow, CoverageSortKey, DisplayFilters, displayed_rows,
};
use crate::plot::common;

/// Options for [`render_coverage_dotplot`] and [`export_coverage_dotplot_png`].
///
/// **DPI is deliberately not a field.** It follows the established pattern: the
/// export function takes `dpi: u32` separately and the pixel data does not
/// depend on it — only the `pHYs` metadata changes. Putting DPI here would
/// imply it affects rendering.
#[derive(Debug, Clone)]
pub struct CoverageDotplotOpts {
    pub width_px: u32,
    pub height_px: u32,
    /// The display-filter chain, passed WHOLE rather than as separate `top_n` /
    /// `min_hit_count` / `min_entry_size` / `sort_key` fields. The renderer
    /// obtains its rows by calling `coverage::displayed_rows(result, filters)`,
    /// so "the rows drawn always equal the rows in the table" is structural
    /// rather than a convention two call sites must remember.
    pub filters: DisplayFilters,
    /// `"Pathway"` / `"Module"`.
    pub mode_label: String,
    /// Species code, or `"<Level> / <Group>"`.
    pub target_label: String,
    /// `|D|`.
    pub detected_total: usize,
    /// `Some((selected, total, threshold))` when a metadata `.csv` was
    /// supplied; `None` when it was not.
    pub group_record: Option<(usize, usize, f64)>,
}

/// Which quantity the X axis and the colour ramp carry.
///
/// The two channels always swap TOGETHER, driven by `filters.sort_key`: the
/// axis you sorted by is the axis you read the ordering off, so putting the
/// sort key anywhere but X leaves a chart whose row order looks arbitrary.
/// Marker size stays `entry_size` on both — it is the denominator behind both
/// quantities, and moving it would leave the reader with no fixed reference.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Encoding {
    /// X = coverage %, colour = hits. The default.
    CoverageOnX,
    /// X = hits, colour = coverage %.
    HitsOnX,
}

impl Encoding {
    fn from_sort_key(key: CoverageSortKey) -> Self {
        match key {
            CoverageSortKey::Hits => Self::HitsOnX,
            // `EntrySize` / `EntryId` are not offered by the UI today; either
            // one leaves the default encoding rather than inventing a third.
            _ => Self::CoverageOnX,
        }
    }

    fn x_desc(self) -> &'static str {
        match self {
            Self::CoverageOnX => "Coverage (% of entry's compounds detected)",
            Self::HitsOnX => "Hits (compounds detected in this entry)",
        }
    }

    fn colorbar_title(self) -> &'static str {
        match self {
            Self::CoverageOnX => "Hits",
            Self::HitsOnX => "Coverage %",
        }
    }

    /// The X value for one row, in the axis's own units.
    fn x_of(self, row: &CoverageRow) -> f64 {
        match self {
            Self::CoverageOnX => row.coverage * 100.0,
            Self::HitsOnX => row.hits as f64,
        }
    }

    /// The colour value for one row, in the colorbar's own units.
    fn color_of(self, row: &CoverageRow) -> f64 {
        match self {
            Self::CoverageOnX => row.hits as f64,
            Self::HitsOnX => row.coverage * 100.0,
        }
    }
}

/// Y-axis label budget per wrapped line, and the line cap. Same values and the
/// same shared `wrap_label` the enrichment plot uses, so entry names break
/// identically on both charts.
const CHARS_PER_LINE: usize = 42;
const MAX_LABEL_LINES: usize = 4;
const LINE_HEIGHT_PX_BASE: f64 = 16.0;
const LABEL_FONT: f64 = 14.0;
const TITLE_FONT: f64 = 20.0;

/// Base marker radius at scale 1.0, shared by the chart dots and the legend's
/// reference dots so a given `entry_size` renders at the SAME size in both.
const DOT_BASE_RADIUS: f64 = 10.0;

/// Shown when no row survives the filters. A valid image, never an error: an
/// over-tight filter is a thing the user did, not a failure of the renderer.
const EMPTY_TEXT: &str = "No entries match the current filters.";

/// The strip's no-statistics note. Pinned by the `coverage-dot-plot` spec.
const NO_STATS_NOTE: &str = "Descriptive coverage — no statistical test";

/// Render the coverage dot plot into an RGBA buffer of `width * height * 4`
/// bytes. The SAME function backs the in-window preview and the PNG export.
pub fn render_coverage_dotplot(
    result: &CoverageResult,
    opts: &CoverageDotplotOpts,
) -> Result<Vec<u8>> {
    common::ensure_font_registered();
    let w = opts.width_px;
    let h = opts.height_px;
    let pixel_count = (w as usize) * (h as usize);
    let mut rgb = vec![0u8; pixel_count * 3];

    // Width-based scale, matching the enrichment plot: the height auto-sizes to
    // the row count, so keying off `min(w, h)` would shrink every font on a
    // sparse result.
    let scale = common::design_scale_by_width(w);
    let legend_w = ((105.0 * scale).round() as u32).min(w * 3 / 10).max(120);
    let chart_w_boundary = (w - legend_w) as i32;

    // THE shared chain — not a re-implementation of it.
    let rows = displayed_rows(result, opts.filters);
    let encoding = Encoding::from_sort_key(opts.filters.sort_key);

    {
        let root = BitMapBackend::with_buffer(&mut rgb, (w, h)).into_drawing_area();
        root.fill(&WHITE)
            .map_err(|e| anyhow!("fill background: {e:?}"))?;

        let (chart_area, legend_area) = root.split_horizontally(chart_w_boundary);

        if rows.is_empty() {
            draw_empty_placeholder(&chart_area, scale)?;
        } else {
            draw_chart(&chart_area, &rows, encoding, scale)?;
        }
        draw_legend(&legend_area, &rows, encoding, scale)?;
        draw_annotation_strip(&root, result, opts, rows.len(), scale)?;

        root.present().map_err(|e| anyhow!("present: {e:?}"))?;
    }

    Ok(common::rgb_to_rgba(&rgb, pixel_count))
}

/// Render and save as PNG, embedding the requested DPI in the pHYs chunk via
/// the same shared encoder the volcano and enrichment plots use.
pub fn export_coverage_dotplot_png(
    result: &CoverageResult,
    opts: &CoverageDotplotOpts,
    dpi: u32,
    out: &Path,
) -> Result<()> {
    let buffer = render_coverage_dotplot(result, opts)?;
    common::encode_png(&buffer, opts.width_px, opts.height_px, dpi, out)
}

fn draw_chart<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    rows: &[&CoverageRow],
    encoding: Encoding,
    scale: f64,
) -> Result<()>
where
    <DB as DrawingBackend>::ErrorType: 'static,
{
    let sp = |v: f64| common::sp(v, scale);
    let su = |v: f64| common::su(v, scale);

    let (_, ah_h) = area.dim_in_pixel();
    let strip_h = su(110.0) as i32;
    let split_y = (ah_h as i32 - strip_h).max(1);
    let (chart_main, _strip_area) = area.split_vertically(split_y);

    let max_x = rows
        .iter()
        .map(|r| encoding.x_of(r))
        .fold(0.0_f64, f64::max);
    let x_max = x_bound(encoding, max_x);

    let max_entry_size = rows.iter().map(|r| r.entry_size).max().unwrap_or(1) as f64;
    let (c_lo, c_hi) = color_span(rows, encoding);

    // Row 0 at the TOP, matching the table's order. Cartesian Y grows upward,
    // so row i sits at y = (n - 1 - i) + 0.5.
    let n = rows.len();
    let row_center = |i: usize| -> f64 { (n - 1 - i) as f64 + 0.5 };

    let mut chart = ChartBuilder::on(&chart_main)
        .margin(sp(5.0))
        .x_label_area_size(su(80.0))
        .y_label_area_size(su(286.0))
        .build_cartesian_2d(0.0..x_max, 0f64..(n as f64))
        .map_err(|e| anyhow!("chart build: {e:?}"))?;

    chart
        .configure_mesh()
        .x_desc(encoding.x_desc())
        .y_desc("KEGG entry")
        .label_style(("sans-serif", su(LABEL_FONT)).into_font().color(&BLACK))
        .axis_desc_style(("sans-serif", su(TITLE_FONT)).into_font().color(&BLACK))
        .x_label_formatter(&|v: &f64| match encoding {
            Encoding::CoverageOnX => format!("{v:.0} %"),
            Encoding::HitsOnX => format!("{v:.0}"),
        })
        .disable_y_mesh()
        .y_labels(1)
        .y_label_formatter(&|_| String::new())
        .light_line_style(RGBAColor(0, 0, 0, 0.0))
        .bold_line_style(common::GRID_LIGHT)
        .axis_style(BLACK)
        .draw()
        .map_err(|e| anyhow!("mesh draw: {e:?}"))?;

    // Hand-drawn horizontal gridlines through each row centre (a continuous-axis
    // mesh would place them BETWEEN rows).
    {
        let (px_x, _) = chart.plotting_area().get_pixel_range();
        for i in 0..n {
            let (_, gy) = chart.backend_coord(&(0.0_f64, row_center(i)));
            chart_main
                .draw(&PathElement::new(
                    vec![(px_x.start, gy), (px_x.end, gy)],
                    ShapeStyle::from(&common::GRID_LIGHT).stroke_width(1),
                ))
                .map_err(|e| anyhow!("y gridline: {e:?}"))?;
        }
    }

    // Hand-drawn Y labels with the shared word-wrap.
    let (row0_x, row0_y) = chart.backend_coord(&(0.0_f64, row_center(0)));
    let row_step_px = if n >= 2 {
        let (_, row1_y) = chart.backend_coord(&(0.0_f64, row_center(1)));
        (row1_y - row0_y).unsigned_abs() as f64
    } else {
        let (_, py_range) = chart.plotting_area().get_pixel_range();
        (py_range.end - py_range.start).max(1) as f64
    };
    let line_height_px = sp(LINE_HEIGHT_PX_BASE);
    let max_lines = {
        let lh = (line_height_px as f64).max(1.0);
        ((row_step_px / lh).floor() as usize).clamp(1, MAX_LABEL_LINES)
    };
    let label_right_x = row0_x - sp(8.0);
    let label_style = ("sans-serif", su(LABEL_FONT))
        .into_font()
        .color(&BLACK)
        .pos(Pos::new(HPos::Right, VPos::Center));
    for (i, row) in rows.iter().enumerate() {
        let (_, cy) = chart.backend_coord(&(0.0_f64, row_center(i)));
        let lines = common::wrap_label(&row.entry_name, CHARS_PER_LINE, max_lines);
        let n_lines = lines.len() as f64;
        for (li, line) in lines.iter().enumerate() {
            let dy = ((li as f64 - (n_lines - 1.0) / 2.0) * line_height_px as f64).round() as i32;
            chart_main
                .draw_text(line, &label_style, (label_right_x, cy + dy))
                .map_err(|e| anyhow!("y label draw: {e:?}"))?;
        }
    }

    // NO reference line: there is no null expectation on this route.

    // Dots: radius from `entry_size` on both encodings; X and colour swap.
    let base_radius = (DOT_BASE_RADIUS * scale).max(2.0);
    for (i, row) in rows.iter().enumerate() {
        let r = (base_radius * (row.entry_size as f64 / max_entry_size.max(1.0)).sqrt()).max(2.0)
            as i32;
        let color = common::t_to_color(normalise(encoding.color_of(row), c_lo, c_hi));
        chart
            .draw_series(std::iter::once(Circle::new(
                (encoding.x_of(row), row_center(i)),
                r,
                color.filled(),
            )))
            .map_err(|e| anyhow!("dot draw: {e:?}"))?;
    }

    Ok(())
}

/// The X upper bound: the smallest round value at or above the data max.
///
/// Percent is capped at 100 — the quantity cannot exceed it. Hit counts have no
/// such ceiling, so they get a decade-scaled ladder instead of a fixed one.
fn x_bound(encoding: Encoding, max_x: f64) -> f64 {
    let m = if max_x.is_finite() {
        max_x.max(0.0)
    } else {
        0.0
    };
    match encoding {
        Encoding::CoverageOnX => {
            for bound in [5.0, 10.0, 20.0, 25.0, 50.0, 75.0, 100.0] {
                if m <= bound {
                    return bound;
                }
            }
            100.0
        }
        Encoding::HitsOnX => {
            // Round up to 1 / 2 / 5 × the enclosing power of ten, so the axis
            // ends on a value the eye reads as round at any magnitude.
            if m <= 5.0 {
                return 5.0;
            }
            let decade = 10f64.powf(m.log10().floor());
            for mult in [1.0, 2.0, 5.0, 10.0] {
                let bound = mult * decade;
                if m <= bound {
                    return bound;
                }
            }
            10.0 * decade
        }
    }
}

/// The span the colour ramp is normalised over, in the colour channel's own
/// units.
///
/// Returns `(lo, hi)` with `hi > lo` guaranteed, so a catalogue where every
/// displayed row shares one colour value still gets a defined `t` rather than a
/// divide-by-zero.
fn color_span(rows: &[&CoverageRow], encoding: Encoding) -> (f64, f64) {
    let lo = rows
        .iter()
        .map(|r| encoding.color_of(r))
        .fold(f64::INFINITY, f64::min);
    let hi = rows
        .iter()
        .map(|r| encoding.color_of(r))
        .fold(f64::NEG_INFINITY, f64::max);
    if !lo.is_finite() || !hi.is_finite() {
        return (0.0, 1.0);
    }
    if hi > lo { (lo, hi) } else { (lo, lo + 1.0) }
}

fn normalise(v: f64, lo: f64, hi: f64) -> f64 {
    let span = (hi - lo).max(1e-12);
    ((v - lo) / span).clamp(0.0, 1.0)
}

fn draw_empty_placeholder<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    scale: f64,
) -> Result<()>
where
    <DB as DrawingBackend>::ErrorType: 'static,
{
    let sp = |v: f64| common::sp(v, scale);
    let su = |v: f64| common::su(v, scale);
    let (_, h) = area.dim_in_pixel();
    area.draw_text(
        EMPTY_TEXT,
        &("sans-serif", su(TITLE_FONT)).into_font().color(&BLACK),
        (sp(20.0), (h as i32 / 2) - sp(20.0)),
    )
    .map_err(|e| anyhow!("empty placeholder draw: {e:?}"))?;
    Ok(())
}

fn draw_legend<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    rows: &[&CoverageRow],
    encoding: Encoding,
    scale: f64,
) -> Result<()>
where
    <DB as DrawingBackend>::ErrorType: 'static,
{
    let sp = |v: f64| common::sp(v, scale);
    let su = |v: f64| common::su(v, scale);
    let font = ("sans-serif", su(LABEL_FONT)).into_font().color(&BLACK);

    // Colorbar for whichever quantity the encoding puts on colour. Hit-count
    // ticks are INTEGERS — they are counts, and a `3.5` tick would invite
    // reading them as a rate; coverage ticks carry one decimal and a `%`.
    area.draw_text(encoding.colorbar_title(), &font, (sp(6.0), sp(28.0)))
        .map_err(|e| anyhow!("colorbar title: {e:?}"))?;

    let bar_left = sp(6.0);
    let bar_right = sp(40.0);
    let bar_top = sp(60.0);
    let bar_height = sp(180.0);
    let steps = 40;
    for i in 0..steps {
        let t_centre = 1.0 - ((i as f64 + 0.5) / steps as f64);
        let color = common::t_to_color(t_centre);
        let y0 = bar_top + ((bar_height as f64) * (i as f64 / steps as f64)) as i32;
        let y1 = bar_top + ((bar_height as f64) * (i + 1) as f64 / steps as f64) as i32;
        area.draw(&Rectangle::new(
            [(bar_left, y0), (bar_right, y1)],
            color.filled(),
        ))
        .map_err(|e| anyhow!("colorbar draw: {e:?}"))?;
    }

    let (lo, hi) = color_span(rows, encoding);
    for value in colorbar_ticks(lo, hi, encoding) {
        let t = normalise(value, lo, hi);
        let y = bar_top + bar_height - (t * bar_height as f64).round() as i32;
        area.draw_text(
            &format_color_tick(value, encoding),
            &font,
            (bar_right + sp(4.0), y - sp(5.0)),
        )
        .map_err(|e| anyhow!("colorbar tick: {e:?}"))?;
    }

    // Entry-size legend below the colorbar.
    let title_y = bar_top + bar_height + sp(25.0);
    area.draw_text("Entry size", &font, (sp(6.0), title_y))
        .map_err(|e| anyhow!("size title: {e:?}"))?;

    let base_radius = (DOT_BASE_RADIUS * scale).max(2.0);
    let max_size = rows.iter().map(|r| r.entry_size).max().unwrap_or(1).max(1);
    let mut y = title_y + sp(36.0);
    for size in size_legend_values(rows) {
        let r = (base_radius * (size as f64 / max_size as f64).sqrt()).max(2.0) as i32;
        area.draw(&Circle::new(
            (sp(14.0), y),
            r,
            RGBColor(80, 80, 80).filled(),
        ))
        .map_err(|e| anyhow!("size dot: {e:?}"))?;
        area.draw_text(&size.to_string(), &font, (sp(28.0), y - sp(5.0)))
            .map_err(|e| anyhow!("size label: {e:?}"))?;
        y += sp(40.0);
    }
    Ok(())
}

/// Colorbar ticks spanning `lo..=hi`, at most five so a wide span does not crowd
/// the bar.
///
/// On the hit-count channel the ticks are INTEGERS — stepping through a count in
/// fractions would misrepresent it — so the step is integer-rounded and
/// duplicates are dropped, which is what lets a 3-wide span produce 3 ticks
/// rather than 5 repeats.
fn colorbar_ticks(lo: f64, hi: f64, encoding: Encoding) -> Vec<f64> {
    if hi <= lo {
        return vec![lo];
    }
    match encoding {
        Encoding::CoverageOnX => {
            // Colour = hits: integer ticks.
            let (l, h) = (lo.round() as i64, hi.round() as i64);
            let span = (h - l).max(1);
            let n = (span + 1).min(5);
            let mut out: Vec<f64> = (0..n)
                .map(|i| (l + (span * i) / (n - 1).max(1)) as f64)
                .collect();
            out.dedup_by(|a, b| (*a - *b).abs() < 1e-9);
            out
        }
        Encoding::HitsOnX => {
            // Colour = coverage %: five evenly-spaced fractional ticks.
            (0..5).map(|i| lo + (hi - lo) * i as f64 / 4.0).collect()
        }
    }
}

fn format_color_tick(value: f64, encoding: Encoding) -> String {
    match encoding {
        Encoding::CoverageOnX => format!("{value:.0}"),
        Encoding::HitsOnX => format!("{value:.0} %"),
    }
}

/// At least three representative sizes when the data allows, so a reader can
/// calibrate the radius scale.
fn size_legend_values(rows: &[&CoverageRow]) -> Vec<usize> {
    let Some(max) = rows.iter().map(|r| r.entry_size).max() else {
        return vec![];
    };
    let min = rows
        .iter()
        .map(|r| r.entry_size)
        .min()
        .unwrap_or(max)
        .max(1);
    if min == max {
        vec![max]
    } else if max - min >= 3 {
        vec![min, (min + max) / 2, max]
    } else {
        vec![min, max]
    }
}

fn draw_annotation_strip<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    result: &CoverageResult,
    opts: &CoverageDotplotOpts,
    displayed: usize,
    scale: f64,
) -> Result<()>
where
    <DB as DrawingBackend>::ErrorType: 'static,
{
    let sp = |v: f64| common::sp(v, scale);
    let su = |v: f64| common::su(v, scale);

    let lines = build_strip(result, opts, displayed);
    let (_, rh) = root.dim_in_pixel();
    let font = ("sans-serif", su(LABEL_FONT)).into_font().color(&BLACK);
    let line_h = sp(LINE_HEIGHT_PX_BASE);
    let n = lines.len() as i32;
    for (i, line) in lines.iter().enumerate() {
        let from_bottom = n - 1 - i as i32;
        let y = (rh as i32) - sp(24.0) - from_bottom * line_h;
        root.draw_text(line, &font, (sp(10.0), y))
            .map_err(|e| anyhow!("annotation strip draw: {e:?}"))?;
    }
    Ok(())
}

/// The annotation strip, as plain strings.
///
/// A testable seam, exactly as the enrichment plot exposes for its own strip:
/// the no-statistical-vocabulary guarantee has to be asserted against the
/// ASSEMBLED STRINGS, because text is unrecoverable from an RGBA buffer.
///
/// `entries_total` is the whole catalogue as computed, before any display
/// filter — so `20 / 318` tells the reader how much of the catalogue they are
/// looking at, not how much survived a filter they can already see the value of.
pub(crate) fn build_strip(
    result: &CoverageResult,
    opts: &CoverageDotplotOpts,
    displayed: usize,
) -> Vec<String> {
    let line1 = format!(
        "{} · {}  ·  {} detected KEGG compounds",
        opts.mode_label, opts.target_label, opts.detected_total
    );

    // The compact group record. A full group list would crowd an already-dense
    // strip; the CSV carries the complete one. But a plot whose group selection
    // is invisible would let two visually different charts of the same dataset
    // look equally authoritative.
    let line2 = match opts.group_record {
        Some((selected, total, threshold)) => format!(
            "groups: {selected} of {total} · detected in >= {:.0}%",
            threshold * 100.0
        ),
        None => "groups: none supplied".to_string(),
    };

    let line3 = format!(
        "Showing {displayed} of {} entries  ·  min entry size {}  ·  min hits {}  ·  top {}",
        result.entries_total,
        opts.filters.min_entry_size,
        opts.filters.min_hit_count,
        opts.filters.top_n
    );

    vec![line1, line2, line3, NO_STATS_NOTE.to_string()]
}

#[cfg(test)]
mod tests {
    use super::*;

    fn row(id: &str, entry_size: usize, hits: usize) -> CoverageRow {
        CoverageRow {
            entry_id: id.to_string(),
            entry_name: format!("Entry {id}"),
            entry_size,
            hits,
            coverage: if entry_size == 0 {
                0.0
            } else {
                hits as f64 / entry_size as f64
            },
            share: 0.0,
            hit_compounds: vec![],
        }
    }

    fn result_of(rows: Vec<CoverageRow>, entries_total: usize) -> CoverageResult {
        CoverageResult {
            entries_without_compounds: rows.iter().filter(|r| r.entry_size == 0).count(),
            rows,
            detected_total: 391,
            entries_total,
            detected_in_entries: 264,
        }
    }

    fn opts() -> CoverageDotplotOpts {
        CoverageDotplotOpts {
            width_px: 800,
            height_px: 600,
            filters: DisplayFilters {
                min_entry_size: 3,
                min_hit_count: 1,
                sort_key: crate::coverage::CoverageSortKey::Coverage,
                top_n: 20,
            },
            mode_label: "Pathway".into(),
            target_label: "hsa".into(),
            detected_total: 391,
            group_record: Some((2, 3, 0.5)),
        }
    }

    /// The buffer is exactly `w * h * 4` bytes, and something was actually drawn.
    ///
    /// This deliberately does NOT compare two renders byte-for-byte. "Preview and
    /// export are the same image" holds because `export_coverage_dotplot_png` calls
    /// `render_coverage_dotplot` — one code path, true by construction — not because
    /// the rasterizer is reproducible. Glyph rasterization runs through plotters →
    /// font-kit → the platform font stack, which is not guaranteed to return
    /// identical bytes for repeated calls; two consecutive renders were observed to
    /// differ on Windows. Owner: the `coverage-dot-plot` capability.
    #[test]
    fn buffer_size_matches_requested_dimensions() {
        let result = result_of(vec![row("a", 42, 18), row("b", 10, 3)], 318);
        let o = CoverageDotplotOpts {
            width_px: 1200,
            height_px: 900,
            ..opts()
        };
        let buf = render_coverage_dotplot(&result, &o).expect("renders");
        assert_eq!(buf.len(), 1200 * 900 * 4);
        // Smoke check: not all zeros, so the renderer drew something rather than
        // silently handing back the blank destination buffer.
        assert!(
            buf.iter().any(|&b| b != 0),
            "buffer is all zeros — nothing rendered"
        );
    }

    /// Two renders of one input agree, compared by digest.
    ///
    /// This is the assertion `correct-coverage-dotplot-reproducibility-claim`
    /// had to remove: it failed reproducibly on Windows while passing on macOS
    /// and Linux, because glyphs were rasterized by the host font stack. With
    /// `plot-typography` the font is embedded and rasterized in pure Rust, so
    /// the property should now hold everywhere — **this test is what verifies
    /// that, and Windows CI is the only place it can be verified.**
    ///
    /// Compared by digest, never `assert_eq!` on the buffers: formatting two
    /// multi-megabyte `Vec<u8>`s into a panic message is what made the original
    /// failure undiagnosable and truncated the CI log.
    #[test]
    fn repeated_renders_agree_by_digest() {
        use std::collections::hash_map::DefaultHasher;
        use std::hash::{Hash, Hasher};

        fn digest(buf: &[u8]) -> u64 {
            let mut h = DefaultHasher::new();
            buf.hash(&mut h);
            h.finish()
        }

        let result = result_of(vec![row("a", 42, 18), row("b", 10, 3)], 318);
        let o = opts();
        let first = digest(&render_coverage_dotplot(&result, &o).expect("first render"));
        let second = digest(&render_coverage_dotplot(&result, &o).expect("second render"));
        assert_eq!(
            first, second,
            "two renders of identical input produced different pixels"
        );
    }

    /// Preview size and export size go through the one renderer: the same
    /// `CoverageResult` renders at both, each to its own exact byte count, without
    /// panicking. Mirrors the volcano renderer's equivalent test.
    #[test]
    fn preview_and_export_sizes_share_one_code_path() {
        let result = result_of(vec![row("a", 42, 18), row("b", 10, 3)], 318);
        let preview = CoverageDotplotOpts {
            width_px: 800,
            height_px: 600,
            ..opts()
        };
        let export = CoverageDotplotOpts {
            width_px: 2400,
            height_px: 1800,
            ..opts()
        };
        let p = render_coverage_dotplot(&result, &preview).expect("preview renders");
        let e = render_coverage_dotplot(&result, &export).expect("export renders");
        assert_eq!(p.len(), 800 * 600 * 4);
        assert_eq!(e.len(), 2400 * 1800 * 4);
        assert!(p.iter().any(|&b| b != 0), "preview buffer is all zeros");
        assert!(e.iter().any(|&b| b != 0), "export buffer is all zeros");
    }

    /// Zero surviving rows is a valid image with a message, never an error and
    /// never a panic. An over-tight filter is something the user did.
    #[test]
    fn an_empty_selection_renders_a_message_not_an_error() {
        let result = result_of(vec![row("a", 42, 1)], 318);
        let mut o = opts();
        o.filters.min_hit_count = 99;
        let buf = render_coverage_dotplot(&result, &o).expect("empty selection still renders");
        assert_eq!(buf.len(), 800 * 600 * 4);
        // The strip still reports the truth: nothing displayed.
        let strip = build_strip(&result, &o, 0);
        assert!(strip[2].starts_with("Showing 0 of 318 entries"));
    }

    /// **The rows drawn equal the rows the table shows** — both go through
    /// `displayed_rows` with the same `DisplayFilters`, so this is a structural
    /// property, not two implementations agreeing.
    #[test]
    fn plot_rows_equal_the_shared_chain_output() {
        let result = result_of(
            vec![
                row("small", 2, 2),   // dropped: entry_size < 3
                row("nohits", 10, 0), // dropped: hits < 1
                row("best", 4, 4),    // coverage 1.00
                row("mid", 10, 5),    // coverage 0.50
            ],
            318,
        );
        let o = opts();
        let rows = displayed_rows(&result, o.filters);
        assert_eq!(
            rows.iter().map(|r| &r.entry_id).collect::<Vec<_>>(),
            vec!["best", "mid"]
        );
        // And the renderer accepts exactly that selection without complaint.
        assert!(render_coverage_dotplot(&result, &o).is_ok());
    }

    /// No statistical vocabulary reaches ANY assembled string. Asserted against
    /// the strings, because text cannot be read back out of an RGBA buffer.
    #[test]
    fn no_statistical_vocabulary_reaches_the_charts_text() {
        let result = result_of(vec![row("a", 42, 18)], 318);
        let o = opts();
        let mut strings = build_strip(&result, &o, 1);
        // Every other author-controlled string on the chart.
        strings.extend([
            "Coverage (% of entry's compounds detected)".to_string(),
            "KEGG entry".to_string(),
            "Hits".to_string(),
            "Entry size".to_string(),
            EMPTY_TEXT.to_string(),
        ]);
        for s in &strings {
            let lower = s.to_lowercase();
            for banned in ["enrichment", "fdr", "q-value", "p-value"] {
                assert!(
                    !lower.contains(banned),
                    "chart text {s:?} contains {banned:?}"
                );
            }
        }
    }

    /// The strip records the run context, the filters, and the note.
    #[test]
    fn the_strip_records_context_filters_and_the_note() {
        let result = result_of(vec![row("a", 42, 18)], 318);
        let strip = build_strip(&result, &opts(), 20);
        assert_eq!(strip[0], "Pathway · hsa  ·  391 detected KEGG compounds");
        assert_eq!(strip[1], "groups: 2 of 3 · detected in >= 50%");
        assert_eq!(
            strip[2],
            "Showing 20 of 318 entries  ·  min entry size 3  ·  min hits 1  ·  top 20"
        );
        assert_eq!(strip[3], "Descriptive coverage — no statistical test");
    }

    /// A run with no metadata `.csv` says so, rather than omitting the term and
    /// leaving the reader to assume every sample was included.
    #[test]
    fn a_no_csv_run_says_so_on_the_strip() {
        let result = result_of(vec![row("a", 42, 18)], 318);
        let o = CoverageDotplotOpts {
            group_record: None,
            ..opts()
        };
        assert_eq!(build_strip(&result, &o, 1)[1], "groups: none supplied");
    }

    /// The percent axis is capped at 100 — the quantity cannot exceed it.
    #[test]
    fn the_percent_x_bound_rounds_up_and_caps_at_100() {
        let e = Encoding::CoverageOnX;
        assert_eq!(x_bound(e, 0.0), 5.0);
        assert_eq!(x_bound(e, 4.2), 5.0);
        assert_eq!(x_bound(e, 42.9), 50.0);
        assert_eq!(x_bound(e, 100.0), 100.0);
        assert_eq!(x_bound(e, f64::NAN), 5.0);
    }

    /// A hit count has no ceiling, so its axis gets a decade-scaled ladder.
    #[test]
    fn the_hits_x_bound_scales_with_magnitude() {
        let e = Encoding::HitsOnX;
        assert_eq!(x_bound(e, 3.0), 5.0);
        assert_eq!(x_bound(e, 18.0), 20.0);
        assert_eq!(x_bound(e, 140.0), 200.0);
        assert_eq!(x_bound(e, 640.0), 1000.0);
        assert_eq!(x_bound(e, f64::NAN), 5.0);
    }

    /// A single colour value across every row still yields a defined `t` rather
    /// than a divide-by-zero, on either encoding.
    #[test]
    fn a_degenerate_color_span_does_not_divide_by_zero() {
        let a = row("a", 10, 5);
        let b = row("b", 20, 5);
        let rows = vec![&a, &b];
        for e in [Encoding::CoverageOnX, Encoding::HitsOnX] {
            let (lo, hi) = color_span(&rows, e);
            assert!(hi > lo, "{e:?}");
            assert!(normalise(e.color_of(&a), lo, hi).is_finite());
        }
        assert_eq!(color_span(&[], Encoding::CoverageOnX), (0.0, 1.0));
    }

    /// Hit-count colorbar ticks are integers, at most five, spanning the data.
    #[test]
    fn colorbar_ticks_are_integers_when_colour_carries_hits() {
        let e = Encoding::CoverageOnX;
        let ticks = colorbar_ticks(1.0, 18.0, e);
        assert_eq!(ticks.first(), Some(&1.0));
        assert_eq!(ticks.last(), Some(&18.0));
        assert!(ticks.len() <= 5);
        assert!(ticks.iter().all(|t| t.fract() == 0.0), "{ticks:?}");
        assert_eq!(colorbar_ticks(4.0, 4.0, e), vec![4.0]);
        assert_eq!(colorbar_ticks(1.0, 2.0, e), vec![1.0, 2.0]);
    }

    /// **Item 6.** The sort key drives the encoding: `Coverage` keeps X =
    /// coverage % / colour = hits; `Hits` swaps both channels together. Marker
    /// size stays `entry_size` on each, because it is the denominator behind
    /// both quantities.
    #[test]
    fn the_sort_key_selects_the_encoding() {
        assert_eq!(
            Encoding::from_sort_key(CoverageSortKey::Coverage),
            Encoding::CoverageOnX
        );
        assert_eq!(
            Encoding::from_sort_key(CoverageSortKey::Hits),
            Encoding::HitsOnX
        );
        // The two keys the UI does not offer fall back to the default rather
        // than inventing a third encoding.
        for key in [CoverageSortKey::EntrySize, CoverageSortKey::EntryId] {
            assert_eq!(Encoding::from_sort_key(key), Encoding::CoverageOnX);
        }

        let r = row("a", 42, 18); // coverage 18/42 ≈ 42.857 %
        assert!((Encoding::CoverageOnX.x_of(&r) - 100.0 * 18.0 / 42.0).abs() < 1e-9);
        assert_eq!(Encoding::CoverageOnX.color_of(&r), 18.0);
        assert_eq!(Encoding::HitsOnX.x_of(&r), 18.0);
        assert!((Encoding::HitsOnX.color_of(&r) - 100.0 * 18.0 / 42.0).abs() < 1e-9);
    }

    /// The axis title and colorbar title swap with the encoding, so a reader
    /// can never mistake which quantity they are looking at.
    #[test]
    fn the_axis_and_colorbar_titles_swap_with_the_encoding() {
        assert_eq!(Encoding::CoverageOnX.colorbar_title(), "Hits");
        assert_eq!(Encoding::HitsOnX.colorbar_title(), "Coverage %");
        assert!(Encoding::CoverageOnX.x_desc().starts_with("Coverage"));
        assert!(Encoding::HitsOnX.x_desc().starts_with("Hits"));
    }

    /// Both encodings render at the requested size without error.
    #[test]
    fn both_encodings_render() {
        let result = result_of(vec![row("a", 42, 18), row("b", 10, 3)], 318);
        for key in [CoverageSortKey::Coverage, CoverageSortKey::Hits] {
            let mut o = opts();
            o.filters.sort_key = key;
            let buf = render_coverage_dotplot(&result, &o).expect("renders");
            assert_eq!(buf.len(), 800 * 600 * 4);
        }
    }

    /// The exported PNG declares its DPI in a `pHYs` chunk, through the same
    /// shared encoder the volcano and enrichment plots use.
    #[test]
    fn export_writes_phys_chunk_with_dpi() {
        let result = result_of(vec![row("a", 42, 18), row("b", 10, 3)], 318);
        let o = CoverageDotplotOpts {
            width_px: 200,
            height_px: 200,
            ..opts()
        };
        let tmp = tempfile::Builder::new()
            .suffix(".png")
            .tempfile()
            .expect("tempfile");
        export_coverage_dotplot_png(&result, &o, 300, tmp.path()).expect("export");
        let bytes = std::fs::read(tmp.path()).expect("read png");

        let mut idx = 8; // skip signature
        let mut found: Option<(u32, u32, u8)> = None;
        while idx + 8 < bytes.len() {
            let len = u32::from_be_bytes(bytes[idx..idx + 4].try_into().unwrap()) as usize;
            let ty = &bytes[idx + 4..idx + 8];
            if ty == b"pHYs" {
                let payload = &bytes[idx + 8..idx + 8 + len];
                found = Some((
                    u32::from_be_bytes(payload[0..4].try_into().unwrap()),
                    u32::from_be_bytes(payload[4..8].try_into().unwrap()),
                    payload[8],
                ));
                break;
            }
            idx += 8 + len + 4;
        }
        let (xppu, yppu, unit) = found.expect("pHYs chunk not found");
        let expected = (300.0_f64 / 0.0254).round() as u32;
        assert_eq!(xppu, expected);
        assert_eq!(yppu, expected);
        assert_eq!(unit, 1);
    }

    /// The DPI changes only the metadata, never a pixel — which is why it is a
    /// parameter of the export rather than a field of the opts.
    #[test]
    fn dpi_does_not_change_the_pixels() {
        let result = result_of(vec![row("a", 42, 18)], 318);
        let o = CoverageDotplotOpts {
            width_px: 200,
            height_px: 200,
            ..opts()
        };
        // `render_coverage_dotplot` takes no DPI at all, so the two exports can
        // differ only in the chunk. Assert the shared buffer directly.
        let a = render_coverage_dotplot(&result, &o).expect("renders");
        let b = render_coverage_dotplot(&result, &o).expect("renders");
        assert_eq!(a, b);
    }

    /// Three representative sizes when the range allows, fewer when it does not.
    #[test]
    fn the_size_legend_offers_representative_values() {
        let a = row("a", 4, 1);
        let b = row("b", 42, 1);
        assert_eq!(size_legend_values(&[&a, &b]), vec![4, 23, 42]);
        assert_eq!(size_legend_values(&[&a]), vec![4]);
        assert_eq!(size_legend_values(&[]), Vec::<usize>::new());
    }
}
