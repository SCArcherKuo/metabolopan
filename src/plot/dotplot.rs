//! Enrichment dot plot renderer. Same chrome conventions as the volcano
//! plot (white background, black axes, no minor gridlines, scaled
//! baseline, pHYs DPI metadata) but with categorical Y axis (entry
//! names — pathways or modules depending on mode) and size-and-color
//! encoded markers.
//!
//! theme: scope-excluded. FDR colour ramp = biological convention, not brand palette.
//!
//! The renderer is the single source of truth for both the in-window
//! preview and the on-disk PNG export. The chart-agnostic scaffold shared
//! with the volcano renderer (RGBA expansion, PNG + pHYs encoding, the
//! design-scale baseline, `sp` / `su`, the major-gridline grey) now lives in
//! `crate::plot::common` — the D17 deferral was lifted by
//! `extract-plot-common-scaffold`.

use anyhow::{Result, anyhow};
use plotters::prelude::*;
use plotters::style::text_anchor::{HPos, Pos, VPos};
use std::path::Path;

use crate::enrichment::types::EnrichmentResult;
use crate::plot::common;

/// Options for `render_dotplot` and `export_dotplot_png`.
#[derive(Debug, Clone, Copy)]
pub struct DotplotOpts {
    pub width_px: u32,
    pub height_px: u32,
    pub fdr_threshold: f64,
    pub top_n: usize,
    /// Live display filter on hit count, read from settings at draw time — NOT
    /// from the result. Storing it per-row froze this filter to the run while
    /// the threshold and `top_n` stayed live; Owner: the `enrichment-dot-plot` capability.
    pub min_hit_count: usize,
    /// Correction method that produced `result.rows[*].fdr`. Affects only the
    /// chrome: the colorbar title (`-log10(FDR (BH))` / `-log10(FDR (BY))` /
    /// `-log10(p-value)`), the annotation strip's `Significance:` line, and the
    /// empty-state placeholder's metric noun. Does NOT alter marker geometry,
    /// the colorbar gradient, axes, or any other pixel.
    ///
    /// Callers MUST pass `result.fdr_method` — the method that fed the
    /// `run_ora` call which produced the result being plotted — never a live
    /// setting, which describes the NEXT run.
    pub fdr_method: crate::dam::fdr::FdrMethod,
    /// Singular lowercase noun for the entity being enriched
    /// (`"pathway"` in pathway mode, `"module"` in module mode). Drives
    /// the Y-axis title (`KEGG <entry_label>`) and the empty-state
    /// placeholder (`No <entry_label>s passed <metric> < ...`); previously
    /// hardcoded to "pathways" which read as a factual error in module-
    /// mode publications. Callers obtain this via
    /// `AnalysisMode::entry_label_singular`.
    pub entry_label: &'static str,
}

/// FDR floor used when computing -log10. BH/BY may underflow exactly to 0.0
/// on huge `|t|`; -log10(0) = +infinity would crash the colorbar lerp.
/// 1e-300 is below f64 subnormal noise but well above f64::MIN_POSITIVE.
const FDR_FLOOR: f64 = 1e-300;

/// Default colorbar span (in -log10 units) when the displayed rows
/// collapse to a single -log10 value (e.g. only one row, or every row
/// shares an FDR). Gives the gradient a 2-decade visual range so the
/// colorbar doesn't degenerate to a single colour.
const COLORBAR_FALLBACK_SPAN: f64 = 2.0;

const REF_LINE: RGBColor = RGBColor(150, 150, 150);

/// Horizontal character budget per wrapped Y-axis label line. Tuned so a
/// ~42-char string at our 14 pt sans-serif fills the widened 286-px
/// `y_label_area_size` (the name column is intentionally wide — the plot
/// and legend shift right to give entry names room). Long names that
/// exceed this on a single line are wrapped at word boundaries (see
/// `wrap_label`).
const CHARS_PER_LINE: usize = 42;

/// Maximum wrapped lines per Y-axis label. Four lines lets even the
/// longest KEGG module names render in full; anything longer gets
/// ellipsis-truncated on the last line. `draw_chart` further caps each
/// row at `floor(row_step_px / line_height_px)` so short rows (large
/// top_n) auto-degrade to fewer lines instead of overlapping neighbours.
const MAX_LABEL_LINES: usize = 4;

/// Pixel height of a single label line at scale 1.0. 14 pt font + leading
/// ≈ 16 px (the standard ~1.15 line-height ratio), tracking `LABEL_FONT`.
const LINE_HEIGHT_PX_BASE: f64 = 16.0;

/// Y-axis label / tick / annotation-strip / legend font (the "small"
/// tier).
const LABEL_FONT: f64 = 14.0;

/// Axis-title / empty-placeholder font (the "large" tier; the X/Y axis
/// titles are the primary use). The legend titles deliberately use
/// `LABEL_FONT`, not this, so the whole legend renders at the small size.
const TITLE_FONT: f64 = 20.0;

/// Base marker radius at scale 1.0, shared by the chart dots and the
/// legend's Hits reference dots so a given hit count renders at the SAME
/// size in both (`radius = DOT_BASE_RADIUS * scale * sqrt(hits/max_hits)`,
/// floored at 2 px). Both sites MUST use this constant — drift here is
/// exactly the chart-vs-legend size mismatch this constant prevents.
const DOT_BASE_RADIUS: f64 = 10.0;

/// Render the dot plot into an RGBA buffer of `width * height * 4` bytes.
pub fn render_dotplot(result: &EnrichmentResult, opts: &DotplotOpts) -> Result<Vec<u8>> {
    common::ensure_font_registered();
    let w = opts.width_px;
    let h = opts.height_px;
    let pixel_count = (w as usize) * (h as usize);
    let mut rgb = vec![0u8; pixel_count * 3];

    // Width-based scale: the dot plot's height auto-sizes to the displayed-row
    // count, so keying the font/element scale off `min(w, h)` (the shared
    // `design_scale`) shrank every font on sparse results. Scaling by the fixed
    // width keeps fonts constant per (width, dpi), independent of entry count.
    let scale = common::design_scale_by_width(w);

    // Reserve a fixed-width legend on the right. Trimmed from 140 → 105 so
    // the main plotting area reclaims ~10% width; the narrow colorbar + the
    // `-log10(FDR (…))` title still fit (verified against the widest title).
    //
    // Widened 105 → 120 by `embed-plot-font`: DejaVu Sans sets
    // `-log10(FDR (BH))` wider than the host face this was tuned against, and
    // at the default 1050 px the title's right margin had collapsed from 17 px
    // to 1 px — not yet clipped, but one glyph away from it. Measured, not
    // eyeballed; the coverage dot plot keeps 105 because its widest colorbar
    // title is `Coverage %`, which leaves ~55 px to spare.
    let legend_w = ((120.0 * scale).round() as u32).min(w * 3 / 10).max(120);
    let chart_w_boundary = (w - legend_w) as i32;

    // Rows the dot plot draws, ordered top-to-bottom: the most-significant
    // `top_n` are selected, then arranged by descending fold enrichment so the
    // largest-fold-enrichment entry sits at the top of the y axis.
    let filtered =
        select_and_order_rows(result, opts.fdr_threshold, opts.min_hit_count, opts.top_n);

    // -log10 colour-mapping span. Computed ONCE so chart dots and the
    // legend colorbar use the same normalisation — without this they
    // would silently disagree.
    let (nl_threshold, max_nl) = neg_log10_span(&filtered, opts.fdr_threshold);

    {
        let root = BitMapBackend::with_buffer(&mut rgb, (w, h)).into_drawing_area();
        root.fill(&WHITE)
            .map_err(|e| anyhow!("fill background: {e:?}"))?;

        let (chart_area, legend_area) = root.split_horizontally(chart_w_boundary);

        if filtered.is_empty() {
            draw_empty_placeholder(
                &chart_area,
                opts.entry_label,
                opts.fdr_method,
                opts.fdr_threshold,
                scale,
            )?;
        } else {
            draw_chart(
                &chart_area,
                &filtered,
                opts.entry_label,
                nl_threshold,
                max_nl,
                scale,
            )?;
        }
        draw_legend(&legend_area, &filtered, opts, nl_threshold, max_nl, scale)?;
        draw_annotation_strip(&root, result, opts.min_hit_count, opts.entry_label, scale)?;

        root.present().map_err(|e| anyhow!("present: {e:?}"))?;
    }

    Ok(common::rgb_to_rgba(&rgb, pixel_count))
}

/// Render and save as PNG, embedding the requested DPI in the pHYs chunk.
pub fn export_dotplot_png(
    result: &EnrichmentResult,
    opts: &DotplotOpts,
    dpi: u32,
    out: &Path,
) -> Result<()> {
    let buffer = render_dotplot(result, opts)?;
    common::encode_png(&buffer, opts.width_px, opts.height_px, dpi, out)
}

/// Select the rows the dot plot draws and order them top-to-bottom.
///
/// Selection (significance) is unchanged from the ORA sort: rows arrive
/// FDR-ascending from `enrichment::ora`, so filtering on
/// `hits >= min_hit_count && fdr < threshold` and truncating to `top_n` keeps the
/// **most significant** `top_n` entries. The kept rows are then re-ordered for
/// DISPLAY by descending fold enrichment (`enrichment_ratio`) so the entry with
/// the largest fold enrichment sits at the top of the y axis — matching the
/// clusterProfiler convention of arranging the y axis by the x-axis metric.
/// Ties break by ascending FDR, then ascending `entry_id` (deterministic).
/// Non-finite ratios are handled gracefully by the `unwrap_or(Equal)` fallback
/// (deferring to the tie-break); in practice displayed rows always have a
/// finite, positive ratio — `hits >= 1` ⇒ `K >= 1` and `total >= 1` ⇒
/// `expected > 0`.
fn select_and_order_rows(
    result: &EnrichmentResult,
    fdr_threshold: f64,
    min_hit_count: usize,
    top_n: usize,
) -> Vec<&crate::enrichment::types::EnrichmentRow> {
    let mut rows: Vec<&crate::enrichment::types::EnrichmentRow> = result
        .rows
        .iter()
        .filter(|r| r.hits >= min_hit_count && r.fdr < fdr_threshold)
        .collect();
    if rows.len() > top_n {
        rows.truncate(top_n);
    }
    rows.sort_by(|a, b| {
        b.enrichment_ratio
            .partial_cmp(&a.enrichment_ratio)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| {
                a.fdr
                    .partial_cmp(&b.fdr)
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .then_with(|| a.entry_id.cmp(&b.entry_id))
    });
    rows
}

fn draw_chart<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    rows: &[&crate::enrichment::types::EnrichmentRow],
    entry_label: &str,
    nl_threshold: f64,
    max_nl: f64,
    scale: f64,
) -> Result<()>
where
    <DB as DrawingBackend>::ErrorType: 'static,
{
    let sp = |v: f64| common::sp(v, scale);
    let su = |v: f64| common::su(v, scale);

    // Reserve a band at the bottom for the multi-line annotation block.
    let (ah_w, ah_h) = area.dim_in_pixel();
    let strip_h = su(110.0) as i32;
    let split_y = (ah_h as i32 - strip_h).max(1);
    let (chart_main, _strip_area) = area.split_vertically(split_y);

    let max_ratio = rows
        .iter()
        .map(|r| {
            if r.enrichment_ratio.is_finite() {
                r.enrichment_ratio
            } else {
                0.0
            }
        })
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let x_max = (max_ratio * 1.05).max(1.5);

    let max_hits = rows.iter().map(|r| r.hits).max().unwrap_or(1) as f64;

    // Lay out pathways top-to-bottom; row 0 (largest fold enrichment, per
    // `select_and_order_rows`) at the top, each at the CENTRE of its 1.0-tall
    // band on a CONTINUOUS 0..n Y axis. The
    // earlier `(0..n).into_segmented()` was dropped: plotters' SegmentedCoord
    // reserves a full `plot_height / n` margin between the last category and
    // `Last` (plotters-0.3.7 discrete.rs::map) — rendered as a blank row above
    // the top entry. A continuous axis leaves only a symmetric half-row pad
    // top & bottom. Cartesian Y grows upward, so row i sits at
    // y = (n - 1 - i) + 0.5 (row 0 highest).
    let n = rows.len();
    let row_center = |i: usize| -> f64 { (n - 1 - i) as f64 + 0.5 };

    let mut chart = ChartBuilder::on(&chart_main)
        .margin(sp(5.0))
        .x_label_area_size(su(80.0))
        .y_label_area_size(su(286.0))
        .build_cartesian_2d(0.0..x_max, 0f64..(n as f64))
        .map_err(|e| anyhow!("chart build: {e:?}"))?;

    // X mesh only. The Y gridlines are hand-drawn through the row centres
    // below — a continuous-axis mesh would place them at integer boundaries
    // (BETWEEN rows). Y labels are hand-drawn too, so plotters' own are
    // suppressed via `disable_y_mesh` + an empty formatter.
    chart
        .configure_mesh()
        // "Fold enrichment (observed / expected)" — disambiguates from
        // clusterProfiler's `GeneRatio = k/n` convention, which readers
        // trained on that tool would otherwise infer from a bare
        // "Enrichment ratio" label. The parenthetical matches the
        // formula used internally: `hits / expected = (k·N) / (m·K)`.
        .x_desc("Fold enrichment (observed / expected)")
        .y_desc(y_axis_title(entry_label))
        .label_style(("sans-serif", su(LABEL_FONT)).into_font().color(&BLACK))
        .axis_desc_style(("sans-serif", su(TITLE_FONT)).into_font().color(&BLACK))
        .disable_y_mesh()
        .y_labels(1)
        .y_label_formatter(&|_| String::new())
        .light_line_style(RGBAColor(0, 0, 0, 0.0))
        .bold_line_style(common::GRID_LIGHT)
        .axis_style(BLACK)
        .draw()
        .map_err(|e| anyhow!("mesh draw: {e:?}"))?;

    // Hand-drawn horizontal gridlines through each row centre (preserves the
    // old segmented-mesh look; drawn before labels/dots so they sit behind).
    {
        let (px_x, _px_y) = chart.plotting_area().get_pixel_range();
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

    // Hand-drawn Y labels with word-wrap. Use `chart.backend_coord` for
    // every row centre so the labels share plotters' exact value → pixel
    // translation (the same one the dots use).
    let (row0_x, row0_y) = chart.backend_coord(&(0.0_f64, row_center(0)));
    let row_step_px = if n >= 2 {
        let (_, row1_y) = chart.backend_coord(&(0.0_f64, row_center(1)));
        (row1_y - row0_y).unsigned_abs() as f64
    } else {
        // Single-row chart has no neighbour; estimate from plot height.
        let (_, py_range) = chart.plotting_area().get_pixel_range();
        (py_range.end - py_range.start).max(1) as f64
    };
    let line_height_px = sp(LINE_HEIGHT_PX_BASE);
    // Auto-fallback to single-line + ellipsis when rows are too short
    // to stack 2 lines without overlap (common when top_n is large).
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
        let (_, row_center_y) = chart.backend_coord(&(0.0_f64, row_center(i)));
        let lines = common::wrap_label(&row.entry_name, CHARS_PER_LINE, max_lines);
        let n_lines = lines.len() as f64;
        for (li, line) in lines.iter().enumerate() {
            let line_offset_y =
                ((li as f64 - (n_lines - 1.0) / 2.0) * line_height_px as f64).round() as i32;
            chart_main
                .draw_text(
                    line,
                    &label_style,
                    (label_right_x, row_center_y + line_offset_y),
                )
                .map_err(|e| anyhow!("y label draw: {e:?}"))?;
        }
    }

    // Reference line at X = 1.0 (no enrichment), spanning the full Y range.
    chart
        .draw_series(std::iter::once(PathElement::new(
            vec![(1.0_f64, 0.0_f64), (1.0_f64, n as f64)],
            ShapeStyle::from(&REF_LINE).stroke_width(su(1.0)),
        )))
        .map_err(|e| anyhow!("ref line: {e:?}"))?;

    // Dots.
    let base_radius = (DOT_BASE_RADIUS * scale).max(2.0);
    for (i, row) in rows.iter().enumerate() {
        if !row.enrichment_ratio.is_finite() {
            continue;
        }
        let r = (base_radius * (row.hits as f64 / max_hits.max(1.0)).sqrt()).max(2.0) as i32;
        let color = fdr_to_color(row.fdr, nl_threshold, max_nl);
        chart
            .draw_series(std::iter::once(Circle::new(
                (row.enrichment_ratio, row_center(i)),
                r,
                color.filled(),
            )))
            .map_err(|e| anyhow!("dot draw: {e:?}"))?;
    }

    let _ = ah_w;
    Ok(())
}

fn draw_empty_placeholder<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    entry_label: &str,
    method: crate::dam::fdr::FdrMethod,
    fdr_threshold: f64,
    scale: f64,
) -> Result<()>
where
    <DB as DrawingBackend>::ErrorType: 'static,
{
    let sp = |v: f64| common::sp(v, scale);
    let su = |v: f64| common::su(v, scale);
    let (w, h) = area.dim_in_pixel();
    let text = empty_placeholder_text(entry_label, method, fdr_threshold);
    area.draw_text(
        &text,
        &("sans-serif", su(TITLE_FONT)).into_font().color(&BLACK),
        (sp(20.0), (h as i32 / 2) - sp(20.0)),
    )
    .map_err(|e| anyhow!("empty placeholder draw: {e:?}"))?;
    let _ = w;
    Ok(())
}

fn draw_legend<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    rows: &[&crate::enrichment::types::EnrichmentRow],
    opts: &DotplotOpts,
    nl_threshold: f64,
    max_nl: f64,
    scale: f64,
) -> Result<()>
where
    <DB as DrawingBackend>::ErrorType: 'static,
{
    let sp = |v: f64| common::sp(v, scale);
    let su = |v: f64| common::su(v, scale);

    // Both legend fonts use the small tier (LABEL_FONT) so the whole
    // legend renders compact — the colorbar/Hits titles are deliberately
    // NOT bumped to TITLE_FONT.
    let title_font = ("sans-serif", su(LABEL_FONT)).into_font().color(&BLACK);
    let row_font = ("sans-serif", su(LABEL_FONT)).into_font().color(&BLACK);

    // FDR colorbar on top of the legend area. Title is in -log10 units
    // (clusterProfiler convention) so the gradient direction is
    // unambiguous: top of bar = highest -log10 = most significant = red.
    // For `NoCorrection`, the values shown ARE raw p-values, NOT q —
    // so the title drops the "FDR (…)" wrapper to stay scientifically
    // honest.
    // Title y stays at sp(28) — anchored just above the frozen bar_top
    // (60), NOT halved, so it doesn't drift up away from the colorbar.
    let colorbar_title = colorbar_title(opts.fdr_method);
    area.draw_text(&colorbar_title, &title_font, (sp(6.0), sp(28.0)))
        .map_err(|e| anyhow!("legend title FDR: {e:?}"))?;

    let bar_left = sp(6.0);
    let bar_right = sp(40.0);
    let bar_top = sp(60.0); // frozen — keeps the colorbar's vertical anchor
    let bar_height = sp(180.0); // frozen — keeps the gradient's resolution
    let steps = 40;
    for i in 0..steps {
        // i=0 paints the TOP slice (most significant end, red).
        // Use the centre t for the slice so the gradient is symmetric
        // across each rectangle.
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
    // Tick labels: both endpoints (`max_nl` at the red top, `nl_threshold`
    // at the yellow bottom) plus any integer -log10 values that fit in
    // between without crowding. `colorbar_tick_positions` greedy-filters
    // candidates against pixel proximity to already-accepted ticks so a
    // 180-px bar over [1.30, 7.0] gets 7 labels, a narrow bar gets fewer.
    let span = (max_nl - nl_threshold).max(1e-12);
    let ticks = colorbar_tick_positions(nl_threshold, max_nl, bar_height, sp(14.0));
    for nl in ticks {
        let t_frac = (nl - nl_threshold) / span;
        let pixel_y = bar_top + bar_height - (t_frac * bar_height as f64).round() as i32;
        area.draw_text(
            &format_tick(nl),
            &row_font,
            (bar_right + sp(4.0), pixel_y - sp(5.0)),
        )
        .map_err(|e| anyhow!("colorbar tick label: {e:?}"))?;
    }

    // Hits size legend below the colorbar. Title embeds the per-figure
    // max so readers comparing two side-by-side figures with different
    // max_hits can see that the dot sizes are NOT directly comparable
    // (the radius is normalised per-plot by `max_hits`, so "largest dot"
    // means a different absolute hit count in each figure).
    let hits_title_y = bar_top + bar_height + sp(25.0);
    let max_hits_label = rows.iter().map(|r| r.hits).max();
    area.draw_text(
        &hits_title(max_hits_label),
        &title_font,
        (sp(6.0), hits_title_y),
    )
    .map_err(|e| anyhow!("hits title: {e:?}"))?;

    let (h_min, h_max) = if rows.is_empty() {
        (1, 1)
    } else {
        let max_h = rows.iter().map(|r| r.hits).max().unwrap_or(1);
        let min_h = rows.iter().map(|r| r.hits).min().unwrap_or(1).max(1);
        (min_h, max_h)
    };
    let reference_hits = if h_min == h_max {
        vec![h_min]
    } else if (h_max as i64 - h_min as i64) >= 3 {
        vec![h_min, (h_min + h_max) / 2, h_max]
    } else {
        vec![h_min, h_max]
    };
    // Reference-dot vertical spacing restored to the full-size rhythm —
    // the dots are DOT_BASE_RADIUS again (not Task-2's halved 5), so the
    // halved 20-px step overlapped them. 40-px step clears radius-10 dots.
    let base_radius = (DOT_BASE_RADIUS * scale).max(2.0);
    let mut y_cursor = hits_title_y + sp(36.0);
    let row_step = sp(40.0);
    for hits in reference_hits {
        let r = (base_radius * (hits as f64 / h_max.max(1) as f64).sqrt()).max(2.0) as i32;
        area.draw(&Circle::new(
            (sp(14.0), y_cursor),
            r,
            RGBColor(80, 80, 80).filled(),
        ))
        .map_err(|e| anyhow!("hits dot: {e:?}"))?;
        area.draw_text(
            &format!("{hits}"),
            &row_font,
            (sp(28.0), y_cursor - sp(5.0)),
        )
        .map_err(|e| anyhow!("hits label: {e:?}"))?;
        y_cursor += row_step;
    }
    Ok(())
}

fn draw_annotation_strip<DB: DrawingBackend>(
    root: &DrawingArea<DB, plotters::coord::Shift>,
    result: &EnrichmentResult,
    min_hit_count: usize,
    entry_label: &str,
    scale: f64,
) -> Result<()>
where
    <DB as DrawingBackend>::ErrorType: 'static,
{
    let sp = |v: f64| common::sp(v, scale);
    let su = |v: f64| common::su(v, scale);

    let lines = build_strip(result, min_hit_count, entry_label);
    let (_, rh) = root.dim_in_pixel();
    let font = ("sans-serif", su(LABEL_FONT)).into_font().color(&BLACK);
    let line_h = sp(LINE_HEIGHT_PX_BASE);
    // Stack the lines upward from the bottom margin: the last line keeps a
    // baseline near where the single-line strip sat; earlier lines stack
    // above it, inside the reserved `strip_h` band.
    let n = lines.len() as i32;
    for (i, line) in lines.iter().enumerate() {
        let from_bottom = n - 1 - i as i32; // line 0 is highest
        let y = (rh as i32) - sp(24.0) - from_bottom * line_h;
        root.draw_text(line, &font, (sp(10.0), y))
            .map_err(|e| anyhow!("annotation strip draw: {e:?}"))?;
    }
    Ok(())
}

/// Build the bottom annotation block as plain-language lines (mode-aware
/// via `entry_label`). The bare N / K / m symbols are intentionally
/// dropped, but their numeric values are kept — in particular the tested-
/// entry count (`= result.rows.len()`, the BH/BY divisor) stays auditable
/// from the figure alone. Never calls `FdrMethod::short_label()`, which is
/// shared by the volcano strip / CSV export / UI radios.
fn build_strip(result: &EnrichmentResult, min_hit_count: usize, entry_label: &str) -> Vec<String> {
    // Line 1 — background universe (the measurable, KEGG-mapped metabolome).
    let line1 = format!(
        "Background universe = {} compounds measured and mapped to KEGG",
        result.universe_size
    );

    // Line 2 — foreground compounds + direction in plain language.
    let direction = match result.direction {
        crate::enrichment::types::EnrichmentDirection::Up => "increased",
        crate::enrichment::types::EnrichmentDirection::Down => "decreased",
        crate::enrichment::types::EnrichmentDirection::Both => "both directions",
    };
    let line2 = format!(
        "Compounds of interest = {} differentially abundant ({direction})",
        result.dam_cpd_size
    );

    // Line 3 — tested entities. `tested` is the FDR-family size; `total`
    // adds back the min_entry_size drops. Suffixes are appended in a fixed
    // order: the min_entry clause first, the hits clause always last.
    let entity_plural = {
        let mut chars = entry_label.chars();
        match chars.next() {
            Some(first) => format!("{}{}s", first.to_uppercase(), chars.as_str()),
            None => "Entriess".to_string(),
        }
    };
    let tested = result.rows.len();
    let dropped = result.entries_dropped_by_min_entry_size;
    let mut line3 = format!("{entity_plural} tested = {tested}");
    if dropped > 0 {
        let total = tested + dropped;
        line3.push_str(&format!(
            " of {total}  ·  {dropped} skipped (< {} compounds each)",
            result.min_entry_size
        ));
    } else if result.min_entry_size > 1 {
        line3.push_str(&format!(
            "  ·  each {entry_label} needs ≥ {} compounds",
            result.min_entry_size
        ));
    }
    if min_hit_count > 1 {
        line3.push_str(&format!("; ≥ {min_hit_count} hits required"));
    }

    // Line 4 — significance basis, in plain language (built locally; never
    // routes through FdrMethod::short_label).
    let line4 = match result.fdr_method {
        crate::dam::fdr::FdrMethod::NoCorrection => {
            "Significance: raw p-value (no FDR correction)".to_string()
        }
        crate::dam::fdr::FdrMethod::BenjaminiHochberg => {
            "Significance: FDR-adjusted, Benjamini–Hochberg (BH)".to_string()
        }
        crate::dam::fdr::FdrMethod::BenjaminiYekutieli => {
            "Significance: FDR-adjusted, Benjamini–Yekutieli (BY)".to_string()
        }
    };

    vec![line1, line2, line3, line4]
}

/// Clamped -log10. Returns -log10(FDR_FLOOR) ≈ 300 for FDR ≤ FDR_FLOOR,
/// preventing +infinity when an FDR underflows to exactly 0.0.
fn neg_log10_clamped(fdr: f64) -> f64 {
    if fdr <= FDR_FLOOR {
        -FDR_FLOOR.log10()
    } else {
        -fdr.log10()
    }
}

/// Compute the -log10(FDR) colorbar span `(nl_threshold, max_nl)`. The
/// chart dots and the legend gradient MUST use the same span — this
/// function is the single source of truth. Falls back to a 2-decade
/// visual span when displayed rows collapse to a single -log10 value
/// (or there are no rows), so the colorbar never degenerates.
fn neg_log10_span(
    rows: &[&crate::enrichment::types::EnrichmentRow],
    fdr_threshold: f64,
) -> (f64, f64) {
    let nl_threshold = neg_log10_clamped(fdr_threshold);
    let max_nl = rows
        .iter()
        .map(|r| neg_log10_clamped(r.fdr))
        .fold(nl_threshold, f64::max);
    let max_nl = if (max_nl - nl_threshold) <= 1e-9 {
        nl_threshold + COLORBAR_FALLBACK_SPAN
    } else {
        max_nl
    };
    (nl_threshold, max_nl)
}

/// Compute -log10(FDR) tick positions for the colorbar. Endpoints
/// (`max_nl` at the red end, `nl_threshold` at the yellow end) are
/// always included; interior integer -log10 values are added when they
/// don't crowd an already-accepted tick (greedy filter on pixel
/// distance using `min_pixel_gap_px`). Endpoints have priority so a
/// dense interior never bumps the range markers off the bar.
///
/// Returned positions are sorted ascending so the caller can iterate
/// bottom-of-bar → top-of-bar.
///
/// Example: span `[1.30, 7.0]` on a 180-px bar with `min_pixel_gap=28`
/// yields `{1.30, 3, 4, 5, 6, 7}` — integer "2" is 22 px from the 1.30
/// endpoint (below the gap), so it's dropped in favour of the endpoint
/// already labelling that region.
fn colorbar_tick_positions(
    nl_threshold: f64,
    max_nl: f64,
    bar_height_px: i32,
    min_pixel_gap_px: i32,
) -> Vec<f64> {
    let span = (max_nl - nl_threshold).max(1e-12);
    let pixel_of = |nl: f64| -> f64 { (1.0 - (nl - nl_threshold) / span) * bar_height_px as f64 };

    // Endpoints first → highest priority. Then interior integers.
    let mut candidates: Vec<f64> = vec![max_nl, nl_threshold];
    let lower = (nl_threshold + 1e-9).ceil() as i64;
    let upper = (max_nl - 1e-9).floor() as i64;
    for k in lower..=upper {
        candidates.push(k as f64);
    }

    let mut accepted: Vec<f64> = Vec::new();
    for c in candidates {
        let cy = pixel_of(c);
        let too_close = accepted
            .iter()
            .any(|a| (pixel_of(*a) - cy).abs() < min_pixel_gap_px as f64);
        if !too_close {
            accepted.push(c);
        }
    }
    accepted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    accepted
}

/// Format a colorbar tick label: integer -log10 values drop the
/// trailing `.00` ("3" reads cleaner than "3.00") while non-integer
/// endpoints keep 2 decimals so `1.30` stays informative.
fn format_tick(nl: f64) -> String {
    let rounded = nl.round();
    if (nl - rounded).abs() < 1e-9 {
        format!("{rounded:.0}")
    } else {
        format!("{nl:.2}")
    }
}

/// Map an FDR value to a colorbar colour. `t` is computed on the
/// -log10(FDR) scale so highly-significant entries (FDR ≪ threshold)
/// are visually distinguishable from borderline ones — the linear
/// `fdr / threshold` mapping that preceded this collapsed every
/// `FDR ≤ 1e-3` entry to indistinguishable red.
fn fdr_to_color(fdr: f64, nl_threshold: f64, max_nl: f64) -> RGBColor {
    let nl = neg_log10_clamped(fdr);
    let span = (max_nl - nl_threshold).max(1e-12);
    let t = ((nl - nl_threshold) / span).clamp(0.0, 1.0);
    common::t_to_color(t)
}

/// Empty-state placeholder text. Pluralises the entity by appending "s"
/// — both supported labels ("pathway" / "module") pluralise regularly.
///
/// Names the significance quantity by the method that produced it: this string
/// is drawn INTO the exported PNG, so under `NoCorrection` it must not carry a
/// correction the figure did not have.
fn empty_placeholder_text(
    entry_label: &str,
    method: crate::dam::fdr::FdrMethod,
    fdr_threshold: f64,
) -> String {
    let metric = method.metric_label();
    format!("No {entry_label}s passed {metric} < {fdr_threshold}")
}

/// Y-axis title: sentence case, matches the X-axis "Enrichment ratio"
/// register. Prefixed with "KEGG" so a standalone figure makes the data
/// source obvious without the caption.
fn y_axis_title(entry_label: &str) -> String {
    format!("KEGG {entry_label}")
}

/// Hits-size-legend title. Embeds the per-figure `max_hits` so readers
/// comparing two figures with different max values can see at a glance
/// that the dot-size scale is per-plot (not absolute) and adjust their
/// visual comparison accordingly. Falls back to plain "Hits" when there
/// is no data (empty-state render).
fn hits_title(max_hits: Option<usize>) -> String {
    match max_hits {
        Some(m) => format!("Hits (max={m})"),
        None => "Hits".to_string(),
    }
}

/// Colorbar title. For `NoCorrection` the bar encodes raw p-values, not
/// adjusted q-values, so the title drops the `FDR (…)` wrapper to stay
/// scientifically honest. For BH / BY the wrapper names the method.
fn colorbar_title(method: crate::dam::fdr::FdrMethod) -> String {
    match method {
        crate::dam::fdr::FdrMethod::NoCorrection => "-log10(p-value)".to_string(),
        _ => format!("-log10(FDR ({}))", method.short_label()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::enrichment::types::{EnrichmentDirection, EnrichmentResult, EnrichmentRow};

    fn sample_row(id: &str, name: &str, hits: usize, ratio: f64, fdr: f64) -> EnrichmentRow {
        EnrichmentRow {
            entry_id: id.into(),
            entry_name: name.into(),
            hits,
            total: 10,
            expected: 1.0,
            enrichment_ratio: ratio,
            p_value: fdr / 2.0,
            fdr,
            hit_kegg_ids: vec![],
        }
    }

    fn sample_result(n: usize) -> EnrichmentResult {
        let mut rows = Vec::new();
        for i in 0..n {
            rows.push(sample_row(
                &format!("p{i:02}"),
                &format!("Pathway {i:02}"),
                (i + 1).min(10),
                2.0 + (i as f64) * 0.1,
                0.001 + (i as f64) * 0.005,
            ));
        }
        EnrichmentResult {
            universe_size: 100,
            dam_cpd_size: 20,
            direction: EnrichmentDirection::Both,
            min_hit_count: 1,
            min_entry_size: 1,
            entries_dropped_by_min_entry_size: 0,
            empty_compound_count: 0,
            rows,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        }
    }

    fn result_from(rows: Vec<EnrichmentRow>) -> EnrichmentResult {
        EnrichmentResult {
            universe_size: 100,
            dam_cpd_size: 20,
            direction: EnrichmentDirection::Both,
            min_hit_count: 1,
            min_entry_size: 1,
            entries_dropped_by_min_entry_size: 0,
            empty_compound_count: 0,
            rows,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiHochberg,
        }
    }

    #[test]
    fn select_and_order_rows_orders_by_fold_enrichment_descending() {
        // Rows arrive FDR-ascending (as ORA emits); display must reorder so the
        // largest fold enrichment is first (top of the y axis).
        let result = result_from(vec![
            sample_row("p0", "P0", 3, 2.0, 0.001),
            sample_row("p1", "P1", 3, 5.0, 0.002),
            sample_row("p2", "P2", 3, 3.0, 0.003),
        ]);
        let ids: Vec<&str> = select_and_order_rows(&result, 0.05, 1, 10)
            .iter()
            .map(|r| r.entry_id.as_str())
            .collect();
        assert_eq!(ids, vec!["p1", "p2", "p0"]); // 5.0 > 3.0 > 2.0
    }

    #[test]
    fn select_and_order_rows_breaks_ties_by_fdr_then_entry_id() {
        // Equal fold enrichment → ascending FDR, then ascending entry_id.
        let result = result_from(vec![
            sample_row("pA", "PA", 3, 4.0, 0.001),
            sample_row("pZ", "PZ", 3, 4.0, 0.001), // same ratio + fdr → entry_id breaks
            sample_row("pM", "PM", 3, 4.0, 0.002), // same ratio, higher fdr → sorts last
        ]);
        let ids: Vec<&str> = select_and_order_rows(&result, 0.05, 1, 10)
            .iter()
            .map(|r| r.entry_id.as_str())
            .collect();
        assert_eq!(ids, vec!["pA", "pZ", "pM"]);
    }

    #[test]
    fn select_and_order_rows_keeps_most_significant_top_n_then_orders_by_fc() {
        // 5 significant rows, FDR-ascending; the two highest-FC entries are also
        // the LEAST significant. top_n = 3 must keep the 3 lowest-FDR (p0,p1,p2)
        // — dropping the high-FC p3,p4 — then order the kept 3 by FC descending.
        let result = result_from(vec![
            sample_row("p0", "P0", 3, 2.0, 0.001),
            sample_row("p1", "P1", 3, 2.2, 0.002),
            sample_row("p2", "P2", 3, 2.4, 0.003),
            sample_row("p3", "P3", 3, 9.0, 0.004),
            sample_row("p4", "P4", 3, 9.5, 0.005),
        ]);
        let ids: Vec<&str> = select_and_order_rows(&result, 0.05, 1, 3)
            .iter()
            .map(|r| r.entry_id.as_str())
            .collect();
        // Selection by significance (p3/p4 excluded), display by FC descending.
        assert_eq!(ids, vec!["p2", "p1", "p0"]);
    }

    #[test]
    fn render_returns_rgba_buffer_of_expected_size() {
        let result = sample_result(5);
        let opts = DotplotOpts {
            width_px: 400,
            height_px: 400,
            fdr_threshold: 0.05,
            min_hit_count: 1,
            top_n: 10,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            entry_label: "pathway",
        };
        let buf = render_dotplot(&result, &opts).expect("render");
        assert_eq!(buf.len(), 400 * 400 * 4);
        // Smoke: not entirely white (we drew something).
        let some_non_white = buf
            .chunks(4)
            .any(|p| p[0] != 255 || p[1] != 255 || p[2] != 255);
        assert!(some_non_white);
    }

    #[test]
    fn empty_result_renders_placeholder() {
        let result = EnrichmentResult {
            universe_size: 100,
            dam_cpd_size: 5,
            direction: EnrichmentDirection::Both,
            min_hit_count: 1,
            min_entry_size: 1,
            entries_dropped_by_min_entry_size: 0,
            empty_compound_count: 0,
            rows: vec![],
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let opts = DotplotOpts {
            width_px: 400,
            height_px: 300,
            fdr_threshold: 0.05,
            min_hit_count: 1,
            top_n: 10,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            entry_label: "pathway",
        };
        let buf = render_dotplot(&result, &opts).expect("empty render");
        assert_eq!(buf.len(), 400 * 300 * 4);
    }

    // ---- annotation strip: plain-language, multi-line block ----
    // build_strip now returns the lines (Vec<String>) and takes the entry
    // label so the wording is mode-aware. The N / K / m symbols are gone;
    // their numeric values stay so the BH/BY divisor remains auditable.

    #[test]
    fn build_strip_default_is_four_plain_language_lines() {
        // sample_result: universe=100, dam_cpd=20, direction Both,
        // min_entry_size=1, dropped=0, min_hit=1, method BY.
        let r = sample_result(5);
        let lines = build_strip(&r, 1, "pathway");
        assert_eq!(lines.len(), 4, "got: {lines:?}");
        assert_eq!(
            lines[0],
            "Background universe = 100 compounds measured and mapped to KEGG"
        );
        assert_eq!(
            lines[1],
            "Compounds of interest = 20 differentially abundant (both directions)"
        );
        // Clean run: no fraction, no skipped/needs, no hits clause.
        assert_eq!(lines[2], "Pathways tested = 5");
        assert_eq!(
            lines[3],
            "Significance: FDR-adjusted, Benjamini–Yekutieli (BY)"
        );
    }

    #[test]
    fn build_strip_drops_the_n_k_m_symbols() {
        let joined = build_strip(&sample_result(5), 1, "pathway").join("\n");
        assert!(!joined.contains("N="), "{joined}");
        assert!(!joined.contains("K="), "{joined}");
        assert!(!joined.contains("m:"), "{joined}");
        assert!(!joined.contains("m="), "{joined}");
    }

    #[test]
    fn build_strip_no_correction_significance_is_self_explanatory() {
        let mut r = sample_result(3);
        r.fdr_method = crate::dam::fdr::FdrMethod::NoCorrection;
        let lines = build_strip(&r, 1, "pathway");
        assert_eq!(lines[3], "Significance: raw p-value (no FDR correction)");
        // The cryptic bare method token must be gone in EITHER spelling — the
        // pre-rename `FDR: None` and the post-rename `FDR: NoCorrection`. The
        // strip builds its wording locally and never calls `short_label()`, and
        // this is what keeps that true.
        let joined = lines.join("\n");
        assert!(!joined.contains("FDR: None"));
        assert!(!joined.contains("FDR: NoCorrection"));
    }

    #[test]
    fn build_strip_bh_names_benjamini_hochberg_with_en_dash() {
        let mut r = sample_result(3);
        r.fdr_method = crate::dam::fdr::FdrMethod::BenjaminiHochberg;
        // En-dash U+2013, copied from the spec — an ASCII hyphen would
        // silently never match.
        assert_eq!(
            build_strip(&r, 1, "pathway")[3],
            "Significance: FDR-adjusted, Benjamini–Hochberg (BH)"
        );
    }

    #[test]
    fn build_strip_direction_phrases_are_plain_language() {
        let mut up = sample_result(1);
        up.direction = EnrichmentDirection::Up;
        assert!(build_strip(&up, 1, "pathway")[1].ends_with("(increased)"));

        let mut down = sample_result(1);
        down.direction = EnrichmentDirection::Down;
        assert!(build_strip(&down, 1, "pathway")[1].ends_with("(decreased)"));

        let both = sample_result(1);
        assert!(build_strip(&both, 1, "pathway")[1].ends_with("(both directions)"));
    }

    #[test]
    fn build_strip_module_mode_uses_module_entity_word() {
        let lines = build_strip(&sample_result(4), 1, "module");
        assert!(
            lines[2].starts_with("Modules tested ="),
            "got: {}",
            lines[2]
        );
    }

    #[test]
    fn build_strip_shows_fraction_and_skipped_when_min_entry_dropped() {
        let mut r = sample_result(5);
        r.min_entry_size = 3;
        r.entries_dropped_by_min_entry_size = 7; // 5 tested + 7 dropped = 12
        assert_eq!(
            build_strip(&r, 1, "pathway")[2],
            "Pathways tested = 5 of 12  ·  7 skipped (< 3 compounds each)"
        );
    }

    #[test]
    fn build_strip_documents_min_entry_when_nothing_dropped() {
        let mut r = sample_result(5);
        r.min_entry_size = 3;
        r.entries_dropped_by_min_entry_size = 0;
        assert_eq!(
            build_strip(&r, 1, "pathway")[2],
            "Pathways tested = 5  ·  each pathway needs ≥ 3 compounds"
        );
    }

    #[test]
    fn build_strip_appends_min_hit_clause_last() {
        let mut r = sample_result(5);
        r.min_entry_size = 3;
        r.entries_dropped_by_min_entry_size = 0;
        // The hits clause reports the LIVE `min_hit_count` the figure was drawn
        // with, not whatever the run started from — the filter is applied at
        // draw time, so the strip must describe the same value.
        assert_eq!(
            build_strip(&r, 3, "pathway")[2],
            "Pathways tested = 5  ·  each pathway needs ≥ 3 compounds; ≥ 3 hits required"
        );
    }

    #[test]
    fn build_strip_min_hit_clause_only_when_above_one() {
        let r = sample_result(5);
        assert_eq!(
            build_strip(&r, 3, "pathway")[2],
            "Pathways tested = 5; ≥ 3 hits required"
        );
        // At 1 (the default, i.e. no filtering) the clause is omitted entirely.
        assert_eq!(build_strip(&r, 1, "pathway")[2], "Pathways tested = 5");
    }

    #[test]
    fn export_writes_phys_chunk_with_dpi() {
        let result = sample_result(3);
        let opts = DotplotOpts {
            width_px: 200,
            height_px: 200,
            fdr_threshold: 0.05,
            min_hit_count: 1,
            top_n: 10,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            entry_label: "pathway",
        };
        let tmp = tempfile::Builder::new()
            .suffix(".png")
            .tempfile()
            .expect("tempfile");
        export_dotplot_png(&result, &opts, 300, tmp.path()).expect("export");
        let bytes = std::fs::read(tmp.path()).expect("read png");

        // Locate the pHYs chunk and decode its xppu / yppu fields.
        let mut idx = 8; // skip signature
        let mut found: Option<(u32, u32, u8)> = None;
        while idx + 8 < bytes.len() {
            let len = u32::from_be_bytes(bytes[idx..idx + 4].try_into().unwrap()) as usize;
            let ty = &bytes[idx + 4..idx + 8];
            if ty == b"pHYs" {
                let payload = &bytes[idx + 8..idx + 8 + len];
                let xppu = u32::from_be_bytes(payload[0..4].try_into().unwrap());
                let yppu = u32::from_be_bytes(payload[4..8].try_into().unwrap());
                let unit = payload[8];
                found = Some((xppu, yppu, unit));
                break;
            }
            idx += 8 + len + 4;
        }
        let (xppu, yppu, unit) = found.expect("pHYs chunk not found");
        let expected = ((300.0_f64) / 0.0254).round() as u32;
        assert_eq!(xppu, expected);
        assert_eq!(yppu, expected);
        assert_eq!(unit, 1);
    }

    #[test]
    fn fdr_to_color_at_threshold_is_yellow_end() {
        // FDR == threshold → nl == nl_threshold → t == 0 → YLORRD_9[0].
        let nl_threshold = -0.05_f64.log10();
        let max_nl = nl_threshold + 4.0;
        let color = fdr_to_color(0.05, nl_threshold, max_nl);
        let (r, g, b) = common::YLORRD_9[0];
        assert_eq!(color, RGBColor(r, g, b));
    }

    #[test]
    fn fdr_to_color_at_data_max_is_red_end() {
        // FDR == 10^(-max_nl) → nl == max_nl → t == 1 → YLORRD_9[8].
        let nl_threshold = -0.05_f64.log10();
        let max_nl = nl_threshold + 4.0;
        let fdr_at_max = 10f64.powf(-max_nl);
        let color = fdr_to_color(fdr_at_max, nl_threshold, max_nl);
        let (r, g, b) = common::YLORRD_9[8];
        assert_eq!(color, RGBColor(r, g, b));
    }

    #[test]
    fn fdr_to_color_zero_fdr_clamps_to_red_without_nan() {
        // BH/BY may underflow to 0.0; renderer must not crash or produce
        // a NaN/+inf-driven garbage colour.
        let nl_threshold = -0.05_f64.log10();
        let max_nl = nl_threshold + 4.0;
        let color = fdr_to_color(0.0, nl_threshold, max_nl);
        let (r, g, b) = common::YLORRD_9[8];
        assert_eq!(color, RGBColor(r, g, b));
    }

    #[test]
    fn fdr_to_color_distinguishes_highly_significant_band() {
        // The whole point of the -log10 change: FDR=1e-6 and FDR=1e-3
        // must produce visibly different colours (in the linear-FDR
        // pre-change they were both essentially saturated red).
        let nl_threshold = -0.05_f64.log10();
        let max_nl = 6.0;
        let c_very = fdr_to_color(1e-6, nl_threshold, max_nl);
        let c_mod = fdr_to_color(1e-3, nl_threshold, max_nl);
        let dr = (c_very.0 as i32 - c_mod.0 as i32).abs();
        let dg = (c_very.1 as i32 - c_mod.1 as i32).abs();
        let db = (c_very.2 as i32 - c_mod.2 as i32).abs();
        // Sum of channel deltas exceeds a perceptual-difference floor
        // of ~30 (well above JND for sRGB midtones).
        assert!(
            dr + dg + db > 30,
            "1e-6 vs 1e-3 not visually distinguishable: Δ=({dr},{dg},{db})"
        );
    }

    #[test]
    fn neg_log10_span_empty_rows_falls_back_to_two_decades() {
        let (nl_threshold, max_nl) = neg_log10_span(&[], 0.05);
        let expected_threshold = -0.05_f64.log10();
        assert!((nl_threshold - expected_threshold).abs() < 1e-12);
        assert!((max_nl - (expected_threshold + COLORBAR_FALLBACK_SPAN)).abs() < 1e-12);
    }

    #[test]
    fn neg_log10_span_uses_data_max_when_rows_present() {
        let row = sample_row("p0", "Pathway 0", 1, 2.0, 1e-5);
        let rows = vec![&row];
        let (nl_threshold, max_nl) = neg_log10_span(&rows, 0.05);
        assert!((nl_threshold - 1.30103).abs() < 1e-4);
        assert!((max_nl - 5.0).abs() < 1e-9);
    }

    #[test]
    fn empty_placeholder_text_pluralises_entry_label() {
        assert_eq!(
            empty_placeholder_text(
                "pathway",
                crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
                0.05
            ),
            "No pathways passed FDR < 0.05"
        );
        assert_eq!(
            empty_placeholder_text(
                "module",
                crate::dam::fdr::FdrMethod::BenjaminiHochberg,
                0.05
            ),
            "No modules passed FDR < 0.05"
        );
        // The metric noun follows the method — this string is drawn into the
        // exported PNG, so an uncorrected run must not claim an FDR.
        assert_eq!(
            empty_placeholder_text("pathway", crate::dam::fdr::FdrMethod::NoCorrection, 0.05),
            "No pathways passed p-value < 0.05"
        );
    }

    #[test]
    fn y_axis_title_prefixes_kegg() {
        assert_eq!(y_axis_title("pathway"), "KEGG pathway");
        assert_eq!(y_axis_title("module"), "KEGG module");
    }

    #[test]
    fn hits_title_embeds_per_figure_max_or_falls_back() {
        assert_eq!(hits_title(Some(12)), "Hits (max=12)");
        assert_eq!(hits_title(Some(1)), "Hits (max=1)");
        // Empty-state render has no data → no max to display.
        assert_eq!(hits_title(None), "Hits");
    }

    #[test]
    fn colorbar_title_uses_p_value_for_no_correction_else_fdr_method() {
        use crate::dam::fdr::FdrMethod;
        assert_eq!(
            colorbar_title(FdrMethod::BenjaminiHochberg),
            "-log10(FDR (BH))"
        );
        assert_eq!(
            colorbar_title(FdrMethod::BenjaminiYekutieli),
            "-log10(FDR (BY))"
        );
        // NoCorrection drops the FDR wrapper — the bar encodes raw p, not q.
        assert_eq!(colorbar_title(FdrMethod::NoCorrection), "-log10(p-value)");
    }

    #[test]
    fn analysis_mode_entry_label_singular_maps_both_variants() {
        use crate::app::AnalysisMode;
        assert_eq!(AnalysisMode::Pathway.entry_label_singular(), "pathway");
        assert_eq!(AnalysisMode::Module.entry_label_singular(), "module");
    }

    #[test]
    fn colorbar_ticks_span_1_30_to_7_includes_endpoints_and_fitting_integers() {
        // 180-px bar over span 5.70 → 31.6 px / unit. Integer "2" sits
        // 22 px from the 1.30 endpoint, below the 28-px gap, so it's
        // filtered (the endpoint already labels the bottom region).
        // Remaining integers 3-6 fit comfortably between 1.30 and 7.
        let ticks = colorbar_tick_positions(-0.05_f64.log10(), 7.0, 180, 28);
        let expected = [-0.05_f64.log10(), 3.0, 4.0, 5.0, 6.0, 7.0];
        assert_eq!(ticks.len(), expected.len());
        for (got, want) in ticks.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-9, "{got} vs {want}");
        }
    }

    #[test]
    fn colorbar_ticks_short_bar_thins_interior_integers() {
        // 60-px bar with 28-px gap = ~2 ticks fit. Endpoints survive,
        // interior integers are dropped greedily.
        let ticks = colorbar_tick_positions(-0.05_f64.log10(), 7.0, 60, 28);
        assert!(ticks.len() <= 3, "got {} ticks, expected ≤3", ticks.len());
        // Endpoints always present.
        assert!((ticks.first().unwrap() - (-0.05_f64.log10())).abs() < 1e-9);
        assert!((ticks.last().unwrap() - 7.0).abs() < 1e-9);
    }

    #[test]
    fn colorbar_ticks_endpoint_matches_integer_no_duplicate() {
        // span [2.0, 6.0]: endpoints are integers themselves; interior
        // integers 3,4,5 fill in — total 5 ticks all integer, no dupes.
        let ticks = colorbar_tick_positions(2.0, 6.0, 180, 28);
        let expected = [2.0, 3.0, 4.0, 5.0, 6.0];
        assert_eq!(ticks.len(), expected.len());
        for (got, want) in ticks.iter().zip(expected.iter()) {
            assert!((got - want).abs() < 1e-9);
        }
    }

    #[test]
    fn colorbar_ticks_tiny_span_returns_only_endpoints() {
        // span [1.30, 1.50] — no interior integer exists. Only endpoints.
        let ticks = colorbar_tick_positions(1.30, 1.50, 180, 28);
        assert_eq!(ticks.len(), 2);
        assert!((ticks[0] - 1.30).abs() < 1e-9);
        assert!((ticks[1] - 1.50).abs() < 1e-9);
    }

    #[test]
    fn format_tick_integers_drop_trailing_decimals_others_keep_two() {
        assert_eq!(format_tick(2.0), "2");
        assert_eq!(format_tick(7.0), "7");
        assert_eq!(format_tick(12.0), "12");
        assert_eq!(format_tick(1.30103), "1.30");
        assert_eq!(format_tick(0.0), "0");
    }
}
