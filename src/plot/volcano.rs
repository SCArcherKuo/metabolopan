//! Volcano plot renderer. Single source of truth used for both the in-window preview
//! (typically 800×800) and the 600 DPI PNG export (6000×6000).
//!
//! theme: scope-excluded. UP/DOWN/NS = biological convention, not brand palette.

use anyhow::{Result, anyhow};
use plotters::prelude::*;
use rand::Rng;
use rand::SeedableRng;
use rand::rngs::StdRng;
use std::path::Path;

use crate::dam::run::classify_trend;
use crate::dam::types::{DamMethod, DamResult, FcBasis, Trend};
use crate::plot::common;

/// Y-axis safety cap (well above the −log10 of f64::MIN_POSITIVE ≈ 308).
const Y_ABSOLUTE_CAP: f64 = 320.0;

/// Margin (in log2-FC units) added to the symmetric X range so ±∞ dock points sit at
/// the visible edge without clipping.
const X_PADDING: f64 = 0.5;

/// Per-point jitter magnitude (±) applied to ±∞ X-axis dock positions
/// (in log2-FC units).
const INF_JITTER: f64 = 0.04;

/// Per-point downward jitter magnitude applied to +∞ Y-axis dock positions
/// (in −log10(q) units). Scale-matched to `INF_JITTER`: the X axis pads by
/// `X_PADDING = 0.5` with `±0.04` jitter (≈ 8 % of the padding); the Y axis
/// pads by `+1.0` (`y_max = y_finite_max + 1.0`), so 8 % of that is `0.08`.
/// Downward-only (range `[0, Y_INF_JITTER]` subtracted from `y_max`) so the
/// jittered points stay strictly within the chart range `[0, y_max]` — the
/// X-axis symmetric ±jitter pattern would push half the points above
/// `y_max` and get them clipped by plotters.
const Y_INF_JITTER: f64 = 0.08;

/// Fixed seed for the ±∞ dock jitter RNG. Using a seeded `StdRng` instead of
/// `thread_rng()` makes the jitter deterministic: the same `DamResult` renders
/// to a byte-identical buffer every time (reproducible PNG exports + preview),
/// and the `y_inf_jitter_spreads_…` regression test is no longer flaky. The
/// jitter still spreads multiple q-saturated / ±∞ features across distinct
/// pixels — a fixed seed yields a fixed-but-well-distributed sequence, which is
/// all the anti-collision behaviour needs.
const JITTER_SEED: u64 = 42; // arbitrary fixed value; chosen so the deterministic
// jitter sequence spreads multi-feature ±∞ docks across distinct pixel rows
// (see `y_inf_jitter_spreads_q_saturated_features_across_rows`).

/// Scatter point colours (also used by the legend with full opacity).
const UP_FILL: RGBAColor = RGBAColor(220, 70, 70, 0.5);
const DOWN_FILL: RGBAColor = RGBAColor(70, 110, 220, 0.5);
const NS_FILL: RGBAColor = RGBAColor(155, 155, 165, 0.5);
const UP_SOLID: RGBColor = RGBColor(220, 70, 70);
const DOWN_SOLID: RGBColor = RGBColor(70, 110, 220);
const NS_SOLID: RGBColor = RGBColor(180, 180, 190);

#[derive(Debug, Clone, Copy, Default)]
struct TrendCounts {
    up: usize,
    down: usize,
    ns: usize,
}

#[derive(Debug, Clone, Copy)]
pub struct VolcanoOpts {
    pub width_px: u32,
    pub height_px: u32,
    pub fc_threshold: f64,
    pub fdr_threshold: f64,
    pub delta_threshold: f64,
    /// FDR method that produced `result.p_adjusted`. Affects only the
    /// annotation-strip label (`FDR(BH)` vs `FDR(BY)`); does NOT alter
    /// marker geometry, threshold lines, or any other pixel.
    pub fdr_method: crate::dam::fdr::FdrMethod,
}

/// Render the volcano plot into an RGBA buffer of `width_px * height_px * 4` bytes.
///
/// Internally plotters writes 3 bytes per pixel (RGB) to a scratch buffer; we then
/// expand to RGBA so callers (egui, image PNG export) can use it directly. The right
/// edge of the canvas is reserved for a trend legend (counts of Up / Down / ns).
pub fn render_volcano(result: &DamResult, opts: &VolcanoOpts) -> Result<Vec<u8>> {
    let w = opts.width_px;
    let h = opts.height_px;
    let pixel_count = (w as usize) * (h as usize);
    let mut rgb = vec![0u8; pixel_count * 3];

    // Shared scale (used by both chart and legend so fonts/dots match visually).
    let scale = common::design_scale(w, h);
    // Legend takes a fixed scaled width but no more than 30 % of the canvas.
    let legend_w = ((220.0 * scale).round() as u32).min(w * 3 / 10).max(120);
    let chart_w_boundary = (w - legend_w) as i32;

    {
        let root = BitMapBackend::with_buffer(&mut rgb, (w, h)).into_drawing_area();
        root.fill(&WHITE)
            .map_err(|e| anyhow!("fill background: {e:?}"))?;

        let (chart_area, legend_area) = root.split_horizontally(chart_w_boundary);
        let counts = draw_chart(&chart_area, result, opts, scale)?;
        draw_legend(&legend_area, counts, scale, result.method)?;
        root.present().map_err(|e| anyhow!("present: {e:?}"))?;
    }

    Ok(common::rgb_to_rgba(&rgb, pixel_count))
}

/// Render the volcano plot and save as PNG, embedding the requested `dpi` as a
/// pHYs chunk so the file declares its physical pixel density (pixels-per-meter).
/// Layout tools such as Word and InDesign read this chunk to set the on-page size;
/// callers that don't care can pass any dpi (e.g. 72) — the pixel data is unaffected.
pub fn export_volcano_png(
    result: &DamResult,
    opts: &VolcanoOpts,
    dpi: u32,
    out: &Path,
) -> Result<()> {
    let buffer = render_volcano(result, opts)?;
    common::encode_png(&buffer, opts.width_px, opts.height_px, dpi, out)
}

fn draw_chart<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    result: &DamResult,
    opts: &VolcanoOpts,
    scale: f64,
) -> Result<TrendCounts>
where
    <DB as DrawingBackend>::ErrorType: 'static,
{
    let (xabs_max, y_finite_max, n_pos_inf, n_neg_inf) = scan_ranges(result);
    let x_edge = xabs_max + X_PADDING;
    let y_max = (y_finite_max + 1.0).min(Y_ABSOLUTE_CAP);

    // Every pixel-denominated constant is scaled relative to an 800-px design baseline
    // so the chart stays legible at any export size (e.g. 6000×6000 at 600 DPI).
    let sp = |v: f64| common::sp(v, scale);
    let su = |v: f64| common::su(v, scale);

    // Reserve the bottom of the chart area for the annotation strip so the chart's
    // X label area (tick numbers + axis title) sits above it with no overlap. The
    // boundary's height ≈ two title-character heights, which is how far the X axis
    // title visually moves up compared with drawing the strip on top of the area.
    let (ah_full, _) = (area.dim_in_pixel().1, ());
    let strip_h = su(116.0) as i32;
    let split_y = (ah_full as i32 - strip_h).max(1);
    let (chart_main, strip_area) = area.split_vertically(split_y);

    let mut chart = ChartBuilder::on(&chart_main)
        .margin(sp(5.0))
        .x_label_area_size(su(80.0))
        .y_label_area_size(su(90.0))
        .build_cartesian_2d(-x_edge..x_edge, 0f64..y_max)
        .map_err(|e| anyhow!("chart build: {e:?}"))?;

    chart
        .configure_mesh()
        .x_desc("Log2(Fold Change)")
        .y_desc("-Log10(FDR)")
        .label_style(("sans-serif", su(21.0)).into_font().color(&BLACK))
        .axis_desc_style(("sans-serif", su(36.0)).into_font().color(&BLACK))
        // Hide minor gridlines (light) by making them transparent; keep major
        // gridlines (bold) as a soft gray on the white background.
        .light_line_style(RGBAColor(0, 0, 0, 0.0))
        .bold_line_style(common::GRID_LIGHT)
        .axis_style(BLACK)
        .draw()
        .map_err(|e| anyhow!("mesh draw: {e:?}"))?;

    // Threshold lines first (below points).
    let log2_fc = opts.fc_threshold.log2();
    let neg_log10_fdr = -opts.fdr_threshold.log10();
    let line_stroke = su(1.0);
    chart
        .draw_series(std::iter::once(PathElement::new(
            vec![(-log2_fc, 0.0), (-log2_fc, y_max)],
            ShapeStyle::from(&BLACK).stroke_width(line_stroke),
        )))
        .map_err(|e| anyhow!("vline left: {e:?}"))?;
    chart
        .draw_series(std::iter::once(PathElement::new(
            vec![(log2_fc, 0.0), (log2_fc, y_max)],
            ShapeStyle::from(&BLACK).stroke_width(line_stroke),
        )))
        .map_err(|e| anyhow!("vline right: {e:?}"))?;
    chart
        .draw_series(std::iter::once(PathElement::new(
            vec![(-x_edge, neg_log10_fdr), (x_edge, neg_log10_fdr)],
            ShapeStyle::from(&BLACK).stroke_width(line_stroke),
        )))
        .map_err(|e| anyhow!("hline: {e:?}"))?;

    // Scatter, grouped by trend so colours layer predictably (ns → down → up).
    // Tuple = (x, y, scaled-pixel-radius); BM features carry per-feature radii
    // mapped from |Cliff's δ|, every other path keeps the uniform su(3.0).
    let mut rng = StdRng::seed_from_u64(JITTER_SEED);
    let mut ns_points: Vec<(f64, f64, u32)> = Vec::new();
    let mut up_points: Vec<(f64, f64, u32)> = Vec::new();
    let mut down_points: Vec<(f64, f64, u32)> = Vec::new();
    let mut out_of_range_count: usize = 0;

    for feat in &result.features {
        let trend = classify_trend(
            feat,
            opts.fc_threshold,
            opts.fdr_threshold,
            opts.delta_threshold,
            result.method,
        );
        let x = if feat.log2_fold_change == f64::INFINITY {
            x_edge + rng.gen_range(-INF_JITTER..=INF_JITTER)
        } else if feat.log2_fold_change == f64::NEG_INFINITY {
            -x_edge + rng.gen_range(-INF_JITTER..=INF_JITTER)
        } else if feat.log2_fold_change.is_nan() {
            continue;
        } else {
            feat.log2_fold_change
        };
        // neg_log10_p_adjusted contract (post-2026-05-26):
        //   - NaN  ⇒ genuine "p couldn't be computed" (BM stratified, Welch
        //              n<2). Drop the point from the plot.
        //   - +INF ⇒ q underflow saturation (q == 0 in BH/BY). Dock just
        //              below y_max with per-point downward jitter (in
        //              [0, Y_INF_JITTER]) so multiple saturated features
        //              don't pile at a single pixel — mirrors the ±INF
        //              X-axis jitter convention. The legend tally still
        //              counts the feature.
        //   - finite ⇒ clamp to y_max.
        let y = if feat.neg_log10_p_adjusted.is_nan() {
            continue;
        } else if feat.neg_log10_p_adjusted == f64::INFINITY {
            y_max - rng.gen_range(0.0..=Y_INF_JITTER)
        } else {
            feat.neg_log10_p_adjusted.min(y_max)
        };
        // BM with a valid Cliff's δ maps |δ| ∈ [0, 1] to baseline radius units
        // in [1.5, 3.9] (linear). Other paths keep su(3.0). See design D1/D2/D7.
        let r = match result.method {
            DamMethod::BrunnerMunzel => match feat.effect_size {
                Some(d) => {
                    debug_assert!(d.abs() <= 1.0 + 1e-12, "Cliff's δ out of range: {d}");
                    if d.abs() > 1.0 + 1e-12 {
                        out_of_range_count += 1;
                    }
                    let d_abs = d.abs().clamp(0.0, 1.0);
                    su(1.5 + d_abs * 2.4)
                }
                None => su(3.0),
            },
            DamMethod::Welch | DamMethod::Student => su(3.0),
        };
        match trend {
            Trend::Up => up_points.push((x, y, r)),
            Trend::Down => down_points.push((x, y, r)),
            Trend::NotSignificant => ns_points.push((x, y, r)),
        }
    }

    if out_of_range_count > 0 {
        tracing::error!(
            count = out_of_range_count,
            "BM volcano: features with |δ| > 1 clamped to 1.0 during size mapping"
        );
    }

    let counts = TrendCounts {
        up: up_points.len(),
        down: down_points.len(),
        ns: ns_points.len(),
    };
    chart
        .draw_series(
            ns_points
                .into_iter()
                .map(|(x, y, r)| Circle::new((x, y), r, NS_FILL.filled())),
        )
        .map_err(|e| anyhow!("ns scatter: {e:?}"))?;
    chart
        .draw_series(
            down_points
                .into_iter()
                .map(|(x, y, r)| Circle::new((x, y), r, DOWN_FILL.filled())),
        )
        .map_err(|e| anyhow!("down scatter: {e:?}"))?;
    chart
        .draw_series(
            up_points
                .into_iter()
                .map(|(x, y, r)| Circle::new((x, y), r, UP_FILL.filled())),
        )
        .map_err(|e| anyhow!("up scatter: {e:?}"))?;

    // Annotation strip — drawn into the reserved strip_area below the chart so it
    // never overlaps the X axis tick labels or title.
    let strip = build_strip(result, opts, n_pos_inf, n_neg_inf);
    strip_area
        .draw_text(
            &strip,
            &("sans-serif", su(28.0)).into_font().color(&BLACK),
            (sp(10.0), sp(38.0)),
        )
        .map_err(|e| anyhow!("annotation draw: {e:?}"))?;

    Ok(counts)
}

/// Draw the trend legend (Up / Down / ns counts) on a dedicated drawing area.
/// Title at top, then three rows with a coloured swatch and `Label: N` text.
/// On BM renders, a second `|δ| size` section is appended below with three
/// reference dots at |δ| = 0.0 / 0.5 / 1.0 sized through the same r_units
/// mapping as the scatter (design D4). Welch / Student renders skip the size
/// section so the legend stays pixel-identical to the pre-change output.
fn draw_legend<DB: DrawingBackend>(
    area: &DrawingArea<DB, plotters::coord::Shift>,
    counts: TrendCounts,
    scale: f64,
    method: DamMethod,
) -> Result<()>
where
    <DB as DrawingBackend>::ErrorType: 'static,
{
    let sp = |v: f64| common::sp(v, scale);
    let su = |v: f64| common::su(v, scale);

    let title_font = ("sans-serif", su(32.0)).into_font().color(&BLACK);
    let row_font = ("sans-serif", su(28.0)).into_font().color(&BLACK);

    area.draw_text("DAM trend", &title_font, (sp(12.0), sp(40.0)))
        .map_err(|e| anyhow!("legend title: {e:?}"))?;

    let dot_x = sp(24.0);
    let text_x = sp(54.0);
    let dot_radius = su(5.0) as i32;
    let row_height = sp(56.0);
    let first_row_y = sp(110.0);
    let text_baseline_adjust = sp(14.0);

    let rows: [(&str, RGBColor, usize); 3] = [
        ("Up", UP_SOLID, counts.up),
        ("Down", DOWN_SOLID, counts.down),
        ("ns", NS_SOLID, counts.ns),
    ];
    for (i, (label, colour, count)) in rows.iter().enumerate() {
        let y = first_row_y + (i as i32) * row_height;
        area.draw(&Circle::new((dot_x, y), dot_radius, colour.filled()))
            .map_err(|e| anyhow!("legend dot: {e:?}"))?;
        area.draw_text(
            &format!("{label}: {count}"),
            &row_font,
            (text_x, y - text_baseline_adjust),
        )
        .map_err(|e| anyhow!("legend row: {e:?}"))?;
    }

    if method == DamMethod::BrunnerMunzel {
        area.draw_text("|δ| size", &title_font, (sp(12.0), sp(290.0)))
            .map_err(|e| anyhow!("size legend title: {e:?}"))?;

        let size_first_row_y = sp(350.0);
        let size_dot_colour = RGBColor(120, 120, 120);
        let size_rows: [f64; 3] = [0.0, 0.5, 1.0];
        for (i, d_abs) in size_rows.iter().enumerate() {
            let y = size_first_row_y + (i as i32) * row_height;
            let r = su(1.5 + d_abs * 2.4) as i32;
            area.draw(&Circle::new((dot_x, y), r, size_dot_colour.filled()))
                .map_err(|e| anyhow!("size legend dot: {e:?}"))?;
            area.draw_text(
                &format!("|δ|: {d_abs:.1}"),
                &row_font,
                (text_x, y - text_baseline_adjust),
            )
            .map_err(|e| anyhow!("size legend row: {e:?}"))?;
        }
    }
    Ok(())
}

fn scan_ranges(result: &DamResult) -> (f64, f64, usize, usize) {
    let mut xabs_max: f64 = 0.0;
    let mut y_finite_max: f64 = 0.0;
    let mut pos_inf: usize = 0;
    let mut neg_inf: usize = 0;
    for feat in &result.features {
        let lf = feat.log2_fold_change;
        if lf == f64::INFINITY {
            pos_inf += 1;
        } else if lf == f64::NEG_INFINITY {
            neg_inf += 1;
        } else if lf.is_finite() {
            let a = lf.abs();
            if a > xabs_max {
                xabs_max = a;
            }
        }
        let y = feat.neg_log10_p_adjusted;
        if y.is_finite() && y > y_finite_max {
            y_finite_max = y;
        }
    }
    if xabs_max == 0.0 {
        xabs_max = 2.0; // sensible default for an empty / no-finite-fc result
    }
    if y_finite_max == 0.0 {
        y_finite_max = 1.0;
    }
    (xabs_max, y_finite_max, pos_inf, neg_inf)
}

fn build_strip(
    result: &DamResult,
    opts: &VolcanoOpts,
    n_pos_inf: usize,
    n_neg_inf: usize,
) -> String {
    let method = result.method.display_name();
    let fdr_label = opts.fdr_method.short_label();
    // FC scale (mean / median / arcsinh-mean) is per-run — every feature
    // in a run shares the same `FcBasis`. Sample from the first feature.
    // For an empty result (no feature passed the prefilter / dedup / Unknown
    // filters) the strip still renders; fall back to the method-natural
    // basis so the field doesn't read `?`.
    let fc_basis = result
        .features
        .first()
        .map(|f| f.fc_basis)
        .unwrap_or(match result.method {
            DamMethod::BrunnerMunzel => FcBasis::Median,
            DamMethod::Welch | DamMethod::Student => FcBasis::Mean,
        });
    let basis_label = fc_basis.label();
    let threshold_part = match result.method {
        DamMethod::Welch | DamMethod::Student => format!(
            "FDR({fdr_label})<{fdr}, FC≥{fc}",
            fdr = opts.fdr_threshold,
            fc = opts.fc_threshold
        ),
        DamMethod::BrunnerMunzel => format!(
            "FDR({fdr_label})<{fdr}, FC≥{fc}, |δ|≥{delta}",
            fdr = opts.fdr_threshold,
            fc = opts.fc_threshold,
            delta = opts.delta_threshold
        ),
    };
    format!(
        "Method: {method} | FC: {basis_label} | {threshold_part} | −∞: {n_neg_inf}  +∞: {n_pos_inf}"
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::dam::types::{DamFeature, FcBasis};

    fn feat(log2_fc: f64, neg_log10: f64) -> DamFeature {
        DamFeature {
            alignment_id: "f".into(),
            metabolite_name: "f".into(),
            inchikey: Some("X".into()),
            average_rt_min: None,
            average_mz: None,
            formula: None,
            smiles: None,
            numerator_mean: 0.0,
            denominator_mean: 0.0,
            numerator_median: 0.0,
            denominator_median: 0.0,
            fold_change: 0.0,
            log2_fold_change: log2_fc,
            fc_basis: FcBasis::Mean,
            p_value: 0.01,
            p_adjusted: 0.01,
            neg_log10_p_adjusted: neg_log10,
            effect_size: None,
        }
    }

    fn small_result() -> DamResult {
        DamResult {
            method: DamMethod::Welch,
            numerator: "A".into(),
            denominator: "B".into(),
            features: vec![
                feat(2.0, 4.0),
                feat(-2.0, 4.0),
                feat(0.1, 1.0),
                feat(f64::INFINITY, 3.0),
                feat(f64::NEG_INFINITY, 3.0),
            ],
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            dedup_report: None,
        }
    }

    #[test]
    fn render_returns_rgba_buffer_of_expected_size() {
        let result = small_result();
        let opts = VolcanoOpts {
            width_px: 400,
            height_px: 300,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let buf = render_volcano(&result, &opts).expect("render");
        assert_eq!(buf.len(), 400 * 300 * 4);
        // Smoke check: not all zeros (something drew).
        let any_nonzero = buf.iter().any(|&b| b != 0);
        assert!(any_nonzero, "buffer is all zeros — nothing rendered");
    }

    #[test]
    fn render_handles_empty_result() {
        let result = DamResult {
            method: DamMethod::Welch,
            numerator: "A".into(),
            denominator: "B".into(),
            features: vec![],
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            dedup_report: None,
        };
        let opts = VolcanoOpts {
            width_px: 200,
            height_px: 200,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let buf = render_volcano(&result, &opts).expect("empty render");
        assert_eq!(buf.len(), 200 * 200 * 4);
    }

    #[test]
    fn build_strip_omits_delta_for_welch() {
        let r = small_result();
        let opts = VolcanoOpts {
            width_px: 1,
            height_px: 1,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let s = build_strip(&r, &opts, 2, 1);
        assert!(s.contains("Welch's t-test"));
        assert!(s.contains("FDR(BY)<0.05"));
        assert!(s.contains("FC≥2"));
        assert!(!s.contains("|δ|"));
        assert!(s.contains("−∞: 1"));
        assert!(s.contains("+∞: 2"));
    }

    #[test]
    fn build_strip_includes_delta_for_bm() {
        let mut r = small_result();
        r.method = DamMethod::BrunnerMunzel;
        let opts = VolcanoOpts {
            width_px: 1,
            height_px: 1,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let s = build_strip(&r, &opts, 0, 0);
        assert!(s.contains("Brunner-Munzel"));
        assert!(s.contains("|δ|≥0.33"));
    }

    #[test]
    fn build_strip_no_delta_for_student() {
        let mut r = small_result();
        r.method = DamMethod::Student;
        let opts = VolcanoOpts {
            width_px: 1,
            height_px: 1,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let s = build_strip(&r, &opts, 0, 0);
        assert!(s.contains("Student's t-test"));
        assert!(s.contains("FDR(BY)<0.05"));
        assert!(s.contains("FC≥2"));
        assert!(!s.contains("|δ|"));
    }

    /// `log_transform=false` Welch / Student runs have `FcBasis::Mean`; the
    /// volcano strip MUST surface it so a reader knows the X axis is the raw
    /// mean ratio, not the arcsinh-scale or median-scale variant. Pins the
    /// 2026-05-29 disclosure change (#6 from /code-review xhigh).
    #[test]
    fn build_strip_discloses_fc_basis_mean_for_welch() {
        let r = small_result();
        let opts = VolcanoOpts {
            width_px: 1,
            height_px: 1,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let s = build_strip(&r, &opts, 0, 0);
        assert!(
            s.contains("FC: mean"),
            "Welch / Mean must disclose 'FC: mean'; got `{s}`"
        );
    }

    /// `log_transform=true` Welch / Student runs use the arcsinh-mean basis;
    /// the strip MUST surface the difference vs raw mean so readers know the
    /// X axis is on the variance-stabilised scale.
    #[test]
    fn build_strip_discloses_fc_basis_arcsinh_mean_for_welch_log_on() {
        let mut r = small_result();
        for f in r.features.iter_mut() {
            f.fc_basis = FcBasis::ArcsinhMean;
        }
        let opts = VolcanoOpts {
            width_px: 1,
            height_px: 1,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let s = build_strip(&r, &opts, 0, 0);
        assert!(
            s.contains("FC: arcsinh-mean"),
            "Welch / ArcsinhMean must disclose 'FC: arcsinh-mean'; got `{s}`"
        );
    }

    /// BM uses the median ratio; the strip MUST disclose it. Without this
    /// the BM X axis read identical to Welch's mean-ratio X axis.
    #[test]
    fn build_strip_discloses_fc_basis_median_for_bm() {
        let mut r = small_result();
        r.method = DamMethod::BrunnerMunzel;
        for f in r.features.iter_mut() {
            f.fc_basis = FcBasis::Median;
        }
        let opts = VolcanoOpts {
            width_px: 1,
            height_px: 1,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let s = build_strip(&r, &opts, 0, 0);
        assert!(
            s.contains("FC: median"),
            "BM / Median must disclose 'FC: median'; got `{s}`"
        );
    }

    /// Y-axis saturation jitter regression: many `+INF` features at the same
    /// finite x must occupy strictly more pixel rows than a single one would,
    /// proving the per-point downward jitter actually spreads them. Without
    /// jitter every q-saturated feature would dock at exactly the same y
    /// pixel and the cluster would compress to one circle's worth of rows.
    /// Pins the 2026-05-29 #7 fix (`Y_INF_JITTER`).
    #[test]
    fn y_inf_jitter_spreads_q_saturated_features_across_rows() {
        let opts = VolcanoOpts {
            width_px: 800,
            height_px: 800,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };

        // Anchor finite feature at y=4 so y_max ≈ 5 across both runs (same
        // axis scaling for fair pixel comparison).
        let mut r1_features = vec![feat(0.5, 4.0)];
        r1_features.push(feat(0.5, f64::INFINITY));
        let r1 = DamResult {
            method: DamMethod::Welch,
            numerator: "A".into(),
            denominator: "B".into(),
            features: r1_features,
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            dedup_report: None,
        };

        let mut r20_features = vec![feat(0.5, 4.0)];
        for _ in 0..20 {
            r20_features.push(feat(0.5, f64::INFINITY));
        }
        let r20 = DamResult {
            method: DamMethod::Welch,
            numerator: "A".into(),
            denominator: "B".into(),
            features: r20_features,
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            dedup_report: None,
        };

        let b1 = render_volcano(&r1, &opts).expect("1 INF render");
        let b20 = render_volcano(&r20, &opts).expect("20 INF render");

        // Saturated band: top ~10 % of chart area, narrow x-strip around the
        // common finite x (log2_fc = 0.5). Chart-main is roughly x ∈ [0, 580]
        // and the finite-feature anchor draws at y = 4 which is well below
        // this band, so the only content in the strip is the jittered cluster.
        let count_rows_in_band = |buf: &[u8]| -> usize {
            let img_w = 800usize;
            let mut rows = std::collections::HashSet::<usize>::new();
            for y in 0..80usize {
                for x in 200..500usize {
                    let i = (y * img_w + x) * 4;
                    if !(buf[i] == 255 && buf[i + 1] == 255 && buf[i + 2] == 255) {
                        rows.insert(y);
                        break;
                    }
                }
            }
            rows.len()
        };

        let rows_1 = count_rows_in_band(&b1);
        let rows_20 = count_rows_in_band(&b20);
        assert!(
            rows_20 > rows_1,
            "20 jittered +INF features must occupy strictly more pixel rows than 1 \
             (rows_1={rows_1}, rows_20={rows_20}). Pre-fix all 20 would have docked \
             at exactly y_max and collapsed to the same row span as 1 dot."
        );
    }

    /// Empty `DamResult` (no feature passed prefilter / Unknown / dedup) still
    /// renders a strip; the basis falls back to the method-natural default
    /// (BM → Median, Welch/Student → Mean) so the field never reads `?`.
    #[test]
    fn build_strip_fallback_basis_when_features_empty() {
        // BM empty result.
        let r_bm = DamResult {
            method: DamMethod::BrunnerMunzel,
            numerator: "A".into(),
            denominator: "B".into(),
            features: vec![],
            skipped: 0,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
            dedup_report: None,
        };
        let opts = VolcanoOpts {
            width_px: 1,
            height_px: 1,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let s_bm = build_strip(&r_bm, &opts, 0, 0);
        assert!(s_bm.contains("FC: median"), "empty BM falls back to median");

        // Welch empty result.
        let mut r_welch = r_bm.clone();
        r_welch.method = DamMethod::Welch;
        let s_welch = build_strip(&r_welch, &opts, 0, 0);
        assert!(
            s_welch.contains("FC: mean"),
            "empty Welch falls back to mean"
        );
    }

    #[test]
    fn export_writes_phys_chunk_with_dpi() {
        let result = small_result();
        let opts = VolcanoOpts {
            width_px: 200,
            height_px: 200,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        };
        let tmp = tempfile::Builder::new()
            .suffix(".png")
            .tempfile()
            .expect("tempfile");
        export_volcano_png(&result, &opts, 600, tmp.path()).expect("export");
        let bytes = std::fs::read(tmp.path()).expect("read png");

        // Locate the pHYs chunk and decode its xppu / yppu / unit fields.
        // PNG chunk layout: 4B length, 4B type, payload, 4B CRC.
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
        let expected = ((600.0_f64) / 0.0254).round() as u32;
        assert_eq!(xppu, expected, "xppu should encode 600 DPI");
        assert_eq!(yppu, expected, "yppu should encode 600 DPI");
        assert_eq!(unit, 1, "unit byte should be Meter (1)");
    }

    /// Count RGBA pixels in `[x0, x1) × [y0, y1)` whose RGB is not pure white.
    fn count_non_white_in_box(
        buf: &[u8],
        img_w: usize,
        x0: usize,
        x1: usize,
        y0: usize,
        y1: usize,
    ) -> usize {
        let mut count = 0;
        for y in y0..y1 {
            for x in x0..x1 {
                let i = (y * img_w + x) * 4;
                if !(buf[i] == 255 && buf[i + 1] == 255 && buf[i + 2] == 255) {
                    count += 1;
                }
            }
        }
        count
    }

    fn bm_opts_800() -> VolcanoOpts {
        VolcanoOpts {
            width_px: 800,
            height_px: 800,
            fc_threshold: 2.0,
            fdr_threshold: 0.05,
            delta_threshold: 0.33,
            fdr_method: crate::dam::fdr::FdrMethod::BenjaminiYekutieli,
        }
    }

    #[test]
    fn bm_feature_radius_scales_with_cliffs_delta_magnitude() {
        // Two BM results identical except for one feature's effect_size.
        // |δ|=1.0 maps to r_units=3.9 (~4 px at scale=1.0); |δ|=0.0 maps to
        // r_units=1.5 (~2 px). Both dots land at the same chart coordinates
        // (finite log2_fc avoids ±∞ jitter), so the whole-buffer non-white
        // pixel count diff reflects only the dot-size diff.
        let opts = bm_opts_800();

        let mut r0 = small_result();
        r0.method = DamMethod::BrunnerMunzel;
        let mut f0 = feat(0.5, 1.5);
        f0.effect_size = Some(0.0);
        r0.features = vec![f0];

        let mut r1 = small_result();
        r1.method = DamMethod::BrunnerMunzel;
        let mut f1 = feat(0.5, 1.5);
        f1.effect_size = Some(1.0);
        r1.features = vec![f1];

        let b0 = render_volcano(&r0, &opts).expect("render |δ|=0");
        let b1 = render_volcano(&r1, &opts).expect("render |δ|=1");

        // Bounding box: (2·su(3.9)+2) × (2·su(3.9)+2) = 10 × 10 centred on
        // the mapped pixel of (log2_fc=0.5, neg_log10=1.5). Computing exact
        // chart-to-pixel arithmetic would couple the test to plotters'
        // internal layout, so use the whole chart area — both renders draw
        // identical axes/gridlines/threshold-lines/legend so the diff in
        // non-white pixels in the chart-main region is exactly the dot diff.
        let nw0 = count_non_white_in_box(&b0, 800, 0, 580, 0, 684);
        let nw1 = count_non_white_in_box(&b1, 800, 0, 580, 0, 684);
        assert!(
            nw1 > nw0,
            "|δ|=1 buffer must have strictly more non-white pixels in the chart area than |δ|=0 (nw0={nw0}, nw1={nw1})"
        );
    }

    #[test]
    fn bm_feature_with_none_effect_size_matches_welch_dot_footprint() {
        // BM with effect_size=None must fall back to su(3.0) — identical
        // chart-area rendering to a Welch render of the same feature.
        let opts = bm_opts_800();

        let mut r_bm = small_result();
        r_bm.method = DamMethod::BrunnerMunzel;
        let mut f_bm = feat(0.5, 1.5);
        f_bm.effect_size = None;
        r_bm.features = vec![f_bm];

        let mut r_welch = small_result();
        r_welch.method = DamMethod::Welch;
        r_welch.features = vec![feat(0.5, 1.5)];

        let b_bm = render_volcano(&r_bm, &opts).expect("BM None");
        let b_welch = render_volcano(&r_welch, &opts).expect("Welch");

        // Chart-main area only (excludes the legend column where BM grows
        // a size section, and excludes the annotation strip where the
        // method-name text differs). In this region both renders should be
        // pixel-identical: same axes, gridlines, threshold lines, and one
        // su(3.0) grey NS dot at the same coords.
        let nw_bm = count_non_white_in_box(&b_bm, 800, 0, 580, 0, 684);
        let nw_welch = count_non_white_in_box(&b_welch, 800, 0, 580, 0, 684);
        assert_eq!(
            nw_bm, nw_welch,
            "BM with effect_size=None must render an identical chart-area footprint to Welch (BM={nw_bm}, Welch={nw_welch})"
        );
    }

    #[test]
    fn bm_legend_grows_size_section_welch_legend_does_not() {
        // The |δ| size legend section lives at sp(290..478) × the legend
        // column (x ∈ [580, 800] at scale=1.0). BM renders MUST draw into
        // this band; Welch renders MUST leave it all white.
        let opts = bm_opts_800();

        let mut r_bm = small_result();
        r_bm.method = DamMethod::BrunnerMunzel;
        let b_bm = render_volcano(&r_bm, &opts).expect("BM");
        let b_welch = render_volcano(&small_result(), &opts).expect("Welch");

        let nw_bm = count_non_white_in_box(&b_bm, 800, 580, 800, 290, 478);
        let nw_welch = count_non_white_in_box(&b_welch, 800, 580, 800, 290, 478);

        assert!(
            nw_bm > 0,
            "BM render must draw the |δ| size legend section (nw_bm={nw_bm})"
        );
        assert_eq!(
            nw_welch, 0,
            "Welch render must leave the size-section band all white (nw_welch={nw_welch})"
        );
    }

    #[test]
    #[should_panic(expected = "Cliff's δ out of range")]
    fn bm_render_panics_on_out_of_range_cliffs_delta_in_debug_builds() {
        // Cliff's δ ∈ [-1, 1] by definition; the renderer asserts this so
        // a future cliffs_delta regression surfaces loudly in dev/test.
        let mut r = small_result();
        r.method = DamMethod::BrunnerMunzel;
        let mut f = feat(0.5, 1.5);
        f.effect_size = Some(1.5);
        r.features = vec![f];

        let opts = bm_opts_800();
        let _ = render_volcano(&r, &opts);
    }
}
