//! Shared, chart-agnostic scaffold for the two `plotters` renderers
//! (`volcano.rs`, `dotplot.rs`). Holds ONLY the pieces that were
//! byte-identical across both files: the RGBA buffer expansion, the PNG +
//! pHYs encoder, the 800-px design-scale baseline, the `sp` / `su` scaling
//! functions, and the major-gridline grey.
//!
//! Deliberately NOT here (per `extract-plot-common-scaffold` design D1): the
//! annotation strips, mesh / axis configuration, legends, and colour ramps —
//! those differ in text and geometry between the two charts and stay per-file
//! so the extraction can't move a single pixel.
//!
//! All items are `pub(crate)`: this is internal scaffold, not part of the
//! public plot API (`render_*` / `export_*` / `*Opts`).

use anyhow::{Context, Result, anyhow};
use plotters::style::RGBColor;
use std::fs::File;
use std::io::BufWriter;
use std::path::Path;

/// 800-px design baseline: every pixel-denominated constant in both renderers
/// scales relative to this so charts stay legible from the in-window preview
/// (≈ 800 px) up to a 600-DPI export (6000 px).
pub(crate) const DESIGN_BASELINE_PX: f64 = 800.0;

/// Major-gridline grey (soft, on a white canvas). Shared by both renderers
/// (volcano used this value inline; dotplot had it as a local `GRID_LIGHT`).
pub(crate) const GRID_LIGHT: RGBColor = RGBColor(210, 210, 210);

/// `scale = (min(w, h) / 800).max(0.25)` — the shared design-scale baseline.
pub(crate) fn design_scale(width_px: u32, height_px: u32) -> f64 {
    (width_px.min(height_px) as f64 / DESIGN_BASELINE_PX).max(0.25)
}

/// `scale = (w / 800).max(0.25)` — width-only variant for the dot plot.
///
/// The dot plot's height auto-sizes to the displayed-row count
/// (`App::stage3_autosize_height_in`) while its width is a fixed
/// `inches × dpi` resolution signal. Keying the font/element scale off
/// `min(w, h)` (the shared `design_scale`) made every font shrink whenever a
/// sparse result's auto-height dropped below the width — an 8 px label at 2
/// rows vs 18 px at 20 rows. Scaling by width alone makes the dot-plot scale
/// constant per (width, dpi), independent of how many entries pass. The
/// volcano keeps `design_scale` — its height is a fixed user setting that does
/// not auto-vary, so it never trips this.
pub(crate) fn design_scale_by_width(width_px: u32) -> f64 {
    (width_px as f64 / DESIGN_BASELINE_PX).max(0.25)
}

/// Scaled pixel (signed): `round(v * scale)`.
pub(crate) fn sp(v: f64, scale: f64) -> i32 {
    (v * scale).round() as i32
}

/// Scaled pixel (unsigned, floored at 1): `round(v * scale).max(1)`.
pub(crate) fn su(v: f64, scale: f64) -> u32 {
    (v * scale).round().max(1.0) as u32
}

/// Plotters writes RGB (3 bytes per pixel); egui / `image` PNG export take
/// RGBA. Allocate a fresh RGBA buffer and copy, setting every alpha to 255.
pub(crate) fn rgb_to_rgba(rgb: &[u8], pixel_count: usize) -> Vec<u8> {
    let mut out = vec![0u8; pixel_count * 4];
    for i in 0..pixel_count {
        out[i * 4] = rgb[i * 3];
        out[i * 4 + 1] = rgb[i * 3 + 1];
        out[i * 4 + 2] = rgb[i * 3 + 2];
        out[i * 4 + 3] = 255;
    }
    out
}

/// Write an RGBA `buffer` (`width_px * height_px * 4` bytes) as a PNG at `out`,
/// embedding the requested `dpi` as a pHYs chunk so the file declares its
/// physical pixel density (pixels-per-meter). Layout tools (Word, InDesign)
/// read the chunk to set the on-page size; callers that don't care can pass any
/// dpi — the pixel data is unaffected. Shared by both `export_*_png` paths.
pub(crate) fn encode_png(
    buffer: &[u8],
    width_px: u32,
    height_px: u32,
    dpi: u32,
    out: &Path,
) -> Result<()> {
    let file = File::create(out)
        .with_context(|| format!("failed to create PNG file at {}", out.display()))?;
    let writer = BufWriter::new(file);

    let mut encoder = png::Encoder::new(writer, width_px, height_px);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    // Convert DPI (dots per inch) → pixels per meter. 1 inch = 0.0254 m.
    let ppu = ((dpi as f64) / 0.0254).round() as u32;
    encoder.set_pixel_dims(Some(png::PixelDimensions {
        xppu: ppu,
        yppu: ppu,
        unit: png::Unit::Meter,
    }));
    let mut png_writer = encoder
        .write_header()
        .map_err(|e| anyhow!("PNG header write: {e}"))?;
    png_writer
        .write_image_data(buffer)
        .map_err(|e| anyhow!("PNG data write: {e}"))?;
    Ok(())
}

// ── Shared chart vocabulary (dot plots) ────────────────────────────────────

/// ColorBrewer YlOrRd 9-step sequential palette. Used as a piecewise-
/// linear lookup so the continuous gradient stays close to ColorBrewer's
/// designed perceptual curve — a 2-anchor sRGB lerp between the
/// extremes (which is what this code used to do) produced a muddy
/// brownish mid-tone (~`#dd6d4e`) the curated palette specifically
/// avoids; the curated mid-bin is `#fd8d3c` (a vivid orange).
///
/// Index 0 = palest yellow `#ffffcc` (t=0, FDR at threshold; least
/// significant in the displayed band). Index 8 = deepest red `#800026`
/// (t=1, FDR at the displayed maximum `-log10`; most significant).
///
/// Source: <https://colorbrewer2.org/#type=sequential&scheme=YlOrRd&n=9>
pub(crate) const YLORRD_9: [(u8, u8, u8); 9] = [
    (255, 255, 204), // #ffffcc  palest yellow (t=0)
    (255, 237, 160), // #ffeda0
    (254, 217, 118), // #fed976
    (254, 178, 76),  // #feb24c
    (253, 141, 60),  // #fd8d3c  curated mid-bin
    (252, 78, 42),   // #fc4e2a
    (227, 26, 28),   // #e31a1c
    (189, 0, 38),    // #bd0026
    (128, 0, 38),    // #800026  deepest red (t=1)
];

/// Pure gradient lookup over the ColorBrewer YlOrRd 9-step palette
/// using piecewise-linear interpolation between adjacent anchors.
/// `t=0` returns palest yellow (`YLORRD_9[0]`), `t=1` returns deepest
/// red (`YLORRD_9[8]`), and any `t` in between lerps within whichever
/// of the 8 adjacent-anchor segments it falls into. Out-of-range `t`
/// values are clamped.
pub(crate) fn t_to_color(t: f64) -> RGBColor {
    let t = t.clamp(0.0, 1.0);
    let last = YLORRD_9.len() - 1;
    let pos = t * last as f64;
    let lo = (pos.floor() as usize).min(last);
    let hi = (lo + 1).min(last);
    let frac = pos - lo as f64;
    let a = YLORRD_9[lo];
    let b = YLORRD_9[hi];
    RGBColor(
        lerp(a.0, b.0, frac),
        lerp(a.1, b.1, frac),
        lerp(a.2, b.2, frac),
    )
}

/// Greedy word-wrap a label into up to `max_lines` lines, each at most
/// `chars_per_line` chars. The last kept line is suffixed with `…` only
/// when the original would have produced more than `max_lines` lines.
/// Words longer than `chars_per_line` are hard-truncated with `…` so
/// the layout never overflows the label area.
///
/// Returns at least one line (an empty string for empty input) so the
/// caller can unconditionally iterate.
pub(crate) fn wrap_label(name: &str, chars_per_line: usize, max_lines: usize) -> Vec<String> {
    if max_lines == 0 {
        return Vec::new();
    }
    let trimmed = name.trim();
    if trimmed.is_empty() {
        return vec![String::new()];
    }

    let mut lines: Vec<String> = Vec::new();
    let mut current = String::new();
    for w in trimmed.split_whitespace() {
        let needs_space = !current.is_empty();
        let candidate_len = current.chars().count() + usize::from(needs_space) + w.chars().count();
        if candidate_len <= chars_per_line {
            if needs_space {
                current.push(' ');
            }
            current.push_str(w);
        } else {
            if !current.is_empty() {
                lines.push(std::mem::take(&mut current));
            }
            if w.chars().count() > chars_per_line {
                // Single word exceeds the per-line budget; hard-truncate.
                let kept: String = w.chars().take(chars_per_line.saturating_sub(1)).collect();
                current = kept;
                current.push('…');
            } else {
                current = w.to_string();
            }
        }
    }
    if !current.is_empty() {
        lines.push(current);
    }

    if lines.len() > max_lines {
        lines.truncate(max_lines);
        let last = lines.last_mut().unwrap();
        // Append "…" without exceeding the per-line budget.
        if last.chars().count() + 1 > chars_per_line {
            let kept: String = last
                .chars()
                .take(chars_per_line.saturating_sub(1))
                .collect();
            *last = kept;
        }
        last.push('…');
    }

    lines
}

fn lerp(a: u8, b: u8, t: f64) -> u8 {
    let t = t.clamp(0.0, 1.0);
    (a as f64 + (b as f64 - a as f64) * t).round() as u8
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rgb_to_rgba_expands_with_opaque_alpha() {
        // 2 pixels: red, green.
        let rgb = [255u8, 0, 0, 0, 255, 0];
        let out = rgb_to_rgba(&rgb, 2);
        assert_eq!(out.len(), 2 * 4, "RGBA buffer is pixel_count * 4 bytes");
        assert_eq!(out, [255, 0, 0, 255, 0, 255, 0, 255]);
        // Every 4th byte (alpha) is opaque.
        assert!(out.chunks(4).all(|p| p[3] == 255));
    }

    #[test]
    fn design_scale_baseline_floor_and_growth() {
        assert_eq!(design_scale(800, 800), 1.0, "1.0 at the 800-px baseline");
        // Below 200 px the raw ratio (< 0.25) is floored at 0.25.
        assert_eq!(design_scale(100, 100), 0.25, "floored at 0.25 below 200 px");
        assert_eq!(design_scale(150, 4000), 0.25, "floor keys off min(w, h)");
        assert_eq!(design_scale(6000, 6000), 7.5, "7.5 at 6000 px");
    }

    #[test]
    fn design_scale_by_width_ignores_height_and_keys_off_width() {
        // Same width → same scale regardless of height. This is the whole point:
        // a sparse dot plot (short canvas) renders fonts at the SAME size as a
        // dense one (tall canvas).
        assert_eq!(
            design_scale_by_width(800),
            1.0,
            "1.0 at the 800-px baseline"
        );
        assert_eq!(
            design_scale_by_width(1050),
            1050.0 / 800.0,
            "3.5in × 300dpi default → 1.3125, independent of row count"
        );
        // Height is not a parameter, so a 480-px-tall vs 2100-px-tall canvas at
        // the same 1050-px width both yield 1.3125 (the bug was the opposite).
        assert_eq!(
            design_scale_by_width(100),
            0.25,
            "floored at 0.25 below 200 px"
        );
        assert_eq!(design_scale_by_width(6000), 7.5, "7.5 at 6000 px");
    }

    #[test]
    fn sp_su_round_and_su_floors_at_one() {
        assert_eq!(sp(2.4, 2.0), 5, "round(4.8) = 5");
        assert_eq!(sp(12.0, 1.0), 12);
        assert_eq!(su(2.6, 1.0), 3, "round(2.6) = 3");
        // su never returns 0 even for a sub-half product.
        assert_eq!(su(0.01, 1.0), 1, "su floors at 1");
        assert_eq!(su(0.0, 5.0), 1, "su floors at 1 even for 0");
    }

    #[test]
    fn encode_png_writes_phys_chunk_with_dpi() {
        // A minimal 2×2 RGBA buffer is enough to exercise the encoder.
        let buffer = vec![0u8; 2 * 2 * 4];
        let tmp = tempfile::Builder::new()
            .suffix(".png")
            .tempfile()
            .expect("tempfile");
        encode_png(&buffer, 2, 2, 300, tmp.path()).expect("encode");
        let bytes = std::fs::read(tmp.path()).expect("read png");

        // Walk PNG chunks (8B sig; each chunk: 4B len, 4B type, payload, 4B CRC)
        // and decode the pHYs payload's xppu / yppu / unit.
        let mut idx = 8;
        let mut found: Option<(u32, u32, u8)> = None;
        while idx + 8 < bytes.len() {
            let len = u32::from_be_bytes(bytes[idx..idx + 4].try_into().unwrap()) as usize;
            let ty = &bytes[idx + 4..idx + 8];
            if ty == b"pHYs" {
                let payload = &bytes[idx + 8..idx + 8 + len];
                let xppu = u32::from_be_bytes(payload[0..4].try_into().unwrap());
                let yppu = u32::from_be_bytes(payload[4..8].try_into().unwrap());
                found = Some((xppu, yppu, payload[8]));
                break;
            }
            idx += 8 + len + 4;
        }
        let (xppu, yppu, unit) = found.expect("pHYs chunk not found");
        let expected = ((300.0_f64) / 0.0254).round() as u32;
        assert_eq!(xppu, expected, "xppu encodes 300 DPI");
        assert_eq!(yppu, expected, "yppu encodes 300 DPI");
        assert_eq!(unit, 1, "unit byte is Meter (1)");
    }

    #[test]
    fn wrap_label_short_name_returns_single_line_unchanged() {
        assert_eq!(wrap_label("Glycolysis", 28, 2), vec!["Glycolysis"]);
    }

    #[test]
    fn wrap_label_long_name_breaks_at_word_boundary_into_two_lines() {
        let lines = wrap_label("Pentose and glucuronate interconversions", 28, 2);
        assert_eq!(lines.len(), 2);
        assert_eq!(lines[0], "Pentose and glucuronate");
        assert_eq!(lines[1], "interconversions");
        // No line exceeds the budget.
        assert!(lines.iter().all(|l| l.chars().count() <= 28));
    }

    #[test]
    fn wrap_label_exceeding_max_lines_ellipsis_truncates_last_line() {
        let lines = wrap_label("Amino sugar and nucleotide sugar metabolism extras", 18, 2);
        assert_eq!(lines.len(), 2);
        assert!(lines.last().unwrap().ends_with('…'));
        assert!(lines.iter().all(|l| l.chars().count() <= 18));
    }

    #[test]
    fn wrap_label_single_word_longer_than_budget_hard_truncates() {
        // "Supercalifragilisticexpialidocious" is 34 chars; budget 20.
        let lines = wrap_label("Supercalifragilisticexpialidocious", 20, 2);
        assert_eq!(lines.len(), 1);
        assert!(lines[0].ends_with('…'));
        assert_eq!(lines[0].chars().count(), 20);
    }

    #[test]
    fn wrap_label_max_lines_one_forces_single_line_with_ellipsis() {
        let lines = wrap_label("Pentose and glucuronate interconversions", 28, 1);
        assert_eq!(lines.len(), 1);
        assert!(lines[0].ends_with('…'));
        assert!(lines[0].chars().count() <= 28);
    }

    #[test]
    fn wrap_label_empty_input_returns_single_empty_line() {
        assert_eq!(wrap_label("", 28, 2), vec![String::new()]);
        assert_eq!(wrap_label("   ", 28, 2), vec![String::new()]);
    }

    #[test]
    fn wrap_label_three_lines_for_very_long_name() {
        // A name needing three lines wraps fully at the 26-char budget with
        // no ellipsis (previously truncated at the old 2-line cap).
        let lines = wrap_label(
            "Glycosaminoglycan biosynthesis heparan sulfate keratan route",
            26,
            3,
        );
        assert_eq!(lines.len(), 3, "got: {lines:?}");
        assert!(lines.iter().all(|l| l.chars().count() <= 26));
        assert!(
            !lines.last().unwrap().ends_with('…'),
            "fits in 3 lines, no ellipsis: {lines:?}"
        );
    }

    #[test]
    fn wrap_label_four_lines_for_very_long_name() {
        // MAX_LABEL_LINES is now 4: a name needing four lines wraps fully at
        // the budget with no ellipsis (previously capped at 3).
        let lines = wrap_label(
            "Amino sugar and nucleotide sugar metabolism related glycan degradation",
            20,
            4,
        );
        assert_eq!(lines.len(), 4, "got: {lines:?}");
        assert!(lines.iter().all(|l| l.chars().count() <= 20));
        assert!(
            !lines.last().unwrap().ends_with('…'),
            "fits in 4 lines, no ellipsis: {lines:?}"
        );
    }

    #[test]
    fn ylorrd_9_endpoints_match_colorbrewer() {
        // Locks the palette in place: a future copy-paste mistake on the
        // ColorBrewer hex values would silently shift the whole gradient.
        assert_eq!(YLORRD_9.len(), 9);
        assert_eq!(YLORRD_9[0], (255, 255, 204)); // #ffffcc
        assert_eq!(YLORRD_9[4], (253, 141, 60)); // #fd8d3c (curated mid-bin)
        assert_eq!(YLORRD_9[8], (128, 0, 38)); // #800026
    }

    #[test]
    fn t_to_color_hits_each_colorbrewer_anchor_exactly() {
        // The 9 anchors must land at t = i/8 with zero interpolation
        // error — confirms the piecewise-linear lookup is anchored on
        // the curated bins, not on an off-by-one neighbour.
        for (i, &(r, g, b)) in YLORRD_9.iter().enumerate() {
            let t = i as f64 / 8.0;
            assert_eq!(
                t_to_color(t),
                RGBColor(r, g, b),
                "anchor {i} (t={t}) drifted"
            );
        }
    }

    #[test]
    fn t_to_color_midpoint_close_to_colorbrewer_mid_bin() {
        // The whole point of the LUT: t=0.5 should land near ColorBrewer's
        // curated mid-bin `#fd8d3c` (vivid orange) — NOT the muddy
        // brownish-orange (~`#dd6d4e`) a 2-anchor sRGB lerp between the
        // extremes would have produced. With 9 anchors, t=0.5 lands exactly
        // on the curated mid-bin.
        let mid = t_to_color(0.5);
        let (r, g, b) = YLORRD_9[4];
        assert_eq!(mid, RGBColor(r, g, b));
        // sRGB lerp G channel of #800026 (G=0) ↔ #ffffcc (G=255): (0 + 255) / 2.
        let two_anchor_mid_green = 255 / 2;
        assert!(
            (mid.1 as i32 - two_anchor_mid_green).abs() > 10,
            "ColorBrewer mid-bin should diverge from a naive 2-anchor lerp"
        );
    }
}
