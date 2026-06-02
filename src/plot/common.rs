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
}
