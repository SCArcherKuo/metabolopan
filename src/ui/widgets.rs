//! Shared widget helpers — the `ui-color.md` §2–§5 component layer.
//!
//! egui has ONE global widget theme; `theme::install` makes the DEFAULT look
//! the §2/§3 **Secondary** component. **Primary** styling is opt-in here, via a
//! scoped `Visuals` override (`ui.scope`) that drives egui's own `interact()`
//! state machine — NOT `Button::fill()` (which freezes one fill across all
//! states and kills hover feedback). The scope discards the override on exit,
//! so a Primary widget never leaks styling onto its siblings.

use rfd::FileDialog;
use std::path::PathBuf;

use crate::dam::fdr::FdrMethod;
use crate::theme;
use egui::{Button, RichText, Stroke};

/// §2 **Primary button** — no border, fill `PRIMARY` (idle) → `PRIMARY_HOVER`
/// (hover) → `PRIMARY_ACTIVE` (pressed), white `ON_PRIMARY` label. Disabled
/// falls through to the §2 Primary-disabled look (`DISABLED` bg,
/// `TEXT_DISABLED` label). Use for the forward CTAs, file pickers, and the
/// `Draw …` actions (per `interactive-component-styles`).
pub fn primary_button(ui: &mut egui::Ui, label: &str, enabled: bool) -> egui::Response {
    primary_button_sized(ui, label, enabled, egui::Vec2::ZERO)
}

/// Like [`primary_button`] but with a minimum button size — for wide triggers
/// (e.g. the KEGG-species selector button).
pub fn primary_button_sized(
    ui: &mut egui::Ui,
    label: &str,
    enabled: bool,
    min_size: egui::Vec2,
) -> egui::Response {
    ui.scope(|ui| {
        let w = &mut ui.style_mut().visuals.widgets;
        w.inactive.weak_bg_fill = theme::PRIMARY;
        w.inactive.bg_stroke = Stroke::NONE;
        w.hovered.weak_bg_fill = theme::PRIMARY_HOVER;
        w.hovered.bg_stroke = Stroke::NONE;
        w.active.weak_bg_fill = theme::PRIMARY_ACTIVE;
        w.active.bg_stroke = Stroke::NONE;
        // Disabled buttons render from `noninteractive`.
        w.noninteractive.weak_bg_fill = theme::DISABLED;
        w.noninteractive.bg_stroke = Stroke::NONE;
        let label_color = if enabled {
            theme::ON_PRIMARY
        } else {
            theme::TEXT_DISABLED
        };
        ui.add_enabled(
            enabled,
            Button::new(RichText::new(label).color(label_color)).min_size(min_size),
        )
    })
    .inner
}

/// Inner padding (px) of a segmented tab: 10 horizontal, 4 vertical.
const SEG_PADDING: egui::Vec2 = egui::vec2(10.0, 4.0);

/// Horizontal gap (px) between a [`segmented_tab_with_icon`] leading icon and
/// its label. Visual-polish constant (see `add-stepper-step-icons` design D3).
const ICON_LABEL_GAP: f32 = 4.0;

/// Leading-icon size as a multiple of the label line height. `> 1.0`, so the
/// icon reads clearly larger than the text; the segment (and therefore the row)
/// grows taller to contain it — see `add-stepper-step-icons` design D3. At the
/// default Button text (~16 px line height) this yields a ~48 px icon.
const ICON_SCALE: f32 = 3.0;

/// The segment height a stepper row needs so a `segmented_tab_with_icon` icon
/// (`ICON_SCALE`× the Button line height) fits. Pass this as the `min_height` of
/// EVERY segment in the row (icon-bearing or not) so they share one height.
///
/// Why pin all segments: egui's `ui.horizontal` lays items out in a single pass
/// and cross-axis-centres each within the row height KNOWN SO FAR. A taller
/// segment added later grows the row but never reflows the earlier, shorter
/// segments — leaving everything before it visually shifted up. Giving every
/// segment the same height makes the FIRST segment establish the full row
/// height, so all segments (and the `>` separators) centre identically. At the
/// default Button text (~16 px) this is ~46 px.
pub fn segmented_icon_row_height(ui: &egui::Ui) -> f32 {
    ui.text_style_height(&egui::TextStyle::Button) * ICON_SCALE + SEG_PADDING.y * 2.0
}

/// §4 **Segmented tab** — one segment of a Segmented Control (sits inside a
/// [`segmented_track`]). `primary` picks the level:
/// - **Primary** selected = white (`ON_PRIMARY`) fill + soft drop shadow
///   ("floating slider") + `PRIMARY` text; hover = `FILL_SECONDARY_HOVER` fill +
///   `TEXT`.
/// - **Secondary** selected = flat `FILL_SECONDARY_HOVER` fill + `HEADING` text,
///   no shadow; hover = NO fill, text darkens to `TEXT`.
///
/// Both: unselected = transparent + `TEXT_SECONDARY`; disabled-unselected =
/// transparent + `TEXT_DISABLED`; disabled-selected = `DISABLED_FILL` +
/// `TEXT_DISABLED` (Primary keeps its float shadow). Custom-painted (not
/// `selectable_label`) so the shadow + the text-only Secondary hover are exact.
///
/// Delegates to [`segmented_tab_with_icon`] with `icon = None`, `min_height = 0`
/// (mirrors `primary_button` → `primary_button_sized`), so this label-only path
/// stays byte-identical to the pre-icon implementation.
pub fn segmented_tab(
    ui: &mut egui::Ui,
    selected: bool,
    label: &str,
    enabled: bool,
    primary: bool,
) -> egui::Response {
    segmented_tab_with_icon(ui, selected, label, enabled, primary, None, 0.0)
}

/// As [`segmented_tab`] but with an optional leading `icon` and a `min_height`
/// floor on the segment's allocated height.
///
/// When `icon = Some(id)` a square icon box ([`ICON_SCALE`]× the label line
/// height) is painted to the LEFT of the label: the segment widens by the icon +
/// [`ICON_LABEL_GAP`] and grows taller to contain it (so the icon reads clearly
/// larger than the text). The `[icon][gap][label]` block is centred so the
/// left/right padding stays symmetric (design D3); icon and label are both
/// vertically centred. The icon is painted with a white tint (no recolour) over
/// the full UV. The §4 per-state fill / shadow / text-colour selection is
/// unchanged — the icon is additive.
///
/// `min_height` floors the allocated segment height (content + padding is used
/// when larger). Callers pin every segment in a row to the same `min_height`
/// (see [`segmented_icon_row_height`]) so a row mixing icon and icon-less
/// segments stays vertically aligned. `min_height = 0.0` imposes no floor.
///
/// With `icon = None` and `min_height = 0.0` the layout and painting are
/// byte-identical to the original `segmented_tab`. The icon param is a
/// `TextureId` (not a `TextureHandle`): the caller's cache owns the lifetime.
pub fn segmented_tab_with_icon(
    ui: &mut egui::Ui,
    selected: bool,
    label: &str,
    enabled: bool,
    primary: bool,
    icon: Option<egui::TextureId>,
    min_height: f32,
) -> egui::Response {
    let padding = SEG_PADDING;
    let font_id = egui::TextStyle::Button.resolve(ui.style());
    let galley =
        ui.fonts(|f| f.layout_no_wrap(label.to_owned(), font_id, egui::Color32::PLACEHOLDER));
    // Icon box is `ICON_SCALE`× the label line height; when present, widen the
    // segment by it + gap and keep the content tall enough to hold it. The final
    // height is floored at `min_height` (so a whole row can share one height).
    // With no icon and `min_height = 0`, this reduces to the original
    // `galley.size() + padding*2` — byte-identical.
    let icon_size = galley.size().y * ICON_SCALE;
    let (icon_extra, content_h) = if icon.is_some() {
        (icon_size + ICON_LABEL_GAP, icon_size.max(galley.size().y))
    } else {
        (0.0, galley.size().y)
    };
    let desired = egui::vec2(
        galley.size().x + icon_extra + padding.x * 2.0,
        (content_h + padding.y * 2.0).max(min_height),
    );
    let sense = if enabled {
        egui::Sense::click()
    } else {
        egui::Sense::hover()
    };
    let (rect, response) = ui.allocate_exact_size(desired, sense);
    let radius = egui::CornerRadius::same(4);

    if ui.is_rect_visible(rect) {
        // (fill, text, with_shadow) per the §4 Navigation-Tabs state table.
        let (fill, text_color, with_shadow) = if !enabled {
            if selected {
                (Some(theme::DISABLED_FILL), theme::TEXT_DISABLED, primary)
            } else {
                (None, theme::TEXT_DISABLED, false)
            }
        } else if selected {
            if primary {
                (Some(theme::ON_PRIMARY), theme::PRIMARY, true)
            } else {
                (Some(theme::FILL_SECONDARY_HOVER), theme::HEADING, false)
            }
        } else if response.hovered() {
            if primary {
                (Some(theme::FILL_SECONDARY_HOVER), theme::TEXT, false)
            } else {
                // Secondary hover: no fill — only the text darkens.
                (None, theme::TEXT, false)
            }
        } else {
            (None, theme::TEXT_SECONDARY, false)
        };

        let painter = ui.painter();
        if with_shadow {
            // box-shadow: 0 2px 6px rgba(45,44,42,0.08) — painted BEFORE the
            // fill so it reads as a soft float under the segment.
            let shadow = egui::Shadow {
                offset: [0, 2],
                blur: 6,
                spread: 0,
                color: egui::Color32::from_rgba_unmultiplied(45, 44, 42, 20),
            };
            painter.add(shadow.as_shape(rect, radius));
        }
        if let Some(fill) = fill {
            painter.rect_filled(rect, radius, fill);
        }
        match icon {
            // `[icon][gap][label]` laid out from the left padding edge; because
            // the segment was widened by exactly `icon_size + gap`, this keeps
            // the padding symmetric (design D3).
            Some(texture_id) => {
                let block_left = rect.min.x + padding.x;
                let icon_rect = egui::Rect::from_min_size(
                    egui::pos2(block_left, rect.center().y - icon_size * 0.5),
                    egui::vec2(icon_size, icon_size),
                );
                painter.image(
                    texture_id,
                    icon_rect,
                    egui::Rect::from_min_max(egui::pos2(0.0, 0.0), egui::pos2(1.0, 1.0)),
                    egui::Color32::WHITE,
                );
                let text_pos = egui::pos2(
                    block_left + icon_size + ICON_LABEL_GAP,
                    rect.center().y - galley.size().y * 0.5,
                );
                painter.galley(text_pos, galley, text_color);
            }
            // Byte-identical to the pre-change path.
            None => {
                let text_pos = rect.center() - galley.size() * 0.5;
                painter.galley(text_pos, galley, text_color);
            }
        }
    }

    response
}

/// Wrap a row of [`segmented_tab`]s in the §4 Primary track container — a
/// `FILL_SECONDARY` rounded background that gives the segmented control its
/// "track + floating slider" look. Inner margin is `4 px` horizontal but only
/// `2 px` vertical, so the track hugs the tabs — it sits only a little taller
/// than the tab segments rather than as a tall slab.
pub fn segmented_track<R>(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    egui::Frame::NONE
        .fill(theme::FILL_SECONDARY)
        .inner_margin(egui::Margin::symmetric(4, 2))
        .corner_radius(egui::CornerRadius::same(6))
        .show(ui, |ui| ui.horizontal(|ui| add(ui)).inner)
        .inner
}

/// §3.3 **Primary dropdown** styling — white (`ON_PRIMARY`) trigger bg +
/// `SURFACE` border, with the open panel's selected item on opaque `PRIMARY`.
/// Wrap the `egui::ComboBox` call in `add`; the caller colours the selected
/// list item `ON_PRIMARY` (see the Num/Den + Group selector call sites). Use
/// for the Numerator/Denominator and Module-mode Group dropdowns.
pub fn primary_dropdown<R>(ui: &mut egui::Ui, add: impl FnOnce(&mut egui::Ui) -> R) -> R {
    ui.scope(|ui| {
        let v = &mut ui.style_mut().visuals;
        v.widgets.inactive.weak_bg_fill = theme::ON_PRIMARY;
        v.widgets.inactive.bg_stroke = Stroke::new(1.0_f32, theme::SURFACE);
        v.widgets.open.weak_bg_fill = theme::ON_PRIMARY;
        v.selection.bg_fill = theme::PRIMARY;
        add(ui)
    })
    .inner
}

/// Render an `egui::ProgressBar` per `ui-color.md` §5: a `FILL_SECONDARY` track
/// (via scoped `extreme_bg_color`) with the given indicator `fill`. Running
/// screens pass `PRIMARY` (in-progress) / `SUCCESS` (done); inline fetch /
/// refresh strips pass `SURFACE` (Secondary).
pub fn progress_bar(
    ui: &mut egui::Ui,
    bar: egui::ProgressBar,
    fill: egui::Color32,
) -> egui::Response {
    ui.scope(|ui| {
        ui.style_mut().visuals.extreme_bg_color = theme::FILL_SECONDARY;
        ui.add(bar.fill(fill))
    })
    .inner
}

/// Render a file picker button (styled §2 Primary — file pickers are a core
/// action per `interactive-component-styles`). When clicked and enabled, opens
/// the system dialog filtered by `filter_label` / `filter_ext` and returns the
/// chosen path (if any).
pub fn file_pick_button(
    ui: &mut egui::Ui,
    label: &str,
    filter_label: &str,
    filter_ext: &str,
    enabled: bool,
) -> Option<PathBuf> {
    let response = primary_button(ui, label, enabled);
    if response.clicked() {
        FileDialog::new()
            .add_filter(filter_label, &[filter_ext])
            .pick_file()
    } else {
        None
    }
}

/// The shared **"PNG export size"** control block — a `theme::HEADING` heading
/// plus a one-row `Width (in) / Height (in) / DPI` triple of `DragValue`s
/// (speeds `0.1 / 0.1 / 10`, ranges `1.0..=40.0 / 1.0..=40.0 / 72..=1200`).
/// Rendered byte-identically on the Stage 2 threshold and Stage 3 result screens.
pub(crate) fn png_export_size_controls(
    ui: &mut egui::Ui,
    width_in: &mut f64,
    height_in: &mut f64,
    dpi: &mut u32,
) {
    ui.label(
        RichText::new("PNG export size")
            .strong()
            .color(theme::HEADING),
    );
    ui.horizontal(|ui| {
        ui.label("Width (in):");
        ui.add(egui::DragValue::new(width_in).speed(0.1).range(1.0..=40.0));
        ui.label("Height (in):");
        ui.add(egui::DragValue::new(height_in).speed(0.1).range(1.0..=40.0));
        ui.label("DPI:");
        ui.add(egui::DragValue::new(dpi).speed(10).range(72..=1200));
    });
}

/// Convert inches × DPI into a clamped pixel pair for plot export. The clamp
/// (`64..=20_000` per axis) guards `plotters` / `png` against absurd sizes.
/// Used by both the volcano (Stage 2) and dot-plot (Stage 3) export paths.
pub(crate) fn export_pixels(width_in: f64, height_in: f64, dpi: u32) -> (u32, u32) {
    let w = (width_in * dpi as f64).round().clamp(64.0, 20_000.0) as u32;
    let h = (height_in * dpi as f64).round().clamp(64.0, 20_000.0) as u32;
    (w, h)
}

/// The shared **save-file dialog** scaffold — a single-filter `rfd` save dialog.
/// Returns the chosen path, or `None` if the user cancelled. Callers keep their
/// own `File::create` + error log + `export_*` call after this (those differ per
/// export: PNG vs CSV vs dedup-audit, distinct error messages, `&mut file` vs
/// path-only).
pub(crate) fn save_dialog(
    filter_label: &str,
    filter_ext: &str,
    default_name: &str,
) -> Option<PathBuf> {
    FileDialog::new()
        .add_filter(filter_label, &[filter_ext])
        .set_file_name(default_name)
        .save_file()
}

/// The (FdrMethod, label) options for [`fdr_method_radios`], in render order.
/// BH + BY always; the `No correction (raw p-values)` (`NoCorrection`) variant only when
/// `include_none` — Stage 3 setup exposes it, Stage 2 setup hides it. Split out
/// as a pure fn so the gating + verbatim label strings are unit-testable.
fn fdr_radio_options(include_none: bool) -> Vec<(FdrMethod, &'static str)> {
    let mut options = vec![
        (
            FdrMethod::BenjaminiHochberg,
            "Benjamini–Hochberg (BH) procedure",
        ),
        (
            FdrMethod::BenjaminiYekutieli,
            "Benjamini–Yekutieli (BY) procedure",
        ),
    ];
    if include_none {
        options.push((FdrMethod::NoCorrection, "No correction (raw p-values)"));
    }
    options
}

/// The shared **FDR-method radio group** — a `FDR correction:` label plus the
/// BH and BY radios (verbatim label strings), and — only when `include_none` —
/// the `No correction (raw p-values)` radio (`FdrMethod::NoCorrection`). Stage 2 setup
/// calls it with `include_none = false` (the `None` variant is Stage-3-only);
/// Stage 3 setup with `include_none = true`. The per-screen grey sub-hint below
/// the radios stays at the call site (the two hint texts differ).
pub(crate) fn fdr_method_radios(ui: &mut egui::Ui, method: &mut FdrMethod, include_none: bool) {
    ui.label("FDR correction:");
    for (value, label) in fdr_radio_options(include_none) {
        ui.radio_value(method, value, label);
    }
}

/// A `- `-prefixed section header → a bold bullet title (`• <text>`) in the
/// darkest ink (`HEADING`), matching the `ui-design.md` `- **Header**`
/// convention. NB: egui's default fonts carry no bold *weight* (`.strong()` only
/// shifts colour, and in this theme `strong_text_color == TEXT`, so it's a
/// no-op), so "bold" is the `HEADING` ink — the app's established emphasis —
/// here paired with the `•` bullet marker. Promoted from `data_tab.rs` so other
/// screens can reuse the convention.
pub(crate) fn section_header(ui: &mut egui::Ui, text: impl std::fmt::Display) {
    ui.label(RichText::new(format!("• {text}")).color(theme::HEADING));
}

/// A `Key: value` data line with the `Key:` label emphasised in `HEADING` ink
/// (the app's "bold") and the value in `value_color`, matching the
/// `ui-design.md` `**Key:**` convention. Rendered as one wrapping `LayoutJob`
/// label so it wraps and spaces naturally. When `text` has no `": "` the whole
/// string is emphasised (e.g. `Groups:`). Promoted from `data_tab.rs`.
pub(crate) fn kv_line_colored(ui: &mut egui::Ui, text: &str, value_color: egui::Color32) {
    let font = egui::TextStyle::Body.resolve(ui.style());
    match text.split_once(": ") {
        Some((key, value)) => {
            let mut job = egui::text::LayoutJob::default();
            job.append(
                &format!("{key}: "),
                0.0,
                egui::text::TextFormat {
                    font_id: font.clone(),
                    color: theme::HEADING,
                    ..Default::default()
                },
            );
            job.append(
                value,
                0.0,
                egui::text::TextFormat {
                    font_id: font,
                    color: value_color,
                    ..Default::default()
                },
            );
            ui.label(job);
        }
        None => {
            ui.label(RichText::new(text).color(theme::HEADING));
        }
    }
}

/// [`kv_line_colored`] with the value in body `TEXT`.
pub(crate) fn kv_line(ui: &mut egui::Ui, text: &str) {
    kv_line_colored(ui, text, theme::TEXT);
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The scope-based Primary helpers MUST restore their `Visuals` override so
    /// a sibling widget rendered afterwards is unaffected (Secondary). Render a
    /// frame and assert the inactive button fill is unchanged after they return.
    /// (`segmented_tab` paints directly and never mutates `Visuals`, so it has
    /// no leak surface.)
    #[test]
    fn primary_helpers_do_not_leak_style() {
        let ctx = egui::Context::default();
        theme::install(&ctx);
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                let before = ui.style().visuals.widgets.inactive.weak_bg_fill;
                assert_eq!(before, theme::FILL_SECONDARY, "default is Secondary");

                let _ = primary_button(ui, "Start", true);
                assert_eq!(
                    ui.style().visuals.widgets.inactive.weak_bg_fill,
                    before,
                    "primary_button must not leak its override"
                );

                let _ = primary_dropdown(ui, |ui| ui.label("x"));
                assert_eq!(
                    ui.style().visuals.widgets.inactive.weak_bg_fill,
                    before,
                    "primary_dropdown must not leak its override"
                );
            });
        });
    }

    /// The icon-less `segmented_tab` path must allocate identically to
    /// `segmented_tab_with_icon(.., None, 0.0)` (it delegates); supplying an icon
    /// must WIDEN the segment (icon box + gap) and — because `ICON_SCALE > 1` —
    /// also make it TALLER. A `min_height` floor must raise the segment to at
    /// least that height (the mechanism that keeps a stepper row aligned).
    #[test]
    fn segmented_tab_icon_none_matches_plain_and_icon_grows() {
        let ctx = egui::Context::default();
        theme::install(&ctx);
        let mut plain = egui::Rect::NOTHING;
        let mut none = egui::Rect::NOTHING;
        let mut icon = egui::Rect::NOTHING;
        let mut pinned = egui::Rect::NOTHING;
        let pin_h = 100.0_f32;
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            let tex = ctx.load_texture(
                "seg_icon",
                egui::ColorImage::from_rgba_unmultiplied([1, 1], &[0, 0, 0, 255]),
                egui::TextureOptions::LINEAR,
            );
            egui::CentralPanel::default().show(ctx, |ui| {
                plain = segmented_tab(ui, false, "Label", true, true).rect;
                none = segmented_tab_with_icon(ui, false, "Label", true, true, None, 0.0).rect;
                icon = segmented_tab_with_icon(ui, false, "Label", true, true, Some(tex.id()), 0.0)
                    .rect;
                pinned = segmented_tab_with_icon(ui, false, "Label", true, true, None, pin_h).rect;
            });
        });
        assert_eq!(
            plain.size(),
            none.size(),
            "icon-less delegation must match plain segmented_tab allocation"
        );
        assert!(
            icon.width() > none.width(),
            "an icon must widen the segment ({} vs {})",
            icon.width(),
            none.width()
        );
        assert!(
            icon.height() > none.height(),
            "an icon (ICON_SCALE > 1) must heighten the segment ({} vs {})",
            icon.height(),
            none.height()
        );
        assert!(
            pinned.height() >= pin_h,
            "min_height must floor the segment height ({} vs {})",
            pinned.height(),
            pin_h
        );
    }

    /// `export_pixels` reproduces the prior inches × DPI numerics and hits the
    /// `64` / `20_000` clamp boundaries (byte-equal to the two deleted copies).
    #[test]
    fn export_pixels_reproduces_prior_numerics() {
        assert_eq!(export_pixels(3.5, 2.2, 300), (1050, 660));
        assert_eq!(export_pixels(10.0, 8.0, 600), (6000, 4800));
        // Below the 64-px floor → clamped up.
        assert_eq!(export_pixels(0.01, 0.01, 72), (64, 64));
        // Above the 20_000-px ceiling → clamped down.
        assert_eq!(export_pixels(40.0, 40.0, 1200), (20_000, 20_000));
    }

    /// `fdr_radio_options` gates the `NoCorrection` variant on `include_none`
    /// and pins all three captions verbatim, so the BH/BY pair cannot drift
    /// between the Stage 2 and Stage 3 screens.
    ///
    /// These captions are NOT the strings any exporter reads — the CSV tag line
    /// renders `FdrMethod::short_label()` and the dot-plot annotation strip
    /// builds its wording locally. The third radio is where the two visibly
    /// differ: its caption is `No correction (raw p-values)`, its short label is
    /// `NoCorrection`. An earlier version of this comment called the captions a
    /// contract the exporters depend on; they never were.
    #[test]
    fn fdr_radio_options_pins_captions_and_gates_no_correction() {
        let without = fdr_radio_options(false);
        assert_eq!(without.len(), 2);
        assert!(
            !without.iter().any(|(m, _)| *m == FdrMethod::NoCorrection),
            "Stage 2 (include_none = false) must NOT offer NoCorrection"
        );
        let with = fdr_radio_options(true);
        assert_eq!(with.len(), 3);
        assert_eq!(with[2].0, FdrMethod::NoCorrection);
        // Verbatim label strings.
        assert_eq!(without[0].1, "Benjamini–Hochberg (BH) procedure");
        assert_eq!(without[1].1, "Benjamini–Yekutieli (BY) procedure");
        assert_eq!(with[2].1, "No correction (raw p-values)");
    }

    /// `png_export_size_controls` renders without panicking (smoke), reusing the
    /// `Context::default()` + `ctx.run` pattern. With no simulated drag input the
    /// three passed `&mut` scalars are left unchanged.
    #[test]
    fn png_export_size_controls_smoke() {
        let ctx = egui::Context::default();
        theme::install(&ctx);
        let mut w = 6.0_f64;
        let mut h = 4.0_f64;
        let mut dpi = 300_u32;
        let _ = ctx.run(egui::RawInput::default(), |ctx| {
            egui::CentralPanel::default().show(ctx, |ui| {
                png_export_size_controls(ui, &mut w, &mut h, &mut dpi);
            });
        });
        assert_eq!((w, h, dpi), (6.0, 4.0, 300));
    }
}
