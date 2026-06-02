//! Metabolopan colour palette + global egui Visuals install.
//!
//! Single source of truth for every UI colour decision in the binary. Token
//! contract + install timing are pinned by `openspec/specs/theme/spec.md`
//! (capability `theme`) and the component system by
//! `openspec/specs/interactive-component-styles/spec.md` — keep the specs and
//! this file in lockstep.
//!
//! The DEFAULT widget appearance installed here is the `ui-color.md` §2/§3
//! **Secondary** component (light fill + `SURFACE` border brightening to
//! `PRIMARY` on hover). **Primary** styling is NOT global — it is opted into
//! per call site by the helpers in `src/ui/widgets.rs` via scoped `Visuals`
//! overrides.
//!
//! All contrast ratios in the doc-comments are computed via WCAG 2.x relative
//! luminance against the actual `BACKGROUND` value `(235, 230, 218)`, NOT
//! against pure white. The `aa_contrast_against_background` test below enforces
//! two tiers: ≥ 4.5:1 (AA normal) for `HEADING`/`TEXT`/`TEXT_SECONDARY`/`ERROR`,
//! and ≥ 3:1 (AA-Large) for `WARNING`/`LINK`/`SUCCESS`.

use egui::{Color32, Stroke, ThemePreference, Visuals};

// ─── Palette tokens ───────────────────────────────────────────────────────
//
// Source — designer-curated water-ink (水墨) palette + the ui-color.md §1–§5
// Primary/Secondary component system + a secondary layout background. 20 public
// tokens.

/// `#EBE6DA` 宣紙白 — primary app background.
pub const BACKGROUND: Color32 = Color32::from_rgb(235, 230, 218);

/// `#CED8D9` 雲影灰 (Cloud Shadow) — secondary layout background (ui-color.md
/// §1.1 "Surface / 次要背景"). Fill of the page chrome panels (top stepper +
/// bottom Data/Log panel) so they read as a distinct surface stacked over
/// `BACKGROUND`. Decorative surface, NOT a text token — no AA requirement.
/// Distinct from `SURFACE` (`#A4C3C9`), which ui-color.md designates the default
/// Border, not the secondary background.
pub const BACKGROUND_SECONDARY: Color32 = Color32::from_rgb(206, 216, 217);

/// `#A4C3C9` 晨霧水藍 — default component border, Secondary-button ACTIVE fill,
/// Pill-tab hover fill, Secondary-progress indicator. ~1.5:1 vs `BACKGROUND`.
pub const SURFACE: Color32 = Color32::from_rgb(164, 195, 201);

/// `#48879F` 晴空藍 — Primary action: Primary-button/Pill default fill, hover
/// border on Secondary widgets, selection/focus stroke, Primary-progress
/// indicator. ~3.22:1 — sub-AA for body text; graphical / CTA use only.
pub const PRIMARY: Color32 = Color32::from_rgb(72, 135, 159);

/// `#396F85` 晴藍墨 — Primary-button HOVER fill (ui-color.md §2). Shares its
/// exact hex with `LINK` (deliberate dual-role: button-hover vs hyperlink).
pub const PRIMARY_HOVER: Color32 = Color32::from_rgb(57, 111, 133);

/// `#285263` 深潭藍 — Primary-button ACTIVE (pressed) fill (§2).
pub const PRIMARY_ACTIVE: Color32 = Color32::from_rgb(40, 82, 99);

/// `#FFFFFF` — text/mark ON Primary fills (Primary-button label, selected
/// Pill-tab label, Primary-dropdown selected item). White.
pub const ON_PRIMARY: Color32 = Color32::from_rgb(255, 255, 255);

/// `#8B735F` 枯木褐 — structural lines ONLY: `widgets.noninteractive.bg_stroke`
/// (Frame borders, `ui.separator()`, indent guides) AND `window_stroke` (modal
/// outlines). NEVER interactive-widget borders (those use `SURFACE`/`PRIMARY`)
/// and never text. ~3.58:1 — sub-AA acceptable for lines/icons.
pub const ACCENT: Color32 = Color32::from_rgb(139, 115, 95);

/// `#E8F0F2` 淡霧藍 — Secondary-component DEFAULT fill: idle Secondary
/// button/input/dropdown bg, unselected Secondary checkbox/radio bg, progress
/// track.
pub const FILL_SECONDARY: Color32 = Color32::from_rgb(232, 240, 242);

/// `#DCE7E9` 深霧藍 — Secondary-component HOVER/SELECTED fill: hovered Secondary
/// button/input, selected Secondary checkbox/radio bg, Secondary-dropdown
/// selected item bg.
pub const FILL_SECONDARY_HOVER: Color32 = Color32::from_rgb(220, 231, 233);

/// `#B5B2AB` 灰質 — disabled strong surface: disabled Primary-button bg,
/// disabled Secondary border, disabled-selected Pill bg.
pub const DISABLED: Color32 = Color32::from_rgb(181, 178, 171);

/// `#E5E3DF` 淺灰質 — disabled light fill: disabled Secondary button/input bg.
pub const DISABLED_FILL: Color32 = Color32::from_rgb(229, 227, 223);

/// `#1A1A19` 焦墨黑 (Charcoal Ink) — top-level headings (H1/H2) + extreme
/// emphasis. Deepest ink. Applied at `ui.heading()` call sites + `.strong()`
/// section headers via `theme::HEADING`.
pub const HEADING: Color32 = Color32::from_rgb(26, 26, 25);

/// `#2D2C2A` 沉墨灰 (Deep Ink) — body text (`override_text_color`), H3/H4,
/// general button/form labels and marks; log pane INFO. ~11.21:1.
pub const TEXT: Color32 = Color32::from_rgb(45, 44, 42);

/// `#5C5A56` 淡墨灰 (Washed Ink) — secondary copy, captions, **dates &
/// timestamps** (the `Cached … ago` line and PubChem/KEGG-conv/modules
/// fetched-date spans), Line-tab unselected label; log pane DEBUG. ~5.5:1.
/// Renamed from `TEXT_MUTED`.
pub const TEXT_SECONDARY: Color32 = Color32::from_rgb(92, 90, 86);

/// `#8F8D88` 霜灰色 (Frost Gray) — disabled-state text + input placeholders +
/// future (not-yet-reached) stepper-pill labels; log pane TRACE. Renamed from
/// `TEXT_DIM`. Intentionally sub-AA (WCAG disabled-text exemption).
pub const TEXT_DISABLED: Color32 = Color32::from_rgb(143, 141, 136);

/// `#A64B44` 印泥紅 (Seal Red) — semantic error + destructive/delete actions
/// (labels, modal title, log pane ERROR, `Visuals::error_fg_color`,
/// Primary-progress error indicator). AA normal.
pub const ERROR: Color32 = Color32::from_rgb(166, 75, 68);

/// `#9C6A21` 赭石橘 (Ochre Amber) — semantic warning / caution / confirm-needed
/// (⚠ icons, log pane WARN, `warn_fg_color`, unsaved-changes prompts,
/// Unassigned-group emphasis). Distinct from `ERROR`. AA-Large tier.
pub const WARNING: Color32 = Color32::from_rgb(156, 106, 33);

/// `#396F85` 晴藍墨 (Link Blue) — egui hyperlink widgets via
/// `Visuals::hyperlink_color`. Same hex as `PRIMARY_HOVER`. AA-Large tier.
pub const LINK: Color32 = Color32::from_rgb(57, 111, 133);

/// `#658A4E` 初葉綠 (Spring Leaf) — success / add actions / Primary-progress
/// completed indicator / picked-file basenames on Stage 1. AA-Large tier.
pub const SUCCESS: Color32 = Color32::from_rgb(101, 138, 78);

// ─── Designer-pinned blend values (used inside `install` only) ────────────
//
// Off-white surface tints with no semantic role outside the Visuals install.

const FAINT_BG: Color32 = Color32::from_rgb(214, 218, 215);
const EXTREME_BG: Color32 = Color32::from_rgb(244, 240, 232);

// ─── Install ──────────────────────────────────────────────────────────────

/// Button padding: widened HORIZONTALLY (12 px, vs egui's default 4) so buttons
/// aren't cramped, but the VERTICAL pad stays at egui's default (1 px) on
/// purpose. egui seeds a horizontal row's cross-axis band to
/// `spacing.interact_size.y` (= 18), and a button taller than that band
/// overflows it and ends up vertically staggered against same-row labels. The
/// only way to grow the band is to raise `interact_size.y` — which ALSO raises
/// radio / checkbox row height and spreads them. So with both "rows centre-align"
/// and "radios keep their tight spacing" required, buttons may grow wider but
/// not taller than the ~18 px band. (`button_padding.y = 1` keeps a button at the
/// 18 px `interact_size.y` floor — the exact height that aligned cleanly before.)
pub const BUTTON_PADDING: egui::Vec2 = egui::vec2(12.0, 1.0);

/// Install the metabolopan palette as the global egui `Visuals` + wider button
/// padding.
///
/// MUST be called once from inside the `eframe::run_native` creation-context
/// callback in `src/main.rs`, BEFORE `App::new` returns.
///
/// Order: pin Light mode → build from `Visuals::light()` → override the
/// documented fields so the DEFAULT widget look is the §2/§3 Secondary
/// component → `ctx.set_visuals(...)` → widen `spacing.button_padding.x` so
/// buttons aren't cramped (leaving `.y` + `interact_size` at egui defaults, so
/// same-row labels stay vertically centred and radios keep their tight spacing).
pub fn install(ctx: &egui::Context) {
    ctx.set_theme(ThemePreference::Light);
    let visuals = build_palette_visuals();
    ctx.set_visuals(visuals);
    // Spacing lives on `Style`, not `Visuals`, so this is applied after
    // `set_visuals` (which only touches visuals) and does not affect the
    // field-by-field Visuals baseline test. ONLY `button_padding` — NOT
    // `interact_size.y`, which would raise radio/checkbox row height AND the
    // horizontal row band (mis-centring same-row labels against taller buttons).
    ctx.style_mut(|style| {
        style.spacing.button_padding = BUTTON_PADDING;
    });
}

/// Build the `Visuals` value `install` writes. Extracted so the test module
/// can construct an `expected_palette_visuals` baseline for full-field
/// `PartialEq` comparison.
fn build_palette_visuals() -> Visuals {
    let mut v = Visuals::light();

    // ── Fills ──
    v.panel_fill = BACKGROUND;
    v.window_fill = BACKGROUND;
    v.faint_bg_color = FAINT_BG;
    v.extreme_bg_color = EXTREME_BG;
    v.code_bg_color = FAINT_BG;

    // ── Widget state: bg_fill (checkbox / radio mark backgrounds — §3.2
    //    Secondary: selection fills the muted 深霧藍, NEVER PRIMARY) ──
    v.widgets.noninteractive.bg_fill = BACKGROUND;
    v.widgets.inactive.bg_fill = FILL_SECONDARY;
    v.widgets.hovered.bg_fill = FILL_SECONDARY_HOVER;
    v.widgets.active.bg_fill = FILL_SECONDARY_HOVER;
    v.widgets.open.bg_fill = FILL_SECONDARY_HOVER;

    // ── Widget state: weak_bg_fill (button / dropdown-trigger backgrounds —
    //    §2/§3.3 Secondary: 淡霧藍 idle → 深霧藍 hover → 晨霧水藍 pressed) ──
    v.widgets.noninteractive.weak_bg_fill = BACKGROUND;
    v.widgets.inactive.weak_bg_fill = FILL_SECONDARY;
    v.widgets.hovered.weak_bg_fill = FILL_SECONDARY_HOVER;
    v.widgets.active.weak_bg_fill = SURFACE;
    v.widgets.open.weak_bg_fill = FILL_SECONDARY_HOVER;

    // ── Widget state: strokes ──
    //
    // `noninteractive.bg_stroke` stays ACCENT — it paints Frame card borders,
    // `ui.separator()` lines, AND CollapsingHeader indent guides (the unified
    // "ACCENT line vocabulary", design D4), kept distinct from interactive
    // borders. The interactive `bg_stroke`s realize the §2/§3 Secondary 1px
    // border: 晨霧水藍 idle, brightening to 晴空藍 on hover / press.
    v.widgets.noninteractive.bg_stroke = Stroke::new(1.0, ACCENT);
    v.widgets.inactive.bg_stroke = Stroke::new(1.0, SURFACE);
    v.widgets.hovered.bg_stroke = Stroke::new(1.0, PRIMARY);
    v.widgets.active.bg_stroke = Stroke::new(1.0, PRIMARY);
    v.widgets.noninteractive.fg_stroke = Stroke::new(1.0, TEXT);
    v.widgets.inactive.fg_stroke = Stroke::new(1.0, TEXT);
    v.widgets.hovered.fg_stroke = Stroke::new(1.0, TEXT);
    // `active.fg_stroke = TEXT` (dark mark/label) so pressed Secondary buttons
    // and checked checkboxes/radios read with the §2/§3 Secondary dark `#2D2C2A`
    // mark. Primary widgets get their white label from the helper's
    // `RichText::color(ON_PRIMARY)`, NOT from this field.
    v.widgets.active.fg_stroke = Stroke::new(1.0, TEXT);

    // ── Selection / focus ──
    //
    // MUST use `from_rgba_unmultiplied` — egui's ecolor crate gamma-decodes the
    // RGB, multiplies by the linear alpha, and gamma-encodes back to a
    // premultiplied sRGB triple. Passing raw bytes via `from_rgba_premultiplied`
    // would skip the linear-space pre-darkening and render an over-saturated
    // halo. The `pill_tab` helper locally overrides `selection.bg_fill` to
    // OPAQUE `PRIMARY` inside its scope; this translucent value is the
    // text-selection halo.
    v.selection.bg_fill = Color32::from_rgba_unmultiplied(
        PRIMARY.r(),
        PRIMARY.g(),
        PRIMARY.b(),
        96, // ~37.6% opacity
    );
    v.selection.stroke = Stroke::new(1.0, PRIMARY);

    // ── Semantic / global text ──
    v.hyperlink_color = LINK;
    v.override_text_color = Some(TEXT);
    v.warn_fg_color = WARNING;
    v.error_fg_color = ERROR;

    // ── Window chrome (modals reuse this — no new modal tokens) ──
    v.window_stroke = Stroke::new(1.0, ACCENT);

    v
}

// ─── Tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// Independently-built `Visuals` baseline mirroring the `theme/spec.md`
    /// override list field-by-field. INTENTIONALLY duplicates the field writes
    /// from `build_palette_visuals` — that duplication is the drift detector:
    /// dropping a single line from `build_palette_visuals` trips
    /// `install_writes_full_palette_via_expected_baseline` via `PartialEq`.
    fn expected_palette_visuals() -> Visuals {
        let mut v = Visuals::light();
        // Fills
        v.panel_fill = BACKGROUND;
        v.window_fill = BACKGROUND;
        v.faint_bg_color = Color32::from_rgb(214, 218, 215);
        v.extreme_bg_color = Color32::from_rgb(244, 240, 232);
        v.code_bg_color = Color32::from_rgb(214, 218, 215);
        // bg_fill per widget state (Secondary checkbox/radio)
        v.widgets.noninteractive.bg_fill = BACKGROUND;
        v.widgets.inactive.bg_fill = FILL_SECONDARY;
        v.widgets.hovered.bg_fill = FILL_SECONDARY_HOVER;
        v.widgets.active.bg_fill = FILL_SECONDARY_HOVER;
        v.widgets.open.bg_fill = FILL_SECONDARY_HOVER;
        // weak_bg_fill per widget state (Secondary button/dropdown)
        v.widgets.noninteractive.weak_bg_fill = BACKGROUND;
        v.widgets.inactive.weak_bg_fill = FILL_SECONDARY;
        v.widgets.hovered.weak_bg_fill = FILL_SECONDARY_HOVER;
        v.widgets.active.weak_bg_fill = SURFACE;
        v.widgets.open.weak_bg_fill = FILL_SECONDARY_HOVER;
        // Strokes
        v.widgets.noninteractive.bg_stroke = Stroke::new(1.0, ACCENT);
        v.widgets.inactive.bg_stroke = Stroke::new(1.0, SURFACE);
        v.widgets.hovered.bg_stroke = Stroke::new(1.0, PRIMARY);
        v.widgets.active.bg_stroke = Stroke::new(1.0, PRIMARY);
        v.widgets.noninteractive.fg_stroke = Stroke::new(1.0, TEXT);
        v.widgets.inactive.fg_stroke = Stroke::new(1.0, TEXT);
        v.widgets.hovered.fg_stroke = Stroke::new(1.0, TEXT);
        v.widgets.active.fg_stroke = Stroke::new(1.0, TEXT);
        // Selection / focus
        v.selection.bg_fill =
            Color32::from_rgba_unmultiplied(PRIMARY.r(), PRIMARY.g(), PRIMARY.b(), 96);
        v.selection.stroke = Stroke::new(1.0, PRIMARY);
        // Semantic / global text
        v.hyperlink_color = LINK;
        v.override_text_color = Some(TEXT);
        v.warn_fg_color = WARNING;
        v.error_fg_color = ERROR;
        // Window chrome
        v.window_stroke = Stroke::new(1.0, ACCENT);
        v
    }

    /// (a) Every public token has full alpha.
    #[test]
    fn all_tokens_have_full_alpha() {
        for c in [
            BACKGROUND,
            BACKGROUND_SECONDARY,
            SURFACE,
            PRIMARY,
            PRIMARY_HOVER,
            PRIMARY_ACTIVE,
            ON_PRIMARY,
            ACCENT,
            FILL_SECONDARY,
            FILL_SECONDARY_HOVER,
            DISABLED,
            DISABLED_FILL,
            HEADING,
            TEXT,
            TEXT_SECONDARY,
            TEXT_DISABLED,
            ERROR,
            WARNING,
            LINK,
            SUCCESS,
        ] {
            assert_eq!(c.a(), 255, "token {c:?} must have full alpha");
        }
    }

    /// All 20 public tokens pin to the EXACT RGB triples documented in the
    /// `theme/spec.md` token table.
    #[test]
    fn all_tokens_have_documented_rgb_triples() {
        assert_eq!(BACKGROUND, Color32::from_rgb(235, 230, 218));
        assert_eq!(BACKGROUND_SECONDARY, Color32::from_rgb(206, 216, 217));
        assert_eq!(SURFACE, Color32::from_rgb(164, 195, 201));
        assert_eq!(PRIMARY, Color32::from_rgb(72, 135, 159));
        assert_eq!(PRIMARY_HOVER, Color32::from_rgb(57, 111, 133));
        assert_eq!(PRIMARY_ACTIVE, Color32::from_rgb(40, 82, 99));
        assert_eq!(ON_PRIMARY, Color32::from_rgb(255, 255, 255));
        assert_eq!(ACCENT, Color32::from_rgb(139, 115, 95));
        assert_eq!(FILL_SECONDARY, Color32::from_rgb(232, 240, 242));
        assert_eq!(FILL_SECONDARY_HOVER, Color32::from_rgb(220, 231, 233));
        assert_eq!(DISABLED, Color32::from_rgb(181, 178, 171));
        assert_eq!(DISABLED_FILL, Color32::from_rgb(229, 227, 223));
        assert_eq!(HEADING, Color32::from_rgb(26, 26, 25));
        assert_eq!(TEXT, Color32::from_rgb(45, 44, 42));
        assert_eq!(TEXT_SECONDARY, Color32::from_rgb(92, 90, 86));
        assert_eq!(TEXT_DISABLED, Color32::from_rgb(143, 141, 136));
        assert_eq!(ERROR, Color32::from_rgb(166, 75, 68));
        assert_eq!(WARNING, Color32::from_rgb(156, 106, 33));
        assert_eq!(LINK, Color32::from_rgb(57, 111, 133));
        assert_eq!(SUCCESS, Color32::from_rgb(101, 138, 78));
    }

    /// `PRIMARY_HOVER` and `LINK` intentionally share `#396F85` (button-hover
    /// vs hyperlink). Pin the share so a future tweak to one without the other
    /// trips loudly rather than silently diverging.
    #[test]
    fn primary_hover_equals_link() {
        assert_eq!(PRIMARY_HOVER, LINK);
    }

    /// (b1) `install` writes the full palette field-by-field — proven via
    /// `PartialEq` against an INDEPENDENTLY-built baseline.
    #[test]
    fn install_writes_full_palette_via_expected_baseline() {
        let ctx = egui::Context::default();
        install(&ctx);
        let expected = expected_palette_visuals();
        assert_eq!(
            ctx.style().visuals,
            expected,
            "post-install Visuals must match the independently-built expected_palette_visuals \
             field-by-field per theme/spec.md"
        );
    }

    /// (b2) `install` is idempotent under intervening mutation.
    #[test]
    fn install_overwrites_intervening_mutation() {
        let ctx = egui::Context::default();
        install(&ctx);
        let expected = expected_palette_visuals();

        ctx.set_visuals(Visuals::default());
        assert_ne!(
            ctx.style().visuals,
            expected,
            "sanity check: intervening mutation must actually differ from palette"
        );

        install(&ctx);
        assert_eq!(
            ctx.style().visuals,
            expected,
            "second install must reset visuals to canonical palette state"
        );
    }

    /// (c) Designer-pinned blend RGBs + the Secondary-default fills landed
    /// verbatim in the install output.
    #[test]
    fn install_writes_pinned_blend_rgbs() {
        let ctx = egui::Context::default();
        install(&ctx);
        let v = &ctx.style().visuals;

        assert_eq!(v.faint_bg_color, Color32::from_rgb(214, 218, 215));
        assert_eq!(v.extreme_bg_color, Color32::from_rgb(244, 240, 232));
        assert_eq!(v.code_bg_color, Color32::from_rgb(214, 218, 215));
        assert_eq!(v.widgets.inactive.bg_fill, FILL_SECONDARY);
        assert_eq!(v.widgets.inactive.weak_bg_fill, FILL_SECONDARY);
    }

    /// (d) `install` pins the egui context's theme preference to Light.
    #[test]
    fn install_pins_light_theme_preference() {
        let ctx = egui::Context::default();
        install(&ctx);
        let pref = ctx.options(|o| o.theme_preference);
        assert_eq!(pref, ThemePreference::Light);
    }

    /// (e) `install` widens `button_padding` so buttons are not cramped
    /// (`Style.spacing`, separate from the Visuals baseline above) WITHOUT
    /// raising `interact_size.y` — so radios / checkboxes keep egui's default
    /// (tight) vertical spacing.
    #[test]
    fn install_applies_roomier_button_padding_only() {
        let ctx = egui::Context::default();
        install(&ctx);
        let spacing = ctx.style().spacing.clone();
        assert_eq!(spacing.button_padding, BUTTON_PADDING);
        // interact_size.y is left at egui's default (NOT raised), so stacked
        // radios/checkboxes are not spread apart.
        assert_eq!(
            spacing.interact_size.y,
            egui::Style::default().spacing.interact_size.y
        );
    }

    /// Spot-check the palette-critical Visuals fields, incl. the Secondary
    /// default widget surfaces. The full-field PartialEq above covers
    /// everything; these give a regression an obvious named diagnostic.
    #[test]
    fn install_writes_palette_critical_fields() {
        let ctx = egui::Context::default();
        install(&ctx);
        let v = &ctx.style().visuals;
        assert!(!v.dark_mode);
        assert_eq!(v.panel_fill, BACKGROUND);
        assert_eq!(v.window_fill, BACKGROUND);
        assert_eq!(v.override_text_color, Some(TEXT));
        assert_eq!(v.hyperlink_color, LINK);
        assert_eq!(v.warn_fg_color, WARNING);
        assert_eq!(v.error_fg_color, ERROR);
        assert_eq!(v.window_stroke.color, ACCENT);
        assert_eq!(v.widgets.noninteractive.bg_stroke.color, ACCENT);
        // Secondary button surface (weak_bg_fill) across the interactive states.
        assert_eq!(v.widgets.noninteractive.weak_bg_fill, BACKGROUND);
        assert_eq!(v.widgets.inactive.weak_bg_fill, FILL_SECONDARY);
        assert_eq!(v.widgets.hovered.weak_bg_fill, FILL_SECONDARY_HOVER);
        assert_eq!(v.widgets.active.weak_bg_fill, SURFACE);
        assert_eq!(v.widgets.open.weak_bg_fill, FILL_SECONDARY_HOVER);
        // Secondary 1px border idle → Primary on hover/press.
        assert_eq!(v.widgets.inactive.bg_stroke.color, SURFACE);
        assert_eq!(v.widgets.hovered.bg_stroke.color, PRIMARY);
        assert_eq!(v.widgets.active.bg_stroke.color, PRIMARY);
        // Dark mark/label on pressed Secondary widgets + checked boxes.
        assert_eq!(v.widgets.active.fg_stroke.color, TEXT);
    }

    // ─── AA contrast ─────────────────────────────────────────────────────

    /// WCAG 2.x relative luminance via `f64::powf`.
    fn wcag_relative_luminance(c: Color32) -> f64 {
        fn channel(c: u8) -> f64 {
            let v = c as f64 / 255.0;
            if v <= 0.040_45 {
                v / 12.92
            } else {
                ((v + 0.055) / 1.055).powf(2.4)
            }
        }
        0.2126 * channel(c.r()) + 0.7152 * channel(c.g()) + 0.0722 * channel(c.b())
    }

    fn contrast_ratio(a: Color32, b: Color32) -> f64 {
        let la = wcag_relative_luminance(a);
        let lb = wcag_relative_luminance(b);
        let (lo, hi) = if la < lb { (la, lb) } else { (lb, la) };
        (hi + 0.05) / (lo + 0.05)
    }

    /// AA gate vs `BACKGROUND`. Tier 1 ≥ 4.5:1 (headings/body/secondary/error);
    /// Tier 2 ≥ 3:1 (warning/link/success). `TEXT_DISABLED` / `PRIMARY` /
    /// `ACCENT` + the surface/fill component tokens are excluded — sub-AA,
    /// disabled, decorative, or non-text by design.
    #[test]
    fn aa_contrast_against_background() {
        // Tier 1 — AA normal (4.5:1).
        for (name, token) in [
            ("HEADING", HEADING),
            ("TEXT", TEXT),
            ("TEXT_SECONDARY", TEXT_SECONDARY),
            ("ERROR", ERROR),
        ] {
            let ratio = contrast_ratio(token, BACKGROUND);
            assert!(
                ratio >= 4.5,
                "{name} must clear AA 4.5:1 against BACKGROUND — got {ratio:.4}"
            );
        }
        // Tier 2 — AA-Large (3:1). Designer-pinned water-ink hues sit just below
        // 4.5:1 on the warm-light background (WARNING ≈ 3.75:1) — documented.
        for (name, token) in [("WARNING", WARNING), ("LINK", LINK), ("SUCCESS", SUCCESS)] {
            let ratio = contrast_ratio(token, BACKGROUND);
            assert!(
                ratio >= 3.0,
                "{name} must clear AA-Large 3:1 against BACKGROUND — got {ratio:.4}"
            );
        }
    }
}
