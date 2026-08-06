//! The pre-stepper `Choose your analysis` screen — the first thing the user
//! sees once the organism roster has loaded.
//!
//! It is **not a numbered stage**: the heading carries no `Stage 0 —` prefix and
//! the stepper is not rendered here (there is nothing to navigate back to yet).
//! Its whole job is to write `app.settings.analysis_route` and move on, which is
//! why picking a card advances immediately rather than arming a separate
//! confirm button.
//!
//! Owner: the `analysis-route-selection` and `interactive-component-styles`
//! capabilities.

use egui::RichText;

use crate::app::{AnalysisRoute, App, AppState};
use crate::theme;
use crate::ui::widgets::primary_button_sized;

/// One route option: the clickable title plus the two supporting lines.
struct RouteCard {
    route: AnalysisRoute,
    title: &'static str,
    summary: &'static str,
    requirements: &'static str,
    /// Extra grey note rendered under the requirements line. Only the coverage
    /// card has one.
    hint: Option<&'static str>,
}

/// Height of a route-card button. Taller than a default button so the two CTAs
/// carry the weight of the only decision on the screen.
const BUTTON_HEIGHT: f32 = 32.0;

const CARDS: [RouteCard; 2] = [
    RouteCard {
        route: AnalysisRoute::DamEnrichment,
        title: "Differential analysis + enrichment",
        summary: "Compare two sample groups, then test which KEGG pathways or modules are \
                  enriched among the significantly changed metabolites.",
        requirements: "Needs a group .csv with at least 2 groups of at least 2 samples.",
        hint: None,
    },
    RouteCard {
        route: AnalysisRoute::KeggCoverage,
        title: "KEGG coverage survey",
        summary: "Map every detected metabolite onto KEGG pathways or modules and report how \
                  completely each one is covered. No statistical test.",
        requirements: "Needs only an MS-DIAL .txt. A group .csv is optional.",
        hint: Some(
            "The PubChem and KEGG caches are shared between the two analyses, so running a \
             coverage survey first makes a later enrichment run faster.",
        ),
    },
];

pub fn show(ui: &mut egui::Ui, app: &mut App) {
    egui::ScrollArea::both()
        .auto_shrink([false, false])
        .show(ui, |ui| {
            // No stage number and no em-dash stage prefix: this screen sits
            // outside the numbered pipeline.
            ui.heading(RichText::new("Choose your analysis").color(theme::HEADING));
            ui.add_space(12.0);

            let mut picked: Option<AnalysisRoute> = None;
            for (i, card) in CARDS.iter().enumerate() {
                if i > 0 {
                    ui.add_space(10.0);
                }
                if render_card(ui, card) {
                    picked = Some(card.route);
                }
            }

            if let Some(route) = picked {
                choose_route(app, route);
            }
        });
}

/// Render one route card; `true` when the user clicked its button.
///
/// **The title is a §2 `primary_button`, not a §4 segmented tab.** A segmented
/// control expresses *the currently selected option among several* — and this
/// screen has no currently selected option, because it is a prompt and
/// pre-selecting one would bias the choice. Both segments therefore rendered in
/// the unselected state, which is transparent with no border: nothing on the
/// card said "this is clickable" until the pointer was already on it. Two
/// filled Primary CTAs say it at rest, which is what a screen whose only job is
/// "pick one of these two actions" needs.
///
/// This is the one screen with TWO Primary CTAs. The usual "one forward CTA per
/// screen" rule exists so the eye knows where the single next step is; here
/// there are genuinely two next steps and neither is the default.
///
/// The button is stretched to the card width so the two read as equals — a
/// button auto-sized to its label would make `KEGG coverage survey` visibly
/// smaller than `Differential analysis + enrichment` and imply a hierarchy the
/// screen does not have.
fn render_card(ui: &mut egui::Ui, card: &RouteCard) -> bool {
    egui::Frame::NONE
        .fill(theme::ON_PRIMARY)
        .stroke(egui::Stroke::new(1.0_f32, theme::SURFACE))
        .inner_margin(egui::Margin::symmetric(10, 8))
        .corner_radius(egui::CornerRadius::same(6))
        .show(ui, |ui| {
            let width = ui.available_width();
            let clicked =
                primary_button_sized(ui, card.title, true, egui::vec2(width, BUTTON_HEIGHT))
                    .clicked();
            ui.add_space(4.0);
            ui.label(RichText::new(card.summary).color(theme::TEXT_SECONDARY));
            ui.label(RichText::new(card.requirements).color(theme::TEXT_SECONDARY));
            if let Some(hint) = card.hint {
                ui.add_space(2.0);
                ui.label(RichText::new(hint).small().color(theme::TEXT_SECONDARY));
            }
            clicked
        })
        .inner
}

/// Record the chosen route and advance to Stage 1.
///
/// This is one of exactly three paths that may write `analysis_route` (the
/// others being the Stage 1 `Change analysis type` escape and a loaded
/// snapshot). Nothing else is touched: `settings`, `inputs`, and `cache` carry
/// over untouched, so a user who takes the Stage 1 escape and re-picks the SAME
/// card lands back exactly where they left.
fn choose_route(app: &mut App, route: AnalysisRoute) {
    app.settings.analysis_route = route;
    tracing::info!(?route, "analysis route chosen");
    app.state = AppState::Stage1Input {
        slot1_mode: None,
        slot2_revealed: false,
        slot2_mode: None,
        error: None,
    };
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Both routes are offered, exactly once each, in the documented order.
    #[test]
    fn both_routes_are_offered_once_each() {
        assert_eq!(CARDS.len(), 2);
        assert_eq!(CARDS[0].route, AnalysisRoute::DamEnrichment);
        assert_eq!(CARDS[1].route, AnalysisRoute::KeggCoverage);
    }

    /// The shared-cache hint belongs to the coverage card only — it is the
    /// argument for trying the cheaper route first, which makes no sense on the
    /// enrichment card.
    #[test]
    fn only_the_coverage_card_carries_the_shared_cache_hint() {
        assert!(CARDS[0].hint.is_none());
        let hint = CARDS[1].hint.expect("coverage card carries the hint");
        assert!(hint.contains("PubChem"));
        assert!(hint.contains("KEGG"));
    }

    /// The card copy is user-facing contract text, pinned so an edit is a
    /// deliberate act rather than a drive-by rewording.
    #[test]
    fn card_copy_matches_the_spec() {
        assert_eq!(CARDS[0].title, "Differential analysis + enrichment");
        assert_eq!(
            CARDS[0].requirements,
            "Needs a group .csv with at least 2 groups of at least 2 samples."
        );
        assert_eq!(CARDS[1].title, "KEGG coverage survey");
        assert_eq!(
            CARDS[1].requirements,
            "Needs only an MS-DIAL .txt. A group .csv is optional."
        );
        assert!(CARDS[1].summary.ends_with("No statistical test."));
    }
}
