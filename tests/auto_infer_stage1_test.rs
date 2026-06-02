//! Composition test for `auto-infer-stage1-ion-mode` §4.1.
//!
//! The full end-to-end path is `parse_msdial_txt → infer_polarity →
//! decide_slot1_mode_on_file_load → AppState::Stage1Input.slot1_mode`. The
//! helpers (`decide_slot1_mode_on_file_load` / `decide_slot2_mode_on_slot1_change`
//! / `opposite_mode`) are private to `src/ui/stage1_input.rs` so we cannot call
//! them directly from an integration test (per task 4.1's fallback clause).
//!
//! What this test covers instead — the two non-private halves of the
//! composition, run on the actual `data/double-mode/` fixtures, so a future
//! regression in either half (parser changes, polarity-detection threshold
//! drift, fixture replacement) is caught with a clear error message rooted at
//! the real file paths used by manual smoke test §5.4:
//!
//! 1. `parse_msdial_txt(POS_fixture) → table` succeeds and `infer_polarity(table)`
//!    returns `AdductPolarityInference::Positive`.
//! 2. Same for the NEG fixture → `Negative`.
//!
//! The `decide_slot1_mode_on_file_load` helper is a single deterministic match
//! over the three `AdductPolarityInference` variants (unit-tested exhaustively
//! in `src/ui/stage1_input.rs::tests::decide_slot1_mode_on_file_load_maps_each_variant`),
//! so the composition `parse → infer → helper → slot1_mode` is closed by these
//! two halves plus the unit test. The actual `AppState` mutation in
//! `apply_picked_file_to_slot` is exercised by the manual smoke tests §5.1–§5.6.
//!
//! NB: `tests/data_msdial_test.rs::infers_polarity_from_bundled_double_mode_{pos,neg}_fixture`
//! already covers the same two halves; this file restates the assertion in the
//! context of the auto-fill feature so a future change archive can grep for
//! `auto-infer-stage1-ion-mode` and find every test that protects the wiring.

use metabolopan::data::{AdductPolarityInference, infer_polarity, parse_msdial_txt};
use std::path::Path;

#[test]
fn pos_fixture_infers_positive_so_slot1_mode_auto_fills_positive() {
    let table = parse_msdial_txt(Path::new("data/double-mode/data-positive.txt"))
        .expect("bundled POS fixture must parse");
    let inferred = infer_polarity(&table);
    assert_eq!(
        inferred,
        AdductPolarityInference::Positive,
        "POS fixture must infer Positive — auto-fill would write slot1_mode = Some(Positive)"
    );
}

#[test]
fn neg_fixture_infers_negative_so_slot1_mode_auto_fills_negative() {
    let table = parse_msdial_txt(Path::new("data/double-mode/data-negative.txt"))
        .expect("bundled NEG fixture must parse");
    let inferred = infer_polarity(&table);
    assert_eq!(
        inferred,
        AdductPolarityInference::Negative,
        "NEG fixture must infer Negative — auto-fill would write slot1_mode = Some(Negative); \
         slot2_mode then auto-fills to Some(Positive) via D3 trigger #3 / opposite_mode"
    );
}
