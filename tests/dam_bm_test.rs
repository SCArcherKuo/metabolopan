//! Brunner-Munzel + Cliff's δ synthetic checks via the public API.

use metabolopan::dam::brunner_munzel::{brunner_munzel_two_tailed, cliffs_delta};

#[test]
fn bm_one_overlap_small_p() {
    let a = [1.0, 2.0, 3.0, 4.0, 5.0];
    let b = [4.5, 6.0, 7.0, 8.0, 9.0];
    let p = brunner_munzel_two_tailed(&a, &b);
    assert!(p.is_finite() && p < 0.01, "got {p}");
}

#[test]
fn bm_perfectly_stratified_nan() {
    let p = brunner_munzel_two_tailed(&[1.0, 2.0, 3.0, 4.0, 5.0], &[6.0, 7.0, 8.0, 9.0, 10.0]);
    assert!(p.is_nan());
}

#[test]
fn cliffs_delta_extremes() {
    assert_eq!(cliffs_delta(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]), -1.0);
    assert_eq!(cliffs_delta(&[4.0, 5.0, 6.0], &[1.0, 2.0, 3.0]), 1.0);
}
