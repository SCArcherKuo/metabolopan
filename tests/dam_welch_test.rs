//! Welch's t-test synthetic checks via the public API.

use metabolopan::dam::welch::welch_t_two_tailed;

#[test]
fn welch_well_separated_small_p() {
    let a = [10.0, 12.0, 11.0, 13.0];
    let b = [1.0, 2.0, 1.5, 0.8];
    let p = welch_t_two_tailed(&a, &b);
    assert!(p.is_finite() && p < 0.01, "expected p < 0.01, got {p}");
}

#[test]
fn welch_nan_aware() {
    let a = [10.0, f64::NAN, 11.0, 13.0];
    let b = [1.0, 2.0, f64::NAN, 0.8];
    let p = welch_t_two_tailed(&a, &b);
    assert!(
        p.is_finite() && p < 0.05,
        "expected NaN-aware Welch to still give small p, got {p}"
    );
}

#[test]
fn welch_too_few_samples_nan() {
    assert!(welch_t_two_tailed(&[1.0], &[2.0, 3.0, 4.0]).is_nan());
}
