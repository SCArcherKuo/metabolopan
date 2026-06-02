//! Optional pre-test transform shared by the parametric DAM paths (Welch and Student).
//!
//! The user-controlled `Log transformation` checkbox (Stage 2 setup) gates
//! whether `arcsinh_in_place` runs on the union of numerator + denominator
//! values before Welch / Student. NaN cells pass through untouched.
//!
//! An earlier `pareto_scale_in_place` helper that ran unconditionally after
//! `arcsinh_in_place` was removed in `add-log-transform-and-scaling` once
//! empirical verification at f64 precision confirmed Pareto scaling is
//! bit-equivalent for univariate t-statistics (per-feature linear rescaling
//! cancels in `t = (m_a − m_b) / sqrt(v_a/na + v_b/nb)`). See design D1 of
//! that change for the verification.

/// Apply `arcsinh(x)` to every non-NaN cell in place; NaN cells are left untouched.
pub fn arcsinh_in_place(values: &mut [f64]) {
    for v in values.iter_mut() {
        if !v.is_nan() {
            *v = v.asinh();
        }
    }
}
