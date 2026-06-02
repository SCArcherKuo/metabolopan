//! Student's t-test synthetic checks via the public API + cross-method
//! divergence regression test (guards against accidentally aliasing Welch and
//! Student).

use metabolopan::dam::student::student_t_two_tailed;
use metabolopan::dam::welch::welch_t_two_tailed;

#[test]
fn student_well_separated_small_p() {
    let a = [10.0, 12.0, 11.0, 13.0];
    let b = [1.0, 2.0, 1.5, 0.8];
    let p = student_t_two_tailed(&a, &b);
    assert!(p.is_finite() && p < 0.01, "expected p < 0.01, got {p}");
}

#[test]
fn student_nan_aware() {
    let a = [10.0, f64::NAN, 11.0, 13.0];
    let b = [1.0, 2.0, f64::NAN, 0.8];
    let p = student_t_two_tailed(&a, &b);
    assert!(
        p.is_finite() && p < 0.05,
        "expected NaN-aware Student to still give small p, got {p}"
    );
}

#[test]
fn student_too_few_samples_nan() {
    assert!(student_t_two_tailed(&[1.0], &[2.0, 3.0, 4.0]).is_nan());
}

/// Cross-method divergence: on a deliberately unequal-variance input both
/// tests must produce finite p-values, but Welch (Satterthwaite df, unequal
/// variances) and Student (pooled variance, df = n+m-2) must NOT collapse to
/// the same number. Guards against future refactors silently aliasing one
/// method to the other.
#[test]
fn student_and_welch_diverge_on_unequal_variances() {
    let a = [1.0, 2.0, 3.0];
    let b = [10.0, 11.0, 30.0];
    let p_welch = welch_t_two_tailed(&a, &b);
    let p_student = student_t_two_tailed(&a, &b);
    assert!(
        p_welch.is_finite(),
        "expected finite Welch p, got {p_welch}"
    );
    assert!(
        p_student.is_finite(),
        "expected finite Student p, got {p_student}"
    );
    assert!(
        (p_welch - p_student).abs() > 1e-3,
        "Welch and Student should diverge on unequal variances, \
         got p_welch={p_welch}, p_student={p_student}"
    );
}
