//! Integration-level FDR check (the unit tests in src/dam/fdr.rs cover the golden
//! vector; this test only confirms the symbols are reachable from the public API,
//! covering both the legacy `benjamini_yekutieli` entry point and the new
//! `adjust_pvalues(p, method)` dispatcher.)

use metabolopan::dam::fdr::{FdrMethod, adjust_pvalues, benjamini_hochberg, benjamini_yekutieli};

#[test]
fn by_golden_vector_matches_r() {
    // R: p.adjust(c(0.001, 0.01, 0.03, 0.05, 0.20), method='BY')
    //  = c(0.011416, 0.057082, 0.114164, 0.142705, 0.456657)
    let q = benjamini_yekutieli(&[0.001, 0.01, 0.03, 0.05, 0.20]);
    let expected = [0.011416, 0.057082, 0.114164, 0.142705, 0.456657];
    for (a, e) in q.iter().zip(expected.iter()) {
        assert!((a - e).abs() < 1e-4, "got {a}, expected {e}");
    }
}

#[test]
fn by_handles_all_one_input() {
    // All large p values should produce all 1.0 (capped).
    let q = benjamini_yekutieli(&[0.9, 0.95, 0.99]);
    assert!(q.iter().all(|v| (*v - 1.0).abs() < 1e-9));
}

#[test]
fn bh_golden_vector_matches_r() {
    // R: p.adjust(c(0.001, 0.01, 0.03, 0.05, 0.20), method='BH')
    //  = c(0.005, 0.025, 0.05, 0.0625, 0.20)
    let q = benjamini_hochberg(&[0.001, 0.01, 0.03, 0.05, 0.20]);
    let expected = [0.005, 0.025, 0.05, 0.0625, 0.20];
    for (a, e) in q.iter().zip(expected.iter()) {
        assert!((a - e).abs() < 1e-9, "got {a}, expected {e}");
    }
}

#[test]
fn dispatcher_via_public_api_matches_underlying_functions() {
    let p = [0.001, 0.01, 0.03, 0.05, 0.20];
    assert_eq!(
        adjust_pvalues(&p, FdrMethod::BenjaminiHochberg),
        benjamini_hochberg(&p)
    );
    assert_eq!(
        adjust_pvalues(&p, FdrMethod::BenjaminiYekutieli),
        benjamini_yekutieli(&p)
    );
}

#[test]
fn cross_method_differs_via_public_api() {
    // Guards against accidental dispatcher aliasing at the public API surface —
    // at least one position must differ by > 1e-12. (BY's c(m) factor inflates
    // the q-values relative to BH on this input.)
    let p = [0.001, 0.01, 0.03, 0.05, 0.20];
    let bh = adjust_pvalues(&p, FdrMethod::BenjaminiHochberg);
    let by = adjust_pvalues(&p, FdrMethod::BenjaminiYekutieli);
    let any_differ = bh.iter().zip(by.iter()).any(|(a, b)| (a - b).abs() > 1e-12);
    assert!(
        any_differ,
        "BH and BY must produce different outputs through adjust_pvalues; bh={bh:?} by={by:?}"
    );
}
