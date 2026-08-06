//! MS-DIAL adduct classification for InChIKey deduplication.
//!
//! Pure functions, no I/O, no `tracing` events (per the
//! `msdial-deduplication` pure-function contract). Classifies an
//! `Adduct type` string into one of four classes that the cascade uses
//! at level 2:
//!
//! | Rank | Class       | Examples                                            |
//! |------|-------------|-----------------------------------------------------|
//! |  0   | `Primary`   | `[M+H]+`, `[M+Na]+`, `[M+NH4]+`, `[M+K]+`, `[M-H]-`, `[M+Cl]-` |
//! |  1   | `NonPrimary`| `[M+FA-H]-`, `[M+H-H2O]+`, missing / unknown        |
//! |  2   | `Dimer`     | `[2M+H]+`, `[2M-H]-`, `[3M-H]-` (any leading n>1)   |
//! |  3   | `Isotope`   | `[M+1]+`, `[M+2]-`, or `isotope_tracking_weight_number > 0` |
//!
//! Within `Primary`, `[M+H]+` and `[M-H]-` carry sub-rank 0 (preferred);
//! the other allowlist entries carry sub-rank 1. Owner: the `msdial-deduplication` capability.

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum AdductClass {
    Primary,
    NonPrimary,
    Dimer,
    Isotope,
}

impl AdductClass {
    /// Variant-name label for the dedup audit CSV's adduct-class cell.
    /// Exhaustive (no wildcard) so a new `AdductClass` variant forces a
    /// compile error HERE rather than a silent CSV mislabel
    /// (`move-labels-onto-types`). The `Primary:{sub}` sub-rank composition
    /// stays in the exporter — it joins a different field, not this label.
    pub fn label(&self) -> &'static str {
        match self {
            AdductClass::Primary => "Primary",
            AdductClass::NonPrimary => "NonPrimary",
            AdductClass::Dimer => "Dimer",
            AdductClass::Isotope => "Isotope",
        }
    }
}

const PRIMARY_ALLOWLIST: &[&str] = &[
    "[M+H]+", "[M+Na]+", "[M+NH4]+", "[M+K]+", "[M-H]-", "[M+Cl]-",
];

/// Classify an MS-DIAL `Adduct type` string given the companion
/// `Isotope tracking weight number` cell. Short-circuit order matches
/// the `msdial-deduplication` spec's "Adduct classification SHALL be a
/// pure function" requirement.
pub fn classify(adduct: Option<&str>, isotope_weight: Option<i32>) -> AdductClass {
    // (1) isotope weight > 0 always wins.
    if isotope_weight.unwrap_or(0) > 0 {
        return AdductClass::Isotope;
    }
    let Some(s) = adduct else {
        // (2) missing adduct string => NonPrimary (decision D5: missing
        // information is not equivalent to being a dimer or isotope).
        return AdductClass::NonPrimary;
    };
    // (3) isotope-string detection: matches `\[M\+\d+\]` (e.g. "[M+1]+", "[M+10]-").
    if matches_isotope_pattern(s) {
        return AdductClass::Isotope;
    }
    // (4) dimer/multimer detection: matches `\[(\d+)M` with captured n > 1.
    if leading_multimer_n(s).is_some_and(|n| n > 1) {
        return AdductClass::Dimer;
    }
    // (5) primary allowlist hit.
    if PRIMARY_ALLOWLIST.contains(&s) {
        return AdductClass::Primary;
    }
    // (6) everything else is NonPrimary.
    AdductClass::NonPrimary
}

/// Sub-rank within the `Primary` class. Returns `0` for `[M+H]+` and
/// `[M-H]-` (most-preferred primary adducts), `1` for any other primary
/// allowlist entry. Callers MUST invoke only when `classify` returned
/// `AdductClass::Primary`; on other inputs the return value is unspecified
/// (the function does not panic but its return value carries no meaning).
pub fn primary_subrank(adduct: &str) -> u8 {
    if adduct == "[M+H]+" || adduct == "[M-H]-" {
        0
    } else {
        1
    }
}

/// True iff the input matches `\[M\+\d+\]` somewhere (e.g. "[M+1]+",
/// "[M+2]-", "[M+10]+"). Byte-scan rather than the `regex` crate per
/// design D11 (no new dependency).
fn matches_isotope_pattern(s: &str) -> bool {
    // Look for the literal "[M+" then at least one ASCII digit then "]".
    let bytes = s.as_bytes();
    let mut i = 0;
    while i + 3 < bytes.len() {
        if &bytes[i..i + 3] == b"[M+" {
            let mut j = i + 3;
            let digit_start = j;
            while j < bytes.len() && bytes[j].is_ascii_digit() {
                j += 1;
            }
            if j > digit_start && j < bytes.len() && bytes[j] == b']' {
                return true;
            }
        }
        i += 1;
    }
    false
}

/// If `s` starts with `[<n>M` where `<n>` is an ASCII integer, return
/// `Some(n)`. Returns `None` for `[M…` (no leading digit), `[2X…`
/// (not "M" after digits), or empty / non-bracketed input.
fn leading_multimer_n(s: &str) -> Option<u32> {
    let bytes = s.as_bytes();
    if bytes.is_empty() || bytes[0] != b'[' {
        return None;
    }
    let mut j = 1;
    let digit_start = j;
    while j < bytes.len() && bytes[j].is_ascii_digit() {
        j += 1;
    }
    if j == digit_start || j >= bytes.len() || bytes[j] != b'M' {
        return None;
    }
    std::str::from_utf8(&bytes[digit_start..j])
        .ok()
        .and_then(|s| s.parse::<u32>().ok())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn adduct_class_label_matches_variant_names() {
        assert_eq!(AdductClass::Primary.label(), "Primary");
        assert_eq!(AdductClass::NonPrimary.label(), "NonPrimary");
        assert_eq!(AdductClass::Dimer.label(), "Dimer");
        assert_eq!(AdductClass::Isotope.label(), "Isotope");
    }

    // ── Primary class + sub-rank ──

    #[test]
    fn primary_m_plus_h_positive() {
        assert_eq!(classify(Some("[M+H]+"), None), AdductClass::Primary);
        assert_eq!(primary_subrank("[M+H]+"), 0);
    }

    #[test]
    fn primary_m_minus_h_negative() {
        assert_eq!(classify(Some("[M-H]-"), None), AdductClass::Primary);
        assert_eq!(primary_subrank("[M-H]-"), 0);
    }

    #[test]
    fn primary_sub_rank_one_for_m_plus_na() {
        assert_eq!(classify(Some("[M+Na]+"), None), AdductClass::Primary);
        assert_eq!(primary_subrank("[M+Na]+"), 1);
    }

    #[test]
    fn primary_sub_rank_one_for_nh4_k_cl() {
        for a in &["[M+NH4]+", "[M+K]+", "[M+Cl]-"] {
            assert_eq!(classify(Some(*a), None), AdductClass::Primary);
            assert_eq!(primary_subrank(a), 1);
        }
    }

    // ── NonPrimary ──

    #[test]
    fn nonprimary_m_fa_h_negative() {
        assert_eq!(classify(Some("[M+FA-H]-"), None), AdductClass::NonPrimary);
    }

    #[test]
    fn nonprimary_m_h_h2o_positive() {
        assert_eq!(classify(Some("[M+H-H2O]+"), None), AdductClass::NonPrimary);
    }

    // ── Dimer ──

    #[test]
    fn dimer_2m_plus_h_positive() {
        assert_eq!(classify(Some("[2M+H]+"), None), AdductClass::Dimer);
    }

    #[test]
    fn dimer_2m_minus_h_negative() {
        assert_eq!(classify(Some("[2M-H]-"), None), AdductClass::Dimer);
    }

    #[test]
    fn dimer_higher_multimer_3m_minus_h() {
        // Any n > 1 multimer classifies as Dimer.
        assert_eq!(classify(Some("[3M-H]-"), None), AdductClass::Dimer);
    }

    // ── Isotope ──

    #[test]
    fn isotope_via_weight_number_overrides_primary() {
        // isotope weight 1 wins even when the adduct string is primary.
        assert_eq!(classify(Some("[M+H]+"), Some(1)), AdductClass::Isotope);
    }

    #[test]
    fn isotope_via_adduct_string_m_plus_1() {
        assert_eq!(classify(Some("[M+1]+"), None), AdductClass::Isotope);
    }

    #[test]
    fn isotope_via_adduct_string_m_plus_2() {
        assert_eq!(classify(Some("[M+2]-"), None), AdductClass::Isotope);
    }

    // ── Edge cases ──

    #[test]
    fn missing_adduct_is_nonprimary_not_isotope() {
        // Decision D5: missing information is not equivalent to being a
        // dimer or isotope; a feature whose adduct cell was blank in the
        // input file might still be a valid primary adduct that someone
        // forgot to annotate.
        assert_eq!(classify(None, None), AdductClass::NonPrimary);
    }

    #[test]
    fn isotope_weight_overrides_dimer_pattern() {
        // [2M+H]+ would normally be Dimer, but isotope_weight=2 short-
        // circuits at level (1).
        assert_eq!(classify(Some("[2M+H]+"), Some(2)), AdductClass::Isotope);
    }

    #[test]
    fn isotope_weight_zero_does_not_override() {
        // isotope_weight=0 means M0 monoisotopic; primary classification
        // proceeds normally.
        assert_eq!(classify(Some("[M+H]+"), Some(0)), AdductClass::Primary);
    }

    // ── helper byte-scan correctness ──

    #[test]
    fn isotope_pattern_helper_matches_m_plus_double_digit() {
        // Anchored regex `\[M\+\d+\]` is satisfied by `[M+10]+` too.
        assert!(matches_isotope_pattern("[M+10]+"));
    }

    #[test]
    fn isotope_pattern_helper_rejects_non_isotope() {
        assert!(!matches_isotope_pattern("[M+H]+"));
        assert!(!matches_isotope_pattern("[2M+H]+"));
        assert!(!matches_isotope_pattern("[M-H]-"));
    }

    #[test]
    fn leading_multimer_extracts_n() {
        assert_eq!(leading_multimer_n("[2M+H]+"), Some(2));
        assert_eq!(leading_multimer_n("[3M-H]-"), Some(3));
        assert_eq!(leading_multimer_n("[10M+H]+"), Some(10));
        assert_eq!(leading_multimer_n("[M+H]+"), None);
        assert_eq!(leading_multimer_n("M+H"), None);
        assert_eq!(leading_multimer_n(""), None);
    }
}
