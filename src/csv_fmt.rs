//! The single shared CSV float formatter used by BOTH the DAM exporter
//! (`dam::export`) and the enrichment exporter (`enrichment::export`). Kept in
//! one place so the two exporters can never drift on how they render `NaN` /
//! `±∞` / finite floats into CSV cells.

/// Format an `f64` for a CSV cell: `NaN → ""` (empty), `+∞ → "inf"`,
/// `-∞ → "-inf"`, finite → the default `{v}` rendering.
pub fn fmt_csv_f64(v: f64) -> String {
    if v.is_nan() {
        String::new()
    } else if v == f64::INFINITY {
        "inf".to_string()
    } else if v == f64::NEG_INFINITY {
        "-inf".to_string()
    } else {
        format!("{v}")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn formats_all_four_branches() {
        assert_eq!(fmt_csv_f64(f64::NAN), "");
        assert_eq!(fmt_csv_f64(f64::INFINITY), "inf");
        assert_eq!(fmt_csv_f64(f64::NEG_INFINITY), "-inf");
        assert_eq!(fmt_csv_f64(0.5), "0.5");
    }
}
