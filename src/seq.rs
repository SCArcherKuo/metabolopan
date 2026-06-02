//! Small sequence helpers shared across modules.

use std::collections::HashSet;
use std::hash::Hash;

/// Deduplicate `input` keeping the FIRST occurrence of each element and
/// preserving input order. Generic over any hashable, equatable, cloneable
/// element; both production callers (`kegg::conv`, `pubchem::resolve`) pass
/// `&[String]`. Extracted from two byte-identical copies by
/// `move-labels-onto-types`.
///
/// NOTE: this is distinct from the `sort(); dedup()` and count-only
/// deduplications in `stage3` — those are deliberately NOT order-preserving
/// and stay separate (see that change's design D4).
pub fn dedupe_preserve_order<T: Eq + Hash + Clone>(input: &[T]) -> Vec<T> {
    let mut seen = HashSet::new();
    let mut out = Vec::with_capacity(input.len());
    for k in input {
        if seen.insert(k.clone()) {
            out.push(k.clone());
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dedupe_preserves_first_occurrence_order() {
        let input = vec![
            "A".to_string(),
            "B".to_string(),
            "A".to_string(),
            "C".to_string(),
            "B".to_string(),
        ];
        let deduped = dedupe_preserve_order(&input);
        assert_eq!(
            deduped,
            vec!["A".to_string(), "B".to_string(), "C".to_string()]
        );
    }
}
