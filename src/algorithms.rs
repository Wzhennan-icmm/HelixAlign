//! Exact pairwise maximal matches, with optional reference/query uniqueness.
//!
//! Coordinates are zero based. A zero byte is a hard sequence boundary and is
//! never part of a match. Callers can join reference records with zero bytes;
//! queries must be processed one FASTA record and one orientation at a time.

use crate::suffix_array::{Match, SparseSuffixArray};

/// Match types corresponding to mummer -mum, -mumreference and -maxmatch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchType {
    MUM,
    MAM,
    MEM,
}

/// Find matches maximal at both ends and unique in both input sequences.
///
/// Query occurrences include overlapping occurrences. For a multi-FASTA query,
/// call this once per record: MUMmer tests uniqueness within that query record,
/// while reference uniqueness is across all records in the reference index.
/// A zero minimum length is invalid and returns an empty result.
pub fn find_mums(reference: &SparseSuffixArray, query: &[u8], min_len: usize) -> Vec<Match> {
    let candidates = find_mams(reference, query, min_len);
    if candidates.is_empty() {
        return candidates;
    }
    // A full query index counts every occurrence, including overlapping ones.
    let query_index = SparseSuffixArray::new(query, 1)
        .expect("a full suffix array with sampling rate 1 is always supported");
    candidates
        .into_iter()
        .filter(|m| is_unique(&query_index, &query[m.query_pos..m.query_pos + m.len]))
        .collect()
}

/// Find matches maximal at both ends and unique in the reference.
///
/// The same matched string may occur more than once in the query. Uniqueness
/// is tested on the completed maximal match, not on its shorter seed.
/// A zero minimum length is invalid and returns an empty result.
pub fn find_mams(reference: &SparseSuffixArray, query: &[u8], min_len: usize) -> Vec<Match> {
    find_mems(reference, query, min_len)
        .into_iter()
        .filter(|m| is_unique(reference, &query[m.query_pos..m.query_pos + m.len]))
        .collect()
}

fn is_unique(index: &SparseSuffixArray, pattern: &[u8]) -> bool {
    matches!(index.search(pattern), Some((start, end)) if start == end)
}

/// Find all exact matches maximal for their particular pair of positions.
///
/// At each query position, look up an exact minimum-length seed, reject pairs
/// that can extend left, and extend every remaining pair to its right boundary.
/// This preserves overlapping matches on different diagonals: containment in
/// the reference and query intervals alone does not make a match redundant.
/// A zero minimum length is invalid and returns an empty result.
pub fn find_mems(reference: &SparseSuffixArray, query: &[u8], min_len: usize) -> Vec<Match> {
    let reference_seq = reference.sequence();
    if min_len == 0 || min_len > query.len() || min_len > reference_seq.len() {
        return Vec::new();
    }

    let mut matches = Vec::new();
    for query_pos in 0..=query.len() - min_len {
        let seed = &query[query_pos..query_pos + min_len];
        if seed.contains(&0) {
            continue;
        }
        let Some((start, end)) = reference.search(seed) else {
            continue;
        };
        for &ref_pos in &reference.suffix_array()[start..=end] {
            if ref_pos > 0
                && query_pos > 0
                && reference_seq[ref_pos - 1] != 0
                && reference_seq[ref_pos - 1] == query[query_pos - 1]
            {
                continue;
            }
            let mut len = min_len;
            while ref_pos + len < reference_seq.len()
                && query_pos + len < query.len()
                && reference_seq[ref_pos + len] != 0
                && reference_seq[ref_pos + len] == query[query_pos + len]
            {
                len += 1;
            }
            matches.push(Match::new(ref_pos, query_pos, len));
        }
    }
    matches.sort_by_key(|m| (m.ref_pos, m.query_pos, m.len));
    matches
}

pub fn run_mummer_algorithm(
    reference: &SparseSuffixArray,
    query: &[u8],
    algorithm: MatchType,
    min_len: usize,
) -> Vec<Match> {
    match algorithm {
        MatchType::MUM => find_mums(reference, query, min_len),
        MatchType::MAM => find_mams(reference, query, min_len),
        MatchType::MEM => find_mems(reference, query, min_len),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn index(sequence: &[u8]) -> SparseSuffixArray {
        SparseSuffixArray::new(sequence, 1).unwrap()
    }

    #[test]
    fn query_repeat_is_mam_but_not_mum() {
        let reference = index(b"ACG");
        let expected = vec![Match::new(0, 0, 3), Match::new(0, 3, 3)];
        assert_eq!(find_mems(&reference, b"ACGACG", 3), expected);
        assert_eq!(find_mams(&reference, b"ACGACG", 3), expected);
        assert!(find_mums(&reference, b"ACGACG", 3).is_empty());
    }

    #[test]
    fn overlapping_query_occurrences_are_not_unique() {
        let reference = index(b"AAA");
        assert_eq!(
            find_mams(&reference, b"AAAA", 3),
            vec![Match::new(0, 0, 3), Match::new(0, 1, 3)]
        );
        assert!(find_mums(&reference, b"AAAA", 3).is_empty());
    }

    #[test]
    fn overlapping_reference_occurrences_are_not_unique() {
        let reference = index(b"AAAA");
        assert_eq!(
            find_mems(&reference, b"AAA", 3),
            vec![Match::new(0, 0, 3), Match::new(1, 0, 3)]
        );
        assert!(find_mams(&reference, b"AAA", 3).is_empty());
        assert!(find_mums(&reference, b"AAA", 3).is_empty());
    }

    #[test]
    fn both_sides_unique_only_the_maximal_match_is_reported() {
        let reference = index(b"ACGT");
        let expected = vec![Match::new(0, 0, 4)];
        assert_eq!(find_mems(&reference, b"ACGT", 2), expected);
        assert_eq!(find_mams(&reference, b"ACGT", 2), expected);
        assert_eq!(find_mums(&reference, b"ACGT", 2), expected);
    }

    #[test]
    fn uniqueness_is_checked_after_extension() {
        let reference = index(b"ACACG");
        // The seed AC is repeated, but the maximal string ACG is unique.
        assert_eq!(find_mums(&reference, b"ACG", 2), vec![Match::new(2, 0, 3)]);
    }

    #[test]
    fn different_diagonals_are_not_removed_by_interval_containment() {
        let reference = index(b"ATAT");
        assert_eq!(
            find_mems(&reference, b"ATAT", 2),
            vec![
                Match::new(0, 0, 4),
                Match::new(0, 2, 2),
                Match::new(2, 0, 2),
            ]
        );
        assert_eq!(find_mums(&reference, b"ATAT", 2), vec![Match::new(0, 0, 4)]);
    }

    #[test]
    fn reference_records_are_boundaries_and_share_uniqueness_scope() {
        let reference = index(b"AC\0GT");
        assert!(find_mems(&reference, b"ACGT", 4).is_empty());
        assert_eq!(
            find_mems(&reference, b"ACGT", 2),
            vec![Match::new(0, 0, 2), Match::new(3, 2, 2)]
        );
        let repeated_records = index(b"ACG\0ACG");
        assert_eq!(
            find_mems(&repeated_records, b"ACG", 3),
            vec![Match::new(0, 0, 3), Match::new(4, 0, 3)]
        );
        assert!(find_mams(&repeated_records, b"ACG", 3).is_empty());
        assert!(find_mums(&repeated_records, b"ACG", 3).is_empty());
    }

    #[test]
    fn equal_barriers_do_not_extend_a_match() {
        let reference = index(b"AC\0GT");
        assert_eq!(
            find_mems(&reference, b"AC\0GT", 2),
            vec![Match::new(0, 0, 2), Match::new(3, 3, 2)]
        );
        assert!(find_mems(&reference, b"AC\0GT", 3).is_empty());
    }

    #[test]
    fn empty_short_and_absent_inputs_return_no_matches() {
        for kind in [MatchType::MEM, MatchType::MAM, MatchType::MUM] {
            assert!(run_mummer_algorithm(&index(b""), b"A", kind, 1).is_empty());
            assert!(run_mummer_algorithm(&index(b"A"), b"", kind, 1).is_empty());
            assert!(run_mummer_algorithm(&index(b"A"), b"C", kind, 1).is_empty());
            assert!(run_mummer_algorithm(&index(b"A"), b"AA", kind, 2).is_empty());
            assert!(run_mummer_algorithm(&index(b"AC"), b"AC", kind, 0).is_empty());
            assert!(run_mummer_algorithm(&index(b"AC"), b"AC", kind, usize::MAX).is_empty());
        }
    }

    // This oracle enumerates every pair of starts and every possible length.
    // It uses neither suffix arrays, seed extension, nor production helpers.
    fn oracle(reference: &[u8], query: &[u8], min_len: usize, kind: MatchType) -> Vec<Match> {
        let mut result = Vec::new();
        for r in 0..reference.len() {
            for q in 0..query.len() {
                for len in min_len..=(reference.len() - r).min(query.len() - q) {
                    let pattern = &reference[r..r + len];
                    if pattern.contains(&0) || pattern != &query[q..q + len] {
                        continue;
                    }
                    let left_maximal = r == 0
                        || q == 0
                        || reference[r - 1] == 0
                        || reference[r - 1] != query[q - 1];
                    let right_maximal = r + len == reference.len()
                        || q + len == query.len()
                        || reference[r + len] == 0
                        || reference[r + len] != query[q + len];
                    if !left_maximal || !right_maximal {
                        continue;
                    }
                    if kind != MatchType::MEM
                        && reference.windows(len).filter(|s| *s == pattern).count() != 1
                    {
                        continue;
                    }
                    if kind == MatchType::MUM
                        && query.windows(len).filter(|s| *s == pattern).count() != 1
                    {
                        continue;
                    }
                    result.push(Match::new(r, q, len));
                }
            }
        }
        result
    }

    fn words(alphabet: &[u8], max_len: usize) -> Vec<Vec<u8>> {
        let mut all = vec![Vec::new()];
        let mut level = vec![Vec::new()];
        for _ in 0..max_len {
            level = level
                .iter()
                .flat_map(|prefix| {
                    alphabet.iter().map(move |&base| {
                        let mut next = prefix.clone();
                        next.push(base);
                        next
                    })
                })
                .collect();
            all.extend(level.iter().cloned());
        }
        all
    }

    fn check_exhaustive(alphabet: &[u8], max_len: usize) {
        let sequences = words(alphabet, max_len);
        for reference in &sequences {
            let sa = index(reference);
            for query in &sequences {
                for min_len in 1..=max_len + 1 {
                    for kind in [MatchType::MEM, MatchType::MAM, MatchType::MUM] {
                        assert_eq!(
                            run_mummer_algorithm(&sa, query, kind, min_len),
                            oracle(reference, query, min_len, kind),
                            "reference={reference:?}, query={query:?}, min_len={min_len}, kind={kind:?}"
                        );
                    }
                }
            }
        }
    }

    #[test]
    fn exhaustive_binary_sequences_match_independent_oracle() {
        // 63 sequences, 3,969 pairs, 6 thresholds, 3 match modes.
        check_exhaustive(b"AC", 5);
    }

    #[test]
    fn exhaustive_sequences_with_boundaries_match_independent_oracle() {
        // 40 sequences, 1,600 pairs, 4 thresholds, 3 match modes.
        check_exhaustive(b"AC\0", 3);
    }
}
