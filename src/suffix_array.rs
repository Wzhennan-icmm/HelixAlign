//! Full suffix array used by the exact-match reference implementation.
//!
//! The historical public name is retained for compatibility. Sparse sampling
//! is not implemented, so construction explicitly rejects sampling rates other
//! than one. Suffix sorting and LCP construction are deliberately simple and
//! are not yet suitable as a performance baseline for large genomes.

use std::cmp::Ordering;

/// A match found between one reference and one oriented query sequence.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Match {
    pub ref_pos: usize,
    pub query_pos: usize,
    pub len: usize,
}

impl Match {
    pub fn new(ref_pos: usize, query_pos: usize, len: usize) -> Self {
        Self {
            ref_pos,
            query_pos,
            len,
        }
    }
}

/// A full suffix array (the legacy type name is retained).
pub struct SparseSuffixArray {
    sequence: Vec<u8>,
    suffix_array: Vec<usize>,
    lcp_array: Vec<usize>,
    k: usize,
}

impl SparseSuffixArray {
    /// Construct a full suffix array; only sampling rate `k = 1` is supported.
    /// The sequence need not end in a sentinel. Empty sequences are supported.
    pub fn new(sequence: &[u8], k: usize) -> Result<Self, String> {
        if k != 1 {
            return Err(
                "Only sampling rate k = 1 is supported; sparse sampling is not implemented"
                    .to_string(),
            );
        }
        let mut suffix_indices: Vec<usize> = (0..sequence.len()).collect();
        suffix_indices.sort_by(|&i, &j| sequence[i..].cmp(&sequence[j..]));
        let lcp_array = Self::compute_lcp_array(sequence, &suffix_indices);
        Ok(Self {
            sequence: sequence.to_vec(),
            suffix_array: suffix_indices,
            lcp_array,
            k,
        })
    }

    fn compute_lcp_array(sequence: &[u8], suffix_array: &[usize]) -> Vec<usize> {
        let mut lcp = vec![0; suffix_array.len()];
        for i in 1..suffix_array.len() {
            lcp[i] = sequence[suffix_array[i - 1]..]
                .iter()
                .zip(sequence[suffix_array[i]..].iter())
                .take_while(|(a, b)| a == b)
                .count();
        }
        lcp
    }

    // Equal means that the entire pattern is a prefix of this suffix. A suffix
    // shorter than the pattern is ordered before it, never reported as a hit.
    fn compare_prefix(&self, suffix_pos: usize, pattern: &[u8]) -> Ordering {
        let suffix = &self.sequence[suffix_pos..];
        suffix[..suffix.len().min(pattern.len())].cmp(pattern)
    }

    /// Return the inclusive suffix-array interval for all exact occurrences.
    /// Empty patterns have no interval. The two half-open partition searches
    /// allow empty results at either end, unlike clamping to a valid suffix.
    pub fn search(&self, pattern: &[u8]) -> Option<(usize, usize)> {
        if pattern.is_empty() {
            return None;
        }
        let start = self
            .suffix_array
            .partition_point(|&pos| self.compare_prefix(pos, pattern) == Ordering::Less);
        let end = self
            .suffix_array
            .partition_point(|&pos| self.compare_prefix(pos, pattern) != Ordering::Greater);
        if start < end {
            Some((start, end - 1))
        } else {
            None
        }
    }

    /// Return every exact occurrence, including overlapping occurrences.
    pub fn find_matches(&self, pattern: &[u8]) -> Vec<Match> {
        match self.search(pattern) {
            Some((start, end)) => self.suffix_array[start..=end]
                .iter()
                .map(|&ref_pos| Match::new(ref_pos, 0, pattern.len()))
                .collect(),
            None => Vec::new(),
        }
    }

    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    pub fn suffix_array(&self) -> &[usize] {
        &self.suffix_array
    }

    pub fn lcp_array(&self) -> &[usize] {
        &self.lcp_array
    }

    pub fn sampling_rate(&self) -> usize {
        self.k
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn positions(sa: &SparseSuffixArray, pattern: &[u8]) -> Vec<usize> {
        let mut result: Vec<_> = sa.find_matches(pattern).iter().map(|m| m.ref_pos).collect();
        result.sort_unstable();
        result
    }

    #[test]
    fn suffix_array_and_lcp_are_exact() {
        let sa = SparseSuffixArray::new(b"banana$", 1).unwrap();
        assert_eq!(sa.suffix_array(), &[6, 5, 3, 1, 0, 4, 2]);
        assert_eq!(sa.lcp_array(), &[0, 0, 1, 3, 0, 0, 2]);
        assert_eq!(positions(&sa, b"ana"), vec![1, 3]);
    }

    #[test]
    fn nonexistent_patterns_do_not_produce_a_clamped_hit() {
        let sa = SparseSuffixArray::new(b"A", 1).unwrap();
        for pattern in [b"C".as_slice(), b"AA", b"\0", b"Z", b""] {
            assert_eq!(sa.search(pattern), None, "pattern={pattern:?}");
            assert!(sa.find_matches(pattern).is_empty());
        }
        assert_eq!(positions(&sa, b"A"), vec![0]);
    }

    #[test]
    fn suffix_end_and_overlap_counts_are_exact_without_a_sentinel() {
        let sa = SparseSuffixArray::new(b"AAAA", 1).unwrap();
        assert_eq!(positions(&sa, b"A"), vec![0, 1, 2, 3]);
        assert_eq!(positions(&sa, b"AA"), vec![0, 1, 2]);
        assert_eq!(positions(&sa, b"AAA"), vec![0, 1]);
        assert_eq!(positions(&sa, b"AAAA"), vec![0]);
        assert!(positions(&sa, b"AAAAA").is_empty());
        assert!(positions(&sa, b"AAAC").is_empty());
    }

    #[test]
    fn empty_index_and_unsupported_sampling_are_explicit() {
        let sa = SparseSuffixArray::new(b"", 1).unwrap();
        assert!(sa.suffix_array().is_empty());
        assert!(sa.lcp_array().is_empty());
        assert_eq!(sa.search(b"A"), None);
        assert_eq!(sa.search(b""), None);
        assert_eq!(sa.sampling_rate(), 1);
        assert!(SparseSuffixArray::new(b"AC", 0).is_err());
        assert!(SparseSuffixArray::new(b"AC", 2).is_err());
    }

    fn binary_words(max_len: usize) -> Vec<Vec<u8>> {
        let mut result = Vec::new();
        for len in 0..=max_len {
            for value in 0..(1usize << len) {
                result.push(
                    (0..len)
                        .map(|bit| if value & (1 << bit) == 0 { b'A' } else { b'C' })
                        .collect(),
                );
            }
        }
        result
    }

    #[test]
    fn exhaustive_lookup_matches_direct_substring_counting() {
        let references = binary_words(6);
        let mut patterns = binary_words(5);
        patterns.extend([
            b"G".to_vec(),
            b"ACG".to_vec(),
            b"AAAAAAA".to_vec(),
            vec![0],
            vec![255],
        ]);
        for reference in references {
            let sa = SparseSuffixArray::new(&reference, 1).unwrap();
            for pattern in &patterns {
                let expected: Vec<_> = if pattern.is_empty() {
                    Vec::new()
                } else {
                    reference
                        .windows(pattern.len())
                        .enumerate()
                        .filter_map(|(position, candidate)| {
                            (candidate == pattern).then_some(position)
                        })
                        .collect()
                };
                assert_eq!(
                    positions(&sa, pattern),
                    expected,
                    "reference={reference:?}, pattern={pattern:?}"
                );
            }
        }
    }
}
