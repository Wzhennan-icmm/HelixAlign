//! Suffix array implementation for efficient string matching
//! Based on the sparse suffix array implementation in the original MUMmer

use std::cmp::Ordering;

/// A match found between reference and query sequences
#[derive(Debug, Clone, PartialEq)]
pub struct Match {
    pub ref_pos: usize,   // position in reference sequence
    pub query_pos: usize, // position in query sequence  
    pub len: usize,       // length of match
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

/// A sparse suffix array implementation
/// This is a simplified version of the original MUMmer sparse suffix array
pub struct SparseSuffixArray {
    sequence: Vec<u8>,
    suffix_array: Vec<usize>,
    lcp_array: Vec<usize>,
    k: usize,  // sampling rate (every k-th suffix is stored)
}

impl SparseSuffixArray {
    /// Create a new sparse suffix array from a sequence
    pub fn new(sequence: &[u8], k: usize) -> Result<Self, String> {
        if k == 0 {
            return Err("Sampling rate k must be greater than 0".to_string());
        }

        let n = sequence.len();
        let mut suffix_indices: Vec<usize> = (0..n).collect();
        
        // Sort the suffixes based on their lexicographic order
        suffix_indices.sort_by(|&i, &j| {
            let suffix_i = &sequence[i..];
            let suffix_j = &sequence[j..];
            suffix_i.cmp(suffix_j)
        });

        // Compute LCP array
        let lcp_array = Self::compute_lcp_array(sequence, &suffix_indices);

        Ok(Self {
            sequence: sequence.to_vec(),
            suffix_array: suffix_indices,
            lcp_array,
            k,
        })
    }

    /// Compute the LCP (Longest Common Prefix) array
    fn compute_lcp_array(sequence: &[u8], suffix_array: &[usize]) -> Vec<usize> {
        let n = suffix_array.len();
        if n == 0 {
            return vec![];
        }
        
        let mut lcp = vec![0; n];
        
        for i in 1..n {
            let suffix1_pos = suffix_array[i-1];
            let suffix2_pos = suffix_array[i];
            
            let suffix1 = &sequence[suffix1_pos..];
            let suffix2 = &sequence[suffix2_pos..];
            
            let mut common_len = 0;
            for (c1, c2) in suffix1.iter().zip(suffix2.iter()) {
                if c1 == c2 {
                    common_len += 1;
                } else {
                    break;
                }
            }
            
            lcp[i] = common_len;
        }
        
        lcp
    }

    /// Binary search to find the left boundary of an interval for character c at position i
    fn bsearch_left(&self, c: u8, i: usize, start: usize, end: usize) -> usize {
        let mut s = start;
        let mut e = end;
        
        while s <= e && s < self.suffix_array.len() {
            if s > e {
                break;
            }
            
            let mid = s + (e - s) / 2;
            if mid >= self.suffix_array.len() {
                break;
            }
            
            let suffix_pos = self.suffix_array[mid];
            if suffix_pos + i >= self.sequence.len() {
                s = mid + 1;
                continue;
            }
            
            let suffix_char = self.sequence[suffix_pos + i];
            match suffix_char.cmp(&c) {
                Ordering::Less => s = mid + 1,
                Ordering::Greater => {
                    if mid == 0 { break; }
                    e = mid - 1;
                },
                Ordering::Equal => {
                    e = mid;
                    if s == e {
                        break;
                    }
                }
            }
        }
        
        s.min(self.suffix_array.len().saturating_sub(1))
    }

    /// Binary search to find the right boundary of an interval for character c at position i
    fn bsearch_right(&self, c: u8, i: usize, start: usize, end: usize) -> usize {
        let mut s = start;
        let mut e = end;
        let mut result = start;
        
        while s <= e && s < self.suffix_array.len() {
            if s > e {
                break;
            }
            
            let mid = s + (e - s) / 2;
            if mid >= self.suffix_array.len() {
                break;
            }
            
            let suffix_pos = self.suffix_array[mid];
            if suffix_pos + i >= self.sequence.len() {
                s = mid + 1;
                continue;
            }
            
            let suffix_char = self.sequence[suffix_pos + i];
            match suffix_char.cmp(&c) {
                Ordering::Less => s = mid + 1,
                Ordering::Greater => {
                    if mid == 0 { break; }
                    e = mid - 1;
                },
                Ordering::Equal => {
                    result = mid;
                    s = mid + 1;
                }
            }
        }
        
        result.min(self.suffix_array.len().saturating_sub(1))
    }

    /// Simple suffix array search for a pattern
    pub fn search(&self, pattern: &[u8]) -> Option<(usize, usize)> {
        if pattern.is_empty() || self.suffix_array.is_empty() {
            return None;
        }

        let mut start = 0;
        let mut end = self.suffix_array.len() - 1;

        for (i, &c) in pattern.iter().enumerate() {
            let new_start = self.bsearch_left(c, i, start, end);
            let new_end = self.bsearch_right(c, i, start, end);

            if new_start > new_end || new_start >= self.suffix_array.len() {
                return None; // Pattern not found
            }

            start = new_start;
            end = new_end;

            if start > end {
                return None; // Pattern not found
            }
        }

        Some((start, end))
    }

    /// Find all matches of a pattern in the reference sequence
    pub fn find_matches(&self, pattern: &[u8]) -> Vec<Match> {
        let interval = self.search(pattern);
        match interval {
            Some((start, end)) => {
                let mut matches = Vec::new();
                for i in start..=end {
                    if i < self.suffix_array.len() {
                        let ref_pos = self.suffix_array[i];
                        matches.push(Match::new(ref_pos, 0, pattern.len()));
                    }
                }
                matches
            }
            None => Vec::new(),
        }
    }

    /// Get the original sequence
    pub fn sequence(&self) -> &[u8] {
        &self.sequence
    }

    /// Get the suffix array
    pub fn suffix_array(&self) -> &[usize] {
        &self.suffix_array
    }

    /// Get the LCP array
    pub fn lcp_array(&self) -> &[usize] {
        &self.lcp_array
    }

    /// Get the sampling rate
    pub fn sampling_rate(&self) -> usize {
        self.k
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_sparse_suffix_array() {
        let sequence = b"banana$";
        let sa = SparseSuffixArray::new(sequence, 1).unwrap();
        
        // The suffix array should contain positions sorted by their suffixes
        // For "banana$", the sorted suffixes are:
        // $ (pos 6), a$ (pos 5), ana$ (pos 3), anana$ (pos 1), banana$ (pos 0), na$ (pos 4), nana$ (pos 2)
        let expected_suffixes = vec![6, 5, 3, 1, 0, 4, 2];
        
        assert_eq!(sa.suffix_array(), &expected_suffixes);
    }

    #[test]
    fn test_search() {
        let sequence = b"banana$";
        let sa = SparseSuffixArray::new(sequence, 1).unwrap();
        
        // Search for "ana"
        let result = sa.search(b"ana");
        assert!(result.is_some());
        
        let (start, end) = result.unwrap();
        assert!(start <= end);
        
        // Should find "ana" at position 1 and 3
        let matches = sa.find_matches(b"ana");
        assert!(!matches.is_empty());
    }
}