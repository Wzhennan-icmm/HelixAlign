//! DNA/RNA sequence handling utilities for MUMmer

use std::fmt;

/// A DNA sequence that can be efficiently processed
#[derive(Debug, Clone, PartialEq)]
pub struct DnaSequence {
    pub sequence: Vec<u8>,
    pub description: String,
}

impl DnaSequence {
    /// Create a new DNA sequence from a string
    pub fn new(seq: &str, description: String) -> Self {
        let sequence = seq.as_bytes().to_vec();
        Self { sequence, description }
    }

    /// Get the length of the sequence
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    /// Check if the sequence is empty
    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }

    /// Get a character at a specific position (0-indexed)
    pub fn get(&self, index: usize) -> Option<u8> {
        self.sequence.get(index).copied()
    }

    /// Convert character to nucleotide code (A=0, C=1, G=2, T=3)
    pub fn char_to_code(c: u8) -> Option<u8> {
        match c {
            b'A' | b'a' => Some(0),
            b'C' | b'c' => Some(1),
            b'G' | b'g' => Some(2),
            b'T' | b't' => Some(3),
            _ => None,
        }
    }

    /// Convert nucleotide code back to character
    pub fn code_to_char(code: u8) -> Option<u8> {
        match code {
            0 => Some(b'A'),
            1 => Some(b'C'),
            2 => Some(b'G'),
            3 => Some(b'T'),
            _ => None,
        }
    }

    /// Get a substring as a new DnaSequence
    pub fn substring(&self, start: usize, end: usize) -> Option<Self> {
        if start <= end && end <= self.sequence.len() {
            let sub_seq = self.sequence[start..end].to_vec();
            Some(DnaSequence {
                sequence: sub_seq,
                description: format!("{}[{}..{}]", self.description, start, end),
            })
        } else {
            None
        }
    }

    /// Reverse complement of the DNA sequence
    pub fn reverse_complement(&self) -> Self {
        let mut complement = Vec::with_capacity(self.sequence.len());
        for &base in self.sequence.iter().rev() {
            let comp_base = match base {
                b'A' | b'a' => b'T',
                b'T' | b't' => b'A',
                b'G' | b'g' => b'C',
                b'C' | b'c' => b'G',
                _ => base, // Keep non-standard bases as is
            };
            complement.push(comp_base);
        }
        DnaSequence {
            sequence: complement,
            description: format!("reverse complement of {}", self.description),
        }
    }
}

impl fmt::Display for DnaSequence {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}: {}", self.description, String::from_utf8_lossy(&self.sequence))
    }
}

/// A collection of DNA sequences
#[derive(Debug, Clone)]
pub struct SequenceCollection {
    pub sequences: Vec<DnaSequence>,
}

impl SequenceCollection {
    pub fn new() -> Self {
        Self {
            sequences: Vec::new(),
        }
    }

    pub fn add_sequence(&mut self, sequence: DnaSequence) {
        self.sequences.push(sequence);
    }

    pub fn len(&self) -> usize {
        self.sequences.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequences.is_empty()
    }

    pub fn total_length(&self) -> usize {
        self.sequences.iter().map(|s| s.len()).sum()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_dna_sequence_creation() {
        let seq = DnaSequence::new("ATCG", "test_sequence".to_string());
        assert_eq!(seq.len(), 4);
        assert_eq!(seq.get(0), Some(b'A'));
        assert_eq!(seq.get(1), Some(b'T'));
        assert_eq!(seq.get(2), Some(b'C'));
        assert_eq!(seq.get(3), Some(b'G'));
    }

    #[test]
    fn test_reverse_complement() {
        let seq = DnaSequence::new("ATCG", "test_sequence".to_string());
        let rev_comp = seq.reverse_complement();
        assert_eq!(String::from_utf8_lossy(&rev_comp.sequence), "CGAT");
    }

    #[test]
    fn test_char_to_code() {
        assert_eq!(DnaSequence::char_to_code(b'A'), Some(0));
        assert_eq!(DnaSequence::char_to_code(b'a'), Some(0));
        assert_eq!(DnaSequence::char_to_code(b'C'), Some(1));
        assert_eq!(DnaSequence::char_to_code(b'c'), Some(1));
        assert_eq!(DnaSequence::char_to_code(b'G'), Some(2));
        assert_eq!(DnaSequence::char_to_code(b'g'), Some(2));
        assert_eq!(DnaSequence::char_to_code(b'T'), Some(3));
        assert_eq!(DnaSequence::char_to_code(b't'), Some(3));
    }
}