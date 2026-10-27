//! Contig-aware, strand-aware exact seeds (not gapped nucmer alignments).

use crate::{FastaRecord, Match, MatchType, SparseSuffixArray, run_mummer_algorithm};

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum Strand {
    Forward,
    Reverse,
}

impl Strand {
    pub fn symbol(self) -> &'static str {
        match self {
            Self::Forward => "+",
            Self::Reverse => "-",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum StrandMode {
    Forward,
    Reverse,
    Both,
}

/// All positions are zero-based and relative to the original forward records.
/// `query_pos` is the left coordinate even for a reverse-strand match.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RecordMatch {
    pub reference_index: usize,
    pub ref_pos: usize,
    pub query_pos: usize,
    pub len: usize,
    pub strand: Strand,
}

pub(crate) fn matching_dna(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .map(|b| match b.to_ascii_uppercase() {
            base @ (b'A' | b'C' | b'G' | b'T') => base,
            _ => 0,
        })
        .collect()
}

pub(crate) fn reverse_matching_dna(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => 0,
        })
        .collect()
}

/// Reference occurrence counts cover every reference record. Query occurrence
/// counts cover only the current query record and orientation, as in mummer.
pub struct ExactMatcher {
    records: Vec<FastaRecord>,
    starts: Vec<usize>,
    suffix_array: SparseSuffixArray,
}

impl ExactMatcher {
    pub fn new(records: Vec<FastaRecord>) -> Result<Self, String> {
        if records.is_empty() || records.iter().any(|r| r.sequence.is_empty()) {
            return Err("Reference must contain nonempty FASTA records".into());
        }
        let mut ids = std::collections::HashSet::new();
        let mut sequence = Vec::new();
        let mut starts = Vec::new();
        for record in &records {
            if !ids.insert(&record.id) {
                return Err(format!("Duplicate reference ID: {}", record.id));
            }
            starts.push(sequence.len());
            sequence.extend(matching_dna(&record.sequence));
            // NUL never matches, including another NUL. It is not a biological base.
            sequence.push(0);
        }
        let suffix_array = SparseSuffixArray::new(&sequence, 1)?;
        Ok(Self {
            records,
            starts,
            suffix_array,
        })
    }

    pub fn references(&self) -> &[FastaRecord] {
        &self.records
    }

    pub fn find_matches(
        &self,
        query: &FastaRecord,
        algorithm: MatchType,
        min_len: usize,
        strands: StrandMode,
    ) -> Result<Vec<RecordMatch>, String> {
        if min_len == 0 {
            return Err("Minimum match length must be greater than zero".into());
        }
        if query.sequence.is_empty() {
            return Err("Query record is empty".into());
        }
        let forward = matching_dna(&query.sequence);
        let mut matches = Vec::new();
        if strands != StrandMode::Reverse {
            for m in run_mummer_algorithm(&self.suffix_array, &forward, algorithm, min_len) {
                matches.push(self.local_match(m, Strand::Forward, forward.len())?);
            }
        }
        if strands != StrandMode::Forward {
            let reverse = reverse_matching_dna(&forward);
            for m in run_mummer_algorithm(&self.suffix_array, &reverse, algorithm, min_len) {
                matches.push(self.local_match(m, Strand::Reverse, forward.len())?);
            }
        }
        matches.sort_by_key(|m| (m.reference_index, m.ref_pos, m.query_pos, m.len, m.strand));
        Ok(matches)
    }

    fn local_match(
        &self,
        m: Match,
        strand: Strand,
        query_len: usize,
    ) -> Result<RecordMatch, String> {
        let index = self
            .starts
            .partition_point(|&start| start <= m.ref_pos)
            .checked_sub(1)
            .ok_or("Internal match before first reference record")?;
        let ref_pos = m.ref_pos - self.starts[index];
        if ref_pos + m.len > self.records[index].sequence.len() || m.query_pos + m.len > query_len {
            return Err("Internal match crosses a FASTA record boundary".into());
        }
        let query_pos = if strand == Strand::Forward {
            m.query_pos
        } else {
            query_len - m.query_pos - m.len
        };
        Ok(RecordMatch {
            reference_index: index,
            ref_pos,
            query_pos,
            len: m.len,
            strand,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse_fasta_records;

    fn check(
        reference: &str,
        query: &str,
        mode: MatchType,
        min_len: usize,
        strands: StrandMode,
    ) -> Vec<Vec<RecordMatch>> {
        let matcher = ExactMatcher::new(parse_fasta_records(reference).unwrap()).unwrap();
        parse_fasta_records(query)
            .unwrap()
            .iter()
            .map(|q| matcher.find_matches(q, mode, min_len, strands).unwrap())
            .collect()
    }

    #[test]
    fn no_reference_or_query_junction_matches() {
        for mode in [MatchType::MEM, MatchType::MAM, MatchType::MUM] {
            assert!(
                check(
                    ">r1\nAAAACCCC\n>r2\nGGGGTTTT\n",
                    ">q\nCCCCGGGG\n",
                    mode,
                    8,
                    StrandMode::Both
                )[0]
                .is_empty()
            );
            assert!(
                check(
                    ">r\nCCCCGGGG\n",
                    ">q1\nAAAACCCC\n>q2\nGGGGTTTT\n",
                    mode,
                    8,
                    StrandMode::Both
                )
                .iter()
                .all(Vec::is_empty)
            );
        }
    }

    #[test]
    fn reference_uniqueness_is_global_query_uniqueness_is_record_local() {
        assert!(
            check(
                ">r1\nACGTTGCA\n>r2\nACGTTGCA\n",
                ">q\nACGTTGCA\n",
                MatchType::MUM,
                8,
                StrandMode::Forward
            )[0]
            .is_empty()
        );
        assert!(
            check(
                ">r1\nACGTTGCA\n>r2\nACGTTGCA\n",
                ">q\nACGTTGCA\n",
                MatchType::MAM,
                8,
                StrandMode::Forward
            )[0]
            .is_empty()
        );
        let matches = check(
            ">r\nACGTTGCA\n",
            ">q1\nACGTTGCA\n>q2\nACGTTGCA\n",
            MatchType::MUM,
            8,
            StrandMode::Forward,
        );
        assert_eq!(matches.iter().map(Vec::len).collect::<Vec<_>>(), vec![1, 1]);
    }

    #[test]
    fn reverse_coordinates_and_palindrome_keep_orientation() {
        let matches = check(
            ">r\nAACCGGTA\n",
            ">q\nCCTACCGGTTG\n",
            MatchType::MUM,
            8,
            StrandMode::Both,
        );
        assert_eq!(
            matches[0],
            vec![RecordMatch {
                reference_index: 0,
                ref_pos: 0,
                query_pos: 2,
                len: 8,
                strand: Strand::Reverse
            }]
        );
        let palindrome = check(
            ">r\nACGT\n",
            ">q\nACGT\n",
            MatchType::MUM,
            4,
            StrandMode::Both,
        );
        assert_eq!(
            palindrome[0].iter().map(|m| m.strand).collect::<Vec<_>>(),
            vec![Strand::Forward, Strand::Reverse]
        );
    }

    #[test]
    fn ambiguous_bases_are_barriers_and_contig_coordinates_are_local() {
        let matches = check(
            ">r1\nAAAA\n>r2\nacgtNNNNtgca\n",
            ">q\nACGTNNNNTGCA\n",
            MatchType::MUM,
            4,
            StrandMode::Forward,
        );
        assert_eq!(
            matches[0],
            vec![
                RecordMatch {
                    reference_index: 1,
                    ref_pos: 0,
                    query_pos: 0,
                    len: 4,
                    strand: Strand::Forward
                },
                RecordMatch {
                    reference_index: 1,
                    ref_pos: 8,
                    query_pos: 8,
                    len: 4,
                    strand: Strand::Forward
                },
            ]
        );
    }
}
