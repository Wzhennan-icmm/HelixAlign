//! Legacy naming for exact seed matching only.
//!
//! This module does not cluster anchors, perform gapped extension, or implement
//! nucmer. Prefer ExactMatcher for explicit FASTA-record handling. Legacy
//! clustering fields retain their defaults for construction compatibility;
//! changing them is rejected instead of silently pretending to apply them.
use crate::{ExactMatcher, FastaRecord, MatchType, RecordMatch, StrandMode};
use rayon::prelude::*;

#[derive(Debug, Clone)]
pub struct NucmerOptions {
    pub match_type: MatchType,
    pub min_len: usize,
    pub forward_only: bool,
    pub reverse_only: bool,
    // Unsupported legacy nucmer settings. Only their defaults are accepted.
    pub break_len: usize,
    pub min_cluster: usize,
    pub diag_diff: usize,
    pub diag_factor: f64,
    pub max_gap: usize,
    pub extend: bool,
    pub optimize: bool,
    pub simplify: bool,
    pub banding: bool,
    pub use_extent: bool,
    pub to_seqend: bool,
    pub do_delta: bool,
    pub do_shadows: bool,
}

impl Default for NucmerOptions {
    fn default() -> Self {
        Self {
            match_type: MatchType::MAM,
            min_len: 20,
            forward_only: false,
            reverse_only: false,
            break_len: 200,
            min_cluster: 65,
            diag_diff: 5,
            diag_factor: 0.12,
            max_gap: 90,
            extend: true,
            optimize: true,
            simplify: true,
            banding: false,
            use_extent: false,
            to_seqend: false,
            do_delta: true,
            do_shadows: false,
        }
    }
}

/// Single-reference compatibility wrapper for exact matching, not nucmer.
/// Results now include strand; the previous Vec<Match> discarded that fact.
pub struct NucmerAligner {
    matcher: ExactMatcher,
    options: NucmerOptions,
}

impl NucmerAligner {
    pub fn new(reference: &[u8], options: NucmerOptions) -> Result<Self, String> {
        if options.min_len == 0 {
            return Err("Minimum match length must be positive".into());
        }
        if options.forward_only && options.reverse_only {
            return Err("forward_only and reverse_only are mutually exclusive".into());
        }
        let defaults = NucmerOptions::default();
        if options.break_len != defaults.break_len
            || options.min_cluster != defaults.min_cluster
            || options.diag_diff != defaults.diag_diff
            || options.diag_factor != defaults.diag_factor
            || options.max_gap != defaults.max_gap
            || options.extend != defaults.extend
            || options.optimize != defaults.optimize
            || options.simplify != defaults.simplify
            || options.banding != defaults.banding
            || options.use_extent != defaults.use_extent
            || options.to_seqend != defaults.to_seqend
            || options.do_delta != defaults.do_delta
            || options.do_shadows != defaults.do_shadows
        {
            return Err("nucmer clustering, extension and delta options are not implemented; use ExactMatcher for exact seeds".into());
        }
        let matcher = ExactMatcher::new(vec![FastaRecord {
            id: "reference".into(),
            sequence: reference.to_vec(),
        }])?;
        Ok(Self { matcher, options })
    }

    /// Return exact matches with explicit strand and zero-based, forward-record
    /// left coordinates. Query uniqueness is local to this query sequence.
    pub fn align(&self, query: &[u8]) -> Vec<RecordMatch> {
        if query.is_empty() {
            return Vec::new();
        }
        let strands = if self.options.forward_only {
            StrandMode::Forward
        } else if self.options.reverse_only {
            StrandMode::Reverse
        } else {
            StrandMode::Both
        };
        let query = FastaRecord {
            id: "query".into(),
            sequence: query.to_vec(),
        };
        self.matcher
            .find_matches(
                &query,
                self.options.match_type,
                self.options.min_len,
                strands,
            )
            .expect("validated nonempty query and positive minimum length")
    }

    /// Use a local thread pool so repeated calls honor their requested size.
    pub fn align_parallel(
        &self,
        queries: &[Vec<u8>],
        num_threads: Option<usize>,
    ) -> Result<Vec<Vec<RecordMatch>>, String> {
        if num_threads == Some(0) {
            return Err("Thread count must be positive".into());
        }
        let available = std::thread::available_parallelism().map_or(1, usize::from);
        let threads = num_threads.unwrap_or_else(|| queries.len().min(available).max(1));
        let pool = rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build()
            .map_err(|e| e.to_string())?;
        Ok(pool.install(|| queries.par_iter().map(|query| self.align(query)).collect()))
    }
}

pub fn align_multiple_sequences_parallel(
    reference: &[u8],
    queries: &[Vec<u8>],
    options: NucmerOptions,
    num_threads: Option<usize>,
) -> Result<Vec<Vec<RecordMatch>>, String> {
    NucmerAligner::new(reference, options)?.align_parallel(queries, num_threads)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::Strand;
    #[test]
    fn legacy_wrapper_preserves_reverse_orientation() {
        let options = NucmerOptions {
            min_len: 8,
            match_type: MatchType::MUM,
            ..Default::default()
        };
        let aligner = NucmerAligner::new(b"AACCGGTA", options).unwrap();
        assert_eq!(
            aligner.align(b"CCTACCGGTTG"),
            vec![RecordMatch {
                reference_index: 0,
                ref_pos: 0,
                query_pos: 2,
                len: 8,
                strand: Strand::Reverse,
            }]
        );
    }
    #[test]
    fn legacy_wrapper_rejects_contradictory_or_inert_settings() {
        for options in [
            NucmerOptions {
                min_len: 0,
                ..Default::default()
            },
            NucmerOptions {
                forward_only: true,
                reverse_only: true,
                ..Default::default()
            },
            NucmerOptions {
                min_cluster: 99,
                ..Default::default()
            },
        ] {
            assert!(NucmerAligner::new(b"ACGT", options).is_err());
        }
    }
}
