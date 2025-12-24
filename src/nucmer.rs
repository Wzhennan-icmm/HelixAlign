use rayon::prelude::*;
use crate::{SparseSuffixArray, run_mummer_algorithm, MatchType, Match, DnaSequence};
use indicatif::{ProgressBar, ProgressStyle};

#[derive(Debug, Clone)]
pub struct NucmerOptions {
    pub match_type: MatchType,
    pub min_len: usize,
    pub forward_only: bool,
    pub reverse_only: bool,
    // Additional nucmer parameters
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
            match_type: MatchType::MAM,  // Default to MAM (MUMREFERENCE equivalent)
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

pub struct NucmerAligner {
    reference_sa: SparseSuffixArray,
    options: NucmerOptions,
}

impl NucmerAligner {
    pub fn new(reference: &[u8], options: NucmerOptions) -> Result<Self, String> {
        let reference_sa = SparseSuffixArray::new(reference, 1)?;
        
        Ok(Self {
            reference_sa,
            options,
        })
    }

    pub fn align(&self, query: &[u8]) -> Vec<Match> {
        let mut all_matches = Vec::new();

        // Forward alignment
        if !self.options.reverse_only {
            let forward_matches = run_mummer_algorithm(
                &self.reference_sa,
                query,
                self.options.match_type.clone(),
                self.options.min_len
            );
            all_matches.extend(forward_matches);
        }

        // Reverse complement alignment
        if !self.options.forward_only {
            // Create a sequence object to use the reverse_complement method
            let query_seq = DnaSequence::new(std::str::from_utf8(query).unwrap_or(""), "query".to_string());
            let rev_query_seq = query_seq.reverse_complement();
            let rev_query = rev_query_seq.sequence;
            
            let reverse_matches = run_mummer_algorithm(
                &self.reference_sa,
                &rev_query,
                self.options.match_type.clone(),
                self.options.min_len
            );
            
            // Adjust reverse matches to original query coordinates
            let adjusted_reverse_matches: Vec<Match> = reverse_matches
                .into_iter()
                .map(|mut m| {
                    // Convert reverse query position back to original query position
                    m.query_pos = query.len() - m.query_pos - m.len;
                    m
                })
                .collect();
                
            all_matches.extend(adjusted_reverse_matches);
        }

        all_matches
    }

    // Parallel version of align that processes multiple query sequences in parallel with progress bar
    pub fn align_parallel(&self, queries: &[Vec<u8>], num_threads: Option<usize>) -> Vec<Vec<Match>> {
        if let Some(threads) = num_threads {
            rayon::ThreadPoolBuilder::new()
                .num_threads(threads)
                .build_global()
                .ok(); // Ignore errors if global pool is already initialized
        }

        // Create progress bar
        let pb = ProgressBar::new(queries.len() as u64);
        pb.set_style(ProgressStyle::default_bar()
            .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
            .unwrap()
            .progress_chars("#>-"));

        let results: Vec<Vec<Match>> = queries
            .par_iter()
            .map(|query| {
                let result = self.align(query);
                pb.inc(1);
                result
            })
            .collect();

        pb.finish_with_message("Alignment completed");
        results
    }
}

// Function to align multiple query sequences in parallel with progress bar
pub fn align_multiple_sequences_parallel(
    reference: &[u8],
    queries: &[Vec<u8>],
    options: NucmerOptions,
    num_threads: Option<usize>,
) -> Result<Vec<Vec<Match>>, String> {
    if let Some(threads) = num_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok(); // Ignore errors if global pool is already initialized
    }

    let aligner = NucmerAligner::new(reference, options)?;
    
    // Create progress bar
    let pb = ProgressBar::new(queries.len() as u64);
    pb.set_style(ProgressStyle::default_bar()
        .template("{spinner:.green} [{elapsed_precise}] [{bar:40.cyan/blue}] {pos}/{len} ({eta})")
        .unwrap()
        .progress_chars("#>-"));

    let results: Vec<Vec<Match>> = queries
        .par_iter()
        .map(|query| {
            let result = aligner.align(query);
            pb.inc(1);
            result
        })
        .collect();

    pb.finish_with_message("Alignment completed");
    
    Ok(results)
}
