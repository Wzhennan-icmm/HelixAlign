//! Rust implementation of HelixAlign - a versatile alignment tool for DNA and protein sequences.
//! This is a command-line tool for finding maximal matches between sequences.

use std::env;
use std::fs;

use helixalign::{SparseSuffixArray, run_mummer_algorithm, MatchType, NucmerAligner, NucmerOptions, Match, parse_fasta, GenomicStats, align_multiple_sequences_parallel, OutputFormat, print_matches_in_format};

fn main() {
    let args: Vec<String> = env::args().collect();
    let program_name = &args[0];
    
    // Check if running as nucmer by program name
    if program_name.contains("nucmer") {
        run_nucmer(args);
        return;
    }
    
    // Otherwise run standard mummer functionality
    if args.len() < 3 {
        print_usage(&args[0]);
        return;
    }
    
    // Parse command line arguments
    let mut min_len = 20;
    let mut algorithm = MatchType::MAM; // Default to MAM (Maximal Almost-Unique Match)
    let mut reference_file = "";
    let mut query_files = Vec::new();
    let mut show_stats = false;
    let mut num_threads: Option<usize> = None;
    let mut output_format = OutputFormat::Default;
    
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-mum" => algorithm = MatchType::MUM,
            "-mumreference" | "-mumcand" => algorithm = MatchType::MAM,  // Same as MAM
            "-maxmatch" => algorithm = MatchType::MEM,
            "-l" => {
                if i + 1 < args.len() {
                    min_len = args[i + 1].parse().expect("Invalid minimum length");
                    i += 1;
                } else {
                    eprintln!("Error: -l requires a value");
                    return;
                }
            }
            "-t" | "--threads" => {
                if i + 1 < args.len() {
                    num_threads = Some(args[i + 1].parse().expect("Invalid thread count"));
                    i += 1;
                } else {
                    eprintln!("Error: -t requires a value");
                    return;
                }
            }
            "-f" | "--format" => {
                if i + 1 < args.len() {
                    output_format = OutputFormat::from_str(&args[i + 1]).unwrap_or(OutputFormat::Default);
                    i += 1;
                } else {
                    eprintln!("Error: -f requires a format (delta, paf, sam)");
                    return;
                }
            }
            "-stats" | "--stats" => {
                show_stats = true;
            }
            arg if !arg.starts_with('-') => {
                if reference_file.is_empty() {
                    reference_file = arg;
                } else {
                    query_files.push(arg.to_string());
                }
            }
            _ => {
                eprintln!("Unknown option: {}", args[i]);
                print_usage(&args[0]);
                return;
            }
        }
        i += 1;
    }
    
    if query_files.is_empty() {
        // If no query files provided, treat the second argument as the only query file
        if args.len() >= 3 && !args[2].starts_with('-') {
            query_files.push(args[2].clone());
        } else {
            eprintln!("Error: No query file provided");
            print_usage(&args[0]);
            return;
        }
    }
    
    // Calculate and print statistics if requested
    if show_stats {
        let ref_sequences = parse_fasta(reference_file);
        let ref_stats = GenomicStats::new(&ref_sequences);
        ref_stats.print_stats("Reference");
        
        for query_file in &query_files {
            let query_sequences = parse_fasta(query_file);
            let query_stats = GenomicStats::new(&query_sequences);
            query_stats.print_stats("Query");
        }
    }
    
    // Set number of threads if specified
    if let Some(threads) = num_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok(); // Ignore errors if global pool is already initialized
    }
    
    // Process each query file
    for query_file in query_files {
        // Read reference and query sequences
        let reference_seq = read_fasta_file(reference_file);
        let query_seq = read_fasta_file(&query_file);
        
        // Create suffix array for reference
        let reference_sa = SparseSuffixArray::new(&reference_seq, 1)
            .expect("Could not create suffix array");
        
        // Find matches - clone algorithm to avoid move error
        let matches = run_mummer_algorithm(&reference_sa, &query_seq, algorithm.clone(), min_len);
        
        // Print matches in the specified format
        print_matches_in_format(&matches, &query_file, &output_format, &reference_seq, &query_seq);
    }
}

fn run_nucmer(args: Vec<String>) {
    if args.len() < 3 {
        print_nucmer_usage(&args[0]);
        return;
    }
    
    // Parse command line arguments for nucmer
    let mut min_len = 20;
    let mut algorithm = MatchType::MAM;
    let mut forward_only = false;
    let mut reverse_only = false;
    let mut reference_file = "";
    let mut query_files = Vec::new();
    let mut show_stats = false;
    let mut num_threads: Option<usize> = None;
    let mut output_format = OutputFormat::Default;
    let mut break_len = 200;
    let mut min_cluster = 65;
    let mut diag_diff = 5;
    let mut diag_factor = 0.12;
    let mut max_gap = 90;
    let mut extend = true;
    let mut optimize = true;
    let mut simplify = true;
    let mut banding = false;
    
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "-mum" => algorithm = MatchType::MUM,
            "-mumreference" | "-mumcand" => algorithm = MatchType::MAM,
            "-maxmatch" => algorithm = MatchType::MEM,
            "-b" | "--breaklen" => {
                if i + 1 < args.len() {
                    break_len = args[i + 1].parse().expect("Invalid break length");
                    i += 1;
                } else {
                    eprintln!("Error: -b requires a value");
                    return;
                }
            }
            "-c" | "--mincluster" => {
                if i + 1 < args.len() {
                    min_cluster = args[i + 1].parse().expect("Invalid minimum cluster length");
                    i += 1;
                } else {
                    eprintln!("Error: -c requires a value");
                    return;
                }
            }
            "-D" | "--diagdiff" => {
                if i + 1 < args.len() {
                    diag_diff = args[i + 1].parse().expect("Invalid diagonal difference");
                    i += 1;
                } else {
                    eprintln!("Error: -D requires a value");
                    return;
                }
            }
            "-d" | "--diagfactor" => {
                if i + 1 < args.len() {
                    diag_factor = args[i + 1].parse().expect("Invalid diagonal factor");
                    i += 1;
                } else {
                    eprintln!("Error: -d requires a value");
                    return;
                }
            }
            "-noextend" => extend = false,
            "-f" | "--forward" => forward_only = true,
            "-g" | "--maxgap" => {
                if i + 1 < args.len() {
                    max_gap = args[i + 1].parse().expect("Invalid max gap");
                    i += 1;
                } else {
                    eprintln!("Error: -g requires a value");
                    return;
                }
            }
            "-l" | "--minmatch" => {
                if i + 1 < args.len() {
                    min_len = args[i + 1].parse().expect("Invalid minimum match length");
                    i += 1;
                } else {
                    eprintln!("Error: -l requires a value");
                    return;
                }
            }
            "-L" | "--minalign" => {
                // For now, we'll just parse this but not implement the functionality
                if i + 1 < args.len() {
                    let _min_align = args[i + 1].parse::<usize>().expect("Invalid minimum alignment length");
                    i += 1;
                } else {
                    eprintln!("Error: -L requires a value");
                    return;
                }
            }
            "-nooptimize" => optimize = false,
            "-r" | "--reverse" => reverse_only = true,
            "-nosimplify" => simplify = false,
            "-banded" => banding = true,
            "-t" | "--threads" => {
                if i + 1 < args.len() {
                    num_threads = Some(args[i + 1].parse().expect("Invalid thread count"));
                    i += 1;
                } else {
                    eprintln!("Error: -t requires a value");
                    return;
                }
            }
            "-f" | "--format" => {
                if i + 1 < args.len() {
                    output_format = OutputFormat::from_str(&args[i + 1]).unwrap_or(OutputFormat::Default);
                    i += 1;
                } else {
                    eprintln!("Error: -f requires a format (delta, paf, sam)");
                    return;
                }
            }
            "-stats" | "--stats" => {
                show_stats = true;
            }
            arg if !arg.starts_with('-') => {
                if reference_file.is_empty() {
                    reference_file = arg;
                } else {
                    query_files.push(arg.to_string());
                }
            }
            _ => {
                eprintln!("Unknown option: {}", args[i]);
                print_nucmer_usage(&args[0]);
                return;
            }
        }
        i += 1;
    }
    
    if query_files.is_empty() {
        // If no query files provided, treat the second argument as the only query file
        if args.len() >= 3 && !args[2].starts_with('-') {
            query_files.push(args[2].clone());
        } else {
            eprintln!("Error: No query file provided");
            print_nucmer_usage(&args[0]);
            return;
        }
    }
    
    // Calculate and print statistics if requested
    if show_stats {
        let ref_sequences = parse_fasta(reference_file);
        let ref_stats = GenomicStats::new(&ref_sequences);
        ref_stats.print_stats("Reference");
        
        for query_file in &query_files {
            let query_sequences = parse_fasta(query_file);
            let query_stats = GenomicStats::new(&query_sequences);
            query_stats.print_stats("Query");
        }
    }
    
    // Set number of threads if specified
    if let Some(threads) = num_threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads)
            .build_global()
            .ok(); // Ignore errors if global pool is already initialized
    }
    
    // Create nucmer aligner with options
    let options = NucmerOptions {
        match_type: algorithm,
        min_len,
        forward_only,
        reverse_only,
        break_len,
        min_cluster,
        diag_diff,
        diag_factor,
        max_gap,
        extend,
        optimize,
        simplify,
        banding,
        use_extent: false,  // Not implemented yet
        to_seqend: !optimize,  // Inverse of optimize
        do_delta: true,      // Always true for nucmer
        do_shadows: !simplify,  // Inverse of simplify
    };
    
    let reference_seq = read_fasta_file(reference_file);
    
    // Process all query files in parallel
    let query_sequences: Vec<Vec<u8>> = query_files
        .iter()
        .map(|f| read_fasta_file(f))
        .collect();
    
    // Align all queries in parallel with progress bar
    let all_matches = align_multiple_sequences_parallel(
        &reference_seq,
        &query_sequences,
        options,
        num_threads,
    ).expect("Could not perform alignments");
    
    // Print matches for each query file in the specified format
    for (i, matches) in all_matches.iter().enumerate() {
        print_matches_in_format(matches, &query_files[i], &output_format, &reference_seq, &query_sequences[i]);
    }
}

fn read_fasta_file(filename: &str) -> Vec<u8> {
    let content = fs::read_to_string(filename)
        .expect("Could not read file");
    
    let mut sequence = Vec::new();
    for line in content.lines() {
        if !line.starts_with('>') {
            sequence.extend_from_slice(line.as_bytes());
        }
    }
    
    // Convert to uppercase and validate DNA sequence
    for base in &mut sequence {
        *base = match *base {
            b'a' | b'A' => b'A',
            b'c' | b'C' => b'C',
            b'g' | b'G' => b'G',
            b't' | b'T' => b'T',
            b'n' | b'N' => b'N',
            _ => b'N', // Default to N for non-standard bases
        };
    }
    
    sequence
}

fn print_usage(program: &str) {
    println!("Usage: {} [options] <reference-file> <query file1> [query file2] ...", program);
    println!("Options:");
    println!("  -mum           compute maximal matches that are unique in both sequences");
    println!("  -mumreference  compute maximal matches that are unique in the reference sequence (default)");
    println!("  -mumcand       same as -mumreference");
    println!("  -maxmatch      compute all maximal matches regardless of their uniqueness");
    println!("  -l <n>         set the minimum length of a match (default: 20)");
    println!("  -t, --threads <n>  number of threads to use (default: all available cores)");
    println!("  -f, --format <format>  output format (default, delta, paf, sam)");
    println!("  -stats         show reference and query sequence statistics (N50, N90, etc.)");
    println!();
    println!("Example:");
    println!("  {} -maxmatch -l 20 -t 4 -f paf reference.fa query.fa", program);
}

fn print_nucmer_usage(program: &str) {
    println!("Usage: {} [options] <reference-file> <query file1> [query file2] ...", program);
    println!("Options:");
    println!("  -mum           compute maximal matches that are unique in both sequences");
    println!("  -mumreference  compute maximal matches that are unique in the reference sequence (default)");
    println!("  -mumcand       same as -mumreference");
    println!("  -maxmatch      compute all maximal matches regardless of their uniqueness");
    println!("  -b, --breaklen <n>      set the distance an alignment extension will attempt to extend poor scoring regions before giving up (default: 200)");
    println!("  -c, --mincluster <n>    sets the minimum length of a cluster of matches (default: 65)");
    println!("  -D, --diagdiff <n>      set the maximum diagonal difference between two adjacent anchors in a cluster (default: 5)");
    println!("  -d, --diagfactor <f>    set the maximum diagonal difference between two adjacent anchors in a cluster as a differential fraction of the gap length (default: 0.12)");
    println!("  -noextend                do not perform cluster extension step");
    println!("  -f, --forward           use only the forward strand of the Query sequences");
    println!("  -g, --maxgap <n>        set the maximum gap between two adjacent matches in a cluster (default: 90)");
    println!("  -l, --minmatch <n>      set the minimum length of a single exact match (default: 20)");
    println!("  -L, --minalign <n>      minimum length of an alignment, after clustering and extension");
    println!("  -nooptimize              no alignment score optimization");
    println!("  -r, --reverse           use only the reverse complement of the Query sequences");
    println!("  -nosimplify              don't simplify alignments by removing shadowed clusters");
    println!("  -banded                  enforce absolute banding of dynamic programming matrix based on diagdiff parameter");
    println!("  -t, --threads <n>       number of threads to use (default: all available cores)");
    println!("  -f, --format <format>   output format (default, delta, paf, sam)");
    println!("  -stats                   show reference and query sequence statistics (N50, N90, etc.)");
    println!();
    println!("Example:");
    println!("  {} -maxmatch -l 20 -t 4 -f sam reference.fa query.fa", program);
}
