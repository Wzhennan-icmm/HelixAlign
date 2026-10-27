//! CLI for record-aware DNA exact matching, not a full nucmer pipeline.
use helixalign::{
    ExactMatcher, FastaRecord, GenomicStats, MatchType, OutputFormat, StrandMode,
    read_fasta_records, write_record_matches,
};
use rayon::prelude::*;
use std::collections::HashSet;
use std::io::{self, Write};
use std::path::Path;

struct Options {
    min_len: usize,
    algorithm: MatchType,
    strands: StrandMode,
    threads: Option<usize>,
    format: OutputFormat,
    stats: bool,
    files: Vec<String>,
}

fn main() {
    if let Err(message) = run() {
        eprintln!("helixalign: {message}");
        std::process::exit(2);
    }
}

fn option_value<'a>(args: &'a [String], i: &mut usize) -> Result<&'a str, String> {
    let option = &args[*i];
    *i += 1;
    args.get(*i)
        .map(String::as_str)
        .ok_or_else(|| format!("{option} requires a value"))
}

fn positive(value: &str, name: &str) -> Result<usize, String> {
    value
        .parse::<usize>()
        .ok()
        .filter(|&n| n > 0)
        .ok_or_else(|| format!("{name} must be a positive integer"))
}

fn parse_options(args: &[String]) -> Result<Options, String> {
    let mut options = Options {
        min_len: 20,
        algorithm: MatchType::MAM,
        strands: StrandMode::Forward,
        threads: None,
        format: OutputFormat::Default,
        stats: false,
        files: Vec::new(),
    };
    let mut algorithm = None;
    let mut strands = None;
    let mut positional_only = false;
    let mut i = 1;
    while i < args.len() {
        let arg = args[i].as_str();
        if positional_only {
            options.files.push(args[i].clone());
            i += 1;
            continue;
        }
        match arg {
            "--" => positional_only = true,
            "-mum" | "-mumreference" | "-mumcand" | "-maxmatch" => {
                let selected = match arg {
                    "-mum" => MatchType::MUM,
                    "-maxmatch" => MatchType::MEM,
                    _ => MatchType::MAM,
                };
                if algorithm.is_some_and(|a| a != selected) {
                    return Err("Choose only one matching mode".into());
                }
                algorithm = Some(selected);
                options.algorithm = selected;
            }
            "-b" | "--both" | "-r" | "--reverse" | "--forward" => {
                let selected = match arg {
                    "-b" | "--both" => StrandMode::Both,
                    "-r" | "--reverse" => StrandMode::Reverse,
                    _ => StrandMode::Forward,
                };
                if strands.is_some_and(|s| s != selected) {
                    return Err("Choose only one strand mode".into());
                }
                strands = Some(selected);
                options.strands = selected;
            }
            "-l" | "--minmatch" => {
                options.min_len = positive(option_value(args, &mut i)?, "Minimum match length")?
            }
            "-t" | "--threads" => {
                options.threads = Some(positive(option_value(args, &mut i)?, "Thread count")?)
            }
            "-f" | "--format" | "--output-format" => {
                let value = option_value(args, &mut i)?;
                options.format = OutputFormat::from_str(value)
                    .ok_or_else(|| format!("Unknown output format: {value}"))?;
                options.format.validate_exact_output()?;
            }
            "-stats" | "--stats" => options.stats = true,
            // DNA-only semantics and reference identifiers are always enabled.
            "-n" | "-F" => {}
            value if value.starts_with('-') => {
                return Err(format!(
                    "Unknown or unsupported option: {value}. See --help."
                ));
            }
            _ => options.files.push(args[i].clone()),
        }
        i += 1;
    }
    if options.files.len() < 2 {
        return Err("Expected a reference FASTA and at least one query FASTA. See --help.".into());
    }
    Ok(options)
}

fn run() -> Result<(), String> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() == 2 && matches!(args[1].as_str(), "-h" | "--help" | "-help") {
        print_usage();
        return Ok(());
    }
    if args.len() == 2 && matches!(args[1].as_str(), "--version" | "-version") {
        println!("helixalign {}", env!("CARGO_PKG_VERSION"));
        return Ok(());
    }
    if Path::new(&args[0])
        .file_name()
        .and_then(|s| s.to_str())
        .is_some_and(|s| s == "nucmer" || s == "nucmer-rs")
    {
        return Err("This binary computes exact matches only; nucmer clustering and gapped extension are not implemented. Invoke helixalign directly.".into());
    }
    let options = parse_options(&args)?;
    let references = read_fasta_records(&options.files[0])?;
    let mut queries = Vec::new();
    let mut ids = HashSet::new();
    for filename in &options.files[1..] {
        for record in read_fasta_records(filename)? {
            if !ids.insert(record.id.clone()) {
                return Err(format!(
                    "Duplicate query ID across input files: {}",
                    record.id
                ));
            }
            queries.push(record);
        }
    }
    if options.stats {
        print_stats("Reference", &references);
        print_stats("Query", &queries);
    }
    let matcher = ExactMatcher::new(references)?;
    // Reuse one reference index. A local pool honors each call's thread setting.
    let available = std::thread::available_parallelism().map_or(1, usize::from);
    let threads = options
        .threads
        .unwrap_or_else(|| queries.len().min(available).max(1));
    let pool_builder = rayon::ThreadPoolBuilder::new().num_threads(threads);
    let pool = pool_builder
        .build()
        .map_err(|e| format!("Cannot create thread pool: {e}"))?;
    let results: Vec<_> = pool.install(|| {
        queries
            .par_iter()
            .map(|query| {
                matcher.find_matches(query, options.algorithm, options.min_len, options.strands)
            })
            .collect()
    });
    let results: Vec<_> = results.into_iter().collect::<Result<_, _>>()?;
    let stdout = io::stdout();
    let mut output = io::BufWriter::new(stdout.lock());
    if options.format == OutputFormat::Matches {
        writeln!(
            output,
            "#ref_id\tquery_id\tstrand\tref_start\tquery_start\tlength"
        )
        .map_err(|e| e.to_string())?;
    }
    for (query, matches) in queries.iter().zip(results.iter()) {
        write_record_matches(
            &mut output,
            matches,
            matcher.references(),
            query,
            options.format,
        )
        .map_err(|e| e.to_string())?;
    }
    output.flush().map_err(|e| e.to_string())
}

fn print_stats(label: &str, records: &[FastaRecord]) {
    let sequences: Vec<_> = records.iter().map(|r| r.sequence.clone()).collect();
    let stats = GenomicStats::new(&sequences);
    eprintln!(
        "{label}: records={} bases={} N50={} N90={} GC={:.2}%",
        stats.num_sequences, stats.total_length, stats.n50, stats.n90, stats.gc_content
    );
}

fn print_usage() {
    println!(
        "HelixAlign: DNA maximal exact matches (not full nucmer alignments)\n\
Usage: helixalign [options] reference.fa query.fa [query2.fa ...]\n\n\
  -mum                 unique across reference records and within this query record\n\
  -mumreference        unique across reference records (default); -mumcand is an alias\n\
  -maxmatch            all maximal exact matches (MEMs)\n\
  -l, --minmatch N     positive minimum match length (default: 20)\n\
  --forward            query forward strand only (default)\n\
  -r, --reverse        query reverse complement only\n\
  -b, --both           both query orientations, separately\n\
  -t, --threads N      positive number of query-processing threads\n\
  -f, --format FMT     default, matches (canonical TSV), or paf\n\
  --output-format FMT  alias for --format\n\
  -stats, --stats      record statistics on stderr\n\
  -n, -F              accepted; ACGT-only matching and reference IDs are always enabled\n\n\
FASTA record boundaries and non-ACGT IUPAC bases are nonmatching barriers.\n\
TSV positions are zero-based starts in the original forward records; strand is explicit.\n\
PAF rows describe exact seeds, with MAPQ 255 (unavailable).\n\
SAM, delta, proteins, sparse sampling and nucmer clustering/extension are unsupported."
    );
}

#[cfg(test)]
mod tests {
    use super::*;
    fn args(values: &[&str]) -> Vec<String> {
        values.iter().map(|s| s.to_string()).collect()
    }
    #[test]
    fn rejects_invalid_parameters_and_unimplemented_formats() {
        for flags in [
            vec!["-l", "0"],
            vec!["-l", "bad"],
            vec!["-t", "0"],
            vec!["-l"],
            vec!["-b", "-r"],
            vec!["-mum", "-maxmatch"],
            vec!["--format", "delta"],
            vec!["--format", "sam"],
            vec!["--format", "unknown"],
            vec!["-c", "10"],
        ] {
            let mut values = vec!["helixalign", "ref.fa", "query.fa"];
            values.extend(flags);
            assert!(
                parse_options(&args(&values)).is_err(),
                "accepted {values:?}"
            );
        }
    }
}
