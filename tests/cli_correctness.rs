//! End-to-end tests of the executable, using hand-calculated exact matches.
//! These exercise parsing, file failures, record identity and serialization;
//! they intentionally do not call the production matching API as an oracle.

use std::fs;
use std::path::PathBuf;
use std::process::{Command, Output};
use std::sync::atomic::{AtomicU64, Ordering};

const TSV_HEADER: &str = "#ref_id\tquery_id\tstrand\tref_start\tquery_start\tlength\n";
static NEXT_CASE: AtomicU64 = AtomicU64::new(0);

struct Case {
    dir: PathBuf,
}

impl Case {
    fn new() -> Self {
        let nanos = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        let dir = std::env::temp_dir().join(format!(
            "helixalign-cli-{}-{nanos}-{}",
            std::process::id(),
            NEXT_CASE.fetch_add(1, Ordering::Relaxed)
        ));
        fs::create_dir(&dir).unwrap();
        Self { dir }
    }

    fn write(&self, name: &str, text: &str) {
        fs::write(self.dir.join(name), text).unwrap();
    }

    fn run(&self, flags: &[&str], files: &[&str]) -> Output {
        Command::new(env!("CARGO_BIN_EXE_helixalign"))
            .current_dir(&self.dir)
            .args(flags)
            .args(files)
            .output()
            .unwrap()
    }
}

impl Drop for Case {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.dir);
    }
}

fn stdout(output: &Output) -> &str {
    std::str::from_utf8(&output.stdout).unwrap()
}

fn stderr(output: &Output) -> &str {
    std::str::from_utf8(&output.stderr).unwrap()
}

fn succeeds(output: &Output) {
    assert!(output.status.success(), "{}", stderr(output));
}

fn fails_cleanly(output: &Output, message: &str) {
    assert_eq!(output.status.code(), Some(2), "{}", stderr(output));
    assert!(
        output.stdout.is_empty(),
        "unexpected output: {}",
        stdout(output)
    );
    assert!(stderr(output).contains("helixalign:"), "{}", stderr(output));
    assert!(stderr(output).contains(message), "{}", stderr(output));
}

fn records_case() -> Case {
    let case = Case::new();
    case.write(
        "ref.fa",
        ">decoy\nGGGGCCCC\n>target description\nAACCGGTA\n",
    );
    case.write("forward.fa", ">q_forward\nTTAACCGGTACC\n");
    case.write("reverse.fa", ">q_reverse\nCCTACCGGTTG\n");
    case
}

#[test]
fn mode_flags_distinguish_mum_from_reference_only_uniqueness() {
    let case = Case::new();
    case.write("ref.fa", ">r\nACG\n");
    case.write("query.fa", ">q\nACGACG\n");
    let rows = "r\tq\t+\t0\t0\t3\nr\tq\t+\t0\t3\t3\n";
    for mode in ["-mumreference", "-mumcand", "-maxmatch"] {
        let output = case.run(&[mode, "-l", "3", "-f", "matches"], &["ref.fa", "query.fa"]);
        succeeds(&output);
        assert_eq!(stdout(&output), format!("{TSV_HEADER}{rows}"));
        assert!(output.stderr.is_empty());
    }
    let output = case.run(
        &["-mum", "-l", "3", "-f", "matches"],
        &["ref.fa", "query.fa"],
    );
    succeeds(&output);
    assert_eq!(stdout(&output), TSV_HEADER);
}

#[test]
fn mem_output_keeps_maximal_matches_on_different_diagonals() {
    let case = Case::new();
    case.write("ref.fa", ">r\nATAT\n");
    case.write("query.fa", ">q\nATAT\n");
    let output = case.run(
        &["-maxmatch", "-l", "2", "-f", "matches"],
        &["ref.fa", "query.fa"],
    );
    succeeds(&output);
    assert_eq!(
        stdout(&output),
        format!("{TSV_HEADER}r\tq\t+\t0\t0\t4\nr\tq\t+\t0\t2\t2\nr\tq\t+\t2\t0\t2\n")
    );
}

#[test]
fn multiple_files_keep_record_local_coordinates_and_strand_in_tsv_and_paf() {
    let case = records_case();
    let files = ["ref.fa", "forward.fa", "reverse.fa"];
    let tsv = case.run(&["-mum", "-b", "-l", "8", "-f", "matches"], &files);
    succeeds(&tsv);
    assert_eq!(
        stdout(&tsv),
        format!("{TSV_HEADER}target\tq_forward\t+\t0\t2\t8\ntarget\tq_reverse\t-\t0\t2\t8\n")
    );
    let paf = case.run(&["-mum", "-b", "-l", "8", "-f", "paf"], &files);
    succeeds(&paf);
    // PAF uses original query coordinates, half-open intervals and unknown MAPQ.
    assert_eq!(
        stdout(&paf),
        concat!(
            "q_forward\t12\t2\t10\t+\ttarget\t8\t0\t8\t8\t8\t255\n",
            "q_reverse\t11\t2\t10\t-\ttarget\t8\t0\t8\t8\t8\t255\n"
        )
    );
}

#[test]
fn forward_is_default_and_reverse_and_both_flags_select_the_expected_hits() {
    let case = records_case();
    let files = ["ref.fa", "forward.fa", "reverse.fa"];
    let forward = "target\tq_forward\t+\t0\t2\t8\n";
    let reverse = "target\tq_reverse\t-\t0\t2\t8\n";
    for (extra, rows) in [
        (vec![], forward.to_string()),
        (vec!["--forward"], forward.to_string()),
        (vec!["-r"], reverse.to_string()),
        (vec!["-b"], format!("{forward}{reverse}")),
    ] {
        let mut flags = vec!["-l", "8", "-f", "matches"];
        flags.extend(extra);
        let output = case.run(&flags, &files);
        succeeds(&output);
        assert_eq!(stdout(&output), format!("{TSV_HEADER}{rows}"));
    }
}

#[test]
fn neither_reference_nor_query_records_can_create_junction_hits() {
    let case = Case::new();
    case.write("split.fa", ">first\nAAAA\n>second\nCCCC\n");
    case.write("joined.fa", ">joined\nAAAACCCC\n");
    for mode in ["-mum", "-mumreference", "-maxmatch"] {
        for files in [["split.fa", "joined.fa"], ["joined.fa", "split.fa"]] {
            let output = case.run(&[mode, "-b", "-l", "8", "-f", "matches"], &files);
            succeeds(&output);
            assert_eq!(stdout(&output), TSV_HEADER);
        }
    }
}

#[test]
fn thread_counts_preserve_order_and_stats_only_use_stderr() {
    let case = records_case();
    let files = ["ref.fa", "forward.fa", "reverse.fa"];
    let serial = case.run(&["-b", "-l", "8", "-f", "matches", "-t", "1"], &files);
    let parallel = case.run(
        &["-b", "-l", "8", "-f", "matches", "-t", "2", "--stats"],
        &files,
    );
    succeeds(&serial);
    succeeds(&parallel);
    assert_eq!(stdout(&parallel), stdout(&serial));
    assert!(serial.stderr.is_empty());
    assert!(stderr(&parallel).contains("Reference: records=2 bases=16"));
    assert!(stderr(&parallel).contains("Query: records=2 bases=23"));
    assert!(!stdout(&parallel).contains("Reference:"));
    assert!(!stdout(&parallel).contains("Query:"));
}

#[test]
fn malformed_empty_duplicate_and_missing_files_fail_without_partial_matches() {
    let case = Case::new();
    case.write("ref.fa", ">r\nACGT\n");
    for (bad_fasta, message) in [
        ("", "contains no records"),
        ("ACGT\n", "before first FASTA header"),
        (">\nACGT\n", "Missing FASTA identifier"),
        (">q\n", "empty"),
        (">q\n>another\nACGT\n", "Empty FASTA record"),
        (">q\nACGT\n>q\nACGT\n", "Duplicate FASTA identifier"),
        (">q\nAC?T\n", "Invalid DNA FASTA byte"),
    ] {
        case.write("bad.fa", bad_fasta);
        let output = case.run(&["-l", "4", "-f", "matches"], &["ref.fa", "bad.fa"]);
        fails_cleanly(&output, message);
    }
    let missing_query = case.run(&["-l", "4"], &["ref.fa", "missing.fa"]);
    fails_cleanly(&missing_query, "Cannot read missing.fa");
    let missing_reference = case.run(&["-l", "4"], &["missing.fa", "ref.fa"]);
    fails_cleanly(&missing_reference, "Cannot read missing.fa");
    case.write("first.fa", ">same\nACGT\n");
    case.write("second.fa", ">same\nACGT\n");
    let duplicated = case.run(
        &["-l", "4", "-f", "matches"],
        &["ref.fa", "first.fa", "second.fa"],
    );
    fails_cleanly(&duplicated, "Duplicate query ID across input files");
    // A good earlier query must not be emitted if a later file is unreadable.
    let late_failure = case.run(
        &["-l", "4", "-f", "matches"],
        &["ref.fa", "first.fa", "missing.fa"],
    );
    fails_cleanly(&late_failure, "Cannot read missing.fa");
}

#[test]
fn invalid_modes_lengths_threads_and_unimplemented_formats_have_nonzero_exit() {
    let case = records_case();
    for (flags, message) in [
        (vec!["-l", "0"], "positive integer"),
        (vec!["-l", "invalid"], "positive integer"),
        (vec!["-t", "0"], "positive integer"),
        (vec!["-b", "-r"], "one strand mode"),
        (vec!["-mum", "-maxmatch"], "one matching mode"),
        (vec!["-f", "sam"], "not implemented correctly"),
        (vec!["-f", "delta"], "not implemented correctly"),
        (vec!["-f", "unknown"], "Unknown output format"),
        (vec!["-c", "10"], "unsupported option"),
    ] {
        let output = case.run(&flags, &["ref.fa", "forward.fa"]);
        fails_cleanly(&output, message);
    }
}

#[test]
fn end_of_options_allows_filenames_starting_with_a_dash() {
    let case = Case::new();
    case.write("-ref.fa", ">r\nACG\n");
    case.write("-query.fa", ">q\nACG\n");
    let output = case.run(
        &["-mum", "-l", "3", "-f", "matches", "--"],
        &["-ref.fa", "-query.fa"],
    );
    succeeds(&output);
    assert_eq!(stdout(&output), format!("{TSV_HEADER}r\tq\t+\t0\t0\t3\n"));
}

#[cfg(unix)]
#[test]
fn nucmer_alias_fails_instead_of_claiming_a_gapped_alignment() {
    let case = records_case();
    let alias = case.dir.join("nucmer");
    std::os::unix::fs::symlink(env!("CARGO_BIN_EXE_helixalign"), &alias).unwrap();
    let output = Command::new(alias)
        .current_dir(&case.dir)
        .args(["-l", "8", "ref.fa", "forward.fa"])
        .output()
        .unwrap();
    fails_cleanly(&output, "exact matches only");
}
