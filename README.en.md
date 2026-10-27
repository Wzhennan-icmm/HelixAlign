# HelixAlign: exact-match usage

This branch implements DNA maximal exact seeds, with explicit record and strand
identity. It does not implement full NUCmer clustering, gapped extension, protein
matching or GPU kernels. The independent audit is in
[audit/CORRECTNESS_AUDIT_20261027.zh.md](audit/CORRECTNESS_AUDIT_20261027.zh.md).

## Definitions

| Mode | Switch | Uniqueness after maximal extension |
|---|---|---|
| MUM | -mum | Once across all reference records and once within the current query record |
| MAM | -mumreference, -mumcand | Once across all reference records |
| MEM | -maxmatch | No uniqueness filter |

Every match is exact and cannot extend either left or right at its own pair of
positions. Overlapping occurrences count toward uniqueness. A shorter MEM on a
different diagonal is not removed simply because its intervals are contained
within a longer match's intervals. Query records and orientations are handled
separately, matching official mummer's scope.

## Build

~~~bash
cargo test --locked
cargo build --release --locked
~~~

The binary is target/release/helixalign. Cargo.lock is committed. The manifest
requires Rust >=1.85 (edition 2024); the audit ran on Rust 1.99.1 and did not
separately verify every supported older compiler.

## Commands

~~~bash
# Reference-unique maximal matches, query forward only (defaults).
target/release/helixalign -l 20 reference.fa query.fa

# Both-side unique matches, both orientations, canonical TSV.
target/release/helixalign -mum -l 20 -b --output-format matches reference.fa query.fa

# All MEMs, reverse only, PAF.
target/release/helixalign -maxmatch -l 20 -r --format paf reference.fa query.fa

# Process query records with an explicitly sized local thread pool.
target/release/helixalign -mumreference -l 20 -b -t 4 -f matches reference.fa q1.fa q2.fa
~~~

- -l/--minmatch must be a positive integer.
- Default strand is forward; --forward is explicit forward, -r/--reverse is
  reverse complement, and -b/--both requests both.
- -f/--format/--output-format accepts default, matches/tsv, or paf.
- -stats/--stats writes sequence statistics to stderr, keeping TSV/PAF clean.
- -n and -F are accepted compatibility flags. ACGT-only matching and reference
  identifiers are always enabled.
- -- ends option parsing and permits paths beginning with a hyphen.
- Unsupported switches and invalid options return exit status 2.

FASTA IDs are the first header token. Empty records, missing headers, illegal
characters and duplicate IDs are errors. Query IDs must also be unique across
multiple query input files. IUPAC letters are retained in coordinates but
non-ACGT letters never match. Wrapped lines, ASCII whitespace, CRLF and case
normalization are supported. No reference or query match crosses a record.

## Output

Canonical matches TSV has a comment header and six columns:

~~~text
#ref_id query_id strand ref_start query_start length
r       q        -      0         2           8
~~~

The actual separator is TAB. All starts are zero-based positions on the
original forward record. Even for the reverse strand, query_start is the
leftmost original-query position; the aligned interval is
[query_start, query_start + length). The default human format uses one-based
left positions and an explicit strand label.

PAF contains real record names, per-record lengths, zero-based half-open
intervals, and strand. MAPQ is 255 (unavailable); rows are exact seeds, not
confidence-calibrated complete read alignments.

SAM/delta output is refused because the former writers did not represent valid
alignment semantics. A binary invoked as nucmer or nucmer-rs is also refused.
No clustering or gapped alignment is inferred from exact seeds.

## Rust API migration

Prefer FastaRecord + ExactMatcher + RecordMatch + write_record_matches.
Core functions find_mums, find_mams, find_mems and run_mummer_algorithm remain
available for a single oriented query byte sequence; byte zero is a hard
nonmatching barrier. Core matching is byte-based: use ExactMatcher for DNA
normalization and record metadata.

SparseSuffixArray is a historical type name. It stores every suffix and now
rejects k != 1. Empty indices are supported. The core matching functions return
an empty vector for invalid min_len=0; the checked ExactMatcher and CLI reject it.

NucmerAligner is retained as a single-record exact-match wrapper, with no
NUCmer clustering or extension. Its align method now returns Vec<RecordMatch>
so direction is retained; align_parallel returns a Result and uses a local
thread pool. Changing legacy inert nucmer settings is rejected.

The old print_matches_in_format API was removed; callers must migrate to
write_record_matches with record context. These changes deliberately expose
incorrect old assumptions and are not a fully compatible API update.

## Validation and scalability

See [REPRODUCE.md](audit/REPRODUCE.md) for the 39 Rust tests and the three-way
differential runner. The finite corpus uses official MUMmer 4.0.1, min_len >=2,
-n and no -c. Known reference edge defects are documented; one-base matching
is covered by an independent exhaustive Rust oracle.

Suffix sorting still compares full suffixes, and adjacent LCP values are
calculated naively. Repeats can make construction and matching expensive.
This is a supplemental small-input correctness implementation. Use a mature,
version-pinned official CPU implementation with equivalent semantics as the
primary GPU performance denominator.
