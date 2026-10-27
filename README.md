# HelixAlign

HelixAlign is a Rust implementation of **DNA maximal exact matching** with MUM,
MAM and MEM modes. It is an independent small-input correctness implementation,
not a GPU aligner or a complete NUCMER replacement.

The 2026-10-27 audit confirmed correctness defects in the original revision
4835117. This branch repairs suffix lookup, pairwise maximality, both-side MUM
uniqueness, FASTA boundaries and strand-aware output.

- [English usage and limitations](README.en.md)
- [中文使用说明与兼容性变化](README.zh.md)
- [独立正确性审核报告与实际结果](audit/CORRECTNESS_AUDIT_20261027.zh.md)
- [Reproduce the audit](audit/REPRODUCE.md)
- [Hand-checkable FASTA fixtures](tests/fixtures/mummer_correctness/README.md)

## Build and use

Requires Rust edition-2024 support. The recorded run used Rust 1.99.1.

~~~bash
cargo test --locked
cargo build --release --locked
target/release/helixalign -mum -l 20 -b --output-format matches reference.fa query.fa
~~~

Canonical TSV columns are ref_id, query_id, strand, ref_start, query_start,
length. Starts are zero-based on the original forward records, including
reverse-strand matches.

## Correctness evidence and limits

The final Rust suite passes 39 tests. For 24 designed and 120 deterministic
random FASTA pairs, all 1,296 candidate comparisons and 1,296 official MUMmer
4.0.1 comparisons agree with an independent brute-force oracle. Detailed
versioned reports and known reference-version caveats are under audit/.

Reference uniqueness is global across reference records; query uniqueness
is local to the current query record and orientation. Only A/C/G/T participate
in matching; other accepted IUPAC symbols and record boundaries stop matches.

This finite suite does not establish whole-genome scalability or complete
MUMmer/NUCmer compatibility. Full suffix sorting and naive LCP construction
remain unsuitable for a mature genome-scale performance baseline. No GPU
speedup is claimed. SAM/delta, proteins, sparse sampling and nucmer
clustering/gapped extension are unsupported and are not silently simulated.
