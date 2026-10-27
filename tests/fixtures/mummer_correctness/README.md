# Exact-match correctness fixtures

These fixtures test the **exact-match seed layer**, not NUCMER clustering,
gapped extension, alignment scoring, or performance. The manifest is
[`cases.json`](cases.json); each case has a small reference/query FASTA pair
and an explicit minimum length. All expected results are independently
computable by enumerating pairs of sequence positions.

## Contract

- Inputs are case-insensitive DNA. Only A, C, G, T may participate in a match;
  any other accepted IUPAC symbol is a barrier (official `mummer -n`).
- A MEM must be exact and cannot extend left or right **at that pair of
  positions**. Containment in another match's coordinate rectangle does not
  make it nonmaximal: another match can be on a different diagonal.
- A MAM is a MEM whose complete matched string occurs once across **all
  reference records**.
- A MUM additionally requires one occurrence within **the current query
  record**, counting overlapping occurrences. Identical strings in separate
  query records do not invalidate MUMs under official MUMmer semantics.
- Query forward and reverse-complement strands are processed separately.
  A palindromic match can therefore occur on both strands.
- No match can cross a FASTA record boundary. Output retains record IDs.
- Tuples are `(ref_id, query_id, strand, ref_start, query_start, length)`.
  Starts are zero-based positions in the original forward records on both
  strands. For an official reverse result `(r1, q_rc1, L)`, normalize as
  `r0 = r1 - 1`, `q0 = query_length - (q_rc1 - 1) - L`.

## Hand-calculated forward examples

In the compact table below, `(r, q, L)` means zero-based reference start,
zero-based query start, and length, on the forward strand.

| Fixture | Reference / query | Minimum | MEM | MAM | MUM |
|---|---|---:|---|---|---|
| `query_repeat` | `ACGT` / `ACGTACGT` | 4 | `(0,0,4)`, `(0,4,4)` | same | empty |
| `both_unique` | `AACCGGTA` / same | 8 | `(0,0,8)` | same | same |
| `reference_repeat` | `ACGTACGT` / `ACGT` | 4 | `(0,0,4)`, `(4,0,4)` | empty | empty |
| `overlapping_query` | `AA` / `AAAA` | 2 | `(0,0,2)`, `(0,1,2)`, `(0,2,2)` | same | empty |
| `homopolymer` | `AAAA` / same | 2 | `(0,0,4)`, `(0,1,3)`, `(0,2,2)`, `(1,0,3)`, `(2,0,2)` | `(0,0,4)` | `(0,0,4)` |
| `nested_diagonals` | `ACAC` / same | 2 | `(0,0,4)`, `(0,2,2)`, `(2,0,2)` | `(0,0,4)` | `(0,0,4)` |
| `no_match` | `AA` / `CC` | 2 | empty | empty | empty |
| `lexicographic_mismatch` | `CAT` / `AA` | 2 | empty | empty | empty |
| `suffix_overrun` | `AC` / `ACG` | 2 | `(0,0,2)` | same | same |
| `reference_boundary` | `r1:AAAA`, `r2:CCCC` / `AACC` | 4 | empty | empty | empty |
| `query_boundary` | `AACC` / `q1:AA`, `q2:CC` | 4 | empty | empty | empty |
| `reference_global_unique` | `r1:ACGT`, `r2:ACGT` / `ACGT` | 4 | one full match in each reference | empty | empty |
| `query_record_unique` | `ACGT` / `q1:ACGT`, `q2:ACGT` | 4 | one full match in each query | same | same |

`reverse_offset` has reference `AACCGGTA` and query `CCTACCGGTTG`, with
minimum 8. Its sole match is `("r", "q", "-", 0, 2, 8)` in each mode.
`palindrome` (`ACGT` / `ACGT`, minimum 4) has both
`("r", "q", "+", 0, 0, 4)` and `("r", "q", "-", 0, 0, 4)` when both
strands are requested; these are distinct results.

The remaining fixtures cover mixed case, ambiguity barriers, too-short input,
terminal suffixes, internal mismatched flanks, overlapping tandem copies, and
several reference/query records with mixed directions.

## Differential runner

From the repository root, build HelixAlign and run:

```sh
cargo build --release
python3 scripts/differential_mummer.py \
  --mummer /absolute/path/to/official/mummer \
  --candidate target/release/helixalign \
  --random-cases 120 --seed 20261027 \
  --report audit/results/differential.json
```

The runner executes every fixture and every generated case under all nine
combinations of `MUM/MAM/MEM` and `forward/reverse/both`. Its independent
brute-force oracle does not call either binary. It compares complete tuples
as multisets, so duplicate output is also a failure. It additionally checks
returned tuples for exactness, left/right maximality, uniqueness, boundaries,
length threshold, and valid strand.

Use `--baseline /path/to/original/helixalign` to test the unmodified upstream
binary. The baseline is invoked only in its supported forward mode. Its
flattened coordinates are mapped back to records; boundary-crossing and
out-of-bounds matches are retained as invalid results, rather than discarded.
Baseline failures document regressions and do not by themselves cause the
runner to exit nonzero. Any official-reference or patched-candidate mismatch
does cause a nonzero exit status.

The JSON report preserves all generated sequences and the random seed, plus
the false positives, false negatives, and invalidity reasons. Commands within
the report identify temporary `ref.fa` and `query.fa` inputs; recreate those
inputs from the named case in the report when reproducing an individual check.

## Official MUMmer version caveats

The audit uses the unmodified official **MUMmer 4.0.1** release. All differential
cases use minimum length at least 2. Its `-l 1` behavior disagrees with the
mathematical oracle (including omissions), so length-one correctness is tested
separately against the brute oracle in Rust. This restriction must not be
silently presented as full equivalence at every possible parameter value.

The runner uses `-F -n`, explicit `-mum` / `-mumreference` / `-maxmatch`, and
`-b` / `-r` for requested directions. It deliberately omits `-c`: a verified
MUMmer 4.0.1 reverse-coordinate conversion case truncates the returned length
when `-c` is requested. Reverse coordinates are converted explicitly above.
