# Independent audit of legacy alignment/output claims

Reviewed original HelixAlign commit `483511789035539b3cc265356af6b70e3b8f8cb9`, preserved at `HelixAlign-original/`. This document records **static code evidence**, not throughput measurements or external format-reader executions. Line numbers below refer to that original version.

## Findings

| Severity | Evidence in original code | Consequence |
|---|---|---|
| P1 | `src/suffix_array.rs:8–12`: `Match` contains only `ref_pos`, `query_pos`, `len`. `src/nucmer.rs:94–104` converts reverse-complement coordinates and merges them with forward matches. | Reverse orientation is irreversibly lost from returned data. Downstream consumers cannot infer the strand reliably; identical coordinates on both strands are possible. |
| P1 | `src/nucmer.rs:66–107`: `align()` only calls `run_mummer_algorithm` on forward and reverse-complement query and returns these matches. The clustering/extension fields at lines 12–24 are never read by matching. | This is an exact-anchor finder, not full NUCmer. No anchor clustering, chain extraction, dynamic-programming/gapped extension, mismatch alignment, shadow filtering, or delta generation is implemented. |
| P1 | `src/output_format.rs:38–58`: outputs `NUCMER` twice, then `> Reference Query`; each match row ends with reference length, query length, match length; no terminating zero is written. | The claimed delta output does not follow the delta grammar. It lacks the source-path header, sequence IDs and lengths in the sequence header, correct three error-count fields, and delta terminators. It cannot be advertised as nucmer-compatible. |
| P1 | `src/output_format.rs:72`: PAF strand hardcoded `+`; lines 67, 74 use filename and generic `reference`; line 81 hardcodes MAPQ 60. | Reverse mappings are mislabeled; contig identity is unavailable; claimed mapping confidence is uncalibrated. Exact segment length itself can legitimately populate PAF matching/block lengths, but confidence should be unavailable (`255`) unless actually estimated. |
| P1 | `src/output_format.rs:98–107`: SAM QNAME is whole filename, FLAG is always 0, RNAME generic, MAPQ 60, CIGAR `<len>M`, SEQ only the matched slice. | Reverse alignments lack FLAG 0x10 and reverse-complemented SEQ. Query position and clipping outside the matched anchor are not represented. Multiple anchors of the same query are all primary and have varying truncated SEQ under the same QNAME. These rows do not faithfully represent alignment of the original query. |
| P1 | `src/output_format.rs:91–93` prints SAM `@HD`/`@SQ` on every call; `src/main.rs:345–346` and 119–132 call the printer per query file. | With multiple query files, headers can appear again after alignment records, violating the required header-before-alignment organization. |
| P2 | `src/main.rs:13–16` chooses NUCmer mode by `argv[0].contains("nucmer")`. Normal mode at 118–132 has only forward matching. | A name/path changes behavior; reverse matching is absent from normal matching path. This is not a reliable explicit mode contract. |
| P2 | `src/main.rs:205` handles `-f` as forward-only; line 247 repeats `-f` as format, shadowed by the earlier branch. Lines 224–228 explicitly parse `--minalign` and discard it. | Documented options can be ignored or interpreted differently. Plain option acceptance is not implementation of their algorithmic effect. |
| P2 | `src/nucmer.rs:57–63` accepts contradictory forward-only and reverse-only options; lines 70 and 81 then skip both branches. | Contradictory configuration silently returns no matches instead of rejecting the configuration. |
| P2 | `src/nucmer.rs:83` converts query bytes via `from_utf8(...).unwrap_or("")`. | Invalid UTF-8 input silently becomes an empty reverse query. A DNA API should validate supported bytes, rather than hide conversion failure. |

### Correct delta outline

The official delta format consists of source filenames, the `NUCMER`/`PROMER` tag, a sequence-pair header containing sequence IDs and lengths, then alignment headers and signed delta entries terminating with zero. An alignment header's last three integers count alignment errors, similarity errors, and stop codons. They are not sequence lengths or exact match length. Reverse alignments use descending query endpoints in original coordinates.

Source specification: https://mummer.sourceforge.net/manual/ , NUCmer output format.

### Relevant output specifications

- PAF: https://github.com/lh3/miniasm/blob/master/PAF.md . Query/target names refer to sequences; coordinates are zero-based half-open; strand is `+` or `-`; MAPQ 255 means unavailable.
- SAM v1.6: https://samtools.github.io/hts-specs/SAMv1.pdf . Header precedes alignment section, FLAG 0x10 denotes reverse-complemented SEQ, clipping can preserve query offset, and primary/secondary/supplementary distinctions must be represented when applicable.

## Smallest safe treatment of legacy public API

The exact-match core's `Match` can reasonably retain coordinates relative to a caller-provided, already oriented sequence. The bug occurs when `NucmerAligner` mixes orientations without preserving that context.

Two explicit solutions are possible:

1. Extend `Match` with an orientation field initialized to forward in `Match::new`, and set it on every reverse match in `NucmerAligner`. Document zero-based original-query left coordinates and include strand in equality/deduplication. This changes public struct literal construction and requires regression tests for forward, reverse, and palindromic matches.
2. Introduce an oriented return type for `NucmerAligner` or replace/deprecate the misleading API in favor of the new record-aware exact matcher. Merely adding a new safe method while leaving the old mixed-strand method silently lossy is not a complete fix.

Either approach should reject `forward_only && reverse_only` and minimum length zero at construction. The public documentation must state that this wrapper performs exact matches only; unused NUCmer parameters should be rejected/removed/deprecated explicitly, rather than presented as working features. Making the CLI fail clearly on unsupported full NUCmer/delta/SAM modes is a suitable scope boundary for this correctness patch.

## Performance and scale claims are not established

The following are direct static observations:

- `src/suffix_array.rs:40–48` creates all suffix positions and sorts by comparing full suffix byte slices. This is a full suffix array, not a sampled one.
- `k` is accepted and stored (`35–58`) and returned by a getter (`232–235`), but does not select every k-th suffix. Changing it cannot provide the documented sparse-memory benefit in the original implementation.
- `src/suffix_array.rs:70–86` recalculates each adjacent LCP from its first byte. For a long homopolymer, this alone requires a quadratic number of byte comparisons. The LCP array is not used by `search()` (`171–197`).
- All original match modes enumerate query starts and every remaining substring length (`src/algorithms.rs:20–26`, `49–55`, `77–83`). For query length m, the number of pattern searches is quadratic in m before accounting for comparison length, repeated occurrences, or candidate storage.
- Every candidate is materialized, and redundancy filtering scans previously retained candidates (`algorithms.rs:108–123`). This adds potentially quadratic work in the candidate count and avoids none of the upfront candidate generation.
- On a 64-bit target, the reference sequence and two `Vec<usize>` arrays alone have approximately `N + 8N + 8N = 17N` payload bytes, before caller copies, sorting temporaries, query indexes, candidates and allocator overhead. This is a static payload estimate, not measured peak RSS.
- Normal CLI rereads/reindexes reference for each query file (`main.rs:118–129`). Setting Rayon global thread count does not parallelize this loop or the matching core.
- Legacy NUCmer parallelism is across supplied query vectors only (`nucmer.rs:126–133`, `163–170`), and normal CLI currently assembles one vector per input file. Thread-pool construction errors are ignored (`113–116`, `148–151`), so a requested thread count need not take effect after earlier global initialization.
- The tracked source inventory has no CUDA/OpenCL/GPU kernels or GPU dependency configuration. `Cargo.toml:9–12` lists clap, Rayon and indicatif. HelixAlign is a separate CPU implementation.

`README.en.md:35`, `128–135`, `168–172` and corresponding Chinese sections claim NUCmer feature parity, sparse memory savings, faster construction and better threading without accompanying reproducible benchmark data in the reviewed repository. These claims cannot establish either feature equivalence or speed advantage.

### Implication for GPU development

After exact-match correctness is repaired and independently tested, HelixAlign may serve as an additional small-input CPU correctness implementation. It must not be the sole performance denominator. A GPU speedup against this simple CPU implementation would confound parallel acceleration with weaker CPU algorithms and different workloads. Compare equivalent exact-match stages against pinned official `mummer` with matching input, strand, alphabet, minimum length and uniqueness semantics. Benchmark the complete `nucmer` pipeline only after clustering and gapped extension are actually implemented and validated.

## Independent confirmation of upstream reference caveats

For official MUMmer v4.0.1 at commit `eb734606f2d516f42a0e0dce7a116bfb88ec1ebf`:

- `include/mummer/sparseSA.hpp:172–174` truncates match length to `max - q`; lines 567 and 725 pass a flipped coordinate when `-c` is used. This explains the directly observed reverse example's incorrect length 2 instead of 8. Without `-c`, the example reports the correct length and RC coordinate.
- `include/mummer/sparseSA.hpp:559` skips MAM traversal depth <=1, inherited by MUM filtering; line 641 skips MEM's minimum-length interval depth <=1 before extending it. This explains both missed one-base matches and MEM's no-output behavior for `-l 1`, even on an identical longer string.
- Full commands/results are in `REFERENCE_SEMANTICS.md`, `semantics_results.json` and `single_base_results.json`.
- These observations concern the specific pinned official executable/options and should not be generalized into a claim that all MUMmer matching is incorrect. Differential comparisons should run without `-c`, with minimum length >=2; the independent brute-force oracle handles excluded edges.
