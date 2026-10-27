# Official MUMmer reference: build and semantics audit

Date: 2026-10-27 (Asia/Shanghai).

## Executed reference

- Upstream: https://github.com/mummer4/mummer
- Release: **v4.0.1**, `mummer -version` reports `4.0.1`.
- Immutable commit: `eb734606f2d516f42a0e0dce7a116bfb88ec1ebf`.
- Release tarball: https://github.com/mummer4/mummer/releases/download/v4.0.1/mummer-4.0.1.tar.gz
- Tarball SHA256: `bc20ae2701a0b2e323e4e515b7cfa18a0f0cb34a4ff5844b289b2de0154e3d3e`.
- Compiler: GCC/G++ 13.3.0, Linux x86_64.
- Built the `mummer` target and its library from unmodified upstream source; did not run the complete upstream package test suite.
- Executable in this workspace: `/workspace/scratch/4f9aebc1eabf/helixalign-audit/reference/mummer-4.0.1/mummer`. It is a libtool wrapper; retain the build tree.

Reproduction:

```bash
curl -L --fail https://github.com/mummer4/mummer/releases/download/v4.0.1/mummer-4.0.1.tar.gz -o mummer-4.0.1.tar.gz
sha256sum mummer-4.0.1.tar.gz
tar --no-same-owner -xzf mummer-4.0.1.tar.gz
cd mummer-4.0.1
./configure --prefix="$PWD/../install"
make -j4 mummer
./mummer -version
```

Initial extraction in this container used ordinary `tar -xzf` and emitted ownership warnings; files were extracted and configure/build succeeded. The command above avoids those ownership warnings. Build logs and machine-readable provenance are retained under `reference/`.

## Contract used for comparison

| HelixAlign concept | Official option | Uniqueness |
|---|---|---|
| MUM | `-mum` | Exactly one occurrence across all reference records, exactly one occurrence within the current query record |
| MAM | `-mumreference` | Exactly one occurrence across all reference records; any count in current query record |
| MEM | `-maxmatch` | No uniqueness restriction |

All reported matches must be exact, left-maximal, and right-maximal for their specific pair of positions. Maximality is pairwise: a shorter MEM may coexist with a longer MEM at the same query start if they match different reference occurrences with different extension limits. Both maximal and maximum must not be conflated.

Official documentation explicitly describes reference-global and query-record-local uniqueness: https://mummer.sourceforge.net/manual/ (section 5.1, program options). The current source loops over query records and runs matching separately per record and strand: https://github.com/mummer4/mummer/blob/eb734606f2d516f42a0e0dce7a116bfb88ec1ebf/src/essaMEM/mummer.cpp .

Recommended invocation for nucleotide differential tests:

```bash
mummer -mum        -n -F -l 8 -b ref.fa qry.fa
mummer -mumreference -n -F -l 8 -b ref.fa qry.fa
mummer -maxmatch   -n -F -l 8 -b ref.fa qry.fa
```

- Use explicit match mode; documentation prose about defaults is inconsistent, whereas the executable usage says `-mumreference` is default.
- `-n` restricts matches to A/C/G/T, case-insensitive. Ambiguous characters are barriers, not wildcards or matching `N` strings.
- With no `-b`/`-r`, only query forward strand is used. `-b` adds reverse complement; `-r` selects only reverse complement.
- Reference coordinates are one-based and local to the named reference record.
- Default reverse coordinates are one-based positions in the reverse-complemented query record.
- Normalize a reverse match to zero-based original-query left coordinate by `q_original_start = query_length - (q_rc_1based - 1) - length`.
- Query uniqueness is tested within the current orientation, not by pooling occurrences across both strands. A unique palindrome can be reported once for each strand.
- Multi-reference sequences are separated by an internal sentinel. Do not concatenate biological records without barriers.
- `-F` forces reference ID in each match row; headers identify query and optionally `Reverse`.

## Empirically verified examples

All tests below used `-n -F -b`; positions shown are one-based. Complete stdout, stderr, commands and FASTA files are under `reference/semantics_cases/` and `reference/semantics_results.json`.

| Case, minimum length | Reference | Query | Verified forward result |
|---|---|---|---|
| Query repeat, 8 | `ACGTTGCA` | `ACGTTGCATACGTTGCA` | MUM empty; MAM/MEM `(1,1,8)` and `(1,10,8)` |
| Both unique, 8 | `TACGTTGCAG` | `AACGTTGCAC` | All modes `(2,2,8)` |
| Reference repeat, 8 | `ACGTTGCATACGTTGCA` | `ACGTTGCA` | MUM/MAM empty; MEM `(1,1,8)` and `(10,1,8)` |
| Repeated across reference records, 8 | r1=`ACGTTGCA`, r2=`ACGTTGCA` | `ACGTTGCA` | MUM/MAM empty; MEM one match on each reference record |
| Repeated across query records, 8 | `ACGTTGCA` | q1=`ACGTTGCA`, q2=`ACGTTGCA` | Every mode reports one match for each query record |
| Reference junction, 8 | r1=`AAAACCCC`, r2=`GGGGTTTT` | `CCCCGGGG` | All modes empty |
| Query junction, 8 | `CCCCGGGG` | q1=`AAAACCCC`, q2=`GGGGTTTT` | All modes empty |
| Overlapping query repeat, 3 | `AAA` | `AAAA` | MUM empty; MAM/MEM `(1,1,3)` and `(1,2,3)` |
| Palindrome, 4 | `ACGT` | `ACGT` | All modes report `(1,1,4)` on both strands |
| Ambiguity barrier, 4 | `ACGTNNNNTGCA` | `ACGTNNNNTGCA` | Forward `(1,1,4)` and `(9,9,4)`; no match spans Ns |

For reverse-only positive example, ref=`AACCGGTA`, query=`CCTACCGGTTG`, minimum 8: all modes report only reverse `(r=1, q_rc=2, length=8)`. Original-query zero-based left coordinate is `11 - 1 - 8 = 2`.

## Important defects in the reference executable discovered during setup

The official implementation is a comparison reference, not an infallible mathematical oracle. Two edge defects were directly reproduced in unmodified v4.0.1:

### 1. Avoid `-c`: reverse lengths can be incorrectly truncated

For the reverse example above, adding `-c` changes the row from `(1,2,8)` to `(1,10,2)` in **all three modes**. The expected length remains 8. Current source passes a converted reverse coordinate into `make_match()`, which then truncates length using `query_length - converted_coordinate`:

- https://github.com/mummer4/mummer/blob/eb734606f2d516f42a0e0dce7a116bfb88ec1ebf/include/mummer/sparseSA.hpp
- `make_match()` uses `min(length, query_length - query_position)`.
- `findMAM_each()` and `find_Lmaximal()` call it after `flip_forward` conversion.

This is a real observed output defect, not a claim that all reverse matching is invalid. Run without `-c` and normalize coordinates in the test harness.

### 2. Do not use official v4.0.1 as the oracle for minimum length 1

`mummer -n -F -l 1` on ref=`A`, query=`A` emits no matches for any of the three modes. On ref=`ACGT`, query=`ACGT`, `-mum`/`-mumreference` emit length 4 while `-maxmatch -l 1` emits no matches. Source has early skip conditions for traversal depth `<= 1`.

Use minimum length **at least 2** for official differential comparisons. Verify one-base edge semantics against the independent brute-force enumerator instead. Complete observed results are in `reference/single_base_results.json`.

These upstream defects were not patched in the reference build. They are documented so the differential test does not incorrectly reward reproducing an upstream defect.
