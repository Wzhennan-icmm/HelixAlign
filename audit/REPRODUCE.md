# Reproducing the 2026-10-27 correctness audit

This tests the exact-match layer, not NUCMER gapped alignment or GPU performance.
Run from this repository's root. Requirements: Rust with edition-2024 support
(the recorded run used 1.99.1), Python 3.10+, GCC/G++, make, curl, tar, and Git.
The Python oracle uses the standard library only. Cargo.lock is committed.

## 1. Rust regression suite and candidate

```bash
cargo test --locked
cargo build --release --locked
```

The recorded final suite has 39 tests: 28 library, 1 CLI parser unit test,
and 10 real CLI subprocess integration tests. Exhaustive loops within those
tests perform 90,642 algorithm comparisons and 8,636 suffix-lookup comparisons
against independent direct-enumeration oracles; these are not additional
independent biological samples.

## 2. Unmodified official reference

```bash
bash scripts/build_reference_mummer.sh /tmp/helixalign-mummer-reference
```

This downloads official v4.0.1, verifies its SHA256, and builds only mummer.
It refuses an existing source directory so that an edited reference is not
silently reused. The executable is a libtool wrapper; retain its build tree.

## 3. Three-way finite-corpus differential gate

```bash
python3 scripts/differential_mummer.py \
  --mummer /tmp/helixalign-mummer-reference/mummer-4.0.1/mummer \
  --candidate target/release/helixalign \
  --random-cases 120 --seed 20261027 \
  --report audit/results/rerun-differential.json
```

Expected: 1,296/1,296 passes for the official executable and separately
1,296/1,296 for the candidate, each compared with the independent Python oracle.
Each complete result is a multiset of `(ref_id, query_id, strand, ref_start0,
query_start0, length)`. Exactness, both-end maximality, uniqueness, record
boundaries, length thresholds, direction, and duplicate tuples are checked.

Both official/candidate mismatches and execution failures make the gate fail.
The original baseline's known mismatches are diagnostic and do not make it
fail when the optional --baseline argument is supplied. See the summary's
per-engine counters rather than interpreting exit 0 as a passing old baseline.

## 4. Reproduce the unmodified HelixAlign failures

Use a separate directory that does not yet exist:

```bash
git worktree add --detach /tmp/helixalign-original-20261027 \
  483511789035539b3cc265356af6b70e3b8f8cb9
cp Cargo.lock /tmp/helixalign-original-20261027/Cargo.lock
cargo test --locked --manifest-path /tmp/helixalign-original-20261027/Cargo.toml
cargo build --release --locked --manifest-path /tmp/helixalign-original-20261027/Cargo.toml
python3 scripts/differential_mummer.py \
  --mummer /tmp/helixalign-mummer-reference/mummer-4.0.1/mummer \
  --baseline /tmp/helixalign-original-20261027/target/release/helixalign \
  --candidate target/release/helixalign \
  --random-cases 120 --seed 20261027 \
  --report audit/results/rerun-with-original-differential.json
```

The old implementation's original 8 unit tests pass. Its forward-only
differential checks fail 298/432, with 4,885 extra and 665 missing tuples in
this corpus. Reverse tests are not attributed to the old ordinary CLI, which
does not support that mode. These are corpus totals, not genome-wide rates.

## 5. Hand-checkable minimal reproduction commands

```bash
target/release/helixalign -mum -l 4 --output-format matches \
  tests/fixtures/mummer_correctness/query_repeat.ref.fa \
  tests/fixtures/mummer_correctness/query_repeat.query.fa
# Only the TSV header: no MUM exists.

target/release/helixalign -mumreference -l 4 --output-format matches \
  tests/fixtures/mummer_correctness/query_repeat.ref.fa \
  tests/fixtures/mummer_correctness/query_repeat.query.fa
# r q + 0 0 4; r q + 0 4 4

target/release/helixalign -maxmatch -l 2 --output-format matches \
  tests/fixtures/mummer_correctness/nested_diagonals.ref.fa \
  tests/fixtures/mummer_correctness/nested_diagonals.query.fa
# r q + 0 0 4; r q + 0 2 2; r q + 2 0 2
```

The fixture README contains additional manual expectations. Full inputs,
commands and mismatch tuples are in the compressed machine reports:

```bash
gzip -dc audit/results/original-differential.json.gz > /tmp/original-differential.json
gzip -dc audit/results/fixed-differential.json.gz > /tmp/fixed-differential.json
```

Temporary input paths in recorded commands no longer exist. Their contents
are preserved in each report's cases list and in failing check records;
the runner recreates them automatically. Machine-specific executable paths
are provenance, not hard-coded requirements for rerunning elsewhere.

## Official version caveats

The differential corpus uses minimum lengths >=2 and does not pass -c.
Unmodified MUMmer 4.0.1's length-one and reverse -c edge behavior disagreed
with directly calculated results; see REFERENCE_SEMANTICS.md and the saved
reference semantic outputs. One-base HelixAlign behavior is instead covered
by the independent Rust exhaustive oracle. Do not claim equivalence for all
MUMmer switches, proteins, or every possible input from this finite suite.
