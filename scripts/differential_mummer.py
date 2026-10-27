#!/usr/bin/env python3
"""Compare exact-match tuples with a brute-force oracle and official MUMmer.

This oracle deliberately does not use a suffix array, seed index, or containment
filter. It enumerates every pair of starts, extends the exact match, and checks
the left flank. Uniqueness counts every occurrence, including overlapping ones.
The tested contract is case-insensitive ACGT DNA (official mummer -n), reference
uniqueness across all reference records, query uniqueness within its own record,
and independently processed strands. All normalized coordinates are zero-based
starts on the original forward sequences, even for reverse-strand matches.

Official MUMmer 4.0.1 has known -l 1 and -c issues; this script uses -l >= 2 and
converts reverse coordinates itself without -c. These restrictions do not alter
the oracle, which supports min_len=1 for independent unit tests.
"""

from __future__ import annotations

import argparse
from collections import Counter
from dataclasses import dataclass
import hashlib
import json
from pathlib import Path
import random
import re
import subprocess
import sys
import tempfile


BASES = frozenset("ACGT")
MODES = {"MUM": "-mum", "MAM": "-mumreference", "MEM": "-maxmatch"}
STRANDS = {"forward": (), "reverse": ("-r",), "both": ("-b",)}
COMPLEMENT = str.maketrans("ACGTRYMKSWBDHVN", "TGCAYRKMSWVHDBN")
Match = tuple[str, str, str, int, int, int]
Records = list[tuple[str, str]]


@dataclass
class Case:
    name: str
    reference: Records
    query: Records
    min_len: int
    description: str = ""


def revcomp(sequence: str) -> str:
    return sequence.upper().translate(COMPLEMENT)[::-1]


def read_fasta(path: Path) -> Records:
    records: Records = []
    name = None
    chunks: list[str] = []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records.append((name, "".join(chunks).upper()))
            name, chunks = line[1:].split()[0], []
        elif line.strip():
            if name is None:
                raise ValueError(f"Sequence before FASTA header: {path}")
            chunks.append("".join(line.split()))
    if name is not None:
        records.append((name, "".join(chunks).upper()))
    if len({name for name, _ in records}) != len(records):
        raise ValueError(f"Duplicate FASTA identifiers: {path}")
    return records


def write_fasta(path: Path, records: Records) -> None:
    path.write_text("".join(f">{name}\n{sequence}\n" for name, sequence in records))


def occurrences(pattern: str, records: Records) -> int:
    return sum(
        sequence[i : i + len(pattern)] == pattern
        for _, raw_sequence in records
        for sequence in (raw_sequence.upper(),)
        for i in range(len(sequence) - len(pattern) + 1)
    )


def brute_force(case: Case, mode: str, strand_setting: str) -> list[Match]:
    if case.min_len < 1:
        raise ValueError("Minimum match length must be positive")
    directions = ("+", "-") if strand_setting == "both" else (
        "+" if strand_setting == "forward" else "-",
    )
    result: list[Match] = []
    refs = [(name, seq.upper()) for name, seq in case.reference]
    for qname, original_query in case.query:
        original_query = original_query.upper()
        for strand in directions:
            query = original_query if strand == "+" else revcomp(original_query)
            for rname, reference in refs:
                for r in range(len(reference)):
                    for q in range(len(query)):
                        if reference[r] not in BASES or reference[r] != query[q]:
                            continue
                        # A match with an equal preceding base is not left-maximal.
                        if r and q and reference[r - 1] in BASES and reference[r - 1] == query[q - 1]:
                            continue
                        length = 0
                        while (r + length < len(reference) and q + length < len(query)
                               and reference[r + length] in BASES
                               and reference[r + length] == query[q + length]):
                            length += 1
                        if length < case.min_len:
                            continue
                        pattern = reference[r : r + length]
                        if mode in ("MUM", "MAM") and occurrences(pattern, refs) != 1:
                            continue
                        if mode == "MUM" and occurrences(pattern, [(qname, query)]) != 1:
                            continue
                        original_q = q if strand == "+" else len(query) - q - length
                        result.append((rname, qname, strand, r, original_q, length))
    return sorted(result)


def validate_matches(matches: list[Match], case: Case, mode: str) -> list[dict]:
    refs = {name: seq.upper() for name, seq in case.reference}
    queries = {name: seq.upper() for name, seq in case.query}
    invalid = []
    for match in matches:
        rname, qname, strand, r, q0, length = match
        reasons = []
        if rname not in refs or qname not in queries:
            reasons.append("unknown_record_or_out_of_range_start")
        elif strand not in ("+", "-"):
            reasons.append("invalid_strand")
        else:
            reference = refs[rname]
            original_query = queries[qname]
            query = original_query if strand == "+" else revcomp(original_query)
            q = q0 if strand == "+" else len(query) - q0 - length
            if length < case.min_len:
                reasons.append("below_minimum_length")
            if min(r, q0, q) < 0 or r + length > len(reference) or q0 + length > len(query):
                reasons.append("crosses_record_boundary_or_out_of_bounds")
            else:
                pattern = reference[r : r + length]
                if not set(pattern) <= BASES:
                    reasons.append("contains_non_ACGT")
                if pattern != query[q : q + length]:
                    reasons.append("not_exact")
                if r and q and reference[r - 1] in BASES and reference[r - 1] == query[q - 1]:
                    reasons.append("not_left_maximal")
                if (r + length < len(reference) and q + length < len(query)
                        and reference[r + length] in BASES
                        and reference[r + length] == query[q + length]):
                    reasons.append("not_right_maximal")
                if mode in ("MUM", "MAM") and occurrences(pattern, case.reference) != 1:
                    reasons.append("not_unique_in_reference_collection")
                if mode == "MUM" and occurrences(pattern, [(qname, query)]) != 1:
                    reasons.append("not_unique_in_query_record")
        if reasons:
            invalid.append({"match": match, "reasons": reasons})
    for match, count in Counter(matches).items():
        if count > 1:
            invalid.append({"match": match, "reasons": ["duplicate_tuple"], "count": count})
    return invalid


def parse_mummer(output: str, case: Case) -> list[Match]:
    qlengths = {name: len(seq) for name, seq in case.query}
    qname = None
    strand = "+"
    matches = []
    for line in output.splitlines():
        if not line.strip():
            continue
        fields = line.split()
        if line.startswith(">"):
            qname = fields[1]
            strand = "-" if "Reverse" in fields[2:] else "+"
            continue
        if qname is None or len(fields) != 4:
            raise ValueError(f"Unexpected official -F output: {line!r}")
        rname, r, q, length = fields
        r, q, length = int(r) - 1, int(q) - 1, int(length)
        if strand == "-":
            q = qlengths[qname] - q - length
        matches.append((rname, qname, strand, r, q, length))
    return sorted(matches)


def parse_candidate(output: str, _case: Case) -> list[Match]:
    matches = []
    for line in output.splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        fields = line.split()
        if fields[0] == "ref_id":
            continue
        if len(fields) != 6:
            raise ValueError(f"Unexpected candidate matches output: {line!r}")
        matches.append((fields[0], fields[1], fields[2], *(int(x) for x in fields[3:])))
    return sorted(matches)


def flattened_coordinate(position: int, records: Records) -> tuple[str, int]:
    offset = 0
    for name, sequence in records:
        if offset <= position < offset + len(sequence):
            return name, position - offset
        offset += len(sequence)
    return "<out-of-range>", position


def parse_legacy(output: str, case: Case) -> list[Match]:
    matches = []
    for line in output.splitlines():
        if not line.strip() or line.startswith(">"):
            continue
        parsed = re.fullmatch(r"\s*Ref:\s+(\d+)\s+Query:\s+(\d+)\s+Len:\s+(\d+)\s*", line)
        if not parsed:
            raise ValueError(f"Unexpected original output: {line!r}")
        r, q, length = (int(x) for x in parsed.groups())
        rname, r = flattened_coordinate(r - 1, case.reference)
        qname, q = flattened_coordinate(q - 1, case.query)
        matches.append((rname, qname, "+", r, q, length))
    return sorted(matches)


def run_engine(command: list[str], parser, case: Case, timeout: float) -> tuple[list[Match], dict]:
    completed = subprocess.run(command, text=True, capture_output=True, timeout=timeout)
    details = {"command": command, "exit_code": completed.returncode}
    if completed.returncode:
        raise RuntimeError(f"exit={completed.returncode}: {completed.stderr[-2000:]}")
    matches = parser(completed.stdout, case)
    return matches, details


def compare(actual: list[Match], expected: list[Match], case: Case, mode: str) -> dict:
    ac, ex = Counter(actual), Counter(expected)
    extra = sorted((ac - ex).elements())
    missing = sorted((ex - ac).elements())
    invalid = validate_matches(actual, case, mode)
    return {"passed": not (extra or missing or invalid), "expected_count": len(expected),
            "actual_count": len(actual), "false_positives": extra, "false_negatives": missing,
            "invalid_matches": invalid}


def fixture_cases(directory: Path) -> list[Case]:
    return [Case(item["name"], read_fasta(directory / item["reference"]),
                 read_fasta(directory / item["query"]), item["min_len"], item["description"])
            for item in json.loads((directory / "cases.json").read_text())]


def random_cases(count: int, seed: int) -> list[Case]:
    rng = random.Random(seed)
    cases = []
    def sequence() -> str:
        length = rng.randint(1, 22)
        if rng.random() < 0.45:
            unit = "".join(rng.choices("ACGT", k=rng.randint(1, 4)))
            text = (unit * (length + 1))[:length]
        else:
            text = "".join(rng.choices("ACGTN", weights=[5, 5, 5, 5, 1], k=length))
        return text.lower() if rng.random() < 0.15 else text
    for i in range(count):
        refs = [(f"r{j}", sequence()) for j in range(rng.randint(1, 3))]
        queries = []
        for j in range(rng.randint(1, 3)):
            text = sequence()
            if rng.random() < 0.6:
                source = rng.choice(refs)[1]
                text = source if rng.random() < 0.5 else revcomp(source)
                if rng.random() < 0.35:
                    text = "A" + text + "C"
            queries.append((f"q{j}", text))
        if len(refs) > 1 and rng.random() < 0.2:
            refs[1] = (refs[1][0], refs[0][1])
        cases.append(Case(f"random_{i:04d}", refs, queries, rng.randint(2, 6), f"seed={seed}"))
    return cases


def self_check() -> None:
    # Hand-calculated assertions guard the oracle independently of either binary.
    repeat = Case("repeat", [("r", "ACGT")], [("q", "ACGTACGT")], 4)
    assert brute_force(repeat, "MUM", "forward") == []
    expected = [("r", "q", "+", 0, 0, 4), ("r", "q", "+", 0, 4, 4)]
    assert brute_force(repeat, "MAM", "forward") == expected
    assert brute_force(repeat, "MEM", "forward") == expected
    nested = Case("nested", [("r", "ACAC")], [("q", "ACAC")], 2)
    assert brute_force(nested, "MEM", "forward") == [
        ("r", "q", "+", 0, 0, 4), ("r", "q", "+", 0, 2, 2), ("r", "q", "+", 2, 0, 2)]
    homo = Case("homo", [("r", "AAAA")], [("q", "AAAA")], 2)
    assert len(brute_force(homo, "MEM", "forward")) == 5
    assert brute_force(homo, "MUM", "forward") == [("r", "q", "+", 0, 0, 4)]
    single = Case("single", [("r", "A")], [("q", "A")], 1)
    assert brute_force(single, "MUM", "forward") == [("r", "q", "+", 0, 0, 1)]
    assert revcomp("ARYKMSWBDHVN") == "NBDHVWSKMRYT"


def file_sha256(path: str) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def binary_identity(binary: str) -> dict:
    path = Path(binary).resolve()
    identity = {"path": str(path), "sha256": file_sha256(str(path))}
    # An uninstalled official build uses a libtool shell launcher; retain both
    # the entry point and the actual executable hashes for reproducibility.
    compiled = path.parent / ".libs" / path.name
    if compiled.is_file():
        identity["libtool_executable"] = {"path": str(compiled), "sha256": file_sha256(str(compiled))}
    return identity


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mummer", required=True, help="Official mummer executable; use v4.0.1")
    parser.add_argument("--candidate", help="Patched HelixAlign binary with --output-format matches")
    parser.add_argument("--baseline", help="Original HelixAlign binary (forward legacy output only)")
    parser.add_argument("--fixtures", type=Path, default=Path(__file__).resolve().parents[1] / "tests/fixtures/mummer_correctness")
    parser.add_argument("--random-cases", type=int, default=120)
    parser.add_argument("--seed", type=int, default=20261027)
    parser.add_argument("--report", type=Path, help="Write complete structured JSON report here")
    parser.add_argument("--timeout", type=float, default=15)
    args = parser.parse_args()
    self_check()
    fixtures = fixture_cases(args.fixtures)
    cases = fixtures + random_cases(args.random_cases, args.seed)
    report = {"schema_version": 1, "contract": "ACGT only; per-record query; global reference uniqueness; 0-based forward coordinates",
              "seed": args.seed, "fixture_cases": len(fixtures), "random_cases": args.random_cases,
              "official_restrictions": "MUMmer 4.0.1: min_len>=2; do not use -c; normalize reverse coordinates in Python",
              "baseline_restrictions": "Original default CLI only supports forward matching; reverse cases not executed",
              "binaries": {}, "checks": [], "summary": {},
              "cases": [{"name": c.name, "description": c.description, "reference": c.reference,
                         "query": c.query, "min_len": c.min_len} for c in cases]}
    for engine in ("mummer", "candidate", "baseline"):
        binary = getattr(args, engine)
        if binary:
            report["binaries"][engine] = binary_identity(binary)
    for engine in ("official", "candidate", "baseline"):
        if engine == "official" or getattr(args, engine):
            report["summary"][engine] = {"checks": 0, "passed": 0, "failed": 0, "errors": 0,
                                         "false_positives": 0, "false_negatives": 0, "invalid_matches": 0}
    with tempfile.TemporaryDirectory(prefix="helixalign-diff-") as tmp:
        ref_path, query_path = Path(tmp) / "ref.fa", Path(tmp) / "query.fa"
        for case in cases:
            write_fasta(ref_path, case.reference)
            write_fasta(query_path, case.query)
            for mode, flag in MODES.items():
                for strands, strand_flags in STRANDS.items():
                    expected = brute_force(case, mode, strands)
                    engines = [("official", [args.mummer, "-F", "-n", flag, "-l", str(case.min_len), *strand_flags], parse_mummer)]
                    if args.candidate:
                        engines.append(("candidate", [args.candidate, flag, "-l", str(case.min_len), "--output-format", "matches", *strand_flags], parse_candidate))
                    if args.baseline and strands == "forward":
                        engines.append(("baseline", [args.baseline, flag, "-l", str(case.min_len)], parse_legacy))
                    for engine, command, output_parser in engines:
                        check = {"case": case.name, "mode": mode, "strands": strands, "engine": engine}
                        summary = report["summary"][engine]
                        summary["checks"] += 1
                        try:
                            actual, details = run_engine(command + [str(ref_path), str(query_path)], output_parser, case, args.timeout)
                            check.update(compare(actual, expected, case, mode))
                            # Fixture tuples remain inspectable even when the check passes.
                            if case in fixtures or not check["passed"]:
                                check.update(expected=expected, actual=actual, invocation=details)
                        except (ValueError, RuntimeError, subprocess.TimeoutExpired) as error:
                            check.update(passed=False, error=str(error), expected=expected, invocation=command)
                            summary["errors"] += 1
                        summary["passed" if check["passed"] else "failed"] += 1
                        for key in ("false_positives", "false_negatives", "invalid_matches"):
                            summary[key] += len(check.get(key, []))
                        if not check["passed"]:
                            check.update(reference=case.reference, query=case.query, min_len=case.min_len)
                        report["checks"].append(check)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({key: report[key] for key in ("seed", "fixture_cases", "random_cases", "summary")}, indent=2))
    # Baseline failures document the regression; official/candidate failures fail the gate.
    return int(any(summary["failed"] for engine, summary in report["summary"].items() if engine != "baseline"))


if __name__ == "__main__":
    sys.exit(main())
