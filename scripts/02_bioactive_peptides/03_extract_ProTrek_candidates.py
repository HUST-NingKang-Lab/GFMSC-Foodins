#!/usr/bin/env python3
"""Export unique ProTrek-positive GFMSC candidates as an activity FASTA."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Iterator, TextIO, Tuple


def read_fasta(handle: TextIO) -> Iterator[Tuple[str, str]]:
    header = None
    sequence = []
    for raw_line in handle:
        line = raw_line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if header is not None:
                yield header.split()[0], "".join(sequence)
            header = line[1:]
            sequence = []
        else:
            sequence.append(line)
    if header is not None:
        yield header.split()[0], "".join(sequence)


def write_fasta_record(handle: TextIO, sequence_id: str, sequence: str) -> None:
    handle.write(f">{sequence_id}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("protrek_csv", type=Path)
    parser.add_argument("gfmsc_fasta", type=Path)
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument(
        "--summary-csv",
        type=Path,
        help="Optional table containing the best reference match per candidate",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    best_matches: Dict[str, Tuple[str, float]] = {}

    with args.protrek_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"query_id", "reference_id", "score"}
        if not required.issubset(reader.fieldnames or []):
            raise SystemExit(
                f"ProTrek CSV must contain columns: {', '.join(sorted(required))}"
            )

        for row in reader:
            query_id = row["query_id"]
            score = float(row["score"])
            if query_id not in best_matches or score > best_matches[query_id][1]:
                best_matches[query_id] = (row["reference_id"], score)

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    found = set()
    with args.gfmsc_fasta.open() as input_handle, args.output_fasta.open("w") as output_handle:
        for sequence_id, sequence in read_fasta(input_handle):
            if sequence_id in best_matches:
                write_fasta_record(output_handle, sequence_id, sequence)
                found.add(sequence_id)

    missing = sorted(set(best_matches) - found)
    if missing:
        raise SystemExit(
            f"{len(missing)} ProTrek-positive IDs were absent from {args.gfmsc_fasta}; "
            f"first missing ID: {missing[0]}"
        )

    if args.summary_csv:
        args.summary_csv.parent.mkdir(parents=True, exist_ok=True)
        with args.summary_csv.open("w", newline="") as handle:
            writer = csv.writer(handle)
            writer.writerow(["query_id", "best_reference_id", "best_score"])
            for query_id in sorted(best_matches):
                reference_id, score = best_matches[query_id]
                writer.writerow([query_id, reference_id, score])

    print(f"Wrote {len(found)} unique candidates to {args.output_fasta}")


if __name__ == "__main__":
    main()

