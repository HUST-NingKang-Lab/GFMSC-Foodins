#!/usr/bin/env python3
"""Merge final activity-specific FASTA files into one peptide table."""

from __future__ import annotations

import argparse
import csv
from collections import OrderedDict
from pathlib import Path
from typing import Iterator, TextIO, Tuple


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


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge LABEL=FASTA inputs into one row per peptide ID."
    )
    parser.add_argument("output_csv", type=Path)
    parser.add_argument(
        "activity_fastas",
        nargs="+",
        metavar="LABEL=FASTA",
        help="For example AMP=AMP_1178.fa ACP=ACP_1204.fa",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    peptides = OrderedDict()

    for specification in args.activity_fastas:
        if "=" not in specification:
            raise SystemExit(f"Invalid input specification: {specification}")
        activity, fasta_name = specification.split("=", 1)
        fasta_path = Path(fasta_name)

        with fasta_path.open() as handle:
            for peptide_id, sequence in read_fasta(handle):
                if peptide_id not in peptides:
                    peptides[peptide_id] = {
                        "sequence": sequence,
                        "activities": [],
                    }
                elif peptides[peptide_id]["sequence"] != sequence:
                    raise ValueError(
                        f"Conflicting sequences found for peptide ID {peptide_id}"
                    )

                if activity not in peptides[peptide_id]["activities"]:
                    peptides[peptide_id]["activities"].append(activity)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["peptide_id", "sequence", "predicted_activities"])
        for peptide_id, record in peptides.items():
            writer.writerow(
                [peptide_id, record["sequence"], ";".join(record["activities"])]
            )

    print(f"Wrote {len(peptides)} unique peptides to {args.output_csv}")


if __name__ == "__main__":
    main()

