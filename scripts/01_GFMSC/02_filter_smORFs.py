#!/usr/bin/env python3
"""Filter Pyrodigal proteins to the final 5-50-aa smORF definition."""

from __future__ import annotations

import argparse
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
                yield header, "".join(sequence)
            header = line[1:]
            sequence = []
        else:
            sequence.append(line)
    if header is not None:
        yield header, "".join(sequence)


def write_record(handle: TextIO, record_id: str, sequence: str) -> None:
    handle.write(f">{record_id}\n")
    for start in range(0, len(sequence), 80):
        handle.write(sequence[start : start + 80] + "\n")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_dir", type=Path, help="Directory of Pyrodigal protein FASTA files")
    parser.add_argument("output_fasta", type=Path)
    parser.add_argument("--min-length", type=int, default=5)
    parser.add_argument("--max-length", type=int, default=50)
    parser.add_argument(
        "--allow-non-atg",
        action="store_true",
        help="Do not require start_type=ATG in the Pyrodigal header",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    fasta_files = sorted(
        path
        for pattern in ("*.fa", "*.faa", "*.fasta")
        for path in args.input_dir.glob(pattern)
    )
    if not fasta_files:
        raise SystemExit(f"No protein FASTA files found in {args.input_dir}")

    args.output_fasta.parent.mkdir(parents=True, exist_ok=True)
    retained = 0
    with args.output_fasta.open("w") as output_handle:
        for fasta_path in fasta_files:
            sample = fasta_path.name
            for suffix in (".proteins.faa", ".fasta", ".faa", ".fa"):
                if sample.endswith(suffix):
                    sample = sample[: -len(suffix)]
                    break

            with fasta_path.open() as input_handle:
                for header, sequence in read_fasta(input_handle):
                    if not args.allow_non_atg and "start_type=ATG" not in header:
                        continue

                    sequence = sequence.rstrip("*")
                    if sequence.startswith("M"):
                        sequence = sequence[1:]
                    if not args.min_length <= len(sequence) <= args.max_length:
                        continue

                    original_id = header.split()[0]
                    write_record(output_handle, f"{sample}|{original_id}", sequence)
                    retained += 1

    print(f"Retained {retained} smORFs in {args.output_fasta}")


if __name__ == "__main__":
    main()

