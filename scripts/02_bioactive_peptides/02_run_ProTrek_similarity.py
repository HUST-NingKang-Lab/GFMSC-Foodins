#!/usr/bin/env python3
"""Identify candidate bioactive peptides by ProTrek sequence similarity."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path
from typing import Iterator, Sequence, TextIO, Tuple


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
        description=(
            "Encode GFMSC and reference peptides with ProTrek-650M and retain "
            "sequence pairs above the specified similarity-score threshold."
        )
    )
    parser.add_argument("query_fasta", type=Path, help="GFMSC representative peptides")
    parser.add_argument("reference_fasta", type=Path, help="Reference peptides for one activity")
    parser.add_argument("output_csv", type=Path)
    parser.add_argument(
        "--protrek-root",
        required=True,
        type=Path,
        help="Root directory of the cloned ProTrek repository",
    )
    parser.add_argument(
        "--weights-dir",
        default="weights/ProTrek_650M",
        help="Path relative to --protrek-root, or an absolute path",
    )
    parser.add_argument("--score-threshold", type=float, default=48.0)
    parser.add_argument("--batch-size", type=int, default=64)
    parser.add_argument(
        "--strict-greater",
        action="store_true",
        help="Use score > threshold instead of the recorded score >= threshold",
    )
    return parser.parse_args()


def batches(values: Sequence[str], size: int) -> Iterator[Sequence[str]]:
    for start in range(0, len(values), size):
        yield values[start : start + size]


def main() -> None:
    args = parse_args()
    protrek_root = args.protrek_root.resolve()
    if not protrek_root.is_dir():
        raise SystemExit(f"ProTrek directory not found: {protrek_root}")

    sys.path.insert(0, str(protrek_root))

    try:
        import torch
        from tqdm import tqdm
        from model.ProTrek.protrek_trimodal_model import ProTrekTrimodalModel
    except ImportError as exc:
        raise SystemExit(
            "ProTrek or one of its Python dependencies could not be imported. "
            "Activate the ProTrek environment and verify --protrek-root."
        ) from exc

    weights_dir = Path(args.weights_dir)
    if not weights_dir.is_absolute():
        weights_dir = protrek_root / weights_dir

    model_config = {
        "protein_config": str(weights_dir / "esm2_t33_650M_UR50D"),
        "text_config": str(
            weights_dir / "BiomedNLP-PubMedBERT-base-uncased-abstract-fulltext"
        ),
        "structure_config": str(weights_dir / "foldseek_t30_150M"),
        "from_checkpoint": str(weights_dir / "ProTrek_650M.pt"),
    }

    with args.query_fasta.open() as handle:
        query_records = list(read_fasta(handle))
    with args.reference_fasta.open() as handle:
        reference_records = list(read_fasta(handle))

    if not query_records or not reference_records:
        raise SystemExit("The query and reference FASTA files must both contain sequences.")

    query_ids, query_sequences = zip(*query_records)
    reference_ids, reference_sequences = zip(*reference_records)

    device = "cuda" if torch.cuda.is_available() else "cpu"
    model = ProTrekTrimodalModel(**model_config).eval().to(device)

    reference_embeddings = []
    with torch.no_grad():
        for reference_batch in tqdm(
            list(batches(reference_sequences, args.batch_size)),
            desc="Encoding reference peptides",
        ):
            reference_embeddings.append(model.get_protein_repr(list(reference_batch)))
    reference_embeddings = torch.cat(reference_embeddings, dim=0)

    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    retained_pairs = 0
    with args.output_csv.open("w", newline="") as output_handle:
        writer = csv.writer(output_handle)
        writer.writerow(["query_id", "reference_id", "score"])

        with torch.no_grad():
            for start in tqdm(
                range(0, len(query_sequences), args.batch_size),
                desc="Comparing GFMSC peptides",
            ):
                query_batch = query_sequences[start : start + args.batch_size]
                query_embeddings = model.get_protein_repr(list(query_batch))
                scores = query_embeddings @ reference_embeddings.T / model.temperature

                if args.strict_greater:
                    matched_rows, matched_columns = torch.where(
                        scores > args.score_threshold
                    )
                else:
                    matched_rows, matched_columns = torch.where(
                        scores >= args.score_threshold
                    )

                matches = []
                for row_tensor, column_tensor in zip(matched_rows, matched_columns):
                    row = int(row_tensor.item())
                    column = int(column_tensor.item())
                    matches.append(
                        (
                            query_ids[start + row],
                            reference_ids[column],
                            float(scores[row, column].item()),
                        )
                    )

                matches.sort(key=lambda item: item[2], reverse=True)
                writer.writerows(matches)
                retained_pairs += len(matches)

    operator = ">" if args.strict_greater else ">="
    print(
        f"Wrote {retained_pairs} query-reference pairs with scores "
        f"{operator} {args.score_threshold} to {args.output_csv}"
    )


if __name__ == "__main__":
    main()

