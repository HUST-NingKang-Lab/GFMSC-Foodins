#!/usr/bin/env python3
"""Calculate pooled residue frequencies using the original workflow."""
import argparse
import csv
from collections import Counter
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output_csv', type=Path)
    parser.add_argument('inputs', nargs='+', metavar='GROUP=FASTA')
    args = parser.parse_args()
    from Bio import SeqIO

    alphabet = 'ACDEFGHIKLMNPQRSTVWY'
    rows = []
    for spec in args.inputs:
        group, path = spec.split('=', 1)
        counts = Counter()
        for record in SeqIO.parse(path, 'fasta'):
            seq = str(record.seq).upper()
            if not seq:
                raise ValueError(f'{path}: empty sequence {record.id}')
            # The original analysis counted every residue in the denominator and
            # reported the frequencies of the 20 standard amino acids.
            counts.update(seq)
        total = sum(counts.values())
        if total == 0:
            raise ValueError(f'Empty FASTA: {path}')
        rows.append([group] + [counts[aa] / total for aa in alphabet])
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['Group'] + list(alphabet))
        writer.writerows(rows)


if __name__ == '__main__':
    main()
