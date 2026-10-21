#!/usr/bin/env python3
"""Calculate the eight physicochemical properties used in the original workflow."""
import argparse
import csv
import re
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('output_csv', type=Path)
    parser.add_argument('inputs', nargs='+', metavar='GROUP=FASTA')
    args = parser.parse_args()
    from Bio import SeqIO
    from Bio.SeqUtils.ProtParam import ProteinAnalysis
    from pyteomics import electrochem
    import peptides

    rows = []
    cleaned_records = 0
    skipped_records = 0
    for spec in args.inputs:
        group, path = spec.split('=', 1)
        seen = set()
        count = 0
        for record in SeqIO.parse(path, 'fasta'):
            raw_seq = str(record.seq).upper()
            # Preserve the original workflow: remove non-standard residues before
            # calculating properties and skip a record only if nothing remains.
            seq = re.sub(r'[^ACDEFGHIKLMNPQRSTVWY]', '', raw_seq)
            if seq != raw_seq:
                cleaned_records += 1
            if not seq:
                skipped_records += 1
                continue
            if record.id in seen:
                raise ValueError(f'Duplicate ID in {group}: {record.id}')
            seen.add(record.id)
            pep = peptides.Peptide(seq)
            rows.append([
                record.id, group, electrochem.charge(seq, pH=7.0),
                ProteinAnalysis(seq).gravy(), pep.boman(), electrochem.pI(seq),
                pep.aliphatic_index(), pep.instability_index(),
                pep.hydrophobic_moment(), pep.energy_cost(),
            ])
            count += 1
        if not count:
            raise ValueError(f'Empty FASTA: {path}')
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.output_csv.open('w', newline='') as handle:
        writer = csv.writer(handle)
        writer.writerow(['ID', 'Group', 'Net Charge', 'Hydrophobicity', 'Boman Index',
                         'Isoelectric Point', 'Aliphatic Index', 'Instability Index',
                         'Amphiphilicity', 'Energy_cost'])
        writer.writerows(rows)
    print(f'Wrote eight properties for {len(rows)} records to {args.output_csv}')
    if cleaned_records:
        print(f'Removed non-standard residues from {cleaned_records} records')
    if skipped_records:
        print(f'Skipped {skipped_records} records with no standard residues')


if __name__ == '__main__':
    main()
