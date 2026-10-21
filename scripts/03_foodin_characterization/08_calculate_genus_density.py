#!/usr/bin/env python3
"""Foodin gene occurrences per assembled Gbp, using the supplied taxon tables."""
import argparse
from pathlib import Path


def main():
    p = argparse.ArgumentParser()
    p.add_argument('occurrences_csv', help='Columns: level, source; one row per gene occurrence')
    p.add_argument('assembled_bp_tsv', help='Columns: level, name, nbps; non-overlapping bp totals')
    p.add_argument('genus_phylum_tsv', help='Columns: genus, phylum; one phylum per genus')
    p.add_argument('output_tsv', type=Path)
    args = p.parse_args()
    import pandas as pd

    occurrences = pd.read_csv(args.occurrences_csv)
    bp = pd.read_csv(args.assembled_bp_tsv, sep='\t')
    taxonomy = pd.read_csv(args.genus_phylum_tsv, sep='\t', dtype=str)
    for table, cols in [(occurrences, ['level', 'source']), (bp, ['level', 'name', 'nbps']),
                        (taxonomy, ['genus', 'phylum'])]:
        if not set(cols).issubset(table.columns):
            raise ValueError(f'Required columns: {cols}')
    occurrences = occurrences.loc[occurrences.level.isin(['genus', 'species'])].copy()
    bp = bp.loc[bp.level.isin(['genus', 'species'])].copy()
    if occurrences.source.isna().any() or bp.name.isna().any():
        raise ValueError('Missing genus/species names')
    # Preserve the original first-token genus extraction for these two input tables.
    occurrences['genus'] = occurrences.source.str.split().str[0]
    bp['genus'] = bp['name'].str.split().str[0]
    bp['nbps'] = pd.to_numeric(bp['nbps'], errors='raise')
    if bp['nbps'].isna().any() or (bp['nbps'] <= 0).any():
        raise ValueError('Assembled base-pair counts must be positive')
    counts = occurrences.groupby('genus').size().rename('amp_genes')
    totals = bp.groupby('genus')['nbps'].sum()
    if not counts.index.difference(totals.index).empty:
        raise ValueError('Some Foodin-bearing genera have no assembled-bp denominator')
    result = totals.to_frame().join(counts).fillna({'amp_genes': 0})
    result['amp_genes'] = result.amp_genes.astype(int)
    result['amps_per_Gbp'] = result.amp_genes * 1e9 / result.nbps
    result = result.loc[result.amp_genes > 0].reset_index()
    taxonomy['genus'] = taxonomy.genus.str.replace(r'^g__', '', regex=True)
    taxonomy = taxonomy[['genus', 'phylum']].drop_duplicates()
    if taxonomy.genus.duplicated().any():
        raise ValueError('Ambiguous genus-to-phylum mapping')
    result = result.merge(taxonomy, on='genus', how='left', validate='one_to_one')
    if result.phylum.isna().any():
        raise ValueError('Phylum annotation missing for a Foodin-bearing genus')
    args.output_tsv.parent.mkdir(parents=True, exist_ok=True)
    result[['genus', 'phylum', 'amp_genes', 'nbps', 'amps_per_Gbp']].to_csv(
        args.output_tsv, sep='\t', index=False)


if __name__ == '__main__':
    main()
