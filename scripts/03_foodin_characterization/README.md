# Foodin characterization scripts

These scripts support the physicochemical, predicted-toxicity and taxonomic
analyses presented in Figures 2 and 3.

## Scripts

1. `01_calculate_physicochemical_properties.py`: calculate the eight peptide
   physicochemical properties produced by the original workflow; six are used
   in Figure 2c.
2. `02_plot_physicochemical_properties.R`: compare property distributions across
   Foodins and three reference AMP collections.
3. `03_calculate_amino_acid_composition.py`: calculate pooled amino-acid frequencies.
4. `04_plot_amino_acid_composition.R`: compare amino-acid composition across datasets.
5. `05_run_ToxinPred.sh`: predict peptide toxicity with ToxinPred 3.
6. `06_summarize_predicted_toxicity.R`: summarize predicted toxicity by dataset.
7. `07_summarize_taxonomic_distribution.R`: summarize assignment ranks and the
   ten genera with the highest Foodin counts.
8. `08_calculate_genus_density.py`: calculate Foodin occurrences per assembled Gbp.
9. `09_plot_genus_density.R`: compare genus-specific Foodin density across phyla.
10. `10_cluster_Foodin_occurrence.R`: cluster binary sample-presence profiles
    using Bray–Curtis dissimilarity and average linkage.

## Usage

Scripts accept input and output paths as command-line arguments. Python scripts
provide `--help`; R and shell scripts display usage instructions when run without
arguments. For physicochemical calculations, non-standard residue symbols are
removed before calculation, following the original workflow. For pooled
amino-acid composition, all residue symbols remain in the denominator while the
20 standard amino-acid frequencies are reported.

The density workflow requires occurrence counts, non-overlapping assembled-base
totals and a genus-to-phylum mapping. Occurrence clustering starts from a prepared
binary peptide-by-sample matrix. The associated occurrence tree and iTOL
annotation files are provided under `data/source_data/figure3/itol`.

## Requirements

Python: Biopython, pyteomics, peptides and pandas. R: ggplot2, readr, tidyr,
dplyr, ggridges, scales, vegan and ape. External software: ToxinPred 3.
