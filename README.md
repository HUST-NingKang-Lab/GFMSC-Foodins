# GFMSC-Foodins analysis code

This repository contains the maintained analysis scripts supporting the
manuscript **“Mining global food microbiomes unveils a rich reservoir of
bioactive peptides.”**

## Repository structure

- `data`: peptide catalogs, metadata and figure-level source data. See
  [data/README.md](data/README.md) for file descriptions.
- `scripts/01_GFMSC`: construction and annotation of the Global Food Microbial
  Small Open Reading Frame Catalog.
- `scripts/02_bioactive_peptides`: bioactive-peptide prediction and sequence
  comparisons.
- `scripts/03_foodin_characterization`: physicochemical, predicted-toxicity and
  taxonomic analyses.
- `scripts/04_kefir_ecology`: Kefir abundance, association and functional
  analyses.
- `scripts/05_in_vitro_validation`: antibacterial, cell-viability and
  haemolysis analyses.
- `scripts/06_in_vivo_validation`: mouse infection, bacterial-burden and
  cytokine analyses.

## Usage

See [scripts/README.md](scripts/README.md) and the README in each numbered
directory for the required inputs, dependencies and command-line interfaces.
Scripts accept explicit input and output paths and should be run in numerical
order within each module when steps depend on earlier outputs.

The repository includes the principal peptide catalogs, figure-level source
data and maintained analysis workflows used in the manuscript. External tools
and databases should be obtained from their respective providers as described
in the module READMEs.

## Contact

Kouyi Zhou, zhoukouyi@hust.edu.cn
