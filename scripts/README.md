# Analysis scripts

The numbered directories contain the maintained analysis workflows for the
Foodin manuscript.

## Modules

1. `01_GFMSC`: smORF prediction, filtering, clustering, taxonomy and novelty.
2. `02_bioactive_peptides`: bioactive-peptide prediction and sequence comparisons.
3. `03_foodin_characterization`: physicochemical, toxicity and taxonomic analyses.
4. `04_kefir_ecology`: Kefir abundance, association and functional analyses.
5. `05_in_vitro_validation`: antibacterial, cell-viability and haemolysis analyses.
6. `06_in_vivo_validation`: mouse infection, bacterial-burden and cytokine analyses.

## Usage

Each directory has a README describing its scripts, inputs and dependencies.
Scripts accept input and output paths rather than relying on a fixed project
location. Large input datasets, external databases and third-party software are
not bundled here.
