# Data

This directory contains the peptide catalogs, sample metadata and source data
distributed with the Foodin analysis code. The existing directory structure is
retained so that peptide resources, metadata and figure-level data can be found
independently of the numbered script modules.

## Directory contents

### `metadata`

`sample_metadata.csv` contains 10,434 records. Of these, 10,112 records describe
MAGs and 322 describe fermented-grain metagenomes for which `MAG_id` is recorded
as `/` because a MAG identifier is not applicable. The remaining fields record
the data source, sample identifier, food category, food type and fermentation
status.

### `peptide_catalogs`

| Directory or file | Records | Description |
| --- | ---: | --- |
| `Foodins/AMP_1178.fa` | 1,178 | Predicted antimicrobial peptides termed Foodins |
| `other_bioactive_peptides/ACP_1204.fa` | 1,204 | Predicted anticancer peptides |
| `other_bioactive_peptides/ADP_80.fa` | 80 | Predicted anti-diabetic peptides |
| `other_bioactive_peptides/AHP_284.fa` | 284 | Predicted anti-hypertensive peptides |
| `other_bioactive_peptides/AIP_288.fa` | 288 | Predicted anti-inflammatory peptides |
| `other_bioactive_peptides/AOP_227.fa` | 227 | Predicted antioxidant peptides |

The six activity FASTA files are multi-label collections and are not mutually
exclusive. Their 3,261 activity assignments correspond to 2,905 unique peptide
sequences. Eighteen Foodins also occur in the ACP collection and four Foodins
also occur in the AIP collection. These overlaps are intentional and should not
be removed when reproducing the activity-class results.

`reference_peptides` contains the human microbiome-derived AMP collection,
DBAASP and the reference sets used for the five other predicted activities.
Some supplied reference records contain ambiguous or non-standard residue
symbols. File names ending in `_uniq` preserve the names of the analysis inputs
and do not imply that every sequence string is unique across records.

### `source_data/figure2`

- `AA_frequence.csv` contains pooled amino-acid frequencies for Foodins and the
  three AMP reference datasets. The legacy filename is retained for consistency
  with the original analysis.
- `Toxinpred3` contains the ToxinPred 3 output for Foodins, DBAASP and the human
  microbiome-derived AMP collection. The AMPSphere result is provided as
  `AMPSphere.csv.gz` and can be read directly by software with gzip support or
  decompressed before use.

### `source_data/figure3`

- `AMP_taxrank_3341.csv` and `AMP_genus.csv` contain taxonomic assignments.
- `AMP_density_genus.tsv` contains genus-level Foodin counts and assembled-base
  normalized density values.
- `itol` contains the occurrence tree and the iTOL annotation files used for
  Figure 3c.

### `source_data/figure4`

- `kefir_abundance.csv` contains relative abundances for 3,622 SGBs across 284
  Kefir samples. Abundance values are percentage points on a 0–100 scale;
  prevalence values in the derived summary are fractions on a 0–1 scale.
- `mag_sgb.csv` maps 10,899 MAGs to SGBs.
- `sgb_summary_with_amp.csv` contains SGB taxonomy, abundance, prevalence,
  Foodin-positive MAG statistics and the 284 sample abundance columns.
- `COG_category_count_matrix.csv` and `CAZy_count_matrix.csv` contain functional
  counts for the same 10,112 MAGs.
- `itol` contains the phylogenetic tree and annotation files used for Figure 4a.

The COG and CAZy matrices have identical MAG coverage. The additional 787 MAGs
in `mag_sgb.csv` do not have corresponding COG or CAZy count rows.

### `source_data/figure5`

The three CSV files contain the processed MIC calls, cell-viability summaries
and haemolysis summaries used for Figure 5. Peptide concentrations are in
μmol L⁻¹. These files contain summarized values rather than the original
replicate-level plate readings.

### `source_data/figure6`

`Figure6_CFU_mouse_level_source_data.csv` contains bacterial burden for six mice
per group, reported as CFU g⁻¹ and log10 CFU g⁻¹. The ELISA source-data file
contains TNF-α, IL-6 and IL-10 concentrations for three mice per group in
pg ml⁻¹.
