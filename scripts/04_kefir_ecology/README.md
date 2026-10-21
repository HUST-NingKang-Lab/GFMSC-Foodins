# Kefir ecology analysis scripts

These scripts support the Kefir microbial ecology and functional analyses
presented in Figure 4.

## Scripts

1. `01_prepare_kefir_metadata.R`: calculate SGB abundance, prevalence and
   Foodin-positive MAG fractions.
2. `02_plot_foodin_abundance.R`: plot abundance against the fraction of
   Foodin-positive MAGs.
3. `03_build_association_network.R`: perform CLR transformation, Spearman tests,
   FDR filtering and positive-edge module detection.
4. `04_plot_COG_profiles.R`: compare COG profiles of selected Kefir species.
5. `05_plot_CAZy_profiles.R`: compare CAZy broad-class and family profiles.
6. `06_calculate_ODS_FDI.R`: calculate ODS and FDI for positively associated SGBs.

`common.R` provides shared input validation and table-reading functions used by
the numbered scripts; it is sourced automatically and is not run separately.

## Usage

Scripts accept input and output paths as command-line arguments. Run step 1
first; its `sgb_summary_with_amp.csv` and `sample_ids.txt` outputs are used by
later steps. Inputs require consistent `MAG_id` and `SGB_id` values. COG and
CAZy tables contain one row per MAG and numeric count columns; ODS/FDI requires
CAZy family counts rather than broad-class totals.

The network uses prevalence ≥0.05, mean abundance >1e-5, a CLR pseudocount of
1e-6, |rho| ≥0.30 and FDR <0.05. The Figure 4a phylogeny and its iTOL annotation
files are provided under `data/source_data/figure4/itol`.

## Requirements

R packages: dplyr, readr, tidyr, tibble, stringr, ggplot2, ggrepel, scales,
pheatmap, purrr, igraph, ggraph, ggforce, ggnewscale and RColorBrewer.
