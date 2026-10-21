# In vitro validation scripts

These scripts support the antibacterial, cell-viability and haemolysis analyses
presented in Figure 5b–d and the associated supplementary figures.

## Scripts

1. `01_plot_in_vitro_validation.R`: apply blank correction, calculate growth
   inhibition and MIC, summarize safety measurements and generate Figure 5b–d.
2. `02_plot_experimental_screening.R`: summarize antibacterial breadth, MIC
   distributions and experimental prioritization.
3. `03_plot_safety_profiles.R`: plot replicate-resolved cell-viability and
   haemolysis profiles.

## Usage

Scripts accept input and output paths as command-line arguments. Run step 1
first; its `mic_call_source_data.csv` output is required by step 2. Antibacterial
inputs consist of two 64 μmol L⁻¹ screening tables and two concentration-series
tables. Cell-viability and haemolysis inputs contain peptide, concentration,
replicate group and response columns.

All peptide concentrations are reported in μmol L⁻¹, following the manuscript.

## Requirements

R packages: ggplot2, patchwork, dplyr, tidyr, readr, stringr, scales, ggrepel
and ragg.
