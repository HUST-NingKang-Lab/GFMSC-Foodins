# In vivo validation scripts

This module supports the murine deep-thigh infection analysis presented in
Figure 6.

## Scripts

1. `01_analyze_in_vivo_validation.R`: calculate mouse-level bacterial burden
   and cytokine measurements, perform one-way ANOVA with Dunnett comparisons,
   and generate Figure 6 and its source-data tables.

## Usage

The script accepts the current CFU workbook, ELISA workbook and output directory
as command-line arguments. Two technical plating replicates are averaged for
each of six mice per group. ELISA concentrations represent three mice per group,
with each mouse calculated from duplicate wells.

The included groups are infection control, Foodin_3, Foodin_21, Foodin_39,
Foodin_40 and polymyxin B. Bacterial burden is analysed as log10 CFU g⁻¹ muscle;
serum TNF-α, IL-6 and IL-10 are reported in pg ml⁻¹.

## Requirements

R packages: readxl, dplyr, tidyr, ggplot2, patchwork, emmeans and ragg.
