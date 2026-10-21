#!/usr/bin/env Rscript
# Adapted from scripts/figure4; methodological changes are listed in README.md.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) stop("Usage: Rscript 04_plot_COG_profiles.R counts.csv summary.csv mag_sgb.csv outdir")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "common.R"))
cog_file <- normalizePath(args[1], mustWork=TRUE)
sgb_file <- normalizePath(args[2], mustWork=TRUE)
map_file <- normalizePath(args[3], mustWork=TRUE)
dir.create(args[4], recursive=TRUE, showWarnings=FALSE)
setwd(normalizePath(args[4]))
.Last <- function() writeLines(c(paste("Command:", paste(commandArgs(), collapse=" ")), capture.output(sessionInfo())), "run_session.txt")

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(pheatmap)

cog_df <- read_counts(cog_file)
sgb_df <- read_summary(sgb_file)
mag_sgb <- read_mapping(map_file)

target_species <- c(
  "Lactobacillus_kefiranofaciens",
  "Lentilactobacillus_kefiri",
  "Lactococcus_lactis",
  "Lactococcus_cremoris",
  "Lactobacillus_helveticus",
  "Streptococcus_thermophilus"
)

stopifnot("MAG_id" %in% colnames(cog_df))
stopifnot(all(c("SGB_id", "species") %in% colnames(sgb_df)))
stopifnot(all(c("MAG_id", "SGB_id") %in% colnames(mag_sgb)))

target_sgb_info <- sgb_df %>%
  filter(species %in% target_species) %>%
  select(SGB_id, species, mean_abundance, prevalence, amp_mag_fraction, sgb_category) %>%
  distinct()

if (nrow(target_sgb_info) == 0) {
  stop("No target species found in sgb_summary_with_amp.csv")
}

target_sgb_info <- target_sgb_info %>%
  mutate(species = factor(species, levels = target_species)) %>%
  arrange(species)

missing_species <- setdiff(target_species, unique(as.character(target_sgb_info$species)))
if (length(missing_species) > 0) {
  message("These target species were NOT found in sgb_summary_with_amp.csv:")
  print(missing_species)
}

cog_cols <- setdiff(colnames(cog_df), "MAG_id")

cog_df <- cog_df %>%
  mutate(across(all_of(cog_cols), as.numeric))

use_cog_cols <- setdiff(cog_cols, intersect(c("R", "S"), cog_cols))

mag_cog_prop <- cog_df %>%
  mutate(total_cog = rowSums(across(all_of(use_cog_cols)), na.rm = TRUE)) %>%
  filter(total_cog > 0) %>%
  mutate(across(all_of(use_cog_cols), ~ .x / total_cog)) %>%
  select(MAG_id, all_of(use_cog_cols))

sgb_cog_prop <- mag_cog_prop %>%
  inner_join(mag_sgb %>% select(MAG_id, SGB_id), by = "MAG_id") %>%
  inner_join(target_sgb_info %>% select(SGB_id, species), by = "SGB_id") %>%
  group_by(SGB_id, species) %>%
  summarise(across(all_of(use_cog_cols), mean, na.rm = TRUE), .groups = "drop")

species_cog_prop <- sgb_cog_prop %>%
  group_by(species) %>%
  summarise(across(all_of(use_cog_cols), mean, na.rm = TRUE), .groups = "drop") %>%
  mutate(species = factor(species, levels = target_species)) %>%
  arrange(species)

min_max_prop <- 0.05

keep_cols <- use_cog_cols[
  sapply(use_cog_cols, function(x) max(species_cog_prop[[x]], na.rm = TRUE) >= min_max_prop)
]

if(!length(keep_cols)) stop("No COG categories meet the 0.05 threshold")
species_cog_prop2 <- species_cog_prop %>%
  select(species, all_of(keep_cols))

write_csv(species_cog_prop2, "selected_species_COG_proportion_table_mean.csv")

mat_prop <- species_cog_prop2 %>%
  column_to_rownames("species") %>%
  as.matrix()

row_zscore <- function(x) {
  if (length(x)<2 || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  } else {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  }
}

mat_z <- t(apply(mat_prop, 1, row_zscore))

anno_row <- target_sgb_info %>%
  group_by(species) %>%
  summarise(
    mean_abundance = mean(mean_abundance, na.rm = TRUE),
    prevalence = mean(prevalence, na.rm = TRUE),
    amp_mag_fraction = mean(amp_mag_fraction, na.rm = TRUE),
    sgb_category = paste(unique(sgb_category), collapse = ";"),
    .groups = "drop"
  ) %>%
  mutate(species = factor(species, levels = target_species)) %>%
  arrange(species) %>%
  column_to_rownames("species")
zscore_colors <- colorRampPalette(c("#7FABD1", "white", "#EC6E66"))(100)
pdf("COG_heatmap_selected_species_row_zscore_mean.pdf", width = 10, height = 4.8)

pheatmap(
  mat_z,
  cluster_rows = FALSE,
  cluster_cols = ncol(mat_prop) > 1,
  annotation_row = anno_row[, c("sgb_category"), drop = FALSE],
  color = zscore_colors,
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = "grey85",
  main = "COG category profile of selected species (row z-score, mean aggregation)"
)

dev.off()

pdf("COG_heatmap_selected_species_raw_proportion_mean.pdf", width = 10, height = 4.8)

pheatmap(
  mat_prop,
  cluster_rows = FALSE,
  cluster_cols = ncol(mat_prop) > 1,
  annotation_row = anno_row[, c("sgb_category"), drop = FALSE],
  color = zscore_colors,
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = "grey85",
  main = "COG category profile of selected species (raw proportion, mean aggregation)"
)

dev.off()

long_df <- species_cog_prop2 %>%
  pivot_longer(cols = -species, names_to = "COG_category", values_to = "proportion")

write_csv(long_df, "selected_species_COG_proportion_long_mean.csv")

cat("✅ Done!\n")
cat("输出文件：\n")
cat("1) selected_species_COG_proportion_table_mean.csv\n")
cat("2) selected_species_COG_proportion_long_mean.csv\n")
cat("3) COG_heatmap_selected_species_row_zscore_mean.pdf\n")
cat("4) COG_heatmap_selected_species_raw_proportion_mean.pdf\n")
