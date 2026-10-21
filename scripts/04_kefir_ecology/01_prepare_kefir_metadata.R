#!/usr/bin/env Rscript
# Adapted from scripts/figure4; methodological changes are listed in README.md.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) stop("Usage: Rscript 01_prepare_kefir_metadata.R abundance.csv foodin_mags.csv mag_sgb.csv outdir")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "common.R"))
abundance_file <- normalizePath(args[1], mustWork=TRUE)
mag_amp_file <- normalizePath(args[2], mustWork=TRUE)
mag_sgb_file <- normalizePath(args[3], mustWork=TRUE)
dir.create(args[4], recursive=TRUE, showWarnings=FALSE)
setwd(normalizePath(args[4]))
.Last <- function() writeLines(c(paste("Command:", paste(commandArgs(), collapse=" ")), capture.output(sessionInfo())), "run_session.txt")


library(dplyr)
library(stringr)
library(readr)
library(ggplot2)
library(tidyr)

ab <- read_csv(abundance_file, show_col_types = FALSE)
mag_amp <- read_csv(mag_amp_file, show_col_types = FALSE)
mag_sgb <- read_mapping(mag_sgb_file)

stopifnot("taxonomy" %in% colnames(ab))
stopifnot("MAG_id" %in% colnames(mag_amp))
stopifnot(all(c("MAG_id", "SGB_id") %in% colnames(mag_sgb)))

extract_rank <- function(x, prefix){
  out <- str_match(x, paste0(prefix, "([^|]+)"))[,2]
  return(out)
}

ab2 <- ab %>%
  mutate(
    order   = extract_rank(taxonomy, "o__"),
    family  = extract_rank(taxonomy, "f__"),
    genus   = extract_rank(taxonomy, "g__"),
    species = extract_rank(taxonomy, "s__"),
    SGB_id  = extract_rank(taxonomy, "t__")
  )

meta_cols <- c("taxonomy", "order", "family", "genus", "species", "SGB_id")
sample_cols <- setdiff(colnames(ab2), meta_cols)

check_id(ab2, "SGB_id")
check_numeric(ab2, sample_cols)
if(length(sample_cols)<3) stop("At least three samples required")
if(length(setdiff(mag_amp$MAG_id, mag_sgb$MAG_id))) stop("Foodin MAG IDs missing from mapping")
writeLines(sample_cols, "sample_ids.txt")
sgb_label <- ab2 %>%
  select(SGB_id, taxonomy, order, family, genus, species) %>%
  distinct()

mag_label <- mag_sgb %>%
  distinct(MAG_id, SGB_id) %>%
  mutate(has_amp = MAG_id %in% mag_amp$MAG_id) %>%
  left_join(
    sgb_label %>% select(SGB_id, order, family, genus, species),
    by = "SGB_id"
  ) %>%
  arrange(species, SGB_id, MAG_id)

sgb_abundance <- ab2 %>%
  rowwise() %>%
  mutate(
    mean_abundance   = mean(c_across(all_of(sample_cols)), na.rm = TRUE),
    median_abundance = median(c_across(all_of(sample_cols)), na.rm = TRUE),
    prevalence       = mean(c_across(all_of(sample_cols)) > 0, na.rm = TRUE)
  ) %>%
  ungroup()

sgb_amp <- mag_label %>%
  group_by(SGB_id) %>%
  summarise(
    n_mag = n_distinct(MAG_id),
    n_amp_mag = sum(has_amp, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    amp_mag_fraction = ifelse(n_mag > 0, n_amp_mag / n_mag, 0)
  )

sgb_summary <- sgb_abundance %>%
  left_join(sgb_amp, by = "SGB_id") %>%
  mutate(
    n_mag = replace_na(n_mag, 0L),
    n_amp_mag = replace_na(n_amp_mag, 0L),
    amp_mag_fraction = replace_na(amp_mag_fraction, 0),
    has_any_amp_mag = n_amp_mag > 0,
    sgb_category = case_when(
      amp_mag_fraction >= 0.1 & mean_abundance >= 0.1 ~ "High_AMP_High_Abundance",
      amp_mag_fraction >= 0.1 & mean_abundance < 0.1 ~ "High_AMP_Low_Abundance",
      amp_mag_fraction < 0.1 & mean_abundance >= 0.1  ~ "Low_AMP_High_Abundance",
      amp_mag_fraction < 0.1 & mean_abundance < 0.1   ~ "Low_AMP_Low_Abundance",
      TRUE ~ "Unknown"
    )
  ) %>%
  relocate(all_of(sample_cols), .after = last_col())

write_csv(sgb_label, "sgb_label_table.csv")
write_csv(mag_label, "mag_label_table.csv")
write_csv(sgb_summary, "sgb_summary_with_amp.csv")
write_csv(sgb_summary %>% filter(n_mag == 0), "sgbs_without_mapped_MAGs.csv")
