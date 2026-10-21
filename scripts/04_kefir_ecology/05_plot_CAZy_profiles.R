#!/usr/bin/env Rscript
# Adapted from scripts/figure4; methodological changes are listed in README.md.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 4) stop("Usage: Rscript 05_plot_CAZy_profiles.R counts.csv summary.csv mag_sgb.csv outdir")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "common.R"))
cazy_file <- normalizePath(args[1], mustWork=TRUE)
sgb_file <- normalizePath(args[2], mustWork=TRUE)
map_file <- normalizePath(args[3], mustWork=TRUE)
dir.create(args[4], recursive=TRUE, showWarnings=FALSE)
setwd(normalizePath(args[4]))
.Last <- function() writeLines(c(paste("Command:", paste(commandArgs(), collapse=" ")), capture.output(sessionInfo())), "run_session.txt")

library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(pheatmap)
library(ggplot2)

cazy_df <- read_counts(cazy_file)
sgb_df  <- read_summary(sgb_file)
mag_sgb <- read_mapping(map_file)

target_species <- c(
  "Lactobacillus_kefiranofaciens",
  "Lentilactobacillus_kefiri",
  "Lactococcus_lactis",
  "Lactococcus_cremoris",
  "Lactobacillus_helveticus",
  "Streptococcus_thermophilus"
)

stopifnot("MAG_id" %in% colnames(cazy_df))
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

cazy_cols <- setdiff(colnames(cazy_df), "MAG_id")

cazy_df <- cazy_df %>%
  mutate(across(all_of(cazy_cols), as.numeric))

broad_classes <- c("GH", "GT", "CE", "CBM", "PL", "AA")
is_broad_matrix <- all(cazy_cols %in% broad_classes)
if(!is_broad_matrix && any(!grepl("^(GH|GT|CE|CBM|PL|AA)[0-9]+(_[0-9]+)?$", cazy_cols))) stop("Use exclusively CAZy families or exclusively broad classes")

if (is_broad_matrix) {
  message("Detected CAZy broad-class matrix.")

  mag_broad <- cazy_df %>%
    select(MAG_id, all_of(cazy_cols))

  for (cls in broad_classes) {
    if (!cls %in% colnames(mag_broad)) {
      mag_broad[[cls]] <- 0
    }
  }

  mag_broad <- mag_broad %>%
    select(MAG_id, all_of(broad_classes))

  mag_family_prop <- NULL

} else {
  message("Detected CAZy family-level matrix.")

  cazy_long <- cazy_df %>%
    pivot_longer(
      cols = -MAG_id,
      names_to = "CAZy_family",
      values_to = "count"
    ) %>%
    mutate(
      CAZy_class = str_extract(CAZy_family, "^(GH|GT|CE|CBM|PL|AA)")
    ) %>%
    filter(!is.na(CAZy_class))

  mag_broad <- cazy_long %>%
    group_by(MAG_id, CAZy_class) %>%
    summarise(count = sum(count, na.rm = TRUE), .groups = "drop") %>%
    tidyr::complete(
      MAG_id,
      CAZy_class = broad_classes,
      fill = list(count = 0)
    ) %>%
    pivot_wider(names_from = CAZy_class, values_from = count) %>%
    select(MAG_id, all_of(broad_classes))

  mag_family_prop <- cazy_long %>%
    group_by(MAG_id) %>%
    mutate(total_CAZy = sum(count, na.rm = TRUE)) %>%
    ungroup() %>%
    filter(total_CAZy > 0) %>%
    mutate(prop = count / total_CAZy)
}

mag_broad_prop <- mag_broad %>%
  mutate(total_CAZy = rowSums(across(all_of(broad_classes)), na.rm = TRUE)) %>%
  filter(total_CAZy > 0) %>%
  mutate(
    GH  = GH  / total_CAZy,
    GT  = GT  / total_CAZy,
    CE  = CE  / total_CAZy,
    CBM = CBM / total_CAZy,
    PL  = PL  / total_CAZy,
    AA  = AA  / total_CAZy
  )

sgb_broad_prop <- mag_broad_prop %>%
  inner_join(mag_sgb %>% select(MAG_id, SGB_id), by = "MAG_id") %>%
  inner_join(target_sgb_info %>% select(SGB_id, species), by = "SGB_id") %>%
  group_by(SGB_id, species) %>%
  summarise(
    GH  = mean(GH, na.rm = TRUE),
    GT  = mean(GT, na.rm = TRUE),
    CE  = mean(CE, na.rm = TRUE),
    CBM = mean(CBM, na.rm = TRUE),
    PL  = mean(PL, na.rm = TRUE),
    AA  = mean(AA, na.rm = TRUE),
    total_CAZy = mean(total_CAZy, na.rm = TRUE),
    .groups = "drop"
  )

species_broad_prop <- sgb_broad_prop %>%
  group_by(species) %>%
  summarise(
    GH  = mean(GH, na.rm = TRUE),
    GT  = mean(GT, na.rm = TRUE),
    CE  = mean(CE, na.rm = TRUE),
    CBM = mean(CBM, na.rm = TRUE),
    PL  = mean(PL, na.rm = TRUE),
    AA  = mean(AA, na.rm = TRUE),
    total_CAZy = mean(total_CAZy, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(species = factor(species, levels = target_species)) %>%
  arrange(species)

write_csv(species_broad_prop, "selected_species_CAZy_broad_proportion_table_mean.csv")

mat_broad_prop <- species_broad_prop %>%
  select(species, GH, GT, CE, CBM, PL, AA) %>%
  column_to_rownames("species") %>%
  as.matrix()

row_zscore <- function(x) {
  if (length(x)<2 || sd(x, na.rm = TRUE) == 0) {
    return(rep(0, length(x)))
  } else {
    return((x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  }
}

mat_broad_z <- t(apply(mat_broad_prop, 1, row_zscore))

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

pdf("CAZy_broad_heatmap_selected_species_row_zscore_mean.pdf", width = 7.5, height = 4.8)

pheatmap(
  mat_broad_z,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = anno_row[, c("sgb_category"), drop = FALSE],
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = "grey85",
  main = "CAZy broad-class profile of selected species (row z-score, mean aggregation)"
)

dev.off()

pdf("CAZy_broad_heatmap_selected_species_raw_proportion_mean.pdf", width = 7.5, height = 4.8)

pheatmap(
  mat_broad_prop,
  cluster_rows = FALSE,
  cluster_cols = TRUE,
  annotation_row = anno_row[, c("sgb_category"), drop = FALSE],
  fontsize_row = 10,
  fontsize_col = 10,
  border_color = "grey85",
  main = "CAZy broad-class profile of selected species (raw proportion, mean aggregation)"
)

dev.off()

p_load <- species_broad_prop %>%
  ggplot(aes(x = species, y = total_CAZy, fill = species)) +
  geom_col(width = 0.7) +
  theme_classic(base_size = 13) +
  labs(
    x = NULL,
    y = "Mean total CAZy count per species",
    title = "Total CAZy load of selected species"
  ) +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", hjust = 0.5)
  )

ggsave("CAZy_total_load_selected_species_mean.pdf", p_load, width = 7, height = 4)

if (!is_broad_matrix) {

  sgb_family_prop <- mag_family_prop %>%
    inner_join(mag_sgb %>% select(MAG_id, SGB_id), by = "MAG_id") %>%
    inner_join(target_sgb_info %>% select(SGB_id, species), by = "SGB_id") %>%
    group_by(SGB_id, species, CAZy_family) %>%
    summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop")

  species_family_prop <- sgb_family_prop %>%
    group_by(species, CAZy_family) %>%
    summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop")

  top_n_family <- 12

  top_families <- species_family_prop %>%
    group_by(CAZy_family) %>%
    summarise(total_prop = sum(prop, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total_prop)) %>%
    slice(1:top_n_family) %>%
    pull(CAZy_family)

  species_family_prop2 <- species_family_prop %>%
    filter(CAZy_family %in% top_families) %>%
    mutate(species = factor(species, levels = target_species)) %>%
    arrange(species)

  write_csv(species_family_prop2, "selected_species_CAZy_top_family_table_mean.csv")

  mat_family_prop <- species_family_prop2 %>%
    pivot_wider(names_from = CAZy_family, values_from = prop, values_fill = 0) %>%
    arrange(species) %>%
    column_to_rownames("species") %>%
    as.matrix()

  mat_family_z <- t(apply(mat_family_prop, 1, row_zscore))
zscore_colors <- colorRampPalette(c("#7FABD1", "white", "#EC6E66"))(100)

  pdf("CAZy_top_family_heatmap_selected_species_row_zscore_mean.pdf", width = 10, height = 4.8)

  pheatmap(
    mat_family_z,
    cluster_rows = FALSE,
    cluster_cols = ncol(mat_family_prop)>1,
    annotation_row = anno_row[, c("sgb_category"), drop = FALSE],
    color = zscore_colors,
    fontsize_row = 10,
    fontsize_col = 9,
    border_color = "grey85",
    main = "Top CAZy family profile of selected species (row z-score, mean aggregation)"
  )

  dev.off()

  pdf("CAZy_top_family_heatmap_selected_species_raw_proportion_mean.pdf", width = 10, height = 4.8)

  pheatmap(
    mat_family_prop,
    cluster_rows = FALSE,
    cluster_cols = ncol(mat_family_prop)>1,
    annotation_row = anno_row[, c("sgb_category"), drop = FALSE],
    fontsize_row = 10,
    fontsize_col = 9,
    border_color = "grey85",
    main = "Top CAZy family profile of selected species (raw proportion, mean aggregation)"
  )

  dev.off()
}

cat("✅ Done!\n")
cat("输出文件：\n")
cat("1) selected_species_CAZy_broad_proportion_table_mean.csv\n")
cat("2) CAZy_broad_heatmap_selected_species_row_zscore_mean.pdf\n")
cat("3) CAZy_broad_heatmap_selected_species_raw_proportion_mean.pdf\n")
cat("4) CAZy_total_load_selected_species_mean.pdf\n")
if (!is_broad_matrix) {
  cat("5) selected_species_CAZy_top_family_table_mean.csv\n")
  cat("6) CAZy_top_family_heatmap_selected_species_row_zscore_mean.pdf\n")
  cat("7) CAZy_top_family_heatmap_selected_species_raw_proportion_mean.pdf\n")
}
