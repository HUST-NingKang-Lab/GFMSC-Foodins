#!/usr/bin/env Rscript
# Adapted from scripts/figure4; methodological changes are listed in README.md.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 5) stop("Usage: Rscript 06_calculate_ODS_FDI.R summary.csv mag_sgb.csv CAZy_family_counts.csv sample_ids.txt outdir")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "common.R"))
sgb_file <- normalizePath(args[1], mustWork=TRUE)
map_file <- normalizePath(args[2], mustWork=TRUE)
cazy_file <- normalizePath(args[3], mustWork=TRUE)
sample_file <- normalizePath(args[4], mustWork=TRUE)
dir.create(args[5], recursive=TRUE, showWarnings=FALSE)
setwd(normalizePath(args[5]))
.Last <- function() writeLines(c(paste("Command:", paste(commandArgs(), collapse=" ")), capture.output(sessionInfo())), "run_session.txt")


library(dplyr)
library(readr)
library(tidyr)
library(tibble)
library(stringr)
library(ggplot2)
library(ggrepel)

prevalence_cutoff <- 0.05
pseudocount <- 1e-6
rho_cutoff <- 0.30
fdr_cutoff <- 0.05
min_shared_samples <- 5

focal_species <- c(
  "s__Lactobacillus_kefiranofaciens",
  "s__Lentilactobacillus_kefiri",
  "s__Lactococcus_lactis",
  "s__Lactococcus_cremoris",
  "s__Lactobacillus_helveticus",
  "s__Streptococcus_thermophilus"
)
focal_species <- unique(c(focal_species, gsub("s__", "", focal_species)))

sgb_df  <- read_summary(sgb_file)
mag_sgb <- read_mapping(map_file)
cazy_df <- read_counts(cazy_file)

stopifnot(all(c("SGB_id", "species", "prevalence", "amp_mag_fraction") %in% colnames(sgb_df)))
stopifnot(all(c("MAG_id", "SGB_id") %in% colnames(mag_sgb)))
stopifnot("MAG_id" %in% colnames(cazy_df))

meta_cols <- c(
  "taxonomy", "order", "family", "genus", "species", "SGB_id",
  "mean_abundance", "median_abundance", "prevalence", "n_mag",
  "n_amp_mag", "amp_mag_fraction", "has_any_amp_mag", "sgb_category"
)

sample_cols <- read_samples(sample_file, sgb_df)

sgb_df <- sgb_df %>% mutate(across(all_of(sample_cols), as.numeric))

sgb_mat <- sgb_df %>%
  filter(prevalence >= prevalence_cutoff) %>%
  select(SGB_id, all_of(sample_cols)) %>%
  column_to_rownames("SGB_id") %>%
  as.matrix()

sgb_meta <- sgb_df %>%
  filter(SGB_id %in% rownames(sgb_mat)) %>%
  select(SGB_id, species, prevalence, amp_mag_fraction, sgb_category) %>%
  mutate(
    amp_group = case_when(
      str_detect(sgb_category, "High_AMP") ~ "AMP-high",
      str_detect(sgb_category, "Low_AMP")  ~ "AMP-low",
      TRUE ~ NA_character_
    )
  )

if(anyNA(sgb_meta$amp_group)) stop("Unrecognized sgb_category; use the summary from step 01")

clr_transform <- function(mat, pseudocount = 1e-6) {
  log_mat <- log(mat + pseudocount)
  clr_mat <- t(apply(log_mat, 2, function(x) x - mean(x)))
  return(clr_mat)
}

clr_mat <- clr_transform(sgb_mat, pseudocount = pseudocount)

sgb_names <- colnames(clr_mat)
if(length(sgb_names)<2) stop("Fewer than two eligible SGBs")
pair_idx <- combn(length(sgb_names), 2)
pair_list <- vector("list", ncol(pair_idx))

for (k in seq_len(ncol(pair_idx))) {
  i <- pair_idx[1, k]; j <- pair_idx[2, k]
  ct <- suppressWarnings(cor.test(clr_mat[, i], clr_mat[, j], method = "spearman", exact = FALSE))
  pair_list[[k]] <- tibble(sgb1 = sgb_names[i], sgb2 = sgb_names[j], rho = unname(ct$estimate), p_value = ct$p.value)
}

all_pairs <- bind_rows(pair_list) %>% mutate(fdr = p.adjust(p_value, method = "BH"))
write_csv(all_pairs, "all_pair_correlations.csv")
pos_pairs <- all_pairs %>%
  filter(rho >= rho_cutoff, fdr < fdr_cutoff)

cat("Retained positive SGB pairs:", nrow(pos_pairs), "\n")

cazy_cols <- setdiff(colnames(cazy_df), "MAG_id")
if(any(!grepl("^(GH|GT|CE|CBM|PL|AA)[0-9]+(_[0-9]+)?$", cazy_cols))) stop("FDI requires CAZy family counts, not broad classes")
cazy_long <- cazy_df %>% pivot_longer(cols = all_of(cazy_cols), names_to = "CAZy_family", values_to = "count")

mag_family_prop <- cazy_long %>%
  group_by(MAG_id) %>%
  mutate(total = sum(count)) %>%
  filter(total > 0) %>%
  mutate(prop = count / total) %>%
  ungroup()

sgb_family_mat <- mag_family_prop %>%
  inner_join(mag_sgb, by = "MAG_id") %>%
  filter(SGB_id %in% rownames(sgb_mat)) %>%
  group_by(SGB_id, CAZy_family) %>%
  summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = CAZy_family, values_from = prop, values_fill = 0) %>%
  column_to_rownames("SGB_id") %>%
  as.matrix()

sgb_family_mat <- sgb_family_mat[, colSums(sgb_family_mat > 0) >= 2, drop = FALSE]

write.csv(sgb_family_mat, "sgb_CAZy_family_proportions.csv", row.names=TRUE)
if(!ncol(sgb_family_mat)) stop("No CAZy families present in at least two eligible SGBs")

cosine_distance <- function(a, b) {
  if (sum(a) == 0 || sum(b) == 0) return(NA_real_)
  1 - (sum(a * b) / (sqrt(sum(a^2)) * sqrt(sum(b^2))))
}

presence_mat <- (sgb_mat > 0) * 1

pair_metrics <- pos_pairs %>%
  filter(sgb1 %in% rownames(sgb_family_mat), sgb2 %in% rownames(sgb_family_mat)) %>%
  left_join(sgb_meta, by = c("sgb1" = "SGB_id")) %>%
  rename(sp1 = species, amp1 = amp_mag_fraction, prev1 = prevalence, group1 = amp_group) %>%
  left_join(sgb_meta, by = c("sgb2" = "SGB_id")) %>%
  rename(sp2 = species, amp2 = amp_mag_fraction, prev2 = prevalence, group2 = amp_group) %>%
  rowwise() %>%
  mutate(
    builder = case_when(
      amp1 < amp2 ~ sgb1,
      amp2 < amp1 ~ sgb2,
      prev1 >= prev2 ~ sgb1,
      TRUE ~ sgb2
    ),
    partner = ifelse(builder == sgb1, sgb2, sgb1),
    builder_sp = ifelse(builder == sgb1, sp1, sp2),
    partner_sp = ifelse(partner == sgb1, sp1, sp2),
    builder_group = ifelse(builder == sgb1, group1, group2),
    partner_group = ifelse(partner == sgb1, group1, group2)
  ) %>%
  mutate(
    n_builder = sum(presence_mat[builder, ] == 1),
    n11 = sum(presence_mat[builder, ] == 1 & presence_mat[partner, ] == 1),
    ODS = ifelse(n_builder > 0, n11 / n_builder, NA_real_),
    FDI = cosine_distance(sgb_family_mat[builder, ], sgb_family_mat[partner, ])
  ) %>%
  ungroup() %>%
  filter(!is.na(ODS), !is.na(FDI), n11 >= min_shared_samples)

pair_metrics <- pair_metrics %>%
  mutate(
    pair_type = case_when(
      builder_group == "AMP-low"  & partner_group == "AMP-high" ~ "AMP-low / AMP-high",
      builder_group == "AMP-high" & partner_group == "AMP-high" ~ "AMP-high / AMP-high",
      TRUE ~ "other"
    ),
    pair_label = paste(gsub("s__", "", builder_sp), gsub("s__", "", partner_sp), sep = " | "),
    focal_pair = builder_sp %in% focal_species & partner_sp %in% focal_species
  )

write_csv(pair_metrics, "kefir_positive_pairs_ODS_FDI_table.csv")

p <- ggplot(pair_metrics, aes(x = ODS, y = FDI, color = pair_type)) +
  geom_point(size = 3.2, alpha = 0.85) +
  geom_text_repel(
    data = pair_metrics %>% filter(focal_pair),
    aes(label = pair_label),
    size = 3.3,
    box.padding = 0.4,
    point.padding = 0.25,
    max.overlaps = Inf,
    seed = 123,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = c(
      "AMP-low / AMP-high" = "#D95F02",
      "AMP-high / AMP-high" = "#E6AB02",
      "other" = "grey65"
    )
  ) +
  theme_classic(base_size = 13) +
  labs(
    x = "ODS (conditional co-occurrence dependency)",
    y = "FDI (functional divergence index; CAZy family)",
    color = "Pair type",
    title = "Positively associated SGB pairs in Kefir"
  ) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("kefir_ODS_FDI_pair_plot_SGB_2.pdf", p, width = 8.5, height = 6.5)

top_pairs <- pair_metrics %>%
  mutate(rank_score = ODS * FDI) %>%
  arrange(desc(rank_score)) %>%
  select(builder, partner, builder_sp, partner_sp, rho, ODS, FDI, pair_type, n11, rank_score)

write_csv(top_pairs, "kefir_top_pairs_by_ODSxFDI_SGB.csv")

cat("✅ Done! Output saved.\n")
