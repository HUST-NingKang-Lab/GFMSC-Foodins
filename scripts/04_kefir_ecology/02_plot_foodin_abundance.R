#!/usr/bin/env Rscript
# Adapted from scripts/figure4; methodological changes are listed in README.md.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) stop("Usage: Rscript 02_plot_foodin_abundance.R summary.csv mag_sgb.csv outdir")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "common.R"))
summary_file <- normalizePath(args[1], mustWork=TRUE)
map_file <- normalizePath(args[2], mustWork=TRUE)
dir.create(args[3], recursive=TRUE, showWarnings=FALSE)
setwd(normalizePath(args[3]))
.Last <- function() writeLines(c(paste("Command:", paste(commandArgs(), collapse=" ")), capture.output(sessionInfo())), "run_session.txt")

library(dplyr)
library(ggplot2)
library(ggrepel)
library(scales)
library(readr)
sgb_summary <- read_summary(summary_file)
mag_sgb <- read_mapping(map_file)
valid_sgbs <- unique(mag_sgb$SGB_id)

plot_df <- sgb_summary %>%
  filter(mean_abundance > 0) %>%
  filter(SGB_id %in% valid_sgbs) %>%
  mutate(
    species_label = paste0(SGB_id, "\n", gsub("_", " ", species)),
    amp_group = ifelse(n_amp_mag > 0, "Foodin-positive", "Foodin-negative")
  )

target_sgbs <- c(
  "SGB8002"
)

label_df <- plot_df %>%
  filter(SGB_id %in% target_sgbs) %>%
  distinct(SGB_id, .keep_all = TRUE)
label_df_pos <- label_df %>% filter(amp_mag_fraction > 0)
label_df_zero <- label_df %>% filter(amp_mag_fraction == 0)

p <- ggplot(plot_df, aes(x = mean_abundance, y = amp_mag_fraction)) +
  geom_point(
    aes(size = prevalence, fill = amp_group),
    shape = 21,
    color = "grey20",
    stroke = 0.35,
    alpha = 0.9
  ) +
  geom_text_repel(
    data = label_df_pos,
    aes(label = species_label),
    size = 3.5,
    seed = 123,
    box.padding = 0.5,
    point.padding = 0.35,
    force = 3,
    max.iter = 20000,
    min.segment.length = 0,
    segment.color = "grey60",
    segment.size = 0.3,
    max.overlaps = Inf
  ) +
  geom_text_repel(
    data = label_df_zero,
    aes(label = species_label),
    nudge_y = 0.05,
    size = 3.5,
    seed = 123,
    box.padding = 0.55,
    point.padding = 0.35,
    force = 4,
    max.iter = 20000,
    min.segment.length = 0,
    segment.color = "grey60",
    segment.size = 0.3,
    max.overlaps = Inf
  ) +
  scale_fill_manual(
    values = c(
      "Foodin-positive" = "#EC6E66",
      "Foodin-negative" = "#7FABD1"
    ),
    name = "Foodin status"
  ) +
  scale_x_log10(
    labels = label_scientific(digits = 2),
    expand = expansion(mult = c(0.08, 0.25))
  ) +
  scale_y_continuous(
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  scale_size_continuous(
    range = c(3, 10),
    name = "Prevalence"
  ) +
  coord_cartesian(clip = "off") +
  labs(
    x = "Mean abundance (log scale)",
    y = "Foodin-positive MAG fraction",
    title = "SGB level: abundance vs Foodin-bearing MAG fraction"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold", size = 15),
    axis.text = element_text(size = 12),
    legend.title = element_text(face = "bold", size = 13),
    legend.text = element_text(size = 12),
    plot.margin = margin(10, 50, 10, 10)
  )

ggsave(
  "kefir_foodin_bubble.pdf",
  p,
  width = 10,
  height = 6.8
)
cat("气泡图已保存为 kefir_foodin_bubble.pdf\n")
write_csv(plot_df, "bubble_source_data.csv")
