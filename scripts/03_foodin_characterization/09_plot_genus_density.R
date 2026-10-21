#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(readr); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('Usage: Rscript 09_plot_genus_density.R density.tsv output_dir')
d <- read_tsv(args[1], show_col_types = FALSE)
if (!all(c('genus', 'phylum', 'amps_per_Gbp') %in% names(d))) stop('Missing density columns')
if (anyDuplicated(d$genus) || anyNA(d$phylum) || any(!is.finite(d$amps_per_Gbp))) stop('Invalid genus density data')
# Original figure includes positive-density genera and phyla with >=4 genera.
d <- d %>% filter(amps_per_Gbp > 0) %>% add_count(phylum) %>% filter(n > 3)
if (!nrow(d)) stop('No phylum has at least four positive-density genera')
palette <- c('#91CCC0', '#7FABD1', '#BD7795', '#F7AC53', '#EC6E66', '#B5CE4E', '#C7C1DE', '#7C7979',
             '#963B79', '#97D0C5', '#52AADC', '#F39865', '#EEB6D4', '#D75B4E', '#C89736', '#2D8875')
phyla <- sort(unique(d$phylum))
if (length(phyla) > length(palette)) stop('More phyla than available palette colors')
labels <- d %>% group_by(phylum) %>% summarize(y = max(amps_per_Gbp) * 1.2, label = paste0('n=', n()), .groups = 'drop')
p <- ggplot(d, aes(phylum, amps_per_Gbp)) +
  geom_point(aes(color = phylum), position = position_jitter(width = 0.25, height = 0, seed = 42), size = 1.5, alpha = 0.8) +
  geom_boxplot(fill = NA, color = 'gray40', outlier.shape = NA, width = 0.6) +
  geom_text(data = labels, aes(phylum, y, label = label), inherit.aes = FALSE, size = 4) +
  scale_color_manual(values = setNames(palette[seq_along(phyla)], phyla)) + scale_y_log10() +
  labs(x = NULL, y = 'Genus-specific Foodin density\n(genes per assembled Gbp)') + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none')
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
write_tsv(d, file.path(args[2], 'genus_density_display_data.tsv'))
ggsave(file.path(args[2], 'genus_density.pdf'), plot = p, width = 9, height = 6, useDingbats = FALSE)
