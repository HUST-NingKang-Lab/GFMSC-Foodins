#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(readr); library(tidyr); library(dplyr); library(ggridges)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('Usage: Rscript 02_plot_physicochemical_properties.R features.csv output_dir')
d <- read_csv(args[1], show_col_types = FALSE)
features <- c('Aliphatic Index', 'Amphiphilicity', 'Boman Index', 'Hydrophobicity', 'Isoelectric Point', 'Net Charge')
if (!all(c('Group', features) %in% names(d))) stop('Missing Group or one of the six properties')
d$Group[d$Group %in% c('c_AMP', 'cAMP(predicted)')] <- 'Foodins'
groups <- c('Foodins', 'DBAASP', 'Human', 'AMPSphere')
if (anyNA(d$Group) || !all(d$Group %in% groups)) stop('Unrecognized dataset label')
long <- d %>% pivot_longer(all_of(features), names_to = 'Feature', values_to = 'Value')
if (any(!is.finite(long$Value))) stop('Properties must be finite numeric values')
# Retain the original per-group 1st-99th percentile trimming for display only.
shown <- long %>% group_by(Feature, Group) %>%
  filter(Value >= quantile(Value, 0.01), Value <= quantile(Value, 0.99)) %>% ungroup()
shown$Group <- factor(shown$Group, levels = rev(groups))
shown$Feature <- factor(shown$Feature, levels = features)
colors <- c(Foodins = '#BD7795', DBAASP = '#91CCC0', Human = '#7FABD1', AMPSphere = '#F7AC53')
p <- ggplot(shown, aes(Value, Group, fill = Group, color = Group)) +
  geom_density_ridges(scale = 1.5, rel_min_height = 0.01, alpha = 0.6) +
  facet_wrap(~Feature, scales = 'free_x', nrow = 2) +
  scale_fill_manual(values = colors) + scale_color_manual(values = colors) +
  labs(x = 'Feature value', y = NULL) + theme_classic() +
  theme(strip.background = element_blank(), strip.text = element_text(face = 'bold'),
        axis.text = element_text(color = 'black'), legend.position = 'none')
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
write_csv(long, file.path(args[2], 'physicochemical_six_properties.csv'))
write_csv(shown, file.path(args[2], 'physicochemical_display_data.csv'))
ggsave(file.path(args[2], 'physicochemical_six_properties.pdf'), plot = p, width = 12, height = 8, useDingbats = FALSE)
