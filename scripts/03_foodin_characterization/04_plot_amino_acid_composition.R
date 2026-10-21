#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(readr); library(tidyr); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('Usage: Rscript 04_plot_amino_acid_composition.R frequencies.csv output.pdf')
d <- read_csv(args[1], show_col_types = FALSE)
names(d)[1] <- 'Group'
d$Group[d$Group %in% c('c_AMP', 'cAMP(predicted)')] <- 'Foodins'
aa <- strsplit('ACDEFGHIKLMNPQRSTVWY', '')[[1]]
if (!all(aa %in% names(d))) stop('Missing amino acid frequency columns')
groups <- c('Foodins', 'DBAASP', 'Human', 'AMPSphere')
if (anyNA(d$Group) || !all(d$Group %in% groups)) stop('Unrecognized dataset label')
long <- pivot_longer(d, all_of(aa), names_to = 'AminoAcid', values_to = 'Frequency')
if (anyNA(long$Frequency) || any(long$Frequency < 0 | long$Frequency > 1)) stop('Invalid frequencies')
long$Group <- factor(long$Group, levels = groups)
colors <- c(Foodins = '#BD7795', DBAASP = '#91CCC0', Human = '#7FABD1', AMPSphere = '#F7AC53')
p <- ggplot(long, aes(AminoAcid, Frequency, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  scale_fill_manual(values = colors) + scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(x = 'Amino acid', y = 'Pooled residue frequency') + theme_classic() +
  theme(legend.position = 'top', axis.text = element_text(color = 'black'))
dir.create(dirname(args[2]), recursive = TRUE, showWarnings = FALSE)
ggsave(args[2], plot = p, width = 12, height = 6, useDingbats = FALSE)
