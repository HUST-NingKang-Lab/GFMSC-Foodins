#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(readr); library(dplyr)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop('Usage: Rscript 07_summarize_taxonomic_distribution.R rank.csv genus.csv output_dir')
tax <- read_csv(args[1], show_col_types = FALSE)
gen <- read_csv(args[2], show_col_types = FALSE)
if (!'rank' %in% names(tax) || !'genus' %in% names(gen)) stop('Inputs require rank and genus columns, respectively')
rank_levels <- c('no rank', 'superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
if (anyNA(tax$rank) || !all(tax$rank %in% rank_levels)) stop('Unexpected/missing ranks; resolve the input annotations first')
rank_data <- tax %>% count(rank, name = 'Count') %>% mutate(Percentage = Count / sum(Count) * 100)
rank_data$rank <- factor(rank_data$rank, levels = rank_levels)
colors <- c('no rank' = '#666666', superkingdom = '#91CCC0', phylum = '#7FABD1', class = '#F7AC53',
            order = '#C7C1DE', family = '#EC6E66', genus = '#2D8875', species = '#BD7795')
p_rank <- ggplot(rank_data, aes('', Percentage, fill = rank)) + geom_col(width = 0.5) +
  scale_fill_manual(values = colors) + scale_y_continuous(limits = c(0, 100.1), expand = c(0, 0)) +
  labs(x = NULL, y = 'Foodin records (%)', fill = 'Taxonomic rank') + theme_classic() +
  theme(axis.line.x = element_blank(), axis.ticks.x = element_blank(), axis.text.x = element_blank())
genus_data <- gen %>% filter(!is.na(genus), genus != '') %>% count(genus, name = 'count') %>% arrange(desc(count), genus)
top <- slice_head(genus_data, n = 10)
p_genus <- ggplot(top, aes(reorder(genus, count), count)) + geom_col(fill = '#91CCC0', width = 0.7) +
  geom_text(aes(label = count), vjust = -0.5, size = 3.5) + labs(x = 'Genus', y = 'Foodin records') +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dir.create(args[3], recursive = TRUE, showWarnings = FALSE)
write_csv(rank_data, file.path(args[3], 'taxonomic_rank_summary.csv'))
write_csv(genus_data, file.path(args[3], 'genus_counts.csv'))
ggsave(file.path(args[3], 'taxonomic_ranks.pdf'), plot = p_rank, width = 5, height = 8, useDingbats = FALSE)
ggsave(file.path(args[3], 'top10_genera.pdf'), plot = p_genus, width = 8, height = 6, useDingbats = FALSE)
