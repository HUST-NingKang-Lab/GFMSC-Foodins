#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(ggplot2); library(readr); library(dplyr); library(tidyr)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) stop('Usage: Rscript 06_summarize_predicted_toxicity.R output_dir Group=result.csv [Group=result.csv ...]')
groups <- c('Foodins', 'DBAASP', 'Human', 'AMPSphere')
tables <- lapply(args[-1], function(spec) {
  pos <- regexpr('=', spec, fixed = TRUE)[1]
  if (pos < 2) stop('Expected Group=result.csv')
  group <- substr(spec, 1, pos - 1)
  if (group %in% c('c_AMP', 'cAMP(predicted)')) group <- 'Foodins'
  if (!group %in% groups) stop('Unrecognized dataset label')
  d <- read_csv(substr(spec, pos + 1, nchar(spec)), col_types = cols(.default = col_character()))
  if (!'Prediction' %in% names(d) || !nrow(d) || anyNA(d$Prediction)) stop('Missing/empty Prediction data')
  label <- tolower(gsub('[^a-z]', '', tolower(d$Prediction)))
  if (!all(label %in% c('toxic', 'nontoxic', 'toxin', 'nontoxin'))) stop('Unrecognized Prediction labels')
  data.frame(Group = group, Prediction = ifelse(label %in% c('toxic', 'toxin'), 'Toxic', 'Non-Toxic'))
})
d <- bind_rows(tables) %>% count(Group, Prediction, name = 'n') %>%
  complete(Group, Prediction = c('Non-Toxic', 'Toxic'), fill = list(n = 0)) %>%
  group_by(Group) %>% mutate(Percentage = n / sum(n)) %>% ungroup()
d$Group <- factor(d$Group, levels = groups)
d$Prediction <- factor(d$Prediction, levels = c('Non-Toxic', 'Toxic'))
p <- ggplot(d, aes(Group, Percentage, fill = Prediction)) + geom_col(width = 0.7) +
  geom_text(aes(label = scales::percent(Percentage, accuracy = 0.1)), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values = c('Non-Toxic' = '#91CCC0', Toxic = '#7FABD1')) +
  scale_y_continuous(labels = scales::percent) + labs(x = NULL, y = 'Proportion', fill = 'Prediction') + theme_classic()
dir.create(args[1], recursive = TRUE, showWarnings = FALSE)
write_csv(d, file.path(args[1], 'predicted_toxicity_summary.csv'))
ggsave(file.path(args[1], 'predicted_toxicity.pdf'), plot = p, width = 8, height = 6, useDingbats = FALSE)
