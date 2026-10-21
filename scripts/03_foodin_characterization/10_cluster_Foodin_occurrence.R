#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(vegan); library(ape)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop('Usage: Rscript 10_cluster_Foodin_occurrence.R presence.csv output_dir')
raw <- read.csv(args[1], check.names = FALSE, stringsAsFactors = FALSE)
if (ncol(raw) < 3 || nrow(raw) < 2) stop('Need >=2 peptides and >=1 sample')
# Original file convention: column 1 = Sequence_ID, column 2 = Cluster_ID.
ids <- as.character(raw[[2]])
if (anyNA(ids) || any(!nzchar(ids)) || anyDuplicated(ids)) stop('Column 2 must contain unique non-empty tree labels')
x <- as.matrix(raw[, -c(1, 2), drop = FALSE])
if (!is.numeric(x) || anyNA(x) || !all(x %in% c(0, 1))) stop('Sample columns must be binary 0/1 without missing values')
if (any(rowSums(x) == 0)) stop('All-zero rows are not valid for this occurrence clustering')
rownames(x) <- ids
hc <- hclust(vegan::vegdist(x, method = 'bray'), method = 'average')
dir.create(args[2], recursive = TRUE, showWarnings = FALSE)
ape::write.tree(ape::as.phylo(hc), file = file.path(args[2], 'Foodin_occurrence_tree.nwk'))
saveRDS(hc, file.path(args[2], 'Foodin_occurrence_clustering.rds'))
# Base-R dendrogram, as in the source script; this is not a sequence phylogeny.
pdf(file.path(args[2], 'Foodin_occurrence_dendrogram.pdf'), width = 15, height = 40, useDingbats = FALSE)
plot(hc, hang = -1, cex = 0.3, main = 'Foodin clustering based on sample presence', xlab = '', sub = '')
dev.off()
