#!/usr/bin/env Rscript
# Adapted from scripts/figure4; methodological changes are listed in README.md.
args <- commandArgs(trailingOnly=TRUE)
if(length(args) != 3) stop("Usage: Rscript 03_build_association_network.R summary.csv sample_ids.txt outdir")
script_file <- sub("^--file=", "", grep("^--file=", commandArgs(), value=TRUE)[1])
source(file.path(dirname(normalizePath(script_file)), "common.R"))
summary_file <- normalizePath(args[1], mustWork=TRUE)
sample_file <- normalizePath(args[2], mustWork=TRUE)
dir.create(args[3], recursive=TRUE, showWarnings=FALSE)
setwd(normalizePath(args[3]))
.Last <- function() writeLines(c(paste("Command:", paste(commandArgs(), collapse=" ")), capture.output(sessionInfo())), "run_session.txt")


library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(igraph)
library(ggraph)
library(ggplot2)
library(ggrepel)
library(ggforce)
library(ggnewscale)
library(tibble)
library(RColorBrewer)
library(grid)

df <- read_summary(summary_file)

meta_cols <- c(
  "order", "family", "genus", "species",
  "mean_abundance", "median_abundance", "prevalence",
  "n_sgb", "n_mag", "n_amp_mag", "amp_mag_fraction", "has_any_amp_mag", "sgb_category"
)

sample_cols <- read_samples(sample_file, df)

df <- df %>%
  mutate(across(all_of(sample_cols), as.numeric))

prev_cutoff  <- 0.05
abund_cutoff <- 1e-5

sp <- df %>%
  filter(prevalence >= prev_cutoff, mean_abundance > abund_cutoff) %>%
  mutate(
    species_label = paste(SGB_id, gsub("_", " ", species)),
    amp_group = ifelse(n_amp_mag > 0, "Foodin-positive", "Foodin-negative")
  )

cat("Filtered SGB number:", nrow(sp), "\n")

mat <- as.matrix(sp[, sample_cols])
if(nrow(sp)<2) stop("Fewer than two eligible SGBs")
rownames(mat) <- sp$SGB_id

mat_pseudo <- mat + 1e-6 

mat_clr <- apply(mat_pseudo, 2, function(x) {
  log(x) - mean(log(x))
})

mat_use <- mat_clr 
cat("CLR transformation applied successfully.\n")

pair_list <- combn(rownames(mat_use), 2, simplify = FALSE)

pairwise_spearman <- function(sp1, sp2, mat_input) {
  x <- as.numeric(mat_input[sp1, ])
  y <- as.numeric(mat_input[sp2, ])

  ct <- suppressWarnings(cor.test(x, y, method = "spearman", exact = FALSE))

  tibble(
    from = sp1,
    to = sp2,
    rho = unname(ct$estimate),
    p_value = ct$p.value
  )
}

cat("Calculating pairwise Spearman correlations (this may take a moment)...\n")
edge_raw <- map_dfr(pair_list, ~pairwise_spearman(.x[1], .x[2], mat_use))

edge_raw <- edge_raw %>%
  mutate(
    fdr = p.adjust(p_value, method = "BH"),
    association = case_when(
      rho > 0 ~ "positive",
      rho < 0 ~ "negative",
      TRUE ~ "neutral"
    )
  )

rho_cutoff <- 0.3
fdr_cutoff <- 0.05

pos_edges <- edge_raw %>%
  filter(
    association == "positive",
    rho >= rho_cutoff,
    fdr < fdr_cutoff
  ) %>%
  arrange(fdr)

neg_edges <- edge_raw %>%
  filter(
    association == "negative",
    rho <= -rho_cutoff,
    fdr < fdr_cutoff
  ) %>%
  arrange(fdr)

cat("Significant positive edges:", nrow(pos_edges), "\n")
cat("Significant negative edges:", nrow(neg_edges), "\n")

write_csv(edge_raw, "all_pair_correlations.csv")
write_csv(bind_rows(pos_edges, neg_edges), "significant_edges.csv")
write_csv(sp, "eligible_sgbs.csv")
if(nrow(pos_edges)+nrow(neg_edges)==0) {
  message("No significant edges at the prespecified thresholds; tables exported.")
  quit(save="no", status=0)
}

nodes <- sp %>%
  transmute(
    name = SGB_id,
    species_label = species_label,
    order = order,
    family = family,
    genus = genus,
    mean_abundance = mean_abundance,
    prevalence = prevalence,
    has_amp = as.logical(has_any_amp_mag),
    amp_group = amp_group
  )

all_significant_edges <- bind_rows(pos_edges, neg_edges)

g <- graph_from_data_frame(
  d = all_significant_edges %>% select(from, to, rho, fdr, association),
  vertices = nodes,
  directed = FALSE
)

g <- delete_vertices(g, V(g)[degree(g) == 0])

cat("Graph node number:", gorder(g), "\n")
cat("Graph edge number:", gsize(g), "\n")

set.seed(123)
g_pos_only <- delete_edges(g, E(g)[E(g)$association == "negative"])
cl_louvain <- cluster_louvain(g_pos_only, weights = E(g_pos_only)$rho)

V(g)$module <- membership(cl_louvain)[V(g)$name]

cat("Module number:", length(unique(V(g)$module)), "\n")
cat("Modularity:", modularity(cl_louvain), "\n")

ends_mat <- ends(g, E(g), names = FALSE)
edge_from_module <- V(g)$module[ends_mat[, 1]]
edge_to_module   <- V(g)$module[ends_mat[, 2]]

E(g)$edge_type <- case_when(
  E(g)$association == "negative" ~ "negative association",
  edge_from_module == edge_to_module ~ "within-module",
  TRUE ~ "between-module"
)

set.seed(123)
lay <- create_layout(g, layout = "stress")
lay_df <- lay %>% as_tibble()
write_csv(igraph::as_data_frame(g, what="vertices"), "network_nodes.csv")
write_csv(igraph::as_data_frame(g, what="edges"), "network_edges.csv")
write_csv(as.data.frame(lay)[, c("name", "x", "y", "module")], "network_layout.csv")

module_ids <- sort(unique(lay_df$module))
module_pal <- setNames(
  colorRampPalette(brewer.pal(8, "Pastel2"))(length(module_ids)),
  module_ids
)

hull_df <- lay_df %>%
  group_by(module) %>%
  mutate(module_size = n()) %>%
  ungroup() %>%
  filter(module_size >= 3)

label_df_full <- lay_df
label_df_simple <- lay_df %>%
  group_by(module) %>%
  slice_max(order_by = mean_abundance, n = 2, with_ties = FALSE) %>%
  ungroup()

make_network_plot <- function(layout_obj, hull_data, label_data, plot_title, plot_subtitle, out_file,
                              width = 12, height = 9, label_size = 3.5, 
                              point_size_range = c(3, 9), edge_width_range = c(0.5, 2.5)) {

  p <- ggraph(layout_obj) +
    ggforce::geom_mark_hull(
      data = hull_data,
      aes(x = x, y = y, group = module, fill = factor(module), colour = factor(module)),
      expand = unit(3, "mm"), concavity = 4, alpha = 0.22, linewidth = 0.7, show.legend = FALSE
    ) +
    scale_fill_manual(values = module_pal, guide = "none") +
    scale_colour_manual(values = module_pal, guide = "none") +
    ggnewscale::new_scale_colour() +
    ggnewscale::new_scale_fill() +

    geom_edge_link(
      aes(width = abs(rho), alpha = edge_type, colour = edge_type, linetype = edge_type),
      show.legend = TRUE
    ) +
    scale_edge_colour_manual(
      values = c(
        "within-module" = "grey40",
        "between-module" = "grey82",
        "negative association" = "tomato"
      ),
      name = "Edge type"
    ) +
    scale_edge_linetype_manual(
      values = c(
        "within-module" = "solid",
        "between-module" = "solid",
        "negative association" = "dashed"
      ),
      guide = "none"
    ) +
    scale_edge_alpha_manual(
      values = c(
        "within-module" = 0.65,
        "between-module" = 0.4,
        "negative association" = 0.8
      ),
      guide = "none"
    ) +
    scale_edge_width_continuous(range = edge_width_range, name = "Spearman |rho|") +

    geom_node_point(
      aes(size = prevalence, fill = amp_group),
      shape = 21, color = "grey20", stroke = 0.25, alpha = 0.95
    ) +
    scale_fill_manual(values = c("Foodin-positive" = "#EC6E66", "Foodin-negative" = "#7FABD1"), name = "Foodin status") +
    scale_size_continuous(range = point_size_range, name = "Prevalence") +

    geom_text_repel(
      data = label_data,
      aes(x = x, y = y, label = species_label),
      size = label_size, seed = 123, box.padding = 0.5, point.padding = 0.22, force = 2.8,
      max.iter = 30000, max.overlaps = Inf, min.segment.length = 0, segment.color = "grey50", segment.size = 0.22
    ) +
    labs(title = plot_title, subtitle = plot_subtitle, x = NULL, y = NULL) +
    theme_void(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12),
      legend.box = "vertical",
      plot.margin = margin(12, 20, 12, 20)
    )

  ggsave(out_file, p, width = width, height = height)
  print(p)
}

make_network_plot(lay, hull_df, label_df_full, "SGB abundance-correlation network", "Full-label version; CLR + Spearman + FDR < 0.05", "correlation_network_full_labels.pdf", width = 15, height = 11, label_size = 3.0)
make_network_plot(lay, hull_df, label_df_simple, "SGB abundance-correlation network", "Simplified-label version; CLR + Spearman + FDR < 0.05", "correlation_network_simple_labels.pdf", width = 12, height = 9, label_size = 3.8)
