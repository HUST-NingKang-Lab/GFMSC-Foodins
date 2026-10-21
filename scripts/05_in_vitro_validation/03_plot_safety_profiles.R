#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(scales)
  library(patchwork)
  library(ragg)
})

# Replicate-resolved safety profiling of the toxicity-tested Foodins.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste(
    "Usage: Rscript 03_plot_safety_profiles.R",
    "cell_viability.csv haemolysis.csv output_dir"
  ))
}

activity_file <- normalizePath(args[1], mustWork = TRUE)
haemolysis_file <- normalizePath(args[2], mustWork = TRUE)
output_dir <- args[3]
panel_dir <- file.path(output_dir, "individual_panels")
source_dir <- file.path(output_dir, "source_data")

dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(activity_file), file.exists(haemolysis_file))

pal <- c(
  teal = "#2D8875",
  teal_light = "#91CCC0",
  coral = "#EC6E66",
  coral_light = "#F4B0AA",
  grey = "#7C7979",
  light_grey = "#F1F1F1",
  black = "#252525"
)

concentration_levels <- c(8, 16, 32, 64, 128)
safe_foodins <- c("Food3", "Food21", "Food26", "Food39", "Food40", "Food46")

format_peptide <- function(x) {
  ifelse(x == "LL37", "LL37", paste("Foodin", str_extract(x, "[0-9]+")))
}

activity_raw <- read_csv(activity_file, show_col_types = FALSE) %>%
  transmute(
    peptide,
    concentration = as.numeric(concentration),
    replicate = group,
    raw_value = as.numeric(activity),
    # Values above one are capped at 100%, as in the main-figure analysis.
    display_value = pmin(pmax(raw_value, 0), 1)
  )

haemolysis_raw <- read_csv(haemolysis_file, show_col_types = FALSE) %>%
  transmute(
    peptide,
    concentration = as.numeric(concentration),
    replicate = group,
    raw_value = as.numeric(Haemolytic),
    display_value = pmax(raw_value, 0)
  )

# Keep LL37 first, followed by Foodins in numerical order.
peptide_order <- union(activity_raw$peptide, haemolysis_raw$peptide)
peptide_order <- c(
  "LL37",
  peptide_order[peptide_order != "LL37"][
    order(as.numeric(str_extract(peptide_order[peptide_order != "LL37"], "[0-9]+")))
  ]
)
peptide_labels <- setNames(format_peptide(peptide_order), peptide_order)

prepare_endpoint <- function(dat) {
  dat %>%
    mutate(
      peptide = factor(peptide, levels = peptide_order),
      concentration = factor(concentration, levels = concentration_levels),
      replicate = factor(replicate, levels = c("group1", "group2", "group3"))
    )
}

activity_plot_data <- prepare_endpoint(activity_raw)
haemolysis_plot_data <- prepare_endpoint(haemolysis_raw)

activity_mean <- activity_plot_data %>%
  group_by(peptide, concentration) %>%
  summarise(
    mean_value = mean(display_value, na.rm = TRUE),
    sd_value = sd(display_value, na.rm = TRUE),
    n = sum(!is.na(display_value)),
    .groups = "drop"
  )

haemolysis_mean <- haemolysis_plot_data %>%
  group_by(peptide, concentration) %>%
  summarise(
    mean_value = mean(display_value, na.rm = TRUE),
    sd_value = sd(display_value, na.rm = TRUE),
    n = sum(!is.na(display_value)),
    .groups = "drop"
  )

# Data-integrity checks: 20 Foodins plus LL37, five concentrations,
# three replicates for Foodins and two replicates for LL37.
count_conditions <- function(dat) {
  dat %>% count(peptide, concentration, name = "n_replicates")
}

activity_n <- count_conditions(activity_plot_data)
haemolysis_n <- count_conditions(haemolysis_plot_data)

stopifnot(
  n_distinct(activity_plot_data$peptide) == 21,
  n_distinct(haemolysis_plot_data$peptide) == 21,
  n_distinct(filter(activity_plot_data, peptide != "LL37")$peptide) == 20,
  all(filter(activity_n, peptide != "LL37")$n_replicates == 3),
  all(filter(haemolysis_n, peptide != "LL37")$n_replicates == 3),
  all(filter(activity_n, peptide == "LL37")$n_replicates == 2),
  all(filter(haemolysis_n, peptide == "LL37")$n_replicates == 2),
  all(concentration_levels %in% activity_plot_data$concentration),
  all(concentration_levels %in% haemolysis_plot_data$concentration)
)

write_csv(activity_raw, file.path(source_dir, "FigS9a_cell_viability_replicates.csv"))
write_csv(activity_mean, file.path(source_dir, "FigS9a_cell_viability_means.csv"))
write_csv(haemolysis_raw, file.path(source_dir, "FigS9b_haemolysis_replicates.csv"))
write_csv(haemolysis_mean, file.path(source_dir, "FigS9b_haemolysis_means.csv"))

theme_supp <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base_size + 0.8,
                                margin = margin(b = 4)),
      plot.tag = element_text(face = "bold", size = base_size + 1.8),
      panel.grid = element_blank(),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.4),
      legend.key.width = unit(7, "mm"),
      legend.key.height = unit(3.2, "mm"),
      legend.spacing.x = unit(1.5, "mm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

# Extend the manuscript palette to all 21 peptides while retaining its hues.
base_foodin_palette <- c(
  "#91CCC0", "#7FABD1", "#F7AC53", "#EC6E66", "#B5CE4E",
  "#BD7795", "#963B79", "#97D0C5", "#F39865", "#52AADC",
  "#C7C1DE", "#EEB6D4", "#C89736", "#2D8875", "#D75B4E",
  "#7C7979"
)
foodin_colours <- grDevices::colorRampPalette(base_foodin_palette)(length(peptide_order))
names(foodin_colours) <- peptide_order
foodin_colours["LL37"] <- "#252525"

make_endpoint_plot <- function(mean_dat, title, y_label, upper_limit) {
  mean_dat <- mean_dat %>%
    mutate(
      peptide = factor(peptide, levels = peptide_order),
      concentration_numeric = as.numeric(as.character(concentration)),
      ymin = pmax(mean_value - sd_value, 0),
      ymax = pmin(mean_value + sd_value, upper_limit)
    )

  ggplot(mean_dat, aes(x = concentration_numeric, y = mean_value,
                       colour = peptide, group = peptide)) +
    geom_vline(xintercept = 64, linetype = "dashed", linewidth = 0.35,
               colour = pal[["grey"]]) +
    geom_hline(yintercept = 0.5, linetype = "longdash", linewidth = 0.35,
               colour = pal[["grey"]]) +
    geom_errorbar(aes(ymin = ymin, ymax = ymax), width = 0.12,
                  linewidth = 0.16, alpha = 0.42) +
    geom_line(linewidth = 0.45, alpha = 0.82) +
    geom_point(size = 1.35, alpha = 0.9) +
    scale_x_continuous(
      trans = "log2",
      breaks = concentration_levels,
      labels = concentration_levels,
      expand = expansion(mult = c(0.015, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, upper_limit),
      breaks = c(0, 0.5, 1),
      labels = percent_format(accuracy = 1),
      expand = expansion(mult = c(0, 0.04)),
      oob = squish
    ) +
    scale_colour_manual(values = foodin_colours, labels = peptide_labels,
                        name = NULL, drop = FALSE) +
    labs(title = title, x = expression(Concentration~(mu*mol~L^{-1})),
         y = y_label) +
    guides(colour = guide_legend(
      nrow = 3, byrow = TRUE,
      override.aes = list(linewidth = 0.55, size = 1.5, alpha = 1)
    )) +
    theme_supp() +
    theme(legend.position = "bottom")
}

p_a <- make_endpoint_plot(
  activity_mean,
  title = "Cell-viability profiles",
  y_label = "Cell viability",
  upper_limit = 1.05
)

p_b <- make_endpoint_plot(
  haemolysis_mean,
  title = "Haemolysis profiles",
  y_label = "Haemolysis",
  upper_limit = 1.10
)

fig_s9 <- (
  (p_a + theme(legend.position = "none")) /
    p_b +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
)

save_pdf <- function(plot, path, width_mm, height_mm) {
  ggsave(
    filename = path,
    plot = plot,
    device = grDevices::pdf,
    width = width_mm,
    height = height_mm,
    units = "mm",
    family = "Helvetica",
    useDingbats = FALSE,
    bg = "white",
    limitsize = FALSE
  )
}

save_png <- function(plot, path, width_mm, height_mm) {
  ggsave(
    filename = path,
    plot = plot,
    device = ragg::agg_png,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    bg = "white",
    limitsize = FALSE
  )
}

save_tiff <- function(plot, path, width_mm, height_mm) {
  ggsave(
    filename = path,
    plot = plot,
    device = ragg::agg_tiff,
    width = width_mm,
    height = height_mm,
    units = "mm",
    dpi = 600,
    compression = "lzw",
    bg = "white",
    limitsize = FALSE
  )
}

stem <- file.path(output_dir, "FigS9_replicate_resolved_safety")
save_pdf(fig_s9, paste0(stem, ".pdf"), 183, 170)
save_png(fig_s9, paste0(stem, ".png"), 183, 170)
save_tiff(fig_s9, paste0(stem, ".tiff"), 183, 170)

# Individual editable PDF panels for downstream assembly.
save_pdf(p_a, file.path(panel_dir, "FigS9a_cell_viability_profiles.pdf"), 183, 92)
save_pdf(p_b, file.path(panel_dir, "FigS9b_haemolysis_profiles.pdf"), 183, 92)

legend_text <- paste0(
  "Supplementary Fig. 9 | Concentration-dependent safety profiles of toxicity-tested Foodins. ",
  "(a) Cell viability and (b) haemolysis measured across five peptide concentrations. ",
  "Points and solid connecting lines show the mean of replicate measurements, and error bars indicate the standard deviation. ",
  "Vertical dashed lines denote 64 μmol L⁻¹, the concentration used for safety prioritization, and horizontal dashed lines indicate the predefined threshold of 0.5. ",
  "Cell-viability values above 1 were capped at 1 for visualization. ",
  "Each condition included three replicate measurements for the 20 Foodins and two replicate measurements for LL37."
)
writeLines(legend_text, file.path(output_dir, "FigS9_legend.txt"))

cat("Fig. S9 completed.\n")
cat("Foodins:", n_distinct(filter(activity_plot_data, peptide != "LL37")$peptide), "\n")
cat("Concentrations:", paste(concentration_levels, collapse = ", "), "μmol L⁻¹\n")
cat("Safety-prioritized at 64 μmol L⁻¹:", paste(format_peptide(safe_foodins), collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n")
