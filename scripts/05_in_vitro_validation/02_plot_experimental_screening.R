#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(scales)
  library(patchwork)
  library(ggrepel)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) {
  stop(paste(
    "Usage: Rscript 02_plot_experimental_screening.R",
    "mic_calls.csv cell_viability.csv haemolysis.csv output_dir"
  ))
}

mic_file <- normalizePath(args[1], mustWork = TRUE)
activity_file <- normalizePath(args[2], mustWork = TRUE)
haemolysis_file <- normalizePath(args[3], mustWork = TRUE)
output_dir <- args[4]
panel_dir <- file.path(output_dir, "individual_panels")
source_dir <- file.path(output_dir, "source_data")
dir.create(panel_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(source_dir, recursive = TRUE, showWarnings = FALSE)

stopifnot(file.exists(mic_file), file.exists(activity_file),
          file.exists(haemolysis_file))

pal <- c(
  teal = "#2D8875", mint = "#91CCC0", blue = "#7FABD1",
  orange = "#F7AC53", coral = "#EC6E66", lilac = "#C7C1DE",
  grey = "#7C7979", light_grey = "#E8E8E8", black = "#252525"
)

theme_supp <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", size = base_size + 0.4),
      plot.tag = element_text(face = "bold", size = base_size + 1.3),
      strip.background = element_blank(),
      panel.grid = element_blank(),
      legend.title = element_text(size = base_size - 0.2),
      legend.text = element_text(size = base_size - 0.5),
      plot.margin = margin(5, 7, 5, 5)
    )
}

mic <- read_csv(mic_file, show_col_types = FALSE) %>%
  filter(peptide != "LL37") %>%
  mutate(active = is.finite(MIC) & MIC <= 64)

activity <- read_csv(activity_file, show_col_types = FALSE)
haemolysis <- read_csv(haemolysis_file, show_col_types = FALSE)

# Determine the six candidates meeting both pre-defined thresholds at 64 μmol L⁻¹.
safety_64 <- activity %>%
  filter(peptide != "LL37", concentration == 64) %>%
  group_by(peptide) %>%
  summarise(cell_viability = mean(activity, na.rm = TRUE), .groups = "drop") %>%
  inner_join(
    haemolysis %>%
      filter(peptide != "LL37", concentration == 64) %>%
      group_by(peptide) %>%
      summarise(haemolysis = mean(Haemolytic, na.rm = TRUE), .groups = "drop"),
    by = "peptide"
  ) %>%
  filter(cell_viability >= 0.5, haemolysis < 0.5)

safe_foodins <- safety_64$peptide
in_vivo_foodins <- c("Food3", "Food21", "Food39", "Food40")
toxicity_foodins <- setdiff(unique(activity$peptide), "LL37")

breadth <- mic %>%
  group_by(peptide) %>%
  summarise(
    pathogens_inhibited = sum(active),
    median_mic = ifelse(any(active), median(MIC[active]), NA_real_),
    .groups = "drop"
  ) %>%
  mutate(
    display_name = str_replace(peptide, "^Food", "Foodin_"),
    selection = case_when(
      peptide %in% in_vivo_foodins ~ "In vivo selected",
      peptide %in% setdiff(safe_foodins, in_vivo_foodins) ~
        "Additional in vitro prioritized",
      peptide %in% toxicity_foodins ~ "Other toxicity-tested",
      TRUE ~ "Other synthesized Foodins"
    ),
    selection = factor(
      selection,
      levels = c("Other synthesized Foodins", "Other toxicity-tested",
                 "Additional in vitro prioritized", "In vivo selected")
    )
  )

# Integrity checks against the manuscript claims.
stopifnot(
  nrow(breadth) == 48,
  sum(breadth$pathogens_inhibited >= 1) == 42,
  sum(breadth$pathogens_inhibited >= 7) == 26,
  sum(breadth$pathogens_inhibited == 8) == 15,
  length(toxicity_foodins) == 20,
  length(safe_foodins) == 6,
  all(in_vivo_foodins %in% safe_foodins)
)

breadth_counts <- breadth %>%
  count(pathogens_inhibited, name = "n_foodins") %>%
  complete(pathogens_inhibited = 0:8, fill = list(n_foodins = 0)) %>%
  mutate(
    activity_class = case_when(
      pathogens_inhibited == 0 ~ "Inactive",
      pathogens_inhibited >= 7 ~ "Broad activity",
      TRUE ~ "Activity against 1-6 pathogens"
    ),
    activity_class = factor(
      activity_class,
      levels = c("Inactive", "Activity against 1-6 pathogens", "Broad activity")
    )
  )

pathogen_labels <- c(
  "A. baumannii ATCC19606 (-)" = "A. baumannii ATCC 19606",
  "E. coli ATCC-25922 (-)" = "E. coli ATCC 25922",
  "E. coli BL21 Ampicillin-resistance (-)" = "E. coli BL21 AmpR",
  "E. faecium ATCC 700221-VRE (+)" = "E. faecium ATCC 700221-VRE",
  "K. pneumoniae ATCC 13883(-)" = "K. pneumoniae ATCC 13883",
  "P. aeruginosa ATCC 27853 (-)" = "P. aeruginosa ATCC 27853",
  "P. aeruginosa PAO1 (-)" = "P. aeruginosa PAO1",
  "S. aureus ATCC-25923 (+)" = "S. aureus ATCC 25923"
)

mic_levels <- c("4", "8", "16", "32", "64", ">64")
mic_distribution <- mic %>%
  mutate(
    pathogen = recode(bacteria, !!!pathogen_labels),
    mic_class = case_when(
      !is.finite(MIC) | MIC > 64 ~ ">64",
      MIC <= 4 ~ "4", MIC <= 8 ~ "8", MIC <= 16 ~ "16",
      MIC <= 32 ~ "32", TRUE ~ "64"
    ),
    mic_class = factor(mic_class, levels = mic_levels),
    pathogen = factor(pathogen, levels = unname(pathogen_labels))
  ) %>%
  count(pathogen, mic_class, name = "n") %>%
  complete(pathogen, mic_class, fill = list(n = 0)) %>%
  group_by(pathogen) %>%
  mutate(proportion = n / sum(n)) %>%
  ungroup()

funnel <- tibble(
  stage = c("Synthesized Foodins", "Toxicity tested",
            "In vitro prioritized", "In vivo selected"),
  n = c(48, length(toxicity_foodins), length(safe_foodins),
        length(in_vivo_foodins))
) %>%
  mutate(
    stage = factor(stage, levels = rev(stage)),
    fill_group = factor(as.character(stage), levels = levels(stage))
  )

write_csv(breadth, file.path(source_dir, "FigS8a_c_Foodin_screening_summary.csv"))
write_csv(mic_distribution,
          file.path(source_dir, "FigS8b_MIC_distribution_by_pathogen.csv"))
write_csv(safety_64, file.path(source_dir, "FigS8_safety_classification_64micromol_per_L.csv"))
write_csv(funnel, file.path(source_dir, "FigS8d_selection_stages.csv"))

p_a <- ggplot(breadth_counts,
              aes(pathogens_inhibited, n_foodins, fill = activity_class)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = ifelse(n_foodins > 0, n_foodins, "")),
            vjust = -0.35, size = 2.4) +
  scale_fill_manual(values = c(
    "Inactive" = pal[["grey"]],
    "Activity against 1-6 pathogens" = pal[["mint"]],
    "Broad activity" = pal[["coral"]]
  )) +
  scale_x_continuous(breaks = 0:8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(title = "Antibacterial spectrum", x = "Pathogens inhibited",
       y = "Number of Foodins", fill = NULL) +
  theme_supp() +
  theme(legend.position = "bottom")

p_b <- ggplot(mic_distribution,
              aes(pathogen, proportion, fill = mic_class)) +
  geom_col(width = 0.72, colour = "black", linewidth = 0.15) +
  coord_flip() +
  scale_y_continuous(labels = percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0))) +
  scale_fill_manual(
    values = c("4" = "#075A52", "8" = pal[["teal"]],
               "16" = pal[["mint"]], "32" = pal[["blue"]],
               "64" = pal[["lilac"]], ">64" = pal[["light_grey"]]),
    drop = FALSE
  ) +
  labs(title = "MIC distributions across pathogens", x = NULL,
       y = "Proportion of Foodins",
       fill = expression(MIC~(mu*mol~L^{-1}))) +
  theme_supp() +
  theme(legend.position = "bottom", legend.key.width = unit(5, "mm"))

selection_colours <- c(
  "Other synthesized Foodins" = pal[["light_grey"]],
  "Other toxicity-tested" = pal[["teal"]],
  "Additional in vitro prioritized" = pal[["orange"]],
  "In vivo selected" = pal[["coral"]]
)

p_c <- ggplot(filter(breadth, pathogens_inhibited > 0),
              aes(pathogens_inhibited, median_mic, colour = selection,
                  shape = selection)) +
  geom_point(size = 2.1, alpha = 0.9,
             position = position_jitter(width = 0.10, height = 0,
                                        seed = 20260903)) +
  geom_text_repel(
    data = filter(breadth, peptide %in% safe_foodins),
    aes(label = display_name), size = 2.1, show.legend = FALSE,
    seed = 20260903, min.segment.length = 0,
    box.padding = 0.25, max.overlaps = Inf
  ) +
  scale_colour_manual(values = selection_colours, drop = FALSE) +
  scale_shape_manual(values = c(16, 17, 15, 18), drop = FALSE) +
  scale_x_continuous(breaks = 1:8, limits = c(0.7, 8.5)) +
  scale_y_continuous(trans = "log2", breaks = c(4, 8, 16, 32, 64),
                     limits = c(3.6, 80), expand = expansion(mult = c(0, 0))) +
  labs(title = "Breadth-potency landscape", x = "Pathogens inhibited",
       y = expression(Median~MIC~(mu*mol~L^{-1})), colour = NULL,
       shape = NULL) +
  guides(
    colour = guide_legend(nrow = 2, byrow = TRUE),
    shape = guide_legend(nrow = 2, byrow = TRUE)
  ) +
  theme_supp() +
  theme(
    legend.position = "bottom",
    legend.justification = "left",
    legend.box.just = "left",
    legend.margin = margin(t = 2),
    legend.spacing.x = unit(1.5, "mm")
  )

p_d <- ggplot(funnel, aes(stage, n, fill = fill_group)) +
  geom_col(width = 0.68, colour = "black", linewidth = 0.25) +
  geom_text(aes(label = n), hjust = -0.25, size = 2.7) +
  coord_flip() +
  scale_fill_manual(values = c(
    "In vivo selected" = pal[["coral"]],
    "In vitro prioritized" = pal[["orange"]],
    "Toxicity tested" = pal[["blue"]],
    "Synthesized Foodins" = pal[["teal"]]
  )) +
  scale_y_continuous(limits = c(0, 54), expand = expansion(mult = c(0, 0))) +
  labs(title = "Experimental prioritization", x = NULL,
       y = "Number of Foodins") +
  theme_supp() +
  theme(legend.position = "none")

fig_s8 <- (p_a | p_b) / (p_c | p_d) +
  plot_layout(widths = c(1, 1.35), heights = c(1, 1)) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")

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
    bg = "white"
  )
}

save_png <- function(plot, path, width_mm, height_mm) {
  agg_png(path, width = width_mm, height = height_mm, units = "mm",
          res = 600, background = "white")
  print(plot)
  dev.off()
}

save_tiff <- function(plot, path, width_mm, height_mm) {
  agg_tiff(path, width = width_mm, height = height_mm, units = "mm",
           res = 600, compression = "lzw", background = "white")
  print(plot)
  dev.off()
}

stem <- file.path(output_dir, "FigS8_experimental_screening")
save_pdf(fig_s8, paste0(stem, ".pdf"), 183, 142)
save_png(fig_s8, paste0(stem, ".png"), 183, 142)
save_tiff(fig_s8, paste0(stem, ".tiff"), 183, 142)

# Vector PDF exports of individual panels for downstream editing.
save_pdf(p_a, file.path(panel_dir, "FigS8a_antibacterial_spectrum.pdf"), 82, 65)
save_pdf(p_b, file.path(panel_dir, "FigS8b_MIC_distribution.pdf"), 101, 65)
save_pdf(p_c, file.path(panel_dir, "FigS8c_breadth_potency.pdf"), 82, 68)
save_pdf(p_d, file.path(panel_dir, "FigS8d_selection_funnel.pdf"), 101, 68)

legend_text <- paste0(
  "Supplementary Fig. 8 | Experimental screening and prioritization of Foodins. ",
  "(a) Distribution of the number of pathogens inhibited by each of the 48 synthesized Foodins. ",
  "A Foodin was considered active when its MIC was no greater than 64 μmol L⁻¹. ",
  "(b) MIC distributions of the 48 Foodins across the eight tested pathogens. ",
  "(c) Relationship between antibacterial breadth and median MIC among active Foodins. ",
  "The six candidates meeting the predefined cell-viability and haemolysis thresholds at 64 μmol L⁻¹ are labelled. ",
  "(d) Stepwise experimental prioritization from 48 synthesized Foodins to 20 toxicity-tested candidates, ",
  "six in vitro-prioritized candidates and four candidates selected for in vivo evaluation."
)
writeLines(legend_text, file.path(output_dir, "FigS8_legend.txt"))

cat("Fig. S8 completed.\n")
cat("Active against >=1 pathogen:", sum(breadth$pathogens_inhibited >= 1), "\n")
cat("Active against >=7 pathogens:", sum(breadth$pathogens_inhibited >= 7), "\n")
cat("Active against all 8 pathogens:", sum(breadth$pathogens_inhibited == 8), "\n")
cat("Safety-prioritized at 64 μmol L⁻¹:", paste(safe_foodins, collapse = ", "), "\n")
cat("Output directory:", output_dir, "\n")
