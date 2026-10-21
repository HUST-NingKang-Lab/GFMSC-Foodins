#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(emmeans)
  library(ragg)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop(paste(
    "Usage: Rscript 01_analyze_in_vivo_validation.R",
    "cfu_workbook.xlsx elisa_workbook.xlsx output_dir"
  ))
}

cfu_file <- normalizePath(args[1], mustWork = TRUE)
elisa_file <- normalizePath(args[2], mustWork = TRUE)
out_dir <- args[3]
processed_dir <- file.path(out_dir, "source_data")
dir.create(processed_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Internal group identifiers follow the raw workbooks; display labels use final Foodin names.
group_levels <- c(
  "Infection control", "Peptide 5", "Peptide 6",
  "Peptide 7", "Peptide 8", "Polymyxin B"
)

group_labels <- c(
  "Infection\ncontrol", "Foodin_3", "Foodin_21",
  "Foodin_39", "Foodin_40", "Polymyxin B"
)
names(group_labels) <- group_levels

# Colors are drawn directly from the palette used in the other manuscript figures.
group_colors <- c(
  "Infection control" = "#7C7979",
  "Peptide 5" = "#91CCC0",
  "Peptide 6" = "#7FABD1",
  "Peptide 7" = "#F7AC53",
  "Peptide 8" = "#EC6E66",
  "Polymyxin B" = "#2D8875"
)

theme_nature <- function(base_size = 7, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.text = element_text(colour = "black", size = base_size - 0.4),
      axis.title = element_text(colour = "black", size = base_size),
      plot.title = element_text(face = "bold", size = base_size + 0.5),
      plot.tag = element_text(face = "bold", size = base_size + 1.5),
      plot.tag.position = c(0, 1),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )
}

theme_set(theme_nature())

mean_sd <- function(data, value) {
  data %>%
    group_by(group) %>%
    summarise(
      n = sum(!is.na(.data[[value]])),
      mean = mean(.data[[value]], na.rm = TRUE),
      sd = sd(.data[[value]], na.rm = TRUE),
      .groups = "drop"
    )
}

p_to_symbol <- function(p) {
  case_when(
    p < 0.0001 ~ "****",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

dunnett_table <- function(data, response, endpoint) {
  model <- aov(reformulate("group", response = response), data = data)
  omnibus <- summary(model)[[1]][["Pr(>F)"]][1]
  emm <- emmeans(model, ~ group)
  cmp <- as.data.frame(
    summary(contrast(emm, method = "trt.vs.ctrl", ref = 1, adjust = "dunnett"))
  )

  cmp %>%
    transmute(
      endpoint = endpoint,
      comparison = as.character(contrast),
      estimate = estimate,
      standard_error = SE,
      degrees_freedom = df,
      statistic = t.ratio,
      adjusted_p_value = p.value,
      significance = p_to_symbol(p.value),
      omnibus_anova_p_value = omnibus
    )
}

# -----------------------------------------------------------------------------
# Bacterial burden: average two plated technical replicates for each mouse.
# -----------------------------------------------------------------------------
cfu_raw <- read_excel(cfu_file, col_names = FALSE, .name_repair = "minimal")
first_cfu_row <- which(cfu_raw[[2]] == "单1-1")[1]
cfu_rows <- first_cfu_row + 0:11

cfu_columns <- tibble(
  group = group_levels,
  plate_column = c(2, 22, 26, 30, 34, 38)
)

cfu_plate <- bind_rows(lapply(seq_len(nrow(cfu_columns)), function(i) {
  col <- cfu_columns$plate_column[i]
  tibble(
    group = cfu_columns$group[i],
    mouse = rep(seq_len(6), each = 2),
    technical_replicate = rep(1:2, 6),
    plate_id = as.character(cfu_raw[[col]][cfu_rows]),
    colony_count = as.numeric(cfu_raw[[col + 1]][cfu_rows]),
    dilution_factor = as.numeric(cfu_raw[[col + 2]][cfu_rows]),
    workbook_cfu_g = as.numeric(cfu_raw[[col + 3]][cfu_rows]),
    cfu_g = colony_count * dilution_factor * 100
  )
}))

stopifnot(isTRUE(all.equal(cfu_plate$cfu_g, cfu_plate$workbook_cfu_g)))

cfu_mouse <- cfu_plate %>%
  group_by(group, mouse) %>%
  summarise(
    cfu_g = mean(cfu_g),
    technical_replicates = n(),
    .groups = "drop"
  ) %>%
  mutate(
    group = factor(group, levels = group_levels),
    log10_cfu_g = log10(cfu_g)
  )

cfu_summary <- mean_sd(cfu_mouse, "log10_cfu_g") %>%
  left_join(
    cfu_mouse %>% group_by(group) %>% summarise(mean_cfu_g = mean(cfu_g), .groups = "drop"),
    by = "group"
  ) %>%
  mutate(
    sem = sd / sqrt(n),
    fold_reduction_vs_infection = first(mean_cfu_g) / mean_cfu_g,
    geometric_mean_cfu_g = 10^mean,
    log10_reduction_vs_infection = first(mean) - mean,
    geometric_fold_reduction_vs_infection = 10^log10_reduction_vs_infection
  )

cfu_stats <- dunnett_table(cfu_mouse, "log10_cfu_g", "Muscle bacterial burden")

# -----------------------------------------------------------------------------
# ELISA: each concentration is calculated from duplicate wells for one mouse.
# -----------------------------------------------------------------------------
elisa_group_numbers <- c(1, 6, 7, 8, 9, 10)

elisa <- bind_rows(lapply(c("TNF-α", "IL-6", "IL-10"), function(cytokine) {
  sheet <- read_excel(
    elisa_file,
    sheet = cytokine,
    col_names = FALSE,
    .name_repair = "minimal"
  )
  first_sample_row <- which(sheet[[2]] == "1-1")[1]

  bind_rows(lapply(seq_along(elisa_group_numbers), function(i) {
    original_group <- elisa_group_numbers[i]
    rows <- first_sample_row + (original_group - 1) * 3 + 0:2
    tibble(
      cytokine = cytokine,
      group = group_levels[i],
      mouse = seq_len(3),
      concentration_pg_ml = as.numeric(sheet[[7]][rows])
    )
  }))
})) %>%
  mutate(
    group = factor(group, levels = group_levels),
    cytokine = factor(cytokine, levels = c("TNF-α", "IL-6", "IL-10"))
  )

elisa_summary <- elisa %>%
  group_by(cytokine) %>%
  group_modify(~ mean_sd(.x, "concentration_pg_ml")) %>%
  ungroup() %>%
  group_by(cytokine) %>%
  mutate(percent_reduction_vs_infection = 100 * (1 - mean / first(mean))) %>%
  ungroup()

elisa_stats <- bind_rows(lapply(levels(elisa$cytokine), function(cytokine_name) {
  cytokine_data <- filter(elisa, .data$cytokine == .env$cytokine_name)
  stopifnot(nrow(cytokine_data) == 18)
  dunnett_table(
    cytokine_data,
    "concentration_pg_ml",
    cytokine_name
  )
}))

stopifnot(
  nrow(cfu_mouse) == 36,
  all(cfu_summary$n == 6),
  nrow(elisa) == 54,
  all(table(elisa$cytokine, elisa$group) == 3),
  nrow(elisa_stats) == 15
)

# Save traceable source data before plotting.
write.csv(cfu_plate, file.path(processed_dir, "Figure6_CFU_plate_level_source_data.csv"), row.names = FALSE)
write.csv(cfu_mouse, file.path(processed_dir, "Figure6_CFU_mouse_level_source_data.csv"), row.names = FALSE)
write.csv(cfu_summary, file.path(processed_dir, "Figure6_CFU_summary.csv"), row.names = FALSE)
write.csv(elisa, file.path(processed_dir, "Figure6_ELISA_mouse_level_source_data.csv"), row.names = FALSE)
write.csv(elisa_summary, file.path(processed_dir, "Figure6_ELISA_summary.csv"), row.names = FALSE)
write.csv(bind_rows(cfu_stats, elisa_stats), file.path(processed_dir, "Figure6_statistics.csv"), row.names = FALSE)

# -----------------------------------------------------------------------------
# Panel (a): experimental design. Event spacing is schematic rather than temporal.
# -----------------------------------------------------------------------------
timeline <- tibble(
  x = 1:8,
  time = c("Day -3", "Day -1", "0 h", "2 h", "7 h", "12 h", "14 h", "26 h"),
  event = c(
    "Cyclophosphamide\n150 mg kg-1, i.p.",
    "Cyclophosphamide\n100 mg kg-1, i.p.",
    "E. coli challenge\nright thigh, i.m.",
    "Foodin\n10 mg kg-1, i.p.",
    "Foodin\n10 mg kg-1, i.p.",
    "PBS or polymyxin B\n5 mg kg-1, i.p.",
    "Foodin\n10 mg kg-1, i.p.",
    "Harvest\nmuscle + serum"
  ),
  event_class = c("Pretreatment", "Pretreatment", "Infection", rep("Peptide", 2), "Comparator", "Peptide", "Collection"),
  y_event = c(1.02, 2.10, 1.02, 2.10, 1.02, 2.10, 1.02, 2.10)
)

event_colors <- c(
  "Pretreatment" = "#C7C1DE",
  "Infection" = "#D75B4E",
  "Peptide" = "#91CCC0",
  "Comparator" = "#B5CE4E",
  "Collection" = "#7FABD1"
)

group_key <- tibble(
  x = seq(1.25, 7.75, length.out = 6),
  y = 3.45,
  group = factor(group_levels, levels = group_levels),
  label = c("Infection control", "Foodin_3", "Foodin_21", "Foodin_39", "Foodin_40", "Polymyxin B")
)

p_a <- ggplot() +
  annotate("segment", x = 1, xend = 8, y = 1.58, yend = 1.58,
           linewidth = 0.55, colour = "#4B4B4B",
           arrow = arrow(length = unit(2.0, "mm"), type = "closed")) +
  geom_segment(
    data = timeline,
    aes(x = x, xend = x, y = 1.58, yend = y_event),
    linewidth = 0.32,
    colour = "#777777"
  ) +
  geom_point(
    data = timeline,
    aes(x = x, y = 1.58, fill = event_class),
    shape = 21,
    size = 3.2,
    stroke = 0.45,
    colour = "white"
  ) +
  geom_text(
    data = timeline,
    aes(x = x, y = y_event, label = event),
    size = 2.25,
    lineheight = 0.92,
    vjust = ifelse(timeline$y_event > 1.58, -0.10, 1.10)
  ) +
  geom_text(
    data = timeline,
    aes(x = x, y = 1.40, label = time),
    size = 2.25,
    fontface = "bold",
    vjust = 1
  ) +
  geom_tile(
    data = group_key,
    aes(x = x, y = y, fill = group),
    width = 0.94,
    height = 0.40,
    colour = "white",
    linewidth = 0.35
  ) +
  geom_text(
    data = group_key,
    aes(x = x, y = y, label = label),
    size = 2.15,
    fontface = "bold",
    colour = "#202020"
  ) +
  annotate("text", x = 0.95, y = 3.88, label = "Experimental groups", hjust = 0,
           size = 2.55, fontface = "bold") +
  scale_fill_manual(values = c(event_colors, group_colors)) +
  coord_cartesian(xlim = c(0.82, 8.18), ylim = c(0.44, 4.02), clip = "off") +
  labs(tag = "(a)") +
  theme_void(base_family = "Helvetica") +
  theme(
    legend.position = "none",
    plot.tag = element_text(face = "bold", size = 8.5),
    plot.tag.position = c(0, 1),
    plot.margin = margin(5, 7, 0, 7)
  )

# -----------------------------------------------------------------------------
# Panel (b): bacterial burden, the principal efficacy endpoint.
# -----------------------------------------------------------------------------
cfu_annotation <- cfu_stats %>%
  mutate(
    group = factor(sub(" - Infection control$", "", comparison), levels = group_levels),
    y = max(cfu_mouse$log10_cfu_g) + 0.28
  )
cfu_y_min <- floor((min(cfu_mouse$log10_cfu_g) - 0.15) * 2) / 2
cfu_y_max <- max(cfu_mouse$log10_cfu_g) + 0.60
cfu_y_breaks <- seq(ceiling(cfu_y_min), floor(cfu_y_max), by = 1)

p_b <- ggplot(cfu_mouse, aes(x = group, y = log10_cfu_g, colour = group)) +
  geom_point(
    position = position_jitter(width = 0.10, height = 0, seed = 24),
    size = 2.25,
    alpha = 0.86
  ) +
  geom_errorbar(
    data = cfu_summary,
    aes(x = group, ymin = mean - sem, ymax = mean + sem),
    inherit.aes = FALSE,
    width = 0.14,
    linewidth = 0.50,
    colour = "black"
  ) +
  geom_point(
    data = cfu_summary,
    aes(x = group, y = mean),
    inherit.aes = FALSE,
    shape = 18,
    size = 3.0,
    colour = "black"
  ) +
  geom_text(
    data = cfu_annotation,
    aes(x = group, y = y, label = significance),
    inherit.aes = FALSE,
    size = 2.65,
    fontface = "bold"
  ) +
  scale_colour_manual(values = group_colors) +
  scale_x_discrete(labels = group_labels, expand = expansion(add = 0.45)) +
  scale_y_continuous(
    limits = c(cfu_y_min, cfu_y_max),
    breaks = cfu_y_breaks,
    expand = expansion(mult = c(0, 0))
  ) +
  labs(
    tag = "(b)",
    title = "Muscle bacterial burden",
    x = NULL,
    y = expression(log[10]~"CFU g"^{-1}~"muscle")
  ) +
  theme_nature(base_size = 7.1) +
  theme(
    axis.text.x = element_text(angle = 26, hjust = 1, vjust = 1, size = 6.3),
    plot.title = element_text(margin = margin(b = 4))
  )

# -----------------------------------------------------------------------------
# Panels (c-e): serum cytokine responses, shown with individual mice and mean ± SD.
# -----------------------------------------------------------------------------
cytokine_plot <- function(cytokine, tag, title, x_breaks = waiver()) {
  cytokine_name <- cytokine
  data <- filter(elisa, .data$cytokine == .env$cytokine_name)
  summary_data <- filter(elisa_summary, .data$cytokine == .env$cytokine_name)
  stat_data <- filter(elisa_stats, .data$endpoint == .env$cytokine_name) %>%
    mutate(group = factor(sub(" - Infection control$", "", comparison), levels = group_levels))

  upper <- max(summary_data$mean + summary_data$sd) * 1.28
  stat_data <- stat_data %>% mutate(x = upper * 0.985)

  ggplot(data, aes(x = concentration_pg_ml, y = group, colour = group)) +
    geom_errorbar(
      data = summary_data,
      aes(y = group, xmin = mean - sd, xmax = mean + sd),
      inherit.aes = FALSE,
      orientation = "y",
      width = 0.18,
      linewidth = 0.48,
      colour = "black"
    ) +
    geom_point(
      position = position_jitter(width = 0, height = 0.075, seed = 24),
      size = 1.72,
      alpha = 0.82
    ) +
    geom_point(
      data = summary_data,
      aes(x = mean, y = group),
      inherit.aes = FALSE,
      shape = 18,
      size = 2.25,
      colour = "black"
    ) +
    geom_text(
      data = stat_data,
      aes(x = x, y = group, label = significance),
      inherit.aes = FALSE,
      hjust = 1,
      size = 2.35,
      fontface = "bold"
    ) +
    scale_colour_manual(values = group_colors) +
    scale_y_discrete(labels = group_labels, limits = rev(group_levels)) +
    scale_x_continuous(
      limits = c(0, upper),
      breaks = x_breaks,
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      tag = tag,
      title = title,
    x = expression("Serum concentration (pg ml"^{-1}*")"),
      y = NULL
    ) +
    theme_nature(base_size = 6.7) +
    theme(
      axis.text.y = element_text(size = 5.8),
      axis.text.x = element_text(size = 5.9),
      axis.title.x = element_text(size = 6.3),
      plot.title = element_text(margin = margin(b = 2)),
      plot.margin = margin(4, 5, 5, 10)
    )
}

p_c <- cytokine_plot("TNF-α", "(c)", expression("TNF-"*alpha), seq(0, 1200, 300))
p_d <- cytokine_plot("IL-6", "(d)", "IL-6", seq(0, 160, 40))
p_e <- cytokine_plot("IL-10", "(e)", "IL-10", seq(0, 900, 300))

right_column <- p_c / p_d / p_e + plot_layout(heights = c(1, 1, 1))
figure6 <- p_a / (p_b | right_column) +
  plot_layout(heights = c(0.75, 2.25), widths = c(1.03, 1))

base_name <- file.path(out_dir, "Figure6_in_vivo_validation")
width_in <- 183 / 25.4
height_in <- 165 / 25.4

ggsave(
  filename = paste0(base_name, ".pdf"),
  plot = figure6,
  device = "pdf",
  width = width_in,
  height = height_in,
  units = "in",
  family = "Helvetica",
  useDingbats = FALSE
)

ragg::agg_tiff(
  paste0(base_name, ".tiff"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 600,
  compression = "lzw"
)
print(figure6)
dev.off()

ragg::agg_png(
  paste0(base_name, ".png"),
  width = width_in,
  height = height_in,
  units = "in",
  res = 300
)
print(figure6)
dev.off()

# Export each unassembled panel separately for downstream editing.
panel_dir <- file.path(out_dir, "individual_panels")
dir.create(panel_dir, showWarnings = FALSE, recursive = TRUE)

save_individual_panel <- function(plot, filename, width_mm, height_mm) {
  width <- width_mm / 25.4
  height <- height_mm / 25.4
  standalone_plot <- plot + labs(tag = NULL)

  ggsave(
    filename = file.path(panel_dir, paste0(filename, ".pdf")),
    plot = standalone_plot,
    device = "pdf",
    width = width,
    height = height,
    units = "in",
    family = "Helvetica",
    useDingbats = FALSE
  )

  ragg::agg_tiff(
    file.path(panel_dir, paste0(filename, ".tiff")),
    width = width,
    height = height,
    units = "in",
    res = 600,
    compression = "lzw"
  )
  print(standalone_plot)
  dev.off()

  ragg::agg_png(
    file.path(panel_dir, paste0(filename, ".png")),
    width = width,
    height = height,
    units = "in",
    res = 300
  )
  print(standalone_plot)
  dev.off()
}

save_individual_panel(p_a, "Figure6a_experimental_design", 183, 58)
save_individual_panel(p_b, "Figure6b_muscle_bacterial_burden", 90, 115)
save_individual_panel(p_c, "Figure6c_TNF_alpha", 90, 58)
save_individual_panel(p_d, "Figure6d_IL6", 90, 58)
save_individual_panel(p_e, "Figure6e_IL10", 90, 58)

legend_text <- paste0(
  "Fig. 6 | Prioritized Foodin candidates reduce bacterial burden and infection-associated cytokine responses ",
  "in a murine deep-thigh infection model. (a) Experimental design. Neutropenia was induced in female ICR mice ",
  "by intraperitoneal administration of cyclophosphamide 3 and 1 days before intramuscular challenge of the right ",
  "thigh with Escherichia coli. Peptides were administered intraperitoneally at 10 mg kg−1 at 2, 7 and 14 h after ",
  "infection. The infection-control group received PBS, and the antibiotic-control group received polymyxin B ",
  "(5 mg kg−1) at 12 h after infection. Thigh muscle and serum were collected at 26 h. (b) Bacterial burden in ",
  "infected thigh muscle. Each coloured point represents one mouse and is the mean of two technical plating ",
  "replicates (n = 6 mice per group); Foodin_3, Foodin_21, Foodin_39 and Foodin_40 correspond to experimental peptide ",
  "groups 5, 6, 7 and 8, respectively. Black diamonds and error bars indicate the mean +/- s.e.m. of log10-transformed ",
  "values. (c-e) Serum concentrations of TNF-alpha (c), IL-6 (d) and IL-10 (e), measured by ELISA. Each point ",
  "represents one mouse and was calculated from duplicate wells (n = 3 mice per group); black diamonds and error ",
  "bars indicate the mean +/- s.d. Statistical significance was assessed by one-way ANOVA followed by Dunnett's ",
  "multiple-comparisons test against the infection-control group. *P < 0.05, **P < 0.01, ***P < 0.001 and ",
  "****P < 0.0001."
)
writeLines(legend_text, file.path(out_dir, "Figure6_legend.txt"))

cat("Figure written to:", out_dir, "\n")
cat("Individual panels written to:", panel_dir, "\n")
cat("Processed data written to:", processed_dir, "\n")
cat("CFU biological replicates per group:", paste(cfu_summary$n, collapse = ", "), "\n")
cat("ELISA biological replicates per group: 3\n")
