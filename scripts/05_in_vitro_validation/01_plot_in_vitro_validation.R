#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
  library(scales)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
  stop(paste(
    "Usage: Rscript 01_plot_in_vitro_validation.R",
    "screen_1_30.csv screen_31_50.csv dilution_1.csv dilution_2.csv",
    "cell_viability.csv haemolysis.csv output_dir"
  ))
}

screen_1_file <- normalizePath(args[1], mustWork = TRUE)
screen_2_file <- normalizePath(args[2], mustWork = TRUE)
dilution_1_file <- normalizePath(args[3], mustWork = TRUE)
dilution_2_file <- normalizePath(args[4], mustWork = TRUE)
activity_file <- normalizePath(args[5], mustWork = TRUE)
hem_file <- normalizePath(args[6], mustWork = TRUE)
out_dir <- args[7]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

palette_contract <- c(
  neutral_dark = "#262626",
  neutral_mid = "#6F6F6F",
  neutral_light = "#E8E8E8",
  signal_deep = "#075A52",
  signal_teal = "#1E9C89",
  signal_soft = "#A7D3C7",
  accent_red = "#C84A3D",
  accent_orange = "#D88A32",
  soft_yellow = "#E7D98E"
)

theme_nature_contract <- function(base_size = 6.5, base_family = "Helvetica") {
  theme_classic(base_size = base_size, base_family = base_family) +
    theme(
      axis.line = element_line(linewidth = 0.35, colour = "black"),
      axis.ticks = element_line(linewidth = 0.35, colour = "black"),
      axis.title = element_text(size = base_size),
      axis.text = element_text(size = base_size - 0.6, colour = "black"),
      legend.title = element_text(size = base_size - 0.3),
      legend.text = element_text(size = base_size - 0.8),
      strip.text = element_text(size = base_size - 0.2, face = "bold"),
      plot.title = element_text(size = base_size + 0.5, face = "bold"),
      panel.grid = element_blank()
    )
}

theme_set(theme_nature_contract())

peptide_order_all <- c("LL37", paste0("Food", 1:50))

read_mic64 <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)
  blank <- raw %>%
    filter(bacteria == "blank") %>%
    pivot_longer(-c(bacteria, group, control), names_to = "peptide", values_to = "blank_od") %>%
    group_by(peptide) %>%
    summarise(blank_mean = mean(blank_od, na.rm = TRUE), .groups = "drop")

  blank_control_mean <- raw %>%
    filter(bacteria == "blank") %>%
    summarise(blank_control_mean = mean(control, na.rm = TRUE)) %>%
    pull(blank_control_mean)

  control <- raw %>%
    filter(bacteria != "blank") %>%
    group_by(bacteria) %>%
    summarise(control_mean = mean(control, na.rm = TRUE), .groups = "drop")

  raw %>%
    filter(bacteria != "blank") %>%
    pivot_longer(-c(bacteria, group, control), names_to = "peptide", values_to = "od") %>%
    left_join(blank, by = "peptide") %>%
    left_join(control, by = "bacteria") %>%
    mutate(
      concentration = 64,
      corrected_sample = pmax(od - blank_mean, 0),
      corrected_control = pmax(control_mean - blank_control_mean, 1e-8),
      inhibition = pmax(0, pmin(1, 1 - corrected_sample / corrected_control)),
      assay = "64 μmol L⁻¹ screen"
    )
}

read_mic_detail <- function(path) {
  raw <- read_csv(path, show_col_types = FALSE)
  blank <- raw %>%
    filter(bacteria == "blank") %>%
    pivot_longer(-c(bacteria, group, control), names_to = "name", values_to = "blank_od") %>%
    mutate(
      peptide = str_replace(name, "_[0-9]+$", ""),
      concentration = as.numeric(str_extract(name, "[0-9]+$"))
    ) %>%
    group_by(peptide, concentration) %>%
    summarise(blank_mean = mean(blank_od, na.rm = TRUE), .groups = "drop")

  blank_control_mean <- raw %>%
    filter(bacteria == "blank") %>%
    summarise(blank_control_mean = mean(control, na.rm = TRUE)) %>%
    pull(blank_control_mean)

  control <- raw %>%
    filter(bacteria != "blank") %>%
    group_by(bacteria) %>%
    summarise(control_mean = mean(control, na.rm = TRUE), .groups = "drop")

  raw %>%
    filter(bacteria != "blank") %>%
    pivot_longer(-c(bacteria, group, control), names_to = "name", values_to = "od") %>%
    mutate(
      peptide = str_replace(name, "_[0-9]+$", ""),
      concentration = as.numeric(str_extract(name, "[0-9]+$"))
    ) %>%
    left_join(blank, by = c("peptide", "concentration")) %>%
    left_join(control, by = "bacteria") %>%
    mutate(
      corrected_sample = pmax(od - blank_mean, 0),
      corrected_control = pmax(control_mean - blank_control_mean, 1e-8),
      inhibition = pmax(0, pmin(1, 1 - corrected_sample / corrected_control)),
      assay = "MIC refinement"
    ) %>%
    select(-name)
}

mic_raw <- bind_rows(
  read_mic64(screen_1_file),
  read_mic64(screen_2_file),
  read_mic_detail(dilution_1_file),
  read_mic_detail(dilution_2_file)
)

mic_mean <- mic_raw %>%
  group_by(bacteria, peptide, concentration) %>%
  summarise(
    inhibition_mean = mean(inhibition, na.rm = TRUE),
    inhibition_sd = sd(inhibition, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

mic_call <- mic_mean %>%
  group_by(bacteria, peptide) %>%
  summarise(
    MIC = ifelse(any(inhibition_mean >= 0.5), min(concentration[inhibition_mean >= 0.5]), Inf),
    max_inhibition = max(inhibition_mean, na.rm = TRUE),
    .groups = "drop"
  )

peptide_summary <- mic_call %>%
  group_by(peptide) %>%
  summarise(
    n_inhibited = sum(is.finite(MIC)),
    median_mic = ifelse(n_inhibited > 0, median(MIC[is.finite(MIC)]), Inf),
    min_mic = min(MIC),
    .groups = "drop"
  )

activity <- read.csv(activity_file, check.names = FALSE) %>%
  group_by(peptide, concentration) %>%
  summarise(
    activity_mean_raw = mean(activity, na.rm = TRUE),
    activity_sd = sd(activity, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    activity_mean = pmin(activity_mean_raw, 1)
  )

hemolysis <- read.csv(hem_file, check.names = FALSE) %>%
  group_by(peptide, concentration) %>%
  summarise(
    hem_mean = mean(Haemolytic, na.rm = TRUE),
    hem_sd = sd(Haemolytic, na.rm = TRUE),
    .groups = "drop"
  )

safety_64 <- activity %>%
  filter(concentration == 64) %>%
  inner_join(hemolysis %>% filter(concentration == 64), by = c("peptide", "concentration")) %>%
  left_join(peptide_summary, by = "peptide") %>%
  mutate(
    n_inhibited = ifelse(is.na(n_inhibited), 0, n_inhibited),
    safety_class = case_when(
      activity_mean >= 0.5 & hem_mean <= 0.5 ~ "Pass both thresholds",
      activity_mean >= 0.5 & hem_mean > 0.5 ~ "Viable but haemolytic",
      activity_mean < 0.5 & hem_mean <= 0.5 ~ "Low viability",
      TRUE ~ "Fail both thresholds"
    )
  )

tox_heatmap <- bind_rows(
  activity %>%
    transmute(peptide, concentration, endpoint = "Cell viability", value = activity_mean),
  hemolysis %>%
    transmute(peptide, concentration, endpoint = "Haemolysis", value = hem_mean)
) %>%
  mutate(
    peptide = factor(peptide, levels = rev(intersect(peptide_order_all, unique(peptide)))),
    concentration = factor(concentration, levels = c(8, 16, 32, 64, 128)),
    value_text = paste0(round(value * 100), "%"),
    text_colour = ifelse(value >= 0.72, "white", "black")
  )

bacteria_order <- c(
  "A. baumannii ATCC19606 (-)",
  "E. coli ATCC-25922 (-)",
  "E. coli BL21 Ampicillin-resistance (-)",
  "E. faecium ATCC 700221-VRE (+)",
  "K. pneumoniae ATCC 13883(-)",
  "P. aeruginosa ATCC 27853 (-)",
  "P. aeruginosa PAO1 (-)",
  "S. aureus ATCC-25923 (+)"
)

bacteria_labels <- c(
  "A. baumannii\nATCC 19606",
  "E. coli\nATCC 25922",
  "E. coli\nBL21 AmpR",
  "E. faecium\nVRE",
  "K. pneumoniae\nATCC 13883",
  "P. aeruginosa\nATCC 27853",
  "P. aeruginosa\nPAO1",
  "S. aureus\nATCC 25923"
)
names(bacteria_labels) <- bacteria_order

mic_order <- intersect(peptide_order_all, unique(mic_call$peptide))
mic_plot_data <- mic_call %>%
  mutate(
    peptide = factor(peptide, levels = rev(mic_order)),
    bacteria = factor(bacteria, levels = bacteria_order),
    MIC_label = case_when(
      MIC == 4 ~ "4",
      MIC == 8 ~ "8",
      MIC == 16 ~ "16",
      MIC == 32 ~ "32",
      MIC == 64 ~ "64",
      TRUE ~ ">64"
    ),
    MIC_label = factor(MIC_label, levels = c("4", "8", "16", "32", "64", ">64")),
    MIC_value = ifelse(is.finite(MIC), MIC, NA_real_),
    MIC_text = ifelse(MIC_label == ">64", "", as.character(MIC_label))
  )

p_a <- ggplot(mic_plot_data, aes(x = bacteria, y = peptide, fill = MIC_value)) +
  geom_tile(width = 0.92, height = 0.92, colour = "white", linewidth = 0.12) +
  geom_text(aes(label = MIC_text), size = 1.55, colour = "white", fontface = "bold") +
  scale_fill_gradientn(
    colours = c("#075A52", "#1E8177", "#58A99B", "#A7D3C7", "#E7D98E"),
    values = rescale(c(4, 8, 16, 32, 64), to = c(0, 1)),
    limits = c(4, 64),
    breaks = c(4, 8, 16, 32, 64),
    na.value = "#ECECEC",
    name = expression(atop(MIC, (mu*mol~L^{-1}))),
    guide = guide_colorbar(
      barheight = grid::unit(115, "mm"),
      barwidth = grid::unit(3.0, "mm"),
      ticks.linewidth = 0.25,
      frame.colour = "black",
      frame.linewidth = 0.25
    )
  ) +
  scale_x_discrete(labels = bacteria_labels, expand = expansion(mult = c(0, 0))) +
  scale_y_discrete(expand = expansion(mult = c(0, 0))) +
  labs(
    tag = "b",
    title = "Blank-corrected MIC across eight pathogens",
    x = NULL,
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 5.3),
    axis.text.y = element_text(size = 4.4),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    plot.title = element_text(margin = margin(b = 3))
  )

safety_label_data <- safety_64 %>%
  filter(activity_mean >= 0.5, hem_mean <= 0.5) %>%
  arrange(desc(n_inhibited), desc(activity_mean), hem_mean) %>%
  slice_head(n = 8)

p_b <- ggplot(safety_64, aes(x = hem_mean * 100, y = activity_mean * 100)) +
  annotate("rect", xmin = -Inf, xmax = 50, ymin = 50, ymax = Inf,
           fill = unname(palette_contract["signal_soft"]), alpha = 0.40) +
  geom_hline(yintercept = 50, linetype = "dashed", linewidth = 0.35, colour = palette_contract["neutral_mid"]) +
  geom_vline(xintercept = 50, linetype = "dashed", linewidth = 0.35, colour = palette_contract["neutral_mid"]) +
  geom_point(aes(fill = safety_class),
             shape = 21, size = 3.2, colour = "white", stroke = 0.45, alpha = 0.95) +
  geom_text(
    data = safety_label_data,
    aes(label = peptide),
    size = 1.95,
    vjust = -0.85,
    colour = palette_contract["neutral_dark"],
    check_overlap = TRUE
  ) +
  scale_fill_manual(
    values = c(
      "Pass both thresholds" = unname(palette_contract["signal_teal"]),
      "Viable but haemolytic" = unname(palette_contract["soft_yellow"]),
      "Low viability" = unname(palette_contract["accent_orange"]),
      "Fail both thresholds" = unname(palette_contract["accent_red"])
    ),
    name = NULL
  ) +
  scale_x_continuous(breaks = c(0, 25, 50, 75), labels = label_number(suffix = "%")) +
  scale_y_continuous(breaks = c(0, 50, 100), labels = label_number(suffix = "%")) +
  coord_cartesian(xlim = c(0, 80), ylim = c(0, 105), clip = "off") +
  labs(
    tag = "c",
    title = expression("Safety window at 64 "*mu*"mol L"^{-1}),
    x = "Haemolysis",
    y = "Cell viability"
  ) +
  guides(fill = "none") +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.title = element_text(margin = margin(b = 3))
  )

p_c <- ggplot(tox_heatmap, aes(x = concentration, y = peptide, fill = value)) +
  geom_tile(width = 0.92, height = 0.92, colour = "white", linewidth = 0.16) +
  geom_text(aes(label = value_text, colour = text_colour), size = 1.45, fontface = "bold") +
  facet_wrap(~endpoint, ncol = 2) +
  scale_fill_gradientn(
    colours = c("#F2F2F2", "#BFDCD4", "#1E9C89", "#075A52"),
    limits = c(0, 1),
    oob = squish,
    breaks = c(0, 0.5, 1),
    labels = c("0%", "50%", "100%"),
    name = "Mean\nvalue",
    guide = guide_colorbar(
      barheight = grid::unit(62, "mm"),
      barwidth = grid::unit(3.0, "mm"),
      ticks.linewidth = 0.25,
      frame.colour = "black",
      frame.linewidth = 0.25
    )
  ) +
  scale_colour_identity() +
  labs(
    tag = "d",
    title = "Toxicity profiles across all tested peptides",
    x = expression(Concentration~(mu*mol~L^{-1})),
    y = NULL
  ) +
  theme(
    axis.text.x = element_text(size = 5.6),
    axis.text.y = element_text(size = 4.8),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    plot.title = element_text(margin = margin(b = 3))
  )

fig <- p_a | (p_b / p_c) +
  plot_layout(widths = c(1.18, 1.0), heights = c(0.62, 1.0)) +
  plot_annotation() &
  theme(
    plot.tag = element_text(size = 8, face = "bold"),
    plot.tag.position = c(0.01, 0.98)
  )

width_mm <- 183
height_mm <- 145
w <- width_mm / 25.4
h <- height_mm / 25.4

pdf_path <- file.path(out_dir, "Figure5b-d_in_vitro_validation.pdf")
tiff_path <- file.path(out_dir, "Figure5b-d_in_vitro_validation.tiff")
png_path <- file.path(out_dir, "Figure5b-d_in_vitro_validation.png")

grDevices::pdf(pdf_path, width = w, height = h, family = "Helvetica", useDingbats = FALSE)
print(fig)
dev.off()

ragg::agg_tiff(tiff_path, width = w, height = h, units = "in", res = 600, compression = "lzw")
print(fig)
dev.off()

ragg::agg_png(png_path, width = w, height = h, units = "in", res = 300)
print(fig)
dev.off()

write_csv(mic_raw, file.path(out_dir, "mic_blank_corrected_replicates_source_data.csv"))
write_csv(mic_mean, file.path(out_dir, "mic_inhibition_source_data.csv"))
write_csv(mic_call, file.path(out_dir, "mic_call_source_data.csv"))
write_csv(safety_64, file.path(out_dir, "toxicity_64_source_data.csv"))
write_csv(tox_heatmap, file.path(out_dir, "toxicity_heatmap_source_data.csv"))
write_csv(activity, file.path(out_dir, "cell_activity_source_data.csv"))
write_csv(hemolysis, file.path(out_dir, "haemolysis_source_data.csv"))

message("Figure written to: ", png_path)
message("Safety-window concentration: 64 μmol L⁻¹")
