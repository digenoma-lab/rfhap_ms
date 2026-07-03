#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# Resolve paths relative to this script so it can be run from any directory.
script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(script_arg) == 1L) {
  normalizePath(sub("^--file=", "", script_arg), mustWork = TRUE)
} else {
  normalizePath("build_supplementary_figure4.R", mustWork = TRUE)
}
script_dir <- dirname(script_path)

input_file <- file.path(script_dir, "Supplementary_Table_ModelPerformance.csv")
pdf_file <- file.path(script_dir, "output", "pdf", "Supplementary_Figure_4.pdf")
png_file <- file.path(script_dir, "output", "png", "Supplementary_Figure_4.png")

dir.create(dirname(pdf_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(png_file), recursive = TRUE, showWarnings = FALSE)

metrics_df <- read.csv(input_file, check.names = TRUE)

# Sensitivity is intentionally excluded because it is mathematically
# equivalent to recall: TP / (TP + FN).
metric_order <- c(
  "F1",
  "Balanced Accuracy",
  "Precision",
  "Recall",
  "Specificity"
)

plot_df <- metrics_df |>
  select(
    Dataset,
    Chemistry,
    Class,
    F1,
    Balanced.Accuracy,
    Precision,
    Recall,
    Specificity
  ) |>
  rename(`Balanced Accuracy` = Balanced.Accuracy) |>
  pivot_longer(
    cols = all_of(metric_order),
    names_to = "Metric",
    values_to = "Value"
  ) |>
  mutate(Metric = factor(Metric, levels = metric_order))

make_metric_plot <- function(metric_name, show_x = FALSE, show_y = FALSE) {
  ggplot(
    filter(plot_df, Metric == metric_name),
    aes(x = Dataset, y = Value, shape = Class, color = Chemistry)
  ) +
    geom_point(size = 3, position = position_dodge(width = 0.45)) +
    scale_y_continuous(
      breaks = seq(0.98, 1.00, by = 0.005),
      labels = scales::label_number(accuracy = 0.001)
    ) +
    coord_cartesian(ylim = c(0.98, 1.00)) +
    labs(
      title = metric_name,
      x = NULL,
      y = if (show_y) "Metric value" else NULL,
      color = "Chemistry",
      shape = "Class"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x = if (show_x) {
        element_text(angle = 30, hjust = 1)
      } else {
        element_blank()
      },
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      plot.margin = margin(6, 7, 6, 7)
    )
}

p_f1 <- make_metric_plot("F1")
p_balanced <- make_metric_plot("Balanced Accuracy")
p_precision <- make_metric_plot("Precision")
p_recall <- make_metric_plot("Recall", show_x = TRUE)
p_specificity <- make_metric_plot("Specificity", show_x = TRUE)

# Twelve equal columns allow a small optical correction. Each panel spans four
# columns, while the bottom row is shifted one column left of the exact center.
figure_design <- "
AAAABBBBCCCC
#DDDDEEEE###
"

plots_grid <-
  p_f1 + p_balanced + p_precision + p_recall + p_specificity +
  plot_layout(
    design = figure_design,
    guides = "collect",
    heights = c(1, 1)
  ) &
  theme(legend.position = "bottom")

shared_y_label <- wrap_elements(
  full = grid::textGrob(
    "Metric value",
    rot = 90,
    gp = grid::gpar(fontsize = 12)
  )
)

p_perf <-
  shared_y_label + plots_grid +
  plot_layout(widths = c(0.035, 1)) +
  plot_annotation(
    title = "RFhap model performance across datasets",
    theme = theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16)
    )
  )

ggsave(
  filename = pdf_file,
  plot = p_perf,
  width = 9,
  height = 6.8,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = png_file,
  plot = p_perf,
  width = 9,
  height = 6.8,
  units = "in",
  dpi = 400,
  bg = "white"
)

message("Created: ", pdf_file)
message("Created: ", png_file)
