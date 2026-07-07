#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_arg <- args[startsWith(args, file_arg)]

if (length(script_arg) == 1) {
  script_path <- normalizePath(sub(file_arg, "", script_arg), mustWork = TRUE)
  project_dir <- dirname(script_path)
} else {
  project_dir <- normalizePath(getwd(), mustWork = TRUE)
}

required_packages <- c(
  "rmarkdown", "dplyr", "tidyr", "readr", "stringr", "forcats",
  "lubridate", "ggplot2", "scales", "patchwork", "ragg"
)

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0) {
  stop(
    "Install the missing packages before rendering: ",
    paste(missing_packages, collapse = ", ")
  )
}

output_dir <- file.path(project_dir, "outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

rmarkdown::render(
  input = file.path(project_dir, "Supplementary_Figure_2.Rmd"),
  output_file = "Supplementary_Figure_2.html",
  output_dir = output_dir,
  knit_root_dir = project_dir,
  clean = TRUE,
  envir = new.env(parent = globalenv())
)

message("Finished. Outputs written to: ", output_dir)
