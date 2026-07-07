# Supplementary Figure 2 v2

This folder contains a self-contained and reproducible reconstruction of
Supplementary Figure 2 from the HG002 Nextflow execution trace.

## Contents

- `Supplementary_Figure_2.Rmd`: complete data-processing and plotting workflow.
- `render_figure.R`: one-command renderer.
- `data/execution_trace_HG002.txt`: source Nextflow trace used by the figure.
- `outputs/`: publication-quality figure files and processed TSV tables.

## Reproduce

```bash
cd suplementary_figure2_v2
Rscript render_figure.R
```

Required R packages:

```r
install.packages(c(
  "rmarkdown", "dplyr", "tidyr", "readr", "stringr", "forcats",
  "lubridate", "ggplot2", "scales", "patchwork", "cowplot", "ragg"
))
```

## Final figure formats

- PDF: vector graphics, recommended for the manuscript.
- PNG: 600 dpi with a white background.
- TIFF: 600 dpi, LZW-compressed, suitable for journal submission systems.

The workflow also writes task-level and process-level TSV files so that every
plotted value can be audited independently of the figure.
