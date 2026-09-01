options(repos = c(CRAN = "https://cloud.r-project.org"))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

cran_packages <- c(
  "renv", "data.table", "ggplot2", "ggrepel", "pheatmap", "readxl",
  "Seurat", "SeuratObject", "Matrix", "patchwork", "dplyr", "coloc",
  "scales", "gridExtra", "cowplot", "ggpubr", "svglite", "ragg"
)
bioc_packages <- c(
  "GEOquery", "limma", "clusterProfiler", "org.Hs.eg.db", "edgeR"
)

missing_cran <- cran_packages[!vapply(cran_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) install.packages(missing_cran)
missing_bioc <- bioc_packages[!vapply(bioc_packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)

message("Dependencies installed. For the exact recorded environment, run renv::restore().")

