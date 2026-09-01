# Shared path bootstrap for every authoritative R script.
# Run scripts from the repository root, or set PIEZO1_PROJECT_ROOT explicitly.

PIEZO1_ROOT <- Sys.getenv("PIEZO1_PROJECT_ROOT", unset = getwd())
PIEZO1_ROOT <- normalizePath(PIEZO1_ROOT, winslash = "/", mustWork = TRUE)
PIEZO1_DATA <- normalizePath(
  Sys.getenv("PIEZO1_DATA_DIR", unset = file.path(PIEZO1_ROOT, "data")),
  winslash = "/", mustWork = FALSE
)
PIEZO1_RESULTS <- normalizePath(
  Sys.getenv("PIEZO1_RESULTS_DIR", unset = file.path(PIEZO1_ROOT, "results")),
  winslash = "/", mustWork = FALSE
)
PIEZO1_FIGURES <- normalizePath(
  Sys.getenv("PIEZO1_FIGURES_DIR", unset = file.path(PIEZO1_ROOT, "figures")),
  winslash = "/", mustWork = FALSE
)
PIEZO1_METADATA <- normalizePath(
  Sys.getenv("PIEZO1_METADATA_DIR", unset = file.path(PIEZO1_ROOT, "metadata")),
  winslash = "/", mustWork = FALSE
)
PIEZO1_REFERENCE <- file.path(PIEZO1_ROOT, "reference")

dir.create(PIEZO1_RESULTS, recursive = TRUE, showWarnings = FALSE)
dir.create(PIEZO1_FIGURES, recursive = TRUE, showWarnings = FALSE)
dir.create(PIEZO1_METADATA, recursive = TRUE, showWarnings = FALSE)
setwd(PIEZO1_ROOT)

data_path <- function(...) file.path(PIEZO1_DATA, ...)
result_path <- function(...) file.path(PIEZO1_RESULTS, ...)
figure_path <- function(...) file.path(PIEZO1_FIGURES, ...)
metadata_path <- function(...) file.path(PIEZO1_METADATA, ...)
reference_path <- function(...) file.path(PIEZO1_REFERENCE, ...)

