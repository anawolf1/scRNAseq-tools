# ------------------------------------------------------------------------------
# setup.R (workflows/03_functional_scores)
#
# Sourced first by this workflow's pipeline script. Loads the shared
# functions and this workflow's functions, config and local overrides;
# checks setting names; loads the packages this workflow needs; resolves
# data, gene-set and results paths and creates the results folder.
#
# This workflow reads population refinement's final checkpoint
# (data/processed/02_population_refined.rds), never an earlier workflow's
# config or decisions.
#
# Usage:
#   source(here::here("workflows", "03_functional_scores", "setup.R"))
# ------------------------------------------------------------------------------

workflow_dir <- here::here("workflows", "03_functional_scores")

source(here::here("shared", "functions", "load_all.R"))
source_folder(file.path(workflow_dir, "functions"))
source(file.path(workflow_dir, "config", "config.R"))

message("scRNAseq-tools version ", trimws(readLines(here::here("VERSION"), warn = FALSE)[1]))

# ------------------------------------------------------------------------------
# Local, dataset-specific values
# ------------------------------------------------------------------------------
default_params <- params
default_steps  <- steps

override_file <- file.path(workflow_dir, "config", "params_override.R")
if (file.exists(override_file)) {
  source(override_file)
} else {
  message(
    "No workflows/03_functional_scores/config/params_override.R found - running with config.R defaults only.\n",
    "Copy config/params_override_example.R to config/params_override.R and ",
    "edit it with your dataset's values."
  )
}

open_sections <- c("cell_sets", "display_order", "palettes")
unknown <- c(
  if (length(setdiff(names(steps), names(default_steps))) > 0)
    paste0("steps$", setdiff(names(steps), names(default_steps))),
  if (length(setdiff(names(params), names(default_params))) > 0)
    paste0("params$", setdiff(names(params), names(default_params)))
)
missing <- character(0)
for (section in setdiff(names(default_params), open_sections)) {
  if (is.list(default_params[[section]]) && !is.null(names(default_params[[section]]))) {
    extra <- setdiff(names(params[[section]]), names(default_params[[section]]))
    if (length(extra) > 0) unknown <- c(unknown, paste0("params$", section, "$", extra))
    dropped <- setdiff(names(default_params[[section]]), names(params[[section]]))
    if (length(dropped) > 0) missing <- c(missing, paste0("params$", section, "$", dropped))
  }
}
if (length(missing) > 0) {
  stop("Setting(s) removed by workflows/03_functional_scores/config/params_override.R: ",
       paste(missing, collapse = ", "),
       ".\nSet single entries (e.g. params$statistics$min_cells_per_sample <- 10) instead of replacing the whole list.")
}
unknown_contexts <- setdiff(names(params$palettes), palette_registry$context)
if (length(unknown_contexts) > 0) unknown <- c(unknown, paste0("params$palettes$", unknown_contexts))
if (length(unknown) > 0) {
  stop("Unknown setting(s) in workflows/03_functional_scores/config/params_override.R: ",
       paste(unknown, collapse = ", "),
       ".\nValid names are in config/config.R (palette contexts: palette_help()).")
}

# ------------------------------------------------------------------------------
# Packages: the core set always, the others only when their section runs
# ------------------------------------------------------------------------------
load_libraries()
load_libraries("writexl", "limma", "babelgene", load_core = FALSE)   # workbooks; symbol/ortholog resolution
if (isTRUE(steps$longitudinal_analysis)) {
  load_libraries("lme4", "emmeans", load_core = FALSE)
}
if (isTRUE(steps$longitudinal_analysis) || isTRUE(steps$cross_sectional_analysis)) {
  load_libraries("ggpubr", load_core = FALSE)   # significance brackets
}

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
data_processed_dir <- params$paths$data_processed_dir %||% here::here("data", "processed")
results_dir        <- params$paths$results_dir        %||% here::here("results")
gene_sets_dir      <- params$gene_sets$dir            %||% here::here("data", "gene_sets")

gene_sets_files <- params$gene_sets$files %||% list.files(gene_sets_dir, pattern = "\\.gmt$", full.names = TRUE)

fig_functional_dir <- file.path(results_dir, "03_functional_scores")
for (d in c(data_processed_dir, fig_functional_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
