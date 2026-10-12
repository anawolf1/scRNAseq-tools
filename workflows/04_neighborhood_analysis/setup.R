# ------------------------------------------------------------------------------
# setup.R (workflows/04_neighborhood_analysis)
#
# Sourced first by this workflow's pipeline script and tool. Loads the
# shared functions and this workflow's functions, config and local
# overrides; checks setting names; loads the packages this workflow needs;
# resolves data and results paths and creates the results folders.
#
# This workflow reads population refinement's final checkpoint and the
# functional scores checkpoint, never an earlier workflow's config or
# decisions.
#
# Usage:
#   source(here::here("workflows", "04_neighborhood_analysis", "setup.R"))
# ------------------------------------------------------------------------------

workflow_dir <- here::here("workflows", "04_neighborhood_analysis")

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
    "No workflows/04_neighborhood_analysis/config/params_override.R found - running with config.R defaults only.\n",
    "Copy config/params_override_example.R to config/params_override.R and ",
    "edit it with your dataset's values."
  )
}

open_sections <- c("display_order", "palettes", "analyses")
unknown <- c(
  if (length(setdiff(names(steps), names(default_steps))) > 0)
    paste0("steps$", setdiff(names(steps), names(default_steps))),
  if (length(setdiff(names(params), names(default_params))) > 0)
    paste0("params$", setdiff(names(params), names(default_params)))
)
missing <- character(0)
# sections with fixed names are checked at every level (e.g. params$da$markers$...)
check_section <- function(default, current, path) {
  if (!is.list(default) || is.null(names(default))) return(invisible(NULL))
  extra <- setdiff(names(current), names(default))
  if (length(extra) > 0) unknown <<- c(unknown, paste0(path, "$", extra))
  dropped <- setdiff(names(default), names(current))
  if (length(dropped) > 0) missing <<- c(missing, paste0(path, "$", dropped))
  for (entry in intersect(names(default), names(current))) {
    check_section(default[[entry]], current[[entry]], paste0(path, "$", entry))
  }
}
for (section in setdiff(names(default_params), open_sections)) {
  check_section(default_params[[section]], params[[section]], paste0("params$", section))
}
if (length(missing) > 0) {
  stop("Setting(s) removed by workflows/04_neighborhood_analysis/config/params_override.R: ",
       paste(missing, collapse = ", "),
       ".\nSet single entries (e.g. params$da$spatial_fdr <- 0.05) instead of replacing the whole list.")
}
unknown_contexts <- setdiff(names(params$palettes), palette_registry$context)
if (length(unknown_contexts) > 0) unknown <- c(unknown, paste0("params$palettes$", unknown_contexts))
if (length(unknown) > 0) {
  stop("Unknown setting(s) in workflows/04_neighborhood_analysis/config/params_override.R: ",
       paste(unknown, collapse = ", "),
       ".\nValid names are in config/config.R (palette contexts: palette_help()).")
}

# ------------------------------------------------------------------------------
# Packages: the core set always, the others only when their section runs
# ------------------------------------------------------------------------------
load_libraries()
load_libraries("miloR", "SingleCellExperiment", "SummarizedExperiment", "writexl", "ggrepel", load_core = FALSE)
if (isTRUE(steps$da_nhood_markers) || isTRUE(steps$nhood_groups)) {
  load_libraries("igraph", load_core = FALSE)
}
if (isTRUE(steps$nhood_groups)) {
  load_libraries("edgeR", "pheatmap", load_core = FALSE)
}
if (isTRUE(steps$milo_de)) {
  load_libraries("miloDE", load_core = FALSE)
}

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
data_processed_dir <- params$paths$data_processed_dir %||% here::here("data", "processed")
results_dir        <- params$paths$results_dir        %||% here::here("results")

input_path <- if (params$input %in% names(checkpoints)) checkpoint_path(params$input, data_processed_dir) else params$input

fig_milo_dir        <- file.path(results_dir, "04_neighborhood_analysis")
fig_milo_da_dir     <- file.path(fig_milo_dir, "01_milo_da")
fig_milo_de_dir     <- file.path(fig_milo_dir, "02_milo_de")
fig_milo_sweeps_dir <- file.path(fig_milo_dir, "exploration")

for (d in c(data_processed_dir, fig_milo_da_dir, fig_milo_de_dir, fig_milo_sweeps_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
