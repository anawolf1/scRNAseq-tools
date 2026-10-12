# ------------------------------------------------------------------------------
# setup.R (workflows/02_population_refinement)
#
# Sourced first by the pipeline script and every tool of this workflow.
# Loads the shared functions and this workflow's functions, config and local
# overrides; checks setting names; loads the packages this workflow needs;
# resolves data/results paths and creates this workflow's results folders.
#
# This workflow reads the QC workflow's checkpoint
# (data/processed/01_qc_filtered.rds), never data/raw/.
#
# Usage:
#   source(here::here("workflows", "02_population_refinement", "setup.R"))
# ------------------------------------------------------------------------------

workflow_dir <- here::here("workflows", "02_population_refinement")

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
    "No workflows/02_population_refinement/config/params_override.R found - running with config.R defaults only.\n",
    "Copy config/params_override_example.R to config/params_override.R and ",
    "edit it with your dataset's values."
  )
}

open_sections <- c("percent_features", "display_order", "palettes")
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
  stop("Setting(s) removed by workflows/02_population_refinement/config/params_override.R: ",
       paste(missing, collapse = ", "),
       ".\nSet single entries (e.g. params$de$padj_cutoff <- 0.01) instead of replacing the whole list.")
}
unknown_contexts <- setdiff(names(params$palettes), palette_registry$context)
if (length(unknown_contexts) > 0) unknown <- c(unknown, paste0("params$palettes$", unknown_contexts))
if (length(unknown) > 0) {
  stop("Unknown setting(s) in workflows/02_population_refinement/config/params_override.R: ",
       paste(unknown, collapse = ", "),
       ".\nValid names are in config/config.R (palette contexts: palette_help()).")
}

# ------------------------------------------------------------------------------
# Packages: the core set always, the others only when their section runs
# ------------------------------------------------------------------------------
load_libraries()
load_libraries("writexl", load_core = FALSE)
if (isTRUE(steps$differential_expression) || isTRUE(steps$exclusion_check)) {
  load_libraries("EnhancedVolcano", load_core = FALSE)
}
if (isTRUE(steps$deg_overlap)) {
  load_libraries("ggvenn", load_core = FALSE)
}
if (isTRUE(steps$exclusion_check)) {
  load_libraries("clusterProfiler", "org.Hs.eg.db", "enrichplot", load_core = FALSE)
}

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
data_processed_dir <- params$paths$data_processed_dir %||% here::here("data", "processed")
results_dir        <- params$paths$results_dir        %||% here::here("results")

fig_refinement_dir           <- file.path(results_dir, "02_population_refinement")
fig_exploration_dir          <- file.path(fig_refinement_dir, "exploration")
fig_clustering_sweeps_dir    <- file.path(fig_exploration_dir, "clustering_sweeps")
fig_candidate_inspection_dir <- file.path(fig_exploration_dir, "candidate_inspection")

for (d in c(data_processed_dir, fig_refinement_dir, fig_clustering_sweeps_dir, fig_candidate_inspection_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
