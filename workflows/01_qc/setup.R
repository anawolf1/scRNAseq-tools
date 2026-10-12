# ------------------------------------------------------------------------------
# setup.R (workflows/01_qc)
#
# Sourced first by every QC pipeline script. Loads the shared functions and
# this workflow's functions, config and local overrides; checks setting
# names; loads the packages this workflow needs; resolves data/results paths
# and the raw data file; creates this workflow's results folders.
#
# Usage:
#   source(here::here("workflows", "01_qc", "setup.R"))
# ------------------------------------------------------------------------------

workflow_dir <- here::here("workflows", "01_qc")

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
    "No workflows/01_qc/config/params_override.R found - running with config.R defaults only.\n",
    "Copy config/params_override_example.R to config/params_override.R and ",
    "edit it with your dataset's values."
  )
}

open_sections <- c("exclusions", "percent_features", "mad_types", "display_order", "palettes")
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
  stop("Setting(s) removed by workflows/01_qc/config/params_override.R: ", paste(missing, collapse = ", "),
       ".\nSet single entries (e.g. params$retained_vs_excluded$padj_cutoff <- 0.01) ",
       "instead of replacing the whole list.")
}
unknown_contexts <- setdiff(names(params$palettes), palette_registry$context)
if (length(unknown_contexts) > 0) unknown <- c(unknown, paste0("params$palettes$", unknown_contexts))
if (length(unknown) > 0) {
  stop("Unknown setting(s) in workflows/01_qc/config/params_override.R: ",
       paste(unknown, collapse = ", "),
       ".\nValid names are in config/config.R (palette contexts: palette_help()).")
}

# ------------------------------------------------------------------------------
# Packages: the core set always, the others only when their section runs
# ------------------------------------------------------------------------------
load_libraries()
if (isTRUE(steps$retained_vs_excluded)) {
  load_libraries("clusterProfiler", "org.Hs.eg.db", "enrichplot", "writexl", load_core = FALSE)
}
if (isTRUE(steps$mad_estimation)) {
  load_libraries("scuttle", load_core = FALSE)
}
if (isTRUE(steps$spatial_qc)) {
  load_libraries("ggnewscale", load_core = FALSE)
}
if (!is.null(params$diagnostic_scatter)) {
  load_libraries("ggrepel", load_core = FALSE)
}

# ------------------------------------------------------------------------------
# Paths
# ------------------------------------------------------------------------------
data_raw_dir       <- params$paths$data_raw_dir       %||% here::here("data", "raw")
data_processed_dir <- params$paths$data_processed_dir %||% here::here("data", "processed")
results_dir        <- params$paths$results_dir        %||% here::here("results")

fig_qc_dir            <- file.path(results_dir, "01_qc")
fig_qc_prefilter_dir  <- file.path(fig_qc_dir, "01_qc_prefilter")
fig_qc_postfilter_dir <- file.path(fig_qc_dir, "02_qc_postfilter")

for (d in c(data_processed_dir, fig_qc_prefilter_dir, fig_qc_postfilter_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ------------------------------------------------------------------------------
# Raw data file: the single .rds in data_raw_dir, unless
# params$paths$raw_data_file is set
# ------------------------------------------------------------------------------
if (!is.null(params$paths$raw_data_file)) {
  raw_data_file <- params$paths$raw_data_file
} else {
  raw_candidates <- list.files(data_raw_dir, pattern = "\\.rds$", full.names = TRUE)
  if (length(raw_candidates) == 0) {
    stop("No .rds file found in ", data_raw_dir, " - copy your dataset there first.")
  }
  if (length(raw_candidates) > 1) {
    stop("Multiple .rds files found in ", data_raw_dir, " - set params$paths$raw_data_file ",
         "in workflows/01_qc/config/params_override.R.")
  }
  raw_data_file <- raw_candidates[1]
}
