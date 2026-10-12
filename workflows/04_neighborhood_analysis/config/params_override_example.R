# ------------------------------------------------------------------------------
# params_override.R (workflows/04_neighborhood_analysis)
#
# EXAMPLE - copy this file to config/params_override.R (that filename is
# gitignored, so your real values never get committed) and edit the values
# below for your own dataset.
#
# Generic defaults live in config/config.R (every valid setting name is
# there); the contrasts and neighborhood parameters live in
# config/decisions.R. A misspelled setting name stops the scripts.
# ------------------------------------------------------------------------------

# --- optional sections: switch on/off -----------------------------------------
# Defaults (all TRUE) in config.R.

# steps$da_nhood_markers <- FALSE
# steps$milo_de          <- FALSE


# --- study description: same values as in the other workflows' overrides ------

params$group_by_vars <- c("Patient", "Timepoint")

params$display_order <- list(
  Timepoint = c("BL", "3M", "FU", "SPR")
)

params$cluster_color_map <- c(
  "Memory"   = "#56B4E9",
  "Effector" = "#F0E442",
  "Naive"    = "#E69F00"
)

params$study_design$repeated_measures_timepoints <- c("BL", "3M", "FU")
params$study_design$cross_sectional_timepoints   <- c("SPR")


# --- marker tables: genes hidden from DA-group marker tables ---------------------

params$da$markers$exclude_gene_patterns <- c("^ENSG", "^RN7SL")


# --- palettes and locations ------------------------------------------------------------

# params$palettes$group_var1 <- "colorbrewer1"
# params$paths$results_dir   <- "/path/to/external/results"
