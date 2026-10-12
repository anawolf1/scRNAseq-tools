# ------------------------------------------------------------------------------
# params_override.R (workflows/03_functional_scores)
#
# EXAMPLE - copy this file to config/params_override.R (that file name is
# gitignored, so dataset values are never committed) and edit it for your
# dataset.
#
# Defaults and the meaning of every setting are in config/config.R; only set
# what differs. Set single entries (params$section$entry <- value) rather
# than replacing a whole list. Which cells each comparison uses goes in
# config/decisions.R.
# ------------------------------------------------------------------------------

# --- optional sections ----------------------------------------------------------
# summaries and plots are on by default; the statistics need the study design below

steps$longitudinal_analysis    <- TRUE
steps$cross_sectional_analysis <- TRUE


# --- study description: the same values as in the other workflows' overrides ----

params$group_by_vars <- c("Patient", "Timepoint")

params$display_order <- list(
  Timepoint = c("BL", "3M", "FU", "SPR")
)

# the same colors as population refinement, by cluster label
params$cluster_color_map <- c(
  "Memory"   = "#56B4E9",
  "Effector" = "#F0E442",
  "Naive"    = "#E69F00"
)


# --- study design -----------------------------------------------------------------
# repeated measures: the same patients over time (tested with an LMM);
# cross-sectional: a separately sampled group (compared with Wilcoxon tests,
# each paired with a repeated timepoint, never two repeated timepoints)

params$study_design$repeated_measures_timepoints <- c("BL", "3M", "FU")
params$study_design$cross_sectional_timepoints   <- c("SPR")
# only if your columns are named differently:
# params$study_design$patient_col   <- "Donor"
# params$study_design$timepoint_col <- "Visit"


# --- gene sets: the species of every set in data/gene_sets/*.gmt ------------------

params$gene_sets$species <- c(
  EXHAUSTION_UP = "human",
  EXHAUSTION_DN = "human",
  MEMORY_MOUSE  = "mouse"
)

# symbols damaged by spreadsheets, or retired, in the .gmt files
params$gene_sets$symbol_fixes <- c(
  "4-SEP" = "SEPT4",
  "MUM1"  = "IRF4"
)


# --- statistics ---------------------------------------------------------------------

# params$statistics$min_cells_per_sample <- 10
# params$statistics$lmer_df <- "satterthwaite"
# params$statistics$cross_sectional_comparisons <- list(c("FU", "SPR"), c("3M", "SPR"))


# --- palettes and locations ------------------------------------------------------------

# params$palettes$group_var1 <- "colorbrewer1"
# params$paths$results_dir   <- "/path/to/external/results"
