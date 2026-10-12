# ------------------------------------------------------------------------------
# decisions.R (workflows/04_neighborhood_analysis)
#
# EXAMPLE - copy this file to config/decisions.R and edit it. Sourced by
# the pipeline and by tools/select_milo_parameters.R.
#
# What is compared, and the neighborhood parameters chosen with
# tools/select_milo_parameters.R.
# ------------------------------------------------------------------------------

# --- what is compared -----------------------------------------------------------------
# Each analysis builds one neighborhood graph from the cells of its
# timepoints and tests its contrasts, c(test, reference), on it. Including a
# timepoint in the graph does not test it. Positive
# logFC = more abundant (DA) or higher expression (DE) in `test`.
# Contrasts between two repeated-measures timepoints (3M vs FU) are tested
# paired by patient; contrasts against SPR are unpaired.
params$analyses <- list(

  # the full landscape: all four timepoints, the same cells the clustering
  # and cluster labels were built on. BL helps define the neighborhoods but
  # is not tested (2 patients).
  primary = list(
    role       = "primary",
    timepoints = c("BL", "3M", "FU", "SPR"),
    contrasts  = list(
      "3M_vs_SPR" = c("3M", "SPR"),
      "FU_vs_SPR" = c("FU", "SPR"),
      "3M_vs_FU"  = c("3M", "FU")
    )
  ),

  # sensitivity to graph construction: the same contrasts without BL cells
  without_BL = list(
    role       = "sensitivity",
    timepoints = c("3M", "FU", "SPR")
  ),

  # sensitivity to graph construction: each contrast on a graph of its two
  # timepoints only
  pairwise = list(role = "pairwise")
)

# --- DA neighborhoods (Part A of tools/select_milo_parameters.R) -----------------------
params$da$k    <- 35
params$da$d    <- 12
params$da$prop <- 0.1
params$da$sweep <- list(
  analysis    = "primary",
  k_tested    = c(20, 25, 30, 35, 40, 50),
  d_tested    = c(10, 12, 15),
  prop_tested = c(0.1, 0.2)
)

# hand-picked neighborhoods of the primary analysis, characterized on their
# own (marker DE, within-sample check, per-sample abundance); numbers as on
# the primary DA graph - they change if k, d, prop or the cells change
params$da$markers$nhood_sets <- list(
  nhood86     = 86,
  nhoods13_40 = c(13, 40)
)

# --- DE neighborhoods (Part B of tools/select_milo_parameters.R) ---------------------
params$de$k <- 8
params$de$sweep <- list(k_tested = c(5, 8, 10, 12, 15), order = 2, prop = 0.2)

# --- neighborhood groups (Part C of tools/select_milo_parameters.R) -----------------
params$da$nhood_groups$resolution <- 1
params$da$nhood_groups$sweep <- list(resolutions_tested = seq(0.2, 2, by = 0.1),
                                     n_boot = 30, frac = 0.8)
