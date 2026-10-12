# ------------------------------------------------------------------------------
# decisions.R (workflows/03_functional_scores)
#
# EXAMPLE - copy this file to config/decisions.R and edit it for your
# dataset. Sourced by 01_functional_scores.R.
#
# Which cells each comparison uses. Scores are computed once for every cell
# (scoring mode "all_cells"), so choosing a cell set only changes which
# cells are compared, never their scores.
# ------------------------------------------------------------------------------

# --- cell sets --------------------------------------------------------------------
# Each set: a label (used in titles) and optional rules, each
# list(column = <metadata column>, values = <values>):
#   keep    - keep only cells with these values
#   exclude - remove cells with these values
# No rule = every cell. Any metadata column works: a protein-gated
# population (e.g. a "Naive" column from index sorting), a cluster label,
# a patient. A value that does not occur in the column stops the run.
# Check a column's values with, e.g., table(data_obj$Naive).

params$cell_sets <- list(
  all_cells = list(
    label = "All cells"
  ),
  non_naive = list(
    label   = "Naive cells excluded",
    exclude = list(column = "Naive", values = "Naive")
  )
  # e.g. by cluster instead:
  # no_naive_cluster = list(label = "Naive cluster excluded",
  #                         exclude = list(column = "cluster_label", values = "Naive"))
)

# --- which cell set(s) each comparison uses -----------------------------------------
# One name, or several to produce the comparison for each (e.g. as a
# sensitivity check: c("non_naive", "all_cells")). Results are grouped by
# cell set: results/03_functional_scores/<cell set>/.

params$analysis_cells <- list(
  cluster_summaries        = "all_cells",   # compares clusters: every cluster, every cell
  program_umaps            = "all_cells",
  timepoint_summaries      = "non_naive",
  trajectory_plot          = "non_naive",
  longitudinal_analysis    = "non_naive",
  cross_sectional_analysis = "non_naive"
)
