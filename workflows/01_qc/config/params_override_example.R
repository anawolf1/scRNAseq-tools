# ------------------------------------------------------------------------------
# params_override.R (workflows/01_qc)
#
# EXAMPLE - copy this file to config/params_override.R (that file name is
# gitignored, so dataset values such as patient IDs are never committed) and
# edit it for your dataset.
#
# Defaults and the meaning of every setting are in config/config.R; only set
# what differs. Set single entries (params$section$entry <- value) rather
# than replacing a whole list. The QC cutoffs you choose after running
# 01_qc_prefilter.R go in config/decisions.R.
# ------------------------------------------------------------------------------

# --- optional sections ----------------------------------------------------------
# switch a section within a script on or off (defaults in config.R)

# steps$mad_estimation       <- FALSE
# steps$retained_vs_excluded <- FALSE   # also off for non-human data
# steps$spatial_qc           <- TRUE    # plate-based protocols only


# --- cells to drop before QC metrics are computed -----------------------------
# metadata column -> values to exclude; any column works

# params$exclusions <- list(
#   Patient = c("P01"),
#   Epitope = c("E1")
# )


# --- percent features -------------------------------------------------------------
# gene-symbol regex per percent.<name> column; config.R has human symbols

# params$percent_features$mt   <- "^mt-"     # mouse
# params$percent_features$ribo <- "^Rp[sl]"  # mouse
# params$percent_features$hb   <- "^HB[^(P)]" # erythrocyte contamination (mouse: "^Hb[ab]-")


# --- MAD-suggested cutoffs -------------------------------------------------------
# direction per QC metric column ("both", "lower", "higher"); NULL = no suggestion

# params$mad_types$percent.ribo <- "higher"
# params$mad_types$percent.hb   <- "higher"
# params$mad_nmads <- 5


# --- grouping and cell counts -------------------------------------------------
# metadata columns to break QC metrics down by (the first two also pick the
# group_var1/group_var2 palettes)
params$group_by_vars <- c("Patient", "Timepoint")

# cell count breakdowns: 1-3 metadata columns each
params$cell_count_vars <- list(
  c("Patient", "Timepoint")
)

# level order in every plot and table: column -> ordered values
params$display_order <- list(
  Timepoint = c("BL", "3M", "FU", "SPR")
)


# --- excluded vs retained cells (steps$retained_vs_excluded) ---------------

# params$retained_vs_excluded$log2fc_cutoff <- 1
# params$retained_vs_excluded$padj_cutoff   <- 0.01


# --- spatial QC column names (steps$spatial_qc) --------------------------------

# params$spatial_qc$plate_col <- "plate"


# --- median scatter per group combination (02_qc_postfilter.R) ----------------
# group_vars[1] labels each point, group_vars[2] (optional) colors it;
# x_metric and y_metric must be QC metric columns

# params$diagnostic_scatter <- list(
#   group_vars = c("Patient", "Timepoint"),
#   x_metric   = "percent.mt",
#   y_metric   = "percent.ribo"
# )


# --- palettes -------------------------------------------------------------------
# palette per plot context; palette_help() lists contexts and defaults

# params$palettes$qc_scatter <- "magma"
# params$palettes$group_var1 <- "dark2"


# --- paths: only for data/results outside this repo ----------------------------

# params$paths$data_raw_dir       <- "/path/to/raw"
# params$paths$data_processed_dir <- "/path/to/processed"
# params$paths$results_dir        <- "/path/to/results"
# params$paths$raw_data_file      <- "/path/to/file.rds"
