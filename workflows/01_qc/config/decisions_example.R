# ------------------------------------------------------------------------------
# decisions.R (workflows/01_qc)
#
# EXAMPLE - copy this file to config/decisions.R and edit it for your
# dataset.
#
# QC cutoffs chosen after reviewing 01_qc_prefilter.R's output (violin
# plots, scatter plot, summary tables, MAD-suggested cutoffs). Sourced by
# 02_qc_postfilter.R.
#
# Only the metrics listed here are filtered; each takes `min`, `max` or both.
# Names are metadata columns (e.g. percent.ribo). Bounds are exclusive:
# cells are kept when value > min and value < max.
# ------------------------------------------------------------------------------

params$thresholds <- list(
  nFeature_RNA = list(min = 600, max = 3500),
  nCount_RNA   = list(max = 1000000),
  percent.mt   = list(max = 50)
)
