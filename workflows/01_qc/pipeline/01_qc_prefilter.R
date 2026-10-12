# ------------------------------------------------------------------------------
# 01_qc_prefilter.R (workflows/01_qc)
#
# Characterizes QC metrics before cell-level filtering, to choose the QC
# cutoffs.
#
# Input:
#   Raw Seurat object (.rds) in data/raw/
#
# Output:
#   results/01_qc/01_qc_prefilter/
#     - QC metric violin plots, overall and per params$group_by_vars column
#     - QC metric summary and median tables
#     - nCount_RNA vs nFeature_RNA scatter plot
#     - cell counts per params$cell_count_vars entry
#     - MAD-suggested cutoffs [optional]
#   data/processed/
#     - 01_qc_prefilter_exclusions.rds (+ _manifest.txt)
#
# Human decision:
#   Choose the QC cutoffs and record them in config/decisions.R (used by
#   02_qc_postfilter.R)
#
# Sections:
#   0. Setup
#   1. Load raw data and apply exclusions
#   2. Compute QC metrics
#   3. Save
#   4. QC assessment
#   5. MAD-suggested cutoffs [optional]
#
# Usage:
#   source(here::here("workflows", "01_qc", "pipeline", "01_qc_prefilter.R"))
#
#   Prerequisites: a Seurat object (.rds) in data/raw/;
#   config/params_override.R copied from params_override_example.R and filled in.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------
source(here::here("workflows", "01_qc", "setup.R"))

# ------------------------------------------------------------------------------
# 1. Load raw data and apply exclusions
# ------------------------------------------------------------------------------
data_raw <- tryCatch(
  readRDS(raw_data_file),
  error = function(e) {
    stop("Could not read '", raw_data_file, "' as an .rds file: ", conditionMessage(e))
  }
)
if (!inherits(data_raw, "Seurat")) {
  stop("'", raw_data_file, "' is not a Seurat object (class: ", paste(class(data_raw), collapse = ", "),
       "). Convert it first, e.g. Seurat::as.Seurat() for a SingleCellExperiment.")
}
print(data_raw)
message("Initial cells: ", ncol(data_raw))

# cells whose value in a metadata column is listed in params$exclusions are
# dropped; a listed value that does not occur stops the run, so a typo
# cannot silently exclude nothing
for (col in names(params$exclusions)) {
  exclude_values <- params$exclusions[[col]]
  if (length(exclude_values) == 0) next

  if (!col %in% colnames(data_raw[[]])) {
    stop("Column '", col, "' not found in metadata. Check params$exclusions.")
  }
  not_found <- setdiff(exclude_values, unique(as.character(data_raw[[col, drop = TRUE]])))
  if (length(not_found) > 0) {
    stop("Value(s) not found in '", col, "': ", paste(not_found, collapse = ", "),
         ". Check params$exclusions.")
  }

  n_before <- ncol(data_raw)
  keep <- !(data_raw[[col, drop = TRUE]] %in% exclude_values)
  data_raw <- subset(data_raw, cells = colnames(data_raw)[keep])
  message("Excluded ", n_before - ncol(data_raw), " cells where ", col,
          " %in% c(", paste(exclude_values, collapse = ", "), ")")
}

# ------------------------------------------------------------------------------
# 2. Compute QC metrics
# ------------------------------------------------------------------------------
for (name in names(params$percent_features)) {
  data_raw[[paste0("percent.", name)]] <-
    Seurat::PercentageFeatureSet(data_raw, pattern = params$percent_features[[name]])
}

qc_info <- build_qc_metric_info(params$percent_features)
qc_metrics <- qc_info$qc_metrics
qc_labels  <- qc_info$qc_labels

# ------------------------------------------------------------------------------
# 3. Save
# ------------------------------------------------------------------------------
write_checkpoint(
  data_raw, "prefilter_exclusions", data_processed_dir,
  workflow = "01_qc", script = "01_qc_prefilter.R",
  settings = list(exclusions       = params$exclusions,
                  percent_features = params$percent_features),
  input_path = raw_data_file
)

# ------------------------------------------------------------------------------
# 4. QC assessment
# ------------------------------------------------------------------------------

# --- violin plot and summary of every QC metric -------------------------------
p_violin <- plot_violin(data_raw, features = qc_metrics, feature_labels = qc_labels,
                        show_points = FALSE, title = "QC metrics (pre-filtering)", ncol = 4)
print(p_violin)
save_plot(p_violin, "qc_violin_before.png", fig_qc_prefilter_dir, width = 7, height = 6)

qc_summary_gt <- gt_generic(summarize_qc_metrics(data_raw, qc_metrics),
                            title = "QC metric summary (pre-filtering)")
print(qc_summary_gt)
save_gt(qc_summary_gt, "qc_summary_before.png", fig_qc_prefilter_dir)

# --- median table and violin plot per grouping column --------------------------
for (group_col in params$group_by_vars) {

  if (!group_col %in% colnames(data_raw[[]])) {
    stop("Column '", group_col, "' not found in metadata. Check params$group_by_vars.")
  }

  medians_gt <- gt_generic(
    format_qc_values(summarize_metrics_by_groups(data_raw[[]], qc_metrics, group_col), qc_metrics),
    title = paste("Median QC metrics by", group_col, "(pre-filtering)")
  )
  print(medians_gt)
  save_gt(medians_gt, paste0("qc_medians_", group_col, "_before.png"), fig_qc_prefilter_dir)

  p_group_violin <- plot_violin(data_raw, features = qc_metrics, feature_labels = qc_labels,
                                group_col = group_col, show_points = FALSE,
                                title = paste("QC metrics by", group_col, "(pre-filtering)"), ncol = 2)
  print(p_group_violin)
  save_plot(p_group_violin, paste0("qc_violin_", group_col, "_before.png"), fig_qc_prefilter_dir,
            width = 10, height = 8)
}

# --- nCount_RNA vs nFeature_RNA, colored by percent.mt if computed -------------
# no threshold lines: no cutoff has been chosen yet
scatter_color <- qc_info$scatter_color_feature
p_scatter <- plot_qc_scatter(
  data_raw,
  color_metric = scatter_color,
  color_label  = if (!is.null(scatter_color)) qc_labels[[scatter_color]],
  title        = "nCount_RNA vs nFeature_RNA (pre-filtering)"
)
print(p_scatter)
save_plot(p_scatter, "qc_scatter_before.png", fig_qc_prefilter_dir, width = 7, height = 6)

# --- cell counts per params$cell_count_vars entry ------------------------------
plot_cell_counts(data_raw, params$cell_count_vars, params$display_order,
                 title = "Cell counts before QC", output_dir = fig_qc_prefilter_dir,
                 file_prefix = "cell_counts_before")

# ------------------------------------------------------------------------------
# 5. MAD-suggested cutoffs
# ------------------------------------------------------------------------------
# advisory only: nothing is filtered here
if (isTRUE(steps$mad_estimation)) {

  unknown_mad <- setdiff(names(params$mad_types), qc_metrics)
  if (length(unknown_mad) > 0) {
    stop("params$mad_types names metric(s) that are not QC metrics: ", paste(unknown_mad, collapse = ", "),
         ". Valid names: ", paste(qc_metrics, collapse = ", "), ".")
  }
  mad_metrics <- Filter(function(metric) !is.null(params$mad_types[[metric]]), qc_metrics)
  if (length(mad_metrics) == 0) {
    stop("steps$mad_estimation is TRUE but params$mad_types sets no metric - set at least one, ",
         "or switch the step off.")
  }

  # counts are log-transformed before the MAD (right-skewed); percentages are
  # not (bounded 0-100, and log(0) is undefined for cells with 0%)
  outliers <- lapply(setNames(mad_metrics, mad_metrics), function(metric) {
    scuttle::isOutlier(data_raw[[metric, drop = TRUE]], nmads = params$mad_nmads,
                       type = params$mad_types[[metric]], log = !startsWith(metric, "percent."))
  })

  # a suggested cutoff, or NA where isOutlier() set no bound on that side
  mad_bound <- function(metric, side) {
    th <- attr(outliers[[metric]], "thresholds")
    value <- if (!is.null(th) && side %in% names(th)) unname(th[[side]]) else NA_real_
    if (is.finite(value)) value else NA_real_
  }

  mad_table <- data.frame(
    metric          = mad_metrics,
    lower           = vapply(mad_metrics, function(m) format_qc_values(mad_bound(m, "lower"), m), character(1)),
    higher          = vapply(mad_metrics, function(m) format_qc_values(mad_bound(m, "higher"), m), character(1)),
    n_cells_flagged = vapply(outliers, sum, integer(1)),
    row.names       = NULL
  )
  mad_gt <- gt_generic(mad_table, title = "MAD-suggested QC cutoffs (advisory, not applied)",
                       subtitle = paste0(params$mad_nmads, " MADs from the median; cells flagged by each cutoff alone"))
  print(mad_gt)
  save_gt(mad_gt, "mad_cutoffs.png", fig_qc_prefilter_dir)

  # the same scatter with the suggested cutoffs drawn; no line where a
  # metric has no suggested cutoff
  mad_range <- function(metric) {
    if (!metric %in% mad_metrics) return(NULL)
    bounds <- c(mad_bound(metric, "lower"), mad_bound(metric, "higher"))
    if (all(is.na(bounds))) NULL else bounds
  }
  p_scatter_mad <- plot_qc_scatter(
    data_raw,
    color_metric = scatter_color,
    color_label  = if (!is.null(scatter_color)) qc_labels[[scatter_color]],
    x_range      = mad_range("nCount_RNA"),
    y_range      = mad_range("nFeature_RNA"),
    title        = "nCount_RNA vs nFeature_RNA (MAD-suggested cutoffs)"
  )
  print(p_scatter_mad)
  save_plot(p_scatter_mad, "qc_scatter_mad.png", fig_qc_prefilter_dir, width = 7, height = 6)

} else {
  message("[skip] mad_estimation (disabled)")
}
