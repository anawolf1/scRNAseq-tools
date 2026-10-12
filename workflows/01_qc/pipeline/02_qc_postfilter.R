# ------------------------------------------------------------------------------
# 02_qc_postfilter.R (workflows/01_qc)
#
# Applies the QC cutoffs chosen after reviewing 01_qc_prefilter.R's output,
# then assesses the result: how many cells were retained, whether QC metrics
# still differ across the study design, and (optionally) whether the
# excluded cells differ from the retained ones in a biologically meaningful
# way.
#
# Input:
#   data/processed/01_qc_prefilter_exclusions.rds (from 01_qc_prefilter.R)
#   config/decisions.R (the chosen cutoffs; see decisions_example.R)
#
# Output:
#   results/01_qc/02_qc_postfilter/
#     - QC metric violin plots, overall and per params$group_by_vars column
#     - QC metric summary and median tables
#     - nCount_RNA vs nFeature_RNA scatter plot with the cutoffs
#     - median scatter per group combination [optional]
#     - cells retained, cells remaining and percent retained
#       (qc_retained_*, qc_cells_remaining_*, qc_pct_retained_*)
#     - excluded vs retained: DE table, GO/KEGG dot plots and
#       qc_de_enrichment.xlsx [optional]
#     - QC by plate/well position [optional]
#   data/processed/
#     - 01_qc_filtered.rds (+ _manifest.txt): input to 02_population_refinement
#
# Sections:
#   0. Setup
#   1. Load the prefiltered data
#   2. Apply the QC cutoffs
#   3. Save
#   4. QC assessment
#   5. Excluded vs retained cells [optional]
#   6. Spatial QC [optional]
#
# Usage:
#   source(here::here("workflows", "01_qc", "pipeline", "02_qc_postfilter.R"))
#
#   Prerequisites: 01_qc_prefilter.R run; config/decisions.R copied from
#   decisions_example.R and filled in.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------
source(here::here("workflows", "01_qc", "setup.R"))

decision_file <- file.path(workflow_dir, "config", "decisions.R")
if (!file.exists(decision_file)) {
  stop("No workflows/01_qc/config/decisions.R found. Run 01_qc_prefilter.R first, look at its ",
       "output, then copy config/decisions_example.R to config/decisions.R and fill in your cutoffs.")
}
source(decision_file)

if (is.null(params$thresholds)) {
  stop("params$thresholds is not set. Check workflows/01_qc/config/decisions.R.")
}

# ------------------------------------------------------------------------------
# 1. Load the prefiltered data
# ------------------------------------------------------------------------------
prefiltered_file <- checkpoint_path("prefilter_exclusions", data_processed_dir)
if (!file.exists(prefiltered_file)) {
  stop("No prefilter_exclusions checkpoint found. Run 01_qc_prefilter.R first.")
}
data_prefiltered <- readRDS(prefiltered_file)
print(data_prefiltered)

qc_info <- build_qc_metric_info(params$percent_features)
qc_metrics <- qc_info$qc_metrics
qc_labels  <- qc_info$qc_labels

# ------------------------------------------------------------------------------
# 2. Apply the QC cutoffs
# ------------------------------------------------------------------------------
# Only metrics listed in params$thresholds are filtered, each by `min`, `max`
# or both. Bounds are exclusive (value > min, value < max): a cell exactly
# on a cutoff is excluded. A cell with a missing value is excluded.
unknown_metrics <- setdiff(names(params$thresholds), colnames(data_prefiltered[[]]))
if (length(unknown_metrics) > 0) {
  stop("params$thresholds names metric(s) not found in the metadata: ",
       paste(unknown_metrics, collapse = ", "), ". Check config/decisions.R ",
       "(names must match metadata columns, e.g. percent.mt).")
}
for (metric in names(params$thresholds)) {
  bounds <- params$thresholds[[metric]]
  if (length(bounds) == 0 || !all(names(bounds) %in% c("min", "max")) ||
      !all(vapply(bounds, function(b) is.numeric(b) && length(b) == 1, logical(1)))) {
    stop("params$thresholds$", metric, " must be a list with a numeric `min` and/or `max`. ",
         "Check config/decisions.R.")
  }
}

qc_retained <- rep(TRUE, ncol(data_prefiltered))
for (metric in names(params$thresholds)) {
  bounds <- params$thresholds[[metric]]
  values <- data_prefiltered[[metric, drop = TRUE]]
  if (!is.null(bounds$min)) qc_retained <- qc_retained & !is.na(values) & values > bounds$min
  if (!is.null(bounds$max)) qc_retained <- qc_retained & !is.na(values) & values < bounds$max
}
if (!any(qc_retained)) {
  stop("No cells pass the cutoffs in config/decisions.R - check them against the distributions ",
       "in results/01_qc/01_qc_prefilter/.")
}

# retained/excluded label on every cell, used in sections 4-6
data_prefiltered$qc <- ifelse(qc_retained, "retained", "excluded")
data_filtered <- subset(data_prefiltered, cells = colnames(data_prefiltered)[qc_retained])

message("QC filtering: ", ncol(data_filtered), " of ", ncol(data_prefiltered), " cells retained (",
        round(100 * ncol(data_filtered) / ncol(data_prefiltered), 1), "%). Cutoffs on: ",
        paste(names(params$thresholds), collapse = ", "))

# ------------------------------------------------------------------------------
# 3. Save
# ------------------------------------------------------------------------------
write_checkpoint(
  data_filtered, "qc_filtered", data_processed_dir,
  workflow = "01_qc", script = "02_qc_postfilter.R",
  settings = list(thresholds = params$thresholds),
  input_path = prefiltered_file
)

# ------------------------------------------------------------------------------
# 4. QC assessment
# ------------------------------------------------------------------------------

# --- violin plot and summary of every QC metric -------------------------------
p_violin <- plot_violin(data_filtered, features = qc_metrics, feature_labels = qc_labels,
                        show_points = FALSE, title = "QC metrics (post-filtering)", ncol = 4)
print(p_violin)
save_plot(p_violin, "qc_violin_after.png", fig_qc_postfilter_dir, width = 7, height = 6)

qc_summary_gt <- gt_generic(summarize_qc_metrics(data_filtered, qc_metrics),
                            title = "QC metric summary (post-filtering)")
print(qc_summary_gt)
save_gt(qc_summary_gt, "qc_summary_after.png", fig_qc_postfilter_dir)

# --- median table and violin plot per grouping column --------------------------
for (group_col in params$group_by_vars) {

  if (!group_col %in% colnames(data_filtered[[]])) {
    stop("Column '", group_col, "' not found in metadata. Check params$group_by_vars.")
  }

  medians_gt <- gt_generic(
    format_qc_values(summarize_metrics_by_groups(data_filtered[[]], qc_metrics, group_col), qc_metrics),
    title = paste("Median QC metrics by", group_col, "(post-filtering)")
  )
  print(medians_gt)
  save_gt(medians_gt, paste0("qc_medians_", group_col, "_after.png"), fig_qc_postfilter_dir)

  p_group_violin <- plot_violin(data_filtered, features = qc_metrics, feature_labels = qc_labels,
                                group_col = group_col, show_points = FALSE,
                                title = paste("QC metrics by", group_col, "(post-filtering)"), ncol = 2)
  print(p_group_violin)
  save_plot(p_group_violin, paste0("qc_violin_", group_col, "_after.png"), fig_qc_postfilter_dir,
            width = 10, height = 8)
}

# --- median scatter per group combination (params$diagnostic_scatter) ---------
# one point per combination of group_vars: group_vars[1] labels each point,
# group_vars[2] (if given) colors it
if (!is.null(params$diagnostic_scatter)) {

  ds <- params$diagnostic_scatter
  if (length(ds$group_vars) < 1 || is.null(ds$x_metric) || is.null(ds$y_metric)) {
    stop("params$diagnostic_scatter needs group_vars (at least 1 column), x_metric and y_metric.")
  }
  missing_metrics <- setdiff(c(ds$x_metric, ds$y_metric), qc_metrics)
  if (length(missing_metrics) > 0) {
    stop("params$diagnostic_scatter: metric(s) not in qc_metrics (check params$percent_features): ",
         paste(missing_metrics, collapse = ", "))
  }
  label_var <- ds$group_vars[1]
  color_var <- if (length(ds$group_vars) >= 2) ds$group_vars[2] else NULL

  medians <- summarize_metrics_by_groups(data_filtered[[]], c(ds$x_metric, ds$y_metric), ds$group_vars)
  medians_gt <- gt_generic(
    format_qc_values(medians, qc_metrics),
    title = paste("Median QC metrics by", paste(ds$group_vars, collapse = " and "), "(post-filtering)")
  )
  print(medians_gt)
  save_gt(medians_gt, paste0("qc_medians_", paste(ds$group_vars, collapse = "_"), "_after.png"),
          fig_qc_postfilter_dir)

  p_medians <- ggplot2::ggplot(medians, ggplot2::aes(x = .data[[ds$x_metric]], y = .data[[ds$y_metric]]))
  if (!is.null(color_var)) {
    color_levels <- ordered_levels(medians[[color_var]], params$display_order[[color_var]])
    p_medians <- p_medians +
      ggplot2::geom_point(ggplot2::aes(color = .data[[color_var]]), size = 3) +
      ggplot2::scale_color_manual(values = resolve_group_colors(color_var, color_levels),
                                  breaks = color_levels)
  } else {
    p_medians <- p_medians + ggplot2::geom_point(size = 3)
  }
  p_medians <- p_medians +
    ggrepel::geom_text_repel(ggplot2::aes(label = .data[[label_var]]), size = 3, show.legend = FALSE) +
    ggplot2::labs(
      title = wrap_label(paste0("Median ", qc_labels[[ds$x_metric]], " vs. ", qc_labels[[ds$y_metric]],
                                ", by ", paste(ds$group_vars, collapse = " and ")), width = 50),
      x = paste("Median", qc_labels[[ds$x_metric]]),
      y = paste("Median", qc_labels[[ds$y_metric]])
    ) +
    ggplot2::theme_classic(base_size = 13)
  print(p_medians)
  save_plot(p_medians, paste0("qc_median_scatter_", paste(ds$group_vars, collapse = "_"), ".png"),
            fig_qc_postfilter_dir, width = 7, height = 6)

} else {
  message("[skip] diagnostic_scatter (not configured)")
}

# --- nCount_RNA vs nFeature_RNA, with the cutoffs ------------------------------
# c(min, max) of a metric's cutoffs (NA where only one is set), or NULL for
# no line
threshold_range <- function(bounds) {
  if (is.null(bounds)) return(NULL)
  c(bounds$min %||% NA_real_, bounds$max %||% NA_real_)
}
scatter_color <- qc_info$scatter_color_feature
p_scatter <- plot_qc_scatter(
  data_filtered,
  color_metric = scatter_color,
  color_label  = if (!is.null(scatter_color)) qc_labels[[scatter_color]],
  x_range      = threshold_range(params$thresholds$nCount_RNA),
  y_range      = threshold_range(params$thresholds$nFeature_RNA),
  title        = "nCount_RNA vs nFeature_RNA (post-filtering)"
)
print(p_scatter)
save_plot(p_scatter, "qc_scatter_after.png", fig_qc_postfilter_dir, width = 7, height = 6)

# --- cells retained per params$cell_count_vars entry ---------------------------
plot_retention(data_prefiltered, data_filtered, params$cell_count_vars, params$display_order,
               title_prefix = "QC filtering", output_dir = fig_qc_postfilter_dir, file_prefix = "qc")

# ------------------------------------------------------------------------------
# 5. Excluded vs retained cells
# ------------------------------------------------------------------------------
# DE (excluded vs retained) and GO/KEGG enrichment of the genes up in each
# group: shows whether the cutoffs removed a biologically distinct population
# rather than low-quality cells
if (isTRUE(steps$retained_vs_excluded) && all(qc_retained)) {
  message("[skip] retained_vs_excluded (no cells were excluded)")
} else if (isTRUE(steps$retained_vs_excluded)) {

  rve <- params$retained_vs_excluded
  run_de_enrichment(
    data_prefiltered, group_col = "qc", ident_1 = "excluded", ident_2 = "retained",
    min_pct = rve$min_pct, log2fc_cutoff = rve$log2fc_cutoff, padj_cutoff = rve$padj_cutoff,
    exclude_gene_patterns = rve$exclude_gene_patterns, n_terms = rve$n_terms,
    label_1 = "Excluded", label_2 = "Retained",
    output_dir = fig_qc_postfilter_dir, file_prefix = "qc"
  )

} else {
  message("[skip] retained_vs_excluded (disabled)")
}

# ------------------------------------------------------------------------------
# 6. Spatial QC
# ------------------------------------------------------------------------------
if (isTRUE(steps$spatial_qc)) {

  sq <- params$spatial_qc

  # --- retained/excluded by well ------------------------------------------------
  p_spatial_qc <- plot_plate_grid(
    data_prefiltered[[]], fill_col = "qc",
    plate_col = sq$plate_col, row_col = sq$row_col, column_col = sq$column_col,
    title = "QC status by well position"
  )
  print(p_spatial_qc)
  save_plot(p_spatial_qc, "spatial_qc_status.png", fig_qc_postfilter_dir, width = 9, height = 7)

  # --- each grouping column by well ---------------------------------------------
  for (group_col in params$group_by_vars) {
    p_spatial_group <- plot_plate_grid(
      data_prefiltered[[]], fill_col = group_col, type = "discrete",
      plate_col = sq$plate_col, row_col = sq$row_col, column_col = sq$column_col,
      title = paste(group_col, "by well position")
    )
    print(p_spatial_group)
    save_plot(p_spatial_group, paste0("spatial_", group_col, ".png"), fig_qc_postfilter_dir,
              width = 16, height = 12)
  }

  # --- each QC metric by well: all cells, then retained cells only ---------------
  # both maps share one color scale (the range over all cells)
  for (metric in qc_metrics) {
    metric_range <- range(data_prefiltered[[metric, drop = TRUE]], na.rm = TRUE)

    p_metric_all <- plot_plate_grid(
      data_prefiltered[[]], fill_col = metric, fill_limits = metric_range,
      plate_col = sq$plate_col, row_col = sq$row_col, column_col = sq$column_col,
      title = paste(qc_labels[[metric]], "by well position (all cells)")
    )
    print(p_metric_all)
    save_plot(p_metric_all, paste0("spatial_", metric, "_all.png"), fig_qc_postfilter_dir,
              width = 16, height = 12)

    p_metric_retained <- plot_plate_grid(
      data_prefiltered[[]], fill_col = metric, fill_limits = metric_range,
      discard_col = "qc", discard_value = "excluded",
      plate_col = sq$plate_col, row_col = sq$row_col, column_col = sq$column_col,
      title = paste(qc_labels[[metric]], "by well position (QC-retained cells)")
    )
    print(p_metric_retained)
    save_plot(p_metric_retained, paste0("spatial_", metric, "_retained.png"), fig_qc_postfilter_dir,
              width = 16, height = 12)
  }

} else {
  message("[skip] spatial_qc (disabled)")
}
