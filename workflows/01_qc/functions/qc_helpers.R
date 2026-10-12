# ------------------------------------------------------------------------------
# qc_helpers.R
#
# Helpers for the QC workflow: summary tables of QC metrics and the QC metric
# scatter plot. Value formatting comes from format_qc_values() (qc_metrics.R).
#
# Functions
# summarize_qc_metrics()          min/quartiles/median/mean/max/missing per QC metric
# summarize_metrics_by_groups()   median of several metrics per group combination
# plot_qc_scatter()               two QC metrics, colored by a third
# ------------------------------------------------------------------------------

#' Summary table of QC metrics
#'
#' One row per metric: min, first quartile, median, mean, third quartile,
#' max and number of missing values, formatted with format_qc_values().
#'
#' @param obj Seurat object.
#' @param metrics Metadata columns to summarize, e.g. qc_metrics from
#'   build_qc_metric_info().
#'
#' @return A data.frame: metric, min, q1, median, mean, q3, max, missing.
#'
#' @examples
#' \dontrun{
#' qc_summary <- summarize_qc_metrics(data_raw, qc_metrics)
#' gt_generic(qc_summary, title = "QC metric summary (pre-filtering)")
#' }
summarize_qc_metrics <- function(obj,
                                 metrics) {
  
  validate_required_cols(obj[[]], metrics, "summarize_qc_metrics",
                         hint = "expects every element of metrics to exist in the object's metadata.")
  
  result <- do.call(
    rbind,
    lapply(metrics, function(metric) {
      x <- obj[[metric, drop = TRUE]]
      stats <- c(min    = min(x, na.rm = TRUE),
                 q1     = unname(stats::quantile(x, 0.25, na.rm = TRUE)),
                 median = stats::median(x, na.rm = TRUE),
                 mean   = mean(x, na.rm = TRUE),
                 q3     = unname(stats::quantile(x, 0.75, na.rm = TRUE)),
                 max    = max(x, na.rm = TRUE))
      data.frame(metric = metric, t(format_qc_values(stats, metric)), missing = sum(is.na(x)),
                 row.names = NULL)
    })
  )
  rownames(result) <- NULL
  result
}
#' Median of several metrics per group combination
#'
#' One row per combination of `group_vars`, one column per metric holding
#' its median (numeric; format with format_qc_values() for display).
#'
#' @param data Cell metadata (a data.frame).
#' @param metrics Numeric columns to summarize.
#' @param group_vars Grouping columns, e.g. c("Patient", "Timepoint").
#'
#' @return A tibble: `group_vars` columns, then one column per metric.
#'
#' @examples
#' \dontrun{
#' medians <- summarize_metrics_by_groups(data_filtered[[]], qc_metrics, c("Patient", "Timepoint"))
#' gt_generic(format_qc_values(medians, qc_metrics), title = "Median QC metrics")
#' }
summarize_metrics_by_groups <- function(data,
                                        metrics,
                                        group_vars) {
  
  validate_required_cols(data, c(metrics, group_vars), "summarize_metrics_by_groups",
                         hint = "expects these columns in the metadata.")
  
  data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(metrics), \(x) stats::median(x, na.rm = TRUE)),
      .groups = "drop"
    )
}
#' QC scatter plot: two metrics colored by a third
#'
#' Scatter plot of one QC metric against another (by default nCount_RNA vs
#' nFeature_RNA), optionally colored by a third metadata column (by default
#' percent.mt). Both axes can be log10-scaled, and dashed lines can mark QC
#' thresholds. Points are drawn in increasing order of color_metric, so
#' high values are on top.
#'
#' @param obj Seurat object.
#' @param x_metric Metadata column for the x-axis. Default "nCount_RNA".
#' @param y_metric Metadata column for the y-axis. Default "nFeature_RNA".
#' @param color_metric Metadata column to color points by, or NULL for no
#'   coloring. Default "percent.mt".
#' @param color_label Legend title for the color scale. Default "% mito. RNA".
#' @param color_limits Length-2 vector with the color scale range. Default
#'   c(0, 100), for a percentage; change it if color_metric is on another
#'   scale.
#' @param x_range Optional c(min, max) threshold for x_metric, drawn as dashed
#'   vertical lines. Use NA for a side without a threshold, e.g. c(NA, 1e6).
#'   NULL (default) draws no lines.
#' @param y_range Optional c(min, max) threshold for y_metric, drawn as dashed
#'   horizontal lines. Same rules as x_range.
#' @param log_scale Logical. Log10-scale both axes. Default TRUE.
#' @param title Plot title. Default "<x_metric> vs <y_metric>".
#' @param title_wrap Maximum characters per line for the title and subtitle.
#'   Default 40.
#' @param color_option viridis option for the color scale. Default
#'   resolve_palette("qc_scatter"), set with params$palettes$qc_scatter.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' # before filtering, with thresholds
#' p_before <- plot_qc_scatter(data_raw, x_range = c(NA, 1e6), y_range = c(600, 3500))
#' save_plot(p_before, "scatter_before.png", fig_dir)
#'
#' # after filtering, without thresholds
#' p_after <- plot_qc_scatter(data_filtered, title = "nCount_RNA vs nFeature_RNA (post-QC)")
#' save_plot(p_after, "scatter_after.png", fig_dir)
#' }
plot_qc_scatter <- function(obj,
                            x_metric      = "nCount_RNA",
                            y_metric      = "nFeature_RNA",
                            color_metric  = "percent.mt",
                            color_label   = "% mito. RNA",
                            color_limits  = c(0, 100),
                            x_range       = NULL,
                            y_range       = NULL,
                            log_scale     = TRUE,
                            title         = NULL,
                            title_wrap    = 40,
                            color_option  = resolve_palette("qc_scatter")) {

  ranges <- list(x_range = x_range, y_range = y_range)
  for (arg in names(ranges)) {
    r <- ranges[[arg]]
    if (!is.null(r) &&
        (!is.numeric(r) || length(r) != 2 ||
         all(is.na(r)) || any(!is.finite(r[!is.na(r)])))) {
      stop("plot_qc_scatter(): `", arg, "` must be a numeric length-2 vector ",
           "with at least one finite value.", call. = FALSE)
    }
  }

  validate_required_cols(obj[[]], c(x_metric, y_metric, color_metric), "plot_qc_scatter",
                         hint = "expects x_metric/y_metric/color_metric to exist in the object's metadata.")

  meta <- obj[[]]
  if (!is.null(color_metric)) {
    meta <- meta[order(meta[[color_metric]]), ]
  }

  has_thresholds <- !is.null(x_range) || !is.null(y_range)

  if (!is.null(color_metric)) {
    p <- ggplot2::ggplot(
      meta,
      ggplot2::aes(x = .data[[x_metric]], y = .data[[y_metric]], color = .data[[color_metric]])
    ) +
      ggplot2::geom_point(size = 0.6, alpha = 0.6) +
      ggplot2::scale_color_viridis_c(option = color_option, name = color_label, limits = color_limits)
  } else {
    p <- ggplot2::ggplot(meta, ggplot2::aes(x = .data[[x_metric]], y = .data[[y_metric]])) +
      ggplot2::geom_point(size = 0.6, alpha = 0.6, color = "#4D4D4D")
  }

  title    <- wrap_label(title %||% paste(x_metric, "vs", y_metric), width = title_wrap)
  subtitle <- wrap_label(
    paste0("n = ", ncol(obj), " cells", if (has_thresholds) "; dashed lines = QC thresholds" else ""),
    width = title_wrap
  )

  p <- p +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::labs(x = x_metric, y = y_metric, title = title, subtitle = subtitle)

  if (log_scale) {
    p <- p +
      ggplot2::scale_x_log10(name = paste0(x_metric, " (log10)"), breaks = scales::breaks_log(n = 6)) +
      ggplot2::scale_y_log10(name = paste0(y_metric, " (log10)"), breaks = scales::breaks_log(n = 6))
  }

  if (!is.null(x_range)) {
    p <- p + ggplot2::geom_vline(xintercept = x_range[!is.na(x_range)],
                                 linetype = "dashed", color = "#666666", linewidth = 0.4)
  }
  if (!is.null(y_range)) {
    p <- p + ggplot2::geom_hline(yintercept = y_range[!is.na(y_range)],
                                 linetype = "dashed", color = "#666666", linewidth = 0.4)
  }

  p
}
