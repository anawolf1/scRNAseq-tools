# ------------------------------------------------------------------------------
# qc_metrics.R
#
# QC metric names and display labels, derived from params$percent_features,
# and the display format of QC values. Used by the QC workflow and by
# population refinement.
#
# Functions
# build_qc_metric_info()   QC metric names, display labels and default scatter color
# format_qc_values()       QC values (a vector or a table's metric columns) as text
# ------------------------------------------------------------------------------

#' QC metric names, display labels and default scatter color
#'
#' The metrics are nFeature_RNA, nCount_RNA and one percent.<name> column per
#' entry of `percent_features`. Known names (mt, ribo, hb) get a descriptive
#' label, `labels` can set or override any of them, and any other name is
#' labelled "Percent <name>".
#'
#' @param percent_features Named list of gene-symbol regex patterns, e.g.
#'   params$percent_features (list(mt = "^MT-", ribo = "^RP[SL]")). An empty
#'   list gives only nFeature_RNA and nCount_RNA.
#' @param labels Optional named character vector of display labels for
#'   percent features, e.g. c(tcr = "TCR V genes (%)"). Overrides the
#'   built-in labels.
#'
#' @return A list: `qc_metrics` (character vector of metadata column names),
#'   `qc_labels` (named character vector of display labels, e.g. for
#'   plot_violin()'s feature_labels) and `scatter_color_feature` (the
#'   column to color QC scatter plots by: "percent.mt" if mt is in
#'   percent_features, else NULL).
#'
#' @examples
#' \dontrun{
#' qc_info <- build_qc_metric_info(params$percent_features)
#' qc_info <- build_qc_metric_info(list(mt = "^MT-", tcr = "^TR[AB]V"),
#'                                 labels = c(tcr = "TCR V genes (%)"))
#' }
build_qc_metric_info <- function(percent_features,
                                 labels = NULL) {

  if (!is.list(percent_features) ||
      (length(percent_features) > 0 &&
       (is.null(names(percent_features)) ||
        anyNA(names(percent_features)) ||
        any(!nzchar(names(percent_features))) ||
        anyDuplicated(names(percent_features))))) {
    stop("build_qc_metric_info(): `percent_features` must be a named list with unique, non-empty names.",
         call. = FALSE)
  }
  if (!is.null(labels) &&
      (!is.character(labels) ||
       is.null(names(labels)) ||
       anyNA(names(labels)) ||
       any(!nzchar(names(labels))) ||
       anyDuplicated(names(labels)) ||
       anyNA(labels))) {
    stop("build_qc_metric_info(): `labels` must be a named character vector with unique, ",
         "non-empty names and no missing values.",
         call. = FALSE)
  }

  percent_names <- names(percent_features) %||% character(0)

  percent_label_lookup <- c(
    mt   = "Mitochondrial genes (%)",
    ribo = "Ribosomal genes (%)",
    hb   = "Hemoglobin genes (%)"
  )
  if (!is.null(labels)) {
    percent_label_lookup[names(labels)] <- labels
  }
  percent_labels <- ifelse(
    percent_names %in% names(percent_label_lookup),
    percent_label_lookup[percent_names],
    paste0("Percent ", percent_names)
  )

  qc_metrics <- c("nFeature_RNA", "nCount_RNA", paste0("percent.", percent_names))
  qc_labels <- c(
    nFeature_RNA = "Genes detected",
    nCount_RNA   = "Total counts",
    setNames(percent_labels, paste0("percent.", percent_names))
  )

  list(
    qc_metrics            = qc_metrics,
    qc_labels             = qc_labels,
    scatter_color_feature = if ("mt" %in% percent_names) "percent.mt" else NULL
  )
}
#' QC values as text: one decimal for percent.*, none otherwise
#'
#' Takes either one metric's values or a whole table: for a data.frame, every
#' column named in `metrics` is formatted and the other columns are left as
#' they are.
#'
#' @param x Numeric vector, or a data.frame/tibble.
#' @param metrics The metric's column name for a vector (e.g. "percent.mt");
#'   for a table, the QC metric columns to format (columns not in `x` are
#'   ignored).
#'
#' @return Character vector (NA where `x` is NA), or `x` with those columns
#'   as text.
#'
#' @examples
#' format_qc_values(c(12.345, NA), "percent.mt")  # "12.3" NA
#' format_qc_values(2531.6, "nFeature_RNA")       # "2532"
#' \dontrun{
#' format_qc_values(summarize_metrics_by_groups(meta, qc_metrics, "Patient"), qc_metrics)
#' }
format_qc_values <- function(x,
                             metrics) {
  if (is.data.frame(x)) {
    for (metric in intersect(metrics, names(x))) {
      x[[metric]] <- format_qc_values(x[[metric]], metric)
    }
    return(x)
  }
  fmt <- if (startsWith(metrics, "percent.")) "%.1f" else "%.0f"
  ifelse(is.na(x), NA_character_, sprintf(fmt, x))
}
