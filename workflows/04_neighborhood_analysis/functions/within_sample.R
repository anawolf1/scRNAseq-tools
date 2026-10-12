# ------------------------------------------------------------------------------
# within_sample.R
#
# Within-sample comparison of a cell set with the other cells of the same
# sample (patient x timepoint). Differences between samples (e.g.
# sequencing depth) cancel, and the sample, not the cell, is the unit of the
# test. Used by characterize_cell_set() (milo_analysis.R) to check whether
# a set's markers are high in its cells across samples or in one sample only.
#
# Functions
# within_sample_shift()        per sample x value: median (or mean) inside vs outside a cell set
# summarize_shift()            per value: consistency across samples and signed-rank test
# plot_within_sample_shift()   paired inside/outside values, one line per sample
# ------------------------------------------------------------------------------

#' Median of each value inside vs outside a cell set, within each sample
#'
#' @param meta Cell metadata (rows = cells, rownames = cell names).
#' @param cells Cell names of the set (e.g. a meta-neighborhood).
#' @param sample_col Sample column (one value per patient x timepoint).
#' @param value_cols Numeric columns to compare (e.g. gene expression).
#' @param transform Function applied before summarizing. Default identity.
#' @param min_cells Minimum non-missing cells required both inside and
#'   outside the set within a sample. Default 3.
#' @param summary_fun Per-sample summary. Default median; use mean for 0/1
#'   indicators (gives the proportion of positive cells).
#'
#' @return Long data.frame: sample, marker, n_in, n_out, stat_in,
#'   stat_out, diff (= stat_in - stat_out, on the transformed scale).
#'   Samples failing min_cells for a marker are omitted for that marker.
within_sample_shift <- function(meta,
                                cells,
                                sample_col,
                                value_cols,
                                transform = identity,
                                min_cells = 3,
                                summary_fun = stats::median) {

  meta <- as.data.frame(meta)
  validate_required_cols(meta, c(sample_col, value_cols), "within_sample_shift")
  in_set <- rownames(meta) %in% cells

  rows <- list()
  for (s in unique(as.character(meta[[sample_col]]))) {
    in_sample <- as.character(meta[[sample_col]]) == s
    for (marker in value_cols) {
      v <- transform(as.numeric(meta[[marker]]))
      inside  <- v[in_sample & in_set & is.finite(v)]
      outside <- v[in_sample & !in_set & is.finite(v)]
      if (length(inside) < min_cells || length(outside) < min_cells) next
      rows[[paste(s, marker)]] <- data.frame(
        sample = s, marker = marker,
        n_in = length(inside), n_out = length(outside),
        stat_in = summary_fun(inside), stat_out = summary_fun(outside),
        diff = summary_fun(inside) - summary_fun(outside)
      )
    }
  }
  if (length(rows) == 0) {
    return(data.frame(sample = character(0), marker = character(0), n_in = integer(0), n_out = integer(0),
                      stat_in = numeric(0), stat_out = numeric(0), diff = numeric(0)))
  }
  dplyr::bind_rows(rows)
}
#' Consistency of within-sample shifts across samples, per marker
#'
#' @param shift_df within_sample_shift() output.
#'
#' @return Data.frame per marker: n_samples, n_higher_inside, median_diff,
#'   p_value (Wilcoxon signed-rank on the per-sample differences; NA with
#'   fewer than 3 samples) and p_adj (BH across markers).
summarize_shift <- function(shift_df) {

  validate_required_cols(shift_df, c("marker", "diff"), "summarize_shift")
  if (nrow(shift_df) == 0) return(data.frame())

  shift_df %>%
    dplyr::group_by(marker) %>%
    dplyr::summarise(
      n_samples       = dplyr::n(),
      n_higher_inside = sum(diff > 0),
      median_diff     = stats::median(diff),
      p_value         = if (dplyr::n() >= 3) suppressWarnings(stats::wilcox.test(diff, mu = 0)$p.value) else NA_real_,
      .groups = "drop"
    ) %>%
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH"),
                  significance = p_to_stars(p_adj))
}
#' Paired inside/outside medians per sample, one panel per marker
#'
#' @param shift_df within_sample_shift() output, with a `timepoint` column
#'   added for the colors.
#' @param timepoint_colors Named colors for the timepoints.
#' @param inside_label Label of the cell set (e.g. "G1").
#' @param y_label Y-axis label.
#' @param title,subtitle Plot title/subtitle.
#'
#' @return A ggplot object.
plot_within_sample_shift <- function(shift_df,
                                     timepoint_colors,
                                     inside_label = "in group",
                                     y_label = "Per-sample median (transformed)",
                                     title = NULL,
                                     subtitle = NULL) {

  validate_required_cols(shift_df, c("sample", "marker", "stat_in", "stat_out", "timepoint"),
                         "plot_within_sample_shift")
  long_df <- shift_df %>%
    tidyr::pivot_longer(c(stat_in, stat_out), names_to = "where", values_to = "median") %>%
    dplyr::mutate(where = factor(ifelse(where == "stat_in", inside_label, "rest of sample"),
                                 levels = c("rest of sample", inside_label)))

  ggplot2::ggplot(long_df, ggplot2::aes(x = where, y = median, group = sample, color = timepoint)) +
    ggplot2::geom_line(alpha = 0.7) +
    ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(values = timepoint_colors) +
    ggplot2::facet_wrap(~ marker, scales = "free_y") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = y_label, color = NULL) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(face = "bold"))
}
