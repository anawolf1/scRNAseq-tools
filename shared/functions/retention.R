# ------------------------------------------------------------------------------
# retention.R
#
# Cell counts per group, and cells kept by a filtering or exclusion step:
# counts remaining and percent retained. Tables and heatmaps are printed
# and saved.
#
# Functions
# plot_cell_counts()   cell counts per breakdown (table, and heatmap for 2-3 columns)
# plot_retention()     cells remaining and percent retained, per breakdown
# ------------------------------------------------------------------------------

#' Cell counts per breakdown
#'
#' For each entry of `cell_count_vars`:
#'   - 1 column: a table of cells per level
#'   - 2 columns: a row x column table (with totals) and a heatmap
#'   - 3 columns: one table with a section per level of the third column,
#'     and a heatmap faceted by it
#' Levels follow `display_order` where given (see ordered_levels()). Cells
#' with a missing (NA) value in a breakdown column are not counted in that
#' breakdown.
#'
#' @param obj Seurat object.
#' @param cell_count_vars List of character vectors (1-3 metadata columns
#'   each), e.g. list("Patient", c("Patient", "Timepoint")).
#' @param display_order Named list (column -> ordered levels), e.g.
#'   params$display_order.
#' @param title Title of every table and heatmap.
#' @param output_dir Directory to save into.
#' @param file_prefix Prefix for every file name: <prefix>_<columns>.png and
#'   <prefix>_heatmap_<columns>.png.
#'
#' @return Invisibly, NULL.
#'
#' @examples
#' \dontrun{
#' plot_cell_counts(data_raw, params$cell_count_vars, params$display_order,
#'                  title = "Cell counts before QC", output_dir = fig_dir,
#'                  file_prefix = "cell_counts_before")
#' }
plot_cell_counts <- function(obj,
                             cell_count_vars,
                             display_order = list(),
                             title,
                             output_dir,
                             file_prefix) {

  for (vars in cell_count_vars) {
    if (!length(vars) %in% 1:3) {
      stop("plot_cell_counts(): each cell_count_vars entry needs 1-3 columns; got ", length(vars), ".",
           call. = FALSE)
    }
    validate_required_cols(obj[[]], vars, "plot_cell_counts",
                           hint = "expects every cell_count_vars column to exist in the object's metadata.")
    label <- paste(vars, collapse = "_")
    values <- lapply(vars, function(var) {
      x <- obj[[var, drop = TRUE]]
      factor(x, levels = ordered_levels(x, display_order[[var]]))
    })

    if (length(vars) == 1) {
      counts_df <- setNames(as.data.frame(table(values[[1]]), stringsAsFactors = FALSE),
                            c(vars[1], "n_cells"))
      counts_gt <- gt_generic(counts_df, title = title, subtitle = vars[1])
      print(counts_gt)
      save_gt(counts_gt, paste0(file_prefix, "_", label, ".png"), output_dir)
      next
    }

    counts <- do.call(table, values)
    subtitle <- if (length(vars) == 2) {
      paste(vars[1], "x", vars[2])
    } else {
      paste0(vars[1], " x ", vars[2], ", by ", vars[3])
    }
    counts_gt <- if (length(vars) == 2) {
      gt_contingency(counts, title = title, subtitle = subtitle)
    } else {
      gt_contingency_faceted(counts, title = title, subtitle = subtitle)
    }
    print(counts_gt)
    save_gt(counts_gt, paste0(file_prefix, "_", label, ".png"), output_dir)

    p_counts <- plot_count_heatmap(
      value_matrix = counts, row_col = vars[1], col_col = vars[2],
      facet_col = if (length(vars) == 3) vars[3],
      title = title, subtitle = subtitle,
      row_order = levels(values[[1]]), col_order = levels(values[[2]]),
      facet_order = if (length(vars) == 3) levels(values[[3]])
    )
    print(p_counts)
    save_plot(p_counts, paste0(file_prefix, "_heatmap_", label, ".png"), output_dir,
              width = if (length(vars) == 3) 9 else 6, height = 5)
  }

  invisible(NULL)
}
#' Cells remaining and percent retained after a filtering step, per breakdown
#'
#' For each entry of `cell_count_vars`:
#'   - 1 column: one table with n_before, n_after and pct_retained per level
#'   - 2 columns: cells remaining (plot_cell_counts()) and percent retained
#'     (table + heatmap, diverging around the overall retention, so groups
#'     that lost more than average stand out)
#'   - 3 columns: cells remaining only (plot_cell_counts())
#' Levels follow `display_order` where given (see ordered_levels()). Cells
#' with a missing (NA) value in a breakdown column are not counted in that
#' breakdown.
#'
#' @param obj_before,obj_after Seurat objects before and after the step
#'   (`obj_after`'s cells must be a subset of `obj_before`'s).
#' @param cell_count_vars List of character vectors (1-3 metadata columns
#'   each), e.g. list(c("Patient", "Timepoint")).
#' @param display_order Named list (column -> ordered levels), e.g.
#'   params$display_order.
#' @param title_prefix Text starting every title (e.g. "QC filtering" or a
#'   round's label).
#' @param output_dir Directory to save into.
#' @param file_prefix Prefix for every file name: <prefix>_retained_*,
#'   <prefix>_cells_remaining_* and <prefix>_pct_retained_*.
#'
#' @return Invisibly, the overall percent of cells retained.
#'
#' @examples
#' \dontrun{
#' plot_retention(data_prefiltered, data_filtered, params$cell_count_vars,
#'                params$display_order, title_prefix = "QC filtering",
#'                output_dir = fig_dir, file_prefix = "qc")
#' }
plot_retention <- function(obj_before,
                           obj_after,
                           cell_count_vars,
                           display_order = list(),
                           title_prefix,
                           output_dir,
                           file_prefix) {

  if (ncol(obj_before) == 0) {
    stop("plot_retention(): obj_before contains no cells.", call. = FALSE)
  }
  if (!all(colnames(obj_after) %in% colnames(obj_before))) {
    stop("plot_retention(): obj_after contains cells that are not in obj_before.", call. = FALSE)
  }
  overall <- 100 * ncol(obj_after) / ncol(obj_before)

  # a column's values in `obj` as a factor with the levels of `obj_before`,
  # ordered by display_order
  as_ordered <- function(obj, var) {
    lv <- ordered_levels(obj_before[[var, drop = TRUE]], display_order[[var]])
    factor(obj[[var, drop = TRUE]], levels = lv)
  }

  for (vars in cell_count_vars) {
    if (!length(vars) %in% 1:3) {
      stop("plot_retention(): each cell_count_vars entry needs 1-3 columns; got ", length(vars), ".",
           call. = FALSE)
    }
    validate_required_cols(obj_before[[]], vars, "plot_retention")
    label <- paste(vars, collapse = "_")

    # --- 1 column: before, after and percent retained in one table -------------
    if (length(vars) == 1) {
      before <- table(as_ordered(obj_before, vars[1]))
      n_before <- as.integer(before)
      n_after  <- as.integer(table(as_ordered(obj_after, vars[1])))
      counts_df <- data.frame(
        level        = names(before),
        n_before     = n_before,
        n_after      = n_after,
        pct_retained = ifelse(n_before > 0, round(100 * n_after / n_before, 0), NA_real_)
      )
      names(counts_df)[1] <- vars[1]
      counts_gt <- gt_generic(counts_df, title = paste(title_prefix, "— cells retained"), subtitle = vars[1])
      print(counts_gt)
      save_gt(counts_gt, paste0(file_prefix, "_retained_", label, ".png"), output_dir)
      next
    }

    # --- 2-3 columns: cells remaining -------------------------------------------
    plot_cell_counts(obj_after, list(vars), display_order,
                     title = paste(title_prefix, "— cells remaining"),
                     output_dir = output_dir, file_prefix = paste0(file_prefix, "_cells_remaining"))
    if (length(vars) == 3) next

    # --- 2 columns: percent retained ---------------------------------------------
    row_var <- vars[1]
    col_var <- vars[2]
    before <- table(as_ordered(obj_before, row_var), as_ordered(obj_before, col_var))
    after  <- table(as_ordered(obj_after, row_var), as_ordered(obj_after, col_var))

    pct <- round(100 * after / before, 0)
    pct[is.nan(pct)] <- NA
    pct_gt <- gt_contingency(pct, title = paste(title_prefix, "— percent retained"),
                             subtitle = paste(row_var, "x", col_var), total_row = FALSE, total_col = FALSE)
    print(pct_gt)
    save_gt(pct_gt, paste0(file_prefix, "_pct_retained_", label, ".png"), output_dir)

    p_pct <- plot_count_heatmap(
      value_matrix = pct, row_col = row_var, col_col = col_var,
      title = paste(title_prefix, "— percent retained"), subtitle = paste(row_var, "x", col_var),
      row_order = rownames(pct), col_order = colnames(pct),
      label_suffix = "%", diverging = TRUE,
      # kept strictly inside 0-100, which the heatmap requires
      midpoint = min(max(overall, 0.5), 99.5), limits = c(0, 100)
    )
    print(p_pct)
    save_plot(p_pct, paste0(file_prefix, "_pct_retained_heatmap_", label, ".png"),
              output_dir, width = 6, height = 5)
  }

  invisible(overall)
}
