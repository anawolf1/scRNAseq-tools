# ------------------------------------------------------------------------------
# plots_generic.R
#
# ggplot builders used across workflows. None of them save to disk; call
# save_plot() separately, e.g. save_plot(p, "violin.png", fig_dir).
#
# Functions
# wrap_label()                wrap long title/subtitle text onto several lines
# needs_white_text()          white or dark text, from a background color's luminance
# ordered_levels()            levels of a column, in display order where one is given
# plot_count_heatmap()        count/value heatmap of two metadata columns, optionally faceted
# plot_violin()               violin plot of one or more metrics, optionally grouped
# plot_umap_by_group()        UMAP colored by a categorical variable, combined and per level
# plot_plate_grid()           a metadata column by well position, one panel per plate
# plot_sweep_tile()           one sweep metric across a two-parameter grid, values in the tiles
# ------------------------------------------------------------------------------

#' Wrap long title/subtitle text onto several lines
#'
#' ggplot2 does not wrap titles; in plots with coord_fixed() a long title can
#' overflow or distort the panel.
#'
#' @param text Character or NULL.
#' @param width Maximum characters per line.
#'
#' @return Character with line breaks, or NULL if `text` is NULL.
wrap_label <- function(text,
                       width = 40) {
  if (is.null(text)) return(NULL)
  paste(strwrap(text, width = width), collapse = "\n")
}
#' White or dark text, from a background color's luminance
#'
#' Uses perceived luminance (ITU-R BT.709 weights; an approximation, without
#' the sRGB linearization of the full WCAG definition).
#'
#' @param colors Colors (hex or any value grDevices::col2rgb() accepts).
#' @param threshold Luminance (0-1) below which text is white. Default 0.5.
#'
#' @return Logical vector, one per color: TRUE = white text.
needs_white_text <- function(colors,
                             threshold = 0.5) {
  
  # rows: red/green/blue; columns: colors
  rgb_vals  <- grDevices::col2rgb(colors) / 255
  luminance <- 0.2126 * rgb_vals["red", ] + 0.7152 * rgb_vals["green", ] + 0.0722 * rgb_vals["blue", ]
  luminance < threshold
}
#' Levels of a column, in display order where one is given
#'
#' The levels present in `values` (NA dropped): first those listed in
#' `display_order`, in that order, then any others alphabetically. Levels in
#' `display_order` that do not occur are left out.
#'
#' @param values Vector of values (character or factor).
#' @param display_order Optional preferred order, e.g.
#'   params$display_order$Timepoint. NULL (default) sorts alphabetically.
#'
#' @return Character vector of levels.
#'
#' @examples
#' ordered_levels(c("FU", "BL", "3M", "X"), c("BL", "3M", "FU", "SPR"))  # "BL" "3M" "FU" "X"
ordered_levels <- function(values,
                           display_order = NULL) {
  present <- unique(as.character(values[!is.na(values)]))
  listed  <- intersect(as.character(display_order), present)
  c(listed, sort(setdiff(present, listed)))
}
#' Count/value heatmap of two metadata columns, optionally faceted by a third
#'
#' Tabulates two categorical metadata columns (or takes a precomputed
#' matrix, e.g. percent retained) and draws a tile heatmap with the value
#' printed in each tile. Text is white or dark depending on each tile's fill.
#'
#' @param obj Seurat object. Ignored if `value_matrix` is given.
#' @param row_col Metadata column for rows (e.g. "Patient").
#' @param col_col Metadata column for columns (e.g. "Timepoint").
#' @param value_matrix Optional precomputed table/matrix to plot instead of
#'   tabulating `obj`; its dimnames are the categories. A 3-way table works
#'   together with `facet_col`. Default NULL.
#' @param facet_col Optional third metadata column to facet by.
#' @param title,subtitle Plot title/subtitle.
#' @param title_wrap Maximum characters per title/subtitle line. Default 40.
#' @param row_order,col_order,facet_order Optional level orders. Default NULL
#'   (table() order).
#' @param label_suffix Text appended to each printed value (e.g. "%").
#'   Default "".
#' @param label_size Font size of the printed values. Default 4.5; lower it
#'   for many rows/columns.
#' @param diverging Logical. TRUE: low/mid/high scale around `midpoint`
#'   (values with a meaningful center); FALSE (default): low-to-high scale.
#' @param low_color,high_color,mid_color Scale colors. NULL (default) takes
#'   them from the "count_heatmap" (diverging) or "count_heatmap_sequential"
#'   palette; each can be set on its own.
#' @param midpoint Center of the diverging scale; required when
#'   `diverging = TRUE`. With `limits`, it must lie strictly inside them;
#'   without, the scale extends equally on both sides of it.
#' @param limits Optional length-2 range of the color scale.
#' @param text_threshold Luminance (0-1) below which tile text is white (see
#'   needs_white_text()). Default 0.5.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_count_heatmap(data_raw, "Patient", "Timepoint",
#'                    title = "Cells before QC", col_order = c("BL", "3M", "FU", "SPR"))
#'
#' # a 3-way table, faceted by its third dimension
#' tbl3 <- table(meta$Patient, meta$Timepoint, meta$Sex)
#' plot_count_heatmap(value_matrix = tbl3, row_col = "Patient", col_col = "Timepoint",
#'                    facet_col = "Sex", title = "Cells by sex")
#' }
plot_count_heatmap <- function(obj = NULL,
                               row_col,
                               col_col,
                               value_matrix = NULL,
                               facet_col = NULL,
                               title = NULL,
                               subtitle = NULL,
                               title_wrap = 40,
                               row_order = NULL,
                               col_order = NULL,
                               facet_order = NULL,
                               label_suffix = "",
                               label_size = 4.5,
                               diverging = FALSE,
                               low_color = NULL,
                               high_color = NULL,
                               mid_color = NULL,
                               midpoint = NULL,
                               limits = NULL,
                               text_threshold = 0.5) {
  
  if (diverging &&
      (!is.numeric(midpoint) || length(midpoint) != 1 || !is.finite(midpoint))) {
    stop("plot_count_heatmap(): `midpoint` must be a single finite numeric value when diverging = TRUE.",
         call. = FALSE)
  }
  if (!is.null(limits) &&
      (!is.numeric(limits) || length(limits) != 2 ||
       any(!is.finite(limits)) || limits[1] >= limits[2])) {
    stop("plot_count_heatmap(): `limits` must be a finite increasing numeric length-2 vector.",
         call. = FALSE)
  }
  if (diverging && !is.null(limits) &&
      (midpoint <= limits[1] || midpoint >= limits[2])) {
    stop("plot_count_heatmap(): midpoint (", midpoint, ") must lie strictly inside limits (",
         limits[1], ", ", limits[2], ").", call. = FALSE)
  }
  if (!is.numeric(text_threshold) || length(text_threshold) != 1 ||
      !is.finite(text_threshold) || text_threshold < 0 || text_threshold > 1) {
    stop("plot_count_heatmap(): `text_threshold` must be a single value between 0 and 1.",
         call. = FALSE)
  }
  
  if (diverging) {
    div_pal    <- resolve_palette("count_heatmap")
    low_color  <- low_color  %||% div_pal$low
    mid_color  <- mid_color  %||% div_pal$mid
    high_color <- high_color %||% div_pal$high
  } else {
    seq_pal    <- resolve_palette("count_heatmap_sequential")
    low_color  <- low_color  %||% seq_pal$low
    high_color <- high_color %||% seq_pal$high
  }
  
  if (!is.null(value_matrix)) {
    # precomputed values: 2D without facet_col, 3D with it
    validate_table_dims(
      value_matrix, expected_dims = if (is.null(facet_col)) 2 else 3, fn_name = "plot_count_heatmap",
      hint = if (is.null(facet_col)) "Pass facet_col for a 3-way value_matrix." else "Omit facet_col for a plain 2D value_matrix."
    )
    df <- as.data.frame(as.table(as.array(value_matrix)))
  } else {
    validate_required_cols(obj[[]], c(row_col, col_col, facet_col), "plot_count_heatmap",
                           hint = "expects row_col/col_col/facet_col to exist in the object's metadata.")
    if (!is.null(facet_col)) {
      df <- as.data.frame(table(
        obj[[row_col, drop = TRUE]], obj[[col_col, drop = TRUE]], obj[[facet_col, drop = TRUE]]
      ))
    } else {
      df <- as.data.frame(table(obj[[row_col, drop = TRUE]], obj[[col_col, drop = TRUE]]))
    }
  }
  # table() output gets the metadata column names; the last column is the
  # value (a count, or the precomputed value from value_matrix)
  names(df) <- if (!is.null(facet_col)) c(row_col, col_col, facet_col, "value") else c(row_col, col_col, "value")
  
  if (!is.null(row_order)) df[[row_col]] <- factor(df[[row_col]], levels = row_order)
  if (!is.null(col_order)) df[[col_col]] <- factor(df[[col_col]], levels = col_order)
  if (!is.null(facet_col) && !is.null(facet_order)) {
    df[[facet_col]] <- factor(df[[facet_col]], levels = facet_order)
  }
  
  df$label <- ifelse(is.na(df$value), "", paste0(df$value, label_suffix))
  
  # tile colors, computed here to choose the text color per tile; the fill
  # scale below uses the same colors, stops and range
  if (diverging) {
    if (!is.null(limits)) {
      value_range <- limits
    } else {
      extent <- suppressWarnings(max(abs(df$value - midpoint), na.rm = TRUE))
      if (!is.finite(extent) || extent == 0) extent <- 1
      value_range <- midpoint + c(-extent, extent)
    }
    fill_colors <- c(low_color, mid_color, high_color)
    fill_stops  <- scales::rescale(c(value_range[1], midpoint, value_range[2]))
  } else {
    value_range <- limits %||% range(df$value, na.rm = TRUE)
    fill_colors <- c(low_color, high_color)
    fill_stops  <- c(0, 1)
  }
  color_fun <- function(x) {
    scales::gradient_n_pal(fill_colors, fill_stops)(scales::rescale(x, from = value_range))
  }
  
  df$text_white <- FALSE
  has_value <- !is.na(df$value)
  df$text_white[has_value] <- needs_white_text(color_fun(df$value[has_value]), threshold = text_threshold)
  
  title    <- wrap_label(title, width = title_wrap)
  subtitle <- wrap_label(subtitle, width = title_wrap)
  
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data[[col_col]], y = .data[[row_col]], fill = value)) +
    ggplot2::geom_tile(color = "#FFFFFF", linewidth = 1.2, na.rm = FALSE) +
    ggplot2::geom_text(ggplot2::aes(label = label, color = text_white), size = label_size, fontface = "bold") +
    ggplot2::scale_color_manual(values = c("TRUE" = "#FFFFFF", "FALSE" = "#333333"), guide = "none") +
    ggplot2::scale_fill_gradientn(colors = fill_colors, values = fill_stops,
                                  limits = value_range, na.value = "#E5E5E5") +
    ggplot2::labs(title = title, subtitle = subtitle) +
    # square tiles
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 16),
      plot.subtitle   = ggplot2::element_text(color = "#666666"),
      axis.title      = ggplot2::element_blank(),
      axis.text       = ggplot2::element_text(face = "bold", size = 12),
      panel.grid      = ggplot2::element_blank(),
      legend.position = "none"
    )
  
  if (!is.null(facet_col)) {
    p <- p + ggplot2::facet_wrap(stats::as.formula(paste("~", facet_col)))
  }
  
  p
}
#' Violin plot of one or more metrics, optionally grouped
#'
#' One violin per feature, or, with `group_col`, one violin per group level
#' and one panel per feature. Optional boxplot, jittered points and
#' threshold lines (drawn only in their feature's panel).
#'
#' @param obj Seurat object.
#' @param features Metadata columns to plot.
#' @param group_col Optional metadata column to group by (e.g. "Patient").
#' @param title,subtitle Plot title/subtitle. Default NULL.
#' @param title_wrap Maximum characters per title/subtitle line. Default 40.
#' @param fill_color Violin fill without `group_col`. Default "#B3B3B3".
#' @param group_colors Named colors, one per level of `group_col`. NULL
#'   (default) uses resolve_group_colors().
#' @param group_order Levels of `group_col` to show, in order; cells of other
#'   levels are left out. NULL (default) shows every level, ordered by
#'   params$display_order (see ordered_levels()).
#' @param feature_labels Optional named vector of panel labels, e.g.
#'   c(nFeature_RNA = "Genes detected"). Default NULL (column names).
#' @param thresholds Optional named list, one entry per feature, each with
#'   `min` and/or `max` (the shape of params$thresholds). Drawn as dashed lines.
#' @param jitter_size,jitter_width,jitter_alpha Jittered point size, width
#'   and transparency. Defaults 0.15, 0.08, 0.25.
#' @param show_points Logical. Show jittered cells. Default TRUE.
#' @param show_boxplot Logical. Show a narrow boxplot. Default TRUE.
#' @param ncol Panel columns when several features are plotted. NULL
#'   (default) = one row.
#' @param legend_position Legend position with `group_col`. Default "right".
#' @param seed Seed for the jitter, so the figure is identical on every run.
#'   Default 123.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_violin(data_raw, features = qc_metrics, title = "QC metrics (pre-filtering)", ncol = 4)
#' plot_violin(data_raw, features = qc_metrics, group_col = "Patient")
#' plot_violin(data_raw, features = qc_metrics, thresholds = params$thresholds, ncol = 4)
#' }
plot_violin <- function(obj,
                        features,
                        group_col        = NULL,
                        title            = NULL,
                        subtitle         = NULL,
                        title_wrap       = 40,
                        fill_color       = "#B3B3B3",
                        group_colors     = NULL,
                        group_order      = NULL,
                        feature_labels   = NULL,
                        thresholds       = NULL,
                        jitter_size      = 0.15,
                        jitter_width     = 0.08,
                        jitter_alpha     = 0.25,
                        show_points      = TRUE,
                        show_boxplot     = TRUE,
                        ncol             = NULL,
                        legend_position  = "right",
                        seed             = 123) {
  
  meta <- obj[[]]
  
  if (!is.character(features) || length(features) == 0 ||
      anyNA(features) || anyDuplicated(features)) {
    stop("plot_violin(): `features` must be a non-empty character vector of unique, non-missing column names.",
         call. = FALSE)
  }
  
  validate_required_cols(meta, c(features, group_col), "plot_violin",
                         hint = "expects every element of features (and group_col, if set) to exist in the object's metadata.")
  
  plot_data <- meta[, c(features, group_col), drop = FALSE] %>%
    tibble::rownames_to_column("cell") %>%
    tidyr::pivot_longer(cols = dplyr::all_of(features), names_to = "feature", values_to = "value")
  plot_data$feature <- factor(plot_data$feature, levels = features)
  
  if (!is.null(group_col)) {
    group_levels <- group_order %||% ordered_levels(plot_data[[group_col]], params$display_order[[group_col]])
    plot_data[[group_col]] <- factor(plot_data[[group_col]], levels = group_levels)
    plot_data <- plot_data[!is.na(plot_data[[group_col]]), , drop = FALSE]
  
    group_colors <- group_colors %||% resolve_group_colors(group_col, group_levels)
    x_aes    <- ggplot2::sym(group_col)
    fill_aes <- ggplot2::sym(group_col)
  } else {
    x_aes    <- ggplot2::sym("feature")
    fill_aes <- NULL
  }
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = !!x_aes, y = value, fill = !!fill_aes))
  
  if (is.null(group_col)) {
    p <- p + ggplot2::geom_violin(fill = fill_color, color = "#000000", linewidth = 0.4, trim = TRUE)
  } else {
    p <- p +
      ggplot2::geom_violin(color = "#000000", linewidth = 0.4, trim = TRUE, position = "identity") +
      ggplot2::scale_fill_manual(values = group_colors)
  }
  
  if (show_boxplot) {
    p <- p + ggplot2::geom_boxplot(
      width = 0.15, fill = "#FFFFFF", color = "#000000",
      linewidth = 0.35, outlier.shape = NA, position = "identity"
    )
  }
  
  if (show_points) {
    p <- p + ggplot2::geom_point(
      position = ggplot2::position_jitter(width = jitter_width, height = 0, seed = seed),
      size = jitter_size, alpha = jitter_alpha
    )
  }
  
  if (!is.null(thresholds)) {
    threshold_rows <- lapply(names(thresholds), function(feature_name) {
      values <- unlist(thresholds[[feature_name]][c("min", "max")], use.names = FALSE)
      if (length(values) == 0) return(NULL)
      data.frame(feature = feature_name, value = values)
    })
    threshold_df <- do.call(rbind, threshold_rows)
  
    if (!is.null(threshold_df) && nrow(threshold_df) > 0) {
      threshold_df$feature <- factor(threshold_df$feature, levels = features)
      threshold_df <- threshold_df[!is.na(threshold_df$feature), ]
      p <- p + ggplot2::geom_hline(
        data = threshold_df, ggplot2::aes(yintercept = value),
        linetype = "dashed", linewidth = 0.5, inherit.aes = FALSE
      )
    }
  }
  
  title    <- wrap_label(title, width = title_wrap)
  subtitle <- wrap_label(subtitle, width = title_wrap)
  
  p <- p +
    ggplot2::labs(x = NULL, y = NULL, title = title, subtitle = subtitle) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      axis.text.x     = ggplot2::element_text(angle = 45, hjust = 1, color = "#000000"),
      axis.text.y     = ggplot2::element_text(color = "#000000"),
      axis.title      = ggplot2::element_text(face = "bold"),
      axis.line       = ggplot2::element_line(linewidth = 0.5),
      axis.ticks      = ggplot2::element_line(linewidth = 0.5),
      panel.grid      = ggplot2::element_blank(),
      legend.position = if (is.null(group_col)) "none" else legend_position,
      legend.title    = ggplot2::element_blank()
    )
  
  if (length(features) > 1) {
    p <- p + ggplot2::facet_wrap(
      ~ feature,
      scales   = if (is.null(group_col)) "free" else "free_y",
      ncol     = ncol %||% length(features),
      labeller = if (!is.null(feature_labels)) ggplot2::as_labeller(feature_labels) else ggplot2::label_value
    )
    # without groups, the panel label already names the feature
    if (is.null(group_col)) {
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
    }
  }
  
  p
}
#' UMAP colored by a categorical variable, combined and per level
#'
#' Two views of the same embedding: all cells colored by `group_col`, and
#' one panel per level with only that level colored over a grey background
#' (shows where each level sits without the others competing for attention).
#'
#' @param obj Seurat object with the reduction computed.
#' @param group_col Categorical metadata column.
#' @param reduction Reduction to plot. Default "umap".
#' @param group_order Levels of `group_col` to show, in order; cells of other
#'   levels are left out. NULL (default) shows every level, ordered by
#'   params$display_order (see ordered_levels()).
#' @param title Title of the combined plot; the per-level plot gets the same
#'   title plus "(one panel per level)". Default "UMAP by <group_col>".
#' @param subtitle Subtitle of the combined plot.
#' @param title_wrap Maximum characters per title/subtitle line. Default 40.
#' @param group_colors Named colors, one per level. NULL (default) uses
#'   resolve_group_colors().
#' @param background_color Color of the other cells in the per-level
#'   panels. Default light grey.
#' @param point_size Point size; background points are drawn smaller.
#'   Default 1.2.
#'
#' @return A list of two ggplot objects: `combined` and `faceted`.
#'
#' @examples
#' \dontrun{
#' umap_patient <- plot_umap_by_group(data_final, "Patient")
#' print(umap_patient$combined)
#' print(umap_patient$faceted)
#' }
plot_umap_by_group <- function(obj,
                               group_col,
                               reduction = "umap",
                               group_order = NULL,
                               title = NULL,
                               subtitle = NULL,
                               title_wrap = 40,
                               group_colors = NULL,
                               background_color = "#D9D9D9",
                               point_size = 1.2) {
  
  if (!reduction %in% Seurat::Reductions(obj)) {
    stop("plot_umap_by_group(): reduction '", reduction, "' not found in object. ",
         "Available reductions: ", paste(Seurat::Reductions(obj), collapse = ", "),
         call. = FALSE)
  }
  validate_required_cols(obj[[]], group_col, "plot_umap_by_group",
                         hint = "expects group_col to exist in the object's metadata.")
  
  title <- title %||% paste("UMAP by", group_col)
  
  embed <- Seurat::Embeddings(obj, reduction)[, 1:2]
  colnames(embed) <- c("Dim1", "Dim2")
  df <- as.data.frame(embed)
  df[[group_col]] <- obj[[group_col, drop = TRUE]]
  
  levels_use <- group_order %||% ordered_levels(df[[group_col]], params$display_order[[group_col]])
  df[[group_col]] <- factor(as.character(df[[group_col]]), levels = levels_use)
  df <- df[!is.na(df[[group_col]]), , drop = FALSE]
  
  group_colors <- group_colors %||% resolve_group_colors(group_col, levels_use)
  
  combined_title    <- wrap_label(title, width = title_wrap)
  combined_subtitle <- wrap_label(subtitle, width = title_wrap)
  faceted_title     <- wrap_label(paste(title, "(one panel per level)"), width = title_wrap)
  
  p_combined <- ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2, color = .data[[group_col]])) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_manual(values = group_colors, name = group_col) +
    ggplot2::labs(title = combined_title, subtitle = combined_subtitle,
                  x = paste0(reduction, "_1"), y = paste0(reduction, "_2")) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(color = "#666666", size = 10),
                   axis.line = ggplot2::element_line(color = "#4D4D4D"),
                   axis.ticks = ggplot2::element_line(color = "#4D4D4D"))
  
  # one copy of the cells per level, so a single ggplot can facet
  df_facet <- dplyr::bind_rows(lapply(levels_use, function(lvl) {
    d <- df
    d$facet_level <- lvl
    d$highlight <- as.character(d[[group_col]]) == lvl
    d
  }))
  
  p_faceted <- ggplot2::ggplot() +
    ggplot2::geom_point(data = df_facet[!df_facet$highlight, ],
                        ggplot2::aes(x = Dim1, y = Dim2), color = background_color, size = point_size * 0.7) +
    ggplot2::geom_point(data = df_facet[df_facet$highlight, ],
                        ggplot2::aes(x = Dim1, y = Dim2, color = facet_level), size = point_size) +
    ggplot2::scale_color_manual(values = group_colors, guide = "none") +
    ggplot2::facet_wrap(~ factor(facet_level, levels = levels_use)) +
    ggplot2::labs(title = faceted_title,
                  x = paste0(reduction, "_1"), y = paste0(reduction, "_2")) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold"),
                   strip.text = ggplot2::element_text(face = "bold"),
                   axis.line = ggplot2::element_line(color = "#4D4D4D"),
                   axis.ticks = ggplot2::element_line(color = "#4D4D4D"))
  
  list(combined = p_combined, faceted = p_faceted)
}
#' A metadata column by well position, one panel per plate
#'
#' Fixed-grid tile plot of every well, faceted by plate. Categorical columns
#' (QC status, patient, cluster) use a discrete scale, numeric columns (QC
#' metrics) a viridis gradient. Categories keep their factor order (e.g.
#' clusters 0, 1, ..., 10), or follow params$display_order (see
#' ordered_levels()). Wells with no row in `data` are shown as "not sorted";
#' with `discard_col`, cells equal to `discard_value` are shown as
#' "discarded". A kept well with an NA `fill_col` gets the scale's NA color.
#' Expects at most one row per plate/row/column.
#'
#' @param data Cell metadata with plate, row and column columns.
#' @param fill_col Column to color wells by.
#' @param plate_col,row_col,column_col Plate, row and column columns.
#'   `row_col` must hold numbers (the row axis runs from min to max);
#'   `column_col` is sorted numerically if all its values are numbers,
#'   alphabetically otherwise.
#' @param type "discrete" or "continuous". NULL (default) detects it from
#'   `fill_col`.
#' @param fill_values Colors for a discrete `fill_col`, named by category or
#'   in category order. NULL (default): the "qc_status" palette when every
#'   category is one of its names (retained/excluded), otherwise
#'   resolve_group_colors(). For clusters pass the cluster palette, e.g.
#'   resolve_palette("cluster", n = 5).
#' @param fill_limits Optional range of the continuous color scale.
#' @param discard_col Optional column marking cells to grey out (e.g. "qc").
#' @param discard_value Value of `discard_col` marking a discarded cell.
#'   Default "excluded".
#' @param not_sorted_color,discarded_color Greys for empty and discarded wells.
#' @param title,subtitle Plot title/subtitle.
#' @param title_wrap Maximum characters per title/subtitle line. Default 40.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_plate_grid(data_prefiltered[[]], fill_col = "qc", title = "QC status by well position")
#' plot_plate_grid(data_final[[]], fill_col = "seurat_clusters",
#'                 fill_values = resolve_palette("cluster", n = 5), title = "Clusters by well position")
#' plot_plate_grid(data_prefiltered[[]], fill_col = "Patient", title = "Patient by well position")
#' }
plot_plate_grid <- function(data,
                            fill_col,
                            plate_col = "Plate",
                            row_col = "Row",
                            column_col = "Column",
                            type = NULL,
                            fill_values = NULL,
                            fill_limits = NULL,
                            discard_col = NULL,
                            discard_value = "excluded",
                            not_sorted_color = "#E5E5E5",
                            discarded_color = "#999999",
                            title = NULL,
                            subtitle = "Fixed grid, all panels same dimensions; light grey = not sorted, dark grey = discarded",
                            title_wrap = 40) {
  
  validate_required_cols(data, c(plate_col, row_col, column_col, fill_col, discard_col),
                         "plot_plate_grid",
                         hint = "expects plate_col/row_col/column_col/fill_col (and discard_col, if set) to exist in data.")
  
  df <- data %>%
    dplyr::rename(Plate    = dplyr::all_of(plate_col),
                  Row      = dplyr::all_of(row_col),
                  Column   = dplyr::all_of(column_col),
                  fill_var = dplyr::all_of(fill_col))
  
  if (!is.null(discard_col)) {
    df <- df %>% dplyr::rename(discard_var = dplyr::all_of(discard_col))
  } else {
    df$discard_var <- FALSE
  }
  
  if (is.null(type)) {
    type <- if (is.numeric(df$fill_var)) "continuous" else "discrete"
  } else if (!type %in% c("discrete", "continuous")) {
    stop("plot_plate_grid(): `type` must be either 'discrete' or 'continuous'.", call. = FALSE)
  }
  
  row_numeric <- suppressWarnings(as.numeric(as.character(df$Row)))
  if (anyNA(row_numeric)) {
    stop("plot_plate_grid(): row_col must contain values coercible to numeric.", call. = FALSE)
  }
  
  all_rows <- as.character(min(row_numeric):max(row_numeric))
  # numeric column labels in numeric order (1, 2, 10), others alphabetically
  col_values  <- unique(as.character(df$Column))
  col_numeric <- suppressWarnings(as.numeric(col_values))
  all_cols    <- if (anyNA(col_numeric)) sort(col_values) else col_values[order(col_numeric)]
  df$Row    <- factor(as.character(df$Row), levels = all_rows)
  df$Column <- factor(as.character(df$Column), levels = all_cols)
  
  dup <- duplicated(df[, c("Plate", "Row", "Column")])
  if (any(dup)) {
    stop("plot_plate_grid(): ", sum(dup), " duplicated plate/row/column combination(s); ",
         "expects at most one cell per well.", call. = FALSE)
  }
  df$present <- TRUE
  
  full_grid <- tidyr::expand_grid(
    Plate  = unique(df$Plate),
    Row    = factor(all_rows, levels = all_rows),
    Column = factor(all_cols, levels = all_cols)
  )
  
  tile_data <- full_grid %>%
    dplyr::left_join(df %>% dplyr::select(Plate, Row, Column, fill_var, discard_var, present),
                     by = c("Plate", "Row", "Column")) %>%
    dplyr::mutate(well_status = dplyr::case_when(
      is.na(present)                 ~ "not sorted",
      discard_var %in% discard_value ~ "discarded",
      TRUE                           ~ "kept"
    ))
  
  p <- ggplot2::ggplot() +
    ggplot2::geom_tile(
      data = dplyr::filter(tile_data, well_status != "kept"),
      ggplot2::aes(x = Column, y = Row, fill = well_status), color = "#FFFFFF"
    ) +
    ggplot2::scale_fill_manual(values = c("not sorted" = not_sorted_color, "discarded" = discarded_color),
                               name = NULL) +
    ggnewscale::new_scale_fill()
  
  if (type == "discrete") {
    # category order: factor levels, numeric order for numbers, else display order
    fill_levels <- if (is.factor(df$fill_var)) {
      intersect(levels(df$fill_var), as.character(df$fill_var))
    } else if (is.numeric(df$fill_var)) {
      as.character(sort(unique(df$fill_var)))
    } else {
      ordered_levels(df$fill_var, params$display_order[[fill_col]])
    }
    if (is.null(fill_values)) {
      status_colors <- resolve_palette("qc_status")
      fill_values <- if (all(fill_levels %in% names(status_colors))) {
        status_colors
      } else {
        resolve_group_colors(fill_col, fill_levels)
      }
    }
    if (is.null(names(fill_values))) fill_values <- setNames(fill_values[seq_along(fill_levels)], fill_levels)
    tile_data$fill_var <- factor(as.character(tile_data$fill_var), levels = fill_levels)
    p <- p +
      ggplot2::geom_tile(
        data = dplyr::filter(tile_data, well_status == "kept"),
        ggplot2::aes(x = Column, y = Row, fill = fill_var), color = "#FFFFFF"
      ) +
      ggplot2::scale_fill_manual(values = fill_values, breaks = fill_levels, name = fill_col)
  } else {
    p <- p +
      ggplot2::geom_tile(
        data = dplyr::filter(tile_data, well_status == "kept"),
        ggplot2::aes(x = Column, y = Row, fill = fill_var), color = "#FFFFFF"
      ) +
      ggplot2::scale_fill_viridis_c(option = resolve_palette("qc_metric"), limits = fill_limits,
                                    name = fill_col)
  }
  
  p +
    ggplot2::facet_wrap(~ Plate, axes = "all") +
    ggplot2::coord_fixed() +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(size = 7),
                   axis.text.y = ggplot2::element_text(size = 7)) +
    ggplot2::labs(title = wrap_label(title, width = title_wrap),
                  subtitle = wrap_label(subtitle, width = title_wrap))
}
#' One sweep metric across a two-parameter grid
#'
#' One tile per parameter combination, with its value printed in it. Every
#' other parameter must be constant in `df` (filter first), except
#' `facet_col`, which draws one panel per value. Text is white or dark
#' depending on each tile's fill.
#'
#' @param df Sweep summary, one row per parameter combination.
#' @param fill_col Metric to show.
#' @param x_col,y_col Parameters on the axes (e.g. "resolution"/"dims" for
#'   clustering, "d"/"k" for Milo).
#' @param title,subtitle Plot title/subtitle.
#' @param title_wrap Maximum characters per title/subtitle line. Default 50.
#' @param low_color,high_color Colors for the lowest and highest values.
#'   NULL (default) takes them from the "sweep_heatmap" palette; each can be
#'   set on its own.
#' @param limits Optional length-2 range of the color scale, e.g. c(0, 1)
#'   for ARI, so plots of the same metric share one scale. NULL (default)
#'   uses the range of the values.
#' @param facet_col Optional parameter with one panel per value.
#' @param text_threshold Luminance (0-1) below which tile text is white (see
#'   needs_white_text()). Default 0.5.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_sweep_tile(sweep_results, "mean_ari_stability", "resolution", "dims",
#'                 title = "Mean ARI stability", limits = c(0, 1))
#' }
plot_sweep_tile <- function(df,
                            fill_col,
                            x_col,
                            y_col,
                            title,
                            subtitle = NULL,
                            title_wrap = 50,
                            low_color = NULL,
                            high_color = NULL,
                            limits = NULL,
                            facet_col = NULL,
                            text_threshold = 0.5) {
  
  validate_required_cols(df, c(x_col, y_col, fill_col, facet_col), "plot_sweep_tile")
  keys <- df[, c(x_col, y_col, facet_col), drop = FALSE]
  if (anyDuplicated(keys)) {
    stop("plot_sweep_tile(): several rows share the same ", paste(c(x_col, y_col, facet_col), collapse = "/"),
         "; another parameter varies. Filter df to one value of it, or use facet_col.", call. = FALSE)
  }
  
  if (!is.null(limits) &&
      (!is.numeric(limits) || length(limits) != 2 ||
       any(!is.finite(limits)) || limits[1] >= limits[2])) {
    stop("plot_sweep_tile(): `limits` must be a finite increasing numeric length-2 vector.",
         call. = FALSE)
  }
  
  sweep_pal  <- resolve_palette("sweep_heatmap")
  low_color  <- low_color  %||% sweep_pal$low
  high_color <- high_color %||% sweep_pal$high
  
  values <- df[[fill_col]]
  # whole numbers without decimals, others to 2 significant digits
  labels <- if (all(is.na(values) | values == round(values))) {
    ifelse(is.na(values), "", format(round(values), trim = TRUE))
  } else {
    ifelse(is.na(values), "", formatC(values, format = "fg", digits = 2))
  }
  
  # same colors and range as the fill scale below, to choose the text color
  value_range <- limits %||% range(values, na.rm = TRUE)
  fill_colors <- scales::col_numeric(c(low_color, high_color), domain = value_range, na.color = "#E5E5E5")
  text_white  <- rep(FALSE, length(values))
  has_value   <- !is.na(values)
  text_white[has_value] <- needs_white_text(fill_colors(values[has_value]), threshold = text_threshold)
  
  plot_df <- data.frame(x = factor(df[[x_col]], levels = sort(unique(df[[x_col]]))),
                        y = factor(df[[y_col]], levels = sort(unique(df[[y_col]]))),
                        value = values, label = labels, text_white = text_white)
  if (!is.null(facet_col)) {
    facet_values <- df[[facet_col]]
    if (is.numeric(facet_values)) facet_values <- signif(facet_values, 3)
    plot_df$facet <- paste(facet_col, "=", facet_values)
  }
  
  p <- ggplot2::ggplot(plot_df, ggplot2::aes(x = x, y = y, fill = value)) +
    ggplot2::geom_tile(color = "#FFFFFF") +
    ggplot2::geom_text(ggplot2::aes(label = label, color = text_white), size = 3.2) +
    ggplot2::scale_color_manual(values = c("TRUE" = "#FFFFFF", "FALSE" = "#333333"), guide = "none") +
    ggplot2::scale_fill_gradient(low = low_color, high = high_color, limits = limits,
                                 name = fill_col, na.value = "#E5E5E5") +
    ggplot2::labs(title = wrap_label(title, width = title_wrap),
                  subtitle = wrap_label(subtitle, width = title_wrap), x = x_col, y = y_col) +
    ggplot2::theme_classic(base_size = 12)
  if (!is.null(facet_col)) p <- p + ggplot2::facet_wrap(~ facet)
  p
}
