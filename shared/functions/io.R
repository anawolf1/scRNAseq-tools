# ------------------------------------------------------------------------------
# io.R
#
# Shared save functions, used across every domain file: save_plot() for
# ggplot objects, save_gt() for gt tables, save_xlsx() for multi-sheet
# Excel workbooks.
#
# Functions
# save_plot()   generic ggplot saving, with a transparent-background default
# save_gt()     generic gt table saving, with a transparent-background default
# save_xlsx()   generic multi-sheet Excel saving, errors on sheet-name collision after cleaning
# ------------------------------------------------------------------------------

#' Save a ggplot to disk
#'
#' Wrapper around ggsave() with the project's default size and resolution
#' and a transparent background. Creates `output_dir` if needed.
#'
#' @param plot A ggplot object.
#' @param filename File name (with extension) to save as, e.g. "umap.png".
#' @param output_dir Directory to save into.
#' @param width Plot width in inches. Default 7.
#' @param height Plot height in inches. Default 6.
#' @param dpi Resolution. Default 300.
#' @param transparent Logical. If TRUE (default), saves with a transparent
#'   background instead of white. Overrides any panel/plot background set
#'   by the plot's own theme.
#'
#' @return Invisibly, the full path written to.
#'
#' @note Transparency depends on the output format/device. Formats such as
#'   PNG support alpha transparency; JPEG does not and will be flattened.
#'
#' @examples
#' \dontrun{
#' save_plot(p_umap, "umap_by_cluster.png", output_dir = fig_dir)
#' save_plot(p_umap, "umap_white_bg.png", output_dir = fig_dir, transparent = FALSE)
#' }
save_plot <- function(plot,
                      filename,
                      output_dir,
                      width = 7,
                      height = 6,
                      dpi = 300,
                      transparent = TRUE) {

  if (!inherits(plot, "ggplot")) {
    stop("save_plot(): `plot` must be a ggplot object, got: ", paste(class(plot), collapse = "/"),
         call. = FALSE)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (transparent) {
    plot <- plot + ggplot2::theme(
      panel.background      = ggplot2::element_rect(fill = "transparent", color = NA),
      plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.background     = ggplot2::element_rect(fill = "transparent", color = NA),
      legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA)
    )
  }

  out_path <- file.path(output_dir, filename)
  ggplot2::ggsave(
    filename = out_path,
    plot     = plot,
    width    = width,
    height   = height,
    dpi      = dpi,
    bg       = if (transparent) "transparent" else "#FFFFFF"
  )
  invisible(out_path)
}
#' Save a gt table to disk
#'
#' Wrapper around gt::gtsave() with the same interface as save_plot().
#' Creates `output_dir` if needed.
#'
#' @param gt_table A gt table object.
#' @param filename File name (with extension) to save as, e.g. "de_table.png".
#' @param output_dir Directory to save into.
#' @param transparent Logical. If TRUE (default), saves with a transparent
#'   table background instead of white.
#' @param vwidth Viewport width (pixels) for the PNG screenshot, passed to
#'   gtsave(). NULL (default) uses webshot2's default, which cuts off wide
#'   tables (many columns) on the right - pass e.g. 1800 for those.
#' @param vheight Viewport height (pixels), as `vwidth`. NULL uses the
#'   default; pass e.g. 2000 for very tall tables.
#'
#' @return Invisibly, the full path written to.
#'
#' @note PNG output is rendered by webshot2, which needs Chrome or Chromium
#'   installed. Transparency depends on the output format/device. Formats
#'   such as PNG support alpha transparency; JPEG does not and will be
#'   flattened.
#'
#' @examples
#' \dontrun{
#' save_gt(cluster_composition_gt, "cluster_composition.png", output_dir = fig_dir)
#'
#' # a wide table (many columns) needs a wider viewport to avoid truncation
#' save_gt(sweep_results_gt, "sweep_full.png", output_dir = fig_dir, vwidth = 1800)
#'
#' # a tall table (many rows) needs a taller viewport
#' save_gt(de_full_gt, "de_full.png", output_dir = fig_dir, vheight = 2000)
#' }
save_gt <- function(gt_table,
                    filename,
                    output_dir,
                    transparent = TRUE,
                    vwidth = NULL,
                    vheight = NULL) {

  if (!inherits(gt_table, "gt_tbl")) {
    stop("save_gt(): `gt_table` must be a gt table object, got: ", paste(class(gt_table), collapse = "/"),
         call. = FALSE)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  if (transparent) {
    gt_table <- gt_table %>% gt::tab_options(table.background.color = "transparent")
  }

  out_path <- file.path(output_dir, filename)
  gtsave_args <- list(data = gt_table, filename = out_path)
  if (!is.null(vwidth))  gtsave_args$vwidth  <- vwidth
  if (!is.null(vheight)) gtsave_args$vheight <- vheight
  do.call(gt::gtsave, gtsave_args)

  invisible(out_path)
}
#' Save a named list of data.frames as a multi-sheet Excel workbook to disk
#'
#' Wrapper around writexl::write_xlsx() with the same interface as
#' save_plot() and save_gt(). Characters Excel does not allow in sheet names
#' ([ ] : * ? / \) are replaced by `_`, and names are cut to 31 characters.
#' If two names collide after cleaning, save_xlsx() errors rather than
#' silently picking one.
#'
#' @param sheets Named list of data.frames, one per sheet. Names become
#'   sheet names (cleaned to valid Excel names if needed).
#' @param filename File name (with .xlsx extension) to save as.
#' @param output_dir Directory to save into.
#'
#' @return Invisibly, the full path written to.
#'
#' @examples
#' \dontrun{
#' save_xlsx(list(Summary = summary_df, Membership = membership_df),
#'           "deg_venn.xlsx", output_dir = fig_dir)
#' }
save_xlsx <- function(sheets,
                      filename,
                      output_dir) {

  # a data.frame is also a list, so it is rejected explicitly
  if (!is.list(sheets) || is.data.frame(sheets)) {
    stop("save_xlsx(): `sheets` must be a named list of data.frames, got: ",
         paste(class(sheets), collapse = "/"), call. = FALSE)
  }
  if (is.null(names(sheets)) || any(names(sheets) == "")) {
    stop("save_xlsx(): every element of `sheets` needs a name (used as its sheet name).", call. = FALSE)
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  sheet_names <- gsub("[\\[\\]:*?/\\\\]", "_", names(sheets))
  sheet_names <- substr(sheet_names, 1, 31)
  if (anyDuplicated(sheet_names)) {
    collided <- unique(sheet_names[duplicated(sheet_names) | duplicated(sheet_names, fromLast = TRUE)])
    stop(
      "save_xlsx(): sheet name(s) collide after cleaning to valid Excel names ",
      "(forbidden characters replaced by `_`, cut to 31 characters): ",
      paste(collided, collapse = ", "),
      ". Shorten or rename the original name(s) in `sheets` so they stay distinguishable.",
      call. = FALSE
    )
  }
  names(sheets) <- sheet_names

  out_path <- file.path(output_dir, filename)
  writexl::write_xlsx(sheets, out_path)
  invisible(out_path)
}
