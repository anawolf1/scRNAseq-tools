# ------------------------------------------------------------------------------
# tables_generic.R
#
# gt table builders used across workflows. None of them save to disk; call
# save_gt() separately, e.g. save_gt(tbl, "cellcounts.png", output_dir).
#
# Functions
# gt_generic()               any 2D table()/matrix/data.frame -> gt
# gt_contingency()           2-variable contingency table, with totals
# gt_contingency_faceted()   3-variable contingency table, one section per level of the third
# gt_de()                    FindMarkers()/FindAllMarkers() results -> gt
# p_to_stars()               p-values -> "***"/"**"/"*"/"ns" labels
# ------------------------------------------------------------------------------

#' Turn a table, matrix or data.frame into a gt table
#'
#' For summary tables in general (parameter sweeps, compositions, QC
#' summaries). A table or matrix gets its row names as the first column;
#' a data.frame or tibble is used as it is. A 1D table() is not supported.
#'
#' @param x A 2D table()/matrix, data.frame or tibble.
#' @param title Table title.
#' @param subtitle Optional table subtitle.
#' @param row_label Header of the row-name column, for table/matrix input.
#'   Default " ".
#'
#' @return A gt table object.
#'
#' @examples
#' \dontrun{
#' tbl <- gt_generic(sweep_results, title = "Clustering parameter sweep")
#' save_gt(tbl, "clustering_sweep.png", fig_dir)
#' }
gt_generic <- function(x,
                       title,
                       subtitle = NULL,
                       row_label = " ") {

  validate_table_dims(x, expected_dims = 2, fn_name = "gt_generic")

  # data.frames also have a 2D dim(), so they are told apart explicitly
  if (!is.data.frame(x)) {
    df <- as.data.frame.matrix(x)
    df <- tibble::rownames_to_column(df, row_label)
  } else {
    df <- as.data.frame(x)
  }

  gt::gt(df) %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::cols_align(align = "center")
}
#' Format a 2-variable contingency table as a gt table
#'
#' Rows and columns are the two variables (e.g. Patient x Timepoint), with
#' optional row and column totals.
#'
#' @param tbl A base table()/matrix of counts (not a data.frame - use
#'   gt_generic() for those).
#' @param title Table title.
#' @param subtitle Optional table subtitle.
#' @param row_label Header of the row-name column. Default " ".
#' @param total_row Logical. Add a totals row. Default TRUE.
#' @param total_col Logical. Add a totals column. Default TRUE. Set both
#'   FALSE for values that should not be summed (e.g. percent retained).
#'
#' @return A gt table object.
#'
#' @examples
#' \dontrun{
#' gt_contingency(table(meta$Patient, meta$Timepoint), title = "Cells per patient and timepoint")
#' }
gt_contingency <- function(tbl,
                           title,
                           subtitle = NULL,
                           row_label = " ",
                           total_row = TRUE,
                           total_col = TRUE) {

  validate_table_dims(tbl, expected_dims = 2, fn_name = "gt_contingency",
                      hint = "For a 3-variable faceted table, use gt_contingency_faceted() instead.")

  if (is.data.frame(tbl)) {
    stop("gt_contingency(): expects a base table()/matrix, not a data.frame. ",
         "For a data.frame, use gt_generic() instead.",
         call. = FALSE)
  }

  df <- as.data.frame.matrix(tbl)

  if (total_col) {
    df$Total <- rowSums(df, na.rm = TRUE)
  }

  if (total_row) {
    totals <- as.data.frame(as.list(colSums(df, na.rm = TRUE)))
    names(totals) <- names(df)
    rownames(totals) <- "Total"
    df <- rbind(df, totals)
  }

  df <- tibble::rownames_to_column(df, row_label)

  gt::gt(df) %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::cols_align(align = "center")
}
#' Format a 3-variable contingency table as one gt table
#'
#' Row x column counts, split by a third variable: one table with a section
#' per level of the third variable.
#'
#' @param tbl A 3-way table (row x column x group), e.g.
#'   table(row_values, col_values, group_values).
#' @param title Table title.
#' @param subtitle Optional table subtitle.
#' @param row_label Header of the row-variable column. Default " ".
#'
#' @return A gt table object.
#'
#' @examples
#' \dontrun{
#' gt_contingency_faceted(table(meta$Patient, meta$Timepoint, meta$Sex),
#'                        title = "Cells before QC", subtitle = "Patient x Timepoint, by Sex")
#' }
gt_contingency_faceted <- function(tbl,
                                   title,
                                   subtitle = NULL,
                                   row_label = " ") {

  validate_table_dims(tbl, expected_dims = 3, fn_name = "gt_contingency_faceted",
                      hint = "For a 2-variable table, use gt_contingency() instead.")

  df <- as.data.frame(tbl)
  names(df) <- c("row_var", "col_var", "group_var", "n")

  df_wide <- df %>%
    tidyr::pivot_wider(names_from = col_var, values_from = n, values_fill = 0)

  names(df_wide)[names(df_wide) == "row_var"] <- row_label

  gt::gt(df_wide, groupname_col = "group_var", rowname_col = row_label) %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::cols_align(align = "center")
}
#' Format differential expression results as a gt table
#'
#' Fixed decimals for log2FC and detection rates, scientific notation for
#' p-values, readable column labels, and log2FC colored with the
#' "de_log2fc" palette.
#'
#' @param de_df FindMarkers()/FindAllMarkers() results for one comparison,
#'   sorted in the order to show (see Details).
#' @param title Table title.
#' @param pct1_label,pct2_label Column labels for `pct.1` and `pct.2` (e.g.
#'   the names of the two groups).
#' @param subtitle Optional table subtitle.
#' @param n_genes Number of genes (rows) to show. Default 20.
#' @param exclude_gene_patterns Optional regex patterns; matching genes are
#'   removed before the top `n_genes` rows are taken (e.g. "^ENSG" for
#'   unannotated Ensembl IDs). Default NULL.
#' @param gene_from_rownames Logical. TRUE (default) when genes are row names
#'   (FindMarkers()); FALSE when they are in a `gene` column (FindAllMarkers()).
#' @param fc_bound Optional bound of the symmetric color scale,
#'   c(-fc_bound, fc_bound). NULL (default) uses the largest |log2FC| among
#'   the rows shown, so colors are not comparable between tables; pass the
#'   same value to several calls to make them comparable.
#'
#' @return A gt table object.
#'
#' @details
#' `de_df` is shown as given: it is not split by cluster or re-sorted, so
#' subset it to one comparison and sort it first.
#'
#' @examples
#' \dontrun{
#' gt_de(de_results, title = "Top DE genes", pct1_label = "Cured", pct2_label = "Chronic")
#'
#' # one table per cluster, on a shared color scale
#' shared_bound <- max(abs(de_by_cluster$avg_log2FC), na.rm = TRUE)
#' for (cl in unique(de_by_cluster$cluster)) {
#'   gt_de(de_by_cluster[de_by_cluster$cluster == cl, ], title = paste("Cluster", cl),
#'         pct1_label = "Cluster", pct2_label = "Other", gene_from_rownames = FALSE,
#'         fc_bound = shared_bound)
#' }
#' }
gt_de <- function(de_df,
                  title,
                  pct1_label,
                  pct2_label,
                  subtitle = NULL,
                  n_genes = 20,
                  exclude_gene_patterns = NULL,
                  gene_from_rownames = TRUE,
                  fc_bound = NULL) {

  validate_required_cols(de_df, c("avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj"), "gt_de",
                         hint = "expects a FindMarkers()/FindAllMarkers()-shaped data.frame.")

  if (!is.numeric(n_genes) || length(n_genes) != 1 || !is.finite(n_genes) ||
      n_genes < 1 || n_genes != as.integer(n_genes)) {
    stop("gt_de(): `n_genes` must be a single positive integer.", call. = FALSE)
  }
  if (!is.null(fc_bound) &&
      (!is.numeric(fc_bound) || length(fc_bound) != 1 || !is.finite(fc_bound) || fc_bound <= 0)) {
    stop("gt_de(): `fc_bound` must be a single positive, finite number, got: ",
         paste(fc_bound, collapse = ", "), call. = FALSE)
  }

  if (gene_from_rownames && !"gene" %in% names(de_df)) {
    de_df <- tibble::rownames_to_column(de_df, "gene")
  }
  if (!"gene" %in% names(de_df)) {
    stop("gt_de(): no 'gene' column and gene_from_rownames = FALSE. Set gene_from_rownames = TRUE ",
         "if genes are row names (FindMarkers()), or provide a 'gene' column (FindAllMarkers()).",
         call. = FALSE)
  }

  if (length(exclude_gene_patterns) > 0) {
    pattern <- paste(exclude_gene_patterns, collapse = "|")
    keep <- !is.na(de_df$gene) & !grepl(pattern, de_df$gene)
    de_df <- de_df[keep, , drop = FALSE]
  }

  de_df <- head(de_df, n_genes)

  if (nrow(de_df) == 0) {
    stop("gt_de(): de_df has no rows.", call. = FALSE)
  }

  fc_bound <- fc_bound %||% max(abs(de_df$avg_log2FC), na.rm = TRUE)
  if (!is.finite(fc_bound) || fc_bound == 0) {
    fc_bound <- 1
  }

  de_pal <- resolve_palette("de_log2fc")

  de_df %>%
    gt::gt() %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::cols_label(
      gene = "Gene", p_val = "p-value", avg_log2FC = "Log2FC",
      pct.1 = pct1_label, pct.2 = pct2_label, p_val_adj = "Adj. p-value"
    ) %>%
    gt::fmt_number(columns = c(avg_log2FC, pct.1, pct.2), decimals = 3) %>%
    gt::fmt_scientific(columns = c(p_val, p_val_adj), decimals = 2) %>%
    gt::data_color(
      columns = avg_log2FC,
      fn = scales::col_numeric(
        palette = c(de_pal$low, de_pal$mid, de_pal$high),
        domain = c(-fc_bound, fc_bound)
      )
    ) %>%
    gt::cols_align(align = "center")
}
#' p-values as significance stars
#'
#' "***" below 0.001, "**" below 0.01, "*" below 0.05, "ns" otherwise;
#' "" for NA.
#'
#' @param p Numeric vector of p-values.
#'
#' @return Character vector, one label per p-value.
#'
#' @examples
#' p_to_stars(c(0.0001, 0.005, 0.03, 0.2, NA))
p_to_stars <- function(p) {
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  )
}
