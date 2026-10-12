# ------------------------------------------------------------------------------
# milo_analysis.R
#
# Analyses, inputs, tests, summaries and plots of miloR differential
# abundance and miloDE neighborhood differential expression. None of these
# save to disk except plot_annotated_heatmap(), which writes its PNG
# (pheatmap draws to a file); call save_plot()/save_xlsx() for the rest.
#
# Functions
# resolve_milo_analyses()       validate params$analyses; expand "pairwise" into one analysis per contrast
# prepare_milo_input()          load the input, keep timepoints, add sample IDs, validate contrasts
# annotate_da_results()         add each neighborhood's majority label (+ "Mixed")
# da_label_summary()            neighborhoods tested / enriched / depleted per label
# nhood_label_composition()     cells and fraction per label within selected neighborhoods
# plot_nhood_composition()      stacked bars of neighborhood composition
# nhood_mean_scores()           mean of numeric cell-level columns per neighborhood
# group_da_nhoods()             merge overlapping significant neighborhoods with the same direction
# da_group_markers()            marker DE of a DA group vs other cells (Seurat FindMarkers)
# plot_marker_volcano()         volcano plot of da_group_markers() output (ggplot)
# plot_da_group_umap()          UMAP highlighting the cells of each DA group
# characterize_cell_set()       marker DE (vs rest, vs own cluster) + within-sample marker check for a cell set
# add_nhood_labels()            label selected neighborhoods on a miloR neighborhood-graph plot
# plot_nhood_program_heatmap()  z-scored mean program score, neighborhood x program
# plot_program_vs_da()          mean program score vs DA logFC, all neighborhoods
# milode_nhood_summary()        DE genes per miloDE neighborhood
# milode_gene_summary()         neighborhoods in which each gene is DE
# nhood_overlap_fraction()      share of each neighborhood's cells inside a set of other neighborhoods
# contrast_design()             sample table + model formula for one contrast (paired if repeated measures)
# nhood_sample_abundance()      per neighborhood x sample: cells and share of the sample's cells
# nhood_variability()           within-condition variability and top-sample share per neighborhood
# plot_nhood_variability()      logFC vs top-sample share and within-condition CV
# set_sample_abundance()        share of each sample's cells inside each cell set (e.g. meta-neighborhoods)
# plot_sample_abundance()       per-sample shares by timepoint, patients joined across repeated timepoints
# cell_da_calls()               per cell: in enriched / depleted / both / no significant neighborhood
# plot_nhood_correlations()     Spearman correlation between neighborhood-level variables
# group_nhoods_by_overlap()     neighborhood groups from ALL neighborhoods (Leiden on shared cells)
# assign_cells_to_groups()      one group per cell: the group most of its neighborhoods belong to
# group_resolution_stability()  how the grouping changes at nearby Leiden resolutions
# group_resolution_sweep()      bootstrap stability of the grouping across resolutions (tools/)
# partition_ari()               adjusted Rand index of two partitions
# test_set_abundance()          per-sample abundance test of non-overlapping cell sets (edgeR)
# set_composition()             cluster and timepoint composition of cell sets (raw and size-normalized)
# plot_set_composition()        stacked composition bars for many cell sets
# plot_annotated_heatmap()      net-program heatmap of neighborhoods or groups (pheatmap)
# plot_set_da_dots()            DA of groups / clusters as a dot plot, sets x contrasts
# ------------------------------------------------------------------------------

#' Validate and expand the analyses in config/decisions.R
#'
#' Each analysis builds one neighborhood graph on the cells of its
#' `timepoints` and tests its `contrasts` on it. Roles:
#'   - "primary" (exactly one): the main landscape; its contrasts are the
#'     study's contrasts
#'   - "sensitivity": another landscape (e.g. adding a timepoint); contrasts
#'     default to the primary's
#'   - "pairwise": expanded into one analysis per primary contrast, each on a
#'     graph of that contrast's two timepoints only (named
#'     "pairwise_<contrast>")
#' Optional `k` / `d` override params$da$k / params$da$d for that analysis.
#'
#' @param analyses params$analyses.
#'
#' @return Named list of analyses, each with role, timepoints, contrasts
#'   (and k, d if given), primary first.
resolve_milo_analyses <- function(analyses) {

  if (length(analyses) == 0 || is.null(names(analyses)) || any(!nzchar(names(analyses)))) {
    stop("resolve_milo_analyses(): params$analyses is empty or has unnamed entries - see config/decisions_example.R.", call. = FALSE)
  }
  roles <- vapply(analyses, function(a) a$role %||% "", character(1))
  bad_roles <- setdiff(roles, c("primary", "sensitivity", "pairwise"))
  if (length(bad_roles) > 0) {
    stop("resolve_milo_analyses(): params$analyses: role must be \"primary\", \"sensitivity\" or \"pairwise\", got: ",
         paste(bad_roles, collapse = ", "), call. = FALSE)
  }
  if (sum(roles == "primary") != 1) stop("resolve_milo_analyses(): params$analyses needs exactly one analysis with role = \"primary\".", call. = FALSE)

  primary <- analyses[[which(roles == "primary")]]
  if (length(primary$contrasts) == 0 || is.null(names(primary$contrasts))) {
    stop("resolve_milo_analyses(): the primary analysis needs named contrasts, each c(test, reference).", call. = FALSE)
  }

  out <- list()
  for (name in names(analyses)[order(roles != "primary")]) {
    a <- analyses[[name]]
    if (a$role == "pairwise") {
      for (ct in names(primary$contrasts)) {
        out[[paste0("pairwise_", ct)]] <- c(list(role = "pairwise", timepoints = primary$contrasts[[ct]],
                                                 contrasts = primary$contrasts[ct]), a[intersect(names(a), c("k", "d"))])
      }
      next
    }
    a$contrasts <- a$contrasts %||% primary$contrasts
    if (length(a$timepoints) == 0) stop("resolve_milo_analyses(): analysis '", name, "' needs `timepoints`.", call. = FALSE)
    for (ct in names(a$contrasts)) {
      if (!all(a$contrasts[[ct]] %in% a$timepoints)) {
        stop("resolve_milo_analyses(): analysis '", name, "': contrast '", ct, "' uses timepoints not in its `timepoints`.", call. = FALSE)
      }
    }
    out[[name]] <- a
  }
  out
}
#' Load and prepare the object shared by the Milo pipeline and tools
#'
#' Loads the input checkpoint (population refinement's final clustering by
#' default, with its cluster_label column), optionally keeps only the
#' timepoints in `timepoints`, adds `milo_sample` (patient x timepoint
#' replicate ID), and checks that every contrast names two timepoints
#' present in the data.
#'
#' @param input_path Path to the input checkpoint (.rds).
#' @param timepoints Timepoints to keep; NULL = all cells.
#' @param contrasts Named list, each c(test, reference) timepoints.
#' @param patient_col,timepoint_col Study-design columns.
#'
#' @return A list: `obj` (Seurat object), `contrasts` (validated),
#'   `sample_table` (samples and cells per timepoint).
prepare_milo_input <- function(input_path,
                               timepoints,
                               contrasts,
                               patient_col,
                               timepoint_col) {

  if (!file.exists(input_path)) {
    stop("prepare_milo_input(): input not found: ", input_path,
         ". Run population refinement first, or set params$input.", call. = FALSE)
  }
  obj <- readRDS(input_path)
  validate_required_cols(obj[[]], c(patient_col, timepoint_col, "cluster_label"), "prepare_milo_input",
                         hint = "expects the study-design columns and population refinement's cluster_label.")

  if (!is.null(timepoints)) {
    all_timepoints <- unique(as.character(obj[[timepoint_col, drop = TRUE]]))
    not_found <- setdiff(timepoints, all_timepoints)
    if (length(not_found) > 0) {
      stop("prepare_milo_input(): timepoint(s) not found in the data: ",
           paste(not_found, collapse = ", "), call. = FALSE)
    }
    keep <- as.character(obj[[timepoint_col, drop = TRUE]]) %in% timepoints
    obj <- subset(obj, cells = colnames(obj)[keep])
    message("Milo input restricted to ", paste(timepoints, collapse = ", "), ": ",
            sum(keep), " of ", length(keep), " cells")
  }

  obj$cluster_label <- droplevels(as.factor(obj$cluster_label))
  obj$milo_sample <- make_sample_id(obj[[]], c(patient_col, timepoint_col))

  if (length(contrasts) == 0 || is.null(names(contrasts)) || any(!nzchar(names(contrasts))) ||
      anyDuplicated(names(contrasts))) {
    stop("prepare_milo_input(): contrasts need unique names.", call. = FALSE)
  }
  present <- unique(as.character(obj[[timepoint_col, drop = TRUE]]))
  for (contrast_name in names(contrasts)) {
    levels_used <- contrasts[[contrast_name]]
    if (length(levels_used) != 2 || levels_used[1] == levels_used[2] || !all(levels_used %in% present)) {
      stop("prepare_milo_input(): contrast '", contrast_name, "' must be two different ", timepoint_col,
           " values present in the data (found: ", paste(sort(present), collapse = ", "), ").", call. = FALSE)
    }
  }

  sample_table <- obj[[]] %>%
    dplyr::group_by(.data[[timepoint_col]]) %>%
    dplyr::summarise(n_samples = dplyr::n_distinct(milo_sample),
                     n_patients = dplyr::n_distinct(.data[[patient_col]]),
                     n_cells = dplyr::n(),
                     min_cells_per_sample = min(table(milo_sample))) %>%
    dplyr::ungroup()

  list(obj = obj, contrasts = contrasts, sample_table = sample_table)
}
#' Add each neighborhood's majority label to testNhoods() results
#'
#' @param milo Milo object the results were computed on.
#' @param da_res testNhoods() output.
#' @param label_col colData column to annotate by (e.g. "cluster_label").
#' @param mixed_threshold Neighborhoods whose majority label covers less
#'   than this fraction of their cells are labeled "Mixed" (e.g.
#'   params$da$mixed_threshold).
#'
#' @return da_res with columns <label_col>, <label_col>_fraction and
#'   <label_col>_mixed (label, or "Mixed").
annotate_da_results <- function(milo,
                                da_res,
                                label_col,
                                mixed_threshold) {

  out <- miloR::annotateNhoods(milo, da_res, coldata_col = label_col)
  fraction_col <- paste0(label_col, "_fraction")
  mixed_col <- paste0(label_col, "_mixed")
  out[[mixed_col]] <- ifelse(out[[fraction_col]] < mixed_threshold, "Mixed", as.character(out[[label_col]]))
  out
}
#' Neighborhoods tested, enriched and depleted per label
#'
#' @param da_res annotate_da_results() output.
#' @param label_col Column to group by (e.g. "cluster_label_mixed").
#' @param alpha SpatialFDR threshold.
#'
#' @return A data.frame, one row per label.
da_label_summary <- function(da_res,
                             label_col,
                             alpha) {

  validate_required_cols(
    da_res,
    c(label_col, "logFC", "SpatialFDR"),
    "da_label_summary"
  )

  da_res %>%
    dplyr::group_by(.data[[label_col]]) %>%
    dplyr::summarise(
      n_nhoods      = dplyr::n(),
      n_significant = sum(SpatialFDR < alpha, na.rm = TRUE),
      n_enriched    = sum(SpatialFDR < alpha & logFC > 0, na.rm = TRUE),
      n_depleted    = sum(SpatialFDR < alpha & logFC < 0, na.rm = TRUE),
      median_logFC  = stats::median(logFC, na.rm = TRUE)
    ) %>%
    dplyr::ungroup()
}
#' Cells and fraction per label within selected neighborhoods
#'
#' @param membership Cell x neighborhood 0/1 matrix (miloR::nhoods()).
#' @param labels One label per cell (rows of `membership`), e.g.
#'   cluster_label or Timepoint.
#' @param nhood_ids Neighborhood indices to describe. NULL = all.
#'
#' @return Long data.frame: Nhood, label, n_cells, fraction (within the
#'   neighborhood). Labels absent from a neighborhood are omitted.
nhood_label_composition <- function(membership,
                                    labels,
                                    nhood_ids = NULL) {

  counts <- nhood_group_counts(membership, labels)
  if (is.null(nhood_ids)) nhood_ids <- seq_len(nrow(counts))
  counts <- counts[nhood_ids, , drop = FALSE]

  out <- data.frame(
    Nhood   = rep(nhood_ids, times = ncol(counts)),
    label   = rep(colnames(counts), each = nrow(counts)),
    n_cells = as.vector(counts)
  )
  out$fraction <- out$n_cells / stats::ave(out$n_cells, out$Nhood, FUN = sum)
  out[out$n_cells > 0, ]
}
#' Stacked bars of neighborhood composition
#'
#' @param composition nhood_label_composition() output.
#' @param nhood_order Neighborhood IDs in display order (bottom to top).
#' @param label_colors Optional named vector, label -> color.
#' @param nhood_text Optional named vector, Nhood (as character) -> axis label.
#' @param legend_title Legend title.
#' @param title,subtitle Plot title/subtitle.
#'
#' @return A ggplot object.
plot_nhood_composition <- function(composition,
                                   nhood_order,
                                   label_colors = NULL,
                                   nhood_text = NULL,
                                   legend_title = "label",
                                   title = NULL,
                                   subtitle = NULL) {

  validate_required_cols(composition, c("Nhood", "label", "fraction"), "plot_nhood_composition")

  row_ids <- as.character(nhood_order)
  composition$row <- factor(as.character(composition$Nhood), levels = row_ids,
                            labels = if (is.null(nhood_text)) row_ids else unname(nhood_text[row_ids]))
  if (!is.null(label_colors)) {
    composition$label <- factor(composition$label, levels = intersect(names(label_colors), unique(composition$label)))
  }

  p <- ggplot2::ggplot(composition, ggplot2::aes(x = fraction, y = row, fill = label)) +
    ggplot2::geom_col(width = 0.85) +
    ggplot2::scale_x_continuous(labels = scales::percent, expand = c(0, 0)) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Cells in neighborhood", y = "Neighborhood",
                  fill = legend_title) +
    ggplot2::theme_classic(base_size = 11)

  if (!is.null(label_colors)) p <- p + ggplot2::scale_fill_manual(values = label_colors)
  p
}
#' Mean of numeric cell-level columns per neighborhood
#'
#' @param membership Cell x neighborhood 0/1 matrix (miloR::nhoods()).
#' @param values Numeric data.frame or matrix, rows named by cell, e.g.
#'   module-score columns.
#'
#' @return A data.frame: Nhood, n_cells_scored, then one column per input
#'   column. Cells missing from `values` are ignored.
nhood_mean_scores <- function(membership,
                              values) {

  values <- as.matrix(values)
  shared <- intersect(rownames(membership), rownames(values))
  if (length(shared) == 0) {
    stop("nhood_mean_scores(): no cell names shared between the neighborhoods and `values`.", call. = FALSE)
  }

  m <- membership[shared, , drop = FALSE]
  v <- values[shared, , drop = FALSE]
  n <- Matrix::colSums(m)
  sums <- as.matrix(Matrix::crossprod(m, v))

  data.frame(Nhood = seq_len(ncol(m)), n_cells_scored = as.numeric(n),
             sweep(sums, 1, pmax(n, 1), "/"), check.names = FALSE)
}
#' Heatmap of mean program score per neighborhood, z-scored within program
#'
#' Programs have different raw ranges, so each is z-scored across the
#' neighborhoods shown: colors compare neighborhoods within a program,
#' not programs with each other.
#'
#' @param score_df nhood_mean_scores() output, optionally with extra
#'   columns (ignored unless listed in `score_cols`).
#' @param score_cols Score columns to show.
#' @param row_label Optional named vector, Nhood (as character) -> row label.
#' @param row_order Optional Nhood IDs in display order, top to bottom
#'   (e.g. by DA logFC). NULL = by Nhood number.
#' @param title,subtitle Plot title/subtitle.
#'
#' @return A ggplot object.
plot_nhood_program_heatmap <- function(score_df,
                                       score_cols,
                                       row_label = NULL,
                                       row_order = NULL,
                                       title = NULL,
                                       subtitle = NULL) {

  validate_required_cols(score_df, c("Nhood", score_cols), "plot_nhood_program_heatmap")

  long_df <- score_df[, c("Nhood", score_cols)] %>%
    tidyr::pivot_longer(-Nhood, names_to = "score_col", values_to = "score") %>%
    dplyr::group_by(score_col) %>%
    dplyr::mutate(z = if (isTRUE(stats::sd(score, na.rm = TRUE) > 0)) as.numeric(scale(score)) else 0) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(program = sub("Score$", "", score_col))

  long_df$program <- factor(long_df$program, levels = order_programs_up_dn(unique(long_df$program)))
  row_ids <- if (is.null(row_order)) as.character(sort(unique(long_df$Nhood))) else rev(as.character(row_order))
  row_text <- if (is.null(row_label)) row_ids else paste0(row_ids, " (", row_label[row_ids], ")")
  long_df$row <- factor(as.character(long_df$Nhood), levels = row_ids, labels = row_text)

  pal <- resolve_palette("module_score_z")

  ggplot2::ggplot(long_df, ggplot2::aes(x = program, y = row, fill = z)) +
    ggplot2::geom_tile(color = "#FFFFFF", linewidth = 0.3) +
    ggplot2::scale_fill_gradient2(low = pal$low, mid = pal$mid, high = pal$high, midpoint = 0,
                                  name = "z-scored\nmean score") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = "Neighborhood") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
#' Mean program score against DA logFC, across all neighborhoods
#'
#' One panel per program; each point is a neighborhood, colored by
#' whether it is DA-significant. Descriptive only: neighborhoods overlap
#' and share cells, so they are not independent observations and no
#' correlation test is reported.
#'
#' @param da_res testNhoods() results joined with nhood_mean_scores().
#' @param score_cols Score columns to plot.
#' @param alpha SpatialFDR threshold used for coloring.
#' @param test_label,ref_label Contrast levels, for the axis title.
#' @param title,subtitle Plot title/subtitle.
#'
#' @return A ggplot object.
plot_program_vs_da <- function(da_res,
                               score_cols,
                               alpha,
                               test_label,
                               ref_label,
                               title = NULL,
                               subtitle = NULL,
                               y_label = "Mean module score") {

  validate_required_cols(da_res, c("Nhood", "logFC", "SpatialFDR", score_cols), "plot_program_vs_da")

  long_df <- da_res[, c("Nhood", "logFC", "SpatialFDR", score_cols)] %>%
    tidyr::pivot_longer(dplyr::all_of(score_cols), names_to = "score_col", values_to = "score") %>%
    dplyr::mutate(
      program = factor(sub("Score$", "", score_col),
                       levels = order_programs_up_dn(unique(sub("Score$", "", score_cols)))),
      DA = dplyr::case_when(
        SpatialFDR < alpha & logFC > 0 ~ paste("Enriched in", test_label),
        SpatialFDR < alpha & logFC < 0 ~ paste("Enriched in", ref_label),
        TRUE ~ "Not significant"
      )
    )
  long_df$DA <- factor(long_df$DA, levels = c(paste("Enriched in", ref_label), "Not significant",
                                              paste("Enriched in", test_label)))
  pal <- resolve_palette("de_log2fc")

  ggplot2::ggplot(long_df, ggplot2::aes(x = logFC, y = score, color = DA)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#999999") +
    ggplot2::geom_point(size = 1.8, alpha = 0.8) +
    ggplot2::scale_color_manual(values = setNames(c(pal$low, "#BFBFBF", pal$high), levels(long_df$DA)),
                                drop = FALSE, name = NULL) +
    ggplot2::facet_wrap(~ program, scales = "free_y") +
    ggplot2::labs(title = title, subtitle = subtitle,
                  x = paste0("DA logFC (", test_label, " vs ", ref_label, ")"), y = y_label) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "bottom", strip.text = ggplot2::element_text(face = "bold", size = 8))
}
#' Merge overlapping significant neighborhoods with the same direction
#'
#' Same idea as miloR::groupNhoods(), written out so the rule is explicit
#' and does not change between miloR versions: two significant
#' neighborhoods are linked when they change in the same direction and
#' share at least `min_overlap` of the smaller one's cells; groups are the
#' connected components of those links. Groups are numbered by cell count.
#'
#' @param membership Cell x neighborhood 0/1 matrix (miloR::nhoods()).
#' @param da_res testNhoods() output (Nhood, logFC, SpatialFDR).
#' @param alpha SpatialFDR threshold.
#' @param min_overlap Shared cells as a fraction of the smaller
#'   neighborhood. Default 0.3.
#'
#' @return A data.frame, one row per significant neighborhood: Nhood,
#'   da_group ("G1", "G2", ...), direction ("enriched"/"depleted"). Zero
#'   rows if nothing is significant.
group_da_nhoods <- function(membership,
                            da_res,
                            alpha,
                            min_overlap = 0.3) {

  sig <- da_res[which(da_res$SpatialFDR < alpha), c("Nhood", "logFC")]
  if (nrow(sig) == 0) {
    return(data.frame(Nhood = integer(0), da_group = character(0), direction = character(0)))
  }

  m <- membership[, sig$Nhood, drop = FALSE]
  shared <- as.matrix(Matrix::crossprod(m))
  sizes <- diag(shared)
  overlap_frac <- shared / outer(sizes, sizes, pmin)
  same_sign <- outer(sign(sig$logFC), sign(sig$logFC), "==")
  adjacency <- (overlap_frac >= min_overlap) & same_sign
  diag(adjacency) <- FALSE

  g <- igraph::graph_from_adjacency_matrix(adjacency * 1, mode = "undirected")
  component <- igraph::components(g)$membership

  group_cells <- vapply(split(seq_along(component), component), function(idx) {
    sum(Matrix::rowSums(m[, idx, drop = FALSE]) > 0)
  }, numeric(1))
  rank <- setNames(seq_along(group_cells), names(sort(group_cells, decreasing = TRUE)))

  data.frame(
    Nhood     = sig$Nhood,
    da_group  = paste0("G", rank[as.character(component)]),
    direction = ifelse(sig$logFC > 0, "enriched", "depleted")
  )
}
#' Marker DE of a set of cells against other cells
#'
#' Seurat::FindMarkers() with explicit cell vectors. The cells were
#' selected in the same expression space being tested, so p-values rank
#' genes but are not calibrated tests (descriptive use only).
#'
#' @param obj Seurat object.
#' @param cells_1,cells_2 Cell names for the group and the comparison.
#' @param log2fc_cutoff,min_pct FindMarkers() logfc.threshold and min.pct.
#' @param exclude_gene_patterns Regular expressions for genes dropped from
#'   the result (e.g. params$da$markers$exclude_gene_patterns).
#'
#' @return FindMarkers() output with a `gene` column, sorted by adjusted
#'   p-value; zero rows if nothing passes.
da_group_markers <- function(obj,
                             cells_1,
                             cells_2,
                             log2fc_cutoff = 0.25,
                             min_pct = 0.1,
                             exclude_gene_patterns = NULL) {

  de <- tryCatch(
    Seurat::FindMarkers(obj, ident.1 = cells_1, ident.2 = cells_2,
                        logfc.threshold = log2fc_cutoff, min.pct = min_pct, verbose = FALSE),
    error = function(e) {
      message("da_group_markers(): FindMarkers() failed: ", conditionMessage(e))
      NULL
    }
  )
  if (is.null(de) || nrow(de) == 0) {
    return(data.frame(p_val = numeric(0), avg_log2FC = numeric(0), pct.1 = numeric(0),
                      pct.2 = numeric(0), p_val_adj = numeric(0), gene = character(0)))
  }
  de$gene <- rownames(de)
  if (length(exclude_gene_patterns) > 0) {
    de <- de[!grepl(paste(exclude_gene_patterns, collapse = "|"), de$gene), , drop = FALSE]
  }
  de[order(de$p_val_adj, -abs(de$avg_log2FC)), , drop = FALSE]
}
#' Volcano plot of marker DE results
#'
#' ggplot volcano (no EnhancedVolcano dependency in this workflow).
#'
#' @param de da_group_markers() output.
#' @param log2fc_cutoff,padj_cutoff Thresholds for highlighting
#'   (|avg_log2FC|, adjusted p-value). Defaults 0.5, 0.05.
#' @param n_labels Number of top genes labeled (by adjusted p, then |log2FC|).
#' @param title,subtitle Plot title/subtitle.
#'
#' @return A ggplot object.
plot_marker_volcano <- function(de,
                                log2fc_cutoff = 0.5,
                                padj_cutoff = 0.05,
                                n_labels = 20,
                                title = NULL,
                                subtitle = NULL) {

  validate_required_cols(de, c("gene", "avg_log2FC", "p_val_adj"), "plot_marker_volcano")
  pal <- resolve_palette("de_log2fc")

  de$neg_log10_p <- -log10(pmax(de$p_val_adj, .Machine$double.xmin))
  de$status <- ifelse(de$p_val_adj < padj_cutoff & de$avg_log2FC >= log2fc_cutoff, "up",
                      ifelse(de$p_val_adj < padj_cutoff & de$avg_log2FC <= -log2fc_cutoff, "down", "ns"))
  label_df <- utils::head(de[de$status != "ns", ], n_labels)

  p <- ggplot2::ggplot(de, ggplot2::aes(x = avg_log2FC, y = neg_log10_p, color = status)) +
    ggplot2::geom_point(size = 1.4, alpha = 0.75) +
    ggplot2::geom_vline(xintercept = c(-log2fc_cutoff, log2fc_cutoff), linetype = "dashed", color = "#999999") +
    ggplot2::geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "#999999") +
    ggplot2::scale_color_manual(values = c(down = pal$low, ns = "#BFBFBF", up = pal$high), guide = "none") +
    ggplot2::labs(title = title, subtitle = subtitle, x = "avg log2FC", y = "-log10 adjusted p") +
    ggplot2::theme_classic(base_size = 12)

  if (nrow(label_df) > 0) {
    p <- p + ggrepel::geom_text_repel(data = label_df, ggplot2::aes(label = gene), size = 3,
                                      color = "#000000", max.overlaps = Inf, min.segment.length = 0)
  }
  p
}
#' UMAP highlighting the cells of each DA group
#'
#' @param obj Seurat object with a "umap" reduction.
#' @param group_cells Named list, group label -> cell names.
#' @param title,subtitle Plot title/subtitle.
#'
#' @return A ggplot object, one panel per group (group cells colored,
#'   all other cells grey).
plot_da_group_umap <- function(obj,
                               group_cells,
                               title = NULL,
                               subtitle = NULL) {

  embed <- as.data.frame(Seurat::Embeddings(obj, "umap")[, 1:2])
  colnames(embed) <- c("Dim1", "Dim2")
  embed$cell <- rownames(embed)

  df <- dplyr::bind_rows(lapply(names(group_cells), function(g) {
    d <- embed
    d$panel <- g
    d$in_group <- d$cell %in% group_cells[[g]]
    d[order(d$in_group), ]
  }))
  df$panel <- factor(df$panel, levels = names(group_cells))

  ggplot2::ggplot(df, ggplot2::aes(x = Dim1, y = Dim2, color = in_group)) +
    ggplot2::geom_point(size = 0.6) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "#D9D9D9", `TRUE` = "#B22222"), guide = "none") +
    ggplot2::facet_wrap(~ panel) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "umap_1", y = "umap_2") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                   strip.text = ggplot2::element_text(size = 8))
}
#' Marker DE and within-sample marker check for one set of cells
#'
#' Used for meta-neighborhoods and for hand-picked neighborhoods. Runs marker DE of
#' the set against (1) all other cells and (2) the other cells of its
#' dominant cluster, saves a volcano plot and top-marker table for each, and
#' checks the set's top markers within each sample (patient x timepoint):
#' mean expression in the sample's cells inside the set vs its other cells.
#' A state marker should be higher inside the set in most samples; a marker
#' carried by one patient rises in one sample only (within_sample_shift(),
#' within_sample.R).
#'
#' @param obj Seurat object (cells = the Milo input).
#' @param cells Cell names of the set.
#' @param set_name Short name used in file names (e.g. "G1", "nhood86").
#' @param title Plot/table title.
#' @param out_dir Output directory.
#' @param file_suffix Appended to file names (e.g. the contrast), or "".
#' @param cluster_col,sample_col,timepoint_col Metadata columns.
#' @param timepoint_colors Named colors for the timepoints.
#' @param marker_params params$da$markers (min_cells, log2fc_cutoff,
#'   min_pct, padj_cutoff, top_n, within_sample_min_cells).
#' @param exclude_gene_patterns Gene patterns dropped from DE results.
#'
#' @return Named list of data.frames: vs_rest, vs_cluster (when run),
#'   within_sample, within_sample_summary (when markers exist).
characterize_cell_set <- function(obj,
                                  cells,
                                  set_name,
                                  title,
                                  out_dir,
                                  file_suffix = "",
                                  cluster_col = "cluster_label",
                                  sample_col = "milo_sample",
                                  timepoint_col = "Timepoint",
                                  timepoint_colors,
                                  marker_params,
                                  exclude_gene_patterns = NULL) {

  min_cells <- marker_params$min_cells
  sheets <- list()
  tag <- function(prefix, extra = NULL) {
    paste0(paste(c(prefix, set_name, extra, if (nzchar(file_suffix)) file_suffix), collapse = "_"), ".png")
  }

  if (length(cells) < min_cells) {
    message("  ", set_name, ": ", length(cells), " cells - below markers$min_cells (", min_cells, "), skipped")
    return(sheets)
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  meta <- obj[[]]
  dominant <- names(which.max(table(as.character(meta[cells, cluster_col]))))

  comparisons <- list(
    vs_rest    = list(cells_2 = setdiff(colnames(obj), cells), label = "vs all other cells"),
    vs_cluster = list(cells_2 = setdiff(colnames(obj)[as.character(meta[[cluster_col]]) == dominant], cells),
                      label = paste0("vs other ", dominant, " cells"))
  )

  for (cmp in names(comparisons)) {
    cells_2 <- comparisons[[cmp]]$cells_2
    if (length(cells_2) < min_cells) {
      message("  ", set_name, " ", cmp, ": comparison has ", length(cells_2), " cells - skipped")
      next
    }
    de <- da_group_markers(obj, cells, cells_2,
                           log2fc_cutoff = marker_params$log2fc_cutoff,
                           min_pct = marker_params$min_pct,
                           exclude_gene_patterns = exclude_gene_patterns)
    sheets[[cmp]] <- de
    if (nrow(de) == 0) next

    p_volcano <- plot_marker_volcano(
      de, log2fc_cutoff = marker_params$log2fc_cutoff, padj_cutoff = marker_params$padj_cutoff, title = title,
      subtitle = paste0(comparisons[[cmp]]$label, "; descriptive - cells were selected in the same expression space")
    )
    print(p_volcano)
    save_plot(p_volcano, tag("volcano", cmp), out_dir, width = 8, height = 6.5)

    de_gt <- gt_de(de, title = paste(title, comparisons[[cmp]]$label),
                   pct1_label = paste0("% ", set_name), pct2_label = "% comparison",
                   subtitle = "Top markers by adjusted p-value (descriptive)", n_genes = 20,
                   gene_from_rownames = FALSE)
    save_gt(de_gt, tag("markers", cmp), out_dir)
  }

  de_rest <- sheets$vs_rest
  if (!is.null(de_rest) && nrow(de_rest) > 0) {
    top_up <- utils::head(de_rest$gene[de_rest$avg_log2FC > 0 & de_rest$p_val_adj < marker_params$padj_cutoff],
                          marker_params$top_n)
    if (length(top_up) > 0) {
      expr <- Seurat::FetchData(obj, vars = top_up)
      check_meta <- cbind(meta[, c(sample_col, timepoint_col)], expr)
      shift <- within_sample_shift(check_meta, cells, sample_col, top_up,
                                   min_cells = marker_params$within_sample_min_cells, summary_fun = mean)
      if (nrow(shift) > 0) {
        shift$timepoint <- as.character(check_meta[[timepoint_col]])[match(shift$sample, check_meta[[sample_col]])]
        sheets$within_sample <- shift
        sheets$within_sample_summary <- summarize_shift(shift)

        p_within <- plot_within_sample_shift(
          shift, timepoint_colors, inside_label = set_name,
          y_label = "Mean log-normalized expression",
          title = paste0(title, ": top markers within each sample"),
          subtitle = paste0("Each line = one sample (patient x timepoint): its cells in ", set_name,
                            " vs its other cells; a state marker rises in most samples")
        )
        print(p_within)
        save_plot(p_within, tag("within_sample"), out_dir,
                  width = min(15, 3.5 * min(length(top_up), 3) + 2), height = 3.5 * ceiling(length(top_up) / 3) + 1.5)
      }
    }
  }
  sheets
}
#' Label selected neighborhoods on a miloR neighborhood-graph plot
#'
#' plotNhoodGraph()/plotNhoodGraphDA() reorder vertices for drawing, so rows
#' of the plot data are not in neighborhood order. Each vertex is named by
#' its index cell, so neighborhood k is the vertex named nhoodIndex(milo)[[k]].
#'
#' @param p ggraph plot from plotNhoodGraph() / plotNhoodGraphDA().
#' @param milo The Milo object the plot was drawn from.
#' @param nhood_ids Neighborhood numbers to label.
#' @param labels Text per neighborhood (default: its number).
#' @param size Label text size.
#'
#' @return The plot with repelled labels added.
add_nhood_labels <- function(p,
                             milo,
                             nhood_ids,
                             labels = as.character(nhood_ids),
                             size = 3) {

  if (length(nhood_ids) == 0) return(p)
  index_cells <- as.character(unlist(miloR::nhoodIndex(milo)))
  rows <- match(index_cells[nhood_ids], as.character(p$data$name))
  if (anyNA(rows)) stop("add_nhood_labels(): some neighborhoods are not in the plot data.", call. = FALSE)

  label_df <- data.frame(x = p$data$x[rows], y = p$data$y[rows], label = labels)
  p + ggrepel::geom_label_repel(data = label_df, ggplot2::aes(x = x, y = y, label = label),
                                inherit.aes = FALSE, size = size, min.segment.length = 0,
                                max.overlaps = Inf, label.size = 0.15, fill = "#FFFFFF", alpha = 0.85)
}
#' DE genes per miloDE neighborhood
#'
#' @param de_stat miloDE::de_test_neighbourhoods() output (data.frame).
#' @param alpha Threshold on pval_corrected_across_genes (within-neighborhood
#'   BH), as used by miloDE::rank_neighbourhoods_by_DE_magnitude().
#'
#' @return A data.frame, one row per neighborhood: Nhood, Nhood_center,
#'   test_performed, n_genes_tested, n_DE_genes, n_specific_DE_genes.
milode_nhood_summary <- function(de_stat,
                                 alpha) {

  ranked <- miloDE::rank_neighbourhoods_by_DE_magnitude(de_stat, pval.thresh = alpha)
  tested <- de_stat %>%
    dplyr::group_by(Nhood) %>%
    dplyr::summarise(test_performed = any(test_performed),
                     n_genes_tested = sum(test_performed & !is.na(pval)), .groups = "drop")

  ranked[, c("Nhood", "Nhood_center", "n_DE_genes", "n_specific_DE_genes")] %>%
    dplyr::left_join(tested, by = "Nhood") %>%
    dplyr::relocate(test_performed, n_genes_tested, .after = Nhood_center)
}
#' In how many neighborhoods is each gene DE
#'
#' A gene counts as DE in a neighborhood when its within-neighborhood
#' BH-adjusted p-value (pval_corrected_across_genes, which controls the FDR
#' across the ~15,000 genes tested there) is below alpha - the same rule as
#' n_DE_genes in milode_nhood_summary(). miloDE's pval_corrected_across_nhoods
#' only corrects each gene across neighborhoods, not across genes, so on its
#' own it is suited to asking about one pre-chosen gene, not to genome-wide
#' screening; it is reported as n_nhoods_spatial for reference.
#'
#' @param de_stat miloDE::de_test_neighbourhoods() output (data.frame).
#' @param alpha FDR threshold.
#' @param extreme_logfc |logFC| above which a result is flagged: in
#'   pseudobulk DE this almost always means the gene is essentially absent
#'   from one group's samples (on/off), which often reflects individual
#'   (e.g. genotype) or technical differences rather than cell state.
#'
#' @return A data.frame, one row per gene DE in at least one neighborhood:
#'   n_nhoods_tested, n_nhoods_DE, n_up, n_down, mean_logFC_DE,
#'   pct_extreme_logFC, n_nhoods_spatial; sorted by n_nhoods_DE.
milode_gene_summary <- function(de_stat,
                                alpha,
                                extreme_logfc = 5) {

  de_stat %>%
    dplyr::filter(test_performed, !is.na(pval)) %>%
    dplyr::group_by(gene) %>%
    dplyr::summarise(
      n_nhoods_tested   = dplyr::n(),
      n_nhoods_DE       = sum(pval_corrected_across_genes < alpha, na.rm = TRUE),
      n_up              = sum(pval_corrected_across_genes < alpha & logFC > 0, na.rm = TRUE),
      n_down            = sum(pval_corrected_across_genes < alpha & logFC < 0, na.rm = TRUE),
      mean_logFC_DE     = if (any(pval_corrected_across_genes < alpha, na.rm = TRUE)) {
        mean(logFC[which(pval_corrected_across_genes < alpha)])
      } else {
        NA_real_
      },
      pct_extreme_logFC = if (any(pval_corrected_across_genes < alpha, na.rm = TRUE)) {
        round(100 * mean(abs(logFC[which(pval_corrected_across_genes < alpha)]) > extreme_logfc), 1)
      } else {
        NA_real_
      },
      n_nhoods_spatial  = sum(pval_corrected_across_nhoods < alpha, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_nhoods_DE > 0) %>%
    dplyr::arrange(dplyr::desc(n_nhoods_DE), pct_extreme_logFC, dplyr::desc(abs(mean_logFC_DE)))
}
#' Share of each neighborhood's cells lying inside a set of other neighborhoods
#'
#' Links two neighborhood sets built on different graphs (e.g. miloDE
#' neighborhoods and DA-significant miloR neighborhoods) through the
#' cells they share.
#'
#' @param membership Cell x neighborhood matrix of the neighborhoods to
#'   describe (rows named by cell).
#' @param reference_membership Cell x neighborhood matrix of the reference
#'   set (rows named by cell).
#' @param reference_nhoods Column indices of `reference_membership` to use
#'   (e.g. DA-significant neighborhoods).
#'
#' @return Numeric vector, one value per column of `membership`: fraction
#'   of its cells that belong to at least one selected reference
#'   neighborhood (0 if none are selected).
nhood_overlap_fraction <- function(membership,
                                   reference_membership,
                                   reference_nhoods) {

  if (length(reference_nhoods) == 0) return(rep(0, ncol(membership)))

  in_reference <- rownames(reference_membership)[
    Matrix::rowSums(reference_membership[, reference_nhoods, drop = FALSE]) > 0
  ]
  inside <- as.numeric(rownames(membership) %in% in_reference)
  as.numeric(Matrix::crossprod(membership, inside)) / pmax(Matrix::colSums(membership), 1)
}
#' Sample table and model formula for one contrast
#'
#' Unpaired (~ condition) when the two timepoints come from different
#' patients. Paired (~ pair + condition) when both are repeated-measures
#' timepoints: only patients sampled at both are kept, and each patient's
#' own level is absorbed by `pair`, so the condition effect is estimated
#' within patients. The tested coefficient is the last one (condition, test
#' vs reference) in both cases.
#'
#' @param sample_design Data.frame, one row per sample, with `milo_sample`,
#'   the patient and the timepoint column.
#' @param levels_used c(test, reference).
#' @param patient_col,timepoint_col Study-design columns.
#' @param repeated_timepoints Repeated-measures timepoints (study design).
#' @param paired Use a paired design when possible. Default TRUE.
#'
#' @return List: `design_df` (rows = samples in the test), `formula`,
#'   `paired` (logical), `dropped` (samples left out: unpaired patients).
contrast_design <- function(sample_design,
                            levels_used,
                            patient_col,
                            timepoint_col,
                            repeated_timepoints,
                            paired = TRUE) {

  design_df <- sample_design[sample_design[[timepoint_col]] %in% levels_used, , drop = FALSE]
  design_df$condition <- factor(as.character(design_df[[timepoint_col]]), levels = rev(levels_used))
  is_paired <- isTRUE(paired) && all(levels_used %in% repeated_timepoints)
  dropped <- character(0)

  if (is_paired) {
    n_levels <- tapply(as.character(design_df[[timepoint_col]]), design_df[[patient_col]],
                       function(x) length(unique(x)))
    complete <- names(n_levels)[n_levels == 2]
    dropped <- design_df$milo_sample[!design_df[[patient_col]] %in% complete]
    design_df <- design_df[design_df[[patient_col]] %in% complete, , drop = FALSE]
    if (length(complete) < 2) {
      stop("contrast_design(): fewer than 2 patients sampled at both ", paste(levels_used, collapse = " and "),
           " - a paired test is not possible.", call. = FALSE)
    }
    design_df$pair <- factor(as.character(design_df[[patient_col]]))
  }

  list(design_df = design_df,
       formula = if (is_paired) ~ pair + condition else ~ condition,
       paired = is_paired,
       dropped = dropped)
}
#' Cells and share of each sample's cells, per neighborhood and sample
#'
#' The share (cells of the sample in the neighborhood / all cells of the
#' sample in the analysis) is the quantity the DA test compares between
#' conditions, before its TMM normalization. Neighborhoods overlap, so a
#' sample's shares do not sum to 1.
#'
#' @param milo Milo object after countCells().
#' @param design_df Samples to include (from contrast_design()), with
#'   `milo_sample` and `condition`.
#' @param patient_col Patient column in `design_df`.
#' @param nhoods Neighborhood IDs; NULL = all.
#'
#' @return Long data.frame: Nhood, milo_sample, condition, patient, n_cells,
#'   sample_total, share.
nhood_sample_abundance <- function(milo,
                                   design_df,
                                   patient_col,
                                   nhoods = NULL) {

  counts <- as.matrix(miloR::nhoodCounts(milo))
  nhoods <- nhoods %||% seq_len(nrow(counts))
  samples <- intersect(design_df$milo_sample, colnames(counts))
  totals <- table(as.character(SummarizedExperiment::colData(milo)$milo_sample))[samples]

  sub <- counts[nhoods, samples, drop = FALSE]
  data.frame(
    Nhood        = rep(nhoods, times = length(samples)),
    milo_sample  = rep(samples, each = length(nhoods)),
    n_cells      = as.vector(sub),
    sample_total = rep(as.numeric(totals), each = length(nhoods))
  ) %>%
    dplyr::mutate(share = n_cells / sample_total) %>%
    dplyr::left_join(design_df[, c("milo_sample", "condition", patient_col)], by = "milo_sample") %>%
    dplyr::rename(patient = dplyr::all_of(patient_col))
}
#' Within-condition variability and top-sample share per neighborhood
#'
#' Two reasons a neighborhood with a large fold change can stay
#' non-significant: its abundance varies a lot between samples of the same
#' condition (high coefficient of variation), or most of its cells come from
#' a single sample (high top-sample share) - one patient's peculiarity
#' rather than a condition effect.
#'
#' @param abundance nhood_sample_abundance() output for one contrast.
#'
#' @return Data.frame per Nhood: mean share per condition, CV per condition,
#'   max_cv, top_sample_share (largest sample's part of the summed shares)
#'   and top_sample.
nhood_variability <- function(abundance) {

  cv <- function(x) if (length(x) < 2 || mean(x) == 0) NA_real_ else stats::sd(x) / mean(x)
  levels_cond <- levels(abundance$condition)   # c(reference, test)

  per_cond <- abundance %>%
    dplyr::group_by(Nhood, condition) %>%
    dplyr::summarise(mean_share = mean(share), cv = cv(share)) %>%
    dplyr::ungroup()
  wide <- per_cond %>%
    tidyr::pivot_wider(names_from = condition, values_from = c(mean_share, cv))

  top <- abundance %>%
    dplyr::group_by(Nhood) %>%
    dplyr::summarise(top_sample_share = if (sum(share) > 0) max(share) / sum(share) else NA_real_,
                     top_sample = milo_sample[which.max(share)]) %>%
    dplyr::ungroup()

  out <- dplyr::left_join(wide, top, by = "Nhood")
  cv_cols <- paste0("cv_", levels_cond)
  out$max_cv <- apply(out[, intersect(cv_cols, names(out)), drop = FALSE], 1, function(x) suppressWarnings(max(x, na.rm = TRUE)))
  out$max_cv[!is.finite(out$max_cv)] <- NA
  as.data.frame(out)
}
#' Why a neighborhood is (not) significant: logFC vs variability
#'
#' Two panels, one point per neighborhood: logFC against the share of its
#' abundance from the single largest sample, and against the larger
#' within-condition CV. Neighborhoods with a large |logFC| on the right of
#' either panel are driven by one sample or vary strongly within a
#' condition - the usual reason they are not significant.
#'
#' @param da_res DA results joined with nhood_variability().
#' @param alpha SpatialFDR threshold.
#' @param test_label,ref_label Contrast levels.
#' @param n_samples Samples in the contrast (sets the 1 / n reference line:
#'   the top-sample share if all samples contributed equally).
#' @param title,subtitle Plot title and subtitle.
#'
#' @return A patchwork plot.
plot_nhood_variability <- function(da_res,
                                   alpha,
                                   test_label,
                                   ref_label,
                                   n_samples,
                                   title = NULL,
                                   subtitle = NULL) {

  validate_required_cols(da_res, c("logFC", "SpatialFDR", "top_sample_share", "max_cv"), "plot_nhood_variability")
  pal <- resolve_palette("de_log2fc")
  df <- da_res
  df$DA <- factor(dplyr::case_when(
    df$SpatialFDR < alpha & df$logFC > 0 ~ paste("Enriched in", test_label),
    df$SpatialFDR < alpha & df$logFC < 0 ~ paste("Enriched in", ref_label),
    TRUE ~ "Not significant"
  ), levels = c(paste("Enriched in", ref_label), "Not significant", paste("Enriched in", test_label)))
  colors <- setNames(c(pal$low, "#B3B3B3", pal$high), levels(df$DA))

  base <- function(x, x_label) {
    ggplot2::ggplot(df, ggplot2::aes(x = .data[[x]], y = logFC, color = DA)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dotted", color = "#999999") +
      ggplot2::geom_point(size = 1.8, alpha = 0.8) +
      ggplot2::scale_color_manual(values = colors, drop = FALSE, name = NULL) +
      ggplot2::labs(x = x_label, y = paste0("DA logFC (", test_label, " vs ", ref_label, ")")) +
      ggplot2::theme_classic(base_size = 11)
  }
  p1 <- base("top_sample_share", "Share of abundance from the largest sample") +
    ggplot2::geom_vline(xintercept = 1 / n_samples, linetype = "dashed", color = "#808080")
  p2 <- base("max_cv", "Within-condition CV (larger of the two)")

  (p1 + p2) + patchwork::plot_layout(guides = "collect") +
    patchwork::plot_annotation(title = title, subtitle = subtitle) &
    ggplot2::theme(legend.position = "bottom")
}
#' Share of each sample's cells inside each cell set
#'
#' @param cell_sets Named list of cell barcodes (e.g. the cells of each
#'   meta-neighborhood).
#' @param meta Cell metadata with `milo_sample`, patient and timepoint.
#' @param patient_col,timepoint_col Study-design columns.
#'
#' @return Long data.frame: set, milo_sample, patient, timepoint, n_cells,
#'   sample_total, share.
set_sample_abundance <- function(cell_sets,
                                 meta,
                                 patient_col,
                                 timepoint_col) {

  samples <- unique(meta[, c("milo_sample", patient_col, timepoint_col)])
  totals <- table(meta$milo_sample)
  dplyr::bind_rows(lapply(names(cell_sets), function(set_name) {
    in_set <- table(factor(meta[intersect(cell_sets[[set_name]], rownames(meta)), "milo_sample"],
                           levels = samples$milo_sample))
    data.frame(set = set_name, milo_sample = samples$milo_sample,
               patient = samples[[patient_col]], timepoint = as.character(samples[[timepoint_col]]),
               n_cells = as.integer(in_set[samples$milo_sample]),
               sample_total = as.integer(totals[samples$milo_sample]))
  })) %>%
    dplyr::mutate(share = n_cells / sample_total, set = factor(set, levels = names(cell_sets)))
}
#' Per-sample shares of cell sets, by timepoint
#'
#' Boxplot + one point per sample; points of the same patient at repeated
#' timepoints are joined by a line, so paired changes are visible. This is
#' the quantity behind a DA result, shown sample by sample.
#'
#' @param abundance set_sample_abundance() output.
#' @param timepoint_levels Order of timepoints on the x axis.
#' @param repeated_timepoints Timepoints whose samples are joined per patient.
#' @param timepoint_colors Named colors.
#' @param title,subtitle Plot title and subtitle.
#'
#' @return A ggplot, one facet per set.
plot_sample_abundance <- function(abundance,
                                  timepoint_levels,
                                  repeated_timepoints,
                                  timepoint_colors,
                                  title = NULL,
                                  subtitle = NULL) {

  df <- abundance[abundance$timepoint %in% timepoint_levels, ]
  df$timepoint <- factor(df$timepoint, levels = timepoint_levels)
  paired_df <- df[df$timepoint %in% repeated_timepoints, ]

  ggplot2::ggplot(df, ggplot2::aes(x = timepoint, y = 100 * share)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = timepoint), outlier.shape = NA, alpha = 0.35, width = 0.6) +
    ggplot2::geom_line(data = paired_df, ggplot2::aes(group = patient), color = "#8C8C8C", linewidth = 0.4) +
    ggplot2::geom_point(ggplot2::aes(color = timepoint), size = 2) +
    ggplot2::scale_fill_manual(values = timepoint_colors, guide = "none") +
    ggplot2::scale_color_manual(values = timepoint_colors, guide = "none") +
    ggplot2::facet_wrap(~ set, scales = "free_y") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = "% of the sample's cells") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(strip.text = ggplot2::element_text(size = 8))
}
#' Per-cell DA call
#'
#' A cell is "enriched" if it lies in at least one significant neighborhood
#' with logFC > 0 and none with logFC < 0, "depleted" for the reverse,
#' "both" if in significant neighborhoods of both signs, "not DA"
#' otherwise. Cell-level calls make results comparable between analyses
#' whose neighborhoods differ (different cells build the graph).
#'
#' @param membership Cells x neighborhoods (miloR::nhoods()).
#' @param da_res DA results (Nhood, logFC, SpatialFDR).
#' @param alpha SpatialFDR threshold.
#'
#' @return Named character vector (cell -> call).
cell_da_calls <- function(membership,
                          da_res,
                          alpha) {

  up   <- da_res$Nhood[which(da_res$SpatialFDR < alpha & da_res$logFC > 0)]
  down <- da_res$Nhood[which(da_res$SpatialFDR < alpha & da_res$logFC < 0)]
  in_up   <- if (length(up))   Matrix::rowSums(membership[, up, drop = FALSE]) > 0 else rep(FALSE, nrow(membership))
  in_down <- if (length(down)) Matrix::rowSums(membership[, down, drop = FALSE]) > 0 else rep(FALSE, nrow(membership))
  setNames(ifelse(in_up & in_down, "both", ifelse(in_up, "enriched", ifelse(in_down, "depleted", "not DA"))),
           rownames(membership))
}
#' Spearman correlation between neighborhood-level variables
#'
#' E.g. mean QC metrics (rows) against mean program scores and DA logFC
#' (columns), across all neighborhoods. Descriptive: neighborhoods overlap,
#' so they are not independent observations.
#'
#' @param nhood_df One row per neighborhood.
#' @param row_vars,col_vars Columns to correlate.
#' @param title,subtitle Plot title and subtitle.
#'
#' @return A ggplot heatmap with the coefficients printed.
plot_nhood_correlations <- function(nhood_df,
                                    row_vars,
                                    col_vars,
                                    title = NULL,
                                    subtitle = NULL) {

  validate_required_cols(nhood_df, c(row_vars, col_vars), "plot_nhood_correlations")
  rho <- stats::cor(nhood_df[, row_vars, drop = FALSE], nhood_df[, col_vars, drop = FALSE],
                    method = "spearman", use = "pairwise.complete.obs")
  df <- data.frame(row = rep(row_vars, times = length(col_vars)),
                   col = rep(col_vars, each = length(row_vars)),
                   rho = as.vector(rho))
  df$row <- factor(df$row, levels = rev(row_vars))
  df$col <- factor(sub("Score$", "", df$col), levels = sub("Score$", "", col_vars))
  pal <- resolve_palette("de_log2fc")

  ggplot2::ggplot(df, ggplot2::aes(x = col, y = row, fill = rho)) +
    ggplot2::geom_tile(color = "#FFFFFF") +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", rho)), size = 3) +
    ggplot2::scale_fill_gradient2(low = pal$low, mid = pal$mid, high = pal$high, midpoint = 0,
                                  limits = c(-1, 1), name = "Spearman rho") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), panel.grid = ggplot2::element_blank())
}
#' Neighborhood groups from all neighborhoods
#'
#' Groups are communities of the neighborhood-overlap graph: neighborhoods
#' are nodes, edges are weighted by the Jaccard index of their cells
#' (shared / union), and Leiden (modularity) partitions the graph. DA
#' results are not used, so the groups can be tested and described without
#' circularity.
#'
#' @param membership Cells x neighborhoods (miloR::nhoods()).
#' @param resolution Leiden resolution (1 = standard modularity).
#' @param seed Random seed.
#'
#' @return Named character vector (neighborhood index -> group "NG1", ...),
#'   groups numbered by size (number of neighborhoods), largest first.
group_nhoods_by_overlap <- function(membership,
                                    resolution = 1,
                                    seed = 123) {

  m <- methods::as(membership > 0, "dMatrix")
  shared <- Matrix::crossprod(m)
  sizes <- Matrix::diag(shared)
  idx <- Matrix::which(Matrix::triu(shared, k = 1) > 0, arr.ind = TRUE)
  jaccard <- shared[idx] / (sizes[idx[, 1]] + sizes[idx[, 2]] - shared[idx])

  g <- igraph::graph_from_data_frame(
    data.frame(from = idx[, 1], to = idx[, 2], weight = jaccard),
    directed = FALSE, vertices = data.frame(name = seq_len(ncol(m)))
  )
  set.seed(seed)
  communities <- igraph::cluster_leiden(g, objective_function = "modularity", weights = igraph::E(g)$weight,
                                        resolution = resolution, n_iterations = -1)
  raw <- igraph::membership(communities)
  ranking <- names(sort(table(raw), decreasing = TRUE))
  setNames(paste0("NG", match(as.character(raw), ranking)), as.character(seq_len(ncol(m))))
}
#' One group per cell
#'
#' A cell can belong to several neighborhoods, possibly of different
#' groups. It is assigned to the group that most of its neighborhoods
#' belong to (ties: the group listed first, i.e. the larger one), so groups
#' become non-overlapping sets of cells and their counts never double-count
#' a cell. Cells in no neighborhood stay unassigned (NA).
#'
#' @param membership Cells x neighborhoods.
#' @param nhood_group group_nhoods_by_overlap() output.
#'
#' @return Named character vector (cell -> group or NA).
assign_cells_to_groups <- function(membership,
                                   nhood_group) {

  groups <- unique(nhood_group[order(as.numeric(sub("NG", "", nhood_group)))])
  indicator <- Matrix::sparseMatrix(i = seq_along(nhood_group), j = match(nhood_group, groups), x = 1,
                                    dims = c(length(nhood_group), length(groups)))
  votes <- as.matrix((membership > 0) %*% indicator)
  best <- max.col(votes, ties.method = "first")
  assigned <- groups[best]
  assigned[rowSums(votes) == 0] <- NA
  setNames(assigned, rownames(membership))
}
#' How the grouping changes at nearby resolutions
#'
#' For each resolution: the number of groups, and the adjusted Rand index
#' of the cells' group assignment against the chosen resolution (1 =
#' identical partition). A grouping that barely changes nearby is not an
#' artifact of the exact resolution chosen.
#'
#' @param membership Cells x neighborhoods.
#' @param resolutions Resolutions to compare.
#' @param chosen The resolution used.
#' @param seed Random seed.
#'
#' @return Data.frame: resolution, n_groups, ari_vs_chosen.
group_resolution_stability <- function(membership,
                                       resolutions,
                                       chosen,
                                       seed = 123) {

  reference <- assign_cells_to_groups(membership, group_nhoods_by_overlap(membership, chosen, seed))
  dplyr::bind_rows(lapply(sort(unique(c(resolutions, chosen))), function(r) {
    groups <- group_nhoods_by_overlap(membership, r, seed)
    data.frame(resolution = r, n_groups = length(unique(groups)),
               ari_vs_chosen = round(partition_ari(reference, assign_cells_to_groups(membership, groups)), 3),
               chosen = r == chosen)
  }))
}
#' Bootstrap stability of the neighborhood grouping across resolutions
#'
#' For each resolution: group all neighborhoods, then repeatedly regroup a
#' random subset of them (with a new seed each time) and compare the two
#' groupings on the shared neighborhoods (adjusted Rand index). A
#' resolution whose groups come back when a fifth of the neighborhoods is
#' left out describes real structure in the overlap graph. Uses no DA
#' result, so choosing on it does not bias the groups.
#'
#' @param membership Cells x neighborhoods.
#' @param resolutions Resolutions to screen.
#' @param n_boot Subsamples per resolution. Default 30.
#' @param frac Share of neighborhoods kept per subsample. Default 0.8.
#' @param seed Random seed.
#'
#' @return Data.frame: resolution, n_groups (all neighborhoods),
#'   mean_ari, sd_ari, smallest_group (neighborhoods).
group_resolution_sweep <- function(membership,
                                   resolutions,
                                   n_boot = 30,
                                   frac = 0.8,
                                   seed = 123) {

  n <- ncol(membership)
  dplyr::bind_rows(lapply(resolutions, function(r) {
    full <- group_nhoods_by_overlap(membership, r, seed)
    set.seed(seed)
    scores <- vapply(seq_len(n_boot), function(b) {
      keep <- sort(sample(n, round(frac * n)))
      sub <- group_nhoods_by_overlap(membership[, keep, drop = FALSE], r, seed + b)
      partition_ari(full[as.character(keep)], sub)
    }, numeric(1))
    data.frame(resolution = r, n_groups = length(unique(full)),
               mean_ari = round(mean(scores), 3), sd_ari = round(stats::sd(scores), 3),
               smallest_group = min(table(full)))
  }))
}
#' Adjusted Rand index of two partitions
#'
#' Pairs where either value is NA are left out. Two partitions that are
#' both a single group (no pair can disagree) give 1.
#'
#' @param a,b Group labels of the same items.
#'
#' @return The adjusted Rand index.
#'
#' @examples
#' partition_ari(c(1, 1, 2, 2), c("x", "x", "y", "y"))   # 1
partition_ari <- function(a,
                          b) {
  keep <- !is.na(a) & !is.na(b)
  tab <- table(a[keep], b[keep])
  comb2 <- function(x) x * (x - 1) / 2
  index    <- sum(comb2(tab))
  rows     <- sum(comb2(rowSums(tab)))
  cols     <- sum(comb2(colSums(tab)))
  expected <- rows * cols / comb2(sum(tab))
  max_index <- (rows + cols) / 2
  if (max_index == expected) return(1)
  (index - expected) / (max_index - expected)
}
#' Per-sample abundance test of non-overlapping cell sets
#'
#' Counts each set's cells per sample and tests them with the contrast's
#' design (edgeR quasi-likelihood negative binomial GLM), using each
#' sample's total cells as library size. With few sets the TMM assumption
#' (most sets unchanged) does not hold, so total cells are used instead.
#' With so few sets, no dispersion or abundance trend is fitted.”
#' Used for neighborhood groups and, as a comparator, for Seurat clusters.
#'
#' @param cell_set Named vector (cell -> set; NA = not in any set).
#' @param meta Cell metadata with `milo_sample` (rownames = cells).
#' @param design contrast_design() output.
#'
#' @return Data.frame per set: logFC (test vs reference), PValue, FDR (BH
#'   across sets), mean share per condition.
test_set_abundance <- function(cell_set,
                               meta,
                               design) {

  design_df <- design$design_df
  samples <- design_df$milo_sample
  in_samples <- rownames(meta)[meta$milo_sample %in% samples]
  totals <- table(factor(meta[in_samples, "milo_sample"], levels = samples))
  sets <- sort(unique(stats::na.omit(cell_set[in_samples])))
  counts <- as.matrix(table(factor(cell_set[in_samples], levels = sets),
                            factor(meta[in_samples, "milo_sample"], levels = samples)))

  X <- stats::model.matrix(design$formula, data = design_df)
  y <- edgeR::estimateDisp(y, X, trend.method = "none")
  fit <- edgeR::glmQLFit(y, X, robust = TRUE, abundance.trend = FALSE)
  fit <- edgeR::glmQLFit(y, X, robust = TRUE)
  res <- edgeR::topTags(edgeR::glmQLFTest(fit, coef = ncol(X)), n = Inf, sort.by = "none")$table

  share <- sweep(counts, 2, as.numeric(totals), "/")
  levels_cond <- levels(design_df$condition)
  out <- data.frame(set = rownames(res), n_cells = rowSums(counts), logFC = res$logFC,
                    PValue = res$PValue, FDR = res$FDR)
  for (lv in rev(levels_cond)) {
    out[[paste0("mean_pct_", lv)]] <- round(100 * rowMeans(share[, design_df$condition == lv, drop = FALSE]), 2)
  }
  out[order(out$PValue), ]
}
#' Cluster and timepoint composition of cell sets
#'
#' Raw composition = share of the set's cells in each category. For
#' timepoints a size-normalized version is added: each timepoint's cells
#' are divided by that timepoint's total cells before taking shares, so a
#' timepoint contributing many cells overall does not dominate every set
#' (1 / n_timepoints = no enrichment).
#'
#' @param cells_by_set Named list of cell barcodes per set.
#' @param meta Cell metadata (rownames = cells).
#' @param column Metadata column (cluster label or timepoint).
#' @param normalize Add the size-normalized shares. Default FALSE.
#'
#' @return Long data.frame: set, category, fraction, type ("raw" / "normalized").
set_composition <- function(cells_by_set,
                            meta,
                            column,
                            normalize = FALSE) {

  categories <- sort(unique(as.character(meta[[column]])))
  totals <- table(factor(meta[[column]], levels = categories))
  dplyr::bind_rows(lapply(names(cells_by_set), function(set_name) {
    counts <- table(factor(meta[intersect(cells_by_set[[set_name]], rownames(meta)), column], levels = categories))
    out <- data.frame(set = set_name, category = categories, fraction = as.numeric(counts) / sum(counts), type = "raw")
    if (normalize) {
      scaled <- as.numeric(counts) / as.numeric(totals)
      out <- rbind(out, data.frame(set = set_name, category = categories,
                                   fraction = scaled / sum(scaled), type = "normalized"))
    }
    out
  }))
}
#' Stacked composition bars for many cell sets
#'
#' @param composition set_composition() output (one type).
#' @param set_order Order of the sets on the x axis (e.g. the heatmap order).
#' @param colors Named colors per category.
#' @param title,subtitle Plot title and subtitle.
#'
#' @return A ggplot.
plot_set_composition <- function(composition,
                                 set_order,
                                 colors,
                                 title,
                                 subtitle = NULL) {

  composition$set <- factor(composition$set, levels = set_order)
  ggplot2::ggplot(composition, ggplot2::aes(x = set, y = 100 * fraction, fill = category)) +
    ggplot2::geom_col(width = 0.9) +
    ggplot2::scale_fill_manual(values = colors, name = NULL) +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = "% of cells") +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1,
                                                       size = if (length(set_order) > 40) 5 else 8))
}
#' Program heatmap (pheatmap)
#'
#' One question only: which programs characterize each row. Rows are
#' neighborhoods or groups; values are mean net program scores, z-scored
#' per program. If `row_groups` is given, rows are kept together by group
#' (gaps between groups) and ordered by hierarchical clustering (Ward,
#' Euclidean) within each group; otherwise all rows are clustered together.
#' Programs (columns) are clustered too. Written straight to a PNG.
#'
#' @param score_matrix Rows x programs (mean scores).
#' @param annotation Optional data.frame of row annotations (rownames = rows).
#' @param annotation_colors Named list of colors per annotation column.
#' @param title Heatmap title.
#' @param path Output PNG path.
#' @param row_groups Optional named vector (row -> group) for grouped order.
#' @param show_rownames Print row names. Default TRUE.
#'
#' @return Invisibly, the row order (row names, top to bottom).
plot_annotated_heatmap <- function(score_matrix,
                                   annotation = NULL,
                                   annotation_colors = NULL,
                                   title,
                                   path,
                                   row_groups = NULL,
                                   show_rownames = TRUE) {

  z <- scale(score_matrix)
  z[is.na(z)] <- 0
  colnames(z) <- sub("_NETScore$|Score$", "", colnames(z))
  pal <- resolve_palette("de_log2fc")

  if (is.null(row_groups)) {
    row_order <- rownames(z)[stats::hclust(stats::dist(z), method = "ward.D2")$order]
    gaps <- NULL
  } else {
    groups <- unique(row_groups[rownames(z)])
    groups <- groups[order(as.numeric(gsub("\\D", "", groups)))]
    row_order <- unlist(lapply(groups, function(g) {
      rows <- rownames(z)[row_groups[rownames(z)] == g]
      if (length(rows) > 2) rows[stats::hclust(stats::dist(z[rows, , drop = FALSE]), method = "ward.D2")$order] else rows
    }))
    gaps <- utils::head(cumsum(table(factor(row_groups[row_order], levels = groups))), -1)
  }
  limit <- max(abs(z))

  pheatmap::pheatmap(
    z[row_order, , drop = FALSE], cluster_rows = FALSE, cluster_cols = ncol(z) > 2, gaps_row = gaps,
    color = grDevices::colorRampPalette(c(pal$low, pal$mid, pal$high))(100),
    breaks = seq(-limit, limit, length.out = 101),
    annotation_row = if (!is.null(annotation)) annotation[row_order, , drop = FALSE],
    annotation_colors = annotation_colors,
    show_rownames = show_rownames, fontsize_row = if (nrow(z) > 40) 5 else 10,
    angle_col = 45, main = title, filename = path,
    width = 5 + 0.7 * ncol(z) + if (!is.null(annotation)) 1.5 else 0,
    height = max(4, min(22, 2.5 + 0.12 * nrow(z)))
  )
  invisible(row_order)
}
#' DA of non-overlapping cell sets, as a dot plot (sets x contrasts)
#'
#' Color = logFC (red: more abundant in the first timepoint of the
#' contrast); size = -log10 FDR; black outline = FDR below `alpha`. One
#' panel per set type (e.g. neighborhood groups, Seurat clusters), so the
#' two levels can be compared side by side.
#'
#' @param tables Named list (panel) of named lists (contrast) of
#'   test_set_abundance() outputs.
#' @param alpha FDR threshold for the outline.
#' @param title,subtitle Plot title and subtitle.
#'
#' @return A ggplot.
plot_set_da_dots <- function(tables,
                             alpha,
                             title,
                             subtitle = NULL) {

  df <- dplyr::bind_rows(lapply(names(tables), function(panel) {
    dplyr::bind_rows(lapply(names(tables[[panel]]), function(ct) {
      t <- tables[[panel]][[ct]]
      data.frame(panel = panel, contrast = ct, set = t$set, logFC = t$logFC, FDR = t$FDR)
    }))
  }))
  df$contrast <- factor(df$contrast, levels = unique(df$contrast))
  df$sig <- !is.na(df$FDR) & df$FDR < alpha
  df$label <- sprintf("%.2f", df$FDR)
  limit <- max(abs(df$logFC), na.rm = TRUE)
  pal <- resolve_palette("de_log2fc")

  ggplot2::ggplot(df, ggplot2::aes(x = contrast, y = set)) +
    ggplot2::geom_point(ggplot2::aes(fill = logFC, size = -log10(FDR)), shape = 21, color = "#666666", stroke = 0.3) +
    ggplot2::geom_point(data = df[df$sig, , drop = FALSE], ggplot2::aes(fill = logFC, size = -log10(FDR)),
                        shape = 21, color = "#000000", stroke = 1.4) +
    ggplot2::geom_text(ggplot2::aes(label = label), nudge_x = 0.32, size = 3, color = "#4D4D4D") +
    ggplot2::scale_fill_gradient2(low = pal$low, mid = pal$mid, high = pal$high, limits = c(-limit, limit)) +
    ggplot2::scale_size_continuous(range = c(2, 9), name = "-log10 FDR") +
    ggplot2::facet_wrap(~ panel, scales = "free_y") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 30, hjust = 1))
}
