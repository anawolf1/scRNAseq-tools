# ------------------------------------------------------------------------------
# cluster_characterization.R
#
# Characterization of a clustered Seurat object with a UMAP: composition
# and stability, QC metrics, UMAPs by group, cell counts, marker columns,
# and differential expression between clusters. Each section of
# characterize_clustering() is switched on or off by an argument, which
# the pipeline sets from the `steps` toggles; no function reads the config.
#
# Functions
# check_marker_significance()   DE results of given genes in every cluster
# plot_cluster_volcanoes()      volcano plot per cluster (EnhancedVolcano)
# resolve_cluster_display()     cluster display labels and colors
# characterize_clustering()     every characterization section of a round
# pairwise_cluster_de()         DE between every pair of clusters
# compare_deg_venn()            Venn diagram and workbook of DEGs shared between 2-4 conditions
# plot_deg_count_comparison()   bar chart of DEG counts per condition
# write_de_workbook()           one-vs-rest DE results as a workbook, one sheet per cluster
#
# Enrichment (GO/KEGG/curated gene sets) is in the shared enrichment.R.
# ------------------------------------------------------------------------------

#' DE results of given genes in every cluster
#'
#' For genes expected to mark a population (e.g. TRDC for gamma-delta T
#' cells): their one-vs-rest results in each cluster. Only rows present in
#' `de_df` are returned; a gene not tested in a cluster (below
#' FindAllMarkers()'s thresholds) has no row there.
#'
#' @param de_df FindAllMarkers() table (cluster, gene, p_val, avg_log2FC,
#'   pct.1, pct.2, p_val_adj).
#' @param genes Genes to report.
#'
#' @return A data.frame, one row per gene and cluster, values rounded for
#'   display.
#'
#' @examples
#' \dontrun{
#' check_marker_significance(de_results, c("TRDC", "TRGC2"))
#' }
check_marker_significance <- function(de_df,
                                      genes) {
  required_cols <- c("cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
  validate_required_cols(de_df, required_cols, "check_marker_significance")
  
  result <- de_df[de_df$gene %in% genes, required_cols]
  rownames(result) <- NULL
  result$avg_log2FC <- round(result$avg_log2FC, 2)
  result$pct.1      <- round(result$pct.1, 2)
  result$pct.2      <- round(result$pct.2, 2)
  result$p_val      <- formatC(result$p_val, format = "e", digits = 2)
  result$p_val_adj  <- formatC(result$p_val_adj, format = "e", digits = 2)
  result
}
#' Volcano plot per cluster
#'
#' One EnhancedVolcano plot per cluster of a FindAllMarkers() table, or one
#' plot of a FindMarkers() table (no `cluster` column; returned as
#' "comparison"). The `n_labels` most significant up and down genes are
#' labeled. Genes to leave out (e.g. "^ENSG") are removed from `de_results`
#' beforehand. The plots are also printed together.
#'
#' @param de_results DE table with `gene`, `avg_log2FC` and `p_val_adj`
#'   (and `cluster` for FindAllMarkers()).
#' @param subtitle Subtitle of every plot; the numbers of up and down genes
#'   are appended.
#' @param log2fc_cutoff |log2FC| above which a gene is up or down. Default 0.5.
#' @param padj_cutoff Adjusted p-value cutoff. Default 0.05.
#' @param n_labels Up and down genes labeled per plot. Default 50.
#' @param cluster_display_labels Optional named vector, cluster ID -> display
#'   label (e.g. c("1" = "1: Naive")), for the plot titles.
#'
#' @return Invisibly, a named list of plots (names: cluster IDs, or
#'   "comparison"). A plot that fails is left out, with a message.
#'
#' @examples
#' \dontrun{
#' volcanoes <- plot_cluster_volcanoes(de_results, subtitle = "vs all other cells")
#' save_plot(volcanoes[["1"]], "volcano_cluster1.png", fig_dir)
#' }
plot_cluster_volcanoes <- function(de_results,
                                   subtitle,
                                   log2fc_cutoff = 0.5,
                                   padj_cutoff = 0.05,
                                   n_labels = 50,
                                   cluster_display_labels = NULL) {
  
  validate_required_cols(de_results, c("gene", "avg_log2FC", "p_val_adj"), "plot_cluster_volcanoes")
  if (!"cluster" %in% colnames(de_results)) de_results$cluster <- "comparison"
  
  volcano_list <- list()
  for (cluster_id in unique(as.character(de_results$cluster))) {
    cluster_de <- de_results[as.character(de_results$cluster) == cluster_id, , drop = FALSE]
    
    significant <- !is.na(cluster_de$p_val_adj) & cluster_de$p_val_adj < padj_cutoff
    is_up   <- significant & cluster_de$avg_log2FC > log2fc_cutoff
    is_down <- significant & cluster_de$avg_log2FC < -log2fc_cutoff
    top_genes <- function(rows) {
      rows <- cluster_de[rows, , drop = FALSE]
      head(rows$gene[order(rows$p_val_adj, -abs(rows$avg_log2FC))], n_labels)
    }
    label_genes <- unique(c(top_genes(is_up), top_genes(is_down)))
    
    title <- if (cluster_id == "comparison") {
      "Differential expression"
    } else {
      paste("Cluster", if (cluster_id %in% names(cluster_display_labels)) cluster_display_labels[[cluster_id]] else cluster_id)
    }
    
    volcano_list[[cluster_id]] <- tryCatch(
      EnhancedVolcano::EnhancedVolcano(
        cluster_de, lab = cluster_de$gene, selectLab = if (length(label_genes) > 0) label_genes,
        x = "avg_log2FC", y = "p_val_adj", pCutoff = padj_cutoff, FCcutoff = log2fc_cutoff,
        pointSize = 2, labSize = 4, labCol = "black", labFace = "plain", boxedLabels = FALSE,
        drawConnectors = TRUE, widthConnectors = 0.5, colConnectors = "grey30", max.overlaps = Inf,
        title = title, subtitle = paste0(subtitle, " (", sum(is_up), " up, ", sum(is_down), " down)"),
        caption = "", legendPosition = "top", legendLabSize = 10, legendIconSize = 3,
        cutoffLineType = "dashed", cutoffLineCol = "black", cutoffLineWidth = 0.4,
        gridlines.major = FALSE, gridlines.minor = FALSE
      ),
      error = function(e) {
        message("plot_cluster_volcanoes(): no plot for ", cluster_id, ": ", conditionMessage(e))
        NULL
      }
    )
  }
  volcano_list <- Filter(Negate(is.null), volcano_list)
  
  if (length(volcano_list) > 0) {
    tryCatch(print(patchwork::wrap_plots(volcano_list, ncol = 3)), error = function(e) {
      message("plot_cluster_volcanoes(): the plot pane is too small for the combined plots; ",
              "the plots are still returned.")
    })
  }
  invisible(volcano_list)
}
#' Cluster display labels and colors
#'
#' Display labels are "<ID>: <label>" for labeled clusters and the ID
#' otherwise. Colors come from resolve_cluster_colors() (shared
#' palettes.R): a label in `cluster_color_map` keeps its color in every
#' round and workflow.
#'
#' @param cluster_levels Cluster IDs, e.g. levels(obj$seurat_clusters).
#' @param cluster_labels Optional named vector, cluster ID -> label.
#' @param cluster_color_map Optional named vector, label -> color.
#'
#' @return A list: `labels` (cluster ID -> display label) and `colors`
#'   (display label -> color).
#'
#' @examples
#' \dontrun{
#' resolve_cluster_display(c("1", "2", "3"), cluster_labels = c("2" = "Naive", "3" = "CD127-"),
#'                         cluster_color_map = c("Naive" = "#E69F00"))
#' }
resolve_cluster_display <- function(cluster_levels,
                                    cluster_labels = NULL,
                                    cluster_color_map = NULL) {
  
  labels <- if (!is.null(cluster_labels)) unname(cluster_labels[cluster_levels]) else rep(NA, length(cluster_levels))
  display_labels <- setNames(ifelse(is.na(labels), cluster_levels, paste0(cluster_levels, ": ", labels)),
                             cluster_levels)
  
  colors <- unname(resolve_cluster_colors(labels, cluster_color_map))
  list(labels = display_labels, colors = setNames(colors, unname(display_labels)))
}
#' Every characterization section of a round
#'
#' For a clustered object with a UMAP: composition and per-cluster
#' stability, QC metrics by cluster, UMAPs by group, cell counts, marker
#' columns and signature genes, and one-vs-rest DE. Each section is
#' switched on or off with its run_* argument. `title_label` is used in
#' titles (e.g. the round's label, "Lineage screening"), `file_prefix` in
#' file names (e.g. "round_1").
#'
#' @param obj Seurat object, clustered (`seurat_clusters`), with a "umap"
#'   reduction.
#' @param qc_metrics QC metric columns.
#' @param qc_labels Optional display labels for `qc_metrics`.
#' @param group_by_vars Metadata columns for the group UMAPs and the cells
#'   per cluster tables.
#' @param cell_count_vars Cell count breakdowns (see plot_cell_counts()).
#' @param display_order Level order per column (e.g. params$display_order).
#' @param cluster_labels Optional named vector, cluster ID -> label (the
#'   round's cluster_labels); unlabeled clusters show their ID.
#' @param cluster_color_map Optional named vector, label -> color (see
#'   resolve_cluster_display()).
#' @param stability_dims,stability_resolution,stability_algorithm,stability_k_param,stability_prune_snn,stability_reduction
#'   The clustering parameters that made `obj`'s clusters, for
#'   per_cluster_stability(). k_param, prune_snn and reduction default to
#'   20, 1/15 and "pca".
#' @param stability_n_reps,stability_frac Passed to per_cluster_stability()
#'   as n_reps and frac. Defaults 10, 0.8.
#' @param marker_check_cols Optional named list, metadata column -> level:
#'   a UMAP colored by each column and a table of its levels per cluster.
#' @param violin_genes Optional genes for a violin plot by cluster.
#' @param highlight_genes Optional genes whose DE results are reported for
#'   every cluster.
#' @param de_min_pct,de_log2fc_cutoff FindAllMarkers() min.pct and
#'   logfc.threshold. Defaults 0.1, 0.5.
#' @param de_padj_cutoff Adjusted p-value cutoff in the volcanoes. Default 0.05.
#' @param de_n_labels Genes labeled per direction in each volcano. Default 50.
#' @param de_exclude_gene_patterns Regex patterns of genes removed from the
#'   DE results (and so from every DE table, plot and the returned table).
#'   Default NULL.
#' @param run_composition_stability,run_qc_metrics,run_group_umaps,run_cell_counts,run_marker_characterization,run_differential_expression
#'   Logical, one per section. All default TRUE.
#' @param seed Seed for the stability subsampling and clustering. Default 123.
#' @param output_dir Directory for every table and plot.
#' @param file_prefix Prefix of every file name.
#' @param title_label Text starting every title.
#'
#' @return Invisibly, a list: `composition`, `qc_medians`, `de_results` and
#'   `highlight_check` (NULL for a section that did not run).
#'
#' @examples
#' \dontrun{
#' characterize_clustering(
#'   obj, qc_metrics = qc_info$qc_metrics, group_by_vars = params$group_by_vars,
#'   cell_count_vars = params$cell_count_vars, display_order = params$display_order,
#'   stability_dims = 10, stability_resolution = 0.4, stability_algorithm = 1,
#'   output_dir = round_dir, file_prefix = "round_1", title_label = "Lineage screening"
#' )
#' }
characterize_clustering <- function(obj,
                                    qc_metrics,
                                    qc_labels = NULL,
                                    group_by_vars = character(0),
                                    cell_count_vars = list(),
                                    display_order = list(),
                                    cluster_labels = NULL,
                                    cluster_color_map = NULL,
                                    stability_dims,
                                    stability_resolution,
                                    stability_algorithm,
                                    stability_k_param = 20,
                                    stability_prune_snn = 1 / 15,
                                    stability_reduction = "pca",
                                    stability_n_reps = 10,
                                    stability_frac = 0.8,
                                    marker_check_cols = NULL,
                                    violin_genes = NULL,
                                    highlight_genes = NULL,
                                    de_min_pct = 0.1,
                                    de_log2fc_cutoff = 0.5,
                                    de_padj_cutoff = 0.05,
                                    de_n_labels = 50,
                                    de_exclude_gene_patterns = NULL,
                                    run_composition_stability = TRUE,
                                    run_qc_metrics = TRUE,
                                    run_group_umaps = TRUE,
                                    run_cell_counts = TRUE,
                                    run_marker_characterization = TRUE,
                                    run_differential_expression = TRUE,
                                    seed = 123,
                                    output_dir,
                                    file_prefix,
                                    title_label) {
  
  cluster_levels <- levels(obj$seurat_clusters)
  cluster_display <- resolve_cluster_display(cluster_levels, cluster_labels, cluster_color_map)
  # the same colors keyed by cluster ID, for plots grouped by seurat_clusters
  cluster_colors_by_id <- setNames(unname(cluster_display$colors), cluster_levels)
  obj$Cluster <- factor(unname(cluster_display$labels[as.character(obj$seurat_clusters)]),
                        levels = unname(cluster_display$labels))
  
  composition <- NULL
  qc_medians <- NULL
  de_results <- NULL
  highlight_check <- NULL
  
  # --- composition and stability ------------------------------------------------
  if (run_composition_stability) {
    stability <- per_cluster_stability(
      obj, dims = stability_dims, resolution = stability_resolution, algorithm = stability_algorithm,
      k_param = stability_k_param, prune_snn = stability_prune_snn, reduction = stability_reduction,
      n_reps = stability_n_reps, frac = stability_frac, seed = seed
    )
    n_cells <- as.integer(table(obj$seurat_clusters)[cluster_levels])
    composition <- data.frame(
      cluster      = unname(cluster_display$labels),
      n_cells      = n_cells,
      pct_cells    = round(100 * n_cells / sum(n_cells), 1),
      mean_jaccard = round(stability$mean_jaccard[match(cluster_levels, stability$cluster)], 2),
      sd_jaccard   = round(stability$sd_jaccard[match(cluster_levels, stability$cluster)], 2)
    )
    composition_gt <- gt_generic(composition, title = paste(title_label, "— cluster composition and stability"))
    print(composition_gt)
    save_gt(composition_gt, paste0(file_prefix, "_composition_stability.png"), output_dir)
  }
  
  # --- QC metrics by cluster ------------------------------------------------------
  if (run_qc_metrics) {
    validate_required_cols(obj[[]], qc_metrics, "characterize_clustering",
                           hint = "expects the QC workflow's metric columns (check that params$percent_features matches it).")
    qc_medians <- obj[[]] %>%
      dplyr::group_by(cluster = Cluster) %>%
      dplyr::summarise(dplyr::across(dplyr::all_of(qc_metrics), \(x) stats::median(x, na.rm = TRUE)),
                       .groups = "drop")
    qc_medians_gt <- gt_generic(format_qc_values(qc_medians, qc_metrics),
                                title = paste(title_label, "— median QC metrics by cluster"))
    print(qc_medians_gt)
    save_gt(qc_medians_gt, paste0(file_prefix, "_qc_medians.png"), output_dir)
    
    qc_violin <- plot_violin(obj, features = qc_metrics, group_col = "Cluster",
                             group_colors = cluster_display$colors, group_order = levels(obj$Cluster),
                             feature_labels = qc_labels, show_points = FALSE,
                             title = paste(title_label, "— QC metrics by cluster"))
    print(qc_violin)
    save_plot(qc_violin, paste0(file_prefix, "_qc_violin.png"), output_dir, width = 10, height = 6)
    
    qc_feature_plot <- Seurat::FeaturePlot(obj, features = qc_metrics, reduction = "umap", ncol = 2,
                                           pt.size = 1, order = TRUE) +
      patchwork::plot_annotation(title = paste(title_label, "— QC metrics"))
    print(qc_feature_plot)
    save_plot(qc_feature_plot, paste0(file_prefix, "_qc_featureplot.png"), output_dir, width = 10, height = 10)
  }
  
  # --- UMAP by cluster and by each grouping column ----------------------------------
  if (run_group_umaps) {
    umap_cluster <- plot_umap_by_group(obj, "Cluster", group_order = levels(obj$Cluster),
                                       group_colors = cluster_display$colors, subtitle = title_label)
    print(umap_cluster$combined)
    save_plot(umap_cluster$combined, paste0(file_prefix, "_umap_cluster.png"), output_dir)
    
    for (group_col in group_by_vars) {
      umap_group <- plot_umap_by_group(obj, group_col, subtitle = title_label)
      print(umap_group$combined)
      save_plot(umap_group$combined, paste0(file_prefix, "_umap_", group_col, ".png"), output_dir)
      print(umap_group$faceted)
      save_plot(umap_group$faceted, paste0(file_prefix, "_umap_", group_col, "_faceted.png"),
                output_dir, width = 10, height = 8)
    }
  }
  
  # --- cell counts: params$cell_count_vars, and cells per cluster and group -----------
  if (run_cell_counts) {
    plot_cell_counts(obj, cell_count_vars, display_order, title = paste(title_label, "— cell counts"),
                     output_dir = output_dir, file_prefix = paste0(file_prefix, "_cell_counts"))
    plot_cell_counts(obj, lapply(group_by_vars, function(group_col) c("Cluster", group_col)),
                     utils::modifyList(display_order, list(Cluster = levels(obj$Cluster))),
                     title = paste(title_label, "— cells per cluster"),
                     output_dir = output_dir, file_prefix = paste0(file_prefix, "_cluster_counts"))
  }
  
  # --- marker columns and signature genes ---------------------------------------------
  if (run_marker_characterization) {
    validate_required_cols(obj[[]], names(marker_check_cols), "characterize_clustering",
                           hint = "expects every marker_check_cols name to be a metadata column.")
    for (col_name in names(marker_check_cols)) {
      marker_dimplot <- Seurat::DimPlot(obj, reduction = "umap", pt.size = 1.3, group.by = col_name) +
        ggplot2::labs(title = paste(title_label, "—", col_name))
      print(marker_dimplot)
      save_plot(marker_dimplot, paste0(file_prefix, "_", col_name, "_dimplot.png"), output_dir)
      
      # NA is shown as its own column, so percentages are of all cells of the cluster
      counts <- table(obj$Cluster, obj[[col_name, drop = TRUE]], useNA = "ifany")
      pct <- round(100 * prop.table(counts, margin = 1), 1)
      counts_pct <- matrix(paste0(counts, " (", pct, "%)"), nrow = nrow(counts), dimnames = dimnames(counts))
      counts_gt <- gt_generic(counts_pct, title = paste(title_label, "—", col_name, "by cluster"),
                              row_label = "Cluster")
      print(counts_gt)
      save_gt(counts_gt, paste0(file_prefix, "_", col_name, "_by_cluster.png"), output_dir)
    }
    
    if (length(violin_genes) > 0) {
      missing_genes <- setdiff(violin_genes, rownames(obj))
      if (length(missing_genes) > 0) {
        message(file_prefix, ": violin_genes not in the object, left out: ", paste(missing_genes, collapse = ", "))
      }
      # VlnPlot() fails on a gene with the same value in every cell (e.g. a
      # marker of a lineage absent from the object), so those are left out
      present_genes <- intersect(violin_genes, rownames(obj))
      constant_genes <- character(0)
      if (length(present_genes) > 0) {
        expression <- Seurat::FetchData(obj, vars = present_genes)
        constant_genes <- present_genes[vapply(expression, function(x) length(unique(x)) <= 1, logical(1))]
      }
      if (length(constant_genes) > 0) {
        message(file_prefix, ": violin_genes with the same value in every cell, left out: ",
                paste(constant_genes, collapse = ", "))
      }
      present_genes <- setdiff(present_genes, constant_genes)
      
      if (length(present_genes) > 0) {
        gene_violin <- Seurat::VlnPlot(obj, group.by = "seurat_clusters", features = present_genes,
                                       cols = cluster_colors_by_id) +
          patchwork::plot_annotation(title = paste(title_label, "— signature genes by cluster"))
        print(gene_violin)
        save_plot(gene_violin, paste0(file_prefix, "_gene_violin.png"), output_dir, width = 16, height = 14)
      }
    }
  }
  
  # --- one-vs-rest differential expression ---------------------------------------------
  if (run_differential_expression) {
    Seurat::Idents(obj) <- obj$seurat_clusters
    de_results <- safe_find_all_markers(obj, test.use = "wilcox", min.pct = de_min_pct,
                                        logfc.threshold = de_log2fc_cutoff, verbose = FALSE)
    if (length(de_exclude_gene_patterns) > 0) {
      de_results <- de_results[!grepl(paste(de_exclude_gene_patterns, collapse = "|"), de_results$gene), ,
                               drop = FALSE]
    }
    
    volcano_plots <- plot_cluster_volcanoes(
      de_results, subtitle = paste(title_label, "vs all other cells"),
      log2fc_cutoff = de_log2fc_cutoff, padj_cutoff = de_padj_cutoff, n_labels = de_n_labels,
      cluster_display_labels = cluster_display$labels
    )
    for (cluster_id in names(volcano_plots)) {
      save_plot(volcano_plots[[cluster_id]], paste0(file_prefix, "_volcano_cluster", cluster_id, ".png"),
                output_dir, width = 10, height = 8)
    }
    
    for (cluster_id in intersect(cluster_levels, as.character(de_results$cluster))) {
      cluster_de <- de_results[as.character(de_results$cluster) == cluster_id, , drop = FALSE]
      cluster_de <- cluster_de[order(cluster_de$p_val_adj), , drop = FALSE]
      de_table <- gt_de(
        cluster_de, title = paste(title_label, "— cluster", cluster_display$labels[[cluster_id]]),
        subtitle = "Top differentially expressed genes vs all other cells",
        pct1_label = paste("Cluster", cluster_display$labels[[cluster_id]]), pct2_label = "Other cells",
        gene_from_rownames = FALSE
      )
      save_gt(de_table, paste0(file_prefix, "_de_table_cluster", cluster_id, ".png"), output_dir)
    }
    
    if (length(highlight_genes) > 0) {
      highlight_check <- check_marker_significance(de_results, highlight_genes)
      print(highlight_check)
      highlight_gt <- gt_generic(highlight_check, title = paste(title_label, "— highlighted genes"))
      save_gt(highlight_gt, paste0(file_prefix, "_highlight_genes.png"), output_dir)
    }
  }
  
  invisible(list(composition = composition, qc_medians = qc_medians,
                 de_results = de_results, highlight_check = highlight_check))
}
#' DE between every pair of clusters
#'
#' FindMarkers() for every pair of the clusters left after `exclude_labels`
#' (e.g. clusters A, B, C and D without D: A vs B, A vs C, B vs C), with a
#' table and a volcano plot per pair. The cells of excluded clusters are
#' removed first, so they are in no comparison. Clusters with the same label
#' are compared as one group.
#'
#' @param obj Seurat object, clustered (`seurat_clusters`).
#' @param cluster_labels Named vector, cluster ID -> label, for every cluster.
#' @param exclude_labels Labels of clusters to leave out. Default none.
#' @param min_pct,log2fc_cutoff FindMarkers() min.pct and logfc.threshold.
#'   Defaults 0.1, 0.5.
#' @param padj_cutoff Adjusted p-value cutoff in the volcanoes. Default 0.05.
#' @param n_labels Genes labeled per direction in each volcano. Default 50.
#' @param exclude_gene_patterns Regex patterns of genes removed from the
#'   results (e.g. "^ENSG"). Default NULL.
#' @param title_label Text starting every title.
#' @param output_dir Directory for the tables and plots.
#' @param file_prefix Prefix of every file name.
#'
#' @return Invisibly, a list named "<label 1> vs <label 2>", each with
#'   `de_results` (FindMarkers() table, genes as row names, sorted by
#'   adjusted p-value) and `volcano` (NULL when there are no results).
#'
#' @examples
#' \dontrun{
#' pairwise_cluster_de(obj, cluster_labels = c("1" = "Memory", "2" = "Effector", "3" = "Naive"),
#'                     exclude_labels = "Naive", title_label = "Final clustering",
#'                     output_dir = fig_dir, file_prefix = "round_2")
#' }
pairwise_cluster_de <- function(obj,
                                cluster_labels,
                                exclude_labels = character(0),
                                min_pct = 0.1,
                                log2fc_cutoff = 0.5,
                                padj_cutoff = 0.05,
                                n_labels = 50,
                                exclude_gene_patterns = NULL,
                                title_label,
                                output_dir,
                                file_prefix) {
  
  cluster_ids <- unique(as.character(obj$seurat_clusters))
  missing_labels <- setdiff(cluster_ids, names(cluster_labels))
  if (length(missing_labels) > 0) {
    stop("pairwise_cluster_de(): `cluster_labels` has no label for cluster(s): ",
         paste(missing_labels, collapse = ", "), call. = FALSE)
  }
  
  obj$cluster_name <- unname(cluster_labels[as.character(obj$seurat_clusters)])
  obj <- subset(obj, cells = colnames(obj)[!obj$cluster_name %in% exclude_labels])
  remaining <- sort(unique(obj$cluster_name))
  if (length(remaining) < 2) {
    stop("pairwise_cluster_de(): fewer than 2 clusters left without ",
         paste(exclude_labels, collapse = ", "), ".", call. = FALSE)
  }
  
  Seurat::Idents(obj) <- obj$cluster_name
  results <- list()
  
  for (pair in utils::combn(remaining, 2, simplify = FALSE)) {
    pair_label <- paste(pair[1], "vs", pair[2])
    file_label <- paste0(gsub("[^A-Za-z0-9]", "", pair[1]), "_vs_", gsub("[^A-Za-z0-9]", "", pair[2]))
    message("[pairwise DE] ", pair_label)
    
    de_results <- tryCatch(
      Seurat::FindMarkers(obj, ident.1 = pair[1], ident.2 = pair[2], test.use = "wilcox",
                          min.pct = min_pct, logfc.threshold = log2fc_cutoff, verbose = FALSE),
      error = function(e) {
        message("pairwise_cluster_de(): FindMarkers() failed for ", pair_label, ": ", conditionMessage(e))
        data.frame(p_val = numeric(0), avg_log2FC = numeric(0), pct.1 = numeric(0),
                   pct.2 = numeric(0), p_val_adj = numeric(0))
      }
    )
    if (length(exclude_gene_patterns) > 0) {
      de_results <- de_results[!grepl(paste(exclude_gene_patterns, collapse = "|"), rownames(de_results)), ,
                               drop = FALSE]
    }
    de_results <- de_results[order(de_results$p_val_adj), , drop = FALSE]
    
    volcano <- NULL
    if (nrow(de_results) == 0) {
      message("pairwise_cluster_de(): no DE results for ", pair_label, "; no table or volcano.")
    } else {
      de_table <- gt_de(de_results, title = paste(title_label, "—", pair_label),
                        subtitle = "Pairwise differential expression",
                        pct1_label = paste0("% expr. (", pair[1], ")"),
                        pct2_label = paste0("% expr. (", pair[2], ")"))
      print(de_table)
      save_gt(de_table, paste0(file_prefix, "_", file_label, "_de_table.png"), output_dir)
      
      volcano <- plot_cluster_volcanoes(
        data.frame(gene = rownames(de_results), de_results), subtitle = pair_label,
        log2fc_cutoff = log2fc_cutoff, padj_cutoff = padj_cutoff, n_labels = n_labels
      )[["comparison"]]
      if (!is.null(volcano)) {
        save_plot(volcano, paste0(file_prefix, "_", file_label, "_volcano.png"), output_dir, width = 10, height = 8)
      }
    }
    results[[pair_label]] <- list(de_results = de_results, volcano = volcano)
  }
  
  invisible(results)
}
#' DEGs shared between 2-4 conditions: Venn diagram and workbook
#'
#' Which genes are DEGs in every condition compared, and which in only
#' some. The conditions can be clusters (a split of FindAllMarkers()
#' results), cluster pairs (pairwise_cluster_de()) or the same population
#' in alternative rounds (tools/compare_rounds.R). A gene is a DEG of a
#' condition when its adjusted p-value < `padj_cutoff` and its log2FC
#' passes `log2fc_cutoff` in `direction`.
#'
#' Saves the Venn diagram, a summary table (genes per Venn region), a
#' workbook (Summary; Membership: one row per DEG with a TRUE/FALSE column
#' and the log2FC and adjusted p-value per condition; Shared_all; one
#' Only_<condition> sheet per condition), and one plain gene list per
#' condition (.csv, no header) for external Venn tools.
#'
#' @param de_list Named list of 2-4 DE tables, one per condition, with
#'   `avg_log2FC` and `p_val_adj`, and genes in a `gene` column or as row
#'   names. Names are used as set labels, columns and sheet names (31
#'   characters at most in Excel).
#' @param log2fc_cutoff,padj_cutoff DEG cutoffs. Defaults 0.5, 0.05.
#' @param exclude_gene_patterns Regex patterns of genes left out. Default NULL.
#' @param direction "up" (log2FC > cutoff), "down" (log2FC < -cutoff) or
#'   "both" (|log2FC| > cutoff). With "both", a gene up in one condition
#'   and down in another counts as shared; "up" and "down" compare genes
#'   that change the same way. Appended to file names and titles unless
#'   "both". Default "both".
#' @param title_label Text starting every title. Default `file_prefix`.
#' @param output_dir Directory to save into.
#' @param file_prefix Prefix of every file name.
#'
#' @return Invisibly, a list: `gene_sets` (DEGs per condition), `membership`,
#'   `venn` (ggplot), `excel_path` and `csv_paths`.
#'
#' @examples
#' \dontrun{
#' de_by_cluster <- split(de_results, de_results$cluster)
#' compare_deg_venn(de_by_cluster[c("1", "2", "3")], direction = "up",
#'                  output_dir = fig_dir, file_prefix = "round_1_deg_overlap")
#' }
compare_deg_venn <- function(de_list,
                             log2fc_cutoff = 0.5,
                             padj_cutoff = 0.05,
                             exclude_gene_patterns = NULL,
                             direction = c("both", "up", "down"),
                             title_label = file_prefix,
                             output_dir,
                             file_prefix) {
  
  direction <- match.arg(direction)
  if (!is.list(de_list) || length(de_list) < 2 || length(de_list) > 4) {
    stop("compare_deg_venn(): `de_list` must be a list of 2 to 4 DE tables, got ", length(de_list), ".",
         call. = FALSE)
  }
  if (is.null(names(de_list)) || any(names(de_list) == "") || anyDuplicated(names(de_list))) {
    stop("compare_deg_venn(): every element of `de_list` needs a unique, non-empty name.", call. = FALSE)
  }
  
  if (direction != "both") {
    file_prefix <- paste0(file_prefix, "_", direction)
    title_label <- paste0(title_label, " (", direction, ")")
  }
  direction_text <- switch(direction,
                           both = paste0("|log2FC| > ", log2fc_cutoff),
                           up   = paste0("log2FC > ", log2fc_cutoff),
                           down = paste0("log2FC < -", log2fc_cutoff))
  condition_names <- names(de_list)
  
  # --- DEGs per condition, with their log2FC and adjusted p-value ---------------------
  sig_stats <- lapply(setNames(condition_names, condition_names), function(nm) {
    df <- as.data.frame(de_list[[nm]])
    if (!"gene" %in% names(df)) df$gene <- rownames(df)
    validate_required_cols(df, c("avg_log2FC", "p_val_adj"), "compare_deg_venn")
    if (length(exclude_gene_patterns) > 0) {
      df <- df[!grepl(paste(exclude_gene_patterns, collapse = "|"), df$gene), , drop = FALSE]
    }
    passes_fc <- switch(direction,
                        both = abs(df$avg_log2FC) > log2fc_cutoff,
                        up   = df$avg_log2FC > log2fc_cutoff,
                        down = df$avg_log2FC < -log2fc_cutoff)
    df <- df[!is.na(df$p_val_adj) & df$p_val_adj < padj_cutoff & passes_fc,
             c("gene", "avg_log2FC", "p_val_adj"), drop = FALSE]
    df[!duplicated(df$gene), , drop = FALSE]
  })
  gene_sets <- lapply(sig_stats, `[[`, "gene")
  
  empty_sets <- condition_names[lengths(gene_sets) == 0]
  if (length(empty_sets) > 0) {
    message("compare_deg_venn(): no DEGs (", direction_text, ", padj < ", padj_cutoff, ") for: ",
            paste(empty_sets, collapse = ", "), "; shown as empty sets.")
  }
  
  # --- Venn diagram ------------------------------------------------------------------
  venn_plot <- ggvenn::ggvenn(gene_sets, fill_color = resolve_palette("venn", n = length(gene_sets)),
                              stroke_size = 0.5, set_name_size = 4, text_size = 4) +
    ggplot2::labs(title = paste(title_label, "— DEG overlap"),
                  subtitle = paste0(direction_text, ", padj < ", padj_cutoff))
  print(venn_plot)
  save_plot(venn_plot, paste0(file_prefix, "_venn.png"), output_dir, width = 8, height = 8)
  
  # --- membership: one row per gene that is a DEG in at least one condition ----------------
  membership <- data.frame(gene = sort(unique(unlist(gene_sets))), stringsAsFactors = FALSE)
  for (nm in condition_names) membership[[nm]] <- membership$gene %in% gene_sets[[nm]]
  membership$n_conditions <- rowSums(membership[condition_names])
  # log2FC and adjusted p-value where the gene is a DEG, NA elsewhere
  for (nm in condition_names) {
    idx <- match(membership$gene, sig_stats[[nm]]$gene)
    membership[[paste0(nm, "_log2FC")]] <- sig_stats[[nm]]$avg_log2FC[idx]
    membership[[paste0(nm, "_padj")]]   <- sig_stats[[nm]]$p_val_adj[idx]
  }
  membership <- membership[order(-membership$n_conditions, membership$gene), ]
  
  # --- summary: genes per Venn region ---------------------------------------------------------
  region <- apply(membership[condition_names], 1, function(row) paste(condition_names[row], collapse = " & "))
  summary_table <- if (length(region) > 0) {
    setNames(as.data.frame(table(region), stringsAsFactors = FALSE), c("shared_between", "n_genes"))
  } else {
    data.frame(shared_between = character(0), n_genes = integer(0))
  }
  summary_table <- summary_table[order(-summary_table$n_genes), , drop = FALSE]
  
  summary_gt <- gt_generic(summary_table, title = paste(title_label, "— DEG overlap summary"))
  print(summary_gt)
  save_gt(summary_gt, paste0(file_prefix, "_venn_summary.png"), output_dir)
  
  # --- workbook and plain gene lists --------------------------------------------------------------
  only_sheets <- lapply(condition_names, function(nm) {
    only <- membership[membership$n_conditions == 1 & membership[[nm]],
                       c("gene", paste0(nm, "_log2FC"), paste0(nm, "_padj")), drop = FALSE]
    names(only) <- c("gene", "log2FC", "padj")
    only[order(-abs(only$log2FC)), ]
  })
  names(only_sheets) <- paste0("Only_", condition_names)
  shared_all <- membership[membership$n_conditions == length(condition_names),
                           c("gene", paste0(rep(condition_names, each = 2), c("_log2FC", "_padj"))),
                           drop = FALSE]
  
  excel_path <- save_xlsx(c(list(Summary = summary_table, Membership = membership, Shared_all = shared_all),
                            only_sheets),
                          paste0(file_prefix, "_deg_venn.xlsx"), output_dir)
  
  csv_paths <- setNames(file.path(output_dir, paste0(file_prefix, "_", make.names(condition_names), "_genes.csv")),
                        condition_names)
  for (nm in condition_names) {
    utils::write.table(gene_sets[[nm]], csv_paths[[nm]], row.names = FALSE, col.names = FALSE, quote = FALSE)
  }
  
  invisible(list(gene_sets = gene_sets, membership = membership, venn = venn_plot,
                 excel_path = excel_path, csv_paths = csv_paths))
}
#' Bar chart of DEG counts per condition
#'
#' The number of DEGs in each set of compare_deg_venn()$gene_sets, e.g. to
#' see whether a condition with more cells simply has more DEGs. Whether
#' the extra DEGs change the biology is answered by the Venn membership and
#' by enrichment (enrichment.R).
#'
#' @param gene_sets Named list of gene vectors, e.g. compare_deg_venn()$gene_sets.
#' @param output_dir Directory to save into.
#' @param file_prefix Prefix of the file name.
#' @param title Plot title. Default "<file_prefix> - DEG count by condition".
#'
#' @return Invisibly, the ggplot object.
#'
#' @examples
#' \dontrun{
#' venn <- compare_deg_venn(de_list, output_dir = fig_dir, file_prefix = "Naive_by_round")
#' plot_deg_count_comparison(venn$gene_sets, output_dir = fig_dir, file_prefix = "Naive_by_round")
#' }
plot_deg_count_comparison <- function(gene_sets,
                                      output_dir,
                                      file_prefix,
                                      title = paste(file_prefix, "— DEG count by condition")) {
  
  if (!is.list(gene_sets) || is.null(names(gene_sets)) || any(names(gene_sets) == "")) {
    stop("plot_deg_count_comparison(): `gene_sets` must be a named list of gene vectors, ",
         "e.g. compare_deg_venn()$gene_sets.", call. = FALSE)
  }
  
  count_df <- data.frame(condition = factor(names(gene_sets), levels = names(gene_sets)),
                         n_degs = lengths(gene_sets))
  bar_plot <- ggplot2::ggplot(count_df, ggplot2::aes(x = condition, y = n_degs, fill = condition)) +
    ggplot2::geom_col(width = 0.6, color = "#000000") +
    ggplot2::geom_text(ggplot2::aes(label = n_degs), vjust = -0.4, size = 4) +
    ggplot2::scale_fill_manual(values = resolve_palette("venn", n = nrow(count_df))) +
    ggplot2::labs(title = title, x = NULL, y = "Number of DEGs") +
    ggplot2::theme_minimal() +
    ggplot2::theme(legend.position = "none")
  
  print(bar_plot)
  save_plot(bar_plot, paste0(file_prefix, "_deg_counts.png"), output_dir, width = 6, height = 5)
  invisible(bar_plot)
}
#' One-vs-rest DE results as a workbook, one sheet per cluster
#'
#' Sheets in cluster-ID order, named "<ID>_<label>" (or "cluster_<ID>"
#' without a label), each sorted by adjusted p-value.
#'
#' @param de_results FindAllMarkers() table with `cluster` and `gene`, and
#'   optionally `cluster_label`.
#' @param filename File name (.xlsx).
#' @param output_dir Directory to save into.
#'
#' @return Invisibly, the file path.
write_de_workbook <- function(de_results,
                              filename,
                              output_dir) {
  
  validate_required_cols(de_results, c("cluster", "gene", "avg_log2FC", "p_val_adj"), "write_de_workbook")
  ids <- as.character(de_results$cluster)
  id_order <- unique(ids[order(suppressWarnings(as.numeric(ids)), ids)])
  
  sheets <- lapply(id_order, function(id) {
    df <- de_results[ids == id, , drop = FALSE]
    df <- df[order(df$p_val_adj, -abs(df$avg_log2FC)), , drop = FALSE]
    df[, c("gene", setdiff(names(df), "gene")), drop = FALSE]
  })
  sheet_names <- vapply(id_order, function(id) {
    label <- if ("cluster_label" %in% names(de_results)) de_results$cluster_label[ids == id][1] else NA
    if (is.na(label) || label == id) paste0("cluster_", id) else paste0(id, "_", label)
  }, character(1))
  # unique within Excel's 31-character limit
  names(sheets) <- make.unique(substr(sheet_names, 1, 28), sep = "_")
  
  save_xlsx(sheets, filename, output_dir)
}
