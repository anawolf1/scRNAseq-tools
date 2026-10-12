# ------------------------------------------------------------------------------
# clustering_evaluation.R
#
# Evaluation of clusterings: parameter sweeps, stability under
# subsampling, marker support, and a closer look at single candidates.
# Expects objects preprocessed with preprocess_seurat() (clustering.R).
# Sweep heatmaps use plot_sweep_tile() (shared plots_generic.R).
#
# Functions
# empty_de_table()                    empty FindAllMarkers()-shaped table
# safe_find_all_markers()             FindAllMarkers(), or an empty table if it fails
# evaluate_clustering()               score every combination of a parameter grid
# per_cluster_stability()             per-cluster Jaccard stability under subsampling
# count_real_markers()                significant positive markers per cluster
# evaluate_all_candidates_markers()   marker support and marker-column prevalence per candidate
# check_cluster_size_imbalance()      correlation of cluster-size imbalance with ARI stability
# plot_size_imbalance()               ARI stability vs. cluster-size imbalance
# combine_round_sweeps()              several rounds' sweeps, worst value per combination
# inspect_clustering_candidate()      one candidate: clusters, stability, top markers
# check_candidate_qc()                QC metrics by cluster for one candidate
# gt_markers_by_cluster()             top markers of several clusters as one gt table
# ------------------------------------------------------------------------------

#' Empty FindAllMarkers()-shaped table
#'
#' @return A zero-row data.frame with `cluster`, `gene`, `p_val`,
#'   `avg_log2FC`, `pct.1`, `pct.2` and `p_val_adj`.
empty_de_table <- function() {
  data.frame(cluster = character(0), gene = character(0), p_val = numeric(0),
             avg_log2FC = numeric(0), pct.1 = numeric(0), pct.2 = numeric(0), p_val_adj = numeric(0))
}
#' FindAllMarkers(), or an empty table if it fails
#'
#' A failure (e.g. a cluster too small to test) is reported as a message
#' and gives the same empty table as no markers, so a sweep continues.
#'
#' @param obj Seurat object, identities set to the clusters.
#' @param ... Passed to Seurat::FindAllMarkers().
#'
#' @return The FindAllMarkers() table, or empty_de_table().
safe_find_all_markers <- function(obj,
                                  ...) {
  de <- tryCatch(Seurat::FindAllMarkers(obj, ...), error = function(e) {
    message("safe_find_all_markers(): FindAllMarkers() failed: ", conditionMessage(e))
    NULL
  })
  if (is.null(de) || nrow(de) == 0) empty_de_table() else de
}
#' Score every combination of a clustering parameter grid
#'
#' Clusters every combination of dims, resolution, k_param and prune_snn,
#' and scores each on separation (mean silhouette, Calinski-Harabasz),
#' clusters without support (fraction of clusters with fewer than
#' `min_markers_per_cluster` positive markers or fewer than
#' `min_cluster_size` cells) and stability (mean ARI between the clusters
#' of all cells and those of subsamples). For narrowing a grid down to
#' candidates; evaluate_all_candidates_markers() and
#' inspect_clustering_candidate() look at those in more detail.
#'
#' Each subsample is reclustered on the reduction computed from all cells,
#' so the ARI measures the stability of the clustering step, not of the
#' whole preprocessing. Every combination uses the same subsamples (`seed` +
#' repetition), so only the swept parameters differ between combinations.
#'
#' @param obj Seurat object with `reduction` computed.
#' @param dims_list,resolution_list Values of dims and resolution to test.
#' @param k_list,prune_list Values of k_param and prune_snn to test.
#'   Defaults 20 and 1/15.
#' @param algorithm FindClusters() algorithm (1 = Louvain, 4 = Leiden).
#' @param reduction Reduction to cluster on. Default "pca".
#' @param n_reps Subsamples per combination. Default 10; with only 10, a
#'   small ARI difference (e.g. 0.84 vs. 0.86) is not meaningful.
#' @param frac Fraction of cells in each subsample. Default 0.8.
#' @param marker_min_pct,marker_log2fc_cutoff,marker_padj_cutoff Marker
#'   thresholds (FindAllMarkers() min.pct and logfc.threshold, positive
#'   markers only). Defaults 0.01, 0.25, 0.05.
#' @param min_markers_per_cluster A cluster with fewer markers is invalid.
#'   Default 1.
#' @param min_cluster_size A cluster with fewer cells is invalid. Default 10.
#' @param exclude_gene_patterns Regex patterns of genes that never count as
#'   markers (e.g. "^ENSG"). Default NULL.
#' @param seed Seed for clustering and subsampling. Default 123.
#'
#' @return A data.frame, one row per combination: dims, resolution, k,
#'   prune, n_clusters, min/median/max_cluster_size, mean_silhouette,
#'   calinski_harabasz, frac_invalid_clusters, mean_ari_stability,
#'   sd_ari_stability. Scores are NA where fewer than 2 clusters formed.
#'
#' @examples
#' \dontrun{
#' sweep <- evaluate_clustering(obj, dims_list = c(10, 20, 30),
#'                              resolution_list = c(0.4, 0.8, 1.2), algorithm = 4)
#' }
evaluate_clustering <- function(obj,
                                dims_list,
                                resolution_list,
                                k_list = 20,
                                prune_list = 1 / 15,
                                algorithm,
                                reduction = "pca",
                                n_reps = 10,
                                frac = 0.8,
                                marker_min_pct = 0.01,
                                marker_log2fc_cutoff = 0.25,
                                marker_padj_cutoff = 0.05,
                                min_markers_per_cluster = 1,
                                min_cluster_size = 10,
                                exclude_gene_patterns = NULL,
                                seed = 123) {

  if (!is.numeric(dims_list) || length(dims_list) == 0 || any(dims_list < 1) ||
      any(dims_list != as.integer(dims_list))) {
    stop("evaluate_clustering(): `dims_list` must contain positive integers.", call. = FALSE)
  }
  validate_reduction_dims(obj, reduction, max(dims_list), "evaluate_clustering")
  if (!all(resolution_list > 0)) {
    stop("evaluate_clustering(): every resolution must be > 0.", call. = FALSE)
  }

  full_embedding <- Seurat::Embeddings(obj, reduction = reduction)
  results <- list()

  for (dims in dims_list) {

    # the distance matrix depends on dims only
    dist_mat <- stats::dist(full_embedding[, 1:dims])

    for (res in resolution_list) for (k in k_list) for (prune in prune_list) {

      obj_tmp <- cluster_seurat(obj, resolution = res, algorithm = algorithm, reduction = reduction,
                                dims = dims, k_param = k, prune_snn = prune, seed = seed,
                                verbose = FALSE)
      clusters_full <- obj_tmp$seurat_clusters
      cluster_sizes <- table(clusters_full)
      row <- data.frame(
        dims = dims, resolution = res, k = k, prune = prune,
        n_clusters = length(cluster_sizes),
        min_cluster_size = min(cluster_sizes), median_cluster_size = stats::median(cluster_sizes),
        max_cluster_size = max(cluster_sizes),
        mean_silhouette = NA_real_, calinski_harabasz = NA_real_,
        frac_invalid_clusters = NA_real_, mean_ari_stability = NA_real_, sd_ari_stability = NA_real_
      )

      if (length(cluster_sizes) >= 2) {
        sil <- cluster::silhouette(as.integer(clusters_full), dist_mat)
        row$mean_silhouette <- round(mean(sil[, "sil_width"]), 3)
        row$calinski_harabasz <- round(tryCatch(
          fpc::calinhara(full_embedding[, 1:dims], as.integer(clusters_full)),
          error = function(e) NA_real_
        ), 1)

        Seurat::Idents(obj_tmp) <- clusters_full
        markers <- safe_find_all_markers(obj_tmp, min.pct = marker_min_pct,
                                         logfc.threshold = marker_log2fc_cutoff,
                                         only.pos = TRUE, verbose = FALSE)
        n_markers <- count_real_markers(markers, cluster_levels = levels(clusters_full),
                                        log2fc_cutoff = marker_log2fc_cutoff,
                                        padj_cutoff = marker_padj_cutoff,
                                        exclude_gene_patterns = exclude_gene_patterns)$n_real_markers
        # a cluster is invalid when it fails either criterion
        row$frac_invalid_clusters <- round(mean(n_markers < min_markers_per_cluster |
                                                  cluster_sizes[levels(clusters_full)] < min_cluster_size), 3)

        ari <- vapply(seq_len(n_reps), function(rep_i) {
          set.seed(seed + rep_i)
          sub_cells <- sample(colnames(obj_tmp), size = floor(frac * ncol(obj_tmp)))
          obj_sub <- cluster_seurat(subset(obj_tmp, cells = sub_cells), resolution = res,
                                    algorithm = algorithm, reduction = reduction, dims = dims,
                                    k_param = k, prune_snn = prune, seed = seed, verbose = FALSE)
          aricode::ARI(as.character(clusters_full[colnames(obj_sub)]), as.character(obj_sub$seurat_clusters))
        }, numeric(1))
        row$mean_ari_stability <- round(mean(ari), 3)
        row$sd_ari_stability <- round(stats::sd(ari), 3)
      }

      message("dims = ", dims, " | res = ", res, " | k = ", k, " | prune = ", round(prune, 3),
              " | clusters = ", row$n_clusters, " | silhouette = ", row$mean_silhouette,
              " | frac invalid = ", row$frac_invalid_clusters, " | ARI = ", row$mean_ari_stability)
      results[[length(results) + 1]] <- row
    }
  }
  do.call(rbind, results)
}
#' Per-cluster Jaccard stability under subsampling
#'
#' Subsamples the cells, reclusters them with cluster_seurat() and the
#' given parameters, and scores each original cluster by its best Jaccard
#' overlap with a subsample cluster. Complements the global ARI of
#' evaluate_clustering(), which can hide one unstable cluster. The score
#' detects a cluster that splits; two clusters that merge can both still
#' score well, which the ARI shows. The parameters must be those that made
#' `obj`'s clusters.
#'
#' @param obj Seurat object, clustered (`seurat_clusters`), with
#'   `reduction` computed.
#' @param dims,resolution,algorithm,k_param,prune_snn,reduction Clustering
#'   parameters (see cluster_seurat()).
#' @param n_reps Subsamples. Default 10.
#' @param frac Fraction of cells in each subsample. Default 0.8.
#' @param seed Seed for clustering and subsampling. Default 123.
#'
#' @return A data.frame, one row per cluster: cluster, mean_jaccard,
#'   sd_jaccard, n_cells.
#'
#' @examples
#' \dontrun{
#' per_cluster_stability(obj, dims = 12, resolution = 0.5, algorithm = 4)
#' }
per_cluster_stability <- function(obj,
                                  dims,
                                  resolution,
                                  algorithm,
                                  k_param = 20,
                                  prune_snn = 1 / 15,
                                  reduction = "pca",
                                  n_reps = 10,
                                  frac = 0.8,
                                  seed = 123) {

  validate_required_cols(obj[[]], "seurat_clusters", "per_cluster_stability")
  clusters_full <- obj$seurat_clusters
  cluster_levels <- levels(clusters_full)

  jaccard <- matrix(NA_real_, nrow = n_reps, ncol = length(cluster_levels),
                    dimnames = list(NULL, cluster_levels))

  for (i in seq_len(n_reps)) {
    set.seed(seed + i)
    sub_cells <- sample(colnames(obj), size = floor(frac * ncol(obj)))
    sub_clusters <- cluster_seurat(subset(obj, cells = sub_cells), resolution = resolution,
                                   algorithm = algorithm, reduction = reduction, dims = dims,
                                   k_param = k_param, prune_snn = prune_snn, seed = seed,
                                   verbose = FALSE)$seurat_clusters

    for (cluster in cluster_levels) {
      orig_cells <- intersect(names(clusters_full)[clusters_full == cluster], sub_cells)
      if (length(orig_cells) == 0) next
      jaccard[i, cluster] <- max(vapply(levels(sub_clusters), function(sub_cluster) {
        sub_cells_k <- names(sub_clusters)[sub_clusters == sub_cluster]
        length(intersect(orig_cells, sub_cells_k)) / length(union(orig_cells, sub_cells_k))
      }, numeric(1)))
    }
  }

  data.frame(
    cluster      = cluster_levels,
    mean_jaccard = round(colMeans(jaccard, na.rm = TRUE), 3),
    sd_jaccard   = round(apply(jaccard, 2, stats::sd, na.rm = TRUE), 3),
    n_cells      = as.integer(table(clusters_full)[cluster_levels]),
    row.names    = NULL
  )
}
#' Significant positive markers per cluster
#'
#' Counts the markers of each cluster with adjusted p-value <
#' `padj_cutoff` and avg_log2FC >= `log2fc_cutoff` (positive markers, as
#' from FindAllMarkers(only.pos = TRUE)), leaving out `exclude_gene_patterns`.
#'
#' @param markers_df FindAllMarkers() table (cluster, gene, p_val_adj,
#'   avg_log2FC); may have zero rows.
#' @param cluster_levels Clusters to report. NULL (default) uses the
#'   clusters in `markers_df`; pass the object's cluster levels so a cluster
#'   without markers is reported with 0.
#' @param log2fc_cutoff Minimum avg_log2FC; use the logfc.threshold the
#'   table was made with.
#' @param padj_cutoff Adjusted p-value cutoff. Default 0.05.
#' @param exclude_gene_patterns Regex patterns of genes not counted.
#'   Default NULL.
#'
#' @return A data.frame: cluster, n_real_markers.
#'
#' @examples
#' \dontrun{
#' count_real_markers(markers, cluster_levels = levels(obj$seurat_clusters), log2fc_cutoff = 0.25)
#' }
count_real_markers <- function(markers_df,
                               cluster_levels = NULL,
                               log2fc_cutoff,
                               padj_cutoff = 0.05,
                               exclude_gene_patterns = NULL) {

  keep <- markers_df$p_val_adj < padj_cutoff & markers_df$avg_log2FC >= log2fc_cutoff
  if (length(exclude_gene_patterns) > 0) {
    keep <- keep & !grepl(paste(exclude_gene_patterns, collapse = "|"), markers_df$gene)
  }
  all_clusters <- cluster_levels %||% unique(as.character(markers_df$cluster))
  counts <- table(factor(markers_df$cluster[keep], levels = all_clusters))

  data.frame(cluster = names(counts), n_real_markers = as.integer(counts), row.names = NULL)
}
#' Marker support and marker-column prevalence per candidate
#'
#' Reclusters each candidate (e.g. rows of evaluate_clustering()'s output),
#' counts each cluster's markers, and optionally reports how strongly a
#' known population (a level of a metadata column, e.g. a protein gate)
#' concentrates in one cluster: the highest fraction of any cluster made
#' up of that level (all cells of the cluster, NA counted as not the
#' level). Checking known populations this way makes the choice
#' partly confirm the expectation rather than test it; report that it was
#' used.
#'
#' @param obj Seurat object with `reduction` computed.
#' @param candidates Data.frame with dims, resolution, k, prune, n_clusters
#'   and mean_ari_stability columns.
#' @param algorithm FindClusters() algorithm.
#' @param reduction Reduction to cluster on. Default "pca".
#' @param marker_min_pct,marker_log2fc_cutoff,marker_padj_cutoff Marker
#'   thresholds. Defaults 0.01, 0.25, 0.05.
#' @param marker_count_threshold Clusters with fewer markers are counted in
#'   `n_clusters_below_threshold`. Default 10.
#' @param check_cols Optional named list, metadata column -> level, e.g.
#'   list(Naive = "Naive", CD127 = "CD127-"). Default NULL.
#' @param exclude_gene_patterns Regex patterns of genes not counted as
#'   markers. Default NULL.
#' @param seed Seed for clustering. Default 123.
#'
#' @return Invisibly, a list: `summary` (one row per candidate: n_clusters,
#'   dims, resolution, k, prune, mean_ari_stability, min_real_markers,
#'   median_real_markers, n_clusters_below_threshold and, per `check_cols`
#'   entry, max_<name>_prevalence) and `objects` (the clustered object of
#'   each candidate, named dims<d>_k<k>_prune<p>_res<r>).
#'
#' @examples
#' \dontrun{
#' evaluate_all_candidates_markers(obj, candidates, algorithm = 4,
#'                                 check_cols = list(Naive = "Naive"))
#' }
evaluate_all_candidates_markers <- function(obj,
                                            candidates,
                                            algorithm,
                                            reduction = "pca",
                                            marker_min_pct = 0.01,
                                            marker_log2fc_cutoff = 0.25,
                                            marker_padj_cutoff = 0.05,
                                            marker_count_threshold = 10,
                                            check_cols = NULL,
                                            exclude_gene_patterns = NULL,
                                            seed = 123) {

  validate_required_cols(candidates, c("dims", "resolution", "k", "prune", "n_clusters", "mean_ari_stability"),
                         "evaluate_all_candidates_markers")
  validate_required_cols(obj[[]], names(check_cols), "evaluate_all_candidates_markers",
                         hint = "expects every check_cols name to be a metadata column.")

  summary_rows <- list()
  clustered_objs <- list()

  for (i in seq_len(nrow(candidates))) {
    cand <- candidates[i, ]
    combo_name <- paste0("dims", cand$dims, "_k", cand$k,
                         "_prune", format(cand$prune, trim = TRUE, scientific = FALSE),
                         "_res", cand$resolution)

    obj_tmp <- cluster_seurat(obj, resolution = cand$resolution, algorithm = algorithm,
                              reduction = reduction, dims = cand$dims, k_param = cand$k,
                              prune_snn = cand$prune, seed = seed, verbose = FALSE)
    markers <- safe_find_all_markers(obj_tmp, min.pct = marker_min_pct,
                                     logfc.threshold = marker_log2fc_cutoff,
                                     only.pos = TRUE, verbose = FALSE)
    n_markers <- count_real_markers(markers, cluster_levels = levels(obj_tmp$seurat_clusters),
                                    log2fc_cutoff = marker_log2fc_cutoff,
                                    padj_cutoff = marker_padj_cutoff,
                                    exclude_gene_patterns = exclude_gene_patterns)$n_real_markers

    row <- data.frame(
      n_clusters = cand$n_clusters, dims = cand$dims, resolution = cand$resolution,
      k = cand$k, prune = cand$prune, mean_ari_stability = cand$mean_ari_stability,
      min_real_markers = min(n_markers), median_real_markers = stats::median(n_markers),
      n_clusters_below_threshold = sum(n_markers < marker_count_threshold)
    )
    for (col_name in names(check_cols)) {
      # fraction of all cells of each cluster; a cell with NA does not count as the level
      is_level <- obj_tmp[[col_name, drop = TRUE]] %in% check_cols[[col_name]]
      prevalence <- tapply(is_level, obj_tmp$seurat_clusters, mean)
      row[[paste0("max_", col_name, "_prevalence")]] <- round(max(prevalence, na.rm = TRUE), 3)
    }

    summary_rows[[i]] <- row
    clustered_objs[[combo_name]] <- obj_tmp
  }

  summary_table <- do.call(rbind, summary_rows)
  summary_table <- summary_table[order(summary_table$n_clusters, -summary_table$min_real_markers), ]
  invisible(list(summary = summary_table, objects = clustered_objs))
}
#' Correlation of cluster-size imbalance with ARI stability across a sweep
#'
#' One large cluster with a few small ones is stable under subsampling
#' almost by construction, so a high ARI can reflect imbalance rather than
#' a good partition. Reports the correlations of the size ratio (largest /
#' smallest cluster) and the number of clusters with the ARI, and the
#' partial correlation of the size ratio with the ARI given the number of
#' clusters. Descriptive; not used to select candidates.
#'
#' @param sweep_df evaluate_clustering() output.
#'
#' @return A one-row data.frame: n_rows, corr_ratio_ari, corr_nclusters_ari,
#'   corr_nclusters_ratio, partial_corr_ratio_ari (NA with fewer than 3
#'   rows with an ARI, or a single number of clusters).
#'
#' @examples
#' \dontrun{
#' gt_generic(check_cluster_size_imbalance(sweep_results), title = "Size imbalance vs. ARI")
#' }
check_cluster_size_imbalance <- function(sweep_df) {

  diag_rows <- sweep_df[!is.na(sweep_df$mean_ari_stability), ]
  size_ratio <- diag_rows$max_cluster_size / diag_rows$min_cluster_size

  if (nrow(diag_rows) < 3 || length(unique(diag_rows$n_clusters)) < 2) {
    return(data.frame(n_rows = nrow(diag_rows), corr_ratio_ari = NA_real_, corr_nclusters_ari = NA_real_,
                      corr_nclusters_ratio = NA_real_, partial_corr_ratio_ari = NA_real_))
  }

  r_ratio_ari       <- stats::cor(size_ratio, diag_rows$mean_ari_stability)
  r_nclusters_ari   <- stats::cor(diag_rows$n_clusters, diag_rows$mean_ari_stability)
  r_nclusters_ratio <- stats::cor(diag_rows$n_clusters, size_ratio)
  partial_r <- (r_ratio_ari - r_nclusters_ari * r_nclusters_ratio) /
    sqrt((1 - r_nclusters_ari^2) * (1 - r_nclusters_ratio^2))

  data.frame(
    n_rows                 = nrow(diag_rows),
    corr_ratio_ari         = round(r_ratio_ari, 3),
    corr_nclusters_ari     = round(r_nclusters_ari, 3),
    corr_nclusters_ratio   = round(r_nclusters_ratio, 3),
    partial_corr_ratio_ari = round(partial_r, 3)
  )
}
#' ARI stability vs. cluster-size imbalance
#'
#' The relationship of check_cluster_size_imbalance() as a scatter plot,
#' colored by the number of clusters.
#'
#' @param sweep_df evaluate_clustering() output.
#' @param title,subtitle Plot title/subtitle. Default NULL.
#'
#' @return A ggplot object, or NULL (with a message) with fewer than 3 rows
#'   with an ARI or a single number of clusters.
#'
#' @examples
#' \dontrun{
#' p <- plot_size_imbalance(sweep_results, title = "ARI stability vs. cluster-size imbalance")
#' }
plot_size_imbalance <- function(sweep_df,
                                title = NULL,
                                subtitle = NULL) {

  diag_rows <- sweep_df[!is.na(sweep_df$mean_ari_stability), ]
  diag_rows$size_ratio <- diag_rows$max_cluster_size / diag_rows$min_cluster_size

  if (nrow(diag_rows) < 3 || length(unique(diag_rows$n_clusters)) < 2) {
    message("plot_size_imbalance(): fewer than 3 rows with an ARI, or a single number of clusters; no plot.")
    return(NULL)
  }

  ggplot2::ggplot(diag_rows, ggplot2::aes(x = size_ratio, y = mean_ari_stability, color = n_clusters)) +
    ggplot2::geom_point(size = 3, alpha = 0.85) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_color_viridis_c(option = resolve_palette("qc_scatter"), name = "n_clusters") +
    ggplot2::labs(x = "Largest / smallest cluster (log scale)", y = "Mean ARI stability",
                  title = title, subtitle = subtitle) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(color = "#666666"),
                   panel.grid    = ggplot2::element_blank(),
                   legend.title  = ggplot2::element_text(size = 9))
}
#' Several rounds' sweeps, worst value per combination
#'
#' For choosing one dims/resolution for rounds that share parameters (e.g.
#' alternative exclusions): their sweeps, saved by
#' tools/select_clustering_parameters.R, merged on the combinations all of
#' them tested, with the worst value of each metric across the rounds.
#'
#' @param sweeps Named list of evaluate_clustering() tables, one per round.
#'
#' @return A data.frame, one row per combination (dims, resolution, k,
#'   prune): per round mean_ari_stability, frac_invalid_clusters and
#'   n_clusters (suffixed with the round name); min_ari (NA if any round's
#'   is NA), max_invalid_frac, and same_n_clusters (TRUE when every round
#'   has the same number of clusters). Sorted by min_ari, highest first.
combine_round_sweeps <- function(sweeps) {
  param_cols  <- c("dims", "resolution", "k", "prune")
  metric_cols <- c("mean_ari_stability", "frac_invalid_clusters", "n_clusters")

  tables <- Map(function(df, round_name) {
    validate_required_cols(df, c(param_cols, metric_cols), "combine_round_sweeps")
    df <- df[, c(param_cols, metric_cols)]
    names(df) <- c(param_cols, paste0(metric_cols, "_", round_name))
    df
  }, sweeps, names(sweeps))
  combined <- Reduce(function(a, b) merge(a, b, by = param_cols), tables)
  if (nrow(combined) == 0) {
    stop("combine_round_sweeps(): the sweeps share no parameter combination; sweep every round on ",
         "the same grid.", call. = FALSE)
  }

  per_round <- function(metric) as.matrix(combined[, paste0(metric, "_", names(sweeps)), drop = FALSE])
  combined$min_ari          <- apply(per_round("mean_ari_stability"), 1, min)
  combined$max_invalid_frac <- apply(per_round("frac_invalid_clusters"), 1, max)
  combined$same_n_clusters  <- apply(per_round("n_clusters"), 1, function(x) length(unique(x)) == 1)
  combined[order(-combined$min_ari, combined$max_invalid_frac), , drop = FALSE]
}
#' One candidate: clusters, stability and top markers
#'
#' For a short list of candidates picked from a sweep (see
#' tools/inspect_clustering_candidates.R): clusters with the given
#' parameters, computes per-cluster stability, and returns the top positive
#' markers of each cluster.
#'
#' @param obj Seurat object with `reduction` computed.
#' @param dims,resolution,algorithm,k_param,prune_snn,reduction Clustering
#'   parameters (see cluster_seurat()). Use the sweep's values, so the
#'   candidate reproduces its sweep result.
#' @param n_reps,frac Passed to per_cluster_stability(). Defaults 10, 0.8.
#' @param marker_min_pct,marker_log2fc_cutoff,marker_padj_cutoff Marker
#'   thresholds. Defaults 0.1, 0.5, 0.05.
#' @param exclude_gene_patterns Regex patterns of genes left out of the
#'   markers. Default NULL.
#' @param top_n_markers Markers kept per cluster, by avg_log2FC. Default 15.
#' @param seed Seed for clustering and subsampling. Default 123.
#'
#' @return A list: `obj` (clustered), `stability` (per_cluster_stability()),
#'   `markers` (top markers per cluster, for gt_markers_by_cluster()), and
#'   the parameters used (`dims`, `resolution`, `algorithm`, `k_param`,
#'   `prune_snn`).
#'
#' @examples
#' \dontrun{
#' result <- inspect_clustering_candidate(obj, dims = 15, resolution = 0.4, algorithm = 4)
#' gt_markers_by_cluster(result$markers, title = "dims15_res0.4")
#' }
inspect_clustering_candidate <- function(obj,
                                         dims,
                                         resolution,
                                         algorithm,
                                         k_param = 20,
                                         prune_snn = 1 / 15,
                                         reduction = "pca",
                                         n_reps = 10,
                                         frac = 0.8,
                                         marker_min_pct = 0.1,
                                         marker_log2fc_cutoff = 0.5,
                                         marker_padj_cutoff = 0.05,
                                         exclude_gene_patterns = NULL,
                                         top_n_markers = 15,
                                         seed = 123) {

  trial_obj <- cluster_seurat(obj, resolution = resolution, algorithm = algorithm, reduction = reduction,
                              dims = dims, k_param = k_param, prune_snn = prune_snn, seed = seed,
                              verbose = FALSE)
  stability <- per_cluster_stability(trial_obj, dims = dims, resolution = resolution, algorithm = algorithm,
                                     k_param = k_param, prune_snn = prune_snn, reduction = reduction,
                                     n_reps = n_reps, frac = frac, seed = seed)

  Seurat::Idents(trial_obj) <- trial_obj$seurat_clusters
  markers <- safe_find_all_markers(trial_obj, test.use = "wilcox", min.pct = marker_min_pct,
                                   logfc.threshold = marker_log2fc_cutoff, only.pos = TRUE, verbose = FALSE)
  markers <- markers[markers$p_val_adj < marker_padj_cutoff, , drop = FALSE]
  if (length(exclude_gene_patterns) > 0) {
    markers <- markers[!grepl(paste(exclude_gene_patterns, collapse = "|"), markers$gene), , drop = FALSE]
  }
  markers <- markers %>%
    dplyr::group_by(cluster) %>%
    dplyr::slice_max(order_by = avg_log2FC, n = top_n_markers, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(cluster, gene, avg_log2FC, p_val_adj) %>%
    dplyr::arrange(cluster, dplyr::desc(avg_log2FC)) %>%
    dplyr::mutate(cluster = as.character(cluster))

  list(obj = trial_obj, stability = stability, markers = markers, dims = dims, resolution = resolution,
       algorithm = algorithm, k_param = k_param, prune_snn = prune_snn)
}
#' QC metrics by cluster for one candidate
#'
#' Median of each QC metric per cluster, and a violin plot: a first look at
#' whether a cluster differs mainly in sequencing depth or mitochondrial
#' content, i.e. may be technical. Descriptive.
#'
#' @param obj Seurat object, clustered (`seurat_clusters`).
#' @param qc_metrics QC metric columns, e.g.
#'   build_qc_metric_info(params$percent_features)$qc_metrics.
#' @param label Candidate label for the titles, e.g. "dims10_res0.4".
#'
#' @return A list: `table` (n_cells and the median of each metric per
#'   cluster, formatted with format_qc_values()) and `plot` (violin plot).
#'
#' @examples
#' \dontrun{
#' qc <- check_candidate_qc(inspection_results[["dims8_res0.6"]]$obj,
#'                          qc_metrics = qc_info$qc_metrics, label = "dims8_res0.6")
#' }
check_candidate_qc <- function(obj,
                               qc_metrics,
                               label = "candidate") {

  validate_required_cols(obj[[]], c("seurat_clusters", qc_metrics), "check_candidate_qc")

  qc_table <- obj[[]] %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(n_cells = dplyr::n(),
                     dplyr::across(dplyr::all_of(qc_metrics), \(x) stats::median(x, na.rm = TRUE)),
                     .groups = "drop")

  cluster_levels <- levels(obj$seurat_clusters)
  qc_violin <- plot_violin(
    obj, features = qc_metrics, group_col = "seurat_clusters",
    group_colors = setNames(resolve_palette("cluster", n = length(cluster_levels)), cluster_levels),
    show_points = FALSE, title = paste(label, "— QC metrics by cluster"), ncol = 2
  )

  list(table = format_qc_values(qc_table, qc_metrics), plot = qc_violin)
}
#' Top markers of several clusters as one gt table
#'
#' One section per cluster; shows the rows given, e.g.
#' inspect_clustering_candidate()$markers.
#'
#' @param markers_df Data.frame with cluster, gene, avg_log2FC and p_val_adj.
#' @param title Table title.
#' @param subtitle Table subtitle.
#'
#' @return A gt table object.
#'
#' @examples
#' \dontrun{
#' gt_markers_by_cluster(result$markers, title = "dims15_res0.4")
#' }
gt_markers_by_cluster <- function(markers_df,
                                  title,
                                  subtitle = "Top markers per cluster, by average log2 fold change") {

  validate_required_cols(markers_df, c("cluster", "gene", "avg_log2FC", "p_val_adj"), "gt_markers_by_cluster")
  if (nrow(markers_df) == 0) {
    stop("gt_markers_by_cluster(): no markers to show; none passed the marker thresholds.", call. = FALSE)
  }

  markers_df %>%
    gt::gt(groupname_col = "cluster") %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::cols_label(gene = "Gene", avg_log2FC = "Log2FC", p_val_adj = "Adj. p-value") %>%
    gt::fmt_number(columns = avg_log2FC, decimals = 3) %>%
    gt::fmt_scientific(columns = p_val_adj, decimals = 2) %>%
    gt::tab_options(
      table.font.size            = gt::px(12),
      heading.title.font.size    = gt::px(17),
      heading.subtitle.font.size = gt::px(11),
      data_row.padding           = gt::px(5),
      row_group.padding          = gt::px(7)
    ) %>%
    gt::opt_row_striping() %>%
    gt::cols_align(align = "center")
}
