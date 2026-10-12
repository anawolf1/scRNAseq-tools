# ------------------------------------------------------------------------------
# select_clustering_parameters.R (workflows/02_population_refinement/tools)
#
# Screens clustering parameters (dims x resolution x k_param x prune_snn)
# for one round, on the preprocessing the round will use. Run by hand, once
# per round (or per subset to check). It records nothing: put the chosen
# values, and the printed `sweep` record, in config/decisions.R.
#
# Staged sweeps are faster than one large grid:
#   1. dims x resolution, k_param and prune_snn fixed; pick a pair
#   2. that pair, several k_param
#   3. if needed, several prune_snn
#
# Output (results/02_population_refinement/exploration/clustering_sweeps/):
#   sweep_<label>_full.png, _candidates_markers.png, metric heatmaps, size
#   imbalance diagnostics, and sweep_<label>.rds (the sweep table, which
#   inspect_clustering_candidates.R can read)
#
# Sections:
#   0. Setup
#   1. Load and preprocess
#   2. Parameter sweep
#   3. Candidate filtering
#   4. Marker support
#   5. Sweep heatmaps
#   6. Cluster size imbalance
#   7. One choice for several rounds [optional]
#
# Usage:
#   source(here::here("workflows", "02_population_refinement", "tools", "select_clustering_parameters.R"))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------
source(here::here("workflows", "02_population_refinement", "setup.R"))
# the rounds, if written yet: cluster labels of round inputs, and
# sweep_for_round below
decision_file <- file.path(workflow_dir, "config", "decisions.R")
if (file.exists(decision_file)) source(decision_file)
load_libraries("cluster", "fpc", "aricode", load_core = FALSE)   # silhouette, Calinski-Harabasz, ARI

# ------------------------------------------------------------------------------
# 1. Load and preprocess (edit per use)
# ------------------------------------------------------------------------------
# what to sweep: a checkpoint name (config.R `checkpoints`, e.g.
# "qc_filtered"), a round ("round_1"), or a file path
sweep_checkpoint <- "qc_filtered"
sweep_algorithm  <- 1   # 1 = Louvain, 4 = Leiden

# [optional] a subset of cells, e.g. one cluster, to see whether it splits
# further; any metadata column. NULL = all cells.
subset_column <- NULL   # e.g. "seurat_clusters" or "cluster_label"
subset_values <- NULL   # e.g. "2"

# [optional] remove params$clustering$exclude_gene_patterns genes first:
# TRUE when the round being screened for has exclude_genes = TRUE
apply_gene_exclusion <- FALSE

# [optional] screen for a round written in config/decisions.R (e.g. an
# alternative exclusion): uses its input, its steps and its exclude_genes
# instead of the settings above; its dims/resolution are not needed yet
sweep_for_round <- NULL   # e.g. "round_2_half"

sweep_steps <- NULL
if (!is.null(sweep_for_round)) {
  target_round <- params$population_refinement$rounds[[sweep_for_round]]
  if (is.null(target_round)) stop("'", sweep_for_round, "' is not a round in config/decisions.R.")
  sweep_checkpoint     <- target_round$input
  sweep_steps          <- target_round$steps
  apply_gene_exclusion <- isTRUE(target_round$exclude_genes)
  message("Screening for ", sweep_for_round, ": input ", sweep_checkpoint, ", ",
          length(sweep_steps), " step(s), gene exclusion ", apply_gene_exclusion)
}

# label of the output files; a subset gets its own, so it does not
# overwrite the sweep of all cells
sweep_label <- paste(c(sweep_for_round %||% basename(sweep_checkpoint),
                       if (!is.null(subset_column)) c(subset_column, paste(subset_values, collapse = "-"))),
                     collapse = "_")
sweep_label <- gsub("[^A-Za-z0-9_.-]+", "_", sub("\\.rds$", "", sweep_label))

# preprocessed with the round's settings (params$clustering)
obj <- prepare_clustering_input(
  sweep_checkpoint, data_processed_dir,
  subset_column         = subset_column,
  subset_values         = subset_values,
  steps                 = sweep_steps,
  exclude_gene_patterns = if (apply_gene_exclusion) params$clustering$exclude_gene_patterns,
  nfeatures             = params$clustering$nfeatures,
  vars_to_regress       = params$clustering$vars_to_regress,
  seed                  = params$clustering$random_seed,
  rounds                = params$population_refinement$rounds
)
print(Seurat::ElbowPlot(obj, ndims = 50))

# look at the elbow plot before setting dims_list below

# ------------------------------------------------------------------------------
# 2. Parameter sweep (edit per use)
# ------------------------------------------------------------------------------
# every combination is clustered n_reps + 1 times: keep the grid small
dims_list       <- c(6, 8, 10, 12, 15, 18, 20)
resolution_list <- seq(0.3, 1.0, by = 0.1)
k_list          <- params$clustering$k_param
prune_list      <- params$clustering$prune_snn

ps <- params$parameter_selection
sweep_results <- evaluate_clustering(
  obj,
  dims_list               = dims_list,
  resolution_list         = resolution_list,
  k_list                  = k_list,
  prune_list              = prune_list,
  algorithm               = sweep_algorithm,
  n_reps                  = ps$stability_n_reps,
  frac                    = ps$stability_frac,
  marker_min_pct          = ps$marker_min_pct,
  marker_log2fc_cutoff    = ps$marker_log2fc_cutoff,
  marker_padj_cutoff      = ps$marker_padj_cutoff,
  min_markers_per_cluster = ps$min_markers_per_cluster,
  min_cluster_size        = ps$min_cluster_size,
  exclude_gene_patterns   = params$de$exclude_gene_patterns,
  seed                    = params$clustering$random_seed
)
print(sweep_results)
saveRDS(sweep_results, file.path(fig_clustering_sweeps_dir, paste0("sweep_", sweep_label, ".rds")))

sweep_results_gt <- gt_generic(sweep_results, title = paste("Clustering parameter sweep:", sweep_label))
print(sweep_results_gt)
save_gt(sweep_results_gt, paste0("sweep_", sweep_label, "_full.png"), fig_clustering_sweeps_dir, vwidth = 1800)

# ------------------------------------------------------------------------------
# 3. Candidate filtering
# ------------------------------------------------------------------------------
# combinations worth a marker-support check (params$parameter_selection;
# min_ari_stability NULL = no ARI filter)
ari_ok <- if (is.null(ps$min_ari_stability)) {
  TRUE
} else {
  !is.na(sweep_results$mean_ari_stability) & sweep_results$mean_ari_stability >= ps$min_ari_stability
}
# frac_invalid_clusters is NA where fewer than 2 clusters formed
invalid_ok <- !is.na(sweep_results$frac_invalid_clusters) &
  sweep_results$frac_invalid_clusters <= ps$max_invalid_frac

candidates <- sweep_results[ari_ok & invalid_ok, ]
message(nrow(candidates), " of ", nrow(sweep_results), " combinations pass candidate filtering.")
print(candidates)

# ------------------------------------------------------------------------------
# 4. Marker support
# ------------------------------------------------------------------------------
marker_check <- NULL
if (nrow(candidates) > 0) {
  marker_check <- evaluate_all_candidates_markers(
    obj,
    candidates             = candidates,
    algorithm              = sweep_algorithm,
    marker_min_pct         = ps$marker_min_pct,
    marker_log2fc_cutoff   = ps$marker_log2fc_cutoff,
    marker_padj_cutoff     = ps$marker_padj_cutoff,
    marker_count_threshold = ps$marker_count_threshold,
    check_cols             = ps$marker_check_cols,
    exclude_gene_patterns  = params$de$exclude_gene_patterns,
    seed                   = params$clustering$random_seed
  )
  marker_check_gt <- gt_generic(marker_check$summary, title = paste("Candidates - marker support:", sweep_label))
  print(marker_check_gt)
  save_gt(marker_check_gt, paste0("sweep_", sweep_label, "_candidates_markers.png"),
          fig_clustering_sweeps_dir, vwidth = 1800)
} else {
  message("No candidates passed filtering; marker support not checked. ",
          "Loosen params$parameter_selection or widen the grid.")
}

# ------------------------------------------------------------------------------
# 5. Sweep heatmaps
# ------------------------------------------------------------------------------
# one metric per plot, dims x resolution; a third varying parameter (k or
# prune) gives one panel per value
params_varying <- function(df) c("dims", "resolution", "k", "prune")[
  vapply(c("dims", "resolution", "k", "prune"), function(p) length(unique(df[[p]])) > 1, logical(1))]

sweep_heatmaps <- function(df, metrics, file_stub, title_stub) {
  varying <- params_varying(df)
  facet <- setdiff(varying, c("dims", "resolution"))
  if (length(facet) > 1) {
    message("No heatmaps for ", file_stub, ": k and prune both vary; sweep them in separate stages.")
    return(invisible(NULL))
  }
  for (m in metrics) {
    if (!m$col %in% names(df)) next
    p <- plot_sweep_tile(df, m$col, x_col = "resolution", y_col = "dims",
                         title = paste0(m$title, " \u2014 ", title_stub),
                         subtitle = paste("algorithm =", sweep_algorithm),
                         limits = m$limits,
                         facet_col = if (length(facet) == 1) facet)
    print(p)
    save_plot(p, paste0("sweep_", sweep_label, "_", file_stub, "_", m$col, ".png"),
              fig_clustering_sweeps_dir, width = if (length(facet) == 1) 12 else 7, height = 5)
  }
}

# every heatmap in the "sweep_heatmap" palette; shares on a fixed 0-1 scale
sweep_heatmaps(sweep_results, list(
  list(col = "mean_ari_stability",    title = "Mean ARI stability",           limits = c(0, 1)),
  list(col = "n_clusters",            title = "Number of clusters"),
  list(col = "frac_invalid_clusters", title = "Fraction of invalid clusters", limits = c(0, 1)),
  list(col = "mean_silhouette",       title = "Mean silhouette width"),
  list(col = "min_cluster_size",      title = "Smallest cluster (cells)")
), "all", sweep_label)

if (!is.null(marker_check)) {
  check_metrics <- list(
    list(col = "min_real_markers",           title = "Markers of the weakest cluster"),
    list(col = "median_real_markers",        title = "Median markers per cluster"),
    list(col = "n_clusters_below_threshold", title = "Clusters with too few markers")
  )
  for (prevalence_col in grep("^max_.*_prevalence$", names(marker_check$summary), value = TRUE)) {
    check_metrics[[length(check_metrics) + 1]] <- list(
      col = prevalence_col, limits = c(0, 1),
      title = paste("Highest cluster share of", sub("^max_(.*)_prevalence$", "\\1", prevalence_col)))
  }
  sweep_heatmaps(marker_check$summary, check_metrics, "candidates", paste(sweep_label, "(candidates)"))
}

# ------------------------------------------------------------------------------
# 6. Cluster size imbalance
# ------------------------------------------------------------------------------
# descriptive, not used to filter. If the partial correlation (given the
# number of clusters) stays close to the raw one, size imbalance raises the
# ARI: a candidate with a high ARI and a high size ratio is less convincing.
imbalance <- check_cluster_size_imbalance(sweep_results)

if (!is.na(imbalance$corr_ratio_ari)) {
  imbalance_gt <- gt_generic(imbalance, title = paste("Cluster size imbalance vs. ARI stability:", sweep_label))
  print(imbalance_gt)
  save_gt(imbalance_gt, paste0("sweep_", sweep_label, "_imbalance_corr.png"), fig_clustering_sweeps_dir)

  p_imbalance <- plot_size_imbalance(sweep_results,
                                     title = paste("ARI stability vs. cluster size imbalance:", sweep_label))
  print(p_imbalance)
  save_plot(p_imbalance, paste0("sweep_", sweep_label, "_ari_vs_imbalance.png"),
            fig_clustering_sweeps_dir, width = 7, height = 5)
} else {
  message("No imbalance diagnostic for ", sweep_label,
          ": fewer than 3 rows with an ARI, or a single number of clusters.")
}

# ready-to-paste `sweep` record for this round's entry in config/decisions.R
cat("\n# paste into this round's `sweep` field in config/decisions.R:\n")
cat("sweep = list(\n")
cat("  input = \"", sweep_checkpoint, "\", algorithm = ", sweep_algorithm, ",\n", sep = "")
if (!is.null(sweep_for_round)) cat("  with_round_steps = TRUE,\n")
if (!is.null(subset_column)) {
  cat("  subset_column = \"", subset_column, "\", subset_values = c(\"",
      paste(subset_values, collapse = "\", \""), "\"),\n", sep = "")
}
cat("  gene_exclusion = ", apply_gene_exclusion, ",\n", sep = "")
cat("  dims_tested = c(", paste(dims_list, collapse = ", "), "),\n", sep = "")
cat("  resolution_tested = c(", paste(resolution_list, collapse = ", "), "),\n", sep = "")
cat("  k_tested = c(", paste(k_list, collapse = ", "), "), prune_tested = c(",
    paste(signif(prune_list, 4), collapse = ", "), ")\n", sep = "")
cat("),\n")

# ------------------------------------------------------------------------------
# 7. One choice for several rounds [optional] (edit per use)
# ------------------------------------------------------------------------------
# Rounds that should share dims/resolution (e.g. alternative exclusions):
# sweep each first (sweep_for_round, the same grid), then list them here.
# Every combination gets the worst value across the rounds: lowest ARI,
# highest fraction of invalid clusters, and whether the number of clusters
# is the same in all. Decide the rule before looking, e.g. the highest
# min_ari among combinations with same_n_clusters and max_invalid_frac <=
# params$parameter_selection$max_invalid_frac. To combine sweeps already
# saved, run section 0 and this section only.
joint_rounds <- NULL   # e.g. c("round_2", "round_2_half", "round_2_third")

if (length(joint_rounds) > 0) {
  sweep_files <- file.path(fig_clustering_sweeps_dir, paste0("sweep_", joint_rounds, ".rds"))
  missing_sweeps <- joint_rounds[!file.exists(sweep_files)]
  if (length(missing_sweeps) > 0) {
    stop("No saved sweep for: ", paste(missing_sweeps, collapse = ", "),
         ". Run this tool with sweep_for_round set to each of them first.")
  }
  joint <- combine_round_sweeps(setNames(lapply(sweep_files, readRDS), joint_rounds))
  joint_label <- paste0("joint_", paste(joint_rounds, collapse = "-"))
  print(head(joint, 10))

  joint_gt <- gt_generic(joint, title = paste("One choice for:", paste(joint_rounds, collapse = ", ")))
  save_gt(joint_gt, paste0("sweep_", joint_label, ".png"), fig_clustering_sweeps_dir, vwidth = 1800)

  for (m in list(
    list(col = "min_ari",          title = "Lowest ARI stability across rounds", limits = c(0, 1)),
    list(col = "max_invalid_frac", title = "Highest fraction of invalid clusters across rounds",
         limits = c(0, 1))
  )) {
    p <- plot_sweep_tile(joint, m$col, x_col = "resolution", y_col = "dims", title = m$title,
                         subtitle = paste(joint_rounds, collapse = ", "),
                         limits = m$limits)
    print(p)
    save_plot(p, paste0("sweep_", joint_label, "_", m$col, ".png"), fig_clustering_sweeps_dir, width = 7, height = 5)
  }
}

# ------------------------------------------------------------------------------
# Next steps (by hand):
#   1. Compare the ARI and marker-support heatmaps and the candidates table.
#   2. Optionally compare a few candidates with inspect_clustering_candidates.R
#      (it can read sweep_<label>.rds).
#   3. Record the chosen dims/resolution (or k_param, prune_snn) and the
#      `sweep` record above in the round's entry in config/decisions.R; with
#      section 7, the same dims/resolution in every round listed there.
# ------------------------------------------------------------------------------
