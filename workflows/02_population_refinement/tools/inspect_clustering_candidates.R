# ------------------------------------------------------------------------------
# inspect_clustering_candidates.R (workflows/02_population_refinement/tools)
#
# After select_clustering_parameters.R: reclusters a few chosen candidates
# and compares their per-cluster stability and top markers; can run a
# round's full characterization on a finalist and save cell lists for a
# round's exclude_cells step.
#
# Each candidate is reclustered, with per-cluster stability (repeated
# subsampling and reclustering) and FindAllMarkers(): 3-6 candidates run in
# reasonable time, more than ~10 take long.
#
# The load settings (input, subset, gene exclusion, steps) must be those of
# the sweep the candidates came from, or the clusterings will differ.
#
# Output (results/02_population_refinement/exploration/candidate_inspection/),
# if save_outputs is TRUE: stability and marker tables per candidate, and
# the full characterization of a finalist. Cell lists (section 5) go to
# data/processed/cell_lists/.
#
# Sections:
#   0. Setup
#   1. Load and preprocess
#   2. Candidates
#   3. Inspect each candidate
#   4. Full characterization [optional]
#
# Usage:
#   source(here::here("workflows", "02_population_refinement", "tools", "inspect_clustering_candidates.R"))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------
source(here::here("workflows", "02_population_refinement", "setup.R"))
decision_file <- file.path(workflow_dir, "config", "decisions.R")
if (file.exists(decision_file)) source(decision_file)

# ------------------------------------------------------------------------------
# 1. Load and preprocess (edit per use; the same values as the sweep)
# ------------------------------------------------------------------------------
inspect_checkpoint   <- "qc_filtered"
inspect_algorithm    <- 1       # 1 = Louvain, 4 = Leiden
subset_column        <- NULL    # e.g. "seurat_clusters", to inspect within one cluster
subset_values        <- NULL    # e.g. "2"
apply_gene_exclusion <- FALSE   # TRUE if the sweep used it

# [optional] for a round written in config/decisions.R (e.g. an alternative
# exclusion): uses its input, steps and exclude_genes instead of the
# settings above
inspect_for_round <- NULL   # e.g. "round_2_half"

inspect_steps <- NULL
if (!is.null(inspect_for_round)) {
  target_round <- params$population_refinement$rounds[[inspect_for_round]]
  if (is.null(target_round)) stop("'", inspect_for_round, "' is not a round in config/decisions.R.")
  inspect_checkpoint   <- target_round$input
  inspect_steps        <- target_round$steps
  apply_gene_exclusion <- isTRUE(target_round$exclude_genes)
}

obj <- prepare_clustering_input(
  inspect_checkpoint, data_processed_dir,
  subset_column         = subset_column,
  subset_values         = subset_values,
  steps                 = inspect_steps,
  exclude_gene_patterns = if (apply_gene_exclusion) params$clustering$exclude_gene_patterns,
  nfeatures             = params$clustering$nfeatures,
  vars_to_regress       = params$clustering$vars_to_regress,
  seed                  = params$clustering$random_seed,
  rounds                = params$population_refinement$rounds
)

# ------------------------------------------------------------------------------
# 2. Candidates (edit per use)
# ------------------------------------------------------------------------------
# one row per candidate; k and prune columns are optional (default
# params$clustering$k_param and $prune_snn). Note why each was picked.
candidates <- data.frame(
  dims       = c(12,  12,  15),
  resolution = c(0.5, 0.6, 0.5)
)
# or rows of a saved sweep table, e.g.:
#   sweep_results <- readRDS(file.path(fig_clustering_sweeps_dir, "sweep_qc_filtered.rds"))
#   candidates <- sweep_results[sweep_results$dims == 12 & dplyr::near(sweep_results$resolution, 0.5), ]

top_n_markers <- 15     # markers per cluster in the tables
save_outputs  <- TRUE   # FALSE: display only

if (!is.data.frame(candidates) || nrow(candidates) == 0 ||
    !all(c("dims", "resolution") %in% names(candidates))) {
  stop("`candidates` must be a non-empty data.frame with dims and resolution columns.")
}
candidates$k     <- candidates$k     %||% params$clustering$k_param
candidates$prune <- candidates$prune %||% params$clustering$prune_snn
print(candidates)

# ------------------------------------------------------------------------------
# 3. Inspect each candidate
# ------------------------------------------------------------------------------
inspection_results <- list()

for (i in seq_len(nrow(candidates))) {
  cand <- candidates[i, ]
  label <- paste0("dims", cand$dims, "_res", cand$resolution,
                  if (cand$k != params$clustering$k_param) paste0("_k", cand$k),
                  if (!isTRUE(all.equal(cand$prune, params$clustering$prune_snn))) paste0("_prune", signif(cand$prune, 3)))
  message("[inspecting] ", label)

  result <- inspect_clustering_candidate(
    obj, dims = cand$dims, resolution = cand$resolution, algorithm = inspect_algorithm,
    k_param = cand$k, prune_snn = cand$prune,
    n_reps = params$characterization$stability_n_reps, frac = params$characterization$stability_frac,
    marker_min_pct = params$de$min_pct, marker_log2fc_cutoff = params$de$log2fc_cutoff,
    marker_padj_cutoff = params$de$padj_cutoff, exclude_gene_patterns = params$de$exclude_gene_patterns,
    top_n_markers = top_n_markers, seed = params$clustering$random_seed
  )
  inspection_results[[label]] <- result

  stability_gt <- gt_generic(result$stability, title = paste(label, "— per-cluster stability (Jaccard)"))
  print(stability_gt)
  if (save_outputs) save_gt(stability_gt, paste0("inspect_", label, "_stability.png"), fig_candidate_inspection_dir)

  if (nrow(result$markers) == 0) {
    message("[", label, "] no markers passed the thresholds; no marker table.")
  } else {
    markers_gt <- gt_markers_by_cluster(result$markers, title = paste("Top cluster markers —", label))
    print(markers_gt)
    if (save_outputs) save_gt(markers_gt, paste0("inspect_", label, "_markers.png"), fig_candidate_inspection_dir)
  }
}
message("Inspected ", length(inspection_results), " candidate(s): ", paste(names(inspection_results), collapse = ", "))

# ------------------------------------------------------------------------------
# 4. Full characterization (optional, for one or two finalists)
# ------------------------------------------------------------------------------
# the characterize_clustering() sections a round gets, on the candidate's
# clustered object (not reclustered). Always saved.
full_characterization_labels <- character(0)   # e.g. c("dims12_res0.5")

qc_info <- build_qc_metric_info(params$percent_features)

for (label in full_characterization_labels) {
  candidate <- inspection_results[[label]]
  if (is.null(candidate)) {
    stop("'", label, "' is not an inspected candidate: ", paste(names(inspection_results), collapse = ", "))
  }
  candidate_obj <- run_umap(candidate$obj, dims = candidate$dims, seed = params$clustering$random_seed)

  characterize_clustering(
    candidate_obj,
    qc_metrics               = qc_info$qc_metrics,
    qc_labels                = qc_info$qc_labels,
    group_by_vars            = params$group_by_vars,
    cell_count_vars          = params$cell_count_vars,
    display_order            = params$display_order,
    stability_dims           = candidate$dims,
    stability_resolution     = candidate$resolution,
    stability_algorithm      = candidate$algorithm,
    stability_k_param        = candidate$k_param,
    stability_prune_snn      = candidate$prune_snn,
    stability_n_reps         = params$characterization$stability_n_reps,
    stability_frac           = params$characterization$stability_frac,
    marker_check_cols        = params$parameter_selection$marker_check_cols,
    violin_genes             = params$characterization$violin_genes,
    highlight_genes          = params$characterization$highlight_genes,
    de_min_pct               = params$de$min_pct,
    de_log2fc_cutoff         = params$de$log2fc_cutoff,
    de_padj_cutoff           = params$de$padj_cutoff,
    de_n_labels              = params$de$n_labels,
    de_exclude_gene_patterns = params$de$exclude_gene_patterns,
    seed                     = params$clustering$random_seed,
    output_dir               = fig_candidate_inspection_dir,
    file_prefix              = paste0("inspect_", label),
    title_label              = label
  )
}

# ------------------------------------------------------------------------------
# 5. Save cell lists (optional, edit per use)
# ------------------------------------------------------------------------------
# cells to remove in a round's exclude_cells step: list name -> the candidate
# and its cluster ID(s). Saved to data/processed/cell_lists/<name>.txt, which
# a round refers to as cells_file = "<name>.txt". Check the cluster sizes
# first (table(inspection_results[[<candidate>]]$obj$seurat_clusters)).
cell_lists_to_save <- list()   # e.g. list(low_quality_half = list(candidate = "dims8_res0.4", clusters = "1"))

for (list_name in names(cell_lists_to_save)) {
  spec <- cell_lists_to_save[[list_name]]
  candidate <- inspection_results[[spec$candidate]]
  if (is.null(candidate)) {
    stop("'", spec$candidate, "' is not an inspected candidate: ", paste(names(inspection_results), collapse = ", "))
  }
  unknown_clusters <- setdiff(spec$clusters, levels(candidate$obj$seurat_clusters))
  if (length(unknown_clusters) > 0) {
    stop("Cluster(s) ", paste(unknown_clusters, collapse = ", "), " not in ", spec$candidate, ".")
  }
  save_cell_list(Seurat::WhichCells(candidate$obj, idents = spec$clusters), list_name, data_processed_dir)
}

# ------------------------------------------------------------------------------
# Next step (by hand):
#   Compare the stability and marker tables of the candidates, and record
#   the chosen dims/resolution in the round's entry in config/decisions.R.
# ------------------------------------------------------------------------------
