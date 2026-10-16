# ------------------------------------------------------------------------------
# 01_population_refinement.R (workflows/02_population_refinement)
#
# Runs the rounds recorded in config/decisions.R. Each round:
#   load input -> steps (remove cells/genes) -> (gene exclusion) ->
#   preprocess -> (batch correction) -> cluster -> UMAP -> characterize ->
#   DEG overlap -> exclusion check -> checkpoint
# See README.md ("How rounds work"), including how to compare alternative
# exclusions as parallel rounds.
#
# A round's input is a checkpoint (e.g. the QC workflow's) or an earlier
# round. A computed round is skipped; if its settings, its input file or
# the labels of the round it builds on changed since, the script stops
# instead of reusing it (force_rerun recomputes).
#
# Input:
#   data/processed/01_qc_filtered.rds (from 01_qc)
#   config/decisions.R (the rounds; see decisions_example.R)
#
# Output:
#   results/02_population_refinement/<round>_<label>/
#     - characterization of the round's clustering, titled with the round's
#       label
#     - <prefix>_de_results.csv and _de_one_vs_rest.xlsx: one-vs-rest DE
#     - DEG overlap Venns and workbooks [optional]
#     - exclusion_check/: excluded vs retained DE, GO/KEGG and workbook,
#       cells retained per group (rounds that remove cells) [optional]
#     - pairwise_de/: DE between every pair of clusters (rounds with
#       `pairwise_de`)
#     - exclusion_steps.png, round_summary.png/.csv
#   data/processed/
#     - 02_population_refinement_<round>.rds (+ _manifest.txt), and
#       _config.rds: the settings used to detect a changed round
#     - 02_population_refined.rds (+ _manifest.txt): the final round, with a
#       cluster_label column; input to the downstream workflows
#
# Human decision:
#   Per round: dims, resolution and algorithm (tools/), exclusion steps and
#   cluster labels, recorded in config/decisions.R; the final round.
#
# Sections:
#   0. Setup
#   1. Check the rounds
#   2. Run the rounds
#   3. Pairwise cluster DE
#   4. Final checkpoint
#
# Usage:
#   source(here::here("workflows", "02_population_refinement", "pipeline", "01_population_refinement.R"))
#
#   Prerequisites: 01_qc run; config/decisions.R copied from
#   decisions_example.R with at least one round.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------
source(here::here("workflows", "02_population_refinement", "setup.R"))

decision_file <- file.path(workflow_dir, "config", "decisions.R")
if (!file.exists(decision_file)) {
  stop("No workflows/02_population_refinement/config/decisions.R found. Choose round 1's ",
       "parameters with tools/select_clustering_parameters.R, then copy ",
       "config/decisions_example.R to config/decisions.R and record the round.")
}
source(decision_file)

# ------------------------------------------------------------------------------
# 1. Check the rounds
# ------------------------------------------------------------------------------
rounds <- params$population_refinement$rounds
if (length(rounds) == 0 || is.null(names(rounds)) || any(!nzchar(names(rounds)))) {
  stop("params$population_refinement$rounds is empty or has unnamed rounds. Check config/decisions.R ",
       "(see decisions_example.R).")
}
round_names <- names(rounds)

round_fields <- c("label", "rationale", "input", "algorithm", "dims", "resolution", "k_param", "prune_snn",
                  "exclude_genes", "batch_correction", "steps", "cluster_labels", "pairwise_de")
for (i in seq_along(round_names)) {
  missing_fields <- setdiff(c("input", "algorithm", "dims", "resolution"), names(rounds[[i]]))
  if (length(missing_fields) > 0) {
    stop("Round '", round_names[i], "' is missing: ", paste(missing_fields, collapse = ", "), ".")
  }
  # fields ending in "sweep" are documentation
  unknown_fields <- setdiff(grep("sweep$", names(rounds[[i]]), value = TRUE, invert = TRUE), round_fields)
  if (length(unknown_fields) > 0) {
    stop("Round '", round_names[i], "' has unknown field(s): ", paste(unknown_fields, collapse = ", "),
         ". Valid fields: ", paste(round_fields, collapse = ", "), " (see decisions_example.R).")
  }
  # a round can only build on a round listed before it
  source_round <- input_round_name(rounds[[i]]$input, round_names)
  if (!is.null(source_round) && match(source_round, round_names) >= i) {
    stop("Round '", round_names[i], "' has input '", rounds[[i]]$input, "', which is not an earlier round.")
  }
}

final_round <- params$population_refinement$final_round
if (!is.null(final_round) && !final_round %in% round_names) {
  stop("params$population_refinement$final_round '", final_round, "' is not a round: ",
       paste(round_names, collapse = ", "))
}
if (any(vapply(rounds, function(r) !is.null(r$pairwise_de), logical(1)))) {
  load_libraries("EnhancedVolcano", load_core = FALSE)
}

force_rerun <- isTRUE(params$population_refinement$force_rerun)

# the settings written to a round's manifest
round_settings <- function(round) {
  round[intersect(c("label", "input", "algorithm", "dims", "resolution", "k_param", "prune_snn",
                    "exclude_genes", "batch_correction", "steps", "preprocessing"), names(round))]
}
qc_info <- build_qc_metric_info(params$percent_features)

# ------------------------------------------------------------------------------
# 2. Run the rounds
# ------------------------------------------------------------------------------
for (round_name in round_names) {

  round <- rounds[[round_name]]
  has_steps <- length(round$steps) > 0

  # --- complete the round's record --------------------------------------------------
  # everything that changes the result is recorded with the checkpoint, so a
  # later change is detected: the clustering defaults, the preprocessing
  # settings, the input file's checksum, a cells_file's checksum and the
  # labels of the input round (steps may select on them)
  round$k_param   <- round$k_param   %||% params$clustering$k_param
  round$prune_snn <- round$prune_snn %||% params$clustering$prune_snn
  round$preprocessing <- list(
    vars_to_regress       = params$clustering$vars_to_regress,
    nfeatures             = params$clustering$nfeatures,
    random_seed           = params$clustering$random_seed,
    exclude_gene_patterns = if (isTRUE(round$exclude_genes)) params$clustering$exclude_gene_patterns
  )
  for (j in seq_along(round$steps)) {
    if (!is.null(round$steps[[j]]$cells_file)) {
      cells_path <- resolve_cells_file(round$steps[[j]]$cells_file, data_processed_dir)
      round$steps[[j]]$cells_md5 <- unname(tools::md5sum(cells_path))
    }
  }
  input_path <- refinement_input_path(round$input, data_processed_dir)
  round$input_md5 <- unname(tools::md5sum(input_path))
  source_round <- input_round_name(round$input, round_names)
  round$input_labels <- if (!is.null(source_round)) rounds[[source_round]]$cluster_labels %||% character(0)

  out_path    <- refinement_checkpoint_path(round_name, data_processed_dir)
  config_path <- refinement_config_path(round_name, data_processed_dir)

  # --- skip a computed round -------------------------------------------------------------
  if (file.exists(out_path) && !force_rerun) {
    unchanged <- round_config_unchanged(round, round_name, data_processed_dir)
    if (!unchanged) {
      stop("Round '", round_name, "' changed since its checkpoint was made (",
           paste(attr(unchanged, "changed_fields"), collapse = ", "), "). Its checkpoint, and every ",
           "round built on it, would mix old and new settings. Set ",
           "params$population_refinement$force_rerun <- TRUE in params_override.R to recompute, ",
           "or revert the change. (input_md5 changes when the input file is rewritten, e.g. by re-running QC.)")
    }
    message("[skip] ", round_name, " already computed with the same settings and input.")
    next
  }

  message("[run] ", round_name, if (!is.null(round$label)) paste0(": ", round$label))
  if (!is.null(round$rationale)) message("  rationale: ", round$rationale)

  round_dir <- round_results_dir(round_name, round, fig_refinement_dir)
  dir.create(round_dir, recursive = TRUE, showWarnings = FALSE)

  # titles use the round's label, file names the round's name
  round_title <- round$label %||% round_name
  title_label <- round_title
  file_prefix <- round_name

  # --- load the input -------------------------------------------------------------------
  data_obj <- readRDS(input_path)
  if (!is.null(source_round)) data_obj <- add_cluster_label_column(data_obj, round$input_labels)
  n_cells_before <- ncol(data_obj)
  n_genes_before <- nrow(data_obj)

  # a malformed step stops the round before preprocessing and clustering
  validate_refinement_steps(data_obj, rounds[[round_name]]$steps, data_processed_dir)

  # the input as loaded, all genes: the baseline of the exclusion check
  input_obj <- if (isTRUE(steps$exclusion_check) && has_steps) data_obj else NULL

  # --- steps: on the input, with its own metadata ------------------------------------------
  if (has_steps) {
    step_result <- apply_refinement_steps(data_obj, round$steps, dir = data_processed_dir)
    data_obj <- step_result$obj
    message("  removed ", length(step_result$removed_cells), " cells and ",
            length(step_result$removed_genes), " genes before clustering")
    audit_gt <- gt_generic(step_result$audit, title = paste("Exclusion steps:", round_title))
    print(audit_gt)
    save_gt(audit_gt, "exclusion_steps.png", round_dir)
  }

  # --- gene exclusion (exclude_genes = TRUE) ----------------------------------------------
  if (isTRUE(round$exclude_genes)) {
    data_obj <- exclude_genes(data_obj, genes = paste(params$clustering$exclude_gene_patterns, collapse = "|"),
                              pattern = TRUE)$obj
  }

  # --- preprocess, from scratch: earlier reductions and graphs no longer apply ---------------
  data_obj <- preprocess_seurat(data_obj, nfeatures = params$clustering$nfeatures,
                                vars_to_regress = params$clustering$vars_to_regress,
                                show_elbow = FALSE, seed = params$clustering$random_seed)$obj

  # --- batch correction -----------------------------------------------------------------------
  cluster_reduction <- "pca"
  if (!is.null(round$batch_correction)) {
    correction <- apply_batch_correction(data_obj, method = round$batch_correction$method,
                                         args = round$batch_correction, seed = params$clustering$random_seed)
    data_obj <- correction$obj
    cluster_reduction <- correction$reduction
  }

  # --- cluster and UMAP ---------------------------------------------------------------------------
  data_obj <- cluster_seurat(data_obj, resolution = round$resolution, algorithm = round$algorithm,
                             reduction = cluster_reduction, dims = round$dims, k_param = round$k_param,
                             prune_snn = round$prune_snn, seed = params$clustering$random_seed)
  data_obj <- run_umap(data_obj, dims = round$dims, reduction = cluster_reduction,
                       seed = params$clustering$random_seed)
  data_obj <- add_cluster_label_column(data_obj, round$cluster_labels)

  # --- characterize ---------------------------------------------------------------------------------
  characterization <- characterize_clustering(
    data_obj,
    qc_metrics                  = qc_info$qc_metrics,
    qc_labels                   = qc_info$qc_labels,
    group_by_vars               = params$group_by_vars,
    cell_count_vars             = params$cell_count_vars,
    display_order               = params$display_order,
    cluster_labels              = round$cluster_labels,
    cluster_color_map           = params$cluster_color_map,
    stability_dims              = round$dims,
    stability_resolution        = round$resolution,
    stability_algorithm         = round$algorithm,
    stability_k_param           = round$k_param,
    stability_prune_snn         = round$prune_snn,
    stability_reduction         = cluster_reduction,
    stability_n_reps            = params$characterization$stability_n_reps,
    stability_frac              = params$characterization$stability_frac,
    marker_check_cols           = params$parameter_selection$marker_check_cols,
    violin_genes                = params$characterization$violin_genes,
    highlight_genes             = params$characterization$highlight_genes,
    de_min_pct                  = params$de$min_pct,
    de_log2fc_cutoff            = params$de$log2fc_cutoff,
    de_padj_cutoff              = params$de$padj_cutoff,
    de_n_labels                 = params$de$n_labels,
    de_exclude_gene_patterns    = params$de$exclude_gene_patterns,
    run_composition_stability   = isTRUE(steps$composition_stability),
    run_qc_metrics              = isTRUE(steps$qc_metrics),
    run_group_umaps             = isTRUE(steps$group_umaps),
    run_cell_counts             = isTRUE(steps$cell_counts),
    run_marker_characterization = isTRUE(steps$marker_characterization),
    run_differential_expression = isTRUE(steps$differential_expression),
    seed                        = params$clustering$random_seed,
    output_dir                  = round_dir,
    file_prefix                 = file_prefix,
    title_label                 = title_label
  )

  # every tested gene (the PNG tables show the top genes); the .csv is also
  # read by tools/compare_rounds.R
  de_results <- characterization$de_results
  if (!is.null(de_results) && nrow(de_results) > 0) {
    de_results$cluster_label <- as.character(data_obj$cluster_label[match(
      as.character(de_results$cluster), as.character(data_obj$seurat_clusters))])
    utils::write.csv(de_results, file.path(round_dir, paste0(file_prefix, "_de_results.csv")), row.names = FALSE)
    write_de_workbook(de_results, paste0(file_prefix, "_de_one_vs_rest.xlsx"), round_dir)
  }

  # --- DEG overlap between clusters ------------------------------------------------------------
  # genes up in each cluster vs the rest, then genes down, so a shared gene
  # changes the same way in both clusters
  if (isTRUE(steps$deg_overlap)) {
    if (is.null(de_results) || nrow(de_results) == 0) {
      message("[skip] ", round_name, " / deg_overlap (needs differential_expression results)")
    } else {
      cluster_ids <- levels(data_obj$seurat_clusters)
      id_to_label <- setNames(as.character(data_obj$cluster_label[match(cluster_ids, data_obj$seurat_clusters)]),
                              cluster_ids)
      chosen <- params$deg_overlap$clusters
      chosen_ids <- if (is.null(chosen)) cluster_ids else cluster_ids[cluster_ids %in% chosen | id_to_label %in% chosen]

      if (length(chosen_ids) < 2 || length(chosen_ids) > 4) {
        message("[skip] ", round_name, " / deg_overlap (", length(chosen_ids), " clusters selected; ",
                "a Venn needs 2-4: set params$deg_overlap$clusters)")
      } else {
        de_by_cluster <- split(de_results, as.character(de_results$cluster))
        de_list <- lapply(chosen_ids, function(id) de_by_cluster[[id]] %||% de_results[0, ])
        names(de_list) <- make.unique(unname(id_to_label[chosen_ids]))
        for (direction in c("up", "down")) {
          compare_deg_venn(de_list, log2fc_cutoff = params$deg_overlap$log2fc_cutoff,
                           padj_cutoff = params$deg_overlap$padj_cutoff, direction = direction,
                           title_label = title_label, output_dir = round_dir,
                           file_prefix = paste0(file_prefix, "_deg_overlap"))
        }
      }
    }
  }

  # --- exclusion check: what was removed, and from which groups -------------------------------------
  # excluded vs retained cells of the input as loaded (all genes)
  if (!is.null(input_obj) && ncol(data_obj) < ncol(input_obj)) {
    check_dir <- file.path(round_dir, "exclusion_check")
    check_prefix <- paste0(round_name, "_excluded_vs_retained")
    input_obj$exclusion_status <- ifelse(colnames(input_obj) %in% colnames(data_obj), "retained", "excluded")
    ec <- params$exclusion_check

    exclusion_de <- run_de_enrichment(
      input_obj, group_col = "exclusion_status", ident_1 = "excluded", ident_2 = "retained",
      min_pct = ec$min_pct, log2fc_cutoff = ec$log2fc_cutoff, padj_cutoff = ec$padj_cutoff,
      exclude_gene_patterns = params$de$exclude_gene_patterns, n_terms = ec$n_terms,
      label_1 = "Excluded", label_2 = "Retained", output_dir = check_dir, file_prefix = check_prefix
    )
    if (nrow(exclusion_de$de_tested) > 0) {
      volcano <- plot_cluster_volcanoes(
        data.frame(gene = rownames(exclusion_de$de_tested), exclusion_de$de_tested),
        subtitle = paste0(round_title, ": excluded (", sum(input_obj$exclusion_status == "excluded"),
                          ") vs retained (", sum(input_obj$exclusion_status == "retained"), ") cells"),
        log2fc_cutoff = ec$log2fc_cutoff, padj_cutoff = ec$padj_cutoff, n_labels = params$de$n_labels
      )[["comparison"]]
      if (!is.null(volcano)) save_plot(volcano, paste0(check_prefix, "_volcano.png"), check_dir, width = 10, height = 8)
    }

    plot_retention(input_obj, data_obj, params$cell_count_vars, params$display_order,
                   title_prefix = round_title, output_dir = check_dir, file_prefix = round_name)
  }
  rm(input_obj)

  n_cells_after <- ncol(data_obj)
  n_genes_after <- nrow(data_obj)
  message(round_name, ": ", n_cells_before, " -> ", n_cells_after, " cells, ",
          n_genes_before, " -> ", n_genes_after, " genes")

  round_summary <- data.frame(
    round = round_name, label = round$label %||% "", input = round$input,
    algorithm = round$algorithm, dims = round$dims, resolution = round$resolution,
    k_param = round$k_param, prune_snn = signif(round$prune_snn, 3),
    n_clusters = nlevels(data_obj$seurat_clusters),
    n_cells_input = n_cells_before, n_cells_final = n_cells_after,
    pct_cells_kept = round(100 * n_cells_after / n_cells_before, 1),
    n_genes_input = n_genes_before, n_genes_final = n_genes_after
  )
  utils::write.csv(round_summary, file.path(round_dir, "round_summary.csv"), row.names = FALSE)
  summary_gt <- gt_generic(round_summary, title = paste("Round summary:", round_title),
                           subtitle = if (has_steps) "after this round's exclusion steps" else "no exclusion steps")
  print(summary_gt)
  save_gt(summary_gt, "round_summary.png", round_dir, vwidth = 1400)

  # --- checkpoint and its recorded settings ---------------------------------------------------------------
  write_checkpoint(data_obj, file = basename(out_path), dir = data_processed_dir,
                   workflow = "02_population_refinement", script = "01_population_refinement.R",
                   settings = round_settings(round), input_path = input_path)
  saveRDS(round, config_path)
}

# ------------------------------------------------------------------------------
# 3. Pairwise cluster DE
# ------------------------------------------------------------------------------
# For rounds with `pairwise_de`. Runs from the saved checkpoints, so it can
# be added to a computed round; skipped when the workbook exists, unless
# force_rerun.
for (round_name in round_names[vapply(rounds, function(r) !is.null(r$pairwise_de), logical(1))]) {
  round <- rounds[[round_name]]
  pairwise_dir <- file.path(round_results_dir(round_name, round, fig_refinement_dir), "pairwise_de")
  workbook <- paste0(round_name, "_pairwise_de.xlsx")
  if (file.exists(file.path(pairwise_dir, workbook)) && !force_rerun) {
    message("[skip] ", round_name, " / pairwise DE (already done)")
    next
  }

  obj <- load_refinement_input(round_name, data_processed_dir)
  obj <- add_cluster_label_column(obj, round$cluster_labels)
  labels <- setNames(as.character(obj$cluster_label), as.character(obj$seurat_clusters))
  labels <- labels[!duplicated(names(labels))]

  exclude <- if (is.list(round$pairwise_de)) round$pairwise_de$exclude else character(0)
  unknown <- setdiff(exclude, labels)
  if (length(unknown) > 0) {
    stop("Round '", round_name, "': pairwise_de$exclude has label(s) not in this round: ",
         paste(unknown, collapse = ", "), ".")
  }
  message("[run] ", round_name, " / pairwise DE",
          if (length(exclude) > 0) paste0(" (without ", paste(exclude, collapse = ", "), ")"))

  pairwise <- pairwise_cluster_de(
    obj, cluster_labels = labels, exclude_labels = exclude,
    min_pct = params$de$min_pct, log2fc_cutoff = params$de$log2fc_cutoff,
    padj_cutoff = params$de$padj_cutoff, n_labels = params$de$n_labels,
    exclude_gene_patterns = params$de$exclude_gene_patterns,
    title_label = round$label %||% round_name, output_dir = pairwise_dir, file_prefix = round_name
  )

  # workbook: a summary sheet (DEGs per pair), then one sheet per pair
  summary_rows <- lapply(names(pairwise), function(pair) {
    de <- pairwise[[pair]]$de_results
    deg <- de[!is.na(de$p_val_adj) & de$p_val_adj < params$de$padj_cutoff &
                abs(de$avg_log2FC) > params$de$log2fc_cutoff, , drop = FALSE]
    parts <- strsplit(pair, " vs ", fixed = TRUE)[[1]]
    data.frame(comparison = pair, first = parts[1], second = parts[2], n_DE = nrow(deg),
               higher_in_first = sum(deg$avg_log2FC > 0), higher_in_second = sum(deg$avg_log2FC < 0))
  })
  pair_sheets <- lapply(pairwise, function(x) data.frame(gene = rownames(x$de_results), x$de_results, row.names = NULL))
  # unique within Excel's 31-character limit
  names(pair_sheets) <- make.unique(substr(names(pair_sheets), 1, 28), sep = "_")
  save_xlsx(c(list(summary = do.call(rbind, summary_rows)), pair_sheets), workbook, pairwise_dir)
  rm(obj, pairwise)
}

# ------------------------------------------------------------------------------
# 4. Final checkpoint
# ------------------------------------------------------------------------------
# The final round under a fixed name, so downstream workflows need only this
# file and its manifest. Labels come from config/decisions.R, so a label
# edited after the round ran is included.
if (!is.null(final_round)) {
  final_path <- refinement_checkpoint_path(final_round, data_processed_dir)
  final_obj <- add_cluster_label_column(readRDS(final_path), rounds[[final_round]]$cluster_labels)

  write_checkpoint(
    final_obj, "population_refined", data_processed_dir,
    workflow = "02_population_refinement", script = "01_population_refinement.R",
    settings = list(
      final_round    = final_round,
      rounds         = paste(round_names, collapse = ", "),
      cluster_labels = rounds[[final_round]]$cluster_labels,
      algorithm      = rounds[[final_round]]$algorithm,
      dims           = rounds[[final_round]]$dims,
      resolution     = rounds[[final_round]]$resolution,
      k_param        = rounds[[final_round]]$k_param %||% params$clustering$k_param,
      prune_snn      = rounds[[final_round]]$prune_snn %||% params$clustering$prune_snn
    ),
    input_path = final_path
  )
  message("Final clustering (", final_round, "): ", ncol(final_obj), " cells, ",
          nlevels(final_obj$seurat_clusters), " clusters")
} else {
  message("[info] params$population_refinement$final_round is not set; no final checkpoint written.")
}
