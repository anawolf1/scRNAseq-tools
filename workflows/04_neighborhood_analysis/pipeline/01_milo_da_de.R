# ------------------------------------------------------------------------------
# 01_milo_da_de.R (workflows/04_neighborhood_analysis)
#
# Differential abundance (miloR) and neighborhood differential expression
# (miloDE). Runs the analyses in config/decisions.R. Each analysis builds one
# set of neighborhoods on the cells of its timepoints and tests each of its
# contrasts (c(test, reference)) on those shared neighborhoods:
#   - primary: the main neighborhood landscape (e.g. BL + 3M + FU + SPR),
#     fully characterized: diagnostics, neighborhood groups,
#     meta-neighborhoods, markers, miloDE
#   - sensitivity: the same contrasts on another landscape (e.g. without
#     BL), DA results only
#   - pairwise: each primary contrast on a landscape built from its two
#     timepoints only
# Analyses are compared per cell (neighborhood IDs differ between them).
#
# Input:
#   data/processed/02_population_refined.rds (from 02_population_refinement)
#   data/processed/03_functional_scores.rds (from 03_functional_scores; steps$nhood_program_scores)
#   config/decisions.R (the analyses and parameters; see decisions_example.R)
#
# Output:
#   results/04_neighborhood_analysis/01_milo_da/<analysis>/
#     - neighborhood size histogram, sample table
#     - <contrast>/ DA graph, beeswarm, p-value histogram, per-cluster summary,
#       results workbook; primary only: variability diagnostic, per-sample
#       abundance of the most extreme neighborhoods, cells of significant
#       neighborhoods on UMAP, compositions, programs and QC per
#       neighborhood, meta-neighborhoods (markers/)
#     - primary only: nhood_groups/ (groups of ALL neighborhoods, Leiden on
#       shared cells: program heatmaps, composition, maps, group- vs
#       cluster-level abundance tests, one-vs-rest DE per group in markers/),
#       meta_nhoods/ (significant neighborhoods merged; descriptive, kept for
#       comparison), chosen_nhoods/
#   results/04_neighborhood_analysis/01_milo_da/robustness/
#     - per contrast: cell-level DA calls compared between analyses
#   results/04_neighborhood_analysis/02_milo_de/<contrast>/ (primary contrasts)
#   data/processed/04_milo_da.rds (primary), 04_milo_da_<analysis>.rds,
#     04_milo_de_<contrast>.rds (each with a _manifest.txt)
#
# Statistical unit:
#   Samples (patient x timepoint). Each test uses only its two timepoints'
#   samples. Contrasts between two repeated-measures timepoints (e.g. 3M vs
#   FU) use a paired design (~ patient + condition, patients sampled at both);
#   contrasts against a separately sampled group are unpaired (~ condition).
#
# Sections:
#   0. Setup
#   1. Load the object
#   2. miloR: each analysis (neighborhoods, tests, diagnostics)   [milo_da]
#   3. Robustness: cell-level calls across analyses               [milo_da]
#   4. miloDE: each primary contrast                              [milo_de]
#
# Usage:
#   source(here::here("workflows", "04_neighborhood_analysis", "pipeline", "01_milo_da_de.R"))
#
#   Prerequisites: 02_population_refinement run (and 03_functional_scores
#   for the program scores); parameters chosen with
#   tools/select_milo_parameters.R and recorded in config/decisions.R.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------
source(here::here("workflows", "04_neighborhood_analysis", "setup.R"))

decision_file <- file.path(workflow_dir, "config", "decisions.R")
if (!file.exists(decision_file)) {
  stop("No workflows/04_neighborhood_analysis/config/decisions.R found. Copy config/decisions_example.R ",
       "to config/decisions.R, set the analyses, and choose k/d with tools/select_milo_parameters.R.")
}
source(decision_file)

patient_col   <- params$study_design$patient_col
timepoint_col <- params$study_design$timepoint_col
repeated_timepoints <- params$study_design$repeated_measures_timepoints

run_da      <- isTRUE(steps$milo_da)
run_markers <- isTRUE(steps$da_nhood_markers)
run_scores  <- isTRUE(steps$nhood_program_scores)
run_de      <- isTRUE(steps$milo_de)

# analyses: validated and expanded (pairwise -> one analysis per primary contrast)
analyses <- resolve_milo_analyses(params$analyses)
primary_name <- names(analyses)[vapply(analyses, function(a) a$role == "primary", logical(1))]
primary <- analyses[[primary_name]]
contrasts <- primary$contrasts

if (run_da && (is.null(params$da$k) || is.null(params$da$d))) {
  stop("params$da$k and params$da$d are not set in config/decisions.R - choose them with ",
       "tools/select_milo_parameters.R (Part A).")
}
de_nhoods <- params$de$neighborhoods
if (!de_nhoods %in% c("da", "milode")) {
  stop("params$de$neighborhoods must be \"da\" or \"milode\".")
}
if (run_de && de_nhoods == "da" && !run_da) {
  stop("params$de$neighborhoods = \"da\" reuses the DA neighborhoods, so steps$milo_da must be on.")
}
de_d <- params$de$d %||% params$da$d
if (run_de && de_nhoods == "milode" && (is.null(params$de$k) || is.null(de_d))) {
  stop("params$de$k (and params$de$d or params$da$d) are not set in config/decisions.R - choose them ",
       "with tools/select_milo_parameters.R (Part B).")
}

da_alpha <- params$da$spatial_fdr
de_alpha <- params$de$alpha
lfc_pal  <- resolve_palette("de_log2fc")
marker_exclude_patterns <- params$da$markers$exclude_gene_patterns
qc_cols_nhood <- params$da$nhood_qc_cols

# ------------------------------------------------------------------------------
# 1. Load the object
# ------------------------------------------------------------------------------
# every timepoint any analysis uses; each analysis then keeps its own
all_timepoints <- unique(unlist(lapply(analyses, `[[`, "timepoints")))
all_contrasts  <- unique(unlist(lapply(analyses, `[[`, "contrasts"), recursive = FALSE))
names(all_contrasts) <- vapply(all_contrasts, paste, character(1), collapse = "_vs_")
prepared <- prepare_milo_input(input_path, all_timepoints, all_contrasts, patient_col, timepoint_col)
full_obj <- prepared$obj
message("Input: ", basename(input_path), " (", ncol(full_obj), " cells of ", paste(all_timepoints, collapse = ", "), ")")

timepoint_levels <- intersect(c(repeated_timepoints, params$study_design$cross_sectional_timepoints), all_timepoints)
if (length(timepoint_levels) == 0) timepoint_levels <- all_timepoints
timepoint_colors <- resolve_group_colors(timepoint_col, timepoint_levels)

cluster_colors <- resolve_cluster_colors(levels(full_obj$cluster_label), params$cluster_color_map)

# functional scores, if requested and available
score_values <- NULL
score_cols   <- character(0)
if (run_scores) {
  score_path <- checkpoint_path("functional_scores", data_processed_dir)
  if (file.exists(score_path)) {
    check_checkpoint_input(score_path, input_path,
                           rerun = "Re-run 03_functional_scores, or switch off steps$nhood_program_scores.")
    scored <- readRDS(score_path)
    score_cols <- grep("Score$", colnames(scored[[]]), value = TRUE)
    score_values <- scored[[]][, score_cols, drop = FALSE]
    n_shared <- length(intersect(colnames(full_obj), rownames(score_values)))
    message("Functional scores: ", n_shared, " of ", ncol(full_obj), " cells found in ", basename(score_path))
    if (n_shared == 0) stop("No cells shared with ", score_path, " - were both built from the same input?")
    rm(scored)
  } else {
    message("[skip] nhood_program_scores (no functional-scores checkpoint at ", score_path,
            "; run 03_functional_scores first)")
    run_scores <- FALSE
  }
}
qc_cols_nhood <- intersect(qc_cols_nhood, colnames(full_obj[[]]))

# ------------------------------------------------------------------------------
# 2. miloR: each analysis
# ------------------------------------------------------------------------------
cell_calls <- list()     # analysis -> contrast -> named vector (cell -> DA call)
da_results_all <- list() # primary analysis: contrast -> results
milo_da <- NULL

for (analysis_name in if (run_da) names(analyses) else character(0)) {
  analysis <- analyses[[analysis_name]]
  is_primary <- analysis$role == "primary"
  analysis_dir <- file.path(fig_milo_da_dir, analysis_name)
  dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
  message("\n[analysis] ", analysis_name, " (", analysis$role, "): graph on ",
          paste(analysis$timepoints, collapse = " + "), "; contrasts ", paste(names(analysis$contrasts), collapse = ", "))

  data_obj <- subset(full_obj, cells = colnames(full_obj)[as.character(full_obj[[timepoint_col, drop = TRUE]]) %in% analysis$timepoints])
  analysis_timepoints <- intersect(timepoint_levels, analysis$timepoints)
  p_umap_timepoint <- plot_umap_by_group(data_obj, timepoint_col, group_order = analysis_timepoints,
                                         title = paste("Cells by", timepoint_col))$combined

  sample_table <- prepared$sample_table[prepared$sample_table[[timepoint_col]] %in% analysis$timepoints, ]
  sample_gt <- gt_generic(sample_table, title = paste("Samples and cells per timepoint \u2014", analysis_name))
  print(sample_gt)
  save_gt(sample_gt, "sample_table.png", analysis_dir)

  # --- neighborhoods + counts ---------------------------------------------------------
  refinement_scheme <- params$da$refinement_scheme
  k_used <- analysis$k %||% params$da$k
  d_used <- analysis$d %||% params$da$d
  milo_obj <- seurat_to_milo(data_obj)
  milo_obj <- build_milo_nhoods(milo_obj, k = k_used, d = d_used, prop = params$da$prop,
                                refinement_scheme = refinement_scheme, seed = params$da$seed)
  milo_obj <- miloR::countCells(milo_obj, samples = "milo_sample",
                                meta.data = as.data.frame(SummarizedExperiment::colData(milo_obj)))
  milo_obj <- miloR::buildNhoodGraph(milo_obj)
  da_membership <- miloR::nhoods(milo_obj)
  da_meta <- as.data.frame(SummarizedExperiment::colData(milo_obj))
  message("  ", ncol(da_membership), " neighborhoods (k = ", k_used, ", d = ", d_used, ")")

  composition <- nhood_composition(
    da_membership, SummarizedExperiment::colData(milo_obj),
    sample_col = "milo_sample", patient_col = patient_col, group_col = timepoint_col,
    cluster_col = "cluster_label", contrasts = analysis$contrasts
  )
  p_nhood_size <- plot_nhood_size_distributions(
    data.frame(combo = paste0("k", k_used, "_d", d_used), size = composition$size),
    target = params$da$size_multiplier * dplyr::n_distinct(data_obj$milo_sample),
    title = paste("Neighborhood size \u2014", analysis_name),
    subtitle = paste0("Dashed line: ", params$da$size_multiplier, " x all samples")
  )
  print(p_nhood_size)
  save_plot(p_nhood_size, "nhood_size_histogram.png", analysis_dir, width = 7, height = 5)

  sample_design <- unique(da_meta[, c("milo_sample", patient_col, timepoint_col)])
  rownames(sample_design) <- sample_design$milo_sample

  nhood_scores <- if (is_primary && run_scores) nhood_mean_scores(da_membership, score_values) else NULL
  nhood_qc <- if (is_primary && length(qc_cols_nhood) > 0) {
    nhood_mean_scores(da_membership, da_meta[, qc_cols_nhood, drop = FALSE])[, c("Nhood", qc_cols_nhood)]
  } else NULL

  analysis_results <- list()
  analysis_groups  <- list()   # primary: meta-neighborhoods defined in each contrast
  cell_calls[[analysis_name]] <- list()

  #-------------------------------------------------------------------------------
  # 2a. Test each contrast
  #-------------------------------------------------------------------------------
  for (contrast_name in names(analysis$contrasts)) {
    levels_used <- analysis$contrasts[[contrast_name]]   # c(test, reference)
    out_dir <- file.path(analysis_dir, contrast_name)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("  [DA] ", contrast_name)

    design <- contrast_design(sample_design, levels_used, patient_col, timepoint_col,
                              repeated_timepoints, paired = params$da$paired_design)
    design_df <- design$design_df
    if (length(design$dropped) > 0) {
      message("    paired design: left out ", length(design$dropped), " sample(s) of patients not sampled at both: ",
              paste(design$dropped, collapse = ", "))
    }

    da_res <- miloR::testNhoods(milo_obj, design = design$formula, design.df = design_df,
                                fdr.weighting = milo_fdr_weighting(refinement_scheme),
                                norm.method = params$da$norm_method, reduced.dim = "PCA")
    da_res <- annotate_da_results(milo_obj, da_res, "cluster_label",
                                  mixed_threshold = params$da$mixed_threshold)
    da_res <- dplyr::left_join(
      da_res,
      composition[, c("Nhood", "size", "patient_purity", paste0(contrast_name, c("_cells", "_n_test", "_n_ref")))],
      by = "Nhood"
    )
    abundance <- nhood_sample_abundance(milo_obj, design_df, patient_col)
    da_res <- dplyr::left_join(da_res, nhood_variability(abundance), by = "Nhood")
    if (!is.null(nhood_scores)) da_res <- dplyr::left_join(da_res, nhood_scores, by = "Nhood")
    if (!is.null(nhood_qc)) da_res <- dplyr::left_join(da_res, nhood_qc, by = "Nhood")
    da_res <- da_res[order(da_res$SpatialFDR), ]
    analysis_results[[contrast_name]] <- da_res
    cell_calls[[analysis_name]][[contrast_name]] <- cell_da_calls(da_membership, da_res, da_alpha)

    n_sig <- sum(da_res$SpatialFDR < da_alpha, na.rm = TRUE)
    message("    ", n_sig, " of ", nrow(da_res), " neighborhoods with SpatialFDR < ", da_alpha,
            " (", sum(da_res$SpatialFDR < da_alpha & da_res$logFC > 0, na.rm = TRUE), " enriched in ",
            levels_used[1], ", ", sum(da_res$SpatialFDR < da_alpha & da_res$logFC < 0, na.rm = TRUE),
            " enriched in ", levels_used[2], ")", if (design$paired) "; paired design" else "")

    subtitle_da <- paste0(nrow(design_df), " samples (", sum(design_df$condition == levels_used[1]), " ",
                          levels_used[1], ", ", sum(design_df$condition == levels_used[2]), " ",
                          levels_used[2], if (design$paired) ", paired by patient" else "",
                          "); graph: ", paste(analysis$timepoints, collapse = " + "),
                          "; colored if SpatialFDR < ", da_alpha,
                          "; logFC > 0 = more abundant in ", levels_used[1])

    # --- neighborhood graph ---------------------------------------------------
    p_graph <- miloR::plotNhoodGraphDA(milo_obj, da_res, alpha = da_alpha, layout = "UMAP") +
      ggplot2::scale_fill_gradient2(low = lfc_pal$low, mid = lfc_pal$mid, high = lfc_pal$high,
                                    midpoint = 0, name = paste0("logFC\n", levels_used[1], " vs ", levels_used[2])) +
      ggplot2::labs(title = paste0("Differential abundance \u2014 ", levels_used[1], " vs ", levels_used[2]),
                    subtitle = wrap_label(subtitle_da, 70))
    if (isTRUE(params$da$label_nhoods)) {
      sig_ids <- da_res$Nhood[which(da_res$SpatialFDR < da_alpha)]
      p_graph <- add_nhood_labels(p_graph, milo_obj, sig_ids)
    }
    p_graph_combined <- p_umap_timepoint + p_graph
    print(p_graph_combined)
    save_plot(p_graph_combined, paste0("da_graph_", contrast_name, ".png"), out_dir, width = 15, height = 7)

    # --- beeswarm by cluster -----------------------------------------------------
    # plotDAbeeswarm() fails when no neighborhood is significant (nothing to
    # color), so that case gets a plain grey version
    p_bee <- if (n_sig > 0) {
      miloR::plotDAbeeswarm(da_res, group.by = "cluster_label_mixed", alpha = da_alpha)
    } else {
      ggplot2::ggplot(da_res, ggplot2::aes(x = logFC, y = cluster_label_mixed)) +
        ggplot2::geom_point(position = ggplot2::position_jitter(height = 0.2, width = 0, seed = params$da$seed),
                            color = "#808080", size = 1.5) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dotted", color = "#999999") +
        ggplot2::labs(x = "Log Fold Change", y = "Cluster") +
        ggplot2::theme_bw(base_size = 12)
    }
    p_bee <- p_bee +
      ggplot2::labs(title = paste0("Neighborhood logFC by cluster \u2014 ", levels_used[1], " vs ", levels_used[2]),
                    subtitle = wrap_label(paste0(if (n_sig == 0) "No neighborhood significant; " else "",
                                                 "Colored if SpatialFDR < ", da_alpha,
                                                 "; 'Mixed' = no label >= ",
                                                 params$da$mixed_threshold, " of cells"), 70))
    print(p_bee)
    save_plot(p_bee, paste0("da_beeswarm_", contrast_name, ".png"), out_dir, width = 8, height = 6)

    # --- p-value histogram (test calibration) ----------------------------------------
    p_hist <- ggplot2::ggplot(da_res, ggplot2::aes(x = PValue)) +
      ggplot2::geom_histogram(breaks = seq(0, 1, 0.05), fill = "#666666", color = "#FFFFFF") +
      ggplot2::labs(title = paste0("Uncorrected p-values \u2014 ", contrast_name),
                    subtitle = "Roughly flat with a peak near 0 is expected; other shapes suggest a poorly fitting model",
                    x = "PValue", y = "Neighborhoods") +
      ggplot2::theme_classic(base_size = 12)
    print(p_hist)
    save_plot(p_hist, paste0("da_pvalue_histogram_", contrast_name, ".png"), out_dir, width = 7, height = 5)

    # --- per-cluster summary --------------------------------------------------------
    label_summary <- da_label_summary(da_res, "cluster_label_mixed", da_alpha)
    label_gt <- gt_generic(label_summary,
                           title = paste0("DA neighborhoods per cluster \u2014 ", levels_used[1], " vs ", levels_used[2]),
                           subtitle = paste0("Significant: SpatialFDR < ", da_alpha, "; enriched = more abundant in ", levels_used[1]))
    print(label_gt)
    save_gt(label_gt, paste0("da_cluster_summary_", contrast_name, ".png"), out_dir)

    # --- why (not) significant: within-condition variability ------------------------------
    if (is_primary) {
      p_var <- plot_nhood_variability(
        da_res, da_alpha, levels_used[1], levels_used[2], n_samples = nrow(design_df),
        title = paste0("Why is a neighborhood (not) significant? \u2014 ", levels_used[1], " vs ", levels_used[2]),
        subtitle = wrap_label(paste0("Left: share of the neighborhood's (size-normalized) abundance from its largest ",
                                     "sample; dashed = all samples equal (1/", nrow(design_df), "). Right: variation ",
                                     "between samples of the same condition. Large |logFC| with high values on ",
                                     "either axis = driven by few samples"), 110)
      )
      print(p_var)
      save_plot(p_var, paste0("da_variability_", contrast_name, ".png"), out_dir, width = 11, height = 5.5)

      # per-sample abundance of the most extreme neighborhoods, significant or not
      n_show <- params$da$n_nhoods_to_show
      extreme <- utils::head(da_res$Nhood[order(-abs(da_res$logFC))], n_show)
      extreme_label <- setNames(
        sprintf("Nhood %d: %s, logFC %.1f, SpatialFDR %.2f", da_res$Nhood, da_res$cluster_label_mixed,
                da_res$logFC, da_res$SpatialFDR),
        as.character(da_res$Nhood))[as.character(extreme)]
      extreme_abundance <- abundance[abundance$Nhood %in% extreme, ]
      extreme_abundance$set <- factor(extreme_label[as.character(extreme_abundance$Nhood)], levels = extreme_label)
      extreme_abundance$timepoint <- as.character(extreme_abundance$condition)
      p_extreme <- plot_sample_abundance(
        extreme_abundance, timepoint_levels = rev(levels(design_df$condition)),
        repeated_timepoints = if (design$paired) levels_used else character(0),
        timepoint_colors = timepoint_colors,
        title = paste0("Per-sample abundance of the ", length(extreme), " neighborhoods with the largest |logFC| \u2014 ",
                       contrast_name),
        subtitle = "Each point = one sample: % of its cells in the neighborhood; this is what the DA test compares"
      )
      print(p_extreme)
      save_plot(p_extreme, paste0("da_extreme_nhoods_abundance_", contrast_name, ".png"), out_dir,
                width = 12, height = 3 * ceiling(length(extreme) / 4) + 1)

      # which cells make up each significant neighborhood
      sig_ids <- da_res$Nhood[which(da_res$SpatialFDR < da_alpha)]
      if (length(sig_ids) > 0) {
        shown <- utils::head(sig_ids[order(-abs(da_res$logFC[match(sig_ids, da_res$Nhood)]))], n_show)
        nhood_cells <- lapply(shown, function(nh) rownames(da_membership)[da_membership[, nh] > 0])
        names(nhood_cells) <- sprintf("Nhood %d (logFC %.1f)", shown, da_res$logFC[match(shown, da_res$Nhood)])
        p_members <- plot_da_group_umap(
          data_obj, nhood_cells,
          title = paste0("Cells of significant neighborhoods \u2014 ", contrast_name),
          subtitle = paste0("Up to ", n_show, " with the largest |logFC|; red = member cells")
        )
        print(p_members)
        save_plot(p_members, paste0("da_sig_nhood_cells_umap_", contrast_name, ".png"), out_dir,
                  width = min(15, 3.5 * min(length(shown), 4) + 1), height = 3.5 * ceiling(length(shown) / 4) + 1)
      }
    }

    # --- program scores of significant neighborhoods ---------------------------------
    sig_res <- da_res[which(da_res$SpatialFDR < da_alpha), ]
    if (!is.null(nhood_scores) && nrow(sig_res) >= 2) {
      row_labels <- setNames(paste0(sig_res$cluster_label_mixed, ", ",
                                    ifelse(sig_res$logFC > 0, levels_used[1], levels_used[2]), "-enriched"),
                             as.character(sig_res$Nhood))
      p_programs <- plot_nhood_program_heatmap(
        sig_res, score_cols, row_label = row_labels,
        row_order = sig_res$Nhood[order(sig_res$logFC, decreasing = TRUE)],
        title = paste0("Programs in DA neighborhoods \u2014 ", levels_used[1], " vs ", levels_used[2]),
        subtitle = paste0("Ordered by DA logFC (top = most enriched in ", levels_used[1],
                          "); mean score z-scored within program across the neighborhoods shown; ",
                          "'Mixed' = no cluster >= ", params$da$mixed_threshold, " of cells")
      )
      print(p_programs)
      save_plot(p_programs, paste0("da_program_heatmap_", contrast_name, ".png"), out_dir,
                width = max(7, 0.6 * length(score_cols) + 4), height = max(4, 0.3 * nrow(sig_res) + 2))
    } else if (run_scores) {
      message("  program heatmap skipped (fewer than 2 significant neighborhoods)")
    }

    # --- programs across all neighborhoods vs DA logFC --------------------------------
    if (!is.null(nhood_scores)) {
      p_prog_da <- plot_program_vs_da(
        da_res, score_cols, da_alpha, levels_used[1], levels_used[2],
        title = paste0("Module scores across neighborhoods \u2014 ", levels_used[1], " vs ", levels_used[2]),
        subtitle = "Each point = one neighborhood (mean score of its cells); descriptive, neighborhoods overlap"
      )
      print(p_prog_da)
      n_prog <- length(score_cols)
      save_plot(p_prog_da, paste0("da_programs_vs_logfc_", contrast_name, ".png"), out_dir,
                width = min(16, 3.2 * min(n_prog, 4) + 1), height = 3 * ceiling(n_prog / 4) + 1.5)
    }

    # --- composition of significant neighborhoods ---------------------------------
    composition_sheets <- list()
    if (is_primary && nrow(sig_res) > 0) {
      nhood_order <- sig_res$Nhood[order(sig_res$logFC)]
      nhood_text <- setNames(sprintf("%d (logFC %.2f)", sig_res$Nhood, sig_res$logFC), as.character(sig_res$Nhood))
      plot_height <- max(4, 0.25 * nrow(sig_res) + 2)

      composition_specs <- list(
        cluster   = list(labels = da_meta$cluster_label, colors = cluster_colors, legend = "Cluster"),
        timepoint = list(labels = da_meta[[timepoint_col]], colors = timepoint_colors, legend = timepoint_col)
      )
      for (spec_name in names(composition_specs)) {
        spec <- composition_specs[[spec_name]]
        comp <- nhood_label_composition(da_membership, spec$labels, sig_res$Nhood)
        p_comp <- plot_nhood_composition(
          comp, nhood_order, label_colors = spec$colors, nhood_text = nhood_text, legend_title = spec$legend,
          title = paste0(spec$legend, " composition of DA neighborhoods \u2014 ", levels_used[1], " vs ", levels_used[2]),
          subtitle = paste0("SpatialFDR < ", da_alpha, "; ordered by logFC (top = most enriched in ", levels_used[1], ")")
        )
        print(p_comp)
        save_plot(p_comp, paste0("da_composition_", spec_name, "_", contrast_name, ".png"), out_dir,
                  width = 8, height = plot_height)
        composition_sheets[[paste0("sig_", spec_name, "_composition")]] <- comp %>%
          dplyr::select(Nhood, label, fraction) %>%
          dplyr::mutate(fraction = round(fraction, 3)) %>%
          tidyr::pivot_wider(names_from = label, values_from = fraction, values_fill = 0) %>%
          dplyr::left_join(sig_res[, c("Nhood", "logFC", "SpatialFDR")], by = "Nhood") %>%
          dplyr::relocate(Nhood, logFC, SpatialFDR) %>%
          dplyr::arrange(logFC)
      }
    }

    # --- QC per neighborhood: vs DA logFC and vs program scores --------------------------
    if (!is.null(nhood_qc)) {
      p_qc_da <- plot_program_vs_da(
        da_res, qc_cols_nhood, da_alpha, levels_used[1], levels_used[2], y_label = "Mean per neighborhood",
        title = paste0("QC metrics across neighborhoods \u2014 ", levels_used[1], " vs ", levels_used[2]),
        subtitle = "A DA signal that follows a QC metric may be technical; descriptive, neighborhoods overlap"
      )
      print(p_qc_da)
      save_plot(p_qc_da, paste0("da_qc_vs_logfc_", contrast_name, ".png"), out_dir,
                width = min(16, 3.2 * min(length(qc_cols_nhood), 4) + 1), height = 3 * ceiling(length(qc_cols_nhood) / 4) + 1.5)
    }

    # --- meta-neighborhoods + marker DE------------------------------------------------------
    if (is_primary && run_markers && nrow(sig_res) > 0) {
      marker_params <- params$da$markers
      marker_dir <- file.path(out_dir, "markers")
      dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)

      groups <- group_da_nhoods(da_membership, da_res, da_alpha, marker_params$min_overlap)
      da_res <- dplyr::left_join(da_res, groups[, c("Nhood", "da_group")], by = "Nhood")
      sig_res <- dplyr::left_join(sig_res, groups[, c("Nhood", "da_group")], by = "Nhood")

      group_cells <- lapply(split(groups$Nhood, groups$da_group), function(nh) {
        rownames(da_membership)[Matrix::rowSums(da_membership[, nh, drop = FALSE]) > 0]
      })
      group_cells <- group_cells[order(as.integer(sub("G", "", names(group_cells))))]
      group_info <- dplyr::bind_rows(lapply(names(group_cells), function(g) {
        nh <- groups$Nhood[groups$da_group == g]
        data.frame(
          da_group  = g,
          direction = paste0(ifelse(groups$direction[groups$da_group == g][1] == "enriched",
                                    levels_used[1], levels_used[2]), "-enriched"),
          n_nhoods  = length(nh),
          nhoods    = paste(nh, collapse = ", "),
          n_cells   = length(group_cells[[g]]),
          dominant_cluster = names(which.max(table(da_meta[group_cells[[g]], "cluster_label"]))),
          n_patients = dplyr::n_distinct(da_meta[group_cells[[g]], patient_col])
        )
      }))
      group_titles <- setNames(paste0(group_info$da_group, ": ", group_info$direction, ", ",
                                      group_info$dominant_cluster, " (", group_info$n_cells, " cells)"),
                               group_info$da_group)
      message("    meta-neighborhoods: ", paste(group_titles, collapse = "; "))
      analysis_groups[[contrast_name]] <- list(cells = group_cells, titles = group_titles, info = group_info,
                                               nhoods = split(groups$Nhood, groups$da_group))

      p_groups <- plot_da_group_umap(
        data_obj, setNames(group_cells, group_titles[names(group_cells)]),
        title = paste0("Cells of meta-neighborhoods \u2014 ", levels_used[1], " vs ", levels_used[2]),
        subtitle = "Meta-neighborhood = significant neighborhoods with the same direction that share cells, merged"
      )
      print(p_groups)
      save_plot(p_groups, paste0("da_groups_umap_", contrast_name, ".png"), marker_dir,
                width = min(15, 4 * length(group_cells) + 1), height = 4.5 * ceiling(length(group_cells) / 3))

      marker_sheets <- list(groups = group_info,
                            cells = dplyr::bind_rows(lapply(names(group_cells), function(g)
                              data.frame(da_group = g, cell = group_cells[[g]]))))

      # contrast_cells_only = TRUE: characterize only the cells of this contrast's
      # two timepoints (the cells the DA test counted), against other cells of
      # those timepoints; FALSE: all cells of the region (the state as a whole)
      marker_obj <- if (isTRUE(marker_params$contrast_cells_only)) {
        subset(data_obj, cells = colnames(data_obj)[as.character(data_obj[[timepoint_col, drop = TRUE]]) %in% levels_used])
      } else {
        data_obj
      }
      for (g in names(group_cells)) {
        sheets_g <- characterize_cell_set(
          marker_obj, intersect(group_cells[[g]], colnames(marker_obj)), set_name = g, title = group_titles[[g]],
          out_dir = marker_dir,
          file_suffix = contrast_name, sample_col = "milo_sample", timepoint_col = timepoint_col,
          timepoint_colors = timepoint_colors, marker_params = marker_params,
          exclude_gene_patterns = marker_exclude_patterns
        )
        if (length(sheets_g) > 0) names(sheets_g) <- paste0(g, "_", names(sheets_g))
        marker_sheets <- c(marker_sheets, sheets_g)
      }
      save_xlsx(marker_sheets, paste0("markers_", contrast_name, ".xlsx"), marker_dir)
    }

    save_xlsx(c(list(all_nhoods = da_res, significant = sig_res, cluster_summary = label_summary,
                     sample_abundance = abundance),
                composition_sheets, list(design = design_df)),
              paste0("milo_da_", contrast_name, ".xlsx"), out_dir)
  }

  #-------------------------------------------------------------------------------
  # 2b. Primary only: meta-neighborhoods across contrasts, chosen neighborhoods
  #-------------------------------------------------------------------------------
  if (is_primary) {
    da_results_all <- analysis_results
    milo_da <- milo_obj

    # every meta-neighborhood (defined in one contrast) followed across all
    # contrasts: its neighborhoods' results in each, and each sample's share of
    # its cells - do the same cells differ at 3M and at FU?
    if (length(analysis_groups) > 0) {
      meta_dir <- file.path(analysis_dir, "meta_nhoods")
      dir.create(meta_dir, recursive = TRUE, showWarnings = FALSE)
      meta_cells <- list()
      meta_rows <- list()
      for (defined_in in names(analysis_groups)) {
        grp <- analysis_groups[[defined_in]]
        for (g in names(grp$cells)) {
          key <- paste0(defined_in, " ", g)
          meta_cells[[paste0(key, ": ", sub("^G[0-9]+: ", "", grp$titles[[g]]))]] <- grp$cells[[g]]
          nh <- grp$nhoods[[g]]
          row <- data.frame(meta_nhood = key, defined_in = defined_in, title = grp$titles[[g]],
                            nhoods = paste(nh, collapse = ", "), n_cells = length(grp$cells[[g]]))
          for (ct in names(analysis_results)) {
            r <- analysis_results[[ct]][analysis_results[[ct]]$Nhood %in% nh, ]
            row[[paste0(ct, "_median_logFC")]] <- round(stats::median(r$logFC), 2)
            row[[paste0(ct, "_n_significant")]] <- paste0(sum(r$SpatialFDR < da_alpha, na.rm = TRUE), "/", nrow(r))
          }
          meta_rows[[key]] <- row
        }
      }
      meta_table <- dplyr::bind_rows(meta_rows)
      meta_gt <- gt_generic(
        meta_table, title = "Meta-neighborhoods across contrasts",
        subtitle = paste0("Each defined in one contrast; its neighborhoods' median logFC and number significant ",
                          "(SpatialFDR < ", da_alpha, ") in every contrast")
      )
      print(meta_gt)
      save_gt(meta_gt, "meta_nhoods_across_contrasts.png", meta_dir, vwidth = 1600)

      meta_abundance <- set_sample_abundance(meta_cells, da_meta, patient_col, timepoint_col)
      p_meta_abundance <- plot_sample_abundance(
        meta_abundance, analysis_timepoints, repeated_timepoints, timepoint_colors,
        title = "Per-sample abundance of each meta-neighborhood",
        subtitle = "Each point = one sample: % of its cells in the meta-neighborhood; lines join a patient's repeated samples"
      )
      print(p_meta_abundance)
      save_plot(p_meta_abundance, "meta_nhoods_sample_abundance.png", meta_dir,
                width = min(16, 4 * min(length(meta_cells), 4) + 1), height = 3.5 * ceiling(length(meta_cells) / 4) + 1)
      save_xlsx(list(meta_nhoods = meta_table, sample_abundance = meta_abundance),
                "meta_nhoods.xlsx", meta_dir)
    }

    nhood_sets <- params$da$markers$nhood_sets
    if (run_markers && length(nhood_sets) > 0) {
      chosen_dir <- file.path(analysis_dir, "chosen_nhoods")
      marker_params <- params$da$markers
      set_cells <- list()
      set_rows <- list()
      chosen_sheets <- list()

      for (set_name in names(nhood_sets)) {
        ids <- as.integer(nhood_sets[[set_name]])
        if (any(is.na(ids)) || any(ids < 1 | ids > ncol(da_membership))) {
          warning("nhood_sets$", set_name, ": neighborhood IDs must be between 1 and ", ncol(da_membership), " - skipped")
          next
        }
        cells <- rownames(da_membership)[Matrix::rowSums(da_membership[, ids, drop = FALSE]) > 0]
        set_cells[[set_name]] <- cells
        tp_share <- round(100 * prop.table(table(as.character(da_meta[cells, timepoint_col]))))
        da_text <- vapply(names(analysis_results), function(ct) {
          r <- analysis_results[[ct]][analysis_results[[ct]]$Nhood %in% ids, ]
          paste0(ct, ": ", paste0(r$Nhood, " (logFC ", round(r$logFC, 2), ", SpatialFDR ", signif(r$SpatialFDR, 2), ")",
                                  collapse = "; "))
        }, character(1))
        set_rows[[set_name]] <- data.frame(
          set = set_name, nhoods = paste(ids, collapse = ", "), n_cells = length(cells),
          dominant_cluster = names(which.max(table(as.character(da_meta[cells, "cluster_label"])))),
          timepoints = paste0(names(tp_share), " ", tp_share, "%", collapse = ", "),
          n_patients = dplyr::n_distinct(da_meta[cells, patient_col]),
          da_results = paste(da_text, collapse = " | ")
        )
        title <- paste0(set_name, ": neighborhood", if (length(ids) > 1) "s " else " ",
                        paste(ids, collapse = ", "), " (", length(cells), " cells)")
        sheets_s <- characterize_cell_set(
          data_obj, cells, set_name = set_name, title = title, out_dir = chosen_dir,
          sample_col = "milo_sample", timepoint_col = timepoint_col, timepoint_colors = timepoint_colors,
          marker_params = marker_params, exclude_gene_patterns = marker_exclude_patterns
        )
        if (length(sheets_s) > 0) names(sheets_s) <- paste0(set_name, "_", names(sheets_s))
        chosen_sheets <- c(chosen_sheets, sheets_s)
      }

      if (length(set_cells) > 0) {
        p_sets <- plot_da_group_umap(data_obj, set_cells, title = "Cells of the chosen neighborhoods")
        print(p_sets)
        save_plot(p_sets, "chosen_nhoods_umap.png", chosen_dir,
                  width = min(15, 4 * length(set_cells) + 1), height = 4.5 * ceiling(length(set_cells) / 3))
        set_abundance <- set_sample_abundance(set_cells, da_meta, patient_col, timepoint_col)
        p_set_abundance <- plot_sample_abundance(
          set_abundance, analysis_timepoints, repeated_timepoints, timepoint_colors,
          title = "Per-sample abundance of the chosen neighborhoods",
          subtitle = "Each point = one sample: % of its cells in the set; lines join a patient's repeated samples"
        )
        print(p_set_abundance)
        save_plot(p_set_abundance, "chosen_nhoods_sample_abundance.png", chosen_dir,
                  width = min(15, 4 * length(set_cells) + 1), height = 4 * ceiling(length(set_cells) / 3))
        save_xlsx(c(list(sets = dplyr::bind_rows(set_rows), sample_abundance = set_abundance), chosen_sheets),
                  "chosen_nhoods_markers.xlsx", chosen_dir)
      }
    }


  }

  #-------------------------------------------------------------------------------
  # 2c. Primary only: neighborhood groups from ALL neighborhoods
  #-------------------------------------------------------------------------------
  # groups = communities of the neighborhood-overlap graph, defined without
  # DA results; each cell then belongs to one group, so groups can be tested
  # per sample and described without circularity
  if (is_primary && isTRUE(steps$nhood_groups)) {
    group_dir <- file.path(analysis_dir, "nhood_groups")
    dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)
    group_params <- params$da$nhood_groups

    nhood_group <- group_nhoods_by_overlap(da_membership, group_params$resolution, group_params$seed)
    cell_group  <- assign_cells_to_groups(da_membership, nhood_group)
    group_levels <- unique(nhood_group[order(as.numeric(sub("NG", "", nhood_group)))])
    group_cells <- split(names(cell_group)[!is.na(cell_group)], factor(cell_group[!is.na(cell_group)], levels = group_levels))
    message("  neighborhood groups: ", length(group_levels), " (Leiden resolution ", group_params$resolution,
            "); ", sum(!is.na(cell_group)), " of ", length(cell_group), " cells assigned")

    # how much the grouping depends on the exact resolution
    stability <- group_resolution_stability(da_membership, group_params$resolutions_checked,
                                            group_params$resolution, group_params$seed)
    stability_gt <- gt_generic(stability, title = "Neighborhood groups at nearby Leiden resolutions",
                               subtitle = "ARI of the cells' group assignment against the resolution used (1 = identical)")
    print(stability_gt)
    save_gt(stability_gt, "group_resolution_stability.png", group_dir)

    # --- what each neighborhood and group is made of ------------------------------------
    da_meta$nhood_group <- cell_group[rownames(da_meta)]
    nhood_cells <- lapply(seq_len(ncol(da_membership)), function(nh) rownames(da_membership)[da_membership[, nh] > 0])
    names(nhood_cells) <- as.character(seq_len(ncol(da_membership)))
    timepoint_tp <- factor(da_meta[[timepoint_col]], levels = analysis_timepoints)
    da_meta[[timepoint_col]] <- timepoint_tp

    majority_label <- function(cells) {
      tab <- sort(table(as.character(da_meta[cells, "cluster_label"])), decreasing = TRUE)
      if (length(tab) == 0 || tab[1] / sum(tab) < params$da$mixed_threshold) "Mixed" else names(tab)[1]
    }
    comp_tp_nhood  <- set_composition(nhood_cells, da_meta, timepoint_col, normalize = TRUE)
    comp_cl_nhood  <- set_composition(nhood_cells, da_meta, "cluster_label")
    comp_tp_group  <- set_composition(group_cells, da_meta, timepoint_col, normalize = TRUE)
    comp_cl_group  <- set_composition(group_cells, da_meta, "cluster_label")

    # --- group-level DA (and the cluster-level comparator), per contrast --------------------
    group_da <- list()
    cluster_da <- list()
    for (contrast_name in names(analysis$contrasts)) {
      design <- contrast_design(sample_design, analysis$contrasts[[contrast_name]], patient_col, timepoint_col,
                                repeated_timepoints, paired = params$da$paired_design)
      group_da[[contrast_name]] <- test_set_abundance(cell_group, da_meta, design)
      cluster_da[[contrast_name]] <- test_set_abundance(
        setNames(as.character(da_meta$cluster_label), rownames(da_meta)), da_meta, design)
    }
    for (what in c("group", "cluster")) {
      tables <- if (what == "group") group_da else cluster_da
      combined <- dplyr::bind_rows(lapply(names(tables), function(ct) cbind(contrast = ct, tables[[ct]][, 1:5])))
      combined[, c("logFC", "PValue", "FDR")] <- signif(combined[, c("logFC", "PValue", "FDR")], 3)
      tab_gt <- gt_generic(
        combined,
        title = if (what == "group") "Abundance test per neighborhood group" else "Abundance test per Seurat cluster (comparator)",
        subtitle = paste0("Cells counted per sample, same designs as the neighborhood tests; BH across ",
                          what, "s within each contrast; logFC > 0 = more abundant in the first group")
      )
      print(tab_gt)
      save_gt(tab_gt, paste0(what, "_abundance_tests.png"), group_dir)
    }

    # --- DA of groups and of Seurat clusters, side by side --------------------------------------
    p_set_da <- plot_set_da_dots(
      list("Neighborhood groups" = group_da, "Seurat clusters (comparator)" = cluster_da), da_alpha,
      title = "Abundance tests: neighborhood groups vs Seurat clusters",
      subtitle = paste0("Cells per sample, same designs as the neighborhood tests. Color = logFC (red = more in the ",
                        "first timepoint); size = -log10 FDR; black outline = FDR < ", da_alpha, "; number = FDR")
    )
    print(p_set_da)
    save_plot(p_set_da, "group_vs_cluster_da.png", group_dir, width = 12, height = 5)

    # --- heatmaps: one question - which programs characterize each row ----------------------
    # rows grouped by Leiden group (gaps), hclust within groups; columns hclust
    group_palette <- setNames(resolve_palette("group_umap", n = length(group_levels)), group_levels)
    label_palette <- c(cluster_colors, Mixed = "#B3B3B3")
    group_label <- vapply(group_cells, majority_label, character(1))   # majority cluster of each group
    if (length(score_cols) > 0) {
      net_cols <- grep("_NETScore$", score_cols, value = TRUE)
      if (length(net_cols) == 0) net_cols <- score_cols

      # neighborhoods: annotated only with their cluster label and group
      nhood_label <- vapply(nhood_cells, majority_label, character(1))
      annotation <- data.frame(group = nhood_group[names(nhood_cells)], cluster = nhood_label,
                               row.names = names(nhood_cells))
      ann_colors <- list(group = group_palette,
                         cluster = label_palette[intersect(names(label_palette), unique(nhood_label))])
      nhood_scores_mat <- as.matrix(nhood_scores[match(as.numeric(rownames(annotation)), nhood_scores$Nhood), net_cols])
      rownames(nhood_scores_mat) <- rownames(annotation)
      nhood_order <- plot_annotated_heatmap(
        nhood_scores_mat, annotation, ann_colors, row_groups = setNames(annotation$group, rownames(annotation)),
        title = paste0("Net program scores of all neighborhoods (z-scored per program)\n",
                       "cluster = majority cluster if >= ", params$da$mixed_threshold * 100, "% of cells, else Mixed"),
        path = file.path(group_dir, "heatmap_nhoods_net_programs.png")
      )

      # groups: mean over each group's cells
      group_scores <- t(vapply(group_cells, function(cells) {
        colMeans(score_values[intersect(cells, rownames(score_values)), net_cols, drop = FALSE])
      }, numeric(length(net_cols))))
      if (length(net_cols) == 1) group_scores <- t(group_scores)
      colnames(group_scores) <- net_cols
      rownames(group_scores) <- paste0(names(group_cells), " (", lengths(group_cells), " cells)")
      g_annotation <- data.frame(cluster = group_label, row.names = rownames(group_scores))
      group_order <- plot_annotated_heatmap(
        group_scores, g_annotation,
        list(cluster = label_palette[intersect(names(label_palette), unique(group_label))]),
        title = "Net program scores of neighborhood groups (z-scored per program)",
        path = file.path(group_dir, "heatmap_groups_net_programs.png")
      )
      group_order <- sub(" \\(.*$", "", group_order)
    } else {
      nhood_order <- names(nhood_cells)[order(as.numeric(sub("NG", "", nhood_group[names(nhood_cells)])))]
      group_order <- group_levels
      message("  [skip] program heatmaps (no functional scores - steps$nhood_program_scores)")
    }

    # --- what distinguishes each group: one-vs-rest DE, as for the Seurat clusters --------------
    marker_params <- params$da$markers
    group_obj <- subset(data_obj, cells = names(cell_group)[!is.na(cell_group)])
    Seurat::Idents(group_obj) <- factor(cell_group[colnames(group_obj)], levels = group_levels)
    group_markers <- tryCatch(
      Seurat::FindAllMarkers(group_obj, logfc.threshold = marker_params$log2fc_cutoff,
                             min.pct = marker_params$min_pct, verbose = FALSE),
      error = function(e) { message("  group markers failed: ", conditionMessage(e)); NULL }
    )
    if (!is.null(group_markers) && nrow(group_markers) > 0) {
      if (length(marker_exclude_patterns) > 0) {
        group_markers <- group_markers[!grepl(paste(marker_exclude_patterns, collapse = "|"), group_markers$gene), ]
      }
      marker_dir <- file.path(group_dir, "markers")
      dir.create(marker_dir, recursive = TRUE, showWarnings = FALSE)
      marker_sheets <- list()
      for (g in group_levels) {
        de_g <- group_markers[group_markers$cluster == g, , drop = FALSE]
        de_g <- de_g[order(de_g$p_val_adj, -abs(de_g$avg_log2FC)), , drop = FALSE]
        marker_sheets[[paste0(g, "_", group_label[[g]])]] <- de_g[, c("gene", setdiff(names(de_g), "gene"))]
        if (nrow(de_g) == 0) next
        p_volcano <- plot_marker_volcano(
          de_g, log2fc_cutoff = marker_params$log2fc_cutoff, padj_cutoff = marker_params$padj_cutoff,
          title = paste0(g, " (mostly ", group_label[[g]], ", ", length(group_cells[[g]]), " cells) vs all other grouped cells"),
          subtitle = "One-vs-rest DE (Seurat FindAllMarkers); descriptive"
        )
        print(p_volcano)
        save_plot(p_volcano, paste0("volcano_", g, ".png"), marker_dir, width = 8, height = 7)
        rownames(de_g) <- make.unique(de_g$gene)
        de_gt <- gt_de(
          de_g,
          title = paste0(g, " (mostly ", group_label[[g]], ")"),
          subtitle = "Top genes vs all other grouped cells",
          pct1_label = paste("%", g),
          pct2_label = "% other groups",
          gene_from_rownames = FALSE
        )
        print(de_gt)
        save_gt(de_gt, paste0("de_table_", g, ".png"), marker_dir)
      }
      names(marker_sheets) <- make.unique(substr(names(marker_sheets), 1, 28), sep = "_")
      save_xlsx(marker_sheets, "nhood_group_markers.xlsx", marker_dir)
    }

    # --- DE between every pair of groups (as pairwise_de for Seurat clusters) ------------------
    if (isTRUE(params$da$nhood_groups$pairwise_de) && length(group_levels) > 1) {
      pair_dir <- file.path(group_dir, "markers", "pairwise")
      dir.create(pair_dir, recursive = TRUE, showWarnings = FALSE)
      pair_sheets <- list()
      summary_rows <- list()
      for (pair in utils::combn(group_levels, 2, simplify = FALSE)) {
        pair_name <- paste(pair[1], "vs", pair[2])
        de_p <- tryCatch(
          Seurat::FindMarkers(group_obj, ident.1 = pair[1], ident.2 = pair[2],
                              logfc.threshold = marker_params$log2fc_cutoff,
                              min.pct = marker_params$min_pct, verbose = FALSE),
          error = function(e) { message("  ", pair_name, " skipped: ", conditionMessage(e)); NULL }
        )
        if (is.null(de_p) || nrow(de_p) == 0) next
        de_p$gene <- rownames(de_p)
        if (length(marker_exclude_patterns) > 0) {
          de_p <- de_p[!grepl(paste(marker_exclude_patterns, collapse = "|"), de_p$gene), ]
        }
        de_p <- de_p[order(de_p$p_val_adj, -abs(de_p$avg_log2FC)), ]
        sig <- de_p[de_p$p_val_adj < marker_params$padj_cutoff & abs(de_p$avg_log2FC) > marker_params$log2fc_cutoff, ]
        summary_rows[[pair_name]] <- data.frame(
          comparison = pair_name, n_DE = nrow(sig),
          higher_in_first = sum(sig$avg_log2FC > 0), higher_in_second = sum(sig$avg_log2FC < 0),
          first = paste0(pair[1], " (", group_label[[pair[1]]], ")"),
          second = paste0(pair[2], " (", group_label[[pair[2]]], ")"))
        pair_sheets[[pair_name]] <- de_p[, c("gene", setdiff(names(de_p), "gene"))]
        p_pair <- plot_marker_volcano(
          de_p, log2fc_cutoff = marker_params$log2fc_cutoff, padj_cutoff = marker_params$padj_cutoff,
          title = paste0(pair[1], " (", group_label[[pair[1]]], ") vs ", pair[2], " (", group_label[[pair[2]]], ")"),
          subtitle = "logFC > 0 = higher in the first group; descriptive"
        )
        print(p_pair)
        save_plot(p_pair, paste0("volcano_", pair[1], "_vs_", pair[2], ".png"), pair_dir, width = 8, height = 7)
      }
      if (length(pair_sheets) > 0) {
        names(pair_sheets) <- make.unique(substr(names(pair_sheets), 1, 28), sep = "_")
        save_xlsx(c(list(summary = do.call(rbind, summary_rows)), pair_sheets),
                  "nhood_group_pairwise_de.xlsx", pair_dir)
      }
    }

    # --- composition bars: all neighborhoods (heatmap order) and all groups ------------------
    for (level in c("nhoods", "groups")) {
      comp_cl <- if (level == "nhoods") comp_cl_nhood else comp_cl_group
      comp_tp <- if (level == "nhoods") comp_tp_nhood else comp_tp_group
      ord <- if (level == "nhoods") nhood_order else group_order
      what <- if (level == "nhoods") "neighborhood" else "neighborhood group"
      p_cl <- plot_set_composition(comp_cl, ord, label_palette, title = paste("Cluster composition of each", what))
      p_tp <- plot_set_composition(comp_tp[comp_tp$type == "raw", ], ord, timepoint_colors,
                                   title = paste("Timepoint composition of each", what), subtitle = "Raw: share of its cells")
      p_tpn <- plot_set_composition(comp_tp[comp_tp$type == "normalized", ], ord, timepoint_colors,
                                    title = paste("Timepoint enrichment of each", what),
                                    subtitle = "Normalized by each timepoint's total cells (equal shares = no enrichment)")
      p_comp <- patchwork::wrap_plots(p_cl, p_tp, p_tpn, ncol = 1)
      print(p_comp)
      save_plot(p_comp, paste0("composition_", level, ".png"), group_dir,
                width = if (level == "nhoods") 16 else 8, height = 11)
    }

    # --- where the groups are -------------------------------------------------------------------
    data_obj$nhood_group <- factor(cell_group[colnames(data_obj)], levels = group_levels)
    p_group_umap <- Seurat::DimPlot(data_obj, group.by = "nhood_group", cols = group_palette, na.value = "#E5E5E5") +
      ggplot2::labs(title = "Cells by neighborhood group",
                    subtitle = "Each cell in the group most of its neighborhoods belong to; grey = in no neighborhood")
    print(p_group_umap)
    save_plot(p_group_umap, "nhood_groups_umap.png", group_dir, width = 8, height = 6.5)

    p_group_facets <- plot_da_group_umap(data_obj, group_cells, title = "Cells of each neighborhood group")
    print(p_group_facets)
    save_plot(p_group_facets, "nhood_groups_umap_facets.png", group_dir,
              width = min(16, 3.5 * min(length(group_cells), 4) + 1), height = 3.5 * ceiling(length(group_cells) / 4) + 1)

    graph_res <- analysis_results[[1]]
    graph_res$NhoodGroup <- nhood_group[as.character(graph_res$Nhood)]
    p_group_graph <- tryCatch(
      miloR::plotNhoodGroups(milo_obj, graph_res, layout = "UMAP") +
        ggplot2::labs(title = "Neighborhood graph colored by group"),
      error = function(e) { message("  group graph skipped: ", conditionMessage(e)); NULL }
    )
    if (!is.null(p_group_graph)) {
      print(p_group_graph)
      save_plot(p_group_graph, "nhood_groups_graph.png", group_dir, width = 9, height = 7)
    }

    # --- per-sample abundance of every group -------------------------------------------------------
    group_abundance <- set_sample_abundance(group_cells, da_meta, patient_col, timepoint_col)
    p_group_abundance <- plot_sample_abundance(
      group_abundance, analysis_timepoints, repeated_timepoints, timepoint_colors,
      title = "Per-sample abundance of each neighborhood group",
      subtitle = "Each point = one sample: % of its cells in the group; lines join a patient's repeated samples"
    )
    print(p_group_abundance)
    save_plot(p_group_abundance, "nhood_groups_sample_abundance.png", group_dir,
              width = min(16, 4 * min(length(group_cells), 4) + 1), height = 3.5 * ceiling(length(group_cells) / 4) + 1)

    # --- cell quality per group: is a group defined by quality rather than state? ------------
    group_qc <- NULL
    if (length(qc_cols_nhood) > 0) {
      qc_df <- data.frame(group = data_obj$nhood_group, timepoint = data_obj[[timepoint_col, drop = TRUE]],
                          data_obj[[]][, qc_cols_nhood, drop = FALSE], check.names = FALSE)
      qc_df <- qc_df[!is.na(qc_df$group), ]
      group_qc <- qc_df %>%
        dplyr::group_by(group) %>%
        dplyr::summarise(n_cells = dplyr::n(),
                         dplyr::across(dplyr::all_of(qc_cols_nhood), ~ stats::median(.x, na.rm = TRUE)),
                         .groups = "drop") %>%
        dplyr::mutate(majority_cluster = unname(group_label[as.character(group)]), .after = group)
      group_qc_gt <- gt_generic(
        format_qc_values(group_qc, qc_cols_nhood), title = "Cell quality per neighborhood group",
        subtitle = "Medians over each group's cells; a group with unusual values may reflect quality rather than state"
      )
      print(group_qc_gt)
      save_gt(group_qc_gt, "nhood_groups_qc.png", group_dir)

      p_group_qc <- plot_violin(
        data_obj[, !is.na(data_obj$nhood_group)], features = qc_cols_nhood, group_col = "nhood_group",
        group_colors = group_palette, group_order = group_levels, show_points = FALSE,
        title = "Cell quality per neighborhood group", subtitle = "Cell-level distributions"
      )
      print(p_group_qc)
      save_plot(p_group_qc, "nhood_groups_qc_violin.png", group_dir,
                width = max(8, 3 * length(qc_cols_nhood)), height = 5)
    }

    save_xlsx(c(Filter(Negate(is.null), list(
      nhoods = data.frame(Nhood = as.numeric(names(nhood_group)), group = nhood_group),
      resolution_stability = stability,
      group_qc = group_qc,
      composition_nhoods = rbind(comp_cl_nhood, comp_tp_nhood),
      composition_groups = rbind(comp_cl_group, comp_tp_group),
      sample_abundance = group_abundance)),
      setNames(group_da, paste0("groupDA_", names(group_da))),
      setNames(cluster_da, paste0("clusterDA_", names(cluster_da)))),
      "nhood_groups.xlsx", group_dir)
  }

  # --- programs/QC vs each other, across neighborhoods (primary) -------------------------
  if (is_primary && !is.null(nhood_qc) && length(score_cols) > 0) {
    nhood_table <- dplyr::left_join(nhood_qc, nhood_scores, by = "Nhood")
    p_qc_scores <- plot_nhood_correlations(
      nhood_table, qc_cols_nhood, score_cols,
      title = "QC metrics vs program scores across neighborhoods",
      subtitle = "Spearman rho of neighborhood means; a program tracking a QC metric may partly reflect cell quality"
    )
    print(p_qc_scores)
    save_plot(p_qc_scores, "nhood_qc_vs_programs.png", analysis_dir,
              width = max(6, 0.6 * length(score_cols) + 3), height = 0.5 * length(qc_cols_nhood) + 2.5)
  }

  # the primary checkpoint also keeps the neighborhood groups and each cell's group
  write_checkpoint(
    list(milo = milo_obj, results = analysis_results, analysis = analysis, params = params$da,
         nhood_group = if (is_primary && isTRUE(steps$nhood_groups)) nhood_group,
         cell_group  = if (is_primary && isTRUE(steps$nhood_groups)) cell_group),
    name = if (is_primary) "milo_da",
    file = if (!is_primary) paste0("04_milo_da_", analysis_name, ".rds"),
    dir = data_processed_dir,
    workflow = "04_neighborhood_analysis", script = "01_milo_da_de.R",
    settings = list(analysis          = analysis_name,
                    role              = analysis$role,
                    timepoints        = analysis$timepoints,
                    contrasts         = names(analysis$contrasts),
                    k                 = k_used,
                    d                 = d_used,
                    prop              = params$da$prop,
                    refinement_scheme = refinement_scheme,
                    seed              = params$da$seed,
                    paired_design     = params$da$paired_design,
                    n_cells           = ncol(milo_obj),
                    n_nhoods          = ncol(da_membership)),
    input_path = input_path
  )
}
if (!run_da) message("[skip] milo_da (disabled)")

# ------------------------------------------------------------------------------
# 3. Robustness: cell-level calls across analyses
# ------------------------------------------------------------------------------
# neighborhood IDs differ between analyses, cells do not: each cell's call
# (in an enriched / depleted / no significant neighborhood) is compared
# with the primary analysis, per contrast
if (run_da && length(analyses) > 1) {
  robust_dir <- file.path(fig_milo_da_dir, "robustness")
  dir.create(robust_dir, recursive = TRUE, showWarnings = FALSE)
  call_levels <- c("enriched", "depleted", "both", "not DA")
  robust_sheets <- list()

  for (contrast_name in names(contrasts)) {
    levels_used <- contrasts[[contrast_name]]
    runs <- names(analyses)[vapply(names(analyses), function(a) contrast_name %in% names(cell_calls[[a]]), logical(1))]
    if (length(runs) < 2) next
    ref_calls <- cell_calls[[primary_name]][[contrast_name]]

    agreement <- dplyr::bind_rows(lapply(setdiff(runs, primary_name), function(a) {
      calls <- cell_calls[[a]][[contrast_name]]
      shared <- intersect(names(ref_calls), names(calls))
      tab <- table(primary = factor(ref_calls[shared], call_levels), other = factor(calls[shared], call_levels))
      data.frame(analysis = a, cells_compared = length(shared),
                 primary_DA_cells = sum(ref_calls[shared] != "not DA"),
                 other_DA_cells = sum(calls[shared] != "not DA"),
                 same_call_pct = round(100 * mean(ref_calls[shared] == calls[shared]), 1),
                 primary_enriched_kept_pct = round(100 * tab["enriched", "enriched"] / max(1, sum(tab["enriched", ])), 1),
                 primary_depleted_kept_pct = round(100 * tab["depleted", "depleted"] / max(1, sum(tab["depleted", ])), 1),
                 opposite_calls = tab["enriched", "depleted"] + tab["depleted", "enriched"])
    }))
    agreement_gt <- gt_generic(
      agreement, title = paste0("Robustness of DA calls \u2014 ", levels_used[1], " vs ", levels_used[2]),
      subtitle = paste0("Per cell, compared with the primary analysis (", primary_name, "); ",
                        "kept_pct = % of the primary's enriched/depleted cells called the same")
    )
    print(agreement_gt)
    save_gt(agreement_gt, paste0("robustness_", contrast_name, ".png"), robust_dir)

    # cells colored by their call, one panel per analysis
    call_df <- dplyr::bind_rows(lapply(runs, function(a) {
      calls <- cell_calls[[a]][[contrast_name]]
      data.frame(cell = names(calls), call = calls, analysis = a)
    }))
    embed <- as.data.frame(Seurat::Embeddings(full_obj, "umap")[, 1:2])
    colnames(embed) <- c("Dim1", "Dim2")
    embed$cell <- rownames(embed)
    call_df <- dplyr::left_join(call_df, embed, by = "cell")
    call_df$call <- factor(ifelse(call_df$call == "enriched", paste("Enriched in", levels_used[1]),
                                  ifelse(call_df$call == "depleted", paste("Enriched in", levels_used[2]),
                                         ifelse(call_df$call == "both", "Both", "Not DA"))),
                           levels = c(paste("Enriched in", levels_used[2]), "Not DA", "Both", paste("Enriched in", levels_used[1])))
    call_df$analysis <- factor(call_df$analysis, levels = runs)
    call_df <- call_df[order(call_df$call == "Not DA", decreasing = TRUE), ]
    p_calls <- ggplot2::ggplot(call_df, ggplot2::aes(x = Dim1, y = Dim2, color = call)) +
      ggplot2::geom_point(size = 0.6) +
      ggplot2::scale_color_manual(values = setNames(c(lfc_pal$low, "#D9D9D9", "#7F7F7F", lfc_pal$high), levels(call_df$call)),
                                  drop = FALSE, name = NULL) +
      ggplot2::facet_wrap(~ analysis) +
      ggplot2::labs(title = paste0("DA calls per cell by analysis \u2014 ", levels_used[1], " vs ", levels_used[2]),
                    subtitle = "Each panel: the cells of that analysis, colored by the significant neighborhoods they fall in",
                    x = "umap_1", y = "umap_2") +
      ggplot2::theme_classic(base_size = 11) +
      ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(),
                     legend.position = "bottom") +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3)))
    print(p_calls)
    save_plot(p_calls, paste0("robustness_cells_", contrast_name, ".png"), robust_dir,
              width = min(15, 4.5 * length(runs)), height = 5.5)

    robust_sheets[[contrast_name]] <- agreement
  }
  if (length(robust_sheets) > 0) save_xlsx(robust_sheets, "robustness.xlsx", robust_dir)
}

# ------------------------------------------------------------------------------
# 4. miloDE: each primary contrast
# ------------------------------------------------------------------------------
if (run_de) {
  primary_obj <- subset(full_obj, cells = colnames(full_obj)[as.character(full_obj[[timepoint_col, drop = TRUE]]) %in% primary$timepoints])
  if (de_nhoods == "milode") sce_all <- methods::as(seurat_to_milo(primary_obj), "SingleCellExperiment")

  for (contrast_name in names(contrasts)) {
    levels_used <- contrasts[[contrast_name]]   # c(test, reference)
    out_dir <- file.path(fig_milo_de_dir, contrast_name)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    message("\n[DE] ", contrast_name)

    if (de_nhoods == "da") {
      # same neighborhoods as DA (same Nhood IDs); cells from other timepoints
      # are removed from every neighborhood, so each test sees only this
      # contrast's cells. "reference" < "test" alphabetically, so the last
      # model-matrix column (tested by default) is test vs reference.
      milo_de <- milo_da
      tp <- as.character(milo_de[[timepoint_col]])
      milo_de$milo_condition <- ifelse(tp == levels_used[1], "test",
                                       ifelse(tp == levels_used[2], "reference", "not_in_contrast"))
      membership <- miloR::nhoods(milo_de)
      membership[!tp %in% levels_used, ] <- 0
      miloR::nhoods(milo_de) <- Matrix::drop0(membership)
    } else {
      keep <- as.character(sce_all[[timepoint_col]]) %in% levels_used
      sce_contrast <- sce_all[, keep]
      # reference level first: the last model-matrix column (tested by default) is test vs reference
      sce_contrast$milo_condition <- factor(as.character(sce_contrast[[timepoint_col]]), levels = rev(levels_used))

      milo_de <- assign_milode_nhoods(
        sce_contrast, k = params$de$k, d = de_d, order = params$de$order,
        prop = params$de$prop, filtering = params$de$filtering,
        builder = params$de$builder,
        seed = params$de$seed
      )
    }
    milo_de <- miloR::buildNhoodGraph(milo_de)
    de_membership <- miloR::nhoods(milo_de)
    message("  ", ncol(de_membership), " neighborhoods (", de_nhoods, "), ",
            sum(Matrix::rowSums(de_membership) > 0), " cells from ", paste(levels_used, collapse = "/"))

    # paired (repeated-measures) contrasts: patient as a covariate, only
    # patients sampled at both timepoints
    paired_de <- isTRUE(params$de$paired_design) && all(levels_used %in% repeated_timepoints)
    if (paired_de) {
      patients_tp <- tapply(as.character(milo_de[[timepoint_col]]), as.character(milo_de[[patient_col]]),
                            function(x) sum(unique(x) %in% levels_used))
      complete <- names(patients_tp)[patients_tp == 2]
      # drop cells of unpaired patients from every neighborhood
      keep <- as.numeric(as.character(milo_de[[patient_col]]) %in% complete)
      membership <- Matrix::Diagonal(x = keep) %*% miloR::nhoods(milo_de)
      dimnames(membership) <- dimnames(miloR::nhoods(milo_de))
      miloR::nhoods(milo_de) <- Matrix::drop0(methods::as(membership, "CsparseMatrix"))
      milo_de$milo_pair <- factor(as.character(milo_de[[patient_col]]))
      message("  paired design: ", length(complete), " patients sampled at both timepoints")
    }
    # miloDE fails as a whole (rather than per neighborhood) when no
    # neighborhood can be tested - e.g. a paired design once samples with
    # too few cells are dropped leaves patients without a partner
    de_stat <- tryCatch(
      miloDE::de_test_neighbourhoods(
        milo_de, sample_id = "milo_sample",
        design = if (paired_de) ~ milo_pair + milo_condition else ~ milo_condition,
        covariates = if (paired_de) c("milo_pair", "milo_condition") else "milo_condition",
        min_n_cells_per_sample = params$de$min_n_cells_per_sample,
        min_count = params$de$min_count,
        output_type = "data.frame", verbose = FALSE
      ),
      error = function(e) {
        message("  miloDE failed: ", conditionMessage(e))
        NULL
      }
    )
    if (is.null(de_stat) || !is.data.frame(de_stat) || nrow(de_stat) == 0) {
      message("  [skip] ", contrast_name, ": no neighborhood could be tested",
              if (paired_de) paste0(" with the paired design. Set params$de$paired_design <- FALSE in ",
                                    "decisions.R to test it unpaired") else "", ".")
      next
    }

    nhood_summary <- milode_nhood_summary(de_stat, de_alpha)
    nhood_summary$dominant_cluster <- nhood_composition(
      de_membership, SummarizedExperiment::colData(milo_de),
      sample_col = "milo_sample", patient_col = patient_col, group_col = timepoint_col,
      cluster_col = "cluster_label"
    )$dominant_cluster[nhood_summary$Nhood]
    nhood_summary$size <- Matrix::colSums(de_membership)[nhood_summary$Nhood]

    # link to the DA results of the same contrast through shared cells
    if (de_nhoods == "da" && contrast_name %in% names(da_results_all)) {
      da_res <- da_results_all[[contrast_name]]
      nhood_summary <- dplyr::left_join(
        nhood_summary,
        da_res[, c("Nhood", "logFC", "SpatialFDR")] %>% dplyr::rename(DA_logFC = logFC, DA_SpatialFDR = SpatialFDR),
        by = "Nhood"
      )
    } else if (!is.null(milo_da) && contrast_name %in% names(da_results_all)) {
      da_res <- da_results_all[[contrast_name]]
      sig_up   <- da_res$Nhood[which(da_res$SpatialFDR < da_alpha & da_res$logFC > 0)]
      sig_down <- da_res$Nhood[which(da_res$SpatialFDR < da_alpha & da_res$logFC < 0)]
      da_membership <- miloR::nhoods(milo_da)
      nhood_summary$da_overlap_enriched <- round(nhood_overlap_fraction(de_membership, da_membership, sig_up)[nhood_summary$Nhood], 2)
      nhood_summary$da_overlap_depleted <- round(nhood_overlap_fraction(de_membership, da_membership, sig_down)[nhood_summary$Nhood], 2)
    }

    # programs per DE neighborhood
    if (run_scores) {
      de_scores <- nhood_mean_scores(de_membership, score_values)
      nhood_summary <- dplyr::left_join(nhood_summary, de_scores[, c("Nhood", score_cols)], by = "Nhood")
      if (nrow(nhood_summary) >= 2) {
        row_labels <- setNames(paste0(nhood_summary$dominant_cluster, ", ", nhood_summary$n_DE_genes, " DE genes"),
                               as.character(nhood_summary$Nhood))
        p_de_programs <- plot_nhood_program_heatmap(
          nhood_summary, score_cols, row_label = row_labels,
          title = paste0("Programs in miloDE neighborhoods \u2014 ", levels_used[1], " vs ", levels_used[2]),
          subtitle = "Mean score per neighborhood (both timepoints' cells), z-scored within program"
        )
        print(p_de_programs)
        save_plot(p_de_programs, paste0("milode_program_heatmap_", contrast_name, ".png"), out_dir,
                  width = max(7, 0.6 * length(score_cols) + 4), height = max(4, 0.3 * nrow(nhood_summary) + 2))
      }
    }

    gene_summary <- milode_gene_summary(de_stat, de_alpha, extreme_logfc = params$de$extreme_logfc)
    message("  ", sum(nhood_summary$test_performed), " of ", nrow(nhood_summary), " neighborhoods tested; ",
            nrow(gene_summary), " genes DE in >= 1 neighborhood")

    nhood_gt <- gt_generic(
      nhood_summary,
      title = paste0("miloDE neighborhoods \u2014 ", levels_used[1], " vs ", levels_used[2]),
      subtitle = paste0("n_DE_genes: within-neighborhood BH < ", de_alpha,
                        if (de_nhoods == "da") "; DA_logFC/DA_SpatialFDR: same neighborhood's DA result"
                        else "; da_overlap: fraction of cells in DA-significant neighborhoods of this contrast")
    )
    print(nhood_gt)
    save_gt(nhood_gt, paste0("milode_nhood_summary_", contrast_name, ".png"), out_dir)

    top_genes <- utils::head(gene_summary, params$de$n_top_genes)
    if (nrow(top_genes) > 0) {
      genes_gt <- gt_generic(
        top_genes,
        title = paste0("Genes DE in the most neighborhoods \u2014 ", levels_used[1], " vs ", levels_used[2]),
        subtitle = paste0("Top ", params$de$n_top_genes, "; DE = within-neighborhood BH < ", de_alpha, "; logFC > 0 = higher in ",
                          levels_used[1], "; pct_extreme_logFC = % of DE calls with |logFC| > ",
                          params$de$extreme_logfc, " (on/off)")
      )
      print(genes_gt)
      save_gt(genes_gt, paste0("milode_top_genes_", contrast_name, ".png"), out_dir)
    }

    p_n_de <- miloDE::plot_milo_by_single_metric(milo_de, nhood_summary, colour_by = "n_DE_genes",
                                                 layout = "UMAP", size_range = c(1.5, 4)) +
      ggplot2::scale_fill_viridis_c(name = "n DE genes") +
      ggplot2::labs(title = paste0("DE genes per neighborhood \u2014 ", levels_used[1], " vs ", levels_used[2]),
                    subtitle = paste0("Within-neighborhood BH < ", de_alpha))
    print(p_n_de)
    save_plot(p_n_de, paste0("milode_n_de_genes_", contrast_name, ".png"), out_dir, width = 8, height = 7)

    for (gene in intersect(params$de$genes_to_plot, unique(de_stat$gene))) {
      p_gene <- miloDE::plot_DE_single_gene(milo_de, de_stat, gene = gene, layout = "UMAP", alpha = de_alpha) +
        ggplot2::labs(title = paste0(gene, " \u2014 ", levels_used[1], " vs ", levels_used[2]),
                      subtitle = paste0("Colored if spatially corrected p < ", de_alpha))
      print(p_gene)
      save_plot(p_gene, paste0("milode_gene_", gene, "_", contrast_name, ".png"), out_dir, width = 8, height = 7)
    }
    missing_genes <- setdiff(params$de$genes_to_plot, unique(de_stat$gene))
    if (length(missing_genes) > 0) message("  genes_to_plot not tested: ", paste(missing_genes, collapse = ", "))

    de_significant <- de_stat[which(de_stat$pval_corrected_across_genes < de_alpha), ]
    save_xlsx(list(nhood_summary = nhood_summary, gene_summary = gene_summary,
                   significant_tests = de_significant),
              paste0("milode_", contrast_name, ".xlsx"), out_dir)

    gz <- gzfile(file.path(out_dir, paste0("milode_all_tests_", contrast_name, ".csv.gz")), "w")
    utils::write.csv(de_stat, gz, row.names = FALSE)
    close(gz)

    write_checkpoint(
      list(milo = milo_de, de_stat = de_stat, nhood_summary = nhood_summary, params = params$de),
      file = paste0("04_milo_de_", contrast_name, ".rds"),
      dir = data_processed_dir,
      workflow = "04_neighborhood_analysis", script = "01_milo_da_de.R",
      settings = list(contrast        = contrast_name,
                      levels          = levels_used,
                      neighborhoods  = de_nhoods,
                      k               = if (de_nhoods == "milode") params$de$k else primary$k %||% params$da$k,
                      d               = if (de_nhoods == "milode") de_d else primary$d %||% params$da$d,
                      paired_design   = paired_de,
                      alpha           = de_alpha,
                      n_cells         = sum(Matrix::rowSums(miloR::nhoods(milo_de)) > 0),
                      n_nhoods        = nrow(nhood_summary),
                      n_nhoods_tested = sum(nhood_summary$test_performed)),
      input_path = input_path
    )
  }
} else {
  message("[skip] milo_de (disabled)")
}

message("\nDifferential abundance finished; analyses: ", paste(names(analyses), collapse = ", "))
