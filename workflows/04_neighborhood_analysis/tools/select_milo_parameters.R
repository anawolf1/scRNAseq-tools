# ------------------------------------------------------------------------------
# select_milo_parameters.R (workflows/04_neighborhood_analysis/tools)
#
# Chooses neighborhood parameters for 01_milo_da_de.R, for one analysis in
# config/decisions.R (sweep_analysis; normally the primary). It records
# nothing: put the chosen values in config/decisions.R.
#
#   Part A - DA neighborhoods (params$da$k, params$da$d, params$da$prop):
#     k x d (x prop) sweep. Target: a typical neighborhood holds >=
#     params$da$size_multiplier x the samples compared, counting only cells
#     from those samples, while staying resolved (cluster purity) and not
#     dominated by one patient.
#   Part B - DE neighborhoods (params$de$k): k sweep per contrast, for
#     params$de$neighborhoods = "milode"; with "da" (DE on the DA
#     neighborhoods) it reports their DE testability instead.
#   Part C - neighborhood groups (params$da$nhood_groups$resolution):
#     subsampling stability of the Leiden grouping across resolutions.
#
# No part uses DA or DE results, so the choice cannot favor significant
# findings.
#
# Output: results/04_neighborhood_analysis/exploration/{milo_da_parameters,
#         milo_de_parameters, nhood_group_resolution}/
#
# Usage:
#   source(here::here("workflows", "04_neighborhood_analysis", "tools", "select_milo_parameters.R"))
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Setup (edit per use)
# ------------------------------------------------------------------------------
source(here::here("workflows", "04_neighborhood_analysis", "setup.R"))
decision_file <- file.path(workflow_dir, "config", "decisions.R")
if (!file.exists(decision_file)) {
  stop("No workflows/04_neighborhood_analysis/config/decisions.R found. Copy config/decisions_example.R ",
       "to config/decisions.R and set the analyses first.")
}
source(decision_file)

run_part_a <- TRUE    # DA neighborhoods
run_part_b <- TRUE    # DE neighborhoods
run_part_c <- FALSE   # Leiden resolution of the neighborhood groups (needs params$da$k/d set)

# which analysis in config/decisions.R to screen (its timepoints and
# contrasts); usually "primary" - other analyses then reuse its k/d/prop
sweep_analysis <- "primary"

# Part A candidates. d: PCs of the clustering's PCA - keep the clustering's
# dims unless there is a reason independent of DA results to change it
k_candidates    <- c(20, 25, 30, 35, 40, 50)
d_candidates    <- c(12)
prop_candidates <- c(0.1)   # e.g. c(0.1, 0.2) to check it; one heatmap set per value
high_purity     <- 0.8      # a neighborhood is dominated by one patient above this fraction

# Part B candidates: second-order miloDE neighborhoods grow roughly with
# k^2, so k stays small; first-order miloR neighborhoods need a larger k
de_k_candidates <- if (params$de$builder == "miloR") c(40, 60, 80, 100, 120) else c(5, 8, 10, 12, 15)

# Part C: resolutions screened, subsamples, and the rule fixed beforehand:
# the highest mean ARI; among resolutions within ari_tolerance of it, the
# highest (finest grouping)
group_resolutions <- seq(0.2, 2, by = 0.1)
group_n_boot      <- 30
group_frac        <- 0.8
ari_tolerance     <- 0.05

# Part B: d to use if neither params$de$d nor params$da$d is set yet
de_d_if_unset <- NULL

patient_col   <- params$study_design$patient_col
timepoint_col <- params$study_design$timepoint_col

analyses <- resolve_milo_analyses(params$analyses)
if (!sweep_analysis %in% names(analyses)) {
  stop("sweep_analysis '", sweep_analysis, "' is not one of: ", paste(names(analyses), collapse = ", "))
}
target <- analyses[[sweep_analysis]]
prepared <- prepare_milo_input(input_path, target$timepoints, target$contrasts, patient_col, timepoint_col)
data_obj  <- prepared$obj
contrasts <- prepared$contrasts

# ------------------------------------------------------------------------------
# Part A. DA neighborhoods: k x d sweep
# ------------------------------------------------------------------------------
if (run_part_a) {
  da_dir <- file.path(fig_milo_sweeps_dir, "milo_da_parameters", sweep_analysis)
  dir.create(da_dir, recursive = TRUE, showWarnings = FALSE)

  sample_gt <- gt_generic(prepared$sample_table, title = "Samples and cells per timepoint",
                          subtitle = paste0("Input: ", basename(input_path), " (", ncol(data_obj), " cells)"))
  print(sample_gt)
  save_gt(sample_gt, "sample_table.png", da_dir)

  if (length(Seurat::Stdev(data_obj[["pca"]])) > 0) {
    p_elbow <- Seurat::ElbowPlot(data_obj, ndims = min(50, length(Seurat::Stdev(data_obj[["pca"]]))))
    print(p_elbow)
    save_plot(p_elbow, "pca_elbow.png", da_dir, width = 7, height = 5)
  }
  input_manifest <- read_manifest(input_path)
  if (!is.null(input_manifest$setting.dims)) message("The input clustering used dims = ", input_manifest$setting.dims)

  milo_base <- seurat_to_milo(data_obj)

  # --- candidates ------------------------------------------------------------
  refinement_scheme <- params$da$refinement_scheme
  seed              <- params$da$seed

  n_pcs <- ncol(SingleCellExperiment::reducedDim(milo_base, "PCA"))
  d_candidates <- d_candidates[d_candidates <= n_pcs]
  if (length(d_candidates) == 0) stop("No d candidate fits within the ", n_pcs, " stored PCs.")

  # --- sweep -----------------------------------------------------------------
  sweep <- evaluate_milo_da_grid(
    milo_base,
    k_grid = k_candidates, d_grid = d_candidates, prop_grid = prop_candidates,
    refinement_scheme = refinement_scheme,
    sample_col = "milo_sample", patient_col = patient_col, group_col = timepoint_col,
    cluster_col = "cluster_label", contrasts = contrasts,
    size_multiplier = params$da$size_multiplier, high_purity = high_purity, seed = seed
  )
  summary_df <- sweep$summary %>% dplyr::arrange(prop, k, d)

  if (nrow(summary_df) == 0) stop("Every combination failed - see the messages above.")

  # --- tables + plots --------------------------------------------------------
  summary_gt <- gt_generic(
    summary_df,
    title = paste("Milo DA neighborhoods \u2014 k x d x prop sweep:", sweep_analysis),
    subtitle = paste0("refinement = ", refinement_scheme, "; targets = ", params$da$size_multiplier,
                      " x samples (all cells, and per contrast)")
  )
  print(summary_gt)
  save_gt(summary_gt, "milo_da_sweep_summary.png", da_dir)
  save_xlsx(list(sweep_summary = summary_df, sample_table = prepared$sample_table),
            "milo_da_sweep_summary.xlsx", da_dir)

  # every heatmap in the "sweep_heatmap" palette; percentages and shares on fixed scales
  heatmaps <- list(
    list(col = "median_size", title = "Median neighborhood size (all cells)"),
    list(col = "pct_ge_target", limits = c(0, 100),
         title = paste0("% neighborhoods >= ", params$da$size_multiplier, " x all samples (",
                        summary_df$target_size[1], " cells)")),
    list(col = "pct_high_patient_purity", limits = c(0, 100),
         title = paste0("% neighborhoods > ", 100 * high_purity, "% from one patient")),
    list(col = "median_cluster_purity", limits = c(0, 1),
         title = "Median cluster purity (resolution)")
  )
  for (contrast_name in names(contrasts)) {
    heatmaps[[length(heatmaps) + 1]] <- list(
      col = paste0(contrast_name, "_pct_ge_target"), limits = c(0, 100),
      title = paste0(contrast_name, ": % neighborhoods >= ", summary_df[[paste0(contrast_name, "_target")]][1],
                     " cells from its samples"))
    heatmaps[[length(heatmaps) + 1]] <- list(
      col = paste0(contrast_name, "_pct_testable"), limits = c(0, 100),
      title = paste0(contrast_name, ": % neighborhoods with >= 2 samples per group"))
  }

  for (prop_value in prop_candidates) {
    prop_df <- summary_df[summary_df$prop == prop_value, ]
    prop_suffix <- if (length(prop_candidates) > 1) paste0("_prop", prop_value) else ""
    for (h in heatmaps) {
      if (!h$col %in% names(prop_df) || nrow(prop_df) == 0) next
      p <- plot_sweep_tile(prop_df, h$col, x_col = "d", y_col = "k", title = h$title,
                           subtitle = paste0("prop = ", prop_value), limits = h$limits)
      print(p)
      save_plot(p, paste0("heatmap_", h$col, prop_suffix, ".png"), da_dir, width = 7, height = 5)
    }
  }

  p_sizes <- plot_nhood_size_distributions(
    sweep$sizes, "size", target = summary_df$target_size[1],
    title = "Neighborhood size, all cells",
    subtitle = paste0("Dashed line: ", params$da$size_multiplier, " x all samples")
  )
  print(p_sizes)
  save_plot(p_sizes, "size_histograms_all_cells.png", da_dir, width = 16, height = 12)

  for (contrast_name in names(contrasts)) {
    p_contrast <- plot_nhood_size_distributions(
      sweep$sizes, paste0(contrast_name, "_cells"),
      target = summary_df[[paste0(contrast_name, "_target")]][1],
      title = paste0("Cells from ", contrast_name, " samples per neighborhood"),
      subtitle = paste0("Dashed line: ", params$da$size_multiplier, " x that contrast's samples")
    )
    print(p_contrast)
    save_plot(p_contrast, paste0("size_histograms_", contrast_name, ".png"), da_dir, width = 16, height = 12)
  }

  message(
    "\nRecord in config/decisions.R:\n",
    "  params$da$k <- <k>\n  params$da$d <- <d>\n  params$da$prop <- <prop>\n",
    "  params$da$sweep <- list(analysis = \"", sweep_analysis, "\", k_tested = c(", paste(k_candidates, collapse = ", "), "), ",
    "d_tested = c(", paste(d_candidates, collapse = ", "), "), prop_tested = c(",
    paste(prop_candidates, collapse = ", "), "))"
  )
}

# ------------------------------------------------------------------------------
# Part B. DE neighborhoods
# ------------------------------------------------------------------------------
if (run_part_b) {
  de_dir <- file.path(fig_milo_sweeps_dir, "milo_de_parameters", sweep_analysis)
  dir.create(de_dir, recursive = TRUE, showWarnings = FALSE)
  sce_base <- seurat_to_milo(data_obj)

  if (params$de$neighborhoods == "da") {
    if (is.null(params$da$k) || is.null(params$da$d)) {
      stop("Set params$da$k and params$da$d first (Part A of this tool).")
    }
    milo_da <- build_milo_nhoods(sce_base, k = params$da$k, d = params$da$d,
                                 prop = params$da$prop,
                                 refinement_scheme = params$da$refinement_scheme,
                                 seed = params$da$seed)
    meta <- as.data.frame(SummarizedExperiment::colData(milo_da))

    testability <- dplyr::bind_rows(lapply(c(1, 3, 5, 10), function(min_cells) {
      composition <- nhood_composition(miloR::nhoods(milo_da), meta, "milo_sample", patient_col,
                                       timepoint_col, contrasts = contrasts, min_cells_per_sample = min_cells)
      dplyr::bind_rows(lapply(names(contrasts), function(contrast_name) {
        n_test <- composition[[paste0(contrast_name, "_n_test")]]
        n_ref  <- composition[[paste0(contrast_name, "_n_ref")]]
        data.frame(contrast = contrast_name, min_n_cells_per_sample = min_cells,
                   n_nhoods = nrow(composition),
                   median_samples_passing = stats::median(n_test + n_ref),
                   pct_testable = round(100 * mean(n_test >= 2 & n_ref >= 2), 1))
      }))
    }))

    testability_gt <- gt_generic(
      testability,
      title = "miloDE on the DA neighborhoods \u2014 testability",
      subtitle = paste0("k = ", params$da$k, ", d = ", params$da$d,
                        "; testable = >= 2 samples per group with >= min_n_cells_per_sample cells")
    )
    print(testability_gt)
    save_gt(testability_gt, "milode_on_da_testability.png", de_dir)
    save_xlsx(list(testability = testability), "milode_on_da_testability.xlsx", de_dir)

    message("\nparams$de$neighborhoods = \"da\": no k sweep needed. Use the table to choose ",
            "params$de$min_n_cells_per_sample (currently ", params$de$min_n_cells_per_sample, ").")
  } else {
    builder <- params$de$builder
    d <- params$de$d %||% params$da$d %||% de_d_if_unset
    if (is.null(d)) stop("Set params$da$d (or params$de$d) in config/decisions.R, or de_d_if_unset above.")

    min_n_cells_per_sample <- params$de$min_n_cells_per_sample

    # --- sweep ---------------------------------------------------------------------
    sweep <- evaluate_milode_grid(
      sce_base,
      k_grid = de_k_candidates, d = d,
      order = params$de$order, prop = params$de$prop,
      filtering = params$de$filtering,
      sample_col = "milo_sample", patient_col = patient_col, group_col = timepoint_col,
      cluster_col = "cluster_label", contrasts = contrasts,
      min_n_cells_per_sample = min_n_cells_per_sample,
      builder = builder,
      seed = params$de$seed
    )
    summary_df <- sweep$summary
    if (nrow(summary_df) == 0) stop("Every combination failed - see the messages above.")

    # --- tables and plots ------------------------------------------------------------
    summary_gt <- gt_generic(
      summary_df,
      title = "miloDE neighborhoods \u2014 k sweep per contrast",
      subtitle = paste0("builder = ", builder, ", d = ", d,
                        if (builder == "miloDE") paste0(", order = ", params$de$order) else "",
                        ", prop = ", params$de$prop,
                        "; a sample passes with >= ", min_n_cells_per_sample, " cells in the neighborhood")
    )
    print(summary_gt)
    save_gt(summary_gt, "milode_sweep_summary.png", de_dir)
    save_xlsx(list(sweep_summary = summary_df), "milode_sweep_summary.xlsx", de_dir)

    metric_long <- summary_df %>%
      dplyr::select(contrast, k, dplyr::any_of(c("pct_testable", "pct_all_samples_pass", "median_cluster_purity"))) %>%
      dplyr::mutate(dplyr::across(dplyr::any_of("median_cluster_purity"), ~ .x * 100)) %>%
      tidyr::pivot_longer(-c(contrast, k), names_to = "metric", values_to = "value")

    p_metrics <- ggplot2::ggplot(metric_long, ggplot2::aes(x = k, y = value, color = contrast)) +
      ggplot2::geom_line() +
      ggplot2::geom_point() +
      ggplot2::facet_wrap(~ metric, nrow = 1) +
      ggplot2::labs(title = "miloDE neighborhood testability and resolution by k",
                    subtitle = "median_cluster_purity shown x 100", y = "%") +
      ggplot2::theme_classic(base_size = 12)
    print(p_metrics)
    save_plot(p_metrics, "milode_metrics_by_k.png", de_dir, width = 12, height = 4.5)

    for (contrast_name in names(contrasts)) {
      p_sizes <- plot_nhood_size_distributions(
        sweep$sizes[sweep$sizes$contrast == contrast_name, ], "size",
        title = paste0("miloDE neighborhood size \u2014 ", contrast_name)
      )
      print(p_sizes)
      save_plot(p_sizes, paste0("milode_size_histograms_", contrast_name, ".png"), de_dir, width = 12, height = 7)
    }
    message("\nRecord the chosen k in config/decisions.R: params$de$k <- <k>")
  }
}

# ------------------------------------------------------------------------------
# Part C. Neighborhood groups: Leiden resolution
# ------------------------------------------------------------------------------
if (run_part_c) {
  if (is.null(params$da$k) || is.null(params$da$d)) stop("Set params$da$k and params$da$d first (Part A).")
  group_dir <- file.path(fig_milo_sweeps_dir, "nhood_group_resolution", sweep_analysis)
  dir.create(group_dir, recursive = TRUE, showWarnings = FALSE)

  # the same neighborhoods the pipeline builds for this analysis
  milo_groups <- build_milo_nhoods(seurat_to_milo(data_obj), k = params$da$k, d = params$da$d,
                                   prop = params$da$prop,
                                   refinement_scheme = params$da$refinement_scheme,
                                   seed = params$da$seed)
  sweep_groups <- group_resolution_sweep(miloR::nhoods(milo_groups), group_resolutions,
                                         n_boot = group_n_boot, frac = group_frac,
                                         seed = params$da$nhood_groups$seed)
  best <- max(sweep_groups$mean_ari)
  chosen <- max(sweep_groups$resolution[sweep_groups$mean_ari >= best - ari_tolerance])
  sweep_groups$chosen_by_rule <- sweep_groups$resolution == chosen
  print(sweep_groups)

  sweep_gt <- gt_generic(sweep_groups, title = paste("Neighborhood groups: Leiden resolution \u2014", sweep_analysis),
                         subtitle = paste0("Mean ARI of the grouping when ", 100 * (1 - group_frac),
                                           "% of neighborhoods are left out (", group_n_boot, " subsamples); ",
                                           "rule: highest mean ARI, finest within ", ari_tolerance, " of it"))
  print(sweep_gt)
  save_gt(sweep_gt, "group_resolution_sweep.png", group_dir)
  save_xlsx(list(resolution_sweep = sweep_groups), "group_resolution_sweep.xlsx", group_dir)

  p_res <- ggplot2::ggplot(sweep_groups, ggplot2::aes(x = resolution, y = mean_ari)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = mean_ari - sd_ari, ymax = mean_ari + sd_ari), fill = "#D9D9D9") +
    ggplot2::geom_line() +
    ggplot2::geom_point(ggplot2::aes(color = chosen_by_rule), size = 2.5) +
    ggplot2::geom_text(ggplot2::aes(label = n_groups), vjust = -1, size = 3) +
    ggplot2::scale_color_manual(values = c(`FALSE` = "#000000", `TRUE` = "#B22222"), guide = "none") +
    ggplot2::labs(title = "Stability of the neighborhood grouping by Leiden resolution",
                  subtitle = "Mean ARI \u00b1 SD over subsamples; numbers = groups; red = chosen by the rule",
                  x = "Leiden resolution", y = "Mean ARI (stability)") +
    ggplot2::theme_classic(base_size = 12)
  print(p_res)
  save_plot(p_res, "group_resolution_sweep.png", group_dir, width = 8, height = 5)

  message("\nRecord in config/decisions.R:\n  params$da$nhood_groups$resolution <- ", chosen,
          "\n  params$da$nhood_groups$sweep <- list(resolutions_tested = seq(", min(group_resolutions), ", ",
          max(group_resolutions), ", by = ", diff(group_resolutions[1:2]), "), n_boot = ", group_n_boot,
          ", frac = ", group_frac, ")",
          "\n  (rule: highest mean ARI, finest within ", ari_tolerance, " of it)")
}
