# ------------------------------------------------------------------------------
# evaluate_milo.R
#
# Building Milo objects and evaluating neighborhood-construction parameters
# (k, d, prop) before differential-abundance (miloR) and neighborhood-level
# differential-expression (miloDE) testing - the Milo counterpart of
# clustering_evaluation.R.
#
# Every metric here describes neighborhood size, sample representation and
# composition only. None uses test results, so parameters can be chosen
# without looking at which neighborhoods come out significant.
#
# Functions
# seurat_to_milo()                 Seurat object -> Milo object with PCA/UMAP
# make_sample_id()                 combine metadata columns into one replicate ID
# build_milo_nhoods()              buildGraph() + makeNhoods() with matching k/d
# milo_fdr_weighting()             testNhoods() fdr.weighting matching a refinement scheme
# nhood_group_counts()             cells per neighborhood x level of a metadata column
# nhood_composition()              per-neighborhood size, sample, patient, cluster, contrast metrics
# summarize_nhood_composition()    one-row summary of nhood_composition() output
# evaluate_milo_da_grid()          k x d (x prop) sweep for miloR DA neighborhoods
# filter_nhoods_greedy()           drop redundant neighborhoods (greedy set cover, pure R)
# assign_milode_nhoods()           miloDE::assign_neighbourhoods() with R-side filtering
# evaluate_milode_grid()           k sweep for miloDE neighborhoods, per contrast
# plot_nhood_size_distributions()  neighborhood-size histograms per parameter combination
# ------------------------------------------------------------------------------

#' Convert a Seurat object to a Milo object, keeping PCA and UMAP
#'
#' Seurat::as.SingleCellExperiment() carries reductions over with upper-case
#' names ("PCA", "UMAP"); this checks they arrived, so downstream calls can
#' rely on reduced.dim = "PCA" and layout = "UMAP".
#'
#' @param obj Seurat object with "pca" (required) and "umap" (optional)
#'   reductions.
#' @param require_umap Logical. Stop if no UMAP is present. Default TRUE,
#'   since neighborhood graph plots use it as their layout.
#'
#' @return A Milo object.
seurat_to_milo <- function(obj,
                           require_umap = TRUE) {

  if (!"pca" %in% Seurat::Reductions(obj)) {
    stop("seurat_to_milo(): no 'pca' reduction in the object.", call. = FALSE)
  }
  if (require_umap && !"umap" %in% Seurat::Reductions(obj)) {
    stop("seurat_to_milo(): no 'umap' reduction in the object (needed as the plotting layout).", call. = FALSE)
  }

  sce <- Seurat::as.SingleCellExperiment(obj)

  if (!"PCA" %in% SingleCellExperiment::reducedDimNames(sce)) {
    stop("seurat_to_milo(): PCA was not carried over by as.SingleCellExperiment(). ",
         "Found: ", paste(SingleCellExperiment::reducedDimNames(sce), collapse = ", "), call. = FALSE)
  }

  miloR::Milo(sce)
}
#' Combine metadata columns into one replicate (sample) ID
#'
#' @param meta Data.frame of cell metadata (colData or Seurat meta.data).
#' @param sample_vars Character vector of column names, e.g.
#'   c("Patient", "Timepoint").
#' @param sep Separator. Default "_".
#'
#' @return Character vector, one sample ID per cell.
make_sample_id <- function(meta,
                           sample_vars,
                           sep = "_") {

  meta <- as.data.frame(meta)
  validate_required_cols(meta, sample_vars, "make_sample_id")
  do.call(paste, c(lapply(sample_vars, function(v) as.character(meta[[v]])), sep = sep))
}
#' Build the kNN graph and neighborhoods with one consistent k and d
#'
#' buildGraph() and makeNhoods() both take k and d; giving them different
#' values builds neighborhoods on a different geometry from the graph.
#' This sets both from the same arguments. Index-cell sampling in
#' makeNhoods() is random, so the seed is set immediately before it.
#'
#' @param milo Milo object with a "PCA" reducedDim.
#' @param k Number of nearest neighbors.
#' @param d Number of PCA dimensions (the first d columns of the stored PCA).
#' @param prop Proportion of cells sampled as neighborhood index cells.
#' @param refinement_scheme "graph" or "reduced_dim", passed to makeNhoods().
#'   Use milo_fdr_weighting() to pick the matching testNhoods() weighting.
#' @param seed Seed set immediately before makeNhoods(). Default 123.
#'
#' @return The Milo object with graph and neighborhoods.
build_milo_nhoods <- function(milo,
                              k,
                              d,
                              prop = 0.1,
                              refinement_scheme = "graph",
                              seed = 123) {

  n_pcs <- ncol(SingleCellExperiment::reducedDim(milo, "PCA"))
  if (d > n_pcs) {
    stop("build_milo_nhoods(): d = ", d, " but the PCA only has ", n_pcs, " dimensions.", call. = FALSE)
  }
  if (!refinement_scheme %in% c("graph", "reduced_dim")) {
    stop("build_milo_nhoods(): refinement_scheme must be \"graph\" or \"reduced_dim\".", call. = FALSE)
  }

  milo <- miloR::buildGraph(milo, k = k, d = d, reduced.dim = "PCA")
  set.seed(seed)
  miloR::makeNhoods(milo, prop = prop, k = k, d = d, refined = TRUE,
                    reduced_dims = "PCA", refinement_scheme = refinement_scheme)
}
#' Spatial-FDR weighting that matches a neighborhood refinement scheme
#'
#' testNhoods()' default "k-distance" weighting relies on distances that a
#' graph-refined neighborhood set does not carry; miloR pairs
#' refinement_scheme = "graph" with fdr.weighting = "graph-overlap".
#'
#' @param refinement_scheme "graph" or "reduced_dim".
#'
#' @return "graph-overlap" or "k-distance".
milo_fdr_weighting <- function(refinement_scheme) {
  switch(refinement_scheme,
         graph       = "graph-overlap",
         reduced_dim = "k-distance",
         stop("milo_fdr_weighting(): unknown refinement_scheme '", refinement_scheme, "'.", call. = FALSE))
}
#' Cells per neighborhood for each level of a metadata column
#'
#' @param membership Cell x neighborhood 0/1 matrix, e.g. miloR::nhoods(milo).
#' @param values Vector (one value per cell, same order as rownames of
#'   `membership`) to tabulate.
#'
#' @return Numeric matrix, neighborhoods x levels of `values`.
nhood_group_counts <- function(membership,
                               values) {

  if (length(values) != nrow(membership)) {
    stop("nhood_group_counts(): `values` must have one entry per cell (row of `membership`).", call. = FALSE)
  }
  indicator <- Matrix::t(Matrix::fac2sparse(factor(as.character(values))))
  as.matrix(Matrix::crossprod(membership, indicator))
}
#' Per-neighborhood composition metrics
#'
#' @param membership Cell x neighborhood 0/1 matrix, e.g. miloR::nhoods(milo).
#' @param meta Cell metadata, rows in the same order as `membership`.
#' @param sample_col Column holding the replicate ID (one per Patient x
#'   Timepoint), e.g. from make_sample_id().
#' @param patient_col Column holding the patient ID.
#' @param group_col Column holding the condition tested (e.g. Timepoint).
#' @param cluster_col Optional column used for a cluster-purity metric
#'   (resolution: larger neighborhoods mix more clusters). NULL to skip.
#' @param contrasts Named list, each c(test_level, reference_level) of
#'   `group_col`. NULL to skip contrast metrics.
#' @param min_cells_per_sample A sample counts as represented in a
#'   neighborhood when it contributes at least this many cells. Default 1.
#'
#' @return A data.frame, one row per neighborhood: Nhood, size, n_samples,
#'   n_patients, patient_purity, (dominant_cluster, cluster_purity), and per
#'   contrast <name>_cells (cells from the contrast's samples),
#'   <name>_n_test, <name>_n_ref (represented samples per group).
nhood_composition <- function(membership,
                              meta,
                              sample_col,
                              patient_col,
                              group_col,
                              cluster_col = NULL,
                              contrasts = NULL,
                              min_cells_per_sample = 1) {

  meta <- as.data.frame(meta)
  validate_required_cols(meta, c(sample_col, patient_col, group_col, cluster_col), "nhood_composition")

  size <- Matrix::colSums(membership)
  sample_counts  <- nhood_group_counts(membership, meta[[sample_col]])
  patient_counts <- nhood_group_counts(membership, meta[[patient_col]])

  out <- data.frame(
    Nhood          = seq_len(ncol(membership)),
    size           = as.numeric(size),
    n_samples      = rowSums(sample_counts >= min_cells_per_sample),
    n_patients     = rowSums(patient_counts > 0),
    patient_purity = apply(patient_counts, 1, max) / size
  )

  if (!is.null(cluster_col)) {
    cluster_counts <- nhood_group_counts(membership, meta[[cluster_col]])
    out$dominant_cluster <- colnames(cluster_counts)[max.col(cluster_counts, ties.method = "first")]
    out$cluster_purity   <- apply(cluster_counts, 1, max) / size
  }

  if (!is.null(contrasts)) {
    sample_group <- unique(meta[, c(sample_col, group_col)])
    if (anyDuplicated(sample_group[[sample_col]])) {
      stop("nhood_composition(): some samples map to more than one ", group_col, " value.", call. = FALSE)
    }
    group_of <- setNames(as.character(sample_group[[group_col]]), sample_group[[sample_col]])

    for (contrast_name in names(contrasts)) {
      levels_used <- contrasts[[contrast_name]]
      test_samples <- intersect(colnames(sample_counts), names(group_of)[group_of == levels_used[1]])
      ref_samples  <- intersect(colnames(sample_counts), names(group_of)[group_of == levels_used[2]])

      out[[paste0(contrast_name, "_cells")]]  <- rowSums(sample_counts[, c(test_samples, ref_samples), drop = FALSE])
      out[[paste0(contrast_name, "_n_test")]] <- rowSums(sample_counts[, test_samples, drop = FALSE] >= min_cells_per_sample)
      out[[paste0(contrast_name, "_n_ref")]]  <- rowSums(sample_counts[, ref_samples, drop = FALSE] >= min_cells_per_sample)
    }
  }

  out
}
#' One-row summary of per-neighborhood composition metrics
#'
#' Size targets follow the miloR rule of thumb that the typical
#' neighborhood should hold at least size_multiplier x the number of
#' samples being compared: the whole dataset for `size`, and each
#' contrast's own samples for its `<name>_cells`.
#'
#' @param composition Output of nhood_composition().
#' @param n_samples_total Number of samples in the whole object.
#' @param contrast_n_samples Named integer vector, samples per contrast
#'   (names matching the contrasts used in nhood_composition()).
#' @param pct_cells_covered Percentage of cells in at least one
#'   neighborhood.
#' @param size_multiplier Default 5.
#' @param high_purity Patient-purity threshold above which a neighborhood
#'   counts as dominated by one patient. Default 0.8.
#' @param min_samples_per_group A contrast neighborhood counts as testable
#'   when each group has at least this many represented samples. Default 2.
#'
#' @return A one-row data.frame.
summarize_nhood_composition <- function(composition,
                                        n_samples_total,
                                        contrast_n_samples = NULL,
                                        pct_cells_covered = NA_real_,
                                        size_multiplier = 5,
                                        high_purity = 0.8,
                                        min_samples_per_group = 2) {

  target_all <- size_multiplier * n_samples_total

  out <- data.frame(
    n_nhoods          = nrow(composition),
    pct_cells_covered = round(pct_cells_covered, 1),
    mean_size         = round(mean(composition$size), 1),
    median_size       = stats::median(composition$size),
    target_size       = target_all,
    pct_ge_target     = round(100 * mean(composition$size >= target_all), 1),
    median_n_samples  = stats::median(composition$n_samples),
    pct_high_patient_purity = round(100 * mean(composition$patient_purity > high_purity), 1)
  )

  if ("cluster_purity" %in% names(composition)) {
    out$median_cluster_purity <- round(stats::median(composition$cluster_purity), 2)
  }

  for (contrast_name in names(contrast_n_samples)) {
    cells  <- composition[[paste0(contrast_name, "_cells")]]
    n_test <- composition[[paste0(contrast_name, "_n_test")]]
    n_ref  <- composition[[paste0(contrast_name, "_n_ref")]]
    target <- size_multiplier * contrast_n_samples[[contrast_name]]

    out[[paste0(contrast_name, "_mean_cells")]]    <- round(mean(cells), 1)
    out[[paste0(contrast_name, "_target")]]        <- target
    out[[paste0(contrast_name, "_pct_ge_target")]] <- round(100 * mean(cells >= target), 1)
    out[[paste0(contrast_name, "_pct_testable")]]  <-
      round(100 * mean(n_test >= min_samples_per_group & n_ref >= min_samples_per_group), 1)
  }

  out
}
#' Sweep k x d (x prop) for miloR differential-abundance neighborhoods
#'
#' Builds neighborhoods once per combination on the full object (the same
#' neighborhoods every contrast is later tested on) and summarizes size,
#' sample representation, patient dominance, cluster mixing, and, per
#' contrast, the cells and samples that contrast's test would actually use.
#'
#' @param milo Milo object (seurat_to_milo()) with sample/patient/group
#'   columns in colData.
#' @param k_grid,d_grid,prop_grid Candidate values.
#' @param refinement_scheme Passed to build_milo_nhoods().
#' @param sample_col,patient_col,group_col,cluster_col See nhood_composition().
#' @param contrasts Named list of c(test_level, reference_level).
#' @param size_multiplier,high_purity See summarize_nhood_composition().
#' @param seed Passed to build_milo_nhoods().
#'
#' @return A list: `summary` (one row per combination) and `sizes` (one row
#'   per neighborhood per combination, for histograms).
evaluate_milo_da_grid <- function(milo,
                                  k_grid,
                                  d_grid,
                                  prop_grid = 0.1,
                                  refinement_scheme = "graph",
                                  sample_col,
                                  patient_col,
                                  group_col,
                                  cluster_col = NULL,
                                  contrasts = NULL,
                                  size_multiplier = 5,
                                  high_purity = 0.8,
                                  seed = 123) {

  meta <- as.data.frame(SummarizedExperiment::colData(milo))
  n_samples_total <- dplyr::n_distinct(meta[[sample_col]])
  contrast_n_samples <- vapply(contrasts, function(lv) {
    dplyr::n_distinct(meta[[sample_col]][meta[[group_col]] %in% lv])
  }, numeric(1))

  summaries <- list()
  sizes <- list()

  for (prop in prop_grid) {
    for (k in k_grid) {
      for (d in d_grid) {
        combo <- paste0("k", k, "_d", d, if (length(prop_grid) > 1) paste0("_p", prop) else "")
        message("  ", combo)

        result <- tryCatch({
          milo_combo <- build_milo_nhoods(milo, k = k, d = d, prop = prop,
                                          refinement_scheme = refinement_scheme, seed = seed)
          membership <- miloR::nhoods(milo_combo)
          composition <- nhood_composition(membership, meta, sample_col, patient_col, group_col,
                                           cluster_col, contrasts)
          summary_row <- summarize_nhood_composition(
            composition, n_samples_total, contrast_n_samples,
            pct_cells_covered = 100 * mean(Matrix::rowSums(membership) > 0),
            size_multiplier = size_multiplier, high_purity = high_purity
          )
          list(summary = cbind(data.frame(k = k, d = d, prop = prop), summary_row),
               sizes = cbind(data.frame(combo = combo, k = k, d = d, prop = prop), composition))
        }, error = function(e) {
          message("    FAILED: ", conditionMessage(e))
          NULL
        })

        if (!is.null(result)) {
          summaries[[combo]] <- result$summary
          sizes[[combo]] <- result$sizes
        }
      }
    }
  }

  list(summary = dplyr::bind_rows(summaries), sizes = dplyr::bind_rows(sizes))
}
#' Drop redundant neighborhoods with a greedy set cover, in plain R
#'
#' Same algorithm as miloDE::filter_neighbourhoods(), which calls
#' RcppGreedySetCover::greedySetCover(): repeatedly keep the neighborhood
#' holding the most not-yet-covered cells until every covered cell is
#' covered. Written in R to avoid the compiled set-cover routine, which can
#' abort the R session on some installations. Ties go to the lower
#' neighborhood index, so the result is deterministic; it can differ from
#' miloDE's own choice only where two neighborhoods tie.
#'
#' @param milo Milo object with neighborhoods.
#'
#' @return The Milo object with the kept neighborhoods, their index cells,
#'   and a rebuilt neighborhood graph.
filter_nhoods_greedy <- function(milo) {

  membership <- miloR::nhoods(milo)
  uncovered <- Matrix::rowSums(membership) > 0
  keep <- integer(0)

  while (any(uncovered)) {
    gain <- Matrix::colSums(membership[uncovered, , drop = FALSE])
    gain[keep] <- -1
    best <- which.max(gain)
    if (gain[best] <= 0) break
    keep <- c(keep, best)
    uncovered[membership[, best] > 0] <- FALSE
  }

  keep <- sort(keep)
  index_cells <- miloR::nhoodIndex(milo)
  miloR::nhoods(milo) <- membership[, keep, drop = FALSE]
  miloR::nhoodIndex(milo) <- index_cells[keep]
  miloR::buildNhoodGraph(milo)
}
#' Assign neighborhoods for miloDE testing
#'
#' builder = "miloDE": miloDE::assign_neighbourhoods() (second-order
#' neighborhoods when order = 2) with filtering = FALSE, then, if filtering
#' is requested, filter_nhoods_greedy() instead of miloDE's compiled filter.
#'
#' builder = "miloR": first-order miloR neighborhoods (buildGraph() +
#' makeNhoods(), the same code path as DA), made larger for DE simply by a
#' larger k. `order` is ignored. Use this if miloDE::assign_neighbourhoods()
#' is unstable on your installation.
#'
#' @param sce SingleCellExperiment (or Milo) with a "PCA" reducedDim.
#' @param k,d,prop Neighborhood parameters.
#' @param order miloDE builder only: 1 or 2.
#' @param filtering Logical. Drop redundant neighborhoods (greedy set cover).
#'   Default TRUE.
#' @param builder "miloDE" (default) or "miloR".
#' @param seed Set immediately before index-cell sampling.
#'
#' @return A Milo object with neighborhoods and a neighborhood graph.
assign_milode_nhoods <- function(sce,
                                 k,
                                 d,
                                 order = 2,
                                 prop = 0.2,
                                 filtering = TRUE,
                                 builder = "miloDE",
                                 seed = 123) {

  if (builder == "miloR") {
    milo <- build_milo_nhoods(miloR::Milo(sce), k = k, d = d, prop = prop,
                              refinement_scheme = "graph", seed = seed)
    milo <- miloR::buildNhoodGraph(milo)
  } else if (builder == "miloDE") {
    set.seed(seed)
    milo <- suppressMessages(miloDE::assign_neighbourhoods(
      sce, reducedDim_name = "PCA", k = k, prop = prop, order = order,
      filtering = FALSE, d = d, verbose = FALSE
    ))
  } else {
    stop("assign_milode_nhoods(): builder must be \"miloDE\" or \"miloR\".", call. = FALSE)
  }

  if (isTRUE(filtering)) milo <- filter_nhoods_greedy(milo)
  milo
}
#' Sweep k for miloDE neighborhoods, separately for each contrast
#'
#' miloDE pseudobulks each neighborhood by sample, dropping samples with
#' fewer than min_n_cells_per_sample cells. Each contrast is evaluated on
#' an object holding only that contrast's cells (the same object the
#' pipeline tests on), and a neighborhood counts as testable when both
#' groups keep at least min_samples_per_group samples.
#'
#' @param sce SingleCellExperiment or Milo object with a "PCA" reducedDim
#'   and sample/patient/group columns in colData.
#' @param k_grid Candidate k values.
#' @param d Number of PCA dimensions.
#' @param order,prop,filtering Passed to assign_milode_nhoods().
#' @param sample_col,patient_col,group_col,cluster_col See nhood_composition().
#' @param contrasts Named list of c(test_level, reference_level).
#' @param min_n_cells_per_sample Same threshold the DE test will use.
#' @param min_samples_per_group Default 2.
#' @param builder Passed to assign_milode_nhoods().
#' @param seed Set before each neighborhood assignment.
#'
#' @return A list: `summary` (one row per contrast x k) and `sizes`.
evaluate_milode_grid <- function(sce,
                                 k_grid,
                                 d,
                                 order = 2,
                                 prop = 0.2,
                                 filtering = TRUE,
                                 sample_col,
                                 patient_col,
                                 group_col,
                                 cluster_col = NULL,
                                 contrasts,
                                 min_n_cells_per_sample = 3,
                                 min_samples_per_group = 2,
                                 builder = "miloDE",
                                 seed = 123) {

  # plain SCE: subsetting a Milo object does not subset its neighborhood slots
  sce <- methods::as(sce, "SingleCellExperiment")
  summaries <- list()
  sizes <- list()

  for (contrast_name in names(contrasts)) {
    levels_used <- contrasts[[contrast_name]]
    keep <- as.character(SummarizedExperiment::colData(sce)[[group_col]]) %in% levels_used
    sce_contrast <- sce[, keep]
    meta <- as.data.frame(SummarizedExperiment::colData(sce_contrast))
    n_samples <- dplyr::n_distinct(meta[[sample_col]])

    for (k in k_grid) {
      message("  ", contrast_name, ": k", k)

      result <- tryCatch({
        milo_de <- assign_milode_nhoods(sce_contrast, k = k, d = d, order = order,
                                        prop = prop, filtering = filtering,
                                        builder = builder, seed = seed)
        membership <- miloR::nhoods(milo_de)
        composition <- nhood_composition(membership, meta, sample_col, patient_col, group_col,
                                         cluster_col, contrasts[contrast_name],
                                         min_cells_per_sample = min_n_cells_per_sample)
        n_test <- composition[[paste0(contrast_name, "_n_test")]]
        n_ref  <- composition[[paste0(contrast_name, "_n_ref")]]

        summary_row <- data.frame(
          contrast          = contrast_name,
          k                 = k,
          n_cells           = ncol(sce_contrast),
          n_samples         = n_samples,
          n_nhoods          = nrow(composition),
          pct_cells_covered = round(100 * mean(Matrix::rowSums(membership) > 0), 1),
          median_size       = stats::median(composition$size),
          median_samples_passing = stats::median(n_test + n_ref),
          pct_testable      = round(100 * mean(n_test >= min_samples_per_group & n_ref >= min_samples_per_group), 1),
          pct_all_samples_pass = round(100 * mean((n_test + n_ref) == n_samples), 1)
        )
        if ("cluster_purity" %in% names(composition)) {
          summary_row$median_cluster_purity <- round(stats::median(composition$cluster_purity), 2)
        }
        list(summary = summary_row,
             sizes = cbind(data.frame(contrast = contrast_name, combo = paste0("k", k), k = k), composition))
      }, error = function(e) {
        message("    FAILED: ", conditionMessage(e))
        NULL
      })

      if (!is.null(result)) {
        summaries[[paste(contrast_name, k)]] <- result$summary
        sizes[[paste(contrast_name, k)]] <- result$sizes
      }
    }
  }

  list(summary = dplyr::bind_rows(summaries), sizes = dplyr::bind_rows(sizes))
}
#' Neighborhood-size histograms, one panel per parameter combination
#'
#' @param sizes_df `sizes` element from evaluate_milo_da_grid() or
#'   evaluate_milode_grid().
#' @param size_col Column to plot. Default "size" (all cells); use
#'   "<contrast>_cells" for the cells a contrast's test uses.
#' @param target Optional value drawn as a dashed line (e.g. 5 x samples).
#' @param title,subtitle Plot title/subtitle.
#'
#' @return A ggplot object.
plot_nhood_size_distributions <- function(sizes_df,
                                          size_col = "size",
                                          target = NULL,
                                          title = NULL,
                                          subtitle = NULL) {

  validate_required_cols(sizes_df, c("combo", size_col), "plot_nhood_size_distributions")
  sizes_df$combo <- factor(sizes_df$combo, levels = unique(sizes_df$combo))

  p <- ggplot2::ggplot(sizes_df, ggplot2::aes(x = .data[[size_col]])) +
    ggplot2::geom_histogram(bins = 25, fill = "#666666", color = "#FFFFFF") +
    ggplot2::facet_wrap(~ combo, scales = "free_y") +
    ggplot2::labs(title = title, subtitle = subtitle, x = size_col, y = "Neighborhoods") +
    ggplot2::theme_classic(base_size = 11)

  if (!is.null(target)) {
    p <- p + ggplot2::geom_vline(xintercept = target, linetype = "dashed", color = "#B22222")
  }
  p
}
