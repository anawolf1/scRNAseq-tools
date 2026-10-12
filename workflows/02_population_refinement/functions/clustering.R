# ------------------------------------------------------------------------------
# clustering.R
#
# Seurat preprocessing, clustering and UMAP. Every setting is an argument;
# none of these functions reads the config or saves to disk. Evaluation of
# clusterings (sweeps, stability) is in clustering_evaluation.R.
#
# Functions
# preprocess_seurat()         NormalizeData -> FindVariableFeatures -> ScaleData -> RunPCA
# cluster_seurat()            FindNeighbors -> FindClusters
# run_umap()                  RunUMAP on a computed reduction
# validate_reduction_dims()   check a reduction and a number of its dimensions
# ------------------------------------------------------------------------------

#' Normalize, scale and run PCA
#'
#' LogNormalize, variable features (vst), ScaleData (optionally regressing
#' metadata columns) and PCA on the variable features. The number of PCs to
#' cluster on is chosen separately, e.g. from the elbow plot.
#'
#' @param obj Seurat object.
#' @param nfeatures Number of variable features. Default 2000.
#' @param vars_to_regress Metadata columns regressed out in ScaleData(), e.g.
#'   c("percent.mt"). Default NULL: none.
#' @param show_elbow Logical. TRUE (default): also return an elbow plot.
#' @param ndims PCs shown in the elbow plot. Default 50.
#' @param reset_analysis_state Logical. TRUE (default): remove existing
#'   reductions, graphs and neighbors first, so none is left from an earlier
#'   round (e.g. a harmony reduction from a round with batch correction).
#' @param seed Seed for RunPCA() (passed as `seed.use`). Default 123.
#'
#' @return A list: `obj` (normalized, scaled, with PCA) and `elbow_plot` (a
#'   ggplot, or NULL with show_elbow = FALSE).
#'
#' @examples
#' \dontrun{
#' result <- preprocess_seurat(obj, nfeatures = 2000, vars_to_regress = "percent.mt")
#' print(result$elbow_plot)
#' }
preprocess_seurat <- function(obj,
                              nfeatures = 2000,
                              vars_to_regress = NULL,
                              show_elbow = TRUE,
                              ndims = 50,
                              reset_analysis_state = TRUE,
                              seed = 123) {

  if (!inherits(obj, "Seurat")) {
    stop("preprocess_seurat(): `obj` must be a Seurat object, got: ", paste(class(obj), collapse = "/"),
         call. = FALSE)
  }
  validate_required_cols(obj[[]], vars_to_regress, "preprocess_seurat",
                         hint = "expects vars_to_regress to name metadata columns.")

  if (reset_analysis_state) {
    obj@reductions <- list()
    obj@graphs <- list()
    obj@neighbors <- list()
  }

  obj <- Seurat::NormalizeData(obj, normalization.method = "LogNormalize", scale.factor = 10000)
  obj <- Seurat::FindVariableFeatures(obj, selection.method = "vst", nfeatures = nfeatures)
  obj <- Seurat::ScaleData(obj, vars.to.regress = vars_to_regress)
  obj <- Seurat::RunPCA(obj, features = Seurat::VariableFeatures(obj), seed.use = seed)

  list(obj = obj, elbow_plot = if (show_elbow) Seurat::ElbowPlot(obj, ndims = ndims) else NULL)
}
#' Cluster on a computed reduction
#'
#' FindNeighbors (annoy nearest-neighbor search) and FindClusters with a
#' fixed seed, so the same settings give the same clusters.
#'
#' @param obj Seurat object with the reduction computed.
#' @param resolution FindClusters() resolution.
#' @param algorithm FindClusters() algorithm: 1 = Louvain, 2 = Louvain with
#'   multilevel refinement, 3 = SLM, 4 = Leiden.
#' @param reduction Reduction to cluster on. Default "pca".
#' @param dims Number of dimensions of `reduction` to use. Default 20.
#' @param k_param Nearest neighbors (FindNeighbors() `k.param`). Default 20.
#' @param prune_snn SNN pruning (FindNeighbors() `prune.SNN`). Default 1/15.
#' @param seed Seed for FindClusters(). Default 123.
#' @param verbose Logical. Print the cluster sizes. Default TRUE.
#'
#' @return `obj` with the neighbor graphs and `seurat_clusters`.
#'
#' @examples
#' \dontrun{
#' obj <- cluster_seurat(obj, resolution = 0.8, algorithm = 4, dims = 20)
#' }
cluster_seurat <- function(obj,
                           resolution,
                           algorithm,
                           reduction = "pca",
                           dims = 20,
                           k_param = 20,
                           prune_snn = 1 / 15,
                           seed = 123,
                           verbose = TRUE) {

  validate_reduction_dims(obj, reduction, dims, "cluster_seurat")

  obj <- Seurat::FindNeighbors(obj, dims = 1:dims, reduction = reduction, k.param = k_param,
                               prune.SNN = prune_snn, nn.method = "annoy", compute.SNN = TRUE,
                               verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = resolution, algorithm = algorithm,
                              random.seed = seed, verbose = FALSE)
  if (verbose) print(table(obj$seurat_clusters))
  obj
}
#' UMAP on a computed reduction
#'
#' For display only; clustering does not use it.
#'
#' @param obj Seurat object with the reduction computed.
#' @param dims Number of dimensions of `reduction` to use.
#' @param reduction Reduction to use. Default "pca".
#' @param n_neighbors,min_dist Passed to Seurat::RunUMAP(). Defaults 20 and
#'   0.01 give a more local layout than Seurat's 30 and 0.3, suited to a few
#'   thousand cells.
#' @param seed Seed. Default 123.
#'
#' @return `obj` with a "umap" reduction.
#'
#' @examples
#' \dontrun{
#' obj <- run_umap(obj, dims = 15)
#' }
run_umap <- function(obj,
                     dims,
                     reduction = "pca",
                     n_neighbors = 20,
                     min_dist = 0.01,
                     seed = 123) {

  validate_reduction_dims(obj, reduction, dims, "run_umap")

  Seurat::RunUMAP(obj, dims = 1:dims, reduction = reduction, n.neighbors = n_neighbors,
                  min.dist = min_dist, seed.use = seed, verbose = FALSE)
}
#' Check a reduction and a number of its dimensions
#'
#' @param obj Seurat object.
#' @param reduction Reduction name.
#' @param dims Number of dimensions to use.
#' @param fn_name Calling function's name (for the error message prefix).
#'
#' @return Invisibly TRUE; stops otherwise.
validate_reduction_dims <- function(obj,
                                    reduction,
                                    dims,
                                    fn_name) {
  if (!inherits(obj, "Seurat")) {
    stop(fn_name, "(): `obj` must be a Seurat object, got: ", paste(class(obj), collapse = "/"),
         call. = FALSE)
  }
  if (!is.numeric(dims) || length(dims) != 1 || !is.finite(dims) || dims < 1 || dims != as.integer(dims)) {
    stop(fn_name, "(): `dims` must be a single positive integer.", call. = FALSE)
  }
  if (!reduction %in% Seurat::Reductions(obj)) {
    stop(fn_name, "(): reduction '", reduction, "' not found. Available: ",
         paste(Seurat::Reductions(obj), collapse = ", "), call. = FALSE)
  }
  n_available <- ncol(Seurat::Embeddings(obj, reduction))
  if (dims > n_available) {
    stop(fn_name, "(): dims (", dims, ") exceeds the ", n_available, " dimensions of '",
         reduction, "'.", call. = FALSE)
  }
  invisible(TRUE)
}
