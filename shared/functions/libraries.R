# ------------------------------------------------------------------------------
# libraries.R
#
# Single registry of every package this repo's code uses directly, with what
# each is for and which workflow(s) use it. tools/write_description.R copies
# this list into DESCRIPTION so that renv records it (renv.lock also pins
# their dependencies). `core = TRUE` packages are loaded by default by
# load_libraries(); each workflow's setup.R loads the rest it needs by name.
#
# Usage:
#   (sourced automatically by shared/functions/load_all.R)
#   load_libraries()                              # core only
#   load_libraries("aricode", "fpc")              # core + these
#   load_libraries("harmony", load_core = FALSE)  # only this, skip core
#   library_help()                                # print the whole registry
#   library_help("volcano")                       # search by name, description or workflow
#   library_help("01_qc")                         # packages a workflow uses
# ------------------------------------------------------------------------------

library_registry <- tibble::tribble(
  ~package,               ~core, ~used_by,                                                                   ~description,
  "here",                 TRUE,  "all",                                                                      "Project-root path resolution via `here::here()` for reproducible file paths",
  "rlang",                TRUE,  "all",                                                                      "Tidy-evaluation helpers and the `%||%` operator",
  "ggplot2",              TRUE,  "all",                                                                      "Plotting",
  "Seurat",               TRUE,  "all",                                                                      "Core scRNAseq object and workflow",
  "dplyr",                TRUE,  "all",                                                                      "Data wrangling",
  "tibble",               TRUE,  "all",                                                                      "Modern data frames; tribble() builds the library and palette registries",
  "tidyr",                TRUE,  "all",                                                                      "Create tidy data",
  "gt",                   TRUE,  "all",                                                                      "Formatted tables",
  "webshot2",             TRUE,  "all",                                                                      "Renders gt tables to PNG (required by save_gt())",
  "patchwork",            TRUE,  "all",                                                                      "Combining ggplot panels into one figure",
  "scales",               TRUE,  "all",                                                                      "Color scales for gt_de(), viridis palettes",
  "Matrix",               FALSE, "all",                                                                      "Sparse matrix support (Seurat dependency)",
  "writexl",              FALSE, "all",                                                                      "Writes result tables to Excel",
  "scuttle",              FALSE, "01_qc",                                                                    "MAD-based QC cutoff suggestions (isOutlier)",
  "ggnewscale",           FALSE, "01_qc",                                                                    "Multiple fill scales in one plot (plot_plate_grid, spatial QC)",
  "clusterProfiler",      FALSE, "01_qc, 02_population_refinement",                                          "GO / KEGG enrichment analysis",
  "enrichplot",           FALSE, "01_qc, 02_population_refinement",                                          "Plots for clusterProfiler enrichment results",
  "org.Hs.eg.db",         FALSE, "01_qc, 02_population_refinement, 03_functional_scores",                    "Human gene ID/annotation database",
  "ggrepel",              FALSE, "01_qc, 04_neighborhood_analysis",                                          "Non-overlapping text labels on plots",
  "cluster",              FALSE, "02_population_refinement",                                                 "Silhouette width for clustering evaluation",
  "aricode",              FALSE, "02_population_refinement",                                                 "Adjusted Rand Index / clustering agreement metrics",
  "fpc",                  FALSE, "02_population_refinement",                                                 "Calinski-Harabasz index",
  "EnhancedVolcano",      FALSE, "02_population_refinement",                                                 "Volcano plots for DE results",
  "harmony",              FALSE, "02_population_refinement",                                                 "Batch/sample integration (only if a round requests it)",
  "ggvenn",               FALSE, "02_population_refinement",                                                 "Venn diagrams of DEG overlap across conditions (compare_deg_venn)",
  "limma",                FALSE, "03_functional_scores",                                                     "alias2SymbolTable() - resolves outdated human gene symbols",
  "babelgene",            FALSE, "03_functional_scores",                                                     "orthologs() - curated mouse-to-human ortholog mapping",
  "lme4",                 FALSE, "03_functional_scores",                                                     "Linear mixed-effects models",
  "emmeans",              FALSE, "03_functional_scores",                                                     "Estimated marginal means and post-hoc contrasts",
  "pbkrtest",             FALSE, "03_functional_scores",                                                     "Kenward-Roger degrees of freedom for LMM contrasts (emmeans)",
  "lmerTest",             FALSE, "03_functional_scores",                                                     "Optional: Satterthwaite df for LMM contrasts (params$statistics$lmer_df)",
  "ggpubr",               FALSE, "03_functional_scores",                                                     "stat_pvalue_manual() - draws significance brackets from a precomputed p-value table",
  "miloR",                FALSE, "04_neighborhood_analysis",                                                 "Differential abundance analysis with Milo",
  "miloDE",               FALSE, "04_neighborhood_analysis",                                                 "Differential expression analysis with Milo",
  "SingleCellExperiment", FALSE, "04_neighborhood_analysis",                                                 "Single-cell data container for Bioconductor workflows",
  "SummarizedExperiment", FALSE, "04_neighborhood_analysis",                                                 "colData() access on Milo/SingleCellExperiment objects",
  "edgeR",                FALSE, "04_neighborhood_analysis",                                                 "Abundance tests of neighborhood groups and clusters (also used by miloR)",
  "pheatmap",             FALSE, "04_neighborhood_analysis",                                                 "Annotated, hclust-ordered program heatmaps of neighborhoods and groups",
  "igraph",               FALSE, "04_neighborhood_analysis",                                                 "Neighborhood graphs: grouping significant neighborhoods, Leiden groups of all neighborhoods"
)

#' Load core libraries plus any additional named packages
#'
#' Loads all packages flagged `core = TRUE` in `library_registry` by
#' default, plus any extra packages named in `...`. Set `load_core = FALSE`
#' to load only the packages named in `...`.
#'
#' @param ... Package names (as strings) to load in addition to (or instead
#'   of) the core set. Must match a `package` entry in `library_registry`.
#' @param load_core Logical. If TRUE (default), also loads every package
#'   flagged `core = TRUE`.
#' @param strict Logical. If TRUE (default), stops with an error listing
#'   any requested package that isn't installed, rather than silently
#'   skipping it. Set FALSE only for ad hoc interactive use where continuing
#'   without some packages is fine.
#'
#' @return Invisibly, a character vector of the packages that were loaded.
#'
#' @examples
#' \dontrun{
#' load_libraries()
#' load_libraries("aricode", "fpc")
#' load_libraries("harmony", load_core = FALSE)
#' }
load_libraries <- function(...,
                           load_core = TRUE,
                           strict = TRUE) {

  extra <- c(...)

  unknown <- setdiff(extra, library_registry$package)
  if (length(unknown) > 0) {
    stop("load_libraries(): Unknown package(s) not in library_registry: ", paste(unknown, collapse = ", "),
         "\nAdd them to library_registry in libraries.R first.",
         call. = FALSE)
  }

  core_pkgs <- if (load_core) library_registry$package[library_registry$core] else character(0)
  libs <- unique(c(core_pkgs, extra))

  installed <- vapply(libs, requireNamespace, logical(1), quietly = TRUE)
  missing <- libs[!installed]

  if (length(missing) > 0) {
    if (strict) {
      stop("load_libraries(): Missing package(s): ", paste(missing, collapse = ", "),
           ". Run renv::restore() to install everything pinned in renv.lock, ",
           "or call load_libraries(..., strict = FALSE) to continue without them.",
           call. = FALSE)
    }
    message("Missing packages (not loaded): ", paste(missing, collapse = ", "))
  }

  invisible(lapply(libs[installed], function(pkg) {
    suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  }))

  invisible(libs[installed])
}
#' Look up what a package in the registry is used for
#'
#' Prints the library registry, optionally filtered to entries whose
#' package name, description or used_by (workflow) matches `pattern`.
#'
#' @param pattern Optional regex (case-insensitive) matched against the
#'   package name, description and used_by columns. Can be a single string or
#'   a character vector of several patterns (rows matching ANY pattern are
#'   returned). NULL (default) prints the full registry.
#'
#' @return Invisibly, the (possibly filtered) registry data.frame.
#'
#' @examples
#' \dontrun{
#' library_help()
#' library_help("volcano")
#' library_help(c("volcano", "enrich"))
#' }
library_help <- function(pattern = NULL) {
  df <- library_registry
  if (!is.null(pattern)) {
    pattern <- paste(pattern, collapse = "|")
    match_pkg  <- grepl(pattern, df$package, ignore.case = TRUE)
    match_desc <- grepl(pattern, df$description, ignore.case = TRUE)
    match_used <- grepl(pattern, df$used_by, ignore.case = TRUE)
    df <- df[match_pkg | match_desc | match_used, ]
  }
  df <- df[order(-df$core, df$package), ]
  print(df, n = Inf)
  invisible(df)
}
