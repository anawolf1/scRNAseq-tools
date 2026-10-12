# ------------------------------------------------------------------------------
# config.R (workflows/04_neighborhood_analysis)
#
# Default configuration and optional section toggles for the neighborhood
# analysis workflow: differential abundance (miloR) and neighborhood-level
# differential expression (miloDE) between timepoints. Dataset-specific
# values belong in config/params_override.R; the analyses and the
# neighborhood parameters you choose belong in config/decisions.R. setup.R
# stops on any setting name not defined here.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1. Step toggles
# ------------------------------------------------------------------------------
steps <- list(
  milo_da              = TRUE,   # miloR differential abundance per analysis and contrast
  nhood_program_scores = TRUE,   # functional-score programs per neighborhood (needs 03_functional_scores)
  da_nhood_markers     = TRUE,   # meta-neighborhoods + chosen neighborhoods: marker DE, within-sample check
  nhood_groups         = TRUE,   # groups of ALL neighborhoods (overlap), heatmaps, maps, group/cluster tests
  milo_de              = TRUE    # miloDE neighborhood differential expression, primary contrasts
)

# ------------------------------------------------------------------------------
# 2. Parameters
# ------------------------------------------------------------------------------
params <- list(

  # Checkpoint to analyze: a name in `checkpoints` below or a file path
  input = "population_refined",

  # --- study description (same values as in the other workflows) --------------
  group_by_vars     = character(0),
  display_order     = list(),
  cluster_color_map = NULL,          # fixed color per cluster label (the same map as 02 and 03)
  study_design = list(
    patient_col                  = "Patient",
    timepoint_col                = "Timepoint",
    repeated_measures_timepoints = character(0),
    cross_sectional_timepoints   = character(0)
  ),

  # --- what is compared (set in config/decisions.R) ---------------------------------
  # analyses: which cells build each neighborhood graph, and which contrasts
  # (c(test, reference) timepoints) are tested on it - see resolve_milo_analyses()
  analyses = list(),

  # --- miloR differential abundance ----------------------------------------------------
  da = list(
    k                 = NULL,      # from tools/select_milo_parameters.R (Part A)
    d                 = NULL,      # PCs for the kNN graph; from the same tool
    prop              = 0.1,       # share of cells sampled as neighborhood index cells
    refinement_scheme = "graph",   # "graph" or "reduced_dim"; sets the matching SpatialFDR weighting
    spatial_fdr       = 0.1,       # SpatialFDR threshold for a significant neighborhood
    mixed_threshold   = 0.5,       # majority cluster below this share -> labeled "Mixed" (an annotation rule)
    size_multiplier   = 5,         # target neighborhood size: this x the samples compared
    norm_method       = "TMM",     # testNhoods() normalization: "TMM" (miloR default), "logMS"
                                   # (total cells per sample) or "RLE"
    paired_design     = TRUE,      # contrasts between repeated-measures timepoints: ~ patient + condition
    n_nhoods_to_show  = 12,        # neighborhoods shown in per-sample and membership plots
    nhood_qc_cols     = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.ribo"),
    # groups of all neighborhoods (steps$nhood_groups): Leiden on the
    # neighborhood-overlap graph (Jaccard of shared cells), no DA input.
    # resolution 1 = standard modularity; nearby resolutions are reported to
    # show how much the grouping depends on it
    nhood_groups = list(
      resolution          = 1,        # choose with tools/select_milo_parameters.R (Part C)
      resolutions_checked = c(0.5, 1, 1.5, 2),
      pairwise_de         = TRUE,     # DE between every pair of groups (volcanoes + workbook)
      seed                = 123,
      sweep               = NULL      # record of the resolution screen (documentation only)
    ),
    label_nhoods      = TRUE,      # number the significant neighborhoods on the DA graph
    seed              = 123,
    sweep             = NULL,      # record of the parameter screen (documentation only)
    # marker DE of meta-neighborhoods and chosen neighborhoods (steps$da_nhood_markers)
    markers = list(
      min_overlap           = 0.3,    # significant neighborhoods sharing >= this share of cells
                                      # (same direction) are merged into one meta-neighborhood
      min_cells             = 20,     # groups / comparison sets smaller than this are skipped
      log2fc_cutoff         = 0.25,
      min_pct               = 0.1,
      padj_cutoff           = 0.05,
      top_n                 = 6,      # top markers checked within each sample
      within_sample_min_cells = 3,    # cells needed inside and outside the set in a sample
      contrast_cells_only   = FALSE,  # TRUE: meta-neighborhoods use only the contrast's two timepoints
      exclude_gene_patterns = character(0),   # genes hidden from marker tables (e.g. "^ENSG")
      nhood_sets            = NULL    # hand-picked neighborhoods (set in decisions.R)
    )
  ),

  # --- miloDE neighborhood differential expression -----------------------------------
  de = list(
    # "milode": separate, larger neighborhoods per contrast, built on that
    #           contrast's cells (more power; linked to DA by shared cells)
    # "da":     DE inside the DA neighborhoods (same IDs as DA; smaller)
    neighborhoods         = "milode",
    builder                = "miloDE",   # "milode" only: "miloDE" (second order) or "miloR" (first order)
    k                      = NULL,       # from tools/select_milo_parameters.R (Part B)
    d                      = NULL,       # NULL = da$d
    order                  = 2,          # "milode" + "miloDE" only
    prop                   = 0.2,
    filtering              = TRUE,       # drop redundant neighborhoods
    paired_design          = TRUE,       # contrasts between repeated-measures timepoints: ~ patient + condition
    min_n_cells_per_sample = 3,          # samples with fewer cells in a neighborhood leave its test
    min_count              = 3,          # edgeR filterByExpr() min.count
    alpha                  = 0.1,        # within-neighborhood FDR for calling a gene DE
    seed                   = 123,
    n_top_genes            = 50,         # genes in the top-genes table
    extreme_logfc          = 5,          # |logFC| above this is flagged as on/off (pct_extreme_logFC)
    genes_to_plot          = character(0),
    sweep                  = NULL
  ),

  # Palette per plot context (see shared/functions/palettes.R); empty = defaults
  palettes = list(),

  # Data/results locations; NULL = folders inside this repo
  paths = list(
    data_processed_dir = NULL,   # default: data/processed
    results_dir        = NULL    # default: results
  )
)

# ------------------------------------------------------------------------------
# 3. Checkpoints (in data/processed/)
# ------------------------------------------------------------------------------
#   population_refined - read: final clustering (02_population_refinement)
#   functional_scores  - read, optional: module scores (03_functional_scores)
#   milo_da            - written: primary analysis' Milo object + DA results
# Other analyses write 04_milo_da_<analysis>.rds, miloDE 04_milo_de_<contrast>.rds
# (each with a manifest).
checkpoints <- c(
  population_refined = "02_population_refined.rds",
  functional_scores  = "03_functional_scores.rds",
  milo_da            = "04_milo_da.rds"
)
