# ------------------------------------------------------------------------------
# config.R (workflows/01_qc)
#
# Default configuration and optional section toggles for the QC workflow.
# Dataset-specific values belong in config/params_override.R; the QC cutoffs
# chosen after looking at 01_qc_prefilter.R's output belong in
# config/decisions.R. setup.R stops on any setting name not defined here.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1. Step toggles
# ------------------------------------------------------------------------------
# Optional sections within 01_qc_prefilter.R and 02_qc_postfilter.R
steps <- list(
  mad_estimation       = TRUE,   # suggested cutoffs from median absolute deviation (prefilter)
  retained_vs_excluded = TRUE,   # DE and GO/KEGG enrichment, excluded vs retained cells (postfilter)
  spatial_qc           = FALSE   # QC by plate/well position (postfilter; plate-based protocols only)
)

# ------------------------------------------------------------------------------
# 2. Parameters
# ------------------------------------------------------------------------------
params <- list(

  # Cells dropped before QC metrics are computed: metadata column -> values,
  # e.g. list(Patient = c("P01"))
  exclusions = list(),

  # Percent-feature patterns (regex for Seurat::PercentageFeatureSet()), one
  # percent.<name> column each. Human gene symbols; override for other species.
  percent_features = list(
    mt   = "^MT-",
    ribo = "^RP[SL]"
  ),

  # MAD outlier direction per QC metric column ("both", "lower" or "higher",
  # passed to scuttle::isOutlier()); a metric not listed, or NULL, gets no
  # suggested cutoff
  mad_types = list(
    nFeature_RNA = "both",     # low = debris/empty wells; high = doublets
    nCount_RNA   = "higher",   # high library size = potential doublets
    percent.mt   = "higher"    # high mitochondrial fraction = dying/damaged cells
  ),

  # Number of MADs from the median for a suggested cutoff
  mad_nmads = 3,

  # Metadata columns to break QC metrics down by. The first two also set which
  # palette each column uses (group_var1, group_var2).
  group_by_vars = character(0),

  # Cell count breakdowns: each entry is 1-3 metadata columns
  # (1: table; 2: table + heatmap; 3: table + heatmap faceted by the third)
  cell_count_vars = list(),

  # Level order in plots and tables: column -> ordered levels
  display_order = list(),

  # Excluded vs retained cells (steps$retained_vs_excluded), passed to
  # run_de_enrichment(). GO/KEGG use human annotation (org.Hs.eg.db, KEGG
  # "hsa"): switch the step off for other species.
  retained_vs_excluded = list(
    min_pct               = 0.10,
    log2fc_cutoff         = 0.50,
    padj_cutoff           = 0.05,
    exclude_gene_patterns = c("^ENSG"),
    n_terms               = 10
  ),

  # Plate/well columns (steps$spatial_qc)
  spatial_qc = list(
    plate_col  = "Plate",
    row_col    = "Row",
    column_col = "Column"
  ),

  # Optional scatter of two median QC metrics per group combination, e.g.
  # list(group_vars = c("Patient", "Timepoint"), x_metric = "percent.mt",
  #      y_metric = "percent.ribo"). NULL = not drawn.
  diagnostic_scatter = NULL,

  # QC cutoffs: set in config/decisions.R
  thresholds = NULL,

  # Palette per plot context (see shared/functions/palettes.R); empty = defaults
  palettes = list(),

  # Data/results locations; NULL = folders inside this repo
  paths = list(
    data_raw_dir       = NULL,   # default: data/raw
    data_processed_dir = NULL,   # default: data/processed
    results_dir        = NULL,   # default: results
    raw_data_file      = NULL    # default: the single .rds in data_raw_dir
  )
)

# ------------------------------------------------------------------------------
# 3. Checkpoints written by this workflow (in data/processed/)
# ------------------------------------------------------------------------------
#   01_qc_prefilter.R  -> prefilter_exclusions
#   02_qc_postfilter.R -> qc_filtered  (input to 02_population_refinement)
checkpoints <- c(
  prefilter_exclusions = "01_qc_prefilter_exclusions.rds",
  qc_filtered          = "01_qc_filtered.rds"
)
