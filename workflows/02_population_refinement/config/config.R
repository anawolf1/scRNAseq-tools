# ------------------------------------------------------------------------------
# config.R (workflows/02_population_refinement)
#
# Default configuration and optional section toggles for the population
# refinement workflow. Dataset-specific values belong in
# config/params_override.R; the rounds you decide on (clustering parameters,
# exclusions, cluster labels) belong in config/decisions.R. setup.R stops on
# any setting name not defined here.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1. Step toggles
# ------------------------------------------------------------------------------
# Optional sections of each round's characterization (01_population_refinement.R)
steps <- list(
  composition_stability   = TRUE,   # cluster sizes and per-cluster stability (Jaccard)
  qc_metrics              = TRUE,   # QC medians, violin and feature plots by cluster
  group_umaps             = TRUE,   # UMAP colored by each params$group_by_vars column
  cell_counts             = TRUE,   # cell counts per params$cell_count_vars entry, and per cluster and group
  marker_characterization = TRUE,   # marker-column UMAPs and tables, signature-gene violins
  differential_expression = TRUE,   # one-vs-rest DE per cluster: volcanoes, tables, workbook
  deg_overlap             = TRUE,   # Venn diagrams of DEGs shared between clusters
                                    # (needs differential_expression)
  exclusion_check         = TRUE    # rounds that remove cells: excluded vs retained DE,
                                    # GO/KEGG and percent retained (human annotation)
)

# ------------------------------------------------------------------------------
# 2. Parameters
# ------------------------------------------------------------------------------
params <- list(

  # --- study description (same values as in the QC workflow) -------------------
  # Metadata columns for UMAPs and breakdowns. The first two also set which
  # palette each column uses (group_var1, group_var2).
  group_by_vars = character(0),

  # Cell count breakdowns: each entry is 1-3 metadata columns
  cell_count_vars = list(),

  # Level order in plots and tables: column -> ordered levels
  display_order = list(),

  # Percent-feature patterns: only the names are used here, to find the
  # percent.<name> columns the QC workflow added
  percent_features = list(
    mt   = "^MT-",
    ribo = "^RP[SL]"
  ),

  # --- preprocessing and clustering ---------------------------------------------
  clustering = list(
    vars_to_regress = character(0),   # metadata columns regressed out in ScaleData()
    nfeatures       = 2000,           # variable features used for PCA
    random_seed     = 123,            # seed for PCA, clustering, subsampling and UMAP
    k_param         = 20,             # nearest neighbors, for rounds that do not set k_param
    prune_snn       = 1 / 15,         # SNN pruning, for rounds that do not set prune_snn
    # Genes removed before preprocessing when a round sets exclude_genes = TRUE.
    # Ribosomal and mitochondrial genes mostly reflect technical or stress
    # variation; cell-type-specific patterns (e.g. TCR V(D)J genes) belong in
    # params_override.R.
    exclude_gene_patterns = c("^RPL", "^RPS", "^MT-")
  ),

  # --- differential expression ----------------------------------------------------
  # One-vs-rest and pairwise DE of a round's clusters. exclude_gene_patterns
  # removes genes from every DE result, table and plot (e.g. "^ENSG"); it does
  # not change what clustering uses (clustering$exclude_gene_patterns does).
  de = list(
    exclude_gene_patterns = character(0),
    min_pct               = 0.1,
    log2fc_cutoff         = 0.5,
    padj_cutoff           = 0.05,
    n_labels              = 50     # top up and top down genes labeled in each volcano
  ),

  # --- characterization of each round -----------------------------------------------
  characterization = list(
    violin_genes     = character(0),   # signature genes for the per-cluster violin plot
    highlight_genes  = character(0),   # genes whose DE results are reported for every cluster
    stability_n_reps = 10,             # subsampling repeats for per-cluster stability
    stability_frac   = 0.8             # fraction of cells kept in each subsample
  ),

  # Fixed color per cluster label (label -> hex), so a population keeps its
  # color across rounds; NULL = colors by position
  cluster_color_map = NULL,

  # --- DEG overlap between clusters (steps$deg_overlap) --------------------------------
  # One Venn for genes up in each cluster vs the rest, one for genes down.
  # clusters: 2-4 cluster labels or IDs; NULL = every cluster (skipped with a
  # message if there are more than 4).
  deg_overlap = list(
    clusters      = NULL,
    log2fc_cutoff = 0.5,
    padj_cutoff   = 0.05
  ),

  # --- excluded vs retained cells (steps$exclusion_check) ------------------------------
  # Passed to run_de_enrichment(); genes are excluded with de$exclude_gene_patterns
  exclusion_check = list(
    min_pct       = 0.1,
    log2fc_cutoff = 0.5,
    padj_cutoff   = 0.05,
    n_terms       = 10
  ),

  # --- clustering parameter screening (tools/) -----------------------------------------
  # Marker thresholds are looser than `de`: they are applied to many
  # candidate clusterings, to flag clusters without any marker support.
  parameter_selection = list(
    min_ari_stability       = NULL,   # candidate filter: minimum mean ARI; NULL = none
    max_invalid_frac        = 0,      # candidate filter: maximum fraction of invalid clusters
    min_markers_per_cluster = 1,      # a cluster with fewer markers is invalid
    min_cluster_size        = 10,     # a cluster with fewer cells is invalid
    marker_count_threshold  = 10,     # candidates table: clusters with fewer markers are counted
    marker_min_pct          = 0.01,
    marker_log2fc_cutoff    = 0.25,
    marker_padj_cutoff      = 0.05,
    stability_n_reps        = 10,     # subsampling repeats for the ARI stability
    stability_frac          = 0.8,    # fraction of cells kept in each subsample
    marker_check_cols       = NULL    # metadata column -> level expected in one cluster,
                                      # e.g. list(Naive = "Naive", CD127 = "CD127-")
  ),

  # --- round engine ------------------------------------------------------------------------
  # rounds and final_round are set in config/decisions.R
  population_refinement = list(
    rounds      = NULL,
    final_round = NULL,
    force_rerun = FALSE   # TRUE = recompute every round, even if its checkpoint is current
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
#   qc_filtered        - read: written by the QC workflow (01_qc)
#   population_refined - written: the final round's result, with a
#                        cluster_label column; input to downstream workflows
# Each round also writes 02_population_refinement_<round>.rds (see
# functions/population_refinement.R).
checkpoints <- c(
  qc_filtered        = "01_qc_filtered.rds",
  population_refined = "02_population_refined.rds"
)
