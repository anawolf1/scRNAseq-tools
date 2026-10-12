# ------------------------------------------------------------------------------
# config.R (workflows/03_functional_scores)
#
# Default configuration and optional section toggles for the functional
# scores workflow: module scores of curated gene programs, summarized by
# cluster and timepoint, with timepoint effects tested on patient means.
# Dataset-specific values belong in config/params_override.R; which cells
# each comparison uses belongs in config/decisions.R. setup.R stops on any
# setting name not defined here.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 1. Step toggles
# ------------------------------------------------------------------------------
# Comparisons produced; the cells each one uses are set in config/decisions.R
# (params$analysis_cells)
steps <- list(
  cluster_summaries        = TRUE,    # score table and violins by cluster (cells, descriptive)
  program_umaps            = TRUE,    # UMAP of each program's z-scored score
  timepoint_summaries      = TRUE,    # score table and violins by timepoint (cells, descriptive)
  trajectory_plot          = TRUE,    # per-patient scores across the repeated timepoints
  longitudinal_analysis    = FALSE,   # LMM: score ~ timepoint + (1 | patient), repeated timepoints
  cross_sectional_analysis = FALSE    # Wilcoxon: repeated timepoints vs cross-sectional group(s)
)

# ------------------------------------------------------------------------------
# 2. Parameters
# ------------------------------------------------------------------------------
params <- list(

  # Checkpoint to score: a name in `checkpoints` below or a file path
  input = "population_refined",

  # --- study description (same values as in the other workflows) --------------
  # The first two group_by_vars columns set the palettes (group_var1, group_var2)
  group_by_vars = character(0),

  # Level order in plots and tables: column -> ordered levels
  display_order = list(),

  # Fixed color per cluster label (label -> hex), the same map as population
  # refinement; NULL = colors by position
  cluster_color_map = NULL,

  # Patient and timepoint columns, and which timepoints are repeated
  # measures (the same patients) or a separately sampled cross-sectional group
  study_design = list(
    patient_col                  = "Patient",
    timepoint_col                = "Timepoint",
    repeated_measures_timepoints = character(0),   # e.g. c("BL", "3M", "FU")
    cross_sectional_timepoints   = character(0)    # e.g. c("SPR")
  ),

  # --- gene sets --------------------------------------------------------------------
  gene_sets = list(
    dir          = NULL,          # NULL = data/gene_sets
    files        = NULL,          # NULL = every .gmt file in `dir`
    # species of each gene set ("human" or "mouse"), required for every set:
    # human sets get current symbols, mouse sets are mapped to human orthologs
    species      = character(0),
    # replacements of .gmt symbols before resolution (matched ignoring case),
    # e.g. c("4-SEP" = "SEPT4") for symbols damaged by spreadsheets
    symbol_fixes = character(0)
  ),

  # --- scoring ------------------------------------------------------------------------
  scoring = list(
    # "all_cells":  score every cell once; a cell has the same score in every
    #               cell set, so differences between sets come only from
    #               which cells are compared
    # "per_subset": score each cell set separately (control genes binned on
    #               that set); scores are not comparable between sets
    mode       = "all_cells",
    seed       = 123,
    net_scores = TRUE   # <base>_NET = <base>_UP - <base>_DN for every paired program
  ),

  # --- statistics (on patient x timepoint means, never on cells) -------------------------
  statistics = list(
    min_cells_per_sample        = 0,                 # samples with fewer cells are left out
    min_patients_lmm            = 3,                 # patients needed to fit a program's LMM
    min_patients_wilcox         = 2,                 # patients needed in each group for a Wilcoxon test
    # df method of the LMM contrasts: "kenward-roger" (pbkrtest),
    # "satterthwaite" (lmerTest) or "asymptotic" (anticonservative with few patients)
    lmer_df                     = "kenward-roger",
    # list of c(group1, group2) timepoint pairs for the Wilcoxon tests; NULL =
    # every repeated timepoint vs every cross-sectional group
    cross_sectional_comparisons = NULL
  ),

  # --- cells used by each comparison (set in config/decisions.R) ----------------------
  # default: every cell, for every comparison
  cell_sets = list(
    all_cells = list(label = "All cells")
  ),
  analysis_cells = list(
    cluster_summaries        = "all_cells",
    program_umaps            = "all_cells",
    timepoint_summaries      = "all_cells",
    trajectory_plot          = "all_cells",
    longitudinal_analysis    = "all_cells",
    cross_sectional_analysis = "all_cells"
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
#   population_refined - read: written by population refinement (02)
#   functional_scores  - written: every cell with its scores (scoring mode
#                        "all_cells"); input to later workflows
#   term2gene_used     - written: the resolved gene sets, genes in the object
checkpoints <- c(
  population_refined = "02_population_refined.rds",
  functional_scores  = "03_functional_scores.rds",
  term2gene_used     = "03_term2gene_used.rds"
)
