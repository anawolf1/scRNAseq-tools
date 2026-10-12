# ------------------------------------------------------------------------------
# params_override.R (workflows/02_population_refinement)
#
# EXAMPLE - copy this file to config/params_override.R (that file name is
# gitignored, so dataset values are never committed) and edit it for your
# dataset.
#
# Defaults and the meaning of every setting are in config/config.R; only set
# what differs. Set single entries (params$section$entry <- value) rather
# than replacing a whole list. The rounds you decide on go in
# config/decisions.R.
# ------------------------------------------------------------------------------

# --- optional sections ------------------------------------------------------------
# switch a section of each round's characterization on or off (defaults in config.R)

# steps$marker_characterization <- FALSE
# steps$deg_overlap             <- FALSE
# steps$exclusion_check         <- FALSE   # also off for non-human data (GO/KEGG)


# --- study description: the same values as in the QC workflow's override ------------

params$group_by_vars <- c("Patient", "Timepoint")

params$cell_count_vars <- list(
  c("Patient", "Timepoint")
)

params$display_order <- list(
  Timepoint = c("BL", "3M", "FU", "SPR")
)

# only if the QC workflow used other percent_features names
# params$percent_features$hb <- "^HB[^(P)]"


# --- preprocessing and clustering -----------------------------------------------------
# metadata columns regressed out in ScaleData()

# params$clustering$vars_to_regress <- c("nFeature_RNA", "percent.mt")

# genes removed before preprocessing in rounds with exclude_genes = TRUE, on top
# of the ribosomal/mitochondrial defaults - e.g. TCR V(D)J genes, which would
# otherwise separate T cells by clonotype rather than state

# params$clustering$exclude_gene_patterns <- c(
#   params$clustering$exclude_gene_patterns, "^TR[ABDG][VDJ]"
# )


# --- differential expression -------------------------------------------------------------
# genes removed from every DE result, table and plot (not from clustering)

# params$de$exclude_gene_patterns <- c("^ENSG", "^RN7SL")
# params$de$log2fc_cutoff <- 1


# --- characterization ---------------------------------------------------------------------
# signature genes for the per-cluster violin plot (e.g. lineage markers, to
# spot contaminating populations)

# params$characterization$violin_genes <- c(
#   "CD3E", "CD8A",    # T cells
#   "CD14", "LYZ",     # monocytes
#   "MS4A1", "CD79A"   # B cells
# )

# genes whose DE results are reported for every cluster
# params$characterization$highlight_genes <- c("TRDC", "TRGC2")

# fixed color per cluster label, so a population keeps its color across
# rounds even when its cluster number changes
# params$cluster_color_map <- c(
#   "Naive"  = "#E69F00",
#   "Memory" = "#56B4E9"
# )


# --- DEG overlap (steps$deg_overlap) --------------------------------------------------------
# 2-4 cluster labels or IDs; NULL = every cluster

# params$deg_overlap$clusters <- c("Naive", "Memory", "Effector")


# --- parameter screening (tools/) -------------------------------------------------------------
# metadata column -> level expected to concentrate in one cluster, e.g. a
# protein-level gate. Use it to confirm known populations are resolved, and
# report that you checked (see evaluate_all_candidates_markers()).

# params$parameter_selection$marker_check_cols <- list(
#   Naive = "Naive",
#   CD127 = "CD127-"
# )
# params$parameter_selection$min_ari_stability <- 0.8


# --- palettes ---------------------------------------------------------------------------------
# palette per plot context; palette_help() lists contexts and defaults

# params$palettes$cluster    <- "dark2"
# params$palettes$group_var1 <- "colorbrewer1"


# --- round engine -------------------------------------------------------------------------------
# TRUE recomputes every round, even if its checkpoint is current

# params$population_refinement$force_rerun <- TRUE


# --- paths: only for data/results outside this repo ----------------------------------------------

# params$paths$data_processed_dir <- "/path/to/processed"
# params$paths$results_dir        <- "/path/to/results"
