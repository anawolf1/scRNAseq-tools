# ------------------------------------------------------------------------------
# decisions.R (workflows/02_population_refinement)
#
# EXAMPLE - copy this file to config/decisions.R and edit it for your
# dataset. Sourced by 01_population_refinement.R and the tools.
#
# The history of population-refinement rounds for this dataset. A round is
# one cycle of:
#   load input -> steps (remove cells/genes) -> preprocess ->
#   (batch correction) -> cluster -> UMAP -> characterize -> save
#
# A round finds populations; the next round removes them: steps run on the
# round's input, before clustering.
#
# Append a new round for each cycle; do not edit a past round's values
# afterwards - if a round needs redoing, add a new one. This file plus git
# history is the record of your choices. Parameters are chosen with the
# scripts in tools/.
#
# Required fields per round: input, algorithm, dims, resolution. A field
# not listed below stops the pipeline; fields ending in "sweep" are
# documentation only.
# ------------------------------------------------------------------------------

params$population_refinement$rounds <- list(

  # --- round 1: screen the QC-filtered cells ------------------------------------
  round_1 = list(
    label      = "Lineage screening",   # title of every figure and table of the round;
                                        # also names the results folder
    rationale  = "Find populations that do not belong before the final clustering",   # optional note
    input      = "qc_filtered",         # checkpoint name (config.R `checkpoints`), an
                                        # earlier round's name, or a file path
    algorithm  = 1,                     # 1 = Louvain, 4 = Leiden (Seurat::FindClusters())
    dims       = 10,                    # PCs used for clustering and UMAP
    resolution = 0.6,

    # optional - default to params$clustering$k_param and $prune_snn
    # k_param   = 20,
    # prune_snn = 1 / 15,

    # optional - remove params$clustering$exclude_gene_patterns genes before
    # preprocessing (e.g. ribosomal, mitochondrial, TCR V(D)J genes)
    # exclude_genes = TRUE,

    # optional - batch correction before clustering; `method` selects the
    # implementation in apply_batch_correction(), the other fields are its arguments
    # batch_correction = list(method = "harmony", group.by = "Patient"),

    # documentation only: the screening behind dims/resolution above (printed
    # by tools/select_clustering_parameters.R)
    sweep = list(
      input = "qc_filtered", algorithm = 1,
      dims_tested = c(6, 8, 10, 12, 15, 18, 20),
      resolution_tested = seq(0.3, 1.0, by = 0.1),
      k_tested = 20, prune_tested = 1 / 15
    ),

    # no exclusion here: this round only clusters and characterizes, and the
    # next round decides what to remove (see round_2)
    steps = list(),

    # cluster number -> name, for THIS round's clustering only (cluster
    # numbers are not stable across rounds). Fill in after looking at the
    # round's results; used in plot labels, the cluster_label column and
    # exclusion by label.
    cluster_labels = c("0" = "Memory",
                       "1" = "Effector",
                       "2" = "Low quality",
                       "3" = "Naive")
  ),

  # --- round 2: remove the low-quality cluster, then cluster what is left ---------
  round_2 = list(
    label         = "Final clustering",
    input         = "round_1",
    algorithm     = 4,
    dims          = 12,
    resolution    = 0.5,
    exclude_genes = TRUE,
    sweep = list(
      input = "round_1", algorithm = 4, with_round_steps = TRUE,
      dims_tested = c(6, 8, 10, 12, 15, 18, 20),
      resolution_tested = seq(0.3, 1.0, by = 0.1)
    ),

    # cells or genes to remove from the input before clustering, applied in
    # order, selecting on the input's metadata (here round 1's cluster_label).
    # expect_n_cells (optional, recommended for cluster-based steps) stops
    # the round if the step would remove a different number of cells than
    # when you decided on it, e.g. because cluster numbers changed.
    steps = list(
      list(type = "exclude_metadata", column = "cluster_label", values = "Low quality",
           expect_n_cells = 58)

      # specific cells, inline or from a file saved with save_cell_list()
      # list(type = "exclude_cells", cells_file = "low_quality_half.txt"),
      # list(type = "exclude_cells", cell_ids = c("CELL_A1", "CELL_B7")),

      # genes; pattern = TRUE treats `genes` as one regular expression
      # list(type = "exclude_genes", genes = c("XIST", "MALAT1")),
      # list(type = "exclude_genes", genes = "^RPL|^RPS", pattern = TRUE)
    ),
    cluster_labels = c("1" = "Memory",
                       "2" = "Effector",
                       "3" = "Naive"),

    # optional - DE between every pair of clusters (volcanoes, tables and a
    # workbook in the round's pairwise_de/ folder). TRUE for all clusters, or
    # leave some out by label. Can be added after the round has run.
    pairwise_de = list(exclude = c("Effector"))
  )

  # --- comparing alternatives: parallel rounds from the same input ----------------
  # e.g. remove all / half / a third of the low-quality cluster, all clustered with
  # the same parameters so only the removed cells differ; compare with
  # tools/compare_rounds.R, keep one as final_round, and leave the others here
  # as the record of what you tried.
  # round_2_half = list(
  #   label = "Remove half of the low-quality cluster", input = "round_1", algorithm = 4,
  #   dims = 12, resolution = 0.5, exclude_genes = TRUE,
  #   steps = list(list(type = "exclude_cells", cells_file = "low_quality_half.txt")),
  #   cluster_labels = c(...)
  # )

)

# the round whose result downstream workflows use: saved as
# data/processed/02_population_refined.rds with a cluster_label column
params$population_refinement$final_round <- "round_2"
