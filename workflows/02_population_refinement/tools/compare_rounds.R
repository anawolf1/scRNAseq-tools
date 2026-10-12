# ------------------------------------------------------------------------------
# compare_rounds.R (workflows/02_population_refinement/tools)
#
# Compares rounds run as alternatives, typically alternative exclusions
# (e.g. removing all, half or a third of a low-quality cluster), each a
# round from the same input. See README.md, "Comparing alternative
# exclusions".
#
# Reads what the pipeline saved for each round (no reclustering):
#   - a summary table: cells kept, clusters, clustering parameters
#   - cluster correspondence: the labels of the cells two rounds share,
#     each round against the first one listed
#   - per population (cluster label) in 2-4 rounds: Venns of its DEGs
#     across rounds (genes up vs the rest, then down) and DEG counts, i.e.
#     whether the population's signature depends on the exclusion
#   - optional GO/KEGG per round and population
#
# Each round's exclusion check (excluded vs retained, percent retained) is
# in its own results folder, under exclusion_check/.
#
# Output: results/02_population_refinement/exploration/round_comparison/<label>/
#
# Sections:
#   0. Setup
#   1. Rounds to compare
#   2. Load what each round saved
#   3. Summary and cluster correspondence
#   4. DEGs per population
#
# Usage:
#   source(here::here("workflows", "02_population_refinement", "tools", "compare_rounds.R"))
#   Prerequisite: the rounds have been run by 01_population_refinement.R.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------
source(here::here("workflows", "02_population_refinement", "setup.R"))
decision_file <- file.path(workflow_dir, "config", "decisions.R")
if (!file.exists(decision_file)) {
  stop("No workflows/02_population_refinement/config/decisions.R found; the rounds to compare are defined there.")
}
source(decision_file)

# ------------------------------------------------------------------------------
# 1. Rounds to compare (edit per use)
# ------------------------------------------------------------------------------
rounds_to_compare <- c("round_2_full", "round_2_half", "round_2_third")
comparison_label  <- "low_quality_exclusion"   # names the output folder
run_enrichment    <- FALSE                     # GO/KEGG per round and population (human annotation)

rounds <- params$population_refinement$rounds
unknown_rounds <- setdiff(rounds_to_compare, names(rounds))
if (length(unknown_rounds) > 0) {
  stop("Not rounds in config/decisions.R: ", paste(unknown_rounds, collapse = ", "))
}
if (length(rounds_to_compare) < 2) stop("Choose at least 2 rounds to compare.")
if (run_enrichment) load_libraries("clusterProfiler", "org.Hs.eg.db", "enrichplot", load_core = FALSE)

out_dir <- file.path(fig_exploration_dir, "round_comparison", comparison_label)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# short display names: the round's label if set
display_names <- setNames(
  make.unique(vapply(rounds_to_compare, function(r) rounds[[r]]$label %||% r, character(1), USE.NAMES = FALSE)),
  rounds_to_compare
)

# ------------------------------------------------------------------------------
# 2. Load what each round saved
# ------------------------------------------------------------------------------
round_data <- lapply(rounds_to_compare, function(r) {
  checkpoint <- refinement_checkpoint_path(r, data_processed_dir)
  if (!file.exists(checkpoint)) stop("Round '", r, "' has no checkpoint; run the pipeline first.")
  obj <- readRDS(checkpoint)
  # the current labels from config/decisions.R (they may postdate the run)
  obj <- add_cluster_label_column(obj, rounds[[r]]$cluster_labels)

  # enrichment background: the genes detected in at least params$de$min_pct
  # of the round's cells, as in the DE tests
  universe <- NULL
  if (run_enrichment) {
    detected <- Matrix::rowMeans(Seurat::GetAssayData(obj, layer = "counts") > 0) >= params$de$min_pct
    universe <- rownames(obj)[detected]
    if (length(params$de$exclude_gene_patterns) > 0) {
      universe <- universe[!grepl(paste(params$de$exclude_gene_patterns, collapse = "|"), universe)]
    }
  }

  r_dir <- round_results_dir(r, rounds[[r]], fig_refinement_dir)
  summary_file <- file.path(r_dir, "round_summary.csv")
  de_file <- list.files(r_dir, pattern = "_de_results\\.csv$", full.names = TRUE)

  list(
    labels   = setNames(as.character(obj$cluster_label), colnames(obj)),
    universe = universe,
    summary  = if (file.exists(summary_file)) utils::read.csv(summary_file) else NULL,
    de       = if (length(de_file) == 1) {
      de <- utils::read.csv(de_file)
      # the current labels, as for the checkpoint above
      ids <- as.character(de$cluster)
      de$cluster_label <- ifelse(ids %in% names(rounds[[r]]$cluster_labels),
                                 unname(rounds[[r]]$cluster_labels[ids]), ids)
      de
    } else NULL
  )
})
names(round_data) <- rounds_to_compare

missing_de <- rounds_to_compare[vapply(round_data, function(x) is.null(x$de), logical(1))]
if (length(missing_de) > 0) {
  message("[info] no saved DE table for: ", paste(missing_de, collapse = ", "),
          " (steps$differential_expression was off); left out of the DEG comparison.")
}

# ------------------------------------------------------------------------------
# 3. Summary and cluster correspondence
# ------------------------------------------------------------------------------
summaries <- do.call(rbind, lapply(round_data, `[[`, "summary"))
if (!is.null(summaries)) {
  summary_gt <- gt_generic(summaries, title = paste("Round comparison:", comparison_label))
  print(summary_gt)
  save_gt(summary_gt, "summary.png", out_dir)
}

# the label of each shared cell in both rounds: whether populations stay
# the same when different cells are removed
reference <- rounds_to_compare[1]
correspondence_sheets <- list()
for (r in rounds_to_compare[-1]) {
  shared_cells <- intersect(names(round_data[[reference]]$labels), names(round_data[[r]]$labels))
  crosstab <- table(round_data[[reference]]$labels[shared_cells], round_data[[r]]$labels[shared_cells],
                    dnn = c(display_names[[reference]], display_names[[r]]))
  crosstab_gt <- gt_contingency(
    crosstab, title = paste("Cluster correspondence:", display_names[[reference]], "vs", display_names[[r]]),
    subtitle = paste(length(shared_cells), "cells present in both rounds")
  )
  print(crosstab_gt)
  save_gt(crosstab_gt, paste0("correspondence_", reference, "_vs_", r, ".png"), out_dir)
  correspondence_sheets[[paste0("vs_", r)]] <- as.data.frame.matrix(crosstab)
}

# ------------------------------------------------------------------------------
# 4. DEGs per population
# ------------------------------------------------------------------------------
# a population is a cluster label, compared across the rounds that have it
de_rounds <- setdiff(rounds_to_compare, missing_de)
populations <- sort(unique(unlist(lapply(round_data[de_rounds], function(x) x$de$cluster_label))))

for (population in populations) {
  de_list <- lapply(de_rounds, function(r) {
    de <- round_data[[r]]$de
    de[!is.na(de$cluster_label) & de$cluster_label == population, , drop = FALSE]
  })
  names(de_list) <- display_names[de_rounds]
  de_list <- Filter(function(df) nrow(df) > 0, de_list)

  if (length(de_list) < 2 || length(de_list) > 4) {
    message("[skip] ", population, ": in ", length(de_list), " round(s) with DE results; a Venn needs 2-4.")
    next
  }

  pop_prefix <- paste0(gsub("[^A-Za-z0-9]+", "", population), "_by_round")
  for (direction in c("up", "down")) {
    venn <- compare_deg_venn(
      de_list,
      log2fc_cutoff         = params$deg_overlap$log2fc_cutoff,
      padj_cutoff           = params$deg_overlap$padj_cutoff,
      exclude_gene_patterns = params$de$exclude_gene_patterns,
      direction             = direction,
      title_label           = paste(population, "by round"),
      output_dir            = out_dir,
      file_prefix           = pop_prefix
    )
    plot_deg_count_comparison(venn$gene_sets, output_dir = out_dir,
                              file_prefix = paste0(pop_prefix, "_", direction))

    if (run_enrichment && direction == "up") {
      for (round_display in names(venn$gene_sets)) {
        r <- names(display_names)[display_names == round_display]
        enrich_gene_list(
          venn$gene_sets[[round_display]], label = paste(population, "-", round_display),
          output_dir = out_dir, file_prefix = paste0(pop_prefix, "_", gsub("[^A-Za-z0-9]+", "", round_display)),
          universe = round_data[[r]]$universe
        )
      }
    }
  }
}

if (length(correspondence_sheets) > 0 || !is.null(summaries)) {
  save_xlsx(c(if (!is.null(summaries)) list(summary = summaries), correspondence_sheets),
            "round_comparison.xlsx", out_dir)
}
message("Round comparison saved in ", out_dir)

# ------------------------------------------------------------------------------
# Next step (by hand):
#   Decide which alternative to keep, note why in its `rationale`, and set
#   params$population_refinement$final_round to it (or build the next round
#   on it). The other rounds stay in config/decisions.R as the record of
#   what was compared.
# ------------------------------------------------------------------------------
