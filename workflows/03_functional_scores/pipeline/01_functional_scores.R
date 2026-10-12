# ------------------------------------------------------------------------------
# 01_functional_scores.R (workflows/03_functional_scores)
#
# Scores every cell for curated functional programs (module scores from
# data/gene_sets/*.gmt), then runs each enabled comparison on the cell
# set(s) chosen for it in config/decisions.R (e.g. clusters compared on all
# cells, timepoints with naive cells excluded):
#   - cluster and timepoint summaries, program UMAPs (cells, descriptive)
#   - per-patient trajectories across the repeated timepoints
#   - longitudinal analysis: LMM on the repeated timepoints
#   - cross-sectional analysis: Wilcoxon, repeated timepoints vs a
#     separately sampled group
#
# Statistical unit: every test uses patient means (one value per patient x
# timepoint), never cells. Cell-level violins are descriptive.
#
# Input:
#   data/processed/02_population_refined.rds (from 02_population_refinement)
#   data/gene_sets/*.gmt
#   config/decisions.R (cell sets per comparison; see decisions_example.R)
#
# Output:
#   results/03_functional_scores/
#     - functional_scores_gene_coverage.png
#     - <cell set>/: the comparisons that use this cell set, cells per
#       sample, and functional_scores_summary_<cell set>.xlsx with every table
#   data/processed/
#     - 03_term2gene_used.rds (+ _manifest.txt): the gene sets scored
#     - 03_functional_scores.rds (+ _manifest.txt): every cell with its
#       scores (scoring mode "all_cells"); input to later workflows. In
#       "per_subset" mode, 03_functional_scores_<cell set>.rds per cell set
#
# Sections:
#   0. Setup
#   1. Load the object
#   2. Gene sets
#   3. Module scores
#   4. Comparisons per cell set
#   5. Checkpoint
#
# Usage:
#   source(here::here("workflows", "03_functional_scores", "pipeline", "01_functional_scores.R"))
#
#   Prerequisites: 02_population_refinement run; .gmt files in
#   data/gene_sets/ with their species in config/params_override.R;
#   config/decisions.R copied from decisions_example.R.
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# 0. Setup
# ------------------------------------------------------------------------------
source(here::here("workflows", "03_functional_scores", "setup.R"))

decision_file <- file.path(workflow_dir, "config", "decisions.R")
if (file.exists(decision_file)) {
  source(decision_file)
} else {
  message("No workflows/03_functional_scores/config/decisions.R found - every comparison uses all cells. ",
          "Copy config/decisions_example.R to config/decisions.R to choose cell sets.")
}

# --- cell sets and the comparisons that use them ------------------------------------
cell_sets <- params$cell_sets
if (length(cell_sets) == 0 || is.null(names(cell_sets)) || any(!nzchar(names(cell_sets))) ||
    anyDuplicated(names(cell_sets))) {
  stop("params$cell_sets needs at least one uniquely named cell set. Check config/decisions.R.")
}
unknown_steps <- setdiff(names(params$analysis_cells), names(steps))
if (length(unknown_steps) > 0) {
  stop("params$analysis_cells names unknown comparison(s): ", paste(unknown_steps, collapse = ", "),
       ". Valid: ", paste(names(steps), collapse = ", "))
}
unknown_sets <- setdiff(unlist(params$analysis_cells), names(cell_sets))
if (length(unknown_sets) > 0) {
  stop("params$analysis_cells uses cell set(s) not in params$cell_sets: ", paste(unknown_sets, collapse = ", "))
}
enabled_steps <- names(steps)[vapply(steps, isTRUE, logical(1))]
missing_assignment <- setdiff(enabled_steps, names(params$analysis_cells))
if (length(missing_assignment) > 0) {
  stop("Enabled comparison(s) with no cell set in params$analysis_cells: ", paste(missing_assignment, collapse = ", "))
}
# the cell sets needed, in the order they are defined
sets_to_run <- intersect(names(cell_sets), unique(unlist(params$analysis_cells[enabled_steps])))

# --- scoring and statistics settings ---------------------------------------------------
scoring_mode <- params$scoring$mode
if (!scoring_mode %in% c("all_cells", "per_subset")) {
  stop("params$scoring$mode must be \"all_cells\" or \"per_subset\", got: ", scoring_mode)
}
stats_settings <- params$statistics
lmer_df <- stats_settings$lmer_df
lmer_df_pkg <- c("kenward-roger" = "pbkrtest", "satterthwaite" = "lmerTest", "asymptotic" = NA)
if (!lmer_df %in% names(lmer_df_pkg)) {
  stop("params$statistics$lmer_df must be one of: ", paste(names(lmer_df_pkg), collapse = ", "))
}

# --- study design ------------------------------------------------------------------------
patient_col   <- params$study_design$patient_col
timepoint_col <- params$study_design$timepoint_col
rm_timepoints <- params$study_design$repeated_measures_timepoints
cs_timepoints <- params$study_design$cross_sectional_timepoints
design_timepoints <- c(rm_timepoints, cs_timepoints)

if (anyDuplicated(design_timepoints)) {
  stop("A timepoint is in both repeated_measures_timepoints and cross_sectional_timepoints.")
}
if (isTRUE(steps$longitudinal_analysis) && length(rm_timepoints) < 2) {
  stop("steps$longitudinal_analysis is on but params$study_design$repeated_measures_timepoints has fewer ",
       "than 2 timepoints. Set it in params_override.R or switch the step off.")
}
if (isTRUE(steps$longitudinal_analysis) && !is.na(lmer_df_pkg[[lmer_df]])) {
  load_libraries(lmer_df_pkg[[lmer_df]], load_core = FALSE)
}
if (isTRUE(steps$cross_sectional_analysis) && length(cs_timepoints) == 0) {
  stop("steps$cross_sectional_analysis is on but params$study_design$cross_sectional_timepoints is empty. ",
       "Set it in params_override.R or switch the step off.")
}

# Wilcoxon pairs: unpaired tests, so two repeated timepoints (the same
# patients) cannot be a pair
cross_sectional_comparisons <- stats_settings$cross_sectional_comparisons %||%
  unlist(lapply(rm_timepoints, function(rm) lapply(cs_timepoints, function(cs) c(rm, cs))), recursive = FALSE)
if (isTRUE(steps$cross_sectional_analysis)) {
  for (pair in cross_sectional_comparisons) {
    if (length(pair) != 2 || !all(pair %in% design_timepoints)) {
      stop("Each params$statistics$cross_sectional_comparisons entry must be two study-design timepoints, got: ",
           paste(pair, collapse = ", "))
    }
    if (all(pair %in% rm_timepoints)) {
      stop("cross_sectional_comparisons: ", paste(pair, collapse = " vs "), " are both repeated timepoints ",
           "(the same patients); the unpaired Wilcoxon test does not apply. The longitudinal analysis compares them.")
    }
  }
}

# ------------------------------------------------------------------------------
# 1. Load the object
# ------------------------------------------------------------------------------
input_path <- if (params$input %in% names(checkpoints)) checkpoint_path(params$input, data_processed_dir) else params$input
if (!file.exists(input_path)) {
  stop("Input not found: ", input_path, ". Run 02_population_refinement first (it writes ",
       checkpoints[["population_refined"]], "), or set params$input.")
}
data_obj <- readRDS(input_path)
input_manifest <- read_manifest(input_path)
message("Input: ", basename(input_path), " (", ncol(data_obj), " cells, ", nrow(data_obj), " genes",
        if (!is.null(input_manifest[["setting.final_round"]])) paste0(", round ", input_manifest[["setting.final_round"]]),
        ")")

validate_required_cols(data_obj[[]], c(patient_col, timepoint_col, "cluster_label"), "01_functional_scores",
                       hint = "expects the study-design columns (params$study_design) and population refinement's cluster_label.")
data_obj$cluster_label <- droplevels(as.factor(data_obj$cluster_label))

# timepoints in study-design order (or display order); values outside the
# study design are left out of the timepoint-based outputs
observed_timepoints <- unique(as.character(data_obj[[timepoint_col, drop = TRUE]]))
timepoint_levels <- if (length(design_timepoints) > 0) {
  design_timepoints
} else {
  ordered_levels(observed_timepoints, params$display_order[[timepoint_col]])
}
not_in_design <- setdiff(observed_timepoints, c(timepoint_levels, NA))
if (length(not_in_design) > 0) {
  message("[info] timepoint value(s) not in the study design, left out of timepoint-based outputs: ",
          paste(not_in_design, collapse = ", "))
}
data_obj[[timepoint_col]] <- factor(as.character(data_obj[[timepoint_col, drop = TRUE]]), levels = timepoint_levels)

# ------------------------------------------------------------------------------
# 2. Gene sets
# ------------------------------------------------------------------------------
if (length(gene_sets_files) == 0) {
  stop("No .gmt file found in ", gene_sets_dir, "; add the curated gene sets there first.")
}
gene_sets <- prepare_gene_sets(read_gmt(gene_sets_files), species = params$gene_sets$species,
                               symbol_fixes = params$gene_sets$symbol_fixes, object_genes = rownames(data_obj))
coverage <- gene_sets$coverage
term2gene_used <- gene_sets$term2gene[gene_sets$term2gene$gene %in% rownames(data_obj), ]

coverage_gt <- gt_generic(coverage, title = "Functional module scores — gene set coverage",
                          subtitle = "Genes per program: in the .gmt, after alias/ortholog resolution, found in the data")
print(coverage_gt)
save_gt(coverage_gt, "functional_scores_gene_coverage.png", fig_functional_dir, vwidth = 1600)

write_checkpoint(term2gene_used, "term2gene_used", data_processed_dir,
                 workflow = "03_functional_scores", script = "01_functional_scores.R",
                 settings = list(gene_set_files = basename(gene_sets_files),
                                 species = params$gene_sets$species,
                                 symbol_fixes = params$gene_sets$symbol_fixes),
                 input_path = input_path)

# ------------------------------------------------------------------------------
# 3. Module scores
# ------------------------------------------------------------------------------
# programs with at least one gene in the data, UP before DN; then net scores
program_names <- order_programs_up_dn(unique(term2gene_used$term))
net_bases <- if (isTRUE(params$scoring$net_scores)) net_program_bases(program_names) else character(0)

score_cell_set <- function(obj) {
  obj <- add_module_scores(obj, term2gene_used, seed = params$scoring$seed)
  if (length(net_bases) > 0) obj <- add_net_scores(obj, net_bases)
  obj
}
scored_obj <- if (scoring_mode == "all_cells") score_cell_set(data_obj) else data_obj

program_names <- c(program_names, paste0(net_bases, "_NET"))
score_cols    <- paste0(program_names, "Score")
program_labels <- setNames(vapply(program_names, function(x) wrap_label(gsub("_", " ", x), width = 20), character(1)),
                           program_names)
score_labels <- setNames(unname(program_labels), score_cols)

# the UP/DN programs and the net scores are drawn as separate figures
net_cols <- paste0(net_bases, "_NETScore")
violin_sets <- Filter(function(v) length(v$cols) > 0, list(
  list(cols = setdiff(score_cols, net_cols), file = "functional_scores", what = "Functional module scores",
       note = NULL),
  list(cols = net_cols, file = "functional_scores_net", what = "Net scores (UP − DN)",
       note = "net = UP score − DN score per cell; higher = closer to the UP signature")
))
message("Scoring mode: ", scoring_mode, "; programs: ", paste(program_names, collapse = ", "))

# ------------------------------------------------------------------------------
# 4. Comparisons per cell set
# ------------------------------------------------------------------------------
for (set_name in sets_to_run) {

  spec  <- cell_sets[[set_name]]
  label <- spec$label %||% set_name
  tag   <- paste0("_", set_name)
  out_dir <- file.path(fig_functional_dir, set_name)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  set_steps <- enabled_steps[vapply(enabled_steps, function(st) set_name %in% params$analysis_cells[[st]], logical(1))]
  runs_here <- function(step) step %in% set_steps
  message("\n[cell set] ", set_name, ": ", label, " - ", paste(set_steps, collapse = ", "))

  # --- select cells ------------------------------------------------------------------
  selection <- select_cell_set(scored_obj[[]], spec, set_name)
  obj <- subset(scored_obj, cells = colnames(scored_obj)[selection$keep])
  obj$cluster_label <- droplevels(obj$cluster_label)
  if (scoring_mode == "per_subset") obj <- score_cell_set(obj)
  message("  ", ncol(obj), " of ", ncol(scored_obj), " cells",
          if (nzchar(selection$rule_text)) paste0(" (", selection$rule_text, ")"))

  settings <- data.frame(
    setting = c("cell set", "label", "input", "cells kept", "cells in input", "rule", "comparisons",
                "scoring", "min cells per sample", "LMM df method"),
    value   = c(set_name, label, basename(input_path), ncol(obj), ncol(scored_obj),
                if (nzchar(selection$rule_text)) selection$rule_text else "none",
                paste(set_steps, collapse = ", "), scoring_mode, stats_settings$min_cells_per_sample, lmer_df)
  )
  subtitle_base <- paste0(label, " (n = ", ncol(obj), " cells)")

  # --- patient x timepoint means: the unit of every test ----------------------------------
  patient_scores   <- NULL
  cells_per_sample <- NULL
  if (any(c("trajectory_plot", "longitudinal_analysis", "cross_sectional_analysis") %in% set_steps)) {
    cells_per_sample <- as.data.frame(table(obj[[patient_col, drop = TRUE]], obj[[timepoint_col, drop = TRUE]]),
                                      stringsAsFactors = FALSE)
    names(cells_per_sample) <- c(patient_col, timepoint_col, "n_cells")
    cells_per_sample <- cells_per_sample[cells_per_sample$n_cells > 0, ]
    cells_per_sample$included <- cells_per_sample$n_cells >= stats_settings$min_cells_per_sample
    cells_per_sample <- cells_per_sample[order(match(cells_per_sample[[timepoint_col]], timepoint_levels),
                                               cells_per_sample[[patient_col]]), ]

    cells_per_sample_gt <- gt_generic(
      cells_per_sample, title = paste("Cells per sample —", label),
      subtitle = paste0("Samples with fewer than ", stats_settings$min_cells_per_sample,
                        " cells are left out of patient-level summaries and tests")
    )
    print(cells_per_sample_gt)
    save_gt(cells_per_sample_gt, paste0("cells_per_sample", tag, ".png"), out_dir)
    if (any(!cells_per_sample$included)) {
      message("  ", sum(!cells_per_sample$included), " sample(s) below min_cells_per_sample left out")
    }

    # the patient and timepoint columns are named Patient and Timepoint from here on
    patient_scores <- summarize_scores_by(obj, score_cols, c(patient_col, timepoint_col)) %>%
      dplyr::rename(Patient = dplyr::all_of(patient_col), Timepoint = dplyr::all_of(timepoint_col)) %>%
      dplyr::filter(!is.na(Timepoint), n >= stats_settings$min_cells_per_sample) %>%
      dplyr::mutate(program = factor(program, levels = program_names)) %>%
      dplyr::arrange(program, Timepoint, Patient) %>%
      dplyr::mutate(program = as.character(program))

    patient_scores_gt <- gt_module_scores(
      patient_scores, title = paste("Functional module scores by patient x timepoint —", label),
      subtitle = "Per-sample mean ± SD across cells; these means are the unit of every test"
    )
    print(patient_scores_gt)
    save_gt(patient_scores_gt, paste0("functional_scores_by_patient_timepoint", tag, ".png"), out_dir)
  }

  # --- trajectories ---------------------------------------------------------------------------
  if (runs_here("trajectory_plot")) {
    if (length(rm_timepoints) == 0) {
      message("  [skip] trajectory_plot (no repeated_measures_timepoints in the study design)")
    } else {
      if (length(cs_timepoints) > 1) {
        message("  trajectory plots: more than one cross-sectional group; none drawn as a reference")
      }
      for (prog in program_names) {
        p_trajectory <- plot_score_trajectories(
          patient_scores[patient_scores$program == prog, ],
          timepoint_levels      = rm_timepoints,
          cross_sectional_group = if (length(cs_timepoints) == 1) cs_timepoints,
          program_labels        = program_labels[prog],
          patient_col           = patient_col,
          title    = paste("Module score trajectory —", gsub("_", " ", prog)),
          subtitle = paste0(subtitle_base, "; patient means, ", paste(rm_timepoints, collapse = " → "))
        )
        print(p_trajectory)
        save_plot(p_trajectory, paste0("functional_score_trajectory_", gsub("[^A-Za-z0-9]+", "_", prog), tag, ".png"),
                  out_dir, width = 8, height = 6)
      }
    }
  }

  # --- cluster summaries ---------------------------------------------------------------------
  cluster_scores <- NULL
  if (runs_here("cluster_summaries")) {
    cluster_scores <- summarize_scores_by(obj, score_cols, "cluster_label") %>%
      dplyr::mutate(program = factor(program, levels = program_names)) %>%
      dplyr::arrange(program, cluster_label) %>%
      dplyr::mutate(program = as.character(program))

    cluster_scores_gt <- gt_module_scores_wide(cluster_scores, group_var = "cluster_label",
                                               title = paste("Functional module scores by cluster —", label),
                                               subtitle = "Mean and SD across cells")
    print(cluster_scores_gt)
    save_gt(cluster_scores_gt, paste0("functional_scores_by_cluster_table", tag, ".png"), out_dir)

    cluster_levels <- levels(obj$cluster_label)
    for (vs in violin_sets) {
      p_by_cluster <- plot_violin(
        obj, features = vs$cols, group_col = "cluster_label",
        group_colors = resolve_cluster_colors(cluster_levels, params$cluster_color_map), group_order = cluster_levels,
        feature_labels = score_labels, show_points = FALSE,
        title    = paste(vs$what, "by cluster —", label),
        subtitle = paste(c("Cell-level distributions (descriptive)", vs$note), collapse = "; ")
      )
      print(p_by_cluster)
      save_plot(p_by_cluster, paste0(vs$file, "_by_cluster", tag, ".png"), out_dir,
                width = max(8, 2.5 * length(vs$cols)), height = 5)
    }
  }

  # --- timepoint summaries ---------------------------------------------------------------------
  timepoint_scores <- NULL
  if (runs_here("timepoint_summaries")) {
    timepoint_obj <- subset(obj, cells = colnames(obj)[!is.na(obj[[timepoint_col, drop = TRUE]])])
    timepoint_scores <- summarize_scores_by(timepoint_obj, score_cols, timepoint_col) %>%
      dplyr::rename(Timepoint = dplyr::all_of(timepoint_col)) %>%
      dplyr::mutate(program = factor(program, levels = program_names)) %>%
      dplyr::arrange(program, Timepoint) %>%
      dplyr::mutate(program = as.character(program))

    timepoint_scores_gt <- gt_module_scores_wide(timepoint_scores, group_var = "Timepoint",
                                                 title = paste("Functional module scores by timepoint —", label),
                                                 subtitle = "Mean and SD across cells")
    print(timepoint_scores_gt)
    save_gt(timepoint_scores_gt, paste0("functional_scores_by_timepoint_table", tag, ".png"), out_dir)

    present_timepoints <- intersect(timepoint_levels, unique(as.character(timepoint_obj[[timepoint_col, drop = TRUE]])))
    for (vs in violin_sets) {
      p_by_timepoint <- plot_violin(
        timepoint_obj, features = vs$cols, group_col = timepoint_col, group_order = present_timepoints,
        feature_labels = score_labels, show_points = FALSE,
        title    = paste(vs$what, "by timepoint —", label),
        subtitle = paste(c("Cell-level distributions (descriptive; tests use patient means)", vs$note), collapse = "; ")
      )
      print(p_by_timepoint)
      save_plot(p_by_timepoint, paste0(vs$file, "_by_timepoint", tag, ".png"), out_dir,
                width = max(8, 2.5 * length(vs$cols)), height = 5)
    }
  }

  # --- program UMAPs ------------------------------------------------------------------------------
  if (runs_here("program_umaps")) {
    for (vs in violin_sets) {
      n_umap_cols <- min(4, length(vs$cols))
      p_umap <- plot_program_umap_grid(
        obj, vs$cols, program_labels = program_labels, ncol = n_umap_cols,
        title = paste(vs$what, "on UMAP —", label),
        subtitle = paste(c("Scores z-scored within each program across the cells of this set", vs$note),
                         collapse = "; ")
      )
      print(p_umap)
      save_plot(p_umap, paste0(vs$file, "_umap", tag, ".png"), out_dir,
                width = 3.5 * n_umap_cols + 1.5, height = 3.5 * ceiling(length(vs$cols) / n_umap_cols) + 1)
    }
  }

  # panel layout of the significance plots
  n_facet_cols <- min(5, length(program_names))
  facet_height <- 4 * ceiling(length(program_names) / n_facet_cols) + 1

  # --- longitudinal analysis (LMM) ---------------------------------------------------------
  lmm_coefficients <- NULL
  lrt_table        <- NULL
  contrast_table   <- NULL
  if (runs_here("longitudinal_analysis")) {
    rm_scores <- patient_scores[patient_scores$Timepoint %in% rm_timepoints, ]
    fits <- list()
    for (prog in program_names) {
      prog_data <- rm_scores[rm_scores$program == prog, ]
      if (dplyr::n_distinct(prog_data$Patient) < stats_settings$min_patients_lmm ||
          dplyr::n_distinct(prog_data$Timepoint) < 2) {
        message("  [skip] LMM for ", prog, ": fewer than ", stats_settings$min_patients_lmm,
                " patients or 2 timepoints")
        next
      }
      fits[[prog]] <- tryCatch(
        fit_timepoint_lmm(prog_data, rm_timepoints, lmer_df, program = prog),
        error = function(e) {
          message("  LMM failed for ", prog, ": ", conditionMessage(e))
          NULL
        }
      )
    }
    fits <- Filter(Negate(is.null), fits)

    if (length(fits) == 0) {
      message("  no LMM could be fitted; no longitudinal outputs")
    } else {
      lmm_coefficients <- dplyr::bind_rows(lapply(fits, `[[`, "coef"))
      lmm_gt <- gt_lmm(lmm_coefficients, title = paste("LMM fixed effects (Timepoint) —", label),
                       subtitle = paste0("score ~ Timepoint + (1 | Patient), REML; reference: ", rm_timepoints[1]))
      print(lmm_gt)
      save_gt(lmm_gt, paste0("functional_scores_lmm_coefficients", tag, ".png"), out_dir)

      lrt_table <- dplyr::bind_rows(lapply(fits, `[[`, "lrt")) %>%
        dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH"), significance = p_to_stars(p_adj))
      lrt_gt <- gt_generic(lrt_table, title = paste("Overall Timepoint effect (likelihood-ratio test) —", label),
                           subtitle = "score ~ Timepoint + (1 | Patient) vs. score ~ 1 + (1 | Patient), ML fits; BH across programs")
      print(lrt_gt)
      save_gt(lrt_gt, paste0("functional_scores_lmm_timepoint_lrt", tag, ".png"), out_dir)
      if (any(lrt_table$singular_fit)) {
        message("  singular LMM fit (patient variance estimated at ~0) for: ",
                paste(lrt_table$program[lrt_table$singular_fit], collapse = ", "))
      }

      contrast_table <- dplyr::bind_rows(lapply(fits, `[[`, "contrasts")) %>%
        dplyr::mutate(significance = p_to_stars(p_adj))
      contrast_gt <- gt_generic(contrast_table, title = paste("Pairwise Timepoint contrasts (LMM) —", label),
                                subtitle = paste0("emmeans, ", lmer_df, " df; BH within each program's ",
                                                  choose(length(rm_timepoints), 2), " comparisons"))
      print(contrast_gt)
      save_gt(contrast_gt, paste0("functional_scores_lmm_contrasts", tag, ".png"), out_dir)

      p_longitudinal <- plot_group_significance(
        rm_scores, contrast_table, group_var = "Timepoint", group_levels = rm_timepoints,
        group_colors = resolve_group_colors(timepoint_col, rm_timepoints),
        program_order = program_names, program_labels = program_labels,
        title = paste("Longitudinal module scores —", label),
        subtitle = paste0("Patient means; LMM contrasts (", lmer_df, " df), BH within program"),
        ncol = n_facet_cols, seed = params$scoring$seed
      )
      print(p_longitudinal)
      save_plot(p_longitudinal, paste0("functional_scores_longitudinal_significance", tag, ".png"), out_dir,
                width = 4 * n_facet_cols, height = facet_height)
    }
  }

  # --- cross-sectional analysis (Wilcoxon) ----------------------------------------------------
  cross_sectional_table <- NULL
  if (runs_here("cross_sectional_analysis")) {
    cross_sectional_table <- compare_groups_wilcox(patient_scores, cross_sectional_comparisons, program_names,
                                                   min_patients = stats_settings$min_patients_wilcox)
    if (anyNA(cross_sectional_table$p_value)) {
      message("  cross-sectional: some tests not run (fewer than ", stats_settings$min_patients_wilcox,
              " patients in a group)")
    }
    cross_gt <- gt_generic(cross_sectional_table,
                           title = paste("Cross-sectional comparisons (Wilcoxon rank-sum) —", label),
                           subtitle = "Patient means, unpaired; BH across programs within each comparison; normal approximation with ties")
    print(cross_gt)
    save_gt(cross_gt, paste0("functional_scores_cross_sectional", tag, ".png"), out_dir)

    tested <- cross_sectional_table[!is.na(cross_sectional_table$p_adj), ]
    if (nrow(tested) > 0) {
      plot_levels <- intersect(timepoint_levels, unique(unlist(cross_sectional_comparisons)))
      p_cross <- plot_group_significance(
        patient_scores[patient_scores$Timepoint %in% plot_levels, ], tested,
        group_var = "Timepoint", group_levels = plot_levels,
        group_colors = resolve_group_colors(timepoint_col, plot_levels),
        program_order = program_names, program_labels = program_labels,
        title = paste("Cross-sectional module scores —", label),
        subtitle = "Patient means; Wilcoxon rank-sum, BH within comparison",
        ncol = n_facet_cols, seed = params$scoring$seed
      )
      print(p_cross)
      save_plot(p_cross, paste0("functional_scores_cross_sectional_plot", tag, ".png"), out_dir,
                width = 4 * n_facet_cols, height = facet_height)
    }
  }

  # --- workbook ---------------------------------------------------------------------------------
  sheets <- Filter(Negate(is.null), list(
    settings             = settings,
    gene_coverage        = coverage,
    cells_per_sample     = cells_per_sample,
    by_patient_timepoint = patient_scores,
    by_cluster           = cluster_scores,
    by_timepoint         = timepoint_scores,
    lmm_coefficients     = lmm_coefficients,
    lmm_timepoint_lrt    = lrt_table,
    lmm_contrasts        = contrast_table,
    cross_sectional      = cross_sectional_table
  ))
  save_xlsx(sheets, paste0("functional_scores_summary", tag, ".xlsx"), out_dir)

  # per-subset scores differ between cell sets, so each set gets its own file
  if (scoring_mode == "per_subset") {
    write_checkpoint(obj, file = paste0("03_functional_scores", tag, ".rds"), dir = data_processed_dir,
                     workflow = "03_functional_scores", script = "01_functional_scores.R",
                     settings = list(cell_set = set_name, rule = selection$rule_text,
                                     programs = paste(program_names, collapse = ", "),
                                     scoring = scoring_mode, seed = params$scoring$seed),
                     input_path = input_path)
  }
}

# ------------------------------------------------------------------------------
# 5. Checkpoint
# ------------------------------------------------------------------------------
# every cell with its scores; later workflows choose their own cells. In
# "per_subset" mode the scores exist per cell set only (section 4).
if (scoring_mode == "all_cells") {
  write_checkpoint(
    scored_obj, "functional_scores", data_processed_dir,
    workflow = "03_functional_scores", script = "01_functional_scores.R",
    settings = list(programs = paste(program_names, collapse = ", "), scoring = scoring_mode,
                    seed = params$scoring$seed),
    input_path = input_path
  )
}

message("\nFunctional scores finished for cell set(s): ", paste(sets_to_run, collapse = ", "))
