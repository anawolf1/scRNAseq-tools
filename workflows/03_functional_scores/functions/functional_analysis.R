# ------------------------------------------------------------------------------
# functional_analysis.R
#
# Functional program scores: gene sets from .gmt files resolved to current
# human symbols, cell sets, module scores and net scores, summaries per
# group, the cross-sectional test, and the plots and tables of the scores.
# None of these save to disk; call save_plot() or save_gt() separately.
# Program ordering (order_programs_up_dn()) is in the shared gene_programs.R.
#
# Functions
# resolve_human_aliases()       outdated human gene symbols -> current HGNC symbols
# resolve_mouse_orthologs()     mouse gene symbols -> human orthologs
# prepare_gene_sets()           .gmt gene sets as human symbols, with a coverage table
# select_cell_set()             cells kept by a cell-set rule (keep/exclude by metadata)
# add_module_scores()           one module score per gene set
# net_program_bases()           programs with both an _UP and a _DN gene set
# add_net_scores()              net scores, <base>_UP - <base>_DN
# summarize_scores_by()         mean, median, SD and number of cells per group
# compare_groups_wilcox()       Wilcoxon rank-sum tests of patient means between groups
# plot_score_trajectories()     per-patient scores across repeated timepoints
# plot_program_umap_grid()      UMAP colored by each program's z-scored score
# parse_vs_comparison()         "<group1>_vs_<group2>" labels -> their two groups
# plot_group_significance()     score boxplots per program, with significance brackets
# gt_module_scores()            score summary as a gt table grouped by program
# gt_module_scores_wide()       score summary as a wide gt table, one spanner per group
# ------------------------------------------------------------------------------

#' Outdated human gene symbols -> current HGNC symbols
#'
#' limma::alias2SymbolTable() replaces aliases and retired symbols (e.g.
#' SEPP1 -> SELENOP); current and unknown symbols are kept. Two inputs can
#' map to the same symbol, so the result can contain duplicates.
#'
#' @param genes Human gene symbols.
#'
#' @return Character vector of the same length as `genes`.
#'
#' @examples
#' \dontrun{
#' resolve_human_aliases(c("SEPP1", "CD8A"))   # "SELENOP" "CD8A"
#' }
resolve_human_aliases <- function(genes) {
  dplyr::coalesce(limma::alias2SymbolTable(genes, species = "Hs"), genes)
}
#' Mouse gene symbols -> human orthologs
#'
#' babelgene::orthologs() (curated consensus of HGNC, HCOP and other
#' sources), so orthologs whose symbols differ by more than letter case are
#' found. A mouse gene with several human orthologs gives all of them;
#' genes without a human ortholog (e.g. Gm models, RIKEN clones, H2 or Klra
#' genes) are dropped.
#'
#' @param genes Mouse gene symbols.
#'
#' @return Unique human gene symbols.
#'
#' @examples
#' \dontrun{
#' resolve_mouse_orthologs(c("Tox", "Havcr2", "Gm12345"))
#' }
resolve_mouse_orthologs <- function(genes) {
  orth <- babelgene::orthologs(genes = genes, species = "mouse", human = FALSE)
  lookup <- split(orth$human_symbol, orth$symbol)
  mapped <- unlist(lapply(genes, function(g) lookup[[g]]), use.names = FALSE)
  unique(mapped[!is.na(mapped) & nzchar(mapped)])
}
#' .gmt gene sets as human symbols, with a coverage table
#'
#' Applies `symbol_fixes` (matched ignoring case), then resolves each gene
#' set by its species: human sets with resolve_human_aliases(), mouse sets
#' with resolve_mouse_orthologs(). A set with no gene left is dropped with a
#' warning.
#'
#' @param term2gene_raw term/gene table from read_gmt().
#' @param species Named vector, gene-set name -> "human" or "mouse"; needed
#'   for every set.
#' @param symbol_fixes Named vector, symbol in the .gmt -> replacement, e.g.
#'   c("4-SEP" = "SEPT4") for symbols damaged by spreadsheets.
#' @param object_genes Genes of the object to be scored (rownames(obj)).
#'
#' @return A list: `term2gene` (resolved term/gene table) and `coverage` (one
#'   row per gene set: species, genes in the .gmt, after resolution, found
#'   in the object, the fraction found, and the missing genes).
#'
#' @examples
#' \dontrun{
#' gene_sets <- prepare_gene_sets(read_gmt(gene_sets_files), params$gene_sets$species,
#'                                params$gene_sets$symbol_fixes, rownames(data_obj))
#' }
prepare_gene_sets <- function(term2gene_raw,
                              species,
                              symbol_fixes = character(0),
                              object_genes) {

  for (from in names(symbol_fixes)) {
    hit <- toupper(term2gene_raw$gene) == toupper(from)
    if (any(hit)) {
      term2gene_raw$gene[hit] <- symbol_fixes[[from]]
      message("Gene symbol fix: ", from, " -> ", symbol_fixes[[from]], " (", sum(hit), " entries)")
    }
  }

  terms <- unique(term2gene_raw$term)
  missing_species <- setdiff(terms, names(species))
  if (length(missing_species) > 0) {
    stop("prepare_gene_sets(): no species set for gene set(s): ", paste(missing_species, collapse = ", "),
         ". Add them to params$gene_sets$species.", call. = FALSE)
  }
  bad_species <- setdiff(unique(species[terms]), c("human", "mouse"))
  if (length(bad_species) > 0) {
    stop("prepare_gene_sets(): species must be \"human\" or \"mouse\", got: ",
         paste(bad_species, collapse = ", "), call. = FALSE)
  }

  resolved <- lapply(setNames(terms, terms), function(term_name) {
    raw_genes <- unique(term2gene_raw$gene[term2gene_raw$term == term_name])
    if (species[[term_name]] == "human") unique(resolve_human_aliases(raw_genes)) else resolve_mouse_orthologs(raw_genes)
  })
  empty <- terms[lengths(resolved) == 0]
  if (length(empty) > 0) {
    warning("prepare_gene_sets(): no genes left after alias/ortholog resolution, dropped: ",
            paste(empty, collapse = ", "), call. = FALSE)
  }

  term2gene <- dplyr::bind_rows(lapply(setdiff(terms, empty), function(term_name) {
    data.frame(term = term_name, gene = resolved[[term_name]], stringsAsFactors = FALSE)
  }))

  coverage <- dplyr::bind_rows(lapply(terms, function(term_name) {
    genes   <- resolved[[term_name]]
    missing <- setdiff(genes, object_genes)
    n_gmt   <- dplyr::n_distinct(term2gene_raw$gene[term2gene_raw$term == term_name])
    data.frame(
      term                  = term_name,
      species               = species[[term_name]],
      n_genes_gmt           = n_gmt,
      n_genes_resolved      = length(genes),
      n_genes_found         = length(genes) - length(missing),
      fraction_found_of_gmt = round((length(genes) - length(missing)) / n_gmt, 3),
      missing_genes         = paste(missing, collapse = ", ")
    )
  }))

  list(term2gene = term2gene, coverage = coverage)
}
#' Cells kept by a cell-set rule
#'
#' A cell set (config/decisions.R) has an optional `label` and optional
#' `keep` and `exclude` rules, each list(column = <metadata column>, values =
#' <values>). `keep` keeps only matching cells, `exclude` removes them; both
#' can be given; no rule keeps every cell. A value that does not occur in
#' the column stops the run, so a renamed label or a typo cannot match
#' nothing.
#'
#' @param meta Cell metadata (rows = cells).
#' @param spec The cell set's definition.
#' @param set_name Its name, for messages.
#'
#' @return A list: `keep` (logical, one per cell) and `rule_text` (the rules
#'   as text, "" if none).
select_cell_set <- function(meta,
                            spec,
                            set_name) {

  keep <- rep(TRUE, nrow(meta))
  rule_text <- character(0)

  for (rule_type in c("keep", "exclude")) {
    rule <- spec[[rule_type]]
    if (is.null(rule)) next

    if (is.null(rule$column) || length(rule$values) == 0) {
      stop("select_cell_set(): cell set '", set_name, "': `", rule_type, "` needs `column` and `values`.",
           call. = FALSE)
    }
    if (!rule$column %in% colnames(meta)) {
      stop("select_cell_set(): cell set '", set_name, "': column '", rule$column, "' not found in the metadata.",
           call. = FALSE)
    }
    column_values <- as.character(meta[[rule$column]])
    not_found <- setdiff(as.character(rule$values), unique(column_values))
    if (length(not_found) > 0) {
      stop("select_cell_set(): cell set '", set_name, "': value(s) not found in '", rule$column, "': ",
           paste(not_found, collapse = ", "), ". Values present: ",
           paste(utils::head(sort(unique(column_values)), 20), collapse = ", "), call. = FALSE)
    }

    hit <- column_values %in% as.character(rule$values)
    keep <- keep & (if (rule_type == "keep") hit else !hit)
    rule_text <- c(rule_text, paste0(rule_type, " ", rule$column, " in c(", paste(rule$values, collapse = ", "), ")"))
  }

  if (!any(keep)) stop("select_cell_set(): cell set '", set_name, "' keeps no cells.", call. = FALSE)
  list(keep = keep, rule_text = paste(rule_text, collapse = "; "))
}
#' One module score per gene set
#'
#' Seurat::AddModuleScore() once per gene set, stored as "<term>Score"
#' (AddModuleScore() itself names it "<term>1"). Pass the genes present in
#' the object (e.g. prepare_gene_sets()$term2gene filtered to rownames(obj)).
#'
#' @param obj Seurat object.
#' @param term2gene term/gene table; every gene must be in `obj`.
#' @param seed Seed for AddModuleScore(). Default 123.
#'
#' @return `obj` with one "<term>Score" metadata column per gene set.
#'
#' @examples
#' \dontrun{
#' obj <- add_module_scores(obj, term2gene_used, seed = 123)
#' }
add_module_scores <- function(obj,
                              term2gene,
                              seed = 123) {

  validate_required_cols(term2gene, c("term", "gene"), "add_module_scores")
  not_in_object <- setdiff(term2gene$gene, rownames(obj))
  if (length(not_in_object) > 0) {
    stop("add_module_scores(): ", length(not_in_object), " gene(s) not in the object, e.g. ",
         paste(utils::head(not_in_object, 5), collapse = ", "), "; pass only genes present.", call. = FALSE)
  }

  for (term_name in unique(term2gene$term)) {
    genes <- unique(term2gene$gene[term2gene$term == term_name])
    obj <- Seurat::AddModuleScore(obj, features = list(genes), name = term_name, seed = seed)

    generated_col <- paste0(term_name, "1")
    if (!generated_col %in% colnames(obj[[]])) {
      stop("add_module_scores(): AddModuleScore() did not create '", generated_col, "'; gene-set names ",
           "should contain only letters, digits, '_' and '.'.", call. = FALSE)
    }
    obj[[paste0(term_name, "Score")]] <- obj[[generated_col, drop = TRUE]]
    obj[[generated_col]] <- NULL
  }
  obj
}
#' Programs with both an _UP and a _DN gene set
#'
#' @param programs Program names (without "Score").
#'
#' @return The bases with both halves, e.g. "BB_EXH" for BB_EXH_UP and BB_EXH_DN.
#'
#' @examples
#' net_program_bases(c("BB_EXH_UP", "BB_EXH_DN", "MEM_UP"))   # "BB_EXH"
net_program_bases <- function(programs) {
  up_bases <- sub("_UP$", "", programs[grepl("_UP$", programs)])
  up_bases[paste0(up_bases, "_DN") %in% programs]
}
#' Net scores, <base>_UP - <base>_DN
#'
#' The per-cell difference of the two halves of a paired program. A shift
#' that raises both halves together (e.g. more genes detected in some
#' samples) cancels in the difference.
#'
#' @param obj Seurat object with "<base>_UPScore" and "<base>_DNScore".
#' @param bases Program bases, e.g. from net_program_bases().
#'
#' @return `obj` with one "<base>_NETScore" column per base; positive values
#'   lean toward the UP signature.
add_net_scores <- function(obj,
                           bases) {
  for (base in bases) {
    cols <- paste0(base, c("_UPScore", "_DNScore"))
    validate_required_cols(obj[[]], cols, "add_net_scores")
    obj[[paste0(base, "_NETScore")]] <- obj[[cols[1], drop = TRUE]] - obj[[cols[2], drop = TRUE]]
  }
  obj
}
#' Mean, median, SD and number of cells per group
#'
#' @param obj Seurat object with score columns.
#' @param score_cols Score columns, e.g. c("BB_EXH_UPScore", "BB_EXH_DNScore").
#' @param group_vars Grouping columns, e.g. c("Patient", "Timepoint").
#'
#' @return A data.frame in long format, one row per group and program:
#'   `group_vars`, `program` (score column without "Score"), `score` (the
#'   mean), `median`, `sd`, `n` (cells).
#'
#' @examples
#' \dontrun{
#' patient_scores <- summarize_scores_by(obj, score_cols, c("Patient", "Timepoint"))
#' }
summarize_scores_by <- function(obj,
                                score_cols,
                                group_vars) {

  meta <- obj[[]]
  validate_required_cols(meta, c(score_cols, group_vars), "summarize_scores_by")

  meta %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(group_vars))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(score_cols),
                    list(mean = ~ mean(.x, na.rm = TRUE),
                         median = ~ stats::median(.x, na.rm = TRUE),
                         sd = ~ stats::sd(.x, na.rm = TRUE)),
                    .names = "{.col}__{.fn}"),
      n = dplyr::n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_longer(cols = -c(dplyr::all_of(group_vars), n),
                        names_to = c("program", ".value"), names_sep = "__") %>%
    dplyr::mutate(program = sub("Score$", "", program)) %>%
    dplyr::rename(score = mean)
}
#' Wilcoxon rank-sum tests of patient means between groups
#'
#' For each pair of groups and each program, an unpaired two-sided
#' Wilcoxon rank-sum test of the patient means. Unpaired, so the two groups
#' must contain different patients (e.g. a repeated timepoint vs a
#' separately sampled group). With ties, wilcox.test() uses the normal
#' approximation. p-values are BH-adjusted across programs within each
#' comparison.
#'
#' @param scores_long Patient means: `Timepoint`, `program` and `score`
#'   (e.g. summarize_scores_by() output with the timepoint column renamed).
#' @param comparisons List of c(group1, group2) pairs of Timepoint values.
#' @param programs Programs to test.
#' @param min_patients Patients needed in each group for a test; fewer
#'   gives NA. Default 2.
#'
#' @return A data.frame: comparison ("<group1>_vs_<group2>"), program,
#'   n_group1, n_group2, median_group1, median_group2, p_value, p_adj,
#'   significance.
#'
#' @examples
#' \dontrun{
#' compare_groups_wilcox(patient_scores, list(c("FU", "SPR")), program_names)
#' }
compare_groups_wilcox <- function(scores_long,
                                  comparisons,
                                  programs,
                                  min_patients = 2) {

  validate_required_cols(scores_long, c("Timepoint", "program", "score"), "compare_groups_wilcox")

  rows <- list()
  for (pair in comparisons) {
    for (prog in programs) {
      x <- scores_long$score[scores_long$program == prog & scores_long$Timepoint == pair[1]]
      y <- scores_long$score[scores_long$program == prog & scores_long$Timepoint == pair[2]]
      rows[[length(rows) + 1]] <- data.frame(
        comparison    = paste0(pair[1], "_vs_", pair[2]),
        program       = prog,
        n_group1      = length(x),
        n_group2      = length(y),
        median_group1 = if (length(x) > 0) stats::median(x) else NA_real_,
        median_group2 = if (length(y) > 0) stats::median(y) else NA_real_,
        p_value       = if (length(x) >= min_patients && length(y) >= min_patients) {
          suppressWarnings(stats::wilcox.test(x, y)$p.value)
        } else {
          NA_real_
        }
      )
    }
  }

  dplyr::bind_rows(rows) %>%
    dplyr::group_by(comparison) %>%
    dplyr::mutate(p_adj = stats::p.adjust(p_value, method = "BH")) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(significance = p_to_stars(p_adj))
}
#' Per-patient scores across repeated timepoints
#'
#' One line per patient, the mean across patients in bold, one panel per
#' program. Optionally one cross-sectional group (e.g. "SPR") is shown to
#' the right, after a dotted line, as unconnected points (its patients are
#' different) with its mean as a diamond.
#'
#' @param scores_long Patient means: `Patient`, `Timepoint`, `program`,
#'   `score` (e.g. from summarize_scores_by()).
#' @param timepoint_levels The repeated timepoints, in order: the trajectory.
#' @param cross_sectional_group Optional single Timepoint value shown as a
#'   reference group. Default NULL.
#' @param program_labels Optional named vector, program -> panel label.
#' @param patient_col The dataset's patient column, for the patient colors
#'   (resolve_group_colors()). Default "Patient".
#' @param title,subtitle Plot title/subtitle. Default NULL.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_score_trajectories(patient_scores, timepoint_levels = c("BL", "3M", "FU"),
#'                         cross_sectional_group = "SPR")
#' }
plot_score_trajectories <- function(scores_long,
                                    timepoint_levels,
                                    cross_sectional_group = NULL,
                                    program_labels = NULL,
                                    patient_col = "Patient",
                                    title = NULL,
                                    subtitle = NULL) {

  validate_required_cols(scores_long, c("Patient", "Timepoint", "program", "score"), "plot_score_trajectories")
  if (!is.character(timepoint_levels) || length(timepoint_levels) == 0 ||
      anyNA(timepoint_levels) || anyDuplicated(timepoint_levels)) {
    stop("plot_score_trajectories(): `timepoint_levels` must be unique, non-missing timepoints.", call. = FALSE)
  }
  if (!is.null(cross_sectional_group) &&
      (length(cross_sectional_group) != 1 || cross_sectional_group %in% timepoint_levels)) {
    stop("plot_score_trajectories(): `cross_sectional_group` must be one timepoint not in `timepoint_levels`.",
         call. = FALSE)
  }

  all_levels <- c(timepoint_levels, cross_sectional_group)
  df <- scores_long[scores_long$Timepoint %in% all_levels, ]
  if (nrow(df) == 0) {
    stop("plot_score_trajectories(): none of the timepoints are in `scores_long`.", call. = FALSE)
  }
  df$Timepoint <- factor(as.character(df$Timepoint), levels = all_levels)

  traj_df <- df[df$Timepoint %in% timepoint_levels, ]
  ref_df  <- df[!df$Timepoint %in% timepoint_levels, ]
  traj_df$Patient <- factor(as.character(traj_df$Patient), levels = sort(unique(as.character(traj_df$Patient))))
  patient_colors <- resolve_group_colors(patient_col, levels(traj_df$Patient))

  mean_df <- traj_df %>%
    dplyr::group_by(program, Timepoint) %>%
    dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")

  p <- ggplot2::ggplot(traj_df, ggplot2::aes(x = Timepoint, y = score, group = Patient)) +
    ggplot2::geom_line(ggplot2::aes(color = Patient), linewidth = 0.7, alpha = 0.45) +
    ggplot2::geom_point(ggplot2::aes(color = Patient), size = 1.8, alpha = 0.7) +
    ggplot2::geom_line(data = mean_df, ggplot2::aes(x = Timepoint, y = mean_score, group = 1),
                       inherit.aes = FALSE, color = "#000000", linewidth = 1.5) +
    ggplot2::geom_point(data = mean_df, ggplot2::aes(x = Timepoint, y = mean_score),
                        inherit.aes = FALSE, color = "#000000", size = 3) +
    ggplot2::scale_color_manual(values = patient_colors) +
    ggplot2::scale_x_discrete(drop = FALSE)

  if (nrow(ref_df) > 0) {
    ref_mean_df <- ref_df %>%
      dplyr::group_by(program, Timepoint) %>%
      dplyr::summarise(mean_score = mean(score, na.rm = TRUE), .groups = "drop")
    p <- p +
      ggplot2::geom_vline(xintercept = length(timepoint_levels) + 0.5, linetype = "dotted", color = "#999999") +
      ggplot2::geom_point(data = ref_df, ggplot2::aes(x = Timepoint, y = score),
                          inherit.aes = FALSE, color = "#808080", size = 1.8, alpha = 0.7) +
      ggplot2::geom_point(data = ref_mean_df, ggplot2::aes(x = Timepoint, y = mean_score),
                          inherit.aes = FALSE, color = "#B22222", shape = 18, size = 4)
  }

  p +
    ggplot2::facet_wrap(~ program, scales = "free_y", ncol = 2,
                        labeller = if (!is.null(program_labels)) ggplot2::as_labeller(program_labels) else ggplot2::label_value) +
    ggplot2::labs(title = title, subtitle = subtitle, x = "Timepoint", y = "Mean module score", color = "Patient") +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title    = ggplot2::element_text(face = "bold", size = 15),
                   plot.subtitle = ggplot2::element_text(color = "#666666"),
                   strip.text    = ggplot2::element_text(face = "bold", size = 10))
}
#' UMAP colored by each program's z-scored score
#'
#' One panel per program. Scores are z-scored within each program across
#' the cells shown, so programs with different ranges share one color scale;
#' a color shows a cell's score relative to the other cells, not an absolute
#' level. Higher-scoring cells are drawn on top. A program with the same
#' score in every cell has no z-score and is drawn grey, with a warning.
#'
#' @param obj Seurat object with the score columns and `reduction`.
#' @param score_cols Score columns, one panel each.
#' @param reduction Reduction to plot. Default "umap".
#' @param program_labels Optional named vector, program -> panel label.
#' @param nrow,ncol Passed to facet_wrap().
#' @param point_size Point size. Default 0.8.
#' @param title,subtitle Plot title/subtitle. Default NULL.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_program_umap_grid(obj, c("BB_EXH_UPScore", "BB_EXH_DNScore"), ncol = 2)
#' }
plot_program_umap_grid <- function(obj,
                                   score_cols,
                                   reduction = "umap",
                                   program_labels = NULL,
                                   nrow = NULL,
                                   ncol = NULL,
                                   point_size = 0.8,
                                   title = NULL,
                                   subtitle = NULL) {

  if (!reduction %in% Seurat::Reductions(obj)) {
    stop("plot_program_umap_grid(): reduction '", reduction, "' not found. Available: ",
         paste(Seurat::Reductions(obj), collapse = ", "), call. = FALSE)
  }
  validate_required_cols(obj[[]], score_cols, "plot_program_umap_grid")

  embed <- Seurat::Embeddings(obj, reduction)[, 1:2]
  colnames(embed) <- c("Dim1", "Dim2")
  df <- cbind(as.data.frame(embed), obj[[]][, score_cols, drop = FALSE])
  program_of <- setNames(sub("Score$", "", score_cols), score_cols)

  df_long <- df %>%
    tidyr::pivot_longer(cols = dplyr::all_of(score_cols), names_to = "score_col", values_to = "score") %>%
    dplyr::mutate(program = program_of[score_col]) %>%
    dplyr::group_by(program) %>%
    dplyr::mutate(score_sd = stats::sd(score, na.rm = TRUE),
                  score_z = if (is.finite(score_sd[1]) && score_sd[1] > 0) as.numeric(scale(score)) else NA_real_) %>%
    dplyr::ungroup()

  constant_programs <- unique(df_long$program[!is.finite(df_long$score_sd) | df_long$score_sd == 0])
  if (length(constant_programs) > 0) {
    warning("plot_program_umap_grid(): no variation, drawn grey: ", paste(constant_programs, collapse = ", "),
            call. = FALSE)
  }
  if (!any(is.finite(df_long$score_z))) {
    stop("plot_program_umap_grid(): no program varies across these cells.", call. = FALSE)
  }
  z_limit <- max(abs(df_long$score_z), na.rm = TRUE)

  df_long$program <- factor(df_long$program, levels = unique(program_of))
  df_long <- df_long[order(df_long$program, df_long$score_z), ]
  module_score_pal <- resolve_palette("module_score_z")

  ggplot2::ggplot(df_long, ggplot2::aes(x = Dim1, y = Dim2, color = score_z)) +
    ggplot2::geom_point(size = point_size) +
    ggplot2::scale_color_gradient2(low = module_score_pal$low, mid = module_score_pal$mid,
                                   high = module_score_pal$high, midpoint = 0,
                                   limits = c(-z_limit, z_limit), na.value = "#D9D9D9",
                                   name = "Module score\n(z, within program)") +
    ggplot2::facet_wrap(~ program, nrow = nrow, ncol = ncol,
                        labeller = if (!is.null(program_labels)) ggplot2::as_labeller(program_labels) else ggplot2::label_value) +
    ggplot2::labs(title = title, subtitle = subtitle, x = paste0(reduction, "_1"), y = paste0(reduction, "_2")) +
    ggplot2::theme_classic(base_size = 13) +
    ggplot2::theme(plot.title        = ggplot2::element_text(face = "bold"),
                   plot.subtitle     = ggplot2::element_text(color = "#666666", size = 10),
                   strip.text        = ggplot2::element_text(face = "bold"),
                   axis.line         = ggplot2::element_line(color = "#4D4D4D"),
                   axis.ticks        = ggplot2::element_line(color = "#4D4D4D"),
                   axis.text         = ggplot2::element_blank(),
                   axis.ticks.length = grid::unit(0, "pt"))
}
#' "<group1>_vs_<group2>" labels -> their two groups
#'
#' Each label is matched against every ordered pair of `valid_levels` and
#' must match exactly one, so levels containing "_" or "_vs_" are split
#' correctly, and a typo stops with an error.
#'
#' @param comparison "<group1>_vs_<group2>" labels.
#' @param valid_levels The group levels the labels are built from.
#'
#' @return A data.frame with `group1` and `group2`, one row per label, in order.
#'
#' @examples
#' parse_vs_comparison(c("BL_vs_3M", "FU_vs_SPR"), c("BL", "3M", "FU", "SPR"))
parse_vs_comparison <- function(comparison,
                                valid_levels) {

  pairs <- expand.grid(group1 = valid_levels, group2 = valid_levels, stringsAsFactors = FALSE)
  pairs <- pairs[pairs$group1 != pairs$group2, ]
  pairs$candidate <- paste0(pairs$group1, "_vs_", pairs$group2)

  out <- lapply(comparison, function(cmp) {
    hit <- pairs[pairs$candidate == cmp, c("group1", "group2")]
    if (nrow(hit) != 1) {
      stop("parse_vs_comparison(): could not match '", cmp, "' to two of the levels: ",
           paste(valid_levels, collapse = ", "), call. = FALSE)
    }
    hit
  })
  do.call(rbind, out)
}
#' Score boxplots per program, with significance brackets
#'
#' One panel per program, one box per group, one point per unit (e.g. per
#' patient). Brackets show the adjusted p-values of `comparisons_table`; no
#' test is run here.
#'
#' @param scores_long Data.frame with `group_var`, `program` and `score`,
#'   one row per point.
#' @param comparisons_table Test results with `comparison`
#'   ("<group1>_vs_<group2>", levels of `group_levels`), `program` and `p_adj`.
#' @param group_var Column of `scores_long` with the groups (e.g. "Timepoint").
#' @param group_levels Group order on the x axis; at least 2.
#' @param group_colors Named colors, one per group. NULL (default) uses
#'   resolve_group_colors(group_var, group_levels).
#' @param program_order Optional panel order; other programs follow.
#' @param program_labels Optional named vector, program -> panel label.
#' @param title,subtitle Plot title/subtitle. Default NULL.
#' @param bracket_gap Vertical gap between stacked brackets, as a fraction
#'   of the panel's score range. Default 0.12.
#' @param nrow,ncol Passed to facet_wrap().
#' @param seed Seed for the jitter, so the figure is identical on every run.
#'   Default 123.
#'
#' @return A ggplot object.
#'
#' @examples
#' \dontrun{
#' plot_group_significance(patient_scores, contrast_table, group_var = "Timepoint",
#'                         group_levels = c("BL", "3M", "FU"))
#' }
plot_group_significance <- function(scores_long,
                                    comparisons_table,
                                    group_var,
                                    group_levels,
                                    group_colors = NULL,
                                    program_order = NULL,
                                    program_labels = NULL,
                                    title = NULL,
                                    subtitle = NULL,
                                    bracket_gap = 0.12,
                                    nrow = NULL,
                                    ncol = NULL,
                                    seed = 123) {

  validate_required_cols(scores_long, c(group_var, "program", "score"), "plot_group_significance")
  validate_required_cols(comparisons_table, c("comparison", "program", "p_adj"), "plot_group_significance",
                         hint = "expects a pairwise-comparison table.")
  if (!is.character(group_levels) || length(group_levels) < 2 || anyNA(group_levels) || anyDuplicated(group_levels)) {
    stop("plot_group_significance(): `group_levels` must be at least 2 unique, non-missing levels.", call. = FALSE)
  }

  plot_data <- scores_long
  plot_data[[group_var]] <- factor(plot_data[[group_var]], levels = group_levels)
  plot_data <- plot_data[!is.na(plot_data[[group_var]]), ]
  if (nrow(plot_data) == 0) {
    stop("plot_group_significance(): none of `group_levels` are in `scores_long$", group_var, "`.", call. = FALSE)
  }

  all_programs  <- unique(plot_data$program)
  program_order <- c(program_order, setdiff(all_programs, program_order))
  plot_data$program <- factor(plot_data$program, levels = program_order)
  group_colors <- group_colors %||% resolve_group_colors(group_var, group_levels)

  program_ranges <- plot_data %>%
    dplyr::group_by(program) %>%
    dplyr::summarise(max_score = max(score, na.rm = TRUE),
                     score_range = diff(range(score, na.rm = TRUE)), .groups = "drop")

  parsed <- parse_vs_comparison(comparisons_table$comparison, group_levels)
  bracket_df <- comparisons_table %>%
    dplyr::mutate(group1      = parsed$group1,
                  group2      = parsed$group2,
                  p.adj.label = paste0("p=", signif(p_adj, 2)),
                  program     = factor(program, levels = program_order)) %>%
    dplyr::filter(!is.na(program)) %>%
    dplyr::left_join(program_ranges, by = "program") %>%
    dplyr::group_by(program) %>%
    dplyr::arrange(group1, group2, .by_group = TRUE) %>%
    dplyr::mutate(y.position = max_score + score_range * bracket_gap * dplyr::row_number()) %>%
    dplyr::ungroup()

  # fill is set on the boxplot layer only: stat_pvalue_manual() would
  # otherwise inherit it, and the bracket table has no group column
  ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[group_var]], y = score)) +
    ggplot2::geom_boxplot(ggplot2::aes(fill = .data[[group_var]]), outlier.shape = NA, width = 0.5, alpha = 0.7) +
    ggplot2::geom_point(position = ggplot2::position_jitter(width = 0.08, height = 0, seed = seed),
                        size = 1, alpha = 0.4) +
    ggplot2::scale_fill_manual(values = group_colors) +
    ggplot2::facet_wrap(~ program, scales = "free_y", nrow = nrow, ncol = ncol, dir = "h",
                        labeller = if (!is.null(program_labels)) ggplot2::as_labeller(program_labels) else ggplot2::label_value) +
    ggpubr::stat_pvalue_manual(bracket_df, label = "p.adj.label", xmin = "group1", xmax = "group2",
                               y.position = "y.position", tip.length = 0.01, size = 3) +
    ggplot2::labs(title = title, subtitle = subtitle, y = "Score", x = NULL) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = "none")
}
#' Score summary as a gt table grouped by program
#'
#' @param scores_summary summarize_scores_by() output.
#' @param title,subtitle Table title/subtitle.
#'
#' @return A gt table, one row group per program.
gt_module_scores <- function(scores_summary,
                             title,
                             subtitle = "Mean ± SD module score per group") {

  validate_required_cols(scores_summary, c("program", "score", "median", "sd", "n"), "gt_module_scores")

  scores_summary %>%
    gt::gt(groupname_col = "program") %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::cols_label(score = "Mean", median = "Median", sd = "SD", n = "N cells") %>%
    gt::fmt_number(columns = c(score, median, sd), decimals = 3) %>%
    gt::cols_align(align = "center")
}
#' Score summary as a wide gt table, one spanner per group
#'
#' One row per program; per group a spanner with the mean and SD. A program
#' with no cells in a group shows an em dash.
#'
#' @param scores_summary summarize_scores_by() output for one grouping
#'   column (e.g. Timepoint or cluster_label).
#' @param group_var That grouping column.
#' @param title,subtitle Table title/subtitle.
#'
#' @return A gt table.
#'
#' @examples
#' \dontrun{
#' gt_module_scores_wide(summarize_scores_by(obj, score_cols, "Timepoint"),
#'                       group_var = "Timepoint", title = "Module scores by timepoint")
#' }
gt_module_scores_wide <- function(scores_summary,
                                  group_var,
                                  title,
                                  subtitle = NULL) {

  validate_required_cols(scores_summary, c(group_var, "program", "score", "sd"), "gt_module_scores_wide")
  group_levels <- unique(as.character(scores_summary[[group_var]]))

  wide_df <- scores_summary %>%
    dplyr::select(program, group_val = dplyr::all_of(group_var), score, sd) %>%
    tidyr::pivot_wider(names_from = group_val, values_from = c(score, sd), names_glue = "{group_val}_{.value}")

  label_list <- list(program = "Program")
  for (lvl in group_levels) {
    label_list[[paste0(lvl, "_score")]] <- "Mean"
    label_list[[paste0(lvl, "_sd")]]    <- "SD"
  }

  tbl <- wide_df %>%
    gt::gt() %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::fmt_number(columns = -program, decimals = 3) %>%
    gt::sub_missing(columns = dplyr::everything(), missing_text = "—") %>%
    gt::cols_label(.list = label_list)
  for (lvl in group_levels) {
    tbl <- tbl %>% gt::tab_spanner(label = lvl, columns = c(paste0(lvl, "_score"), paste0(lvl, "_sd")))
  }
  tbl %>% gt::cols_align(align = "center")
}
