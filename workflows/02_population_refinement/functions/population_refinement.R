# ------------------------------------------------------------------------------
# population_refinement.R
#
# The round engine of pipeline/01_population_refinement.R: checkpoint paths,
# change detection, round inputs, and the steps that remove cells or genes.
# A step is one of three general operations (exclude_metadata,
# exclude_cells, exclude_genes); any metadata column (patient, plate,
# cluster, ...) can be used.
#
# Functions
# refinement_checkpoint_path()   path to a round's checkpoint
# refinement_config_path()       path to a round's recorded settings
# round_results_dir()            a round's results folder (round name + label)
# round_config_unchanged()       compare a round's settings with those recorded with its checkpoint
# refinement_input_path()        a round's input as a file path
# load_refinement_input()        load a round's input
# prepare_clustering_input()     load, subset, exclude and preprocess, for tools/
# input_round_name()             the round an input refers to, or NULL
# add_cluster_label_column()     cluster IDs -> labels, as a metadata column
# exclude_metadata()             remove cells whose metadata column has given values
# exclude_cells()                remove cells by barcode
# exclude_genes()                remove genes
# apply_batch_correction()       batch correction by method name
# validate_refinement_steps()    check a round's steps before anything is computed
# apply_refinement_steps()       run a round's steps on its input, with an audit table
# resolve_cells_file()           path of an exclude_cells step's cells_file
# resolve_step_cells()           barcodes of an exclude_cells step
# save_cell_list()               save barcodes to data/processed/cell_lists/
# ------------------------------------------------------------------------------

#' Path to a round's checkpoint
#'
#' @param round_name Round name, e.g. "round_1".
#' @param dir Directory the checkpoint lives in.
#'
#' @return File path.
refinement_checkpoint_path <- function(round_name,
                                       dir) {
  file.path(dir, paste0("02_population_refinement_", round_name, ".rds"))
}
#' Path to a round's recorded settings
#'
#' Saved next to each checkpoint, so a later run can tell whether the
#' round's settings changed since (see round_config_unchanged()).
#'
#' @param round_name Round name, e.g. "round_1".
#' @param dir Directory the checkpoint lives in.
#'
#' @return File path.
refinement_config_path <- function(round_name,
                                   dir) {
  file.path(dir, paste0("02_population_refinement_", round_name, "_config.rds"))
}
#' A round's results folder
#'
#' <base_dir>/<round>_<label>, with the label in lower case and
#' non-alphanumeric characters replaced by "_"; <base_dir>/<round> for a
#' round without a label. Used by the pipeline (to write) and by
#' tools/compare_rounds.R (to read).
#'
#' @param round_name Round name.
#' @param round The round's entry in config/decisions.R.
#' @param base_dir fig_refinement_dir.
#'
#' @return Folder path.
round_results_dir <- function(round_name,
                              round,
                              base_dir) {
  slug <- if (!is.null(round$label)) {
    paste0("_", gsub("^_|_$", "", tolower(gsub("[^a-zA-Z0-9]+", "_", round$label))))
  } else {
    ""
  }
  file.path(base_dir, paste0(round_name, slug))
}
#' Compare a round's settings with those recorded with its checkpoint
#'
#' Only the fields that change the result are compared: input,
#' input_md5 (checksum of the input file), input_labels (cluster labels of
#' the input round, which steps may select on), algorithm, dims,
#' resolution, k_param, prune_snn, batch_correction, steps, preprocessing
#' and exclude_genes. label, rationale, cluster_labels and the sweep
#' records are not. Field order does not matter.
#'
#' @param round The round's current settings, as completed by the pipeline.
#' @param round_name Round name.
#' @param dir Directory the checkpoint lives in.
#'
#' @return TRUE if unchanged, otherwise FALSE with attribute
#'   "changed_fields" naming what differs ("recorded settings" when the
#'   settings file is missing).
round_config_unchanged <- function(round,
                                   round_name,
                                   dir) {
  config_path <- refinement_config_path(round_name, dir)
  if (!file.exists(config_path)) {
    return(structure(FALSE, changed_fields = "recorded settings"))
  }

  result_fields <- c(
    "input", "input_md5", "input_labels", "algorithm", "dims", "resolution",
    "k_param", "prune_snn", "batch_correction", "steps", "preprocessing", "exclude_genes"
  )
  keep_result_fields <- function(r) r[sort(intersect(names(r), result_fields))]

  current  <- keep_result_fields(round)
  recorded <- keep_result_fields(readRDS(config_path))
  if (identical(current, recorded)) return(TRUE)

  fields  <- union(names(current), names(recorded))
  changed <- fields[!vapply(fields, function(f) identical(current[[f]], recorded[[f]]), logical(1))]
  structure(FALSE, changed_fields = changed)
}
#' A round's input as a file path
#'
#' A round's `input` is a checkpoint name from config.R's `checkpoints`
#' (e.g. "qc_filtered"), an earlier round's name (e.g. "round_1"), or a
#' file path (absolute, or relative to `dir`).
#'
#' @param name Checkpoint name, round name or file path.
#' @param dir Directory checkpoints live in; also the base of a relative path.
#'
#' @return File path. Stops if no file exists.
refinement_input_path <- function(name,
                                  dir) {
  path <- if (name %in% names(checkpoints)) {
    checkpoint_path(name, dir)
  } else {
    refinement_checkpoint_path(name, dir)
  }
  if (file.exists(path)) return(path)

  literal_path <- if (file.exists(name)) name else file.path(dir, name)
  if (!file.exists(literal_path)) {
    stop("refinement_input_path(): no file for '", name, "' (looked for ", path, " and ",
         literal_path, "). For a round name, check that the round has been run.", call. = FALSE)
  }
  literal_path
}
#' Load a round's input
#'
#' @inheritParams refinement_input_path
#'
#' @return The loaded object.
load_refinement_input <- function(name,
                                  dir) {
  readRDS(refinement_input_path(name, dir))
}
#' Load and preprocess an object for clustering-parameter screening
#'
#' The first step of tools/select_clustering_parameters.R and
#' tools/inspect_clustering_candidates.R: load an input, optionally keep a
#' subset of cells, run a round's steps, remove genes, and
#' preprocess as a round does. Parameters are only transferable to a round
#' that is preprocessed the same way, so pass the round's settings.
#'
#' A subset gets its own PCA: structure specific to those cells (e.g.
#' within one cluster) is then not diluted by the variation of all cells.
#'
#' @param input Checkpoint name, round name or file path (see
#'   refinement_input_path()).
#' @param dir Directory checkpoints live in.
#' @param subset_column,subset_values Keep only cells whose `subset_column`
#'   is in `subset_values`. Default NULL: all cells.
#' @param steps A round's steps. Default NULL.
#' @param exclude_gene_patterns Regex patterns of genes removed before
#'   preprocessing (params$clustering$exclude_gene_patterns for a round with
#'   exclude_genes = TRUE). Default NULL.
#' @param nfeatures,vars_to_regress,seed Passed to preprocess_seurat().
#' @param rounds params$population_refinement$rounds: when `input` is a
#'   round, its current cluster_labels are added. Default NULL.
#'
#' @return The preprocessed Seurat object (PCA computed, not clustered).
prepare_clustering_input <- function(input,
                                     dir,
                                     subset_column = NULL,
                                     subset_values = NULL,
                                     steps = NULL,
                                     exclude_gene_patterns = NULL,
                                     nfeatures = 2000,
                                     vars_to_regress = NULL,
                                     seed = 123,
                                     rounds = NULL) {

  obj <- load_refinement_input(input, dir)
  source_round <- input_round_name(input, names(rounds))
  if (!is.null(source_round)) obj <- add_cluster_label_column(obj, rounds[[source_round]]$cluster_labels)

  if (!is.null(subset_column)) {
    validate_required_cols(obj[[]], subset_column, "prepare_clustering_input")
    keep <- obj[[subset_column, drop = TRUE]] %in% subset_values
    if (!any(keep)) {
      stop("prepare_clustering_input(): no cells with ", subset_column, " in c(",
           paste(subset_values, collapse = ", "), ").", call. = FALSE)
    }
    obj <- subset(obj, cells = colnames(obj)[keep])
    message("Subset to ", subset_column, " in c(", paste(subset_values, collapse = ", "),
            "): ", ncol(obj), " cells.")
  }

  if (length(steps) > 0) {
    validate_refinement_steps(obj, steps, dir)
    step_result <- apply_refinement_steps(obj, steps, dir = dir)
    obj <- step_result$obj
    message("Steps removed ", length(step_result$removed_cells), " cells: ", ncol(obj), " remain.")
  }

  if (length(exclude_gene_patterns) > 0) {
    obj <- exclude_genes(obj, genes = paste(exclude_gene_patterns, collapse = "|"), pattern = TRUE)$obj
  }

  preprocess_seurat(obj, nfeatures = nfeatures, vars_to_regress = vars_to_regress,
                    show_elbow = FALSE, seed = seed)$obj
}
#' The round an input refers to
#'
#' @param name A round's `input`.
#' @param round_names Names of all rounds in config/decisions.R.
#'
#' @return The round name for an input like "round_1"; NULL for a
#'   checkpoint name or file path.
input_round_name <- function(name,
                             round_names) {
  if (name %in% round_names) name else NULL
}
#' Cluster IDs -> labels, as a metadata column
#'
#' @param obj Seurat object with a cluster ID column.
#' @param cluster_labels Named character vector (cluster ID -> label), e.g. a
#'   round's `cluster_labels`; NULL uses the IDs.
#' @param cluster_col Metadata column with the cluster IDs. Default
#'   "seurat_clusters".
#'
#' @return `obj` with a factor column `cluster_label`, levels in cluster-ID
#'   order. Clusters without a label keep their ID.
add_cluster_label_column <- function(obj,
                                     cluster_labels = NULL,
                                     cluster_col = "seurat_clusters") {

  validate_required_cols(obj[[]], cluster_col, "add_cluster_label_column")
  ids <- as.character(obj[[cluster_col, drop = TRUE]])
  labels <- ids
  if (!is.null(cluster_labels)) {
    mapped <- unname(cluster_labels[ids])
    labels[!is.na(mapped)] <- mapped[!is.na(mapped)]
  }
  id_order <- order(suppressWarnings(as.numeric(ids)), ids)
  obj$cluster_label <- factor(labels, levels = unique(labels[id_order]))
  obj
}
#' Remove cells whose metadata column has given values
#'
#' Any metadata column works: a cluster (`"seurat_clusters", "1"`), a
#' patient (`"Patient", "P01"`), a plate, and so on. Cells with NA in
#' `column` are kept unless `values` includes NA.
#'
#' @param obj Seurat object.
#' @param column Metadata column.
#' @param values Values of `column` to remove.
#'
#' @return A list: `obj` (without those cells) and `removed_cells` (their
#'   barcodes).
#'
#' @examples
#' \dontrun{
#' result <- exclude_metadata(obj, "seurat_clusters", "1")
#' result <- exclude_metadata(obj, "Patient", c("P01", "P02"))
#' }
exclude_metadata <- function(obj,
                             column,
                             values) {
  validate_required_cols(obj[[]], column, "exclude_metadata")
  match_idx <- obj[[column, drop = TRUE]] %in% values
  removed_cells <- colnames(obj)[match_idx]
  list(obj = subset(obj, cells = colnames(obj)[!match_idx]), removed_cells = removed_cells)
}
#' Remove cells by barcode
#'
#' For cells that no metadata column identifies, e.g. one subcluster of a
#' population, saved with save_cell_list().
#'
#' @param obj Seurat object.
#' @param cell_ids Barcodes to remove.
#'
#' @return A list: `obj` (without those cells) and `removed_cells` (the
#'   barcodes of `cell_ids` present in `obj`).
#'
#' @examples
#' \dontrun{
#' result <- exclude_cells(obj, c("CELL_A1", "CELL_B7"))
#' }
exclude_cells <- function(obj,
                          cell_ids) {
  removed_cells <- intersect(cell_ids, colnames(obj))
  list(obj = subset(obj, cells = setdiff(colnames(obj), removed_cells)), removed_cells = removed_cells)
}
#' Remove genes
#'
#' @param obj Seurat object.
#' @param genes Gene names, or one regular expression with `pattern = TRUE`.
#' @param pattern Logical. TRUE: `genes` is a regular expression matched
#'   against gene names (e.g. "^MT-|^RP[SL]"). Default FALSE.
#'
#' @return A list: `obj` (without those genes) and `removed_genes`.
#'
#' @examples
#' \dontrun{
#' result <- exclude_genes(obj, c("XIST", "MALAT1"))
#' result <- exclude_genes(obj, "^TR[ABDG][VDJ]", pattern = TRUE)
#' }
exclude_genes <- function(obj,
                          genes,
                          pattern = FALSE) {

  if (isTRUE(pattern) && length(genes) != 1) {
    stop("exclude_genes(): with `pattern = TRUE`, `genes` must be a single regular expression.",
         call. = FALSE)
  }

  all_genes <- rownames(obj)
  removed_genes <- if (pattern) all_genes[grepl(genes, all_genes)] else intersect(genes, all_genes)
  list(obj = subset(obj, features = setdiff(all_genes, removed_genes)), removed_genes = removed_genes)
}
#' Batch correction by method name
#'
#' Implemented: "harmony" (harmony::RunHarmony(), with `group.by` as its
#' group.by.vars; other harmony arguments are not passed). Other methods
#' are added as further branches.
#'
#' @param obj Seurat object with PCA computed.
#' @param method Method name, e.g. "harmony".
#' @param args Named list of the method's arguments, e.g.
#'   list(group.by = "Plate") for harmony.
#' @param seed Seed. Default 123.
#'
#' @return A list: `obj` (corrected) and `reduction` (the name of the
#'   corrected reduction, to cluster on).
apply_batch_correction <- function(obj,
                                   method,
                                   args = list(),
                                   seed = 123) {
  if (!is.character(method) || length(method) != 1 || is.na(method)) {
    stop("apply_batch_correction(): `method` must be a single method name.", call. = FALSE)
  }

  if (method == "harmony") {
    if (is.null(args$group.by)) {
      stop("apply_batch_correction(): harmony needs `group.by` (the batch column).", call. = FALSE)
    }
    validate_required_cols(obj[[]], args$group.by, "apply_batch_correction")
    set.seed(seed)
    obj <- harmony::RunHarmony(obj, group.by.vars = args$group.by)
    return(list(obj = obj, reduction = "harmony"))
  }

  stop("apply_batch_correction(): unsupported method '", method, "'.", call. = FALSE)
}
#' Check a round's steps before anything is computed
#'
#' Checks each step's type and fields (a field the step type does not use
#' stops, so a typo is not ignored), that a metadata column exists in the
#' input and that a cells_file exists, so a mistake stops the round before
#' preprocessing and clustering.
#'
#' @param obj The object the steps run on (the round's input).
#' @param steps The round's steps. NULL or an empty list passes.
#' @param dir data_processed_dir, where cells_file lists are looked up.
#'
#' @return Invisibly TRUE; stops otherwise.
validate_refinement_steps <- function(obj,
                                      steps,
                                      dir) {

  for (i in seq_along(steps)) {
    step <- steps[[i]]
    prefix <- paste0("validate_refinement_steps(): step ", i, " (", step$type %||% "no type", "): ")

    if (is.null(step$type) || !step$type %in% c("exclude_metadata", "exclude_cells", "exclude_genes")) {
      stop(prefix, "type must be \"exclude_metadata\", \"exclude_cells\" or \"exclude_genes\".", call. = FALSE)
    }

    step_fields <- c("type", switch(step$type,
      exclude_metadata = c("column", "values", "expect_n_cells"),
      exclude_cells    = c("cells_file", "cell_ids", "expect_n_cells"),
      exclude_genes    = c("genes", "pattern")
    ))
    unknown_fields <- setdiff(names(step), step_fields)
    if (length(unknown_fields) > 0) {
      stop(prefix, "unknown field(s): ", paste(unknown_fields, collapse = ", "), ". Valid fields: ",
           paste(step_fields, collapse = ", "), " (steps always run on the input, before clustering).",
           call. = FALSE)
    }

    if (!is.null(step$expect_n_cells)) {
      if (step$type == "exclude_genes") {
        stop(prefix, "`expect_n_cells` applies to steps that remove cells.", call. = FALSE)
      }
      if (!is.numeric(step$expect_n_cells) || length(step$expect_n_cells) != 1 ||
          !is.finite(step$expect_n_cells) || step$expect_n_cells < 0 ||
          step$expect_n_cells != round(step$expect_n_cells)) {
        stop(prefix, "`expect_n_cells` must be a single whole number.", call. = FALSE)
      }
    }

    if (step$type == "exclude_metadata") {
      if (!is.character(step$column) || length(step$column) != 1 || length(step$values) == 0) {
        stop(prefix, "needs one `column` name and non-empty `values`.", call. = FALSE)
      }
      if (!step$column %in% colnames(obj[[]])) {
        stop(prefix, "column '", step$column, "' not found in the input's metadata.", call. = FALSE)
      }
    }

    if (step$type == "exclude_cells") {
      has_ids  <- is.character(step$cell_ids) && length(step$cell_ids) > 0
      has_file <- is.character(step$cells_file) && length(step$cells_file) == 1
      if (has_ids == has_file) {
        stop(prefix, "needs exactly one of `cell_ids` (barcodes) or `cells_file` (a file of barcodes).",
             call. = FALSE)
      }
      if (has_file) resolve_cells_file(step$cells_file, dir)
    }

    if (step$type == "exclude_genes") {
      if (!is.character(step$genes) || length(step$genes) == 0) {
        stop(prefix, "needs a non-empty character vector `genes`.", call. = FALSE)
      }
      if (isTRUE(step$pattern) && length(step$genes) != 1) {
        stop(prefix, "with `pattern = TRUE`, `genes` must be a single regular expression.", call. = FALSE)
      }
    }
  }

  invisible(TRUE)
}
#' Run a round's steps on its input, in order
#'
#' Each step is a list with `type` ("exclude_metadata", "exclude_cells" or
#' "exclude_genes") and that function's arguments; `cells_file` can replace
#' `cell_ids` (see resolve_step_cells()). Each step runs on the previous
#' step's result. Steps run on the round's input, before preprocessing and
#' clustering, and select on the input's metadata (e.g. the previous
#' round's cluster_label).
#'
#' `expect_n_cells` (optional) is the number of cells the step should
#' remove; the round stops if it removes a different number. Cluster
#' numbers only stay the same while every input and setting does, so this
#' is recommended for cluster-based steps.
#'
#' @param obj Seurat object.
#' @param steps The round's steps. NULL or an empty list does nothing.
#' @param dir data_processed_dir, where cells_file lists are looked up.
#'
#' @return A list: `obj` (after the steps), `removed_cells` and
#'   `removed_genes` (over all steps run), and `audit` (one row per step
#'   run: type, target, cells and genes removed; NULL if there are none).
#'
#' @examples
#' \dontrun{
#' steps <- list(
#'   list(type = "exclude_metadata", column = "Patient", values = "P01"),
#'   list(type = "exclude_metadata", column = "cluster_label", values = "Low quality", expect_n_cells = 58)
#' )
#' result <- apply_refinement_steps(obj, steps, dir = data_processed_dir)
#' }
apply_refinement_steps <- function(obj,
                                   steps,
                                   dir) {
  removed_cells <- character(0)
  removed_genes <- character(0)
  audit_rows <- list()

  for (i in seq_along(steps)) {
    step <- steps[[i]]
    step_result <- switch(step$type,
      exclude_metadata = exclude_metadata(obj, step$column, step$values),
      exclude_cells    = exclude_cells(obj, resolve_step_cells(step, dir)),
      exclude_genes    = exclude_genes(obj, step$genes, pattern = isTRUE(step$pattern))
    )
    n_cells_removed <- length(step_result$removed_cells)
    if (!is.null(step$expect_n_cells) && n_cells_removed != step$expect_n_cells) {
      stop("apply_refinement_steps(): step ", i, " (", step$type, ") would remove ", n_cells_removed,
           " cells, but expect_n_cells = ", step$expect_n_cells, ". Cluster numbers or the input ",
           "may have changed since this step was decided.", call. = FALSE)
    }
    obj <- step_result$obj

    target <- switch(step$type,
      exclude_metadata = paste0(step$column, " = ", paste(step$values, collapse = ", ")),
      exclude_cells    = step$cells_file %||% paste(length(step$cell_ids), "cell IDs"),
      exclude_genes    = paste(step$genes, collapse = ", ")
    )
    audit_rows[[length(audit_rows) + 1]] <- data.frame(
      step = i, type = step$type, target = target,
      n_cells_removed = n_cells_removed,
      n_genes_removed = length(step_result$removed_genes)
    )
    removed_cells <- c(removed_cells, step_result$removed_cells)
    removed_genes <- c(removed_genes, step_result$removed_genes)
  }

  list(
    obj           = obj,
    removed_cells = unique(removed_cells),
    removed_genes = unique(removed_genes),
    audit         = if (length(audit_rows) > 0) do.call(rbind, audit_rows) else NULL
  )
}
#' Path of an exclude_cells step's cells_file
#'
#' Looked up in `dir`/cell_lists/ (where save_cell_list() writes), then in
#' `dir`, then as given (an absolute path, or relative to the working
#' directory).
#'
#' @param cells_file File name or path.
#' @param dir data_processed_dir.
#'
#' @return File path. Stops if the file is not found.
resolve_cells_file <- function(cells_file,
                               dir) {
  candidates <- c(file.path(dir, "cell_lists", cells_file), file.path(dir, cells_file), cells_file)
  path <- candidates[file.exists(candidates)][1]
  if (is.na(path)) {
    stop("resolve_cells_file(): cells_file '", cells_file, "' not found in ",
         file.path(dir, "cell_lists"), ", ", dir, " or as a path.", call. = FALSE)
  }
  path
}
#' Barcodes of an exclude_cells step
#'
#' From the step's `cell_ids`, or from its `cells_file` (one barcode per
#' line, e.g. written by save_cell_list(); see resolve_cells_file()).
#'
#' @param step The step.
#' @param dir data_processed_dir.
#'
#' @return Character vector of barcodes.
resolve_step_cells <- function(step,
                               dir) {
  if (!is.null(step$cell_ids)) return(as.character(step$cell_ids))
  cells <- trimws(readLines(resolve_cells_file(step$cells_file, dir), warn = FALSE))
  cells[nzchar(cells)]
}
#' Save barcodes for an exclude_cells step
#'
#' Writes one barcode per line to `dir`/cell_lists/<name>.txt, so a round
#' refers to the list by name (`cells_file = "<name>.txt"`).
#'
#' @param cells Barcodes (e.g. from Seurat::WhichCells()).
#' @param name File name without extension, e.g. "low_quality_half".
#' @param dir data_processed_dir.
#'
#' @return Invisibly, the file path.
#'
#' @examples
#' \dontrun{
#' cells <- Seurat::WhichCells(inspection_results[["dims8_res0.4"]]$obj, idents = "1")
#' save_cell_list(cells, "low_quality_half", data_processed_dir)
#' }
save_cell_list <- function(cells,
                           name,
                           dir) {
  if (length(cells) == 0) stop("save_cell_list(): no cells given.", call. = FALSE)
  out_dir <- file.path(dir, "cell_lists")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(out_dir, paste0(name, ".txt"))
  writeLines(as.character(cells), path)
  message("Saved ", length(cells), " cells to ", path)
  invisible(path)
}
