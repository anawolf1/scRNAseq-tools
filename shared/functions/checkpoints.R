# ------------------------------------------------------------------------------
# checkpoints.R
#
# Saving and locating pipeline checkpoints. Each workflow's config/config.R
# defines `checkpoints` (name -> file name); files live in data/processed/,
# which is the only link between workflows.
#
# Every checkpoint is written with a small plain-text manifest next to it
# (<file>_manifest.txt): which workflow and script made it, when, with which
# repo version, its settings, its input and the input's checksum, and the
# file's md5 checksum; for Seurat objects also the cell and gene counts,
# assays, reductions and metadata columns. A downstream workflow can read the
# manifest to know what it was given without reading the upstream workflow's
# config.
#
# Functions
# checkpoint_path()         build the full path to a registered pipeline checkpoint
# write_checkpoint()        save a checkpoint + its manifest
# read_manifest()           read a checkpoint's manifest as a named list
# check_checkpoint_input()  stop if a checkpoint's input changed since it was made
# ------------------------------------------------------------------------------

#' Build the full path to a pipeline checkpoint file
#'
#' @param name One of the names in `checkpoints` (e.g. "qc_filtered").
#' @param dir Directory the checkpoint lives in.
#'
#' @return Full file path (character).
#'
#' @examples
#' \dontrun{
#' obj <- readRDS(checkpoint_path("qc_filtered", data_processed_dir))
#' }
checkpoint_path <- function(name, dir) {
  if (!name %in% names(checkpoints)) {
    stop("checkpoint_path(): unknown checkpoint '", name,
         "'. Check `checkpoints` in this workflow's config/config.R.", call. = FALSE)
  }
  file.path(dir, checkpoints[[name]])
}
#' Save a checkpoint and write its manifest
#'
#' @param obj Object to save (usually a Seurat object).
#' @param name One of the names in `checkpoints`. Give this or `file`.
#' @param dir Directory to save into (usually data_processed_dir).
#' @param workflow Workflow folder name, e.g. "01_qc".
#' @param script Script that produced it, e.g. "02_qc_postfilter.R".
#' @param settings Named list of the settings that determined the result
#'   (e.g. the QC thresholds). Values are written with deparse(), one line
#'   each.
#' @param input_path File the checkpoint was made from; its checksum is
#'   recorded so later workflows can detect a changed input.
#' @param file File name in `dir` ending in .rds, for checkpoints whose name
#'   is not fixed (e.g. one per contrast). Give this or `name`.
#'
#' @return Invisibly, the checkpoint path.
#'
#' @examples
#' \dontrun{
#' write_checkpoint(data_raw, "qc_filtered", data_processed_dir,
#'                  workflow = "01_qc", script = "02_qc_postfilter.R",
#'                  settings = list(thresholds = params$thresholds))
#' }
write_checkpoint <- function(obj,
                             name = NULL,
                             dir,
                             workflow,
                             script,
                             settings = list(),
                             input_path = NULL,
                             file = NULL) {

  if (is.null(name) == is.null(file)) {
    stop("write_checkpoint(): give either `name` or `file`.", call. = FALSE)
  }
  path <- if (!is.null(file)) file.path(dir, file) else checkpoint_path(name, dir)
  
  if (!grepl("\\.rds$", path)) {
    stop("write_checkpoint(): the checkpoint file must end in .rds, got: ", basename(path), call. = FALSE)
  }
  saveRDS(obj, path)

  version_file <- here::here("VERSION")
  version <- if (file.exists(version_file)) trimws(readLines(version_file, warn = FALSE)[1]) else "unknown"

  lines <- c(
    paste0("checkpoint: ", basename(path)),
    paste0("workflow: ", workflow),
    paste0("script: ", script),
    paste0("created: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("repo_version: ", version),
    paste0("md5: ", unname(tools::md5sum(path)))
  )
  if (inherits(obj, "Seurat")) {
    lines <- c(lines,
               paste0("n_cells: ", ncol(obj)),
               paste0("n_genes: ", nrow(obj)),
               paste0("assays: ", paste(Seurat::Assays(obj), collapse = ", ")),
               paste0("reductions: ", paste(Seurat::Reductions(obj), collapse = ", ")),
               paste0("metadata: ", paste(colnames(obj[[]]), collapse = ", ")))
  }
  if (!is.null(input_path)) {
    lines <- c(lines,
               paste0("input: ", basename(input_path)),
               paste0("input_md5: ", unname(tools::md5sum(input_path))))
  }
  for (key in names(settings)) {
    lines <- c(lines, paste0("setting.", key, ": ", paste(deparse(settings[[key]], width.cutoff = 500L), collapse = " ")))
  }

  writeLines(lines, sub("\\.rds$", "_manifest.txt", path))
  message("Checkpoint saved: ", path)
  invisible(path)
}
#' Read a checkpoint's manifest
#'
#' @param path Path to the checkpoint (.rds) or to its manifest.
#'
#' @return Named list (key -> value, all character), or NULL if there is no
#'   manifest (e.g. an object not saved with write_checkpoint()).
read_manifest <- function(path) {
  manifest <- if (grepl("_manifest\\.txt$", path)) path else sub("\\.rds$", "_manifest.txt", path)
  if (!file.exists(manifest)) return(NULL)

  lines <- readLines(manifest, warn = FALSE)
  keys <- sub(":.*$", "", lines)
  values <- sub("^[^:]*: ?", "", lines)
  as.list(setNames(values, keys))
}
#' Stop if a checkpoint was made from an earlier version of its input
#'
#' Compares the input checksum recorded in the checkpoint's manifest with the
#' current input file: results built from an older input must not be
#' combined with the current one. Stops if they differ or if the input file
#' is missing. A checkpoint with no recorded checksum (no manifest, or saved
#' without `input_path`) is reported, not stopped on.
#'
#' @param path The checkpoint being read.
#' @param input_path The file it should have been made from.
#' @param rerun What to do about it, appended to the error.
#'
#' @return Invisibly TRUE (current) or NA (nothing recorded).
check_checkpoint_input <- function(path,
                                   input_path,
                                   rerun = NULL) {
  if (!file.exists(input_path)) {
    stop("check_checkpoint_input(): input not found: ", input_path, call. = FALSE)
  }
  recorded <- read_manifest(path)$input_md5
  if (is.null(recorded)) {
    message("[info] ", basename(path), " records no input checksum - cannot confirm it was made ",
            "from the current ", basename(input_path), ".")
    return(invisible(NA))
  }
  if (!identical(recorded, unname(tools::md5sum(input_path)))) {
    stop("check_checkpoint_input(): ", basename(path), " was made from an earlier version of ", basename(input_path),
         ". Combining them would mix results from different inputs. ", rerun %||% "", call. = FALSE)
  }
  invisible(TRUE)
}
