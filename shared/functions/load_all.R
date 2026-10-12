# ------------------------------------------------------------------------------
# load_all.R
#
# Defines source_folder() and uses it to source every file in
# shared/functions/. Each workflow's setup.R then calls source_folder() on
# its own functions/ folder, so there is one loader for the whole repo
# rather than a copy in every workflow.
#
# Uses here::here(): the repo has its own .Rproj at the top level, so
# here() always resolves to the repo root regardless of which script calls it.
#
# Usage:
#   source(here::here("shared", "functions", "load_all.R"))
#   source_folder(here::here("workflows", "01_qc", "functions"))
# ------------------------------------------------------------------------------

#' Source every .R file in a folder
#'
#' Files are sourced in alphabetical order (list.files()' default order is
#' OS-dependent). load_all.R itself is skipped. A missing folder is not an
#' error, so a workflow without its own functions/ folder still works.
#'
#' @param dir Folder to source.
#'
#' @return Invisibly, the files sourced.
source_folder <- function(dir) {
  if (!dir.exists(dir)) return(invisible(character(0)))

  files <- sort(list.files(dir, pattern = "\\.R$", full.names = TRUE), method = "radix")
  files <- files[basename(files) != "load_all.R"]
  for (f in files) {
    tryCatch(
      source(f),
      error = function(e) stop("source_folder(): failed to load ", basename(f), " (", dir, "): ",
                               conditionMessage(e), call. = FALSE)
    )
  }
}

# null-default operator used throughout; built into R >= 4.4, defined here for older R
if (!exists("%||%", envir = baseenv())) {
  `%||%` <- function(x, y) if (is.null(x)) y else x
}

source_folder(here::here("shared", "functions"))
