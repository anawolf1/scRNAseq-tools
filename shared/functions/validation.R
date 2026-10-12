# ------------------------------------------------------------------------------
# validation.R
#
# Generic input-validation helpers shared across function files.
#
# Functions
# validate_table_dims()          dim()-length check
# validate_required_cols()       required-column check
# ------------------------------------------------------------------------------

#' Validate that an object has the expected number of dimensions
#'
#' Used by table and heatmap helpers that need a 2- or 3-way table.
#'
#' @param x The object to check.
#' @param expected_dims Integer, the required length of dim(x).
#' @param fn_name Calling function's name (for the error message prefix).
#' @param hint Optional extra sentence appended to the wrong-dimensions
#'   message, e.g. pointing to a companion function that handles that case.
#'
#' @return Invisibly TRUE if valid; stops otherwise.
validate_table_dims <- function(x,
                                expected_dims,
                                fn_name,
                                hint = NULL) {

  if (is.null(dim(x))) {
    stop(fn_name, "() expects a ", expected_dims, "-dimensional table()/matrix. ",
         "Got an object with no dim() at all (class: ", paste(class(x), collapse = "/"), ").", call. = FALSE)
  }

  if (length(dim(x)) != expected_dims) {
    stop(fn_name, "() expects a ", expected_dims, "-dimensional table()/matrix. ",
         "Got an object with ", length(dim(x)), " dimensions: ", paste(dim(x), collapse = " x "), ".",
         if (!is.null(hint)) paste0(" ", hint) else "",
         call. = FALSE)
  }

  invisible(TRUE)
}
#' Validate that a data.frame has all required columns
#'
#' Used throughout the repository to stop early, with a clear message, when
#' an input lacks a column a function relies on.
#'
#' @param df The data.frame to check.
#' @param required_cols Character vector of column names that must be present.
#' @param fn_name Calling function's name (for the error message prefix).
#' @param hint Optional sentence describing the expected input shape,
#'   inserted between the function name and "Missing column(s): ...".
#'
#' @return Invisibly TRUE if valid; stops otherwise.
validate_required_cols <- function(df,
                                   required_cols,
                                   fn_name,
                                   hint = NULL) {

  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    prefix <- if (!is.null(hint)) paste0(fn_name, "() ", hint, " ") else paste0(fn_name, "(): ")
    stop(prefix, "Missing column(s): ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }

  invisible(TRUE)
}
