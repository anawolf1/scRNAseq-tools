# ------------------------------------------------------------------------------
# gene_programs.R
#
# Helpers for gene programs (gene sets scored as modules), shared by the
# functional scores and neighborhood analysis workflows.
#
# Functions
# order_programs_up_dn()   order paired _UP/_DN programs for tables and faceted plots
# ------------------------------------------------------------------------------

#' Order program names with each _UP directly before its _DN
#'
#' Programs sharing a base name are kept together, _UP first (e.g.
#' BB_EXH_UP, BB_EXH_DN, OTHER_UP, OTHER_DN). Other names follow in their
#' original order; none are dropped.
#'
#' @param programs Character vector of program names.
#'
#' @return The same names, reordered.
#'
#' @examples
#' order_programs_up_dn(c("BB_EXH_DN", "OTHER", "BB_EXH_UP", "OTHER_UP"))
order_programs_up_dn <- function(programs) {
  base   <- sub("_(UP|DN)$", "", programs)
  paired <- unlist(lapply(unique(base),
                          function(b) intersect(paste0(b, c("_UP", "_DN")), programs)))
  c(paired, setdiff(programs, paired))
}
