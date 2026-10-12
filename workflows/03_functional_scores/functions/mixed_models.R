# ------------------------------------------------------------------------------
# mixed_models.R
#
# Linear mixed models of patient means across repeated timepoints: the fit
# with its likelihood-ratio test and pairwise contrasts, fixed effects as
# tables, and emmeans contrast labels in the "<group1>_vs_<group2>" form
# used by plot_group_significance(). None of these save to disk; call
# save_gt() separately.
#
# Functions
# fit_timepoint_lmm()        score ~ Timepoint + (1 | Patient): fixed effects, LRT, pairwise contrasts
# tidy_lmm()                 fixed effects of an lmer() model as a data.frame
# gt_lmm()                   tidy_lmm() output as a gt table, grouped by model
# emmeans_contrast_to_vs()   emmeans "A - B" labels -> "A_vs_B"
# ------------------------------------------------------------------------------

#' score ~ Timepoint + (1 | Patient): fixed effects, LRT and pairwise contrasts
#'
#' Fits the model to one program's patient means. The overall Timepoint
#' effect is a likelihood-ratio test of ML fits against score ~ 1 +
#' (1 | Patient); fixed effects and pairwise contrasts (emmeans,
#' BH-adjusted within the program) use the REML fit.
#'
#' @param data Patient means of one program: `score`, `Timepoint`, `Patient`.
#' @param timepoint_levels Timepoints in order; the first is the reference.
#' @param lmer_df Degrees-of-freedom method for the contrasts, passed to
#'   emmeans as lmer.df: "kenward-roger" (pbkrtest), "satterthwaite"
#'   (lmerTest) or "asymptotic".
#' @param program Program name stored in the results.
#'
#' @return A list: `coef` (tidy_lmm()), `lrt` (one row: program,
#'   n_patients, n_samples, chisq, df, p_value, singular_fit) and
#'   `contrasts` (program, comparison "<a>_vs_<b>", estimate, SE, df,
#'   statistic, p_adj, df_method).
#'
#' @examples
#' \dontrun{
#' fit_timepoint_lmm(prog_data, c("BL", "3M", "FU"), "kenward-roger", program = "BB_EXH_NET")
#' }
fit_timepoint_lmm <- function(data,
                              timepoint_levels,
                              lmer_df,
                              program) {

  validate_required_cols(data, c("score", "Timepoint", "Patient"), "fit_timepoint_lmm")
  data$Timepoint <- factor(as.character(data$Timepoint), levels = timepoint_levels)

  model_reml    <- lme4::lmer(score ~ Timepoint + (1 | Patient), data = data, REML = TRUE)
  model_full_ml <- lme4::lmer(score ~ Timepoint + (1 | Patient), data = data, REML = FALSE)
  model_null_ml <- lme4::lmer(score ~ 1 + (1 | Patient), data = data, REML = FALSE)
  lrt <- stats::anova(model_null_ml, model_full_ml)

  emm <- emmeans::emmeans(model_reml, ~ Timepoint, lmer.df = lmer_df)
  contrasts_df <- as.data.frame(emmeans::contrast(emm, method = "pairwise", adjust = "BH"))
  names(contrasts_df)[names(contrasts_df) %in% c("t.ratio", "z.ratio")] <- "statistic"

  list(
    coef = tidy_lmm(model_reml, model_label = program),
    lrt = data.frame(
      program      = program,
      n_patients   = dplyr::n_distinct(data$Patient),
      n_samples    = nrow(data),
      chisq        = lrt$Chisq[2],
      df           = lrt$Df[2],
      p_value      = lrt$`Pr(>Chisq)`[2],
      singular_fit = lme4::isSingular(model_reml)
    ),
    contrasts = data.frame(
      program    = program,
      comparison = emmeans_contrast_to_vs(contrasts_df$contrast, timepoint_levels),
      estimate   = contrasts_df$estimate,
      SE         = contrasts_df$SE,
      df         = contrasts_df$df,
      statistic  = contrasts_df$statistic,
      p_adj      = contrasts_df$p.value,
      df_method  = lmer_df
    )
  )
}
#' Fixed effects of an lmer() model as a data.frame
#'
#' One row per fixed-effect term, labeled with `model_label` so several
#' models can be combined with dplyr::bind_rows(). Values are not rounded;
#' gt_lmm() formats them for display. There is no p-value column: lme4 does
#' not compute one (contrast p-values come from emmeans).
#'
#' @param model Model fitted with lme4::lmer() or lmerTest::lmer().
#' @param model_label Label stored in the `Model` column.
#'
#' @return A data.frame: Model, Term, Estimate, SE, t_value.
#'
#' @examples
#' \dontrun{
#' coefs <- dplyr::bind_rows(
#'   tidy_lmm(model_exh, "BB_EXH_NET"),
#'   tidy_lmm(model_mem, "BB_MEM_NET")
#' )
#' }
tidy_lmm <- function(model,
                     model_label) {

  if (!inherits(model, "lmerMod")) {
    stop("tidy_lmm(): expects a model fitted with lme4::lmer(), got: ",
         paste(class(model), collapse = "/"), call. = FALSE)
  }

  coefs <- summary(model)$coefficients

  data.frame(
    Model    = model_label,
    Term     = rownames(coefs),
    Estimate = coefs[, "Estimate"],
    SE       = coefs[, "Std. Error"],
    t_value  = coefs[, "t value"],
    row.names = NULL
  )
}
#' tidy_lmm() output as a gt table, grouped by model
#'
#' Estimate, SE and t value rounded to three decimals, one row group per
#' model.
#'
#' @param lmm_df Data.frame from tidy_lmm(), or several of them combined with
#'   dplyr::bind_rows().
#' @param title,subtitle Table title/subtitle.
#'
#' @return A gt table.
#'
#' @examples
#' \dontrun{
#' gt_lmm(coefs, title = "LMM: score ~ Timepoint", subtitle = "Reference: BL")
#' }
gt_lmm <- function(lmm_df,
                   title,
                   subtitle = NULL) {

  validate_required_cols(lmm_df, c("Model", "Term", "Estimate", "SE", "t_value"), "gt_lmm",
                         hint = "expects tidy_lmm() output.")

  lmm_df %>%
    gt::gt(groupname_col = "Model") %>%
    gt::tab_header(title = title, subtitle = subtitle) %>%
    gt::cols_label(t_value = "t value") %>%
    gt::fmt_number(columns = c(Estimate, SE, t_value), decimals = 3) %>%
    gt::cols_align(align = "center")
}
#' emmeans contrast labels as "<group1>_vs_<group2>"
#'
#' emmeans labels pairwise contrasts "A - B", sometimes with a level in
#' parentheses ("BL - (3M)"). This returns "A_vs_B", keeping the direction,
#' and stops if a label does not split into two of `levels`. Levels that
#' contain " - " or parentheses cannot be parsed.
#'
#' @param contrast emmeans contrast labels (the `contrast` column of
#'   as.data.frame(emmeans::contrast(...))).
#' @param levels The factor levels the contrasts were built from.
#'
#' @return Character vector, one "<a>_vs_<b>" label per contrast.
#'
#' @examples
#' emmeans_contrast_to_vs(c("BL - 3M", "BL - (3M)"), c("BL", "3M", "FU"))
emmeans_contrast_to_vs <- function(contrast,
                                   levels) {

  clean <- gsub("[()]", "", as.character(contrast))
  parts <- strsplit(clean, " - ", fixed = TRUE)

  vapply(seq_along(parts), function(i) {
    p <- trimws(parts[[i]])
    if (length(p) != 2 || !all(p %in% levels)) {
      stop("emmeans_contrast_to_vs(): could not split '", contrast[i],
           "' into two of the levels: ", paste(levels, collapse = ", "), call. = FALSE)
    }
    paste0(p[1], "_vs_", p[2])
  }, character(1))
}
