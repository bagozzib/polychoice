## ---------------------------------------------------------------------------
## broom (tidyverse) integration for "choicemodel" fits; registered only when
## broom is installed (Suggests). modelsummary consumes these; texreg is
## served in zzz.R.
## ---------------------------------------------------------------------------

#' Tidy a choice model (broom method)
#'
#' @param x A `"choicemodel"` fit.
#' @param conf.int If `TRUE`, add `conf.low`/`conf.high` (see [confint.choicemodel()]).
#' @param conf.level Confidence level for the interval.
#' @param ... Unused.
#' @return A data frame with one row per parameter: `term`, `estimate`,
#'   `std.error`, `statistic`, `p.value`, and `component` (`"outcome"`,
#'   `"generic"`, `"fixed effect"`, `"inflation"`, `"nest"`, `"scale"`,
#'   `"covariance"`, or `"random intercept"`); inflation terms keep their
#'   `infl_` prefix so both equations' terms stay distinct in tables. The
#'   statistic and the intervals follow [summary.choicemodel()] and
#'   [confint.choicemodel()] (z against 1 for inclusive values and scales,
#'   intervals on the log scale for positive structural parameters, profile
#'   penalized-likelihood intervals for a penalized fit).
#' @examples
#' if (requireNamespace("broom", quietly = TRUE)) {
#'   set.seed(1)
#'   d <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
#'   m <- mnl(y ~ x1 + x2, data = d)
#'   broom::tidy(m, conf.int = TRUE)
#'   broom::glance(m)
#'   head(broom::augment(m))
#' }
#' @seealso [summary.choicemodel()]; `modelsummary::modelsummary()` consumes
#'   these methods, and `texreg::screenreg()` uses the package's `extract()`
#'   method; both carry the fit's `converged`/`boundary`/`ill_conditioned` flags.
#' @family broom methods
#' @exportS3Method broom::tidy
tidy.choicemodel <- function(x, conf.int = FALSE, conf.level = 0.95, ...) {
  s <- summary(x)
  ct <- s$coefficients
  out <- data.frame(term = rownames(ct), estimate = ct[, 1L], std.error = ct[, 2L], statistic = ct[, 3L],
                    p.value = ct[, 4L], component = s$block, row.names = NULL, stringsAsFactors = FALSE)
  if (conf.int && !is.null(x$vcov)) {
    ci <- stats::confint(x, level = conf.level)
    out$conf.low <- ci[, 1L]; out$conf.high <- ci[, 2L]
  }
  out
}

#' Glance at a choice model (broom method)
#'
#' @param x A `"choicemodel"` fit.
#' @param ... Unused.
#' @return A one-row data frame: `logLik`, `AIC`, `BIC`, `df`, `nobs`, for
#'   inflated models `inflate` and `share_inflated`, the fit-quality flags
#'   `converged`, `boundary`, `ill_conditioned`, and `n_modes` (see the section
#'   "Fit-quality flags" of [mnl()]), for fixed-effects fits
#'   `bias_correct` and `jackknife_refused` (`TRUE` when a requested jackknife
#'   correction was refused), and for penalized fits (`firth = TRUE` or a prior)
#'   `penalized = TRUE`: their `logLik` is the ordinary log-likelihood at the
#'   penalized estimates, and `AIC` and `BIC` are not maximum-likelihood
#'   criteria.
#' @family broom methods
#' @exportS3Method broom::glance
glance.choicemodel <- function(x, ...) {
  out <- data.frame(logLik = x$loglik, AIC = -2 * x$loglik + 2 * x$df, BIC = -2 * x$loglik + log(x$nobs_weighted) * x$df,
                    df = x$df, nobs = x$n, row.names = NULL)
  if (x$inflated) { out$inflate <- paste(x$inflate, collapse = ", "); out$share_inflated <- stats::weighted.mean(x$y %in% .pc_kc(x), x$weights) }
  out$converged <- isTRUE(x$converged); out$boundary <- isTRUE(x$boundary); out$ill_conditioned <- isTRUE(x$ill_conditioned)
  out$n_modes <- if (is.null(x$n_modes)) 1L else as.integer(x$n_modes)
  if (!is.null(x$fe)) { out$bias_correct <- if (is.null(x$bias_correct)) "none" else x$bias_correct; out$jackknife_refused <- !is.null(x$jackknife_refused) }
  if (.pc_penalized(x)) out$penalized <- TRUE                       # logLik, AIC, BIC at penalized estimates, not maximum-likelihood ones
  out
}

#' Augment the estimation cases with fitted classes and probabilities (broom method)
#'
#' @param x A `"choicemodel"` fit.
#' @param ... Unused.
#' @return The case-level estimation data (one row per case) with `.fitted`
#'   (the modal alternative), `.prob_<level>` columns, and for inflated models
#'   `.regime` (fitted \eqn{P(s = 1)}).
#' @family broom methods
#' @exportS3Method broom::augment
augment.choicemodel <- function(x, ...) {
  d <- x$design$data_case
  P <- x$fitted.values
  d$.fitted <- factor(x$levels[max.col(P, ties.method = "first")], levels = x$levels)   # from the same probabilities as the .prob_ columns
  colnames(P) <- paste0(".prob_", x$levels)
  d <- cbind(d, as.data.frame(P, row.names = rownames(d)))
  if (x$inflated) d$.regime <- unname(x$phi)
  d
}
