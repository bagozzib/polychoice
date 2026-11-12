## ---------------------------------------------------------------------------
## Priors on the coefficients (`prior =` in mnl() and rologit()): the
## log-likelihood penalized by the log prior density of each coefficient,
## normal, Cauchy or Student t centred at zero, on coefficients scaled as
## Gelman et al. (2008) and arm::bayesglm scale them. Definitions fixed in
## DESIGN.md, section 1, Phase 4c. The penalized machinery of Firth's penalty
## (two starts, scoring, profile intervals, penalized likelihood ratios)
## serves both (firth.R).
## ---------------------------------------------------------------------------

#' Prior distributions for the coefficients of a penalized fit
#'
#' Specifies independent priors centred at zero for the coefficients of
#' [mnl()] and [rologit()] (`prior = coef_prior(...)`): the fit maximizes the
#' log-likelihood plus the log prior density, the posterior mode. Under
#' separation the maximum-likelihood estimates do not exist, and the prior (like
#' Firth's penalty, `firth = TRUE`, which is Jeffreys' prior) decides how large
#' the separated estimates are (Rainey 2016). Jeffreys' prior is a default
#' justified by bias reduction; it depends on the design and is not independent
#' across coefficients, and can move an estimate against the direction of the
#' separation (Beiser-McGrath 2022). Independent priors make the assumption
#' explicit, and fitting under several (a skeptical, an informative and an
#' enthusiastic scale, Firth's penalty, and Gelman et al.'s default) shows how
#' much of a separated estimate comes from the prior.
#'
#' @details
#' The defaults are Gelman et al.'s (2008): Cauchy priors with scale 2.5 for
#' the slopes and 10 for the constants, on coefficients scaled to the input.
#' With `autoscale = TRUE` the scale of coefficient \eqn{k} is divided by
#' \eqn{c_k}: twice the standard deviation of its column when the column takes
#' more than two values, the difference of its two values when it takes two
#' (a dummy), as `arm::bayesglm` does. The column scales are taken over the
#' estimation cases for a case-specific covariate and over the available
#' alternatives for an alternative-specific (generic) one (the spread
#' across all available alternatives, between cases as well as within them;
#' Gelman et al.'s inputs are case-level, and the conditional logit's
#' within-case spread would be the closer analogue), with the frequency
#' weights (a weighted fit and the fit to the rows repeated by their weights
#' get the same prior). The constants (the intercept of the case-specific
#' part, which gives the alternative-specific constants) take
#' `intercept_scale`, unscaled. The covariates are not centred, unlike Gelman et
#' al.'s inputs, so the constants' prior concerns the log odds where every
#' covariate is zero: centre a covariate whose zero lies far from its values
#' (a calendar year) for that prior to mean what it says. With scale 2.5 a
#' normal prior puts 95 percent of its mass on log odds within about 5 per two
#' standard deviations of a covariate; the default Cauchy, with its heavier
#' tails, puts 30 percent beyond 5. The prior's normalizing constants are left
#' out of the penalized log-likelihood (they do not depend on the
#' coefficients).
#'
#' The fit reports the log-likelihood at the posterior mode (`loglik`), the
#' maximized objective (`loglik_penalized`), and the inverse of the curvature
#' of the log posterior at the mode as the covariance (the model's information
#' plus the prior's curvature). `arm::bayesglm` differs for t and Cauchy
#' priors: its approximate EM algorithm reports a related, always positive
#' definite matrix at its own fixed point, and the two agree for the normal
#' prior. For a t or Cauchy prior far out in its tail the exact curvature can
#' fail to be positive definite (a separated coefficient); the prior's negative
#' curvature is then left out, with a warning. Profile penalized-likelihood intervals ([confint()]) and penalized
#' likelihood-ratio tests ([term_test()], [combine_test()], [lr_test()] between fits with the
#' same prior, [rank_test()]) are the recommended inference under separation,
#' as for Firth's penalty.
#'
#' @param family `"cauchy"` (default), `"normal"`, or `"t"` (Student t with
#'   `df` degrees of freedom).
#' @param scale The prior scale of the slopes (default 2.5). A named vector sets
#'   the scale of the coefficients it names (by their names in `coef()`), and
#'   an unnamed element, if any, the default for the others: a separated
#'   coefficient can be given a prior of its own (Rainey's partial prior).
#'   With `autoscale`, a named scale is divided by its column's scale like the
#'   default. A name the model does not have is ignored with a message, so the
#'   restricted model of a penalized likelihood-ratio test ([lr_test()]) can
#'   take the full model's prior.
#' @param df Degrees of freedom of the t prior (`family = "t"` only).
#' @param intercept_scale The prior scale of the constants (default 10).
#' @param autoscale If `TRUE` (default), the slopes' scales are relative to
#'   their columns (see Details); if `FALSE`, they are on the coefficients as
#'   they are. The column scales use the frequency weights, so a fit whose
#'   weights total at most 2 (normalized or survey weights, which are not
#'   frequency weights) is refused under `autoscale = TRUE`: rescale the
#'   weights to counts, or set `autoscale = FALSE`. A scale named for a
#'   coefficient the model lacks is ignored with a message, which suggests a
#'   near match when the name looks like a misspelling.
#' @return An object of class `"coef_prior"`.
#' @references
#' Beiser-McGrath, L.F. (2022). Separation and rare events. *Political Science
#' Research and Methods*, 10, 428-437. \doi{10.1017/psrm.2020.46}
#'
#' Gelman, A., Jakulin, A., Pittau, M.G. and Su, Y.-S. (2008). A weakly
#' informative default prior distribution for logistic and other regression
#' models. *Annals of Applied Statistics*, 2, 1360-1383. \doi{10.1214/08-AOAS191}
#'
#' Rainey, C. (2016). Dealing with separation in logistic regression models.
#' *Political Analysis*, 24, 339-355. \doi{10.1093/pan/mpw014}
#' @examples
#' data(defection)
#' ## the birth year centred, so that the constants' prior concerns typical respondents
#' defection$birth_c <- defection$birth_year - round(mean(defection$birth_year, na.rm = TRUE))
#' f <- outcome ~ economic * indoctrination + age_joined + birth_c + school_age + male
#' mc <- mnl(f, data = defection, ref = "captured", prior = coef_prior())
#' mn <- mnl(f, data = defection, ref = "captured", prior = coef_prior("normal", scale = 4.5))
#' cbind(cauchy = coef(mc), normal = coef(mn))
#' @seealso [mnl()], [rologit()], [separation()]
#' @export
coef_prior <- function(family = c("cauchy", "normal", "t"), scale = 2.5, df = NULL, intercept_scale = 10, autoscale = TRUE) {
  family <- match.arg(family)
  if (family == "t") {
    if (is.null(df) || !is.numeric(df) || length(df) != 1L || !is.finite(df) || df <= 0) stop("A t prior needs 'df', one positive number.", call. = FALSE)
    df <- as.numeric(df)                       # 7L and 7 are one prior (local review 2026-11-06)
  } else {
    if (!is.null(df)) stop("'df' is for family = \"t\" only.", call. = FALSE)
    df <- if (family == "cauchy") 1 else Inf
  }
  if (!is.numeric(scale) || !length(scale) || any(!is.finite(scale)) || any(scale <= 0)) stop("'scale' must be positive numbers.", call. = FALSE)
  nm <- names(scale)
  if (length(scale) > 1L && (is.null(nm) || sum(!nzchar(nm)) > 1L))
    stop("'scale' is one number, or a named vector (with at most one unnamed element, the default).", call. = FALSE)
  if (!is.null(nm) && anyDuplicated(nm[nzchar(nm)])) stop("'scale' names a coefficient twice.", call. = FALSE)
  if (!is.numeric(intercept_scale) || length(intercept_scale) != 1L || !is.finite(intercept_scale) || intercept_scale <= 0)
    stop("'intercept_scale' must be one positive number.", call. = FALSE)
  if (!is.logical(autoscale) || length(autoscale) != 1L || is.na(autoscale)) stop("'autoscale' must be TRUE or FALSE.", call. = FALSE)
  structure(list(family = family, scale = scale, df = df, intercept_scale = intercept_scale, autoscale = autoscale), class = "coef_prior")
}

#' @export
print.coef_prior <- function(x, ...) {
  cat(.pc_penalty_label(list(prior = x)), "; constants: scale ", format(x$intercept_scale), "\n", sep = "")
  invisible(x)
}

## the column scale c_k of a covariate's values (DESIGN.md Phases 4c and 4d)
## over the rows of positive weight: the difference of its two values, twice
## its frequency-weighted standard deviation when it takes more (the standard
## deviation of the rows repeated by their weights), 1 when it is constant
.pc_prior_colscale <- function(v, w = rep(1, length(v))) {
  ok <- is.finite(v) & w > 0; v <- v[ok]; w <- w[ok]
  u <- unique(v)
  if (length(u) <= 1L) 1 else if (length(u) == 2L) abs(diff(u)) else {
    m <- sum(w * v) / sum(w)
    2 * sqrt(sum(w * (v - m)^2) / (sum(w) - 1))
  }
}

## The prior resolved on a design: the scale s_k of every coefficient in theta
## units, the constants marked. In mnl() and rologit() designs theta holds the
## case-specific block (p columns of X, one block per non-base alternative,
## D$idx$beta) and then the generic columns of W.
.pc_prior_resolve <- function(prior, D) {
  if (!inherits(prior, "coef_prior")) stop("'prior' must be made by coef_prior().", call. = FALSE)
  ## the autoscaled columns' standard deviations are those of the rows
  ## repeated by their frequency weights, which needs a weight total above 2
  ## (normalized survey weights, summing to 1, are not frequency weights)
  if (isTRUE(prior$autoscale) && sum(D$w) <= 2)
    stop("The weights total ", format(sum(D$w)), ": a prior's autoscaling reads them as frequency weights (the rows repeated), which needs a total above 2. ",
         "Normalized or survey weights are not supported; rescale frequency weights to counts, or give coef_prior(autoscale = FALSE).", call. = FALSE)
  nm <- names(prior$scale)
  default <- if (is.null(nm)) prior$scale[1L] else if (any(!nzchar(nm))) prior$scale[!nzchar(nm)][1L] else 2.5
  over <- if (is.null(nm)) numeric(0) else prior$scale[nzchar(nm)]
  ## a scale named for a coefficient the model lacks is ignored with a message
  ## (DESIGN.md Phase 4d): the restricted model of a penalized likelihood
  ## ratio takes the full model's prior, whose named scales may concern the
  ## coefficients it drops
  bad <- setdiff(names(over), D$par_names)
  if (length(bad)) {
    ## a near match (a case or one- or two-character difference) is likely a
    ## typo, unless the name's term is a variable of the data that this model
    ## leaves out (a restricted model under the full model's prior; round 7,
    ## D-5), where the hint would point to another coefficient
    vars <- names(D$data)
    near <- vapply(bad, function(b) {
      case_only <- tolower(b) %in% tolower(D$par_names)
      term <- if (sub(":.*$", "", b) %in% D$levels) sub("^[^:]+:", "", b) else b
      dropped <- all(strsplit(term, ":", fixed = TRUE)[[1L]] %in% vars)
      d <- utils::adist(tolower(b), tolower(D$par_names))[1L, ]
      if (case_only || (!dropped && min(d) <= 2L)) paste0(b, " (did you mean ", D$par_names[which.min(d)], "?)") else b }, character(1))
    message("'scale' names coefficient(s) the model does not have: ", paste(near, collapse = ", "), "; their scales are ignored (the names are those of coef()).")
    over <- over[setdiff(names(over), bad)]
  }
  s <- rep(default, D$K); const <- rep(FALSE, D$K); cc <- rep(1, D$K)
  if (D$p) {
    xs <- vapply(seq_len(D$p), function(j) .pc_prior_colscale(D$X[, j], D$w), numeric(1))
    isc <- colnames(D$X) == "(Intercept)"
    for (m in seq_len(D$J - 1L)) {
      ix <- D$idx$beta[, m]
      cc[ix] <- ifelse(isc, 1, xs); const[ix] <- isc
    }
  }
  if (D$q) {
    av <- as.vector(t(D$avail)) == 1L
    Wa <- D$W[av, , drop = FALSE]; wa <- rep(D$w, each = D$J)[av]
    ix <- D$p * (D$J - 1L) + seq_len(D$q)
    cc[ix] <- vapply(seq_len(D$q), function(j) .pc_prior_colscale(Wa[, j], wa), numeric(1))
  }
  s[const] <- prior$intercept_scale
  if (length(over)) s[match(names(over), D$par_names)] <- over
  if (!prior$autoscale) cc[!const] <- 1
  sk <- s / ifelse(const, 1, cc)
  list(family = prior$family, df = prior$df, scale = prior$scale, intercept_scale = prior$intercept_scale, autoscale = prior$autoscale,
       s = stats::setNames(sk, D$par_names), column_scale = stats::setNames(ifelse(const, NA_real_, cc), D$par_names), constant = const,
       spec = prior)
}

## The penalized objective of a design with a prior, in the engine's form (as
## .pc_firth_fns): nll the negative of the log-likelihood plus the log prior,
## grad its gradient, scores the per-case scores plus an equal share of the
## prior's gradient per unit of weight (so that, weighted, they sum to zero at
## the mode), `info` the curvature of the scoring (the model's information
## plus the prior's where positive), and `vcov_info` the curvature of the log
## posterior (the model's information plus the prior's), which the covariance
## inverts where it is positive definite (.pc_vcov).
.pc_prior_fns <- function(f, D) {
  pr <- D$prior; s <- unname(pr$s); df <- pr$df; normal <- identical(pr$family, "normal")
  ro <- identical(D$family, "rologit")
  a <- list(D$X, D$W, D$off, D$avail, if (ro) D$str$rank else D$avail, if (ro) 1L else 0L, D$J, D$base - 1L,
            if (ro) as.integer(D$str$split %||% 0L) else 0L)
  w <- D$w; sw <- sum(w)
  memo <- new.env(parent = emptyenv())
  info0 <- function(th) {
    th <- as.numeric(th)
    if (!is.null(memo$th) && identical(memo$th, th)) return(memo$I)
    I <- do.call(pc_firth_info_cpp, c(list(th), a, list(w)))
    memo$th <- th; memo$I <- I
    I
  }
  ## a split design (rank_test) repeats the coefficients: the prior applies to each copy
  sk <- function(th) rep_len(s, length(th))
  lp <- function(th) { z <- th / sk(th); if (normal) -0.5 * sum(z^2) else -(df + 1) / 2 * sum(log1p(z^2 / df)) }
  glp <- function(th) { ss <- sk(th); z <- th / ss; if (normal) -z / ss else -(df + 1) * z / (ss * (df + z^2)) }
  hlp <- function(th) { ss <- sk(th); z <- th / ss; if (normal) 1 / ss^2 else (df + 1) * (df - z^2) / (ss^2 * (df + z^2)^2) }
  list(lli = f$lli,
       nll = function(th) f$nll(th) - lp(th),
       grad = function(th) f$grad(th) - glp(th),
       hess = NULL,
       info = function(th) info0(th) + diag(pmax(hlp(th), 0), length(th)),
       vcov_info = function(th) info0(th) + diag(hlp(th), length(th)),
       scores = function(th) { S <- f$scores(th); S + matrix(glp(th) / sw, nrow(S), ncol(S), byrow = TRUE) },
       probs = f$probs, nll0 = f$nll, penalty = lp)
}
