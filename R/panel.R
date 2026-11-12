## ---------------------------------------------------------------------------
## Panel tools that compose with every estimator: the Mundlak / correlated-
## random-effects device. Cluster-robust standard errors and unit fixed
## effects (`fe`) are built into the estimators.
## ---------------------------------------------------------------------------

#' Mundlak (correlated random effects) device
#'
#' Augments a data frame with the unit-level means of the numeric case-specific
#' covariates in the formula (every part of a multipart formula; in the long
#' layout only variables that are constant within a case are averaged, over
#' the unit's cases), and returns the augmented formula and data. A covariate
#' constant within every unit is not averaged (it is its own unit mean; its
#' coefficient is the between-unit effect), with a message. Fitting
#' [remnl()] on the result implements the Mundlak / correlated-random-effects
#' specification: the coefficients on the original covariates are the
#' within-unit effects while those on the means capture (and test) the
#' correlation between the covariates and the unit effect, without the
#' incidental-parameters bias of unit dummies in short panels. A pooled
#' [mnl()] fit of the augmented formula with `cluster = "unit"` is the
#' population-averaged version: its slopes are attenuated towards zero by the
#' unit-effect variance (in the vignette's example by about a quarter), as in
#' any nonlinear model, and are consistent for average effects, not for the
#' within-unit coefficients.
#'
#' @param formula A model formula, possibly with `|` parts.
#' @param data A data frame.
#' @param unit Column name identifying the unit.
#' @param id Long layout: the case column (so that each case counts once in the
#'   unit mean); `NULL` for the wide layout.
#' @param weights Optional frequency weights (a column name or a vector aligned
#'   to the rows of `data`); the unit means are then weighted means, as the fit
#'   will be.
#' @param suffix Suffix of the added mean columns (default `"_mean"`).
#' @return A list with `formula` (the augmented formula, keeping the `|`
#'   structure), `data`, and `added` (the names of the mean columns).
#' @references Mundlak, Y. (1978). On the pooling of time series and cross
#'   section data. \emph{Econometrica}, 46, 69-85. \doi{10.2307/1913646}
#' @examples
#' ## 40 units of 15 cases; each unit's effect on alternative b is
#' ## correlated with its mean of x1
#' set.seed(1)
#' unit <- rep(1:40, each = 15); u <- rnorm(40)
#' x1 <- rnorm(600) + u[unit]; x2 <- rnorm(600)
#' V <- cbind(0, 0.5 - x1 + 0.8 * x2 + u[unit], -0.3 + 0.6 * x1 - 0.5 * x2)
#' y <- apply(exp(V) / rowSums(exp(V)), 1, function(p) sample(c("a", "b", "c"), 1, prob = p))
#' d <- data.frame(y = factor(y), x1, x2, unit)
#' m <- mundlak(y ~ x1 + x2, d, unit = "unit")
#' m$formula
#' fit <- mnl(m$formula, m$data, cluster = "unit")
#' ## the unit mean of x1 picks up the unit effect in b's equation; without
#' ## the means x1's slope absorbs it. This pooled fit estimates average
#' ## effects: its x1 slope is close to the within-unit -1 of the simulation
#' ## only because the unit mean takes up nearly all of the unit effect here
#' ## (15 cases per unit). remnl() on m$formula estimates the within-unit
#' ## coefficients.
#' coef(fit)[grep("x1", names(coef(fit)))]
#' coef(mnl(y ~ x1 + x2, d, cluster = "unit"))[c("b:x1", "c:x1")]
#' @seealso [mnl()], [imnl()]; the `fe` and `cluster` arguments of the estimators.
#' @family panel tools
#' @export
mundlak <- function(formula, data, unit, id = NULL, weights = NULL, suffix = "_mean") {
  if (!unit %in% names(data)) stop("'unit' must be a column in 'data'.", call. = FALSE)
  if (!is.null(id) && !id %in% names(data)) stop("'id' must be a column in 'data'.", call. = FALSE)
  wv <- if (is.null(weights)) rep(1, nrow(data)) else if (is.character(weights)) { if (!weights %in% names(data)) stop("'weights' must be a column in 'data'.", call. = FALSE); as.numeric(data[[weights]]) } else as.numeric(weights)
  if (length(wv) != nrow(data)) stop("'weights' must have one entry per row of 'data'.", call. = FALSE)
  F <- Formula::Formula(formula)
  np <- length(F)[2L]
  uf <- data[[unit]]
  aug <- data
  added <- character(0); flat <- character(0)
  ## case-level weights for the mean: each case once
  one_per_case <- if (is.null(id)) rep(TRUE, nrow(data)) else !duplicated(data[[id]])
  parts_rhs <- character(np)
  for (i in seq_len(np)) {
    tl <- attr(stats::terms(F, rhs = i), "term.labels")
    vars <- all.vars(stats::formula(F, lhs = 0L, rhs = i))
    for (v in intersect(vars, names(data))) {
      if (!is.numeric(data[[v]]) || v == unit || identical(v, id)) next
      if (!is.null(id)) {
        ## alternative-specific variables (varying within a case) are not averaged
        if (any(stats::ave(data[[v]], data[[id]], FUN = function(z) max(z) - min(z)) > 1e-10)) next
      }
      ## a covariate constant within every unit is its own unit mean: averaging
      ## it would duplicate the column (stress finding F16)
      vv0 <- data[[v]][one_per_case]; uu0 <- uf[one_per_case]; ok0 <- !is.na(vv0)
      if (all(tapply(vv0[ok0], uu0[ok0], function(z) diff(range(z))) <= 1e-10 * max(1, max(abs(vv0[ok0]))))) { flat <- c(flat, v); next }
      nm <- paste0(v, suffix)
      if (!nm %in% names(aug)) {
        ## the unit means are frequency-weighted when weights are given (as the fit will be)
        vv <- data[[v]][one_per_case]; ww <- wv[one_per_case]; ok <- !is.na(vv)
        mu <- tapply(vv[ok] * ww[ok], uf[one_per_case][ok], sum) / tapply(ww[ok], uf[one_per_case][ok], sum)
        aug[[nm]] <- as.numeric(mu[as.character(uf)])
      }
      added <- unique(c(added, nm))
      tl <- c(tl, nm)
    }
    parts_rhs[i] <- if (length(tl)) paste(tl, collapse = " + ") else "1"
    if (attr(stats::terms(F, rhs = i), "intercept") == 0L) parts_rhs[i] <- paste("0 +", parts_rhs[i])
  }
  lhs <- deparse(attr(F, "lhs")[[1L]], width.cutoff = 500L)
  f <- stats::as.formula(paste(lhs, "~", paste(parts_rhs, collapse = " | ")))
  environment(f) <- environment(formula)
  flat <- unique(flat)
  if (length(flat)) message("Not averaged (constant within every unit, so the variable is its own unit mean): ", paste(flat, collapse = ", "), ".")
  if (!length(added)) warning("No numeric case-level covariates found to average.", call. = FALSE)
  list(formula = f, data = aug, added = added)
}
