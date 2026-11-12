## ---------------------------------------------------------------------------
## Quantities of interest: relative risk ratios and first differences in
## predicted probabilities (delta-method or simulation intervals).
## ---------------------------------------------------------------------------

#' Relative risk ratios
#'
#' Exponentiated coefficients of a choice model: for an alternative-specific
#' coefficient \eqn{\beta_{jm}}, \eqn{e^{\beta_{jm}}} is the factor by which a
#' unit increase in covariate \eqn{m} multiplies the odds of alternative
#' \eqn{j} relative to the base (Stata's `mlogit, rrr`); for a generic
#' coefficient it multiplies the odds of any alternative relative to any other
#' per unit of the alternative-specific covariate; for the split equation it
#' multiplies the odds of the choice regime (logit link only).
#'
#' The odds reading holds for the multinomial logit, the inflated multinomial logit (whose
#' outcome stage is a multinomial logit within the choice regime), and the
#' random-intercept logit, where it is conditional on the unit's intercepts
#' (the population-averaged odds ratio is closer to one). It does not hold for
#' the nested logit (within a nest the ratio is \eqn{e^{\beta/\lambda}}, across
#' nests it is not constant), the HEV model (the odds depend on the scales), or
#' the probit, so `rrr()` refuses those fits: use [first_difference()].
#'
#' @param object A `"choicemodel"` fit.
#' @param level Confidence level for the (Wald, exponentiated) interval.
#' @return A data frame with `term`, `block`, `rrr`, `lower`, `upper`, and the
#'   z statistic and p-value of the underlying coefficient.
#' @examples
#' set.seed(1)
#' d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
#' rrr(mnl(y ~ x1 + x2, data = d))
#' @seealso [first_difference()], [summary.choicemodel()]
#' @family quantities of interest
#' @export
rrr <- function(object, level = 0.95) {
  if (!inherits(object, "choicemodel")) stop("'object' must be a fitted polychoice model.", call. = FALSE)
  if (object$family %in% c("nlogit", "hev", "mnp"))
    stop("Exponentiated coefficients are not odds ratios in a ", .pc_title(object), " (the odds depend on ",
         switch(object$family, nlogit = "the inclusive-value parameters", hev = "the scales", mnp = "the error covariance"),
         "); use first_difference() for effects on the probabilities.", call. = FALSE)
  est <- object$coefficients
  v <- if (is.null(object$vcov)) rep(NA_real_, length(est)) else diag(object$vcov)
  v[is.finite(v) & v <= 0] <- NA_real_
  se <- sqrt(pmax(v, 0))
  zq <- stats::qnorm(1 - (1 - level) / 2)
  keep <- object$blocks %in% c("outcome", "generic", "inflation")
  if (object$inflated && object$link != "logit") keep <- keep & object$blocks != "inflation"
  out <- data.frame(term = names(est)[keep], block = object$blocks[keep], rrr = exp(est[keep]),
                    lower = exp(est[keep] - zq * se[keep]), upper = exp(est[keep] + zq * se[keep]),
                    z = est[keep] / se[keep], p.value = 2 * stats::pnorm(-abs(est[keep] / se[keep])), row.names = NULL)
  class(out) <- c("rrr", "data.frame")
  out
}

#' @export
print.rrr <- function(x, digits = 3, ...) {
  y <- as.data.frame(x); y[c("rrr", "lower", "upper", "z")] <- lapply(y[c("rrr", "lower", "upper", "z")], round, digits)
  y$p.value <- format.pval(y$p.value, digits = digits)
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}

## a one-case profile for the wide layout: numeric covariates at their
## (weighted) means; factors, logicals, characters, the fixed-effect unit, and
## numeric variables used as factors in the formula (factor(year)) at their
## (weighted) modes
.pc_profile <- function(object) {
  D <- object$design
  d <- D$data_case
  w <- object$weights
  xl <- c(names(D$tm$x$xlevels), names(D$tm$w$xlevels), names(D$tm$z$xlevels))
  wrapped <- sub("^(as\\.)?factor\\((.*)\\)$", "\\2", xl[grepl("^(as\\.)?factor\\(", xl)])
  categorical <- c(D$fe, wrapped)
  prof <- lapply(names(d), function(nm) {
    v <- d[[nm]]
    if (is.numeric(v) && !nm %in% categorical) sum(w * v) / sum(w)
    else {
      tb <- tapply(w, v, sum); m <- names(tb)[which.max(tb)]
      if (is.factor(v)) factor(m, levels = levels(v)) else if (is.logical(v)) as.logical(m) else if (is.numeric(v)) as.numeric(m) else m
    }
  })
  names(prof) <- names(d)
  as.data.frame(prof, stringsAsFactors = FALSE)
}

## set a variable in a data frame to a value, keeping factor levels
.pc_set_var <- function(d, var, value) {
  v <- d[[var]]
  if (is.factor(v)) {
    if (!value %in% levels(v)) stop("'", value, "' is not a level of ", var, " (", paste(levels(v), collapse = ", "), ").", call. = FALSE)
    d[[var]] <- factor(rep(as.character(value), nrow(d)), levels = levels(v))
  } else if (is.logical(v)) d[[var]] <- rep(as.logical(value), nrow(d))
  else if (is.character(v)) d[[var]] <- rep(as.character(value), nrow(d))
  else d[[var]] <- rep(as.numeric(value), nrow(d))
  d
}

#' First differences in predicted probabilities
#'
#' The change in the predicted probability of every alternative when a covariate
#' moves from `from` to `to`, holding the other covariates at a profile (the
#' estimation cases' weighted means and modes; the default in the wide layout),
#' at a supplied one-case `newdata`, or averaged over the estimation cases
#' (`average = TRUE`, the average first difference / average marginal effect;
#' the default in the long layout). Intervals by the delta method (numerical
#' Jacobian of the quantity with respect to the parameters) or by simulation
#' from the asymptotic normal distribution of the estimates (the practice of
#' Bagozzi and Marchetti 2017; this is not a bootstrap). A fit with an offset
#' carries it along: the estimation offsets for `average = TRUE`, their
#' weighted mean at a profile, and `offset` must be given with `newdata`.
#' For a fixed-effects fit the profile holds the unit constants at the modal
#' unit's values, an arbitrary choice in a balanced panel: `average = TRUE`
#' (the effect averaged over the estimation cases and their units) is the
#' usual quantity there.
#'
#' For inflated models, `type = "prob"` moves the probabilities of the observed
#' alternatives (regime mixture included) and `type = "prob_outcome"` the
#' outcome-stage probabilities, that is, the effects among cases in the
#' choice regime, which is what Bagozzi and Marchetti (2017) report when
#' comparing BIMNL with MNL effects; `type = "regime"` reports the change in
#' \eqn{P(s = 1)} (with several inflated alternatives, the probabilities of the
#' choice regime and of captivity to each). `stage` restricts the covariate change to one equation when
#' the covariate enters both: with `stage = "outcome"` the covariate moves in
#' the outcome equations while the split equation keeps it at `from`, and
#' conversely for `"inflation"`. The two stage effects need not add up to the
#' total (`"both"`): the mixture is nonlinear, so the remainder is an
#' interaction of the two changes.
#'
#' In the long layout the covariate is changed on every row (a case-specific
#' covariate) or only on the rows of one alternative (`alternative =`, for an
#' alternative-specific covariate), and `average = TRUE` is the default since
#' there is no natural one-case profile.
#'
#' @param object A `"choicemodel"` fit.
#' @param var Name of the covariate to change.
#' @param from,to Its two values (a factor level for a factor).
#' @param newdata Optional data frame giving the profile (one case: one row in
#'   the wide layout, that case's rows in the long layout).
#' @param average If `TRUE`, average the first differences over the estimation
#'   cases (weighted by the frequency weights). `NULL` (default) means `FALSE`
#'   in the wide layout (a profile) and `TRUE` in the long layout, where there
#'   is no natural one-case profile.
#' @param alternative Long layout: the alternative whose rows receive the
#'   change (for an alternative-specific covariate); `NULL` changes all rows.
#' @param type `"prob"` (default), `"prob_outcome"`, or `"regime"` (inflated
#'   models).
#' @param stage `"both"` (default), `"outcome"`, or `"inflation"`: which
#'   equation(s) see the change (inflated models).
#' @param ci `"delta"` (default), `"sim"` (simulation from the asymptotic
#'   normal distribution of the estimates), or `"none"`.
#' @param level Confidence level.
#' @param R Number of simulation draws for `ci = "sim"`.
#' @param offset Offset for `newdata` when the fit has one (an n x J matrix in
#'   the wide layout, a vector aligned to the rows in the long layout).
#' @return A data frame of class `"pc_fd"` with one row per alternative (for
#'   `type = "regime"`, one row per regime: one with a single inflated
#'   alternative, the choice regime and each captive regime with several):
#'   `from` and `to` are the predicted
#'   probabilities at the two covariate values, `diff` their difference, and
#'   `lower`, `upper` the interval for `diff`.
#' @references Bagozzi, B.E. and Marchetti, K. (2017). Distinguishing
#'   occasional abstention from routine indifference in models of vote choice.
#'   \emph{Political Science Research and Methods}, 5(2), 277-294.
#'   \doi{10.1017/psrm.2015.42}
#' @examples
#' set.seed(1)
#' d <- rmnl(1000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
#' m <- mnl(y ~ x1 + x2, data = d)
#' first_difference(m, "x1", 0, 1)
#' first_difference(m, "x1", 0, 1, average = TRUE, ci = "sim", R = 500)
#' @seealso [predict.choicemodel()], [rrr()]
#' @family quantities of interest
#' @export
first_difference <- function(object, var, from, to, newdata = NULL, average = NULL, alternative = NULL,
                             type = c("prob", "prob_outcome", "regime"), stage = c("both", "outcome", "inflation"),
                             ci = c("delta", "sim", "none"), level = 0.95, R = 1000, offset = NULL) {
  if (!inherits(object, "choicemodel")) stop("'object' must be a fitted polychoice model.", call. = FALSE)
  type <- match.arg(type); stage <- match.arg(stage); ci <- match.arg(ci)
  D <- object$design
  if (!object$inflated && type != "prob") stop("type = \"", type, "\" applies to inflated models only.", call. = FALSE)
  if (!object$inflated && stage != "both") stop("'stage' applies to inflated models only.", call. = FALSE)
  if (is.null(average)) average <- D$layout == "long" && is.null(newdata)
  base_d <- if (!is.null(newdata)) {
    if (!is.data.frame(newdata)) stop("'newdata' must be a data frame.", call. = FALSE)
    if (D$layout == "wide" && nrow(newdata) != 1L) stop("'newdata' must be a one-row data frame (one case).", call. = FALSE)
    newdata
  } else if (average) D$data else {
    if (D$layout == "long") stop("In the long layout a profile needs 'newdata' (the rows of one case, one per alternative) or average = TRUE.", call. = FALSE)
    .pc_profile(object)
  }
  if (!var %in% names(base_d)) stop("'", var, "' is not a variable of the model data.", call. = FALSE)
  ## the offset of the fit travels with the estimation data (average = TRUE) and
  ## is the weighted mean offset at a profile; for user data it must be supplied
  off_use <- if (!is.null(newdata)) offset
             else if (average) { if (D$layout == "wide") D$off else D$off[cbind(D$ci, D$ji)] }
             else if (isTRUE(D$has_offset)) matrix(colSums(D$off * object$weights) / sum(object$weights), 1L) else NULL
  if (!is.null(alternative)) {
    if (D$layout != "long") stop("'alternative' applies to the long layout.", call. = FALSE)
    if (!alternative %in% D$levels) stop("'alternative' must be one of ", paste(D$levels, collapse = ", "), ".", call. = FALSE)
  }
  set_rows <- function(d, value) {
    if (is.null(alternative)) return(.pc_set_var(d, var, value))
    rows <- as.character(d[[D$alt]]) == alternative
    d2 <- .pc_set_var(d[rows, , drop = FALSE], var, value)
    d[rows, ] <- d2
    d
  }
  d_from <- set_rows(base_d, from); d_to <- set_rows(base_d, to)
  w <- if (average) { if (D$layout == "wide") object$weights else object$weights } else NULL
  ## the two designs are built once (the Jacobian evaluates the quantity many
  ## times); the stage switch takes the outcome design from one frame and the
  ## inflation design from the other
  N_from <- .pc_design_new(object, d_from, offset = off_use); N_to <- .pc_design_new(object, d_to, offset = off_use)
  fns <- .pc_fns(object$design)
  rs <- if (identical(object$family, "remnl")) max(exp(object$theta[D$idx$str])) else NULL     # the prediction rule, held at the estimate
  qfun <- function(th) {
    if (stage == "both") return(.pc_q_N(fns, th, N_to, N_to, type, w, rs) - .pc_q_N(fns, th, N_from, N_from, type, w, rs))
    .pc_q_N(fns, th, if (stage == "inflation") N_from else N_to, if (stage == "outcome") N_from else N_to, type, w, rs) -
      .pc_q_N(fns, th, N_from, N_from, type, w, rs)
  }
  ## the quantities are functions of the internal parameter vector (the
  ## structural transforms are inside fn$probs), so the internal covariance is used
  pf <- .pc_q_N(fns, object$theta, N_from, N_from, type, w, rs)
  pt <- pf + qfun(object$theta)
  comp <- if (type == "regime") { if (length(D$kc) > 1L) c("choice regime", paste0("captive: ", D$levels[D$kc])) else "P(choice regime)" } else D$levels
  lower <- upper <- rep(NA_real_, length(comp)); method <- "none"
  if (ci != "none") {
    if (is.null(object$vcov_internal)) warning("No covariance available (se = \"none\"); intervals are NA.", call. = FALSE)
    else {
      V <- object$vcov_internal; zq <- stats::qnorm(1 - (1 - level) / 2)
      na_msg <- paste0("The covariance has no entries for some parameters the quantity depends on (with unit fixed effects ",
                       "and se = \"bootstrap\" the constants and unit effects have none); use se = \"analytic\", \"robust\", or \"cluster\" for these intervals.")
      if (ci == "delta") {
        G <- numDeriv::jacobian(qfun, object$theta)
        if (any(colSums(abs(G)) > 0 & apply(is.na(V), 2L, any))) stop(na_msg, call. = FALSE)
        V[is.na(V)] <- 0
        qf <- rowSums((G %*% V) * G)
        qa <- rowSums((abs(G) %*% abs(V)) * abs(G))                  # the scale of the quadratic form's terms
        se <- ifelse(qf < -sqrt(.Machine$double.eps) * qa, NA_real_, sqrt(pmax(qf, 0)))   # an indefinite covariance: no interval
        lower <- (pt - pf) - zq * se; upper <- (pt - pf) + zq * se; method <- "delta"
      } else {
        if (anyNA(V)) stop(na_msg, call. = FALSE)
        draws <- MASS::mvrnorm(R, object$theta, V)
        sims <- apply(draws, 1L, qfun)
        if (is.null(dim(sims))) sims <- matrix(sims, nrow = 1L)
        a <- (1 - level) / 2
        lower <- apply(sims, 1L, stats::quantile, a); upper <- apply(sims, 1L, stats::quantile, 1 - a)
        method <- paste0("simulation (R=", R, ")")
      }
    }
  }
  out <- data.frame(alternative = comp, from = pf, to = pt, diff = pt - pf, lower = lower, upper = upper, row.names = NULL)
  attr(out, "var") <- var; attr(out, "values") <- c(from, to); attr(out, "method") <- method; attr(out, "level") <- level
  attr(out, "type") <- type; attr(out, "average") <- average; attr(out, "stage") <- stage
  class(out) <- c("pc_fd", "data.frame")
  out
}

## quantity from prebuilt designs: the outcome part (X, W, offsets,
## availability) from No and the inflation part (Z) from Ni; averaged
## probabilities or the regime probability
.pc_q_N <- function(fns, theta, No, Ni, type, w, rs = NULL) {
  pr <- if (is.null(rs)) fns$probs(theta, X = No$X, W = No$W, Z = Ni$Z, off = No$off, avail = No$avail, fe_off = No$fe_off)
        else fns$probs(theta, X = No$X, W = No$W, Z = Ni$Z, off = No$off, avail = No$avail, fe_off = No$fe_off, rule_sigma = rs)
  if (is.null(w)) w <- rep(1, No$n)
  if (type == "regime") return(if (is.null(pr$regime)) sum(w * pr$phi) / sum(w) else colSums(pr$regime * w) / sum(w))
  P <- if (type == "prob") pr$prob else pr$prob_outcome
  colSums(P * w) / sum(w)
}


#' @export
print.pc_fd <- function(x, digits = 3, ...) {
  v <- attr(x, "values")
  cat(sprintf("First differences: %s from %s to %s%s (%s%s)\n", attr(x, "var"), format(v[1L]), format(v[2L]),
              if (isTRUE(attr(x, "average"))) ", averaged over the estimation cases" else " at the covariate profile",
              if (attr(x, "method") == "none") "no intervals" else sprintf("%s intervals at %g%%", attr(x, "method"), 100 * attr(x, "level")),
              if (attr(x, "stage") != "both") paste0("; ", attr(x, "stage"), " stage only") else ""))
  if (attr(x, "type") != "prob") cat("Quantity:", if (attr(x, "type") == "prob_outcome") "outcome-stage probabilities (choice regime)" else "P(choice regime)", "\n")
  y <- as.data.frame(x); y[-1L] <- lapply(y[-1L], round, digits)
  print.data.frame(y, row.names = FALSE)
  invisible(x)
}
