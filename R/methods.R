## ---------------------------------------------------------------------------
## S3 methods for "choicemodel" fits (all six estimators): print, summary, coef, vcov,
## logLik, nobs, fitted, residuals, predict, confint, update, formula, simulate.
## ---------------------------------------------------------------------------

#' @export
print.choicemodel <- function(x, digits = 4, ...) {
  cat(.pc_title(x), "\n")
  if (!is.null(x$call)) { cat("Call: "); print(x$call) }
  cat("Alternatives: ", .pc_level_list(x$levels), "  (base: ", x$ref, ")\n", sep = "")
  cat("\nCoefficients:\n")
  print(round(x$coefficients, digits))
  cat(sprintf("\nLog-likelihood: %.*f on %d parameters; %d cases%s\n", digits - 1L, x$loglik, x$df, x$n,
              if (x$nobs_weighted != x$n) sprintf(" (weight total %g)", x$nobs_weighted) else ""))
  if (.pc_penalized(x)) .pc_firth_note(x)
  if (!isTRUE(x$converged)) cat("Warning: the optimizer did not converge.\n")
  if (isTRUE(x$boundary)) cat("Warning: the fit is at a boundary (see ?", if (x$inflated) "imnl" else "mnl", ").\n", sep = "")
  if (isTRUE(x$bias_correct == "jackknife")) cat("Slopes are split-panel jackknife bias-corrected; constants, unit effects, log-likelihood, and standard errors are maximum-likelihood (see ?mnl).\n")
  invisible(x)
}

## the lines a penalized fit prints (print and summary): Firth's penalty or a
## prior (the posterior mode)
.pc_firth_note <- function(x) {
  if (isTRUE(x$firth))
    cat(sprintf("Penalized (Firth) estimates: penalized log-likelihood %.3f; the log-likelihood above is the ordinary one at these estimates.\n", x$loglik_penalized))
  else
    cat(sprintf("Posterior-mode estimates under a %s; log-likelihood plus log prior %.3f; the log-likelihood above is the ordinary one at these estimates.\n",
                .pc_penalty_label(x), x$loglik_penalized))
  if (length(x$separation))
    cat("Separation: the maximum-likelihood estimates of ", paste(utils::head(x$separation, 6L), collapse = ", "), if (length(x$separation) > 6L) ", ..." else "",
        " diverge; their penalized estimates are determined largely by the ", if (isTRUE(x$firth)) "penalty" else "prior", ". Use profile intervals (confint())",
        if (isTRUE(x$firth)) "" else " and compare other priors", ".\n", sep = "")
  invisible(x)
}

#' @export
coef.choicemodel <- function(object, ...) object$coefficients

#' @export
vcov.choicemodel <- function(object, ...) object$vcov

#' @export
logLik.choicemodel <- function(object, ...) {
  ll <- structure(object$loglik, df = object$df, nobs = object$nobs_weighted, class = "logLik")
  if (.pc_penalized(object)) attr(ll, "penalized") <- TRUE         # the ordinary log-likelihood at penalized estimates
  ll
}

#' @export
nobs.choicemodel <- function(object, ...) object$n

#' @export
fitted.choicemodel <- function(object, ...) object$fitted.values

#' @export
formula.choicemodel <- function(x, ...) x$formula

#' Residuals of a choice model
#'
#' @param object A `"choicemodel"` fit.
#' @param type `"response"` (the \eqn{n \times J} matrix of indicator minus
#'   fitted probability) or `"pearson"` (response residuals divided by
#'   \eqn{\sqrt{P(1-P)}}).
#' @param ... Unused.
#' @return An \eqn{n \times J} matrix.
#' @export
residuals.choicemodel <- function(object, type = c("response", "pearson"), ...) {
  type <- match.arg(type)
  P <- object$fitted.values
  Y <- matrix(0, nrow(P), ncol(P), dimnames = dimnames(P)); Y[cbind(seq_len(nrow(P)), object$y)] <- 1
  R <- Y - P
  if (type == "pearson") R <- R / sqrt(pmax(P * (1 - P), .Machine$double.eps))
  R
}

#' Summarize a choice model
#'
#' @param object A `"choicemodel"` fit.
#' @param ... Unused.
#' @return An object of class `"summary.choicemodel"` with the coefficient
#'   table (`coefficients`: estimate, standard error, z, p), the parameter
#'   `block` of every row (`"outcome"`, `"generic"`, `"fixed effect"`,
#'   `"inflation"`, `"nest"`, `"scale"`, `"covariance"`, `"random intercept"`),
#'   and the fit statistics. The z statistics of the inclusive-value and scale
#'   blocks test the value 1 (the multinomial logit); those of the multinomial probit
#'   covariance block test the values under independent, identically
#'   distributed errors (variances 1, covariances 0.5 under the normalization
#'   of the first variance to 1); the random-intercept
#'   standard deviations get no z test (their null is on the boundary); a
#'   parameter whose variance is negative (the Hessian is not positive
#'   definite) has standard error `NA`. The printed table of a penalized fit
#'   (`firth = TRUE` or a prior) marks `[sep]` the coefficients that move along
#'   a separating direction of the data ([separation()]). Their Wald tests rest
#'   on a quadratic approximation of the penalized log-likelihood that fails
#'   along that direction and can be far off (Heinze and Schemper 2002); under
#'   a prior they
#'   can have little power (Rainey 2024). Test them by `term_test(type = "lr")` or
#'   [lr_test()]. A
#'   maximum-likelihood fit names its diverging coefficients in a warning at
#'   fit time instead.
#' @seealso [fit_stats()] for the pseudo-R-squared measures.
#' @export
summary.choicemodel <- function(object, ...) {
  est <- object$coefficients
  v <- if (is.null(object$vcov)) rep(NA_real_, length(est)) else diag(object$vcov)
  v[is.finite(v) & v <= 0] <- NA_real_                                           # a negative or exactly zero variance (an indefinite or generalized inverse): no standard error
  se <- sqrt(pmax(v, 0))
  ## the hypothesis of each z test: 1 for inclusive values and scales (the
  ## multinomial logit); the iid values for the multinomial probit covariance ((I + 11') / 2: variances
  ## 1, covariances 0.5); 0 otherwise
  null <- ifelse(object$blocks %in% c("nest", "scale"), 1, 0)
  cv <- object$blocks == "covariance"
  null[cv] <- ifelse(grepl("^var:", names(est)[cv]), 1, 0.5)
  z <- (est - null) / se
  z[object$blocks == "random intercept"] <- NA_real_                               # sd = 0 is on the boundary: use lr_test()
  ct <- cbind(Estimate = est, `Std. Error` = se, `z value` = z, `Pr(>|z|)` = 2 * stats::pnorm(-abs(z)))
  rownames(ct) <- names(est)
  structure(list(call = object$call, title = .pc_title(object), coefficients = ct, block = object$blocks,
                 levels = object$levels, ref = object$ref, inflate = object$inflate, se = object$se,
                 loglik = object$loglik, df = object$df, n = object$n, nobs_weighted = object$nobs_weighted,
                 aic = -2 * object$loglik + 2 * object$df, bic = -2 * object$loglik + log(object$nobs_weighted) * object$df,
                 share_inflated = if (object$inflated) vapply(.pc_kc(object), function(c) stats::weighted.mean(object$y == c, object$weights), numeric(1)) else NULL,
                 mean_phi = if (object$inflated) { ka <- rowSums(object$design$avail[, .pc_kc(object), drop = FALSE]) > 0L
                                                  stats::weighted.mean(object$phi[ka], object$weights[ka]) } else NULL,
                 converged = object$converged, boundary = object$boundary, ill_conditioned = object$ill_conditioned, singular = object$singular,
                 n_modes = object$n_modes, n_starts = if (is.null(object$screen_logliks)) length(object$start_logliks) else length(object$screen_logliks),
                 reported_type = if (is.null(object$search)) NULL else object$search$type[object$search$reported],
                 fe_info = object$fe_info, bias_correct = object$bias_correct, jackknife_refused = object$jackknife_refused,
                 firth = object$firth, prior = object$prior, loglik_penalized = object$loglik_penalized, separation = object$separation), class = "summary.choicemodel")
}

#' @export
print.summary.choicemodel <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  cat(x$title, "\n")
  if (!is.null(x$call)) { cat("Call: "); print(x$call) }
  cat("Alternatives: ", .pc_level_list(x$levels), "  (base: ", x$ref, ")\n", sep = "")
  ct <- x$coefficients
  ## a penalized fit: the separated coefficients' rows are marked (their Wald
  ## tests are unreliable; DESIGN.md Phase 4d)
  sepr <- if (.pc_penalized(x) && length(x$separation)) rownames(ct) %in% x$separation else rep(FALSE, nrow(ct))
  if (any(sepr)) rownames(ct)[sepr] <- paste(rownames(ct)[sepr], "[sep]")
  show_block <- function(lab, title) {
    rows <- x$block == lab
    if (!any(rows)) return(invisible())
    cat("\n", title, ":\n", sep = "")
    stats::printCoefmat(ct[rows, , drop = FALSE], digits = digits, signif.stars = TRUE, na.print = "NA")
  }
  show_block("outcome", "Outcome equations (alternative:term)")
  show_block("generic", "Generic coefficients (alternative-specific covariates)")
  if (any(x$block == "fixed effect")) {
    nfe <- sum(x$block == "fixed effect")
    cat("\n", nfe, " unit fixed-effect constants not shown (summary(x)$coefficients holds them).\n", sep = "")
  }
  show_block("inflation", if (length(x$inflate) > 1L) "Split equations (log P(choice regime) / P(captive to the alternative) = z'delta)" else "Split equation (P(choice regime) = F(z'delta))")
  show_block("nest", "Inclusive-value parameters lambda (z tests are against 1 = multinomial logit; the joint test is lr_test() against the mnl() fit)")
  show_block("scale", "Gumbel scale parameters (base alternative = 1; z tests are against 1; the joint test is lr_test() against the mnl() fit)")
  show_block("covariance", "Covariance of the utility differences from the base (first variance fixed at 1; z tests are against the iid values, variances 1 and covariances 0.5)")
  show_block("random intercept", "Random alternative-specific intercepts: standard deviations (no z test: sd = 0 is on the boundary; use lr_test() against the mnl() fit)")
  cat(sprintf("\nStandard errors: %s. Log-likelihood %.3f on %d parameters, %d cases%s. AIC %.2f, BIC %.2f.\n",
              x$se, x$loglik, x$df, x$n, if (x$nobs_weighted != x$n) sprintf(" (weight total %g)", x$nobs_weighted) else "", x$aic, x$bic))
  if (!is.null(x$share_inflated) && length(x$inflate) > 1L)
    cat(sprintf("Inflated alternatives %s: observed shares %s; mean fitted P(choice regime) %.3f.\n", paste0("'", x$inflate, "'", collapse = ", "),
                paste(sprintf("%.3f", x$share_inflated), collapse = ", "), x$mean_phi))
  else if (!is.null(x$share_inflated))
    cat(sprintf("Inflated alternative '%s': observed share %.3f; mean fitted P(choice regime) %.3f.\n", x$inflate, x$share_inflated, x$mean_phi))
  ## the search's maxima (object$search; a penalized fit's two runs in
  ## object$start_logliks)
  if (!is.null(x$n_modes) && !is.na(x$n_modes)) {
    if (!is.null(x$reported_type) && x$reported_type != "maximum" && x$n_modes > 0L)
      cat(sprintf("Note: the reported fit, the highest end of the search, is not a local maximum; the search from %d starts also reached %d below it (see object$search).\n",
                  x$n_starts, x$n_modes))
    else if (x$n_modes > 1L)
      cat(sprintf("Note: %d distinct local maxima among the ends of the search from %d starts; the highest is reported (see object$%s).\n",
                  x$n_modes, x$n_starts, if (is.null(x$reported_type)) "start_logliks" else "search"))
  }
  if (.pc_penalized(x)) {
    .pc_firth_note(x)
    cat("The z tests are Wald tests; profile penalized-likelihood intervals (confint()) and penalized likelihood-ratio tests (term_test(type = \"lr\"), lr_test()) are recommended. AIC and BIC above are not maximum-likelihood criteria.\n")
    if (any(sepr))
      cat("[sep]: a separated coefficient. Its Wald z and p-value rest on a quadratic approximation that fails along the separating direction",
          if (isTRUE(x$firth)) " and can be far off (Heinze and Schemper 2002)"
          else "; Wald tests under a prior can have little power (Rainey 2024)",
          ". Test it by term_test(type = \"lr\") or lr_test().\n", sep = "")
  }
  if (!isTRUE(x$converged)) cat("Warning: the optimizer did not converge.\n")
  if (isTRUE(x$boundary)) cat("Warning: the fit is at a boundary.\n")
  ## a singular information matrix: the table's rows for the parameters the
  ## data do not identify look like any others, so the footer says so
  if (isTRUE(x$singular)) cat("Warning: the information matrix is singular (a generalized inverse was used): some parameters are not identified by the data, and their estimates and standard errors are not meaningful.\n")
  else if (isTRUE(x$ill_conditioned)) cat("Warning: the information matrix is ill-conditioned.\n")
  if (isTRUE(x$bias_correct == "jackknife")) cat("Slopes are split-panel jackknife bias-corrected; constants, unit effects, log-likelihood, and standard errors are maximum-likelihood (see ?mnl).\n")
  if (!is.null(x$jackknife_refused)) cat("The split-panel jackknife was REFUSED (", x$jackknife_refused, "); the estimates are maximum-likelihood.\n", sep = "")
  invisible(x)
}

#' Confidence intervals for a choice model
#'
#' Wald intervals from the fit's covariance, percentile intervals from the
#' bootstrap replicates when the model was fit with `se = "bootstrap"`, or,
#' for a penalized fit (`firth = TRUE` or a `prior` in [mnl()] and
#' [rologit()]), profile penalized-likelihood intervals. The
#' Wald intervals of the positive structural parameters (inclusive values,
#' scales, standard deviations, variances) are computed on the log scale and
#' transformed back, so they stay positive. With unit fixed effects and the bootstrap, the
#' alternative-specific constants and the unit effects have no replicates and
#' their intervals are `NA`; percentile intervals are refused for a
#' jackknife-corrected fit (the replicates are not corrected).
#'
#' Profile intervals (Heinze and Schemper 2002) hold the values \eqn{c} of a
#' coefficient for which twice the drop of the penalized log-likelihood,
#' maximized over the other parameters with the coefficient held at \eqn{c},
#' is at most the chi-square quantile on one degree of freedom (the full
#' model's penalty throughout). Under separation they are asymmetric and wider
#' than the Wald intervals on the side towards which the maximum-likelihood
#' estimate runs off. A limit not reached within 50 standard errors is `NA`,
#' with a warning. With few cases per coefficient, and under separation, the
#' penalized likelihood with one coefficient held fixed can have several local
#' maxima. Each limit is therefore checked by maximizing afresh at it, from the
#' estimate and from zero; a limit at which a higher maximum turns up (after a
#' second search that takes the highest of the three starts at every point),
#' or at which the maximization does not converge, is `NA` with a warning,
#' since the number found would not be a limit of the interval.
#'
#' @param object A `"choicemodel"` fit.
#' @param parm Parameters (names or positions); default all.
#' @param level Confidence level.
#' @param type `"wald"`, `"percentile"` (bootstrap fits only), or `"profile"`
#'   (penalized fits only). Default: `"profile"` for a penalized fit, `"wald"`
#'   otherwise.
#' @param ... Unused.
#' @return A two-column matrix.
#' @references
#' Heinze, G. and Schemper, M. (2002). A solution to the problem of separation
#' in logistic regression. \emph{Statistics in Medicine}, 21, 2409-2419.
#' \doi{10.1002/sim.1047}
#' @export
confint.choicemodel <- function(object, parm, level = 0.95, type = c("wald", "percentile", "profile"), ...) {
  type <- if (missing(type) && .pc_penalized(object)) "profile" else match.arg(type)
  est <- object$coefficients
  pn <- names(est)
  if (missing(parm)) parm <- pn else if (is.numeric(parm)) parm <- pn[parm]
  a <- (1 - level) / 2
  if (type == "profile") {
    if (!.pc_penalized(object)) stop("Profile intervals are offered for penalized fits (firth = TRUE, or a prior); use type = \"wald\".", call. = FALSE)
    ci <- .pc_firth_profile(object, parm, level)
  } else if (type == "percentile") {
    if (is.null(object$boot)) stop("Percentile intervals need a fit with se = \"bootstrap\".", call. = FALSE)
    if (identical(object$bias_correct, "jackknife") && is.null(object$jackknife_refused))
      stop("The bootstrap replicates are those of the maximum-likelihood estimator, not of the jackknife-corrected one; use type = \"wald\".", call. = FALSE)
    B <- object$boot$replicates[, parm, drop = FALSE]                  # natural scale
    ci <- t(apply(B, 2L, function(b) if (all(is.na(b))) c(NA_real_, NA_real_) else stats::quantile(b, probs = c(a, 1 - a), na.rm = TRUE)))
  } else {
    if (is.null(object$vcov)) stop("No covariance matrix (se = \"none\").", call. = FALSE)
    v <- diag(object$vcov); v[is.finite(v) & v <= 0] <- NA_real_
    se <- sqrt(pmax(v, 0))[parm]
    zq <- stats::qnorm(1 - a)
    ## positive structural parameters (inclusive values, scales, standard
    ## deviations, variances) get the interval on the log scale, transformed back:
    ## a symmetric interval on the natural scale can cross zero
    pos <- object$blocks %in% c("nest", "scale", "random intercept") | grepl("^var:", pn)
    pos <- pos[match(parm, pn)] & is.finite(est[parm]) & est[parm] > 0
    lo <- est[parm] - zq * se; hi <- est[parm] + zq * se
    lo[pos] <- est[parm][pos] * exp(-zq * se[pos] / est[parm][pos]); hi[pos] <- est[parm][pos] * exp(zq * se[pos] / est[parm][pos])
    ci <- cbind(lo, hi)
  }
  dimnames(ci) <- list(parm, paste(format(100 * c(a, 1 - a), trim = TRUE, digits = 3), "%"))
  ci
}

#' Predictions from a choice model
#'
#' @param object A `"choicemodel"` fit.
#' @param newdata Optional data frame (same layout as the fit: for the long
#'   layout it must contain the `id` and `alt` columns, one row per available
#'   alternative; the choice indicator is not needed). It must hold every
#'   covariate of every equation (for an inflated fit the split-equation
#'   variables too) and, for a fixed-effects fit, the unit column (only units
#'   seen in the fit). Factor levels must be among those of the fit. With
#'   `NULL`, the estimation cases.
#' @param type `"prob"` (default): the \eqn{n \times J} matrix of probabilities
#'   of the observed alternatives, \eqn{P^*} (for a random-intercept fit the
#'   population-averaged probabilities); `"class"`: the modal alternative (a
#'   factor); `"link"`: the \eqn{n \times J} matrix of utilities \eqn{V_{ij}}
#'   (\eqn{-\infty} where unavailable); for nested-logit fits
#'   `"prob_conditional"`: the probabilities of the alternatives within their
#'   nests, \eqn{P(j \mid \text{nest})}; for inflated models also
#'   `"prob_outcome"`: the outcome-stage probabilities \eqn{P_{ij}} (the
#'   choice probabilities of a case in the choice regime, which is what Bagozzi and
#'   Marchetti (2017) hold the split equation at when comparing BIMNL and MNL
#'   effects), `"regime"`: \eqn{P(s_i = 1 \mid z_i)}, the probability of
#'   not being an inflated (always-\eqn{k}) case (with several inflated
#'   alternatives the matrix of the regime probabilities: the choice regime,
#'   then captivity to each inflated alternative), and `"posterior"` (estimation
#'   cases only): the probability of the choice regime given the observed
#'   choice, \eqn{\phi_i P_{ik} / P^*_{ik}} for a case that chose \eqn{k} and 1
#'   otherwise, which separates routine from occasional choosers of \eqn{k}
#'   (with several inflated alternatives the matrix of the regime posteriors);
#'   for rank-ordered logit fits also `"ranking"` (estimation cases only):
#'   the probability of each observed ranking. For a rank-ordered fit `"prob"`
#'   and `"class"` describe first choices.
#' @param offset Offset for `newdata` when the fit has one (an \eqn{n \times J}
#'   matrix in the wide layout, a vector aligned to the rows in the long
#'   layout; a long-layout fit whose offset was a column name finds it in
#'   `newdata`). A fit with an offset refuses `newdata` without one.
#' @param ... Unused.
#' @return A matrix (rows named by case), a factor, or a numeric vector.
#' @export
predict.choicemodel <- function(object, newdata = NULL, type = c("prob", "class", "link", "prob_outcome", "regime", "prob_conditional", "posterior", "ranking"),
                                offset = NULL, ...) {
  type <- match.arg(type)
  if (type == "ranking") {
    if (!identical(object$family, "rologit")) stop("type = \"ranking\" applies to rank-ordered logit fits only.", call. = FALSE)
    if (!is.null(newdata)) stop("type = \"ranking\" is the probability of the OBSERVED ranking; it is available for the estimation cases only (newdata = NULL).", call. = FALSE)
    return(stats::setNames(exp(object$lli), object$case_id))
  }
  if (type %in% c("prob_outcome", "regime", "posterior") && !object$inflated) stop("type = \"", type, "\" applies to inflated models only.", call. = FALSE)
  if (type == "prob_conditional" && object$family != "nlogit") stop("type = \"prob_conditional\" applies to nested logit fits only.", call. = FALSE)
  if (type == "posterior" && !is.null(newdata)) stop("type = \"posterior\" is the regime probability given the OBSERVED choice; it is available for the estimation cases only (newdata = NULL).", call. = FALSE)
  D <- object$design
  fn <- .pc_fns(D)
  if (is.null(newdata)) {
    pr <- fn$probs(object$theta); ids <- object$case_id
    if (type == "posterior" && length(D$kc) > 1L) {
      ## several captive alternatives: P(choice regime | y) = pi_0 P_y / P*_y and
      ## P(captive to c | y) = 1{y = c} pi_c / P*_y
      Py <- pr$prob_outcome[cbind(seq_len(D$n), D$y)]; Psy <- pr$prob[cbind(seq_len(D$n), D$y)]
      post <- matrix(0, D$n, 1L + length(D$kc), dimnames = list(ids, c("choice", D$levels[D$kc])))
      post[, 1L] <- pr$phi * Py / Psy
      for (t in seq_along(D$kc)) { it <- D$y == D$kc[t]; post[it, 1L + t] <- pr$regime[it, 1L + t] / Psy[it] }
      return(post)
    }
    if (type == "posterior") {
      ## P(choice regime | y): 1 for a case that chose another alternative, phi P_k / P*_k for a case that chose the inflated one
      post <- rep(1, D$n); ik <- D$y == D$k
      post[ik] <- pr$phi[ik] * pr$prob_outcome[ik, D$k] / pr$prob[ik, D$k]
      return(stats::setNames(post, ids))
    }
  } else {
    N <- .pc_design_new(object, newdata, offset)
    pr <- fn$probs(object$theta, X = N$X, W = N$W, Z = N$Z, off = N$off, avail = N$avail, fe_off = N$fe_off); ids <- N$case_id
    nna <- sum(!stats::complete.cases(pr$prob))
    if (nna) warning(nna, " case(s) of 'newdata' have missing covariate values; their predictions are NA.", call. = FALSE)
  }
  nm <- function(M) { dimnames(M) <- list(ids, D$levels); M }
  switch(type,
         prob = nm(pr$prob),
         prob_outcome = nm(pr$prob_outcome),
         prob_conditional = nm(pr$prob_conditional),
         link = nm(pr$V),
         regime = if (length(D$kc) > 1L) { Rg <- pr$regime; dimnames(Rg) <- list(ids, c("choice", D$levels[D$kc])); Rg } else stats::setNames(pr$phi, ids),
         class = factor(D$levels[vapply(seq_len(nrow(pr$prob)), function(i) if (anyNA(pr$prob[i, ])) NA_integer_ else which.max(pr$prob[i, ]), integer(1))],
                        levels = D$levels))
}

#' Update and refit a choice model
#'
#' Works on the stored call as [stats::update()] does, with multi-part formulas
#' handled by [Formula::update.Formula()]: `update(fit, . ~ . + x3)` changes the
#' first part, `update(fit, . ~ . | . + z2)` the second.
#'
#' @param object A `"choicemodel"` fit.
#' @param formula. Changes to the formula (optional).
#' @param ... Other arguments of the original call to replace.
#' @param evaluate If `FALSE`, return the updated call instead of refitting.
#' @return The refitted model (or the call).
#' @export
update.choicemodel <- function(object, formula., ..., evaluate = TRUE) {
  call <- object$call
  if (is.null(call)) stop("The fit has no stored call.", call. = FALSE)
  extras <- match.call(expand.dots = FALSE)$...
  if (!missing(formula.)) call$formula <- stats::formula(stats::update(object$formula, formula.))
  if (length(extras)) {
    existing <- !is.na(match(names(extras), names(call)))
    for (a in names(extras)[existing]) call[[a]] <- extras[[a]]
    if (any(!existing)) { call <- c(as.list(call), extras[!existing]); call <- as.call(call) }
  }
  if (evaluate) eval(call, parent.frame()) else call
}

#' Simulate outcomes from a fitted choice model
#'
#' Draws `nsim` outcome vectors from the fitted model (holding the covariates
#' fixed): the parametric bootstrap of the model. For every family but the
#' random-intercept one the draw of a case is from its fitted probabilities of
#' the observed alternatives (the GHK probabilities of a multinomial probit case are
#' renormalized to one); for [remnl()] one intercept vector is drawn per unit
#' and the unit's cases are then drawn given it, so the within-unit dependence
#' of the model is reproduced. Frequency weights are left as they are: one
#' draw per row, so for a weighted fit the draws are not the parametric
#' bootstrap of the weighted data (a row of weight w would need w draws;
#' [boot_pvalue()] expands the rows for that reason). For nominal outcomes there is no meaningful residual ordering,
#' so these draws serve test calibration and predictive checks rather than
#' DHARMa-style quantile residuals.
#'
#' @param object A `"choicemodel"` fit.
#' @param nsim Number of simulated outcome vectors.
#' @param seed Optional seed (as in [stats::simulate()]; the caller's
#'   random-number stream is restored afterwards).
#' @param ... Unused.
#' @param ties For a [rologit()] fit: `"untied"` (default), every case ranking
#'   as many positions as it ranked in the data, without ties; or `"observed"`,
#'   each case's observed pattern of rank groups (the same group sizes in the
#'   same order, a tie staying a tie at its positions), as [boot_lr_test()]
#'   draws rankings.
#' @return A data frame with `nsim` factor columns (`sim_1`, ...), one row per
#'   estimation case. For a [rologit()] fit, a list of `nsim` integer matrices
#'   (cases by alternatives) of simulated ranks (0 for the unranked and
#'   unavailable alternatives).
#' @export
simulate.choicemodel <- function(object, nsim = 1, seed = NULL, ..., ties = c("untied", "observed")) {
  ties_given <- !missing(ties); ties <- match.arg(ties)
  if (ties_given && !identical(object$family, "rologit")) message("'ties' applies to rank-ordered fits only and is ignored here.")
  if (!is.null(seed)) {
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(1)
    old <- get(".Random.seed", envir = globalenv()); on.exit(assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  D <- object$design; n <- D$n; J <- D$J
  if (identical(object$family, "rologit")) {
    ## Gumbel utilities rank the available alternatives; each case keeps the
    ## number of positions it ranked in the data
    V <- pc_ro_probs_cpp(object$theta, D$X, D$W, D$off, D$avail, J, D$base - 1L)$V
    L <- rowSums(D$str$rank > 0L & D$avail == 1L)
    grp <- if (ties == "observed") lapply(seq_len(n), function(i) { r <- D$str$rank[i, D$avail[i, ] == 1L]; tabulate(r[r > 0L]) }) else NULL
    out <- lapply(seq_len(nsim), function(s) {
      U <- V + matrix(-log(-log(stats::runif(n * J))), n, J)
      if (ties == "observed") { R <- .pc_ro_coarsen(U, D$avail, grp); dimnames(R) <- list(object$case_id, object$levels); return(R) }
      U[D$avail == 0L] <- -Inf
      pos <- t(apply(-U, 1L, rank, ties.method = "first"))
      R <- matrix(as.integer(pos), n, J); R[pos > L | D$avail == 0L] <- 0L
      dimnames(R) <- list(object$case_id, object$levels)
      R
    })
    names(out) <- paste0("sim_", seq_len(nsim))
    return(out)
  }
  draw <- function(P) {
    P <- P / rowSums(P)                                        # GHK rows need not sum to one exactly
    u <- stats::runif(n); cp <- t(apply(P, 1L, cumsum))
    idx <- rowSums(u > cp) + 1L; idx[idx > J] <- J
    idx
  }
  ## the systematic utilities of a random-intercept fit, once (a one-node rule:
  ## only V is needed, not the population-averaged probabilities)
  V0 <- if (identical(object$family, "remnl"))
    pc_re_probs_cpp(object$theta, D$X, D$W, D$off, D$avail, rep(0L, n), J, D$base - 1L, D$str$uidx, D$str$ustart, matrix(0, 1L, J - 1L), 0)$V
  out <- lapply(seq_len(nsim), function(s) {
    P <- if (identical(object$family, "remnl")) {
      ## the model draws one intercept vector per unit, then the unit's cases given it
      sig <- exp(object$theta[D$idx$str]); G <- nlevels(D$unit); ui <- as.integer(D$unit)
      U <- matrix(stats::rnorm(G * length(sig)), G) * rep(sig, each = G)
      V <- V0
      nb <- setdiff(seq_len(J), D$base)
      V[, nb] <- V[, nb] + U[ui, , drop = FALSE]
      V[D$avail == 0L] <- -Inf
      E <- exp(V - apply(V, 1L, max)); E / rowSums(E)
    } else object$fitted.values
    factor(object$levels[draw(P)], levels = object$levels)
  })
  names(out) <- paste0("sim_", seq_len(nsim))
  out <- as.data.frame(out, stringsAsFactors = FALSE)
  rownames(out) <- object$case_id
  out
}
