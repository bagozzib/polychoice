## ---------------------------------------------------------------------------
## Model fit and model comparison: fit statistics, classification, the Vuong
## test, likelihood-ratio tests, side-by-side comparison.
## ---------------------------------------------------------------------------

## the alternative-specific-constants-only (null) fit on the same cases and
## availability: the McFadden baseline shared by mnl and imnl fits
.pc_null_loglik <- function(object) {
  D <- object$design
  ro <- identical(D$family, "rologit")
  if (is.null(D$fe) && !D$asc && !("(Intercept)" %in% colnames(D$X))) {
    ## no constants in the model: the null is the equiprobable choice over the
    ## available set (for rankings, at every stage)
    if (ro) return(sum(D$w * pc_ro_lli_cpp(numeric(0), matrix(0, D$n, 0L), matrix(0, D$n * D$J, 0L), matrix(0, D$n, D$J), D$avail,
                                           D$str$rank, D$J, D$base - 1L, 0L)))
    return(sum(D$w * log(1 / rowSums(D$avail))))
  }
  D0 <- D
  D0$firth <- NULL; D0$prior <- NULL            # the null of a penalized fit is the maximum-likelihood constants-only fit
  D0$X <- matrix(1, D$n, 1L, dimnames = list(NULL, "(Intercept)")); D0$p <- 1L
  D0$W <- matrix(0, D$n * D$J, 0L); D0$q <- 0L
  D0$Z <- matrix(0, D$n, 0L); D0$r <- 0L; D0$k <- 0L; D0$kc <- integer(0)
  nonbase <- setdiff(seq_len(D$J), D$base)
  D0$par_names <- paste0(D$levels[nonbase], ":(Intercept)"); D0$K <- length(nonbase)
  D0$idx <- list(beta = matrix(seq_len(D0$K), nrow = 1L), gamma = integer(0), delta = integer(0), str = integer(0))
  D0$blocks <- rep("outcome", D0$K); D0$fe_map <- NULL
  if (!ro) { D0$family <- "mnl"; D0$str <- .pc_structural("mnl", NULL, D$levels, D$base) }     # a ranking keeps its rank matrix
  .pc_refit(D0)$loglik
}

#' Fit statistics for a choice model
#'
#' Log-likelihood, information criteria, and the pseudo-R-squared measures in
#' common use for multinomial models, all relative to the same null: the model
#' with alternative-specific constants only (equal choice probabilities over the
#' available alternatives when the model has no constants). Inflated and
#' uninflated fits share this null, so their McFadden measures are comparable.
#'
#' The null fit keeps the fit's availability pattern, including the
#' alternatives that unit fixed effects switched off for never-chosen pairs.
#' For a random-intercept fit `n` is the number of cases (not units), and so
#' is the \eqn{n} of the BIC. For a penalized fit (`firth = TRUE` or a prior)
#' the log-likelihood is the ordinary one at the penalized estimates.
#' The information criteria and R-squared measures built on it are therefore
#' not maximum-likelihood ones (the print says so). The null is the
#' maximum-likelihood constants-only fit, and `loglik_penalized` holds the
#' maximized objective.
#'
#' @param object A `"choicemodel"` fit.
#' @return A one-row data frame of class `"fit_stats"` with `loglik`,
#'   `loglik_null`, `df`, `n` (cases), `AIC`, `BIC` (with the weight total as
#'   \eqn{n}), `mcfadden` (\eqn{1 - \ell/\ell_0}), `mcfadden_adj`
#'   (\eqn{1 - (\ell - K)/\ell_0}), `cox_snell`, `nagelkerke`, `count` (the
#'   weighted share of cases whose modal predicted alternative is the chosen
#'   one), and `count_adj` (the proportional reduction in error relative to
#'   always predicting the most frequent alternative, Long 1997).
#' @references Long, J.S. (1997). \emph{Regression Models for Categorical and
#'   Limited Dependent Variables}. Sage.
#' @examples
#' set.seed(1)
#' d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
#' fit_stats(mnl(y ~ x1 + x2, data = d))
#' @seealso [classification()], [compare_models()], [vuong()]
#' @family model comparison
#' @export
fit_stats <- function(object) {
  if (!inherits(object, "choicemodel")) stop("'object' must be a fitted polychoice model.", call. = FALSE)
  ll <- object$loglik; ll0 <- .pc_null_loglik(object); K <- object$df
  nw <- object$nobs_weighted; w <- object$weights
  P <- object$fitted.values
  modal <- apply(P, 1L, which.max)
  correct <- sum(w * (modal == object$y)) / nw
  shares <- tapply(w, factor(object$y, levels = seq_len(object$design$J)), sum); shares[is.na(shares)] <- 0
  pmax_share <- max(shares) / nw
  out <- data.frame(loglik = ll, loglik_null = ll0, df = K, n = object$n,
                    AIC = -2 * ll + 2 * K, BIC = -2 * ll + log(nw) * K,
                    mcfadden = 1 - ll / ll0, mcfadden_adj = 1 - (ll - K) / ll0,
                    cox_snell = 1 - exp(2 * (ll0 - ll) / nw),
                    nagelkerke = (1 - exp(2 * (ll0 - ll) / nw)) / (1 - exp(2 * ll0 / nw)),
                    count = correct, count_adj = (correct - pmax_share) / (1 - pmax_share))
  if (.pc_penalized(object)) { out$loglik_penalized <- object$loglik_penalized; attr(out, "penalty") <- .pc_penalty_label(object) }
  class(out) <- c("fit_stats", "data.frame")
  out
}

#' @export
print.fit_stats <- function(x, digits = 3, ...) {
  cat(sprintf("Log-likelihood %.*f (null %.*f), %d parameters, %d cases; AIC %.*f, BIC %.*f\n", digits, x$loglik, digits, x$loglik_null,
              x$df, x$n, digits, x$AIC, digits, x$BIC))
  if (!is.null(x$loglik_penalized))
    cat(sprintf("Penalized fit (%s): the log-likelihood is at the penalized estimates (penalized log-likelihood %.*f); AIC, BIC and the R2 measures are not those of maximum likelihood.\n",
                attr(x, "penalty") %||% "Firth", digits, x$loglik_penalized))
  cat(sprintf("McFadden R2 %.*f (adjusted %.*f); Cox-Snell %.*f; Nagelkerke %.*f\n", digits, x$mcfadden, digits, x$mcfadden_adj,
              digits, x$cox_snell, digits, x$nagelkerke))
  cat(sprintf("Count R2 %.*f (adjusted: %.*f proportional reduction in error over the modal alternative)\n", digits, x$count, digits, x$count_adj))
  invisible(x)
}

#' Classification table of a choice model
#'
#' Cross-tabulates the chosen alternative against the modal predicted
#' alternative (weighted by the frequency weights), with the share correctly
#' classified and the proportional reduction in error relative to always
#' predicting the most frequent alternative. Modal classification understates a
#' model's value for minority alternatives: a case with probability 0.4 on a
#' rare alternative is classified to the common one even though 0.4 may be far
#' above its base rate; see the probability-based measures in [fit_stats()].
#'
#' @param object A `"choicemodel"` fit.
#' @param newdata Optional data frame (same layout as the fit) with the
#'   response; default the estimation cases. A hold-out set may lack some
#'   alternative or covariate level; only the response and the predictors are
#'   read from it.
#' @param weights Optional frequency weights for `newdata` (a column name or a
#'   vector aligned to its rows); the estimation cases use the fit's weights.
#' @param offset Offset for `newdata` when the fit has one (see
#'   [predict.choicemodel()]).
#' @return An object of class `"pc_classification"`: a list with `table`
#'   (observed in rows, predicted in columns), `correct`, and `pre`.
#' @examples
#' set.seed(1)
#' d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
#' classification(mnl(y ~ x1 + x2, data = d))
#' @seealso [fit_stats()]
#' @family model comparison
#' @export
classification <- function(object, newdata = NULL, weights = NULL, offset = NULL) {
  if (!inherits(object, "choicemodel")) stop("'object' must be a fitted polychoice model.", call. = FALSE)
  if (is.null(newdata)) {
    ## in-sample: the fitted probabilities (for a jackknife-corrected fit, those
    ## of the maximum-likelihood fit, as in fit_stats())
    obs <- object$levels[object$y]; pred <- object$levels[max.col(object$fitted.values, ties.method = "first")]; w <- object$weights
  } else {
    ## the observed alternatives are read off newdata directly (no estimation
    ## design is built: a hold-out set may lack an alternative or a level)
    if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
    F <- object$design$F
    yraw <- eval(attr(F, "lhs")[[1L]], newdata, environment(F))
    if (object$layout == "wide") {
      obs <- as.character(if (is.factor(yraw)) yraw else yraw)
      bad <- setdiff(unique(obs), object$levels)
      if (length(bad)) stop("'newdata' has response value(s) not among the fit's alternatives: ", paste(bad, collapse = ", "), ".", call. = FALSE)
      w <- if (is.null(weights)) rep(1, nrow(newdata)) else .pc_aligned(weights, newdata, "weights")
    } else if (identical(object$family, "rologit")) {
      ## rankings: the first choice is the alternative of smallest positive rank
      ## (a tie at the top: the first of the tied in level order, as in the fit)
      ids <- newdata[[object$id]]; cases <- unique(ids); ci <- match(ids, cases)
      altf <- factor(newdata[[object$alt]], levels = object$levels)
      if (anyNA(altf)) stop("'newdata' has alternatives that are not in the model.", call. = FALSE)
      rk <- as.numeric(yraw); rk[!is.finite(rk) | rk <= 0] <- Inf
      ord <- order(ci, rk, as.integer(altf))
      firsts <- ord[!duplicated(ci[ord])]
      if (any(!is.finite(rk[firsts]))) stop("Some case(s) in 'newdata' rank no alternative.", call. = FALSE)
      obs <- character(length(cases)); obs[ci[firsts]] <- as.character(altf[firsts])
      w <- if (is.null(weights)) rep(1, length(cases)) else .pc_aligned(weights, newdata, "weights")[match(seq_along(cases), ci)]
    } else {
      ids <- newdata[[object$id]]; cases <- unique(ids); ci <- match(ids, cases)
      ch <- if (is.factor(yraw)) as.character(yraw) else yraw
      ch <- as.numeric(ch)
      nch <- tabulate(ci[ch == 1], nbins = length(cases))
      if (any(nch != 1L)) stop(sum(nch != 1L), " case(s) in 'newdata' do not have exactly one chosen row.", call. = FALSE)
      obs <- character(length(cases)); obs[ci[ch == 1]] <- as.character(newdata[[object$alt]][ch == 1])
      bad <- setdiff(unique(obs), object$levels)
      if (length(bad)) stop("'newdata' has chosen alternative(s) not in the fit: ", paste(bad, collapse = ", "), ".", call. = FALSE)
      w <- if (is.null(weights)) rep(1, length(cases)) else .pc_aligned(weights, newdata, "weights")[match(seq_along(cases), ci)]
    }
    pred <- as.character(predict(object, newdata = newdata, type = "class", offset = offset))
  }
  lev <- object$levels
  tab <- stats::xtabs(w ~ factor(obs, levels = lev) + factor(pred, levels = lev))
  names(dimnames(tab)) <- c("observed", "predicted")
  correct <- sum(diag(tab)) / sum(tab)
  pmax_share <- max(rowSums(tab)) / sum(tab)
  structure(list(table = tab, correct = correct, pre = (correct - pmax_share) / (1 - pmax_share)), class = "pc_classification")
}

#' @export
print.pc_classification <- function(x, digits = 3, ...) {
  print(x$table)
  cat(sprintf("\nCorrectly classified: %.*f; proportional reduction in error over the modal alternative: %.*f\n", digits, x$correct, digits, x$pre))
  invisible(x)
}

## the Vuong statistic with the raw, AIC-, and BIC-corrected numerators; l1,
## l2 per-case log-likelihoods, w frequency weights, ddf = df1 - df2
.pc_vuong_stat <- function(l1, l2, w, ddf) {
  m <- l1 - l2
  nw <- sum(w)
  mbar <- sum(w * m) / nw
  s <- sqrt(sum(w * (m - mbar)^2) / nw)
  if (!is.finite(s) || s <= 0) stop("The per-case log-likelihoods of the two models are identical; the Vuong statistic is undefined.", call. = FALSE)
  num <- c(raw = nw * mbar, AIC = nw * mbar - ddf, BIC = nw * mbar - ddf * log(nw) / 2)
  list(stat = num / (sqrt(nw) * s), mean_diff = mbar, sd_diff = s, n = nw)
}

#' Vuong test for non-nested or overlapping choice models
#'
#' Vuong's (1989) test comparing two fits to the same cases; positive values
#' favour the first model. It is the comparison the inflation literature uses
#' between an inflated model and its uninflated counterpart (Bagozzi and
#' Marchetti 2017; Greene 1994), where the restriction \eqn{\phi = 1} sits at
#' the boundary so that the likelihood-ratio statistic is not chi-square. Its
#' standard normal reference does not hold there: when the uninflated model is
#' true the two models are nested at the boundary, the per-case differences
#' degenerate, and the statistic is not N(0, 1) (Wilson 2015, on the same use
#' for zero-inflated counts; a size study under the MNL null is part of the
#' package's verification, `verification/vuong_size_mc.R` in the development
#' repository). The calibrated alternative is the parametric bootstrap of the
#' likelihood-ratio and Vuong statistics under the restricted model,
#' [boot_lr_test()]. Between strictly non-nested models the normal reference
#' holds; for overlapping ones (sharing some covariates) Vuong's procedure
#' first tests whether the per-case differences have zero variance (the
#' models are then equivalent on these data), which this function reports
#' only as an error when the differences are identical. With
#' frequency weights the statistic equals the one computed on duplicated cases.
#' Reported with the raw numerator and with the AIC and BIC corrections for the
#' difference in parameter counts. When the two fits are practically the same
#' model on the data (one at the other's boundary, such as an inflated fit with
#' \eqn{\phi} at 1), the per-case differences are round-off, every statistic
#' divides by their spread, and the corrected ones run to minus millions: when
#' that spread totals less than 0.01 (\eqn{\sqrt{n}\,\hat\omega < 0.01}) the
#' statistics are `NA` and a note says the models cannot be told apart. When either model is a [remnl()] fit, whose
#' likelihood is per unit, the contributions of both are summed by unit and
#' the statistic (and the \eqn{n} of its BIC correction) is over units, while
#' [fit_stats()] counts cases.
#'
#' The statistics assume independent cases: a fit's `cluster` does not enter them.
#'
#' @param m1,m2 Two `"choicemodel"` fits to the same cases.
#' @return A data frame of class `"vuong"` with one row per correction:
#'   `statistic`, `p_m1_better` (one-sided, m1 preferred), `p_m2_better`, and
#'   `p_two_sided`.
#' @references
#' Bagozzi, B.E. and Marchetti, K. (2017). Distinguishing occasional abstention
#' from routine indifference in models of vote choice. \emph{Political Science
#' Research and Methods}, 5(2), 277-294. \doi{10.1017/psrm.2015.42}
#'
#' Greene, W.H. (1994). Accounting for excess zeros and sample selection in
#' Poisson and negative binomial regression models. Working Paper EC-94-10,
#' Stern School of Business, New York University.
#'
#' Vuong, Q.H. (1989). Likelihood ratio tests for model selection and
#' non-nested hypotheses. \emph{Econometrica}, 57, 307-333. \doi{10.2307/1912557}
#'
#' Wilson, P. (2015). The misuse of the Vuong test for non-nested models to
#' test for zero-inflation. \emph{Economics Letters}, 127, 51-53.
#' \doi{10.1016/j.econlet.2014.12.029}
#' @examples
#' set.seed(2)
#' d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
#'           delta = c(0.5, 1.2), inflate = "a")
#' vuong(imnl(y ~ x1 + x2 | z1, data = d, inflate = "a"), mnl(y ~ x1 + x2, data = d))
#' @seealso [lr_test()], [compare_models()]
#' @family model comparison
#' @export
vuong <- function(m1, m2) {
  if (!inherits(m1, "choicemodel") || !inherits(m2, "choicemodel")) stop("'m1' and 'm2' must be fitted polychoice models.", call. = FALSE)
  .pc_no_firth(m1, "vuong()"); .pc_no_firth(m2, "vuong()")
  if (xor(identical(m1$family, "rologit"), identical(m2$family, "rologit")))
    stop("A rank-ordered logit is fit to rankings and the other model to single choices; their likelihoods are not comparable.", call. = FALSE)
  if (identical(m1$family, "rologit") && !identical(m1$design$str$rank, m2$design$str$rank))
    stop("The two rank-ordered fits use different rankings (a different 'depth'?); their likelihoods are not comparable.", call. = FALSE)
  al <- .pc_vuong_align(m1, m2)
  if (is.null(al))
    stop("The two models must be fit to the same cases (same rows, response, and weights).", call. = FALSE)
  l1 <- al$l1; l2 <- al$l2; w <- al$w
  if (m1$family == "remnl" || m2$family == "remnl") {
    ## the random-intercept likelihood is per unit: aggregate the other model's
    ## contributions over the same units (the case check above aligns them)
    if (m1$family == "remnl" && m2$family == "remnl" && !identical(as.character(m1$design$unit), as.character(m2$design$unit)))
      stop("Two random-intercept models must share the unit variable.", call. = FALSE)
    ui <- as.integer(if (m1$family == "remnl") m1$design$unit else m2$design$unit)
    agg <- function(m, l) if (m$family == "remnl") l else as.vector(rowsum(l, ui))
    l1 <- agg(m1, l1); l2 <- agg(m2, l2)
    w <- as.vector(rowsum(w, ui)) / tabulate(ui)
  }
  vs <- .pc_vuong_stat(l1, l2, w, m1$df - m2$df)
  ## fits that are practically the same model (one at the other's boundary): the
  ## spread of the per-case differences is round-off, and every statistic
  ## divides by it (the corrected ones then run to minus millions)
  spread <- sqrt(vs$n) * vs$sd_diff
  equiv <- spread < 0.01
  if (equiv) vs$stat[] <- NA_real_
  infl_pair <- xor(isTRUE(m1$inflated), isTRUE(m2$inflated)) && identical(m1$family, "mnl") && identical(m2$family, "mnl")
  out <- data.frame(correction = names(vs$stat), statistic = unname(vs$stat),
                    p_m1_better = stats::pnorm(unname(vs$stat), lower.tail = FALSE),
                    p_m2_better = stats::pnorm(unname(vs$stat)),
                    p_two_sided = 2 * stats::pnorm(-abs(unname(vs$stat))), row.names = NULL)
  attr(out, "models") <- c(deparse(substitute(m1))[1L], deparse(substitute(m2))[1L])
  attr(out, "labels") <- c(.pc_title(m1), .pc_title(m2))
  equiv_note <- if (equiv) paste0("The two fits give practically the same per-case log-likelihoods (the spread of their differences totals ",
                                  format(spread, digits = 2), " < 0.01): one is at the other's boundary, the models cannot be told apart on these ",
                                  "data, and the statistics are undefined (NA). ") else ""
  if (equiv) attr(out, "note") <- equiv_note
  if (infl_pair) attr(out, "note") <- paste0(equiv_note, "The inflated multinomial logit and the multinomial logit are nested at the boundary phi = 1, where the normal reference ",
                                             "does not hold for testing the uninflated model; boot_lr_test() calibrates the statistic.")
  class(out) <- c("vuong", "data.frame")
  out
}

#' @export
print.vuong <- function(x, digits = 3, ...) {
  nm <- attr(x, "models"); lb <- attr(x, "labels")
  cat("Vuong test:  m1 =", nm[1L], "[", lb[1L], "]\n             m2 =", nm[2L], "[", lb[2L], "]\n")
  cat("positive statistics favor m1\n\n")
  y <- as.data.frame(x); y[-1L] <- lapply(y[-1L], round, digits)
  print.data.frame(y, row.names = FALSE)
  if (!is.null(attr(x, "note"))) cat("Note:", attr(x, "note"), "\n")
  invisible(x)
}

## The per-case log-likelihood contributions of two fits on the same cases;
## NULL when the cases do not match.
.pc_vuong_align <- function(m1, m2) {
  if (m1$n == m2$n && identical(m1$y, m2$y) && identical(m1$case_id, m2$case_id) && isTRUE(all.equal(m1$weights, m2$weights)))
    return(list(l1 = m1$lli, l2 = m2$lli, w = m1$weights))
  NULL
}

## The nesting checks shared by lr_test() and boot_lr_test(). Returns the kind
## of restriction: "interior" (the chi-square reference applies), "sd" (the
## multinomial logit inside a random-intercept fit: sd = 0 is on the boundary and the
## chi-square p-value is conservative), or "inflation" (the multinomial logit inside the
## inflated multinomial logit: phi = 1 is on the boundary and the split coefficients are
## not identified there, so only the parametric bootstrap applies; accepted
## only when allow_inflation = TRUE).
.pc_nesting <- function(restricted, full, allow_inflation = FALSE) {
  if (!inherits(restricted, "choicemodel") || !inherits(full, "choicemodel")) stop("Both arguments must be fitted polychoice models.", call. = FALSE)
  if (!identical(restricted$fe, full$fe))
    stop("One model has unit fixed effects and the other does not; the fixed effects' -Inf constants and the incidental-parameters asymptotics leave no chi-square reference for this comparison.", call. = FALSE)
  if (restricted$n != full$n || !identical(restricted$y, full$y) || !identical(restricted$case_id, full$case_id))
    stop("The two models must be fit to the same cases (same rows and response).", call. = FALSE)
  if (!isTRUE(all.equal(restricted$weights, full$weights))) stop("The two models must be fit with the same weights.", call. = FALSE)
  if (!identical(restricted$design$avail, full$design$avail))
    stop("The two models are fit with different choice sets (availability); they are not nested.", call. = FALSE)
  if (identical(restricted$family, full$family) && !identical(restricted$design$str[c("nodes", "draws")], full$design$str[c("nodes", "draws")]))
    stop("The two fits use different numerical approximations (quadrature nodes or simulation draws); refit with the same settings.", call. = FALSE)
  if (xor(identical(restricted$family, "rologit"), identical(full$family, "rologit")))
    stop("A rank-ordered logit is fit to rankings and the other model to single choices; the models are not nested.", call. = FALSE)
  if (identical(full$family, "rologit") && !identical(restricted$design$str$rank, full$design$str$rank))
    stop("The two rank-ordered fits use different rankings (a different 'depth'?); they are not nested.", call. = FALSE)
  nonstr <- function(m) names(m$coefficients)[m$blocks %in% c("outcome", "generic", "fixed effect")]
  subset_of <- function(a, b) all(a %in% b)
  if (restricted$inflated != full$inflated) {
    ## one inflated alternative with an intercept-only split: the captive odds
    ## at zero are a regular boundary point (DESIGN.md Phase 4e)
    const_split <- !restricted$inflated && restricted$family == "mnl" && full$family == "mnl" && length(full$inflate) == 1L &&
      identical(colnames(full$design$Z), "(Intercept)")
    if (!(allow_inflation && !restricted$inflated && restricted$family == "mnl" && full$family == "mnl") && !const_split)
      stop("An mnl() vs imnl() comparison is not a standard nested test (phi = 1 is on the boundary of the parameter space and the split coefficients are not identified there); use boot_lr_test(). ",
           "With one inflated alternative and an intercept-only split (y ~ x | 1) the boundary is regular and lr_test() applies.", call. = FALSE)
    if (!subset_of(nonstr(restricted), nonstr(full)))
      stop("The multinomial logit has outcome covariates the inflated model lacks; the models are not nested.", call. = FALSE)
    return(if (allow_inflation) "inflation" else "captive")
  }
  if (restricted$inflated && full$inflated && !identical(restricted$inflate, full$inflate) && all(restricted$inflate %in% full$inflate) &&
      identical(restricted$link, "logit") && identical(full$link, "logit")) {
    ## fewer inflated alternatives: a captive share of zero, on the boundary
    if (!allow_inflation)
      stop("A model with fewer inflated alternatives is nested at the boundary (a captive share of zero, where its split coefficients are not identified); use boot_lr_test().", call. = FALSE)
    if (!subset_of(nonstr(restricted), nonstr(full)) || !identical(colnames(restricted$design$Z), colnames(full$design$Z)))
      stop("The two inflated models differ in their covariates; they are not nested.", call. = FALSE)
    return("inflation")
  }
  if (restricted$inflated && (!identical(restricted$inflate, full$inflate) || !identical(restricted$link, full$link)))
    stop("The two inflated models differ in the inflated alternative or the link; they are not nested.", call. = FALSE)
  if (!identical(restricted$family, full$family)) {
    if (!(restricted$family == "mnl" && full$family %in% c("nlogit", "hev", "remnl")))
      stop("Only the multinomial logit nested in a nested logit, an HEV, or a random-intercept fit can be compared across families; the multinomial probit does not nest the multinomial logit (use vuong()).", call. = FALSE)
    if (!identical(nonstr(restricted), nonstr(full)))
      stop("The two fits must have the same covariates (the restriction is lambda = 1 / scale = 1 / sd = 0).", call. = FALSE)
    return(if (full$family == "remnl") "sd" else "interior")
  }
  ## the same family: the restricted parameters are a subset of the full ones,
  ## with the family's structure shared
  Dr <- restricted$design; Df <- full$design
  if (full$family == "nlogit") {
    if (!identical(Dr$str$nest, Df$str$nest)) stop("The two nested logits have different nests; they are not nested.", call. = FALSE)
    equal_vs_free <- isTRUE(Dr$str$equal) && !isTRUE(Df$str$equal)
    if (!subset_of(nonstr(restricted), nonstr(full)) || (!equal_vs_free && !subset_of(names(restricted$coefficients), names(full$coefficients))))
      stop("The restricted model has parameters the full model lacks; the models are not nested.", call. = FALSE)
    return("interior")
  }
  if (full$family == "remnl" && !identical(as.character(Dr$unit), as.character(Df$unit)))
    stop("The two random-intercept fits must share the unit variable.", call. = FALSE)
  if (!subset_of(names(restricted$coefficients), names(full$coefficients)))
    stop("The restricted model has parameters the full model lacks (", paste(utils::head(setdiff(names(restricted$coefficients), names(full$coefficients)), 3L), collapse = ", "),
         "); the models are not nested.", call. = FALSE)
  "interior"
}

## the restricted estimates embedded in the full model's parameter space at the
## restriction: zero for the added coefficients, lambda = 1 and scale = 1 (zero
## on the log scale), an sd near zero, and the split equation at phi = 1 - 1e-8.
## A refit started there cannot end below the restricted log-likelihood (up to
## the 1e-8 left in the split), so a replicate LR is not under-maximized at
## the restriction.
.pc_embed <- function(theta_r, Dr, Df) {
  th <- stats::setNames(numeric(Df$K), Df$par_names)
  common <- intersect(Dr$par_names, Df$par_names)
  th[common] <- theta_r[match(common, Dr$par_names)]
  if (Df$family == "nlogit" && Dr$family == "nlogit" && isTRUE(Dr$str$equal) && !isTRUE(Df$str$equal))
    th[Df$idx$str] <- theta_r[Dr$idx$str][1L]
  if (Df$family == "remnl" && Dr$family != "remnl") th[Df$idx$str] <- -7
  if (Df$k > 0L) {
    kr <- if (Dr$k > 0L) (if (length(Dr$kc)) Dr$kc else Dr$k) else integer(0)
    kf <- if (length(Df$kc)) Df$kc else Df$k
    nm_f <- function(c, term) if (length(kf) > 1L) paste0("infl_", Df$levels[c], ":", term) else paste0("infl_", term)
    ## a single-inflated fit inside a several-inflated model: its split block is
    ## that alternative's block of the larger model
    if (length(kf) > 1L && length(kr) == 1L)
      for (term in colnames(Dr$Z)) { from <- paste0("infl_", term); to <- nm_f(kr, term); if (to %in% names(th)) th[to] <- theta_r[match(from, Dr$par_names)] }
    ## the alternatives the restricted model does not inflate: captive share 1e-8
    for (c in setdiff(kf, kr)) {
      th[nm_f(c, colnames(Df$Z))] <- 0
      ic <- nm_f(c, "(Intercept)")
      if (ic %in% names(th)) th[ic] <- if (identical(Df$link, "probit")) stats::qnorm(1 - 1e-8) else stats::qlogis(1 - 1e-8)
    }
  }
  unname(th)
}

## the maximum over the engine's own multistart (exactly as for the original
## fit) and any extra starts; NULL when every attempt fails
.pc_best_fit <- function(D, extra = list(), maxit = 1000L, reltol = 1e-10) {
  quiet <- function(expr) tryCatch(suppressWarnings(suppressMessages(expr)), error = function(e) { .pc_pass_limit(e); NULL })
  f <- quiet(.pc_fit_design(D, se = "none", maxit = maxit, reltol = reltol, check = FALSE))
  best <- if (is.null(f)) NULL else list(par = unname(f$theta), value = -f$loglik, converged = isTRUE(f$converged))
  ## a concave log-likelihood (the multinomial logit and the rank-ordered logit) has one
  ## maximum: the further starts matter only when the fit did not converge
  ## (separation, where the statistic takes the supremum anyway). They took
  ## 92 of 114 s of a ranking boot_lr_test() (local review 2026-11-06)
  if (length(extra) && !(.pc_concave(D) && !is.null(best) && best$converged)) {
    fn <- .pc_fns(D); sc <- .pc_scales(D)
    for (s in extra) {
      op <- quiet(.pc_optimize(fn, s, sc, maxit, reltol))
      if (!is.null(op) && is.finite(op$value) && (is.null(best) || op$value < best$value - 1e-10))
        best <- list(par = op$par, value = op$value, converged = isTRUE(op$converged))
    }
  }
  best
}

#' Likelihood-ratio test for nested choice models
#'
#' For a restricted model that is a special case of the full one: the same
#' estimator (both [mnl()] or both [imnl()] with the same inflated alternative
#' and link) on the same cases, with the restricted parameters a subset of the
#' full ones (checked by name; two nested logits must share their nests, and a
#' common inclusive value, `equal_lambda = TRUE`, nests the free ones); or an
#' [mnl()] fit against the [nlogit()], [hev()], or [remnl()] fit with the
#' same covariates, which tests the
#' inclusive-value parameters against 1 (the nested-logit test of IIA), the
#' scale parameters against 1 (homoskedastic errors), or the random-intercept
#' standard deviations against 0.
#'
#' Separation. When the data of an [mnl()] or [rologit()] fit (without unit
#' fixed effects or inflation) are separated, the maximum-likelihood estimates
#' do not exist and the fit stops somewhere along the separating direction;
#' the statistic then compares the suprema of the two log-likelihoods, which
#' [separation()] makes computable: 0 under complete separation, and under
#' quasi-complete separation the value reached by moving along the
#' separating direction and re-maximizing. Likelihood-ratio tests stay valid
#' under separation where Wald tests do not (Rainey 2024), but with few cases
#' or rankings per coefficient the chi-square reference can be far off;
#' [boot_lr_test()] calibrates it. [term_test()] and [combine_test()] use the
#' same suprema.
#'
#' Penalized fits. Two penalized fits (`firth = TRUE`, or the same prior) are
#' compared by the penalized likelihood ratio: the restricted maximum is taken
#' over the full model's penalized log-likelihood (Heinze and Schemper 2001).
#' The priors count as the same when they share the family, degrees of
#' freedom, autoscaling and constants' scale, and the scales of the
#' coefficients both fits have. A scale named for a coefficient only the full
#' model has is ignored by the restricted fit, with a message. A penalized
#' fit against an unpenalized one, or two different penalties, are refused.
#'
#' Boundaries. Two of these restrictions sit on the boundary of the parameter
#' space and the chi-square reference does not hold exactly. For the
#' random-intercept model (\eqn{\sigma = 0}) the null distribution is a mixture
#' of chi-squares with fewer degrees of freedom, so the reported p-value is
#' conservative and the test can be used as such (a message says so). For the
#' inflated model against the multinomial logit (\eqn{\phi = 1}) the split-equation slopes
#' are not identified under the null, so no simple correction exists: the
#' function refuses that pair; [boot_lr_test()] calibrates it by simulation
#' (the normal reference of [vuong()] is not justified there either). The
#' exception is an inflated model with one inflated alternative and an
#' intercept-only split (`y ~ x | 1`, the dogit with a constant captive
#' share): its captive odds at zero are a regular boundary point, and the
#' reference is the 50:50 mixture of chi-square(q) and chi-square(q + 1), q
#' the further restrictions (Self and Liang 1987), as in
#' `captivity_test(type = "lr")`. The
#' multinomial probit does not nest the multinomial logit: use [vuong()]. Two fits of the
#' same family other than the multinomial logit (two HEV, nested-logit, multinomial probit, or
#' random-intercept fits) are accepted when the restricted parameters are a
#' subset of the full ones and the fits use the same quadrature nodes or
#' simulation draws. The function also refuses fits with different choice
#' sets, and a fixed-effects fit against one without fixed effects; two
#' fixed-effects fits are compared with a message when the panels are short,
#' since the incidental-parameters bias shifts the statistic.
#'
#' Clustered data. The statistic assumes independent cases: a fit's `cluster` does not enter it.
#'
#' @param restricted,full Two `"choicemodel"` fits on the same cases.
#' @return A data frame of class `"lr_test"` with the statistic, degrees of
#'   freedom, and p-value; the attribute `"suprema"` is `TRUE` when a
#'   supremum replaced a separated fit's value. For the constant-share
#'   inflated multinomial logit, `df` is q + 1 (the share and the q further
#'   restrictions), the attribute `"mixture"` holds q, and `p.value` is `NA`,
#'   with a warning, when the inflated fit is degenerate (outcome
#'   coefficients diverging, or a saddle).
#' @references
#' Rainey, C. (2024). Hypothesis tests under separation. *Political
#' Analysis*, 32, 172-185. \doi{10.1017/pan.2023.28}
#'
#' Self, S.G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
#' likelihood estimators and likelihood ratio tests under nonstandard
#' conditions. *Journal of the American Statistical Association*, 82,
#' 605-610. \doi{10.1080/01621459.1987.10478472}
#' @examples
#' set.seed(1)
#' d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
#' lr_test(mnl(y ~ x1, data = d), mnl(y ~ x1 + x2, data = d))
#' @seealso [boot_lr_test()], [vuong()], [compare_models()]
#' @family model comparison
#' @export
lr_test <- function(restricted, full) {
  kind <- .pc_nesting(restricted, full)
  if (full$df <= restricted$df) stop("'full' must have more parameters than 'restricted'.", call. = FALSE)
  if (.pc_penalized(restricted) || .pc_penalized(full)) {
    ## penalized fits: the penalized likelihood ratio, the restricted maximum
    ## taken over the full model's penalized log-likelihood (Heinze and Schemper
    ## 2001); the separately penalized restricted fit carries another penalty
    if (!(.pc_penalized(restricted) && .pc_penalized(full)))
      stop("One fit is penalized (firth = TRUE or a prior) and the other is not; their likelihoods are not comparable. Fit both with the same 'firth' or 'prior'.", call. = FALSE)
    if (!.pc_same_penalty(restricted, full))
      stop("The two fits are penalized differently (different penalties or prior scales); their penalized likelihoods are not comparable. Fit both with the same 'firth' or 'prior'.", call. = FALSE)
    extra <- match(setdiff(names(full$coefficients), names(restricted$coefficients)), names(full$coefficients))
    Rm <- matrix(0, length(extra), full$df); Rm[cbind(seq_along(extra), extra)] <- 1
    stat <- .pc_firth_lr(full, Rm, start = .pc_embed(restricted$theta, restricted$design, full$design))
    df <- full$df - restricted$df
    out <- data.frame(statistic = max(stat, 0), df = df, p.value = stats::pchisq(max(stat, 0), df, lower.tail = FALSE))
    attr(out, "penalized") <- TRUE; attr(out, "penalty") <- .pc_penalty_label(full)
    class(out) <- c("lr_test", "data.frame")
    return(out)
  }
  if (kind == "sd")
    message("The restriction sd = 0 is on the boundary of the parameter space: the chi-square p-value is conservative (the null distribution is a mixture of chi-squares with fewer degrees of freedom); boot_lr_test() gives a calibrated p-value.")
  df <- full$df - restricted$df
  if (!is.null(full$fe) && isTRUE(full$fe_info$Tbar < 10))
    message("Both fits have unit fixed effects with few cases per unit (median ", full$fe_info$Tbar, "): the incidental-parameters bias shifts the ",
            "LR statistic, so the chi-square p-value can be far too small in short panels.")
  ## suprema (DESIGN.md Phase 4c): under separation a fit's value is where the optimizer stopped
  sp <- .pc_sup_pair(restricted$design, restricted$theta, restricted$loglik, full$design, full$theta, full$loglik)
  stat <- 2 * (sp$f - sp$r)
  ## an alternative whose probability vanishes at a supremum (never chosen,
  ## its comparisons separated) contributes no coefficients to it: the degrees
  ## of freedom count the parameters that enter the suprema (local review
  ## 2026-11-06, finding 10)
  absent <- character(0)
  if (kind != "captive") {
    ab_f <- .pc_absent_alts(full$design); ab_r <- .pc_absent_alts(restricted$design)
    if (length(ab_f) || length(ab_r)) {
      dfe <- .pc_k_effective(full$design, ab_f) - .pc_k_effective(restricted$design, ab_r)
      if (dfe >= 1L && dfe != df) { df <- dfe; absent <- full$design$levels[union(ab_f, ab_r)] }
    }
  }
  if (kind == "captive") {
    cr <- .pc_captive_rule(stat, full, restricted$loglik, df - 1L)
    if (cr$below) warning("The restricted model has the higher log-likelihood; the fits may not be nested or one did not converge.", call. = FALSE)
    if (cr$degenerate) warning("The inflated fit is degenerate (", cr$why, "); the chi-bar reference does not apply. Use boot_lr_test().", call. = FALSE)
    out <- data.frame(statistic = cr$statistic, df = df, p.value = if (cr$degenerate) NA_real_ else .pc_chibar_p(cr$statistic, df - 1L))
    attr(out, "mixture") <- df - 1L
  } else {
    if (stat < -1e-6) warning("The restricted model has the higher log-likelihood; the fits may not be nested or one did not converge.", call. = FALSE)
    out <- data.frame(statistic = max(stat, 0), df = df, p.value = stats::pchisq(max(stat, 0), df, lower.tail = FALSE))
  }
  if (sp$used) attr(out, "suprema") <- TRUE
  if (length(absent)) attr(out, "absent") <- absent
  class(out) <- c("lr_test", "data.frame")
  out
}

## the statistic of the multinomial logit against the constant-share inflated multinomial
## logit (DESIGN.md Phase 4e), shared by lr_test() and boot_lr_test(): stat =
## 2 (l_full - l_restricted), q the restrictions inside the space.
## * A statistic within max(1e-6, 2e-8 |l|) of zero is the optimizer's
##   shortfall at the boundary and counts as 0.
## * With q = 0 a full fit whose share collapsed (phi -> 1) is the multinomial logit
##   itself, and no warning is given when it ends a little below. Its
##   statistic is not otherwise set to 0: a share below the boundary flag's
##   1e-6 can still carry a real likelihood ratio (local review 2026-11-06,
##   finding 8). With q > 0 it is the multinomial logit with the further coefficients, and
##   its ratio is theirs (round 7, A-1).
## * Degenerate: a positive statistic from a fit whose outcome coefficients
##   diverge or that stopped at a saddle (the split coefficients diverge at a
##   collapse and do not count).
## lr_test(), boot_lr_test() and captivity_test(type = "lr") all use this rule.
.pc_captive_rule <- function(stat, full, l_restricted, q) {
  tol <- max(1e-6, 2e-8 * max(1, abs(l_restricted)))
  collapsed <- q == 0L && isTRUE(full$phi_boundary)
  st <- if (stat < tol) 0 else stat
  dv <- setdiff(full$diverging, grep("^infl_", full$diverging, value = TRUE))
  degenerate <- st > 0 && (length(dv) > 0L || isTRUE(full$saddle))
  list(statistic = st, below = stat < -tol && !collapsed, degenerate = degenerate,
       why = if (length(dv)) paste0("diverging: ", paste(utils::head(dv, 3L), collapse = ", ")) else "a saddle")
}

## the p-value of the 50:50 mixture of chi-square(q) and chi-square(q + 1),
## the reference of a likelihood ratio with one parameter on the boundary of
## its space and q restrictions inside it (Self and Liang 1987; DESIGN.md
## Phase 4e); chi-square(0) is the point mass at zero
.pc_chibar_p <- function(stat, q) {
  if (!is.finite(stat)) return(NA_real_)
  p0 <- if (q == 0L) as.numeric(stat <= 0) else stats::pchisq(stat, q, lower.tail = FALSE)
  0.5 * p0 + 0.5 * stats::pchisq(stat, q + 1L, lower.tail = FALSE)
}

#' @export
print.lr_test <- function(x, digits = 4, ...) {
  cat(if (isTRUE(attr(x, "penalized"))) paste0("Penalized likelihood-ratio test (", attr(x, "penalty") %||% "Firth", ")\n") else "Likelihood-ratio test\n")
  cat(sprintf("  LR = %.*f on %d df, p = %s\n", digits, x$statistic, x$df, format.pval(x$p.value, digits = digits)))
  if (!is.null(q <- attr(x, "mixture")))
    cat(sprintf("  The captive share is on the boundary of its space: p-value from the 50:50 mixture of chi-square(%d) and chi-square(%d).\n", q, q + 1L))
  if (isTRUE(attr(x, "suprema"))) cat(" ", sub("^; ", "", .pc_sup_note), "\n", sep = "")
  if (length(a <- attr(x, "absent")))
    cat("  Never chosen, its probability zero at the suprema: ", paste(a, collapse = ", "), "; its coefficients do not enter them and are not counted in df.\n", sep = "")
  invisible(x)
}

#' Compare several choice models side by side
#'
#' @param ... Named or unnamed `"choicemodel"` fits.
#' @return A data frame of class `"compare_models"` with one row per model:
#'   the model type, inflated alternative, cases, parameters, log-likelihood,
#'   AIC, BIC, McFadden R-squared, and the `converged`, `boundary`, and
#'   `ill_conditioned` flags and `n_modes`; with a penalized fit among them a
#'   `penalized` column, and the print notes that their log-likelihoods and
#'   criteria are at the penalized estimates, not maximum-likelihood ones.
#' @examples
#' set.seed(2)
#' d <- rmnl(1000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
#'           delta = c(0.5, 1.2), inflate = "a")
#' compare_models(mnl = mnl(y ~ x1 + x2, data = d),
#'                imnl = imnl(y ~ x1 + x2 | z1, data = d, inflate = "a"))
#' @seealso [vuong()], [lr_test()], [fit_stats()]
#' @family model comparison
#' @export
compare_models <- function(...) {
  ms <- list(...)
  if (!length(ms)) stop("Give at least one fitted model.", call. = FALSE)
  nm <- names(ms); dn <- vapply(substitute(list(...))[-1L], function(e) deparse(e)[1L], character(1))
  if (is.null(nm)) nm <- dn else nm[nm == ""] <- dn[nm == ""]
  for (i in seq_along(ms)) if (!inherits(ms[[i]], "choicemodel")) stop("Argument ", i, " is not a fitted polychoice model.", call. = FALSE)
  if (length(ms) > 1L) for (i in seq_along(ms)[-1L])
    if (ms[[i]]$n != ms[[1L]]$n || !identical(ms[[i]]$case_id, ms[[1L]]$case_id) || !isTRUE(all.equal(ms[[i]]$weights, ms[[1L]]$weights)))
      warning("Model ", nm[i], " is not fit to the same cases (or weights) as ", nm[1L], "; the information criteria are not comparable.", call. = FALSE)
  rows <- lapply(seq_along(ms), function(i) {
    m <- ms[[i]]
    fs <- fit_stats(m)
    data.frame(model = nm[i], type = .pc_title(m), inflate = if (m$inflated) paste(m$inflate, collapse = ", ") else "", n = m$n, df = m$df,
               loglik = m$loglik, AIC = fs$AIC, BIC = fs$BIC, mcfadden = fs$mcfadden,
               converged = isTRUE(m$converged), boundary = isTRUE(m$boundary), ill_conditioned = isTRUE(m$ill_conditioned),
               n_modes = if (is.null(m$n_modes)) 1L else as.integer(m$n_modes), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  ## penalized fits: their log-likelihoods and criteria are not maximum-likelihood ones
  if (any(vapply(ms, .pc_penalized, logical(1)))) out$penalized <- vapply(ms, .pc_penalized, logical(1))
  class(out) <- c("compare_models", "data.frame")
  out
}

#' @export
print.compare_models <- function(x, digits = 3, ...) {
  y <- as.data.frame(x)
  y$loglik <- round(y$loglik, digits); y$AIC <- round(y$AIC, digits); y$BIC <- round(y$BIC, digits); y$mcfadden <- round(y$mcfadden, digits)
  print.data.frame(y, row.names = FALSE)
  if (!is.null(y$penalized) && any(y$penalized))
    cat("Penalized fits (Firth's penalty or a prior): log-likelihood at the penalized estimates; their AIC, BIC and R2 are not maximum-likelihood ones.\n")
  invisible(x)
}
