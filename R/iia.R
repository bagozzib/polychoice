## ---------------------------------------------------------------------------
## Tests of the independence of irrelevant alternatives (IIA) / misspecification
## of the multinomial logit: Hausman-McFadden (1984), Small-Hsiao (1985), and
## the composite-likelihood Hausman and GMM overidentification tests of Fok and
## Paap (2025). Every test takes a fitted mnl() object, so a test is always on
## the model the user reports. Definitions are fixed in DESIGN.md, section 5.
## ---------------------------------------------------------------------------

.pc_check_mnl <- function(object, what) {
  if (!inherits(object, "choicemodel")) stop("'object' must be a fitted polychoice model.", call. = FALSE)
  if (object$inflated) stop(what, " is a test of the multinomial logit; it does not apply to an inflated fit. Test the mnl() fit of the same specification.", call. = FALSE)
  if (!identical(object$family, "mnl")) stop(what, " is a test of the multinomial logit; it does not apply to a ", object$family, " fit. Test the mnl() fit of the same specification.", call. = FALSE)
  if (!is.null(object$fe)) stop(what, " does not apply to a fixed-effects fit: the unit constants would enter the contrast and the half-sample or resampled fits. Fit the specification without 'fe' or with mundlak(); a cluster = unit then makes the covariance-based tests cluster-robust (the Small-Hsiao and likelihood-ratio statistics take the cases as independent).", call. = FALSE)
  if (object$design$J < 3L) stop(what, " needs at least three alternatives.", call. = FALSE)
  invisible(TRUE)
}

## the design re-expressed with another base alternative (same data), with the
## parameter names / blocks of that parameterization
.pc_set_base <- function(D, newbase) {
  D$base <- newbase
  nonbase <- setdiff(seq_len(D$J), newbase)
  beta_names <- if (D$p) as.vector(vapply(nonbase, function(j) paste0(D$levels[j], ":", colnames(D$X)), character(D$p))) else character(0)
  D$par_names <- c(beta_names, colnames(D$W), if (D$r) paste0("infl_", colnames(D$Z)))
  D$idx$beta <- if (D$p) matrix(seq_len(D$p * (D$J - 1L)), nrow = D$p, dimnames = list(colnames(D$X), D$levels[nonbase])) else matrix(integer(0), 0L, D$J - 1L)
  D
}

## the exact linear map from the parameters under base b to those under base c:
## beta_j^(c) = beta_j - beta_c (j not in {b, c}), beta_b^(c) = -beta_c, gamma unchanged
.pc_rebase_matrix <- function(D, newbase) {
  K <- D$p * (D$J - 1L) + D$q
  A <- matrix(0, K, K)
  old_nb <- setdiff(seq_len(D$J), D$base); new_nb <- setdiff(seq_len(D$J), newbase)
  blk <- function(nb, j) { m <- match(j, nb); (m - 1L) * D$p + seq_len(D$p) }
  for (j in new_nb) {
    rows <- blk(new_nb, j)
    if (j != D$base) A[rows, blk(old_nb, j)] <- diag(D$p)
    A[rows, blk(old_nb, newbase)] <- A[rows, blk(old_nb, newbase)] - diag(D$p)
  }
  if (D$q) { g <- D$p * (D$J - 1L) + seq_len(D$q); A[g, g] <- diag(D$q) }
  A
}

## full-model estimate, covariance and scores expressed under base `newbase`
.pc_rebase_fit <- function(object, newbase, V, S) {
  D <- object$design
  if (newbase == D$base) return(list(theta = object$theta, vcov = V, scores = S, D = D))
  A <- .pc_rebase_matrix(D, newbase)
  Dn <- .pc_set_base(D, newbase)
  theta <- as.vector(A %*% object$theta); names(theta) <- Dn$par_names
  Vn <- A %*% V %*% t(A); dimnames(Vn) <- list(Dn$par_names, Dn$par_names)
  ## scores transform with the inverse transpose: s' = A^{-T} s (theta' = A theta)
  Sn <- if (is.null(S)) NULL else S %*% solve(A)
  list(theta = theta, vcov = Vn, scores = Sn, D = Dn)
}

## resolve `drop` to alternative indices
.pc_resolve_drop <- function(object, drop) {
  lev <- object$levels
  if (is.null(drop)) return(seq_along(lev))
  idx <- if (is.numeric(drop)) as.integer(drop) else match(as.character(drop), lev)
  if (anyNA(idx) || any(idx < 1L | idx > length(lev))) stop("'drop' must name alternatives among: ", paste(lev, collapse = ", "), ".", call. = FALSE)
  unique(idx)
}

## a quadratic form with the usual safeguards: returns the statistic, whether
## the matrix was positive semi-definite, and whether a generalized inverse was used
.pc_quad_form <- function(delta, Vd, ref = NULL) {
  ## with ref, the diagonal of the efficient estimator's covariance, the
  ## contrast is first scaled by the estimators' own standard deviations (the
  ## statistic is unchanged), so that a direction in which the contrast has no
  ## variance relative to the estimators can be recognized: when the two
  ## estimators are the same function of the data (a saturated model, every
  ## estimator reproducing the cell shares) the contrast and its covariance are
  ## both zero up to rounding, and judged against its own largest eigenvalue
  ## that rounding would pass for information. Such directions are dropped
  ## (degrees of freedom reduced); with none left the test is undefined (NA)
  if (!is.null(ref) && length(ref) == length(delta) && all(is.finite(ref) & ref > 0)) {
    sr <- sqrt(ref); delta <- delta / sr; Vd <- Vd / tcrossprod(sr)
    es <- eigen((Vd + t(Vd)) / 2, symmetric = TRUE)
    keep <- abs(es$values) > sqrt(.Machine$double.eps)
    psd0 <- all(es$values >= -sqrt(.Machine$double.eps))
    if (!any(keep)) return(list(stat = NA_real_, psd = psd0, ginv = TRUE, df = 0L, degenerate = TRUE))
    if (!all(keep)) {
      z <- crossprod(es$vectors[, keep, drop = FALSE], delta)
      return(list(stat = sum(z^2 / es$values[keep]), psd = psd0, ginv = TRUE, df = sum(keep), degenerate = FALSE))
    }
  }
  ev <- eigen((Vd + t(Vd)) / 2, symmetric = TRUE, only.values = TRUE)$values
  psd <- all(ev >= -sqrt(.Machine$double.eps) * max(abs(ev), 1))
  inv <- tryCatch(solve(Vd), error = function(e) { .pc_pass_limit(e); NULL })
  ginv_used <- is.null(inv)
  if (ginv_used) inv <- MASS::ginv(Vd)
  ## with a generalized inverse the statistic lives in the column space of Vd:
  ## its degrees of freedom are the rank (a cluster sandwich with fewer
  ## clusters than parameters has rank below the dimension)
  rk <- if (ginv_used) sum(abs(ev) > sqrt(.Machine$double.eps) * max(abs(ev), 1e-300)) else length(delta)
  list(stat = as.numeric(t(delta) %*% inv %*% delta), psd = psd, ginv = ginv_used, df = rk, degenerate = FALSE)
}
## the IIA contrasts compare estimators of the same coefficients; a coefficient
## that diverges (quasi-complete separation) has no estimate, and contrasts
## through it can reject at any size (Kardasheva 2013 in the replication
## stress test: Q1 776 on 19 df with the separating terms, 123 on 54 without)
.pc_separation_warning <- function(object, what,
                                   why = "so the contrasts compare estimates that do not exist and can reject spuriously; drop or merge the separating terms before testing IIA") {
  dv <- tryCatch(.pc_diverging(object$theta, object$design, object$scale), error = function(e) { .pc_pass_limit(e); logical(0) })
  nm <- object$design$par_names[dv %in% TRUE]
  if (length(nm)) warning(what, ": the fit has diverging coefficients (", paste(utils::head(nm, 6), collapse = ", "), if (length(nm) > 6) ", ..." else "",
                          "; quasi-complete separation), ", why, ".", call. = FALSE)
}
## the score tests' form of the warning
.pc_score_sep_why <- "so the score statistic is evaluated where the optimizer stopped, at estimates that do not exist, and its value depends on that point (fitted probabilities near zero weigh heavily); drop or merge the separating terms first"
.pc_degenerate_warning <- function(what) {
  if (length(what)) warning(what, ": the contrast has no variance relative to the estimators (the two estimators coincide, as in a saturated model ",
                            "where every one reproduces the cell shares), so there is nothing to test and the statistic is NA.", call. = FALSE)
}

#' Hausman-McFadden test of IIA
#'
#' For each dropped alternative, compares the multinomial logit estimated on the
#' cases that did not choose it (restricted choice set) with the full-sample
#' estimate of the same parameters: under IIA both are consistent and the
#' restricted one is less efficient, so
#' \eqn{(\hat\beta_r - \hat\beta_f)'(V_r - V_f)^{-1}(\hat\beta_r - \hat\beta_f)}
#' is asymptotically chi-square with as many degrees of freedom as restricted
#' parameters (Hausman and McFadden 1984).
#'
#' @details
#' Dropping the base alternative is handled exactly: the full model is
#' re-expressed with the first remaining alternative as the base through the
#' linear map \eqn{\beta_j^{(c)} = \beta_j - \beta_c}, and the restricted model
#' is estimated with that base, so the compared parameters coincide (the
#' misalignment that produces spurious rejections for the base-drop row in some
#' implementations does not arise).
#'
#' `vcov = "difference"` is the classic form, with both covariances the inverse
#' observed information (the form the Hausman lemma requires; the fit's own
#' robust or cluster covariance is not used here). The difference of covariance
#' matrices need not be positive semi-definite in finite samples; the result
#' then reports `psd = FALSE`, a negative statistic is shown as computed, and
#' its p-value is 1 (Hausman and McFadden 1984, p. 1226, and Long and Freese
#' 2014 read a negative statistic as no evidence against IIA). A positive
#' statistic with `psd = FALSE` is unreliable: in the package's size study the
#' difference form rejected well above the nominal level. It stays the default
#' because it is the form the literature and Stata's `hausman` report, so that
#' published values can be reproduced; the joint form is the one to rely on.
#' `vcov = "joint"` instead estimates the covariance of the difference from the
#' stacked per-case scores of the two fits (the seemingly-unrelated-estimation
#' device of Weesie 1999; Stata's `mlogtest, suest`): \eqn{V_r + V_f - C - C'}
#' with the sandwich covariances and cross-covariance, positive semi-definite
#' by construction and valid with heteroskedasticity and, when the fit has a
#' `cluster`, clustering.
#'
#' @param object A fitted [mnl()] model with at least three alternatives.
#'   A fit with diverging coefficients (separation) draws a warning: the
#'   contrasts then compare estimates that do not exist and can reject
#'   spuriously.
#' @param drop Alternatives to drop, one at a time (labels or indices); default
#'   all of them.
#' @param vcov `"difference"` (default) or `"joint"`; see Details.
#' @return A data frame of class `"iia_test"` with one row per dropped
#'   alternative: `alternative`, `statistic`, `df`, `p.value`, `psd` (was the
#'   covariance of the difference positive semi-definite), and `n_restricted`.
#' @references
#' Hausman, J.A. and McFadden, D. (1984). Specification tests for the
#' multinomial logit model. \emph{Econometrica}, 52, 1219-1240.
#' \doi{10.2307/1910997}
#'
#' Weesie, J. (1999). Seemingly unrelated estimation and the cluster-adjusted
#' sandwich estimator. \emph{Stata Technical Bulletin}, 52, 34-47.
#'
#' Long, J.S. and Freese, J. (2014). \emph{Regression Models for Categorical
#' Dependent Variables Using Stata}, 3rd ed. Stata Press.
#' @examples
#' set.seed(1)
#' d <- rmnl(1200, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
#' m <- mnl(y ~ x1 + x2, data = d)
#' hausman_mcfadden(m)
#' hausman_mcfadden(m, drop = "b", vcov = "joint")
#' @seealso [small_hsiao()], [fp_tests()], [iia_tests()]
#' @family IIA tests
#' @export
hausman_mcfadden <- function(object, drop = NULL, vcov = c("difference", "joint")) {
  .pc_check_mnl(object, "The Hausman-McFadden test")
  .pc_no_firth(object, "The Hausman-McFadden test")
  .pc_separation_warning(object, "The Hausman-McFadden test")
  vcov <- match.arg(vcov)
  D <- object$design
  fn <- .pc_fns(D)
  theta_f <- object$theta
  ## information-based covariance of the full fit (recomputed: the fit's own may be robust)
  Vf_info <- .pc_strip_vcov(.pc_vcov(theta_f, fn, "analytic", NULL, object$scale, D$w))
  Sf <- fn$scores(theta_f)
  drops <- .pc_resolve_drop(object, drop)
  rows <- lapply(drops, function(d) {
    newbase <- if (d == D$base) setdiff(seq_len(D$J), d)[1L] else D$base
    Dr <- .pc_drop_alternative(D, d, newbase)
    fit_r <- .pc_refit(Dr, vcov = TRUE)
    theta_r <- fit_r$coefficients
    rb <- .pc_rebase_fit(object, newbase, Vf_info, Sf)
    sel <- match(names(theta_r), names(rb$theta))
    if (anyNA(sel)) stop("internal: restricted parameters do not match the full model.")
    delta <- theta_r - rb$theta[sel]
    if (vcov == "difference") {
      Vd <- fit_r$vcov - rb$vcov[sel, sel, drop = FALSE]
    } else {
      ## suest: joint sandwich of the two M-estimators over the full set of cases
      Sr_full <- matrix(0, D$n, length(theta_r))
      kept <- which(D$y != d)
      Sr_full[kept, ] <- fit_r$fn$scores(theta_r)
      Hr_inv <- .pc_strip_vcov(.pc_vcov(theta_r, fit_r$fn, "analytic", NULL, fit_r$scale, Dr$w))         # bread of the restricted fit
      Hf_inv <- rb$vcov                                                                  # bread of the (re-based) full fit
      ## Hf_inv[sel, sel] is not the bread of the sub-vector; use the full bread and select after
      wS_r <- Sr_full * D$w; wS_f <- rb$scores * D$w
      if (!is.null(D$cluster)) { wS_r <- as.matrix(rowsum(wS_r, D$cluster)); wS_f <- as.matrix(rowsum(wS_f, D$cluster)); Sr_u <- wS_r; Sf_u <- wS_f }
      else { Sr_u <- Sr_full * sqrt(D$w); Sf_u <- rb$scores * sqrt(D$w) }
      Vr <- Hr_inv %*% crossprod(Sr_u) %*% Hr_inv
      Vf <- (Hf_inv %*% crossprod(Sf_u) %*% Hf_inv)[sel, sel, drop = FALSE]
      C <- (Hr_inv %*% crossprod(Sr_u, Sf_u) %*% Hf_inv)[, sel, drop = FALSE]
      Vd <- Vr + Vf - C - t(C)
    }
    qf <- .pc_quad_form(delta, Vd, ref = diag(rb$vcov)[sel])
    data.frame(alternative = D$levels[d], statistic = qf$stat, df = qf$df,
               p.value = stats::pchisq(max(qf$stat, 0), qf$df, lower.tail = FALSE),
               psd = qf$psd, n_restricted = Dr$n, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  dg <- is.na(out$statistic) & out$df == 0
  .pc_degenerate_warning(if (any(dg)) paste0("Hausman-McFadden, dropping ", paste(out$alternative[dg], collapse = ", ")))
  attr(out, "method") <- paste0("Hausman-McFadden test of IIA (", if (vcov == "difference") "variance difference" else "joint (suest) covariance", ")")
  attr(out, "note") <- if (vcov == "difference" && any(!out$psd)) "psd = FALSE: the covariance difference is not positive semi-definite. A negative statistic is no evidence against IIA (p-value 1); a positive one is unreliable (the difference form over-rejects); use vcov = \"joint\"." else NULL
  class(out) <- c("iia_test", "data.frame")
  out
}

## a stratified half-split of the cases by chosen alternative
.pc_half_split <- function(y, J) {
  A <- integer(0)
  for (j in seq_len(J)) {
    idx <- which(y == j); nj <- length(idx)
    if (nj >= 2L) A <- c(A, sample(idx, floor(nj / 2)))
  }
  list(A = sort(A), B = sort(setdiff(seq_along(y), A)))
}

#' Small-Hsiao test of IIA
#'
#' The likelihood-ratio form of the IIA test with a random half-split of the
#' cases (Small and Hsiao 1985). The full model is estimated on both halves,
#' \eqn{\hat\theta_1} and \eqn{\hat\theta_2}, and combined as
#' \eqn{\hat\theta_{12} = \sqrt{1/2}\,\hat\theta_1 + (1 - \sqrt{1/2})\,\hat\theta_2};
#' the restricted model (alternative \eqn{d} dropped) is estimated on the second
#' half, \eqn{\hat\theta_r}, and the statistic is
#' \eqn{-2[\ell_2^r(\hat\theta_{12}) - \ell_2^r(\hat\theta_r)]} with the
#' second-half restricted log-likelihood \eqn{\ell_2^r}, chi-square on as many
#' degrees of freedom as restricted parameters under IIA. The weights
#' \eqn{\sqrt{1/2}} and \eqn{1 - \sqrt{1/2}} are what make the statistic
#' chi-square; the variant that drops them and also estimates the restricted
#' model on both halves (as in some earlier implementations) is twice a
#' chi-square and rejects far too often.
#'
#' The split is stratified by chosen alternative. Results can depend heavily on
#' the random split: one split may give a large statistic that the next does
#' not. The test is one split chosen before the data are seen: fix `seed` in
#' advance (`small_hsiao(m, seed = )`) and report that split. In the package's
#' size study one split held its level (rejection 0.040 to 0.047 at 0.05),
#' while the median over ten splits rejected in none of 300 samples, so the
#' median is not a test at its level. `splits > 1` repeats the statistic on
#' independent splits and reports every one together with the median and the
#' share of splits rejecting at the 5 percent level, a description of how much
#' the statistic varies (Long and Freese 2014 recommend inspecting this
#' variability; [iia_tests()] reports the median over ten splits);
#' [boot_pvalue()] can calibrate the median. When the pre-specified split's
#' half-sample fits are separated (`separated`), its statistic reflects the
#' separation rather than IIA: report the row as uninformative, or calibrate
#' it with `boot_pvalue(m, small_hsiao, seed = )` and the same seed. A seed
#' chosen after seeing a result is not a test. Dropping the base alternative is
#' handled by re-basing, as in [hausman_mcfadden()]. The test is not
#' frequency-weight exact: the halves are balanced in cases, not in weight. The
#' cases are taken as independent: a fit's `cluster` does not enter the
#' statistic. Not for fixed-effects fits.
#'
#' @param object A fitted [mnl()] model with at least three alternatives.
#'   A fit with diverging coefficients (separation) draws a warning: the
#'   contrasts then compare estimates that do not exist and can reject
#'   spuriously.
#' @param drop Alternatives to drop, one at a time; default all.
#' @param splits Number of independent half-splits (default 1, the test;
#'   more describe how the statistic varies across splits, see Details).
#' @param seed Optional seed for the splits; the caller's random-number stream
#'   is restored afterwards.
#' @return A data frame of class `"iia_test"` with one row per split and
#'   dropped alternative: `split`, `alternative`, `statistic`, `df`, `p.value`,
#'   and `separated` (a fit behind that row has diverging coefficients in an
#'   equation the statistic uses: the quasi-separation of a half-sample that a
#'   rare alternative with binary covariates invites, after which the statistic
#'   may reflect the separation rather than IIA; a divergence confined to the
#'   dropped alternative's own equation does not mark the row);
#'   with `splits > 1` the attribute `"summary"` holds per alternative the
#'   median statistic over the splits, `median_p` (the chi-square p-value of
#'   that median, as [iia_tests()] reports it; descriptive, since the median of
#'   split statistics is not chi-square, and with an even number of splits not
#'   the median of the split p-values), the share of splits rejecting at
#'   5 percent, and the share of separated splits. [boot_pvalue()] calibrates the median (see the IIA vignette).
#' @references
#' Long, J.S. and Freese, J. (2014). \emph{Regression Models for Categorical
#' Dependent Variables Using Stata}, 3rd ed. College Station, TX: Stata Press.
#'
#' Small, K.A. and Hsiao, C. (1985). Multinomial logit specification tests.
#' \emph{International Economic Review}, 26, 619-627. \doi{10.2307/2526707}
#' @examples
#' set.seed(1)
#' d <- rmnl(1200, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
#' m <- mnl(y ~ x1 + x2, data = d)
#' small_hsiao(m, seed = 1)
#' small_hsiao(m, drop = "c", splits = 5, seed = 1)
#' @seealso [hausman_mcfadden()], [fp_tests()], [iia_tests()]
#' @family IIA tests
#' @export
small_hsiao <- function(object, drop = NULL, splits = 1L, seed = NULL) {
  .pc_check_mnl(object, "The Small-Hsiao test")
  .pc_no_firth(object, "The Small-Hsiao test")
  .pc_separation_warning(object, "The Small-Hsiao test")
  splits <- .pc_check_whole(splits, "splits")
  if (!is.null(seed)) {
    ## seed the splits without disturbing the caller's random-number stream
    if (!exists(".Random.seed", envir = globalenv(), inherits = FALSE)) stats::runif(1)
    old <- get(".Random.seed", envir = globalenv()); on.exit(assign(".Random.seed", old, envir = globalenv()))
    set.seed(seed)
  }
  D <- object$design
  if (any(D$w != 1)) warning("The Small-Hsiao halves are balanced in cases, not in weight; with frequency weights the sqrt(1/2) combination is only approximate.", call. = FALSE)
  drops <- .pc_resolve_drop(object, drop)
  wgt <- sqrt(0.5)
  rows <- list()
  for (s in seq_len(splits)) {
    sp <- .pc_half_split(D$y, D$J)
    DA <- .pc_subset_design(D, sp$A); DB <- .pc_subset_design(D, sp$B)
    if (.pc_own_coef(D) && (any(tabulate(DA$y, nbins = D$J) == 0L) || any(tabulate(DB$y, nbins = D$J) == 0L)))
      stop("A half-sample has no case choosing some alternative; the Small-Hsiao test needs every alternative in both halves.", call. = FALSE)
    r1 <- .pc_refit(DA, start = object$theta); r2 <- .pc_refit(DB, start = object$theta)
    th1 <- r1$coefficients; th2 <- r2$coefficients
    th12 <- wgt * th1 + (1 - wgt) * th2
    ## a half-sample fit with diverging coefficients (quasi-separation, likely
    ## with a rare alternative and binary covariates) makes the split's
    ## statistic reflect the separation rather than IIA
    ## the alternatives whose equations diverge in either half (0 for a generic
    ## coefficient, which enters every equation)
    div_alts <- function(dv, Dh) {
      nb <- setdiff(seq_len(Dh$J), Dh$base); a <- integer(0)
      if (Dh$p && any(dv[as.vector(Dh$idx$beta)]))
        a <- nb[unique(which(matrix(dv[as.vector(Dh$idx$beta)], nrow = Dh$p), arr.ind = TRUE)[, 2L])]
      if (length(Dh$idx$gamma) && any(dv[Dh$idx$gamma])) a <- c(a, 0L)
      unique(a)
    }
    sep_alts <- unique(c(div_alts(.pc_diverging(th1, DA, r1$scale), DA),
                         div_alts(.pc_diverging(th2, DB, r2$scale), DB)))
    for (d in drops) {
      newbase <- if (d == D$base) setdiff(seq_len(D$J), d)[1L] else D$base
      DrB <- .pc_drop_alternative(DB, d, newbase)
      fit_r <- .pc_refit(DrB)
      ## theta_12 in the restricted parameterization (re-based when d is the base)
      th12_c <- if (newbase == D$base) th12 else { A <- .pc_rebase_matrix(D, newbase); v <- as.vector(A %*% th12); names(v) <- .pc_set_base(D, newbase)$par_names; v }
      sel <- match(DrB$par_names, names(th12_c))
      if (anyNA(sel)) stop("internal: restricted parameters do not match the full model.")
      ll12 <- -fit_r$fn$nll(unname(th12_c[sel]))
      stat <- -2 * (ll12 - fit_r$loglik)
      ## a row is affected unless every diverging equation is the dropped
      ## alternative's own (dropping the base re-bases every other equation on a
      ## new base, so then any divergence matters), or the restricted fit diverges
      half_sep <- length(sep_alts) > 0L && !(d != D$base && all(sep_alts == d))
      sep <- half_sep || any(.pc_diverging(fit_r$coefficients, DrB, fit_r$scale))
      rows[[length(rows) + 1L]] <- data.frame(split = s, alternative = D$levels[d], statistic = stat, df = DrB$K,
                                              p.value = stats::pchisq(max(stat, 0), DrB$K, lower.tail = FALSE), separated = sep, stringsAsFactors = FALSE)
    }
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, "method") <- "Small-Hsiao test of IIA (stratified half-split, sqrt(1/2) weighting)"
  if (splits > 1L) {
    sm <- do.call(rbind, lapply(split(out, out$alternative), function(g)
      data.frame(alternative = g$alternative[1L], median_statistic = stats::median(g$statistic), df = g$df[1L],
                 median_p = stats::pchisq(stats::median(g$statistic), g$df[1L], lower.tail = FALSE), share_reject_05 = mean(g$p.value < 0.05),
                 share_separated = mean(g$separated), stringsAsFactors = FALSE)))
    rownames(sm) <- NULL
    attr(out, "summary") <- sm[match(D$levels[drops], sm$alternative), ]
  }
  class(out) <- c("iia_test", "data.frame")
  out
}

## the pairwise machinery of Fok and Paap: for each pair (j, m) the cases with
## both alternatives available and y in {j, m}, the difference design
## Delta = D_j - D_m (so V_j - V_m = Delta' theta + o_j - o_m), the binary
## outcome 1{y = j}, and the instruments zeta = (x, w_j - w_m)
.pc_pairs <- function(D, pairs) {
  J <- D$J; n <- D$n; p <- D$p; q <- D$q; K <- p * (J - 1L) + q
  nonbase <- setdiff(seq_len(J), D$base)
  blk <- function(j) if (j == D$base || p == 0L) integer(0) else (match(j, nonbase) - 1L) * p + seq_len(p)
  pl <- if (pairs == "base") lapply(setdiff(seq_len(J), D$base), function(j) c(j, D$base)) else utils::combn(J, 2L, simplify = FALSE)
  lapply(pl, function(pr) {
    j <- pr[1L]; m <- pr[2L]
    cases <- which(D$avail[, j] == 1L & D$avail[, m] == 1L & (D$y == j | D$y == m))
    nc <- length(cases)
    Delta <- matrix(0, nc, K)
    if (p) { if (j != D$base) Delta[, blk(j)] <- D$X[cases, , drop = FALSE]; if (m != D$base) Delta[, blk(m)] <- Delta[, blk(m)] - D$X[cases, , drop = FALSE] }
    Wd <- if (q) D$W[(cases - 1L) * J + j, , drop = FALSE] - D$W[(cases - 1L) * J + m, , drop = FALSE] else matrix(0, nc, 0L)
    if (q) Delta[, p * (J - 1L) + seq_len(q)] <- Wd
    list(j = j, m = m, cases = cases, Delta = Delta, odiff = D$off[cases, j] - D$off[cases, m],
         yj = as.numeric(D$y[cases] == j), zeta = cbind(D$X[cases, , drop = FALSE], Wd), w = D$w[cases])
  })
}

## composite-likelihood estimator over the pairs by Newton (the pairwise logit
## log-likelihood is concave); returns the estimate, the sensitivity H (sum of
## weighted Lambda(1-Lambda) Delta Delta'), and the per-case summed scores
.pc_cl_fit <- function(P, theta0, n, K) {
  theta <- theta0
  for (it in seq_len(50L)) {
    g <- numeric(K); H <- matrix(0, K, K)
    for (pr in P) {
      eta <- as.vector(pr$Delta %*% theta) + pr$odiff; lam <- stats::plogis(eta)
      g <- g + colSums(pr$w * (pr$yj - lam) * pr$Delta)
      H <- H + crossprod(pr$Delta * sqrt(pr$w * lam * (1 - lam)))
    }
    step <- tryCatch(solve(H, g), error = function(e) { .pc_pass_limit(e); MASS::ginv(H) %*% g })
    theta <- theta + as.vector(step)
    if (max(abs(step)) < 1e-10) break
  }
  S <- matrix(0, n, K); H <- matrix(0, K, K)
  for (pr in P) {
    eta <- as.vector(pr$Delta %*% theta) + pr$odiff; lam <- stats::plogis(eta)
    S[pr$cases, ] <- S[pr$cases, ] + (pr$yj - lam) * pr$Delta
    H <- H + crossprod(pr$Delta * sqrt(pr$w * lam * (1 - lam)))
  }
  list(theta = theta, H = H, scores = S, converged = it < 50L)
}

## per-case moment matrix (n x M) and the Jacobian of the mean moment (M x K)
.pc_gmm_moments <- function(P, theta, n, K, nz, w) {
  ## G holds the UNWEIGHTED per-case moments; the weights enter once, in the
  ## mean sum(w g) / N_w and in the covariance sum(w g g') / N_w (the duplicated-
  ## rows forms), never squared
  M <- length(P) * nz; nw <- sum(w)
  G <- matrix(0, n, M); Dm <- matrix(0, M, K)
  for (b in seq_along(P)) {
    pr <- P[[b]]
    eta <- as.vector(pr$Delta %*% theta) + pr$odiff; lam <- stats::plogis(eta)
    cols <- (b - 1L) * nz + seq_len(nz)
    G[pr$cases, cols] <- (pr$yj - lam) * pr$zeta
    Dm[cols, ] <- -crossprod(pr$zeta, pr$w * lam * (1 - lam) * pr$Delta) / nw
  }
  list(G = G, D = Dm, gbar = colSums(G * w) / nw, Omega = crossprod(G * sqrt(w)) / nw, nw = nw)
}

#' Fok-Paap tests of multinomial logit misspecification
#'
#' The two tests of Fok and Paap (2025), built on the property that under the
#' multinomial logit the choice between any two alternatives \eqn{j, m}, given
#' that one of them is chosen, is a binary logit in the utility difference
#' \eqn{V_j - V_m}. Pairs are the \eqn{J - 1} base pairs or all
#' \eqn{J(J-1)/2} pairs.
#'
#' **Composite-likelihood Hausman test (Q1).** The pairwise composite
#' likelihood estimator \eqn{\hat\theta_{CL}} maximizes the sum over pairs and
#' cases of the binary logit log-likelihoods; its covariance is the Godambe
#' sandwich. Under correct specification the MLE and the CL estimator are both
#' consistent and only the MLE is efficient, so
#' \eqn{Q_1 = (\hat\theta_{ML} - \hat\theta_{CL})'(V_{CL} - V_{ML})^{-1}(\hat\theta_{ML} - \hat\theta_{CL})}
#' is chi-square on \eqn{\dim\theta} degrees of freedom (F&P eq. 18).
#' Under the multinomial logit the composite-likelihood estimator is nearly
#' efficient (its variance is within a few percent of the MLE's in most
#' directions), so \eqn{V_{CL} - V_{ML}} is a small difference of two
#' separately estimated matrices, and the sampling noise of the sandwich makes
#' it indefinite in nearly every sample: in the package's simulations
#' (`verification/` in the development repository) it was never positive
#' semi-definite, 44 percent of the statistics were negative, and the
#' rejection rates were 0.033, 0.042, and 0.049 at the 1, 5, and 10 percent
#' levels. The default `vcov = "joint"` therefore estimates the covariance of
#' \eqn{\hat\theta_{ML} - \hat\theta_{CL}} directly, from the stacked
#' influence functions \eqn{V_{ML} s^{ML}_i - H_{CL}^{-1} s^{CL}_i} of the two
#' estimators (summed by cluster when the fit has a `cluster`), as the joint
#' form of [hausman_mcfadden()] does; it is positive semi-definite by
#' construction and its rejection rates in the same simulations were 0.015,
#' 0.062, and 0.112. `vcov = "difference"` is the formula of F&P eq. 18, for
#' comparison with published values; `psd` reports whether its covariance
#' difference was positive semi-definite, and the p-value refers
#' `max(statistic, 0)` to the chi-square (a negative statistic is no evidence
#' against the model).
#'
#' **GMM overidentification test (Q2).** With moment conditions
#' \eqn{E[(1\{y = j\} - \Lambda(V_j - V_m))\, 1\{y \in \{j, m\}\}\, \zeta] = 0}
#' for every pair, instruments \eqn{\zeta = (x, w_j - w_m)}, the parameters are
#' overidentified and Hansen's J statistic at the efficient GMM estimate is
#' chi-square on (number of moments \eqn{-} number of parameters) degrees of
#' freedom (F&P eq. 19). `gmm = "two-step"` (default) is Hansen's two-step
#' estimator: the weight matrix is the inverse empirical covariance of the
#' moments at the MLE (F&P eq. 15), the parameters are re-estimated under it,
#' and the statistic is evaluated at that estimate with its own weight matrix.
#' `gmm = "iterated"` repeats the update until the estimate stops moving, as
#' Fok and Paap do; with many moments relative to the sample the iteration can
#' fail to settle (a warning is issued and the last iterate is used), which is
#' why it is not the default. `gmm = "projected"` evaluates the moments at the
#' MLE and projects out the estimated directions (Newey's conditional-moment /
#' C form), which has the same asymptotic distribution without a GMM
#' optimization. All three have the same degrees of freedom. With case-specific
#' covariates only, the overidentifying restrictions carry information only
#' through the covariates (Fok and Paap 2025, Section 2.2): with few or weak
#' covariates the moment covariance is near singular and the test has little
#' power; its reciprocal condition number is reported.
#'
#' @param object A fitted [mnl()] model with at least three alternatives.
#'   A fit with diverging coefficients (separation) draws a warning: the
#'   contrasts then compare estimates that do not exist and can reject
#'   spuriously.
#' @param pairs `"all"` (default) or `"base"`.
#' @param tests Which tests to compute: `"cl"`, `"gmm"`, or both (default).
#' @param gmm `"two-step"` (default), `"iterated"`, or `"projected"`; see Details.
#' @param vcov The covariance in \eqn{Q_1}: `"joint"` (default) or
#'   `"difference"`; see Details.
#' @param maxit Maximum number of weight-matrix updates for `gmm = "iterated"`.
#' @return A data frame of class `"iia_test"` with one row per test:
#'   `test`, `statistic`, `df`, `p.value`, `psd`, and `pairs`. The attribute
#'   `"details"` holds the CL and GMM estimates and the reciprocal condition
#'   number of the moment covariance.
#' @references
#' Fok, D. and Paap, R. (2025). New misspecification tests for multinomial
#' logit models. \emph{Journal of Choice Modelling}, 54, 100531.
#' \doi{10.1016/j.jocm.2024.100531}
#'
#' Varin, C., Reid, N. and Firth, D. (2011). An overview of composite
#' likelihood methods. \emph{Statistica Sinica}, 21, 5-42.
#'
#' Hansen, L.P. (1982). Large sample properties of generalized method of
#' moments estimators. \emph{Econometrica}, 50, 1029-1054.
#' \doi{10.2307/1912775}
#' @examples
#' set.seed(1)
#' d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
#' m <- mnl(y ~ x1 + x2, data = d)
#' fp_tests(m)
#' fp_tests(m, gmm = "projected")
#' fp_tests(m, gmm = "iterated")
#' @seealso [hausman_mcfadden()], [small_hsiao()], [iia_tests()]
#' @family IIA tests
#' @export
fp_tests <- function(object, pairs = c("all", "base"), tests = c("cl", "gmm"), gmm = c("two-step", "iterated", "projected"),
                     vcov = c("joint", "difference"), maxit = 50L) {
  .pc_check_mnl(object, "The Fok-Paap tests")
  .pc_no_firth(object, "The Fok-Paap tests")
  .pc_separation_warning(object, "The Fok-Paap tests")
  pairs <- match.arg(pairs); gmm <- match.arg(gmm); vcov <- match.arg(vcov)
  tests <- match.arg(tests, several.ok = TRUE)
  D <- object$design
  n <- D$n; K <- D$K
  fn <- .pc_fns(D)
  theta_ml <- object$theta
  V_ml <- .pc_strip_vcov(.pc_vcov(theta_ml, fn, "analytic", NULL, object$scale, D$w))
  P <- .pc_pairs(D, pairs)
  if (any(vapply(P, function(pr) length(pr$cases), integer(1)) < 2L)) stop("Some pair of alternatives has fewer than two cases; the pairwise tests are not computable.", call. = FALSE)
  rows <- list(); details <- list(pairs = pairs)
  if ("cl" %in% tests) {
    cl <- .pc_cl_fit(P, theta_ml, n, K)
    Jm <- crossprod(cl$scores * sqrt(D$w))
    Hinv <- tryCatch(solve(cl$H), error = function(e) { .pc_pass_limit(e); MASS::ginv(cl$H) })
    V_cl <- Hinv %*% Jm %*% Hinv
    Vd <- if (vcov == "difference") V_cl - V_ml else {
      ## the covariance of theta_ML - theta_CL from the stacked influence
      ## functions of the two estimators: positive semi-definite by construction
      Psi <- fn$scores(theta_ml) %*% V_ml - cl$scores %*% Hinv
      crossprod(if (!is.null(D$cluster)) as.matrix(rowsum(Psi * D$w, D$cluster)) else Psi * sqrt(D$w))
    }
    qf <- .pc_quad_form(theta_ml - cl$theta, Vd, ref = diag(V_ml))
    if (isTRUE(qf$degenerate)) .pc_degenerate_warning("The composite-likelihood Hausman test (Q1)")
    if (!cl$converged) warning("The composite-likelihood estimator did not converge in 50 Newton steps; Q1 is unreliable.", call. = FALSE)
    rows[[length(rows) + 1L]] <- data.frame(test = paste0("Fok-Paap CL Hausman (Q1, ", vcov, ")"), statistic = qf$stat, df = qf$df,
                                            p.value = stats::pchisq(max(qf$stat, 0), qf$df, lower.tail = FALSE), psd = qf$psd, pairs = pairs, stringsAsFactors = FALSE)
    details$theta_cl <- stats::setNames(cl$theta, D$par_names); details$vcov_cl <- V_cl; details$cl_converged <- cl$converged
  }
  do_gmm <- "gmm" %in% tests
  if (do_gmm) {
    nz <- D$p + D$q
    M <- length(P) * nz
    df2 <- M - K
    if (df2 < 1L) {
      msg <- paste0("No overidentifying restrictions with ", pairs, " pairs (", M, " moments for ", K, " parameters); the GMM test needs pairs = \"all\"",
                    if (pairs == "all") " and alternative-specific covariates or J > 3" else "", ".")
      if (identical(tests, "gmm")) stop(msg, call. = FALSE)
      warning(msg, " Only the CL Hausman test is reported.", call. = FALSE)
      do_gmm <- FALSE
    }
  }
  if (do_gmm) {
    mom <- .pc_gmm_moments(P, theta_ml, n, K, nz, D$w)
    nw <- mom$nw
    Omega <- mom$Omega
    rc <- rcond(Omega)
    Winv <- function(Om) tryCatch(solve(Om), error = function(e) { .pc_pass_limit(e); MASS::ginv(Om) })
    if (gmm == "projected") {
      W <- Winv(Omega); gbar <- mom$gbar; Dm <- mom$D
      Pm <- W - W %*% Dm %*% solve(crossprod(Dm, W %*% Dm)) %*% t(Dm) %*% W
      stat <- nw * as.numeric(t(gbar) %*% Pm %*% gbar)
      theta_g <- theta_ml; iters <- 0L; conv <- TRUE
    } else {
      theta_g <- theta_ml; conv <- FALSE; obj_old <- Inf
      n_updates <- if (gmm == "two-step") 1L else maxit
      for (iters in seq_len(n_updates)) {
        mom_g <- .pc_gmm_moments(P, theta_g, n, K, nz, D$w)
        W <- Winv(mom_g$Omega)
        ## Gauss-Newton steps on the GMM objective with this weight matrix (the
        ## objective is close to quadratic in theta: a few steps suffice)
        th <- theta_g
        for (inner in seq_len(20L)) {
          mm <- .pc_gmm_moments(P, th, n, K, nz, D$w)
          gbar <- mm$gbar; Dm <- mm$D
          step <- tryCatch(solve(crossprod(Dm, W %*% Dm), crossprod(Dm, W %*% gbar)), error = function(e) { .pc_pass_limit(e); NULL })
          if (is.null(step)) break
          th <- th - as.vector(step)
          if (max(abs(step)) < 1e-8 * max(1, max(abs(th)))) break
        }
        gb <- .pc_gmm_moments(P, th, n, K, nz, D$w)$gbar; obj <- nw * as.numeric(t(gb) %*% W %*% gb)
        ## converged when the estimate or the objective stops moving (iterated GMM
        ## can cycle at round-off level between successive weight matrices)
        if (gmm == "two-step") { theta_g <- th; conv <- TRUE; break }
        if (max(abs(th - theta_g)) < 1e-6 * max(1, max(abs(th))) || abs(obj - obj_old) < 1e-8 * max(1, abs(obj))) { theta_g <- th; conv <- TRUE; break }
        theta_g <- th; obj_old <- obj
      }
      if (!conv) warning("The iterated GMM estimator did not settle within ", maxit, " weight-matrix updates; the last iterate is used. ",
                         "Consider gmm = \"two-step\" or \"projected\".", call. = FALSE)
      mom_g <- .pc_gmm_moments(P, theta_g, n, K, nz, D$w)
      W <- Winv(mom_g$Omega); gbar <- mom_g$gbar
      stat <- nw * as.numeric(t(gbar) %*% W %*% gbar)
      rc <- rcond(mom_g$Omega)
    }
    rows[[length(rows) + 1L]] <- data.frame(test = paste0("Fok-Paap GMM overidentification (Q2, ", gmm, ")"), statistic = stat, df = df2,
                                            p.value = stats::pchisq(max(stat, 0), df2, lower.tail = FALSE), psd = NA, pairs = pairs, stringsAsFactors = FALSE)   # no covariance difference (as in iia_tests())
    details$theta_gmm <- stats::setNames(theta_g, D$par_names); details$gmm_iterations <- iters; details$gmm_converged <- conv
    details$moment_rcond <- rc; details$n_moments <- M
    if (is.finite(rc) && rc < 1e-12) warning("The moment covariance matrix is nearly singular (reciprocal condition number ", format(rc, digits = 2),
                                             "): with few or weak covariates the pairwise overidentifying restrictions carry little information.", call. = FALSE)
  }
  out <- do.call(rbind, rows); rownames(out) <- NULL
  attr(out, "method") <- paste0("Fok-Paap (2025) misspecification tests, ", pairs, " pairs")
  attr(out, "details") <- details
  class(out) <- c("iia_test", "data.frame")
  out
}

#' The battery of IIA tests
#'
#' Runs [hausman_mcfadden()] (both covariance forms), [small_hsiao()], and
#' [fp_tests()] on one fitted [mnl()] model and collects the results in one
#' table. See the individual tests for definitions and caveats; in particular,
#' with case-specific covariates only, every test draws its power from the
#' covariates' nonlinear information, and the tests are not powerful against
#' every departure from IIA.
#'
#' @param object A fitted [mnl()] model with at least three alternatives.
#'   A fit with diverging coefficients (separation) draws a warning: the
#'   contrasts then compare estimates that do not exist and can reject
#'   spuriously.
#' @param drop Alternatives for the Hausman-McFadden and Small-Hsiao tests;
#'   default all.
#' @param splits Number of Small-Hsiao half-splits (default 10; the battery
#'   reports the median statistic over splits with the chi-square p-value of
#'   that median, a descriptive figure).
#' @param seed Optional seed for the Small-Hsiao splits.
#' @param pairs Pairs for the Fok-Paap tests (`"all"` or `"base"`).
#' @return A data frame of class `"iia_battery"` with `test`, `alternative`
#'   (dropped; blank for the pairwise tests), `statistic`, `df`, `p.value`,
#'   `psd` (whether the Hausman-McFadden covariance difference is positive
#'   semidefinite; `NA` where it does not apply), and `share_separated` (for
#'   the Small-Hsiao rows, the share of half-sample splits with a separated
#'   fit; `NA` otherwise).
#' @examples
#' set.seed(1)
#' d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
#' iia_tests(mnl(y ~ x1 + x2, data = d), seed = 1)
#' @seealso [hausman_mcfadden()], [small_hsiao()], [fp_tests()]
#' @family IIA tests
#' @export
iia_tests <- function(object, drop = NULL, splits = 10L, seed = NULL, pairs = c("all", "base")) {
  .pc_check_mnl(object, "The IIA battery")
  .pc_no_firth(object, "The IIA battery")
  pairs <- match.arg(pairs)
  .pc_separation_warning(object, "The IIA battery")
  ## the component tests' own separation warnings would repeat the battery's
  once <- function(expr) withCallingHandlers(expr, warning = function(w) if (grepl("the fit has diverging coefficients", conditionMessage(w), fixed = TRUE)) invokeRestart("muffleWarning"))
  hm <- once(hausman_mcfadden(object, drop = drop, vcov = "difference"))
  hj <- once(hausman_mcfadden(object, drop = drop, vcov = "joint"))
  sh <- once(small_hsiao(object, drop = drop, splits = splits, seed = seed))
  fp <- once(fp_tests(object, pairs = pairs))
  ## the median row's p-value is the chi-square p-value of the median statistic,
  ## so that the row is consistent (descriptive: the median of split statistics
  ## is not chi-square; boot_pvalue() calibrates it)
  ## psd applies to the covariance-difference rows only (NA elsewhere);
  ## share_separated to the Small-Hsiao rows only
  sh_rows <- if (splits > 1L) { s <- attr(sh, "summary"); data.frame(test = "Small-Hsiao (median over splits)", alternative = s$alternative, statistic = s$median_statistic, df = s$df,
                                                                      p.value = stats::pchisq(s$median_statistic, s$df, lower.tail = FALSE), psd = NA, share_separated = s$share_separated) }
            else data.frame(test = "Small-Hsiao", alternative = sh$alternative, statistic = sh$statistic, df = sh$df, p.value = sh$p.value, psd = NA, share_separated = as.numeric(sh$separated))
  q1 <- grepl("Q1", fp$test)
  out <- rbind(data.frame(test = "Hausman-McFadden (difference)", alternative = hm$alternative, statistic = hm$statistic, df = hm$df, p.value = hm$p.value, psd = hm$psd, share_separated = NA_real_),
               data.frame(test = "Hausman-McFadden (joint)", alternative = hj$alternative, statistic = hj$statistic, df = hj$df, p.value = hj$p.value, psd = hj$psd, share_separated = NA_real_),
               sh_rows,
               data.frame(test = fp$test, alternative = "", statistic = fp$statistic, df = fp$df, p.value = fp$p.value, psd = ifelse(q1, fp$psd, NA), share_separated = NA_real_))
  rownames(out) <- NULL
  attr(out, "method") <- "Tests of IIA / misspecification of the multinomial logit"
  class(out) <- c("iia_battery", "iia_test", "data.frame")
  out
}

#' @export
print.iia_test <- function(x, digits = 3, ...) {
  cat(attr(x, "method"), "\n")
  y <- as.data.frame(x)
  if ("statistic" %in% names(y)) y$statistic <- round(y$statistic, digits)
  if ("p.value" %in% names(y)) y$p.value <- format.pval(y$p.value, digits = digits)
  if ("share" %in% names(y)) y$share <- round(y$share, digits)
  print.data.frame(y, row.names = FALSE)
  if (!is.null(attr(x, "note"))) cat("Note:", attr(x, "note"), "\n")
  if (!is.null(attr(x, "summary"))) { cat("\nAcross splits:\n"); s <- attr(x, "summary"); s$median_statistic <- round(s$median_statistic, digits); s$median_p <- round(s$median_p, digits); print.data.frame(s, row.names = FALSE) }
  if (is.null(attr(x, "note")) && !is.null(x$psd) && any(!x$psd, na.rm = TRUE))
    cat("psd = FALSE: the covariance difference was not positive semi-definite for that row. A negative statistic is no evidence against the model (p-value 1); a positive one from a variance-difference form is unreliable (see the test's help).\n")
  if (!is.null(x$separated) && any(x$separated, na.rm = TRUE)) cat("separated = TRUE: a fit behind that row has diverging coefficients (quasi-separation in a half-sample); its statistic may reflect the separation rather than IIA.\n")
  if (!is.null(x$share_separated) && any(x$share_separated > 0, na.rm = TRUE)) cat("share_separated > 0: some Small-Hsiao splits had diverging half-sample fits; see small_hsiao(..., splits = ) for the split-level rows.\n")
  if (!is.null(x$test) && any(x$test == "Small-Hsiao (median over splits)"))
    cat("Small-Hsiao (median over splits): descriptive, not a test at its level (it rejected in none of 300 samples in the package's size study); for the test, run one split with its seed fixed in advance (small_hsiao(m, seed = ), or iia_tests(m, splits = 1, seed = )).\n")
  invisible(x)
}
