## ---------------------------------------------------------------------------
## Further existing specification tests for the multinomial logit: combining
## alternatives (Cramer and Ridder 1991; Anderson 1984 / the Wald form of
## Long and Freese), per-term joint tests, the Fagerland-Hosmer-Bofin
## goodness-of-fit test, the nested-logit score (LM) test of IIA, and a
## parametric-bootstrap calibration usable with any test in the package.
## ---------------------------------------------------------------------------

## selector matrix rows: the coefficients of term columns `cols` in alternative
## j's equation (0 for the base)
.pc_beta_index <- function(object, j, cols) {
  D <- object$design
  if (j == D$base) return(rep(NA_integer_, length(cols)))
  D$idx$beta[cols, D$levels[j]]
}

## the Wald statistic th' V^{-1} th, NA when V is singular (reciprocal
## condition number of the diagonally normalized matrix below 1e-12): a
## coefficient at a boundary (a separated level) has no variance, and the
## generalized inverse the fit falls back on there is not one. The system is
## solved normalized, (th/d)' (V/dd')^{-1} (th/d): with variances on very
## different scales solve() would refuse V itself (an unnormalized condition
## number of 1e-17 on Cao 2025 and Huth et al. 2011) though it is well posed.
.pc_wald <- function(th, V) {
  if (!all(is.finite(V))) return(NA_real_)
  dn <- sqrt(pmax(abs(diag(V)), .Machine$double.eps))
  Vn <- V / tcrossprod(dn)
  rc <- tryCatch(rcond(Vn), error = function(e) { .pc_pass_limit(e); 0 })
  if (!is.finite(rc) || rc < 1e-12) return(NA_real_)
  z <- th / dn
  as.numeric(crossprod(z, solve(Vn, z)))
}
.pc_wald_na_warning <- function(what) {
  if (length(what)) warning("The covariance of the coefficients of ", paste(what, collapse = ", "), " is singular (a coefficient at a boundary, such as a ",
                            "separated level; see the fit's warnings), so ", if (length(what) > 1L) "their" else "its", " Wald statistic is NA; type = \"lr\" refits instead.", call. = FALSE)
}

#' Test whether two alternatives can be combined
#'
#' Tests the hypothesis that two alternatives are indistinguishable with
#' respect to the covariates: all coefficients other than the constants are
#' equal in the two equations (when one of the two is the base, that the
#' other's coefficients are zero). The Wald form uses the fitted covariance;
#' the likelihood-ratio form (Cramer and Ridder 1991) refits the multinomial
#' logit with the two alternatives pooled and adds the binomial constant-only
#' log-likelihood of the split within the pooled alternative:
#' \eqn{\ell_R = \ell_{pooled} + n_j \log(n_j/n_{jk}) + n_k \log(n_k/n_{jk})}.
#' Both are chi-square on the number of non-constant terms (Stata's
#' `mlogtest, combine` and `lrcomb`): the restriction leaves the two constants
#' free, through the pooled equation's constant and the split share. An
#' alternative that is never chosen contributes \eqn{0 \log 0 = 0}.
#'
#' The likelihood-ratio form under separation and penalties:
#' * For a maximum-likelihood fit of separated data it compares the suprema
#'   of the log-likelihoods (see [lr_test()]).
#' * For a penalized fit (`firth = TRUE` or a prior) it is the penalized
#'   likelihood ratio, the restriction imposed on the full model's penalized
#'   log-likelihood (no pooled refit).
#'
#' The likelihood-ratio form assumes independent cases (a fit's `cluster` does not enter it); the Wald form uses the fit's covariance, cluster-robust when the fit has one.
#'
#' @param object A fitted [mnl()] model (for the LR form: case-specific
#'   covariates only, in either layout, and every alternative available; the
#'   Wald form works for any multinomial-logit fit).
#' @param alternatives Two alternatives (labels). `NULL` (default) tests every
#'   pair.
#' @param type `"wald"` (default) or `"lr"`.
#' @return A data frame of class `"iia_test"` with one row per pair:
#'   `pair`, `statistic`, `df`, `p.value`.
#' @references
#' Cramer, J.S. and Ridder, G. (1991). Pooling states in the multinomial logit
#' model. \emph{Journal of Econometrics}, 47, 267-272.
#' \doi{10.1016/0304-4076(91)90102-J}
#'
#' Long, J.S. and Freese, J. (2014). \emph{Regression Models for Categorical
#' Dependent Variables Using Stata}, 3rd ed. Stata Press.
#' @examples
#' set.seed(1)
#' ## alternatives b and c share their slopes
#' d <- rmnl(1000, beta = cbind(b = c(0.5, -1, 0.8), c = c(0.3, -1, 0.8)))
#' m <- mnl(y ~ x1 + x2, data = d)
#' combine_test(m)
#' combine_test(m, c("b", "c"), type = "lr")
#' @seealso [term_test()], [iia_tests()]
#' @family specification tests
#' @export
combine_test <- function(object, alternatives = NULL, type = c("wald", "lr")) {
  .pc_check_mnl(object, "combine_test()")
  type <- match.arg(type)
  D <- object$design
  if (D$p == 0L) stop("The model has no case-specific covariates; there is nothing to combine.", call. = FALSE)
  cols <- setdiff(colnames(D$X), "(Intercept)")
  if (!length(cols)) stop("The model has only constants; the combine test needs covariates.", call. = FALSE)
  pairs <- if (is.null(alternatives)) utils::combn(D$levels, 2L, simplify = FALSE) else {
    if (length(alternatives) != 2L || !all(alternatives %in% D$levels) || alternatives[1L] == alternatives[2L])
      stop("'alternatives' must be two different response levels.", call. = FALSE)
    list(as.character(alternatives))
  }
  firth <- .pc_penalized(object)                       # Firth's penalty or a prior
  ## a penalized fit's LR form restricts the full model directly (below), so
  ## the pooling conditions of the Cramer-Ridder form do not arise
  if (type == "lr" && !firth && D$layout == "long" && D$q > 0L)
    stop("The likelihood-ratio form pools two alternatives and is defined for case-specific covariates only; use type = \"wald\".", call. = FALSE)
  if (type == "lr" && !firth && any(D$avail == 0L))
    stop("The likelihood-ratio form needs every alternative available to every case (pooling two alternatives presumes both were in the choice set); use type = \"wald\".", call. = FALSE)
  ## the maximum-likelihood forms compare suprema (DESIGN.md Phase 4c): under
  ## separation the fits' values are where the optimizer stopped
  sup_full <- if (type == "lr" && !firth) .pc_sup_loglik(D, object$theta, object$loglik) else NULL
  sup_used <- !is.null(sup_full) && sup_full$type %in% c("complete", "quasi-complete")
  rows <- lapply(pairs, function(pr) {
    j <- match(pr[1L], D$levels); k <- match(pr[2L], D$levels)
    ## the restriction: equal non-constant coefficients of j and k (zero when one is the base)
    ij <- .pc_beta_index(object, j, cols); ik <- .pc_beta_index(object, k, cols)
    Rm <- matrix(0, length(cols), length(object$theta))
    for (t in seq_along(cols)) { if (!is.na(ij[t])) Rm[t, ij[t]] <- 1; if (!is.na(ik[t])) Rm[t, ik[t]] <- Rm[t, ik[t]] - 1 }
    if (type == "lr" && firth) {
      ## penalized likelihood ratio: Cramer and Ridder's restriction imposed on
      ## the full model's penalized log-likelihood (Heinze and Schemper 2001)
      stat <- .pc_firth_lr(object, Rm, start = object$theta)
    } else if (type == "wald") {
      if (is.null(object$vcov)) stop("The Wald form needs a covariance matrix (se != \"none\").", call. = FALSE)
      delta <- as.vector(Rm %*% object$theta)
      Vd <- Rm %*% object$vcov_internal %*% t(Rm)
      stat <- .pc_wald(delta, Vd)
    } else {
      ## pooled fit: merge k into j
      Dp <- D
      keep_alt <- setdiff(seq_len(D$J), k)
      yp <- D$y; yp[yp == k] <- j
      Dp$y <- match(yp, keep_alt); Dp$J <- D$J - 1L; Dp$levels <- D$levels[keep_alt]
      Dp$base <- match(if (D$base == k) j else D$base, keep_alt)
      Dp$avail <- D$avail[, keep_alt, drop = FALSE]; Dp$off <- D$off[, keep_alt, drop = FALSE]
      Dp$W <- D$W[rep(keep_alt, D$n) + rep((seq_len(D$n) - 1L) * D$J, each = Dp$J), , drop = FALSE]
      nb <- setdiff(seq_len(Dp$J), Dp$base)
      Dp$par_names <- c(as.vector(vapply(nb, function(a) paste0(Dp$levels[a], ":", colnames(D$X)), character(D$p))), colnames(D$W))
      Dp$idx <- list(beta = matrix(seq_len(D$p * (Dp$J - 1L)), nrow = D$p), gamma = if (D$q) D$p * (Dp$J - 1L) + seq_len(D$q) else integer(0), delta = integer(0), str = integer(0))
      Dp$K <- length(Dp$par_names); Dp$blocks <- c(rep("outcome", D$p * (Dp$J - 1L)), rep("generic", D$q)); Dp$fe_map <- NULL; Dp$fe_off <- NULL
      rp <- .pc_refit_quiet(Dp)
      sp <- .pc_sup_loglik(Dp, rp$coefficients, rp$loglik)
      if (sp$type %in% c("complete", "quasi-complete")) sup_used <<- TRUE
      lp <- sp$value
      ## the split of the pooled choices back into j and k (0 log 0 = 0: an
      ## alternative never chosen contributes nothing)
      nj <- sum(D$w[D$y == j]); nk <- sum(D$w[D$y == k])
      xlx <- function(a) if (a > 0) a * log(a / (nj + nk)) else 0
      lr_restricted <- lp + xlx(nj) + xlx(nk)
      stat <- 2 * (sup_full$value - lr_restricted)
    }
    data.frame(pair = paste(pr, collapse = " & "), statistic = stat, df = length(cols),
               p.value = stats::pchisq(max(stat, 0), length(cols), lower.tail = FALSE), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  if (type == "wald") .pc_wald_na_warning(out$pair[is.na(out$statistic)])
  attr(out, "method") <- paste0("Test of combining alternatives (", if (type == "wald") "Wald" else if (firth) paste0("penalized likelihood ratio, ", .pc_penalty_label(object)) else "likelihood ratio, Cramer-Ridder", ")",
                                if (sup_used) .pc_sup_note else "")
  class(out) <- c("iia_test", "data.frame")
  out
}

#' Joint tests of each term across equations
#'
#' For every term of the model, tests that its coefficients are zero in all
#' equations at once: the Wald form from the fitted covariance, or the
#' likelihood-ratio form by refitting without the term (Stata's `mlogtest,
#' wald` and `lr`). For a generic coefficient on an alternative-specific
#' covariate the test has one degree of freedom. The terms are those of the
#' outcome equations; for an inflated model the split-equation coefficients
#' are not included (test them by [lr_test()] between two [imnl()] fits that
#' differ in the split equation, or from the `infl_` rows of
#' [summary.choicemodel()]).
#'
#' The likelihood-ratio form works as follows:
#' * The model without the term is fitted as the estimators fit a model,
#'   multistart and saddle step included, so that both log-likelihoods are
#'   maxima; a warning says when that fit did not converge or still stops
#'   at a saddle.
#' * For a maximum-likelihood fit of separated data ([mnl()], [rologit()]
#'   without unit fixed effects or inflation), it compares the suprema of the
#'   log-likelihoods, not where the optimizer stopped (see [lr_test()]). The
#'   method line says so.
#' * For a penalized fit (`firth = TRUE` or a prior), it is the penalized
#'   likelihood ratio: the term's coefficients are set to zero in the full
#'   model's penalized log-likelihood.
#'
#' The Wald form of a penalized fit uses the penalized estimates and their
#' covariance; for a separated coefficient it can be far off (see
#' [summary.choicemodel()]).
#'
#' The likelihood-ratio form assumes independent cases (a fit's `cluster` does not enter it); the Wald form uses the fit's covariance, cluster-robust when the fit has one.
#'
#' @param object A fitted choice model.
#' @param type `"wald"` (default) or `"lr"`.
#' @return A data frame of class `"iia_test"` with one row per term:
#'   `term`, `statistic`, `df`, `p.value`.
#' @examples
#' set.seed(1)
#' d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
#' term_test(mnl(y ~ x1 + x2, data = d))
#' @seealso [combine_test()], [lr_test()]
#' @family specification tests
#' @export
term_test <- function(object, type = c("wald", "lr")) {
  if (!inherits(object, "choicemodel")) stop("'object' must be a fitted polychoice model.", call. = FALSE)
  type <- match.arg(type)
  D <- object$design
  ## terms of the case-specific part and of the alternative-specific part
  tx <- if (!is.null(D$tm$x)) D$tm$x$terms else NULL
  tw <- if (!is.null(D$tm$w)) D$tm$w$terms else NULL
  items <- list()
  if (!is.null(tx)) {
    lab <- attr(tx, "term.labels")
    M <- stats::model.matrix(tx, stats::model.frame(tx, D$data_case, na.action = stats::na.pass))
    asg <- attr(M, "assign")
    for (t in seq_along(lab)) {
      cols <- intersect(colnames(M)[asg == t], colnames(D$X))
      if (length(cols)) items[[length(items) + 1L]] <- list(term = lab[t], part = "x", cols = cols)
    }
  }
  if (!is.null(tw)) {
    lab <- attr(tw, "term.labels")
    for (t in seq_along(lab)) {
      ## the columns of a term come from the model matrix's assign attribute
      ## (a prefix match would let `w1` swallow `w10`)
      cols <- if (!is.null(D$tm$w$assign)) D$tm$w$cols[D$tm$w$assign == t]
              else grep(paste0("^", gsub("([][(){}.*+?^$|\\\\])", "\\\\\\1", lab[t]), "($|:)"), colnames(D$W), value = TRUE)
      cols <- intersect(cols, names(object$theta)); if (length(cols)) items[[length(items) + 1L]] <- list(term = lab[t], part = "w", cols = cols) }
  }
  if (!length(items)) stop("No covariate terms to test.", call. = FALSE)
  ## the maximum-likelihood form compares suprema (DESIGN.md Phase 4c): under
  ## separation the fits' values are where the optimizer stopped
  sup_full <- if (type == "lr" && !.pc_penalized(object)) .pc_sup_loglik(D, object$theta, object$loglik) else NULL
  sup_used <- !is.null(sup_full) && sup_full$type %in% c("complete", "quasi-complete")
  ## an alternative whose probability vanishes at the supremum contributes no
  ## coefficients to it (.pc_absent_alts; local review 2026-11-06, finding 10)
  absent <- if (type == "lr" && !.pc_penalized(object)) .pc_absent_alts(D) else integer(0)
  rows <- lapply(items, function(it) {
    idx <- if (it$part == "x") as.vector(D$idx$beta[it$cols, , drop = FALSE]) else match(it$cols, names(object$theta))
    if (type == "wald") {
      if (is.null(object$vcov_internal)) stop("The Wald form needs a covariance matrix (se != \"none\").", call. = FALSE)
      th <- object$theta[idx]; V <- object$vcov_internal[idx, idx, drop = FALSE]
      stat <- .pc_wald(th, V)
    } else if (.pc_penalized(object)) {
      ## penalized likelihood ratio: the term's coefficients set to zero in the
      ## full model's penalized log-likelihood (Heinze and Schemper 2001)
      Rm <- matrix(0, length(idx), length(object$theta)); Rm[cbind(seq_along(idx), idx)] <- 1
      stat <- .pc_firth_lr(object, Rm, start = replace(object$theta, idx, 0))
    } else {
      ## the restricted fit drops the term's columns from the design (no re-evaluation of the call)
      D2 <- D
      if (it$part == "x") D2$X <- D$X[, setdiff(colnames(D$X), it$cols), drop = FALSE]
      else D2$W <- D$W[, setdiff(colnames(D$W), it$cols), drop = FALSE]
      D2 <- .pc_reindex(D2)
      r2 <- .pc_refit_checked(D2)                     # with the saddle step, as the user's fit
      s2 <- .pc_sup_loglik(D2, r2$coefficients, r2$loglik)
      if (s2$type %in% c("complete", "quasi-complete")) sup_used <<- TRUE
      stat <- 2 * (sup_full$value - s2$value)
    }
    dfi <- if (type == "lr" && length(absent) && it$part == "x") length(it$cols) * (D$J - 1L - length(absent)) else length(idx)
    data.frame(term = it$term, statistic = stat, df = dfi, p.value = stats::pchisq(max(stat, 0), dfi, lower.tail = FALSE), stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows); rownames(out) <- NULL
  if (type == "wald") .pc_wald_na_warning(out$term[is.na(out$statistic)])
  attr(out, "method") <- paste0("Joint tests of each term across equations (", if (type == "wald") "Wald" else if (.pc_penalized(object)) paste0("penalized likelihood ratio, ", .pc_penalty_label(object)) else "likelihood ratio", ")",
                                if (sup_used) .pc_sup_note else "",
                                if (length(absent)) paste0("; never chosen, its probability zero at the suprema and its coefficients not counted in df: ", paste(D$levels[absent], collapse = ", ")) else "")
  class(out) <- c("iia_test", "data.frame")
  out
}

#' Hosmer-Lemeshow-type goodness-of-fit test for multinomial models
#'
#' The test of Fagerland, Hosmer and Bofin (2008): cases are sorted by the
#' fitted probability of the base alternative and split into `groups` groups
#' of equal size; the weighted observed and expected counts of every
#' alternative in every group give a Pearson statistic that is approximately
#' chi-square on \eqn{(g - 2)(J - 1)} degrees of freedom under a correct
#' model. Small expected counts (below 5) are reported; the test is a
#' calibration check, not a test of IIA.
#'
#' @param object A fitted choice model with independent single choices: not a
#'   [remnl()] or [rologit()] fit, and not a penalized fit (`firth = TRUE` or a
#'   prior). The cases are taken as independent; a fit's `cluster` does not
#'   enter the statistic.
#' @param groups Number of groups (default 10).
#' @return A list of class `"hl_test"` with `statistic`, `df`, `p.value`,
#'   `observed` and `expected` (group by alternative), and `n_small` (expected
#'   counts below 5).
#' @references Fagerland, M.W., Hosmer, D.W. and Bofin, A.M. (2008).
#'   Multinomial goodness-of-fit tests for logistic regression models.
#'   \emph{Statistics in Medicine}, 27, 4238-4253. \doi{10.1002/sim.3202}
#' @examples
#' set.seed(1)
#' d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
#' hl_test(mnl(y ~ x1 + x2, data = d))
#' @family specification tests
#' @export
hl_test <- function(object, groups = 10L) {
  if (!inherits(object, "choicemodel")) stop("'object' must be a fitted polychoice model.", call. = FALSE)
  if (identical(object$family, "remnl")) stop("hl_test() assumes independent cases; it does not apply to a random-intercept fit (the cases of a unit are dependent).", call. = FALSE)
  if (identical(object$family, "rologit")) stop("hl_test() is a test of single-choice fits; it does not apply to a rank-ordered logit.", call. = FALSE)
  .pc_no_firth(object, "hl_test()", "its chi-square reference is derived for maximum-likelihood fitted probabilities", penalized_forms = FALSE)
  groups <- .pc_check_whole(groups, "groups", 3)
  D <- object$design
  P <- object$fitted.values; w <- object$weights; y <- object$y
  ord <- order(P[, D$base])
  cw <- cumsum(w[ord]); g <- pmin(groups, floor(cw / (sum(w) / groups) + 1e-12) + 1L)
  grp <- integer(D$n); grp[ord] <- g
  O <- E <- matrix(0, groups, D$J, dimnames = list(paste0("g", seq_len(groups)), D$levels))
  for (i in seq_len(D$n)) { O[grp[i], y[i]] <- O[grp[i], y[i]] + w[i]; E[grp[i], ] <- E[grp[i], ] + w[i] * P[i, ] }
  ok <- E > 0
  stat <- sum((O[ok] - E[ok])^2 / E[ok])
  ## the groups that received cases and the alternatives with expected counts
  ## (ties or heavy weights can leave a group empty; local review 2026-11-06)
  g_eff <- sum(rowSums(E) > 0); j_eff <- sum(colSums(E) > 0)
  if (g_eff < 3L) stop("Fewer than three groups received cases; use fewer 'groups' or more data.", call. = FALSE)
  df <- (g_eff - 2L) * (j_eff - 1L)
  structure(list(statistic = stat, df = df, p.value = stats::pchisq(stat, df, lower.tail = FALSE), observed = O, expected = E,
                 n_small = sum(E < 5), groups = groups), class = "hl_test")
}

#' @export
print.hl_test <- function(x, digits = 3, ...) {
  cat("Fagerland-Hosmer-Bofin goodness-of-fit test (", x$groups, " groups by the base alternative's fitted probability)\n", sep = "")
  cat(sprintf("  statistic = %.*f on %d df, p = %s\n", digits, x$statistic, x$df, format.pval(x$p.value, digits = digits)))
  if (x$n_small) cat("  ", x$n_small, " expected count(s) below 5.\n", sep = "")
  invisible(x)
}

#' Nested-logit score test of IIA
#'
#' The Lagrange-multiplier (score) test of the multinomial logit against the
#' nested logit with the given nests, evaluated at the multinomial-logit
#' estimates where every inclusive-value parameter equals 1 (McFadden 1987;
#' Hausman and McFadden 1984, section 4). Only the restricted (MNL) fit is
#' needed: the score of the nested-logit log-likelihood with respect to the
#' inclusive-value parameters at \eqn{\lambda = 1} is combined with the
#' information matrix (by default its outer-product form, McFadden's regression
#' form, positive semidefinite by construction; or the differenced Hessian with
#' `information = "hessian"`). The test is directional, towards the stated nesting
#' structure, and the sign of each nest's score shows whether the data pull
#' its \eqn{\lambda} below 1 (correlated errors within the nest) or above.
#'
#' The statistic assumes independent cases: a fit's `cluster` does not enter it.
#'
#' @param object A fitted [mnl()] model.
#'   A fit with diverging coefficients (separation) draws a warning: the
#'   statistic is then evaluated where the optimizer stopped, at estimates
#'   that do not exist.
#' @param nests A named list assigning every alternative to a nest, as in
#'   [nlogit()].
#' @param equal_lambda If `TRUE`, one shared inclusive-value parameter (one
#'   degree of freedom).
#' @param information `"opg"` (default; the outer product of the scores, the form of
#'   McFadden's regression-based test, positive semidefinite by construction) or
#'   `"hessian"` (the observed information at the restricted estimate, which is
#'   not guaranteed positive definite there; a negative statistic is reported
#'   with a warning and no p-value).
#' @return A data frame of class `"iia_test"` with the statistic, degrees of
#'   freedom, and p-value; the attribute `"scores"` holds the per-nest scores.
#' @references
#' Hausman, J.A. and McFadden, D. (1984). Specification tests for the
#' multinomial logit model. \emph{Econometrica}, 52, 1219-1240.
#' \doi{10.2307/1910997}
#'
#' McFadden, D. (1987). Regression-based specification tests for the
#' multinomial logit model. \emph{Journal of Econometrics}, 34, 63-82.
#' \doi{10.1016/0304-4076(87)90067-4}
#' @examples
#' set.seed(2)
#' d <- rnlogit(2000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)),
#'              nests = list(n1 = c("a", "b"), n2 = c("c", "d")), lambda = c(0.5, 0.5))
#' m <- mnl(y ~ x1 + x2, data = d)
#' nested_lm_test(m, nests = list(n1 = c("a", "b"), n2 = c("c", "d")))
#' @seealso [nlogit()], [lr_test()], [iia_tests()], [boot_pvalue()]
#' @family IIA tests
#' @export
nested_lm_test <- function(object, nests, equal_lambda = FALSE, information = c("opg", "hessian")) {
  .pc_check_mnl(object, "The nested-logit score test")
  .pc_no_firth(object, "The nested-logit score test", "a score test is evaluated at the maximum-likelihood estimates, where the score of the free parameters vanishes")
  .pc_separation_warning(object, "The nested-logit score test", .pc_score_sep_why)
  information <- match.arg(information)
  D <- object$design
  S <- .pc_structural("nlogit", list(nests = nests, equal = isTRUE(equal_lambda)), D$levels, D$base)
  D2 <- D; D2$family <- "nlogit"; D2$str <- S
  D2$par_names <- c(D$par_names, S$names); D2$blocks <- c(D$blocks, S$block)
  D2$idx$str <- D$K + seq_len(S$n); D2$K <- D$K + S$n
  fn <- .pc_fns(D2)
  th0 <- c(object$theta, rep(0, S$n))
  Sc <- fn$scores(th0)
  s <- colSums(D$w * Sc)
  ## observed information: the Hessian of the negative log-likelihood (fn$grad is its gradient), or the outer product of the scores
  I <- if (information == "opg") crossprod(Sc * sqrt(D$w)) else {
    ## differences in the standardized coordinates (as for the fit's own
    ## covariance), so that the steps scale with each parameter's covariate,
    ## mapped back: H = diag(sc) H_s diag(sc)
    sc <- .pc_scales(D2)
    Hs <- .pc_hess_fd(function(ts) fn$grad(ts / sc) / sc, th0 * sc); attr(Hs, "asym") <- NULL
    sc * Hs * rep(sc, each = length(sc))
  }
  ## an alternative whose probability vanishes at the supremum (never chosen,
  ## .pc_absent_alts) leaves the model there: its coefficients, whose scores
  ## and information are round-off, leave the statistic, and a nest it leaves
  ## with one alternative has no inclusive value to test (local review
  ## 2026-11-06)
  df <- S$n; dropped <- character(0); s_full <- s
  ab <- .pc_absent_alts(D)
  if (length(ab)) {
    if (D$base %in% ab) stop("The base alternative ", D$levels[D$base], " is never chosen; refit with another 'ref'.", call. = FALSE)
    live <- vapply(S$nests, function(a) sum(!(match(a, D$levels) %in% ab)), integer(1))
    multi <- lengths(S$nests) > 1L
    lam_keep <- if (S$equal) any(live[multi] > 1L) else live[multi] > 1L
    if (!any(lam_keep)) stop("Once the never-chosen alternative(s) ", paste(D$levels[ab], collapse = ", "), " are taken out (their probability is zero at the ",
                             "maximum), no nest has two alternatives; there is no inclusive value to test.", call. = FALSE)
    nb <- setdiff(seq_len(D$J), D$base)
    gone <- c(if (D$p) as.vector(D$idx$beta[, match(ab, nb), drop = FALSE]), D2$idx$str[!lam_keep])
    keep <- setdiff(seq_len(D2$K), gone)
    dropped <- S$names[!lam_keep]
    s <- s[keep]; I <- I[keep, keep, drop = FALSE]; df <- sum(lam_keep)
  }
  stat <- .pc_wald(s, I)
  if (is.na(stat)) stop("The information matrix at the multinomial-logit estimate is singular (a coefficient at a boundary, such as a separated level; ",
                        "see the fit's warnings), so the LM statistic is not defined.", call. = FALSE)
  if (stat < 0) warning("The observed information at the multinomial-logit estimate is not positive definite in the nest directions (negative statistic); use information = \"opg\".", call. = FALSE)
  out <- data.frame(test = paste0("Nested-logit LM test (", information, ")"), statistic = stat, df = df,
                    p.value = if (stat < 0) NA_real_ else stats::pchisq(stat, df, lower.tail = FALSE), stringsAsFactors = FALSE)
  attr(out, "method") <- paste0("Lagrange-multiplier test of IIA against a nested logit",
                                if (length(ab)) paste0(" (never chosen, with probability zero at the maximum, and taken out: ", paste(D$levels[ab], collapse = ", "),
                                                       if (length(dropped)) paste0("; not tested: ", paste(dropped, collapse = ", "), ", whose nest keeps one alternative") else "", ")") else "")
  attr(out, "scores") <- stats::setNames(s_full[D2$idx$str], S$names)
  class(out) <- c("iia_test", "data.frame")
  out
}

#' Tests for captive (inflated) choosers: score and likelihood ratio
#'
#' Tests the multinomial logit against the dogit model of Gaudry and Dagenais
#' (1979), in which a share of the decision makers is captive to an
#' alternative: \eqn{P^*_{ij} = (P_{ij} + \theta_j) / (1 + \sum_m \theta_m)},
#' \eqn{\theta_j \ge 0}, with the sum over the case's choice set. The
#' multinomial logit is \eqn{\theta = 0}. The joint statistic is the Lagrange
#' multiplier statistic of Tse (1987, eqs. 4-5), with the expected information,
#' extended here to choice sets that vary across cases and to
#' alternative-specific covariates; the one-sided z statistic for each
#' alternative is this package's addition. The score form needs only the
#' [mnl()] fit; the likelihood-ratio form (`type = "lr"`, below) refits. Both
#' need at least three alternatives. The score of case \eqn{i} in
#' \eqn{\theta_j} at
#' \eqn{\theta = 0} is \eqn{1\{y_i = j\}/P_{ij} - 1}: it is positive when the
#' cases that chose \eqn{j} are ones the logit gave a low probability to,
#' which is what captive choosers produce. The expected information at the
#' logit estimates, with the coefficients partialled out, gives a one-sided z
#' statistic for each alternative (captivity cannot be negative) and a joint
#' chi-square statistic over the alternatives tested. A rejection for the
#' base alternative of a vote-choice model is evidence of routine choosers of
#' it, the case the inflated multinomial logit [imnl()] models.
#'
#' The constant captivity parameters make the logit an ordinary boundary
#' point of the dogit (\eqn{\theta = 0}), so the score test is regular there,
#' unlike the comparison with [imnl()], whose split coefficients are not
#' identified at \eqn{\phi = 1}; it can detect covariate-dependent captivity
#' through its average (but see the loss of power below). Tse refers the joint statistic to the
#' chi-square with \eqn{J} degrees of freedom, reads the boundary null as a
#' test for arbitrarily small \eqn{\theta}, and found its size at 0.05
#' acceptable in his simulations (three alternatives, 50 and 300 cases).
#' When the logit's fitted probabilities of the tested alternatives do not vary
#' across cases (no case-specific and no alternative-specific covariates),
#' captivity is not identified and the statistic is undefined.
#'
#' The weights \eqn{1/P_{ij}} make the statistic heavy-tailed when some fitted
#' probabilities are small: a single case that chose an alternative the logit
#' gave a probability of 0.001 adds about 1000 to its score. In the package's
#' simulations under the logit (three and four alternatives, 500 and 2,000
#' cases, smallest fitted probabilities around 0.0005) the one-sided tests
#' rejected in 2.8 to 6.3 percent of samples at the 0.05 level and the joint
#' test in 5.7 to 7.3 percent, but at the 0.01 level they rejected in 1.3 to
#' 3.7 percent: for small significance levels calibrate by simulation with
#' [boot_pvalue()]. Against a baseline captive share of about 15 percent the
#' test of that alternative rejected in 47 percent of samples of 500 cases and
#' 94 percent of samples of 2,000. The expected information is used: the
#' outer product of the scores, which would also allow clustering, sized far
#' worse in the same simulations (one-sided rejection at 0.05 of 0.3 to 2.3
#' percent, joint 16 to 23 percent).
#'
#' The same weights cost power. A case whose fitted probability of an
#' alternative is small and that did not choose it adds -1 to the score but
#' \eqn{1/P_{ij} - 1} to its variance. On [pa2006] (Replication I's logit) the
#' one-sided z of abstention is 0.50, while [imnl()] with a constant split,
#' the test's own alternative, raises the log-likelihood by 10.5 on one
#' parameter, at a captive share of 0.20. The score is 37.3 (a numerical
#' derivative of the log-likelihood agrees), but its expected information
#' (before the coefficients are partialled out) is
#' 10,186, half of it from one respondent whose fitted probability of
#' abstaining is 0.0002 and who did not abstain; and along the captive share
#' the log-likelihood is convex at zero, so the gain lies away from the null,
#' where a local test does not look. A small statistic is therefore no
#' evidence against captivity. (On [fluvax] the two agree: z = 0.63 for
#' "don't know", and the constant-share fit gains 0.41.)
#'
#' `type = "lr"` measures that gain. For each alternative it fits the
#' inflated multinomial logit with that alternative inflated and an
#' intercept-only split equation, which is the dogit with a constant captive
#' share of it. The fit uses the logit's own cases, weights, offsets and
#' choice sets, [imnl()]'s search, and a start at the logit itself. It reports
#' the likelihood ratio against the logit and the captive share. The logit is
#' a captive share of zero, a boundary point of the share's space. The point
#' is a regular one when the other parameters are identified there and the
#' share's score is not a combination of theirs, that is, when the fitted
#' probabilities of the alternative vary in a way its constants and slopes
#' cannot absorb. The reference is then the 50:50 mixture of chi-square(0)
#' and chi-square(1) (Self and Liang 1987), whose p-value is half the
#' chi-square(1) tail (and 1 for a statistic of 0).
#'
#' * On `pa2006` it gives 20.95 for abstention (p = 2.4e-6, a captive share
#'   of 0.196).
#' * There, where one respondent carries half the information, a parametric
#'   bootstrap from the logit (199 replicates) put 0.61 of the statistics at
#'   zero (the mixture puts half there) and exceeded the mixture's 0.10, 0.05
#'   and 0.01 critical values in 0.101, 0.040 and 0.010 of the replicates:
#'   the upper tail agreed with the mixture within Monte Carlo error, and the
#'   reference was, if anything, conservative.
#' * Where few cases dominate the information (the attribute `fits` reports
#'   `top_share`, the largest row's share; with frequency weights a row of
#'   weight w stands for w cases), [boot_lr_test()], which accepts this pair,
#'   calibrates it.
#' * A saturated outcome design (a single factor, say) absorbs a constant
#'   share into the alternative's constants: the ratio is then 0 and the
#'   share is reported as not identified.
#' * A fit at which outcome coefficients diverge, or that stops at a saddle,
#'   keeps its statistic with a p-value of NA and a note, since the reference
#'   does not hold there. The fits' warnings are in `attr(, "fits")`.
#'
#' The form refits the inflated model once per alternative (about 20 seconds
#' each on `pa2006`) and has no joint row: with several captive shares on the
#' boundary the reference is a chi-square mixture whose weights depend on the
#' information. Each row is a separate test, so with every alternative tested
#' (the default) one can reject by chance. The comparison with an inflated model whose split equation
#' has covariates is not of this kind; [boot_lr_test()] calibrates it. Like
#' the other likelihood-ratio and score statistics of the package, both forms
#' assume independent cases and ignore a fit's `cluster`.
#'
#' @param object A fitted [mnl()] model.
#'   A fit with diverging coefficients (separation) draws a warning: the
#'   statistic is then evaluated where the optimizer stopped, at estimates
#'   that do not exist.
#' @param alternatives The alternatives to test (labels); default all.
#' @param type `"score"` (default; computed from the logit fit alone) or
#'   `"lr"` (the likelihood ratio against the constant-share inflated multinomial logit,
#'   one refit per alternative).
#' @return A data frame of class `"iia_test"`. With `type = "score"`, one row
#'   per alternative (`statistic` the one-sided z, `p.value` its upper-tail
#'   p-value; `NA` for an alternative no case chose, which carries no evidence
#'   on captivity to it and is left out of the joint test) and a last row for
#'   the joint test (`statistic` the chi-square, computed on the correlation
#'   scale of the scores, `df`, `p.value`). With `type = "lr"`, one row per alternative
#'   (`statistic` the likelihood ratio, `df` 1, `p.value` from the mixture,
#'   `share` the fitted captive share), and the attribute `fits`, one row per
#'   alternative:
#'   * `alternative`, `loglik` (the inflated fit's), `captive_share`;
#'   * `converged`, `boundary`, `n_modes` (the search's distinct maxima);
#'   * `top_share` (above);
#'   * `warnings` (the inflated fit's, joined), `note`;
#'   * `status`:
#'     - `"logit"`: the fit reached the logit; statistic 0, p-value 1, and
#'       share 0, or NA when the share is not identified;
#'     - `"interior"`: the mixture p-value;
#'     - `"degenerate"`: diverging outcome coefficients or a saddle; the
#'       statistic is kept and the p-value is NA;
#'     - `"failed"`: the fit failed or was refused; statistic and p-value NA.
#'   The attribute `note` collects the fits' notes.
#' @references
#' Gaudry, M.J.I. and Dagenais, M.G. (1979). The dogit model.
#' \emph{Transportation Research Part B}, 13, 105-111.
#' \doi{10.1016/0191-2615(79)90028-6}
#'
#' Self, S.G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
#' likelihood estimators and likelihood ratio tests under nonstandard
#' conditions. \emph{Journal of the American Statistical Association}, 82,
#' 605-610. \doi{10.1080/01621459.1987.10478472}
#'
#' Tse, Y.K. (1987). A diagnostic test for the multinomial logit model.
#' \emph{Journal of Business and Economic Statistics}, 5, 283-286.
#' \doi{10.1080/07350015.1987.10509587}
#' @examples
#' set.seed(6)
#' d <- rmnl(2000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
#'           delta = c(1.5, 0), inflate = "a")
#' m <- mnl(y ~ x1 + x2, data = d)
#' captivity_test(m)
#' \donttest{
#' ## the likelihood ratio against a constant captive share of a (one refit)
#' captivity_test(m, alternatives = "a", type = "lr")
#' }
#' @seealso [imnl()], [boot_lr_test()], [boot_pvalue()]
#' @family specification tests
#' @export
captivity_test <- function(object, alternatives = NULL, type = c("score", "lr")) {
  type <- match.arg(type)
  .pc_check_mnl(object, "captivity_test()")
  .pc_no_firth(object, "captivity_test()", if (type == "score") "a score test is evaluated at the maximum-likelihood estimates, where the score of the free parameters vanishes"
                                           else "the likelihood ratio compares maximum-likelihood fits", penalized_forms = FALSE)
  .pc_separation_warning(object, "captivity_test()", if (type == "score") .pc_score_sep_why else "the logit's log-likelihood is then where its optimizer stopped")
  D <- object$design
  tj <- if (is.null(alternatives)) seq_len(D$J) else match(as.character(alternatives), D$levels)
  if (!length(tj) || anyNA(tj) || anyDuplicated(tj)) stop("'alternatives' must name distinct response categories.", call. = FALSE)
  if (type == "lr") return(.pc_captivity_lr(object, tj))
  ci <- .pc_captivity_info(object, tj)
  s <- ci$s; V <- ci$V; dv <- ci$dv; ok <- ci$ok
  ## an alternative no case chose carries no evidence on captivity to it: its
  ## score is fixed near zero by its vanishing probabilities (local review
  ## 2026-11-06, finding 6), so it is not tested
  never <- tabulate(D$y[D$w > 0], D$J)[tj] == 0L
  ok <- ok & !never
  if (!any(ok)) stop("The captivity directions are not identified at the logit estimates (the fitted probabilities do not vary across cases: the model needs case-specific or alternative-specific covariates, ",
                     "and an alternative must have been chosen to be tested).", call. = FALSE)
  z <- ifelse(ok, s / sqrt(pmax(dv, 0)), NA_real_)
  ## the joint statistic on the correlation scale, with one eigenvalue rule
  ## for its rank and its inverse: the captive directions' information can
  ## differ by orders of magnitude, and a rule relative to the largest raw
  ## eigenvalue would drop well-identified directions (finding 6)
  zk <- z[ok]; Rk <- V[ok, ok, drop = FALSE] / tcrossprod(sqrt(dv[ok]))
  e <- eigen((Rk + t(Rk)) / 2, symmetric = TRUE)
  keep <- e$values > 1e-10 * max(e$values)
  rk <- sum(keep)
  joint <- sum(as.vector(crossprod(e$vectors[, keep, drop = FALSE], zk))^2 / e$values[keep])
  out <- data.frame(test = c(paste0("captive: ", D$levels[tj]), "joint"), statistic = c(z, joint),
                    df = c(rep(NA_integer_, length(tj)), rk),
                    p.value = c(stats::pnorm(z, lower.tail = FALSE), stats::pchisq(joint, rk, lower.tail = FALSE)), stringsAsFactors = FALSE)
  attr(out, "method") <- paste0("Score test of the multinomial logit against captive choosers (dogit; one-sided z per alternative, joint chi-square",
                                if (rk < sum(ok)) paste0(" on the ", rk, " identified directions") else "",
                                if (any(never)) paste0("; never chosen, not tested: ", paste(D$levels[tj][never], collapse = ", ")) else "", ")")
  attr(out, "scores") <- stats::setNames(s, D$levels[tj])
  class(out) <- c("iia_test", "data.frame")
  out
}

## the score of the captive odds of each alternative in tj at the logit, and
## their expected information after the logit's coefficients are partialled
## out (V); `ok` marks the identified shares (V's diagonal not negligible
## against the unadjusted information). Used by both forms of captivity_test().
.pc_captivity_info <- function(object, tj) {
  D <- object$design
  fn <- .pc_fns(D)
  th <- object$theta
  P <- fn$probs(th)$prob
  A <- D$avail == 1L; w <- D$w; n <- D$n
  Y <- matrix(0, n, D$J); Y[cbind(seq_len(n), D$y)] <- 1
  At <- A[, tj, drop = FALSE]; Pt <- P[, tj, drop = FALSE]
  S <- ifelse(At, Y[, tj, drop = FALSE] / ifelse(At, Pt, 1) - 1, 0)       # per-case scores in theta
  s <- colSums(w * S)
  ## expected information: theta-theta, theta-beta (the MNL score of a case that
  ## chose j, D_ij - Dbar_i), and beta-beta (the MNL information)
  Itt <- diag(colSums(w * At / ifelse(At, Pt, 1)), length(tj)) - crossprod(At * sqrt(w))
  Itb <- t(vapply(tj, function(j) {
    Dj <- D; Dj$y <- ifelse(A[, j], j, D$y)
    colSums((w * A[, j]) * .pc_fns(Dj)$scores(th))
  }, numeric(D$K)))
  Ibb <- fn$hess(th)
  ## the projection on the logit's information, solved in the diagonally
  ## normalized form; a direction with no information (a separated coefficient
  ## whose probabilities have underflowed) or a singular matrix stops the test
  dI <- abs(diag(Ibb)); dn <- sqrt(pmax(dI, .Machine$double.eps))
  Ibn <- Ibb / tcrossprod(dn)
  rc <- if (all(is.finite(Ibn))) tryCatch(rcond(Ibn), error = function(e) { .pc_pass_limit(e); 0 }) else 0
  proj <- if (is.finite(rc) && rc >= 1e-12 && all(dI > 1e-12 * max(dI))) tryCatch(solve(Ibn, t(Itb) / dn) / dn, error = function(e) { .pc_pass_limit(e); NULL }) else NULL
  if (is.null(proj))
    stop("The information matrix of the logit fit is singular (a coefficient at a boundary, such as a separated level, or an unidentified ",
         "design; see the fit's warnings); the score test projects on it and needs it invertible.", call. = FALSE)
  V <- Itt - Itb %*% proj
  dv <- diag(V)
  list(s = s, Itt = Itt, V = V, dv = dv, ok = dv > 1e-10 * pmax(diag(Itt), 1))
}

## captivity_test(type = "lr") (DESIGN.md Phase 4e): for each alternative the
## constant-share inflated multinomial logit on the logit's own design, by
## imnl()'s search and once more from the logit embedded at the boundary
## (phi = 1 - plogis(-30)), the higher kept; the likelihood ratio against the
## logit with the 50:50 chi-square(0) / chi-square(1) reference.
## * A fit that reaches the logit itself (its captive share collapses, or the
##   ratio is within max(1e-6, 2e-8 |l|) of zero: the optimizer's shortfall at
##   the boundary) gives 0 and p = 1; a share the singular information leaves
##   unidentified (a saturated outcome design) is NA, with a note.
## * A degenerate fit (outcome coefficients diverging, or a saddle the step
##   could not leave) keeps its statistic but gets p = NA and a note: the
##   chi-bar reference does not apply there (round 6, B-N1, D-new-4).
## * An alternative the inflated model refuses, or whose fit fails, gets NA
##   and its reason. The fits' warnings are kept in `fits$warnings`.
## * `fits$top_share`: the largest single case's share of the alternative's
##   expected information about the share at zero (before the coefficients
##   are partialled out), the dominance that strains the reference (round 6,
##   D-new-1).
.pc_captivity_lr <- function(object, tj) {
  D <- object$design
  cnt <- tabulate(D$y[D$w > 0], D$J)
  P0 <- .pc_fns(D)$probs(object$theta)$prob; A <- D$avail == 1L
  l0 <- object$loglik; tol <- max(1e-6, 2e-8 * max(1, abs(l0)))
  ## which shares are identified: the score form's rule, on the share's
  ## information after the logit's coefficients are partialled out (round 7,
  ## A-2); NULL when the logit's own information is singular
  idf <- tryCatch(.pc_captivity_info(object, tj)$ok, error = function(e) { .pc_pass_limit(e); NULL })
  rows <- lapply(seq_along(tj), function(i) {
    j <- tj[i]; lab <- D$levels[j]
    info <- D$w * A[, j] * (1 / pmax(P0[, j], .Machine$double.xmin) - 1); top <- if (sum(info) > 0) max(info) / sum(info) else NA_real_
    na <- function(why) data.frame(alternative = lab, loglik = NA_real_, captive_share = NA_real_, converged = NA, boundary = NA, n_modes = NA_integer_,
                                   top_share = top, warnings = "", note = why, status = "failed", stringsAsFactors = FALSE)
    if (any(cnt == 0L)) return(na(paste0("not fitted: alternative(s) never chosen (", paste(D$levels[cnt == 0L], collapse = ", "), "), which the inflated model refuses")))
    Di <- .pc_inflate_design(D, j)
    ## a fit with its own warnings, so that a refit that fails leaves the first
    ## fit's intact (round 7, B-6a)
    run <- function(start = NULL) {
      ww <- character(0)
      r <- tryCatch(withCallingHandlers(suppressMessages(.pc_fit_design(Di, se = "analytic", start = start)),
                                        warning = function(x) { ww <<- c(ww, conditionMessage(x)); invokeRestart("muffleWarning") }),
                    error = function(e) { .pc_pass_limit(e); e })
      list(fit = r, w = ww)
    }
    r1 <- run()
    failed <- inherits(r1$fit, "error")
    nm <- if (failed) NA_integer_ else as.integer(r1$fit$n_modes %||% 1L)        # the search's count of distinct local maxima
    ## the logit embedded at the boundary as a further start (as boot_lr_test()
    ## does), refit from where it ends when that is above the search's result by
    ## more than the tolerance in the ratio's units, or when the search failed
    ## (round 7, I-N4 and C-N6)
    te <- numeric(Di$K); te[-Di$idx$delta] <- object$theta; te[Di$idx$delta] <- 30
    oe <- tryCatch(suppressWarnings(.pc_optimize(.pc_fns(Di), te, .pc_scales(Di), 1000L, 1e-10)), error = function(e) { .pc_pass_limit(e); NULL })
    if (!is.null(oe) && is.finite(oe$value) && (failed || 2 * (-oe$value - r1$fit$loglik) > tol)) {
      r2 <- run(oe$par)
      if (!inherits(r2$fit, "error")) r1 <- r2
    }
    if (inherits(r1$fit, "error")) return(na(paste("the inflated fit failed:", conditionMessage(r1$fit))))
    f <- r1$fit
    lr <- 2 * (f$loglik - l0)
    cs <- 1 - stats::plogis(f$coefficients[["infl_(Intercept)"]])
    ## the rule of lr_test() and boot_lr_test() (.pc_captive_rule; local review
    ## 2026-11-06: this function had its own copy, which had drifted)
    cr <- .pc_captive_rule(lr, f, l0, 0L)
    dv <- setdiff(f$diverging, grep("^infl_", f$diverging, value = TRUE))   # the split intercept diverges at a collapse
    status <- if (cr$statistic == 0) "logit" else if (cr$degenerate) "degenerate" else "interior"
    unid <- if (!is.null(idf)) !idf[i] else isTRUE(f$ill_conditioned)
    note <- switch(status,
      logit = if (cr$below) "the inflated fit ended below the logit" else if (unid) "the share is not identified (the outcome design absorbs a constant share)" else "",
      degenerate = paste0("the inflated fit is degenerate (", cr$why, if (length(dv) > 3L) ", ..." else "",
                          "); the chi-bar reference does not apply, see boot_lr_test()"),
      interior = if (!isTRUE(f$converged)) "the inflated fit did not converge" else "")
    data.frame(alternative = lab, loglik = f$loglik, captive_share = if (status == "logit" && unid) NA_real_ else if (status == "logit") 0 else cs,
               converged = isTRUE(f$converged), boundary = isTRUE(f$boundary), n_modes = nm, top_share = top,
               warnings = paste(unique(r1$w), collapse = " | "), note = note, status = status, stringsAsFactors = FALSE)
  })
  fits <- do.call(rbind, rows)
  raw <- 2 * (fits$loglik - l0)
  if (any(fits$note == "the inflated fit ended below the logit", na.rm = TRUE))
    warning("For ", paste(fits$alternative[fits$note == "the inflated fit ended below the logit"], collapse = ", "), " the inflated fit's log-likelihood is below the logit's; its search did not reach the logit at the boundary.", call. = FALSE)
  stat <- ifelse(fits$status == "failed", NA_real_, ifelse(fits$status == "logit", 0, raw))
  pv <- ifelse(fits$status %in% c("failed", "degenerate"), NA_real_, vapply(stat, function(s) if (is.na(s)) NA_real_ else .pc_chibar_p(s, 0L), numeric(1)))
  out <- data.frame(test = paste0("captive: ", fits$alternative), statistic = stat, df = ifelse(is.na(stat), NA_integer_, 1L),
                    p.value = pv, share = fits$captive_share, stringsAsFactors = FALSE)
  attr(out, "method") <- "Likelihood-ratio test of the multinomial logit against a constant captive share of each alternative (the inflated multinomial logit with an intercept-only split; p-value from the 50:50 mixture of chi-square(0) and chi-square(1))"
  attr(out, "fits") <- fits
  if (any(nzchar(fits$note))) attr(out, "note") <- paste(paste0(fits$alternative, ": ", fits$note)[nzchar(fits$note)], collapse = "; ")
  class(out) <- c("iia_test", "data.frame")
  out
}

#' Parametric-bootstrap p-values for a specification test
#'
#' Calibrates any test in the package (or a user function returning a data
#' frame with a `statistic` column) by the parametric bootstrap: `B` outcome
#' vectors are drawn from the fitted model with the covariates held fixed, the
#' model is refit on each (with the fit's own covariance type; a bootstrap
#' covariance is replaced by the analytic one inside the replicates), the test
#' is recomputed, and the p-value of every row is \eqn{(1 + k) / (B' + 1)},
#' with \eqn{B'} the replicates that succeeded and \eqn{k} the number of
#' their statistics at least as large as the observed one. This is the
#' reference distribution that remains valid when the asymptotic chi-square is
#' unreliable (small samples, rare alternatives, Hausman statistics with a
#' non-positive-definite covariance difference).
#'
#' Frequency weights (integers) are expanded to one row per case, for the
#' simulated data and for the observed statistic alike, so that both follow
#' one definition; a test that refits from the data rather than from the
#' fitted design cannot be calibrated with weights (the package's own tests all
#' work from the fitted design; a user-written test that refits with [mnl()]
#' from `object$data` is one that refits from the data), and a test that calls
#' [update()] on the fit re-evaluates the original call on the original data,
#' so it must work from the fitted object (its `data` component carries the
#' simulated response). Replicates whose refit
#' warns (a separated simulated sample, an ill-conditioned information matrix)
#' are kept and counted (`n_flagged`): a few are expected with rare
#' alternatives; when they are many, the statistic itself is unstable on data
#' like these, and the p-value should be read with that in mind. The
#' attribute `p.boot_unflagged` leaves the flagged replicates out. It is a
#' p-value conditional on unflagged samples, whose size has not been studied,
#' and is no substitute for `p.boot`. A replicate in which an alternative is
#' never chosen, and whose maximum-likelihood estimate therefore does not
#' exist, is dropped and counted in `n_failed`, so the bootstrap distribution
#' is conditional on the estimate existing. For [mnl()] without fixed
#' effects, [separation()] decides that (an alternative without a constant
#' can be never chosen and still have an estimate); for the other models an
#' alternative with coefficients of its own must be chosen
#' ([boot_lr_test()] keeps such replicates where the separation check
#' applies). For a penalized fit, whose estimates exist without every
#' alternative chosen, such replicates are kept and refit under the same
#' penalty.
#'
#' The replicates are drawn from the fitted model, so the test's null
#' hypothesis must be that model: the IIA tests, [nested_lm_test()],
#' [captivity_test()], and user-written tests of that kind. [term_test()] and
#' [combine_test()] test a restriction of the fitted model, which the
#' replicates would not satisfy, and are refused; for a term, fit the model
#' without it and use [boot_lr_test()]. The package's tests derived for
#' maximum-likelihood estimates refuse a penalized fit (`firth = TRUE` or a
#' prior), so they cannot be calibrated on one here.
#'
#' Every replicate's outcomes are drawn in the session before any refit runs,
#' so a result under [set.seed()] does not depend on `cores` (for a test with
#' random draws of its own, such as [small_hsiao()], when its `seed` is
#' fixed). The replicates draw cases independently: with clustered data the
#' bootstrap distribution ignores the clustering.
#'
#' @param object A fitted [mnl()] model.
#' @param test A test function taking the fitted [mnl()] model as its first
#'   argument, whose null hypothesis is that model, and returning a data frame
#'   with a `statistic` column (e.g. [hausman_mcfadden()], [small_hsiao()],
#'   [fp_tests()], [nested_lm_test()], [captivity_test()]; not [hl_test()],
#'   which returns a list, and not [term_test()] or [combine_test()], whose
#'   null is a restriction of the model).
#' @param B Number of bootstrap replicates (default 199).
#' @param cores Number of parallel workers.
#' @param ... Further arguments to `test`.
#' @return The test's data frame with the column `p.boot` added (and
#'   `p.value` kept), the attribute `"boot_statistics"` holding the \eqn{B}
#'   replicate statistics per row, `"B"`, `"n_failed"` (replicates whose refit
#'   or test failed, dropped), `"n_flagged"` (replicates whose refit raised
#'   a warning; they are kept, and the count is printed with the method line),
#'   `"flagged"` (which replicates, aligned with the columns of
#'   `"boot_statistics"`), `"flag_reasons"` (a table of the first warning of
#'   each flagged refit), and, when some but not all replicates are flagged,
#'   `"p.boot_unflagged"`, the bootstrap p-values over the unflagged replicates.
#'   When the observed fit raised no such warning, those are the replicates that
#'   resemble the data: a separated simulated sample can produce an extreme
#'   statistic that says nothing about the data at hand.
#' @examples
#' \donttest{
#' set.seed(3)
#' d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
#' m <- mnl(y ~ x1 + x2, data = d)
#' boot_pvalue(m, hausman_mcfadden, B = 99, drop = "b")
#' }
#' @seealso [hausman_mcfadden()], [small_hsiao()], [fp_tests()], [nested_lm_test()]
#' @family IIA tests
#' @export
boot_pvalue <- function(object, test, B = 199L, cores = 1L, ...) {
  .pc_check_mnl(object, "boot_pvalue()")
  ## the replicates are drawn from the fit, which must be the test's null
  ## (round 7, B-1): the null of term_test() and combine_test() is a restriction
  ## of the fit, and replicates from the fit would not satisfy it
  if (identical(test, term_test) || identical(test, combine_test))
    stop("boot_pvalue() draws its replicates from the fitted model, the null hypothesis of the IIA and score tests. The null of ",
         if (identical(test, term_test)) "term_test()" else "combine_test()", " is a restriction of that model, which the replicates would not satisfy, ",
         "so the bootstrap p-value would be wrong. For a term, fit the model without it and use boot_lr_test(restricted, full).", call. = FALSE)
  B <- .pc_check_whole(B, "B", 10)
  obs <- test(object, ...)
  if (!is.data.frame(obs) || !"statistic" %in% names(obs)) stop("'test' must return a data frame with a 'statistic' column.", call. = FALSE)
  D <- object$design
  ## the replicates are refit with the fit's own covariance type (a bootstrap
  ## covariance is replaced by the analytic one inside the replicates)
  se_b <- if (object$se %in% c("analytic", "opg", "robust", "cluster")) object$se else "analytic"
  ## frequency weights: a row with weight w is w cases, so the simulated data
  ## expand the design to one row per case (integer weights only), and the
  ## observed statistic is recomputed on the same expanded design so that both
  ## follow one definition (a split or a grouping of rows would otherwise differ)
  if (any(D$w != 1)) {
    if (any(abs(D$w - round(D$w)) > 1e-8)) stop("boot_pvalue() needs integer frequency weights (every row is replicated 'weight' times in the simulated data).", call. = FALSE)
    D$w[] <- 1
    D <- .pc_subset_design(D, rep(seq_len(D$n), round(object$design$w)))
    oe <- tryCatch(suppressWarnings(suppressMessages(.pc_fit_design(D, se = se_b, start = object$theta, call = object$call))), error = function(e) { .pc_pass_limit(e); NULL })
    if (is.null(oe)) stop("The fit could not be reproduced on the expanded (one row per case) design.", call. = FALSE)
    obs <- test(oe, ...)
  }
  P <- .pc_fns(D)$probs(object$theta)$prob
  cp <- t(apply(P / rowSums(P), 1L, cumsum))
  ## every outcome vector is drawn here, before any test runs: a test that seeds
  ## the random-number stream (small_hsiao(seed = ) through ...) would otherwise
  ## make the replicates dependent
  U <- matrix(stats::runif(B * D$n), B, D$n)
  lhs <- all.vars(attr(D$F, "lhs")[[1L]])
  one <- function(b) {
    yb <- rowSums(U[b, ] > cp) + 1L; yb[yb > D$J] <- D$J
    Db <- D; Db$y <- yb
    ## the data carried by the refit show the simulated response too (a test that
    ## refits from the data, such as a user-written one, must see it)
    if (!is.null(Db$data) && length(lhs) == 1L && lhs %in% names(Db$data)) {
      if (D$layout == "wide") Db$data[[lhs]] <- factor(D$levels[yb], levels = D$levels) else Db$data[[lhs]] <- as.integer(D$ji == yb[D$ci])
    }
    ## a replicate without a maximum-likelihood estimate is dropped (the
    ## conditional null); in the exact check's scope an unchosen alternative
    ## drops it only when the replicate is separated (finding 7)
    if (.pc_unchosen(D, tabulate(yb, nbins = D$J)) &&
        (!.pc_sep_scope(D) || isTRUE(tryCatch(.pc_separation(Db)$separated, error = function(e) { .pc_pass_limit(e); TRUE })))) return(NULL)
    ## the refit's own warnings (a separated simulated sample, an ill-conditioned
    ## information matrix) are counted, not shown: they describe the replicate
    flagged <- FALSE; why <- NA_character_
    ob <- tryCatch(withCallingHandlers(.pc_fit_design(Db, se = se_b, call = object$call),
                                       warning = function(w) { if (!flagged) why <<- conditionMessage(w); flagged <<- TRUE; invokeRestart("muffleWarning") },
                                       message = function(m) invokeRestart("muffleMessage")),
                   error = function(e) { .pc_pass_limit(e); NULL })
    if (is.null(ob)) return(NULL)
    r <- tryCatch(suppressWarnings(suppressMessages(test(ob, ...))), error = function(e) { .pc_pass_limit(e); NULL })
    if (is.null(r) || nrow(r) != nrow(obs)) NULL else list(stat = r$statistic, flagged = flagged, why = why)
  }
  res <- .pc_lapply(seq_len(B), one, cores = cores)
  ok <- Filter(Negate(is.null), res)
  if (length(ok) < B / 2) stop("Most bootstrap replicates failed; the test or the model cannot be recomputed on simulated data.", call. = FALSE)
  Bm <- do.call(cbind, lapply(ok, `[[`, "stat"))
  n_flagged <- sum(vapply(ok, `[[`, logical(1), "flagged"))
  why <- vapply(ok, `[[`, character(1), "why"); why <- why[!is.na(why)]
  ## the first warning of each flagged refit, up to the end of its parenthesis
  ## (or 120 characters), with standalone numbers masked so one cause is one row
  if (length(why)) attr(obs, "flag_reasons") <- sort(table(substr(gsub("(?<![A-Za-z0-9_.:])[-+]?[0-9][0-9.e+-]*", "#", sub("\\): .*", ")", why), perl = TRUE), 1L, 120L)), decreasing = TRUE)
  flg <- vapply(ok, `[[`, logical(1), "flagged")
  ## a replicate whose statistic is NA for a row (a Wald contrast with a
  ## singular covariance, say) is left out of that row's p-value only (local
  ## review 2026-11-06, finding 9); an infinite one counts
  pb <- function(cols) vapply(seq_len(nrow(obs)), function(r) { v <- Bm[r, cols]; v <- v[!is.na(v)]; (1 + sum(v >= obs$statistic[r])) / (length(v) + 1) }, numeric(1))
  obs$p.boot <- pb(rep(TRUE, ncol(Bm)))
  n_na <- rowSums(is.na(Bm))
  if (any(n_na > 0L)) attr(obs, "n_na") <- stats::setNames(n_na, if (!is.null(obs$test)) obs$test else NULL)
  if (n_flagged > 0L && n_flagged < ncol(Bm)) attr(obs, "p.boot_unflagged") <- pb(!flg)
  attr(obs, "boot_statistics") <- Bm; attr(obs, "flagged") <- flg; attr(obs, "B") <- ncol(Bm); attr(obs, "n_failed") <- B - length(ok); attr(obs, "n_flagged") <- n_flagged
  attr(obs, "method") <- paste0(attr(obs, "method"), " -- parametric bootstrap, ", ncol(Bm), " replicates",
                                if (B - length(ok) > 0L) paste0(", ", B - length(ok), " failed") else "",
                                if (n_flagged > 0L) paste0(", ", n_flagged, " with a warning in the refit (attribute \"flag_reasons\")") else "",
                                if (any(n_na > 0L)) paste0(", up to ", max(n_na), " with an NA statistic left out of a row (attribute \"n_na\")") else "")
  obs
}

## the rankings of a bootstrap replicate (DESIGN.md Phase 4c): utilities U
## (n x J) order each case's available alternatives, and the order is cut to
## the case's observed rank groups `grp[[i]]` (their sizes, in order): the
## first grp[[i]][1] alternatives get rank 1, the next grp[[i]][2] rank 2, and
## so on; the rest are unranked (0)
.pc_ro_coarsen <- function(U, avail, grp) {
  U[avail == 0L] <- -Inf
  R <- matrix(0L, nrow(U), ncol(U))
  for (i in seq_len(nrow(U))) {
    g <- grp[[i]]
    if (length(g)) R[i, order(-U[i, ])[seq_len(sum(g))]] <- rep.int(seq_along(g), g)
  }
  R
}

#' Parametric-bootstrap likelihood-ratio (and Vuong) test of nested fits
#'
#' Calibrates the comparison of a restricted fit with a fit that nests it by
#' simulation under the restricted model: `B` outcome vectors are drawn from the
#' restricted fit, both models are refit on each, and the bootstrap
#' distribution of the likelihood-ratio statistic (and, for an [mnl()] versus
#' [imnl()] pair, of the raw Vuong statistic) gives the p-values. This is the
#' calibrated comparison for restrictions on the boundary of the parameter
#' space, where neither the chi-square reference of [lr_test()] nor the normal
#' reference of [vuong()] holds: the inflated model against the multinomial
#' logit (\eqn{\phi = 1}, where the split coefficients are not identified),
#' the random-intercept model against the logit (\eqn{\sigma = 0}). It applies
#' to the pairs [lr_test()] accepts (a subset of the parameters within one
#' family, with the nest structure shared for the nested logit, or the logit
#' with the same covariates inside a nested logit, an HEV, or a
#' random-intercept fit), to the logit inside the inflated multinomial logit,
#' and to an inflated multinomial logit inside one with more inflated
#' alternatives; and to two nested [rologit()] fits of the same rankings,
#' where few rankings per coefficient make the chi-square reference of
#' [lr_test()] unreliable. Penalized fits (`firth = TRUE` or a prior) are
#' refused.
#'
#' @details
#' Each replicate refits both models the way the original fits were obtained
#' (the engine's own starting values and, for the inflated model, its
#' multistart), and in addition from the original estimates and, for the full
#' model, from the restricted replicate estimates embedded in the full model at
#' the restriction (zero for the added coefficients, \eqn{\lambda = 1},
#' scales 1, a standard deviation near zero, or the split equation at
#' \eqn{\phi = 1 - 10^{-8}}). The full refit therefore cannot end below the
#' restricted one, and the replicate statistic is maximized as thoroughly as
#' the observed one, which is computed the same way. A replicate takes about
#' as long as the two original fits; `cores` runs replicates in parallel. The
#' outcomes, and the seeds of the ranking draws, come from the session's
#' generator before the refits run, so a result under [set.seed()] does not
#' depend on `cores`.
#'
#' The outcomes are drawn from the restricted fit's probabilities, case by
#' case, or for a random-intercept restricted fit unit by unit (one intercept
#' vector per unit, then the unit's cases). Rankings are drawn as the order of
#' the restricted fit's utilities plus independent standard Gumbel errors, cut
#' to each case's observed pattern: the same rank groups, with the same sizes in
#' the same order (a tie stays a tie at its positions), and the alternatives
#' below them unranked; the truncation and the ties are thus treated as
#' independent of the utilities.
#'
#' Under separation (in [mnl()] and [rologit()] fits without unit fixed
#' effects or inflation) the statistic compares the suprema of the two
#' log-likelihoods, not the values where the optimizer stopped (see
#' [lr_test()]); replicates whose refit is separated are kept and counted,
#' since separation can be common under the null when the rankings or choices
#' are few or the attributes discrete (a replicate in which an alternative is
#' never chosen, or never ranked
#' above another, is such a replicate). When the restricted fit is itself
#' separated, the outcomes are drawn from the point at which its supremum is
#' approached, where the separated comparisons are near-deterministic (the
#' restricted model's limit), not from where its optimizer stopped. Outside
#' that scope a replicate in which an alternative with coefficients of its own
#' is never chosen fails. Replicates whose refit fails are
#' dropped and counted; the p-value is \eqn{(1 + k) / (B' + 1)} with \eqn{B'}
#' the successful replicates and \eqn{k} the number of replicate statistics at
#' least as large as the observed one. Frequency weights are honoured by
#' expanding the design to one row per case (integer weights only; not with
#' random intercepts, where a weighted unit stands for several units with
#' their own intercepts). Fixed-effects fits are refused (the unit constants
#' would be re-estimated on simulated outcomes that may make them
#' unidentified).
#'
#' @param restricted,full Two `"choicemodel"` fits to the same cases, the
#'   second nesting the first.
#' @param B Number of bootstrap replicates (default 199).
#' @param cores Number of parallel workers.
#' @return A data frame of class `"iia_test"` with one row per statistic
#'   (`LR`, and `Vuong` for an mnl / imnl pair): the observed statistic,
#'   `p.value`, and `p.boot`. For the LR, `p.value` is the chi-square p-value
#'   with `df` degrees of freedom for an interior restriction, the
#'   (conservative) chi-square p-value for \eqn{\sigma = 0}, the 50:50
#'   chi-square mixture's for the mnl / imnl pair with one inflated
#'   alternative and an intercept-only split (with [lr_test()]'s rule: `NA`
#'   when the inflated fit is degenerate, its outcome coefficients diverging
#'   or at a saddle, as the method line then says), and `NA` for any other
#'   mnl / imnl pair, which has no valid asymptotic reference; for
#'   the Vuong statistic `p.value` is `NA`, and `p.boot` is one-sided
#'   in favour of the full model. The attribute `"boot_statistics"` holds
#'   the replicates; `n_failed` (dropped), `n_nonconverged` (replicates in
#'   which either refit did not converge; kept), and `n_negative` (replicates
#'   whose full refit fell short of the restricted one by more than the
#'   embedding's own gap at a boundary restriction; set to 0 and kept) are
#'   attributes too, and the method line reports each that is not zero. A few
#'   non-converged replicates are common for the inflated model on data
#'   simulated without inflation, where the split equation runs towards
#'   \eqn{\phi = 1} and the optimizer stops short of a maximum that is not
#'   attained; `n_negative` guards against an optimizer that ends above its own
#'   start and should be 0. The outcomes are drawn independently case by case
#'   (unit by unit for a random-intercept restricted fit), so the calibration
#'   assumes independent cases; with clustered data the bootstrap p-value is not
#'   calibrated. The split equation must contain an intercept named
#'   `(Intercept)`; a constant column under another name is not recognized.
#' @examples
#' \donttest{
#' set.seed(2)
#' d <- rmnl(400, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
#' m0 <- mnl(y ~ x1 + x2, data = d)
#' m1 <- imnl(y ~ x1 + x2 | 1, data = d)           # a constant captive share; none in the data
#' boot_lr_test(m0, m1, B = 19)
#' }
#' @seealso [lr_test()], [vuong()], [boot_pvalue()]
#' @family model comparison
#' @export
boot_lr_test <- function(restricted, full, B = 199L, cores = 1L) {
  if (!inherits(restricted, "choicemodel") || !inherits(full, "choicemodel")) stop("Both arguments must be fitted polychoice models.", call. = FALSE)
  .pc_no_firth(restricted, "boot_lr_test()"); .pc_no_firth(full, "boot_lr_test()")
  if (!is.null(restricted$fe) || !is.null(full$fe)) stop("boot_lr_test() does not apply to fixed-effects fits.", call. = FALSE)
  ro <- identical(restricted$family, "rologit")
  if (ro != identical(full$family, "rologit")) stop("boot_lr_test() compares two rologit() fits of the same rankings, or two single-choice fits.", call. = FALSE)
  kind <- .pc_nesting(restricted, full, allow_inflation = TRUE)
  if (full$df <= restricted$df) stop("'full' must have more parameters than 'restricted'.", call. = FALSE)
  if (kind == "inflation" && !"(Intercept)" %in% colnames(full$design$Z))
    stop("The inflated model's split equation has no intercept, so phi cannot approach 1 for every case and the model does not nest the multinomial logit; add the intercept.", call. = FALSE)
  B <- .pc_check_whole(B, "B", 10)
  Dr <- restricted$design; Df <- full$design
  if (any(Dr$w != 1)) {
    if (identical(restricted$family, "remnl") || identical(full$family, "remnl"))
      stop("boot_lr_test() does not take frequency weights with random intercepts (a weighted unit stands for several units, each with its own intercepts); refit on the expanded data.", call. = FALSE)
    if (any(abs(Dr$w - round(Dr$w)) > 1e-8)) stop("boot_lr_test() needs integer frequency weights.", call. = FALSE)
    rows <- rep(seq_len(Dr$n), round(Dr$w))
    Dr$w[] <- 1; Df$w[] <- 1
    Dr <- .pc_subset_design(Dr, rows); Df <- .pc_subset_design(Df, rows)
  }
  ddf <- full$df - restricted$df
  ## the observed statistics, with the full fit also tried from the restricted
  ## estimates embedded at the restriction (as in every replicate)
  emb <- .pc_embed(restricted$theta, Dr, Df)
  op <- tryCatch(suppressWarnings(.pc_optimize(.pc_fns(Df), emb, .pc_scales(Df), 1000L, 1e-10)), error = function(e) { .pc_pass_limit(e); NULL })
  ll_full <- full$loglik; th_full <- full$theta
  if (!is.null(op) && is.finite(op$value) && -op$value > ll_full + 1e-6) {
    message("The full model's maximum from the restricted estimates (", format(-op$value, digits = 8), ") is above the supplied fit's (",
            format(full$loglik, digits = 8), "); the higher one is used.")
    ll_full <- -op$value; th_full <- op$par
  }
  ## suprema under separation (DESIGN.md Phase 4c)
  sp_obs <- .pc_sup_pair(Dr, restricted$theta, restricted$loglik, Df, th_full, ll_full)
  lr_obs <- max(2 * (sp_obs$f - sp_obs$r), 0)
  obs <- data.frame(test = "LR (parametric bootstrap under the restricted model)", statistic = lr_obs, df = ddf,
                    p.value = if (kind == "inflation") NA_real_ else stats::pchisq(lr_obs, ddf, lower.tail = FALSE), stringsAsFactors = FALSE)
  ## one inflated alternative with an intercept-only split: the captive odds
  ## at zero are a regular boundary point, and the 50:50 chi-square mixture is
  ## the asymptotic reference, as in lr_test() (DESIGN.md Phase 4e)
  ## (with lr_test()'s rule: the tolerance, a collapsed share, and NA for a
  ## degenerate full fit; the flags are those of the point used, round 7)
  captive <- kind == "inflation" && !restricted$inflated && length(full$inflate) == 1L && identical(colnames(full$design$Z), "(Intercept)")
  captive_degenerate <- FALSE
  if (captive) {
    ## the flags of the point used: the engine's (a fit started there), not a
    ## hand-built copy of its rules (local review 2026-11-06, finding 12)
    fl <- if (identical(th_full, full$theta)) full else
      tryCatch(withCallingHandlers(.pc_fit_design(Df, se = "none", start = th_full),
                                   message = function(m) invokeRestart("muffleMessage"), warning = function(w) invokeRestart("muffleWarning")),
               error = function(e) { .pc_pass_limit(e); full })
    cr <- .pc_captive_rule(2 * (sp_obs$f - sp_obs$r), fl, restricted$loglik, ddf - 1L)
    captive_degenerate <- cr$degenerate
    obs$statistic[1L] <- cr$statistic                 # the statistic the p-values use
    obs$p.value[1L] <- if (cr$degenerate) NA_real_ else .pc_chibar_p(cr$statistic, ddf - 1L)
  }
  if (kind == "inflation") {
    l_full <- .pc_fns(Df)$lli(th_full); l_res <- .pc_fns(Dr)$lli(restricted$theta)
    v_obs <- tryCatch(unname(.pc_vuong_stat(l_full, l_res, Dr$w, ddf)$stat["raw"]), error = function(e) { .pc_pass_limit(e); NA_real_ })
    obs <- rbind(obs, data.frame(test = "Vuong, raw (full vs restricted; parametric bootstrap)", statistic = v_obs, df = NA_integer_,
                                 p.value = NA_real_, stringsAsFactors = FALSE))
  }
  ## the point the outcomes are drawn from: the restricted fit's estimates, or,
  ## when the restricted design is separated, the point its supremum is
  ## approached at, where the separated comparisons are near-deterministic (the
  ## limit of the restricted model; DESIGN.md Phase 4d). Where the optimizer
  ## stopped along the separating direction depends on its tolerance.
  th_draw <- if (isTRUE(sp_obs$r_separated)) sp_obs$r_theta else restricted$theta
  in_scope <- .pc_sep_scope(Dr) && .pc_sep_scope(Df)
  ## the draws: uniforms for the outcomes and, for a random-intercept
  ## restricted fit, standard normals for the unit intercepts, all drawn here
  ## so that the result does not depend on the number of workers
  Pr <- if (ro) NULL else .pc_fns(Dr)$probs(th_draw)
  U <- if (ro) NULL else matrix(stats::runif(B * Dr$n), B, Dr$n)
  ## rankings (DESIGN.md Phase 4c): the restricted utilities plus standard
  ## Gumbel errors order each case's available alternatives, and the order is
  ## coarsened to the case's observed pattern (its rank groups' sizes, in
  ## order; the alternatives below them unranked). One seed per replicate,
  ## drawn here, and the session's generator kinds, so that the result does not
  ## depend on the number of workers (a worker's stream is L'Ecuyer's)
  if (ro) {
    Vr <- pc_ro_probs_cpp(th_draw, Dr$X, Dr$W, Dr$off, Dr$avail, Dr$J, Dr$base - 1L)$V
    grp <- lapply(seq_len(Dr$n), function(i) { r <- Dr$str$rank[i, Dr$avail[i, ] == 1L]; tabulate(r[r > 0L]) })
    seeds <- sample.int(.Machine$integer.max, B)
  }
  rk <- RNGkind()
  draw_rank <- function(b) {
    old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL
    set.seed(seeds[b], kind = rk[1L], normal.kind = rk[2L], sample.kind = rk[3L])
    Gm <- matrix(-log(-log(stats::runif(Dr$n * Dr$J))), Dr$n, Dr$J)
    if (is.null(old)) rm(".Random.seed", envir = globalenv()) else assign(".Random.seed", old, envir = globalenv())
    .pc_ro_coarsen(Vr + Gm, Dr$avail, grp)
  }
  re <- identical(restricted$family, "remnl")
  if (re) { G <- nlevels(Dr$unit); m <- length(Dr$idx$str); Zu <- matrix(stats::rnorm(B * G * m), B, G * m) }
  draw_y <- function(b) {
    P <- if (re) {
      sig <- exp(restricted$theta[Dr$idx$str]); ui <- as.integer(Dr$unit)
      V <- Pr$V; nb <- setdiff(seq_len(Dr$J), Dr$base)
      V[, nb] <- V[, nb] + (matrix(Zu[b, ], G, m) * rep(sig, each = G))[ui, , drop = FALSE]
      V[Dr$avail == 0L] <- -Inf
      E <- exp(V - apply(V, 1L, max)); E / rowSums(E)
    } else Pr$prob
    cp <- t(apply(P / rowSums(P), 1L, cumsum))
    yb <- rowSums(U[b, ] > cp) + 1L; yb[yb > Dr$J] <- Dr$J
    yb
  }
  one <- function(b) {
    if (ro) {
      ## a ranking replicate: the first choice (the first tied alternative in
      ## level order) is the design's y; an alternative never ranked above
      ## another is separation, kept (its supremum enters the statistic)
      Rb <- draw_rank(b)
      yb <- apply(Rb, 1L, function(r) which(r == 1L)[1L])
      Drb <- Dr; Drb$str$rank <- Rb; Drb$y <- yb; Dfb <- Df; Dfb$str$rank <- Rb; Dfb$y <- yb
    } else {
      yb <- draw_y(b)
      cb <- tabulate(yb, nbins = Dr$J)
      ## an alternative never chosen is separation in the separation check's
      ## scope, kept as in rankings (its supremum enters the statistic; DESIGN.md
      ## Phase 4d); elsewhere it leaves a parameter without an estimate and the
      ## replicate fails
      if (!in_scope && (.pc_unchosen(Dr, cb) || .pc_unchosen(Df, cb))) return("an alternative is never chosen")
      Drb <- Dr; Drb$y <- yb; Dfb <- Df; Dfb$y <- yb
    }
    fr <- .pc_best_fit(Drb, extra = unique(list(restricted$theta, th_draw)))
    if (is.null(fr)) return("the restricted refit failed")
    eb <- .pc_embed(fr$par, Drb, Dfb)
    ## at a boundary restriction (sd = 0, phi = 1) the embedding is a point near
    ## the limit, below the restricted maximum by a small gap; a full refit
    ## below the restricted one by no more than that gap is not a failure
    gap <- if (kind == "interior") 0 else max(0, tryCatch(.pc_fns(Dfb)$nll(eb), error = function(e) { .pc_pass_limit(e); Inf }) - fr$value)
    ff <- .pc_best_fit(Dfb, extra = list(eb, th_full))
    if (is.null(ff)) return("the full refit failed")
    spb <- .pc_sup_pair(Drb, fr$par, -fr$value, Dfb, ff$par, -ff$value)
    lr <- 2 * (spb$f - spb$r)
    ## the captive share: the observed statistic's tolerance (.pc_captive_rule)
    lrs <- if (captive && lr < max(1e-6, 2e-8 * max(1, abs(fr$value)))) 0 else max(lr, 0)
    out <- c(lr = lrs, negative = lr < -(1e-6 + 2 * gap), converged = fr$converged && ff$converged, separated = spb$used)
    if (kind == "inflation")
      out <- c(out, vuong = tryCatch(unname(.pc_vuong_stat(.pc_fns(Dfb)$lli(ff$par), .pc_fns(Drb)$lli(fr$par), Drb$w, ddf)$stat["raw"]),
                                     error = function(e) { .pc_pass_limit(e); NA_real_ }))
    out
  }
  res <- .pc_lapply(seq_len(B), one, cores = cores)
  why <- unlist(Filter(is.character, res))
  ok <- Filter(Negate(is.character), res)
  if (length(ok) < B / 2) stop("Most bootstrap replicates failed", if (length(why)) paste0(" (most often: ", names(sort(table(why), decreasing = TRUE))[1L], ")") else "",
                               "; the models cannot be refit on simulated data.", call. = FALSE)
  Bm <- do.call(cbind, ok)
  stat_rows <- c("lr", if (kind == "inflation") "vuong")
  obs$p.boot <- vapply(seq_along(stat_rows), function(r) {
    s <- Bm[stat_rows[r], ]; s <- s[is.finite(s)]
    if (!is.finite(obs$statistic[r]) || !length(s)) NA_real_ else (1 + sum(s >= obs$statistic[r])) / (length(s) + 1)
  }, numeric(1))
  ## a separated refit cannot meet the optimizer's tolerance along the separated
  ## direction; its supremum enters the statistic, so it is counted apart
  n_failed <- B - length(ok); n_neg <- sum(Bm["negative", ]); n_sep <- sum(Bm["separated", ] == 1)
  n_nc <- sum(!Bm["converged", ] & Bm["separated", ] == 0)
  attr(obs, "boot_statistics") <- Bm[stat_rows, , drop = FALSE]; attr(obs, "B") <- ncol(Bm)
  attr(obs, "n_failed") <- n_failed; attr(obs, "n_nonconverged") <- n_nc; attr(obs, "n_negative") <- n_neg; attr(obs, "n_separated") <- n_sep
  attr(obs, "method") <- paste0("Parametric bootstrap under the restricted model, ", ncol(Bm), " replicates",
                                if (ro) " (rankings drawn with each case's observed rank groups)" else "",
                                if (n_failed > 0L) paste0(", ", n_failed, " failed") else "",
                                if (n_sep > 0L) paste0(", ", n_sep, " with a separated refit (its log-likelihood supremum used)") else "",
                                if (n_nc > 0L) paste0(", ", n_nc, " with a refit that did not converge") else "",
                                if (sp_obs$used) "; the observed statistic uses log-likelihood suprema (the data are separated)" else "",
                                if (isTRUE(sp_obs$r_separated)) "; the draws come from the separated restricted model's limit (the point its supremum is approached at)" else "",
                                if (n_neg > 0L) paste0(", ", n_neg, " with a full refit below the restricted one (set to 0)") else "",
                                if (kind == "sd") "; the chi-square p-value is conservative (sd = 0 is on the boundary)" else "",
                                if (captive && captive_degenerate) "; the full fit is degenerate (outcome coefficients diverging, or a saddle), so the LR has no asymptotic reference (p.value NA)"
                                else if (captive) "; the LR's asymptotic p-value is the 50:50 chi-square mixture (a constant captive share)"
                                else if (kind == "inflation") "; no valid asymptotic reference for the LR (phi = 1 is on the boundary)" else "")
  class(obs) <- c("iia_test", "data.frame")
  obs
}
