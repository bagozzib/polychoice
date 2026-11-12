## ---------------------------------------------------------------------------
## Split-panel (half-panel) jackknife for the slopes of a fixed-effects fit
## (Dhaene and Jochmans 2015), with the assumption checks the method requires.
## Mirrors underdisp:::.fe_jackknife, adapted to a vector of slopes and to the
## multinomial outcome (DESIGN.md, phase 2).
##
## The maximum-likelihood unit constants carry an incidental-parameters bias of
## order 1/T that leaks into the slopes. The correction refits the model on the
## first and second halves of every unit's cases (rows split in the order
## supplied, which should be temporal order) and reports
##   theta_jack = 2 * theta_full - (theta_half1 + theta_half2) / 2
## for the slopes (every parameter that is not a unit dummy or an
## alternative-specific constant; the constants refer to a reference unit that
## may differ between halves and are nuisance here).
##
## VALIDITY GATE. The identity requires the two half-panels to estimate the same
## pseudo-true parameter, differing only through the 1/T term (time
## homogeneity). Two detectors:
##   * asymmetric heterogeneity (a regime change): the half-panel slope vectors
##     disagree. Placebo-calibrated: five random splits of every unit's cases
##     into halves of the time split's sizes, ignoring the time order, have
##     the time split's structure (all units, half-panels of the same length),
##     so under time homogeneity their distances and the time split's share a
##     distribution; their mean Euclidean distance measures that noise, and the
##     time-split distance must not exceed kappa = 2.5 times it (calibrated by
##     verification/jackknife_gate_mc.R under the ledger's selection rule). A
##     regime change at mid-panel separates the time halves but is balanced in
##     every random split. The splits are drawn from fixed seeds with the
##     caller's random-number stream restored, so the gate is reproducible;
##   * symmetric heterogeneity (a smooth trend): both halves differ from the
##     full fit the same way, which the first detector cannot see. Detected on
##     the full fit: the residuals y_ij - P_ij of some alternative correlate with
##     the within-unit time order (Bonferroni over alternatives at 0.01). The
##     cure is a time control in the formula, which the message says.
## The gate REFUSES the correction (the caller returns the ML fit with a
## warning) on either detector or on a failed refit; it WARNS when the median
## panel length is below 6 or the implied correction is large.
## ---------------------------------------------------------------------------

## the design of a subset of cases with the fixed-effect dummies regenerated on
## the subset (never-chosen pairs and empty dummies handled as in a fresh fit)
.pc_half_design <- function(D, cases) {
  S <- .pc_subset_design(D, cases)
  if (is.null(D$fe_map) || !nrow(D$fe_map)) return(S)
  S$W <- S$W[, !colnames(S$W) %in% D$fe_map$name, drop = FALSE]
  if (!is.null(S$fe_off)) S$avail[S$fe_off] <- 1L                  # back to the choice sets in the data
  ## the fixed-effects rules and the single-alternative drop, repeated until
  ## nothing changes (as in .pc_design)
  repeat {
    S$unit <- droplevels(S$unit)
    fx <- .pc_fe_columns(S$unit, S$y, S$J, S$base, S$levels, S$n, S$avail, D$fe, asc = !isFALSE(D$asc), build = FALSE)
    single <- rowSums(fx$avail) < 2L
    if (!any(single)) break
    if (all(single)) stop("No case of the subsample has two or more available alternatives.", call. = FALSE)
    S <- .pc_subset_design(S, which(!single))
  }
  fx <- .pc_fe_columns(S$unit, S$y, S$J, S$base, S$levels, S$n, S$avail, D$fe, asc = !isFALSE(D$asc))   # the dummies, once
  S$avail <- fx$avail; S$fe_map <- fx$map; S$fe_off <- fx$off; S$W <- cbind(S$W, fx$W)
  .pc_subsample_check(.pc_reindex(S))
}

## the identification checks of .pc_design() for a subsample design (half
## panels, placebo splits, bootstrap replicates): every alternative chosen in a
## case where another is available, no column zero on every case, and the
## scale-free identification check of the choice design. A subsample that
## fails stops with the reason, which the callers report (a refused jackknife,
## a failed replicate) instead of refitting a flat direction that would stay
## at its starting value. An alternative with coefficients of its own that is
## never chosen in the subsample (in a ranking: never ranked above another
## available alternative) fails a maximum-likelihood replicate, whose estimate
## does not exist (DESIGN.md, Phase 4b: a separated resample is dropped, the
## author's decision of 2026-11-05); a penalized replicate (Firth's penalty or
## a prior) is finite and kept.
.pc_subsample_check <- function(S) {
  pen <- .pc_penalized(S)
  if (identical(S$family, "rologit") && !pen && S$p > 0L) {
    gone <- .pc_rank_wins(S$str$rank, S$avail) == 0L
    if (any(gone)) stop("in the subsample ", paste(S$levels[gone], collapse = ", "), " is ranked above no other available alternative", call. = FALSE)
  }
  gone <- if (identical(S$family, "rologit") || pen || !.pc_own_coef(S)) FALSE else tabulate(S$y[rowSums(S$avail) >= 2L], nbins = S$J) == 0L
  ## in the exact check's scope (mnl() without fixed effects) the estimate can
  ## exist although an alternative is never chosen (no constant, a covariate
  ## of both signs): the linear program decides (local review 2026-11-06,
  ## finding 7), as it does for the fit
  if (any(gone) && .pc_sep_scope(S) && !isTRUE(tryCatch(.pc_separation(S)$separated, error = function(e) { .pc_pass_limit(e); TRUE }))) gone[] <- FALSE
  if (any(gone)) stop("in the subsample ", paste(S$levels[gone], collapse = ", "), " is chosen in no case where another alternative is available", call. = FALSE)
  if (identical(S$family, "nlogit") && !.pc_own_coef(S)) {
    gn <- .pc_unchosen_nests(S$str, S$levels, tabulate(S$y[rowSums(S$avail) >= 2L], nbins = S$J))
    if (length(gn)) stop("in the subsample no alternative of nest(s) ", paste(gn, collapse = ", "), " is chosen in a case where another alternative is available", call. = FALSE)
  }
  zz <- if (!is.null(S$Z) && ncol(S$Z)) colnames(S$Z)[colSums(S$Z != 0) == 0] else character(0)
  zc <- c(colnames(S$X)[colSums(S$X != 0) == 0], colnames(S$W)[colSums(S$W != 0) == 0],
          if (length(zz)) paste0("infl_", zz))                       # (paste0 of a zero-length vector is not empty)
  if (length(zc)) stop("in the subsample the column(s) ", paste(zc, collapse = ", "), " are zero on every case", call. = FALSE)
  .pc_identification_check(S)
  S
}

## the number of cases in which each alternative is ranked above another
## available alternative (the stages it can win)
.pc_rank_wins <- function(R, avail) {
  Ra <- R * (avail == 1L)
  mx <- apply(Ra, 1L, max)
  unr <- rowSums(avail == 1L & R == 0L) > 0L
  colSums(avail == 1L & R > 0L & (R < mx | unr))
}

## returns list(par = corrected slopes, halves, noise, warn, gate) or
## list(refused = reason, gate). `full = TRUE` (the gate Monte Carlo) computes
## every statistic regardless of the early refusals; `kappa` and `alpha` are the
## gate constants (see verification/jackknife_gate_mc.R for their calibration)
.pc_jackknife <- function(fit, D, maxit = 1000L, reltol = 1e-10, kappa = 2.5, alpha = 0.01, full = FALSE) {
  iu <- as.integer(D$unit)
  ord <- stats::ave(seq_len(D$n), iu, FUN = seq_along)
  len <- stats::ave(seq_len(D$n), iu, FUN = length)
  slopes <- D$par_names[!(D$blocks %in% "fixed effect") & !grepl(":\\(Intercept\\)$", D$par_names)]
  if (!length(slopes)) return(list(refused = "the model has no slope coefficients to correct"))
  warn <- character(0)
  if (stats::median(tabulate(iu)) < 6)
    warn <- c(warn, "median panel length is below 6; half-panel fits are short and the correction may be unreliable")
  ## symmetric heterogeneity: residual time trend on the full fit, per alternative
  tt <- ord - stats::ave(ord, iu)
  R <- (outer(D$y, seq_len(D$J), "==") * 1) - fit$fitted.values
  ## the frequency-weighted correlation over the cases where j is available
  ## (elsewhere the residual is 0 by construction), tested on sum(w) - 2 df
  pv <- vapply(seq_len(D$J), function(j) {
    a <- D$avail[, j] == 1L
    if (sum(a) < 3L) return(1)
    r <- R[a, j]; t0 <- tt[a]; w <- D$w[a]; N <- sum(w)
    mr <- sum(w * r) / N; mt <- sum(w * t0) / N
    vr <- sum(w * (r - mr)^2); vt <- sum(w * (t0 - mt)^2)
    if (vr <= 0 || vt <= 0 || N <= 3) return(1)
    rho <- sum(w * (r - mr) * (t0 - mt)) / sqrt(vr * vt)
    tstat <- rho * sqrt((N - 2) / max(1 - rho^2, 1e-300))
    2 * stats::pt(-abs(tstat), N - 2) }, numeric(1))
  gate <- list(trend_p = stats::setNames(pv, D$levels), dist = NA_real_, noise = NA_real_, placebo = rep(NA_real_, 5L), size = NA_real_, corr = NA_real_)
  trend_refusal <- sprintf(
    "the full fit's residuals for alternative %s trend with the within-unit time order (p = %.2g), so the model omits time structure and the split-panel jackknife's time-homogeneity requirement fails; include a smooth time control (a trend in the period) in the formula and re-fit (period dummies are not identified within the half panels)",
    D$levels[which.min(pv)], min(pv))
  if (min(pv) < alpha / D$J && !full) return(list(refused = trend_refusal, gate = gate))
  last_err <- NULL
  refit <- function(cases) tryCatch({
    r <- .pc_refit(.pc_half_design(D, cases), maxit = maxit, reltol = reltol)
    if (!r$converged) NULL else r$coefficients[slopes]
  }, error = function(e) { .pc_pass_limit(e); last_err <<- conditionMessage(e); NULL })
  ## the temporal halves; a unit with an odd number of cases has two possible
  ## splits (floor and ceiling of T/2), both are used and averaged (Dhaene and
  ## Jochmans 2015, section 3)
  splits_h <- if (any(len %% 2 == 1)) list(ord <= floor(len / 2), ord <= ceiling(len / 2)) else list(ord <= len / 2)
  halves <- lapply(splits_h, function(h1) list(p1 = refit(which(h1)), p2 = refit(which(!h1))))
  bad <- vapply(halves, function(hh) is.null(hh$p1) || is.null(hh$p2) || anyNA(hh$p1) || anyNA(hh$p2), logical(1))
  if (any(bad)) return(list(refused = paste0("a half-panel refit failed or did not converge",
                                             if (!is.null(last_err)) paste0(" (", sub("[.]$", "", last_err), ")") else ""), gate = gate))
  p1 <- Reduce(`+`, lapply(halves, `[[`, "p1")) / length(halves)
  p2 <- Reduce(`+`, lapply(halves, `[[`, "p2")) / length(halves)
  ## asymmetric heterogeneity: the time-split distance against placebo splits,
  ## each a random half (of the time split's size) of every unit's cases
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  placebo_h1 <- lapply(seq_len(5L), function(s) {
    set.seed(7919L * s)
    stats::ave(stats::runif(D$n), iu, FUN = rank) <= floor(len / 2)
  })
  if (had_seed) assign(".Random.seed", old_seed, envir = globalenv()) else rm(".Random.seed", envir = globalenv())
  placebo <- vapply(placebo_h1, function(h1) {
    q1 <- refit(which(h1)); q2 <- refit(which(!h1))
    if (is.null(q1) || is.null(q2) || anyNA(q1) || anyNA(q2)) NA_real_ else sqrt(sum((q1 - q2)^2)) }, numeric(1))
  noise <- mean(placebo, na.rm = TRUE)
  gate$placebo <- placebo; gate$noise <- noise
  if (!is.finite(noise) || sum(is.finite(placebo)) < 3L)
    return(list(refused = "the placebo-split fits failed, so the time-homogeneity check cannot be performed", gate = gate))
  theta_full <- fit$theta[slopes]
  size <- max(1, sqrt(sum(theta_full^2)))
  dis <- mean(vapply(halves, function(hh) sqrt(sum((hh$p1 - hh$p2)^2)), numeric(1)))
  jk <- 2 * theta_full - (p1 + p2) / 2
  corr <- sqrt(sum((theta_full - (p1 + p2) / 2)^2))
  gate$dist <- dis; gate$size <- size; gate$corr <- corr
  out <- list(par = jk, halves = cbind(first = p1, second = p2), noise = noise, placebo = placebo, warn = warn, gate = gate)
  if (min(pv) < alpha / D$J) { out$refused <- trend_refusal; return(out) }
  if (dis > max(0.05 * size, kappa * noise)) {
    out$refused <- sprintf(
      "the half-panel slope estimates disagree (distance %.3f) far beyond the noise of random within-unit splits (mean %.3f); the panel appears time-heterogeneous and the split-panel jackknife is not valid here (Dhaene and Jochmans 2015)",
      dis, noise)
    return(out)
  }
  if (corr > 0.15 * size)
    out$warn <- c(warn, sprintf("the implied correction (norm %.3f) is large relative to the validated 1/T magnitudes; verify the panel is time-homogeneous before relying on it", corr))
  out
}

## apply the correction to a fitted object (called by .pc_fit)
.pc_apply_jackknife <- function(fit, D, maxit = 1000L, reltol = 1e-10) {
  jk <- .pc_jackknife(fit, D, maxit, reltol)
  if (!is.null(jk$refused)) {
    warning("bias_correct = \"jackknife\" REFUSED: ", jk$refused, ". Returning the uncorrected maximum-likelihood estimates.", call. = FALSE)
    fit$bias_correct <- "none"
    fit$jackknife_refused <- jk$refused
    fit$jackknife <- list(gate = jk$gate)
    return(fit)
  }
  for (w in jk$warn) warning("bias_correct = \"jackknife\": ", w, ".", call. = FALSE)
  fit$uncorrected <- list(coefficients = fit$coefficients, theta = fit$theta)
  fit$coefficients[names(jk$par)] <- jk$par
  fit$theta[names(jk$par)] <- jk$par
  fit$jackknife <- list(halves = jk$halves, placebo_noise = jk$noise, placebo = jk$placebo, gate = jk$gate)
  fit$bias_correct <- "jackknife"
  fit$jackknife_refused <- NULL
  fit
}
