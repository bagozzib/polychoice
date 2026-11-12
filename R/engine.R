## ---------------------------------------------------------------------------
## The shared fitting engine behind every estimator: design, starting values,
## column-scaled BFGS with the analytic gradient and an exact-Hessian Newton
## polish, analytic / robust / cluster / bootstrap covariances, one object
## class ("choicemodel") so every method is written once.
## ---------------------------------------------------------------------------

## column scales for the optimizer, one per parameter: case-specific columns by
## their (weighted) standard deviation across cases; alternative-specific
## columns by the standard deviation of their within-case deviations over the
## available alternatives (only within-case variation is identified, so the
## between-case component would inflate the scale); split columns across cases
.pc_scales <- function(D) {
  sx <- if (D$p) .pc_colscale(D$X, D$w) else numeric(0)
  sw <- numeric(0)
  if (D$q) {
    av <- as.vector(t(D$avail)) == 1L
    case <- rep(seq_len(D$n), each = D$J)[av]
    Wa <- D$W[av, , drop = FALSE]
    Wd <- Wa - apply(Wa, 2L, function(v) stats::ave(v, case))
    sw <- .pc_colscale(Wd)
  }
  sz <- if (D$r) rep(.pc_colscale(D$Z, D$w), max(1L, length(D$kc))) else numeric(0)     # one split block per inflated alternative
  c(rep(sx, D$J - 1L), sw, sz, rep(1, length(D$idx$str)))      # structural parameters: scale 1
}

## BFGS on the standardized parameters ts = theta * scale (the effect of a
## one-standard-deviation change in each column, so every coordinate has a
## comparable curvature), then a damped Newton polish with the analytic
## Hessian; returns the parameters on their natural scale. `polish = FALSE`
## runs the BFGS stage only (the screening pass of the multistart).
## The search holds the structural log-parameters (log scales, log inclusive
## values, the probit's log variances, log standard deviations) within
## +/- .pc_str_bound, the value the boundary step moves them to: past the
## boundary flag (|log| > 7) the likelihood is usually flat, and there the HEV
## integrals are costly (on Palik et al. 2025 two thirds of the search's time
## went to points past the flag, stress finding F22). A point beyond the bound
## is evaluated at the bound, with a zero gradient in that coordinate, and the
## optimizer returns the bounded point; the fit's diagnostics and standard
## errors are computed there from the likelihood itself. Where the
## log-likelihood still rises at the bound (on Palik by 0.67 or more from 8 to
## 12 on the log scale), the fit's boundary warning says so
## (.pc_rising_at_bound).
## Under hev() and mnprobit() the search (BFGS and the Newton polish) runs in
## the coordinates of .pc_search_map (fn$search), along which the likelihood's
## ridge towards an infinite error scale is straight; convergence is judged in
## the model's own coordinates, as under every other family.
.pc_str_bound <- 8
.pc_optimize <- function(fn, start, scale, maxit = 1000L, reltol = 1e-10, polish = TRUE, newton_first = FALSE) {
  ## no parameters (a restricted model without constants that drops its only
  ## term): the objective at the empty vector
  if (!length(start)) { v <- fn$nll(numeric(0)); return(list(par = numeric(0), value = v, converged = is.finite(v), max_grad = 0, bfgs_code = NA_integer_, newton_it = 0L)) }
  sr <- fn$search
  s <- .pc_std(fn, scale)
  sb <- if (is.null(sr)) s else .pc_std(list(nll = function(u) fn$nll(sr$from(u)), grad = function(u) sr$grad(u, fn$grad(sr$from(u))), bound = fn$bound), scale)
  ts0 <- start * scale
  ## a globally concave log-likelihood (the plain MNL) is solved by damped Newton
  ## directly; BFGS is the fallback when that does not converge
  if (newton_first) {
    nw <- .pc_newton(s$fs, s$gs, s$hs, ts0, 50L)
    if (nw$converged) return(list(par = s$bd(nw$ts / scale), value = nw$val, converged = TRUE, max_grad = nw$max_grad, bfgs_code = NA_integer_, newton_it = nw$it))
  }
  ## a penalized objective (fn$info, firth = TRUE): Fisher scoring, the damped
  ## Newton iteration with the model's information as its curvature (positive
  ## definite, so every step is an ascent direction), as Bull, Mak and Greenwood
  ## (2002), logistf and coxphf compute the penalized estimates. Under
  ## separation BFGS's first steps follow the steep separated direction far out,
  ## onto a flat stretch where it stalls (a ranking design where the maximum is
  ## at 4.6 and BFGS stopped at 34.7); the step cap and step-halving keep the
  ## scoring near the maximum. BFGS and the polish continue from its endpoint
  ## only when it stops short.
  if (!is.null(fn$info)) {
    hsI <- function(ts) fn$info(ts / scale) / tcrossprod(scale)
    nw <- .pc_newton(s$fs, s$gs, hsI, ts0, 200L)
    if (nw$converged) return(list(par = nw$ts / scale, value = nw$val, converged = TRUE, max_grad = nw$max_grad, bfgs_code = NA_integer_, newton_it = nw$it))
    if (is.finite(nw$val)) { ts0 <- nw$ts; start <- ts0 / scale }
  }
  ub0 <- if (is.null(sr)) ts0 else sr$to(start) * scale
  op <- tryCatch(stats::optim(ub0, sb$fs, sb$gs, method = "BFGS", control = list(maxit = maxit, reltol = reltol)),
                 error = function(e) { .pc_pass_limit(e); NULL })
  if (is.null(op) || !is.finite(op$value)) op <- list(par = ub0, value = sb$fs(ub0), convergence = 1L)
  nw <- .pc_newton(sb$fs, sb$gs, sb$hs, op$par, if (polish) 25L else 0L, val0 = op$value)
  if (!is.null(sr)) {
    ## the search's gradient weighs an alternative's coefficients by its
    ## utility's standard deviation, so its tolerance is not the model's: where
    ## the point misses the model's tolerance, Newton iterations in the model's
    ## coordinates finish the polish
    ts <- sr$from(sb$bd(nw$ts / scale)) * scale; it <- nw$it
    g <- s$gs(ts); mg <- if (all(is.finite(g))) max(abs(g)) else Inf
    nw <- if (polish && !(mg < .pc_gtol(nw$val))) .pc_newton(s$fs, s$gs, s$hs, ts, 25L, val0 = nw$val)
          else list(ts = ts, val = nw$val, it = 0L, max_grad = mg, converged = is.finite(nw$val) && mg < .pc_gtol(nw$val))
    nw$it <- nw$it + it
  }
  list(par = s$bd(nw$ts / scale), value = nw$val, converged = nw$converged, max_grad = nw$max_grad,
       bfgs_code = op$convergence, newton_it = nw$it)
}

## the convergence tolerance on the largest scaled gradient, relative to the
## objective
.pc_gtol <- function(val) 1e-5 * max(1, abs(val))

## the objective, gradient and Hessian on the standardized parameters, held
## within fn$bound when the family has one
.pc_std <- function(fn, scale) {
  bd <- if (is.null(fn$bound)) identity else fn$bound
  fs <- function(ts) fn$nll(bd(ts / scale))
  gs <- if (is.null(fn$bound)) function(ts) fn$grad(ts / scale) / scale
        else function(ts) { th <- ts / scale; tb <- bd(th); g <- fn$grad(tb) / scale; g[th != tb] <- 0; g }
  ## analytic Hessian when the family has one, otherwise central differences of
  ## the analytic gradient in the standardized coordinates. A coordinate past
  ## the bound has the row and column of a fixed parameter (a unit diagonal),
  ## so the Newton step leaves it in place and solves for the others
  hs <- function(ts) {
    H <- if (is.null(fn$hess)) .pc_hess_fd(gs, ts) else fn$hess(bd(ts / scale)) / tcrossprod(scale)
    if (!is.null(fn$bound)) { th <- ts / scale; z <- th != bd(th); H[z, ] <- 0; H[, z] <- 0; diag(H)[z] <- 1 }
    H
  }
  list(bd = bd, fs = fs, gs = gs, hs = hs)
}

## damped Newton iterations on the standardized objective: a full step when it
## improves the objective, halved up to twelve times otherwise, a gradient step
## when the Hessian direction is not a descent direction
.pc_newton <- function(fs, gs, hs, ts, maxit, val0 = NULL) {
  val <- if (is.null(val0)) fs(ts) else val0
  it <- 0L
  for (i in seq_len(maxit)) {
    g <- gs(ts)
    if (!all(is.finite(g))) break
    if (max(abs(g)) < 1e-9 * max(1, abs(val))) break
    H <- hs(ts)
    step <- tryCatch(solve(H, g), error = function(e) { .pc_pass_limit(e); NULL })
    if (is.null(step) || !all(is.finite(step)) || sum(g * step) <= 0) step <- g / max(abs(g))   # not a descent direction: a unit gradient step
    ## no more than two standardized units per iteration: along a flat direction
    ## (a log scale running off) an unbounded step would walk the parameter far
    ## past where the objective stopped changing
    smax <- max(abs(step)); if (smax > 2) step <- step * (2 / smax)
    ## accept the longest step that does not increase the objective (steps near
    ## the optimum move the parameters without changing the objective in double
    ## precision while the gradient still shrinks)
    lam <- 1; accepted <- FALSE
    for (h in seq_len(12L)) {
      tn <- ts - lam * step; vn <- fs(tn)
      if (is.finite(vn) && vn <= val + 1e-14 * max(1, abs(val))) { ts <- tn; val <- min(val, vn); accepted <- TRUE; break }
      lam <- lam / 2
    }
    it <- i
    if (!accepted) break
  }
  g <- gs(ts); mg <- if (all(is.finite(g))) max(abs(g)) else Inf
  list(ts = ts, val = val, it = it, max_grad = mg, converged = is.finite(val) && is.finite(mg) && mg < .pc_gtol(val))
}

## the structural parameters estimated on the log scale (log scales, log
## inclusive values, log standard deviations, the probit's log variances),
## whose boundary (0 or infinity) lies at minus or plus infinity
.pc_str_log <- function(D) D$idx$str[D$blocks[D$idx$str] %in% c("scale", "nest", "random intercept") | grepl("^var:", D$par_names[D$idx$str])]

## A structural log-parameter can run towards its boundary with the
## log-likelihood still rising (a degenerate nest, a scale going to zero): the
## optimizer then stalls short of the |log| > 7 boundary flag, the gradient
## along the ridge never quite vanishing, and the fit would be reported as not
## converged. The parameters beyond 3 in magnitude (a ratio above 20) are
## tried at 8 on their own side (the search's bound, .pc_str_bound: a ratio
## near 3000, past the flag) with the others re-optimized, first all together
## and, when that lowers the log-likelihood, one at a time, largest first,
## each kept with those already kept when the log-likelihood does not fall
## (DESIGN.md Phase 4d: one parameter at an interior maximum near e^4 would
## otherwise block another running off). The result is the boundary fit the
## diagnostics then flag; a genuine interior maximum at a large ratio loses
## the comparison and is kept. Under hev() and mnprobit() the move holds the
## search's coordinates (.pc_search_map), so an alternative's coefficients
## follow its scale, and the re-optimization searches in them.
.pc_boundary_step <- function(fn, D, op, scale, maxit, reltol) {
  str_log <- .pc_str_log(D)
  if (!length(str_log) || is.null(op) || !is.finite(op$value)) return(op)
  th <- op$par
  cand <- str_log[abs(th[str_log]) > 3 & abs(th[str_log]) <= 7]
  if (!length(cand)) return(op)
  r <- .pc_boundary_try(fn, op, cand, scale, maxit, reltol)
  if (!is.null(r) || length(cand) == 1L) return(if (is.null(r)) op else r)
  best <- op; kept <- integer(0)
  for (c in cand[order(-abs(th[cand]))]) {
    ## a candidate the earlier moves brought back inside (|log| <= 3) is no
    ## longer one; its bound side is taken where it now is (round 6, I-6)
    if (abs(best$par[c]) <= 3) next
    r <- .pc_boundary_try(fn, best, c(kept, c), scale, maxit, reltol)
    if (!is.null(r)) { best <- r; kept <- c(kept, c) }
  }
  best
}

## one try of the boundary step: the parameters `fixed` moved to the bound on
## their side of the point op$par, the others re-optimized; the result when
## the log-likelihood does not fall, NULL otherwise
.pc_boundary_try <- function(fn, op, fixed, scale, maxit, reltol) {
  th <- op$par
  tf <- sign(th[fixed]) * .pc_str_bound
  free <- setdiff(seq_along(th), fixed)
  sr <- fn$search
  t0 <- th; t0[fixed] <- tf
  if (!is.null(sr)) { u0 <- sr$to(th); u0[fixed] <- tf; t0 <- sr$from(u0) }
  full <- function(b) { t <- t0; t[free] <- b; t }
  fnr <- list(nll = function(b) fn$nll(full(b)), grad = function(b) fn$grad(full(b))[free],
              hess = if (is.null(fn$hess)) NULL else function(b) fn$hess(full(b))[free, free, drop = FALSE],
              bound = if (is.null(fn$bound)) NULL else function(b) fn$bound(full(b))[free])
  if (!is.null(sr)) {
    ur <- function(v) { u <- u0; u[free] <- v; u }
    fnr$search <- list(to = function(b) sr$to(full(b))[free], from = function(v) sr$from(ur(v))[free],
                       grad = function(v, gb) { g <- numeric(length(t0)); g[free] <- gb; sr$grad(ur(v), g)[free] })
  }
  r <- tryCatch(.pc_optimize(fnr, t0[free], scale[free], maxit, reltol), error = function(e) { .pc_pass_limit(e); NULL })
  if (is.null(r) || !is.finite(r$value) || r$value > op$value + 1e-10 * max(1, abs(op$value))) return(NULL)
  if (!isTRUE(r$converged)) {
    ## near the boundary the remaining directions can be very unequally curved
    ## (a degenerate nest: the within-nest coefficient differences shrink with
    ## lambda, their curvature growing like 1 / lambda^2, so the Hessian's
    ## diagonal spans ten orders of magnitude), and a scaled gradient above the
    ## usual tolerance then moves the log-likelihood by less than rounding: the
    ## fit is converged when a Newton step would lower the objective by less
    ## than 1e-9 of it (the Newton decrement)
    sf <- scale[free]; gs <- function(ts) fnr$grad(ts / sf) / sf
    ts <- r$par * sf; g <- gs(ts)
    decrement <- function(H) {
      if (is.null(H) || !all(is.finite(H))) return(NA_real_)
      H <- (H + t(H)) / 2; attr(H, "asym") <- NULL
      st <- tryCatch(solve(H, g), error = function(e) { .pc_pass_limit(e); NULL })
      if (is.null(st) || !all(is.finite(st))) NA_real_ else sum(g * st) / 2
    }
    thr <- 1e-9 * max(1, abs(r$value))
    dec <- decrement(tryCatch(if (is.null(fnr$hess)) .pc_hess_fd(gs, ts) else fnr$hess(r$par) / tcrossprod(sf), error = function(e) { .pc_pass_limit(e); NULL }))
    ## central differences of the gradient are too noisy at that curvature
    ## (their asymmetry reaches 1e6, with spurious negative eigenvalues), so an
    ## inconclusive decrement is recomputed from a Richardson-extrapolated Jacobian
    if (!(is.finite(dec) && dec >= 0 && dec < thr))
      dec <- decrement(tryCatch(numDeriv::jacobian(gs, ts), error = function(e) { .pc_pass_limit(e); NULL }))
    if (is.finite(dec) && dec >= 0 && dec < thr) { r$converged <- TRUE; r$max_grad <- max(abs(g)) }
  }
  list(par = if (is.null(fn$bound)) full(r$par) else fn$bound(full(r$par)), value = r$value, converged = r$converged, max_grad = r$max_grad, bfgs_code = r$bfgs_code,
       newton_it = r$newton_it, boundary_fixed = fixed)
}

## the uninflated design of an inflated one (for starting values and the
## baseline comparison)
.pc_strip_inflation <- function(D) {
  D0 <- D
  D0$k <- 0L; D0$kc <- integer(0); D0$Z <- matrix(0, D$n, 0L); D0$r <- 0L
  D0$par_names <- D$par_names[-D$idx$delta]; D0$blocks <- D$blocks[-D$idx$delta]
  D0$idx$delta <- integer(0); D0$K <- length(D0$par_names)
  D0
}

## the multinomial-logit design underlying a nested logit / HEV / probit design
## (same beta / gamma blocks, no structural block): the starting fit
.pc_as_mnl <- function(D) {
  D0 <- D
  if (length(D$idx$str)) { D0$par_names <- D$par_names[-D$idx$str]; D0$blocks <- D$blocks[-D$idx$str] }
  D0$idx$str <- integer(0); D0$family <- "mnl"; D0$str <- .pc_structural("mnl", NULL, D$levels, D$base)
  D0$K <- length(D0$par_names)
  D0
}

## Starting values. MNL: zero (the log-likelihood is globally concave). imnl:
## the inflated likelihood can have several local maxima (two on each bundled
## data set, and many separation-induced ones on others), so a grid of fourteen
## starts is run and the best kept:
## beta from the MNL fit or from an MNL that down-weights the inflated-category
## cases by a first-pass P(regime 1); delta from a binary fit of 1{y != k} on z
## (full, halved, doubled, tripled, or intercept-only slopes) with shifts of
## the split intercept.
## The log-likelihoods of the screened starts (object$screen_logliks) and of the
## polished fits (object$start_logliks) are stored; the polished fits' ends are
## grouped into distinct points and classified (object$search,
## .pc_search_record), and the number of distinct local maxima among them is
## reported (object$n_modes).
.pc_starts <- function(D, fn, scale, maxit, reltol) {
  if (D$family == "rologit") return(list(numeric(D$K)))           # globally concave, as the MNL
  if (D$family != "mnl") {
    ## beta / gamma from the multinomial logit on the same design, structural
    ## parameters at the MNL (lambda = 1, theta = 1) and one perturbed value
    D0 <- .pc_as_mnl(D)
    b0 <- .pc_optimize(.pc_fns(D0), numeric(D0$K), scale[-D$idx$str], maxit, reltol, newton_first = TRUE)$par
    ## probit utilities are on a different scale: the logit's error differences
    ## have variance pi^2/3, the probit's first difference variance 1
    if (D$family == "mnp") b0 <- b0 / sqrt(pi^2 / 3)
    ns <- D$str$n
    s_alt <- switch(D$family, nlogit = rep(log(0.7), ns), hev = rep(-0.3, ns), mnp = D$str$start_alt, remnl = rep(log(0.5), ns), numeric(ns))
    return(list(c(b0, D$str$start %||% numeric(ns)), c(b0, s_alt)))
  }
  if (D$k == 0L) return(list(numeric(D$K)))
  D0 <- .pc_strip_inflation(D)
  fn0 <- .pc_fns(D0)
  b0 <- .pc_optimize(fn0, numeric(D0$K), scale[-D$idx$delta], maxit, reltol, newton_first = TRUE)$par
  if (length(D$kc) > 1L) {
    ## several captive alternatives: for each, the single-alternative start (a
    ## logit of 1{y != c} on z, which counts every chooser of c as captive and so
    ## overstates captivity), at full and halved slopes, the intercepts shifted
    ## together towards less captivity
    ic <- match("(Intercept)", colnames(D$Z))
    st <- lapply(D$kc, function(c) {
      notc <- as.numeric(D$y != c)
      d_int <- numeric(D$r); if (!is.na(ic)) d_int[ic] <- stats::qlogis(min(max(stats::weighted.mean(notc, D$w), 0.05), 0.95))
      d0 <- tryCatch(suppressWarnings(stats::glm.fit(D$Z, notc, weights = D$w, family = stats::binomial()))$coefficients, error = function(e) { .pc_pass_limit(e); NULL })
      if (is.null(d0) || anyNA(d0) || any(!is.finite(d0))) d0 <- d_int
      list(d0 = d0, d_int = d_int)
    })
    shift <- function(dd, s) { if (!is.na(ic)) dd[ic] <- dd[ic] + s else dd <- dd + s / max(1, D$r); dd }
    mk <- function(f, s, which = "d0") c(b0, unlist(lapply(st, function(x) shift(f * x[[which]], s))))
    return(list(mk(1, 0), mk(1, 1), mk(1, 2), mk(1, 3), mk(0.5, 1), mk(0.5, 2), mk(1, 1, "d_int"), mk(1, 2.5, "d_int")))
  }
  notk <- as.numeric(D$y != D$k)
  fam <- if (D$link == "logit") stats::binomial("logit") else stats::binomial("probit")
  share <- min(max(stats::weighted.mean(notk, D$w), 0.05), 0.95)
  ic <- match("(Intercept)", colnames(D$Z))
  d_int <- numeric(D$r); if (!is.na(ic)) d_int[ic] <- fam$linkfun(share)
  d0 <- tryCatch(suppressWarnings(stats::glm.fit(D$Z, notk, weights = D$w, family = fam))$coefficients, error = function(e) { .pc_pass_limit(e); NULL })
  if (is.null(d0) || anyNA(d0) || any(!is.finite(d0))) d0 <- d_int
  shift <- function(dd, s) { if (!is.na(ic)) dd[ic] <- dd[ic] + s else dd <- dd + s / max(1, D$r); dd }
  ## an EM-flavoured beta start: the inflated-category cases weighted by the first-pass P(regime 1)
  phi0 <- fam$linkinv(as.vector(D$Z %*% d0))
  Dw <- D0; Dw$w <- D$w * ifelse(D$y == D$k, pmin(pmax(phi0, 0.02), 1), 1)
  bw <- tryCatch(.pc_optimize(.pc_fns(Dw), b0, scale[-D$idx$delta], maxit, reltol, newton_first = TRUE)$par, error = function(e) { .pc_pass_limit(e); b0 })
  ## sharpened split slopes (doubled and tripled, the intercept shifted towards
  ## the choice regime) reach the maxima at which the split separates a few
  ## cases into the inflated regime: on Glasgow (2001) every one of the ten
  ## starts above polished to -1462.76 while maxima at -1442 to -1453 exist
  list(c(b0, shift(d0, 0)), c(b0, shift(d0, 1)), c(b0, shift(d0, 2)), c(b0, shift(d0, -1)),
       c(b0, shift(0.5 * d0, 0)), c(b0, shift(0.5 * d0, 1)),
       c(b0, shift(d_int, 0)), c(b0, shift(d_int, 1.5)),
       c(bw, shift(d0, 0)), c(bw, shift(d0, 1)),
       c(b0, shift(2 * d0, 1)), c(b0, shift(2 * d0, 2)), c(b0, shift(3 * d0, 1)), c(b0, shift(3 * d0, 2)))
}

## covariance of the estimates: analytic (inverse observed information), robust
## (sandwich over cases), or cluster (sandwich over clusters), computed in the
## standardized parameterization ts = theta * scale and mapped back
## (V_theta = V_s / (scale scale')). With frequency weights the meat is
## sum_i w_i s_i s_i' (a weight of w equals w copies of the case), the bread the
## Hessian of the weighted log-likelihood. The reciprocal condition number is
## that of the diagonally normalized Hessian (a scale-free measure).
.pc_vcov <- function(theta, fn, type, cluster, scale, w) {
  ## a penalized fit (firth = TRUE) inverts the model's information at its
  ## estimates (fn$info), not the Hessian of the penalized objective; a fit
  ## with a prior inverts the curvature of the log posterior (fn$vcov_info:
  ## the information plus the prior's), or, where that is not positive
  ## definite, the scoring's curvature (the prior's negative part left out;
  ## attribute prior_clipped, DESIGN.md Phase 4c)
  clipped <- FALSE
  Hs <- tryCatch(if (!is.null(fn$vcov_info)) {
                   Hv <- fn$vcov_info(theta)
                   dv <- sqrt(pmax(abs(diag(Hv)), .Machine$double.xmin))
                   if (!all(is.finite(Hv)) || inherits(tryCatch(chol(Hv / tcrossprod(dv)), error = function(e) { .pc_pass_limit(e); e }), "error")) { clipped <- TRUE; Hv <- fn$info(theta) }
                   Hv / tcrossprod(scale)
                 } else if (!is.null(fn$info)) fn$info(theta) / tcrossprod(scale)
                 else if (is.null(fn$hess)) .pc_hess_fd(function(ts) fn$grad(ts / scale) / scale, theta * scale)
                 else fn$hess(theta) / tcrossprod(scale), error = function(e) { .pc_pass_limit(e); NULL })
  if (is.null(Hs) || !all(is.finite(Hs))) return(NULL)
  asym <- attr(Hs, "asym"); attr(Hs, "asym") <- NULL
  dn <- sqrt(pmax(abs(diag(Hs)), .Machine$double.eps))
  Hn <- Hs / tcrossprod(dn)
  noise <- if (is.null(asym)) NULL else asym / tcrossprod(dn)                 # normalized asymmetry, per entry
  rc <- tryCatch(rcond(Hn), error = function(e) { .pc_pass_limit(e); NA_real_ })
  bread <- tryCatch(solve(Hs), error = function(e) { .pc_pass_limit(e); NULL })
  used_ginv <- is.null(bread)
  if (used_ginv) bread <- MASS::ginv(Hs)
  Vs <- if (type == "analytic") bread else {
    S <- fn$scores(theta)
    if (anyNA(S)) return(NULL)
    S <- sweep(S, 2L, scale, "/")
    meat <- if (type == "cluster") {
      if (is.null(cluster)) stop("se = \"cluster\" requires 'cluster'.", call. = FALSE)
      crossprod(as.matrix(rowsum(S * w, group = cluster)))
    } else crossprod(S * sqrt(w))
    if (type == "opg") tryCatch(solve(meat), error = function(e) { .pc_pass_limit(e); MASS::ginv(meat) }) else bread %*% meat %*% bread
  }
  V <- Vs / tcrossprod(scale)
  attr(V, "rcond") <- rc; attr(V, "ginv") <- used_ginv; attr(V, "hn") <- Hn; attr(V, "noise") <- noise; attr(V, "hdiag") <- diag(Hs); attr(V, "hrow") <- apply(abs(Hs), 1L, max)
  if (clipped) attr(V, "prior_clipped") <- TRUE
  bd <- diag(bread); attr(V, "info_se") <- ifelse(bd > 0, sqrt(pmax(bd, 0)), NA_real_) / scale     # inverse-information se, whatever the type
  V
}

## the inverse-information standard errors at theta (NULL when the information
## is singular and a generalized inverse would stand in): the input of the
## separation rule's Wald criterion (.pc_diverge)
.pc_info_se <- function(theta, D, scale) {
  V <- tryCatch(.pc_vcov(theta, .pc_fns(D), "analytic", NULL, scale, .pc_obs(D)$w), error = function(e) { .pc_pass_limit(e); NULL })
  if (is.null(V) || isTRUE(attr(V, "ginv"))) NULL else attr(V, "info_se")
}

## the saddle check: the smallest eigenvalue of the normalized Hessian of the
## negative log-likelihood over the parameters not at a boundary (`keep`; a
## parameter at the boundary, a scale or a coefficient run off to infinity,
## has no curvature to speak of, and its diagonal can carry numerical noise of
## either sign). An eigenvalue below -tol_eig makes the point a saddle only if
## the log-likelihood rises along its eigenvector: on a curved ridge of exact
## non-identification (a saturated design, as the single binary covariate of
## DiGiuseppe et al. 2024) a straight line leaves the ridge, the Hessian has a
## slightly negative eigenvalue there, and the log-likelihood falls along it
## on both sides. The probe steps are 1e-3 to 1 normalized units (a unit is
## one conditional standard error); a gain below the optimizer's relative
## tolerance on the log-likelihood is no ascent. When the eigenvalue predicts
## a gain the smallest step would show and the log-likelihood curves down
## there instead, the Hessian contradicts the likelihood (as the corrupted
## differenced Hessian of the HEV quadrature error, ledger 2026-10-14).
## Returns the eigenvalue and "none", "saddle" or "inconsistent", and at a
## saddle the probe point where the objective fell most (`point`).
.pc_saddle <- function(nll, theta, scale, Hn, hdiag, keep, value, reltol, tol_eig) {
  out <- list(type = "none", min_eig = NA_real_)
  if (is.null(Hn) || !any(keep)) return(out)
  e <- tryCatch(eigen(Hn[keep, keep, drop = FALSE], symmetric = TRUE), error = function(e) { .pc_pass_limit(e); NULL })
  if (is.null(e) || !all(is.finite(e$values))) return(out)
  out$min_eig <- min(e$values)
  if (out$min_eig >= -tol_eig) return(out)
  v <- numeric(length(theta)); v[keep] <- e$vectors[, which.min(e$values)]
  v <- v / sqrt(pmax(abs(hdiag), .Machine$double.eps)) / scale          # one normalized unit, on the scale of theta
  tol_ll <- reltol * max(1, abs(value))
  steps <- c(1e-3, 1e-2, 1e-1, 1)
  f <- vapply(c(steps, -steps), function(s) tryCatch(nll(theta + s * v), error = function(e) { .pc_pass_limit(e); NA_real_ }), numeric(1))
  if (any(is.finite(f) & value - f > tol_ll)) {
    ## the probe point where the objective fell most: the saddle step's start
    out$type <- "saddle"; b <- which.min(ifelse(is.finite(f), f, Inf)); out$point <- theta + c(steps, -steps)[b] * v
    return(out)
  }
  c2 <- (f[1L] + f[5L] - 2 * value) / steps[1L]^2                         # the likelihood's own curvature at the smallest step
  if (is.finite(c2) && c2 > 0 && 0.5 * abs(out$min_eig) * steps[1L]^2 > tol_ll) out$type <- "inconsistent"
  out
}

## drop the diagnostic attributes .pc_vcov() attaches
.pc_strip_vcov <- function(V) {
  if (!is.null(V)) attr(V, "rcond") <- attr(V, "ginv") <- attr(V, "hn") <- attr(V, "noise") <- attr(V, "hdiag") <- attr(V, "hrow") <- attr(V, "info_se") <- attr(V, "prior_clipped") <- NULL
  V
}

## the observation level of a design: cases, or units for the random-intercept
## family (whose likelihood contributions and scores are per unit)
.pc_obs <- function(D) {
  if (D$family == "remnl") list(w = D$str$wu, cluster = D$str$ucluster) else list(w = D$w, cluster = D$cluster)
}

## Nonparametric bootstrap over cases (or over clusters / units). Frequency
## weights are honoured: without clusters a replicate draws multinomial weights
## with the weight total as its size, which is the bootstrap of the duplicated
## rows (exactly so for integer weights); for random intercepts without a user
## cluster the same is done over units (a unit of weight w is w independent
## units); with a user cluster the clusters are drawn and carry their rows'
## weights (the duplicated rows stay inside their cluster). With unit fixed effects the dummy
## columns are regenerated on every resample (never-chosen rule, dropped units),
## so a replicate may have fewer dummies than the fit; replicate vectors are
## aligned to the fit's parameters by name (NA where a dummy is absent). For the
## inflated model every replicate also runs the two deterministic challengers of
## the main fit so that a resample does not settle on the lower mode unseen.
## frequency weights under se = "bootstrap": a replicate draws as many cases
## as the weights total, so normalized or survey weights (a total of at most 2)
## would give a sample of two (round 7, G-N2), as for the priors'
## autoscaling. The check applies where the draw uses the weights: cases
## without a cluster, unit fixed effects or random intercepts, and remnl()
## without a cluster, which draws its units by their weights (local review
## 2026-11-06, finding 4). Run before the fit, so that the refusal does not
## wait for it, and again in .pc_boot().
.pc_boot_weight_check <- function(D) {
  if (!any(D$w != 1)) return(invisible())
  has_fe <- !is.null(D$fe_map) && nrow(D$fe_map) > 0L
  remnl_w <- D$family == "remnl" && is.null(D$cluster)
  plain <- is.null(D$cluster) && !has_fe && D$family != "remnl"
  if ((plain || remnl_w) && sum(if (remnl_w) D$str$wu else D$w) <= 2)
    stop("se = \"bootstrap\" draws as many cases as the frequency weights total, and these total at most 2 (normalized or survey weights are not ",
         "frequency weights). Rescale the weights to counts, or use another 'se'.", call. = FALSE)
  invisible()
}

.pc_boot <- function(D, theta, scale, nboot, cores, maxit, reltol) {
  n <- D$n
  cl <- D$cluster
  if (D$family == "remnl" && is.null(cl)) cl <- as.character(D$unit)        # units are the resampling blocks
  weighted <- any(D$w != 1)
  has_fe <- !is.null(D$fe_map) && nrow(D$fe_map) > 0L
  ## with unit fixed effects the units are resampled (a case bootstrap would
  ## give each unit a random number of cases and add a second incidental-
  ## parameters bias to the replicates); a unit drawn twice is two units
  if (has_fe && is.null(cl)) cl <- as.character(D$unit)
  ## a unit drawn twice becomes two units only when the clusters nest the units
  fe_nested <- has_fe && all(lengths(lapply(split(cl, D$unit), unique)) <= 1L)
  if (has_fe && !fe_nested)
    stop("With unit fixed effects the bootstrap resamples units; a 'cluster' that does not nest the units would resample cases within units ",
         "(a random number of cases per unit and a second incidental-parameters bias). Use the unit, or a clustering of whole units.", call. = FALSE)
  remnl_w <- D$family == "remnl" && is.null(D$cluster) && weighted
  .pc_boot_weight_check(D)                      # also before the fit (.pc_estimate)
  rows_by_cl <- if (!is.null(cl)) split(seq_len(n), cl) else NULL
  draw <- function() {
    if (remnl_w) {
      ui <- as.integer(D$unit); wu <- D$str$wu
      cnt <- as.vector(stats::rmultinom(1L, size = max(2L, round(sum(wu))), prob = wu))
      return(list(idx = seq_len(n), block = NULL, w = cnt[ui]))
    }
    if (is.null(cl)) {
      if (weighted) return(list(idx = seq_len(n), block = NULL, w = as.vector(stats::rmultinom(1L, size = max(2L, round(sum(D$w))), prob = D$w))))
      return(list(idx = sample.int(n, n, replace = TRUE), block = NULL, w = NULL))
    }
    u <- unique(cl); g <- sample(u, length(u), replace = TRUE)
    parts <- rows_by_cl[g]
    list(idx = unlist(parts, use.names = FALSE), block = rep(seq_along(g), lengths(parts)), w = NULL)
  }
  one <- function(b) {
    dr <- draw(); idx <- dr$idx
    ## the replicate design; a resample whose fixed effects cannot be built (an
    ## alternative left unconnected to the base, no usable case) is a failed
    ## replicate, not an error of the fit
    Db <- tryCatch({
      if (!is.null(dr$w)) {
        ## the replicate weights travel with their rows through the subsetting;
        ## the rows drawn are checked as an unweighted resample is (a weighted
        ## fit and its rows repeated by their weights bootstrap alike)
        Dw <- D; Dw$w <- dr$w; idx <- which(dr$w > 0)
        if (has_fe) .pc_half_design(Dw, idx) else .pc_subsample_check(.pc_subset_design(Dw, idx))
      } else if (has_fe) {
        ## a unit drawn twice is two units, each with its own fixed effects
        S <- .pc_subset_design(D, idx)
        if (fe_nested) S$unit <- factor(paste(dr$block, as.character(D$unit[idx]), sep = "/"))
        .pc_half_design(S, seq_len(S$n))
      } else .pc_subsample_check(.pc_subset_design(D, idx))
    }, error = function(e) { .pc_pass_limit(e); conditionMessage(e) })
    if (is.character(Db)) return(list(failed = Db))                 # the reason, tabulated below
    if (D$family == "remnl" && !is.null(dr$block)) {
      ## a unit drawn twice is two units in the resample
      Db$unit <- factor(paste(dr$block, as.character(D$unit[idx]), sep = "/"))
      Db <- .pc_re_index(Db)
    }
    fnb <- .pc_fns(Db)
    same <- identical(Db$par_names, D$par_names)
    th0 <- if (same) theta else { v <- theta[match(Db$par_names, D$par_names)]; v[is.na(v)] <- 0; v }
    scb <- if (same) scale else .pc_scales(Db)
    op <- tryCatch(.pc_optimize(fnb, th0, scb, maxit, reltol, newton_first = .pc_concave(D)), error = function(e) { .pc_pass_limit(e); NULL })
    if (is.null(op) || !is.finite(op$value)) return(list(failed = "the refit failed"))
    if (D$k > 0L) for (f in c(0.5, 1.5)) {
      chal <- th0; chal[Db$idx$delta] <- chal[Db$idx$delta] * f
      ch <- tryCatch(.pc_optimize(fnb, chal, scb, maxit, reltol), error = function(e) { .pc_pass_limit(e); NULL })
      if (!is.null(ch) && is.finite(ch$value) && ch$value < op$value - 1e-8) op <- ch
    }
    op <- tryCatch(.pc_boundary_step(fnb, Db, op, scb, maxit, reltol), error = function(e) { .pc_pass_limit(e); op })
    ## by name; a replicate may carry a dummy the fit has not (its reference unit
    ## can differ), which has no place in the fit's vector
    par <- if (same) op$par else { v <- rep(NA_real_, D$K); mi <- match(Db$par_names, D$par_names); v[mi[!is.na(mi)]] <- op$par[!is.na(mi)]; v }
    list(par = par, converged = op$converged, loglik = -op$value)
  }
  res <- .pc_lapply(seq_len(nboot), one, cores = cores)
  failed <- vapply(res, function(r) if (is.null(r$failed)) NA_character_ else r$failed, character(1))
  okres <- res[is.na(failed)]
  if (length(okres) < max(2L, nboot / 2))
    stop("Fewer than half of the ", nboot, " bootstrap replicates are usable", if (any(!is.na(failed))) paste0(" (most often: ", names(sort(table(failed), decreasing = TRUE))[1L], ")") else "",
         "; the bootstrap is not informative on these data.", call. = FALSE)
  ## the causes of the failed replicates, with the numbers masked so that one
  ## cause is one row
  fail_reasons <- if (any(!is.na(failed))) sort(table(gsub("(?<![A-Za-z0-9_.:])[-+]?[0-9][0-9.e+-]*", "#", failed[!is.na(failed)], perl = TRUE)), decreasing = TRUE) else NULL
  if (length(okres) < 2L) stop("The bootstrap produced fewer than two usable replicates.", call. = FALSE)
  Bi <- do.call(rbind, lapply(okres, `[[`, "par"))
  ## with unit fixed effects the constants and the unit dummies are
  ## normalizations: the reference unit of a constant can change between
  ## resamples, and a dummy is absent when its unit or choice is, so only the
  ## remaining parameters (slopes, split and structural parameters) get
  ## bootstrap replicates; the others are NA
  if (has_fe) Bi[, D$blocks == "fixed effect" | grepl(":\\(Intercept\\)$", D$par_names)] <- NA_real_
  B <- matrix(apply(Bi, 1L, .pc_natural, D = D), ncol = D$K, byrow = TRUE)
  colnames(B) <- colnames(Bi) <- D$par_names
  list(replicates = B, replicates_internal = Bi, logliks = vapply(okres, `[[`, numeric(1), "loglik"),
       n_failed = nboot - length(okres), n_nonconverged = sum(!vapply(okres, `[[`, logical(1), "converged")),
       fail_reasons = fail_reasons, cluster = !is.null(cl), weighted = weighted, nboot = nboot)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

## .pc_fit(): the engine. Arguments are documented on the exported wrappers.
.pc_fit <- function(formula, data, inflated = FALSE, inflate = NULL, ref = NULL, link = "logit",
                    id = NULL, alt = NULL, weights = NULL, offset = NULL, cluster = NULL, fe = NULL,
                    se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"), nboot = 200L, cores = 1L,
                    start = NULL, maxit = 1000L, reltol = 1e-10, call = NULL, se_missing = FALSE,
                    family = "mnl", str = NULL, bias_correct = "none", crosscheck = FALSE, firth = FALSE, prior = NULL) {
  se <- match.arg(se)
  .pc_check_whole(maxit, "maxit", 0); .pc_check_whole(cores, "cores")
  if (!is.logical(firth) || length(firth) != 1L || is.na(firth)) stop("'firth' must be TRUE or FALSE.", call. = FALSE)
  if (firth && !is.null(prior)) stop("Give either firth = TRUE (Jeffreys' prior) or a 'prior', not both.", call. = FALSE)
  if (!is.null(prior) && !inherits(prior, "coef_prior")) stop("'prior' must be made by coef_prior().", call. = FALSE)
  if ((firth || !is.null(prior)) && !is.null(fe))
    stop(if (firth) "firth = TRUE" else "A prior", " is not available with unit fixed effects ('fe'): the penalty would act on every unit constant, and its interplay with the ",
         "never-chosen rule and the jackknife correction is not established. Use cluster = unit, or mundlak() with remnl().", call. = FALSE)
  if (bias_correct == "jackknife" && is.null(fe)) stop("bias_correct = \"jackknife\" needs unit fixed effects ('fe').", call. = FALSE)
  if (!is.null(cluster) && se_missing) se <- "cluster"
  if (!is.null(cluster) && !se %in% c("cluster", "bootstrap"))
    warning("'cluster' is ignored unless se = \"cluster\" (or se = \"bootstrap\", where it defines the resampling blocks).", call. = FALSE)
  if (se == "cluster" && is.null(cluster)) stop("se = \"cluster\" requires 'cluster'.", call. = FALSE)
  if (se == "bootstrap") .pc_check_whole(nboot, "nboot", 2)
  D <- .pc_design(formula, data, inflated = inflated, inflate = inflate, ref = ref, link = link, id = id, alt = alt,
                  weights = weights, offset = offset, cluster = cluster, fe = fe, family = family, str = str)
  if (D$K == 0L) stop("The model has no parameters.", call. = FALSE)
  D$crosscheck <- isTRUE(crosscheck)            # kept in the design, so the parametric-bootstrap refits search as the fit did
  if (firth) D$firth <- TRUE                    # likewise: every refit of the design is penalized
  ## Jeffreys' penalty is half the log-determinant of the information, which
  ## the weights scale: read as frequency weights (the rows repeated) they
  ## must total more than 2, as a prior's autoscaling requires; normalized
  ## weights would let the penalty swamp the likelihood (local review
  ## 2026-11-06, finding 14)
  if (firth && sum(D$w) <= 2)
    stop("The weights total ", format(sum(D$w)), ": firth = TRUE reads them as frequency weights (the rows repeated), which needs a total above 2. ",
         "Normalized or survey weights are not supported; rescale frequency weights to counts.", call. = FALSE)
  if (!is.null(prior)) D$prior <- .pc_prior_resolve(prior, D)
  if (!is.null(cluster) && se == "cluster" && length(unique(D$cluster)) < 2L) stop("'cluster' is constant on the estimation rows.", call. = FALSE)
  if (se == "bootstrap") .pc_boot_weight_check(D)
  fit <- .pc_fit_design(D, se = se, nboot = nboot, cores = cores, start = start, maxit = maxit, reltol = reltol, call = call)
  if (bias_correct == "jackknife") fit <- .pc_apply_jackknife(fit, D, maxit, reltol)
  fit
}

## The structural log-parameters held at the search's bound along which the
## log-likelihood still rises: its gradient points out of the bound by more
## than the convergence tolerance (fn$grad is the gradient of the negative
## log-likelihood). There the likelihood is not flat, and the reported
## log-likelihood is below the model's supremum (stress finding F22: on Palik
## et al. 2025 by 0.67 or more). Returns their positions in theta.
.pc_rising_at_bound <- function(fn, theta, scale, str_log, value) {
  at <- str_log[abs(theta[str_log]) >= .pc_str_bound - 1e-8]
  if (!length(at)) return(integer())
  g <- tryCatch(fn$grad(theta), error = function(e) { .pc_pass_limit(e); NULL })
  if (is.null(g)) return(integer())
  up <- -sign(theta[at]) * g[at] / scale[at]
  at[is.finite(up) & up > .pc_gtol(value)]
}

## The warning of a fit whose structural log-parameters `at` are past the
## boundary flag, of which those named in `rising` are at the search's bound
## with the log-likelihood still rising (.pc_rising_at_bound). The probit does
## not nest the logit, so its comparison with mnl() is vuong(), not lr_test().
## With four or more alternatives the probit's probabilities are GHK
## simulations, and near a singular covariance one draw's term can switch
## between about 0 and 1 over a small change in the parameters (on a simulated
## design at a log Cholesky diagonal of -8 the derivative along a coefficient
## rises in 110 of 200 steps of 1e-4, where the exact likelihood is concave;
## stress finding F24), so there the simulated likelihood is not smooth and its
## gradient does not measure the distance to the maximum.
.pc_str_boundary_msg <- function(D, at, rising = character()) {
  cmp <- if (D$family == "mnp") "vuong()" else "lr_test()"
  paste0("The fit is at a boundary: ", paste(at, collapse = ", "), if (length(at) > 1L) " are at " else " is at ",
         switch(D$family,
                hev = "zero or infinity (the alternative's utility is estimated as deterministic relative to the others)",
                nlogit = "zero or infinity (a degenerate nest)",
                remnl = "zero or infinity (at zero the model reduces to the multinomial logit in that alternative)",
                "zero or infinity (a singular covariance)"),
         if (length(rising))
           paste0("; the log-likelihood still rises past the search's bound (+/-8 on the log scale) in ", paste(rising, collapse = ", "),
                  ", so the reported log-likelihood is below the model's supremum and a comparison with mnl() by ", cmp,
                  " or information criteria understates this model; the standard errors are unreliable.")
         else paste0("; the likelihood is flat in that direction and the standard errors are unreliable. Compare with mnl() by ",
                     cmp, " / information criteria."),
         if (D$family == "mnp" && D$J > 3L)
           " Near a singular covariance the GHK simulated likelihood is not smooth (one draw's contribution can jump over a small change in the parameters), so the gradient and the convergence check are unreliable there too."
         else "")
}

## .pc_fit_design(): estimation, covariances, and the object from a finished
## design. Also the entry point of the parametric bootstrap (boot_pvalue()),
## which replaces D$y and refits without re-parsing the data.
.pc_fit_design <- function(D, se = "analytic", nboot = 200L, cores = 1L, start = NULL, maxit = 1000L, reltol = 1e-10, call = NULL,
                           check = TRUE, .saddle = 0L, .inherit = NULL) {
  fn <- .pc_fns(D)
  scale <- .pc_scales(D)
  if (!is.null(start)) {
    if (length(start) != D$K) stop("'start' must have length ", D$K, " (the number of parameters).", call. = FALSE)
    starts <- list(as.numeric(start))
  } else starts <- .pc_starts(D, fn, scale, maxit, reltol)
  screen_vals <- NULL; crosscheck <- NULL
  if (length(starts) == 1L) {
    op <- .pc_optimize(fn, starts[[1L]], scale, maxit, reltol, newton_first = .pc_concave(D)); vals <- op$value
    ends <- list(op); best <- 1L
  } else {
    ms <- .pc_multistart(fn, scale, starts, maxit, reltol)
    if (isTRUE(D$crosscheck) && D$family == "mnl" && D$k > 0L && length(D$kc) == 1L) ms <- .pc_crosscheck(D, fn, scale, ms, maxit, reltol)
    if (D$k > 0L) ms <- .pc_challenge(fn, D, scale, ms, maxit, reltol)
    op <- ms$op; vals <- ms$vals; conv <- ms$conv; screen_vals <- ms$screen_vals; crosscheck <- ms$crosscheck
    ends <- ms$ends; best <- ms$best
  }
  ## a penalized fit (Firth's penalty or a prior) without a user start: a
  ## second run from the maximum-likelihood estimate, the higher maximum kept
  ## (.pc_firth_second); `firth` below means penalized either way
  firth <- .pc_penalized(D); firth_ml <- NULL; firth_modes <- NULL
  if (firth && is.null(start)) {
    f2 <- .pc_firth_second(D, fn, scale, maxit, reltol, op)
    op <- f2$op; vals <- f2$vals; firth_ml <- f2$ml; firth_modes <- f2$modes
  }
  op <- .pc_boundary_step(fn, D, op, scale, maxit, reltol)
  theta <- op$par; names(theta) <- D$par_names
  coefficients <- .pc_natural(theta, D); names(coefficients) <- D$par_names
  obs <- .pc_obs(D)
  ## the information-based diagnostics' inputs, before anything is reported, so
  ## that a fit stopped at a saddle can be moved off it first (below)
  V <- NULL; boot <- NULL; Vn <- NULL
  if (se %in% c("analytic", "opg", "robust", "cluster")) V <- .pc_vcov(theta, fn, se, obs$cluster, scale, obs$w)
  ## the conditioning and saddle diagnostics come from the information matrix
  ## whatever the standard errors are (computed here for the bootstrap and for
  ## se = "none"; internal refits pass check = FALSE)
  Vd <- if (!is.null(V) && !is.null(attr(V, "hn"))) V else if (check) .pc_vcov(theta, fn, "analytic", NULL, scale, obs$w) else NULL
  rc <- if (is.null(Vd) || is.null(attr(Vd, "rcond"))) NA_real_ else attr(Vd, "rcond")
  hn <- if (is.null(Vd)) NULL else attr(Vd, "hn")
  noise_m <- if (is.null(Vd)) NULL else attr(Vd, "noise")
  hdiag <- if (is.null(Vd)) NULL else attr(Vd, "hdiag"); hrow <- if (is.null(Vd)) NULL else attr(Vd, "hrow")
  used_ginv <- !is.null(Vd) && isTRUE(attr(Vd, "ginv"))
  ## (1) separation / boundary: a coefficient whose standardized
  ## magnitude (a one-SD effect) exceeds 15, or a 0/1 column (a dummy, a unit
  ## constant) whose own coefficient exceeds 15: rare dummies have a tiny SD and
  ## the one-SD criterion would miss them, while an effect-over-the-range
  ## criterion misfires on wide continuous covariates (travel minutes); and
  ## whose Wald |z| is below 1.96 (the likelihood flat along it), which spares
  ## large but precisely estimated coefficients; for the inflated model a split
  ## equation with phi at 1 (or 0) everywhere is the MNL boundary
  ## In mnl() and rologit() without fixed effects or inflation the exact check
  ## replaces that rule (.pc_separation: a linear program on the data's
  ## comparisons, DESIGN.md Phase 4b): the diverging coefficients are those
  ## along the separating direction, and an overlapping data set has none,
  ## however large its estimates. A penalized fit's estimates stay finite; the
  ## check names the separated coefficients of its data (`separation`)
  se_info <- if (!is.null(Vd) && !used_ginv) attr(Vd, "info_se") else NULL
  sepchk <- if (check) .pc_separation_fit(D, if (firth && !is.null(firth_ml) && all(is.finite(firth_ml$par))) firth_ml$par else theta) else NULL
  diverge <- if (firth) rep(FALSE, length(theta)) else if (!is.null(sepchk)) D$par_names %in% sepchk$terms else .pc_diverge(theta, D, scale, se = se_info)
  str_log <- .pc_str_log(D)
  ## (3) a negative eigenvalue of the normalized Hessian over the interior
  ## parameters along which the log-likelihood rises: a saddle, not a maximum
  ## (.pc_saddle). A differenced Hessian is only as accurate as its noise
  ## (measured by its asymmetry), so an eigenvalue within ten times that noise
  ## of zero is not evidence of a saddle
  ## a parameter with no curvature and no coupling (its whole row of the raw
  ## Hessian below 1e-8 of the largest diagonal entry: an inflation coefficient
  ## on a covariate whose cases all choose the inflated alternative with
  ## certainty) is a flat, unidentified direction, reported by the
  ## singular-information warning; normalized by the root of its near-zero
  ## diagonal, its round-off would read as an eigenvalue of -1, a spurious saddle
  ## a penalized fit: the saddle check is on the Hessian of the penalized
  ## objective (central differences of its gradient), not on the information
  if (firth && check) {
    Hp <- tryCatch(.pc_hess_fd(function(ts) fn$grad(ts / scale) / scale, theta * scale), error = function(e) { .pc_pass_limit(e); NULL })
    if (!is.null(Hp) && all(is.finite(Hp))) {
      dn <- sqrt(pmax(abs(diag(Hp)), .Machine$double.eps))
      noise_m <- attr(Hp, "asym") / tcrossprod(dn); attr(Hp, "asym") <- NULL
      hn <- Hp / tcrossprod(dn); hdiag <- diag(Hp); hrow <- apply(abs(Hp), 1L, max)
    }
  }
  flat <- if (is.null(hrow) || length(hrow) != length(theta)) rep(FALSE, length(theta)) else hrow <= 1e-8 * max(abs(hdiag))
  at_bound <- diverge | flat | (seq_along(theta) %in% str_log[abs(theta[str_log]) > 7])
  noise <- if (is.null(noise_m) || !any(!at_bound)) 0 else max(noise_m[!at_bound, !at_bound, drop = FALSE])   # over the parameters checked
  sad <- .pc_saddle(fn$nll, theta, scale, hn, hdiag, !at_bound, op$value, reltol, max(1e-8, 10 * noise))
  ## the saddle step (DESIGN.md Phase 4d): at a saddle the search restarts from
  ## the probe point of .pc_saddle where the objective fell most along the
  ## eigenvector of the negative eigenvalue; when the refit there ends higher,
  ## the fit is redone from that maximum (at most three times), the
  ## multistart's record kept. A refit that ends no higher leaves the saddle
  ## reported as before.
  if (check && sad$type == "saddle" && !is.null(sad$point) && .saddle < 3L) {
    esc <- tryCatch(.pc_optimize(fn, sad$point, scale, maxit, reltol), error = function(e) { .pc_pass_limit(e); NULL })
    if (!is.null(esc) && is.finite(esc$value)) esc <- tryCatch(.pc_boundary_step(fn, D, esc, scale, maxit, reltol), error = function(e) { .pc_pass_limit(e); esc })
    if (!is.null(esc) && is.finite(esc$value) && esc$value < op$value - 1e-8 * max(1, abs(op$value))) {
      message("The fit stopped at a saddle (", if (firth) "penalized " else "", "log-likelihood ", format(-op$value, digits = 10),
              "); a step along the direction in which it rises and a refit reach ", format(-esc$value, digits = 10), ".")
      ## the search's ends are kept (the saddle among them, which the record
      ## of the refit classifies as one)
      inh <- if (is.null(.inherit)) list(vals = vals, screen_vals = screen_vals, crosscheck = crosscheck, n_starts = length(starts),
                                         ends = ends, firth_modes = firth_modes) else .inherit
      return(.pc_fit_design(D, se = se, nboot = nboot, cores = cores, start = esc$par, maxit = maxit, reltol = reltol, call = call, check = check,
                            .saddle = .saddle + 1L, .inherit = inh))
    }
  }
  if (!is.null(.inherit)) {
    ## a refit from a saddle step: the record of the original search, with
    ## the maximum the step reached (the reported one)
    vals <- c(.inherit$vals, op$value); screen_vals <- .inherit$screen_vals; crosscheck <- .inherit$crosscheck
    firth_modes <- .inherit$firth_modes
    best <- length(.inherit$ends) + best; ends <- c(.inherit$ends, ends)
  }
  if (D$family == "nlogit" && any(coefficients[D$idx$str] > 1 + 1e-8))
    message("Inclusive-value parameter(s) above 1 (", paste(D$par_names[D$idx$str][coefficients[D$idx$str] > 1 + 1e-8], collapse = ", "),
            "): outside the range consistent with random-utility maximization for all covariate values; the nesting may be wrong or unnecessary (lr_test() against the mnl() fit).")
  n_starts <- if (!is.null(.inherit)) .inherit$n_starts else length(starts)
  if (!is.null(firth_modes))
    warning("The penalized likelihood has two local maxima (penalized log-likelihoods ", formatC(firth_modes[1L], digits = 6, format = "f"),
            " from the start at zero and ", formatC(firth_modes[2L], digits = 6, format = "f"),
            " from the maximum-likelihood estimate); the higher is reported. See object$start_logliks.", call. = FALSE)
  if (se %in% c("analytic", "opg", "robust", "cluster")) {
    if (is.null(V)) warning("The covariance matrix could not be computed (non-finite Hessian); standard errors are NA.", call. = FALSE)
  } else if (se == "bootstrap") {
    boot <- .pc_boot(D, theta, scale, nboot, cores, maxit, reltol)
    V <- stats::cov(boot$replicates_internal)
    Vn <- stats::cov(boot$replicates)
    if (!is.null(D$fe_map) && nrow(D$fe_map))
      message("With unit fixed effects the bootstrap resamples units and covers the slopes: the alternative-specific constants and the unit effects have no bootstrap standard errors (NA).")
    if (boot$n_failed) warning(boot$n_failed, " of ", nboot, " bootstrap replicates failed and were dropped (most often: ", names(boot$fail_reasons)[1L],
                               "; object$boot$fail_reasons tabulates the causes)",
                               if (!is.null(D$fe_map) && nrow(D$fe_map)) "; the standard errors are conditional on resamples in which the fixed effects can be built" else "", ".", call. = FALSE)
    if (boot$n_nonconverged) warning(boot$n_nonconverged, " of ", nboot - boot$n_failed, " bootstrap replicates did not converge; their estimates enter the bootstrap covariance (see object$boot$logliks).", call. = FALSE)
  }
  if (!is.null(V)) { V <- .pc_strip_vcov(V); dimnames(V) <- list(D$par_names, D$par_names) }
  ## the reported covariance is on the natural scale (delta method through the
  ## structural transform; identity for the multinomial logit)
  if (!is.null(V) && is.null(Vn)) {
    Vn <- if (length(D$idx$str)) { Jm <- numDeriv::jacobian(function(t) .pc_natural(t, D), theta); Jm %*% V %*% t(Jm) } else V
  }
  if (!is.null(Vn)) dimnames(Vn) <- list(D$par_names, D$par_names)
  if (D$family == "remnl") {
    pru <- .pc_re_pred_rule(max(exp(theta[D$idx$str])), D$str$nodes, D$J - 1L)
    if (pru$capped && pru$err > 1e-6)
      warning("The population-averaged probabilities (fitted(), predict()) use a quadrature rule capped at ", length(pru$logw), " points; at these ",
              "standard deviations their error may reach about ", format(pru$err, digits = 2), " (the likelihood itself is not affected).", call. = FALSE)
  }
  ## refits used only for their maximum (check = FALSE: boot_lr_test()) skip the
  ## random-intercept model's fine prediction rule
  pr <- if (!check && D$family == "remnl") fn$probs(theta, coarse = TRUE) else fn$probs(theta)
  P <- pr$prob; colnames(P) <- D$levels; rownames(P) <- D$case_id
  Pout <- if (is.null(pr$prob_outcome)) P else pr$prob_outcome; colnames(Pout) <- D$levels; rownames(Pout) <- D$case_id
  if (D$family == "hev") pc_hev_unsettled(TRUE)
  lli <- fn$lli(theta)
  if (D$family == "hev" && (nu <- pc_hev_unsettled(TRUE)) > 0L)
    warning("At the estimate the HEV integral of ", nu, " case(s) stopped at its panel or depth limit before reaching its 1e-10 tolerance; ",
            "their log-likelihood contributions are less accurate (a scale ratio far outside the data's range).", call. = FALSE)
  ## diagnostics: (1) separation and boundaries (the diverging coefficients,
  ## computed above)
  separation <- if (!firth || !check) character(0) else if (!is.null(sepchk)) sepchk$terms
                else .pc_firth_separation(D, scale, maxit, reltol, theta = theta, se = se_info, op = firth_ml)
  if (length(separation))
    message("The maximum-likelihood estimates of ", paste(utils::head(separation, 6L), collapse = ", "), if (length(separation) > 6L) ", ..." else "",
            " diverge (separation); their penalized estimates are determined largely by the penalty. Profile penalized-likelihood intervals (confint()) ",
            "and penalized likelihood-ratio tests are the recommended inference.")
  pb <- .pc_phi_boundary(D, pr); phi_boundary <- pb$phi_boundary; collapsed <- pb$collapsed
  boundary <- any(diverge) || phi_boundary
  ## a log scale or log inclusive-value parameter beyond +/- 7 (a ratio above
  ## 1100) is the family's boundary: a nearly deterministic alternative (HEV) or
  ## a degenerate nest; the likelihood is usually flat there and the information
  ## singular. Where it still rises at the search's bound, the warning says so
  str_boundary <- length(str_log) > 0L && any(abs(theta[str_log]) > 7)
  boundary <- boundary || str_boundary
  rising <- if (str_boundary) .pc_rising_at_bound(fn, theta, scale, str_log, op$value) else integer()
  ## under separation the log-likelihood rises without bound along the
  ## diverging coefficients, so their gradient never vanishes and the
  ## optimizer's tolerance cannot be met; the fit is converged when the
  ## remaining parameters are (the divergence is reported on its own)
  if (!isTRUE(op$converged) && any(diverge)) {
    rest <- !diverge & !(seq_along(theta) %in% str_log[abs(theta[str_log]) > 7])
    gr <- tryCatch(fn$grad(theta) / scale, error = function(e) { .pc_pass_limit(e); NULL })
    if (any(rest) && !is.null(gr) && all(is.finite(gr[rest])) && max(abs(gr[rest])) < .pc_gtol(op$value)) {
      op$converged <- TRUE; op$max_grad <- max(abs(gr[rest]))
    }
  }
  ## (3) the saddle check (computed above, with the saddle step)
  ill <- !is.null(Vd) && (used_ginv || (is.finite(rc) && rc < 1e-9))
  if (sad$type == "saddle") {
    op$converged <- FALSE
    warning("The Hessian at the solution is not positive definite (smallest normalized eigenvalue ", format(sad$min_eig, digits = 2),
            "): the point is a saddle, not a maximum, and the standard errors are not valid; ", if (.saddle >= 3L) "three steps off saddles have been taken (the limit)" else "a step off it along the direction of ascent did not end higher",
            ". Supply other starting values ('start').", call. = FALSE)
  } else if (!op$converged) warning("The optimizer did not converge (largest scaled gradient ", format(op$max_grad, digits = 3), "); inspect the fit.", call. = FALSE)
  if (sad$type == "inconsistent")
    warning("The Hessian at the solution has a negative eigenvalue (smallest normalized eigenvalue ", format(sad$min_eig, digits = 2),
            ") that the log-likelihood does not bear out (it falls on both sides along that direction): the Hessian is inaccurate there and the standard errors are not valid.", call. = FALSE)
  if (boundary) warning(if (str_boundary) .pc_str_boundary_msg(D, D$par_names[str_log][abs(theta[str_log]) > 7], D$par_names[rising])
                        else if (phi_boundary && length(collapsed))
                          paste0("The fit is at a boundary: captivity to ", paste(collapsed, collapse = ", "), " has collapsed (its share is below 1e-6 for every case); ",
                                 "the fit is the model without ", if (length(collapsed) > 1L) "them" else "it", "; compare the two by boot_lr_test().")
                        else if (phi_boundary) "The fit is at a boundary: the split equation puts every case in one regime (the inflated model collapses onto the MNL); compare with mnl() by boot_lr_test(), or by lr_test() when the split has only an intercept."
                        else if (isTRUE(sepchk$separated))
                          paste0("The data are ", if (sepchk$type == "complete") "completely" else "quasi-completely", " separated (see separation()): the maximum-likelihood estimates of ",
                                 paste(utils::head(D$par_names[diverge], 6L), collapse = ", "), if (sum(diverge) > 6L) ", ..." else "",
                                 " do not exist (the log-likelihood rises along a direction in them without reaching a maximum); the reported values are where the optimizer stopped, ",
                                 "and they and their standard errors are not meaningful. firth = TRUE gives finite estimates.")
                        else paste0("Coefficient(s) diverge (", paste(utils::head(D$par_names[diverge], 6L), collapse = ", "), if (sum(diverge) > 6L) ", ..." else "", "): ",
                                    if (any(D$blocks[diverge] == "fixed effect")) "a unit constant at infinity (a unit with a nearly constant choice) or " else "",
                                    "possible quasi-complete separation; the estimates and standard errors of those parameters are unreliable."), call. = FALSE)
  if (ill) warning(if (used_ginv) "The information matrix is singular (a generalized inverse was used); some parameters are not identified by the data."
                   else paste0("The information matrix is ill-conditioned (reciprocal condition number ", format(rc, digits = 2), "); the model may be weakly identified."), call. = FALSE)
  ## the record of the search's ends and the number of distinct local maxima
  ## among them (.pc_search_record). A penalized fit counts the maxima of its
  ## two runs (.pc_firth_second); refits used only for their maximum
  ## (check = FALSE) keep no record
  search <- NULL
  if (firth) n_modes <- if (is.null(start) || !is.null(.inherit)) 1L + !is.null(firth_modes) else 1L
  else if (!check) n_modes <- NA_integer_
  else {
    rep_type <- if (boundary) "boundary" else if (sad$type == "saddle") "saddle" else if (!isTRUE(op$converged)) "not converged" else "maximum"
    rep_detail <- switch(rep_type, boundary = .pc_boundary_detail(D, diverge, str_log[abs(theta[str_log]) > 7], pb),
                         saddle = sprintf("smallest normalized eigenvalue %.2g", sad$min_eig),
                         `not converged` = sprintf("largest scaled gradient %.2g", op$max_grad),
                         maximum = if (ill) "ill-conditioned information" else "")
    search <- .pc_search_record(D, fn, scale, obs$w, ends, best, op, list(type = rep_type, detail = rep_detail, se = se_info), reltol)
    n_modes <- sum(search$table$type == "maximum")
    if (D$k > 0L) .pc_search_message(search$table, n_modes, n_starts)
  }
  obj <- structure(list(
    call = call, formula = D$F, layout = D$layout, levels = D$levels, ref = D$levels[D$base], base = D$base,
    family = D$family, inflated = D$k > 0L, inflate = if (D$k > 0L) D$levels[D$kc] else NULL, k = D$k, link = if (D$k > 0L) D$link else NULL,
    coefficients = coefficients, theta = theta, vcov = Vn, vcov_internal = V, se = se, nboot = if (se == "bootstrap") nboot else NULL, boot = boot,
    loglik = if (firth) -fn$nll0(theta) else -op$value, lli = lli, n = D$n, nobs_weighted = sum(D$w), df = D$K,
    fitted.values = P, prob_outcome = Pout, phi = if (D$k > 0L) stats::setNames(pr$phi, D$case_id) else NULL,
    regime = if (length(D$kc) > 1L) { Rg <- pr$regime; dimnames(Rg) <- list(D$case_id, c("choice", D$levels[D$kc])); Rg } else NULL,
    y = D$y, weights = D$w, cluster = D$cluster, data = D$data, case_id = D$case_id,
    converged = op$converged, max_grad = op$max_grad, boundary = boundary, ill_conditioned = ill, singular = used_ginv, rcond = rc,
    diverging = D$par_names[diverge], phi_boundary = isTRUE(phi_boundary), saddle = identical(sad$type, "saddle"),
    start_logliks = -vals, screen_logliks = if (is.null(screen_vals)) NULL else -screen_vals, n_modes = n_modes,
    search = search$table, search_theta = search$theta, crosscheck = crosscheck,
    par_names = D$par_names, blocks = D$blocks, idx = D$idx, scale = scale,
    fe = D$fe, fe_info = D$fe_info, id = D$id, alt = D$alt, n_dropped_rows = D$n_dropped_rows,
    design = D), class = c(if (D$k > 0L) "imnl" else if (D$family == "mnp") "mnprobit" else D$family, "choicemodel"))
  ## a penalized fit: `loglik` is the ordinary log-likelihood at the penalized
  ## estimates, `loglik_penalized` the maximized objective
  if (firth) {
    if (isTRUE(D$firth)) obj$firth <- TRUE else obj$prior <- D$prior
    obj$loglik_penalized <- -op$value; obj$separation <- separation
  }
  if (isTRUE(attr(V, "prior_clipped")) || isTRUE(attr(Vd, "prior_clipped")))
    warning("The curvature of the log posterior is not positive definite at the mode (a t or Cauchy prior far out in its tail, as on a separated coefficient): ",
            "the covariance leaves out the prior's negative curvature, and the standard errors understate the uncertainty there. Use profile intervals (confint()).", call. = FALSE)
  obj
}

## the separation rule of the fit diagnostics (see .pc_fit_design): a
## coefficient whose one-SD effect exceeds 15, or a 0/1 column whose own
## coefficient exceeds 15, and (with standard errors) whose |z| is below 1.96
.pc_diverge <- function(theta, D, scale, se = NULL) {
  is01 <- function(M) if (ncol(M)) apply(M, 2L, function(v) all(v == 0 | v == 1)) else logical(0)
  avv <- as.vector(t(D$avail)) == 1L
  binary <- c(rep(if (D$p) is01(D$X) else logical(0), D$J - 1L), if (D$q) is01(D$W[avv, , drop = FALSE]) else logical(0),
              if (D$r) rep(is01(D$Z), max(1L, length(D$kc))) else logical(0), rep(FALSE, length(D$idx$str)))
  big <- abs(theta * scale) > 15 | (binary & abs(theta) > 15)
  ## along a diverging coefficient the likelihood flattens, so its Wald
  ## statistic tends to zero (Hauck and Donner 1977); a large coefficient
  ## estimated with precision is not diverging (the square and cube of a raw
  ## count of years, as in Liou 2024, have one-SD effects of 15-20 at |z| > 4,
  ## their size the collinearity of the polynomial). With the inverse-information
  ## standard errors at hand (not when the information is singular), a large
  ## coefficient counts as diverging only when its |z| is below 1.96.
  if (!is.null(se) && length(se) == length(theta)) big <- big & !(is.finite(se) & se > 0 & abs(theta) / se >= stats::qnorm(0.975))
  big
}

## multistart in two stages: every start is screened by a short BFGS run and
## then polished to convergence from its screened point (starts whose screens
## end at the same point are polished once). Polishing only the best few
## screened starts made the result depend on the screen: a start that screens
## lower can polish higher, so adding starts could lower the maximum found (on
## Huth et al. 2011 and Koyama 2025 in the replication stress test).
.pc_multistart <- function(fn, scale, starts, maxit, reltol) {
  screen <- lapply(starts, function(s) .pc_optimize(fn, s, scale, min(maxit, 150L), max(reltol, 1e-6), polish = FALSE))
  sv <- vapply(screen, `[[`, numeric(1), "value")
  keep <- is.finite(sv)
  for (i in which(keep)) for (j in which(keep[seq_len(i - 1L)]))
    if (abs(sv[i] - sv[j]) < 1e-6 && max(abs((screen[[i]]$par - screen[[j]]$par) * scale)) < 1e-4) { keep[i] <- FALSE; break }
  if (!any(keep)) keep[] <- TRUE
  fits <- lapply(screen[keep], function(s) .pc_optimize(fn, s$par, scale, maxit, reltol))
  vals <- vapply(fits, `[[`, numeric(1), "value")
  ## `ends`: every polished fit, for the record of the search's ends
  ## (.pc_search_record); `best` indexes the one reported
  list(op = fits[[which.min(vals)]], vals = vals, conv = vapply(fits, function(f) isTRUE(f$converged), logical(1)), screen_vals = sv,
       ends = fits, best = which.min(vals))
}

## challengers of an inflated fit: the incumbent's split coefficients halved,
## enlarged by half, and doubled, repeated from a new incumbent while one of
## them improves (at most three rounds). They start from the incumbent alone,
## so where the split coefficients diverge (many flat plateaus) the plateau
## reached depends on which polished fit is best (stress finding F19; the
## author declined challenging the best three, 2026-10-21)
.pc_challenge <- function(fn, D, scale, ms, maxit, reltol) {
  for (round in 1:3) {
    improved <- FALSE
    for (f in c(0.5, 1.5, 2)) {
      s2 <- ms$op$par; s2[D$idx$delta] <- s2[D$idx$delta] * f
      ch <- .pc_optimize(fn, s2, scale, maxit, reltol)
      better <- is.finite(ch$value) && ch$value < ms$op$value - 1e-8
      if (isTRUE(ch$converged) || better) { ms$vals <- c(ms$vals, ch$value); ms$conv <- c(ms$conv, isTRUE(ch$converged)); ms$ends <- c(ms$ends, list(ch)) }
      if (better) { ms$op <- ch; ms$best <- length(ms$ends); improved <- TRUE }
    }
    if (!improved) break
  }
  ms
}

## the split-link cross-check of imnl() (one inflated alternative): the same
## search under the other link, its split coefficients rescaled to this link
## and polished here; the higher maximum is kept. The logistic and normal
## densities at zero are 1/4 and 1/sqrt(2 pi), so a logit coefficient is about
## 1.6 times the probit one where the split probabilities are near one half.
## The two links' surfaces differ: on Glasgow (2001) the logit's own search
## stops 4.9 units below the maximum it reaches from the probit's estimates.
.pc_crosscheck <- function(D, fn, scale, ms, maxit, reltol) {
  D2 <- D; D2$link <- if (D$link == "logit") "probit" else "logit"; D2$link_code <- if (D2$link == "logit") 1L else 0L
  fn2 <- .pc_fns(D2)
  m2 <- .pc_challenge(fn2, D2, scale, .pc_multistart(fn2, scale, .pc_starts(D2, fn2, scale, maxit, reltol), maxit, reltol), maxit, reltol)
  s <- m2$op$par; s[D$idx$delta] <- s[D$idx$delta] * (if (D$link == "logit") 1.6 else 1 / 1.6)
  cx <- .pc_optimize(fn, s, scale, maxit, reltol)
  improved <- is.finite(cx$value) && cx$value < ms$op$value - 1e-8
  if (isTRUE(cx$converged) || improved) { ms$vals <- c(ms$vals, cx$value); ms$conv <- c(ms$conv, isTRUE(cx$converged)); ms$ends <- c(ms$ends, list(cx)) }
  ms$crosscheck <- list(link = D2$link, loglik_other_link = -m2$op$value, loglik_own_search = -ms$op$value, loglik_from_other = -cx$value, improved = improved)
  if (improved) { ms$op <- cx; ms$best <- length(ms$ends) }
  ms
}

## the split equation of an inflated fit at one regime, over the cases where
## the inflated alternative can be chosen (elsewhere phi is 1 by convention:
## regime 0 is impossible without k in the choice set); with several captive
## alternatives a captive share below 1e-6 wherever the alternative is
## available has collapsed (the model without it), and a choice regime of share
## below 1e-6 everywhere is the degenerate limit. `pr` holds the fit's
## probabilities (fn$probs)
.pc_phi_boundary <- function(D, pr) {
  ka <- if (D$k > 0L) D$avail[, D$k] == 1L else logical(0)
  phi_boundary <- D$k > 0L && any(ka) && (all(pr$phi[ka] > 1 - 1e-6) || all(pr$phi[ka] < 1e-6))
  collapsed <- character(0)
  if (length(D$kc) > 1L) {
    col <- vapply(seq_along(D$kc), function(t) { kt <- D$avail[, D$kc[t]] == 1L; any(kt) && all(pr$regime[kt, 1L + t] < 1e-6) }, logical(1))
    collapsed <- D$levels[D$kc][col]
    anyc <- rowSums(D$avail[, D$kc, drop = FALSE]) > 0L
    phi_boundary <- any(col) || (any(anyc) && all(pr$phi[anyc] < 1e-6))
  }
  list(phi_boundary = phi_boundary, collapsed = collapsed)
}

## what puts a point at a boundary, in words for the record of the search's
## ends ("" when nothing does): structural log-parameters past the flag
## (positions `str_out`), a collapsed split equation (.pc_phi_boundary), and
## diverging coefficients
.pc_boundary_detail <- function(D, diverge, str_out, pb) {
  parts <- character(0)
  if (length(str_out)) parts <- c(parts, paste0(paste(D$par_names[str_out], collapse = ", "), " at zero or infinity"))
  if (length(pb$collapsed)) parts <- c(parts, paste0("captivity to ", paste(pb$collapsed, collapse = ", "), " collapsed"))
  else if (isTRUE(pb$phi_boundary)) parts <- c(parts, "the split equation puts every case in one regime")
  if (any(diverge)) parts <- c(parts, paste0("diverging: ", paste(utils::head(D$par_names[diverge], 4L), collapse = ", "), if (sum(diverge) > 4L) ", ..." else ""))
  paste(parts, collapse = "; ")
}

## The type of an end of the search other than the reported fit, by the rules
## the reported fit's flags follow (.pc_fit_design): "boundary" (a diverging
## coefficient by the exact check of separation() or the size-and-Wald rule,
## a structural log-parameter beyond 7, or a split equation at one regime),
## "saddle" (.pc_saddle on the normalized Hessian over the parameters with
## curvature), "not converged", or "maximum". Returns the type, its details,
## and the inverse-information standard errors (NULL when the information is
## singular), which the grouping of the ends uses.
.pc_end_type <- function(D, fn, scale, w, end, reltol) {
  theta <- end$par
  Vd <- tryCatch(.pc_vcov(theta, fn, "analytic", NULL, scale, w), error = function(e) { .pc_pass_limit(e); NULL })
  ginv <- !is.null(Vd) && isTRUE(attr(Vd, "ginv"))
  se <- if (is.null(Vd) || ginv) NULL else attr(Vd, "info_se")
  sepchk <- .pc_separation_fit(D, theta)
  diverge <- if (!is.null(sepchk)) D$par_names %in% sepchk$terms else .pc_diverge(theta, D, scale, se = se)
  str_log <- .pc_str_log(D)
  pb <- if (D$k > 0L) .pc_phi_boundary(D, fn$probs(theta)) else list(phi_boundary = FALSE, collapsed = character(0))
  bd <- .pc_boundary_detail(D, diverge, str_log[abs(theta[str_log]) > 7], pb)
  if (nzchar(bd)) return(list(type = "boundary", detail = bd, se = se))
  if (is.null(Vd)) return(list(type = "not converged", detail = "the Hessian could not be computed", se = NULL))
  hdiag <- attr(Vd, "hdiag"); hrow <- attr(Vd, "hrow"); noise_m <- attr(Vd, "noise")
  keep <- !(hrow <= 1e-8 * max(abs(hdiag)))                           # parameters with curvature (as `flat` in .pc_fit_design)
  noise <- if (is.null(noise_m) || !any(keep)) 0 else max(noise_m[keep, keep, drop = FALSE])
  sad <- .pc_saddle(fn$nll, theta, scale, attr(Vd, "hn"), hdiag, keep, end$value, reltol, max(1e-8, 10 * noise))
  if (sad$type == "saddle") return(list(type = "saddle", detail = sprintf("smallest normalized eigenvalue %.2g", sad$min_eig), se = se))
  if (!isTRUE(end$converged)) return(list(type = "not converged", detail = sprintf("largest scaled gradient %.2g", end$max_grad), se = se))
  rc <- attr(Vd, "rcond")
  list(type = "maximum", detail = if (ginv || (is.finite(rc) && rc < 1e-9)) "ill-conditioned information" else "", se = se)
}

## The record of the search's ends (object$search and object$search_theta):
## the reported fit and every polished fit of the search (the multistart's,
## the challengers' and the cross-check's that are kept, and before a saddle
## step the original search's), grouped into distinct points. Two ends are the
## same point when their log-likelihoods differ by less than 0.01 and either
## their estimates by less than a tenth of the first's inverse-information
## standard error (where it has none, by less than 0.01 in the standardized
## coordinates theta * scale) or every fitted probability by less than 1e-4.
## The second is the data's own test: ends at other estimates with the same
## probabilities lie on a ridge of non-identification (the saturated designs
## of DiGiuseppe et al. 2024, where the estimates differ by several units and
## the probabilities agree to 1e-8), and the point's `detail` says how many
## ends joined it that way. The end the reported fit came from (`best`, which
## the boundary step can have moved) is the reported point. The ends are taken
## in decreasing order of log-likelihood, each point classified at its
## highest end (.pc_end_type; the reported point by its own flags, `rep`).
## Returns `table` (loglik, type, reached = the number of ends at the point,
## reported, detail), ordered by log-likelihood, and `theta`, the estimates
## at each point on the internal scale of `start`, one row per row of `table`.
.pc_search_record <- function(D, fn, scale, w, ends, best, op, rep, reltol) {
  pts <- list(list(par = op$par, ll = -op$value, se = rep$se, type = rep$type, detail = rep$detail, n = 0L, ridge = 0L, reported = TRUE))
  same_est <- function(par, p) {
    d <- abs(par - p$par)
    ok <- if (is.null(p$se)) rep(FALSE, length(d)) else is.finite(p$se) & p$se > 0
    all(d[ok] < 0.1 * p$se[ok]) && all(d[!ok] * scale[!ok] < 0.01)
  }
  ## the fitted probabilities, computed once per end or point and only where the
  ## estimates differ (list(NULL) when they cannot be computed)
  probs <- function(par) list(tryCatch(fn$probs(par)$prob, error = function(e) { .pc_pass_limit(e); NULL }))
  same_prob <- function(P1, P2) !is.null(P1) && !is.null(P2) && identical(dim(P1), dim(P2)) && identical(is.na(P1), is.na(P2)) &&
    all(abs(P1 - P2) < 1e-4, na.rm = TRUE)
  vals <- vapply(ends, function(e) if (is.finite(e$value)) e$value else Inf, numeric(1))
  for (i in order(vals)) {
    e <- ends[[i]]
    if (i == best) { pts[[1L]]$n <- pts[[1L]]$n + 1L; next }
    if (!is.finite(e$value)) next
    near <- which(vapply(pts, function(p) abs(-e$value - p$ll) < 0.01, logical(1)))
    hit <- 0L; by_prob <- FALSE; pe <- NULL
    for (j in near) if (same_est(e$par, pts[[j]])) { hit <- j; break }
    if (!hit) for (j in near) {
      if (is.null(pe)) pe <- probs(e$par)
      if (is.null(pts[[j]]$prob)) pts[[j]]$prob <- probs(pts[[j]]$par)
      if (same_prob(pe[[1L]], pts[[j]]$prob[[1L]])) { hit <- j; by_prob <- TRUE; break }
    }
    if (hit) { pts[[hit]]$n <- pts[[hit]]$n + 1L; pts[[hit]]$ridge <- pts[[hit]]$ridge + by_prob; next }
    ty <- .pc_end_type(D, fn, scale, w, e, reltol)
    pts[[length(pts) + 1L]] <- list(par = e$par, ll = -e$value, se = ty$se, type = ty$type, detail = ty$detail, n = 1L, ridge = 0L, reported = FALSE, prob = pe)
  }
  o <- order(-vapply(pts, `[[`, numeric(1), "ll"), !vapply(pts, `[[`, logical(1), "reported"))
  pts <- pts[o]
  detail <- vapply(pts, function(p) {
    if (p$ridge == 0L) return(p$detail)
    r <- sprintf("%d %s at other estimates with the same fitted probabilities: the estimates are not identified", p$ridge, if (p$ridge == 1L) "end" else "ends")
    if (nzchar(p$detail)) paste0(p$detail, "; ", r) else r
  }, character(1))
  table <- data.frame(loglik = vapply(pts, `[[`, numeric(1), "ll"), type = vapply(pts, `[[`, character(1), "type"),
                      reached = vapply(pts, `[[`, integer(1), "n"), reported = vapply(pts, `[[`, logical(1), "reported"),
                      detail = detail, stringsAsFactors = FALSE)
  theta <- matrix(unlist(lapply(pts, `[[`, "par")), nrow = length(pts), byrow = TRUE, dimnames = list(NULL, D$par_names))
  list(table = table, theta = theta)
}

## the message of an inflated fit about the search's maxima: when it reached
## more than one, or when the reported fit (the highest end) is not one but a
## lower end is
.pc_search_message <- function(tab, n_modes, n_starts) {
  mx <- tab$loglik[tab$type == "maximum"]
  fmt <- function(v) { s <- sprintf("%.1f", v); if (anyDuplicated(s)) s <- sprintf("%.3f", v); paste(s, collapse = ", ") }
  rtype <- tab$type[tab$reported]
  other <- sum(tab$type != "maximum" & !tab$reported)
  if (rtype == "maximum" && n_modes > 1L)
    message("The inflated likelihood has ", n_modes, " distinct local maxima among the ends of the search from ", n_starts, " starts (log-likelihoods about ",
            fmt(mx), "); the highest is reported.",
            if (other) paste0(" The search also ended at ", other, if (other == 1L) " point that is" else " points that are",
                              " not a maximum (a boundary, a saddle, or no convergence).") else "",
            " See object$search.")
  else if (rtype != "maximum" && n_modes > 0L)
    message("The reported fit, the highest end of the search, is ",
            switch(rtype, boundary = "at a boundary", saddle = "a saddle", "not converged"), " (see the warning); the search also reached ",
            n_modes, if (n_modes == 1L) " local maximum" else " distinct local maxima", " below it (log-likelihood", if (n_modes > 1L) "s", " about ", fmt(mx),
            "). object$search lists the ends of the search, and object$search_theta holds their estimates (for 'start').")
  invisible()
}

## a restricted fit whose maximum enters a likelihood-ratio statistic beside a
## user's fit (term_test(type = "lr")): fitted as a user's fit is, with the
## saddle step, so that both sides of the statistic are maxima (local review
## 2026-11-06, finding 1). Its messages are dropped; a warning is kept when it
## did not converge (unless its design is separated, where the statistic
## takes the supremum) or when it still stops at a saddle. Returns the
## internal parameter vector and the log-likelihood, as .pc_refit() does.
.pc_refit_checked <- function(D, maxit = 1000L, reltol = 1e-10) {
  if (D$K == 0L) return(.pc_refit_quiet(D, maxit = maxit, reltol = reltol))   # no parameters: evaluated, not fitted
  f <- withCallingHandlers(.pc_fit_design(D, se = "none", maxit = maxit, reltol = reltol),
                           message = function(m) invokeRestart("muffleMessage"), warning = function(w) invokeRestart("muffleWarning"))
  if (!isTRUE(f$converged)) {
    sep <- .pc_sep_scope(D) && isTRUE(tryCatch(.pc_separation(D, f$theta)$separated, error = function(e) { .pc_pass_limit(e); FALSE }))
    if (!sep) warning("A restricted refit did not converge (largest scaled gradient ", format(f$max_grad, digits = 3), ").", call. = FALSE)
  }
  if (isTRUE(f$saddle)) warning("A restricted refit stopped at a saddle that the saddle step did not leave; the likelihood-ratio statistic may be too large.", call. = FALSE)
  list(coefficients = f$theta, loglik = f$loglik, converged = f$converged)
}

## refit the same specification on another design (bootstrap, tests): returns
## the parameter vector, the log-likelihood, and optionally the covariance
.pc_refit <- function(D, start = NULL, maxit = 1000L, reltol = 1e-10, vcov = FALSE, se = "analytic") {
  fn <- .pc_fns(D)
  scale <- .pc_scales(D)
  starts <- if (!is.null(start)) list(start) else .pc_starts(D, fn, scale, maxit, reltol)
  fits <- lapply(starts, function(s) .pc_optimize(fn, s, scale, maxit, reltol, newton_first = .pc_concave(D)))
  op <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  op <- .pc_boundary_step(fn, D, op, scale, maxit, reltol)
  theta <- op$par; names(theta) <- D$par_names
  if (!op$converged) warning("A refit did not converge (largest scaled gradient ", format(op$max_grad, digits = 3), ").", call. = FALSE)
  out <- list(coefficients = theta, loglik = -op$value, converged = op$converged, fn = fn, scale = scale)
  if (vcov) {
    obs <- .pc_obs(D)
    V <- .pc_vcov(theta, fn, se, obs$cluster, scale, obs$w)
    if (isTRUE(attr(V, "ginv"))) warning("A refit's information matrix is singular (a generalized inverse was used); some of its parameters are not identified by the data.", call. = FALSE)
    if (!is.null(V)) { V <- .pc_strip_vcov(V); dimnames(V) <- list(D$par_names, D$par_names) }
    out$vcov <- V
  }
  out
}
