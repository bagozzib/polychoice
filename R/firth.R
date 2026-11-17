## ---------------------------------------------------------------------------
## Firth's penalty (`firth = TRUE` in mnl() and rologit()): the log-likelihood
## penalized by half the log-determinant of the information, Jeffreys' prior.
## Definitions fixed in DESIGN.md, section 1, Phase 4; kernel src/firth.cpp.
## ---------------------------------------------------------------------------

## the logit families whose maximum-likelihood log-likelihood is globally
## concave, solved by damped Newton steps first; the penalized objective's
## concavity is not known, so a penalized fit starts with BFGS
.pc_concave <- function(D) D$family %in% c("mnl", "rologit") && D$k == 0L && !.pc_penalized(D)

## The penalized objective of a design, in the engine's form: nll is the
## negative penalized log-likelihood -(l + 1/2 log|I|), grad its gradient (the
## gradient of the negative log-likelihood minus the penalty's), scores the
## per-case adjusted scores (the score plus the case's share of the penalty's
## gradient; with the weights they sum to zero at the estimate), and `info` the
## information I itself, which the covariance inverts (the model's information
## at the penalized estimates, as brglm2, logistf and coxphf report). The
## Hessian of the penalized objective is left to central differences of its
## gradient. `nll0` keeps the unpenalized negative log-likelihood. I is
## memoised at the last parameter vector: the optimizer asks for the objective
## and the gradient at the same point. The log-determinant and the inverse are
## those of the equilibrated information I_s = I / (d d'), d = sqrt(diag(I)):
## log|I| = log|I_s| + 2 sum log d and I^{-1} = I_s^{-1} / (d d'). In theta
## units I inherits the covariates' scales (Poast's 2007 IOC bids: 0.3 to
## 1.7e8, a reciprocal condition number of 1e-18 where I_s has 0.02), solve()
## refuses it, and a generalized inverse would truncate it.
.pc_firth_fns <- function(f, D) {
  ro <- identical(D$family, "rologit")
  a <- list(D$X, D$W, D$off, D$avail, if (ro) D$str$rank else D$avail, if (ro) 1L else 0L, D$J, D$base - 1L,
            if (ro) as.integer(D$str$split %||% 0L) else 0L)
  w <- D$w
  memo <- new.env(parent = emptyenv())
  info <- function(th) {
    th <- as.numeric(th)
    if (!is.null(memo$th) && identical(memo$th, th)) return(memo$I)
    I <- do.call(pc_firth_info_cpp, c(list(th), a, list(w)))
    memo$th <- th; memo$I <- I
    I
  }
  eq <- function(I) { d <- sqrt(diag(I)); d[!(is.finite(d) & d > 0)] <- 1; list(Is = I / tcrossprod(d), d = d) }
  logdet <- function(th) {
    I <- info(th)
    if (!all(is.finite(I)) || !all(diag(I) > 0)) return(-Inf)       # a coefficient without information: I is singular
    e <- eq(I)
    dt <- tryCatch(determinant(e$Is, logarithm = TRUE), error = function(e) { .pc_pass_limit(e); NULL })
    if (is.null(dt) || dt$sign <= 0 || !is.finite(dt$modulus)) -Inf else as.numeric(dt$modulus) + 2 * sum(log(e$d))
  }
  adj <- function(th) {
    e <- eq(info(th))
    Ii <- .pc_inverse(e$Is)$inv / tcrossprod(e$d)
    do.call(pc_firth_adj_cpp, c(list(as.numeric(th)), a, list(Ii)))
  }
  list(lli = f$lli,
       nll = function(th) f$nll(th) - 0.5 * logdet(th),
       grad = function(th) f$grad(th) - colSums(w * adj(th)),
       hess = NULL, info = info,
       scores = function(th) f$scores(th) + adj(th),
       probs = f$probs, nll0 = f$nll, penalty = function(th) 0.5 * logdet(th))
}

## The second start of a penalized fit without a user start (DESIGN.md, Phase
## 4b). The penalized likelihood can have two local maxima (a high-leverage
## covariate: on the research notes' separated binary design the run from
## zero finds the lower one, on their overlapping one a run from the
## maximum-likelihood estimate does), so the fit also runs from the endpoint
## of the maximum-likelihood fit of the same design (concave, Newton first),
## halved towards zero until the penalized objective is finite there (far out
## along a separating direction the information underflows), and keeps the
## higher maximum. `op` is the run from zero. Returns the kept run, the
## objective values of the runs, the maximum-likelihood fit, and `modes`, the
## two penalized log-likelihoods when both runs converged to distinct maxima
## (points more than 1e-3 apart in the standardized coordinates whose values
## differ by more than 1e-6 (1 + |l*|)).
.pc_firth_second <- function(D, fn, scale, maxit, reltol, op) {
  D0 <- D; D0$firth <- FALSE; D0$prior <- NULL
  ml <- tryCatch(suppressWarnings(.pc_optimize(.pc_fns(D0), numeric(D0$K), scale, maxit, reltol, newton_first = TRUE)),
                 error = function(e) { .pc_pass_limit(e); NULL })
  out <- list(op = op, vals = op$value, ml = ml, modes = NULL)
  if (is.null(ml) || !all(is.finite(ml$par))) return(out)
  s <- ml$par
  for (h in 0:60) {
    v <- tryCatch(fn$nll(s), error = function(e) { .pc_pass_limit(e); NA_real_ })
    if (is.finite(v)) break
    s <- s / 2
  }
  if (!is.finite(v)) return(out)
  op2 <- tryCatch(.pc_optimize(fn, s, scale, maxit, reltol), error = function(e) { .pc_pass_limit(e); NULL })
  ## a run from far out along a separating direction can stop short of a
  ## maximum (the information is close to singular there). It is then run
  ## again from the point of the segment from that endpoint to zero where
  ## the penalized objective is highest, and the better run is kept
  if (!is.null(op2) && is.finite(op2$value) && !isTRUE(op2$converged)) {
    tt <- 2^-(1:30)
    vv <- vapply(tt, function(t) tryCatch(fn$nll(t * s), error = function(e) { .pc_pass_limit(e); NA_real_ }), numeric(1))
    if (any(is.finite(vv))) {
      op3 <- tryCatch(.pc_optimize(fn, tt[which.min(vv)] * s, scale, maxit, reltol), error = function(e) { .pc_pass_limit(e); NULL })
      if (!is.null(op3) && is.finite(op3$value) && op3$value <= op2$value) op2 <- op3
    }
  }
  if (is.null(op2) || !is.finite(op2$value)) return(out)
  out$vals <- c(op$value, op2$value)
  if (isTRUE(op$converged) && isTRUE(op2$converged) && max(abs((op$par - op2$par) * scale)) > 1e-3 &&
      abs(op$value - op2$value) > 1e-6 * (1 + abs(op$value)))
    out$modes <- -c(op$value, op2$value)
  if (!is.finite(op$value) || op2$value < op$value) out$op <- op2
  out
}

## the separation report of a penalized fit when the exact check
## (.pc_separation) is unavailable (a solver failure): the maximum-likelihood
## fit of the same design (`op`, from .pc_firth_second, or refitted: concave,
## Newton first) and the engine's separation rule on
## it, keeping the coefficients whose maximum-likelihood estimate also lies
## more than 3 of the penalized fit's standard errors from the penalized
## estimate (theta, se). Along a separating direction the maximum-likelihood
## estimate runs off while the penalized one stays put; a large constant that
## the rule flags because it is imprecise (the constants of a model with an
## uncentered year of birth, as in Cook, Niehaus and Zuhlke's models, where
## the two estimates differ by 0.1 standard errors) is not separated. Returns
## the names of the separated coefficients.
.pc_firth_separation <- function(D, scale, maxit, reltol, theta = NULL, se = NULL, op = NULL) {
  D0 <- D; D0$firth <- FALSE; D0$prior <- NULL
  if (is.null(op)) op <- tryCatch(suppressWarnings(.pc_optimize(.pc_fns(D0), numeric(D0$K), scale, maxit, reltol, newton_first = TRUE)),
                                  error = function(e) { .pc_pass_limit(e); NULL })
  if (is.null(op) || !all(is.finite(op$par))) return(character(0))
  dv <- .pc_diverge(op$par, D0, scale, se = .pc_info_se(op$par, D0, scale))
  if (!is.null(theta) && length(se) == length(theta))
    dv <- dv & !(is.finite(se) & abs(op$par - theta) <= 3 * se)
  D$par_names[dv]
}

## The maximum of a penalized objective (the full model's penalty) under the
## linear restriction Rm theta = 0: in the standardized coordinates
## ts = theta * scale the restriction is (Rm / scale) ts = 0, whose solutions
## are ts = N phi with N an orthonormal basis of that null space (from the QR
## decomposition), so the search runs over phi with unit scales. Returns the
## restricted maximum's objective value (negative penalized log-likelihood),
## its point, and convergence. The search runs from `start` projected onto the
## restriction and from zero (the fit's own start), and keeps the higher
## maximum: a start built from the full estimates can sit far from the
## restricted maximum (zeroing the year-of-birth coefficients of a model whose
## constants absorb a year near 1980 leaves those constants at 165, where the
## restricted maximum has them near 0).
.pc_firth_restricted <- function(fn, scale, Rm, start, maxit = 1000L, reltol = 1e-10) {
  Rs <- sweep(Rm, 2L, scale, "/")
  qrR <- qr(t(Rs))
  N <- qr.Q(qrR, complete = TRUE)[, -seq_len(qrR$rank), drop = FALSE]
  th_of <- function(phi) as.vector(N %*% phi) / scale
  fnr <- list(nll = function(phi) fn$nll(th_of(phi)),
              grad = function(phi) as.vector(crossprod(N, fn$grad(th_of(phi)) / scale)),
              info = function(phi) crossprod(N, (fn$info(th_of(phi)) / tcrossprod(scale)) %*% N))     # the scoring's curvature in phi
  run <- function(phi0) tryCatch(.pc_optimize(fnr, phi0, rep(1, ncol(N)), maxit, reltol), error = function(e) { .pc_pass_limit(e); NULL })
  fits <- Filter(function(op) !is.null(op) && is.finite(op$value), list(run(as.vector(crossprod(N, start * scale))), run(numeric(ncol(N)))))
  if (!length(fits)) return(list(value = NA_real_, theta = rep(NA_real_, length(scale)), converged = FALSE))
  op <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  list(value = op$value, theta = th_of(op$par), converged = isTRUE(op$converged))
}

## the penalized likelihood-ratio statistic 2 [l*(theta_hat) - max_{Rm theta = 0} l*]
## of a penalized fit (Heinze and Schemper 2001), restricted from `start`
.pc_firth_lr <- function(object, Rm, start, fn = .pc_fns(object$design), scale = object$scale) {
  rr <- .pc_firth_restricted(fn, scale, Rm, start)
  if (!rr$converged) warning("The restricted penalized fit did not converge; the penalized likelihood-ratio statistic may be inaccurate.", call. = FALSE)
  2 * (object$loglik_penalized + rr$value)
}

## Profile penalized-likelihood intervals (Heinze and Schemper 2002): the values
## c of a coefficient with 2 [l*(theta_hat) - max_{theta_r = c} l*(theta)] at most
## qchisq(level, 1), the full model's penalty in both terms. Each limit is
## bracketed by stepping out from the estimate in multiples of the Wald
## standard error (1, 2, 4, ... up to 50) and found by uniroot(); NA with a
## warning when no limit lies within 50 standard errors.
##
## The restricted maximum at a point is first taken along a path of warm
## starts from the estimate. The restricted penalized likelihood can have
## several local maxima (few cases per coefficient, separation), and a path
## can follow a lower one, so a limit is checked by maximizing there afresh
## from the estimate and from zero. Where that finds a higher maximum the
## side is searched again with every point's maximum the highest of the
## three starts. A limit is reported when the highest restricted maximum
## found at it converged and reproduces the quantile; otherwise it is NA with a
## warning, since a number there would be a point on some path and not a
## limit (DESIGN.md Phase 4g).
.pc_firth_profile <- function(object, parm, level) {
  D <- object$design; fn <- .pc_fns(D); scale <- object$scale
  th <- object$theta; pn <- names(object$coefficients)
  q <- stats::qchisq(level, 1)
  se <- if (is.null(object$vcov_internal)) rep(NA_real_, length(th)) else sqrt(pmax(diag(object$vcov_internal), 0))
  missed <- character(0); unverified <- character(0)
  one <- function(nm) {
    r <- match(nm, pn); free <- setdiff(seq_along(th), r)
    warm <- th
    ## the restricted maximizations at cval: from the previous point's solution
    ## (`warm`), and with `cold` also from the estimate and from zero. Far out
    ## along a separated coefficient the information underflows, the penalty's
    ## log-determinant is -Inf and the penalized log-likelihood -Inf: that
    ## point is beyond the limit (a large positive value), not a failure of
    ## the search (DESIGN.md Phase 4d)
    solve_at <- function(cval, starts = c("warm", "cold")) {
      full <- function(b) { t <- numeric(length(th)); t[free] <- b; t[r] <- cval; t }
      fnr <- list(nll = function(b) fn$nll(full(b)), grad = function(b) fn$grad(full(b))[free],
                  info = function(b) fn$info(full(b))[free, free, drop = FALSE])
      s0 <- c(if ("warm" %in% starts) list(warm[free]), if ("cold" %in% starts) list(th[free], numeric(length(free))))
      ops <- lapply(s0, function(s) tryCatch(.pc_optimize(fnr, s, scale[free], 1000L, 1e-10), error = function(e) { .pc_pass_limit(e); list(value = Inf, converged = FALSE) }))
      vals <- vapply(ops, function(o) if (is.finite(o$value)) o$value else Inf, numeric(1))
      op <- ops[[which.min(vals)]]
      if ("warm" %in% starts && is.finite(op$value)) warm <<- full(op$par)
      v <- 2 * (object$loglik_penalized + min(vals)) - q
      list(v = if (!is.finite(v)) 1e10 else v, converged = isTRUE(op$converged))
    }
    step <- if (is.finite(se[r]) && se[r] > 0) se[r] else max(1, abs(th[r]))
    ## one side's limit: list(root, at) with `at` the deviance and convergence
    ## of the restricted maximum there, or NULL when no limit is bracketed
    lim <- function(side, starts) {
      dev <- function(cval) solve_at(cval, starts)$v
      lo <- th[r]; g_lo <- -q; warm <<- th
      for (k in c(1, 2, 4, 8, 16, 32, 50)) {
        hi <- th[r] + side * k * step
        g_hi <- tryCatch(dev(hi), error = function(e) { .pc_pass_limit(e); NA_real_ })
        if (!is.finite(g_hi)) break
        if (g_hi > 0) {
          ## the root to 1e-6 of a standard error, in the coefficient's own
          ## units (an absolute tolerance would exceed the standard error of a
          ## covariate measured in large units; local review 2026-11-06,
          ## finding 13)
          root <- tryCatch(stats::uniroot(function(cv) dev(cv), sort(c(lo, hi)), f.lower = if (side > 0) g_lo else g_hi,
                                          f.upper = if (side > 0) g_hi else g_lo, tol = 1e-6 * step)$root,
                           error = function(e) { .pc_pass_limit(e); NA_real_ })
          if (!is.finite(root)) break
          return(list(root = root, at = tryCatch(solve_at(root, starts), error = function(e) { .pc_pass_limit(e); NULL })))
        }
        lo <- hi; g_lo <- g_hi
      }
      NULL
    }
    ## pass 1 follows the warm path; at its limit the maximum is taken again
    ## from all three starts. A higher maximum there (the path was on a lower
    ## branch) sends the side to pass 2, where every point's maximum is the
    ## highest of the three starts. A limit is reported when the highest
    ## maximum found at it converged and reproduces the quantile; pass 2 has
    ## no further check, since its maxima are already the highest found
    side_limit <- function(side) {
      l <- lim(side, "warm")
      if (is.null(l)) { missed <<- c(missed, nm); return(NA_real_) }
      at <- tryCatch(solve_at(l$root, c("warm", "cold")), error = function(e) { .pc_pass_limit(e); NULL })
      if (!is.null(at) && isTRUE(at$converged) && abs(at$v) <= 1e-3) return(l$root)
      if (!is.null(at) && at$v < -1e-3) {
        l <- lim(side, c("warm", "cold"))
        if (!is.null(l) && !is.null(l$at) && isTRUE(l$at$converged) && abs(l$at$v) <= 1e-3) return(l$root)
      }
      unverified <<- c(unverified, nm)
      NA_real_
    }
    c(side_limit(-1), side_limit(1))
  }
  ci <- t(vapply(parm, one, numeric(2)))
  if (length(missed)) warning("No profile limit within 50 standard errors for ", paste(unique(missed), collapse = ", "),
                              ": the penalized likelihood is too flat there (NA).", call. = FALSE)
  if (length(unverified)) warning("A profile limit of ", paste(unique(unverified), collapse = ", "), " could not be verified and is not reported (NA): ",
                                  "at the point found, the restricted penalized maximum did not converge or does not reproduce the limit ",
                                  "(the restricted penalized likelihood has several local maxima).", call. = FALSE)
  ci
}

## a penalized design or fit: Firth's penalty or a prior (DESIGN.md Phase 4c).
## The penalized machinery (two starts, scoring, profile intervals, penalized
## likelihood ratios) serves both; only the penalty differs
.pc_penalized <- function(x) isTRUE(x$firth) || !is.null(x$prior)

## the penalty's name for printed output: "Firth", or e.g. "Cauchy(0, 2.5) prior"
.pc_penalty_label <- function(x) {
  if (isTRUE(x$firth)) return("Firth")
  pr <- x$prior
  if (is.null(pr)) return(NULL)
  fam <- switch(pr$family, normal = "normal", cauchy = "Cauchy", t = paste0("t(", format(pr$df), ")"))
  sc <- if (length(pr$scale) == 1L && is.null(names(pr$scale))) paste0("scale ", format(pr$scale)) else "scales set by coefficient"
  paste0(fam, " prior (", sc, if (isTRUE(pr$autoscale)) ", autoscaled" else "", ")")
}

## whether two penalized fits carry the same penalty (DESIGN.md Phase 4d):
## both Firth's, or priors of the same family, degrees of freedom, autoscaling
## and constants' scale whose scales agree on the coefficients the fits share
## (a scale set by name for a coefficient only the full model has is that
## coefficient's own, and the restricted model, which lacks it, ignores it)
.pc_same_penalty <- function(a, b) {
  if (isTRUE(a$firth) || isTRUE(b$firth)) return(isTRUE(a$firth) && isTRUE(b$firth))
  pa <- a$prior; pb <- b$prior
  if (is.null(pa) || is.null(pb)) return(FALSE)
  if (!identical(pa$family, pb$family) || !isTRUE(all.equal(as.numeric(pa$df), as.numeric(pb$df))) || !identical(isTRUE(pa$autoscale), isTRUE(pb$autoscale)) ||
      !isTRUE(all.equal(pa$intercept_scale, pb$intercept_scale))) return(FALSE)
  sh <- intersect(names(pa$s), names(pb$s))
  length(sh) > 0L && isTRUE(all.equal(unname(pa$s[sh]), unname(pb$s[sh]), tolerance = 1e-10))
}

## the functions that do not take a penalized fit
.pc_no_firth <- function(object, what, why = "its reference distribution is derived for maximum-likelihood estimates", penalized_forms = TRUE) {
  if (.pc_penalized(object))
    stop(what, " does not take a penalized fit (", if (isTRUE(object$firth)) "firth = TRUE" else "a prior", "): ", why,
         ". Use the maximum-likelihood fit", if (penalized_forms) ", or the penalized likelihood-ratio forms of term_test(), combine_test() and lr_test()" else "", ".", call. = FALSE)
  invisible(TRUE)
}
