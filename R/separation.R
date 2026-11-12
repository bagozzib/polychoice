## ---------------------------------------------------------------------------
## Separation in the conditional, multinomial and rank-ordered logit: an exact
## check by linear programming (DESIGN.md, section 1, Phase 4b). A finite
## maximum-likelihood estimate exists iff the structure vectors (the
## differences of the design rows of the alternatives each observation orders)
## overlap: no direction v != 0 with s'v >= 0 for all of them (Albert and
## Anderson 1984; Speckman, Lee and Sun 2009; Sablica, Hornik and Rusch 2026,
## Theorem 1, rankings and ties included).
## ---------------------------------------------------------------------------

## in scope: mnl() and rologit() designs without unit fixed effects or
## inflation, and not rank_test()'s split designs
.pc_sep_scope <- function(D) {
  D$family %in% c("mnl", "rologit") && D$k == 0L && (is.null(D$fe_map) || !nrow(D$fe_map)) &&
    !(identical(D$family, "rologit") && !identical(as.integer(D$str$split %||% 0L), 0L)) &&
    D$K == D$p * (D$J - 1L) + D$q
}

## the structure vectors of the comparisons P (case, above, below), in theta
## coordinates: the design row of the alternative above minus that of the one
## below, where the case-specific block of non-base alternative m is X_i for
## m and zero otherwise, followed by the generic columns W_ij
.pc_sep_struct <- function(D, P) {
  p <- D$p; q <- D$q; J <- D$J; nb <- setdiff(seq_len(J), D$base)
  S <- matrix(0, nrow(P), D$K)
  if (p) for (m in seq_along(nb)) {
    cols <- (m - 1L) * p + seq_len(p)
    a <- P[, 2L] == nb[m]; if (any(a)) S[a, cols] <- D$X[P[a, 1L], , drop = FALSE]
    b <- P[, 3L] == nb[m]; if (any(b)) S[b, cols] <- S[b, cols, drop = FALSE] - D$X[P[b, 1L], , drop = FALSE]
  }
  if (q) {
    Wa <- D$W[(P[, 1L] - 1L) * J + P[, 2L], , drop = FALSE]; Wb <- D$W[(P[, 1L] - 1L) * J + P[, 3L], , drop = FALSE]
    Sw <- Wa - Wb
    ## a difference within rounding is zero: 0.1 + 0.2 against 0.3, or a large
    ## level whose tie differs in the last digit (round 7, D-1). Rounding is
    ## measured against the larger of the operands and the column's typical
    ## size (the median of its nonzero absolute values): after centring, a
    ## tie's operands can be tiny while its rounding error is that of the
    ## uncentred values (local review 2026-11-06, finding 3), and a median,
    ## unlike the column's largest entry, is not moved by an outlier. The
    ## case-specific blocks copy design values and need no such rule.
    cs <- apply(abs(D$W), 2L, function(v) { v <- v[v > 0 & is.finite(v)]; if (length(v)) stats::median(v) else 0 })
    thr <- pmax(abs(Wa), abs(Wb), matrix(cs, nrow(Wa), q, byrow = TRUE))
    Sw[abs(Sw) <= 64 * .Machine$double.eps * thr] <- 0
    S[, p * (J - 1L) + seq_len(q)] <- Sw
  }
  S
}

## the ordered pairs (case, above, below) the data report: in mnl() the choice
## above every other available alternative; in rologit() each alternative of a
## rank group above each of the next group (ties as groups, the available
## unranked alternatives last); cases with zero weight excluded
.pc_sep_pairs <- function(D) {
  A <- D$avail == 1L; A[D$w <= 0, ] <- FALSE; keep <- which(rowSums(A) > 0L)
  if (D$family == "mnl") {
    ix <- which(A, arr.ind = TRUE)
    ix <- ix[ix[, 2L] != D$y[ix[, 1L]], , drop = FALSE]
    P <- cbind(ix[, 1L], D$y[ix[, 1L]], ix[, 2L])
    return(P[order(P[, 1L], P[, 3L]), , drop = FALSE])
  } else {
    R <- D$str$rank
    out <- lapply(keep, function(i) {
      av <- which(A[i, ]); r <- R[i, av]; lev <- sort(unique(r[r > 0L]))
      grp <- lapply(lev, function(g) av[r == g]); if (any(r == 0L)) grp <- c(grp, list(av[r == 0L]))
      if (length(grp) < 2L) return(NULL)
      do.call(rbind, lapply(seq_len(length(grp) - 1L), function(g) cbind(i, rep(grp[[g]], each = length(grp[[g + 1L]])), rep(grp[[g + 1L]], length(grp[[g]])))))
    })
  }
  P <- do.call(rbind, out)
  if (is.null(P)) matrix(0L, 0L, 3L) else P
}

## The linear program, in K variables: maximize the sum of s_i'v over the
## comparisons not yet found separable, subject to s_i'v >= 0 for every
## comparison and -1 <= v_k <= 1; the comparisons it makes positive (above
## tol) are separable, and the program is solved again over the rest until it
## adds none. A separating direction positive on a remaining comparison would
## make the objective positive, so every separable comparison is found; the
## sum of the solutions is positive on all of them (the cone of separating
## directions is closed under addition).
## A re-solve after the objective changes starts from the previous basis; where
## lp_solve fails from there (status 3, "unbounded", on Kardasheva 2013 in the
## stress battery, where a fresh model solves) the program is rebuilt and
## solved afresh once.
.pc_sep_lp <- function(S, tol = 1e-9) {
  N <- nrow(S); K <- ncol(S)
  mk <- function() {
    lp <- lpSolveAPI::make.lp(N, K)
    for (k in seq_len(K)) { nz <- which(S[, k] != 0); if (length(nz)) lpSolveAPI::set.column(lp, k, S[nz, k], nz) }
    lpSolveAPI::set.constr.type(lp, rep(">=", N))
    lpSolveAPI::set.rhs(lp, numeric(N))
    lpSolveAPI::set.bounds(lp, lower = rep(-1, K), upper = rep(1, K))
    lpSolveAPI::lp.control(lp, sense = "max")
    lp
  }
  lp <- mk()
  sep <- logical(N); v <- numeric(K); it <- 0L
  repeat {
    obj <- colSums(S[!sep, , drop = FALSE])
    lpSolveAPI::set.objfn(lp, obj)
    status <- solve(lp); it <- it + 1L
    if (status != 0L && it > 1L) { lp <- mk(); lpSolveAPI::set.objfn(lp, obj); status <- solve(lp) }
    if (status != 0L) stop("the linear program of the separation check did not solve (lp_solve status ", status, ")", call. = FALSE)
    if (lpSolveAPI::get.objective(lp) <= tol) break
    vi <- lpSolveAPI::get.variables(lp)
    new <- !sep & drop(S %*% vi) > tol
    if (!any(new)) break
    sep <- sep | new; v <- v + vi
    if (all(sep)) break
  }
  list(sep = sep, v = v, iterations = it)
}

## the utilities (n x J) at theta of an in-scope design: the case-specific
## block of each non-base alternative, the generic columns and the offsets;
## -Inf where an alternative is unavailable
.pc_sep_utilities <- function(theta, D) {
  n <- D$n; J <- D$J; p <- D$p; q <- D$q; nb <- setdiff(seq_len(J), D$base)
  V <- matrix(0, n, J)
  if (p) V[, nb] <- D$X %*% matrix(theta[seq_len(p * (J - 1L))], p, J - 1L)
  if (q) V <- V + matrix(drop(D$W %*% theta[p * (J - 1L) + seq_len(q)]), n, J, byrow = TRUE)
  V <- V + D$off
  V[D$avail != 1L] <- -Inf
  V
}

## Positive weights u on the comparisons P whose combination S'u is the
## gradient of the log-likelihood at theta. In mnl(), w_i p_il on the choice
## over l. In a ranking without ties the gradient of stage t is the sum over
## the alternatives l still unranked of p_tl (D_{a_t} - D_l), and D_{a_t} -
## D_l is the sum of the consecutive comparisons between them; with T_s the
## sum of exp(V) over the alternatives left at stage s, the comparison of
## a_r over a_{r+1} collects w_i T_{r+1} sum_{t <= r} 1/T_t, and that of the
## last ranked a_m over an unranked x collects w_i exp(V_x) sum_{t <= m} 1/T_t.
## NULL with ties (a tie group's gradient is not a positive combination of
## the comparisons).
.pc_sep_weights <- function(D, theta, P) {
  V <- .pc_sep_utilities(theta, D)
  E <- exp(V - apply(V, 1L, max))
  if (D$family == "mnl") return(D$w[P[, 1L]] * (E / rowSums(E))[P[, c(1L, 3L)]])
  A <- D$avail == 1L; R <- D$str$rank * A
  if (any(apply(R, 1L, function(r) anyDuplicated(r[r > 0L]) > 0L))) return(NULL)
  mr <- max(R)
  Tm <- vapply(seq_len(mr + 1L), function(s) rowSums(E * (A & (R >= s | R == 0L))), numeric(D$n))
  Tm <- matrix(Tm, D$n)
  C <- t(apply(1 / Tm[, seq_len(mr), drop = FALSE], 1L, cumsum)); if (mr == 1L) C <- t(C)
  i <- P[, 1L]; r <- R[P[, c(1L, 2L)]]; below <- R[P[, c(1L, 3L)]]
  D$w[i] * C[cbind(i, r)] * ifelse(below > 0L, Tm[cbind(i, r + 1L)], E[P[, c(1L, 3L)]])
}

## Overlap certified without the linear program (Stiemke's lemma: no v has
## Sv >= 0 with Sv != 0 iff some u > 0 has S'u = 0): the weights u of the
## fitted probabilities, corrected by the least-norm change that makes S'u
## vanish (at a maximum-likelihood estimate S'u is the gradient, near zero),
## twice; certified when every corrected weight stays above 1e-12 of the
## largest and the remaining S'u is at most 1e-9 of sum |S|'u. Failing that,
## the linear program decides.
## The correction is the residual of u on the columns of S by a pivoted QR (a
## least-squares projection that does not square S's condition number, as the
## normal equations would).
.pc_sep_certify <- function(S, u, rs = 1) {
  if (is.null(u) || length(u) != nrow(S) || !all(is.finite(u)) || any(u <= 0)) return(FALSE)
  qs <- tryCatch(qr(S), error = function(e) { .pc_pass_limit(e); NULL })
  if (is.null(qs)) return(FALSE)
  us <- u
  for (k in 1:2) {
    us <- tryCatch(qr.resid(qs, us), error = function(e) { .pc_pass_limit(e); NULL })
    if (is.null(us) || !all(is.finite(us))) return(FALSE)
  }
  ## positivity in the rows' own units (rs: the row scales) and in the scaled
  ## ones: the residual's round-off does not shrink with a row's scale, so the
  ## first test alone could pass a zero weight on a badly scaled row (round 7,
  ## I-N2); failing either, the linear program decides
  if (!all(is.finite(us)) || min(us / rs) <= 1e-12 * max(u / rs) || min(us) <= 1e-12 * max(u)) return(FALSE)
  max(abs(crossprod(S, us))) <= 1e-9 * max(crossprod(abs(S), us))
}

## the distinct rows of S (sorting, exact) and the index of every row's
## distinct row
.pc_sep_unique <- function(S) {
  n <- nrow(S)
  if (n < 2L) return(list(rows = seq_len(n), group = seq_len(n)))
  o <- do.call(order, unname(as.data.frame(S)))
  dup <- c(FALSE, rowSums(S[o[-1L], , drop = FALSE] != S[o[-n], , drop = FALSE]) == 0)
  g <- integer(n); g[o] <- cumsum(!dup)
  list(rows = o[!dup], group = g)
}

## The check of a design: overlap, or the separation type, the separated terms,
## the direction (unit length, theta coordinates), the cases with a separated
## comparison, and the counts. A comparison whose structure vector is zero, or
## round-off after the column scaling (two
## available alternatives with the same design row in a case) is favoured by no
## direction: it is left out of the program and the certificate but counted
## (`n_zero`), and with one the data are at most quasi-completely separated,
## since its term stays below zero along every direction (DESIGN.md Phase 4d).
## Columns are scaled by their largest absolute value and then rows by theirs
## (neither changes which comparisons a direction favours, and the program's
## tolerance becomes relative to each row); with `theta` (any finite point: the
## fit's estimates) overlap is first sought by the certificate of
## .pc_sep_certify, and otherwise the program is solved on the distinct rows.
## The separated terms are the coefficients not orthogonal to the span of the
## cone of separating directions, which is the null space of the comparisons
## no separating direction favours (singular values at most 1e-9 of the
## largest): every coefficient that moves along some separating direction,
## whichever vertex the solver returns.
.pc_separation <- function(D, theta = NULL, tol = 1e-9) {
  res <- list(separated = FALSE, type = "overlap", terms = character(0), direction = NULL, cases = D$case_id[integer(0)],
              n_comparisons = 0L, n_separated = 0L, n_zero = 0L, method = "linear program")
  P <- .pc_sep_pairs(D)
  if (!nrow(P)) return(res)
  S <- .pc_sep_struct(D, P)
  res$n_comparisons <- nrow(S)
  sc <- apply(abs(S), 2L, max); sc[sc == 0] <- 1
  S <- S / rep(sc, each = nrow(S))
  ## a zero comparison: its structure vector vanishes (.pc_sep_struct() sets a
  ## difference within rounding of its operands, or of the column's typical
  ## size, to zero, so two design rows equal but for their last digits give
  ## one); it is set aside before the rows are scaled. The rule is not
  ## relative to a column's largest entry, so a small genuine difference in a
  ## column with an outlier stays (round 7, D-1; round 6 used 1e-10 of the
  ## column's largest entry)
  rs <- if (ncol(S)) apply(abs(S), 1L, max) else rep(0, nrow(S))      # no parameters: every comparison is zero
  nz <- rs > 0
  res$n_zero <- sum(!nz)
  if (!all(nz)) { S <- S[nz, , drop = FALSE]; P <- P[nz, , drop = FALSE]; rs <- rs[nz] }
  if (!nrow(S)) return(res)
  S <- S / rs
  if (!is.null(theta) && length(theta) == D$K && all(is.finite(theta))) {
    w <- tryCatch(.pc_sep_weights(D, theta, P), error = function(e) { .pc_pass_limit(e); NULL })
    if (.pc_sep_certify(S, if (!is.null(w)) w * rs, rs)) { res$method <- "certificate"; return(res) }   # the weights of the rescaled rows
  }
  u <- .pc_sep_unique(S)
  lp <- .pc_sep_lp(S[u$rows, , drop = FALSE], tol)
  if (!any(lp$sep)) return(res)
  sep <- lp$sep[u$group]
  E <- S[u$rows[!lp$sep], , drop = FALSE]
  K <- ncol(S)
  Vn <- if (!nrow(E)) diag(K) else {
    s <- svd(E, nu = 0L, nv = K)
    rk <- sum(s$d > 1e-9 * s$d[1L])
    s$v[, seq_len(K) > rk, drop = FALSE]
  }
  vt <- lp$v / sc; vt <- vt / sqrt(sum(vt^2))
  res$separated <- TRUE
  res$type <- if (all(sep) && res$n_zero == 0L) "complete" else "quasi-complete"
  res$terms <- D$par_names[rowSums(Vn^2) > 1e-10]
  res$direction <- stats::setNames(vt, D$par_names)
  res$cases <- D$case_id[sort(unique(P[sep, 1L]))]
  res$n_separated <- sum(sep)
  res$sep_pairs <- P[sep, , drop = FALSE]                    # (case, above, below) of the separated comparisons
  res
}

## The alternatives whose probability vanishes at the supremum of a
## maximum-likelihood fit in the separation check's scope: never chosen (in a
## ranking, never ranked above another available alternative) and below the
## chosen one in separated comparisons only, so that along the separating
## direction their utility falls behind in every case. Their coefficients do
## not enter the supremum, which is the likelihood of the model without them,
## so a likelihood ratio of suprema does not count them in its degrees of
## freedom (local review 2026-11-06, finding 10). An alternative never chosen
## whose estimate exists (no constant, a covariate of both signs) is not
## absent. Returns their indices (possibly empty).
.pc_absent_alts <- function(D) {
  if (!.pc_sep_scope(D) || .pc_penalized(D)) return(integer(0))
  live <- D$w > 0
  cnt <- if (identical(D$family, "rologit")) .pc_rank_wins(D$str$rank * live, D$avail) else tabulate(D$y[live & rowSums(D$avail) >= 2L], nbins = D$J)
  cand <- which(cnt == 0L & colSums(D$avail[live, , drop = FALSE] == 1L) > 0L)
  if (!length(cand)) return(integer(0))
  sp <- tryCatch(.pc_separation(D), error = function(e) { .pc_pass_limit(e); NULL })
  if (is.null(sp) || !isTRUE(sp$separated) || is.null(sp$sep_pairs)) return(integer(0))
  P <- .pc_sep_pairs(D)
  cand[vapply(cand, function(j) { nb <- sum(P[, 3L] == j); nb > 0L && sum(sp$sep_pairs[, 3L] == j) == nb }, logical(1))]
}

## the number of parameters of a design that enter its supremum: those of an
## absent alternative's case-specific block do not (.pc_absent_alts)
.pc_k_effective <- function(D, absent = .pc_absent_alts(D)) D$K - length(absent) * D$p

## the diverging coefficients of a design fitted at theta (logical, in
## D$par_names order), as the IIA tests judge their fits and half-sample
## refits: the exact check where it applies, the size-and-Wald rule
## (.pc_diverge) elsewhere or when the solver fails
.pc_diverging <- function(theta, D, scale) {
  if (.pc_sep_scope(D)) {
    s <- tryCatch(.pc_separation(D, theta), error = function(e) { .pc_pass_limit(e); NULL })
    if (!is.null(s)) return(D$par_names %in% s$terms)
  }
  .pc_diverge(theta, D, scale, se = .pc_info_se(theta, D, scale))
}

## The supremum of the maximum-likelihood log-likelihood of a design, from a
## point theta where a fit stopped with value ll (DESIGN.md, section 1, Phases
## 4c and 4d). In scope of the check:
## * overlap: ll;
## * complete separation (every comparison separated, none with a zero
##   structure vector): 0, since along the reported direction every reported
##   comparison is favoured, so every choice probability, and a tie group's
##   exact probability, tends to one;
## * quasi-complete separation: the value reached by moving along the
##   direction until every separated comparison's utility difference is at
##   least 40 and re-maximizing there, repeated until the gain is below 1e-10
##   (at most five rounds). Along the direction the separated comparisons'
##   terms vanish and the others do not change. In a multinomial design the
##   value is then within the number of separated comparisons times e^-40 (and
##   the optimizer's tolerance) of the supremum; in a ranking the move bounds
##   the consecutive comparisons only, and the re-maximization and the rounds
##   take up the rest.
## Out of scope (or when the check fails): ll. Returns the value, the
## separation type ("overlap", "complete", "quasi-complete"; NA out of scope),
## and `theta`, the point reached: the fit's point under overlap or out of
## scope, the moved and re-maximized point under quasi-complete separation,
## the moved point under complete separation (where the separated
## comparisons are near-deterministic; the parametric bootstrap draws from it).
## The penalty of a penalized design is ignored: this is the likelihood's
## supremum.
.pc_sup_loglik <- function(D, theta, ll) {
  D0 <- D; D0$firth <- NULL; D0$prior <- NULL
  out <- list(value = ll, type = NA_character_, theta = theta)
  if (!.pc_sep_scope(D0) || length(theta) != D0$K || !all(is.finite(theta))) return(out)
  ## a check that fails leaves the value where the optimizer stopped, and says
  ## so: paired with another fit's supremum it would bias the ratio (local
  ## review 2026-11-06, finding 11)
  sep <- tryCatch(.pc_separation(D0, theta), error = function(e) {
    .pc_pass_limit(e)
    warning("The separation check failed (", conditionMessage(e), "); a likelihood-ratio statistic uses this fit's log-likelihood where its optimizer stopped, ",
            "which is below the supremum if its data are separated.", call. = FALSE)
    NULL
  })
  if (is.null(sep)) return(out)
  if (!sep$separated) { out$type <- "overlap"; return(out) }
  P <- .pc_sep_pairs(D0); S <- .pc_sep_struct(D0, P)
  v <- unname(sep$direction)
  m <- drop(S %*% v)
  keep <- m > 1e-9 * max(abs(m))                          # the separated comparisons (the direction favours them)
  S <- S[keep, , drop = FALSE]; m <- m[keep]
  od <- D0$off[P[keep, c(1L, 2L), drop = FALSE]] - D0$off[P[keep, c(1L, 3L), drop = FALSE]]
  push <- function(th) th + max(0, (40 - (drop(S %*% th) + od)) / m) * v
  if (sep$type == "complete") return(list(value = 0, type = "complete", theta = push(theta)))
  fn <- .pc_fns(D0); scale <- .pc_scales(D0)
  best <- ll; th <- theta
  for (round in 1:5) {
    th1 <- push(th)
    op <- tryCatch(suppressWarnings(.pc_optimize(fn, th1, scale, 1000L, 1e-12, newton_first = TRUE)), error = function(e) { .pc_pass_limit(e); NULL })
    val <- if (!is.null(op) && is.finite(op$value)) -op$value else tryCatch(-fn$nll(th1), error = function(e) { .pc_pass_limit(e); NA_real_ })
    if (!is.finite(val) || val - best < 1e-10) {
      ## no further gain: the moved point is kept when it is as good (the
      ## separated comparisons are then near-deterministic there)
      if (is.finite(val) && val >= best - 1e-8 * max(1, abs(best))) { best <- max(best, val); th <- if (!is.null(op) && is.finite(op$value)) op$par else th1 }
      break
    }
    best <- val; th <- if (!is.null(op) && is.finite(op$value)) op$par else th1
  }
  list(value = best, type = "quasi-complete", theta = th)
}

## the suprema of two nested fits' log-likelihoods (restricted r, full f), each
## from its design, point, and value: both by .pc_sup_loglik when both designs
## are in the check's scope, else both values as they are (mixing a supremum
## with an out-of-scope fit's value could invert the comparison). `used` is
## TRUE when either supremum differs from its fit's value (separation);
## `r_separated` when the restricted design is separated, and `r_theta` the
## point its supremum is approached at (.pc_sup_loglik)
.pc_sup_pair <- function(Dr, thr, llr, Df, thf, llf) {
  if (!(.pc_sep_scope(Dr) && .pc_sep_scope(Df))) return(list(r = llr, f = llf, used = FALSE, r_separated = FALSE, r_theta = thr))
  sr <- .pc_sup_loglik(Dr, thr, llr); sf <- .pc_sup_loglik(Df, thf, llf)
  sepd <- c("complete", "quasi-complete")
  list(r = sr$value, f = sf$value, used = any(c(sr$type, sf$type) %in% sepd), r_separated = isTRUE(sr$type %in% sepd), r_theta = sr$theta)
}

## the method-line note of a statistic computed from suprema
.pc_sup_note <- "; log-likelihood suprema under separation (where the optimizer stopped is not the supremum)"

## a refit whose maximum enters a likelihood-ratio statistic: its warning that
## it did not converge is kept unless its design is separated (there the
## optimizer cannot meet its tolerance along the separated direction, and the
## statistic takes the supremum, .pc_sup_loglik)
.pc_refit_quiet <- function(D, ...) {
  w <- character(0)
  r <- withCallingHandlers(.pc_refit(D, ...), warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  if (length(w)) {
    sep <- .pc_sep_scope(D) && isTRUE(tryCatch(.pc_separation(D, r$coefficients)$separated, error = function(e) { .pc_pass_limit(e); FALSE }))
    if (!sep) for (m in w) warning(m, call. = FALSE)
  }
  r
}

## the check inside a fit: NULL when out of scope; on a solver failure a
## warning and NULL (the fit then uses the size-and-Wald rule)
.pc_separation_fit <- function(D, theta = NULL) {
  if (!.pc_sep_scope(D)) return(NULL)
  tryCatch(.pc_separation(D, theta), error = function(e) {
    .pc_pass_limit(e)
    warning("The separation check failed (", conditionMessage(e), "); diverging coefficients are judged by their size and Wald statistics instead.", call. = FALSE)
    NULL
  })
}

#' Separation in multinomial, conditional and rank-ordered logit data
#'
#' Checks whether the data of an [mnl()] or [rologit()] fit are separated, so
#' that the maximum-likelihood estimates do not exist. Each observation
#' reports that some alternatives are preferred to others: the choice to every
#' other available alternative, or, in a ranking, each alternative of a rank
#' group to each of the next (ties as groups, the unranked alternatives last).
#' A finite maximum-likelihood estimate exists if and only if no direction of
#' the coefficients weakly favours every one of these comparisons and strictly
#' favours at least one (overlap; Albert and Anderson 1984 for logistic
#' regression, Speckman, Lee and Sun 2009 for choice models, Sablica, Hornik
#' and Rusch 2026 for rankings). Along such a direction the log-likelihood
#' rises towards its supremum and the estimates run off. Overlap is first
#' sought by a certificate: positive weights on the comparisons that balance
#' them exactly (Stiemke's lemma), built from the fitted probabilities at the
#' fit's estimates, where the gradient of the log-likelihood is such a
#' weighted sum (not for tied rankings). Otherwise the check solves a linear
#' program (lp_solve, through \pkg{lpSolveAPI}). Both are exact up to
#' floating-point tolerances.
#'
#' The maximum-likelihood fits of [mnl()] and [rologit()] (without unit fixed
#' effects or inflation) run this check and name the diverging coefficients in
#' a warning; penalized fits (`firth = TRUE` or a prior) store them in
#' `separation`, and `summary()` marks them `[sep]`. The other estimators
#' judge divergence by the size of a coefficient and its Wald statistic.
#'
#' A comparison between two alternatives with the same design row in a case
#' (possible without alternative-specific constants), or with rows that
#' differ only by rounding (each attribute's difference within 64 machine
#' epsilons of the larger of its two values and the attribute's typical
#' size, the median of its nonzero absolute values; such as attributes
#' computed as 0.1 + 0.2 and 0.3, centred or not),
#' is favoured by no direction. It is counted (`n_zero`), and with one the data are at most
#' quasi-completely separated, since that choice's probability stays below one
#' along every direction.
#'
#' @param object A fit from [mnl()] or [rologit()], without unit fixed effects
#'   or inflation (maximum likelihood or penalized: separation is a property
#'   of the data and the model, not of the fit).
#' @return An object of class `polychoice_separation`: `separated` (logical);
#'   `type` (`"overlap"`, `"quasi-complete"`, or `"complete"`: every comparison
#'   separated); `terms`, the coefficients that move along some separating
#'   direction (a term can have a zero coordinate in the reported `direction`
#'   and still move along another); `direction`, a separating direction (unit
#'   length; along it every comparison some separating direction can favour is
#'   favoured); `cases`, the identifiers of the cases with a separated
#'   comparison; `n_comparisons` (all the comparisons the data report),
#'   `n_separated`, `n_zero` (the comparisons with identical design rows, up
#'   to rounding); and
#'   `method`, `"certificate"` or `"linear program"`.
#' @references
#' Albert, A. and Anderson, J.A. (1984). On the existence of maximum
#' likelihood estimates in logistic regression models. *Biometrika*, 71, 1-10.
#' \doi{10.1093/biomet/71.1.1}
#'
#' Sablica, L., Hornik, K. and Rusch, T. (2026). The geometry of existence and
#' uniqueness of maximum likelihood estimation in categorical response models.
#' arXiv:2610.08100.
#'
#' Speckman, P.L., Lee, J. and Sun, D. (2009). Existence of the MLE and
#' propriety of posteriors for a general multinomial choice model.
#' *Statistica Sinica*, 19, 731-748.
#' @examples
#' data(defection)
#' ## Cook, Niehaus and Zuhlke's Model 3: no captured combatant joined for
#' ## economic reasons without being indoctrinated (an empty cell)
#' m3 <- mnl(outcome ~ economic * indoctrination + age_joined + birth_year + school_age + male,
#'           data = defection, ref = "captured", firth = TRUE)
#' separation(m3)
#' @export
separation <- function(object) {
  if (!inherits(object, "choicemodel") || is.null(object$design))
    stop("separation() takes a fit from mnl() or rologit().", call. = FALSE)
  D <- object$design
  if (!.pc_sep_scope(D))
    stop("separation() covers mnl() and rologit() fits without unit fixed effects or inflation; ",
         "the other estimators judge divergence by the size and Wald statistic of a coefficient.", call. = FALSE)
  structure(.pc_separation(D, object$theta), class = "polychoice_separation")
}

#' @export
print.polychoice_separation <- function(x, ...) {
  if (!x$separated) {
    cat("Overlap: no separating direction among the ", x$n_comparisons, " comparisons the data report;\n",
        "the maximum-likelihood estimates exist.\n", sep = "")
    if (isTRUE(x$n_zero > 0L)) cat(x$n_zero, " of them between alternatives with identical design rows.\n", sep = "")
    return(invisible(x))
  }
  cat(if (x$type == "complete") "Complete" else "Quasi-complete", " separation: ", x$n_separated, " of the ", x$n_comparisons,
      " comparisons the data report are separated, in ", length(x$cases), " case(s);\n",
      "the maximum-likelihood estimates do not exist.\n", sep = "")
  if (isTRUE(x$n_zero > 0L))
    cat(x$n_zero, " comparison(s) between alternatives with identical design rows, which no direction can separate.\n", sep = "")
  cat("Coefficients moving along some separating direction (the reported direction's coordinates):\n")
  print(round(x$direction[x$terms], 4))
  invisible(x)
}
