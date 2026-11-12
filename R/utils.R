## ---------------------------------------------------------------------------
## Internal utilities shared by the design builder, the engine, the tests, and
## the methods: argument alignment, rank checks, finite-difference Hessians,
## the parallel replicate driver used by the bootstrap.
## ---------------------------------------------------------------------------

## resolve a column-name-or-vector argument (weights / offset / cluster) to a
## vector with one entry per row of `data`
.pc_aligned <- function(x, data, what) {
  if (is.null(x)) return(NULL)
  v <- if (is.character(x) && length(x) == 1L) {
    if (!x %in% names(data)) stop("'", what, "' names a column (", x, ") that is not in 'data'.", call. = FALSE)
    data[[x]]
  } else x
  if (NROW(v) != nrow(data))
    stop("'", what, "' must be a column name in 'data' or a vector with one entry per row of 'data' (", nrow(data), ").", call. = FALSE)
  v
}

## a count-valued argument (nboot, cores, maxit, splits) is one whole number >= lower
.pc_check_whole <- function(x, name, lower = 1) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x != floor(x) || x < lower)
    stop("'", name, "' must be a whole number of at least ", lower, ".", call. = FALSE)
  invisible(as.integer(x))
}

## whether the alternatives carry coefficients of their own (constants,
## case-specific slopes, a split equation, scales, a covariance), which are not
## identified for an alternative that is never chosen; a conditional logit with
## generic coefficients only does not, and neither does a nested logit with
## generic coefficients only, whose inclusive-value parameters belong to the
## nests (see .pc_unchosen_nests and .pc_design)
.pc_own_coef <- function(D) !(D$family %in% c("mnl", "nlogit")) || D$k > 0L || D$p > 0L

## the nests of a nested logit none of whose alternatives is chosen (`cnt`: the
## choices of each alternative): the nest's inclusive-value parameter then
## enters only the denominators, which fall with it, so it runs to zero. A
## parameter shared by every nest (equal_lambda) is lost only when no nest of
## two or more alternatives has a chosen one.
.pc_unchosen_nests <- function(str, levels, cnt) {
  multi <- lengths(str$nests) > 1L
  none <- multi & vapply(str$nests, function(a) all(cnt[match(a, levels)] == 0L), logical(1))
  if (isTRUE(str$equal) && any(multi & !none)) return(character(0))
  names(str$nests)[none]
}

## whether the choice counts `cnt` leave a parameter unidentified, or without
## a maximum-likelihood estimate: the rule of .pc_subsample_check, applied to
## the simulated responses of the bootstrap tests (a penalized design keeps
## them: its estimates are finite)
.pc_unchosen <- function(D, cnt)
  if (.pc_penalized(D)) FALSE else if (.pc_own_coef(D)) any(cnt == 0L) else identical(D$family, "nlogit") && length(.pc_unchosen_nests(D$str, D$levels, cnt)) > 0L

## A time limit set by the caller (setTimeLimit(), or a timeout built on it)
## fires as an error at the next interrupt check, wherever that falls, and R
## clears the limit as it fires. The package's own error handlers (a failed
## start, a singular solve, a failed replicate) pass it on: swallowed, it would
## leave the rest of the fit running with no limit at all. (Inside compiled
## code the limit arrives as an interrupt, which no error handler catches.)
.pc_pass_limit <- function(e) {
  lim <- gettext(c("reached elapsed time limit", "reached CPU time limit", "reached session elapsed time limit", "reached session CPU time limit"), domain = "R")
  if (conditionMessage(e) %in% lim) stop(e)
  invisible()
}

## rank check for a design block: refuse an aliased column (an empty factor
## level, an exact linear combination) by name rather than letting a flat
## likelihood direction reach the optimizer
.pc_rank_check <- function(M, what) {
  if (ncol(M) == 0L || nrow(M) == 0L) return(invisible(TRUE))
  qx <- qr(M)
  if (qx$rank < ncol(M)) {
    bad <- colnames(M)[qx$pivot[(qx$rank + 1L):ncol(M)]]
    stop("The ", what, " design matrix is rank deficient: ", paste(bad, collapse = ", "),
         " is a linear combination of the other columns (an empty factor level, an exact ",
         "linear combination); drop it or another member of the dependent set.", call. = FALSE)
  }
  invisible(TRUE)
}

## central finite-difference Hessian of a gradient function at th (symmetrized)
.pc_hess_fd <- function(grad, th, h = 1e-5) {
  K <- length(th)
  H <- matrix(0, K, K)
  for (j in seq_len(K)) {
    d <- h * max(1, abs(th[j]))
    tp <- th; tp[j] <- tp[j] + d
    tm <- th; tm[j] <- tm[j] - d
    H[, j] <- (grad(tp) - grad(tm)) / (2 * d)
  }
  ## the asymmetry of the differenced Hessian measures its numerical noise
  ## (the exact Hessian is symmetric); kept for the saddle check's tolerance
  structure((H + t(H)) / 2, asym = abs(H - t(H)))
}

## column scales for the optimizer: the (weighted) standard deviation of each
## design column, so the fit and the differenced Hessians do not depend on the
## units of measurement; constant columns (the intercept, dummies that are
## constant) keep scale 1
.pc_colscale <- function(M, w = NULL) {
  if (ncol(M) == 0L) return(numeric(0))
  wv <- if (is.null(w)) rep(1, nrow(M)) else as.numeric(w)
  sw <- sum(wv)
  s <- apply(M, 2L, function(v) {
    m <- sum(wv * v) / sw
    r <- if (sw > 0) sqrt(sum(wv * (v - m)^2) / sw) else NA_real_     # invariant to rescaling the weights
    if (!is.finite(r) || r <= 1e-10 * max(abs(v))) 1 else r          # relative: a covariate on a tiny scale is still scaled
  })
  unname(s)
}

## Parallel replicate driver for the bootstrap: `cores = 1` runs lapply();
## `cores > 1` runs a PSOCK cluster that loads this package from the calling
## session's library path and draws from L'Ecuyer streams seeded from the
## calling session's RNG, so set.seed() before the call makes a parallel run
## reproducible (the parallel and serial streams differ).
.pc_lapply <- function(X, FUN, cores = 1L) {
  cores <- as.integer(cores)[1L]
  if (is.na(cores) || cores < 1L) cores <- 1L
  if (cores == 1L || length(X) < 2L) return(lapply(X, FUN))
  cores <- min(cores, length(X))
  cl <- parallel::makeCluster(cores)
  on.exit(parallel::stopCluster(cl), add = TRUE)
  parallel::clusterCall(cl, function(lp, pkg) { .libPaths(lp); suppressPackageStartupMessages(library(pkg, character.only = TRUE)); NULL },
                        .libPaths(), utils::packageName())
  seed <- sample.int(.Machine$integer.max, 1L)
  parallel::clusterSetRNGStream(cl, iseed = seed)
  parallel::parLapply(cl, X, FUN)
}

## log-sum-exp over the rows of a matrix with -Inf for unavailable entries
.pc_rowlse <- function(V) {
  m <- apply(V, 1L, max)
  m + log(rowSums(exp(V - m)))
}

## the inflated (captive) alternatives of a fit (level indices)
.pc_kc <- function(x) if (length(x$design$kc)) x$design$kc else x$k

## the alternatives for print(): every label up to 12, else the first 8 and
## the count (unlabelled alternatives, such as potential cabinets, can number
## hundreds)
.pc_level_list <- function(lv) if (length(lv) <= 12L) paste(lv, collapse = ", ") else
  paste0(paste(lv[1:8], collapse = ", "), ", ... (", length(lv), " alternatives)")

## a model title for print / compare / texreg
.pc_title <- function(x) {
  cond <- x$layout == "long" && x$design$q > 0L
  fam <- if (is.null(x$family)) "mnl" else x$family
  if (fam == "nlogit") return(paste0("Nested logit (", length(x$design$str$nests), " nests", if (isTRUE(x$design$str$equal)) ", equal inclusive-value parameters" else "", ")"))
  if (fam == "hev") return("Heteroskedastic extreme value model")
  if (fam == "rologit") return(paste0("Rank-ordered logit", if (!is.null(x$design$str$depth)) paste0(" (first ", x$design$str$depth, " positions)") else ""))
  if (fam == "remnl") return(paste0(if (cond) "Conditional logit" else "Multinomial logit", " with random alternative-specific intercepts (", x$design$str$G, " units, ", x$design$str$nodes, "-point Gauss-Hermite per dimension)"))
  if (fam == "mnp") return(paste0("Multinomial probit (", if (x$design$J == 2L) "binary probit" else if (x$design$J == 3L) "exact bivariate normal" else paste0("GHK, ", x$design$str$draws, " draws"), ")"))
  if (x$inflated) sprintf("Inflated %s (inflated alternative%s: %s, %s split)", if (cond) "conditional logit" else "multinomial logit",
                          if (length(x$inflate) > 1L) "s" else "", paste(x$inflate, collapse = ", "), x$link)
  else if (cond) "Conditional logit" else "Multinomial logit"
}
