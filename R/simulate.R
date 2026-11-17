## ---------------------------------------------------------------------------
## Data generators for the models in the package: multinomial logit (optionally
## inflated), nested logit, multinomial probit, HEV. Wide layout by default;
## with `gamma` (generic coefficients on alternative-specific covariates) the
## data come in the long layout, which is where the structural parameters of
## the nested logit, HEV, and multinomial probit are well identified (case-specific
## covariates identify them only through the functional form; ledger).
## ---------------------------------------------------------------------------

#' Simulate multinomial-logit (optionally inflated) choice data
#'
#' Draws a data set from the model of [mnl()] / [imnl()]: standard normal
#' case-specific covariates `x1`, `x2`, ... (or those supplied in `x`),
#' alternative-specific coefficients `beta` with the first level as the base,
#' optionally alternative-specific covariates `w1`, `w2`, ... with generic
#' coefficients `gamma`, and, when `delta` is given, a split equation on `z1`,
#' `z2`, ... (or `z`) that puts cases in the captive regime, where they always choose the
#' inflated alternative.
#'
#' @param n Number of cases.
#' @param beta A \eqn{p \times (J-1)} matrix: rows are the intercept and the
#'   \eqn{p - 1} covariates, columns the non-base alternatives (column names
#'   become their labels; the base is `"a"` unless `levels` is given).
#' @param delta Optional split-equation coefficients (intercept first) on
#'   `z`; `NULL` for a plain MNL. With several inflated alternatives, a matrix
#'   with one column per inflated alternative: the coefficients of
#'   \eqn{\log(\pi_{i0} / \pi_{ic})}, the log-odds of the choice regime against
#'   captivity to that alternative (see [imnl()]).
#' @param inflate The inflated alternative (a label); default the base. Several
#'   labels draw several captive alternatives (logit split only); the output
#'   then also has the column `regime` (`"choice"` or the captive alternative).
#' @param link Link of the split equation.
#' @param gamma Optional generic coefficients on alternative-specific
#'   covariates `w1`, `w2`, ... (standard normal draws for every case and
#'   alternative). When given, the data are returned in the long layout: one
#'   row per case and alternative with columns `id`, `alt`, `chosen` (0/1), the
#'   `w`'s, the `x`'s (and, for inflated draws of [rmnl()], the `z`'s and `s`), ready for
#'   `chosen ~ w1 | x1 + x2` with `id = "id", alt = "alt"`.
#' @param x Optional data frame of case-specific covariates (`n` rows) to use
#'   instead of standard normal draws.
#' @param z Optional data frame of split-equation covariates (`n` rows) to use instead of
#'   standard normal draws.
#' @param levels Optional labels of all \eqn{J} alternatives (base first).
#' @return A data frame: in the wide layout `y` (factor) and the covariates,
#'   plus the latent regime `s` (1 = choice regime) for inflated draws; in
#'   the long layout as described under `gamma`. The true parameters are stored
#'   in the attribute `"truth"`.
#' @examples
#' set.seed(1)
#' d <- rmnl(1000, beta = cbind(b = c(0.2, -1, 0.5), c = c(-0.5, 0.8, -0.3)))
#' table(d$y)
#' di <- rmnl(1000, beta = cbind(b = c(0.2, -1, 0.5), c = c(-0.5, 0.8, -0.3)),
#'            delta = c(0.5, 1), inflate = "a")
#' table(di$y, di$s)
#' dl <- rmnl(500, beta = cbind(b = c(0.2, -1, 0.5), c = c(-0.5, 0.8, -0.3)), gamma = 0.8)
#' head(dl)
#' @seealso [mnl()], [imnl()], [simulate.choicemodel()]
#' @export
rmnl <- function(n, beta, delta = NULL, inflate = NULL, link = c("logit", "probit"), gamma = NULL, x = NULL, z = NULL, levels = NULL) {
  link <- match.arg(link)
  beta <- as.matrix(beta)
  J <- ncol(beta) + 1L
  levels <- .pc_sim_levels(beta, levels)
  sx <- .pc_sim_x(n, beta, x)
  V <- cbind(0, sx$X %*% beta)
  sw <- .pc_sim_w(n, J, gamma); V <- V + sw$V
  P <- exp(V - apply(V, 1L, max)); P <- P / rowSums(P)
  y <- apply(P, 1L, function(pp) sample.int(J, 1L, prob = pp))
  truth <- list(beta = beta, levels = levels, gamma = gamma)
  extra <- NULL
  if (!is.null(delta)) {
    r <- if (is.matrix(delta)) nrow(delta) else length(delta)
    if (is.null(z)) {
      Z <- cbind(1, matrix(stats::rnorm(n * (r - 1L)), n))
      zdf <- as.data.frame(Z[, -1L, drop = FALSE]); names(zdf) <- if (r > 1L) paste0("z", seq_len(r - 1L)) else character(0)
    } else {
      zdf <- as.data.frame(z); if (nrow(zdf) != n || ncol(zdf) != r - 1L) stop("'z' must have n rows and length(delta) - 1 columns.")
      Z <- cbind(1, as.matrix(zdf))
    }
    k <- if (is.null(inflate)) 1L else match(as.character(inflate), levels)
    if (anyNA(k) || anyDuplicated(k)) stop("'inflate' must name distinct levels.")
    if (length(k) == 1L) {
      eta <- as.vector(Z %*% delta)
      phi <- if (link == "logit") stats::plogis(eta) else stats::pnorm(eta)
      s <- stats::rbinom(n, 1L, phi)
      y[s == 0L] <- k
      extra <- cbind(zdf, s = s)
    } else {
      ## several captive alternatives: regimes by a multinomial logit, log(pi_0 / pi_c) = z' delta_c
      if (link != "logit") stop("Several inflated alternatives take the logit split.")
      if (ncol(delta) != length(k)) stop("With several inflated alternatives 'delta' is a matrix with one column per inflated alternative.")
      A <- cbind(0, -(Z %*% delta)); Pr <- exp(A - apply(A, 1L, max)); Pr <- Pr / rowSums(Pr)
      reg <- apply(Pr, 1L, function(pp) sample.int(length(pp), 1L, prob = pp))
      y[reg > 1L] <- k[reg[reg > 1L] - 1L]
      extra <- cbind(zdf, s = as.integer(reg == 1L), regime = factor(c("choice", levels[k])[reg], levels = c("choice", levels[k])))
    }
    truth$delta <- delta; truth$inflate <- levels[k]; truth$link <- link
  }
  out <- .pc_sim_out(y, levels, sx$xdf, sw, extra)
  attr(out, "truth") <- truth
  out
}

## shared helpers for the generators: labels, case-specific covariates,
## alternative-specific covariates with generic coefficients, output layout
.pc_sim_levels <- function(beta, levels) {
  J <- ncol(beta) + 1L
  if (is.null(levels)) levels <- c("a", if (!is.null(colnames(beta))) colnames(beta) else letters[1L + seq_len(J - 1L)])
  if (length(levels) != J) stop("'levels' must have ncol(beta) + 1 labels.")
  levels
}
.pc_sim_x <- function(n, beta, x) {
  p <- nrow(beta)
  if (is.null(x)) {
    X <- cbind(1, matrix(stats::rnorm(n * (p - 1L)), n))
    xdf <- as.data.frame(X[, -1L, drop = FALSE]); names(xdf) <- if (p > 1L) paste0("x", seq_len(p - 1L)) else character(0)
  } else {
    xdf <- as.data.frame(x); if (nrow(xdf) != n || ncol(xdf) != p - 1L) stop("'x' must have n rows and nrow(beta) - 1 columns.")
    X <- cbind(1, as.matrix(xdf))
  }
  list(X = X, xdf = xdf)
}
.pc_sim_w <- function(n, J, gamma) {
  if (is.null(gamma)) return(list(V = matrix(0, n, J), W = NULL))
  q <- length(gamma)
  W <- array(stats::rnorm(n * J * q), c(n, J, q))
  V <- matrix(0, n, J)
  for (k in seq_len(q)) V <- V + W[, , k] * gamma[k]
  list(V = V, W = W)
}
.pc_sim_out <- function(y, levels, xdf, sw, extra = NULL) {
  n <- length(y); J <- length(levels)
  if (is.null(sw$W)) {
    out <- data.frame(y = factor(levels[y], levels = levels), xdf)
    if (!is.null(extra)) out <- cbind(out, extra)
    return(out)
  }
  out <- data.frame(id = rep(seq_len(n), each = J), alt = factor(rep(levels, n), levels = levels),
                    chosen = as.integer(rep(y, each = J) == rep(seq_len(J), n)))
  for (k in seq_len(dim(sw$W)[3L])) out[[paste0("w", k)]] <- as.vector(t(sw$W[, , k]))     # case-major, alternatives in level order
  rows <- rep(seq_len(n), each = J)
  if (ncol(xdf)) out <- cbind(out, xdf[rows, , drop = FALSE])
  if (!is.null(extra)) out <- cbind(out, extra[rows, , drop = FALSE])
  rownames(out) <- NULL
  out
}

#' Simulate nested-logit choice data
#'
#' Draws a data set from the nested logit of [nlogit()]: utilities from `beta`
#' (and `gamma`) as in [rmnl()], nests as in [nlogit()], and inclusive-value
#' parameters `lambda` (one per nest, in the order of `nests`; a singleton
#' nest's value is irrelevant). With case-specific covariates only, `lambda`
#' is identified just through the curvature of the inclusive values and
#' recovered imprecisely; alternative-specific covariates (`gamma`) identify
#' it sharply.
#'
#' @inheritParams rmnl
#' @param nests A named list assigning every level to a nest.
#' @param lambda Inclusive-value parameters, one per nest.
#' @return A data frame as in [rmnl()]; the truth in the attribute `"truth"`.
#' @examples
#' set.seed(1)
#' d <- rnlogit(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
#'              nests = list(n1 = c("a", "b"), n2 = "c"), lambda = c(0.5, 1))
#' table(d$y)
#' @seealso [nlogit()], [rmnl()]
#' @export
rnlogit <- function(n, beta, nests, lambda, gamma = NULL, x = NULL, levels = NULL) {
  beta <- as.matrix(beta); J <- ncol(beta) + 1L
  levels <- .pc_sim_levels(beta, levels)
  if (length(lambda) != length(nests)) stop("'lambda' needs one value per nest.")
  nest <- integer(J); for (m in seq_along(nests)) nest[match(nests[[m]], levels)] <- m
  if (any(nest == 0L)) stop("'nests' must assign every level.")
  sx <- .pc_sim_x(n, beta, x)
  V <- cbind(0, sx$X %*% beta)
  sw <- .pc_sim_w(n, J, gamma); V <- V + sw$V
  P <- matrix(0, n, J)
  lam <- lambda[nest]
  I <- sapply(seq_along(nests), function(m) { mem <- which(nest == m); a <- V[, mem, drop = FALSE] / lambda[m]; mx <- apply(a, 1L, max); mx + log(rowSums(exp(a - mx))) })
  B <- sweep(I, 2L, lambda, "*"); Bm <- apply(B, 1L, max); Pn <- exp(B - Bm) / rowSums(exp(B - Bm))
  for (j in seq_len(J)) P[, j] <- exp(V[, j] / lam[j] - I[, nest[j]]) * Pn[, nest[j]]
  y <- apply(P, 1L, function(pp) sample.int(J, 1L, prob = pp))
  out <- .pc_sim_out(y, levels, sx$xdf, sw)
  attr(out, "truth") <- list(beta = beta, levels = levels, nests = nests, lambda = lambda, gamma = gamma)
  out
}

#' Simulate multinomial probit choice data
#'
#' Draws a data set from the multinomial probit of [mnprobit()]: utilities
#' from `beta` (and `gamma`) as in [rmnl()] plus multivariate normal errors
#' whose differences from the base alternative have covariance `sigma`
#' (\eqn{(J-1) \times (J-1)}, first element 1 to match the estimator's
#' normalization). With case-specific covariates only the covariance is
#' fragilely identified (Keane 1992); alternative-specific covariates
#' (`gamma`) identify it.
#'
#' @inheritParams rmnl
#' @param sigma Covariance of the error differences from the base, a
#'   \eqn{(J-1) \times (J-1)} matrix.
#' @return A data frame as in [rmnl()]; the truth in the attribute `"truth"`.
#' @examples
#' set.seed(1)
#' d <- rmnp(500, beta = cbind(b = c(0.3, -0.6, 0.5), c = c(-0.2, 0.4, -0.3)),
#'           sigma = matrix(c(1, 0.5, 0.5, 1), 2), gamma = 1)
#' head(d)
#' @seealso [mnprobit()], [rmnl()]
#' @export
rmnp <- function(n, beta, sigma, gamma = NULL, x = NULL, levels = NULL) {
  beta <- as.matrix(beta); J <- ncol(beta) + 1L; m <- J - 1L
  levels <- .pc_sim_levels(beta, levels)
  sigma <- as.matrix(sigma)
  if (!all(dim(sigma) == m)) stop("'sigma' must be a (J-1) x (J-1) matrix.")
  sx <- .pc_sim_x(n, beta, x)
  V <- cbind(0, sx$X %*% beta)
  sw <- .pc_sim_w(n, J, gamma); V <- V + sw$V
  E <- MASS::mvrnorm(n, rep(0, m), sigma)                   # error differences eta_k = eps_k - eps_base
  U <- V; U[, -1L] <- U[, -1L] + E
  y <- apply(U, 1L, which.max)
  out <- .pc_sim_out(y, levels, sx$xdf, sw)
  attr(out, "truth") <- list(beta = beta, levels = levels, sigma = sigma, gamma = gamma)
  out
}

#' Simulate heteroskedastic extreme value choice data
#'
#' Draws a data set from the HEV model of [hev()]: utilities from `beta` (and
#' `gamma`) as in [rmnl()] plus independent Gumbel errors with scales `scale`
#' (one per level, in level order; the base alternative's scale should be 1 to
#' match the estimator's normalization).
#'
#' @inheritParams rmnl
#' @param scale Gumbel scale parameters, one per level.
#' @return A data frame as in [rmnl()]; the truth in the attribute `"truth"`.
#' @examples
#' set.seed(1)
#' d <- rhev(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), scale = c(1, 0.6, 1.5))
#' table(d$y)
#' @seealso [hev()], [rmnl()]
#' @export
rhev <- function(n, beta, scale, gamma = NULL, x = NULL, levels = NULL) {
  beta <- as.matrix(beta); J <- ncol(beta) + 1L
  levels <- .pc_sim_levels(beta, levels)
  if (length(scale) != J) stop("'scale' needs one value per level.")
  sx <- .pc_sim_x(n, beta, x)
  V <- cbind(0, sx$X %*% beta)
  sw <- .pc_sim_w(n, J, gamma); V <- V + sw$V
  E <- matrix(-log(-log(stats::runif(n * J))), n, J)             # standard Gumbel draws
  U <- V + sweep(E, 2L, scale, "*")
  y <- apply(U, 1L, which.max)
  out <- .pc_sim_out(y, levels, sx$xdf, sw)
  attr(out, "truth") <- list(beta = beta, levels = levels, scale = scale, gamma = gamma)
  out
}

#' Simulate rank-ordered choice data
#'
#' Draws rankings from the rank-ordered logit of [rologit()]: utilities from
#' `beta` (and `gamma`) as in [rmnl()] plus independent standard Gumbel errors;
#' each case ranks its alternatives by utility. `depth` keeps only the first
#' positions (the rest unranked, rank 0), and `ties` merges each ranked
#' position into the one above it with that probability, independently of the
#' order within the merged group, which is the coarsening under which the
#' exact tie likelihood of [rologit()] is the likelihood of the reported
#' ranking.
#'
#' @inheritParams rmnl
#' @param gamma Optional generic coefficients on alternative-specific
#'   covariates `w1`, `w2`, ... (standard normal draws for every case and
#'   alternative). The data are in the long layout with or without them (see
#'   Value).
#' @param depth Optional number of ranked positions (default: every
#'   alternative is ranked).
#' @param ties Probability that a ranked position is tied with the position
#'   above it (default 0).
#' @return A data frame in the long layout, one row per case and alternative:
#'   `id`, `alt`, `rank` (1 = most preferred; equal values are ties; 0 =
#'   unranked), the `w`'s and the `x`'s, ready for `rank ~ w1 | x1 + x2` with
#'   `id = "id", alt = "alt"`. The true parameters are in the attribute
#'   `"truth"`.
#' @examples
#' set.seed(1)
#' d <- rrologit(300, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6), d = c(0.2, 0.4)), depth = 2)
#' head(d, 8)
#' @seealso [rologit()], [rmnl()]
#' @export
rrologit <- function(n, beta, gamma = NULL, x = NULL, levels = NULL, depth = NULL, ties = 0) {
  beta <- as.matrix(beta); J <- ncol(beta) + 1L
  levels <- .pc_sim_levels(beta, levels)
  if (!is.null(depth)) depth <- .pc_check_whole(depth, "depth", 1)
  if (!is.numeric(ties) || length(ties) != 1L || ties < 0 || ties >= 1) stop("'ties' must be a probability in [0, 1).")
  sx <- .pc_sim_x(n, beta, x)
  V <- cbind(0, sx$X %*% beta)
  sw <- .pc_sim_w(n, J, gamma); V <- V + sw$V
  U <- V + matrix(-log(-log(stats::runif(n * J))), n, J)
  ## the position of every alternative (1 = highest utility), then the group of
  ## every position after the random merges, then the cut at `depth`
  pos <- t(apply(-U, 1L, rank, ties.method = "first"))
  merge <- matrix(stats::runif(n * J) < ties, n, J); merge[, 1L] <- FALSE
  grp <- t(apply(!merge, 1L, cumsum))                          # group index of every position
  rk <- matrix(grp[cbind(rep(seq_len(n), J), as.vector(pos))], n, J)
  if (!is.null(depth)) {
    start <- t(apply(grp, 1L, function(g) match(g, g)))        # the position at which each position's group starts
    gstart <- matrix(start[cbind(rep(seq_len(n), J), as.vector(pos))], n, J)
    rk[gstart > depth] <- 0L
  }
  out <- data.frame(id = rep(seq_len(n), each = J), alt = factor(rep(levels, n), levels = levels), rank = as.vector(t(rk)))
  if (!is.null(sw$W)) for (k in seq_len(dim(sw$W)[3L])) out[[paste0("w", k)]] <- as.vector(t(sw$W[, , k]))
  rows <- rep(seq_len(n), each = J)
  if (ncol(sx$xdf)) out <- cbind(out, sx$xdf[rows, , drop = FALSE])
  rownames(out) <- NULL
  attr(out, "truth") <- list(beta = beta, levels = levels, gamma = gamma, depth = depth, ties = ties)
  out
}
