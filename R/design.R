## ---------------------------------------------------------------------------
## From formula + data to the choice design the C++ engine consumes.
##
## Two layouts, one grammar. The right-hand side has up to three parts
## separated by `|` (package Formula); the LAST part is the inflation equation
## of imnl():
##   wide (one row per case; the response is the chosen alternative):
##       y ~ x                 mnl()
##       y ~ x | z             imnl()
##   long (one row per case-alternative pair; `id` and `alt` name the case and
##   alternative columns; the response is the 0/1 choice indicator):
##       choice ~ w            mnl(): alternative-specific covariates only
##       choice ~ w | x        mnl()
##       choice ~ w | x | z    imnl()  (use `| 1 |` when there are no x)
## x: case-specific covariates, alternative-specific coefficients beta_j (the
## intercept of x, or of w in the long layout, gives the alternative-specific
## constants). w: alternative-specific covariates, generic coefficient gamma.
## z: split-equation covariates, delta.
##
## Unit fixed effects (`fe`) are alternative-specific constants per unit,
## implemented as alternative-specific dummy columns in W (.pc_fe_columns): a
## unit that never chooses alternative j has that constant at -Inf, which is
## the likelihood with j unavailable to the unit, so those rows get avail = 0
## and no dummy (imnl() refuses fixed effects). The common
## constant of alternative j refers to a reference unit: the first unit (in
## level order) that chose both j and the base, or else one that links j to an
## alternative already anchored (breadth-first); a unit that never chose the
## base omits one dummy. The rules and the drop of cases left with a single
## alternative are repeated until nothing changes.
##
## The design D returned by .pc_design() holds: X (n x p), W ((n*J) x q, row
## (i-1)*J + j), Z (n x r), off (n x J), avail (n x J, 0/1), y (1..J), w
## (weights), cluster, the terms/xlevels/contrasts of every part for
## predict(), the estimation data, and the parameter names and index blocks.
## ---------------------------------------------------------------------------

## split the Formula into part indices for the layout / inflation combination
.pc_parts <- function(formula, inflated, layout) {
  F <- Formula::Formula(formula)
  if (length(F)[1L] != 1L) stop("The formula must have exactly one response.", call. = FALSE)
  np <- length(F)[2L]
  if (layout == "wide") {
    if (inflated) {
      if (np != 2L) stop("imnl() with one row per case takes `y ~ x | z` (outcome covariates | inflation covariates).", call. = FALSE)
      parts <- list(x = 1L, w = NULL, z = 2L)
    } else {
      if (np != 1L) stop("mnl() with one row per case takes `y ~ x`; a `|` part is the split equation of imnl(), ",
                         "or the case-specific part of a long layout (give `id` and `alt`).", call. = FALSE)
      parts <- list(x = 1L, w = NULL, z = NULL)
    }
  } else {
    if (inflated) {
      if (np != 3L) stop("imnl() with a long layout takes `choice ~ w | x | z` (alternative-specific | case-specific | ",
                         "inflation covariates); write `| 1 |` for an empty case-specific part.", call. = FALSE)
      parts <- list(w = 1L, x = 2L, z = 3L)
    } else {
      if (!np %in% 1:2) stop("mnl() with a long layout takes `choice ~ w` or `choice ~ w | x`.", call. = FALSE)
      parts <- list(w = 1L, x = if (np == 2L) 2L else NULL, z = NULL)
    }
  }
  parts$F <- F
  parts
}

## model matrix of one part on data frame `d`; at fit time (tm NULL) the terms
## with predvars, the factor levels, and the contrasts are recorded; at
## prediction time they are reused so that poly()/ns()/factor columns match
.pc_part_matrix <- function(F, rhs, d, tm = NULL) {
  if (is.null(tm)) {
    tt <- stats::terms(F, rhs = rhs)
    mf <- stats::model.frame(tt, d, na.action = stats::na.pass)
    tt <- attr(mf, "terms")
    M <- stats::model.matrix(tt, mf)
    fac <- vapply(mf, is.factor, logical(1))
    if (attr(tt, "response") > 0L) fac[attr(tt, "response")] <- FALSE   # the response is not a predictor level set
    xlev <- if (any(fac)) lapply(mf[fac], levels) else NULL
    list(M = M, tm = list(terms = tt, xlevels = xlev, contrasts = attr(M, "contrasts")))
  } else {
    tt <- stats::delete.response(tm$terms)
    mf <- stats::model.frame(tt, d, na.action = stats::na.pass, xlev = tm$xlevels)
    M <- stats::model.matrix(tt, mf, contrasts.arg = tm$contrasts)
    list(M = M, tm = tm)
  }
}

## the per-alternative unit dummies for `fe` (see the file header): returns the
## columns to append to W, a map (column name, alternative, unit), the
## availability mask with never-chosen unit-alternative pairs switched off, and
## `off`, the cells this rule switched off (as opposed to alternatives absent
## from a case's choice set in the data; the inflated kernel treats the two
## differently, see .pc_kernel_avail()).
##
## Normalization. A unit's constants are identified only through the
## alternatives it chooses: a unit that chose the base pins its constants
## against the base; a unit that never chose the base (the base is then
## unavailable to it, the limiting likelihood) pins only their differences, so
## one of its dummies is omitted. The common constant of alternative j refers
## to a reference unit: the first unit that chose both j and the base, or else
## (breadth first) a unit that chose j and an alternative already anchored; a
## reference unit that never chose the base omits the dummy of the alternative
## it anchors through. An alternative not connected to the base through the
## units' choices has no identified constant, which is refused by name.
.pc_fe_columns <- function(unit, y, J, base, levels, n, avail, fe, asc = TRUE, build = TRUE) {
  G <- levels(unit)
  ui <- as.integer(unit)
  ## a case with a single available alternative in the data (a captive choice)
  ## says nothing about the unit's preferences: it does not count as a choice
  informative <- rowSums(avail) >= 2L
  chosen <- matrix(FALSE, length(G), J)
  for (j in seq_len(J)) chosen[, j] <- tabulate(ui[y == j & informative], nbins = length(G)) > 0L
  avail_in <- avail
  no_base <- which(!chosen[, base])
  if (length(no_base)) avail[ui %in% no_base, base] <- 0L
  refg <- integer(J); own_ref <- rep(NA_integer_, length(G))
  multi <- rowSums(chosen) >= 2L                 # a unit that chose one alternative links nothing
  if (asc) {
  for (j in seq_len(J)) if (j != base) {
    cand <- which(chosen[, j] & chosen[, base] & multi); if (length(cand)) refg[j] <- cand[1L]
  }
  anchored <- c(base, which(refg > 0L))
  repeat {
    progress <- FALSE
    for (g in which(multi)) {
      S <- which(chosen[g, ]); new <- setdiff(S, anchored); old <- intersect(S, anchored)
      if (!length(new) || !length(old)) next
      l <- if (base %in% old) base else if (!is.na(own_ref[g]) && own_ref[g] %in% old) own_ref[g] else old[1L]
      refg[new] <- g; anchored <- c(anchored, new); progress <- TRUE
      if (l != base) own_ref[g] <- l
    }
    if (!progress) break
  }
  loose <- setdiff(which(colSums(chosen[multi, , drop = FALSE]) > 0L), anchored)
  if (length(loose))
    stop("With unit fixed effects the constant(s) of ", paste(levels[loose], collapse = ", "),
         " are not identified relative to the base '", levels[base], "': no unit that chose ",
         if (length(loose) > 1L) "them" else "it", " also chose the base or an alternative linked to it. ",
         "Drop the fixed effects or the cases of those alternatives.", call. = FALSE)
  }
  ## without alternative-specific constants there is nothing to anchor: every
  ## unit keeps all its dummies (a unit that never chose the base drops one,
  ## its own normalization)
  for (g in no_base) if (is.na(own_ref[g])) { S <- setdiff(which(chosen[g, ]), base); if (length(S)) own_ref[g] <- S[1L] }
  cols <- list(); map <- list()
  dropped <- 0L
  for (j in seq_len(J)) {
    if (j == base) next
    never <- which(!chosen[, j])
    if (length(never)) { avail[ui %in% never, j] <- 0L; dropped <- dropped + length(never) }
    if (build) for (g in which(chosen[, j])) {          # build = FALSE: the rules only (the fixed-point loops)
      if (g == refg[j]) next
      if (!is.na(own_ref[g]) && own_ref[g] == j) next
      col <- numeric(n * J)
      col[(which(ui == g) - 1L) * J + j] <- 1
      cols[[length(cols) + 1L]] <- col
      map[[length(map) + 1L]] <- data.frame(name = paste0(levels[j], ":", fe, G[g]), alt = j, unit = G[g], stringsAsFactors = FALSE)
    }
  }
  Wfe <- if (length(cols)) do.call(cbind, cols) else matrix(0, n * J, 0L)
  map <- if (length(map)) do.call(rbind, map) else data.frame(name = character(0), alt = integer(0), unit = character(0))
  colnames(Wfe) <- map$name
  list(W = Wfe, map = map, avail = avail, off = avail_in == 1L & avail == 0L, dropped_pairs = dropped,
       no_base = length(no_base), no_base_units = G[no_base])
}

## rebuild the fixed-effect dummy columns from the map for a vector of unit labels
.pc_fe_rebuild <- function(map, u, J, n) {
  Wfe <- matrix(0, n * J, nrow(map), dimnames = list(NULL, map$name))
  for (cc in seq_len(nrow(map))) Wfe[(which(u == map$unit[cc]) - 1L) * J + map$alt[cc], cc] <- 1
  Wfe
}

## Exact identification check: the parameters enter only through the within-case
## differences of the stacked design rows D_ij = (e_j (x) x_i, w_ij) over the
## available alternatives, so the Gram matrix of the within-case demeaned rows
## must have full rank. Catches an alternative-specific column that is constant
## within every case (belongs in the case-specific part), alternative dummies
## duplicating the constants, and any other aliasing, and names the columns.
.pc_identification_check <- function(D) {
  n <- D$n; J <- D$J; p <- D$p; q <- D$q; K0 <- p * (J - 1L) + q
  if (K0 == 0L) return(invisible(TRUE))
  A <- D$avail; Ji <- rowSums(A)
  nonbase <- setdiff(seq_len(J), D$base)
  Wbar <- if (q) { Wb <- matrix(0, n, q); for (j in seq_len(J)) Wb <- Wb + A[, j] * D$W[(seq_len(n) - 1L) * J + j, , drop = FALSE]; Wb / Ji } else NULL
  G <- matrix(0, K0, K0)
  for (j in seq_len(J)) {
    Dj <- matrix(0, n, K0)
    if (p) for (m in seq_along(nonbase)) {
      coefm <- (j == nonbase[m]) - A[, nonbase[m]] / Ji        # e_j minus the within-case mean of e_m
      Dj[, (m - 1L) * p + seq_len(p)] <- D$X * coefm
    }
    if (q) Dj[, p * (J - 1L) + seq_len(q)] <- D$W[(seq_len(n) - 1L) * J + j, , drop = FALSE] - Wbar
    Dj <- Dj * A[, j]
    G <- G + crossprod(Dj)
  }
  nm <- D$par_names[seq_len(K0)]
  dg <- diag(G)
  ## a column with no within-case variation first: its within-case sum of
  ## squares is zero up to round-off relative to its own raw sum of squares
  ## (relative to the largest column it would refuse any covariate on a scale
  ## 1e6 times another's)
  raw <- numeric(K0)
  if (p) for (m in seq_along(nonbase)) raw[(m - 1L) * p + seq_len(p)] <- colSums(D$X^2 * A[, nonbase[m]])
  if (q) {
    ## the alternative-specific columns against their own centred sum of squares
    ## (a shifted column, w + 1e9, varies within cases just as w does)
    Wm <- Reduce(`+`, lapply(seq_len(J), function(j) colSums(D$W[(seq_len(n) - 1L) * J + j, , drop = FALSE] * A[, j]))) / sum(A)
    for (j in seq_len(J)) raw[p * (J - 1L) + seq_len(q)] <- raw[p * (J - 1L) + seq_len(q)] +
      colSums(sweep(D$W[(seq_len(n) - 1L) * J + j, , drop = FALSE], 2L, Wm)^2 * A[, j])
  }
  zero <- dg <= 1e-12 * raw
  if (any(zero))
    stop("The choice design is not identified: column(s) ", paste(nm[zero], collapse = ", "),
         " do not vary within cases across the available alternatives (a case-specific variable belongs in the ",
         "case-specific part; a variable constant within units is absorbed by unit fixed effects).", call. = FALSE)
  sc <- sqrt(dg)
  qx <- qr(G / tcrossprod(sc), tol = 1e-9)
  if (qx$rank < K0) {
    bad <- nm[qx$pivot[(qx$rank + 1L):K0]]
    stop("The choice design is not identified: column(s) ", paste(bad, collapse = ", "),
         " are linear combinations of the other columns (for example alternative dummies next to the alternative-specific constants).", call. = FALSE)
  }
  invisible(TRUE)
}

## The rank matrix of rologit() from the long-layout rank column `rk` (rows
## indexed by case ci and alternative ji): per case the dense rank of every
## ranked alternative (1 = the most preferred group; equal values are ties), 0
## for the unranked and the unavailable. With `depth`, a group that starts
## below that position becomes unranked (a group straddling the cut stays
## whole). A tie group with alternatives ranked below it is summed over its
## orders, at most 7! of them; a larger one is refused. `top` is the first
## choice of every case: the first alternative, in level order, of its top group.
.pc_rank_matrix <- function(rk, ci, ji, n, J, depth = NULL) {
  pos <- rk > 0
  dense <- integer(length(rk))
  if (any(pos)) dense[pos] <- as.integer(stats::ave(rk[pos], ci[pos], FUN = function(v) match(v, sort(unique(v)))))
  if (!is.null(depth)) {
    start <- stats::ave(dense, ci, FUN = function(v) vapply(v, function(r) sum(v > 0L & v < r), numeric(1)) + 1)
    dense[pos & start > depth] <- 0L
  }
  big <- tapply(dense, ci, function(v) { tb <- tabulate(v[v > 0L]); any(vapply(seq_along(tb), function(r) tb[r] > 7L && (any(v > r) || any(v == 0L)), logical(1))) })
  if (any(big)) stop(sum(big), " case(s) tie more than 7 alternatives at a rank above other alternatives; their likelihood sums over every order of the tied ",
                     "alternatives (more than 5040). Break the ties or use 'depth' to rank only the first positions.", call. = FALSE)
  R <- matrix(0L, n, J); R[cbind(ci, ji)] <- dense
  list(R = R, top = max.col(R == 1L, ties.method = "first"))
}

## ---------------------------------------------------------------------------
## Model families. Every family shares the design (X, W, Z, off, avail), the
## beta / gamma (/ delta) blocks, and the engine; a family adds a block of
## structural parameters at the end of theta and its own C++ kernel:
##   mnl     none (the inflated model's delta block sits before the structural block)
##   nlogit  log inclusive-value parameters, one per non-singleton nest (or one
##           shared when equal = TRUE); natural scale lambda = exp(ell)
##   hev     log Gumbel scale parameters of the non-base alternatives (the base
##           scale is 1); natural scale theta = exp(t)
##   mnp     the free Cholesky elements of the covariance of the utility
##           differences from the base ((J-1) x (J-1), element (1,1) fixed at
##           1): natural scale = the covariance elements themselves
## .pc_structural() validates the family's specification and returns the
## structural block (names, internal-to-natural transform, family-specific
## maps); .pc_natural() maps the full internal vector to the reported one.
## ---------------------------------------------------------------------------
.pc_structural <- function(family, str, levels, base) {
  J <- length(levels)
  if (family == "mnl") return(list(names = character(0), block = character(0), n = 0L))
  if (family == "rologit") return(list(names = character(0), block = character(0), n = 0L, depth = str$depth))
  if (family == "nlogit") {
    nests <- str$nests
    if (!is.list(nests) || is.null(names(nests)) || any(names(nests) == "")) stop("'nests' must be a named list of alternative labels, e.g. list(fast = c(\"air\", \"train\"), slow = c(\"bus\", \"car\")).", call. = FALSE)
    all_alt <- unlist(nests, use.names = FALSE)
    if (anyDuplicated(all_alt)) stop("'nests': alternative(s) ", paste(unique(all_alt[duplicated(all_alt)]), collapse = ", "), " appear in more than one nest.", call. = FALSE)
    missing <- setdiff(levels, all_alt); extra <- setdiff(all_alt, levels)
    if (length(extra)) stop("'nests' names alternative(s) not in the response: ", paste(extra, collapse = ", "), ".", call. = FALSE)
    if (length(missing)) stop("'nests' must assign every alternative; missing: ", paste(missing, collapse = ", "), ".", call. = FALSE)
    if (length(nests) < 2L) stop("'nests' must define at least two nests (one nest is the multinomial logit).", call. = FALSE)
    nest <- integer(J)
    for (m in seq_along(nests)) nest[match(nests[[m]], levels)] <- m - 1L
    size <- vapply(nests, length, integer(1))
    if (all(size == 1L)) stop("Every nest is a singleton; the model is the multinomial logit.", call. = FALSE)
    equal <- isTRUE(str$equal)
    lmap <- rep(-1L, length(nests))
    if (equal) { lmap[size > 1L] <- 0L; names <- "lambda" }
    else { lmap[size > 1L] <- seq_len(sum(size > 1L)) - 1L; names <- paste0("lambda:", names(nests)[size > 1L]) }
    if (any(size == 1L)) message("Singleton nest(s) ", paste(names(nests)[size == 1L], collapse = ", "), ": their inclusive-value parameter does not enter the likelihood and is fixed at 1.")
    list(names = names, block = rep("nest", length(names)), n = length(names), nest = nest, lmap = lmap, nests = nests, equal = equal,
         transform = function(v) exp(v))
  } else if (family == "hev") {
    nb <- setdiff(seq_len(J), base)
    list(names = paste0("scale:", levels[nb]), block = rep("scale", J - 1L), n = J - 1L, transform = function(v) exp(v))
  } else if (family == "mnp") {
    ## the GHK simulator's Halton sequences use the first 20 primes, one per
    ## dimension of the differenced utilities
    if (J > 21L) stop("mnprobit() handles at most 21 alternatives (the GHK simulator's 20 Halton dimensions); this model has ", J, ".", call. = FALSE)
    m <- J - 1L
    nb <- setdiff(seq_len(J), base)
    ## free Cholesky elements of the m x m differenced covariance, element (1,1) = 1:
    ## row-major lower triangle without (1,1); diagonal elements on the log scale
    idx <- which(lower.tri(matrix(0, m, m), diag = TRUE), arr.ind = TRUE)
    idx <- idx[order(idx[, 1L], idx[, 2L]), , drop = FALSE]
    idx <- idx[!(idx[, 1L] == 1L & idx[, 2L] == 1L), , drop = FALSE]
    diag_el <- idx[, 1L] == idx[, 2L]
    cov_names <- apply(idx, 1L, function(r) paste0("cov:", levels[nb][r[1L]], ",", levels[nb][r[2L]]))
    cov_names[diag_el] <- paste0("var:", levels[nb][idx[diag_el, 1L]])
    draws <- if (is.null(str$draws)) 200L else .pc_check_whole(str$draws, "draws", 20)
    ## starting values: the covariance of the differences under iid errors,
    ## (I + 11')/2 (unit first variance), and independent differences (L = I)
    L0 <- t(chol((diag(m) + matrix(1, m, m)) / 2))
    start <- ifelse(diag_el, log(L0[idx]), L0[idx])
    list(names = cov_names, block = rep("covariance", nrow(idx)), n = nrow(idx), chol_idx = idx, diag_el = diag_el, m = m,
         draws = draws, start = start, start_alt = numeric(nrow(idx)),
         transform = function(v) {
           L <- matrix(0, m, m); L[1L, 1L] <- 1
           L[idx] <- ifelse(diag_el, exp(v), v)
           S <- L %*% t(L)
           S[idx]
         })
  } else if (family == "remnl") {
    nb <- setdiff(seq_len(J), base)
    if (is.null(str$unit) || !is.character(str$unit) || length(str$unit) != 1L) stop("'unit' must name a column in 'data'.", call. = FALSE)
    list(names = paste0("sd:", levels[nb]), block = rep("random intercept", J - 1L), n = J - 1L, unit = str$unit,
         nodes = .pc_check_whole(if (is.null(str$nodes)) 12L else str$nodes, "nodes", 2), transform = function(v) exp(v))
  } else stop("internal: unknown family ", family)
}

## probabilists' Gauss-Hermite rule (nodes and weights for expectations under
## N(0, 1)) by the Golub-Welsch eigenvalue method
.pc_gauss_hermite <- function(K) {
  if (K == 1L) return(list(nodes = 0, weights = 1))
  i <- seq_len(K - 1L); b <- sqrt(i / 2)
  Jm <- matrix(0, K, K); Jm[cbind(i, i + 1L)] <- b; Jm[cbind(i + 1L, i)] <- b
  e <- eigen(Jm, symmetric = TRUE)
  w <- e$vectors[1L, ]^2
  list(nodes = e$values * sqrt(2), weights = w / sum(w))
}

## unit indexing for the random-intercept family: cases grouped by unit, unit
## weights (frequency weights replicate units, so they must be constant within
## one), unit clusters, rebuilt whenever the cases of a design change
.pc_re_index <- function(D) {
  D$unit <- droplevels(D$unit)
  u <- as.integer(D$unit); G <- nlevels(D$unit)
  D$str$uidx <- order(u) - 1L
  D$str$ustart <- c(0L, cumsum(tabulate(u, nbins = G)))
  first <- match(seq_len(G), u)
  wu <- D$w[first]
  if (any(abs(D$w - wu[u]) > 1e-12)) stop("'weights' must be constant within a unit when units share random intercepts (they replicate units).", call. = FALSE)
  D$str$wu <- wu
  if (!is.null(D$cluster)) {
    cu <- D$cluster[first]
    if (any(D$cluster != cu[u])) stop("'cluster' must not split units that share random intercepts.", call. = FALSE)
    D$str$ucluster <- cu
  } else D$str$ucluster <- NULL
  D$str$G <- G
  D
}

## the full internal parameter vector mapped to the reported (natural) scale
.pc_natural <- function(theta, D) {
  out <- theta
  if (length(D$idx$str)) out[D$idx$str] <- D$str$transform(theta[D$idx$str])
  out
}

## the main builder
.pc_design <- function(formula, data, inflated = FALSE, inflate = NULL, ref = NULL,
                       link = c("logit", "probit"), id = NULL, alt = NULL, weights = NULL,
                       offset = NULL, cluster = NULL, fe = NULL, family = "mnl", str = NULL) {
  link <- match.arg(link)
  if (!is.data.frame(data)) data <- as.data.frame(data)
  if (nrow(data) == 0L) stop("'data' has no rows.", call. = FALSE)
  if (xor(is.null(id), is.null(alt))) stop("A long layout needs both 'id' (case column) and 'alt' (alternative column).", call. = FALSE)
  layout <- if (is.null(id)) "wide" else "long"
  P <- .pc_parts(formula, inflated, layout); F <- P$F
  for (nm in c(id, alt, fe)) if (!is.character(nm) || length(nm) != 1L || !nm %in% names(data))
    stop("'", nm, "' must name a column in 'data'.", call. = FALSE)

  ## ---- aligned auxiliary vectors, then complete cases over every model variable
  wv <- .pc_aligned(weights, data, "weights")
  cv <- .pc_aligned(cluster, data, "cluster")
  ov <- if (layout == "long") .pc_aligned(offset, data, "offset") else NULL
  if (layout == "wide" && !is.null(offset) && (!is.matrix(offset) || nrow(offset) != nrow(data)))
    stop("With one row per case, 'offset' must be an n x J matrix (one column per alternative in level order) or NULL.", call. = FALSE)
  vars <- all.vars(F)
  if ("." %in% vars) stop("A '.' in the formula is not supported; list the covariates explicitly.", call. = FALSE)
  miss <- setdiff(vars, names(data))
  if (length(miss)) stop("Variable(s) not in 'data': ", paste(miss, collapse = ", "), ".", call. = FALSE)
  mf_all <- stats::model.frame(F, data = data, na.action = stats::na.pass)
  ok <- stats::complete.cases(mf_all)
  if (!is.null(wv)) ok <- ok & !is.na(wv)
  if (!is.null(cv)) ok <- ok & !is.na(cv)
  if (!is.null(ov)) ok <- ok & !is.na(ov)
  if (!is.null(fe)) ok <- ok & !is.na(data[[fe]])
  if (family == "remnl") {
    if (is.null(str$unit) || !is.character(str$unit) || length(str$unit) != 1L || !str$unit %in% names(data))
      stop("'unit' must name a column in 'data'.", call. = FALSE)
    ok <- ok & !is.na(data[[str$unit]])
  }
  if (layout == "long") ok <- ok & !is.na(data[[id]]) & !is.na(data[[alt]])
  if (layout == "wide" && !is.null(offset)) ok <- ok & stats::complete.cases(offset)
  ## in the long layout a missing value on one row drops the whole case: dropping
  ## the row alone would silently remove that alternative from the choice set
  n_dropped_cases <- 0L
  if (layout == "long") {
    bad <- unique(data[[id]][!ok & !is.na(data[[id]])])
    if (length(bad)) { n_dropped_cases <- length(bad); ok <- ok & !(data[[id]] %in% bad) }
  }
  n_dropped <- sum(!ok)
  if (!any(ok)) stop("No complete cases: every row has a missing value in a model variable (", paste(vars, collapse = ", "), ").", call. = FALSE)
  rows <- which(ok)
  if (n_dropped_cases)
    message(n_dropped_cases, " case(s) dropped entirely because a model variable is missing on one of their rows (missing values do not change choice sets).")
  d <- data[rows, , drop = FALSE]
  if (!is.null(wv)) wv <- wv[rows]
  if (!is.null(cv)) cv <- cv[rows]
  if (!is.null(ov)) ov <- ov[rows]
  if (layout == "wide" && !is.null(offset)) offset <- offset[rows, , drop = FALSE]
  ## character covariates become factors here, so that their levels are recorded
  ## for prediction (model.matrix would otherwise take the levels present in
  ## whatever data it sees); the response, the id / alt / unit columns are not touched
  lhs_vars <- all.vars(attr(F, "lhs")[[1L]])
  keep_as_is <- c(lhs_vars, id, alt, fe, if (family == "remnl") str$unit)
  for (nm in intersect(vars, names(d))) {
    if (is.character(d[[nm]]) && !nm %in% keep_as_is) d[[nm]] <- factor(d[[nm]])
    if (is.factor(d[[nm]])) d[[nm]] <- droplevels(d[[nm]])
  }

  ## ---- response, alternatives, cases
  yraw <- eval(attr(F, "lhs")[[1L]], d, environment(F))
  if (family == "rologit") {
    ## the rank-ordered logit: the response is the rank; cases whose ranking
    ## carries no information (nothing ranked, or every alternative tied) drop
    ## before anything is indexed, so that no alternative survives only in them
    if (layout != "long") stop("rologit() takes the long layout (one row per case and alternative, with 'id' and 'alt'); the response is the rank.", call. = FALSE)
    rk <- if (is.factor(yraw) || is.logical(yraw)) NA_real_ else suppressWarnings(as.numeric(yraw))
    if (any(!is.finite(rk)) || any(rk < 0))
      stop("The response of rologit() is the rank of each alternative: non-negative numbers, smaller preferred, equal numbers for ties, 0 for an alternative left unranked.", call. = FALSE)
    idv <- as.character(d[[id]])
    info <- vapply(split(rk, idv), function(r) { r1 <- r[r > 0]; length(r1) > 0L && (any(r == 0) || length(unique(r1)) > 1L) }, logical(1))
    if (any(!info)) {
      keep <- idv %in% names(info)[info]
      if (!any(keep)) stop("No case ranks one alternative above another; the rankings carry no information.", call. = FALSE)
      message(sum(!info), " case(s) whose ranking carries no information (no alternative ranked above another) are dropped.")
      d <- d[keep, , drop = FALSE]; rows <- rows[keep]; yraw <- yraw[keep]
      if (!is.null(wv)) wv <- wv[keep]
      if (!is.null(cv)) cv <- cv[keep]
      if (!is.null(ov)) ov <- ov[keep]
      for (nm in intersect(vars, names(d))) if (is.factor(d[[nm]])) d[[nm]] <- droplevels(d[[nm]])
    }
  }
  if (layout == "wide") {
    if (is.character(yraw))
      stop("The response is a character vector; its categories would be ordered alphabetically (",
           paste(sort(unique(yraw)), collapse = ", "), "). Supply a factor with the intended level order, ",
           "e.g. factor(y, levels = c(...)), or integer codes.", call. = FALSE)
    if (is.logical(yraw)) yraw <- factor(yraw, levels = c(FALSE, TRUE))
    yf <- if (is.factor(yraw)) droplevels(yraw) else factor(yraw, levels = sort(unique(yraw)))
    levels <- levels(yf); J <- length(levels)
    y <- as.integer(yf); n <- length(y)
    case_id <- rownames(d)
    ci <- seq_len(n); ji <- y
  } else {
    altf <- d[[alt]]
    if (is.character(altf))
      stop("The alternative column '", alt, "' is a character vector; its categories would be ordered alphabetically (",
           paste(sort(unique(altf)), collapse = ", "), "), which ",
           if (is.null(ref)) "decides the base alternative and " else "", "sets the order of the coefficient blocks and of the prediction columns. ",
           "Supply a factor with the intended level order.", call. = FALSE)
    altf <- if (is.factor(altf)) droplevels(altf) else factor(altf, levels = sort(unique(altf)))
    levels <- levels(altf); J <- length(levels)
    ji <- as.integer(altf)
    cases <- unique(d[[id]]); ci <- match(d[[id]], cases); n <- length(cases)
    case_id <- as.character(cases)
    if (anyDuplicated(cbind(ci, ji))) stop("Some cases have more than one row for the same alternative ('", id, "' x '", alt, "' must be unique).", call. = FALSE)
    if (family == "rologit") {
      ro <- .pc_rank_matrix(as.numeric(yraw), ci, ji, n, J, str$depth)
      Rk <- ro$R; y <- ro$top
    } else {
      ch <- yraw
      if (is.factor(ch)) ch <- as.character(ch)
      if (is.character(ch) || !all(as.numeric(ch) %in% c(0, 1)))
        stop("The response of a long layout must be the 0/1 (or logical) choice indicator.", call. = FALSE)
      ch <- as.numeric(ch)
      nch <- tabulate(ci[ch == 1], nbins = n)
      if (any(nch != 1L)) {
        bad <- case_id[nch != 1L]
        stop(length(bad), " case(s) do not have exactly one chosen alternative (e.g. ", paste(utils::head(bad, 3L), collapse = ", "),
             "); a case whose chosen row had a missing covariate loses that row -- drop such cases from 'data'.", call. = FALSE)
      }
      y <- integer(n); y[ci[ch == 1]] <- ji[ch == 1]
    }
  }
  if (family != "rologit") Rk <- NULL
  if (J < 2L) stop("The response has a single category; a choice model needs at least 2.", call. = FALSE)
  if (J == 2L && inflated) stop("An inflated model with two alternatives is a split-population binary model; imnl() needs at least 3 alternatives.", call. = FALSE)
  base <- if (is.null(ref)) 1L else match(as.character(ref), levels)
  if (is.na(base)) stop("'ref' must be one of the response categories (", paste(levels, collapse = ", "), ").", call. = FALSE)
  ## the inflated (captive) alternatives: k the first, kc all of them; several
  ## take the logit split (a multinomial logit over the regimes, captive_ll.cpp)
  k <- 0L; kc <- integer(0)
  if (inflated) {
    kc <- if (is.null(inflate)) base else match(as.character(inflate), levels)
    if (!length(kc) || anyNA(kc)) stop("'inflate' must be one or more of the response categories (", paste(levels, collapse = ", "), ").", call. = FALSE)
    if (anyDuplicated(kc)) stop("'inflate' names an alternative more than once.", call. = FALSE)
    if (length(kc) > 1L && link != "logit")
      stop("Several inflated alternatives take the logit split equation (the regimes follow a multinomial logit); link = \"probit\" applies to one inflated alternative.", call. = FALSE)
    k <- kc[1L]
  }
  ## an alternative with coefficients of its own must be chosen somewhere (for
  ## rankings: ranked above another available alternative somewhere); checked
  ## once the design is built (below) and again once the cases are final
  cnt <- tabulate(y, nbins = J)

  ## ---- case-level vectors (long: one value per case, constant within case)
  first_row <- match(seq_len(n), ci)
  case_const <- function(v, what) {
    if (is.null(v)) return(NULL)
    if (layout == "wide") return(v)
    vv <- v[first_row]
    if (any(v != vv[ci])) stop("'", what, "' must be constant within a case in a long layout.", call. = FALSE)
    vv
  }
  w <- case_const(wv, "weights")
  if (is.null(w)) w <- rep(1, n) else {
    w <- as.numeric(w)
    if (any(!is.finite(w)) || any(w < 0)) stop("'weights' must be finite and non-negative (frequency weights).", call. = FALSE)
    if (any(abs(w - round(w)) > 1e-8))
      message("Non-integer weights are treated as frequency weights (the weight total is the sample size of the standard errors, ",
              "information criteria, and tests); survey (probability) weights are not supported, see ?mnl.")
  }
  ## a case of weight zero carries nothing: it is dropped below (its
  ## availability row is cleared before the fixed-effects rules, so it counts as
  ## nobody's choice)
  zero_w <- w == 0
  if (all(zero_w)) stop("Every case has weight zero.", call. = FALSE)
  cntw <- tabulate(y[!zero_w], nbins = J)                     # checked with cnt once the design is built
  clus <- case_const(cv, "cluster")
  if (!is.null(clus)) clus <- as.character(clus)
  unit <- if (!is.null(fe)) factor(case_const(d[[fe]], fe)) else NULL

  ## ---- design matrices
  dcase <- d[first_row, , drop = FALSE]
  X <- matrix(0, n, 0L); Z <- matrix(0, n, 0L); W <- matrix(0, n * J, 0L)
  tm <- list(x = NULL, w = NULL, z = NULL)
  asc <- TRUE
  if (!is.null(P$w)) {
    pw <- .pc_part_matrix(F, P$w, d); tm$w <- pw$tm
    Wr <- pw$M
    asc <- "(Intercept)" %in% colnames(Wr)
    tm$w$assign <- attr(Wr, "assign")[colnames(Wr) != "(Intercept)"]       # term of every column, for term_test()
    Wr <- Wr[, colnames(Wr) != "(Intercept)", drop = FALSE]
    tm$w$cols <- colnames(Wr)
    if (ncol(Wr)) {
      W <- matrix(0, n * J, ncol(Wr), dimnames = list(NULL, colnames(Wr)))
      W[(ci - 1L) * J + ji, ] <- Wr
    }
  }
  if (!is.null(P$x)) {
    px <- .pc_part_matrix(F, P$x, dcase); tm$x <- px$tm
    X <- px$M
    if (layout == "wide") asc <- "(Intercept)" %in% colnames(X)
    if (layout == "long") {
      Xall <- .pc_part_matrix(F, P$x, d, tm$x)$M
      if (ncol(Xall) && any(abs(Xall - X[ci, , drop = FALSE]) > 1e-8))
        stop("The case-specific part of the formula contains a variable that varies within a case; ",
             "alternative-specific variables belong in the first part.", call. = FALSE)
      has_int <- "(Intercept)" %in% colnames(X)
      ## the constants follow the first part only; `| 0` (mlogit's way of
      ## dropping them) leaves them in, which a user from mlogit would not expect
      if (asc && !has_int)
        message("The alternative-specific constants are in the model: in the long layout they are switched off by `0 +` (or `- 1`) in the first ",
                "part of the formula, and the second part's intercept has no effect (`| 0` does not drop them, unlike mlogit).")
      if (asc && !has_int) X <- cbind("(Intercept)" = 1, X)
      if (!asc && has_int) X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
  } else if (layout == "long" && asc) {
    X <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  }
  if (inflated) {
    pz <- .pc_part_matrix(F, P$z, dcase); tm$z <- pz$tm
    Z <- pz$M
    if (ncol(Z) == 0L) stop("The split equation has no terms; it needs at least an intercept.", call. = FALSE)
  }
  ## constants, case-specific slopes, a split equation or structural parameters
  ## belong to particular alternatives and are not identified for one never
  ## chosen; a conditional logit with generic coefficients only (no constants:
  ## unlabelled alternatives such as potential cabinets) is identified from the
  ## choices among the available alternatives, and an alternative available
  ## but never chosen informs it like any other. So is a nested logit with
  ## generic coefficients only, except for a nest none of whose alternatives
  ## is chosen (checked once the nests are built, below). In mnl() without
  ## unit fixed effects or inflation, and in rologit(), such an alternative is
  ## not refused: its constant and slopes run to minus infinity, a separation
  ## that the exact check reports (and firth = TRUE keeps finite); the other
  ## models' own parameters of the alternative (scales, nests, variances, the
  ## split equation, unit constants) would not be identified at all
  own <- !(family %in% c("mnl", "nlogit")) || inflated || ncol(X) > 0L
  refuse_unchosen <- own && family != "rologit" && !(family == "mnl" && !inflated && is.null(fe))
  if (refuse_unchosen) {
    if (any(cnt == 0L)) stop("Alternative(s) never chosen on the estimation rows: ", paste(levels[cnt == 0L], collapse = ", "),
                             "; their coefficients are not identified. Drop the level (droplevels) or merge it.", call. = FALSE)
    if (any(cntw == 0L)) stop("Alternative(s) chosen only by cases of weight zero: ", paste(levels[cntw == 0L], collapse = ", "),
                              "; their coefficients are not identified.", call. = FALSE)
  }
  avail <- matrix(0L, n, J)
  if (layout == "wide") avail[] <- 1L else avail[cbind(ci, ji)] <- 1L
  if (any(zero_w)) { avail[zero_w, ] <- 0L; message(sum(zero_w), " case(s) of weight zero are dropped.") }
  off <- matrix(0, n, J)
  if (layout == "wide" && !is.null(offset)) {
    if (ncol(offset) != J) stop("'offset' must have one column per alternative (", J, "), in level order.", call. = FALSE)
    off[] <- as.numeric(offset)
  }
  if (layout == "long" && !is.null(ov)) off[cbind(ci, ji)] <- as.numeric(ov)
  if (family == "hev" && J == 2L)
    warning("With two alternatives the HEV scale is identified only through the functional form of the binary choice probability; treat its estimate with caution.", call. = FALSE)

  ## ---- unit fixed effects and cases left with a single alternative. The
  ## rules switch off the alternatives a unit never chose, which can leave a
  ## case with one alternative; dropping it can change which alternatives the
  ## unit chose, so rule and drop are repeated until nothing changes.
  fe_info <- NULL; fe_map <- NULL; fx <- NULL; fe_off <- NULL
  avail_data <- avail; W0 <- W; n_single <- 0L; n_dropped_fe <- 0L
  repeat {
    if (!is.null(fe)) {
      if (nlevels(unit) < 2L) stop(if (n_dropped_fe) "After dropping the cases with a single available alternative, fewer than two units remain for 'fe'."
                                   else "'fe' has a single unit; nothing to absorb.", call. = FALSE)
      fx <- .pc_fe_columns(unit, y, J, base, levels, n, avail_data, fe, asc = asc, build = FALSE)
      avail <- fx$avail; fe_off <- fx$off
    } else avail <- avail_data
    single <- rowSums(avail) < 2L
    if (!any(single)) break
    keep <- !single
    n_single <- n_single + sum(single & !zero_w); n_dropped_fe <- n_dropped_fe + sum(single)
    X <- X[keep, , drop = FALSE]; Z <- Z[keep, , drop = FALSE]
    W0 <- W0[rep(keep, each = J), , drop = FALSE]
    off <- off[keep, , drop = FALSE]; avail_data <- avail_data[keep, , drop = FALSE]
    y <- y[keep]; w <- w[keep]; zero_w <- zero_w[keep]; if (!is.null(clus)) clus <- clus[keep]
    if (!is.null(Rk)) Rk <- Rk[keep, , drop = FALSE]
    case_id <- case_id[keep]; dcase <- dcase[keep, , drop = FALSE]
    if (layout == "long") { kr <- keep[ci]; d <- d[kr, , drop = FALSE]; ci <- match(ci[kr], which(keep)); ji <- ji[kr]; rows <- rows[kr] }
    else { d <- d[keep, , drop = FALSE]; rows <- rows[keep]; ci <- seq_len(sum(keep)); ji <- y }
    n <- sum(keep)
    if (!is.null(unit)) unit <- droplevels(unit[keep])
    if (n == 0L) stop("No case is left with two or more available alternatives.", call. = FALSE)
  }
  if (!is.null(fe)) fx <- .pc_fe_columns(unit, y, J, base, levels, n, avail_data, fe, asc = asc)   # the dummies, once, at the fixed point
  W <- if (!is.null(fe)) cbind(W0, fx$W) else W0
  if (layout == "long") first_row <- match(seq_len(n), ci)       # case_const() below must see the retained cases
  if (!is.null(fe)) fe_map <- fx$map else fe_off <- NULL
  if (!is.null(fe) && (fx$dropped_pairs || fx$no_base))
    message(fx$dropped_pairs, " unit-alternative pair(s) are never chosen; those alternatives are treated as unavailable ",
            "to those units (their fixed effect is -Inf)",
            if (fx$no_base) paste0("; ", fx$no_base, " unit(s) never choose the base alternative, which is unavailable to them") else "",
            if (n_single) paste0("; ", n_single, " case(s) left with a single alternative are dropped") else "", ".")
  else if (n_single)
    message(n_single, " case(s) with a single available alternative carry no information and are dropped.")
  ## every alternative must be chosen in some case where it competes with
  ## another (captive choices and units that always choose it identify nothing),
  ## and no column may be zero on every retained case (a level present only in
  ## dropped cases); refused where an alternative never chosen is (above). In
  ## a ranking an alternative never ranked above another available one is the
  ## same separation, and without alternative-specific coefficients it
  ## restricts nothing (DESIGN.md, Phase 4b)
  if (refuse_unchosen) {
    cnt2 <- tabulate(y[rowSums(avail) >= 2L], nbins = J)
    if (any(cnt2 == 0L))
      stop("Alternative(s) ", paste(levels[cnt2 == 0L], collapse = ", "), " are chosen only in cases where no other alternative is available",
           if (!is.null(fe)) " (with unit fixed effects: only by units that always choose them)" else "",
           "; their coefficients are not identified.", call. = FALSE)
  }
  zc <- c(colnames(X)[colSums(X != 0) == 0 & colnames(X) != "(Intercept)"], colnames(W0)[colSums(W0 != 0) == 0], colnames(Z)[colSums(Z != 0) == 0 & colnames(Z) != "(Intercept)"])
  if (length(zc)) stop("Column(s) ", paste(zc, collapse = ", "), " are zero on every case of the estimation sample (a variable that never departs from zero ",
                       "there, or a factor level present only in cases that were dropped), so their coefficients are not identified; drop the variable or merge the level.", call. = FALSE)
  if (!is.null(fe)) {
    ## a case-specific covariate constant within every unit, or an
    ## alternative-specific one constant within every unit and alternative, is
    ## a combination of the unit dummies (and the constants): the fixed effects
    ## absorb it (stress finding F15: the party's incumbent candidate with
    ## constituency effects, Arzheimer and Evans 2012)
    ui <- as.integer(unit)
    within_flat <- function(v, g) all(tapply(v, g, function(z) diff(range(z))) <= 1e-10 * max(1, max(abs(v))))
    xs <- setdiff(colnames(X), "(Intercept)")
    flat_x <- xs[vapply(xs, function(cn) within_flat(X[, cn], ui), logical(1))]
    avr <- as.vector(t(avail)) == 1L
    cell <- paste(rep(ui, each = J), rep(seq_len(J), n))[avr]
    flat_w <- colnames(W0)[vapply(seq_len(ncol(W0)), function(c) within_flat(W0[avr, c], cell), logical(1))]
    if (length(c(flat_x, flat_w)))
      stop("With unit fixed effects ('", fe, "'), ", paste(c(flat_x, flat_w), collapse = ", "), if (length(c(flat_x, flat_w)) > 1L) " do" else " does",
           " not vary within units", if (length(flat_w)) " (for each alternative)" else "", ", so the fixed effects absorb ",
           if (length(c(flat_x, flat_w)) > 1L) "them" else "it", "; drop ", if (length(c(flat_x, flat_w)) > 1L) "them" else "it",
           " or the fixed effects.", call. = FALSE)
    Tbar <- stats::median(table(unit))
    fe_info <- list(fe = fe, units = levels(unit), n_units = nlevels(unit), unavailable_pairs = fx$dropped_pairs,
                    no_base_units = fx$no_base_units, cases_dropped = n_dropped_fe, Tbar = Tbar)
    if (Tbar < 10)
      warning(sprintf("Unit fixed effects with few cases per unit (median %d): maximum-likelihood unit constants carry incidental-parameters bias in nonlinear models%s", as.integer(Tbar),
                      if (family == "mnl") "; consider bias_correct = \"jackknife\", or remnl() on the formula that mundlak() returns."
                      else "; consider remnl() on the formula that mundlak() returns."), call. = FALSE)
  }

  ## ---- parameter names, index blocks, identification
  p <- ncol(X); q <- ncol(W); r <- ncol(Z)
  nonbase <- setdiff(seq_len(J), base)
  S <- .pc_structural(family, str, levels, base)
  if (family == "nlogit" && !own) {
    bad <- .pc_unchosen_nests(S, levels, cnt)
    if (length(bad)) stop("No alternative of nest(s) ", paste(bad, collapse = ", "), " is chosen on the estimation rows; the inclusive-value parameter is not identified (it runs to zero). Merge the nest into another or drop its alternatives.", call. = FALSE)
    bad <- .pc_unchosen_nests(S, levels, cntw)
    if (length(bad)) stop("The alternatives of nest(s) ", paste(bad, collapse = ", "), " are chosen only by cases of weight zero; the inclusive-value parameter is not identified.", call. = FALSE)
  }
  beta_names <- if (p) as.vector(vapply(nonbase, function(j) paste0(levels[j], ":", colnames(X)), character(p))) else character(0)
  dn <- .pc_delta_names(Z, kc, levels); rd <- length(dn)          # one split block per inflated alternative
  par_names <- c(beta_names, colnames(W), dn, S$names)
  idx_beta <- if (p) matrix(seq_len(p * (J - 1L)), nrow = p, dimnames = list(colnames(X), levels[nonbase])) else matrix(integer(0), 0L, J - 1L)
  idx_gamma <- if (q) p * (J - 1L) + seq_len(q) else integer(0)
  idx_delta <- if (rd) p * (J - 1L) + q + seq_len(rd) else integer(0)
  idx_str <- if (S$n) p * (J - 1L) + q + rd + seq_len(S$n) else integer(0)
  blocks <- c(rep("outcome", p * (J - 1L)), rep("generic", q), rep("inflation", rd), S$block)
  if (!is.null(fe_map) && nrow(fe_map)) blocks[idx_gamma][colnames(W) %in% fe_map$name] <- "fixed effect"
  if (p) {
    const <- vapply(seq_len(p), function(j) all(X[, j] == X[1L, j]), logical(1)) & colnames(X) != "(Intercept)"
    if (any(const)) stop("Constant covariate(s) in the case-specific part: ", paste(colnames(X)[const], collapse = ", "),
                         "; with the alternative-specific constants they are not identified.", call. = FALSE)
  }
  if (r) .pc_rank_check(Z, "inflation")
  D <- list(layout = layout, F = F, parts = P[c("x", "w", "z")], n = n, J = J, levels = levels, base = base, k = k, kc = kc,
            link = link, link_code = if (link == "logit") 1L else 0L, family = family, str = S,
            y = y, X = X, W = W, Z = Z, off = off, avail = avail, w = w, cluster = clus,
            case_id = case_id, ci = ci, ji = ji, rows = rows, data = d, data_case = dcase,
            id = id, alt = alt, fe = fe, fe_info = fe_info, fe_map = fe_map, fe_off = fe_off, unit = unit, asc = asc,
            tm = tm, p = p, q = q, r = r, K = length(par_names), par_names = par_names, blocks = blocks,
            idx = list(beta = idx_beta, gamma = idx_gamma, delta = idx_delta, str = idx_str),
            n_dropped_rows = n_dropped, n_dropped_cases = n_dropped_cases, has_offset = !is.null(offset),
            offset_arg = if (layout == "long" && is.character(offset)) offset else NULL)
  .pc_identification_check(D)
  if (family == "rologit") { D$str$rank <- Rk; D$str$split <- 0L }
  if (family == "remnl") {
    ## the units sharing random intercepts, the product Gauss-Hermite grid
    unit <- factor(case_const(d[[S$unit]], S$unit))
    if (nlevels(unit) < 2L) stop("'unit' has a single level; random intercepts need several units.", call. = FALSE)
    if (nlevels(unit) == n) warning("Every unit has a single case: the random-intercept variances are identified only through the functional form.", call. = FALSE)
    D$unit <- unit; D$re <- S$unit
    D <- .pc_re_index(D)
    gg <- .pc_gh_grid(S$nodes, J - 1L); grid <- gg$Z
    D$str$Z <- gg$Z; D$str$logw <- gg$logw
    if (nrow(grid) > 20000L) warning("The product Gauss-Hermite grid has ", nrow(grid), " points; reduce 'nodes' or the number of alternatives if the fit is slow.", call. = FALSE)
  }
  if (family == "nlogit") {
    ## a nest whose members are never jointly available to any case, or an
    ## alternative never chosen, cannot identify its lambda; alternatives are
    ## already checked, nests are checked here
    for (m in seq_along(S$nests)) if (S$lmap[m] >= 0L) {
      mem <- match(S$nests[[m]], levels)
      if (!any(rowSums(D$avail[, mem, drop = FALSE]) >= 2L))
        stop("Nest '", names(S$nests)[m], "' never has two of its alternatives available to the same case; its inclusive-value parameter is not identified.", call. = FALSE)
    }
  }
  D
}

## the closures the engine optimizes: all take the full INTERNAL parameter
## vector theta; `hess` is NULL for a family without an analytic Hessian (the
## engine then differences the gradient); `bound` holds the structural
## log-parameters within the search's bound and `search` maps to the search's
## coordinates (.pc_optimize)
.pc_fns <- function(D) {
  f <- .pc_fns_family(D)
  if (isTRUE(D$firth)) return(.pc_firth_fns(f, D))              # mnl() and rologit() only: no structural parameters
  if (!is.null(D$prior)) return(.pc_prior_fns(f, D))            # likewise (DESIGN.md Phase 4c)
  sl <- .pc_str_log(D)
  if (length(sl)) f$bound <- function(th) { th[sl] <- pmin(pmax(th[sl], -.pc_str_bound), .pc_str_bound); th }
  f$search <- .pc_search_map(D)
  f
}

## The search's coordinates under the families with an error scale per
## alternative (hev, mnprobit). Where an alternative's error scale runs to
## infinity, its coefficients grow in proportion (only their ratio to the
## scale is identified in the limit), a ridge that is curved in the
## coefficients and that the search followed for thousands of evaluations
## (stress finding F23). The search therefore measures each non-base
## alternative's coefficients in units of the standard deviation of its
## utility difference from the base, relative to that standard deviation at
## the iid start: sqrt((1 + theta_j^2) / 2) under hev(), where
## Var(theta_j e_j - e_b) is proportional to 1 + theta_j^2, and sqrt(Sigma_jj)
## under mnprobit(), one at the start (I + 11') / 2. Towards a zero scale the
## search's coordinates are the coefficients (up to a constant), towards an
## infinite one their ratio to the scale; the map is smooth and one to one,
## so the maximum is the same. `to` and `from` map the parameters, `grad` the
## gradient at from(u) into the search's coordinates; the structural
## parameters are the same in both.
.pc_search_map <- function(D) {
  if (!(D$family %in% c("hev", "mnp")) || !D$p) return(NULL)
  is <- D$idx$str; m <- D$J - 1L
  if (D$family == "hev") {
    sd <- function(th) sqrt((1 + exp(2 * th[is])) / 2)
    dsd <- function(th) diag(exp(2 * th[is]) / (2 * sd(th)), m)          # d sd_j / d log theta_j
  } else {
    ci <- D$str$chol_idx; de <- D$str$diag_el
    chol_L <- function(th) { L <- matrix(0, m, m); L[1L, 1L] <- 1; v <- th[is]; L[ci] <- ifelse(de, exp(v), v); L }
    sd <- function(th) sqrt(rowSums(chol_L(th)^2))
    dsd <- function(th) {                                                 # d sd_r / d (free element of row r of L)
      L <- chol_L(th); e <- L[ci]
      M <- matrix(0, m, length(is)); M[cbind(ci[, 1L], seq_along(is))] <- ifelse(de, e^2, e) / sqrt(rowSums(L^2))[ci[, 1L]]
      M
    }
  }
  ib <- as.vector(D$idx$beta); p <- D$p
  list(to = function(th) { th[ib] <- th[ib] / rep(sd(th), each = p); th },
       from = function(u) { u[ib] <- u[ib] * rep(sd(u), each = p); u },
       grad = function(u, g) {
         gb <- matrix(g[ib], p); g[ib] <- g[ib] * rep(sd(u), each = p)
         g[is] <- g[is] + as.vector(crossprod(dsd(u), colSums(gb * matrix(u[ib], p))))
         g
       })
}
.pc_fns_family <- function(D) {
  w <- D$w
  if (D$family == "mnl" && length(D$kc) > 1L) {
    a <- list(D$X, D$W, D$Z, D$off, D$avail, D$y - 1L, D$J, D$base - 1L, D$kc - 1L)
    return(list(
      lli    = function(th) do.call(pc_cap_lli_cpp, c(list(th), a)),
      nll    = function(th) -sum(w * do.call(pc_cap_lli_cpp, c(list(th), a))),
      grad   = function(th) -do.call(pc_cap_grad_cpp, c(list(th), a[1:6], list(w), a[7:9])),
      hess   = function(th) -do.call(pc_cap_hess_cpp, c(list(th), a[1:6], list(w), a[7:9])),
      scores = function(th) do.call(pc_cap_scores_cpp, c(list(th), a)),
      probs  = function(th, X = D$X, W = D$W, Z = D$Z, off = D$off, avail = D$avail, fe_off = NULL)
        pc_cap_probs_cpp(th, X, W, Z, off, avail, rep(0L, nrow(X)), D$J, D$base - 1L, D$kc - 1L)))
  }
  if (D$family == "mnl") {
    ka <- .pc_kernel_avail(D$avail, D$off, D$fe_off, D$k)
    a <- list(D$X, D$W, D$Z, ka$off, ka$avail, D$y - 1L, D$J, D$base - 1L, D$k - 1L, D$link_code)
    return(list(
      lli    = function(th) do.call(pc_lli_cpp, c(list(th), a)),
      nll    = function(th) -sum(w * do.call(pc_lli_cpp, c(list(th), a))),
      grad   = function(th) -do.call(pc_grad_cpp, c(list(th), a[1:6], list(w), a[7:10])),   # gradient of nll
      hess   = function(th) -do.call(pc_hess_cpp, c(list(th), a[1:6], list(w), a[7:10])),   # Hessian of nll
      scores = function(th) do.call(pc_scores_cpp, c(list(th), a)),
      probs  = function(th, X = D$X, W = D$W, Z = D$Z, off = D$off, avail = D$avail, fe_off = D$fe_off) {
        kp <- .pc_kernel_avail(avail, off, fe_off, D$k)
        pc_probs_cpp(th, X, W, Z, kp$off, kp$avail, rep(0L, nrow(X)), D$J, D$base - 1L, D$k - 1L, D$link_code) }))
  }
  if (D$family == "rologit") {
    a <- list(D$X, D$W, D$off, D$avail, D$str$rank, D$J, D$base - 1L, D$str$split %||% 0L)
    return(list(
      lli    = function(th) do.call(pc_ro_lli_cpp, c(list(th), a)),                         # one entry per ranking
      nll    = function(th) -sum(w * do.call(pc_ro_lli_cpp, c(list(th), a))),
      grad   = function(th) -do.call(pc_ro_grad_cpp, c(list(th), a[1:5], list(w), a[6:8])),
      hess   = function(th) -do.call(pc_ro_hess_cpp, c(list(th), a[1:5], list(w), a[6:8])),
      scores = function(th) do.call(pc_ro_scores_cpp, c(list(th), a)),
      ## first-choice probabilities: the multinomial logit over the choice set
      probs  = function(th, X = D$X, W = D$W, Z = NULL, off = D$off, avail = D$avail, fe_off = NULL)
        pc_ro_probs_cpp(th, X, W, off, avail, D$J, D$base - 1L)))
  }
  if (D$family == "nlogit") {
    nest <- D$str$nest; lmap <- D$str$lmap
    a <- list(D$X, D$W, D$off, D$avail, D$y - 1L, D$J, D$base - 1L, nest, lmap)
    return(list(
      lli    = function(th) do.call(pc_nl_lli_cpp, c(list(th), a)),
      nll    = function(th) -sum(w * do.call(pc_nl_lli_cpp, c(list(th), a))),
      grad   = function(th) -do.call(pc_nl_grad_cpp, c(list(th), a[1:5], list(w), a[6:9])),
      hess   = NULL,
      scores = function(th) do.call(pc_nl_scores_cpp, c(list(th), a)),
      probs  = function(th, X = D$X, W = D$W, Z = NULL, off = D$off, avail = D$avail, fe_off = NULL)
        pc_nl_probs_cpp(th, X, W, off, avail, rep(0L, nrow(X)), D$J, D$base - 1L, nest, lmap)))
  }
  if (D$family == "hev") {
    a <- list(D$X, D$W, D$off, D$avail, D$y - 1L, D$J, D$base - 1L)
    return(list(
      lli    = function(th) do.call(pc_hev_lli_cpp, c(list(th), a)),
      nll    = function(th) -sum(w * do.call(pc_hev_lli_cpp, c(list(th), a))),
      grad   = function(th) -do.call(pc_hev_grad_cpp, c(list(th), a[1:5], list(w), a[6:7])),
      hess   = NULL,
      scores = function(th) do.call(pc_hev_scores_cpp, c(list(th), a)),
      probs  = function(th, X = D$X, W = D$W, Z = NULL, off = D$off, avail = D$avail, fe_off = NULL)
        pc_hev_probs_cpp(th, X, W, off, avail, rep(0L, nrow(X)), D$J, D$base - 1L)))
  }
  if (D$family == "remnl") {
    S <- D$str; wu <- S$wu
    a <- list(D$X, D$W, D$off, D$avail, D$y - 1L, D$J, D$base - 1L, S$uidx, S$ustart, S$Z, S$logw)
    return(list(
      lli    = function(th) do.call(pc_re_lli_cpp, c(list(th), a)),                       # one entry per unit
      nll    = function(th) -sum(wu * do.call(pc_re_lli_cpp, c(list(th), a))),
      grad   = function(th) -do.call(pc_re_grad_cpp, c(list(th), a[1:5], list(wu), a[6:11])),
      hess   = NULL,
      scores = function(th) do.call(pc_re_scores_cpp, c(list(th), a)),                    # one row per unit
      probs  = function(th, X = D$X, W = D$W, Z = NULL, off = D$off, avail = D$avail, fe_off = NULL, coarse = FALSE, rule_sigma = NULL) {
        ## coarse: the likelihood's own rule (refits whose probabilities are
        ## discarded); rule_sigma: hold the rule chosen at the estimate (numerical
        ## derivatives must not see the rule change between nearby parameters)
        gp <- if (coarse) list(Z = S$Z, logw = S$logw)
              else .pc_re_pred_rule(if (is.null(rule_sigma)) max(exp(th[D$idx$str])) else rule_sigma, S$nodes, D$J - 1L, S)
        pc_re_probs_cpp(th, X, W, off, avail, rep(0L, nrow(X)), D$J, D$base - 1L, S$uidx, S$ustart, gp$Z, gp$logw) }))
  }
  if (D$family == "mnp") {
    a <- list(D$X, D$W, D$off, D$avail, D$y - 1L, D$J, D$base - 1L, D$str$draws)
    return(list(
      lli    = function(th) do.call(pc_mnp_lli_cpp, c(list(th), a)),
      nll    = function(th) -sum(w * do.call(pc_mnp_lli_cpp, c(list(th), a))),
      grad   = function(th) -do.call(pc_mnp_grad_cpp, c(list(th), a[1:5], list(w), a[6:8])),
      hess   = NULL,
      scores = function(th) do.call(pc_mnp_scores_cpp, c(list(th), a)),
      probs  = function(th, X = D$X, W = D$W, Z = NULL, off = D$off, avail = D$avail, fe_off = NULL)
        pc_mnp_probs_cpp(th, X, W, off, avail, rep(0L, nrow(X)), D$J, D$base - 1L, D$str$draws)))
  }
  stop("internal: unknown family ", D$family)
}

## The inflated kernel reads avail[i, k] == 0 as "k is not in the choice set
## of case i": regime 0 (always k) is impossible there and the case is in the
## choice regime with certainty. The fixed-effects rule also switches k off for
## a unit that never chose it, but that is the limit of the unit's k constant
## at -Inf, where regime 0 stays possible (the case keeps its phi_i factor and
## P*_k = 1 - phi_i). Those cells are passed to the kernel as available with
## utility -Inf (probability 0 in the choice regime).
.pc_kernel_avail <- function(avail, off, fe_off, k) {
  if (k == 0L || is.null(fe_off)) return(list(avail = avail, off = off))
  r <- fe_off[, k] & avail[, k] == 0L
  if (!any(r)) return(list(avail = avail, off = off))
  avail[r, k] <- 1L; off[r, k] <- -Inf
  list(avail = avail, off = off)
}

## The rule for the population-averaged probabilities of a random-intercept
## fit, an integral over N(0, sigma^2) that no data narrow. The logistic's
## poles lie at distance pi / sigma in the standard normal variable z, so a
## fixed rule's error grows with sigma: Gauss-Hermite with Q points per
## dimension is off by about exp(-2 pi sqrt(Q) / sigma), the trapezoidal rule
## with step h by about exp(-2 pi^2 / (sigma h)). Up to 7 sigma^2 <= 128 points
## (about 1e-6) the rule is Gauss-Hermite with max(nodes, 7 sigma^2) points,
## taken from the ladder 2^k so that it does not change with every small step
## in sigma (numerical derivatives); beyond, the trapezoidal rule in z on
## [-9, 9] with a power-of-two number of points at least 18 sigma (step at most
## 1 / sigma: error about 1e-8), which needs no eigenproblem and grows only
## linearly with sigma. The product over the dimensions is capped at 5e4
## points (or nodes^m, if larger); `capped` and `err` (the error bound at the
## cap) say when it binds.
.pc_re_pred_rule <- function(sigma, nodes, m, S = NULL) {
  qmax <- max(nodes, floor((5e4)^(1 / m)))                  # never coarser than the likelihood rule
  need <- 10 * sigma^2
  if (need <= 128) {
    q <- max(nodes, 2L^ceiling(log2(max(need, 2))))
    capped <- q > qmax; q <- min(q, qmax)
    if (!capped && !is.null(S) && q == nodes) return(list(Z = S$Z, logw = S$logw, capped = FALSE, err = 0))
    g <- .pc_gh_grid(q, m)
    g$capped <- capped; g$err <- if (capped) exp(-2 * pi * sqrt(q) / sigma) else 0
    return(g)
  }
  q <- 2L^ceiling(log2(18 * sigma + 1))
  capped <- q > qmax; q <- min(q, qmax)
  if (capped && exp(-2 * pi * sqrt(q) / sigma) < exp(-2 * pi^2 * (q - 1) / (18 * sigma))) {
    ## at the cap Gauss-Hermite on q points beats the trapezoid on q points
    g <- .pc_gh_grid(q, m); g$capped <- TRUE; g$err <- exp(-2 * pi * sqrt(q) / sigma)
    return(g)
  }
  z <- seq(-9, 9, length.out = q); lw <- stats::dnorm(z, log = TRUE); lw <- lw - log(sum(exp(lw)))
  grid <- as.matrix(expand.grid(rep(list(seq_len(q)), m)))
  list(Z = matrix(z[grid], nrow(grid)), logw = rowSums(matrix(lw[grid], nrow(grid))), capped = capped,
       err = if (capped) exp(-2 * pi^2 * (q - 1) / (18 * sigma)) else 0)
}

## Gauss-Hermite rules already computed, by size (the eigenproblem is O(q^3))
.pc_gh_cache <- new.env(parent = emptyenv())

## the product probabilists' Gauss-Hermite grid with q points in each of m
## dimensions: nodes (rows of Z) and log weights
.pc_gh_grid <- function(q, m) {
  key <- as.character(q)
  gh <- .pc_gh_cache[[key]]
  if (is.null(gh)) { gh <- .pc_gauss_hermite(q); assign(key, gh, envir = .pc_gh_cache) }
  grid <- as.matrix(expand.grid(rep(list(seq_along(gh$nodes)), m)))
  list(Z = matrix(gh$nodes[grid], nrow(grid)), logw = rowSums(matrix(log(gh$weights[grid]), nrow(grid))))
}

## a design restricted to a subset of cases (bootstrap, half-samples, restricted samples)
.pc_subset_design <- function(D, cases) {
  J <- D$J
  out <- D
  out$n <- length(cases)
  out$X <- D$X[cases, , drop = FALSE]; out$Z <- D$Z[cases, , drop = FALSE]
  out$W <- D$W[as.vector(outer(seq_len(J), (cases - 1L) * J, "+")), , drop = FALSE]
  out$off <- D$off[cases, , drop = FALSE]; out$avail <- D$avail[cases, , drop = FALSE]
  if (!is.null(D$fe_off)) out$fe_off <- D$fe_off[cases, , drop = FALSE]
  out$y <- D$y[cases]; out$w <- D$w[cases]
  if (!is.null(D$cluster)) out$cluster <- D$cluster[cases]
  out$case_id <- D$case_id[cases]
  out$data_case <- D$data_case[cases, , drop = FALSE]
  if (!is.null(D$unit)) out$unit <- D$unit[cases]
  if (D$family == "rologit") out$str$rank <- D$str$rank[cases, , drop = FALSE]
  out$data <- NULL; out$ci <- NULL; out$ji <- NULL; out$rows <- NULL
  if (D$family == "remnl") out <- .pc_re_index(out)
  out
}

## a design with alternative `drop` removed from the choice set (cases that
## chose it are removed too); the base moves to `newbase` (a level index of the
## original design) when the dropped alternative was the base
.pc_drop_alternative <- function(D, drop, newbase = NULL) {
  J <- D$J
  keep_cases <- which(D$y != drop)
  S <- .pc_subset_design(D, keep_cases)
  keep_alt <- setdiff(seq_len(J), drop)
  newJ <- J - 1L
  S$J <- newJ; S$levels <- D$levels[keep_alt]
  S$W <- S$W[rep(keep_alt, S$n) + rep((seq_len(S$n) - 1L) * J, each = newJ), , drop = FALSE]
  S$off <- S$off[, keep_alt, drop = FALSE]; S$avail <- S$avail[, keep_alt, drop = FALSE]
  if (!is.null(S$fe_off)) S$fe_off <- S$fe_off[, keep_alt, drop = FALSE]
  S$y <- match(S$y, keep_alt)
  base_old <- if (D$base == drop) { if (is.null(newbase)) keep_alt[1L] else newbase } else D$base
  S$base <- match(base_old, keep_alt)
  S$k <- if (D$k > 0L) match(D$k, keep_alt) else 0L
  if (D$k > 0L && is.na(S$k)) stop("internal: cannot drop the inflated alternative.")
  ## fixed-effect dummies for the dropped alternative disappear
  if (!is.null(S$fe_map) && nrow(S$fe_map)) {
    gone <- S$fe_map$alt == drop
    if (any(gone)) { S$W <- S$W[, !colnames(S$W) %in% S$fe_map$name[gone], drop = FALSE]; S$fe_map <- S$fe_map[!gone, , drop = FALSE] }
    S$fe_map$alt <- match(S$fe_map$alt, keep_alt)
  }
  S$q <- ncol(S$W)
  nonbase <- setdiff(seq_len(newJ), S$base)
  beta_names <- if (S$p) as.vector(vapply(nonbase, function(j) paste0(S$levels[j], ":", colnames(S$X)), character(S$p))) else character(0)
  S$par_names <- c(beta_names, colnames(S$W), if (S$r) paste0("infl_", colnames(S$Z)))
  S$idx <- list(beta = if (S$p) matrix(seq_len(S$p * (newJ - 1L)), nrow = S$p, dimnames = list(colnames(S$X), S$levels[nonbase])) else matrix(integer(0), 0L, newJ - 1L),
                gamma = if (S$q) S$p * (newJ - 1L) + seq_len(S$q) else integer(0),
                delta = if (S$r) S$p * (newJ - 1L) + S$q + seq_len(S$r) else integer(0),
                str = integer(0))
  S$str <- .pc_structural("mnl", NULL, S$levels, S$base)
  S$blocks <- c(rep("outcome", S$p * (newJ - 1L)), rep("generic", S$q), rep("inflation", S$r))
  if (!is.null(S$fe_map) && nrow(S$fe_map)) S$blocks[S$idx$gamma][colnames(S$W) %in% S$fe_map$name] <- "fixed effect"
  S$K <- length(S$par_names)
  S
}

## designs on new data for predict(): X, W, Z, off, avail aligned to the fit
.pc_design_new <- function(object, newdata, offset = NULL) {
  D <- object$design
  if (!is.data.frame(newdata)) newdata <- as.data.frame(newdata)
  J <- D$J; F <- D$F; P <- D$parts
  if (D$layout == "wide") {
    n <- nrow(newdata); ci <- seq_len(n); ji <- NULL
    dcase <- newdata; dall <- newdata
    avail <- matrix(1L, n, J)
    case_id <- rownames(newdata)
  } else {
    for (nm in c(D$id, D$alt)) if (!nm %in% names(newdata)) stop("'newdata' must contain the '", nm, "' column.", call. = FALSE)
    altf <- factor(newdata[[D$alt]], levels = D$levels)
    if (anyNA(altf)) stop("'newdata' has alternatives that are not in the model: ", paste(unique(as.character(newdata[[D$alt]])[is.na(altf)]), collapse = ", "), ".", call. = FALSE)
    ji <- as.integer(altf)
    cases <- unique(newdata[[D$id]]); ci <- match(newdata[[D$id]], cases); n <- length(cases)
    if (anyDuplicated(cbind(ci, ji))) stop("'newdata' has more than one row for the same case and alternative.", call. = FALSE)
    first_row <- match(seq_len(n), ci)
    dcase <- newdata[first_row, , drop = FALSE]; dall <- newdata
    avail <- matrix(0L, n, J); avail[cbind(ci, ji)] <- 1L
    case_id <- as.character(cases)
  }
  X <- matrix(0, n, 0L); W <- matrix(0, n * J, 0L); Z <- matrix(0, n, 0L)
  if (!is.null(D$tm$w)) {
    Wr <- .pc_part_matrix(F, P$w, dall, D$tm$w)$M
    Wr <- Wr[, colnames(Wr) != "(Intercept)", drop = FALSE]
    if (ncol(Wr)) { W <- matrix(0, n * J, ncol(Wr), dimnames = list(NULL, colnames(Wr))); W[(ci - 1L) * J + ji, ] <- Wr }
  }
  if (!is.null(D$tm$x)) {
    X <- .pc_part_matrix(F, P$x, dcase, D$tm$x)$M
    if (D$layout == "long") {
      has_int <- "(Intercept)" %in% colnames(X)
      if (D$asc && !has_int) X <- cbind("(Intercept)" = 1, X)
      if (!D$asc && has_int) X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
  } else if (D$layout == "long" && D$asc) X <- matrix(1, n, 1L, dimnames = list(NULL, "(Intercept)"))
  if (!is.null(D$tm$z)) Z <- .pc_part_matrix(F, P$z, dcase, D$tm$z)$M
  if (!is.null(D$fe)) {
    if (!D$fe %in% names(dcase)) stop("'newdata' must contain the fixed-effect column '", D$fe, "'.", call. = FALSE)
    u <- as.character(dcase[[D$fe]])
    unseen <- setdiff(unique(u), D$fe_info$units)
    if (length(unseen)) stop("'newdata' has unit(s) not in the fit: ", paste(unseen, collapse = ", "), ".", call. = FALSE)
    W <- cbind(W, .pc_fe_rebuild(D$fe_map, u, J, n))
    ## an alternative never chosen by a unit in the fit stays unavailable to it,
    ## and so does the base for a unit that never chose it
    ufit <- as.character(D$unit); av0 <- avail
    for (j in seq_len(J)) if (j != D$base) avail[!(u %in% unique(ufit[D$avail[, j] == 1L])), j] <- 0L
    if (length(D$fe_info$no_base_units)) avail[u %in% D$fe_info$no_base_units, D$base] <- 0L
    fe_off <- av0 == 1L & avail == 0L
    if (any(rowSums(avail) == 0L))
      stop("'newdata' has case(s) none of whose alternatives was ever chosen by their unit in the fit: ",
           paste(utils::head(case_id[rowSums(avail) == 0L], 5L), collapse = ", "), ".", call. = FALSE)
  } else fe_off <- NULL
  if (ncol(X) != D$p || ncol(W) != D$q || ncol(Z) != D$r) stop("internal: the new design does not match the fit.")
  off <- matrix(0, n, J)
  if (!is.null(offset)) {
    if (D$layout == "wide") { if (!is.matrix(offset) || nrow(offset) != n || ncol(offset) != J) stop("'offset' must be an n x J matrix.", call. = FALSE); off[] <- offset }
    else { ov <- .pc_aligned(offset, newdata, "offset"); off[cbind(ci, ji)] <- as.numeric(ov) }
  } else if (!is.null(D$offset_arg) && D$offset_arg %in% names(newdata)) {
    off[cbind(ci, ji)] <- as.numeric(newdata[[D$offset_arg]])
  } else if (!is.null(D$offset_arg)) stop("The fit has an offset ('", D$offset_arg, "'); 'newdata' must contain it or 'offset' must be given.", call. = FALSE)
  else if (isTRUE(D$has_offset)) stop("The fit has an offset; supply 'offset' for 'newdata' (an n x J matrix in the wide layout, a vector aligned to its rows in the long layout).", call. = FALSE)
  list(X = X, W = W, Z = Z, off = off, avail = avail, fe_off = fe_off, n = n, ci = ci, ji = ji, case_id = case_id)
}

## the constant-share inflated multinomial logit of alternative k on a multinomial logit fit's own
## design (DESIGN.md Phase 4e): the split block, one intercept, added after the
## outcome and generic blocks; the same cases, weights, offsets and choice sets
.pc_inflate_design <- function(D, k) {
  pos <- D$p * (D$J - 1L) + D$q
  Di <- D
  Di$k <- k; Di$kc <- k; Di$link <- "logit"; Di$link_code <- 1L
  Di$Z <- matrix(1, D$n, 1L, dimnames = list(NULL, "(Intercept)")); Di$r <- 1L
  Di$par_names <- append(D$par_names, "infl_(Intercept)", after = pos)
  Di$blocks <- append(D$blocks, "inflation", after = pos)
  Di$idx$delta <- pos + 1L
  if (length(D$idx$str)) Di$idx$str <- D$idx$str + 1L
  Di$K <- D$K + 1L
  Di$tm$z <- NULL; Di$crosscheck <- FALSE
  Di
}

## the names of the split-equation coefficients: infl_<term> with one inflated
## alternative, infl_<alternative>:<term> per alternative with several
.pc_delta_names <- function(Z, kc, levels) {
  if (!ncol(Z) || !length(kc)) return(character(0))
  if (length(kc) == 1L) return(paste0("infl_", colnames(Z)))
  as.vector(vapply(kc, function(c) paste0("infl_", levels[c], ":", colnames(Z)), character(ncol(Z))))
}

## rebuild the parameter names, index blocks, and block labels of a design whose
## X, W, or Z columns were changed (the structural block is kept)
.pc_reindex <- function(D) {
  p <- ncol(D$X); q <- ncol(D$W); r <- if (is.null(D$Z)) 0L else ncol(D$Z); J <- D$J
  nonbase <- setdiff(seq_len(J), D$base)
  nb <- p * (J - 1L)
  beta_names <- if (p) as.vector(vapply(nonbase, function(j) paste0(D$levels[j], ":", colnames(D$X)), character(p))) else character(0)
  ns <- if (is.null(D$str)) 0L else D$str$n
  D$p <- p; D$q <- q; D$r <- r
  dn <- if (r) .pc_delta_names(D$Z, D$kc, D$levels) else character(0); rd <- length(dn)
  D$idx <- list(beta = if (p) matrix(seq_len(nb), nrow = p, dimnames = list(colnames(D$X), D$levels[nonbase])) else matrix(integer(0), 0L, J - 1L),
                gamma = if (q) nb + seq_len(q) else integer(0),
                delta = if (rd) nb + q + seq_len(rd) else integer(0),
                str = if (ns) nb + q + rd + seq_len(ns) else integer(0))
  D$par_names <- c(beta_names, colnames(D$W), dn, if (ns) D$str$names)
  D$blocks <- c(rep("outcome", nb), rep("generic", q), rep("inflation", rd), if (ns) D$str$block)
  if (!is.null(D$fe_map) && nrow(D$fe_map)) D$blocks[D$idx$gamma][colnames(D$W) %in% D$fe_map$name] <- "fixed effect"
  D$K <- length(D$par_names)
  D
}
