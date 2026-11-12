## the exact separation check, the second start of penalized fits, and the
## ranking refusal removed (DESIGN.md, section 1, Phase 4b)

## two ranking data sets on which the size-and-Wald rule erred (research-note
## checks of 2026-11-05): in the first, z is 1 only for an alternative tied at
## the top, a quasi-complete separation the rule missed (z stopped at 11.4,
## standard error 2e4); the second overlaps, with large finite estimates the
## rule flagged
sep_data <- function(k) {
  if (k == 1L) {
    J <- 3L
    rank <- c(2,1,0, 0,1,1, 1,2,0, 1,0,1, 1,1,0, 0,1,1, 1,0,1, 1,0,2, 1,0,1, 1,0,1)
    x <- c(0.876452, 0.496041, -1.848763, -0.066377, 0.31387, -0.350491, 0.826792, 1.268712, 0.257375, 2.009813,
           0.217542, 0.509293, 1.833373, 0.186765, 0.732526, 0.94163, 0.667248, 0.051365, 0.160492, -0.335944,
           -0.31281, 1.068216, -2.096234, -1.069518, 0.635528, 0.638724, -0.105694, -0.277311, -1.136835, -0.781315)
    z <- replace(numeric(30), 27, 1)
  } else {
    J <- 4L
    rank <- c(0,0,1,0, 0,0,0,1, 1,0,0,0, 0,0,0,1, 0,0,1,0, 0,0,0,1, 0,0,0,1, 0,1,0,0, 0,0,1,0, 0,1,0,0)
    x <- c(1.49884, -0.251995, 1.626804, -0.459435, 0.381983, 0.064831, -2.073914, -0.430614, 2.838549, -2.250825,
           -0.140266, -2.526229, -0.189105, 0.391992, -0.020409, 1.232824, 2.061712, -0.013801, -1.29628, -0.980172,
           -0.272961, 0.477893, 1.366371, 0.492323, 0.727865, -0.120372, 0.833485, 1.798347, -1.497155, -0.238547,
           1.553789, -0.447106, -0.800764, 1.460403, 1.294819, 1.447115, 0.912047, 1.287701, -0.369469, -1.34036)
    z <- c(0,1,1,0, 1,0,0,0, 1,0,0,1, 0,0,0,0, 0,0,0,0, 0,0,0,1, 0,0,0,0, 0,1,0,1, 0,0,0,0, 0,0,0,0)
  }
  data.frame(case = rep(1:10, each = J), alt = factor(rep(seq_len(J), 10)), x = x, z = z, rank = rank)
}

test_that("separation(): the two data sets the size-and-Wald rule got wrong", {
  d1 <- sep_data(1L)
  expect_warning(f1 <- rologit(rank ~ x + z, data = d1, id = "case", alt = "alt"), "quasi-completely separated.*estimates of z do not exist")
  expect_true(f1$boundary)
  s1 <- separation(f1)
  expect_true(s1$separated); expect_identical(s1$type, "quasi-complete"); expect_identical(s1$terms, "z")
  expect_identical(s1$n_separated, 1L); expect_equal(unname(s1$direction["z"]), 1)
  expect_output(print(s1), "Quasi-complete separation: 1 of the 20 comparisons")
  ## the log-likelihood rises along the reported direction from the estimate
  nll <- .pc_fns(f1$design)$nll
  v <- vapply(c(0, 1, 4, 16), function(t) -nll(f1$theta + t * s1$direction), numeric(1))
  expect_true(all(diff(v) >= -1e-12)); expect_gt(v[4], v[1])
  ## overlap with large finite estimates: no boundary
  d2 <- sep_data(2L)
  expect_no_warning(f2 <- rologit(rank ~ x + z, data = d2, id = "case", alt = "alt"))
  expect_false(f2$boundary); expect_true(f2$converged)
  expect_gt(max(abs(f2$theta)), 15)
  s2 <- separation(f2)
  expect_false(s2$separated); expect_identical(s2$type, "overlap"); expect_length(s2$terms, 0L)
  expect_output(print(s2), "Overlap")
})

test_that("separation() agrees with detectseparation on binary data through mnl() and rologit()", {
  skip_on_cran()
  skip_if_not_installed("detectseparation")
  set.seed(20261105)
  agree_m <- agree_r <- n_sep <- 0L; tot <- 0L
  for (r in 1:80) {
    n <- sample(c(8L, 15L, 30L), 1L); x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.2); y <- rbinom(n, 1, plogis(-1 + 1.5 * x1 + 2 * x2))
    if (length(unique(y)) < 2L || length(unique(x2)) < 2L) next
    ds <- suppressWarnings(stats::glm(y ~ x1 + x2, family = stats::binomial, method = detectseparation::detect_separation))$outcome
    fm <- suppressWarnings(mnl(y ~ x1 + x2, data = data.frame(y = factor(y), x1 = x1, x2 = x2)))
    dl <- data.frame(case = rep(seq_len(n), each = 2L), alt = factor(rep(1:2, n)), x1 = as.vector(rbind(0, x1)), x2 = as.vector(rbind(0, x2)),
                     a2 = rep(0:1, n), rank = as.vector(rbind(2L - (y == 0), 2L - (y == 1))))
    fr <- suppressWarnings(rologit(rank ~ 0 + a2 + x1 + x2, data = dl, id = "case", alt = "alt"))
    tot <- tot + 1L; n_sep <- n_sep + ds
    agree_m <- agree_m + (separation(fm)$separated == ds); agree_r <- agree_r + (separation(fr)$separated == ds)
  }
  expect_gt(tot, 50L); expect_gt(n_sep, 5L); expect_lt(n_sep, tot - 5L)
  expect_identical(agree_m, tot); expect_identical(agree_r, tot)
})

test_that("rologit(): an alternative never ranked above another is separation, not refused", {
  set.seed(3)
  dl <- do.call(rbind, lapply(1:60, function(i) {
    w <- rnorm(4L); u <- c(0.5, 0, -0.3, -Inf) + 0.8 * w + stats::rlogis(4L)   # alternative 4 always last
    data.frame(case = i, alt = factor(1:4), w = w, rank = rank(-u))
  }))
  expect_warning(f0 <- rologit(rank ~ w, data = dl, id = "case", alt = "alt"), "quasi-completely separated.*4:\\(Intercept\\) do not exist")
  expect_true(f0$boundary)
  expect_identical(separation(f0)$terms, "4:(Intercept)")
  expect_message(f1 <- rologit(rank ~ w, data = dl, id = "case", alt = "alt", firth = TRUE), "4:\\(Intercept\\) diverge")
  expect_true(f1$converged); expect_true(all(abs(coef(f1)) < 10))
  expect_identical(f1$separation, "4:(Intercept)")
  ## without alternative-specific coefficients it restricts nothing
  expect_no_warning(f2 <- rologit(rank ~ 0 + w, data = dl, id = "case", alt = "alt"))
  expect_false(f2$boundary); expect_false(separation(f2)$separated)
})

test_that("firth = TRUE: the second start finds the higher of two penalized maxima", {
  ## binary logits without intercept (research notes, section 9): on the
  ## separated design the run from zero ends at the lower maximum, on the
  ## overlapping one the run from the maximum-likelihood estimate does
  bin2rank <- function(x, n, s) {
    y <- unlist(Map(function(nn, ss) c(rep(1, ss), rep(0, nn - ss)), n, s)); xx <- rep(x, n)
    do.call(rbind, lapply(seq_along(y), function(i)
      data.frame(case = i, alt = factor(1:2), x = c(xx[i], 0), rank = if (y[i] == 1) c(1L, 2L) else c(2L, 1L))))
  }
  dA <- bin2rank(c(1, 5), c(3, 1), c(0, 0))
  expect_warning(fA <- suppressMessages(rologit(rank ~ 0 + x, data = dA, id = "case", alt = "alt", firth = TRUE)),
                 "two local maxima \\(penalized log-likelihoods -1.051245 from the start at zero and -0.955465")
  expect_equal(unname(fA$theta), -1.922670, tolerance = 1e-6 / 1.92)
  expect_equal(fA$loglik_penalized, -0.95546465, tolerance = 1e-8)
  expect_identical(fA$n_modes, 2L); expect_identical(fA$separation, "x")
  dC <- bin2rank(c(1, 10), c(5, 1), c(1, 0))
  expect_warning(fC <- rologit(rank ~ 0 + x, data = dC, id = "case", alt = "alt", firth = TRUE),
                 "-2.064269 from the start at zero and -2.568373 from the maximum-likelihood estimate")
  expect_equal(unname(fC$theta), -0.157299, tolerance = 1e-6 / 0.157)
  expect_equal(fC$start_logliks, c(-2.06426868, -2.56837275), tolerance = 1e-8)
  expect_length(fC$separation, 0L)
  ## a user start is a single run
  expect_no_warning(fS <- rologit(rank ~ 0 + x, data = dC, id = "case", alt = "alt", firth = TRUE, start = -1.2))
  expect_equal(unname(fS$theta), -1.090573, tolerance = 1e-6)
})

test_that("separation() refuses what it does not cover", {
  expect_error(separation(1), "takes a fit from mnl\\(\\) or rologit\\(\\)")
  m <- suppressWarnings(mnl(y ~ x1, data = data.frame(y = factor(rep(1:3, 10)), x1 = rnorm(30))))
  m$design$k <- 1L
  expect_error(separation(m), "without unit fixed effects or inflation")
})

test_that("bootstrap: a resample without a chooser of an alternative fails a maximum-likelihood replicate and is kept by a penalized one", {
  set.seed(11)
  n <- 40L
  d <- data.frame(x = rnorm(n), y = factor(c(rep("c", 2L), sample(c("a", "b"), n - 2L, TRUE)), levels = c("a", "b", "c")))
  set.seed(5); m0 <- suppressWarnings(mnl(y ~ x, data = d, se = "bootstrap", nboot = 30L))
  expect_gt(m0$boot$n_failed, 0L)
  expect_match(names(m0$boot$fail_reasons)[1L], "is chosen in no case where another alternative is available")
  set.seed(5); m1 <- suppressMessages(suppressWarnings(mnl(y ~ x, data = d, se = "bootstrap", nboot = 30L, firth = TRUE)))
  expect_identical(m1$boot$n_failed, 0L)
  ## rankings: the replicate check runs (a ranked-last alternative fails maximum-likelihood replicates only)
  dl <- do.call(rbind, lapply(1:30, function(i) {
    w <- rnorm(3L); data.frame(case = i, alt = factor(1:3), w = w, rank = if (i == 1L) c(2L, 3L, 1L) else c(rank(-(w[1:2] + stats::rlogis(2L))), 3L))
  }))
  set.seed(5); r0 <- suppressWarnings(rologit(rank ~ w, data = dl, id = "case", alt = "alt", se = "bootstrap", nboot = 20L))
  expect_false(any(grepl("could not find", names(r0$boot$fail_reasons))))
  expect_gt(r0$boot$n_failed, 0L)
  set.seed(5); r1 <- suppressMessages(suppressWarnings(rologit(rank ~ w, data = dl, id = "case", alt = "alt", se = "bootstrap", nboot = 20L, firth = TRUE)))
  expect_identical(r1$boot$n_failed, 0L)
})

test_that("the overlap certificate: its weights give the gradient, and it agrees with the program", {
  ## S'u is the gradient of the log-likelihood at any theta: mnl() with weights
  ## and varying choice sets (long layout), rologit() with partial rankings
  set.seed(8)
  n <- 120L; J <- 4L
  dl <- data.frame(id = rep(seq_len(n), each = J), alt = factor(rep(letters[1:J], n)), w = rnorm(n * J), x = rep(rnorm(n), each = J))
  dl <- dl[!(dl$alt == "d" & dl$id %% 3L == 0L), ]
  u <- 0.7 * dl$w + ifelse(dl$alt == "b", 0.4 * dl$x, 0) + stats::rlogis(nrow(dl))
  dl$chosen <- as.integer(stats::ave(u, dl$id, FUN = function(v) v == max(v)))
  wt <- stats::runif(n, 0.5, 2)
  dl$wt <- wt[dl$id]
  m <- mnl(chosen ~ w | x, data = dl, id = "id", alt = "alt", weights = "wt")
  dl$rank <- stats::ave(-u, dl$id, FUN = function(v) { r <- rank(v); ifelse(r <= 2, r, 0) })     # top-2 partial rankings
  r <- rologit(rank ~ w | x, data = dl, id = "id", alt = "alt", weights = "wt")
  for (f in list(m, r)) {
    D <- f$design
    P <- .pc_sep_pairs(D); S <- .pc_sep_struct(D, P)
    th <- f$theta + stats::rnorm(D$K, sd = 0.3)
    w <- .pc_sep_weights(D, th, P)
    expect_true(all(w > 0))
    expect_equal(drop(crossprod(S, w)), -.pc_fns(D)$grad(th), tolerance = 1e-10)
    ## overlapping data: certified at the estimates, and the program agrees
    s1 <- separation(f); s2 <- .pc_separation(D)
    expect_identical(s1$method, "certificate"); expect_false(s1$separated)
    expect_identical(s2$method, "linear program"); expect_false(s2$separated)
  }
  ## ties: no certificate, the program decides
  dl$rank[dl$rank == 2] <- 1
  rt <- rologit(rank ~ w | x, data = dl, id = "id", alt = "alt")
  expect_null(.pc_sep_weights(rt$design, rt$theta, .pc_sep_pairs(rt$design)))
  expect_identical(separation(rt)$method, "linear program")
  ## separated data are never certified
  expect_identical(separation(suppressWarnings(rologit(rank ~ x + z, data = sep_data(1L), id = "case", alt = "alt")))$method, "linear program")
})
