## Regressions from the replication stress test (verification/stress/FINDINGS.md).

test_that("a nest whose inclusive value runs to zero is a converged boundary fit, not a failed one (F1)", {
  skip_on_cran()
  data(pa2006)
  f <- vote ~ canvass + phone + democrat + republican + age + female + voted2004
  w <- character(0)
  m <- withCallingHandlers(suppressMessages(nlogit(f, data = pa2006, nests = list(base = "abstain", rest = c("democrat", "republican_other", "nonresponse")))),
                           warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_true(m$converged)
  expect_true(m$boundary)
  expect_lt(m$theta[m$idx$str], -7)
  ## the profile log-likelihood rises monotonically to -2123.9129 as lambda -> 0
  ## (-2123.934 at lambda = 0.01, -2123.915 at 0.001): the boundary fit attains it
  expect_gt(m$loglik, -2123.9136)
  expect_true(any(grepl("lambda:rest is at zero or infinity", w)))
  expect_false(any(grepl("did not converge", w)))
})

test_that("Wald and score tests on a singular information report it instead of failing in solve() (F4)", {
  set.seed(3)
  d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d$g <- factor(sample(c("u", "v", "w"), 600, replace = TRUE))
  d$y[d$g != "u" & d$y == "c"] <- "a"                    # c never occurs in two of the groups
  m <- suppressWarnings(mnl(y ~ x1 + x2 + g, data = d))
  ## separated coefficients pushed until their probabilities underflow: no information
  m2 <- m; m2$theta[grep("^c:g", names(m$theta))] <- -60
  expect_error(suppressWarnings(captivity_test(m2)), "information matrix of the multinomial logit fit is singular")   # after its separation warning (DESIGN.md Phase 4d)
  ## a rank-deficient covariance block: the term's Wald statistic is NA, the others stand
  m3 <- m; jj <- grep(":g", names(m$theta)); m3$vcov_internal[jj, jj] <- 1
  expect_warning(tt <- term_test(m3), "covariance of the coefficients of g is singular")
  expect_true(is.na(tt$statistic[tt$term == "g"]))
  expect_equal(tt$statistic[tt$term == "x1"], term_test(m)$statistic[term_test(m)$term == "x1"])
  m4 <- m; m4$vcov_internal[] <- 1
  expect_warning(ct <- combine_test(m4), "is singular")
  expect_true(all(is.na(ct$statistic)))
})

test_that("a flat inflation coefficient is not reported as a saddle (F5)", {
  ## every case with m = 1 chooses the inflated alternative: the m coefficients
  ## diverge and infl_m is flat; before the fix this seed was reported as a saddle
  B <- cbind(b = c(-1.5, 0.4, -0.3), c = c(-2.5, -0.5, 0.2), d = c(-1.8, 0.3, 0.6), e = c(-1.6, -0.2, -0.4))
  set.seed(119)
  d <- rmnl(500, beta = B, delta = c(1.5, 0.5), inflate = "a")
  d$m <- stats::rbinom(500, 1, 0.05); d$y[d$m == 1] <- "a"; d$z2 <- stats::rnorm(500)
  w <- character(0)
  fit <- withCallingHandlers(suppressMessages(imnl(y ~ x1 + x2 + m | z1 + z2 + x1 + m, data = d, inflate = "a", link = "probit")),
                             warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("saddle", w)))
  expect_true(any(grepl("diverge", w)))
  expect_true(fit$converged)
})

test_that("a flat direction of an unidentified inflated model is not reported as a saddle (F18)", {
  skip_on_cran()
  ## one binary covariate and three alternatives: the MNL reproduces the 2 x 3
  ## table, the inflated model has two unidentified directions, and its fits sit
  ## on a curved ridge where the Hessian can have a slightly negative
  ## eigenvalue; before the fix these three seeds were reported as saddles
  sim <- function(seed, n = 1500) {
    set.seed(seed); x <- stats::rbinom(n, 1, 0.5)
    p <- cbind(0.25 + 0.1 * x, 0.35 - 0.05 * x, 0.40 - 0.05 * x)
    data.frame(y = factor(apply(p, 1, function(pr) sample(0:2, 1, prob = pr))), x = x)
  }
  for (s in c(3, 16, 37)) {
    d <- sim(s); w <- character(0)
    fit <- withCallingHandlers(suppressMessages(imnl(y ~ x | x, data = d, inflate = "2", ref = "2", link = "probit")),
                               warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    expect_false(any(grepl("saddle|does not bear out", w)))
    expect_true(fit$converged)
    expect_equal(fit$loglik, suppressWarnings(mnl(y ~ x, data = d, ref = "2"))$loglik, tolerance = 1e-6)
  }
  ## the rule: an ascent along the negative eigenvector is a saddle; a negative
  ## eigenvalue the likelihood contradicts is an inaccurate Hessian; one too
  ## small for the likelihood to rise is neither
  chk <- function(f, Hn) polychoice:::.pc_saddle(f, c(0, 0), c(1, 1), Hn, c(1, 1), c(TRUE, TRUE), 0, 1e-10, 1e-8)$type
  expect_equal(chk(function(t) t[1]^2 - t[2]^2, diag(c(1, -1))), "saddle")
  expect_equal(chk(function(t) t[1]^2 + t[2]^2, diag(c(1, -1))), "inconsistent")
  expect_equal(chk(function(t) t[1]^2 + 10 * t[2]^4 - 5e-7 * t[2]^2, diag(c(1, -1e-6))), "none")
})

test_that("the IIA contrasts of a saturated logit are undefined, not significant (F7)", {
  ## one binary covariate: every estimator reproduces the 2 x 3 table of cell
  ## shares, so the contrasts and their covariances are zero up to rounding
  set.seed(7)
  d <- data.frame(t = rep(0:1, each = 400))
  d$y <- factor(ifelse(stats::runif(800) < ifelse(d$t == 1, 0.5, 0.3), "b", ifelse(stats::runif(800) < 0.6, "a", "c")))
  m <- mnl(y ~ t, data = d)
  expect_warning(hj <- hausman_mcfadden(m, vcov = "joint"), "contrast has no variance")
  expect_true(all(is.na(hj$statistic)))
  expect_warning(fp <- fp_tests(m, tests = "cl"), "contrast has no variance")
  expect_true(is.na(fp$statistic[1]))
})

test_that("the IIA battery warns once when the fit is separated (F8)", {
  set.seed(3)
  d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d$g <- factor(sample(c("u", "v", "w"), 600, replace = TRUE))
  d$y[d$g != "u" & d$y == "c"] <- "a"
  m <- suppressWarnings(mnl(y ~ x1 + x2 + g, data = d))
  w <- character(0)
  r <- withCallingHandlers(iia_tests(m, splits = 2, seed = 1), warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_equal(sum(grepl("diverging coefficients", w)), 1L)
  expect_true(any(grepl("c:gv", w)))
})

test_that("Q1 leaves the diverging directions out of its comparison, wherever the fit stopped", {
  ## the separated design of the test above: no case with g = v or w chooses c,
  ## so c:gv and c:gw each run off on their own (two directions)
  set.seed(3)
  d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d$g <- factor(sample(c("u", "v", "w"), 600, replace = TRUE))
  d$y[d$g != "u" & d$y == "c"] <- "a"
  m <- suppressWarnings(mnl(y ~ x1 + x2 + g, data = d)); K <- length(coef(m))
  expect_warning(f <- fp_tests(m, tests = "cl"), "Q1 leaves the diverging directions out of its comparison")
  det <- attr(f, "details")
  expect_setequal(det$q1_left_out, c("c:gv", "c:gw")); expect_identical(det$q1_left_out_rank, 2L)
  expect_identical(as.integer(f$df), K - 2L)
  i <- match(det$q1_left_out, names(m$theta)); expect_false(anyNA(i))
  expect_lt(max(abs(det$vcov_cl[i, ])), 1e-8)                      # no variance along the directions left out
  expect_equal(unname(det$theta_cl[i]), unname(m$theta[i]))          # the estimator kept the fit's values there
  expect_output(print(f), "Q1 leaves 2 diverging directions out of its comparison")
  ## the statistic does not depend on how far out the optimizer stopped
  for (s in c(4, 10)) {
    m2 <- m; m2$theta[i] <- m2$theta[i] + s * sign(m2$theta[i])
    f2 <- suppressWarnings(fp_tests(m2, tests = "cl"))
    expect_equal(f2$statistic, f$statistic, tolerance = 1e-3); expect_identical(f2$df, f$df)
  }
  expect_lt(f$statistic, 10 * f$df)                                  # not the statistic of a comparison through them (hundreds or more)
  ## a regular fit is compared on every coefficient
  m0 <- mnl(y ~ x1 + x2, data = d); f0 <- fp_tests(m0, tests = "cl")
  expect_null(attr(f0, "details")$q1_left_out); expect_identical(as.integer(f0$df), length(coef(m0)))
})

test_that("a separated reference cell is one direction: Q1 keeps the estimable sum", {
  ## no case with g = u (the reference level) chooses c: c's intercept falls and
  ## c:gv, c:gw rise, their sums with the intercept staying finite
  set.seed(5)
  d <- rmnl(900, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d$g <- factor(sample(c("u", "v", "w"), 900, replace = TRUE))
  d$y[d$g == "u" & d$y == "c"] <- "a"
  m <- suppressWarnings(mnl(y ~ x1 + x2 + g, data = d)); K <- length(coef(m))
  f <- suppressWarnings(fp_tests(m, tests = "cl")); det <- attr(f, "details")
  expect_setequal(det$q1_left_out, c("c:(Intercept)", "c:gv", "c:gw")); expect_identical(det$q1_left_out_rank, 1L)
  expect_identical(as.integer(f$df), K - 1L)
  ## moving along the diverging direction leaves the statistic where it was
  i <- match(c("c:(Intercept)", "c:gv", "c:gw"), names(m$theta))
  m2 <- m; m2$theta[i] <- m2$theta[i] + c(-8, 8, 8)
  f2 <- suppressWarnings(fp_tests(m2, tests = "cl"))
  expect_equal(f2$statistic, f$statistic, tolerance = 1e-3); expect_identical(f2$df, f$df)
  ## the battery's table says what Q1 did
  b <- suppressWarnings(iia_tests(m, splits = 1, seed = 1))
  expect_identical(attr(b, "q1_left_out")$rank, 1L)
  expect_output(print(b), "Q1 leaves 1 diverging direction out of its comparison")
})

test_that("a Wald statistic on a badly scaled covariance is computed, not refused by solve() (F10)", {
  ## variances 1e-10 and 1e10 with correlation 0.5: well posed (normalized
  ## condition number 1/3), but solve() refuses V itself
  V <- matrix(c(1e-10, 0.5, 0.5, 1e10), 2)
  th <- c(1e-5, 1e5)
  expect_error(solve(V, th))
  Vn <- V / tcrossprod(sqrt(diag(V)))
  z <- th / sqrt(diag(V))
  expect_equal(polychoice:::.pc_wald(th, V), as.numeric(crossprod(z, solve(Vn, z))), tolerance = 1e-12)
  expect_true(is.na(polychoice:::.pc_wald(th, matrix(c(1, 1, 1, 1), 2))))
})

test_that("imnl() polishes every screened start and records the link cross-check", {
  skip_on_cran()
  set.seed(2)
  d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), delta = c(0.5, 1.2), inflate = "a")
  expect_null(suppressMessages(imnl(y ~ x1 + x2 | z1, data = d, inflate = "a"))$crosscheck)    # opt-in
  m <- suppressMessages(imnl(y ~ x1 + x2 | z1, data = d, inflate = "a", crosscheck = TRUE))
  ## one polished fit per distinct screened point (a cap on the number polished
  ## made the maximum depend on the screen), then the cross-check and challengers
  expect_gte(length(m$start_logliks), sum(!duplicated(round(m$screen_logliks, 6))))
  expect_equal(m$crosscheck$link, "probit")
  expect_gte(m$loglik, max(m$crosscheck$loglik_own_search, m$crosscheck$loglik_from_other) - 1e-8)
  mp <- suppressMessages(imnl(y ~ x1 + x2 | z1, data = d, inflate = "a", link = "probit", crosscheck = TRUE))
  expect_equal(mp$crosscheck$link, "logit")
})

test_that("the character-alternative refusal names what the level order decides", {
  dl <- rmnl(100, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6)), gamma = 0.7)
  dl$alt <- as.character(dl$alt)
  expect_error(mnl(chosen ~ w1 | x1, data = dl, id = "id", alt = "alt"), "decides the base alternative")
  e <- tryCatch(mnl(chosen ~ w1 | x1, data = dl, id = "id", alt = "alt", ref = "b"), error = conditionMessage)
  expect_match(e, "sets the order of the coefficient blocks")
  expect_false(grepl("base alternative", e))
})

test_that("a conditional logit with generic coefficients only takes alternatives never chosen (F11)", {
  ## unlabelled alternatives (potential cabinets in Baeck et al. 2025): no
  ## constants, so an alternative that is available but never chosen identifies
  ## nothing of its own and informs the generic coefficient like any other
  set.seed(4)
  n <- 300; J <- 4
  d <- data.frame(id = rep(seq_len(n), each = J), alt = factor(rep(letters[1:J], n)), w = stats::rnorm(n * J))
  u <- 0.8 * d$w - ifelse(d$alt == "d", 50, 0) - log(-log(stats::runif(n * J)))
  d$chosen <- as.integer(stats::ave(u, d$id, FUN = function(v) v == max(v)))
  expect_equal(sum(d$chosen[d$alt == "d"]), 0)
  m <- mnl(chosen ~ 0 + w | 1, data = d, id = "id", alt = "alt")
  ll <- function(b) sum(tapply(seq_len(nrow(d)), d$id, function(r) b * d$w[r][d$chosen[r] == 1] - log(sum(exp(b * d$w[r])))))
  o <- stats::optimize(ll, c(-5, 5), maximum = TRUE, tol = 1e-10)
  expect_equal(unname(coef(m)), o$maximum, tolerance = 1e-6)
  expect_equal(m$loglik, o$objective, tolerance = 1e-8)
  ## with constants the never-chosen alternative's constant runs off: a
  ## separation, reported (DESIGN.md Phase 4b); the other coefficients tend to
  ## the fit without that alternative, and the penalized fit is finite
  expect_warning(mc <- mnl(chosen ~ w | 1, data = d, id = "id", alt = "alt"), "quasi-completely separated.*d:\\(Intercept\\) do not exist")
  expect_identical(separation(mc)$terms, "d:(Intercept)")
  m3 <- mnl(chosen ~ w | 1, data = d[d$alt != "d", ], id = "id", alt = "alt")
  expect_equal(coef(mc)[c("b:(Intercept)", "c:(Intercept)", "w")], coef(m3)[c("b:(Intercept)", "c:(Intercept)", "w")], tolerance = 1e-4)
  expect_message(mf <- mnl(chosen ~ w | 1, data = d, id = "id", alt = "alt", firth = TRUE), "d:\\(Intercept\\) diverge")
  expect_true(mf$converged); expect_true(all(abs(coef(mf)) < 15))
  ## other models keep the refusal: the alternative's own parameters would not be identified at all
  expect_error(hev(chosen ~ w | 1, data = d, id = "id", alt = "alt"), "never chosen")
})

test_that("a nested logit with generic coefficients only takes alternatives never chosen, not a nest never chosen (F20)", {
  set.seed(4)
  n <- 300; J <- 4
  d <- data.frame(id = rep(seq_len(n), each = J), alt = factor(rep(letters[1:J], n)), w = stats::rnorm(n * J))
  u <- 0.8 * d$w - ifelse(d$alt == "d", 50, 0) - log(-log(stats::runif(n * J)))
  d$chosen <- as.integer(stats::ave(u, d$id, FUN = function(v) v == max(v)))
  nst <- list(ab = c("a", "b"), cd = c("c", "d"))
  m <- suppressMessages(suppressWarnings(nlogit(chosen ~ 0 + w | 1, data = d, nests = nst, id = "id", alt = "alt")))
  ## the nested-logit log-likelihood written out, on the log scale:
  ## log P(j) = V_j / l_m - I_m + l_m I_m - log sum_k e^(l_k I_k), I_m = log sum over the nest of e^(V / l_m)
  W <- matrix(d$w, n, J, byrow = TRUE); Y <- matrix(d$chosen, n, J, byrow = TRUE); nest <- c(1, 1, 2, 2)
  lse <- function(M) { mx <- apply(M, 1L, max); mx + log(rowSums(exp(M - mx))) }
  ll <- function(p) {
    lam <- exp(p[2:3]); V <- p[1] * W; L <- matrix(lam[nest], n, J, byrow = TRUE)
    I <- sapply(1:2, function(k) lse(V[, nest == k, drop = FALSE] / lam[k]))
    lp <- V / L - I[, nest] + I[, nest] * L - lse(I * rep(lam, each = n))
    sum(lp * Y)
  }
  o <- stats::optim(c(0.8, 0, 0), function(p) -ll(p), method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
  expect_equal(m$loglik, -o$value, tolerance = 1e-6)
  ## with constants the never-chosen alternative's constant runs off: refused
  expect_error(nlogit(chosen ~ w | 1, data = d, nests = nst, id = "id", alt = "alt"), "never chosen")
  ## a nest with no chosen alternative: its inclusive-value parameter runs to zero
  u2 <- 0.8 * d$w - ifelse(d$alt %in% c("c", "d"), 50, 0) - log(-log(stats::runif(n * J)))
  d$chosen2 <- as.integer(stats::ave(u2, d$id, FUN = function(v) v == max(v)))
  expect_error(nlogit(chosen2 ~ 0 + w | 1, data = d, nests = nst, id = "id", alt = "alt"), "No alternative of nest\\(s\\) cd")
  ## with one parameter for every nest, nest ab still identifies it
  m2 <- suppressMessages(suppressWarnings(nlogit(chosen2 ~ 0 + w | 1, data = d, nests = nst, equal_lambda = TRUE, id = "id", alt = "alt")))
  expect_true(is.finite(m2$loglik))
})

test_that("a time limit set by the caller is not swallowed by the package's error handlers", {
  e <- tryCatch({ setTimeLimit(elapsed = 0.2, transient = TRUE); repeat sqrt(2) }, error = identity)
  setTimeLimit()
  expect_error(polychoice:::.pc_pass_limit(e), conditionMessage(e), fixed = TRUE)
  expect_null(polychoice:::.pc_pass_limit(simpleError("Lapack routine dgesv: system is exactly singular")))
  skip_on_cran()
  ## a bootstrap catches each replicate's errors as a failed replicate; the
  ## limit must stop the fit, not one replicate (after which nothing would)
  set.seed(2)
  n <- 2000; x <- stats::rnorm(n); y <- factor(sample(c("a", "b", "c"), n, TRUE))
  d <- data.frame(y, x)
  t0 <- proc.time()[["elapsed"]]
  expect_error({ setTimeLimit(elapsed = 1, transient = TRUE); mnl(y ~ x, data = d, se = "bootstrap", nboot = 5000, cores = 1) }, "time limit")
  setTimeLimit()
  expect_lt(proc.time()[["elapsed"]] - t0, 20)
})

test_that("a covariate absorbed by the unit fixed effects is named (F15), and mundlak() skips unit-constant covariates (F16)", {
  set.seed(15)
  d <- rmnl(600, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6)), gamma = 0.7)
  d$unit <- rep(1:40, each = 15 * 3)[seq_len(nrow(d))]
  ## an alternative attribute constant within unit and alternative (a constituency's incumbent)
  d$inc <- as.numeric((d$unit + as.integer(d$alt)) %% 3 == 0)
  expect_error(suppressWarnings(mnl(chosen ~ w1 + inc | x1, data = d, id = "id", alt = "alt", fe = "unit")),
               "inc does not vary within units .for each alternative., so the fixed effects absorb it")
  ## a case-specific covariate constant within units
  d$z <- d$unit / 10
  expect_error(suppressWarnings(mnl(chosen ~ w1 | x1 + z, data = d, id = "id", alt = "alt", fe = "unit")), "z does not vary within units")
  expect_message(mu <- mundlak(chosen ~ w1 | x1 + z, d, unit = "unit", id = "id"), "Not averaged.*z")
  expect_equal(mu$added, "x1_mean")
})

test_that("a log-parameter at the search's bound with the log-likelihood still rising is named in the boundary warning (F22)", {
  ## fn$grad is the gradient of the negative log-likelihood: at +8 a negative
  ## entry, at -8 a positive one, means the log-likelihood rises out of the bound
  fn <- list(grad = function(th) c(0.5, -0.4, 0.4, -0.3, -1e-9))
  sl <- 2:5
  ## position 4 rises into the bound, position 5 is within the tolerance (1e-3 at |ll| = 100)
  expect_identical(polychoice:::.pc_rising_at_bound(fn, c(1, 8, -8, -8, 8), rep(1, 5), sl, 100), c(2L, 3L))
  expect_identical(polychoice:::.pc_rising_at_bound(fn, c(1, 7.5, -7.5, 0, 0), rep(1, 5), sl, 100), integer())
  msg <- polychoice:::.pc_str_boundary_msg(list(family = "hev", J = 5L), c("scale:3", "scale:4"), "scale:4")
  expect_match(msg, "still rises past the search's bound (+/-8 on the log scale) in scale:4", fixed = TRUE)
  expect_match(msg, "below the model's supremum", fixed = TRUE)
  expect_match(msg, "scale:3, scale:4 are at zero or infinity", fixed = TRUE)
  expect_false(grepl("flat in that direction", msg))
  msg1 <- polychoice:::.pc_str_boundary_msg(list(family = "hev", J = 5L), "scale:3")
  expect_match(msg1, "flat in that direction", fixed = TRUE)
  expect_match(msg1, "scale:3 is at zero or infinity", fixed = TRUE)
})
