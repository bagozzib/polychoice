## DESIGN.md Phase 4g: the Fok-Paap tests in other units, a singular moment
## covariance, the tests on an ill-conditioned fit, Wald tests on a singular
## fit, the polish to stationarity, and profile limits that do not verify.

sim3 <- function(n, seed) {
  set.seed(seed)
  x1 <- rnorm(n, 50, 10); x2 <- runif(n, 1990, 2010); x3 <- rbinom(n, 1, 0.4)
  eta <- cbind(0, 0.03 * (x1 - 50) + 0.4 * x3, -0.05 * (x2 - 2000) - 0.3 * x3, 0.02 * (x1 - 50) + 0.04 * (x2 - 2000))
  y <- factor(apply(exp(eta) / rowSums(exp(eta)), 1, function(p) sample(c("a", "b", "c", "d"), 1, prob = p)))
  data.frame(y, x1, x2, x3)
}

test_that("Q1 and Q2 do not depend on the units or the centering of the covariates", {
  d <- sim3(700, 11)
  d2 <- transform(d, x1 = x1 * 1000, x2 = x2 / 365)
  d3 <- transform(d, x1 = as.numeric(scale(x1)), x2 = as.numeric(scale(x2)))
  ft <- lapply(list(d, d2, d3), function(dd) suppressWarnings(fp_tests(mnl(y ~ x1 + x2 + x3, data = dd))))
  for (k in 2:3) {
    expect_equal(ft[[k]]$statistic, ft[[1]]$statistic, tolerance = 1e-6)
    expect_identical(ft[[k]]$df, ft[[1]]$df)
    expect_equal(attr(ft[[k]], "details")$moment_rcond, attr(ft[[1]], "details")$moment_rcond, tolerance = 1e-4)
  }
  expect_true(all(is.finite(ft[[1]]$statistic)))
})

test_that("a singular moment covariance: Q2's degrees of freedom are its rank minus the parameters", {
  ## two binary covariates: four cells, so one of the nine pairwise moments is a combination of the others,
  ## and Q2 has the two degrees of freedom of the likelihood-ratio test against the saturated model
  set.seed(8); n <- 900; x1 <- rbinom(n, 1, 0.5); x2 <- rbinom(n, 1, 0.5); eta <- cbind(0, 0.4 * x1 - 0.2 * x2, -0.3 * x1 + 0.3 * x2)
  d <- data.frame(x1, x2, y = factor(apply(exp(eta) / rowSums(exp(eta)), 1, function(p) sample(c("a", "b", "c"), 1, prob = p))))
  m <- mnl(y ~ x1 + x2, data = d)
  w <- character(0); ft <- withCallingHandlers(fp_tests(m), warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  dt <- attr(ft, "details")
  expect_identical(dt$n_moments, 9L); expect_identical(dt$moment_rank, 8L); expect_identical(ft$df[2], 2L)
  expect_true(any(grepl("The moment covariance matrix is singular (rank 8 of 9 moments)", w, fixed = TRUE)))
  ms <- mnl(y ~ x1 * x2, data = d)                              # the saturated model
  expect_lt(abs(ft$statistic[2] - 2 * (ms$loglik - m$loglik)), 0.2)   # asymptotically the same statistic
  ## one binary covariate: the model is saturated, no restriction is left, and Q2 is NA
  m1 <- mnl(y ~ x1, data = d)
  w <- character(0); f1 <- withCallingHandlers(fp_tests(m1), warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_identical(attr(f1, "details")$moment_rank, 4L); expect_true(is.na(f1$statistic[2])); expect_identical(f1$df[2], 0L)
  expect_true(any(grepl("and no overidentifying restriction is left; Q2 is not computed (NA).", w, fixed = TRUE)))
})

test_that("the Fok-Paap tests are not computed on an ill-conditioned information", {
  ## a covariate and a near copy of it: identified, and determined to few digits
  d <- sim3(600, 3); set.seed(9); d$x1b <- d$x1 + 1e-3 * rnorm(600)
  m <- suppressWarnings(mnl(y ~ x1 + x1b + x3, data = d))
  expect_true(m$ill_conditioned); expect_false(m$singular)
  expect_warning(ft <- fp_tests(m), "The Fok-Paap tests are not computed (NA): the fit's information matrix is ill-conditioned", fixed = TRUE)
  expect_true(all(is.na(ft$statistic))); expect_true(all(is.na(ft$p.value)))
  expect_lt(attr(ft, "details")$info_rcond, 1e-9)
})

test_that("a Wald test of parameters a singular fit does not identify is NA", {
  ## the saturated design of test-search.R: the split's coefficients are not identified
  set.seed(5); n <- 400; x <- rbinom(n, 1, 0.5)
  d <- data.frame(x = x, y = factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.5, 0.3, 0.2))))
  m <- suppressWarnings(imnl(y ~ x | x, data = d, inflate = c("a", "b")))
  expect_true(m$singular)                                  # on every platform: the fit is polished to stationarity first
  w <- character(0); tt <- withCallingHandlers(term_test(m, type = "wald"), warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") })
  expect_true(any(is.na(tt$statistic)))
  expect_true(any(grepl("involve parameters the data do not identify", w, fixed = TRUE)))
  ## the bootstrap covariance of the same fit has full rank; the rule looks at the information, not at it
  set.seed(1); mb <- suppressWarnings(imnl(y ~ x | x, data = d, inflate = c("a", "b"), se = "bootstrap", nboot = 10))
  if (isTRUE(mb$singular)) expect_true(all(is.na(suppressWarnings(term_test(mb, type = "wald"))$statistic)))
  ## contrasts of a null direction that is not a coordinate: H identifies 3 theta1 + theta2 only
  fake <- list(singular = TRUE, theta = c(0, 0, 0), design = NULL)
  H <- tcrossprod(c(3, 1, 0)) + diag(c(0, 0, 2))
  with_mocked_bindings(.pc_fns = function(D) list(hess = function(th) H), .package = "polychoice", {
    expect_true(polychoice:::.pc_wald_identified(fake, rbind(c(3, 1, 0))))
    expect_true(polychoice:::.pc_wald_identified(fake, rbind(c(0, 0, 1))))
    expect_false(polychoice:::.pc_wald_identified(fake, rbind(c(1 / 3, 1, 0))))
    expect_false(polychoice:::.pc_wald_identified(fake, rbind(c(1, 0, 0))))
  })
  ## an identified fit is untouched
  m0 <- mnl(y ~ x, data = d)
  expect_true(polychoice:::.pc_wald_identified(m0, diag(length(m0$theta))))
  expect_true(all(is.finite(term_test(m0, type = "wald")$statistic)))
})

test_that("the polish to stationarity: an unidentified model's information is singular after it, an identified one's is not", {
  rc <- function(fn, th) { H <- fn$hess(th); dn <- sqrt(abs(diag(H))); dn[!(dn > 0)] <- 1; rcond(H / tcrossprod(dn)) }
  ## not identified: the same covariate twice in a multinomial logit cannot be fitted through the formula
  ## interface, so the saturated inflated model again, moved off its stationary point
  set.seed(5); n <- 400; x <- rbinom(n, 1, 0.5)
  d <- data.frame(x = x, y = factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.5, 0.3, 0.2))))
  m <- suppressWarnings(imnl(y ~ x | x, data = d, inflate = c("a", "b")))
  fn <- polychoice:::.pc_fns(m$design)
  expect_false(is.null(fn$hess))
  set.seed(1); th <- m$theta + 1e-7 * rnorm(length(m$theta))
  pol <- polychoice:::.pc_stationary_polish(fn, th, m$scale)
  expect_false(is.null(pol)); expect_gt(rc(fn, th), 1e-10)
  expect_lte(pol$value, fn$nll(th) + 1e-9)
  expect_lt(pol$max_grad, max(abs(fn$grad(th) / m$scale)))
  expect_lt(rc(fn, pol$par), 1e-13)
  ## identified and well-conditioned: the polish leaves the conditioning where it is
  d2 <- sim3(500, 2); m2 <- mnl(y ~ x1 + x3, data = d2); fn2 <- polychoice:::.pc_fns(m2$design)
  set.seed(1); th2 <- m2$theta + 1e-7 * rnorm(length(m2$theta))
  pol2 <- polychoice:::.pc_stationary_polish(fn2, th2, m2$scale)
  expect_false(is.null(pol2))
  expect_equal(unname(pol2$par), unname(m2$theta), tolerance = 1e-6)
  expect_gt(rc(fn2, pol2$par), 1e-9)
  expect_false(m2$singular); expect_false(m2$ill_conditioned)
})

test_that("profile limits: verified ones are limits, and one that does not verify is NA with a warning", {
  ## an ordinary penalized fit: every limit reproduces the quantile and nothing is warned
  d <- sim3(150, 4); mf <- mnl(y ~ x1 + x3, data = d, firth = TRUE)
  expect_silent(ci <- confint(mf, c("b:x1", "c:x3")))
  expect_true(all(is.finite(ci)))
  fn <- polychoice:::.pc_fns(mf$design); k <- match("b:x1", names(coef(mf))); free <- setdiff(seq_along(mf$theta), k)
  for (cv in ci["b:x1", ]) {
    full <- function(b) { t <- numeric(length(mf$theta)); t[free] <- b; t[k] <- cv; t }
    o <- optim(mf$theta[free], function(b) fn$nll(full(b)), method = "BFGS", control = list(reltol = 1e-14, maxit = 2000))
    expect_equal(2 * (mf$loglik_penalized + o$value), qchisq(0.95, 1), tolerance = 1e-3)
  }
  ## rankings of four items by four judges with seven coefficients: completely separated, and on the
  ## platform this was written on three limits do not verify
  set.seed(4); J <- 4; nj <- 4
  r <- do.call(rbind, lapply(seq_len(nj), function(i) data.frame(id = i, item = seq_len(J), rank = sample(J), matrix(rnorm(J * 7), J, 7))))
  ms <- suppressWarnings(suppressMessages(rologit(rank ~ 0 + X1 + X2 + X3 + X4 + X5 + X6 + X7, data = r, id = "id", alt = "item", firth = TRUE)))
  w <- character(0); cs <- withCallingHandlers(confint(ms), warning = function(c) { w <<- c(w, conditionMessage(c)); invokeRestart("muffleWarning") },
                                               message = function(c) invokeRestart("muffleMessage"))
  ## whatever this sample does, a limit is a number only when it verifies: no "too narrow" intervals are reported
  expect_false(any(grepl("too narrow", w, fixed = TRUE)))
  ## a missing limit always comes with one of the two warnings, and the unverified warning with a missing limit
  expect_true(!anyNA(cs) || any(grepl("could not be verified|too flat there", w)))
  expect_true(!any(grepl("could not be verified and is not reported (NA)", w, fixed = TRUE)) || anyNA(cs))
  skip_on_os(c("mac", "linux", "solaris"))                # where the search ends is platform-dependent on this design
  expect_true(any(grepl("could not be verified and is not reported (NA)", w, fixed = TRUE))); expect_true(anyNA(cs))
})
