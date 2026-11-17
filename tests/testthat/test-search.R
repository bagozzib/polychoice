## The record of the search's ends (object$search, object$search_theta) and
## the count of distinct local maxima (n_modes): every end classified by the
## rules of the fit-quality flags, the reported fit the highest end.

test_that("a separated logit has no maximum: n_modes is 0 and its one end is a boundary", {
  set.seed(3)
  d <- data.frame(x = rnorm(200)); d$y <- factor(ifelse(d$x > 0.8, "c", sample(c("a", "b"), 200, TRUE)))
  m <- suppressWarnings(mnl(y ~ x, data = d))
  expect_true(m$boundary)
  expect_identical(m$n_modes, 0L)
  expect_identical(m$search$type, "boundary"); expect_true(m$search$reported); expect_identical(m$search$reached, 1L)
  expect_match(m$search$detail, "diverging: c:(Intercept), c:x", fixed = TRUE)
  expect_identical(polychoice:::glance.choicemodel(m)$n_modes, 0L)
  data(fluvax)
  mo <- mnl(vaccinate ~ condition + party, data = fluvax)
  expect_identical(mo$n_modes, 1L); expect_identical(mo$search$type, "maximum"); expect_false(mo$singular)
  expect_identical(unname(mo$search_theta[1, ]), unname(mo$theta))
  ## a penalized fit counts the maxima of its two runs and keeps no record
  mf <- suppressMessages(mnl(y ~ x, data = d, firth = TRUE))
  expect_null(mf$search); expect_identical(mf$n_modes, 1L)
})

test_that("fluvax: an end at which a split coefficient diverges is not counted as a maximum", {
  skip_on_cran()
  data(fluvax)
  f <- vaccinate ~ condition + age + female + educ + income + race + party + knowledge | age + educ + knowledge
  expect_message(mi <- imnl(f, data = fluvax, inflate = "dont_know"), "2 distinct local maxima.*1 point that is not a maximum")
  expect_identical(mi$n_modes, 2L)
  expect_identical(mi$search$type, c("maximum", "boundary", "maximum"))
  expect_identical(mi$search$reported, c(TRUE, FALSE, FALSE))
  expect_equal(mi$search$loglik, c(-1276.86985, -1278.59009, -1279.57696), tolerance = 1e-5 / 1276)
  expect_match(mi$search$detail[2], "infl_knowledge", fixed = TRUE)
  ## every polished end is in one row; the reported row holds the fit's estimates
  expect_identical(sum(mi$search$reached), sum(is.finite(mi$start_logliks)))
  expect_identical(unname(mi$search_theta[1, ]), unname(mi$theta))
  expect_identical(colnames(mi$search_theta), mi$par_names)
  expect_output(print(summary(mi)), "2 distinct local maxima among the ends of the search from 14 starts")
  ## the lower maximum, refitted from its row
  m3 <- imnl(f, data = fluvax, inflate = "dont_know", start = mi$search_theta[3, ])
  expect_equal(m3$loglik, mi$search$loglik[3], tolerance = 1e-9)
  expect_identical(m3$n_modes, 1L)
})

test_that("the record tells maxima apart by their estimates, not by a rounded log-likelihood", {
  skip_on_cran()
  ## the captive design without an exclusion restriction (verification/captive_mc.R,
  ## n = 2000, replication 2): two maxima, -2367.099 and -2368.577
  beta <- cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)); delta <- cbind(a = c(1.0, 0.8), c = c(1.8, -0.7))
  set.seed(20261019 + 5000 * 3 + 2)
  X <- data.frame(x1 = rnorm(2000), x2 = rnorm(2000))
  d <- rmnl(2000, beta = beta, delta = delta, inflate = c("a", "c"), x = X, z = data.frame(z1 = X$x1))
  m <- suppressMessages(imnl(y ~ x1 + x2 | x1, data = d, inflate = c("a", "c")))
  expect_identical(m$n_modes, 2L); expect_identical(m$search$type, c("maximum", "maximum"))
  D <- m$design; fn <- polychoice:::.pc_fns(D); sc <- m$scale
  op1 <- list(par = m$theta, value = -m$loglik, converged = TRUE, max_grad = 0)
  rp <- list(type = "maximum", detail = "", se = sqrt(diag(m$vcov_internal)))
  ## the second maximum given a log-likelihood 0.005 below the first (one
  ## value at 0.1 rounding): its estimates are far apart, so it is a second maximum
  op2 <- list(par = m$search_theta[2, ], value = -m$loglik + 0.005, converged = TRUE, max_grad = 0)
  s <- polychoice:::.pc_search_record(D, fn, sc, D$w, list(op1, op2), 1L, op1, rp, 1e-10)
  expect_identical(s$table$type, c("maximum", "maximum")); expect_identical(s$table$reached, c(1L, 1L))
  ## an end a millionth of a standard error from the fit is the same point
  op3 <- list(par = m$theta + 1e-6 * rp$se, value = -m$loglik + 1e-9, converged = TRUE, max_grad = 0)
  s1 <- polychoice:::.pc_search_record(D, fn, sc, D$w, list(op1, op3), 1L, op1, rp, 1e-10)
  expect_identical(nrow(s1$table), 1L); expect_identical(s1$table$reached, 2L); expect_identical(s1$table$detail, "")
})

test_that("ends on a ridge of non-identification are one maximum: a saturated design", {
  ## one binary covariate, three alternatives, two of them inflated with the
  ## covariate in the split: more parameters than cell probabilities, so the
  ## search ends at different estimates with the same fitted probabilities
  ## (as in DiGiuseppe et al. 2024); grouped by the estimates alone, these were
  ## eleven "distinct" maxima at one log-likelihood
  set.seed(1)
  x <- rbinom(600, 1, 0.5)
  p <- cbind(a = 0.5 - 0.2 * x, b = 0.3 + 0.1 * x, c = 0.2 + 0.1 * x)
  d <- data.frame(y = factor(apply(p, 1, function(pr) sample(c("a", "b", "c"), 1, prob = pr))), x = x)
  ## the information is singular at the maximum itself; where the search
  ## stopped its reciprocal condition number was below 1e-15 on Windows and
  ## 2.5e-12 under R 4.5.3 on Linux, so the fit is polished to stationarity
  ## before the flag is set (DESIGN.md Phase 4g) and is singular on both
  expect_warning(m <- imnl(y ~ x | x, data = d, inflate = c("a", "b")), "not identified by the data")
  expect_true(m$ill_conditioned); expect_lt(m$rcond, 1e-9)
  ## the summary's rows for the unidentified parameters look like the others,
  ## so its footer says that their estimates and standard errors mean nothing
  expect_true(m$singular); {
    expect_output(print(summary(m)), "singular (a generalized inverse was used): some parameters are not identified by the data, and their estimates and standard errors are not meaningful", fixed = TRUE)
    expect_match(m$search$detail, "ends at other estimates with the same fitted probabilities: the estimates are not identified", fixed = TRUE)
  }
  expect_identical(m$n_modes, 1L); expect_identical(m$search$type, "maximum")
  expect_identical(m$search$reached, sum(is.finite(m$start_logliks)))
  ## two of the ends, refitted: other estimates, the same probabilities
  D <- m$design; fn <- polychoice:::.pc_fns(D)
  op1 <- list(par = m$theta, value = -m$loglik, converged = TRUE, max_grad = 0)
  m2 <- suppressWarnings(imnl(y ~ x | x, data = d, inflate = c("a", "b"), start = m$theta + c(rep(0, 4), 0.5, 0, -0.5, 0)))
  expect_equal(m2$loglik, m$loglik, tolerance = 1e-8)
  expect_gt(max(abs(m2$theta - m$theta) * m$scale), 0.01)
  expect_lt(max(abs(fn$probs(m2$theta)$prob - fn$probs(m$theta)$prob)), 1e-4)
  op2 <- list(par = m2$theta, value = -m2$loglik, converged = TRUE, max_grad = 0)
  rp <- list(type = "maximum", detail = "", se = NULL)
  s <- polychoice:::.pc_search_record(D, fn, m$scale, D$w, list(op1, op2), 1L, op1, rp, 1e-10)
  expect_identical(s$table$reached, 2L)
  expect_match(s$table$detail, "^1 end at other estimates with the same fitted probabilities")
})

test_that("a reported fit at a boundary with a lower interior maximum: the message names the maximum", {
  skip_on_cran()
  ## replication 6 of the same design: the highest end diverges, a lower one is a maximum
  beta <- cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)); delta <- cbind(a = c(1.0, 0.8), c = c(1.8, -0.7))
  set.seed(20261019 + 5000 * 3 + 6)
  X <- data.frame(x1 = rnorm(2000), x2 = rnorm(2000))
  d <- rmnl(2000, beta = beta, delta = delta, inflate = c("a", "c"), x = X, z = data.frame(z1 = X$x1))
  expect_message(m <- suppressWarnings(imnl(y ~ x1 + x2 | x1, data = d, inflate = c("a", "c"))),
                 "is at a boundary \\(see the warning\\); the search also reached 1 local maximum below it")
  expect_true(m$boundary)
  expect_identical(m$n_modes, 1L)
  expect_identical(m$search$type[m$search$reported], "boundary")
  expect_identical(sum(m$search$type == "maximum"), 1L)
  expect_lt(max(m$search$loglik[m$search$type == "maximum"]), m$loglik)
  expect_output(print(suppressWarnings(summary(m))), "is not a local maximum; the search from 8 starts also reached 1 below it")
})

test_that("ill-conditioned maxima within 0.01 in log-likelihood are one point", {
  ## the ends' types are set by hand (.pc_end_type mocked), on a small fit: two
  ## ends at different estimates and different fitted probabilities
  set.seed(4); d <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | x1, data = d, inflate = "a")))
  D <- m$design; fn <- polychoice:::.pc_fns(D); sc <- polychoice:::.pc_scales(D)
  par2 <- m$theta; par2[1] <- par2[1] + 0.5
  expect_gt(max(abs(fn$probs(par2)$prob - fn$probs(m$theta)$prob)), 1e-3)
  op1 <- list(par = m$theta, value = -m$loglik, converged = TRUE, max_grad = 0)
  rec <- function(detail, gap, rep_detail = detail) {
    testthat::local_mocked_bindings(.pc_end_type = function(...) list(type = "maximum", detail = detail, se = NULL), .package = "polychoice")
    op2 <- list(par = par2, value = -m$loglik + gap, converged = TRUE, max_grad = 0)
    polychoice:::.pc_search_record(D, fn, sc, D$w, list(op1, op2), 1L, op1, list(type = "maximum", detail = rep_detail, se = NULL), 1e-10)$table
  }
  ill <- "ill-conditioned information"
  s <- rec(ill, 0.005)
  expect_identical(nrow(s), 1L); expect_identical(s$reached, 2L)
  expect_identical(s$detail, "ill-conditioned information; 1 other ill-conditioned end within 0.01 in log-likelihood, not counted apart")
  expect_identical(nrow(rec(ill, 0.02)), 2L)                      # further apart than 0.01: two points
  expect_identical(nrow(rec("", 0.005)), 2L)                      # well-conditioned maxima stay apart
  expect_identical(nrow(rec(ill, 0.005, rep_detail = "")), 2L)    # and so does an ill-conditioned end beside a well-conditioned point
})

test_that("a reported fit at a boundary: ill-conditioned maxima below it are one point", {
  set.seed(4); d <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | x1, data = d, inflate = "a")))
  D <- m$design; fn <- polychoice:::.pc_fns(D); sc <- polychoice:::.pc_scales(D)
  testthat::local_mocked_bindings(.pc_end_type = function(...) list(type = "maximum", detail = polychoice:::.pc_ill_detail, se = NULL), .package = "polychoice")
  op1 <- list(par = m$theta, value = -m$loglik, converged = TRUE, max_grad = 0)
  mk <- function(shift, gap) { p <- m$theta; p[1] <- p[1] + shift; list(par = p, value = -m$loglik + gap, converged = TRUE, max_grad = 0) }
  ends <- list(op1, mk(0.5, 0.020), mk(1.0, 0.024), mk(1.5, 0.027), mk(2.0, 0.050))
  s <- polychoice:::.pc_search_record(D, fn, sc, D$w, ends, 1L, op1, list(type = "boundary", detail = "diverging: x", se = NULL), 1e-10)
  expect_identical(s$table$type, c("boundary", "maximum", "maximum")); expect_identical(s$table$reached, c(1L, 3L, 1L))
  expect_identical(s$table$detail[1], "diverging: x")
  expect_identical(s$table$detail[2], "ill-conditioned information; 2 other ill-conditioned ends within 0.01 in log-likelihood, not counted apart")
  expect_identical(nrow(s$theta), 3L); expect_equal(unname(s$theta[2, 1]), unname(m$theta[1]) + 0.5)   # the highest end's estimates are kept
})

test_that("the inverse of a singular matrix: the rule is on the normalized matrix, at 1e-13", {
  inv <- polychoice:::.pc_inverse
  Q <- qr.Q(qr(matrix(c(2, 1, 0, 1, 3, 1, 0, 1, 4), 3)))
  mat <- function(small) Q %*% diag(c(1, 0.5, small)) %*% t(Q)
  a <- inv(mat(1e-11)); expect_false(a$ginv); expect_equal(a$inv, solve(mat(1e-11)), tolerance = 1e-6)
  b <- inv(mat(1e-15)); expect_true(b$ginv); expect_lt(b$rcond, 1e-13); expect_identical(b$rank, 2L); expect_identical(a$rank, 3L)
  expect_lt(max(abs(b$inv)), 1e3)                                  # the null direction is dropped, not inverted to 1e15
  ## badly scaled but regular: the same matrix in other units is not singular
  s <- c(1, 1e-9, 1e9); M <- mat(0.3) * tcrossprod(s)
  r <- inv(M); expect_false(r$ginv); expect_equal(r$inv * tcrossprod(s), solve(mat(0.3)), tolerance = 1e-8)
  ## and a singular one stays singular in other units
  expect_true(inv(mat(1e-15) * tcrossprod(s))$ginv)
  ## a statistic's weight matrix: with the tolerance sqrt(eps) a direction at 1e-10 is dropped too, once the matrix is singular
  Q4 <- qr.Q(qr(matrix(c(2, 1, 0, 1, 1, 3, 1, 0, 0, 1, 4, 1, 1, 0, 1, 5), 4))); M4 <- Q4 %*% diag(c(1, 0.5, 1e-10, 1e-15)) %*% t(Q4)
  expect_identical(inv(M4)$rank, 3L); expect_identical(inv(M4, gtol = sqrt(.Machine$double.eps))$rank, 2L)
})
