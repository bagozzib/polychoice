## priors on the coefficients (prior = coef_prior(...) in mnl() and rologit();
## DESIGN.md, section 1, Phase 4c)

## the log-likelihood plus the log prior written here, at the package's
## estimate: its value and its gradient in the standardized coordinates
.prior_check <- function(m) {
  D <- m$design; s <- m$prior$s; sc <- m$scale; pr <- m$prior
  lp <- function(th) { z <- th / s; if (pr$family == "normal") -0.5 * sum(z^2) else -(pr$df + 1) / 2 * sum(log1p(z^2 / pr$df)) }
  D0 <- D; D0$prior <- NULL; fn0 <- .pc_fns(D0)
  obj <- function(ts) -fn0$nll(ts / sc) + lp(ts / sc)
  list(value = obj(m$theta * sc), grad = numDeriv::grad(obj, m$theta * sc))
}

test_that("coef_prior(): the posterior mode against the definition (mnl, normal, Cauchy, t, scales by name)", {
  data(defection)
  f <- outcome ~ economic * indoctrination + age_joined + birth_year + school_age + male
  for (pr in list(coef_prior("normal"), coef_prior(), coef_prior("t", df = 3), coef_prior("normal", scale = c(2, "switched:economic" = 0.5)),
                  coef_prior("cauchy", autoscale = FALSE))) {
    m <- suppressWarnings(suppressMessages(mnl(f, data = defection, ref = "captured", prior = pr)))
    expect_true(m$converged)
    chk <- .prior_check(m)
    expect_equal(chk$value, m$loglik_penalized, tolerance = 1e-10)
    expect_lt(max(abs(chk$grad)), 1e-5)
    expect_null(m$firth); expect_false(is.null(m$prior))
    expect_true(all(is.finite(coef(m))))
    ## the data are separated (Cook, Niehaus and Zuhlke's Model 3): reported as for Firth's fits
    expect_true("switched:economic" %in% m$separation)
  }
  ## the named scale is that coefficient's own
  mn <- suppressWarnings(suppressMessages(mnl(f, data = defection, ref = "captured", prior = coef_prior("normal", scale = c(2, "switched:economic" = 0.5)))))
  cs <- mn$prior$column_scale
  expect_equal(unname(mn$prior$s["switched:economic"]), 0.5 / unname(cs["switched:economic"]))
  expect_equal(unname(mn$prior$s["switched:male"]), 2 / unname(cs["switched:male"]))
})

test_that("coef_prior(): the scales of the coefficients (autoscale), the constants, and the checks", {
  set.seed(1)
  n <- 200
  d <- data.frame(y = factor(sample(c("a", "b", "c"), n, TRUE)), dummy = rbinom(n, 1, 0.3), five = 5 * rbinom(n, 1, 0.5), cont = rnorm(n, sd = 3))
  m <- mnl(y ~ dummy + five + cont, data = d, prior = coef_prior("normal", scale = 2.5, intercept_scale = 10))
  s <- m$prior$s
  expect_equal(unname(s[c("b:(Intercept)", "c:(Intercept)")]), c(10, 10))
  expect_equal(unname(s[c("b:dummy", "c:dummy")]), c(2.5, 2.5))                 # two values 0 and 1
  expect_equal(unname(s["b:five"]), 2.5 / 5)                                     # two values 0 and 5
  expect_equal(unname(s["c:cont"]), 2.5 / (2 * sd(d$cont)))                      # twice the standard deviation
  ma <- mnl(y ~ dummy + five + cont, data = d, prior = coef_prior("normal", scale = 2.5, autoscale = FALSE))
  expect_equal(unname(ma$prior$s["c:cont"]), 2.5)
  ## a very wide prior reproduces maximum likelihood on overlapping data
  m0 <- mnl(y ~ dummy + five + cont, data = d)
  mw <- mnl(y ~ dummy + five + cont, data = d, prior = coef_prior("normal", scale = 1e6, intercept_scale = 1e6))
  expect_equal(coef(mw), coef(m0), tolerance = 1e-8)
  ## checks
  expect_error(coef_prior("t"), "needs 'df'")
  expect_error(coef_prior("normal", df = 3), "for family = \"t\" only")
  expect_error(coef_prior(scale = -1), "positive")
  expect_error(coef_prior(scale = c(1, 2)), "named vector")
  expect_message(mnl(y ~ dummy, data = d, prior = coef_prior(scale = c(nope = 1))), "does not have: nope")
  expect_error(mnl(y ~ dummy, data = d, prior = list(family = "normal")), "coef_prior")
  expect_error(mnl(y ~ dummy, data = d, firth = TRUE, prior = coef_prior()), "not both")
  expect_output(print(coef_prior("t", df = 4, scale = 1)), "t\\(4\\) prior \\(scale 1, autoscaled\\); constants: scale 10")
})

test_that("coef_prior(): rologit() with ties, rank_test(), and the penalized inference", {
  set.seed(2)
  dr <- rrologit(200, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6), d = c(0.2, 0.4)), gamma = 0.7, ties = 0.3)
  rp <- suppressMessages(rologit(rank ~ w1 | x1, data = dr, id = "id", alt = "alt", prior = coef_prior("t", df = 5, scale = 1)))
  chk <- .prior_check(rp)
  expect_equal(chk$value, rp$loglik_penalized, tolerance = 1e-10); expect_lt(max(abs(chk$grad)), 1e-5)
  rt <- rank_test(rp)
  expect_match(attr(rt, "method"), "t\\(5\\) prior")
  expect_true(all(is.finite(rt$statistic)))
  ci <- confint(rp, parm = c("w1", "b:x1"))
  expect_true(all(is.finite(ci))); expect_true(all(ci[, 1] < coef(rp)[rownames(ci)] & coef(rp)[rownames(ci)] < ci[, 2]))
  tt <- term_test(rp, type = "lr")
  expect_match(attr(tt, "method"), "penalized likelihood ratio, t\\(5\\) prior")
  expect_output(print(rp), "Posterior-mode estimates under a t\\(5\\) prior")
  expect_true(isTRUE(attr(logLik(rp), "penalized")))
  ## two fits with the same prior compare by the penalized likelihood ratio; different penalties are refused
  r0 <- suppressMessages(rologit(rank ~ w1 | 1, data = dr, id = "id", alt = "alt", prior = coef_prior("t", df = 5, scale = 1)))
  lt <- lr_test(r0, rp)
  expect_true(isTRUE(attr(lt, "penalized"))); expect_output(print(lt), "Penalized likelihood-ratio test \\(t\\(5\\) prior")
  expect_error(lr_test(r0, suppressMessages(rologit(rank ~ w1 | x1, data = dr, id = "id", alt = "alt", firth = TRUE))), "penalized differently")
  expect_error(lr_test(r0, suppressMessages(rologit(rank ~ w1 | x1, data = dr, id = "id", alt = "alt", prior = coef_prior("t", df = 5, scale = 2)))),
               "penalized differently")
  expect_error(lr_test(suppressMessages(rologit(rank ~ w1 | 1, data = dr, id = "id", alt = "alt")), rp), "One fit is penalized")
  expect_error(iia_tests(mnl(y ~ x1, data = rmnl(200, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6))), prior = coef_prior())), "a prior")
})

test_that("coef_prior(): weighted scales, a named scale the restricted model lacks, and the bootstrap (DESIGN.md Phase 4d)", {
  set.seed(3)
  n <- 150
  d <- data.frame(y = factor(sample(c("a", "b", "c"), n, TRUE)), x1 = rnorm(n), x2 = rnorm(n), dm = rbinom(n, 1, 0.4), wt = rep(1:3, length.out = n))
  dd <- d[rep(seq_len(n), d$wt), ]
  pr <- coef_prior("normal")
  mw <- mnl(y ~ x1 + dm, data = d, weights = "wt", prior = pr); md <- mnl(y ~ x1 + dm, data = dd, prior = pr)
  expect_equal(mw$prior$s, md$prior$s, tolerance = 1e-12)             # the frequency-weighted column scales
  expect_equal(coef(mw), coef(md), tolerance = 1e-8)
  ## a scale named for x2 in the full model: the restricted model ignores it (a message), and the penalties agree on the shared coefficients
  pn <- coef_prior("normal", scale = c(2.5, "b:x2" = 1))
  f <- mnl(y ~ x1 + x2, data = d, prior = pn)
  expect_message(mnl(y ~ x1, data = d, prior = pn), "does not have: b:x2")
  r <- suppressMessages(mnl(y ~ x1, data = d, prior = pn))
  lt <- lr_test(r, f)
  expect_true(isTRUE(attr(lt, "penalized"))); expect_true(is.finite(lt$statistic))
  r2 <- mnl(y ~ x1, data = d, prior = coef_prior("normal"))
  expect_equal(lr_test(r2, f)$statistic, lt$statistic, tolerance = 1e-10)
  expect_error(lr_test(mnl(y ~ x1, data = d, prior = coef_prior("normal", scale = 1)), f), "penalized differently")
  ## a prior fit's bootstrap keeps resamples with a never-chosen alternative, as Firth's does
  d3 <- data.frame(y = factor(c(rep("a", 40), rep("b", 40), "c", "c"), levels = c("a", "b", "c")), x = rnorm(82))
  set.seed(4); mb <- suppressWarnings(suppressMessages(mnl(y ~ x, data = d3, prior = coef_prior(), se = "bootstrap", nboot = 20)))
  expect_equal(mb$boot$n_failed, 0)
  expect_output(print(fit_stats(mw)), "Penalized fit \\(normal prior")
})

test_that("coef_prior(): the covariance leaves out the prior's negative curvature where the posterior's is not positive definite", {
  Hpd <- diag(2); Hnot <- matrix(c(1, 0, 0, -0.5), 2)
  fn <- list(vcov_info = function(th) Hnot, info = function(th) Hpd, scores = function(th) matrix(0, 1, 2))
  V <- .pc_vcov(c(0, 0), fn, "analytic", NULL, c(1, 1), 1)
  expect_true(isTRUE(attr(V, "prior_clipped"))); expect_equal(unclass(.pc_strip_vcov(V)), diag(2), ignore_attr = TRUE)
  fn$vcov_info <- function(th) 2 * Hpd
  V2 <- .pc_vcov(c(0, 0), fn, "analytic", NULL, c(1, 1), 1)
  expect_null(attr(V2, "prior_clipped")); expect_equal(unclass(.pc_strip_vcov(V2)), diag(2) / 2, ignore_attr = TRUE)
})

test_that("a prior's autoscaling refuses weights that total 2 or less, and names near matches of absent scales (round 6, I-2 / D-new-7)", {
  set.seed(45)
  d <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d$wt <- 1 / nrow(d)                                                  # normalized weights, total 1
  expect_error(suppressMessages(mnl(y ~ x1 + x2, data = d, weights = "wt", prior = coef_prior())), "total above 2")
  expect_message(mnl(y ~ x1 + x2, data = d, prior = coef_prior(scale = c(5, "b:X1" = 1))), "did you mean b:x1")
  ## a term the model leaves out (a restricted model under the full model's
  ## prior) gets no hint: b:x2 is not a typo of b:x1 (round 7, D-5)
  msgs <- character(0)
  withCallingHandlers(mnl(y ~ x1, data = d, prior = coef_prior(scale = c(5, "b:x2" = 1))),
                      message = function(m) { msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("b:x2", msgs))); expect_false(any(grepl("did you mean", msgs)))
})
