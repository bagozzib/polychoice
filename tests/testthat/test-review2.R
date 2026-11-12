## Regression tests for the defects found by the simulated-referee panel
## (round 2): the anchoring of the constants, the weighted fixed-effects bootstrap, boot_lr_test(), the
## nesting checks, rrr() outside the logit family, the identification check
## with large-scale covariates, the probit z tests, zero weights, the
## diagnostics for every standard-error type, prediction edge cases, and the
## random-intercept predictions at large standard deviations.

test_that("fixed effects: constants are anchored through the units' choices; a disconnected alternative is refused by name", {
  set.seed(12); G <- 60; T <- 6; unit <- rep(seq_len(G), each = T); x1 <- rnorm(G * T)
  y <- sample(1:2, G * T, TRUE); grp <- unit > 40
  y[grp] <- sample(2:3, sum(grp), TRUE)          # c is chosen only by units that never choose the base a
  d <- data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1, unit = factor(unit))
  m <- suppressWarnings(suppressMessages(mnl(y ~ x1, data = d, fe = "unit")))
  expect_true(m$converged)
  fn <- polychoice:::.pc_fns(m$design); expect_lt(max(abs(fn$grad(m$theta))), 1e-4)
  ## units 41-60 choose only c: they drop out, and c is chosen by no unit that also chooses another alternative
  y2 <- y; y2[grp] <- 3L
  d2 <- data.frame(y = factor(c("a", "b", "c")[y2], levels = c("a", "b", "c")), x1, unit = factor(unit))
  expect_error(suppressWarnings(suppressMessages(mnl(y ~ x1, data = d2, fe = "unit"))), "Alternative\\(s\\) c are chosen only")
})

test_that("weighted bootstrap with fixed effects: weights follow their rows; constants and unit effects get no replicates", {
  skip_on_cran()                                    # slow; runs in the full suite (check time, 2026-11-06)
  set.seed(13); G <- 80; T <- 6; unit <- rep(seq_len(G), each = T); x1 <- rnorm(G * T)
  V <- cbind(0, 0.3 + 0.8 * x1 + rnorm(G)[unit], -0.2 - 0.6 * x1 + rnorm(G)[unit]); P <- exp(V) / rowSums(exp(V))
  y <- apply(P, 1, function(p) sample.int(3, 1, prob = p))
  d <- data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1, unit = factor(unit), w = sample(1:3, G * T, TRUE))
  m <- suppressWarnings(suppressMessages(mnl(y ~ x1, data = d, fe = "unit", weights = "w", se = "bootstrap", nboot = 20)))
  expect_equal(m$boot$n_failed, 0L)
  se <- sqrt(diag(vcov(m)))
  expect_true(all(is.finite(se[c("b:x1", "c:x1")])))
  expect_true(all(is.na(se[m$blocks == "fixed effect"]))); expect_true(is.na(se[["b:(Intercept)"]]))
  expect_error(first_difference(m, "x1", 0, 1), "no entries for some parameters")
  ci <- confint(m, type = "percentile"); expect_true(all(is.finite(ci[c("b:x1", "c:x1"), ])))
})

test_that("boot_lr_test(): nesting checked, the restriction embedded, no replicate below the restricted fit", {
  skip_on_cran()                                    # slow; runs in the full suite (check time, 2026-11-06)
  set.seed(14); d <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m0 <- mnl(y ~ x1 + x2, data = d); m1 <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | x1, data = d)))
  b <- boot_lr_test(m0, m1, B = 12)
  expect_equal(nrow(b), 2L); expect_true(is.na(b$p.value[1])); expect_equal(attr(b, "n_negative"), 0L)
  expect_true(all(attr(b, "boot_statistics")["lr", ] >= 0))
  expect_true(is.na(b$p.value[2]))                                                  # no normal reference at the boundary
  expect_error(boot_lr_test(mnl(y ~ x2, data = d), mnl(y ~ x1, data = d), B = 10), "not nested")
  expect_error(boot_lr_test(mnl(y ~ x1, data = d), suppressWarnings(mnprobit(y ~ x1 + x2, data = d)), B = 10), "does not nest")
  ## random intercepts with weights other than 1 are refused
  d$u <- factor(rep(1:100, each = 5)); d$w <- rep(1:2, each = 250)
  expect_error(boot_lr_test(mnl(y ~ x1, data = d, weights = "w"), suppressWarnings(remnl(y ~ x1, data = d, unit = "u", weights = "w", nodes = 6)), B = 10),
               "frequency weights with random intercepts")
})

test_that("lr_test(): a common inclusive value nests the free ones; rrr() refuses families without an odds reading", {
  skip_if_not_installed("AER")
  data("TravelMode", package = "AER"); TM <- TravelMode; TM$chosen <- as.integer(TM$choice == "yes")
  nests <- list(n1 = c("air", "train"), n2 = c("bus", "car"))
  ne <- suppressWarnings(suppressMessages(nlogit(chosen ~ wait + gcost | income, data = TM, id = "individual", alt = "mode", nests = nests, equal_lambda = TRUE)))
  nf <- suppressWarnings(suppressMessages(nlogit(chosen ~ wait + gcost | income, data = TM, id = "individual", alt = "mode", nests = nests)))
  lr <- lr_test(ne, nf); expect_equal(lr$df, 1L); expect_equal(lr$statistic, 2 * (nf$loglik - ne$loglik), tolerance = 1e-10)
  nx <- suppressWarnings(suppressMessages(nlogit(chosen ~ wait + gcost | income, data = TM, id = "individual", alt = "mode",
                                                 nests = list(n1 = c("air", "bus"), n2 = c("train", "car")))))
  expect_error(lr_test(ne, nx), "different nests")
  expect_error(rrr(nf), "not odds ratios")
})

test_that("identification check is scale-free; probit covariances are tested against the iid values; zero weights drop", {
  set.seed(15); d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5))); d$big <- d$x2 * 1e7
  m <- mnl(y ~ x1 + big, data = d)
  expect_equal(unname(coef(m)[c("b:big", "c:big")] * 1e7), unname(coef(mnl(y ~ x1 + x2, data = d))[c("b:x2", "c:x2")]), tolerance = 1e-6)
  set.seed(16); dp <- rmnp(1500, beta = cbind(b = c(0.3, 0.8, -0.5), c = c(-0.2, -0.6, 0.7)), sigma = matrix(c(1, .5, .5, 1), 2), gamma = 1)
  mp <- suppressWarnings(mnprobit(chosen ~ w1 | x1 + x2, data = dp, id = "id", alt = "alt"))
  s <- summary(mp)$coefficients
  expect_equal(unname(s["cov:c,b", "z value"]), unname((s["cov:c,b", "Estimate"] - 0.5) / s["cov:c,b", "Std. Error"]), tolerance = 1e-12)
  expect_equal(unname(s["var:c", "z value"]), unname((s["var:c", "Estimate"] - 1) / s["var:c", "Std. Error"]), tolerance = 1e-12)
  d$w <- c(0, 0, rep(1, nrow(d) - 2))
  expect_message(mz <- mnl(y ~ x1, data = d, weights = "w"), "weight zero")
  expect_equal(mz$n, nrow(d) - 2L)
  expect_equal(coef(mz), coef(mnl(y ~ x1, data = d[-(1:2), ])), tolerance = 1e-8)
})

test_that("diagnostics do not depend on the standard errors; prediction edge cases", {
  set.seed(17); d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  expect_equal(mnl(y ~ x1, data = d, se = "none")$rcond, mnl(y ~ x1, data = d)$rcond, tolerance = 1e-10)
  d$yr <- sample(2001:2004, nrow(d), TRUE)
  mf <- mnl(y ~ x1 + factor(yr), data = d)
  expect_true(all(is.finite(first_difference(mf, "x1", 0, 1)$diff)))       # the profile takes yr at its mode
  nd <- d[1:3, ]; nd$x1[2] <- NA
  expect_warning(predict(mf, newdata = nd, type = "class"), "missing covariate")
  pc <- suppressWarnings(predict(mf, newdata = nd, type = "class"))
  expect_true(is.na(pc[2])); expect_false(anyNA(pc[-2]))
})

test_that("random-intercept predictions stay accurate at large standard deviations", {
  skip_on_cran()
  set.seed(1); G <- 150; T <- 6; n <- G * T; unit <- rep(seq_len(G), each = T); x1 <- rnorm(n); u <- rnorm(G, 0, 5)
  p <- plogis(0.4 + 0.9 * x1 + u[unit])
  d <- data.frame(y = factor(ifelse(runif(n) < p, "b", "a"), levels = c("a", "b")), x1, unit = factor(unit))
  m <- suppressWarnings(remnl(y ~ x1, data = d, unit = "unit"))
  V <- polychoice:::.pc_fns(m$design)$probs(m$theta)$V; sg <- coef(m)[["sd:b"]]
  pa <- vapply(V[, 2] - V[, 1], function(e) integrate(function(z) plogis(e + sg * z) * dnorm(z), -Inf, Inf, rel.tol = 1e-12)$value, numeric(1))
  expect_lt(max(abs(fitted(m)[, "b"] - pa)), 1e-5)               # the fixed 12-point rule was off by 1e-2 here
})
