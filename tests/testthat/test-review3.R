## Regression tests for the defects found by the simulated-referee panel
## (round 3): the inflated model's never-chosen rule under unit fixed effects,
## the fixed-effects rules iterated to a fixed point, fixed effects without
## alternative-specific constants, the fixed-effects bootstrap (units
## resampled, failed replicates counted), the random-intercept prediction
## rule, the scale-free identification and scaling, the comparison guards,
## the Fok-Paap Q1 covariance, and the reporting of the IIA tests.

fe_panel <- function(G, T, seed, phi = NULL, sd = 2) {
  set.seed(seed); unit <- rep(seq_len(G), each = T); n <- G * T
  x1 <- rnorm(n); z1 <- rnorm(n)
  V <- cbind(0, 1 + 0.8 * x1 + rnorm(G, sd = sd)[unit], 0.5 - 0.6 * x1 + rnorm(G, sd = sd)[unit]); P <- exp(V) / rowSums(exp(V))
  y <- apply(P, 1, function(p) sample.int(3, 1, prob = p))
  if (!is.null(phi)) y <- ifelse(runif(n) < plogis(qlogis(phi) + 0.8 * z1), y, 1L)
  data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1, z1, unit = factor(unit))
}

test_that("fixed effects: the rules and the single-alternative drop are iterated; captive and always-chosen alternatives are refused by cause", {
  ## long layout: unit 1's only b choice is in a {b, c} case that becomes single once c (never chosen) is switched off
  set.seed(35); rows <- list()
  for (g in 1:40) for (t in 1:6) {
    cs <- if (g == 1 && t == 1) c("b", "c") else if (g == 1) c("a", "b") else c("a", "b", "c")
    ch <- if (g == 1 && t == 1) "b" else if (g == 1) "a" else sample(cs, 1)
    for (a in cs) rows[[length(rows) + 1]] <- data.frame(id = paste(g, t), alt = a, unit = g, x = rnorm(1), chosen = as.integer(a == ch))
  }
  dl <- do.call(rbind, rows); dl$alt <- factor(dl$alt, levels = c("a", "b", "c")); dl$unit <- factor(dl$unit)
  m <- suppressWarnings(suppressMessages(mnl(chosen ~ x, data = dl, id = "id", alt = "alt", fe = "unit")))
  expect_true(m$converged); expect_false(m$boundary)
  expect_false("b:unit1" %in% names(coef(m)))                          # unit 1 is left choosing a only: it drops out
  ## an alternative chosen only by units that always choose it identifies nothing
  set.seed(12); G <- 60; T <- 6; unit <- rep(seq_len(G), each = T); x1 <- rnorm(G * T)
  y <- sample(1:2, G * T, TRUE); y[unit > 50] <- 3L
  d2 <- data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1, unit = factor(unit))
  expect_error(suppressWarnings(suppressMessages(mnl(y ~ x1, data = d2, fe = "unit"))), "only by units that always choose them")
  ## alternatives chosen only by units unconnected to the base (they choose c and d, never a or b) are refused by name
  y4 <- sample(1:2, G * T, TRUE); grp <- unit > 40; y4[grp] <- sample(3:4, sum(grp), TRUE)
  d4 <- data.frame(y = factor(c("a", "b", "c", "d")[y4], levels = c("a", "b", "c", "d")), x1, unit = factor(unit))
  expect_error(suppressWarnings(suppressMessages(mnl(y ~ x1, data = d4, fe = "unit"))), "constant\\(s\\) of c, d are not identified relative to the base")
  ## a factor level present only in dropped (zero-weight) cases is refused with that cause
  d5 <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d5$f <- factor(ifelse(seq_len(300) <= 3, "rare", sample(c("p", "q"), 300, TRUE))); d5$w <- ifelse(d5$f == "rare", 0, 1)
  expect_error(suppressMessages(mnl(y ~ x1 + f, data = d5, weights = "w")), "zero on every case of the estimation sample")
})

test_that("fixed effects without alternative-specific constants keep every unit's dummies (the same model)", {
  d <- fe_panel(50, 8, 34)
  a1 <- suppressWarnings(suppressMessages(mnl(y ~ x1, data = d, fe = "unit")))
  a0 <- suppressWarnings(suppressMessages(mnl(y ~ x1 - 1, data = d, fe = "unit")))
  expect_equal(a0$loglik, a1$loglik, tolerance = 1e-8); expect_equal(a0$df, a1$df)
  expect_equal(coef(a0)[c("b:x1", "c:x1")], coef(a1)[c("b:x1", "c:x1")], tolerance = 1e-6)
})

test_that("fixed-effects bootstrap: units are resampled, a replicate that cannot be built counts as failed", {
  skip_on_cran()
  set.seed(32); G <- 60; T <- 6; unit <- rep(seq_len(G), each = T); x1 <- rnorm(G * T)
  y <- sample(1:2, G * T, TRUE); y[unit == 1] <- rep(c(2L, 3L), 3); y[unit %in% 2:4] <- 3L
  d <- data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1, unit = factor(unit))
  w <- character(0)
  m <- withCallingHandlers(suppressMessages(mnl(y ~ x1, data = d, fe = "unit", se = "bootstrap", nboot = 20)),
                           warning = function(e) { w <<- c(w, conditionMessage(e)); invokeRestart("muffleWarning") })
  expect_true(m$boot$cluster)                                             # units, not cases
  expect_gt(m$boot$n_failed, 0L); expect_true(any(grepl("bootstrap replicates failed.*fixed effects can be built", w)))
  expect_false(is.null(m$boot$fail_reasons))
  expect_true(is.finite(sqrt(vcov(m)["b:x1", "b:x1"])))
})

test_that("random intercepts: the prediction rule is linear in sigma beyond the Gauss-Hermite range and says when it is capped", {
  r <- polychoice:::.pc_re_pred_rule(30, 12L, 1L)
  expect_lte(length(r$logw), 1024L); expect_false(r$capped)
  eta <- c(-2, 0, 1.5); sig <- 30
  zz <- seq(-12, 12, length.out = 20001); w <- dnorm(zz); w <- w / sum(w)
  fine <- sapply(eta, function(e) sum(w * plogis(e + sig * zz)))
  rule <- sapply(eta, function(e) sum(exp(r$logw) * plogis(e + sig * r$Z[, 1])))
  expect_equal(rule, fine, tolerance = 1e-10)
  g <- polychoice:::.pc_re_pred_rule(1.5, 12L, 1L); expect_equal(length(g$logw), 32L)    # Gauss-Hermite on the 2^k ladder above 10 sigma^2
  c4 <- polychoice:::.pc_re_pred_rule(5, 12L, 3L); expect_true(c4$capped); expect_gt(c4$err, 1e-6)
  expect_lte(length(c4$logw), 5e4)
})

test_that("identification and scaling are scale-free: shifted and tiny-scale covariates", {
  set.seed(37); d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d$x3 <- d$x1 * 1e-12
  m1 <- mnl(y ~ x1 + x2, data = d); m3 <- mnl(y ~ x3 + x2, data = d)
  expect_true(m3$converged)
  expect_equal(unname(coef(m3)[c("b:x3", "c:x3")]) * 1e-12, unname(coef(m1)[c("b:x1", "c:x1")]), tolerance = 1e-6)
  skip_if_not_installed("AER")
  data("TravelMode", package = "AER"); TM <- TravelMode; TM$chosen <- as.integer(TM$choice == "yes"); TM$wait9 <- TM$wait + 1e9
  c1 <- mnl(chosen ~ wait + vcost | income, data = TM, id = "individual", alt = "mode")
  c2 <- mnl(chosen ~ wait9 + vcost | income, data = TM, id = "individual", alt = "mode")
  expect_equal(unname(coef(c2)["wait9"]), unname(coef(c1)["wait"]), tolerance = 1e-4)   # the shift by 1e9 costs digits in the centred covariate: 1e-6 relative on macOS arm64
})

test_that("comparison guards: split equation without intercept, fixed effects on one side, and the Vuong reference at the boundary", {
  set.seed(36); d <- rmnl(400, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m0 <- mnl(y ~ x1 + x2, data = d)
  expect_error(boot_lr_test(m0, suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | x1 - 1, data = d))), B = 10), "no intercept")
  dp <- fe_panel(40, 8, 41)
  expect_error(lr_test(mnl(y ~ x1, data = dp), suppressWarnings(suppressMessages(mnl(y ~ x1, data = dp, fe = "unit")))), "fixed effects")
  v <- vuong(suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | x1, data = d))), m0)
  expect_match(attr(v, "note"), "boot_lr_test")
  expect_output(print(v), "Note:")
})

test_that("Fok-Paap Q1: the joint covariance is positive semi-definite; the variance difference is available", {
  set.seed(1); d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
  m <- mnl(y ~ x1 + x2, data = d)
  fj <- fp_tests(m, tests = "cl"); fd <- fp_tests(m, tests = "cl", vcov = "difference")
  expect_true(fj$psd); expect_gte(fj$statistic, 0)
  expect_match(fj$test, "joint"); expect_match(fd$test, "difference")
  ## the joint covariance equals the outer product of the stacked influence functions
  D <- m$design; fn <- polychoice:::.pc_fns(D)
  V_ml <- polychoice:::.pc_strip_vcov(polychoice:::.pc_vcov(m$theta, fn, "analytic", NULL, m$scale, D$w))
  cl <- polychoice:::.pc_cl_fit(polychoice:::.pc_pairs(D, "all"), m$theta, D$n, D$K)
  Psi <- fn$scores(m$theta) %*% V_ml - cl$scores %*% solve(cl$H); dl <- m$theta - cl$theta
  expect_equal(fj$statistic, as.numeric(t(dl) %*% solve(crossprod(Psi), dl)), tolerance = 1e-8)
})

test_that("reporting: Small-Hsiao median p-value and separated splits, flagged bootstrap reasons, first-difference header, glance flags", {
  set.seed(1); d <- rmnl(1200, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
  m <- mnl(y ~ x1 + x2, data = d)
  sh <- small_hsiao(m, drop = "c", splits = 4, seed = 2)
  s <- attr(sh, "summary")
  expect_equal(s$median_p, pchisq(median(sh$statistic), sh$df[1], lower.tail = FALSE))
  expect_true(is.logical(sh$separated)); expect_true("share_separated" %in% names(s))
  bp <- boot_pvalue(m, hausman_mcfadden, B = 12, drop = "b")
  if (attr(bp, "n_flagged") > 0L) expect_equal(sum(attr(bp, "flag_reasons")), attr(bp, "n_flagged"))
  di <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), delta = c(0.5, 1), inflate = "a")
  mi <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | x1 + z1, data = di, inflate = "a")))
  expect_output(print(first_difference(mi, "x1", 0, 1, stage = "outcome")), "intervals at 95%; outcome stage only")
  g <- polychoice:::glance.choicemodel(mi); expect_true(all(c("n_modes", "converged", "boundary", "ill_conditioned") %in% names(g)))
  gf <- polychoice:::glance.choicemodel(suppressWarnings(suppressMessages(mnl(y ~ x1, data = fe_panel(30, 6, 43), fe = "unit"))))
  expect_equal(gf$bias_correct, "none"); expect_false(gf$jackknife_refused)
})
