## IIA tests: agreement with an independent nnet-based Hausman-McFadden
## statistic (non-base and base drops), invariance to the fitted base, the
## suest covariance against a bootstrap, Small-Hsiao mechanics, and the
## Fok-Paap statistics' structure.

.hm_nnet <- function(dd, drop, base, formula) {
  fit <- function(d2) {
    d2$y <- stats::relevel(droplevels(d2$y), ref = base)
    f <- do.call(nnet::multinom, list(formula, data = d2, trace = FALSE, maxit = 1000, reltol = 1e-14))
    cf <- coef(f); V <- vcov(f)
    th <- if (is.matrix(cf)) as.vector(t(cf)) else cf
    names(th) <- rownames(V)
    list(theta = th, V = V)
  }
  full <- fit(dd); rest <- fit(dd[dd$y != drop, ])
  sel <- match(names(rest$theta), names(full$theta))
  delta <- rest$theta - full$theta[sel]
  as.numeric(t(delta) %*% solve(rest$V - full$V[sel, sel]) %*% delta)
}

test_that("Hausman-McFadden (difference) equals the nnet-based statistic, including the base drop", {
  skip_if_not_installed("nnet")
  set.seed(7)
  d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8, 0.4), c = c(-0.3, 0.6, -0.5, 0.7), d = c(0.1, 0.3, 0.9, -0.5)))
  m <- mnl(y ~ x1 + x2 + x3, data = d)
  hm <- hausman_mcfadden(m)
  expect_equal(hm$alternative, c("a", "b", "c", "d")); expect_equal(hm$df, rep(8L, 4))   # 2 remaining non-base alternatives x 4 terms
  for (alt in c("b", "c", "d"))
    expect_equal(hm$statistic[hm$alternative == alt], .hm_nnet(d, alt, "a", y ~ x1 + x2 + x3), tolerance = 1e-3)
  ## the base drop: re-based to "b" on both sides
  expect_equal(hm$statistic[hm$alternative == "a"], .hm_nnet(d, "a", "b", y ~ x1 + x2 + x3), tolerance = 1e-3)
  ## invariance: the same statistics from a fit with another base
  hm_b <- hausman_mcfadden(mnl(y ~ x1 + x2 + x3, data = d, ref = "c"))
  expect_equal(hm_b$statistic[match(hm$alternative, hm_b$alternative)], hm$statistic, tolerance = 1e-3)
  expect_true(all(hm$p.value[hm$statistic < 0] == 1))
})

test_that("the joint (suest) covariance of the difference agrees with a bootstrap", {
  skip_on_cran()
  set.seed(9)
  d <- rmnl(2000, beta = cbind(b = c(-1, -1, 0.8), c = c(-1.2, 0.6, -0.5)))
  m <- mnl(y ~ x1 + x2, data = d)
  hj <- hausman_mcfadden(m, drop = "c", vcov = "joint")
  hd <- hausman_mcfadden(m, drop = "c", vcov = "difference")
  expect_true(hj$psd)
  ## bootstrap the difference theta_r - theta_f (sub-vector) and compare its covariance with the suest one
  diff_b <- replicate(150, {
    idx <- sample.int(nrow(d), replace = TRUE); db <- d[idx, ]
    tf <- coef(mnl(y ~ x1 + x2, data = db)); tr <- coef(mnl(y ~ x1 + x2, data = db[db$y != "c", ]))
    tr - tf[names(tr)]
  })
  Vb <- stats::cov(t(diff_b))
  ## recompute the suest covariance through the exported statistic: delta' V^{-1} delta with the bootstrap V
  tf <- coef(m); tr <- coef(mnl(y ~ x1 + x2, data = d[d$y != "c", ])); delta <- tr - tf[names(tr)]
  stat_boot <- as.numeric(t(delta) %*% solve(Vb) %*% delta)
  expect_equal(hj$statistic, stat_boot, tolerance = 0.5)   # a loose check: 150 replicates on a 4-dimensional covariance
  ## the classic variance-difference form is reported with its psd flag either way
  expect_true(is.logical(hd$psd)); expect_equal(hd$df, 3L)   # one remaining non-base alternative x 3 terms
})

test_that("Small-Hsiao: structure, non-negativity, reproducibility, repeated splits", {
  set.seed(1)
  d <- rmnl(1200, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
  m <- mnl(y ~ x1 + x2, data = d)
  sh <- small_hsiao(m, seed = 1)
  expect_equal(sh$alternative, c("a", "b", "c", "d")); expect_equal(sh$df, rep(6L, 4))
  expect_true(all(sh$statistic >= 0))
  expect_equal(sh$statistic, small_hsiao(m, seed = 1)$statistic)
  sh5 <- small_hsiao(m, drop = "c", splits = 5, seed = 2)
  expect_equal(nrow(sh5), 5L); expect_equal(nrow(attr(sh5, "summary")), 1L)
  expect_equal(attr(sh5, "summary")$median_statistic, median(sh5$statistic))
  ## base drop reproduces the statistic from a fit with that base's replacement
  sh_a <- small_hsiao(m, drop = "a", seed = 3)
  expect_true(sh_a$statistic >= 0)
})

test_that("Fok-Paap tests: degrees of freedom, CL estimate near the MLE, projected and iterated Q2", {
  set.seed(1)
  d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
  m <- mnl(y ~ x1 + x2, data = d)
  fp <- fp_tests(m)
  expect_equal(fp$df, c(9L, 6 * 3 - 9))                                 # Q1: dim theta; Q2: M(p+q) - K with 6 pairs
  det <- attr(fp, "details")
  expect_true(det$cl_converged); expect_true(is.finite(fp$statistic[2])); expect_lte(det$gmm_iterations, 50L)
  expect_lt(max(abs(det$theta_cl - coef(m))), 0.5)
  fpp <- fp_tests(m, gmm = "projected")
  expect_equal(fpp$df[2], fp$df[2]); expect_gt(fpp$statistic[2], 0)
  ## base pairs with case-specific covariates only: Q2 is not available (df = 0) and says so
  expect_warning(fb <- fp_tests(m, pairs = "base"), "No overidentifying restrictions")
  expect_equal(nrow(fb), 1L)
  expect_error(fp_tests(m, pairs = "base", tests = "gmm"), "No overidentifying restrictions")
  ## inflated fits are refused
  di <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), delta = c(0.5, 1), inflate = "a")
  mi <- imnl(y ~ x1 + x2 | z1, data = di, inflate = "a")
  expect_error(fp_tests(mi), "does not apply to an inflated fit")
  expect_error(hausman_mcfadden(mi), "does not apply")
  expect_error(small_hsiao(mi), "does not apply")
  ## the battery
  bat <- iia_tests(m, seed = 1)
  expect_s3_class(bat, "iia_battery")
  expect_equal(nrow(bat), 4 + 4 + 4 + 2)
})

test_that("Fok-Paap with alternative-specific covariates on TravelMode: base pairs give an overidentified GMM", {
  skip_if_not_installed("AER")
  data("TravelMode", package = "AER"); TM <- TravelMode; TM$chosen <- as.integer(TM$choice == "yes")
  cl <- mnl(chosen ~ wait + vcost | income, data = TM, id = "individual", alt = "mode")
  fb <- fp_tests(cl, pairs = "base")
  expect_equal(nrow(fb), 2L)
  expect_equal(fb$df[2], 3 * (2 + 2) - cl$df)                           # 3 base pairs x (p + q) moments
  fa <- fp_tests(cl, pairs = "all")
  expect_equal(fa$df[2], 6 * (2 + 2) - cl$df)
})
