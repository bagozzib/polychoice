## The further specification tests: combine (Wald = LR to first order on
## poolable alternatives; detection when they differ), per-term tests (Wald vs
## LR agreement, LR equals lr_test), Hosmer-Lemeshow-type fit, the nested-logit
## score test (size on MNL data, power on nested data, direction of the
## scores), and the parametric bootstrap wrapper.

test_that("combine_test(): poolable alternatives are not rejected, distinct ones are; Wald and LR agree", {
  set.seed(1)
  d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(0.3, -1, 0.8)))     # b and c share slopes
  m <- mnl(y ~ x1 + x2, data = d)
  cw <- combine_test(m); cl <- combine_test(m, type = "lr")
  expect_equal(cw$pair, c("a & b", "a & c", "b & c")); expect_equal(cw$df, rep(2L, 3))
  expect_gt(cw$p.value[cw$pair == "b & c"], 0.05); expect_gt(cl$p.value[cl$pair == "b & c"], 0.05)
  expect_lt(cw$p.value[cw$pair == "a & b"], 1e-6); expect_lt(cl$p.value[cl$pair == "a & b"], 1e-6)
  expect_equal(cw$statistic[3], cl$statistic[3], tolerance = 0.3)
  ## the LR form by hand: pooled fit plus the binomial split term
  dp <- d; dp$y <- factor(ifelse(dp$y == "c", "b", as.character(dp$y)), levels = c("a", "b"))
  mp <- mnl(y ~ x1 + x2, data = dp)
  nb <- sum(d$y == "b"); nc <- sum(d$y == "c")
  lr_hand <- 2 * (m$loglik - (mp$loglik + nb * log(nb / (nb + nc)) + nc * log(nc / (nb + nc))))
  expect_equal(cl$statistic[3], lr_hand, tolerance = 1e-6)
  expect_error(combine_test(m, c("a", "zzz")), "two different response levels")
})

test_that("term_test(): Wald and LR forms agree and LR equals lr_test()", {
  set.seed(2)
  d <- rmnl(1200, beta = cbind(b = c(0.5, -1, 0.0), c = c(-0.3, 0.6, 0.0)))     # x2 has no effect
  d$g <- factor(sample(c("u", "v", "w"), nrow(d), TRUE))
  m <- mnl(y ~ x1 + x2 + g, data = d)
  tw <- term_test(m); tl <- term_test(m, type = "lr")
  expect_equal(tw$term, c("x1", "x2", "g")); expect_equal(tw$df, c(2L, 2L, 4L))
  expect_lt(tw$p.value[1], 1e-8); expect_gt(tw$p.value[2], 0.01)
  ## Wald and LR agree for moderate effects; for the strong x1 term the Wald form is the smaller one (Hauck-Donner), both decisive
  expect_equal(tw$statistic[2:3], tl$statistic[2:3], tolerance = 0.25)
  expect_lt(tw$p.value[1], 1e-12); expect_lt(tl$p.value[1], 1e-12); expect_lt(tw$statistic[1], tl$statistic[1])
  expect_equal(tl$statistic[2], lr_test(mnl(y ~ x1 + g, data = d), m)$statistic, tolerance = 1e-6)
})

test_that("hl_test(): calibrated model passes, observed and expected totals agree", {
  set.seed(3)
  d <- rmnl(2000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  h <- hl_test(mnl(y ~ x1 + x2, data = d))
  expect_equal(h$df, 16L); expect_gt(h$p.value, 0.01)
  expect_equal(colSums(h$observed), colSums(h$expected), tolerance = 1e-8)
  expect_equal(sum(h$observed), 2000)
  expect_error(hl_test(mnl(y ~ x1 + x2, data = d), groups = 2), "at least 3")
})

test_that("nested_lm_test(): near-nominal on MNL data, rejects nested data with the right direction", {
  set.seed(4)
  nests <- list(n1 = c("a", "b"), n2 = c("c", "d"))
  d0 <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)))
  t0 <- nested_lm_test(mnl(y ~ x1 + x2, data = d0), nests); expect_equal(t0$df, 2L); expect_gt(t0$p.value, 0.001)
  t0e <- nested_lm_test(mnl(y ~ x1 + x2, data = d0), nests, equal_lambda = TRUE); expect_equal(t0e$df, 1L)
  t0o <- nested_lm_test(mnl(y ~ x1 + x2, data = d0), nests, information = "hessian"); expect_true(is.finite(t0o$statistic))
  ## lambda is weakly identified with case-specific covariates only (ledger): a strong departure from IIA is needed for power
  d1 <- rnlogit(3000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)), nests = nests, lambda = c(0.3, 0.3))
  m1 <- mnl(y ~ x1 + x2, data = d1)
  t1 <- nested_lm_test(m1, nests)
  expect_lt(t1$p.value, 1e-4)
  expect_equal(length(attr(t1, "scores")), 2L)
  t1h <- nested_lm_test(m1, nests, information = "hessian"); expect_equal(t1h$statistic, t1$statistic, tolerance = 0.5)
  ## the LM statistic is close to the LR statistic of the nested-logit fit
  n1 <- nlogit(y ~ x1 + x2, data = d1, nests = nests)
  expect_true(all(coef(n1)[c("lambda:n1", "lambda:n2")] < 1))                   # lambda pulled below 1 in both nests
  expect_equal(t1$statistic, lr_test(m1, n1)$statistic, tolerance = 0.5)
})

test_that("boot_pvalue() calibrates a test and keeps the observed table", {
  skip_on_cran()
  set.seed(5)
  d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
  m <- mnl(y ~ x1 + x2, data = d)
  b <- boot_pvalue(m, hausman_mcfadden, B = 30, drop = "b", vcov = "joint")
  expect_true("p.boot" %in% names(b)); expect_equal(nrow(b), 1L)
  expect_true(b$p.boot >= 0 && b$p.boot <= 1); expect_equal(attr(b, "B"), 30L)
  expect_equal(dim(attr(b, "boot_statistics")), c(1L, 30L))
  b2 <- boot_pvalue(m, nested_lm_test, B = 20, nests = list(n1 = c("a", "b"), n2 = c("c", "d")))
  expect_equal(nrow(b2), 1L); expect_true(is.finite(b2$p.boot))
})
