## DESIGN.md section 1, Phase 4e: the likelihood-ratio test of the logit
## against the inflated multinomial logit with a constant captive share, by
## captivity_test(type = "lr") and lr_test(), with the 50:50 mixture of
## chi-square(0) and chi-square(1) as reference

test_that("the mixture p-value", {
  chibar <- polychoice:::.pc_chibar_p
  expect_equal(chibar(0, 0L), 1)
  expect_equal(chibar(stats::qchisq(0.9, 1), 0L), 0.05, tolerance = 1e-12)
  expect_equal(chibar(3, 2L), 0.5 * stats::pchisq(3, 2, lower.tail = FALSE) + 0.5 * stats::pchisq(3, 3, lower.tail = FALSE))
  expect_true(is.na(chibar(NA_real_, 0L)))
})

test_that("captivity_test(type = \"lr\") is the likelihood ratio of imnl(y ~ x | 1) against the logit, on both layouts and with weights", {
  set.seed(41)
  B <- cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5))
  d <- rmnl(800, beta = B, delta = c(qlogis(0.8), 0), inflate = "a")
  d$w <- sample(1:3, nrow(d), replace = TRUE)
  m <- mnl(y ~ x1 + x2, data = d, weights = "w")
  expect_no_warning(ct <- captivity_test(m, alternatives = c("a", "b"), type = "lr"))
  expect_equal(ct$test, c("captive: a", "captive: b"))
  ma <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | 1, data = d, inflate = "a", weights = "w")))
  expect_equal(ct$statistic[1L], 2 * (ma$loglik - m$loglik), tolerance = 1e-6)
  expect_equal(ct$share[1L], unname(1 - stats::plogis(coef(ma)[["infl_(Intercept)"]])), tolerance = 1e-4)
  expect_equal(ct$p.value[1L], 0.5 * stats::pchisq(ct$statistic[1L], 1, lower.tail = FALSE))
  expect_lt(ct$p.value[1L], 0.001)                                          # a fifth of the cases captive to a
  expect_true(all(ct$p.value[ct$statistic == 0] == 1))
  expect_true(all(attr(ct, "fits")$status %in% c("logit", "interior")))
  expect_true(all(attr(ct, "fits")$top_share > 0 & attr(ct, "fits")$top_share < 1))
  expect_identical(nrow(attr(ct, "fits")), 2L)
  ## lr_test() of the same pair: the same statistic and p-value
  lt <- lr_test(m, ma)
  expect_equal(lt$statistic, ct$statistic[1L], tolerance = 1e-6)
  expect_equal(lt$p.value, ct$p.value[1L], tolerance = 1e-6)
  expect_identical(attr(lt, "mixture"), 0L)
  expect_output(print(lt), "50:50 mixture")
  ## the long layout: the same model, the same statistic
  skip_on_cran()
  dl <- do.call(rbind, lapply(seq_len(nrow(d)), function(i) data.frame(id = i, alt = factor(c("a", "b", "c")), chosen = as.integer(c("a", "b", "c") == d$y[i]),
                                                                        x1 = d$x1[i], x2 = d$x2[i], w = d$w[i])))
  ml <- mnl(chosen ~ 1 | x1 + x2, data = dl, id = "id", alt = "alt", weights = "w")   # the constants follow the first part
  expect_equal(ml$loglik, m$loglik, tolerance = 1e-8)
  expect_no_warning(cl <- captivity_test(ml, alternatives = "a", type = "lr"))
  expect_equal(cl$statistic, ct$statistic[1L], tolerance = 1e-6)
})

test_that("lr_test() refuses a split with covariates and several inflated alternatives; the LR form's refusals", {
  skip_on_cran()
  set.seed(42)
  B <- cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5))
  d <- rmnl(600, beta = B, delta = c(qlogis(0.8), 0.5), inflate = "a")
  m <- mnl(y ~ x1 + x2, data = d)
  mz <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | z1, data = d, inflate = "a")))
  expect_error(lr_test(m, mz), "intercept-only split")
  m2 <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | 1, data = d, inflate = c("a", "c"))))
  expect_error(lr_test(m, m2), "boot_lr_test")
  ## a further outcome covariate (two coefficients inside): the mixture of chi-square(2) and chi-square(3)
  d$x3 <- rnorm(nrow(d))
  m3 <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 + x3 | 1, data = d, inflate = "a")))
  lt <- lr_test(m, m3)
  expect_identical(lt$df, 3L); expect_identical(attr(lt, "mixture"), 2L)
  expect_equal(lt$p.value, polychoice:::.pc_chibar_p(lt$statistic, 2L))
  ## two alternatives; penalized fits
  d2 <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8)))
  expect_error(captivity_test(mnl(y ~ x1 + x2, data = d2), type = "lr"), "three alternatives")
  expect_error(captivity_test(mnl(y ~ x1 + x2, data = d, firth = TRUE), type = "lr"), "penalized")
})

test_that("captivity_test(type = \"lr\"): the logit at the boundary, an unidentified share, a never-chosen alternative (round 6)", {
  skip_on_cran()
  ## a saturated outcome design (one factor): a constant share of an
  ## alternative is absorbed by its constants, so the ratio is 0, p 1, and the
  ## share is not identified (D-new-2); no warning (I-1 / A-new-4)
  set.seed(43)
  d <- data.frame(y = factor(sample(c("a", "b", "c"), 900, TRUE, prob = c(0.3, 0.3, 0.4))), g = factor(sample(1:3, 900, TRUE)))
  m <- mnl(y ~ g, data = d)
  expect_no_warning(ct <- captivity_test(m, alternatives = "a", type = "lr"))
  expect_identical(ct$statistic, 0); expect_identical(ct$p.value, 1)
  expect_true(is.na(ct$share)); expect_match(attr(ct, "note"), "not identified")
  expect_identical(attr(ct, "fits")$status, "logit")
  ## a never-chosen alternative (on offer in the long layout, never chosen):
  ## NA with the reason
  dl <- data.frame(id = rep(1:300, each = 3), alt = factor(rep(c("a", "b", "c"), 300)), x = rep(rnorm(300), each = 3))
  pick <- ifelse(stats::runif(300) < stats::plogis(dl$x[dl$alt == "a"]), "a", "b")
  dl$chosen <- as.integer(as.character(dl$alt) == rep(pick, each = 3))
  m2 <- suppressWarnings(suppressMessages(mnl(chosen ~ 1 | x, data = dl, id = "id", alt = "alt")))
  c2 <- suppressWarnings(captivity_test(m2, alternatives = "a", type = "lr"))
  expect_true(is.na(c2$statistic)); expect_match(attr(c2, "note"), "never chosen")
  ## boot_lr_test() on the constant-split pair reports the mixture p-value too (B-N6 / C-5)
  set.seed(44)
  ds <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  ms <- mnl(y ~ x1 + x2, data = ds)
  mi <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | 1, data = ds, inflate = "a")))
  bt <- suppressWarnings(suppressMessages(boot_lr_test(ms, mi, B = 10)))
  expect_false(is.na(bt$p.value[1L]))
  expect_equal(bt$p.value[1L], polychoice:::.pc_chibar_p(bt$statistic[1L], 0L))
  expect_match(attr(bt, "method"), "chi-square mixture")
})

test_that("lr_test() with further restrictions: a collapsed share keeps the ratio of the further coefficients (round 7, A-1 / I-N1)", {
  skip_on_cran()
  ## no captivity, a strong x3: the constant-share fit with x3 can collapse to
  ## the logit with x3, and its ratio against the logit without x3 is x3's
  set.seed(101)
  d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8, 0.6), c = c(-0.3, 0.6, -0.5, -0.6)))
  m0 <- mnl(y ~ x1 + x2, data = d); m3 <- mnl(y ~ x1 + x2 + x3, data = d)
  m1 <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 + x3 | 1, data = d, inflate = "a")))
  lt <- lr_test(m0, m1)
  expect_identical(lt$df, 3L)
  expect_equal(lt$statistic, 2 * (m1$loglik - m0$loglik), tolerance = 1e-8)
  expect_gte(lt$statistic, 2 * (m3$loglik - m0$loglik) - 1e-6)                 # at least x3's ratio
  expect_gt(lt$statistic, 20); expect_lt(lt$p.value, 1e-4)
  expect_equal(lt$p.value, polychoice:::.pc_chibar_p(lt$statistic, 2L))
})

test_that("a degenerate constant-share fit: NA p-values in captivity_test(), lr_test() and boot_lr_test() (round 7, D-4 / B-2)", {
  skip_on_cran()
  ## every case with x > 1 chooses a, the rest choose at random: the logit's
  ## data overlap, while the constant-share model sends the outcome slopes to
  ## infinity and lets the share absorb the a-choosers below 1 (the pres2004
  ## mechanism; probe r7_captivity_lr_degenerate_construct)
  set.seed(1); n <- 400
  x <- stats::rnorm(n); y <- ifelse(x > 1, "a", sample(c("a", "b", "c"), n, TRUE))
  d <- data.frame(y = factor(y, levels = c("a", "b", "c")), x = x)
  m <- mnl(y ~ x, data = d)
  expect_false(separation(m)$separated)
  ct <- captivity_test(m, alternatives = "a", type = "lr")
  expect_identical(attr(ct, "fits")$status, "degenerate")
  expect_gt(ct$statistic, 10); expect_true(is.na(ct$p.value)); expect_match(attr(ct, "note"), "degenerate")
  mi <- suppressWarnings(suppressMessages(imnl(y ~ x | 1, data = d, inflate = "a")))
  expect_warning(lt <- lr_test(m, mi), "degenerate")
  expect_true(is.na(lt$p.value)); expect_equal(lt$statistic, ct$statistic, tolerance = 1e-6)
  set.seed(2)
  bt <- suppressWarnings(suppressMessages(boot_lr_test(m, mi, B = 10)))
  expect_true(is.na(bt$p.value[1L])); expect_match(attr(bt, "method"), "degenerate")
})

test_that("boot_pvalue() refuses tests whose null is a restriction of the fit (round 7, B-1)", {
  set.seed(7)
  d <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m <- mnl(y ~ x1 + x2, data = d)
  expect_error(boot_pvalue(m, term_test, B = 10), "restriction of that model")
  expect_error(boot_pvalue(m, combine_test, B = 10), "restriction of that model")
})
