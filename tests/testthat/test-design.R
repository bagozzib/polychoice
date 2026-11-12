## Design and engine identities: re-basing, frequency weights, robust and
## cluster covariances by hand, offsets, prediction on new data, fixed
## effects, and the refusals.

test_that("re-basing leaves the fit invariant and maps the coefficients exactly", {
  set.seed(11)
  d <- rmnl(1000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m <- mnl(y ~ x1 + x2, data = d); m2 <- mnl(y ~ x1 + x2, data = d, ref = "b")
  expect_equal(m2$loglik, m$loglik, tolerance = 1e-9)
  expect_equal(unname(fitted(m2)), unname(fitted(m)), tolerance = 1e-6)
  expect_equal(unname(coef(m2)[c("c:(Intercept)", "c:x1", "c:x2")]),
               unname(coef(m)[c("c:(Intercept)", "c:x1", "c:x2")] - coef(m)[c("b:(Intercept)", "b:x1", "b:x2")]), tolerance = 1e-5)
  expect_equal(unname(coef(m2)[c("a:(Intercept)", "a:x1", "a:x2")]), -unname(coef(m)[c("b:(Intercept)", "b:x1", "b:x2")]), tolerance = 1e-5)
  expect_equal(m2$ref, "b"); expect_equal(m$levels, c("a", "b", "c"))
})

test_that("frequency weights equal duplicated cases; offsets shift the utilities", {
  set.seed(2)
  d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  w <- sample(1:3, nrow(d), TRUE)
  mw <- mnl(y ~ x1 + x2, data = d, weights = w)
  md <- mnl(y ~ x1 + x2, data = d[rep(seq_len(nrow(d)), w), ])
  expect_equal(coef(mw), coef(md), tolerance = 1e-7)
  expect_equal(vcov(mw), vcov(md), tolerance = 1e-6)
  expect_equal(mw$loglik, md$loglik, tolerance = 1e-8)
  expect_equal(mw$nobs_weighted, sum(w)); expect_equal(mw$n, nrow(d))
  expect_equal(fit_stats(mw)$BIC, fit_stats(md)$BIC, tolerance = 1e-6)
  expect_equal(vuong(mw, mnl(y ~ x1, data = d, weights = w))$statistic, vuong(md, mnl(y ~ x1, data = d[rep(seq_len(nrow(d)), w), ]))$statistic, tolerance = 1e-6)
  ## an offset equal to a coefficient times a covariate removes that coefficient's need
  m <- mnl(y ~ x1 + x2, data = d)
  off <- cbind(0, coef(m)["b:x2"] * d$x2, coef(m)["c:x2"] * d$x2)
  mo <- mnl(y ~ x1, data = d, offset = off)
  expect_equal(mo$loglik, m$loglik, tolerance = 1e-6)
  expect_equal(unname(coef(mo)), unname(coef(m)[c("b:(Intercept)", "b:x1", "c:(Intercept)", "c:x1")]), tolerance = 1e-4)
})

test_that("robust and cluster covariances equal the hand-built sandwiches", {
  set.seed(3)
  d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d$g <- sample(1:30, nrow(d), TRUE)
  m <- mnl(y ~ x1 + x2, data = d)
  S <- .pc_fns(m$design)$scores(coef(m))
  expect_equal(unname(vcov(mnl(y ~ x1 + x2, data = d, se = "robust"))), unname(vcov(m) %*% crossprod(S) %*% vcov(m)), tolerance = 1e-8)
  mc <- mnl(y ~ x1 + x2, data = d, cluster = "g")
  expect_equal(mc$se, "cluster")
  expect_equal(unname(vcov(mc)), unname(vcov(m) %*% crossprod(rowsum(S, d$g)) %*% vcov(m)), tolerance = 1e-8)
  expect_warning(mnl(y ~ x1 + x2, data = d, cluster = "g", se = "analytic"), "ignored")
  expect_error(mnl(y ~ x1 + x2, data = d, se = "cluster"), "requires 'cluster'")
})

test_that("predict() on new data matches fitted values and handles factors and types", {
  set.seed(4)
  d <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), delta = c(0.5, 1), inflate = "a")
  d$g <- factor(sample(c("u", "v", "w"), nrow(d), TRUE))
  m <- imnl(y ~ x1 + g | z1, data = d, inflate = "a")
  expect_equal(predict(m, newdata = d[1:9, ]), fitted(m)[1:9, ], tolerance = 1e-10)
  expect_equal(unname(predict(m, newdata = d[1:9, ], type = "regime")), unname(m$phi[1:9]), tolerance = 1e-10)
  expect_equal(predict(m, type = "prob_outcome"), m$prob_outcome, tolerance = 1e-12)
  expect_equal(rowSums(predict(m, type = "prob")), rep(1, nrow(d)), ignore_attr = TRUE)
  expect_s3_class(predict(m, type = "class"), "factor")
  expect_equal(dim(predict(m, type = "link")), c(nrow(d), 3L))
  nd <- d[1:3, ]; nd$g <- factor("zzz")
  expect_error(predict(m, newdata = nd))
  expect_error(predict(mnl(y ~ x1, data = d), type = "regime"), "inflated models only")
})

test_that("the long layout is parsed, checked, and predicted correctly", {
  set.seed(5)
  n <- 300; J <- 3
  long <- data.frame(id = rep(seq_len(n), each = J), alt = factor(rep(c("car", "bus", "train"), n), levels = c("car", "bus", "train")),
                     cost = runif(n * J), time = runif(n * J), inc = rep(rnorm(n), each = J))
  V <- -1.5 * long$cost + 0.8 * long$time + 0.5 * long$inc * (long$alt == "bus") + 0.2 * (long$alt == "train")
  Vm <- matrix(V, n, J, byrow = TRUE); P <- exp(Vm) / rowSums(exp(Vm))
  ch <- apply(P, 1, function(p) sample(J, 1, prob = p)); long$chosen <- as.integer(rep(ch, each = J) == rep(1:J, n))
  m <- mnl(chosen ~ cost + time | inc, data = long, id = "id", alt = "alt")
  expect_equal(m$layout, "long"); expect_equal(names(coef(m))[1:4], c("bus:(Intercept)", "bus:inc", "train:(Intercept)", "train:inc"))
  expect_equal(names(coef(m))[5:6], c("cost", "time"))
  expect_equal(unname(predict(m, newdata = long[1:6, ])), unname(fitted(m)[1:2, ]), tolerance = 1e-10)
  ## a case-constant variable in the alternative-specific part is not identified: refused by name
  expect_error(mnl(chosen ~ cost + inc, data = long, id = "id", alt = "alt"), "do not vary within cases")
  ## an alternative-specific variable in the case-specific part: refused
  expect_error(mnl(chosen ~ cost | time, data = long, id = "id", alt = "alt"), "varies within a case")
  ## duplicated (id, alt) and a case with two choices
  bad <- rbind(long, long[1, ]); expect_error(mnl(chosen ~ cost | inc, data = bad, id = "id", alt = "alt"), "more than one row")
  bad2 <- long; bad2$chosen[1:2] <- 1L; expect_error(mnl(chosen ~ cost | inc, data = bad2, id = "id", alt = "alt"), "exactly one chosen")
  expect_error(mnl(chosen ~ cost | inc, data = long, id = "id"), "both 'id'")
  ## weights must be constant within a case
  long$w <- runif(nrow(long)); expect_error(mnl(chosen ~ cost | inc, data = long, id = "id", alt = "alt", weights = "w"), "constant within a case")
})

test_that("refusals are informative: rank deficiency, character response, never-chosen alternative", {
  set.seed(6)
  d <- rmnl(400, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  expect_error(mnl(y ~ x1 + I(2 * x1), data = d), "not identified")
  expect_error(mnl(as.character(y) ~ x1, data = d), "character vector")
  expect_error(mnl(y ~ x1 + I(x1 * 0 + 1), data = d), "Constant covariate")
  d2 <- d; d2$y <- factor(d2$y, levels = c("a", "b", "c", "zzz")); expect_silent(mnl(y ~ x1, data = d2))   # unused level dropped
  expect_error(imnl(y ~ x1 | x1, data = d, inflate = "q"), "'inflate' must be one")
  expect_error(mnl(y ~ x1, data = d, ref = "q"), "'ref' must be one")
  expect_error(mnl(y ~ x1 | x1, data = d), "split equation")
  expect_error(imnl(y ~ x1, data = d), "imnl\\(\\) with one row per case takes")
  dd <- d; dd$x1[1:20] <- NA
  m <- mnl(y ~ x1 + x2, data = dd); expect_equal(m$n, 380L); expect_equal(m$n_dropped_rows, 20L)
})

test_that("unit fixed effects: never-chosen pairs become unavailable and predictions follow", {
  set.seed(8)
  n <- 1600
  dfe <- data.frame(x1 = rnorm(n), unit = factor(rep(1:40, each = 40)))
  ue <- rnorm(40, 0, 1.5)[as.integer(dfe$unit)]
  V <- cbind(0, 0.5 - dfe$x1 + ue, -0.3 + 0.6 * dfe$x1); P <- exp(V) / rowSums(exp(V))
  dfe$y <- factor(c("a", "b", "c")[apply(P, 1, function(p) sample(3, 1, prob = p))], levels = c("a", "b", "c"))
  mfe <- suppressWarnings(mnl(y ~ x1, data = dfe, fe = "unit"))               # a unit with one choice of some alternative is weakly identified
  expect_true(mfe$converged)
  expect_equal(sum(mfe$blocks == "fixed effect"), 2L * 39L)
  expect_equal(unname(coef(mfe)[c("b:x1", "c:x1")]), c(-1, 0.6), tolerance = 0.2)
  expect_equal(unname(predict(mfe, newdata = dfe[1:5, ])), unname(fitted(mfe)[1:5, ]), tolerance = 1e-10)
  ## unit 3 never chooses c: c unavailable to it, its fitted P(c) is zero, no dummy for (c, 3)
  dfe2 <- dfe; dfe2$y[dfe2$unit == "3" & dfe2$y == "c"] <- "a"
  expect_message(mfe2 <- suppressWarnings(mnl(y ~ x1, data = dfe2, fe = "unit")), "never chosen")
  expect_equal(max(fitted(mfe2)[dfe2$unit == "3", "c"]), 0)
  expect_false("c:unit3" %in% names(coef(mfe2)))
  expect_equal(max(predict(mfe2, newdata = dfe2[dfe2$unit == "3", ][1:3, ])[, "c"]), 0)
  expect_error(predict(mfe2, newdata = data.frame(x1 = 0, unit = factor("999"))), "not in the fit")
  ## a short panel warns about incidental parameters
  dshort <- dfe[dfe$unit %in% levels(dfe$unit)[1:20], ]; dshort$unit <- factor(rep(1:100, each = 8))
  w <- testthat::capture_warnings(suppressMessages(mnl(y ~ x1, data = dshort, fe = "unit")))
  expect_true(any(grepl("incidental-parameters", w)))
})

test_that("round 6 (C-3): the overlap print counts comparisons with identical design rows; the long layout's `| 0` message", {
  ## comparisons +1 (case 1), -1 (case 2), 0 (case 3: identical rows), +1 (case 4): overlap with one zero comparison
  d <- data.frame(id = rep(1:4, each = 2), alt = factor(rep(c("p", "q"), 4)), chosen = c(1, 0, 0, 1, 1, 0, 0, 1), x = c(1, 0, 1, 0, 1, 1, 0, 1))
  m <- mnl(chosen ~ 0 + x, data = d, id = "id", alt = "alt")
  s <- separation(m)
  expect_false(s$separated)
  expect_identical(c(s$n_comparisons, s$n_zero), c(4L, 1L))
  expect_output(print(s), "1 of them between alternatives with identical design rows")
  ## `| 0` does not drop the constants in the long layout (unlike mlogit), and says so
  set.seed(3)
  dl <- data.frame(id = rep(1:60, each = 3), alt = factor(rep(c("a", "b", "c"), 60)), w = stats::rnorm(180))
  u <- dl$w - log(-log(stats::runif(180))); dl$chosen <- as.integer(stats::ave(u, dl$id, FUN = function(v) v == max(v)))
  expect_message(m0 <- mnl(chosen ~ w | 0, data = dl, id = "id", alt = "alt"), "does not drop them")
  expect_true(all(c("b:(Intercept)", "c:(Intercept)") %in% names(coef(m0))))
})
