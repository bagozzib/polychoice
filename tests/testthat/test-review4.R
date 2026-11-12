## Regression tests for the defects found by the simulated-referee panel
## (round 4): the random-intercept long layout with dropped cases, the vuong()
## alignment guard, identification checks of subsample designs, the IIA
## battery's columns, the fixed-effects bootstrap guards, the fit-statistic
## null under fixed effects, the comparison messages, and numerical paths.

long_re <- function(G = 60, T = 4, seed = 1) {
  set.seed(seed); rows <- list(); u <- matrix(rnorm(G * 2, sd = 0.8), G)
  for (g in seq_len(G)) for (t in seq_len(T)) {
    x <- rnorm(1); w <- rnorm(3)
    v <- c(0, 0.3 + 0.8 * x + u[g, 1], -0.2 - 0.6 * x + u[g, 2]) + w
    ch <- sample.int(3, 1, prob = exp(v) / sum(exp(v)))
    rows[[length(rows) + 1L]] <- data.frame(id = paste(g, t), alt = c("a", "b", "c"), unit = g, x1 = x, w1 = w, chosen = as.integer(seq_len(3) == ch))
  }
  d <- do.call(rbind, rows); d$alt <- factor(d$alt, levels = c("a", "b", "c")); d$unit <- factor(d$unit); d
}

test_that("remnl(), long layout: a captive case in the first or the last unit is dropped and the units stay aligned", {
  skip_on_cran()
  d <- long_re()
  for (g in c("1", "60")) {
    cid <- d$id[d$unit == g][1]
    dd <- d[!(d$id == cid & d$chosen == 0L), ]              # that case keeps only its chosen row
    m <- suppressWarnings(suppressMessages(remnl(chosen ~ w1 | x1, data = dd, id = "id", alt = "alt", unit = "unit")))
    expect_equal(m$n, 60L * 4L - 1L)
    expect_equal(length(m$design$unit), m$n)
    expect_equal(as.character(m$design$unit), as.character(dd$unit[match(m$case_id, dd$id)]))
    expect_true(m$converged)
  }
})

test_that("vuong(): different case sets are refused; fits that are practically the same model give NA", {
  skip_on_cran()
  set.seed(31); G <- 150; T <- 12; unit <- rep(seq_len(G), each = T); n <- G * T
  x1 <- rnorm(n)
  V <- cbind(0, 1 + 0.8 * x1 + rnorm(G, sd = 3)[unit], 0.5 - 0.6 * x1 + rnorm(G, sd = 3)[unit]); P <- exp(V) / rowSums(exp(V))
  d <- data.frame(y = factor(c("a", "b", "c")[apply(P, 1, function(p) sample.int(3, 1, prob = p))], levels = c("a", "b", "c")),
                  x1, z1 = rnorm(n), unit = factor(unit))
  m0 <- suppressWarnings(suppressMessages(mnl(y ~ x1, data = d, fe = "unit")))
  d$x2 <- rnorm(n); d$x2[1:80] <- NA
  expect_error(vuong(m0, suppressWarnings(suppressMessages(mnl(y ~ x1 + x2, data = d, fe = "unit")))), "same cases")
  expect_error(imnl(y ~ x1 | z1, data = d, fe = "unit"), "does not take unit fixed effects")
  ## an inflated fit at phi = 1 under the MNL null (replication 8 of verification/vuong_size_mc.R):
  ## the per-case differences are round-off, and the corrected statistics used to run to minus millions
  set.seed(20261014 + 8); dv <- rmnl(1000, beta = cbind(b = c(0.2, -1, 0.5), c = c(-0.5, 0.8, -0.3))); dv$z1 <- rnorm(1000)
  mi <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | z1, data = dv, inflate = "a")))
  v <- vuong(mi, mnl(y ~ x1 + x2, data = dv))
  expect_true(all(is.na(v$statistic))); expect_true(all(is.na(v$p_two_sided)))
  expect_match(attr(v, "note"), "cannot be told apart")
})

test_that("subsample designs are identification-checked: period dummies refuse the jackknife by name", {
  skip_on_cran()
  set.seed(41); G <- 80; T <- 6; unit <- rep(seq_len(G), each = T); n <- G * T; t <- rep(seq_len(T), G)
  x1 <- rnorm(n)
  V <- cbind(0, 0.3 + 0.8 * x1 + rnorm(G)[unit], -0.2 - 0.6 * x1 + rnorm(G)[unit]); P <- exp(V) / rowSums(exp(V))
  d <- data.frame(y = factor(c("a", "b", "c")[apply(P, 1, function(p) sample.int(3, 1, prob = p))], levels = c("a", "b", "c")),
                  x1, period = factor(t), unit = factor(unit))
  m <- suppressMessages(withCallingHandlers(mnl(y ~ x1 + period, data = d, fe = "unit", bias_correct = "jackknife"),
                                            warning = function(w) invokeRestart("muffleWarning")))
  expect_identical(m$bias_correct, "none")
  expect_match(m$jackknife_refused, "zero on every case|not identified")
  ## a fixed-effects bootstrap whose cluster cuts across units is refused
  d$cl <- factor(rep(1:6, length.out = n))
  expect_error(suppressWarnings(suppressMessages(mnl(y ~ x1, data = d, fe = "unit", se = "bootstrap", cluster = "cl", nboot = 4))),
               "does not nest the units")
  ## the subsample check passes a sound inflated design (no false zero split column)
  di <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), delta = c(0.5, 1), inflate = "a")
  mi <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | z1, data = di, inflate = "a")))
  expect_true(is.list(polychoice:::.pc_subsample_check(polychoice:::.pc_subset_design(mi$design, seq_len(400)))))
})

test_that("IIA battery: psd only on covariance-difference rows, Small-Hsiao separation carried; per-row separation", {
  set.seed(1); d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
  b <- iia_tests(mnl(y ~ x1 + x2, data = d), seed = 1, splits = 3)
  expect_true(all(is.na(b$psd[grepl("Small-Hsiao|Q2", b$test)])))
  expect_false(anyNA(b$psd[grepl("Hausman-McFadden|Q1", b$test)]))
  expect_true(all(is.finite(b$share_separated[grepl("Small-Hsiao", b$test)])))
  expect_true(all(is.na(b$share_separated[!grepl("Small-Hsiao", b$test)])))
  skip_on_cran()
  data(pa2006, package = "polychoice")
  m <- mnl(vote ~ canvass + phone + democrat + republican + age + female + voted2004, data = pa2006)
  sh <- small_hsiao(m, splits = 1, seed = 1)          # split 1: the nonresponse equation separates in the first half
  expect_true(all(sh$separated[sh$alternative != "nonresponse"]))
  expect_false(sh$separated[sh$alternative == "nonresponse"])
})

test_that("fit statistics under fixed effects: one null for the same model; lr_test() warns for short panels", {
  set.seed(34); G <- 50; T <- 8; unit <- rep(seq_len(G), each = T); x1 <- rnorm(G * T); x2 <- rnorm(G * T)
  V <- cbind(0, 0.3 + 0.8 * x1 + rnorm(G)[unit], -0.2 - 0.6 * x1 + rnorm(G)[unit]); P <- exp(V) / rowSums(exp(V))
  d <- data.frame(y = factor(c("a", "b", "c")[apply(P, 1, function(p) sample.int(3, 1, prob = p))], levels = c("a", "b", "c")), x1, x2, unit = factor(unit))
  a1 <- suppressWarnings(suppressMessages(mnl(y ~ x1, data = d, fe = "unit")))
  a0 <- suppressWarnings(suppressMessages(mnl(y ~ x1 - 1, data = d, fe = "unit")))
  expect_equal(fit_stats(a0)$loglik_null, fit_stats(a1)$loglik_null, tolerance = 1e-8)
  a2 <- suppressWarnings(suppressMessages(mnl(y ~ x1 + x2, data = d, fe = "unit")))
  expect_message(lr_test(a1, a2), "incidental-parameters bias")
})

test_that("numerics: nested LM statistic invariant to covariate units; probit floor gives a finite gradient; remnl J = 3 rule", {
  set.seed(2)
  d <- rnlogit(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)),
               nests = list(n1 = c("a", "b"), n2 = c("c", "d")), lambda = c(0.6, 0.6))
  nests <- list(n1 = c("a", "b"), n2 = c("c", "d"))
  t1 <- nested_lm_test(mnl(y ~ x1 + x2, data = d), nests, information = "hessian")
  d$x1k <- d$x1 * 1e3
  t2 <- nested_lm_test(mnl(y ~ x1k + x2, data = d), nests, information = "hessian")
  expect_equal(t2$statistic, t1$statistic, tolerance = 1e-5)
  ## the probit floor: a case pushed far into a tail keeps a finite gradient
  dp <- rmnp(300, beta = cbind(b = c(0.3, -0.6, 0.5), c = c(-0.2, 0.4, -0.3)), sigma = matrix(c(1, 0.6, 0.6, 1.5), 2))
  mp <- suppressWarnings(mnprobit(y ~ x1 + x2, data = dp))
  fn <- polychoice:::.pc_fns(mp$design); th <- mp$theta; th[2] <- 200
  expect_true(all(is.finite(fn$grad(th))))
  ## remnl at J = 3: the population-averaged probabilities against a fine product trapezoid
  skip_on_cran()
  set.seed(5); G <- 80; T <- 3; unit <- rep(seq_len(G), each = T); x <- rnorm(G * T)
  u <- matrix(rnorm(G * 2, sd = 2.5), G)
  V <- cbind(0, 0.3 + 0.8 * x + u[unit, 1], -0.2 - 0.6 * x + u[unit, 2]); P <- exp(V) / rowSums(exp(V))
  dr <- data.frame(y = factor(c("a", "b", "c")[apply(P, 1, function(p) sample.int(3, 1, prob = p))], levels = c("a", "b", "c")), x, unit = factor(unit))
  mr <- suppressWarnings(remnl(y ~ x, data = dr, unit = "unit", nodes = 8L))
  s <- exp(mr$theta[mr$design$idx$str]); eta <- cbind(0, as.vector(mr$design$X[1:3, ] %*% mr$theta[1:2]), as.vector(mr$design$X[1:3, ] %*% mr$theta[3:4]))
  z <- seq(-8, 8, length.out = 401); wz <- dnorm(z); wz <- wz / sum(wz)
  fine <- t(sapply(1:3, function(i) { acc <- numeric(3)
    for (a in seq_along(z)) { e <- cbind(eta[i, 1], eta[i, 2] + s[1] * z[a], eta[i, 3] + s[2] * z)
      p <- exp(e - apply(e, 1, max)); p <- p / rowSums(p); acc <- acc + wz[a] * colSums(wz * p) }
    acc }))
  expect_equal(unname(fitted(mr)[1:3, ]), fine, tolerance = 1e-7)
})
