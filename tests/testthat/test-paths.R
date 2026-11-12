## Paths that the review panel found no test exercising (round 3, A9 and C6):
## each is run once on a small design and checked for its defining property.

small_panel <- function(G = 40, T = 5, seed = 1) {
  set.seed(seed); unit <- rep(seq_len(G), each = T); n <- G * T; x1 <- rnorm(n); z1 <- rnorm(n)
  V <- cbind(0, 0.3 + 0.8 * x1 + rnorm(G)[unit], -0.2 - 0.6 * x1 + rnorm(G)[unit]); P <- exp(V) / rowSums(exp(V))
  y <- apply(P, 1, function(p) sample.int(3, 1, prob = p)); y <- ifelse(runif(n) < plogis(1.5 + z1), y, 1L)
  data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1, z1, unit = factor(unit), g = factor((unit - 1) %/% 4))
}

test_that("start = reproduces the fit; a cluster bootstrap draws clusters", {
  d <- small_panel()
  m <- mnl(y ~ x1, data = d)
  ms <- mnl(y ~ x1, data = d, start = m$theta)
  expect_equal(ms$loglik, m$loglik, tolerance = 1e-10); expect_equal(coef(ms), coef(m), tolerance = 1e-6)
  set.seed(2); mb <- suppressWarnings(mnl(y ~ x1, data = d, se = "bootstrap", cluster = "g", nboot = 8))
  expect_true(mb$boot$cluster); expect_equal(nrow(mb$boot$replicates) + mb$boot$n_failed, 8L)
})

test_that("boot_pvalue() with a cluster-robust fit and boot_lr_test() for the interior and sd kinds", {
  skip_on_cran()
  d <- small_panel(60, 5, 3)
  mc <- mnl(y ~ x1, data = d, cluster = "g")
  set.seed(3); bp <- boot_pvalue(mc, hausman_mcfadden, B = 10, drop = "b", vcov = "joint")
  expect_true(is.finite(bp$p.boot))
  m0 <- mnl(y ~ x1, data = d)
  m1 <- suppressWarnings(suppressMessages(remnl(y ~ x1, data = d, unit = "unit", nodes = 8)))
  set.seed(4); b1 <- boot_lr_test(m0, m1, B = 10)
  expect_match(attr(b1, "method"), "conservative"); expect_true(all(attr(b1, "boot_statistics")["lr", ] >= 0))
  n1 <- suppressWarnings(suppressMessages(nlogit(y ~ x1, data = d, nests = list(n0 = "a", n1 = c("b", "c")))))
  set.seed(5); b2 <- boot_lr_test(m0, n1, B = 10)
  expect_true(is.finite(b2$p.value[1])); expect_true(is.finite(b2$p.boot[1]))
})

test_that("imnl() refuses unit fixed effects and points to the alternatives", {
  d <- small_panel(40, 12, 6)
  expect_error(imnl(y ~ x1 | z1, data = d, fe = "unit"), "does not take unit fixed effects")
  expect_error(imnl(y ~ x1 | z1, data = d, fe = "unit"), "cluster")
})

test_that("simulate() for the nested logit, HEV, and probit; a weighted random-intercept bootstrap", {
  skip_on_cran()
  d <- small_panel(50, 4, 8)
  for (m in list(suppressWarnings(suppressMessages(nlogit(y ~ x1, data = d, nests = list(n0 = "a", n1 = c("b", "c"))))),
                 suppressWarnings(suppressMessages(hev(y ~ x1, data = d))),
                 suppressWarnings(suppressMessages(mnprobit(y ~ x1, data = d))))) {
    s <- simulate(m, nsim = 2, seed = 1)
    expect_equal(dim(s), c(m$n, 2L)); expect_true(all(unlist(lapply(s, function(v) v %in% m$levels))))
  }
  d$w <- rep(1:2, length.out = nlevels(d$unit))[as.integer(d$unit)]
  set.seed(9); mr <- suppressWarnings(suppressMessages(remnl(y ~ x1, data = d, unit = "unit", nodes = 6, weights = "w", se = "bootstrap", nboot = 4)))
  expect_true(mr$boot$weighted); expect_true(is.finite(sqrt(vcov(mr)["b:x1", "b:x1"])))
})

test_that("first_difference(alternative =) in the long layout and a parallel bootstrap", {
  skip_on_cran(); skip_if_not_installed("AER")
  data("TravelMode", package = "AER"); TM <- TravelMode; TM$chosen <- as.integer(TM$choice == "yes")
  cl <- mnl(chosen ~ wait + vcost | income, data = TM, id = "individual", alt = "mode")
  fd <- first_difference(cl, "wait", 10, 20, alternative = "car")
  expect_lt(fd$diff[fd$alternative == "car"], 0)                             # a longer wait for the car: car less likely
  expect_equal(sum(fd$diff), 0, tolerance = 1e-10)
  d <- small_panel(30, 4, 10)
  set.seed(11); mp <- suppressWarnings(mnl(y ~ x1, data = d, se = "bootstrap", nboot = 4, cores = 2))
  expect_equal(nrow(mp$boot$replicates) + mp$boot$n_failed, 4L)
})
