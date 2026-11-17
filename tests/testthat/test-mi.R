## Multiple imputation: mice::pool() on a list of fits, through tidy() and
## glance(), gives Rubin's rules; the complete-data degrees of freedom are
## passed (mice 3.19 otherwise takes n J - p from the residual matrix).

rubin <- function(fits) {
  Q <- sapply(fits, coef); U <- sapply(fits, function(f) diag(vcov(f))); m <- length(fits)
  list(est = rowMeans(Q), se = sqrt(rowMeans(U) + (1 + 1 / m) * apply(Q, 1, var)))
}

test_that("mice::pool() of a list of mnl() fits equals Rubin's rules computed by hand", {
  skip_if_not_installed("mice"); skip_if_not_installed("broom")
  set.seed(11); d <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d$x1[sample(500, 80)] <- NA; d$x2[sample(500, 60)] <- NA
  imp <- mice::mice(d, m = 5, printFlag = FALSE, seed = 3)
  fits <- lapply(1:5, function(i) mnl(y ~ x1 + x2, data = mice::complete(imp, i)))
  p <- mice::pool(fits, dfcom = nobs(fits[[1]]) - length(coef(fits[[1]])))
  s <- summary(p); r <- rubin(fits); o <- match(names(r$est), as.character(s$term))
  expect_false(anyNA(o))
  expect_equal(s$estimate[o], unname(r$est), tolerance = 1e-10)
  expect_equal(s$std.error[o], unname(r$se), tolerance = 1e-10)
  expect_true(all(s$df <= 494))                                   # Barnard-Rubin: never above the complete-data value
  expect_true("singular" %in% names(broom::glance(fits[[1]])))
})

test_that("mice::pool() keeps both equations' terms of inflated fits apart", {
  skip_on_cran(); skip_if_not_installed("mice"); skip_if_not_installed("broom")
  set.seed(21); n <- 1500; X <- data.frame(x1 = rnorm(n), x2 = rnorm(n)); Z <- data.frame(z1 = rnorm(n))
  d <- rmnl(n, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), delta = cbind(a = c(0.5, 1.5)), inflate = "a", x = X, z = Z)
  d$x2[sample(n, 150)] <- NA
  imp <- mice::mice(d, m = 3, printFlag = FALSE, seed = 5)
  fits <- lapply(1:3, function(i) suppressMessages(imnl(y ~ x1 + x2 | z1, data = mice::complete(imp, i), inflate = "a")))
  expect_true(all(vapply(fits, function(f) f$converged && !f$boundary, logical(1))))
  s <- summary(mice::pool(fits, dfcom = nobs(fits[[1]]) - length(coef(fits[[1]]))))
  r <- rubin(fits); o <- match(names(r$est), as.character(s$term))
  expect_false(anyNA(o)); expect_true(any(startsWith(names(r$est), "infl_")))
  expect_equal(s$estimate[o], unname(r$est), tolerance = 1e-10)
  expect_equal(s$std.error[o], unname(r$se), tolerance = 1e-10)
})
