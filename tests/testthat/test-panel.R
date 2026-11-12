## Fixed-effects panels: the split-panel jackknife and its validity gate.

.sim_panel <- function(G, T, beta, sd_u = 1, shift = 0, trend = 0, seed = 1) {
  set.seed(seed)
  J <- ncol(beta) + 1L
  unit <- rep(seq_len(G), each = T); time <- rep(seq_len(T), G); n <- G * T
  x1 <- rnorm(n); x2 <- rnorm(n)
  U <- matrix(rnorm(G * (J - 1L), 0, sd_u), G)                         # alternative-specific unit effects
  V <- cbind(0, cbind(1, x1, x2) %*% beta + U[unit, ])
  V[, -1] <- V[, -1] + shift * x1 * (time > T / 2)                      # regime change: the x1 slopes jump at mid-panel
  V[, 2] <- V[, 2] + trend * time                                       # a smooth trend not in the formula
  P <- exp(V - apply(V, 1, max)); P <- P / rowSums(P)
  y <- apply(P, 1, function(pp) sample.int(J, 1L, prob = pp))
  data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1 = x1, x2 = x2, unit = factor(unit), time = time)
}
.B <- cbind(b = c(0.3, 0.8, -0.5), c = c(-0.2, -0.6, 0.7))

test_that("bias_correct = 'jackknife' corrects the slopes of a stationary fixed-effects panel and leaves the rest alone", {
  skip_on_cran()
  d <- .sim_panel(300, 8, .B, seed = 11)
  m0 <- suppressWarnings(mnl(y ~ x1 + x2, data = d, fe = "unit"))
  mj <- suppressWarnings(mnl(y ~ x1 + x2, data = d, fe = "unit", bias_correct = "jackknife"))
  expect_equal(mj$bias_correct, "jackknife")
  expect_equal(mj$uncorrected$coefficients, m0$coefficients)
  sl <- c("b:x1", "b:x2", "c:x1", "c:x2"); truth <- c(0.8, -0.5, -0.6, 0.7)
  expect_false(any(mj$coefficients[sl] == m0$coefficients[sl]))
  ## the ML slopes are inflated by the incidental-parameters bias; the correction pulls them in
  expect_gt(mean(abs(m0$coefficients[sl]) / abs(truth)), mean(abs(mj$coefficients[sl]) / abs(truth)))
  expect_lt(sum((mj$coefficients[sl] - truth)^2), sum((m0$coefficients[sl] - truth)^2))
  ## constants, unit dummies, log-likelihood, and standard errors are the ML ones
  expect_equal(mj$coefficients["b:(Intercept)"], m0$coefficients["b:(Intercept)"])
  fe_names <- mj$par_names[mj$blocks == "fixed effect"]
  expect_equal(mj$coefficients[fe_names], m0$coefficients[fe_names])
  expect_equal(mj$loglik, m0$loglik); expect_equal(vcov(mj), vcov(m0))
  expect_equal(dim(mj$jackknife$halves), c(4L, 2L))
  expect_output(print(mj), "jackknife")
  expect_error(mnl(y ~ x1 + x2, data = d, bias_correct = "jackknife"), "fe")
})

test_that("the gate refuses the correction on time-heterogeneous panels", {
  skip_on_cran()
  ## a regime change: the x1 slopes jump in the second half
  d1 <- .sim_panel(300, 8, .B, shift = 1.5, seed = 12)
  w1 <- character(0)
  m1 <- withCallingHandlers(mnl(y ~ x1 + x2, data = d1, fe = "unit", bias_correct = "jackknife"),
                            warning = function(cnd) { w1 <<- c(w1, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("REFUSED", w1))); expect_true(any(grepl("disagree|trend with the within-unit time order", w1)))   # a slope jump also moves the mean residuals over time
  expect_equal(m1$bias_correct, "none"); expect_null(m1$uncorrected)
  ## a smooth trend omitted from the formula: the residual-trend detector
  d2 <- .sim_panel(300, 8, .B, trend = 0.25, seed = 13)
  w2 <- character(0)
  m2 <- withCallingHandlers(mnl(y ~ x1 + x2, data = d2, fe = "unit", bias_correct = "jackknife"),
                            warning = function(cnd) { w2 <<- c(w2, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("REFUSED", w2))); expect_true(any(grepl("trend with the within-unit time order", w2)))
  ## with the trend in the formula the correction goes through
  w3 <- character(0)
  m3 <- withCallingHandlers(mnl(y ~ x1 + x2 + time, data = d2, fe = "unit", bias_correct = "jackknife"),
                            warning = function(cnd) { w3 <<- c(w3, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("REFUSED", w3))); expect_equal(m3$bias_correct, "jackknife")
})
