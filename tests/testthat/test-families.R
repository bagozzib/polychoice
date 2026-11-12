## The non-MNL families: nested logit and HEV. Derivatives of the C++ kernels
## against numDeriv, agreement with mlogit on TravelMode, reduction to the
## MNL, parameter recovery, and the methods on these fits.

test_that("nested-logit kernel: gradient and scores match numDeriv; lambda = 1 is the MNL", {
  set.seed(4)
  n <- 250; J <- 4; p <- 3; q <- 1
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n)); W <- matrix(rnorm(n * J * q), n * J, q)
  off <- matrix(rnorm(n * J, 0, 0.2), n, J); avail <- matrix(1L, n, J); avail[cbind(sample(n, 60), sample(J, 60, TRUE))] <- 0L
  nest <- c(0L, 0L, 1L, 1L); lmap <- c(0L, 1L)
  theta <- c(rnorm(p * (J - 1) + q, 0, 0.6), -0.3, 0.2)
  pr <- pc_nl_probs_cpp(theta, X, W, off, avail, rep(0L, n), J, 0L, nest, lmap)$prob
  y <- apply(pr, 1, function(pp) sample(J, 1, prob = pp)) - 1L; w <- runif(n, 0.5, 2)
  f <- function(th) sum(w * pc_nl_lli_cpp(th, X, W, off, avail, y, J, 0L, nest, lmap))
  g <- pc_nl_grad_cpp(theta, X, W, off, avail, y, w, J, 0L, nest, lmap)
  expect_equal(g, numDeriv::grad(f, theta), tolerance = 1e-6)
  S <- pc_nl_scores_cpp(theta, X, W, off, avail, y, J, 0L, nest, lmap)
  expect_equal(colSums(w * S), g, tolerance = 1e-10)
  expect_equal(rowSums(pr), rep(1, n), tolerance = 1e-12)
  th1 <- theta; th1[(length(th1) - 1):length(th1)] <- 0
  l_nl <- pc_nl_lli_cpp(th1, X, W, off, avail, y, J, 0L, nest, lmap)
  l_mnl <- pc_lli_cpp(th1[1:(length(th1) - 2)], X, W, matrix(0, n, 0), off, avail, y, J, 0L, -1L, 1L)
  expect_equal(l_nl, l_mnl, tolerance = 1e-12)
})

test_that("HEV kernel: gradient and scores match numDeriv; unit scales give the MNL; probabilities sum to one", {
  skip_on_cran()                                    # slow; runs in the full suite (check time, 2026-11-06)
  set.seed(5)
  n <- 200; J <- 4; p <- 3; q <- 1
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n)); W <- matrix(rnorm(n * J * q), n * J, q)
  off <- matrix(0, n, J); avail <- matrix(1L, n, J); avail[cbind(sample(n, 40), sample(J, 40, TRUE))] <- 0L
  theta <- c(rnorm(p * (J - 1) + q, 0, 0.6), 0.3, -0.4, 0.5)
  pr <- pc_hev_probs_cpp(theta, X, W, off, avail, rep(0L, n), J, 0L)$prob
  y <- apply(pr, 1, function(pp) sample(J, 1, prob = pp)) - 1L; w <- runif(n, 0.5, 2)
  f <- function(th) sum(w * pc_hev_lli_cpp(th, X, W, off, avail, y, J, 0L))
  g <- pc_hev_grad_cpp(theta, X, W, off, avail, y, w, J, 0L)
  expect_equal(g, numDeriv::grad(f, theta), tolerance = 1e-6)
  S <- pc_hev_scores_cpp(theta, X, W, off, avail, y, J, 0L)
  expect_equal(colSums(w * S), g, tolerance = 1e-10)
  expect_equal(rowSums(pr), rep(1, n), tolerance = 1e-9)                       # quadrature error only
  th1 <- theta; th1[(length(th1) - 2):length(th1)] <- 0
  l_hev <- pc_hev_lli_cpp(th1, X, W, off, avail, y, J, 0L)
  l_mnl <- pc_lli_cpp(th1[1:(length(th1) - 3)], X, W, matrix(0, n, 0), off, avail, y, J, 0L, -1L, 1L)
  expect_equal(l_hev, l_mnl, tolerance = 1e-9)
  ## extreme scale ratios and utility gaps stay finite and accurate
  th2 <- theta; th2[(length(th2) - 2):length(th2)] <- c(1.5, -1.5, 0); th2[1] <- 6
  pr2 <- pc_hev_probs_cpp(th2, X, W, off, avail, rep(0L, n), J, 0L)$prob
  expect_equal(rowSums(pr2), rep(1, n), tolerance = 1e-8)
  expect_true(all(is.finite(pc_hev_grad_cpp(th2, X, W, off, avail, y, w, J, 0L))))
  ## independent oracle: integrate() on the probability formula, in the
  ## trapezoid regime (ratios <= 20) and the adaptive Gauss-Kronrod regime
  ## (a scale pushed towards zero makes a competitor's factor a near step)
  Vfun <- function(th) { V <- matrix(0, n, J); for (j in 2:J) V[, j] <- X %*% th[(p * (j - 2) + 1):(p * (j - 1))]; V + matrix(W %*% th[p * (J - 1) + 1], n, J, byrow = TRUE) }
  oracle <- function(th, i, j) { V <- Vfun(th); thv <- c(1, exp(th[(length(th) - 2):length(th)])); ks <- setdiff(which(avail[i, ] == 1), j)
    f <- function(ww) sapply(ww, function(w) exp(-sum(exp(-(V[i, j] - V[i, ks]) / thv[ks] - thv[j] / thv[ks] * w)) - w - exp(-w)))
    integrate(f, -Inf, Inf, rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value }
  full <- which(rowSums(avail) == 4)[1:4]
  ## the gradient is checked in every regime, up to a scale ratio of 1.6e5
  ## (before the step breakpoints the sharpest regimes failed here), and the
  ## finite-difference Hessian of the analytic gradient (which the standard
  ## errors use) against numDeriv's Hessian of the log-likelihood
  for (tv in c(0, -6, -9, -12, 8)) {
    th3 <- theta; th3[length(theta)] <- tv
    pr3 <- pc_hev_probs_cpp(th3, X, W, off, avail, rep(0L, n), J, 0L)$prob
    expect_equal(rowSums(pr3), rep(1, n), tolerance = 1e-9)
    for (i in full) expect_equal(unname(pr3[i, ]), sapply(1:4, function(j) oracle(th3, i, j)), tolerance = 1e-9)
    f3 <- function(th) sum(w * pc_hev_lli_cpp(th, X, W, off, avail, y, J, 0L))
    g3 <- function(th) pc_hev_grad_cpp(th, X, W, off, avail, y, w, J, 0L)
    expect_equal(g3(th3), numDeriv::grad(f3, th3), tolerance = 1e-5)
    if (tv %in% c(0, -9)) expect_equal(polychoice:::.pc_hess_fd(g3, th3), numDeriv::hessian(f3, th3), tolerance = 1e-4, ignore_attr = TRUE)
  }
  ## the theta_d -> 0 limit: P(d) = prod_k F((V_d - V_k) / theta_k)
  th4 <- theta; th4[length(theta)] <- -20
  pr4 <- pc_hev_probs_cpp(th4, X, W, off, avail, rep(0L, n), J, 0L)$prob
  V4 <- Vfun(th4); thv <- c(1, exp(th4[(length(th4) - 2):length(th4)]))
  Pd <- exp(-rowSums(sapply(1:3, function(k) avail[, k] * exp(-(V4[, 4] - V4[, k]) / thv[k]))))
  ok <- avail[, 4] == 1
  expect_equal(pr4[ok, 4], Pd[ok], tolerance = 1e-8)
})

test_that("nlogit() reproduces mlogit on TravelMode (equal and free inclusive-value parameters)", {
  skip_if_not_installed("mlogit"); skip_if_not_installed("dfidx"); skip_if_not_installed("AER")
  data("TravelMode", package = "AER"); TM <- TravelMode; TM$chosen <- as.integer(TM$choice == "yes")
  ## mlogit reports BHHH (outer-product) standard errors for the nested logit and HEV: compare like with like
  nests <- list(fast = c("air", "train"), slow = c("bus", "car"))
  m1 <- nlogit(chosen ~ wait + vcost + travel + gcost | income + size, data = TM, id = "individual", alt = "mode", nests = nests, equal_lambda = TRUE, se = "opg")
  expect_true(m1$converged)
  expect_equal(m1$loglik, -168.9346, tolerance = 1e-6)                          # mlogit un.nest.el = TRUE
  expect_equal(unname(coef(m1)["lambda"]), 1.49823, tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(m1)))["lambda"]), 0.66727, tolerance = 1e-3)
  expect_equal(unname(coef(m1)[c("wait", "vcost", "travel", "gcost")]), c(-0.13244, -0.06183, -0.01230, 0.05421), tolerance = 1e-3)
  m2 <- nlogit(chosen ~ wait + vcost + travel + gcost | income + size, data = TM, id = "individual", alt = "mode", nests = nests, se = "opg")
  expect_equal(m2$loglik, -165.2235, tolerance = 1e-6)                          # mlogit un.nest.el = FALSE
  expect_equal(unname(coef(m2)[c("lambda:fast", "lambda:slow")]), c(1.79796, 0.79800), tolerance = 1e-4)
  expect_equal(unname(sqrt(diag(vcov(m2)))[c("lambda:fast", "lambda:slow")]), c(0.84629, 0.35080), tolerance = 1e-3)
  expect_equal(unname(coef(m2)[c("train:(Intercept)", "train:income")]), c(1.12262, -0.08694), tolerance = 1e-3)
  ## the LR test of IIA against the mnl fit
  m0 <- mnl(chosen ~ wait + vcost + travel + gcost | income + size, data = TM, id = "individual", alt = "mode")
  lr <- lr_test(m0, m2); expect_equal(lr$df, 2L); expect_equal(lr$statistic, 2 * (m2$loglik - m0$loglik))
  ## probabilities: within-nest and marginal consistent
  P <- predict(m2); Pc <- predict(m2, type = "prob_conditional")
  expect_equal(rowSums(P), rep(1, nrow(P)), ignore_attr = TRUE, tolerance = 1e-12)
  expect_equal(unname(predict(m2, newdata = TM[1:8, ])), unname(P[1:2, ]), tolerance = 1e-10)
  expect_s3_class(m2, "nlogit")
})

test_that("hev() on TravelMode: the likelihood is exact where mlogit's is not, and the maximum is at the boundary", {
  skip_on_cran()
  skip_if_not_installed("mlogit"); skip_if_not_installed("dfidx"); skip_if_not_installed("AER")
  data("TravelMode", package = "AER"); TM <- TravelMode; TM$chosen <- as.integer(TM$choice == "yes")
  ## mlogit (heterosc = TRUE) reports -172.4258 at scales train 0.83021, bus
  ## 0.56441, car 0.07192 with wait -0.06022, gcost -0.01205 and intercepts /
  ## income / size as below. Its Gauss-Laguerre quadrature is inaccurate at a
  ## scale ratio of 14: the exact log-likelihood at that point is -176.19
  ## (integrate() oracle, ledger), and the likelihood keeps rising as the car
  ## scale goes to zero (a deterministic car utility), so the maximum is at
  ## the boundary, where the information in the car scale vanishes.
  w <- character(0)
  m <- withCallingHandlers(hev(chosen ~ wait + gcost | income + size, data = TM, id = "individual", alt = "mode"),
                           warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("boundary", w))); expect_false(any(grepl("saddle|not positive definite", w)))
  expect_true(m$converged)
  expect_gt(m$loglik, -169.98); expect_lt(m$loglik, -169.96)
  expect_lt(unname(coef(m)["scale:car"]), 0.01)
  expect_true(m$boundary)
  D <- m$design; fn <- polychoice:::.pc_fns(D)
  ## at the boundary every competitor of the car carries a step of width
  ## 1 / r (r about 3000 here); the kernel resolves it smoothly, so the analytic
  ## gradient matches numerical differentiation (an earlier panel layout left
  ## the step in a node-free sliver now and then: the gradient was off by 1.1,
  ## and the finite-difference Hessian showed spurious negative curvature)
  expect_equal(fn$grad(m$theta), numDeriv::grad(fn$nll, m$theta), tolerance = 1e-5)
  ## the kernel at mlogit's solution equals the integrate() oracle and is far
  ## below mlogit's reported value
  ml <- c("train:(Intercept)" = -0.81054, "bus:(Intercept)" = -1.51838, "car:(Intercept)" = -3.15405, wait = -0.06022, gcost = -0.01205,
          "train:income" = -0.04753, "bus:income" = -0.00776, "car:income" = -0.01367, "train:size" = 0.74786, "bus:size" = 0.37309, "car:size" = 0.34308,
          "scale:train" = log(0.83021), "scale:bus" = log(0.56441), "scale:car" = log(0.07192))
  th_ml <- ml[names(m$theta)]
  ll_k <- -fn$nll(th_ml)
  P <- fn$probs(th_ml); V <- P$V; thv <- c(air = 1, train = 0.83021, bus = 0.56441, car = 0.07192)[D$levels]
  ll_or <- 0
  for (i in seq_len(D$n)) { yy <- D$y[i]; ks <- setdiff(which(D$avail[i, ] == 1), yy)
    f <- function(ww) sapply(ww, function(w) exp(-sum(exp(-(V[i, yy] - V[i, ks]) / thv[ks] - thv[yy] / thv[ks] * w)) - w - exp(-w)))
    ll_or <- ll_or + log(integrate(f, -Inf, Inf, rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value) }
  expect_equal(ll_k, ll_or, tolerance = 1e-8)
  expect_lt(ll_k, -176); expect_gt(ll_k, -176.5)
  ## the profile in the car scale is flat at the boundary and lower at mlogit's value
  th <- m$theta; th["scale:car"] <- log(0.07192); expect_lt(-fn$nll(th), m$loglik - 0.1)
  th["scale:car"] <- -12; expect_equal(-fn$nll(th), m$loglik, tolerance = 1e-4)
  m0 <- mnl(chosen ~ wait + gcost | income + size, data = TM, id = "individual", alt = "mode")
  lr <- lr_test(m0, m); expect_equal(lr$df, 3L); expect_gt(lr$statistic, 0)
  expect_error(lr_test(m0, nlogit(chosen ~ wait + gcost | income, data = TM, id = "individual", alt = "mode", nests = list(a = c("air", "train"), b = c("bus", "car")))), "same covariates")
})

test_that("nlogit() and hev() recover their parameters and behave on simulated data", {
  skip_on_cran()
  ## alternative-specific covariates (long layout) identify lambda and the scales sharply (ledger)
  set.seed(3)
  nests <- list(n1 = c("a", "b"), n2 = c("c", "d"))
  B4 <- cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3))
  d <- rnlogit(1500, beta = B4, nests = nests, lambda = c(0.6, 0.8), gamma = 1)
  expect_message(m <- nlogit(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt", nests = list(n1 = c("a", "b"), n2 = "c", n3 = "d")), "Singleton")
  m <- nlogit(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt", nests = nests)
  expect_true(m$converged)
  truth <- c(as.vector(B4), 1, 0.6, 0.8)
  z <- (coef(m) - truth) / sqrt(diag(vcov(m)))
  expect_lt(max(abs(z)), 3.5)
  expect_equal(unname(coef(m)[c("lambda:n1", "lambda:n2")]), c(0.6, 0.8), tolerance = 0.2)
  expect_equal(rowSums(fitted(m)), rep(1, m$n), ignore_attr = TRUE, tolerance = 1e-12)
  fd <- first_difference(m, "x1", 0, 1); expect_equal(nrow(fd), 4L); expect_true(all(is.finite(fd$lower)))
  expect_error(nlogit(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt", nests = list(n1 = c("a", "b"), n2 = c("c"))), "must assign every alternative")
  expect_error(nlogit(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt", nests = list(n1 = c("a", "b", "c", "d"))), "at least two nests")
  ## the wide layout identifies lambda only through the functional form: judge recovery on the SE scale
  dw <- rnlogit(3000, beta = B4, nests = nests, lambda = c(0.6, 0.8))
  mw <- nlogit(y ~ x1 + x2, data = dw, nests = nests)
  zl <- (coef(mw)[c("lambda:n1", "lambda:n2")] - c(0.6, 0.8)) / sqrt(diag(vcov(mw)))[c("lambda:n1", "lambda:n2")]
  expect_lt(max(abs(zl)), 3.5)
  ## HEV
  dh <- rhev(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), scale = c(1, 0.6, 1.5), gamma = 1)
  mh <- hev(chosen ~ w1 | x1 + x2, data = dh, id = "id", alt = "alt")
  expect_true(mh$converged)
  zh <- (coef(mh)[c("scale:b", "scale:c")] - c(0.6, 1.5)) / sqrt(diag(vcov(mh)))[c("scale:b", "scale:c")]
  expect_lt(max(abs(zh)), 3.5)
  expect_equal(unname(coef(mh)[c("scale:b", "scale:c")]), c(0.6, 1.5), tolerance = 0.25)
  expect_equal(rowSums(fitted(mh)), rep(1, mh$n), ignore_attr = TRUE, tolerance = 1e-8)
  s <- summary(mh); expect_true(any(s$block == "scale"))
  expect_error(rrr(mh), "not odds ratios")                                  # exp(coef) has no odds reading under HEV
  cm <- compare_models(mnl = mnl(chosen ~ w1 | x1 + x2, data = dh, id = "id", alt = "alt"), hev = mh); expect_equal(cm$df, c(7L, 9L))
  ## bootstrap replicates are on the natural scale
  set.seed(2)
  sub <- dh[dh$id <= 300, ]
  mb <- suppressWarnings(hev(chosen ~ w1 | x1 + x2, data = sub, id = "id", alt = "alt", se = "bootstrap", nboot = 6))
  expect_true(all(mb$boot$replicates[, "scale:b"] > 0))
  expect_equal(dim(confint(mb, type = "percentile")), c(9L, 2L))
})
