## Analytic gradient, scores, and Hessian of the C++ engine against numDeriv,
## in every layout x inflation x link combination, and on extreme split indices.

.rand_design <- function(n, J, p, q, r, k, link, sparse) {
  X <- if (p > 0) cbind(1, matrix(rnorm(n * (p - 1)), n)) else matrix(0, n, 0)
  W <- if (q > 0) matrix(rnorm(n * J * q), n * J, q) else matrix(0, n * J, 0)
  Z <- if (k >= 0) cbind(1, matrix(rnorm(n * (r - 1)), n)) else matrix(0, n, 0)
  off <- matrix(rnorm(n * J, 0, 0.3), n, J); avail <- matrix(1L, n, J)
  if (sparse) { drop <- cbind(sample(n, n %/% 3), sample(J, n %/% 3, TRUE)); avail[drop] <- 0L }
  K <- p * (J - 1) + q + (if (k >= 0) r else 0)
  theta <- rnorm(K, 0, 0.7)
  pr <- pc_probs_cpp(theta, X, W, Z, off, avail, rep(0L, n), J, 0L, k, link)$prob
  y <- apply(pr, 1, function(pp) sample(J, 1, prob = pp)) - 1L
  list(X = X, W = W, Z = Z, off = off, avail = avail, y = y, w = runif(n, 0.5, 2), theta = theta, J = J, k = k, link = link)
}

test_that("gradient, scores, and Hessian match numerical derivatives", {
  set.seed(11)
  configs <- list(list(300, 3, 3, 0, 0, -1L, 1L, FALSE), list(300, 4, 3, 2, 0, -1L, 1L, TRUE),
                  list(300, 3, 3, 0, 2, 0L, 1L, FALSE), list(300, 3, 3, 0, 2, 0L, 0L, FALSE),
                  list(300, 4, 3, 2, 3, 2L, 1L, TRUE), list(300, 4, 2, 2, 2, 3L, 0L, TRUE), list(200, 5, 0, 3, 0, -1L, 1L, FALSE))
  for (cf in configs) {
    D <- do.call(.rand_design, cf)
    f <- function(th) sum(D$w * pc_lli_cpp(th, D$X, D$W, D$Z, D$off, D$avail, D$y, D$J, 0L, D$k, D$link))
    g <- pc_grad_cpp(D$theta, D$X, D$W, D$Z, D$off, D$avail, D$y, D$w, D$J, 0L, D$k, D$link)
    expect_equal(g, numDeriv::grad(f, D$theta), tolerance = 1e-6)
    S <- pc_scores_cpp(D$theta, D$X, D$W, D$Z, D$off, D$avail, D$y, D$J, 0L, D$k, D$link)
    expect_equal(colSums(D$w * S), g, tolerance = 1e-10)
    H <- pc_hess_cpp(D$theta, D$X, D$W, D$Z, D$off, D$avail, D$y, D$w, D$J, 0L, D$k, D$link)
    Hn <- numDeriv::hessian(f, D$theta)
    expect_equal(H, Hn, tolerance = 1e-5)
    expect_equal(H, t(H), tolerance = 1e-12)
    P <- pc_probs_cpp(D$theta, D$X, D$W, D$Z, D$off, D$avail, D$y, D$J, 0L, D$k, D$link)
    expect_equal(rowSums(P$prob), rep(1, nrow(D$X)), tolerance = 1e-12)
    expect_equal(rowSums(P$prob_outcome), rep(1, nrow(D$X)), tolerance = 1e-12)
  }
})

test_that("extreme split intercepts keep the derivatives finite and exact", {
  set.seed(3)
  n <- 200; J <- 3; X <- cbind(1, rnorm(n)); Z <- cbind(1, rnorm(n)); off <- matrix(0, n, J); avail <- matrix(1L, n, J); W <- matrix(0, n * J, 0)
  for (d0 in c(-12, -4, 4, 12)) {
    theta <- c(0.3, -0.5, -0.2, 0.8, d0, 0.5)
    pr <- pc_probs_cpp(theta, X, W, Z, off, avail, rep(0L, n), J, 0L, 0L, 1L)$prob
    y <- apply(pr, 1, function(pp) sample(J, 1, prob = pp)) - 1L
    f <- function(th) sum(pc_lli_cpp(th, X, W, Z, off, avail, y, J, 0L, 0L, 1L))
    g <- pc_grad_cpp(theta, X, W, Z, off, avail, y, rep(1, n), J, 0L, 0L, 1L)
    H <- pc_hess_cpp(theta, X, W, Z, off, avail, y, rep(1, n), J, 0L, 0L, 1L)
    expect_true(all(is.finite(g)) && all(is.finite(H)))
    expect_equal(g, numDeriv::grad(f, theta), tolerance = 1e-6)
    expect_equal(H, numDeriv::hessian(f, theta), tolerance = 1e-5)
  }
})
