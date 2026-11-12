## The multinomial probit: exact probabilities against mvtnorm, derivatives
## against numDeriv (exact J = 3 and GHK J = 4), parameter recovery, and the
## scale relation with the logit.

.mnp_design <- function(n, J, p, q, seed, sparse = FALSE) {
  set.seed(seed)
  X <- cbind(1, matrix(rnorm(n * (p - 1)), n)); W <- if (q > 0) matrix(rnorm(n * J * q), n * J, q) else matrix(0, n * J, 0)
  off <- matrix(0, n, J); avail <- matrix(1L, n, J)
  if (sparse) avail[cbind(sample(n, n %/% 5), sample(J, n %/% 5, TRUE))] <- 0L
  m <- J - 1; nc <- m * (m + 1) / 2 - 1
  cpar <- c(rnorm(nc, 0, 0.3))
  theta <- c(rnorm(p * (J - 1) + q, 0, 0.5), cpar)
  list(X = X, W = W, off = off, avail = avail, theta = theta, J = J, nc = nc)
}

test_that("three-alternative probabilities equal mvtnorm::pmvnorm and derivatives match numDeriv", {
  skip_if_not_installed("mvtnorm")
  D <- .mnp_design(60, 3, 2, 1, 7)
  n <- 60; J <- 3
  pr <- pc_mnp_probs_cpp(D$theta, D$X, D$W, D$off, D$avail, rep(0L, n), J, 0L, 50L)
  expect_equal(rowSums(pr$prob), rep(1, n), tolerance = 1e-10)
  ## recompute case 3's probabilities from the covariance of the differences
  S <- pc_mnp_sigma_cpp(D$theta[(length(D$theta) - D$nc + 1):length(D$theta)], 2L)
  V <- pr$V[3, ]
  ## P(y = base): eta_k <= V_b - V_k, eta ~ N(0, S)
  p0 <- mvtnorm::pmvnorm(upper = c(V[1] - V[2], V[1] - V[3]), sigma = S)
  ## P(y = alt 2): eps_1 - eps_2 = -eta_1 <= V_2 - V_1 and eps_3 - eps_2 = eta_2 - eta_1 <= V_2 - V_3
  A2 <- rbind(c(-1, 0), c(-1, 1)); p1 <- mvtnorm::pmvnorm(upper = c(V[2] - V[1], V[2] - V[3]), sigma = A2 %*% S %*% t(A2))
  A3 <- rbind(c(0, -1), c(1, -1)); p2 <- mvtnorm::pmvnorm(upper = c(V[3] - V[1], V[3] - V[2]), sigma = A3 %*% S %*% t(A3))
  expect_equal(unname(pr$prob[3, ]), c(p0, p1, p2), tolerance = 1e-8)
  y <- apply(pr$prob, 1, function(pp) sample(J, 1, prob = pp)) - 1L; w <- runif(n, 0.5, 2)
  f <- function(th) sum(w * pc_mnp_lli_cpp(th, D$X, D$W, D$off, D$avail, y, J, 0L, 50L))
  g <- pc_mnp_grad_cpp(D$theta, D$X, D$W, D$off, D$avail, y, w, J, 0L, 50L)
  expect_equal(g, numDeriv::grad(f, D$theta), tolerance = 1e-5)
  Sc <- pc_mnp_scores_cpp(D$theta, D$X, D$W, D$off, D$avail, y, J, 0L, 50L)
  expect_equal(colSums(w * Sc), g, tolerance = 1e-10)
})

test_that("GHK (four alternatives): probabilities near mvtnorm, derivatives match numDeriv with fixed draws", {
  skip_if_not_installed("mvtnorm")
  D <- .mnp_design(40, 4, 2, 1, 8, sparse = TRUE)
  n <- 40; J <- 4; R <- 400L
  pr <- pc_mnp_probs_cpp(D$theta, D$X, D$W, D$off, D$avail, rep(0L, n), J, 0L, R)
  full <- which(rowSums(D$avail) == 4L)
  expect_lt(max(abs(rowSums(pr$prob[full, ]) - 1)), 0.01)
  S <- pc_mnp_sigma_cpp(D$theta[(length(D$theta) - D$nc + 1):length(D$theta)], 3L)
  i <- full[1]; V <- pr$V[i, ]
  p0 <- mvtnorm::pmvnorm(upper = V[1] - V[2:4], sigma = S)
  expect_equal(unname(pr$prob[i, 1]), as.numeric(p0), tolerance = 2e-3)
  y <- apply(pr$prob, 1, function(pp) { pp[D$avail[1, ] == 0] <- 0; sample(J, 1, prob = pmax(pp, 0)) }) - 1L
  for (ii in seq_len(n)) if (D$avail[ii, y[ii] + 1] == 0L) y[ii] <- which(D$avail[ii, ] == 1L)[1] - 1L
  w <- runif(n, 0.5, 2)
  f <- function(th) sum(w * pc_mnp_lli_cpp(th, D$X, D$W, D$off, D$avail, y, J, 0L, R))
  g <- pc_mnp_grad_cpp(D$theta, D$X, D$W, D$off, D$avail, y, w, J, 0L, R)
  expect_equal(g, numDeriv::grad(f, D$theta), tolerance = 1e-5)
})

test_that("mnprobit() recovers the parameters, reports the covariance, and relates to the logit", {
  skip_on_cran()
  set.seed(5)
  sig <- matrix(c(1, 0.6, 0.6, 1.5), 2)
  B3 <- cbind(b = c(0.3, -0.6, 0.5), c = c(-0.2, 0.4, -0.3))
  ## alternative-specific covariates identify the covariance (Keane 1992; see ?mnprobit)
  d <- rmnp(2000, beta = B3, sigma = sig, gamma = 1)
  m <- mnprobit(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
  expect_true(m$converged); expect_s3_class(m, "mnprobit"); expect_false(m$boundary)
  truth <- c(as.vector(B3), 1, 0.6, 1.5)
  z <- (coef(m) - truth) / sqrt(diag(vcov(m)))
  expect_lt(max(abs(z)), 3.5)
  expect_equal(unname(coef(m)[c("cov:c,b", "var:c")]), c(0.6, 1.5), tolerance = 0.3)
  expect_equal(rowSums(fitted(m)), rep(1, m$n), ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(unname(predict(m, newdata = d[d$id <= 5, ])), unname(fitted(m)[1:5, ]), tolerance = 1e-10)
  s <- summary(m); expect_true(any(s$block == "covariance"))
  ml <- mnl(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
  expect_error(lr_test(ml, m), "does not nest")
  v <- vuong(m, ml); expect_true(is.finite(v$statistic[1]))
  ## iid errors: logit coefficients are about sqrt(pi^2 / 3) = 1.8 times the probit ones
  set.seed(6)
  d0 <- rmnp(2000, beta = B3, sigma = (diag(2) + 1) / 2, gamma = 1)
  mp <- mnprobit(chosen ~ w1 | x1 + x2, data = d0, id = "id", alt = "alt")
  ml0 <- mnl(chosen ~ w1 | x1 + x2, data = d0, id = "id", alt = "alt")
  expect_equal(unname(coef(ml0)["w1"] / coef(mp)["w1"]), sqrt(pi^2 / 3), tolerance = 0.1)
  ## with case-specific covariates only the covariance is fragilely identified:
  ## the fit runs towards a singular covariance and says so
  set.seed(5)
  dw <- rmnp(1500, beta = B3, sigma = sig)
  mw <- suppressWarnings(mnprobit(y ~ x1 + x2, data = dw)); expect_s3_class(mw, "mnprobit")   # may or may not converge (ledger)
})

test_that("the boundary warning points the probit to vuong() and, under GHK, says the likelihood is rough there", {
  msg <- function(fam, J) polychoice:::.pc_str_boundary_msg(list(family = fam, J = J), "var:c")
  expect_true(grepl("by vuong()", msg("mnp", 3L), fixed = TRUE))
  expect_false(grepl("GHK", msg("mnp", 3L)))                      # three alternatives: exact bivariate normal
  expect_true(grepl("GHK simulated likelihood is not smooth", msg("mnp", 4L)))
  expect_true(grepl("by lr_test()", msg("hev", 4L), fixed = TRUE))
  expect_false(grepl("GHK", msg("hev", 4L)))
})
