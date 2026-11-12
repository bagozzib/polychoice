## Random alternative-specific intercepts: the Gauss-Hermite rule, kernel
## derivatives, reduction to the MNL, the glmer oracle at J = 2, recovery, and
## the methods.

.sim_re <- function(G, T, beta, sd_u, gamma = NULL, seed = 1) {
  set.seed(seed)
  J <- ncol(beta) + 1L; levels <- c("a", colnames(beta)); n <- G * T
  unit <- rep(seq_len(G), each = T)
  X <- cbind(1, matrix(rnorm(n * (nrow(beta) - 1L)), n))
  U <- sweep(matrix(rnorm(G * (J - 1L)), G), 2L, sd_u, "*")
  V <- cbind(0, X %*% beta + U[unit, , drop = FALSE])
  Wa <- NULL
  if (!is.null(gamma)) { Wa <- matrix(rnorm(n * J), n, J); V <- V + Wa * gamma }
  P <- exp(V - apply(V, 1L, max)); P <- P / rowSums(P)
  y <- apply(P, 1L, function(pp) sample.int(J, 1L, prob = pp))
  xdf <- as.data.frame(X[, -1L, drop = FALSE]); names(xdf) <- paste0("x", seq_len(ncol(xdf)))
  if (is.null(gamma)) return(data.frame(y = factor(levels[y], levels = levels), xdf, unit = factor(unit)))
  rows <- rep(seq_len(n), each = J)
  data.frame(id = rows, alt = factor(rep(levels, n), levels = levels), chosen = as.integer(rep(y, each = J) == rep(seq_len(J), n)),
             w1 = as.vector(t(Wa)), xdf[rows, , drop = FALSE], unit = factor(unit[rows]))
}
.B3 <- cbind(b = c(0.3, 0.8, -0.5), c = c(-0.2, -0.6, 0.7))

test_that("the Gauss-Hermite rule integrates normal moments", {
  gh <- polychoice:::.pc_gauss_hermite(8L)
  expect_equal(sum(gh$weights), 1)
  expect_equal(sum(gh$weights * gh$nodes), 0, tolerance = 1e-12)
  expect_equal(sum(gh$weights * gh$nodes^2), 1, tolerance = 1e-12)
  expect_equal(sum(gh$weights * gh$nodes^4), 3, tolerance = 1e-10)
  expect_equal(sum(gh$weights * exp(gh$nodes)), exp(0.5), tolerance = 1e-8)
})

test_that("random-intercept kernel: gradient and scores match numDeriv; sd -> 0 gives the MNL", {
  d <- .sim_re(60, 4, .B3, sd_u = c(0.8, 1.2), seed = 2)
  ## the analytic gradient is the quadrature of the exact gradient, numDeriv the
  ## derivative of the quadrature (whose adaptive nodes move with theta); the two
  ## agree as the rule is refined (max difference 1.7e-2 at 6 nodes, 1.1e-3 at
  ## 8, 7.7e-7 at 12, 1.9e-8 at 20), so the comparison uses 20 nodes
  D <- polychoice:::.pc_design(y ~ x1 + x2, d, family = "remnl", str = list(unit = "unit", nodes = 20L))
  fn <- polychoice:::.pc_fns(D)
  expect_equal(D$str$G, 60L); expect_equal(dim(D$str$Z), c(400L, 2L)); expect_equal(sum(exp(D$str$logw)), 1)
  set.seed(1); theta <- c(rnorm(6, 0, 0.5), log(0.8), log(1.2))
  expect_equal(fn$grad(theta), numDeriv::grad(fn$nll, theta), tolerance = 1e-6)
  expect_equal(length(fn$lli(theta)), 60L)
  expect_equal(colSums(fn$scores(theta)), -fn$grad(theta), tolerance = 1e-10)
  ## sd -> 0: the unit log-likelihood is the sum of the MNL contributions of its cases
  theta0 <- theta; theta0[7:8] <- -20
  fn0 <- polychoice:::.pc_fns(polychoice:::.pc_as_mnl(D))
  expect_equal(fn$lli(theta0), as.vector(rowsum(fn0$lli(theta[1:6]), as.integer(D$unit))), tolerance = 1e-10)
  expect_equal(fn$probs(theta0)$prob, fn0$probs(theta[1:6])$prob, tolerance = 1e-10)
  ## weights must be constant within units
  d$w <- rep(c(1, 2), length.out = nrow(d))
  expect_error(remnl(y ~ x1 + x2, data = d, unit = "unit", weights = "w"), "constant within a unit")
})

test_that("at J = 2 the model is the binomial GLMM: lme4::glmer with adaptive quadrature agrees", {
  skip_on_cran()                                    # slow; runs in the full suite (check time, 2026-11-06)
  skip_if_not_installed("lme4")
  d <- .sim_re(150, 6, cbind(b = c(0.4, 0.9, -0.6)), sd_u = 1.1, seed = 3)
  m <- remnl(y ~ x1 + x2, data = d, unit = "unit", nodes = 40L)
  g <- lme4::glmer(I(y == "b") ~ x1 + x2 + (1 | unit), data = d, family = binomial, nAGQ = 25L)
  expect_equal(m$loglik, as.numeric(logLik(g)), tolerance = 1e-8)              # 40 nodes: 5e-9 apart in the ledger
  expect_equal(unname(coef(m)[c("b:(Intercept)", "b:x1", "b:x2")]), unname(lme4::fixef(g)), tolerance = 1e-3)
  expect_equal(unname(coef(m)["sd:b"]), sqrt(as.numeric(lme4::VarCorr(g)$unit)[1]), tolerance = 1e-3)
  expect_equal(unname(sqrt(diag(vcov(m)))[1:3]), unname(sqrt(diag(as.matrix(vcov(g))))), tolerance = 0.02)
})

test_that("the adaptive rule stays accurate on long panels with a large sd (where a fixed rule fails)", {
  skip_on_cran(); skip_if_not_installed("lme4")
  ## 150 units x 50 cases, sd 2: the fixed 12-node rule centred at zero was off
  ## by 22 log-likelihood units and 25 percent in sd here; the adaptive rule is
  ## within 1e-4 (verification/remnl_adaptive_check.R)
  set.seed(7); G <- 150; T <- 50; n <- G * T; unit <- rep(seq_len(G), each = T)
  x1 <- rnorm(n); x2 <- rnorm(n); u <- rnorm(G, 0, 2)
  p <- plogis(0.4 + 0.9 * x1 - 0.6 * x2 + u[unit])
  d <- data.frame(y = factor(ifelse(runif(n) < p, "b", "a"), levels = c("a", "b")), x1, x2, unit = factor(unit))
  m <- remnl(y ~ x1 + x2, data = d, unit = "unit")
  g <- lme4::glmer(I(y == "b") ~ x1 + x2 + (1 | unit), data = d, family = binomial, nAGQ = 25L)
  expect_lt(abs(m$loglik - as.numeric(logLik(g))), 1e-3)
  expect_equal(unname(coef(m)["sd:b"]), sqrt(as.numeric(lme4::VarCorr(g)$unit)[1]), tolerance = 1e-4)
})

test_that("remnl() recovers the parameters and behaves on simulated data", {
  skip_on_cran()
  d <- .sim_re(400, 5, .B3, sd_u = c(0.8, 1.2), seed = 4)
  m <- remnl(y ~ x1 + x2, data = d, unit = "unit", nodes = 8L)
  expect_true(m$converged); expect_s3_class(m, "remnl")
  truth <- c(0.3, 0.8, -0.5, -0.2, -0.6, 0.7, 0.8, 1.2)
  z <- (coef(m) - truth) / sqrt(diag(vcov(m))); expect_lt(max(abs(z)), 3.5)
  expect_equal(unname(coef(m)[c("sd:b", "sd:c")]), c(0.8, 1.2), tolerance = 0.25)
  expect_equal(length(m$lli), 400L)
  expect_equal(rowSums(fitted(m)), rep(1, 2000), ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(unname(predict(m, newdata = d[1:5, ])), unname(fitted(m)[1:5, ]), tolerance = 1e-10)
  s <- summary(m); expect_true(any(s$block == "random intercept"))
  expect_output(print(m), "random alternative-specific intercepts")
  m0 <- mnl(y ~ x1 + x2, data = d)
  expect_message(lr <- lr_test(m0, m), "boundary"); expect_equal(lr$df, 2L); expect_gt(lr$statistic, 10)
  v <- vuong(m, m0); expect_true(is.finite(v$statistic[1])); expect_gt(v$statistic[1], 0)
  expect_error(hausman_mcfadden(m), "does not apply")
  mr <- remnl(y ~ x1 + x2, data = d, unit = "unit", se = "robust", nodes = 8L); expect_equal(dim(vcov(mr)), c(8L, 8L))
  set.seed(1)
  mb <- remnl(y ~ x1 + x2, data = d[d$unit %in% levels(d$unit)[1:120], ], unit = "unit", se = "bootstrap", nboot = 5L, nodes = 8L)
  expect_true(all(mb$boot$replicates[, "sd:b"] > 0))
  expect_error(remnl(y ~ x1 + x2, data = d), "unit")
  ## long layout with a generic covariate
  dl <- .sim_re(300, 4, .B3, sd_u = c(0.8, 1.2), gamma = 1, seed = 5)
  ml <- remnl(chosen ~ w1 | x1 + x2, data = dl, id = "id", alt = "alt", unit = "unit", nodes = 8L)
  expect_true(ml$converged)
  zz <- (coef(ml) - c(0.3, 0.8, -0.5, -0.2, -0.6, 0.7, 1, 0.8, 1.2)) / sqrt(diag(vcov(ml))); expect_lt(max(abs(zz)), 3.5)
})
