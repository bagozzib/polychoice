## Oracle agreement: mnl() against nnet::multinom (wide layout) and mlogit
## (long layout, alternative-specific covariates, varying choice sets).

test_that("mnl() reproduces nnet::multinom: coefficients, log-likelihood, covariance, fitted", {
  skip_if_not_installed("nnet")
  set.seed(11)
  d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m <- mnl(y ~ x1 + x2, data = d)
  mm <- nnet::multinom(y ~ x1 + x2, data = d, trace = FALSE, maxit = 500, reltol = 1e-12)
  expect_true(m$converged)
  expect_equal(unname(coef(m)), as.vector(t(coef(mm))), tolerance = 1e-5)
  expect_equal(m$loglik, as.numeric(logLik(mm)), tolerance = 1e-8)
  expect_equal(unname(vcov(m)), unname(vcov(mm)), tolerance = 1e-4)
  expect_equal(unname(fitted(m)), unname(fitted(mm)), tolerance = 1e-5)
  expect_equal(names(coef(m)), c("b:(Intercept)", "b:x1", "b:x2", "c:(Intercept)", "c:x1", "c:x2"))
})

test_that("mnl() reproduces mlogit on TravelMode, with and without varying choice sets", {
  skip_if_not_installed("mlogit"); skip_if_not_installed("dfidx"); skip_if_not_installed("AER")
  data("TravelMode", package = "AER")
  TM <- TravelMode; TM$chosen <- as.integer(TM$choice == "yes")
  cl <- mnl(chosen ~ wait + vcost + travel + gcost | income + size, data = TM, id = "individual", alt = "mode")
  TMd <- dfidx::dfidx(TM, idx = c("individual", "mode"), choice = "choice")
  ml <- mlogit::mlogit(choice ~ wait + vcost + travel + gcost | income + size, data = TMd)
  map <- vapply(names(coef(cl)), function(nm) if (grepl(":", nm)) { p <- strsplit(nm, ":", fixed = TRUE)[[1]]; paste0(p[2], ":", p[1]) } else nm, character(1))
  expect_true(cl$converged)
  expect_equal(cl$loglik, as.numeric(logLik(ml)), tolerance = 1e-7)
  expect_equal(unname(coef(cl)), unname(coef(ml)[map]), tolerance = 1e-5)
  expect_equal(unname(sqrt(diag(vcov(cl)))), unname(sqrt(diag(vcov(ml)))[map]), tolerance = 1e-5)
  ## varying choice sets: remove some unchosen rows
  set.seed(5)
  TM2 <- TM[-sample(which(TM$choice == "no"), 150), ]
  cl2 <- mnl(chosen ~ wait + vcost + travel + gcost | income + size, data = TM2, id = "individual", alt = "mode")
  ml2 <- mlogit::mlogit(choice ~ wait + vcost + travel + gcost | income + size, data = dfidx::dfidx(TM2, idx = c("individual", "mode"), choice = "choice"))
  expect_equal(cl2$loglik, as.numeric(logLik(ml2)), tolerance = 1e-7)
  expect_equal(unname(coef(cl2)), unname(coef(ml2)[map]), tolerance = 1e-5)
  expect_true(any(cl2$design$avail == 0L))
  ## alternative-specific covariates only, no constants
  cl3 <- mnl(chosen ~ 0 + wait + vcost, data = TM, id = "individual", alt = "mode")
  ml3 <- mlogit::mlogit(choice ~ 0 + wait + vcost, data = TMd)
  expect_equal(cl3$loglik, as.numeric(logLik(ml3)), tolerance = 1e-7)
  expect_equal(names(coef(cl3)), c("wait", "vcost"))
})
