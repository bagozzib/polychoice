## Regression tests for the defects found by the simulated-referee panel
## (round 1): prediction with offsets, hold-out classification, the inflated
## model with an unavailable inflated alternative, the fixed-effects base rule,
## missing values in the long layout, character covariates, term bookkeeping,
## frequency weights in the tests and bootstraps, seeds inside boot_pvalue(),
## simulate() for random intercepts, intervals for positive parameters, and
## the guards the reviewers asked for.

test_that("offsets travel into predict(newdata), first_difference(), and classification()", {
  set.seed(1)
  d <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  off <- cbind(0, 0.4 * d$x1, -0.3 * d$x2)
  m <- mnl(y ~ x1 + x2, data = d, offset = off)
  expect_true(m$design$has_offset)
  expect_error(predict(m, newdata = d), "offset")
  expect_equal(unname(predict(m, newdata = d, offset = off)), unname(fitted(m)), tolerance = 1e-12)
  ## average first difference with the offset equals the hand computation
  fd <- first_difference(m, "x1", 0, 1, average = TRUE, ci = "none")
  d0 <- d; d0$x1 <- 0; d1 <- d; d1$x1 <- 1
  hand <- colMeans(predict(m, newdata = d1, offset = off)) - colMeans(predict(m, newdata = d0, offset = off))
  expect_equal(fd$diff, unname(hand), tolerance = 1e-10)
  expect_equal(classification(m, newdata = d, offset = off)$correct, classification(m)$correct)
  ## long layout with a numeric offset vector
  dl <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), gamma = 0.7)
  ov <- 0.2 * dl$w1
  ml <- mnl(chosen ~ w1 | x1 + x2, data = dl, id = "id", alt = "alt", offset = ov)
  expect_error(predict(ml, newdata = dl), "offset")
  expect_equal(unname(predict(ml, newdata = dl, offset = ov)), unname(fitted(ml)), tolerance = 1e-12)
  ## a wide offset with the wrong number of columns is refused
  expect_error(mnl(y ~ x1 + x2, data = d, offset = off[, 1:2]), "one column per alternative")
})

test_that("classification(newdata) reads the observed alternatives off the hold-out set", {
  set.seed(2)
  d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m <- mnl(y ~ x1 + x2, data = d[1:400, ])
  hold <- d[401:600, ]
  cl <- classification(m, newdata = hold)
  expect_equal(sum(cl$table), 200)
  expect_equal(rowSums(cl$table), unname(c(table(hold$y))), ignore_attr = TRUE)
  ## a hold-out set missing one alternative keeps the right labels
  hold2 <- hold[hold$y != "b", ]
  cl2 <- classification(m, newdata = hold2)
  expect_equal(unname(rowSums(cl2$table)["b"]), 0)
  expect_equal(unname(rowSums(cl2$table)["c"]), sum(hold2$y == "c"))
  expect_equal(unname(classification(m, newdata = hold, weights = rep(2, 200))$table), 2 * unname(cl$table))
})

test_that("an inflated alternative unavailable to a case gets no inflation mass", {
  set.seed(3)
  dl <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), delta = c(0.4, 0.9), inflate = "b", gamma = 0.6)
  ## drop alternative b from the choice set of 150 cases (none of them chose b)
  bcases <- unique(dl$id[dl$alt == "b" & dl$chosen == 1])
  drop_ids <- setdiff(unique(dl$id), bcases)[1:150]
  dl2 <- dl[!(dl$id %in% drop_ids & dl$alt == "b"), ]
  m <- imnl(chosen ~ w1 | x1 + x2 | z1, data = dl2, id = "id", alt = "alt", inflate = "b")
  P <- predict(m); phi <- predict(m, type = "regime")
  i_na <- rownames(P) %in% as.character(drop_ids)
  expect_true(all(P[i_na, "b"] == 0))
  expect_equal(rowSums(P), rep(1, nrow(P)), ignore_attr = TRUE, tolerance = 1e-12)
  expect_true(all(phi[i_na] == 1))
  ## those cases contribute the plain MNL probability: the derivatives still match numDeriv
  D <- m$design; fn <- polychoice:::.pc_fns(D)
  expect_equal(fn$grad(m$theta), numDeriv::grad(fn$nll, m$theta), tolerance = 1e-5)
  expect_equal(fn$hess(m$theta), numDeriv::jacobian(fn$grad, m$theta), tolerance = 1e-4)
  post <- predict(m, type = "posterior")
  expect_true(all(post >= 0 & post <= 1)); expect_true(all(post[D$y != D$k] == 1))
})

test_that("fixed effects: a unit that never chooses the base has the base unavailable and one constant fewer; constant choosers drop out", {
  set.seed(4)
  G <- 40; T <- 12; unit <- rep(seq_len(G), each = T); x1 <- rnorm(G * T)
  V <- cbind(0, 0.3 + 0.8 * x1 + rnorm(G)[unit], -0.2 - 0.6 * x1 + rnorm(G)[unit])
  P <- exp(V) / rowSums(exp(V)); y <- apply(P, 1, function(p) sample.int(3, 1, prob = p))
  y[unit == 1] <- rep(2:3, 6)            # unit 1 never chooses the base but chooses both b and c
  y[unit == 2] <- 2L                      # unit 2 always chooses b
  y[unit == 3] <- rep(1:3, 4)             # unit 3 chooses everything: the first such unit anchors the constants
  y[unit == 4] <- rep(1:3, 4)             # unit 4 chooses everything and is not the anchor
  d <- data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1, unit = factor(unit))
  expect_message(m <- suppressWarnings(mnl(y ~ x1, data = d, fe = "unit")), "never choose the base")
  D <- m$design; u <- as.character(D$unit)
  expect_true("1" %in% D$fe_info$no_base_units)
  expect_true(all(D$avail[u == "1", 1] == 0L))
  expect_equal(sum(grepl(":unit1$", D$par_names)), 1L)                               # one of unit 1's two constants is omitted
  expect_equal(sum(grepl(":unit3$", D$par_names)), 0L)                               # the anchoring unit carries no constants
  expect_equal(sum(grepl(":unit4$", D$par_names)), 2L)                               # an ordinary unit keeps both
  expect_false(any(u == "2"))                                                        # a constant chooser drops out
  expect_equal(D$fe_info$cases_dropped, T)
  expect_equal(unname(fitted(m)[u == "1", "a"]), rep(0, sum(u == "1")))
  expect_true(m$converged)
  expect_false(any(abs(m$coefficients[m$blocks == "fixed effect"]) > 15))
  expect_equal(unname(predict(m, newdata = d[d$unit == "1", ])), unname(fitted(m)[u == "1", ]), tolerance = 1e-10)
  ## the scores vanish at the optimum, including the retained constant of unit 1
  fn <- polychoice:::.pc_fns(D); expect_lt(max(abs(fn$grad(m$theta))), 1e-4)
  ## the fitted b-share of unit 1 equals its observed share (the retained constant is a sufficient statistic)
  expect_equal(sum(fitted(m)[u == "1", "c"]), sum(d$y[d$unit == "1"] == "c"), tolerance = 1e-6)
})

test_that("long layout: a missing value drops the whole case; single-alternative cases drop; character alt is refused", {
  set.seed(5)
  dl <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), gamma = 0.7)
  dl$w1[dl$id == 5 & dl$alt == "c"] <- NA                                           # an unchosen row of case 5
  expect_message(m <- mnl(chosen ~ w1 | x1 + x2, data = dl, id = "id", alt = "alt"), "dropped entirely")
  expect_equal(m$n, 299L); expect_false("5" %in% m$case_id)
  dl2 <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), gamma = 0.7)
  keep <- !(dl2$id %in% 1:10 & dl2$alt != "a" & !(dl2$chosen == 1))                # leave cases 1:10 with the chosen row only
  dl2 <- dl2[keep, ]
  expect_message(m2 <- mnl(chosen ~ w1 | x1 + x2, data = dl2, id = "id", alt = "alt"), "single available alternative")
  expect_equal(m2$n, 300L - sum(table(dl2$id) == 1L))
  dl3 <- dl2; dl3$alt <- as.character(dl3$alt)
  expect_error(mnl(chosen ~ w1 | x1 + x2, data = dl3, id = "id", alt = "alt"), "character vector")
})

test_that("character covariates are recorded like factors; a new value is refused at prediction", {
  set.seed(6)
  d <- rmnl(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d$g <- sample(c("u", "v", "w"), nrow(d), TRUE)
  mc <- mnl(y ~ x1 + g, data = d); mf <- mnl(y ~ x1 + g, data = transform(d, g = factor(g)))
  expect_equal(coef(mc), coef(mf), tolerance = 1e-8)
  expect_equal(unname(predict(mc, newdata = d[d$g != "w", ])), unname(fitted(mc)[d$g != "w", ]), tolerance = 1e-10)
  expect_error(predict(mc, newdata = transform(d[1:3, ], g = "x")), "factor .* has new level|new level")
  fd <- first_difference(mc, "g", "u", "w"); expect_equal(nrow(fd), 3L); expect_true(all(is.finite(fd$diff)))
})

test_that("term_test() keeps w1 and w10 apart; combine_test(type = 'lr') refuses varying choice sets", {
  set.seed(7)
  dl <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), gamma = c(0.7, -0.4))
  names(dl)[names(dl) == "w2"] <- "w10"
  m <- mnl(chosen ~ w1 + w10 | x1 + x2, data = dl, id = "id", alt = "alt")
  tt <- term_test(m)
  expect_equal(tt$df[tt$term == "w1"], 1L); expect_equal(tt$df[tt$term == "w10"], 1L)
  dl2 <- dl[!(dl$id %in% 1:20 & dl$alt == "c" & dl$chosen == 0), ]
  m2 <- mnl(chosen ~ 1 | x1 + x2, data = dl2, id = "id", alt = "alt")
  expect_error(combine_test(m2, c("b", "c"), type = "lr"), "available to every case")
})

test_that("frequency weights: Fok-Paap Q1 and Q2 equal the duplicated-rows computation; the weighted bootstrap runs", {
  set.seed(8)
  d <- rmnl(400, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
  w <- sample(1:3, 400, TRUE)
  mw <- mnl(y ~ x1 + x2, data = d, weights = w); md <- mnl(y ~ x1 + x2, data = d[rep(seq_len(400), w), ])
  fw <- fp_tests(mw); fd <- fp_tests(md)
  expect_equal(fw$statistic, fd$statistic, tolerance = 1e-6)
  nw <- nested_lm_test(mw, list(n1 = c("a", "b"), n2 = c("c", "d"))); nd <- nested_lm_test(md, list(n1 = c("a", "b"), n2 = c("c", "d")))
  expect_equal(nw$statistic, nd$statistic, tolerance = 1e-8)
  mb <- mnl(y ~ x1 + x2, data = d, weights = w, se = "bootstrap", nboot = 12)
  expect_true(mb$boot$weighted); expect_equal(dim(mb$boot$replicates), c(12L, 9L))
  expect_true(all(is.finite(sqrt(diag(vcov(mb))))))
})

test_that("boot_pvalue() with a reseeding test draws independent replicates and expands integer weights", {
  set.seed(9)
  d <- rmnl(400, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
  m <- mnl(y ~ x1 + x2, data = d)
  b <- boot_pvalue(m, small_hsiao, B = 15, drop = "b", seed = 1)
  reps <- attr(b, "boot_statistics")
  expect_gt(length(unique(round(reps[1, ], 8))), 8L)                                # replicates differ although the test reseeds
  mw <- mnl(y ~ x1 + x2, data = d, weights = rep(2, 400))
  bw <- boot_pvalue(mw, hausman_mcfadden, B = 10, drop = "b", vcov = "joint")
  expect_true(bw$p.boot >= 0 && bw$p.boot <= 1)
  expect_error(boot_pvalue(mnl(y ~ x1 + x2, data = d, weights = runif(400, 0.5, 2)), hausman_mcfadden, B = 10, drop = "b"), "integer frequency weights")
})

test_that("simulate() for random intercepts reproduces within-unit dependence; positive parameters get log-scale intervals", {
  set.seed(10)
  G <- 150; T <- 6; unit <- rep(seq_len(G), each = T); x1 <- rnorm(G * T); u <- 1.5 * rnorm(G)
  V <- cbind(0, 0.3 + 0.8 * x1 + u[unit], -0.2 - 0.6 * x1 + 1.5 * rnorm(G)[unit])
  P <- exp(V) / rowSums(exp(V)); y <- apply(P, 1, function(p) sample.int(3, 1, prob = p))
  d <- data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1, unit = factor(unit))
  m <- remnl(y ~ x1, data = d, unit = "unit", nodes = 8L)
  s <- simulate(m, nsim = 20, seed = 1)
  ## the share of b within a unit varies across units far more than under independence
  within_var <- mean(sapply(s, function(v) var(tapply(v == "b", d$unit, mean))))
  pb <- mean(fitted(m)[, "b"]); indep_var <- pb * (1 - pb) / T
  expect_gt(within_var, 2 * indep_var)
  ci <- confint(m)
  expect_true(all(ci[grepl("^sd:", rownames(ci)), 1] > 0))
  s2 <- summary(m); expect_true(all(is.na(s2$coefficients[s2$block == "random intercept", "z value"])))
  expect_error(hl_test(m), "independent cases")
})

test_that("guards: fe fits are refused by the IIA tests, lr_test() checks nesting by name, the probit class is mnprobit", {
  set.seed(11)
  G <- 80; T <- 8; unit <- rep(seq_len(G), each = T); x1 <- rnorm(G * T)
  V <- cbind(0, 0.3 + 0.8 * x1 + rnorm(G)[unit], -0.2 - 0.6 * x1 + rnorm(G)[unit])
  P <- exp(V) / rowSums(exp(V)); y <- apply(P, 1, function(p) sample.int(3, 1, prob = p))
  d <- data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1, x2 = rnorm(G * T), unit = factor(unit))
  mfe <- suppressWarnings(suppressMessages(mnl(y ~ x1 + x2, data = d, fe = "unit")))
  expect_error(hausman_mcfadden(mfe), "fixed-effects"); expect_error(small_hsiao(mfe), "fixed-effects")
  expect_error(boot_pvalue(mfe, hausman_mcfadden, B = 10), "fixed-effects")
  expect_error(lr_test(mnl(y ~ x2, data = d), mnl(y ~ x1, data = d)), "not nested")
  expect_error(lr_test(mnl(y ~ x1, data = d), mnl(y ~ x1 + x2, data = d, weights = rep(2, nrow(d)))), "same weights")
  dp <- rmnp(500, beta = cbind(b = c(0.3, -0.6, 0.5), c = c(-0.2, 0.4, -0.3)), sigma = matrix(c(1, 0.5, 0.5, 1), 2), gamma = 1)
  mp <- mnprobit(chosen ~ w1 | x1 + x2, data = dp, id = "id", alt = "alt")
  expect_s3_class(mp, "mnprobit"); expect_false(inherits(mp, "mnp"))
  expect_output(print(mp), "Alternatives: a, b, c  \\(base: a\\)")
})

test_that("boot_lr_test() calibrates mnl against imnl and the nested-logit summary tests lambda against 1", {
  skip_on_cran()
  set.seed(12)
  d <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m0 <- mnl(y ~ x1 + x2, data = d); m1 <- suppressMessages(imnl(y ~ x1 + x2 | x1, data = d))
  b <- boot_lr_test(m0, m1, B = 12)
  expect_equal(b$test[1], "LR (parametric bootstrap under the restricted model)")
  expect_equal(nrow(b), 2L); expect_true(all(b$p.boot >= 0 & b$p.boot <= 1)); expect_equal(ncol(attr(b, "boot_statistics")), attr(b, "B"))
  nests <- list(n1 = c("a", "b"), n2 = "c")
  dn <- rnlogit(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), nests = nests, lambda = c(0.5, 1), gamma = 1)
  mn <- suppressMessages(nlogit(chosen ~ w1 | x1 + x2, data = dn, id = "id", alt = "alt", nests = nests))
  s <- summary(mn); r <- s$coefficients["lambda:n1", ]
  expect_equal(unname(r["z value"]), unname((r["Estimate"] - 1) / r["Std. Error"]), tolerance = 1e-12)
})
