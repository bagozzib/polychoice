## The local review of 2026-11-06 (internal_reviews/local_review_2026-11-06):
## one test per finding that a compact example reaches. Finding 1 (term_test's
## restricted refit with the saddle step) and 15 are tested on the Liou probit
## fixture in test-phase4d.R and checked in the ledger; 11 by a mocked failure.

test_that("combine_test() refuses the same alternative twice (finding 2)", {
  set.seed(1)
  d <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m <- mnl(y ~ x1 + x2, data = d)
  expect_error(combine_test(m, c("b", "b"), type = "lr"), "two different response levels")
  expect_error(combine_test(m, c("b", "b")), "two different response levels")
})

test_that("a tie up to rounding is a tie after centring (finding 3)", {
  set.seed(3); n <- 60; J <- 3
  L <- data.frame(id = rep(1:n, each = J), alt = factor(rep(c("a", "b", "c"), n)))
  L$cost <- runif(n * J); L$cost[1] <- 0.1 + 0.2; L$cost[2] <- 0.3; L$cc <- L$cost - 0.3
  V <- -10 * L$cost
  L$ch <- unlist(tapply(V, L$id, function(v) { z <- integer(J); z[which.max(v)] <- 1L; z })); L$ch[1:3] <- c(1L, 0L, 0L)
  s1 <- separation(suppressWarnings(mnl(ch ~ 0 + cost, data = L, id = "id", alt = "alt")))
  s2 <- separation(suppressWarnings(mnl(ch ~ 0 + cc, data = L, id = "id", alt = "alt")))
  expect_identical(s2$type, s1$type)
  expect_identical(s2$n_zero, s1$n_zero)
  expect_identical(s2$type, "quasi-complete")
})

test_that("remnl()'s bootstrap refuses weights totalling at most 2 (finding 4)", {
  set.seed(4)
  d <- rmnl(300, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6))); d$unit <- rep(1:30, each = 10); d$wt <- 1 / 300
  expect_error(suppressMessages(remnl(y ~ x1, data = d, unit = "unit", weights = "wt", se = "bootstrap", nboot = 5)), "total at most 2")
})

test_that("rank_test() leaves out coefficients a block of stages does not identify (finding 5)", {
  set.seed(5); n <- 200
  d <- rrologit(n, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6), d = c(0.2, 0.4)))
  for (i in unique(d$id)) {                       # b always in the first two positions
    ix <- which(d$id == i); b <- ix[d$alt[ix] == "b"]
    if (d$rank[b] > 2) { sw <- ix[d$rank[ix] == 2]; d$rank[sw] <- d$rank[b]; d$rank[b] <- 2 }
  }
  r0 <- suppressWarnings(rologit(rank ~ 1 | x1, data = d, id = "id", alt = "alt"))
  expect_warning(rt <- rank_test(r0, after = 2), "do not identify b:\\(Intercept\\), b:x1")
  expect_equal(rt$df, c(4L, 4L))
  expect_true(all(is.na(attr(rt, "coefficients")[c("b:(Intercept)", "b:x1"), "difference"])))
  r1 <- suppressWarnings(suppressMessages(rologit(rank ~ 1 | x1, data = d, id = "id", alt = "alt", firth = TRUE)))
  expect_error(rank_test(r1, after = 2), "penalized fit needs every coefficient identified")
})

test_that("captivity_test() leaves a never-chosen alternative out of the joint test (finding 6)", {
  set.seed(1); n <- 200; J <- 4
  dl <- data.frame(id = rep(1:n, each = J), alt = factor(rep(letters[1:J], n)), w = rnorm(n * J), x = rep(rnorm(n), each = J))
  u <- 0.7 * dl$w + ifelse(dl$alt == "d", -Inf, 0) + rlogis(n * J)
  dl$chosen <- as.integer(ave(u, dl$id, FUN = function(v) v == max(v)))
  m <- suppressWarnings(mnl(chosen ~ w | x, data = dl, id = "id", alt = "alt"))
  ct <- suppressWarnings(captivity_test(m))
  abc <- suppressWarnings(captivity_test(m, alternatives = c("a", "b", "c")))
  expect_true(is.na(ct$statistic[ct$test == "captive: d"]))
  expect_equal(ct$statistic[ct$test == "joint"], abc$statistic[abc$test == "joint"], tolerance = 1e-8)
  expect_equal(ct$df[ct$test == "joint"], 3L)
  expect_match(attr(ct, "method"), "never chosen, not tested: d")
})

test_that("a never-chosen alternative whose estimate exists does not fail a resample (finding 7)", {
  set.seed(7); n <- 300
  dl <- data.frame(id = rep(1:n, each = 3), alt = factor(rep(c("a", "b", "c"), n)), x = rep(rnorm(n), each = 3), w = rnorm(3 * n))
  u <- 0.8 * dl$w + ifelse(dl$alt == "b", 1, 0) * dl$x + ifelse(dl$alt == "c", -Inf, 0) + rlogis(3 * n)
  dl$chosen <- as.integer(ave(u, dl$id, FUN = function(v) v == max(v)))
  m <- mnl(chosen ~ 0 + w | x, data = dl, id = "id", alt = "alt")
  expect_identical(separation(m)$type, "overlap")
  set.seed(1)
  mb <- mnl(chosen ~ 0 + w | x, data = dl, id = "id", alt = "alt", se = "bootstrap", nboot = 20)
  expect_true(all(is.finite(sqrt(diag(vcov(mb))))))
})

test_that("a tiny but real captive share keeps its statistic (finding 8)", {
  cr <- polychoice:::.pc_captive_rule(7.72, list(phi_boundary = TRUE, diverging = character(0), saddle = FALSE), -500, 0L)
  expect_equal(cr$statistic, 7.72)
  expect_false(cr$degenerate)
  cr0 <- polychoice:::.pc_captive_rule(-1e-3, list(phi_boundary = TRUE, diverging = character(0), saddle = FALSE), -500, 0L)
  expect_equal(cr0$statistic, 0); expect_false(cr0$below)            # a collapsed share a little below: no warning
})

test_that("an NA replicate statistic leaves only its own row (finding 9)", {
  set.seed(9)
  d <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m <- mnl(y ~ x1 + x2, data = d)
  odd <- sum(m$y == 1L) %% 2L
  tf <- function(o) data.frame(test = c("t1", "t2"), statistic = c(if (sum(o$y == 1L) %% 2L != odd) NA_real_ else 1, 2))
  b <- boot_pvalue(m, tf, B = 30)
  expect_false(anyNA(b$p.boot))
  expect_true(attr(b, "n_na")[["t1"]] > 0L)
  expect_equal(attr(b, "n_na")[["t2"]], 0)
})

test_that("likelihood ratios of suprema do not count a never-chosen alternative's coefficients (finding 10)", {
  set.seed(1)
  d <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.2, 0.3)), gamma = c(1))
  cc <- unique(d$id[d$alt == "c" & d$chosen == 1])
  d$chosen[d$id %in% cc & d$alt == "c"] <- 0L; d$chosen[d$id %in% cc & d$alt == "b"] <- 1L
  m1 <- suppressWarnings(mnl(chosen ~ w1 | x1, data = d, id = "id", alt = "alt"))
  m2 <- suppressWarnings(mnl(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt"))
  dd <- droplevels(d[d$alt != "c", ])
  r1 <- mnl(chosen ~ w1 | x1, data = dd, id = "id", alt = "alt"); r2 <- mnl(chosen ~ w1 | x1 + x2, data = dd, id = "id", alt = "alt")
  lt <- lr_test(m1, m2); lr <- lr_test(r1, r2)
  expect_equal(lt$df, lr$df); expect_equal(lt$statistic, lr$statistic, tolerance = 1e-5)
  expect_identical(attr(lt, "absent"), "c")
  tt <- suppressWarnings(term_test(m2, type = "lr"))
  expect_equal(tt$df[tt$term == "x2"], 2L); expect_equal(tt$df[tt$term == "w1"], 1L)
})

test_that("a failed separation check in a supremum warns (finding 11)", {
  set.seed(11)
  d <- rmnl(200, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m <- mnl(y ~ x1 + x2, data = d)
  local_mocked_bindings(.pc_separation = function(...) stop("lp failed"), .package = "polychoice")
  expect_warning(s <- polychoice:::.pc_sup_loglik(m$design, m$theta, m$loglik), "separation check failed")
  expect_equal(s$value, m$loglik)
})

test_that("profile intervals are equivariant under rescaling a covariate (finding 13)", {
  set.seed(13)
  d <- rmnl(400, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5))); d$x1b <- d$x1 * 1e6
  f1 <- suppressMessages(mnl(y ~ x1 + x2, data = d, firth = TRUE)); f2 <- suppressMessages(mnl(y ~ x1b + x2, data = d, firth = TRUE))
  expect_no_warning(c2 <- confint(f2, "b:x1b"))
  expect_equal(unname(c2 * 1e6), unname(confint(f1, "b:x1")), tolerance = 1e-6)
})

test_that("firth = TRUE refuses weights totalling at most 2 (finding 14)", {
  set.seed(14)
  d <- rmnl(200, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5))); d$wt <- 1 / 200
  expect_error(suppressMessages(mnl(y ~ x1 + x2, data = d, weights = "wt", firth = TRUE)), "total above 2")
})

test_that("two t priors with df 7L and 7 are the same prior (review, cut items)", {
  set.seed(15)
  d <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  a <- mnl(y ~ x1, data = d, prior = coef_prior("t", df = 7L)); b <- mnl(y ~ x1 + x2, data = d, prior = coef_prior("t", df = 7))
  expect_s3_class(lr_test(a, b), "lr_test")
})

test_that("nested_lm_test() takes a never-chosen alternative out (review, cut items)", {
  set.seed(16); n <- 300; J <- 4
  dl <- data.frame(id = rep(1:n, each = J), alt = factor(rep(letters[1:J], n)), w = rnorm(n * J), x1 = rep(rnorm(n), each = J))
  u <- 0.7 * dl$w + ifelse(dl$alt == "d", -Inf, 0) + rlogis(n * J)
  dl$chosen <- as.integer(ave(u, dl$id, FUN = function(v) v == max(v)))
  m <- suppressWarnings(mnl(chosen ~ w | x1, data = dl, id = "id", alt = "alt"))
  nl <- suppressWarnings(nested_lm_test(m, nests = list(n1 = c("a", "b"), n2 = c("c", "d"))))
  expect_equal(nl$df, 1L)
  expect_match(attr(nl, "method"), "not tested: lambda:n2")
  expect_gt(nl$p.value, 0.01)                     # IIA holds in these data
})

test_that("hl_test() counts the groups that received cases (review, cut items)", {
  set.seed(17)
  d <- rmnl(200, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5))); d$wt <- 1; d$wt[1] <- 300
  h <- hl_test(mnl(y ~ x1 + x2, data = d, weights = "wt"))
  g <- sum(rowSums(h$expected) > 0)
  expect_lt(g, 10L)
  expect_equal(h$df, (g - 2L) * 2L)
})

test_that("Firth's penalty runs with more than 64 alternatives (review, cut items)", {
  set.seed(18); n <- 60; J <- 70
  d <- data.frame(id = rep(1:n, each = J), alt = factor(rep(sprintf("a%02d", 1:J), n)), w = rnorm(n * J))
  u <- 0.8 * d$w + rlogis(n * J); d$chosen <- as.integer(ave(u, d$id, FUN = function(v) v == max(v)))
  f <- suppressMessages(mnl(chosen ~ 0 + w, data = d, id = "id", alt = "alt", firth = TRUE))
  expect_true(is.finite(coef(f)[["w"]]))
})
