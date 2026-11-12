## likelihood-ratio statistics from suprema under separation, and the
## parametric bootstrap for rank-ordered fits (DESIGN.md, section 1, Phase 4c)

## eight rankings of three alternatives ordered exactly by 3 x1 - 2 x2, with
## two covariates of noise: completely separated, and still completely
## separated without x3 or without x4
sup_rank_data <- function() {
  set.seed(11)
  n <- 8L; J <- 3L
  d <- data.frame(case = rep(seq_len(n), each = J), alt = rep(seq_len(J), n),
                  x1 = round(rnorm(n * J), 3), x2 = round(rnorm(n * J), 3), x3 = round(rnorm(n * J), 3), x4 = round(rnorm(n * J), 3))
  d$rank <- ave(-(3 * d$x1 - 2 * d$x2), d$case, FUN = rank)
  d
}

test_that("suprema: complete separation gives 0, and term_test() compares suprema (the Poast pattern)", {
  d <- sup_rank_data()
  m <- suppressWarnings(rologit(rank ~ 0 + x1 + x2 + x3 + x4, data = d, id = "case", alt = "alt"))
  expect_identical(separation(m)$type, "complete")
  expect_lt(m$loglik, 0)                                            # where the optimizer stopped
  sc <- .pc_sup_loglik(m$design, m$theta, m$loglik)
  expect_identical(sc[c("value", "type")], list(value = 0, type = "complete"))
  ## the point returned: the fit moved along the direction until every
  ## comparison's utility difference is at least 40, where the log-likelihood is near 0
  expect_gt(-polychoice:::.pc_fns(m$design)$nll(sc$theta), -1e-12)
  expect_no_warning(tt <- term_test(m, type = "lr"))
  for (k in c("x3", "x4")) {
    r <- suppressWarnings(rologit(reformulate(c("0", setdiff(paste0("x", 1:4), k)), "rank"), data = d, id = "case", alt = "alt"))
    expect_identical(separation(r)$type, "complete")
    expect_identical(tt$statistic[tt$term == k], 0)                 # both suprema are 0
    expect_identical(tt$p.value[tt$term == k], 1)
  }
  expect_match(attr(tt, "method"), "suprema under separation")
  ## dropping x1 removes the separation: the statistic is minus twice the restricted maximum
  r1 <- rologit(rank ~ 0 + x2 + x3 + x4, data = d, id = "case", alt = "alt")
  expect_identical(separation(r1)$type, "overlap")
  expect_equal(tt$statistic[tt$term == "x1"], -2 * r1$loglik, tolerance = 1e-7)
  ## lr_test() between two separated fits
  r3 <- suppressWarnings(rologit(rank ~ 0 + x1 + x2 + x4, data = d, id = "case", alt = "alt"))
  lt <- lr_test(r3, m)
  expect_identical(lt$statistic, 0); expect_true(isTRUE(attr(lt, "suprema")))
  expect_output(print(lt), "suprema under separation")
})

test_that("suprema: quasi-complete separation, against the fit of the data without the separated cases", {
  set.seed(1)
  n <- 300; x <- rbinom(n, 1, 0.2); z <- rnorm(n)
  y <- ifelse(x == 1, "c", sample(c("a", "b", "c"), n, TRUE))
  dq <- data.frame(y = factor(y), x = x, z = z)
  mq <- suppressWarnings(mnl(y ~ x + z, data = dq))
  expect_identical(separation(mq)$type, "quasi-complete")
  s <- .pc_sup_loglik(mq$design, mq$theta, mq$loglik)
  ## along the separating direction the x = 1 cases' probabilities tend to one,
  ## and x leaves the rest: the supremum is the x = 0 cases' own maximum
  mh <- mnl(y ~ z, data = dq[dq$x == 0, ])
  expect_identical(s$type, "quasi-complete")
  expect_equal(s$value, mh$loglik, tolerance = 1e-9)
  expect_gte(s$value, mq$loglik - 1e-12)
  ## term_test() and lr_test() use it; on overlapping data nothing changes
  mq0 <- suppressWarnings(mnl(y ~ z, data = dq))
  lt <- lr_test(mq0, mq)
  expect_equal(lt$statistic, 2 * (s$value - .pc_sup_loglik(mq0$design, mq0$theta, mq0$loglik)$value), tolerance = 1e-9)
  set.seed(2); d <- rmnl(400, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m1 <- mnl(y ~ x1 + x2, data = d); m0 <- mnl(y ~ x1, data = d)
  expect_identical(.pc_sup_loglik(m1$design, m1$theta, m1$loglik), list(value = m1$loglik, type = "overlap", theta = m1$theta))
  expect_equal(lr_test(m0, m1)$statistic, 2 * (m1$loglik - m0$loglik))
  expect_null(attr(lr_test(m0, m1), "suprema"))
})

test_that("suprema: out of scope (inflated or fixed-effects designs) the values stay as they are", {
  set.seed(3); d <- rmnl(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  mi <- suppressMessages(suppressWarnings(imnl(y ~ x1 | 1, data = d)))
  expect_identical(.pc_sup_loglik(mi$design, mi$theta, mi$loglik), list(value = mi$loglik, type = NA_character_, theta = mi$theta))
})

test_that("boot_lr_test(): rank-ordered fits, the observed rank groups kept in every replicate", {
  ## the coarsening: each case keeps its group sizes in order. Case 1, groups
  ## (1, 2): the best alternative ranked 1, the next two tied at 2; case 2,
  ## one group of 1 among the three available; case 3, groups (2, 1, 1): the
  ## two best tied at 1, then 2, then 3
  U <- rbind(c(3, 1, 2, 0), c(0, 5, 4, 1), c(1, 2, 3, 4))
  avail <- rbind(c(1L, 1L, 1L, 1L), c(1L, 1L, 1L, 0L), c(1L, 1L, 1L, 1L))
  R <- .pc_ro_coarsen(U, avail, list(c(1L, 2L), 1L, c(2L, 1L, 1L)))
  expect_identical(R, rbind(c(1L, 2L, 2L, 0L), c(0L, 1L, 0L, 0L), c(3L, 2L, 1L, 1L)))
  set.seed(4)
  dr <- rrologit(150, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6), d = c(0.2, 0.4)), gamma = 0.7, depth = 2, ties = 0.2)
  r0 <- rologit(rank ~ w1 | 1, data = dr, id = "id", alt = "alt")
  r1 <- rologit(rank ~ w1 | x1, data = dr, id = "id", alt = "alt")
  set.seed(5); bt <- boot_lr_test(r0, r1, B = 19)
  expect_true(bt$p.boot > 0 && bt$p.boot <= 1)
  expect_equal(bt$statistic, lr_test(r0, r1)$statistic, tolerance = 1e-6)
  expect_match(attr(bt, "method"), "observed rank groups")
  set.seed(5); bt2 <- boot_lr_test(r0, r1, B = 19)
  expect_identical(bt$p.boot, bt2$p.boot)
  expect_error(boot_lr_test(rologit(rank ~ w1 | 1, data = dr, id = "id", alt = "alt", firth = TRUE),
                            rologit(rank ~ w1 | x1, data = dr, id = "id", alt = "alt", firth = TRUE), B = 19), "penalized")
})

test_that("boot_lr_test(): separated rankings, the observed statistic and the replicates from suprema", {
  d <- sup_rank_data()
  r <- suppressWarnings(rologit(rank ~ 0 + x2 + x3 + x4, data = d, id = "case", alt = "alt"))
  m <- suppressWarnings(rologit(rank ~ 0 + x1 + x2 + x3 + x4, data = d, id = "case", alt = "alt"))
  set.seed(6); bt <- suppressMessages(suppressWarnings(boot_lr_test(r, m, B = 19)))
  expect_equal(bt$statistic, -2 * r$loglik, tolerance = 1e-7)       # the full model's supremum is 0
  expect_match(attr(bt, "method"), "the observed statistic uses log-likelihood suprema")
  expect_true(attr(bt, "n_separated") >= 1L)
})

## DESIGN.md Phase 4d: comparisons with identical design rows
test_that("separation: a comparison with a zero structure vector makes the data at most quasi-complete (round 5, B1/D1/I1)", {
  ## B1: case 3's two alternatives share x = 1; the others are separated by x
  d <- data.frame(id = rep(1:3, each = 2), alt = factor(rep(c("p", "q"), 3)), chosen = c(1, 0, 1, 0, 1, 0), x = c(1, 0, 1, 0, 1, 1))
  m <- suppressWarnings(mnl(chosen ~ 0 + x, d, id = "id", alt = "alt"))
  s <- separation(m)
  expect_identical(s$type, "quasi-complete")
  expect_identical(c(s$n_comparisons, s$n_separated, s$n_zero), c(3L, 2L, 1L))
  expect_output(print(s), "identical design rows")
  expect_equal(.pc_sup_loglik(m$design, m$theta, m$loglik)$value, log(1 / 2), tolerance = 1e-9)
  ## D1: J = 3, case 1 ties x on its alternatives a and b; z separates it. The
  ## full model is completely separated (supremum 0); the restricted one keeps
  ## the tie (supremum log 1/2), so the LR of z is 2 log 2
  d2 <- data.frame(id = rep(1:4, each = 3), alt = factor(rep(c("a", "b", "c"), 4)), chosen = c(1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0),
                   x = c(1, 1, 0, 0, 2, 1, 0, 1, 2, 2, 1, 0), z = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
  mf <- suppressWarnings(mnl(chosen ~ 0 + x + z, d2, id = "id", alt = "alt"))
  mr <- suppressWarnings(mnl(chosen ~ 0 + x, d2, id = "id", alt = "alt"))
  expect_identical(separation(mf)$type, "complete"); expect_identical(separation(mr)$type, "quasi-complete")
  expect_equal(.pc_sup_loglik(mr$design, mr$theta, mr$loglik)$value, log(1 / 2), tolerance = 1e-9)
  tt <- term_test(mf, type = "lr")
  expect_equal(tt$statistic[tt$term == "z"], 2 * log(2), tolerance = 1e-8)
  expect_equal(lr_test(mr, mf)$statistic, 2 * log(2), tolerance = 1e-8)
  ## I1: a generic-only ranking, z binary; case 1 ranks its two z = 1
  ## alternatives first and second (a zero comparison), and every other
  ## comparison puts z = 1 above z = 0 (cases 2 and 3 rank their z = 1
  ## alternative first, the rest unranked); the supremum is log 1/2 (case 1's
  ## first stage)
  dr <- data.frame(case = rep(1:3, each = 3), alt = factor(rep(c("a", "b", "c"), 3)), rank = c(1, 2, 3, 0, 1, 0, 0, 0, 1),
                   z = c(1, 1, 0, 0, 1, 0, 0, 0, 1))
  rz <- suppressWarnings(rologit(rank ~ 0 + z, data = dr, id = "case", alt = "alt"))
  expect_identical(separation(rz)$type, "quasi-complete")
  expect_equal(.pc_sup_loglik(rz$design, rz$theta, rz$loglik)$value, log(1 / 2), tolerance = 1e-8)
})

test_that("separation: rows equal but for rounding are zero comparisons (round 6, B-N2 / D-new-3 / I-4)", {
  ## B-N2: case 1's two values differ by 5.6e-17 (0.1 + 0.2 against 0.3); the
  ## other two cases are separated by x
  d <- data.frame(id = rep(1:3, each = 2), alt = factor(rep(c("p", "q"), 3)), chosen = c(1, 0, 1, 0, 1, 0), x = c(0.1 + 0.2, 0.3, 0, 1, 0, 1))
  m <- suppressWarnings(mnl(chosen ~ 0 + x, d, id = "id", alt = "alt"))
  s <- separation(m)
  expect_identical(s$type, "quasi-complete"); expect_identical(s$n_zero, 1L)
  ## D-new-3: the round-5 D1 data with case 1's tie computed in floating point,
  ## in either order: the supremum of the restricted model stays log 1/2 and
  ## the LR of z 2 log 2
  for (x1 in list(c(0.1 + 0.2, 0.3, 0), c(0.3, 0.1 + 0.2, 0))) {
    d2 <- data.frame(id = rep(1:4, each = 3), alt = factor(rep(c("a", "b", "c"), 4)), chosen = c(1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0),
                     x = c(x1, 0, 2, 1, 0, 1, 2, 2, 1, 0), z = c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0))
    mr <- suppressWarnings(mnl(chosen ~ 0 + x, d2, id = "id", alt = "alt"))
    mf <- suppressWarnings(mnl(chosen ~ 0 + x + z, d2, id = "id", alt = "alt"))
    expect_identical(separation(mr)$type, "quasi-complete")
    expect_equal(.pc_sup_loglik(mr$design, mr$theta, mr$loglik)$value, log(1 / 2), tolerance = 1e-8)
    expect_equal(lr_test(mr, mf)$statistic, 2 * log(2), tolerance = 1e-7)
  }
  ## round 7, D-1: the rule is relative to the operands, not to the column's
  ## spread. (a) A large level whose tie differs in the last digit is a zero
  ## comparison: x shifted by 2e6, case 1's alternatives a and b one ulp apart
  d3 <- data.frame(id = rep(1:4, each = 3), alt = factor(rep(c("a", "b", "c"), 4)), chosen = c(1, 0, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0),
                   x = c(1, 1, 0, 0, 2, 1, 0, 1, 2, 2, 1, 0) + 2e6)
  d3$x[2] <- d3$x[1] + 2^-32
  expect_gt(d3$x[2], d3$x[1])
  mr3 <- suppressWarnings(mnl(chosen ~ 0 + x, d3, id = "id", alt = "alt"))
  s3 <- separation(mr3)
  expect_identical(s3$type, "quasi-complete"); expect_identical(s3$n_zero, 1L)
  ## (b) a genuine difference stays, however wide the column: comparisons +1,
  ## -1, +1e10 and +2 overlap (the round-6 rule set the first two aside, as
  ## below 1e-10 of the column's largest, and reported a false separation)
  d4 <- data.frame(id = rep(1:4, each = 2), alt = factor(rep(c("p", "q"), 4)), chosen = c(1, 0, 0, 1, 1, 0, 0, 1), x = c(1, 0, 1, 0, 1e10, 0, 0, 2))
  m4 <- suppressWarnings(mnl(chosen ~ 0 + x, d4, id = "id", alt = "alt"))
  s4 <- separation(m4)
  expect_false(s4$separated); expect_identical(s4$n_zero, 0L)
})

test_that("boot_lr_test(): single-choice replicates with a never-chosen alternative are kept in scope; serial and parallel draws agree", {
  set.seed(7)
  n <- 60; x <- rnorm(n)
  y <- factor(sample(c("a", "b", "c"), n, TRUE, prob = c(0.48, 0.48, 0.04)), levels = c("a", "b", "c"))
  dd <- data.frame(y = y, x = x)
  m0 <- suppressWarnings(mnl(y ~ 1, data = dd)); m1 <- suppressWarnings(mnl(y ~ x, data = dd))
  set.seed(8); bt <- suppressMessages(suppressWarnings(boot_lr_test(m0, m1, B = 19)))
  expect_identical(attr(bt, "n_failed"), 0L)                          # replicates without a "c" are separation, kept
  expect_gte(attr(bt, "n_separated"), 1L)                             # and at least one occurred (round 7, G-N3)
  ## rankings: the replicates' draws do not depend on the number of workers
  skip_on_cran()
  set.seed(4)
  dr <- rrologit(80, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6), d = c(0.2, 0.4)), gamma = 0.7, depth = 2)
  r0 <- rologit(rank ~ w1 | 1, data = dr, id = "id", alt = "alt"); r1 <- rologit(rank ~ w1 | x1, data = dr, id = "id", alt = "alt")
  set.seed(9); b1 <- boot_lr_test(r0, r1, B = 10)
  set.seed(9); b2 <- boot_lr_test(r0, r1, B = 10, cores = 2)
  expect_identical(attr(b1, "boot_statistics"), attr(b2, "boot_statistics"))
})

test_that("boot_lr_test(): a separated restricted fit draws from the limit of the restricted model", {
  d <- sup_rank_data()
  r <- suppressWarnings(rologit(rank ~ 0 + x1 + x2 + x3, data = d, id = "case", alt = "alt"))
  m <- suppressWarnings(rologit(rank ~ 0 + x1 + x2 + x3 + x4, data = d, id = "case", alt = "alt"))
  set.seed(10); bt <- suppressMessages(suppressWarnings(boot_lr_test(r, m, B = 10)))
  expect_match(attr(bt, "method"), "separated restricted model's limit")
})

test_that("combine_test(type = \"lr\"): a never-chosen alternative contributes 0 log 0 = 0", {
  ## the long layout keeps an alternative that no case chooses (the wide layout drops an unused level)
  set.seed(12)
  n <- 120; y <- sample(c("a", "b"), n, TRUE); x <- rnorm(n)
  d <- data.frame(id = rep(seq_len(n), each = 3), alt = factor(rep(c("a", "b", "c"), n)), x = rep(x, each = 3))
  d$chosen <- as.integer(as.character(d$alt) == rep(y, each = 3))
  m <- suppressWarnings(mnl(chosen ~ 1 | x, data = d, id = "id", alt = "alt"))
  expect_no_warning(ct <- suppressMessages(combine_test(m, c("b", "c"), type = "lr")))
  expect_true(is.finite(ct$statistic))
})

test_that("term_test(type = \"lr\"): a restricted model without parameters is evaluated without warnings", {
  set.seed(13)
  d <- data.frame(id = rep(1:40, each = 3), alt = factor(rep(c("a", "b", "c"), 40)), w = rnorm(120))
  d$chosen <- as.integer(ave(d$w + rnorm(120), d$id, FUN = function(v) v == max(v)))
  m <- mnl(chosen ~ 0 + w, data = d, id = "id", alt = "alt")
  expect_no_warning(tt <- term_test(m, type = "lr"))
  expect_equal(tt$statistic, 2 * (m$loglik - sum(-log(3) * 40)), tolerance = 1e-8)
})
