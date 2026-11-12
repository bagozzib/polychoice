## The rank-ordered logit: the exploded-data and first-choice identities, the
## exact tie likelihood against a brute-force sum, derivatives, the J = 2 and
## mlogit oracles, the stage-stability test, and the methods and refusals.

.Bro <- cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3))

## one pseudo-case per stage, the alternatives not yet ranked available (no ties)
.explode <- function(d) {
  do.call(rbind, lapply(split(d, d$id), function(di) {
    di <- di[di$rank > 0 | TRUE, ]
    di <- di[order(ifelse(di$rank > 0, di$rank, Inf)), ]
    L <- sum(di$rank > 0); J <- nrow(di)
    do.call(rbind, lapply(seq_len(min(L, J - 1L)), function(l) {
      s <- di[l:J, ]; s$chosen <- as.integer(seq_len(nrow(s)) == 1L); s$pid <- paste(di$id[1L], l); s
    }))
  }))
}

test_that("rologit(): the exploded data fitted by mnl() and depth = 1 against mnl() on the first choices", {
  set.seed(11)
  d <- rrologit(300, beta = .Bro, gamma = 0.7, depth = 3)
  m <- rologit(rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
  ex <- .explode(d); ex$alt <- factor(ex$alt, levels = levels(d$alt))
  me <- mnl(chosen ~ w1 | x1 + x2, data = ex, id = "pid", alt = "alt")
  expect_equal(unname(coef(me)), unname(coef(m)), tolerance = 1e-8)
  expect_equal(me$loglik, m$loglik, tolerance = 1e-10)
  expect_equal(sqrt(diag(vcov(me))), sqrt(diag(vcov(m))), tolerance = 1e-6, ignore_attr = TRUE)
  m1 <- rologit(rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt", depth = 1)
  d1 <- d; d1$chosen <- as.integer(d1$rank == 1)
  mf <- mnl(chosen ~ w1 | x1 + x2, data = d1, id = "id", alt = "alt")
  expect_equal(unname(coef(mf)), unname(coef(m1)), tolerance = 1e-8)
  expect_equal(mf$loglik, m1$loglik, tolerance = 1e-10)
  expect_s3_class(m, "rologit")
  expect_match(polychoice:::.pc_title(m), "Rank-ordered logit")
})

test_that("rologit(): exact tie likelihood against the sum over orders; derivatives against numDeriv", {
  set.seed(12)
  d <- suppressMessages(rrologit(250, beta = .Bro, gamma = 0.7, ties = 0.35, depth = 3))
  m <- suppressMessages(rologit(rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt"))
  fn <- polychoice:::.pc_fns(m$design)
  set.seed(1); th <- m$theta + rnorm(length(m$theta), 0, 0.1)
  R <- m$design$str$rank; A <- m$design$avail; V <- fn$probs(th)$V
  perms <- function(v) if (length(v) <= 1L) list(v) else do.call(c, lapply(seq_along(v), function(k) lapply(perms(v[-k]), function(p) c(v[k], p))))
  brute <- function(i) {
    r <- R[i, ]; unr <- which(r == 0L & A[i, ] == 1L); ll <- 0
    for (g in sort(unique(r[r > 0L]))) {
      T <- which(r == g); below <- c(which(r > g), unr)
      if (!length(below)) next
      ll <- ll + log(sum(vapply(perms(T), function(p) prod(vapply(seq_along(p), function(s)
        exp(V[i, p[s]]) / (sum(exp(V[i, p[s:length(p)]])) + sum(exp(V[i, below]))), numeric(1))), numeric(1))))
    }
    ll
  }
  tied <- which(apply(R, 1L, function(r) any(duplicated(r[r > 0L]))))
  expect_gt(length(tied), 10L)
  lli <- fn$lli(th)
  for (i in c(tied[1:5], 1:3)) expect_equal(lli[i], brute(i), tolerance = 1e-12)
  g_num <- numDeriv::grad(function(t) sum(fn$lli(t)), th)
  expect_equal(-fn$grad(th), g_num, tolerance = 1e-6)
  H_num <- numDeriv::jacobian(function(t) -fn$grad(t), th)
  expect_equal(-fn$hess(th), H_num, tolerance = 1e-6)
  expect_equal(colSums(fn$scores(th)), -fn$grad(th), tolerance = 1e-10)
  ## a ranking whose alternatives are all tied carries nothing: such cases drop
  dd <- d[d$id %in% unique(d$id)[1:40], ]
  dd$rank[dd$id == dd$id[1]] <- 1
  expect_message(rologit(rank ~ w1 | x1 + x2, data = dd, id = "id", alt = "alt"), "carries no information")
})

test_that("rologit(): J = 2 is the binary logit; mlogit's rank-ordered logit on its Game data", {
  set.seed(13)
  d2 <- rrologit(400, beta = cbind(b = c(0.3, -0.8)))
  m2 <- rologit(rank ~ 1 | x1, data = d2, id = "id", alt = "alt")
  g2 <- glm(I(rank == 1) ~ x1, family = binomial, data = d2[d2$alt == "b", ])
  expect_equal(unname(coef(m2)), unname(coef(g2)), tolerance = 1e-6)
  expect_equal(m2$loglik, as.numeric(logLik(g2)), tolerance = 1e-8)
  skip_if_not_installed("mlogit"); skip_if_not_installed("dfidx")
  data("Game", package = "mlogit", envir = environment())
  alts <- c("GameBoy", "GameCube", "PC", "PlayStation", "PSPortable", "Xbox")
  long <- do.call(rbind, lapply(seq_len(nrow(Game)), function(i)
    data.frame(id = i, alt = alts, rank = unlist(Game[i, paste0("ch.", alts)]), own = unlist(Game[i, paste0("own.", alts)]),
               hours = Game$hours[i], age = Game$age[i])))
  long$alt <- factor(long$alt, levels = alts)
  mr <- rologit(rank ~ own | hours + age, data = long, id = "id", alt = "alt", ref = "PC")
  G <- dfidx::dfidx(Game, varying = 1:12, choice = "ch", ranked = TRUE, idnames = c("chid", "alt"))
  mg <- mlogit::mlogit(ch ~ own | hours + age, G, reflevel = "PC")
  conv <- vapply(names(coef(mr)), function(s) if (s == "own") "own" else { p <- strsplit(s, ":", fixed = TRUE)[[1L]]; paste0(p[2L], ":", p[1L]) }, character(1))
  expect_equal(unname(coef(mr)), unname(coef(mg)[conv]), tolerance = 1e-5)
  expect_equal(mr$loglik, as.numeric(logLik(mg)), tolerance = 1e-8)
})

test_that("rologit(): weights are frequency weights; the bootstrap, term_test(), and the methods", {
  set.seed(14)
  d <- rrologit(200, beta = .Bro, gamma = 0.7, depth = 2)
  wt <- rep(1:2, length.out = 200)
  dw <- d; dw$wt <- wt[dw$id]
  mw <- rologit(rank ~ w1 | x1 + x2, data = dw, id = "id", alt = "alt", weights = "wt")
  dd <- d[rep(seq_len(nrow(d)), wt[d$id]), ]
  dd$id <- paste(dd$id, ave(seq_len(nrow(dd)), dd$id, dd$alt, FUN = seq_along))
  md <- rologit(rank ~ w1 | x1 + x2, data = dd, id = "id", alt = "alt")
  expect_equal(coef(mw), coef(md), tolerance = 1e-8)
  expect_equal(mw$loglik, md$loglik, tolerance = 1e-8)
  expect_equal(vcov(mw), vcov(md), tolerance = 1e-6)
  m <- rologit(rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
  set.seed(2); mb <- rologit(rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt", se = "bootstrap", nboot = 20L)
  expect_equal(dim(mb$boot$replicates), c(20L, 10L))
  tt <- term_test(m, type = "lr")
  expect_equal(nrow(tt), 3L)
  m_nox2 <- rologit(rank ~ w1 | x1, data = d, id = "id", alt = "alt")
  expect_equal(tt$statistic[tt$term == "x2"], 2 * (m$loglik - m_nox2$loglik), tolerance = 1e-6)
  expect_equal(lr_test(m_nox2, m)$statistic, 2 * (m$loglik - m_nox2$loglik))
  expect_equal(sum(log(predict(m, type = "ranking"))), m$loglik, tolerance = 1e-10)
  expect_equal(unname(rowSums(fitted(m))), rep(1, m$n), tolerance = 1e-12)
  fs <- fit_stats(m)
  m0 <- rologit(rank ~ 1, data = d, id = "id", alt = "alt")
  expect_equal(fs$loglik_null, m0$loglik, tolerance = 1e-8)
  expect_equal(classification(m, newdata = d)$correct, classification(m)$correct)
  s <- simulate(m, nsim = 2, seed = 3)
  expect_length(s, 2L); expect_equal(dim(s[[1L]]), c(m$n, 4L))
  expect_equal(unname(rowSums(s[[1L]] > 0)), rowSums(m$design$str$rank > 0))
  fd <- first_difference(m, "x1", from = 0, to = 1)
  expect_s3_class(fd, "pc_fd")
})

test_that("rank_test(): the split fit nests the rologit fit; refusals of the single-choice tools", {
  set.seed(15)
  d <- rrologit(300, beta = .Bro, gamma = 0.7)
  m <- rologit(rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
  rt <- rank_test(m, after = 1)
  expect_equal(rt$df, c(10L, 10L))
  cf <- attr(rt, "coefficients")
  expect_equal(nrow(cf), 10L)
  expect_true(all(rt$statistic >= 0))
  ## the late block alone, refit on the stages after position 1 (pseudo-cases), equals the split fit's late coefficients
  ex <- .explode(d); ex$alt <- factor(ex$alt, levels = levels(d$alt))
  late <- mnl(chosen ~ w1 | x1 + x2, data = ex[!grepl(" 1$", ex$pid), ], id = "pid", alt = "alt")
  expect_equal(unname(cf$late), unname(coef(late)), tolerance = 1e-6)
  expect_error(rank_test(m, after = 3), "No ranking has a stage after position 3")
  expect_error(rank_test(mnl(y ~ x1, data = rmnl(100, beta = cbind(b = c(0, 1))))), "takes a rologit")
  m1 <- mnl(I(rank == 1) ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
  expect_error(vuong(m, m1), "not comparable")
  expect_error(lr_test(m1, m), "not nested")
  expect_error(hl_test(m), "single-choice")
  expect_error(hausman_mcfadden(m, drop = "b"), "multinomial logit")
  ## boot_lr_test() takes two nested rologit() fits since DESIGN.md Phase 4c (test-suprema.R), but not a ranking against a single choice
  expect_error(boot_lr_test(m1, m, B = 10), "two rologit\\(\\) fits")
  expect_error(rologit(rank ~ x1, data = d), "long layout")
  bad <- d; bad$rank[1] <- -1
  expect_error(rologit(rank ~ w1 | x1, data = bad, id = "id", alt = "alt"), "rank of each alternative")
  big <- rrologit(30, beta = cbind(b = c(0, 1), c = c(0, 1), d = c(0, 1), e = c(0, 1), f = c(0, 1), g = c(0, 1), h = c(0, 1), i = c(0, 1)))
  big$rank[big$id == 1] <- c(rep(1, 8), 2)
  expect_error(rologit(rank ~ 1 | x1, data = big, id = "id", alt = "alt"), "tie more than 7")
  md3 <- rologit(rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt", depth = 2)
  expect_error(lr_test(md3, m), "different rankings")
})
