## Firth's penalty (firth = TRUE in mnl() and rologit(); DESIGN.md, section 1,
## Phase 4): the oracles brglm2::brmultinom, logistf::logistf and
## coxphf::coxphf (a stratified Cox model through disjoint time windows), the
## penalized gradient against numDeriv, the tie information against a
## brute-force enumeration, invariances, the penalized likelihood-ratio and
## profile machinery, and the refusals. The reference packages are called with
## `::` and never attached: attaching logistf breaks coxphftest()'s formula
## handling.

.fsim <- function(n, seed) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rbinom(n, 1, 0.3)
  X <- cbind(1, x1, x2); B <- cbind(0, c(0.3, 0.8, -1), c(-0.5, -0.4, 1.5))
  pr <- exp(X %*% B); pr <- pr / rowSums(pr)
  y <- factor(apply(pr, 1, function(p) sample(3, 1, prob = p)), labels = c("a", "b", "c"))
  data.frame(y, x1, x2)
}
## rankings of J alternatives by n cases, alternative attributes x (normal) and z (binary)
.frank <- function(n, J, seed, depth = J - 1L) {
  set.seed(seed)
  d <- do.call(rbind, lapply(seq_len(n), function(i) {
    x <- rnorm(J); z <- rbinom(J, 1, 0.5); u <- 0.8 * x + 0.5 * z - log(-log(runif(J)))
    data.frame(case = i, alt = factor(seq_len(J)), x = x, z = z, rank = rank(-u))
  }))
  d$rank[d$rank > depth] <- 0L
  d
}
## the same rankings as a Cox model with each case on its own time window
.fwindows <- function(d) {
  depth <- max(d$rank)
  d$event <- as.integer(d$rank > 0); d$time <- ifelse(d$event == 1L, d$rank, depth + 1L)
  d$start <- 10 * d$case; d$stop <- d$start + d$time
  d
}

test_that("firth = TRUE: finite under separation, the ML fit's divergence named, invariances, and firth = FALSE unchanged", {
  d <- .fsim(60, 3)
  d$y[d$x2 == 1 & d$y == "c"] <- "a"                      # c never chosen when x2 = 1: c:x2 separated
  expect_message(m <- mnl(y ~ x1 + x2, data = d, firth = TRUE), "c:x2 diverge")
  expect_true(m$converged); expect_true(isTRUE(m$firth)); expect_false(m$boundary)
  expect_identical(m$separation, "c:x2")
  expect_true(all(is.finite(coef(m))) && abs(coef(m)["c:x2"]) < 10)
  ml <- suppressWarnings(mnl(y ~ x1 + x2, data = d))
  expect_lt(m$loglik, ml$loglik)                          # the ordinary log-likelihood at the penalized estimates
  expect_equal(m$loglik_penalized, m$loglik + polychoice:::.pc_fns(m$design)$penalty(m$theta), tolerance = 1e-10)
  ## the penalty is invariant to the base alternative
  mr <- suppressMessages(mnl(y ~ x1 + x2, data = d, ref = "c", firth = TRUE))
  expect_equal(fitted(mr), fitted(m), tolerance = 1e-8)
  expect_equal(mr$loglik_penalized, m$loglik_penalized, tolerance = 1e-10)
  ## frequency weights are replicated cases
  w <- rep(1:2, length.out = nrow(d))
  mw <- suppressMessages(mnl(y ~ x1 + x2, data = cbind(d, w = w), weights = "w", firth = TRUE))
  md <- suppressMessages(mnl(y ~ x1 + x2, data = d[rep(seq_len(nrow(d)), w), ], firth = TRUE))
  expect_equal(coef(mw), coef(md), tolerance = 1e-8)
  expect_equal(vcov(mw), vcov(md), tolerance = 1e-8)
  ## the default is maximum likelihood, bit for bit
  d0 <- .fsim(60, 3)
  expect_identical(coef(mnl(y ~ x1 + x2, data = d0)), coef(mnl(y ~ x1 + x2, data = d0, firth = FALSE)))
  expect_null(mnl(y ~ x1 + x2, data = d0)$firth)
  ## the adjusted scores sum to zero at the estimate
  S <- polychoice:::.pc_fns(m$design)$scores(m$theta)
  expect_lt(max(abs(colSums(m$weights * S))), 1e-4)
})

test_that("firth = TRUE: refusals", {
  d <- .fsim(80, 5)
  m <- mnl(y ~ x1 + x2, data = d, firth = TRUE); ml <- mnl(y ~ x1 + x2, data = d)
  expect_error(hausman_mcfadden(m), "does not take a penalized fit")
  expect_error(small_hsiao(m, seed = 1), "does not take a penalized fit")
  expect_error(fp_tests(m), "does not take a penalized fit")
  expect_error(iia_tests(m), "does not take a penalized fit")
  expect_error(captivity_test(m), "does not take a penalized fit")
  expect_error(nested_lm_test(m, nests = list(n1 = c("b", "c"))), "does not take a penalized fit")
  expect_error(vuong(m, ml), "does not take a penalized fit")
  expect_error(lr_test(mnl(y ~ x1, data = d), m), "One fit is penalized")
  expect_error(mnl(y ~ x1, data = cbind(d, u = rep(1:8, 10)), fe = "u", firth = TRUE), "not available with unit fixed effects")
  expect_error(mnl(y ~ x1, data = d, firth = NA), "TRUE or FALSE")
  expect_error(confint(ml, type = "profile"), "penalized fits")
})

test_that("firth = TRUE: penalized likelihood-ratio and profile machinery", {
  d <- .fsim(80, 5)
  m <- mnl(y ~ x1 + x2, data = d, firth = TRUE)
  m1 <- mnl(y ~ x1, data = d, firth = TRUE)
  lr <- lr_test(m1, m)
  tt <- term_test(m, type = "lr")
  ## lr_test() of two penalized fits restricts the full model's penalized likelihood
  ## (not the separately penalized restricted fit), as term_test() does
  expect_equal(lr$statistic, tt$statistic[tt$term == "x2"], tolerance = 1e-6)
  expect_true(isTRUE(attr(lr, "penalized")))
  expect_false(isTRUE(all.equal(lr$statistic, 2 * (m$loglik_penalized - m1$loglik_penalized))))
  cb <- combine_test(m, c("b", "c"), type = "lr")
  expect_true(is.finite(cb$statistic) && cb$statistic > 0)
  expect_match(attr(cb, "method"), "penalized")
  ## profile interval limits: twice the drop of the profile penalized likelihood is the chi-square quantile
  ci <- confint(m, "b:x1")
  fn <- polychoice:::.pc_fns(m$design); r <- match("b:x1", names(coef(m)))
  for (lim in ci) {
    ## the maximum over the other parameters with theta_r held at the limit
    free <- setdiff(seq_along(m$theta), r)
    f <- function(b) { t <- m$theta; t[free] <- b; t[r] <- lim; fn$nll(t) }
    op <- optim(m$theta[free], f, method = "BFGS", control = list(reltol = 1e-14, maxit = 1000))
    expect_equal(2 * (m$loglik_penalized + op$value), stats::qchisq(0.95, 1), tolerance = 1e-4)
  }
  expect_equal(confint(m, type = "wald"), { s <- sqrt(diag(vcov(m))); cbind(coef(m) - 1.959964 * s, coef(m) + 1.959964 * s) }, tolerance = 1e-6, ignore_attr = TRUE)
  ## reporting
  fs <- fit_stats(m)
  expect_equal(fs$loglik_penalized, m$loglik_penalized)
  expect_true(isTRUE(attr(logLik(m), "penalized")))
  expect_output(print(summary(m)), "Penalized \\(Firth\\) estimates")
  cm <- compare_models(ml = mnl(y ~ x1 + x2, data = d), firth = m)
  expect_identical(cm$penalized, c(FALSE, TRUE))
  skip_if_not_installed("broom")
  expect_true(broom::glance(m)$penalized)
})

test_that("firth = TRUE: mnl() against brglm2::brmultinom (estimates and standard errors)", {
  skip_on_cran()
  skip_if_not_installed("brglm2")
  d <- .fsim(60, 3)
  m <- mnl(y ~ x1 + x2, data = d, firth = TRUE)
  bm <- brglm2::brmultinom(y ~ x1 + x2, data = d, ref = 1, epsilon = 1e-12, maxit = 500)
  expect_equal(unname(coef(m)), as.vector(t(coef(bm))), tolerance = 1e-8)
  expect_equal(unname(sqrt(diag(vcov(m)))), unname(sqrt(diag(vcov(bm)))), tolerance = 1e-8)
  ## weighted
  w <- rep(c(1, 3), length.out = nrow(d))
  mw <- mnl(y ~ x1 + x2, data = cbind(d, w = w), weights = "w", firth = TRUE)
  bw <- brglm2::brmultinom(y ~ x1 + x2, data = cbind(d, w = w), weights = w, ref = 1, epsilon = 1e-12, maxit = 500)
  expect_equal(unname(coef(mw)), as.vector(t(coef(bw))), tolerance = 1e-8)
})

test_that("firth = TRUE: rologit() and the long-layout mnl() against coxphf through disjoint time windows", {
  skip_on_cran()
  skip_if_not_installed("coxphf")
  skip_if_not_installed("survival")
  Surv <- survival::Surv
  for (depth in c(3L, 2L)) {
    d <- .frank(50, 4, 7, depth)
    f <- rologit(rank ~ 0 + x + z, data = d, id = "case", alt = "alt", firth = TRUE)
    dw <- .fwindows(d)
    cf <- coxphf::coxphf(Surv(start, stop, event) ~ x + z, data = dw, firth = TRUE, pl = TRUE, maxit = 200)
    expect_equal(unname(coef(f)), unname(coef(cf)), tolerance = 1e-7)
    expect_equal(unname(vcov(f)), unname(cf$var), tolerance = 1e-7)
    expect_equal(unname(confint(f)), unname(log(cbind(cf$ci.lower, cf$ci.upper))), tolerance = 1e-6)
    ct <- coxphf::coxphftest(Surv(start, stop, event) ~ x + z, data = dw, test = ~ z)
    tt <- term_test(f, type = "lr")
    expect_equal(tt$statistic[tt$term == "z"], 2 * diff(ct$loglik), tolerance = 1e-6)
  }
  ## varying choice sets: some alternatives unavailable to some cases
  d <- .frank(60, 4, 8, 2L)
  d <- d[!(d$case %% 3 == 0 & d$alt == "4" & d$rank == 0), ]
  f <- rologit(rank ~ 0 + x + z, data = d, id = "case", alt = "alt", firth = TRUE)
  cf <- coxphf::coxphf(Surv(start, stop, event) ~ x + z, data = .fwindows(d), firth = TRUE, pl = FALSE, maxit = 200)
  expect_equal(unname(coef(f)), unname(coef(cf)), tolerance = 1e-7)
  ## the conditional logit (one choice per case) is a Cox model with one event per stratum
  d1 <- .frank(60, 4, 9, 1L); d1$chosen <- as.integer(d1$rank == 1L)
  mc <- mnl(chosen ~ 0 + x + z, data = d1, id = "case", alt = "alt", firth = TRUE)
  cf1 <- coxphf::coxphf(Surv(start, stop, event) ~ x + z, data = .fwindows(d1), firth = TRUE, pl = FALSE, maxit = 200)
  expect_equal(unname(coef(mc)), unname(coef(cf1)), tolerance = 1e-7)
  ## a separated ranking design: coxphf needs more iterations and a smaller step
  ds <- .frank(50, 4, 7); ds$x <- ifelse(ds$rank > 0, 5 - ds$rank, 0) + 0.01 * rnorm(nrow(ds))   # x orders every ranking
  fs <- suppressMessages(rologit(rank ~ 0 + x + z, data = ds, id = "case", alt = "alt", firth = TRUE))
  expect_true("x" %in% fs$separation)
  cs <- coxphf::coxphf(Surv(start, stop, event) ~ x + z, data = .fwindows(ds), firth = TRUE, pl = FALSE, maxit = 2000, maxstep = 0.1)
  expect_equal(unname(coef(fs)), unname(coef(cs)), tolerance = 1e-4)
})

test_that("firth = TRUE: two alternatives against logistf::logistf (case- and alternative-specific covariates)", {
  skip_on_cran()
  skip_if_not_installed("logistf")
  d <- .fsim(60, 3)
  d$y2 <- factor(ifelse(d$y == "a", "a", "b"))
  m <- mnl(y2 ~ x1 + x2, data = d, firth = TRUE)
  lf <- logistf::logistf(I(y2 == "b") ~ x1 + x2, data = d)
  expect_equal(unname(coef(m)), unname(coef(lf)), tolerance = 1e-6)
  expect_equal(unname(confint(m)), unname(cbind(lf$ci.lower, lf$ci.upper)), tolerance = 1e-5)
  tt <- term_test(m, type = "lr")                         # one coefficient per term at J = 2
  expect_equal(tt$p.value, unname(lf$prob[-1L]), tolerance = 1e-5)
  ## the conditional logit with J = 2 is the logit on the differenced attributes
  set.seed(4); n <- 70
  w0 <- rnorm(n); w1 <- rnorm(n); ch <- rbinom(n, 1, plogis(0.3 + 0.9 * (w1 - w0)))
  long <- data.frame(id = rep(seq_len(n), each = 2), alt = factor(rep(c("A", "B"), n)), w = as.vector(rbind(w0, w1)), chosen = as.vector(rbind(1 - ch, ch)))
  mc <- mnl(chosen ~ w, data = long, id = "id", alt = "alt", firth = TRUE)
  lc <- logistf::logistf(ch ~ dw, data = data.frame(ch = ch, dw = w1 - w0))
  expect_equal(unname(coef(mc)), unname(coef(lc)), tolerance = 1e-6)   # the constant of B, then w
})

test_that("firth = TRUE: the penalized gradient against numDeriv (weights, offsets, choice sets, ties)", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  d <- .fsim(50, 11); d$w <- rep(c(1, 2, 0.5), length.out = 50)
  off <- matrix(rnorm(150, sd = 0.3), 50, 3)
  m <- mnl(y ~ x1 + x2, data = d, weights = "w", offset = off, firth = TRUE)
  fn <- polychoice:::.pc_fns(m$design); th <- m$theta + 0.1
  expect_equal(fn$grad(th), numDeriv::grad(fn$nll, th), tolerance = 1e-6)
  dl <- .frank(40, 4, 12, 1L); dl$chosen <- as.integer(dl$rank == 1L)
  dl <- dl[!(dl$case %% 4 == 0 & dl$alt == "2" & dl$chosen == 0L), ]
  ml <- mnl(chosen ~ x + z, data = dl, id = "case", alt = "alt", firth = TRUE)
  fl <- polychoice:::.pc_fns(ml$design); tl <- ml$theta - 0.1
  expect_equal(fl$grad(tl), numDeriv::grad(fl$nll, tl), tolerance = 1e-6)
  ## rankings with ties of two and three and partial rankings
  dt <- .frank(40, 5, 13, 4L)
  dt$rank[dt$case %% 2 == 0 & dt$rank == 3] <- 2L
  dt$rank[dt$case %% 3 == 0 & dt$rank %in% 2:4] <- 2L
  ft <- rologit(rank ~ x + z, data = dt, id = "case", alt = "alt", firth = TRUE)
  ff <- polychoice:::.pc_fns(ft$design); tt <- ft$theta + 0.05
  expect_equal(ff$grad(tt), numDeriv::grad(ff$nll, tt), tolerance = 1e-6)
  expect_true(ft$converged)
  ## rank_test() of a penalized fit: both maxima under the split model's penalty
  rt <- rank_test(ft)
  expect_true(all(is.finite(rt$statistic)))
})

test_that("firth = TRUE: covariates on very different scales (the information equilibrated before it is inverted)", {
  ## in theta units the information inherits the covariates' scales (Poast's
  ## 2007 IOC bids: 0.3 to 1.7e8); its reciprocal condition number falls below
  ## machine precision, solve() refuses it, and a generalized inverse in its
  ## place gave a wrong penalty gradient and an unconverged fit
  dr <- .frank(40, 4, 31); db <- dr; db$x <- db$x * 1e9
  fr <- rologit(rank ~ 0 + x + z, data = dr, id = "case", alt = "alt", firth = TRUE)
  fb <- rologit(rank ~ 0 + x + z, data = db, id = "case", alt = "alt", firth = TRUE)
  expect_lt(rcond(polychoice:::.pc_fns(fb$design)$info(fb$theta)), .Machine$double.eps)
  expect_true(fb$converged)
  expect_equal(unname(coef(fb)["x"]) * 1e9, unname(coef(fr)["x"]), tolerance = 1e-6)
  expect_equal(unname(coef(fb)["z"]), unname(coef(fr)["z"]), tolerance = 1e-6)
  expect_equal(fb$loglik, fr$loglik, tolerance = 1e-8)
  expect_equal(fb$loglik_penalized, fr$loglik_penalized + log(1e9), tolerance = 1e-8)   # 1/2 log|I| shifts by log c
  d <- .fsim(60, 5); d2 <- d; d2$x1 <- d2$x1 * 1e9
  m <- mnl(y ~ x1 + x2, data = d, firth = TRUE); m2 <- mnl(y ~ x1 + x2, data = d2, firth = TRUE)
  expect_true(m2$converged)
  expect_equal(fitted(m2), fitted(m), tolerance = 1e-6)
})

test_that("firth = TRUE: the tie information against a brute-force sum over subsets and orders", {
  skip_on_cran()
  skip_if_not_installed("numDeriv")
  dt <- .frank(12, 5, 21, 4L)
  dt$rank[dt$case %% 2 == 0 & dt$rank == 3] <- 2L
  dt$rank[dt$case %% 3 == 0 & dt$rank %in% 2:4] <- 2L
  ft <- rologit(rank ~ x + z, data = dt, id = "case", alt = "alt", firth = TRUE)
  D <- ft$design; th <- ft$theta
  perms <- function(v) if (length(v) <= 1L) list(v) else do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
  logpi <- function(V, S) { B <- setdiff(seq_along(V), S)
    log(sum(vapply(perms(S), function(p) prod(vapply(seq_along(p), function(s) exp(V[p[s]]) / sum(exp(V[c(p[s:length(p)], B)])), 0)), 0))) }
  row <- function(i, j) { v <- numeric(D$K); nb <- setdiff(seq_len(D$J), D$base)
    k <- match(j, nb); if (!is.na(k)) v[D$idx$beta[, k]] <- D$X[i, ]
    v[D$idx$gamma] <- D$W[(i - 1L) * D$J + j, ]; v }
  I <- matrix(0, D$K, D$K)
  for (i in seq_len(D$n)) {
    r <- D$str$rank[i, ]; av <- which(D$avail[i, ] == 1L); ranked <- av[r[av] > 0]
    for (g in sort(unique(r[ranked]))) {
      Tg <- ranked[r[ranked] == g]; below <- c(ranked[r[ranked] > g], av[r[av] == 0])
      if (!length(below)) next
      R <- c(Tg, below); DR <- t(vapply(R, function(j) row(i, j), numeric(D$K))); V <- as.vector(DR %*% th)
      for (S in utils::combn(length(R), length(Tg), simplify = FALSE)) {
        gS <- numDeriv::grad(function(v) logpi(v, S), V)
        I <- I + D$w[i] * exp(logpi(V, S)) * crossprod(DR, tcrossprod(gS) %*% DR)
      }
    }
  }
  expect_equal(polychoice:::.pc_fns(D)$info(th), I, tolerance = 1e-7)
})

test_that("firth = TRUE: robust, cluster and bootstrap standard errors", {
  skip_on_cran()
  d <- .fsim(80, 5); d$cl <- rep(1:20, 4)
  mr <- mnl(y ~ x1 + x2, data = d, se = "robust", firth = TRUE)
  ma <- mnl(y ~ x1 + x2, data = d, firth = TRUE)
  expect_equal(coef(mr), coef(ma))
  expect_true(all(is.finite(sqrt(diag(vcov(mr))))))
  mc <- mnl(y ~ x1 + x2, data = d, cluster = "cl", firth = TRUE)
  expect_true(all(is.finite(sqrt(diag(vcov(mc))))))
  set.seed(2)
  mb <- mnl(y ~ x1 + x2, data = d, se = "bootstrap", nboot = 20, firth = TRUE)
  expect_true(all(is.finite(sqrt(diag(vcov(mb))))))
  expect_true(isTRUE(mb$design$firth))
})

test_that("defection: Cook, Niehaus and Zuhlke's (2018) penalized Models 3-5 reproduce (their replication log, brglm2)", {
  skip_on_cran()                                    # slow; runs in the full suite (check time, 2026-11-06)
  ## demobilized estimate, se; switched estimate, se (Warning_RAP_Comment_RLog.txt in their Dataverse deposit)
  pub <- list(
    m3 = rbind(economic = c(0.9031805751, 1.59797978, 0.69757593, 1.75029952), `economic:indoctrination` = c(-0.6321131557, 1.66418546, 0.77267120, 1.87034678),
               indoctrination = c(0.1568017149, 0.54424642, -1.96551679, 0.71733740), age_joined = c(0.0816491447, 0.04198641, 0.02530272, 0.05299049),
               birth_year = c(0.0001572316, 0.03020921, -0.08442600, 0.03866640), school_age = c(0.0245105956, 0.03386011, 0.08269713, 0.05106702),
               male = c(0.7597087600, 0.31491293, 1.12274475, 0.61665525), `(Intercept)` = c(-0.3780652830, 60.31547373, 165.51298419, 77.21848203)),
    m4 = rbind(ideological = c(-0.236534862, 0.45926956, -1.90341969, 0.98878284), `ideological:abuse` = c(-1.094041361, 1.89519333, 0.01511906, 2.35259547),
               abuse = c(0.318759244, 0.69609056, 0.93396477, 0.96239609), age_joined = c(0.088786614, 0.04217539, 0.04539671, 0.05339326),
               birth_year = c(-0.002376447, 0.03008903, -0.09132592, 0.03977966), school_age = c(0.021315697, 0.03422214, 0.09443064, 0.05216715),
               male = c(0.739919943, 0.31430080, 1.01765641, 0.60457714), `(Intercept)` = c(4.787093556, 60.03859208, 177.57786800, 79.33317491)),
    m5 = rbind(economic = c(0.590380490, 0.51543869, 1.59667822, 0.65389416), `economic:abuse` = c(-2.309670766, 1.48970024, -3.24454028, 2.25454527),
               abuse = c(0.823347256, 0.86893104, 1.80528722, 1.05983176), age_joined = c(0.079200085, 0.04272320, 0.04103121, 0.05367542),
               birth_year = c(-0.003589635, 0.03080472, -0.07210012, 0.03873285), school_age = c(0.021830564, 0.03427289, 0.08227029, 0.05159915),
               male = c(0.758994777, 0.31602762, 1.07425404, 0.60841462), `(Intercept)` = c(7.217791702, 61.47177644, 139.21534137, 77.30933517)))
  rhs <- list(m3 = "economic * indoctrination", m4 = "ideological * abuse", m5 = "economic * abuse")
  for (m in names(pub)) {
    f <- suppressMessages(mnl(stats::as.formula(paste("outcome ~", rhs[[m]], "+ age_joined + birth_year + school_age + male")),
                              data = defection, ref = "captured", firth = TRUE))
    expect_true(f$converged)
    P <- pub[[m]]; tm <- rownames(P); se <- sqrt(diag(vcov(f)))
    ## element by element (testthat's tolerance is a mean over the vector, which
    ## the constants of 60-177 would dominate): estimates within 1e-4 of their
    ## published standard errors, standard errors within 1e-4 relative
    est <- unname(f$coefficients[c(paste0("demobilized:", tm), paste0("switched:", tm))]); pse <- unname(c(P[, 2], P[, 4]))
    expect_lt(max(abs(est - unname(c(P[, 1], P[, 3]))) / pse), 1e-4)
    expect_lt(max(abs(unname(se[c(paste0("demobilized:", tm), paste0("switched:", tm))]) - pse) / pse), 1e-4)
  }
  ## the separated coefficients are those of the empty cells (Cook et al.'s Table 2), not the imprecise constants
  f3 <- suppressMessages(mnl(outcome ~ economic * indoctrination + age_joined + birth_year + school_age + male, data = defection, ref = "captured", firth = TRUE))
  expect_false(any(grepl("Intercept", f3$separation)))
  expect_true(all(c("demobilized:economic", "demobilized:economic:indoctrination") %in% f3$separation))
  ## the restricted maxima of the penalized LR tests start from zero too: zeroing
  ## birth_year from the full estimates leaves the constants near 165
  expect_no_warning(tt <- term_test(f3, type = "lr"))
  expect_lt(tt$statistic[tt$term == "birth_year"], 20)
  ## the profile intervals of the separated coefficients: a warm-started
  ## search point near a limit that stops short is not a failure when the
  ## limit, solved afresh, reproduces the quantile (no warning; the limits of
  ## 2026-11-06, which plain optim() confirms to six decimals of the deviance)
  sel <- c("demobilized:economic", "demobilized:economic:indoctrination")
  expect_no_warning(ci <- confint(f3, sel))
  expect_equal(unname(ci), rbind(c(-1.48854, 5.82850), c(-5.60934, 1.96956)), tolerance = 1e-5)
})
