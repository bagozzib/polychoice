## Firth's penalty (firth = TRUE): bias, interval coverage, and the size of the
## penalized likelihood-ratio test (DESIGN.md, section 1, Phase 4; designs and
## criteria fixed 2026-11-05, before any run).
##   F1  multinomial logit, J = 3, wide, n = 50 and 200 (after Bull, Mak and
##       Greenwood 2002, Table 6, models 9 and 11): x1 = 1{z1 > 0}, x2 = z2,
##       (z1, z2) bivariate normal with correlation 0.752 (cor(x1, x2) = 0.6),
##       and a null covariate x3 ~ N(0, 1).
##   F2  rank-ordered logit, J = 4, n = 100 and 300 cases, complete rankings
##       and the same data cut at depth 1 and 2: constants b 0.5, c 0.25,
##       d -0.25 (a the base), a normal attribute (0.8), a binary attribute with
##       prevalence 0.1 (1.5) and a null normal attribute (0).
##   F3  F2 at n = 100 with ties: in each ranking, independently of the
##       utilities, positions 1-3 are tied with probability 0.1, else positions
##       2-3 with probability 0.3.
## Reported per cell: the share of samples whose maximum-likelihood fit is
## separated (the penalized fit's diagnostic); for the penalized fit
## convergence, mean bias, RMSE and median absolute deviation over all
## samples, and maximum likelihood's over the samples where it converged
## without a boundary; coverage of the 95 percent profile and Wald intervals
## for every coefficient; the size at 0.05 of the penalized LR test of the
## null term; the time per replication.
##   F4  (added 2026-11-06, DESIGN.md Phase 4d, before any run) a design
##       following that in Rainey (2024, section 3.1): n = 1000; a binary s equal to 1 for the
##       five highest of a latent normal correlated 0.06 with three controls
##       z1-z3 ~ N(0, 0.5^2) (coefficients 1); constant -4.1. J = 2 (mnl() on
##       two alternatives, the binary logit) and J = 3 (s in the b equation
##       only; the c equation has constant -1, the controls with coefficient
##       1, and no s), with beta_s in {0, 2, -2}. Reported: the rejection rate
##       at 0.05 of the test of beta_s by Firth's Wald test, the penalized
##       likelihood ratio (lr_test() between the Firth fits with and without
##       s) and the maximum-likelihood ratio from the suprema (lr_test()
##       between the ML fits), and the share of samples separated; also the
##       three rates among the samples where beta_s itself is separated,
##       where it is not, and where no s = 1 case chooses b (added in round
##       6, before any run).
## Criteria: (a) every penalized fit converges; (b) profile coverage within
## [0.93, 0.97] for every coefficient at the larger n of F1 and F2 and in F3;
## (c) the penalized LR test's size at 0.05 within [0.035, 0.065] in every
## cell; (d) at the larger n, every penalized estimate's mean bias within 3
## Monte Carlo standard errors or 5 percent of its true value; (e) in F4 at
## beta_s = 0 the penalized LR and the ML LR from the suprema reject in at
## most 0.08 of samples; Firth's Wald rate is reported, the evidence for
## summary()'s [sep] marks. The summary() rule of DESIGN.md is read from (b)'s
## Wald column (Wald coverage outside [0.92, 0.98] where profile coverage is
## inside).
## Round 6 (2026-11-06, before any run; criteria unchanged): the Firth Wald
## rate is reported overall and among the samples where beta_s is separated.
## Round 7 (2026-11-06, before any run; criteria unchanged): with five s = 1
## cases and a base rate near 0.024, F4 carries almost no information about
## beta_s. No event among the five (probability about 0.88 at beta_s = 0)
## gives an ML ratio of about 0.25, one event about 2.6, and only two or more
## (probability about 0.006) pass 3.84. So criterion (e) can hardly fail, and
## the power at beta_s = -2 is near 0 by construction. F5 adds the same
## design with 40 s = 1 cases (J = 2 and 3, beta_s in {0, 2, -2}), where the
## tests can differ; it is reported, with no criterion.
##   F6  (added 2026-11-06 at the author's decision, before any run) the
##       covariates of design 911 of Rainey's (2024) replication archive
##       (Harvard Dataverse doi:10.7910/DVN/6EYRJG, CC0; verification/data/
##       rainey2024_dgp911.csv: n = 1000, five s = 1, three controls), held
##       fixed across samples; constant -4.1, the controls'
##       coefficients 1; J = 2, beta_s in {0, 2, -2}. Reported as F4, with the
##       share of samples with an empty cell of the s-by-y table.
## Criterion (f): in each F6 cell, the ML LR test's rejection rate at 0.05 is
## within max(0.01, 3 Monte Carlo standard errors of the difference) of the
## rate of the archive's 2,500 samples per cell (output/summarized-
## simulations.rds: 0.0624 / 0.4460 / 0.4336 at beta_s = 0 / 2 / -2; ours is
## the ML LR from the suprema). A failure is audited before release. As first
## written (2026-11-06), (f) also compared Firth's Wald rate; it was narrowed
## before the run, at the author's decision.
## Run from the package root:  Rscript verification/firth_mc.R
## Env: FIRTH_REPS (default 1000), FIRTH_WORKERS (default 6; PSOCK workers).
suppressPackageStartupMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."
log <- file.path(out_dir, "firth_mc.log"); rds <- file.path(out_dir, "firth_mc.rds")
cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
R <- as.integer(Sys.getenv("FIRTH_REPS", "1000")); W <- as.integer(Sys.getenv("FIRTH_WORKERS", "6"))

B9 <- cbind(b = c(-1.4, 1, 1, 0), c = c(-1.4, 1, -1, 0))        # rows: (Intercept), x1, x2, x3
B11 <- cbind(b = c(-1.4, 2, 1, 0), c = c(0, 2, 1, 0))
cells <- list(
  list(id = "F1 model 9, n = 50", fam = "mnl", B = B9, n = 50L, large = FALSE),
  list(id = "F1 model 9, n = 200", fam = "mnl", B = B9, n = 200L, large = TRUE),
  list(id = "F1 model 11, n = 50", fam = "mnl", B = B11, n = 50L, large = FALSE),
  list(id = "F1 model 11, n = 200", fam = "mnl", B = B11, n = 200L, large = TRUE),
  list(id = "F2 complete, n = 100", fam = "ro", n = 100L, depth = 3L, ties = FALSE, large = FALSE),
  list(id = "F2 depth 2, n = 100", fam = "ro", n = 100L, depth = 2L, ties = FALSE, large = FALSE),
  list(id = "F2 depth 1, n = 100", fam = "ro", n = 100L, depth = 1L, ties = FALSE, large = FALSE),
  list(id = "F2 complete, n = 300", fam = "ro", n = 300L, depth = 3L, ties = FALSE, large = TRUE),
  list(id = "F2 depth 2, n = 300", fam = "ro", n = 300L, depth = 2L, ties = FALSE, large = TRUE),
  list(id = "F2 depth 1, n = 300", fam = "ro", n = 300L, depth = 1L, ties = FALSE, large = TRUE),
  list(id = "F3 ties, n = 100", fam = "ro", n = 100L, depth = 3L, ties = TRUE, large = TRUE))
ro_truth <- c(`b:(Intercept)` = 0.5, `c:(Intercept)` = 0.25, `d:(Intercept)` = -0.25, w1 = 0.8, w2 = 1.5, w3 = 0)

gen_mnl <- function(n, B) {
  z1 <- stats::rnorm(n); z2 <- 0.752 * z1 + sqrt(1 - 0.752^2) * stats::rnorm(n)
  X <- cbind(1, x1 = as.integer(z1 > 0), x2 = z2, x3 = stats::rnorm(n))
  V <- cbind(0, X %*% B); P <- exp(V - apply(V, 1L, max)); P <- P / rowSums(P)
  y <- factor(apply(P, 1L, function(p) sample.int(3L, 1L, prob = p)), levels = 1:3, labels = c("a", "b", "c"))
  data.frame(y = y, x1 = X[, 2], x2 = X[, 3], x3 = X[, 4])
}
gen_ro <- function(n, depth, ties) {
  J <- 4L; asc <- c(a = 0, b = 0.5, c = 0.25, d = -0.25)
  d <- do.call(rbind, lapply(seq_len(n), function(i) {
    w1 <- stats::rnorm(J); w2 <- stats::rbinom(J, 1L, 0.1); w3 <- stats::rnorm(J)
    u <- asc + 0.8 * w1 + 1.5 * w2 - log(-log(stats::runif(J)))
    r <- rank(-u)
    if (ties) { v <- stats::runif(1); if (v < 0.1) r[r <= 3] <- 1 else if (v < 0.4) r[r %in% 2:3] <- 2 }
    data.frame(id = i, alt = factor(names(asc), levels = names(asc)), w1 = w1, w2 = w2, w3 = w3, rank = as.integer(r))
  }))
  d$rank[d$rank > depth] <- 0L
  d
}

## one replication: the penalized fit (estimates, Wald and profile intervals,
## penalized LR test of the null term, separation) and the ML fit
one <- function(cell, seed) {
  set.seed(seed)
  out <- tryCatch({
    if (cell$fam == "mnl") {
      d <- gen_mnl(cell$n, cell$B)
      truth <- c(stats::setNames(cell$B[, "b"], paste0("b:", c("(Intercept)", "x1", "x2", "x3"))),
                 stats::setNames(cell$B[, "c"], paste0("c:", c("(Intercept)", "x1", "x2", "x3"))))
      fit <- function(firth) mnl(y ~ x1 + x2 + x3, data = d, firth = firth)
      null_term <- "x3"
    } else {
      d <- gen_ro(cell$n, cell$depth, cell$ties)
      truth <- ro_truth
      fit <- function(firth) rologit(rank ~ w1 + w2 + w3, data = d, id = "id", alt = "alt", firth = firth)
      null_term <- "w3"
    }
    t0 <- Sys.time()
    f <- suppressWarnings(suppressMessages(fit(TRUE)))
    pci <- suppressWarnings(confint(f, type = "profile"))
    wci <- confint(f, type = "wald")
    tt <- suppressWarnings(term_test(f, type = "lr"))
    secs <- as.numeric(Sys.time() - t0, units = "secs")
    m <- suppressWarnings(suppressMessages(fit(FALSE)))
    est <- coef(f)[names(truth)]
    data.frame(par = names(truth), truth = unname(truth), est = unname(est),
               cover_profile = unname(pci[names(truth), 1] <= truth & truth <= pci[names(truth), 2]),
               cover_wald = unname(wci[names(truth), 1] <= truth & truth <= wci[names(truth), 2]),
               ml = unname(coef(m)[names(truth)]), ml_ok = isTRUE(m$converged) && !isTRUE(m$boundary),
               conv = isTRUE(f$converged), separated = length(f$separation) > 0L,
               p_null = tt$p.value[tt$term == null_term], secs = secs, failed = FALSE, row.names = NULL)
  }, error = function(e) data.frame(par = NA_character_, failed = TRUE, msg = conditionMessage(e)))
  out$seed <- seed
  out
}

## F4: a design following that in Rainey (2024, section 3.1); one replication returns the three tests' p-values
gen_rainey <- function(J, beta_s, ns = 5L) {
  n <- 1000L; rho <- 0.06
  Sig <- matrix(rho, 4, 4); diag(Sig) <- 1
  L <- matrix(stats::rnorm(n * 4), n) %*% chol(Sig)
  s <- as.integer(rank(-L[, 1]) <= ns)
  z <- 0.5 * L[, 2:4]
  eta_b <- -4.1 + beta_s * s + rowSums(z)
  V <- if (J == 2L) cbind(0, eta_b) else cbind(0, eta_b, -1 + rowSums(z))
  P <- exp(V - apply(V, 1L, max)); P <- P / rowSums(P)
  y <- factor(apply(P, 1L, function(p) sample.int(J, 1L, prob = p)), levels = seq_len(J), labels = c("a", "b", "c")[seq_len(J)])
  data.frame(y = y, s = s, z1 = z[, 1], z2 = z[, 2], z3 = z[, 3])
}
## F6: the covariates of Rainey's design 911, the same in every sample
X911 <- utils::read.csv(file.path(out_dir, "data", "rainey2024_dgp911.csv"), comment.char = "#")
gen_rainey911 <- function(beta_s) {
  p <- stats::plogis(-4.1 + beta_s * X911$s + X911$z1 + X911$z2 + X911$z3)
  y <- factor(1L + stats::rbinom(nrow(X911), 1L, p), levels = 1:2, labels = c("a", "b"))
  data.frame(y = y, s = X911$s, z1 = X911$z1, z2 = X911$z2, z3 = X911$z3)
}
## s enters the b equation only: in the long layout it is the alternative-
## specific variable sb = s for alternative b (0 for the others), with a
## generic coefficient, so every test is of the one coefficient beta_s
one_rainey <- function(cell, seed) {
  set.seed(seed)
  out <- tryCatch({
    w <- if (isTRUE(cell$dgp911)) gen_rainey911(cell$beta_s) else gen_rainey(cell$J, cell$beta_s, if (is.null(cell$ns)) 5L else cell$ns)
    alts <- levels(w$y)
    d <- data.frame(id = rep(seq_len(nrow(w)), each = length(alts)), alt = factor(rep(alts, nrow(w)), levels = alts))
    d$chosen <- as.integer(w$y[d$id] == d$alt); d$sb <- w$s[d$id] * (d$alt == "b")
    for (v in c("z1", "z2", "z3")) d[[v]] <- w[[v]][d$id]
    f1 <- function(form, firth) suppressWarnings(suppressMessages(mnl(form, data = d, id = "id", alt = "alt", firth = firth)))
    full <- chosen ~ sb | z1 + z2 + z3; rest <- chosen ~ 1 | z1 + z2 + z3
    ff <- f1(full, TRUE); fr <- f1(rest, TRUE); mf <- f1(full, FALSE); mr <- f1(rest, FALSE)
    sm <- summary(ff)$coefficients
    data.frame(p_wald_firth = sm["sb", "Pr(>|z|)"], p_lr_firth = suppressWarnings(lr_test(fr, ff))$p.value,
               p_lr_ml = suppressWarnings(lr_test(mr, mf))$p.value, separated = length(ff$separation) > 0L,
               sep_s = "sb" %in% ff$separation, any_s = any(w$s[w$y == "b"] == 1),
               zero_cell = any(table(factor(w$y == "b", c(FALSE, TRUE)), factor(w$s, 0:1)) == 0),
               conv = isTRUE(ff$converged), failed = FALSE, seed = seed)
  }, error = function(e) data.frame(failed = TRUE, msg = conditionMessage(e), seed = seed))
  out
}
rainey_cells <- lapply(c(2L, 3L), function(J) lapply(c(0, 2, -2), function(b) list(id = sprintf("F4 Rainey, J = %d, beta_s = %+g", J, b), J = J, beta_s = b)))
rainey_cells <- unlist(rainey_cells, recursive = FALSE)
## F5 (round 7): the same design with 40 s = 1 cases, where the tests can differ
rainey_cells <- c(rainey_cells, unlist(lapply(c(2L, 3L), function(J) lapply(c(0, 2, -2), function(b)
  list(id = sprintf("F5 Rainey design, 40 s = 1 cases, J = %d, beta_s = %+g", J, b), J = J, beta_s = b, ns = 40L))), recursive = FALSE))
## F6: Rainey's design 911, with the archive's ML LR rate (2,500 samples per cell)
rainey_cells <- c(rainey_cells, lapply(list(c(0, 0.0624), c(2, 0.4460), c(-2, 0.4336)), function(v)
  list(id = sprintf("F6 Rainey design 911, J = 2, beta_s = %+g", v[1]), J = 2L, beta_s = v[1], dgp911 = TRUE, his_lr_ml = v[2])))

cl <- parallel::makeCluster(W)                              # stopped after the last cell (on.exit does not act at the top level of Rscript)
invisible(parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(polychoice))))
parallel::clusterExport(cl, c("gen_mnl", "gen_ro", "ro_truth", "one", "gen_rainey", "one_rainey", "X911", "gen_rainey911"))
res <- list()
for (k in seq_along(cells)) {
  cell <- cells[[k]]
  t0 <- Sys.time()
  rows <- parallel::parLapply(cl, 20261105L + 100000L * k + seq_len(R), function(s, cell) one(cell, s), cell = cell)
  failed <- vapply(rows, function(x) isTRUE(x$failed[1]), logical(1))
  x <- do.call(rbind, rows[!failed])
  res[[cell$id]] <- list(cell = cell, rows = x, failed = rows[failed])
  saveRDS(res, rds)                                         # checkpoint after every cell
  reps <- x[!duplicated(x$seed), ]
  say("%s: %d of %d replications (%d failed: %s); ML separated in %.3f; penalized not converged %d; %.2f s per replication (wall %.0f s)",
      cell$id, nrow(reps), R, sum(failed), if (any(failed)) paste(unique(vapply(rows[failed], function(z) substr(z$msg[1], 1, 60), "")), collapse = " / ") else "-",
      mean(reps$separated), sum(!reps$conv), mean(reps$secs), as.numeric(Sys.time() - t0, units = "secs"))
  say("  penalized LR test of the null term: size at 0.05 %.3f (MC se %.3f)", mean(reps$p_null < 0.05), sqrt(0.05 * 0.95 / nrow(reps)))
  for (p in unique(x$par)) {
    xp <- x[x$par == p, ]; e <- xp$est - xp$truth; b <- mean(e); mcse <- stats::sd(xp$est) / sqrt(nrow(xp))
    ok <- xp$ml_ok; eml <- xp$ml[ok] - xp$truth[ok]
    say("  %-14s truth %+.2f | penalized bias %+.4f (%+.1f MC se) RMSE %.4f MAD %.4f | ML (%d ok) bias %+.4f RMSE %.4f | coverage profile %.3f Wald %.3f",
        p, xp$truth[1], b, b / mcse, sqrt(mean(e^2)), stats::median(abs(e)), sum(ok), if (any(ok)) mean(eml) else NA, if (any(ok)) sqrt(mean(eml^2)) else NA,
        mean(xp$cover_profile, na.rm = TRUE), mean(xp$cover_wald))
  }
}
## F4: the three tests of beta_s in Rainey's design
res_rainey <- list()
for (k in seq_along(rainey_cells)) {
  cell <- rainey_cells[[k]]
  t0 <- Sys.time()
  rows <- parallel::parLapply(cl, 20261106L + 100000L * k + seq_len(R), function(s, cell) one_rainey(cell, s), cell = cell)
  failed <- vapply(rows, function(x) isTRUE(x$failed[1]), logical(1))
  x <- do.call(rbind, rows[!failed])
  res_rainey[[cell$id]] <- list(cell = cell, rows = x, failed = rows[failed])
  saveRDS(list(main = res, rainey = res_rainey), rds)
  say("%s: %d of %d replications (%d failed); separated %.3f; Firth not converged %d; rejection at 0.05: Firth Wald %.3f, penalized LR %.3f, ML LR (suprema) %.3f (MC se at 0.05 %.3f; wall %.0f s)",
      cell$id, nrow(x), R, sum(failed), mean(x$separated), sum(!x$conv), mean(x$p_wald_firth < 0.05), mean(x$p_lr_firth < 0.05), mean(x$p_lr_ml < 0.05),
      sqrt(0.05 * 0.95 / max(1, nrow(x))), as.numeric(Sys.time() - t0, units = "secs"))
  if (isTRUE(cell$dgp911))
    say("  empty s-by-y cell %.3f; the archive's ML LR rate %.4f", mean(x$zero_cell), cell$his_lr_ml)
  ## conditional on whether beta_s itself is separated (round 6): the [sep]
  ## note in summary() rests on the Wald test's behaviour in each group
  for (g in list(list("beta_s separated", x$sep_s), list("beta_s not separated", !x$sep_s), list("no s = 1 case chooses b", !x$any_s)))
    say("  %-24s %4d samples: Firth Wald %s, penalized LR %s, ML LR %s", g[[1]], sum(g[[2]]),
        if (any(g[[2]])) sprintf("%.3f", mean(x$p_wald_firth[g[[2]]] < 0.05)) else "-",
        if (any(g[[2]])) sprintf("%.3f", mean(x$p_lr_firth[g[[2]]] < 0.05)) else "-",
        if (any(g[[2]])) sprintf("%.3f", mean(x$p_lr_ml[g[[2]]] < 0.05)) else "-")
}
parallel::stopCluster(cl)
## the criteria
say("\nCriteria")
for (nm in names(res_rainey)) {
  x <- res_rainey[[nm]]$rows; cell <- res_rainey[[nm]]$cell
  if (cell$beta_s != 0 || !is.null(cell$ns) || isTRUE(cell$dgp911)) next   # F4 only: F5 is reported, with no criterion; F6 has (f)
  e <- mean(x$p_lr_firth < 0.05) <= 0.08 && mean(x$p_lr_ml < 0.05) <= 0.08
  say("  %-34s (e) %s [penalized LR %.3f, ML LR %.3f; Firth Wald %.3f]", nm, e, mean(x$p_lr_firth < 0.05), mean(x$p_lr_ml < 0.05), mean(x$p_wald_firth < 0.05))
}
for (nm in names(res_rainey)) {
  x <- res_rainey[[nm]]$rows; cell <- res_rainey[[nm]]$cell
  if (!isTRUE(cell$dgp911)) next
  tol <- function(p1, p2) { pb <- (p1 + p2) / 2; max(0.01, 3 * sqrt(pb * (1 - pb) * (1 / nrow(x) + 1 / 2500))) }
  w <- mean(x$p_wald_firth < 0.05, na.rm = TRUE); l <- mean(x$p_lr_ml < 0.05, na.rm = TRUE)
  f <- abs(l - cell$his_lr_ml) <= tol(l, cell$his_lr_ml)          # (f) as amended: the ML LR only
  say("  %-34s (f) %s [ML LR %.3f against %.4f; Firth Wald %.3f, no criterion]", nm, f, l, cell$his_lr_ml, w)
}
for (nm in names(res)) {
  x <- res[[nm]]$rows; cell <- res[[nm]]$cell; reps <- x[!duplicated(x$seed), ]
  a <- all(reps$conv)
  cov <- tapply(x$cover_profile, x$par, mean, na.rm = TRUE); covw <- tapply(x$cover_wald, x$par, mean)
  b <- if (cell$large) all(cov >= 0.93 & cov <= 0.97) else NA
  sz <- mean(reps$p_null < 0.05); cc <- sz >= 0.035 && sz <= 0.065
  dd <- if (cell$large) all(vapply(split(x, x$par), function(xp) { bb <- mean(xp$est - xp$truth); abs(bb) <= 3 * stats::sd(xp$est) / sqrt(nrow(xp)) || abs(bb) <= 0.05 * abs(xp$truth[1]) }, logical(1))) else NA
  wald_rule <- cell$large && any(covw < 0.92 | covw > 0.98) && isTRUE(b)
  say("  %-22s (a) %s  (b) %s  (c) %s [%.3f]  (d) %s  | Wald coverage %.3f-%.3f%s", nm, a, b, cc, sz, dd, min(covw), max(covw),
      if (wald_rule) "  -> summary() rule triggered" else "")
}
say("polychoice %s from %s", as.character(packageVersion("polychoice")), find.package("polychoice"))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n", file = log, append = TRUE)
capture.output(sessionInfo(), file = log, append = TRUE)
