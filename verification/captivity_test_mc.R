## Size and power of captivity_test() (DESIGN.md, Phase 3 A; the likelihood-
## ratio form, Phase 4e). Size: data from the multinomial logit, J = 3 and 4,
## two case-specific covariates, n = 500 and 2000. Power: data from imnl()
## with baseline captivity (alternative a) of about 15 percent that depends on
## a covariate z1, at n = 500 and 2000, the test of alternative a.
## The score form (the default): size at 0.05 and 0.01 for each alternative's
## one-sided z and for the joint chi-square, with the smallest fitted
## probability of the data. The first run (build e4d451767e7a) also studied an
## outer-product form, which sized far worse and was removed (ledger).
## The likelihood-ratio form (type = "lr"), on the same data sets: size at
## 0.10, 0.05 and 0.01 for each alternative, the share of statistics at zero
## (the mixture puts half its mass there), and the inflated fits that did not
## converge; power for alternative a.
## Each replication sets its own seed (20261019 + 10000 k + r), so the results
## do not depend on the number of workers, and the score form's are those of
## the earlier serial runs. Round 6 added two size designs with small fitted
## probabilities (the coefficients tripled; a rare alternative), and logs the
## LR form's raw statistics, statuses and warnings, and the largest single
## case's share of the information.
## Run from the package root:  Rscript verification/captivity_test_mc.R
## Env: CT_REPS_SIZE (default 1000), CT_REPS_POWER (default 500), CT_WORKERS
## (default 6; PSOCK workers load polychoice from R_LIBS).
suppressPackageStartupMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."
log <- file.path(out_dir, "captivity_test_mc.log"); rds <- file.path(out_dir, "captivity_test_mc.rds")
cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
R_size <- as.integer(Sys.getenv("CT_REPS_SIZE", "1000")); R_pow <- as.integer(Sys.getenv("CT_REPS_POWER", "500"))
W <- as.integer(Sys.getenv("CT_WORKERS", "6"))
B3 <- cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5))
B4 <- cbind(B3, d = c(0.2, 0.4, 0.3))
one <- function(d, alternatives = NULL) {
  m <- mnl(y ~ x1 + x2, data = d)
  e <- captivity_test(m, alternatives = alternatives)
  nw <- 0L
  l <- withCallingHandlers(captivity_test(m, alternatives = alternatives, type = "lr"), warning = function(w) { nw <<- nw + 1L; invokeRestart("muffleWarning") })
  fl <- attr(l, "fits")
  list(z_e = e$statistic, p_e = e$p.value, minP = min(fitted(m)), test = e$test,
       lr = l$statistic, p_lr = l$p.value, conv_lr = fl$converged, test_lr = l$test,
       raw_lr = 2 * (fl$loglik - m$loglik), status = fl$status, top = fl$top_share, n_warn = nw)
}
cl <- parallel::makeCluster(W)                              # stopped at the end (on.exit does not act at the top level of Rscript)
invisible(parallel::clusterEvalQ(cl, suppressPackageStartupMessages(library(polychoice))))
run <- function(R, k, gen, alternatives = NULL)
  Filter(Negate(is.null), parallel::parLapply(cl, seq_len(R), function(r, k, gen, alternatives, one) {
    set.seed(20261019 + 10000 * k + r); tryCatch(one(gen(), alternatives), error = function(e) NULL) }, k, gen, alternatives, one))
res <- list(); band <- function(a, R) a + c(-2, 2) * sqrt(a * (1 - a) / R)
## a size design: the score and likelihood-ratio forms on the same data sets.
## For the LR form also the raw statistic 2 (l_inflated - l_logit) before the
## floor, the refits that ended below the logit (raw below -1e-6), the
## degenerate ones (p-value NA), the warnings raised, and the largest single
## case's share of the information (round 6, H-N6, B-N5, D-new-1)
size_design <- function(label, k, gen) {
  t0 <- Sys.time()
  rows <- run(R_size, k, gen)
  tests <- rows[[1L]]$test; tests_lr <- rows[[1L]]$test_lr
  Pe <- do.call(rbind, lapply(rows, `[[`, "p_e")); minP <- vapply(rows, `[[`, numeric(1), "minP")
  Pl <- do.call(rbind, lapply(rows, `[[`, "p_lr")); Ll <- do.call(rbind, lapply(rows, `[[`, "lr")); Cl <- do.call(rbind, lapply(rows, `[[`, "conv_lr"))
  Rl <- do.call(rbind, lapply(rows, `[[`, "raw_lr")); St <- do.call(rbind, lapply(rows, `[[`, "status")); Tp <- do.call(rbind, lapply(rows, `[[`, "top"))
  nw <- vapply(rows, `[[`, integer(1), "n_warn")
  res[[label]] <<- list(p_expected = Pe, z_expected = do.call(rbind, lapply(rows, `[[`, "z_e")), minP = minP, p_lr = Pl, lr = Ll, conv_lr = Cl,
                        raw_lr = Rl, status = St, top_share = Tp, n_warn = nw)
  saveRDS(res, rds)
  b05 <- band(0.05, length(rows)); b10 <- band(0.10, length(rows)); b01 <- band(0.01, length(rows))
  say("%s: %d of %d replications, %.2f s each (%d workers); smallest fitted probability median %.4f (min %.5f); bands at 0.10 [%.3f, %.3f], 0.05 [%.3f, %.3f], 0.01 [%.3f, %.3f]",
      label, length(rows), R_size, as.numeric(Sys.time() - t0, units = "secs") / R_size, W, stats::median(minP), min(minP), b10[1], b10[2], b05[1], b05[2], b01[1], b01[2])
  for (j in seq_along(tests))
    say("  score %-12s rejection at 0.05 / 0.01: %.3f / %.3f", tests[j], mean(Pe[, j] < 0.05, na.rm = TRUE), mean(Pe[, j] < 0.01, na.rm = TRUE))
  for (j in seq_along(tests_lr))
    say("  LR    %-12s rejection at 0.10 / 0.05 / 0.01: %.3f / %.3f / %.3f; statistic 0 in %.3f; not converged %d; below the logit %d; degenerate %d; largest case's information share median %.3f (max %.3f)",
        tests_lr[j], mean(Pl[, j] < 0.10, na.rm = TRUE), mean(Pl[, j] < 0.05, na.rm = TRUE), mean(Pl[, j] < 0.01, na.rm = TRUE),
        mean(Ll[, j] == 0, na.rm = TRUE), sum(!Cl[, j], na.rm = TRUE), sum(Rl[, j] < -1e-6, na.rm = TRUE), sum(St[, j] == "degenerate", na.rm = TRUE),
        stats::median(Tp[, j], na.rm = TRUE), max(Tp[, j], na.rm = TRUE))
  say("  replications with a warning from the LR form: %d", sum(nw > 0L))
}
k <- 0L
for (B in list(B3, B4)) for (n in c(500, 2000)) {
  k <- k + 1L
  size_design(sprintf("size: J = %d, n = %d", ncol(B) + 1L, n), k, local({ B <- B; n <- n; function() rmnl(n, beta = B) }))
}
for (n in c(500, 2000)) {
  k <- k + 1L; label <- sprintf("power: J = 3, n = %d, baseline captivity about 15 percent", n)
  t0 <- Sys.time()
  rows <- run(R_pow, k, local({ n <- n; B3 <- B3; function() rmnl(n, beta = B3, delta = c(stats::qlogis(0.85), 0.5), inflate = "a") }), alternatives = "a")
  Pe <- do.call(rbind, lapply(rows, `[[`, "p_e")); Pl <- do.call(rbind, lapply(rows, `[[`, "p_lr"))
  res[[label]] <- list(p_expected = Pe, p_lr = Pl)
  saveRDS(res, rds)
  say("%s: %d of %d replications, %.2f s each; rejection of captive: a at 0.05: score %.3f, LR %.3f",
      label, length(rows), R_pow, as.numeric(Sys.time() - t0, units = "secs") / R_pow, mean(Pe[, 1L] < 0.05), mean(Pl[, 1L] < 0.05, na.rm = TRUE))
}
## small fitted probabilities (round 6, D-new-1): the coefficients tripled,
## and a rare alternative (about 3 percent): few cases then dominate the
## information of a captive share
k <- k + 1L
size_design("size, small probabilities: J = 3, n = 2000, coefficients x 3", k, local({ B3 <- B3; function() rmnl(2000, beta = 3 * B3) }))
k <- k + 1L
size_design("size, rare alternative: J = 3, n = 2000, c chosen by about 3 percent", k,
            local({ B3 <- B3; function() rmnl(2000, beta = cbind(b = B3[, "b"], c = c(-3.5, 0.6, -0.5))) }))
parallel::stopCluster(cl)
say("polychoice %s from %s", as.character(packageVersion("polychoice")), find.package("polychoice"))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n", file = log, append = TRUE)
capture.output(sessionInfo(), file = log, append = TRUE)
