## Calibration Monte Carlo for the split-panel jackknife gate of mnl(fe = ...).
## For every replication the gate statistics are stored (trend p-values per
## alternative, half-panel distance, placebo noise, corrected slopes) so that
## refusal rates can be computed for any pair of constants (kappa, alpha)
## without refitting. Panels: 150 units, J = 3, T in {6, 10, 20}; three DGPs:
## stationary (the gate should rarely refuse), a regime change in the x1
## slopes at mid-panel (it should refuse), and an omitted smooth time trend
## (it should refuse and the trend detector should be the one that fires).
## A replication counts as refused when the gate would refuse it: a residual
## trend (min p below alpha / J), a failed half-panel or placebo refit, or a
## half-panel distance above max(0.05 * size, kappa * noise). Bias is reported
## signed, for the ML and jackknife slopes on the SAME (accepted) replications.
## Run from the package root:  Rscript verification/jackknife_gate_mc.R
## Env: JK_REPS_NULL (default 200), JK_REPS_ALT (default 40), JK_SEED_BASE
## (default 5000: the selection run; the holdout uses 9000, through
## verification/jackknife_gate_holdout.R), JK_TAG (the file name stem of the
## log and results, default jackknife_gate_mc).
suppressPackageStartupMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."
tag <- Sys.getenv("JK_TAG", "jackknife_gate_mc"); seed_base <- as.integer(Sys.getenv("JK_SEED_BASE", "5000"))
log <- file.path(out_dir, paste0(tag, ".log")); rds <- file.path(out_dir, paste0(tag, ".rds"))
cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
## the grid and the shipped constants: kappa is the smallest grid value whose
## stationary refusal rate is at most 0.05 at every T, alpha = 0.01 (the rule
## of the ledger entry "Jackknife gate, second design", written before this
## design's calibration was run); JK_KAPPA overrides the shipped value
kappa_grid <- c(1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4)
kappa_env <- Sys.getenv("JK_KAPPA", ""); alpha_shipped <- 0.01
kappa_shipped <- if (nzchar(kappa_env)) as.numeric(kappa_env) else 2.5      # selected 2026-10-14 (ledger)
refused <- function(x, alpha = alpha_shipped, kappa = kappa_shipped)
  x$min_trend_p < alpha / 3 | !x$halves_ok | x$n_placebo < 3 | (is.finite(x$dist) & x$dist > pmax(0.05 * x$size, kappa * x$noise))

B <- cbind(b = c(0.3, 0.8, -0.5), c = c(-0.2, -0.6, 0.7)); sl <- c("b:x1", "b:x2", "c:x1", "c:x2"); truth <- c(0.8, -0.5, -0.6, 0.7)
sim <- function(G, T, shift = 0, trend = 0) {
  unit <- rep(seq_len(G), each = T); time <- rep(seq_len(T), G); n <- G * T
  x1 <- rnorm(n); x2 <- rnorm(n); U <- matrix(rnorm(G * 2), G)
  V <- cbind(0, cbind(1, x1, x2) %*% B + U[unit, ]); V[, -1] <- V[, -1] + shift * x1 * (time > T / 2); V[, 2] <- V[, 2] + trend * time
  P <- exp(V - apply(V, 1, max)); P <- P / rowSums(P)
  y <- apply(P, 1, function(pp) sample.int(3L, 1L, prob = pp))
  data.frame(y = factor(c("a", "b", "c")[y], levels = c("a", "b", "c")), x1 = x1, x2 = x2, unit = factor(unit))
}
R0 <- as.integer(Sys.getenv("JK_REPS_NULL", "200")); R1 <- as.integer(Sys.getenv("JK_REPS_ALT", "40")); G <- 150
say("seed base %d, %d stationary and %d alternative replications per cell; polychoice %s from %s; build %s", seed_base, R0, R1,
    as.character(utils::packageVersion("polychoice")), find.package("polychoice"), Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"))
configs <- expand.grid(T = c(6L, 10L, 20L), dgp = c("stationary", "regime", "trend"), stringsAsFactors = FALSE)
res <- list()
for (i in seq_len(nrow(configs))) {
  T <- configs$T[i]; dgp <- configs$dgp[i]; R <- if (dgp == "stationary") R0 else R1
  set.seed(seed_base + 10 * T + match(dgp, c("stationary", "regime", "trend")))
  t0 <- Sys.time()
  rows <- vector("list", R)
  for (r in seq_len(R)) {
    d <- sim(G, T, shift = if (dgp == "regime") 1.5 else 0, trend = if (dgp == "trend") 0.15 else 0)
    m <- tryCatch(suppressWarnings(suppressMessages(mnl(y ~ x1 + x2, data = d, fe = "unit"))), error = function(e) NULL)
    if (is.null(m)) next
    jk <- tryCatch(polychoice:::.pc_jackknife(m, m$design, full = TRUE), error = function(e) NULL)
    if (is.null(jk) || is.null(jk$gate)) next
    g <- jk$gate
    rows[[r]] <- data.frame(T = T, dgp = dgp, rep = r, min_trend_p = min(g$trend_p), dist = g$dist, noise = g$noise, size = g$size, corr = g$corr,
                            placebo_min = suppressWarnings(min(g$placebo, na.rm = TRUE)), placebo_max = suppressWarnings(max(g$placebo, na.rm = TRUE)),
                            n_placebo = sum(is.finite(g$placebo)), halves_ok = !is.null(jk$par),
                            t(setNames(m$coefficients[sl], paste0("ml_", sl))),
                            t(setNames(sqrt(diag(vcov(m)))[sl], paste0("se_", sl))),
                            t(setNames(if (is.null(jk$par)) rep(NA_real_, 4) else jk$par[sl], paste0("jk_", sl))), check.names = FALSE)
  }
  res[[i]] <- do.call(rbind, rows)
  saveRDS(do.call(rbind, res), rds)
  x <- res[[i]]
  ## refusal rates under the shipped constants and the detector that fires
  trend_ref <- x$min_trend_p < alpha_shipped / 3; fail_ref <- !x$halves_ok | x$n_placebo < 3
  dis_ref <- is.finite(x$dist) & x$dist > pmax(0.05 * x$size, kappa_shipped * x$noise)
  acc <- !refused(x)
  bias <- function(cols, rows) colMeans(x[rows, cols, drop = FALSE], na.rm = TRUE) - truth
  say("T = %2d %-10s: %d reps, %.1f s/rep | shipped gate: refused %.3f (trend %.3f, refit failure %.3f, disagreement %.3f) | dist/noise quantiles 50/90/99: %s",
      T, dgp, nrow(x), as.numeric(Sys.time() - t0, units = "secs") / max(1, nrow(x)), mean(!acc), mean(trend_ref), mean(fail_ref), mean(dis_ref),
      paste(sprintf("%.2f", stats::quantile(x$dist / x$noise, c(0.5, 0.9, 0.99), na.rm = TRUE)), collapse = "/"))
  ## root mean squared error and the coverage of 95 percent Wald intervals (the
  ## jackknife slope with the ML standard error, as the package reports it)
  rmse <- function(cols, rows) sqrt(colMeans(sweep(as.matrix(x[rows, cols, drop = FALSE]), 2L, truth)^2, na.rm = TRUE))
  cover <- function(cols, rows) colMeans(abs(sweep(as.matrix(x[rows, cols, drop = FALSE]), 2L, truth)) <= 1.959964 * as.matrix(x[rows, paste0("se_", sl), drop = FALSE]), na.rm = TRUE)
  if (any(acc)) say("        accepted reps: RMSE ML %s | jackknife %s ; Wald coverage ML %s | jackknife %s",
      paste(sprintf("%.3f", rmse(paste0("ml_", sl), acc)), collapse = " "), paste(sprintf("%.3f", rmse(paste0("jk_", sl), acc)), collapse = " "),
      paste(sprintf("%.2f", cover(paste0("ml_", sl), acc)), collapse = " "), paste(sprintf("%.2f", cover(paste0("jk_", sl), acc)), collapse = " "))
  say("        signed bias (b:x1 b:x2 c:x1 c:x2), accepted reps (%d): ML %s | jackknife %s ; all reps ML %s",
      sum(acc), paste(sprintf("%+.3f", bias(paste0("ml_", sl), acc)), collapse = " "), paste(sprintf("%+.3f", bias(paste0("jk_", sl), acc)), collapse = " "),
      paste(sprintf("%+.3f", bias(paste0("ml_", sl), rep(TRUE, nrow(x)))), collapse = " "))
}
all <- do.call(rbind, res); saveRDS(all, rds)
## refusal rates over a grid of constants
say("--- refusal rates by (alpha, kappa): rows T x dgp ---")
for (alpha in c(0.01, 0.001)) for (kappa in kappa_grid) {
  r <- sapply(split(all, list(all$T, all$dgp)), function(x) mean(refused(x, alpha, kappa)))
  say("alpha %.3f kappa %.2f: %s", alpha, kappa, paste(sprintf("%s=%.3f", names(r), r), collapse = "  "))
}
## the selection rule, applied to this run
st <- all[all$dgp == "stationary", ]
ok_k <- vapply(kappa_grid, function(k) all(tapply(refused(st, alpha_shipped, k), st$T, mean) <= 0.05), logical(1))
say("selection rule (alpha %.2f, stationary refusal <= 0.05 at every T): smallest kappa %s; shipped kappa %.1f: stationary refusal %s",
    alpha_shipped, if (any(ok_k)) format(kappa_grid[which(ok_k)[1L]]) else "none", kappa_shipped,
    paste(sprintf("T=%s %.3f", names(tapply(refused(st), st$T, mean)), tapply(refused(st), st$T, mean)), collapse = ", "))
say("done; statistics in %s", rds)
say("JK_REPS_NULL %d, JK_REPS_ALT %d, polychoice %s from %s", R0, R1, as.character(packageVersion("polychoice")), find.package("polychoice"))
capture.output(sessionInfo(), file = log, append = TRUE)
