## Parameter recovery and standard-error calibration of mnprobit() on the
## design of the "Beyond IIA" vignette's example: J = 3, an alternative-
## specific covariate w1 (gamma = 1), two case-specific covariates, error-
## difference covariance [1, 0.6; 0.6, 1.5], n = 1500 and n = 3000 (with three
## alternatives the probabilities are exact bivariate normal integrals, so
## any bias is the estimator's own, not simulation error). The n = 3000 cell
## (added 2026-11-06 at the author's decision, before its run) shows how the
## biases at n = 1500 change with n.
## Reports per sample size and parameter the mean bias and its Monte Carlo
## standard error, the mean standard error against the Monte Carlo standard
## deviation, and the coverage of the 95 percent Wald intervals (MNP_REPS
## replications per size, default 200; seeds 7000 + r at each size).
## Run from the package root: Rscript verification/mnprobit_recovery_mc.R
## (writes verification/mnprobit_recovery_mc.log and _results.rds, a list by
## sample size)
suppressMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."   # run from the package root or from this folder
log <- file.path(out_dir, "mnprobit_recovery_mc.log"); cat("", file = log)
say <- function(...) { msg <- paste0(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
R <- as.integer(Sys.getenv("MNP_REPS", "200")); sizes <- c(1500L, 3000L)
beta <- cbind(b = c(0.3, -0.6, 0.5), c = c(-0.2, 0.4, -0.3)); sigma <- matrix(c(1, 0.6, 0.6, 1.5), 2)
truth <- c(`b:(Intercept)` = 0.3, `b:x1` = -0.6, `b:x2` = 0.5, `c:(Intercept)` = -0.2, `c:x1` = 0.4, `c:x2` = -0.3, w1 = 1, `cov:c,b` = 0.6, `var:c` = 1.5)
t0 <- proc.time()[["elapsed"]]
all_res <- list()
for (n in sizes) {
  res <- vector("list", R)
  for (r in seq_len(R)) {
    set.seed(7000 + r)
    res[[r]] <- tryCatch({
      d <- rmnp(n, beta = beta, sigma = sigma, gamma = 1)
      m <- suppressWarnings(suppressMessages(mnprobit(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")))
      list(est = coef(m), se = sqrt(diag(vcov(m))), converged = isTRUE(m$converged), boundary = isTRUE(m$boundary))
    }, error = function(e) { say("n = ", n, ", rep ", r, " failed: ", conditionMessage(e)); NULL })
    if (r %% 25L == 0L) {
      say("n = ", n, ", rep ", r, " elapsed ", round(proc.time()[["elapsed"]] - t0), " s")
      all_res[[as.character(n)]] <- res; saveRDS(all_res, file.path(out_dir, "mnprobit_recovery_mc_results.rds"))
    }
  }
  all_res[[as.character(n)]] <- res
  saveRDS(all_res, file.path(out_dir, "mnprobit_recovery_mc_results.rds"))
}
for (n in sizes) {
  ok <- Filter(Negate(is.null), all_res[[as.character(n)]])
  nm <- names(ok[[1]]$est)
  say("\nmnprobit recovery: ", length(ok), " of ", R, " replications, n = ", n, "; converged ", sum(vapply(ok, `[[`, logical(1), "converged")),
      ", boundary ", sum(vapply(ok, `[[`, logical(1), "boundary")))
  E <- do.call(rbind, lapply(ok, `[[`, "est")); S <- do.call(rbind, lapply(ok, `[[`, "se"))
  say(sprintf("%-16s %8s %8s %8s %9s %8s %8s %8s", "parameter", "truth", "mean", "bias", "bias/mcse", "MC sd", "mean se", "cover95"))
  for (j in nm) {
    tr <- if (j %in% names(truth)) truth[[j]] else NA_real_
    b <- mean(E[, j]) - tr; mcse <- sd(E[, j]) / sqrt(nrow(E))
    cov95 <- if (is.na(tr)) NA_real_ else mean(abs(E[, j] - tr) <= 1.959964 * S[, j], na.rm = TRUE)
    say(sprintf("%-16s %8.3f %8.3f %+8.3f %9.2f %8.3f %8.3f %8.3f", j, tr, mean(E[, j]), b, b / mcse, sd(E[, j]), mean(S[, j], na.rm = TRUE), cov95))
  }
}
say("seconds: ", round(proc.time()[["elapsed"]] - t0))
say("Session: ", R.version.string, "  polychoice ", as.character(packageVersion("polychoice")), " from ", find.package("polychoice"))
say("build: ", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"))
