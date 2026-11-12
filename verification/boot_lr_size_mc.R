## Size of boot_lr_test() for the logit against the inflated multinomial logit when the
## logit is true (review round 2, reviewer B): the replicate statistic must be
## maximized as thoroughly as the observed one, or the bootstrap p-value is too
## small. Rejection rates of p.boot (LR and raw Vuong) at the 5 and 10 percent
## levels, with the Monte Carlo standard error, and the count of replicates
## whose full refit fell below the restricted one (should be zero).
## Run from the package root:  Rscript verification/boot_lr_size_mc.R
## With B replicates the bootstrap test rejects at level a exactly
## floor(a (B + 1)) / (B + 1) under the null; B = 19 makes that 0.05 at the 5
## percent level and 0.10 at the 10 percent level, and R = 300 gives a Monte
## Carlo standard error of about 0.013 at 0.05.
## Env: BLR_REPS (default 300), BLR_B (default 19).
suppressPackageStartupMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."
log <- file.path(out_dir, "boot_lr_size_mc.log"); rds <- file.path(out_dir, "boot_lr_size_mc.rds")
cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
R <- as.integer(Sys.getenv("BLR_REPS", "300")); B <- as.integer(Sys.getenv("BLR_B", "19")); n <- 500
beta <- cbind(b = c(0.2, -1, 0.5), c = c(-0.5, 0.8, -0.3))
rows <- vector("list", R); t0 <- Sys.time()
for (r in seq_len(R)) {
  set.seed(40000 + r)
  d <- rmnl(n, beta = beta); d$z1 <- rnorm(n)
  b <- tryCatch({
    m0 <- mnl(y ~ x1 + x2, data = d)
    m1 <- suppressWarnings(suppressMessages(imnl(y ~ x1 + x2 | z1, data = d, inflate = "a")))
    suppressWarnings(suppressMessages(boot_lr_test(m0, m1, B = B)))
  }, error = function(e) { say("rep %d failed: %s", r, conditionMessage(e)); NULL })
  if (is.null(b)) next
  bs <- attr(b, "boot_statistics")
  rows[[r]] <- data.frame(rep = r, lr = b$statistic[1], p_lr = b$p.boot[1], vuong = b$statistic[2], p_vuong = b$p.boot[2],
                          B_used = ncol(bs), B_vuong = sum(is.finite(bs["vuong", ])),
                          n_negative = attr(b, "n_negative"), n_failed = attr(b, "n_failed"), n_nonconverged = attr(b, "n_nonconverged"))
  if (r %% 10L == 0L) saveRDS(do.call(rbind, rows), rds)
}
x <- do.call(rbind, rows); saveRDS(x, rds)
mcse <- function(p) sqrt(p * (1 - p) / nrow(x))
say("boot_lr_test(mnl, imnl) under the MNL null: %d of %d replications, n = %d, B = %d, %.0f s per replication",
    nrow(x), R, n, B, as.numeric(Sys.time() - t0, units = "secs") / R)
## the exact level of a bootstrap test with B' usable replicates is
## floor(a (B' + 1)) / (B' + 1); a replication with B' < B (failed replicates)
## has a lower exact level, so the mean over replications is the reference
for (a in c(0.05, 0.10)) { ex_lr <- mean(floor(a * (x$B_used + 1)) / (x$B_used + 1)); ex_v <- mean(floor(a * (x$B_vuong + 1)) / (x$B_vuong + 1))
  say("  level %.2f: LR rejects %.3f (mean exact level %.3f, MC se %.3f); Vuong rejects %.3f (mean exact level %.3f; %d replications without a Vuong p-value)",
      a, mean(x$p_lr <= a), ex_lr, mcse(ex_lr), mean(x$p_vuong <= a, na.rm = TRUE), ex_v, sum(is.na(x$p_vuong))) }
say("  usable bootstrap replicates per replication: min %d, mean %.1f (of B = %d)", min(x$B_used), mean(x$B_used), B)
say("  replicates below the restricted fit: %d in all; failed %d; non-converged %d", sum(x$n_negative), sum(x$n_failed), sum(x$n_nonconverged))
say("polychoice %s from %s", as.character(packageVersion("polychoice")), find.package("polychoice"))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n", file = log, append = TRUE)
capture.output(sessionInfo(), file = log, append = TRUE)
