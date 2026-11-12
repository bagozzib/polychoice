## Size of vuong(imnl, mnl) when the multinomial logit is true (review
## round 1, reviewer G7). The inflated model nests the MNL at the boundary
## phi = 1 (the split-equation intercept at +infinity, where the split slopes
## are not identified), so neither Vuong's normal reference nor a chi-square
## reference for the LR statistic is justified; this run measures how far off
## they are. DGP: rmnl(n = 1000, beta) with an independent z1 that plays no
## role; models mnl(y ~ x1 + x2) and imnl(y ~ x1 + x2 | z1, inflate = "a").
## Recorded per replication: the raw / AIC / BIC Vuong statistics (positive
## favours imnl), the LR statistic, and the imnl boundary flag.
## Run from the package root:  Rscript verification/vuong_size_mc.R
## Env: VUONG_REPS (default 400).
suppressPackageStartupMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."
log <- file.path(out_dir, "vuong_size_mc.log"); rds <- file.path(out_dir, "vuong_size_mc.rds")
cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
R <- as.integer(Sys.getenv("VUONG_REPS", "400")); n <- 1000
beta <- cbind(b = c(0.2, -1, 0.5), c = c(-0.5, 0.8, -0.3))
t0 <- Sys.time()
rows <- vector("list", R)
for (r in seq_len(R)) {
  set.seed(20261014 + r)                                   # one seed per replication: any one can be regenerated alone
  d <- rmnl(n, beta = beta); d$z1 <- rnorm(n)
  m0 <- tryCatch(mnl(y ~ x1 + x2, data = d), error = function(e) { say("rep %d: mnl failed: %s", r, conditionMessage(e)); NULL })
  if (is.null(m0)) next
  w <- character(0)
  m1 <- tryCatch(withCallingHandlers(suppressMessages(imnl(y ~ x1 + x2 | z1, data = d, inflate = "a")),
                                     warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }),
                 error = function(e) { say("rep %d: imnl failed: %s", r, conditionMessage(e)); NULL })
  if (is.null(m1)) next
  v <- tryCatch(vuong(m1, m0), error = function(e) { say("rep %d: vuong failed: %s", r, conditionMessage(e)); NULL })
  if (is.null(v)) next
  rows[[r]] <- data.frame(rep = r, lr = 2 * (m1$loglik - m0$loglik), v_raw = v$statistic[v$correction == "raw"],
                          v_aic = v$statistic[v$correction == "AIC"], v_bic = v$statistic[v$correction == "BIC"],
                          boundary = m1$boundary, converged = m1$converged, n_warn = length(w))
  if (r %% 50L == 0L) saveRDS(do.call(rbind, rows), rds)  # checkpoint
}
x <- do.call(rbind, rows); saveRDS(x, rds)
z <- stats::qnorm(0.975)
say("vuong(imnl, mnl) under the MNL null: %d of %d replications fitted, n = %d, %.1f s per replication", nrow(x), R, n,
    as.numeric(Sys.time() - t0, units = "secs") / R)
## vuong() returns NA when the two fits are practically the same model (the
## inflated fit at phi = 1): those replications are counted, and the rates are
## over the replications with a statistic and over all replications (an NA
## counted as no rejection, as a user would read it)
say("  statistics undefined (fits practically the same model): %d of %d replications", sum(is.na(x$v_raw)), nrow(x))
for (cn in c("v_raw", "v_aic", "v_bic")) { v <- x[[cn]]; f <- is.finite(v)
  say("  %-6s median %+.3f, IQR %.3f | reject two-sided 0.05: %.3f of defined, %.3f of all | favours imnl (one-sided 0.05): %.3f | favours mnl: %.3f", cn,
      stats::median(v[f]), stats::IQR(v[f]), mean(abs(v[f]) > z), sum(abs(v[f]) > z) / length(v), mean(v[f] > stats::qnorm(0.95)), mean(v[f] < -stats::qnorm(0.95))) }
say("  LR: mean %.3f, quantiles 50/90/95: %s; reject vs chi-square(2) %.3f, vs 50:50 chi-square(1)/(2) mixture %.3f",
    mean(x$lr), paste(sprintf("%.2f", stats::quantile(x$lr, c(0.5, 0.9, 0.95))), collapse = "/"),
    mean(stats::pchisq(x$lr, 2, lower.tail = FALSE) < 0.05),
    mean(0.5 * stats::pchisq(x$lr, 1, lower.tail = FALSE) + 0.5 * stats::pchisq(x$lr, 2, lower.tail = FALSE) < 0.05))
say("  imnl at the boundary in %.3f of replications; not converged in %.3f; replications with warnings %.3f",
    mean(x$boundary), mean(!x$converged), mean(x$n_warn > 0))
say("polychoice %s from %s", as.character(packageVersion("polychoice")), find.package("polychoice"))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n", file = log, append = TRUE)
capture.output(sessionInfo(), file = log, append = TRUE)
