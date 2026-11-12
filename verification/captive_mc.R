## Parameter recovery and interval coverage of imnl() with two inflated
## (captive) alternatives (DESIGN.md, Phase 3 A): data drawn by rmnl() from the
## model itself, J = 4 alternatives, captive alternatives a (the base) and c,
## outcome covariates x1, x2. Four designs: n = 2000 and 5000, with an
## exclusion restriction (the split equations on z1, not in the outcome
## equations) and without one (the split equations on x1).
## Criteria fixed before the run (as recovery_coverage_mc): 95 percent Wald
## coverage (analytic) inside 0.95 +/- 2 sqrt(0.95 * 0.05 / R), bias within
## three Monte Carlo standard errors; the shares of boundary, non-converged,
## and multimodal fits are reported.
## Run from the package root:  Rscript verification/captive_mc.R
## Env: CAP_REPS (default 200).
suppressPackageStartupMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."
log <- file.path(out_dir, "captive_mc.log"); rds <- file.path(out_dir, "captive_mc.rds")
cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
R <- as.integer(Sys.getenv("CAP_REPS", "200"))
z <- stats::qnorm(0.975)
beta <- cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3))
delta <- cbind(a = c(1.0, 0.8), c = c(1.8, -0.7))
designs <- list(list(n = 2000, excl = TRUE), list(n = 5000, excl = TRUE), list(n = 2000, excl = FALSE), list(n = 5000, excl = FALSE))
band <- 2 * sqrt(0.95 * 0.05 / R)
res <- list()
for (k in seq_along(designs)) {
  ds <- designs[[k]]; n <- ds$n
  zname <- if (ds$excl) "z1" else "x1"
  truth <- c(stats::setNames(as.vector(beta), paste0(rep(colnames(beta), each = 3), ":", c("(Intercept)", "x1", "x2"))),
             stats::setNames(as.vector(delta), paste0("infl_", rep(colnames(delta), each = 2), ":", c("(Intercept)", zname))))
  label <- sprintf("n = %d, %s", n, if (ds$excl) "split on z1 (exclusion restriction)" else "split on x1 (no exclusion restriction)")
  t0 <- Sys.time(); flags <- c(boundary = 0L, nonconv = 0L, modes = 0L, failed = 0L)
  rows <- lapply(seq_len(R), function(r) {
    set.seed(20261019 + 5000 * k + r)                      # one seed per replication
    tryCatch({
      X <- data.frame(x1 = stats::rnorm(n), x2 = stats::rnorm(n))
      Zd <- if (ds$excl) data.frame(z1 = stats::rnorm(n)) else data.frame(z1 = X$x1)
      d <- rmnl(n, beta = beta, delta = delta, inflate = c("a", "c"), x = X, z = Zd)
      f <- stats::as.formula(paste("y ~ x1 + x2 |", zname))
      m <- suppressWarnings(suppressMessages(imnl(f, data = d, inflate = c("a", "c"))))
      if (isTRUE(m$boundary)) flags["boundary"] <<- flags["boundary"] + 1L
      if (!isTRUE(m$converged)) flags["nonconv"] <<- flags["nonconv"] + 1L
      if (isTRUE(m$n_modes > 1L)) flags["modes"] <<- flags["modes"] + 1L
      est <- coef(m); se <- sqrt(diag(vcov(m)))
      data.frame(rep = r, par = names(est), truth = unname(truth[names(est)]), est = unname(est), se = unname(se),
                 cover = abs(est - truth[names(est)]) <= z * se, boundary = isTRUE(m$boundary), converged = isTRUE(m$converged), row.names = NULL)
    }, error = function(e) { flags["failed"] <<- flags["failed"] + 1L; NULL })
  })
  x <- do.call(rbind, rows)
  res[[label]] <- list(par = x, flags = flags)
  saveRDS(res, rds)                                         # checkpoint after every design
  say("%s: %d of %d replications, %.1f s per replication; boundary %d, not converged %d, several maxima %d, failed %d; coverage band [%.3f, %.3f]",
      label, length(unique(x$rep)), R, as.numeric(Sys.time() - t0, units = "secs") / R, flags["boundary"], flags["nonconv"], flags["modes"], flags["failed"],
      0.95 - band, 0.95 + band)
  for (p in unique(x$par)) {
    xp <- x[x$par == p, ]; b <- mean(xp$est - xp$truth); mcse <- stats::sd(xp$est) / sqrt(nrow(xp))
    say("  %-20s truth %+.3f bias %+.4f (%+.1f MC se) | sd %.4f, mean se %.4f | coverage %.3f",
        p, xp$truth[1], b, b / mcse, stats::sd(xp$est), mean(xp$se), mean(xp$cover))
  }
}
say("polychoice %s from %s", as.character(packageVersion("polychoice")), find.package("polychoice"))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n", file = log, append = TRUE)
capture.output(sessionInfo(), file = log, append = TRUE)
