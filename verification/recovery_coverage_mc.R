## Parameter recovery and interval coverage (DESIGN.md section 6, item 4;
## review round 1, reviewer H7): mnl() and imnl() on data drawn from their own
## models by rmnl(), n = 5000, 95 percent Wald intervals from the analytic
## (inverse observed information) and the robust (sandwich) covariance.
## Targets: bias within two Monte Carlo standard errors, coverage inside the
## binomial band 0.95 +/- 2 sqrt(0.95 * 0.05 / R).
## Run from the package root:  Rscript verification/recovery_coverage_mc.R
## Env: COV_REPS_MNL (default 400), COV_REPS_IMNL (default 200).
suppressPackageStartupMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."
log <- file.path(out_dir, "recovery_coverage_mc.log"); rds <- file.path(out_dir, "recovery_coverage_mc.rds")
cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
R_mnl <- as.integer(Sys.getenv("COV_REPS_MNL", "400")); R_imnl <- as.integer(Sys.getenv("COV_REPS_IMNL", "200")); n <- 5000
beta <- cbind(b = c(0.2, -1, 0.5), c = c(-0.5, 0.8, -0.3)); delta <- c(0.5, 1)
z <- stats::qnorm(0.975)
one <- function(m, truth) {
  fn <- polychoice:::.pc_fns(m$design)
  Vr <- polychoice:::.pc_strip_vcov(polychoice:::.pc_vcov(m$theta, fn, "robust", NULL, m$scale, m$design$w))
  est <- coef(m); sa <- sqrt(diag(vcov(m))); sr <- sqrt(diag(Vr))
  data.frame(par = names(est), truth = unname(truth[names(est)]), est = unname(est), se_a = unname(sa), se_r = unname(sr),
             cover_a = abs(est - truth[names(est)]) <= z * sa, cover_r = abs(est - truth[names(est)]) <= z * sr, row.names = NULL)
}
summarise <- function(x, label, R) {
  band <- 2 * sqrt(0.95 * 0.05 / R)
  say("%s: %d replications, n = %d; coverage band [%.3f, %.3f]", label, length(unique(x$rep)), n, 0.95 - band, 0.95 + band)
  for (p in unique(x$par)) {
    xp <- x[x$par == p, ]; b <- mean(xp$est - xp$truth); mcse <- stats::sd(xp$est) / sqrt(nrow(xp))
    say("  %-18s truth %+.3f bias %+.4f (MC se %.4f, %+.1f se) | sd of estimates %.4f, mean analytic se %.4f, robust se %.4f | coverage analytic %.3f robust %.3f",
        p, xp$truth[1], b, mcse, b / mcse, stats::sd(xp$est), mean(xp$se_a), mean(xp$se_r), mean(xp$cover_a), mean(xp$cover_r))
  }
}
truth_mnl <- c(stats::setNames(beta[, "b"], paste0("b:", c("(Intercept)", "x1", "x2"))),
               stats::setNames(beta[, "c"], paste0("c:", c("(Intercept)", "x1", "x2"))))
## one seed per replication (any one can be regenerated alone); a failed
## replication is dropped and counted, not fatal
t0 <- Sys.time()
xm <- do.call(rbind, lapply(seq_len(R_mnl), function(r) {
  set.seed(20261015 + r)
  tryCatch({ d <- rmnl(n, beta = beta); cbind(rep = r, one(mnl(y ~ x1 + x2, data = d), truth_mnl)) }, error = function(e) NULL) }))
saveRDS(list(mnl = xm), rds)                                 # checkpoint after the first part
say("mnl: %.2f s per replication; %d of %d replications completed", as.numeric(Sys.time() - t0, units = "secs") / R_mnl, length(unique(xm$rep)), R_mnl)
summarise(xm, "mnl()", R_mnl)
truth_imnl <- c(truth_mnl, "infl_(Intercept)" = delta[1], infl_z1 = delta[2])
t0 <- Sys.time(); nb <- 0L
xi <- do.call(rbind, lapply(seq_len(R_imnl), function(r) {
  set.seed(20261016 + r * 7919)
  tryCatch({
    d <- rmnl(n, beta = beta, delta = delta, inflate = "a")
    m <- suppressMessages(imnl(y ~ x1 + x2 | z1, data = d, inflate = "a"))
    if (m$boundary || !m$converged) nb <<- nb + 1L
    cbind(rep = r, one(m, truth_imnl)) }, error = function(e) NULL) }))
say("imnl: %.2f s per replication; %d of %d replications completed; boundary or non-converged fits: %d", as.numeric(Sys.time() - t0, units = "secs") / R_imnl, length(unique(xi$rep)), R_imnl, nb)
summarise(xi, "imnl()", R_imnl)
saveRDS(list(mnl = xm, imnl = xi), rds)
say("polychoice %s from %s", as.character(packageVersion("polychoice")), find.package("polychoice"))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n", file = log, append = TRUE)
capture.output(sessionInfo(), file = log, append = TRUE)
