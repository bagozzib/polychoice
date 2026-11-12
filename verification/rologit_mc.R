## Parameter recovery and interval coverage of rologit(), and the size of
## rank_test() under the model (DESIGN.md, Phase 3 B). Rankings drawn by
## rrologit() from the rank-ordered logit itself: J = 4 and 5 alternatives, an
## alternative-specific covariate w1 (generic coefficient) and two
## case-specific covariates; complete rankings, partial rankings (the first two
## positions), and partial rankings with ties (each position tied with the one
## above with probability 0.25, the first three positions). n = 500 rankings.
## Criteria fixed before the run: 95 percent Wald coverage (analytic and
## robust) inside the binomial band 0.95 +/- 2 sqrt(0.95 * 0.05 / R), bias
## within three Monte Carlo standard errors; size of rank_test(after = 1) at
## 0.05 inside 0.05 +/- 2 sqrt(0.05 * 0.95 / R) for the LR and the Wald forms.
## Run from the package root:  Rscript verification/rologit_mc.R
## Env: RO_REPS (default 300).
suppressPackageStartupMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."
log <- file.path(out_dir, "rologit_mc.log"); rds <- file.path(out_dir, "rologit_mc.rds")
cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
R <- as.integer(Sys.getenv("RO_REPS", "300")); n <- 500
z <- stats::qnorm(0.975)
B4 <- cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3))
B5 <- cbind(B4, e = c(-0.4, -0.2, 0.6))
designs <- list(
  list(label = "J = 4, complete rankings", beta = B4, depth = NULL, ties = 0),
  list(label = "J = 4, first two positions", beta = B4, depth = 2, ties = 0),
  list(label = "J = 4, first three positions, ties 0.25", beta = B4, depth = 3, ties = 0.25),
  list(label = "J = 5, complete rankings", beta = B5, depth = NULL, ties = 0),
  list(label = "J = 5, first three positions, ties 0.25", beta = B5, depth = 3, ties = 0.25))
gamma <- 0.7
truth_of <- function(beta) {
  tr <- unlist(lapply(colnames(beta), function(a) stats::setNames(beta[, a], paste0(a, ":", c("(Intercept)", "x1", "x2")))))
  c(tr, w1 = gamma)
}
one <- function(m, truth) {
  fn <- polychoice:::.pc_fns(m$design)
  Vr <- polychoice:::.pc_strip_vcov(polychoice:::.pc_vcov(m$theta, fn, "robust", NULL, m$scale, m$design$w))
  est <- coef(m); sa <- sqrt(diag(vcov(m))); sr <- sqrt(diag(Vr))
  data.frame(par = names(est), truth = unname(truth[names(est)]), est = unname(est), se_a = unname(sa), se_r = unname(sr),
             cover_a = abs(est - truth[names(est)]) <= z * sa, cover_r = abs(est - truth[names(est)]) <= z * sr, row.names = NULL)
}
band_c <- 2 * sqrt(0.95 * 0.05 / R); band_s <- 2 * sqrt(0.05 * 0.95 / R)
res <- list()
for (k in seq_along(designs)) {
  ds <- designs[[k]]; truth <- truth_of(ds$beta)
  t0 <- Sys.time(); nfail <- 0L; ndrop <- 0L
  rows <- lapply(seq_len(R), function(r) {
    set.seed(20261019 + 1000 * k + r)                      # one seed per replication
    tryCatch({
      d <- rrologit(n, beta = ds$beta, gamma = gamma, depth = ds$depth, ties = ds$ties)
      msg <- character(0)
      m <- withCallingHandlers(rologit(rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt"),
                               message = function(cnd) { msg <<- c(msg, conditionMessage(cnd)); invokeRestart("muffleMessage") })
      if (length(msg)) ndrop <<- ndrop + 1L
      rt <- rank_test(m, after = 1)
      list(par = cbind(rep = r, one(m, truth)), rt = data.frame(rep = r, lr = rt$statistic[1], wald = rt$statistic[2], df = rt$df[1],
                                                                 p_lr = rt$p.value[1], p_wald = rt$p.value[2]),
           conv = isTRUE(m$converged), boundary = isTRUE(m$boundary))
    }, error = function(e) { nfail <<- nfail + 1L; NULL })
  })
  rows <- Filter(Negate(is.null), rows)
  x <- do.call(rbind, lapply(rows, `[[`, "par")); rt <- do.call(rbind, lapply(rows, `[[`, "rt"))
  res[[ds$label]] <- list(par = x, rank_test = rt, n_failed = nfail)
  saveRDS(res, rds)                                         # checkpoint after every design
  say("%s: %d of %d replications (%d failed, %d with dropped cases), n = %d, %.2f s per replication; not converged %d, boundary %d",
      ds$label, length(rows), R, nfail, ndrop, n, as.numeric(Sys.time() - t0, units = "secs") / R,
      sum(!vapply(rows, `[[`, logical(1), "conv")), sum(vapply(rows, `[[`, logical(1), "boundary")))
  say("  coverage band [%.3f, %.3f]; size band [%.3f, %.3f]", 0.95 - band_c, 0.95 + band_c, 0.05 - band_s, 0.05 + band_s)
  for (p in unique(x$par)) {
    xp <- x[x$par == p, ]; b <- mean(xp$est - xp$truth); mcse <- stats::sd(xp$est) / sqrt(nrow(xp))
    say("  %-16s truth %+.3f bias %+.4f (%+.1f MC se) | sd %.4f, mean se analytic %.4f robust %.4f | coverage analytic %.3f robust %.3f",
        p, xp$truth[1], b, b / mcse, stats::sd(xp$est), mean(xp$se_a), mean(xp$se_r), mean(xp$cover_a), mean(xp$cover_r))
  }
  say("  rank_test(after = 1), df %d: rejection at 0.05 LR %.3f, Wald %.3f; mean LR %.2f (df %d)",
      rt$df[1], mean(rt$p_lr < 0.05), mean(rt$p_wald < 0.05), mean(rt$lr), rt$df[1])
}
say("polychoice %s from %s", as.character(packageVersion("polychoice")), find.package("polychoice"))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n", file = log, append = TRUE)
capture.output(sessionInfo(), file = log, append = TRUE)
