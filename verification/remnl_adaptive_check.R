## Accuracy of the adaptive Gauss-Hermite rule of remnl() against the
## lme4::glmer oracle (J = 2: the binomial GLMM with a random intercept,
## nAGQ = 25), over panel length T and random-intercept sd. The fixed rule
## that preceded the adaptive one degraded as T grew (the integrand narrows
## below the node spacing); the adaptive rule should improve instead.
## Run from the package root:  Rscript verification/remnl_adaptive_check.R
suppressPackageStartupMessages({ library(polychoice); library(lme4) })
out_dir <- if (dir.exists("verification")) "verification" else "."
log <- file.path(out_dir, "remnl_adaptive_check.log"); cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
sim <- function(G, T, beta, sd_u, seed) {
  set.seed(seed); n <- G * T; unit <- rep(seq_len(G), each = T)
  X <- cbind(1, matrix(rnorm(n * 2), n)); u <- rnorm(G, 0, sd_u)
  V <- cbind(0, X %*% beta + u[unit]); P <- exp(V - apply(V, 1, max)); P <- P / rowSums(P)
  y <- apply(P, 1, function(pp) sample.int(2L, 1L, prob = pp))
  data.frame(y = factor(c("a", "b")[y], levels = c("a", "b")), x1 = X[, 2], x2 = X[, 3], unit = factor(unit))
}
say("remnl adaptive Gauss-Hermite vs lme4::glmer(nAGQ = 25); G = 150 units, seed 7; dll = remnl - glmer log-likelihood")
for (T in c(6L, 20L, 50L)) for (s in c(1, 2)) {
  d <- sim(150, T, cbind(b = c(0.4, 0.9, -0.6)), s, 7)
  g <- glmer(I(y == "b") ~ x1 + x2 + (1 | unit), data = d, family = binomial, nAGQ = 25L)
  llg <- as.numeric(logLik(g)); sdg <- sqrt(as.numeric(VarCorr(g)$unit)[1])
  cells <- vapply(c(8L, 12L, 40L), function(q) {
    t0 <- Sys.time(); m <- remnl(y ~ x1 + x2, data = d, unit = "unit", nodes = q)
    sprintf("%2d nodes: dll %+.1e sd %+.1e%% (%.1f s)", q, m$loglik - llg, 100 * (coef(m)[["sd:b"]] / sdg - 1),
            as.numeric(Sys.time() - t0, units = "secs"))
  }, character(1))
  say("T = %2d, sd = %.0f: glmer ll %.4f sd %.4f | %s", T, s, llg, sdg, paste(cells, collapse = " | "))
}
say("R %s, polychoice %s, lme4 %s", getRversion(), packageVersion("polychoice"), packageVersion("lme4"))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n", file = log, append = TRUE)
capture.output(sessionInfo(), file = log, append = TRUE)
