## Fok-Paap Q1: the covariance of theta_ML - theta_CL. Under the MNL the
## composite-likelihood estimator is nearly efficient, so the variance
## difference V_CL - V_ML (F&P eq. 18) is a small difference of two separately
## estimated matrices. This check compares, over Q1_REPS replications of the
## IIA size-study design (J = 3, four case-specific covariates, n = 1000):
##   * the Monte Carlo covariance of theta_ML - theta_CL with the averages of
##     the variance difference and of the joint (influence-function) estimate;
##   * the generalized eigenvalues of V_CL against V_ML, per replication and
##     on average (how close CL is to efficient, and how noisy the sandwich is);
##   * the share of negative statistics and the rejection rates of both forms.
## Run from the package root: Rscript verification/q1_covariance_check.R
## (writes verification/q1_covariance_check.log)
suppressMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."   # run from the package root or from this folder
log <- file.path(out_dir, "q1_covariance_check.log"); cat("", file = log)
say <- function(...) { msg <- paste0(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
R <- as.integer(Sys.getenv("Q1_REPS", "1000")); n <- 1000L
beta <- cbind(b = c(0.5, -1, 0.8, 0.4, -0.6), c = c(-0.3, 0.6, -0.5, 0.7, 0.2))
ge <- function(A, B) { Li <- solve(t(chol(B))); eigen(Li %*% A %*% t(Li), symmetric = TRUE, only.values = TRUE)$values }
ev <- function(A) eigen((A + t(A)) / 2, symmetric = TRUE, only.values = TRUE)$values
t0 <- proc.time()[["elapsed"]]
rds <- file.path(out_dir, "q1_covariance_check.rds")
out <- lapply(seq_len(R), function(r) tryCatch({
  set.seed(5000 + r)
  d <- rmnl(n, beta = beta)
  m <- mnl(y ~ x1 + x2 + x3 + x4, data = d)
  fj <- fp_tests(m, tests = "cl"); fd <- fp_tests(m, tests = "cl", vcov = "difference")
  det <- attr(fd, "details")
  D <- m$design; fn <- polychoice:::.pc_fns(D)
  V_ml <- polychoice:::.pc_strip_vcov(polychoice:::.pc_vcov(m$theta, fn, "analytic", NULL, m$scale, D$w))
  cl <- polychoice:::.pc_cl_fit(polychoice:::.pc_pairs(D, "all"), m$theta, D$n, D$K)
  Hinv <- solve(cl$H)
  Psi <- fn$scores(m$theta) %*% V_ml - cl$scores %*% Hinv
  list(dl = m$theta - det$theta_cl, V_cl = det$vcov_cl, V_ml = V_ml, V_j = crossprod(Psi),
       g = range(ge(det$vcov_cl, V_ml)),
       stat = c(joint = fj$statistic, difference = fd$statistic), p = c(joint = fj$p.value, difference = fd$p.value), psd = fd$psd)
}, error = function(e) { say("rep ", r, " failed: ", conditionMessage(e)); NULL }))
failed <- vapply(out, is.null, logical(1)); out <- out[!failed]; R <- length(out)
saveRDS(out, rds)
if (any(failed)) say(sum(failed), " replications failed and are left out")
DL <- do.call(rbind, lapply(out, `[[`, "dl")); K <- ncol(DL)
avg <- function(nm) Reduce(`+`, lapply(out, `[[`, nm)) / R
fmt <- function(v) paste(formatC(v, format = "e", digits = 2), collapse = " ")
say("Fok-Paap Q1 covariance check: ", R, " replications, n = ", n, ", K = ", K, " parameters, seeds 5000 + r")
say("eigenvalues of the Monte Carlo covariance of theta_ML - theta_CL: ", fmt(ev(cov(DL))))
say("eigenvalues of the average joint estimate:                      ", fmt(ev(avg("V_j"))))
say("eigenvalues of the average variance difference V_CL - V_ML:     ", fmt(ev(avg("V_cl") - avg("V_ml"))))
G <- do.call(rbind, lapply(out, `[[`, "g"))
say(sprintf("generalized eigenvalues of V_CL against V_ML: on average %.3f to %.3f; per replication, median smallest %.3f, median largest %.3f",
            min(ge(avg("V_cl"), avg("V_ml"))), max(ge(avg("V_cl"), avg("V_ml"))), median(G[, 1]), median(G[, 2])))
S <- do.call(rbind, lapply(out, `[[`, "stat")); P <- do.call(rbind, lapply(out, `[[`, "p"))
say(sprintf("variance difference positive semi-definite in %.3f of replications; negative statistics: difference %.3f, joint %.3f",
            mean(vapply(out, `[[`, logical(1), "psd")), mean(S[, "difference"] < 0), mean(S[, "joint"] < 0)))
for (a in c(0.01, 0.05, 0.10))
  say(sprintf("rejection rate at %.2f: joint %.3f, difference %.3f (MC standard error %.3f)", a, mean(P[, "joint"] < a), mean(P[, "difference"] < a), sqrt(a * (1 - a) / R)))
say(sprintf("joint statistic: mean %.2f (chi-square %d), variance %.1f (%d)", mean(S[, "joint"]), K, var(S[, "joint"]), 2L * K))
say("seconds: ", round(proc.time()[["elapsed"]] - t0))
say("Session: ", R.version.string, "  polychoice ", as.character(packageVersion("polychoice")), " from ", find.package("polychoice"))
say("build: ", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"))
capture.output(sessionInfo(), file = log, append = TRUE)
