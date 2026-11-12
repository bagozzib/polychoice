## Null rejection rates of the IIA battery (DESIGN.md, section 6, item 5):
## 300 replications of an MNL data-generating process with J = 3, four
## case-specific covariates, n = 1000; every statistic saved per replication.
## Small-Hsiao is recorded for one split and as the battery reports it (the
## median statistic over ten splits, referred to the chi-square); Fok-Paap Q1
## in both covariance forms (joint, the default, and the variance difference).
## Run from the package root: Rscript verification/size_check_iia.R
## (writes verification/size_check_iia.log and size_check_iia_results.csv)
out_dir <- if (dir.exists("verification")) "verification" else "."   # run from the package root or from this folder
sink(file.path(out_dir, "size_check_iia.log"), split = TRUE)
suppressMessages(library(polychoice))
R <- 300L; n <- 1000L
beta <- cbind(b = c(0.5, -1, 0.8, 0.4, -0.6), c = c(-0.3, 0.6, -0.5, 0.7, 0.2))
out <- data.frame(rep = seq_len(R), hm_diff_b = NA_real_, hm_diff_a = NA_real_, hm_joint_b = NA_real_, hm_joint_a = NA_real_,
                  sh_b = NA_real_, sh_a = NA_real_, q1_joint = NA_real_, q2_iter = NA_real_, q2_proj = NA_real_,
                  q2_twostep = NA_real_, sh10_b = NA_real_, sh10_a = NA_real_, q1_diff = NA_real_, hm_diff_b_psd = NA, q1_psd = NA, seconds = NA_real_)
t0 <- proc.time()[["elapsed"]]
for (r in seq_len(R)) {
  t1 <- proc.time()[["elapsed"]]
  set.seed(20261013 + r)                                   # one seed per replication: any one can be regenerated alone
  ok <- tryCatch({
  d <- rmnl(n, beta = beta)
  m <- mnl(y ~ x1 + x2 + x3 + x4, data = d)
  hd <- hausman_mcfadden(m, drop = c("b", "a"))
  hj <- hausman_mcfadden(m, drop = c("b", "a"), vcov = "joint")
  sh <- small_hsiao(m, drop = c("b", "a"))
  fp <- suppressWarnings(fp_tests(m))                                    # Q1 (joint) and the two-step Q2
  fpd <- suppressWarnings(fp_tests(m, tests = "cl", vcov = "difference")) # Q1, variance difference (F&P eq. 18)
  fpi <- suppressWarnings(fp_tests(m, tests = "gmm", gmm = "iterated"))
  fpp <- suppressWarnings(fp_tests(m, tests = "gmm", gmm = "projected"))
  out[r, 2:10] <- c(hd$p.value[1], hd$p.value[2], hj$p.value[1], hj$p.value[2], sh$p.value[1], sh$p.value[2], fp$p.value[1], fpi$p.value[1], fpp$p.value[1])
  out$q2_twostep[r] <- fp$p.value[2]
  s10 <- attr(small_hsiao(m, drop = c("b", "a"), splits = 10L), "summary")
  out$sh10_b[r] <- stats::pchisq(s10$median_statistic[1], s10$df[1], lower.tail = FALSE)
  out$sh10_a[r] <- stats::pchisq(s10$median_statistic[2], s10$df[2], lower.tail = FALSE)
  out$q1_diff[r] <- fpd$p.value[1]
  out$hm_diff_b_psd[r] <- hd$psd[1]; out$q1_psd[r] <- fpd$psd[1]
  TRUE }, error = function(e) { cat("rep", r, "failed:", conditionMessage(e), "\n"); FALSE })
  out$seconds[r] <- proc.time()[["elapsed"]] - t1
  if (r %% 50L == 0L) { cat("rep", r, "elapsed", round(proc.time()[["elapsed"]] - t0), "s\n"); write.csv(out, file.path(out_dir, "size_check_iia_results.csv"), row.names = FALSE) }
}
write.csv(out, file.path(out_dir, "size_check_iia_results.csv"), row.names = FALSE)
cat("\nRejection rates at nominal 0.05 (MC standard error about 0.013) and 0.10:\n")
pv <- out[, c(2:10, which(names(out) %in% c("q2_twostep", "sh10_b", "sh10_a", "q1_diff")))]
cat("replications completed:", sum(stats::complete.cases(pv)), "of", R, "\n")
print(round(rbind(`0.05` = colMeans(pv < 0.05, na.rm = TRUE), `0.10` = colMeans(pv < 0.10, na.rm = TRUE)), 3))
cat("\nShare of replications with a non-PSD covariance difference: HM(diff, drop b)", mean(!out$hm_diff_b_psd), " Q1 (difference)", mean(!out$q1_psd), "\n")
cat("Mean seconds per replication:", round(mean(out$seconds), 2), "\n")
cat("Session:", R.version.string, " polychoice", as.character(packageVersion("polychoice")), "from", find.package("polychoice"), "\n")
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n")
print(sessionInfo())
sink()
