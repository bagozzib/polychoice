## Bagozzi and Marchetti (2017, PSRM 5(2), 277-294), Table 2 (Replication II,
## pres2004) and Table 1 (Replication I, pa2006): which model does each fit
## statistic favour, computed from the fitted probabilities (Vuong, PRE, AIC,
## BIC, LR)? Assertions on the same numbers: tests/testthat/test-published.R.
## Run from the package root:  Rscript verification/table2_check.R
out_dir <- if (dir.exists("verification")) "verification" else "."
sink(file.path(out_dir, "table2_check.log"), split = TRUE)
suppressMessages(library(polychoice)); data(pres2004)
f <- vote ~ educ + female + black + hispanic + pid7 + bush_values + battleground + south + gmb + agegroup + evangelical +
  evangelical:gmb + catholic + catholic:gmb + other_religion + secular + secular:gmb + ideology + mobilization
n0 <- mnl(f, data = pres2004)
n1 <- suppressMessages(imnl(vote ~ educ + female + black + hispanic + pid7 + bush_values + battleground + south + gmb + agegroup + evangelical +
             evangelical:gmb + catholic + catholic:gmb + other_religion + secular + secular:gmb + ideology + mobilization |
             educ + agegroup + mobilization + secular, data = pres2004))
y <- as.integer(pres2004$vote); n <- length(y)
P0 <- fitted(n0); P1 <- fitted(n1)
vuong_of <- function(PA, PB) { m <- log(PA[cbind(seq_len(n), y)] / PB[cbind(seq_len(n), y)]); sqrt(n) * mean(m) / sqrt(mean((m - mean(m))^2)) }
cat(sprintf("Vuong (MNL vs BIMNL; < -1.96 favours BIMNL, > 1.96 favours MNL): %.3f, two-sided p %.3f\n", vuong_of(P0, P1), 2 * pnorm(-abs(vuong_of(P0, P1)))))
## PRE as in the script: share correctly classified by the modal category, relative to the modal baseline
modal <- function(P) apply(P, 1, which.max)
pre <- function(P) { correct <- mean(modal(P) == y); pm <- max(table(y)) / n; c(correct = correct, pre = (correct - pm) / (1 - pm)) }
cat("MNL   correct / PRE:", round(pre(P0), 4), "\n")
cat("BIMNL correct / PRE:", round(pre(P1), 4), "\n")
## AIC / BIC / LR as printed
cat(sprintf("AIC: MNL %.2f, BIMNL %.2f (favours %s)\n", AIC(n0), AIC(n1), if (AIC(n1) < AIC(n0)) "BIMNL" else "MNL"))
cat(sprintf("BIC: MNL %.2f, BIMNL %.2f (favours %s)\n", BIC(n0), BIC(n1), if (BIC(n1) < BIC(n0)) "BIMNL" else "MNL"))
lr <- 2 * (n1$loglik - n0$loglik); cat(sprintf("LR = %.3f on 5 df, p = %.4f (favours %s at 0.05)\n", lr, pchisq(lr, 5, lower.tail = FALSE), if (pchisq(lr, 5, lower.tail = FALSE) < 0.05) "BIMNL" else "MNL"))
## Replication I for completeness
data(pa2006)
m0 <- mnl(vote ~ canvass + phone + democrat + republican + canvass:democrat + canvass:republican + phone:democrat + phone:republican +
            household + age + female + voted2004 + district156, data = pa2006)
m1 <- suppressMessages(imnl(vote ~ canvass + phone + democrat + republican + canvass:democrat + canvass:republican + phone:democrat + phone:republican +
             household + age + female + voted2004 + district156 | household + age + female + voted2004 + district156, data = pa2006))
y <- as.integer(pa2006$vote); n <- length(y)
cat(sprintf("\npa2006: Vuong (MNL vs BIMNL) %.3f; PRE MNL %.4f BIMNL %.4f; AIC favours %s; BIC favours %s; LR p = %.2e\n",
            vuong_of(fitted(m0), fitted(m1)), pre(fitted(m0))["pre"], pre(fitted(m1))["pre"],
            if (AIC(m1) < AIC(m0)) "BIMNL" else "MNL", if (BIC(m1) < BIC(m0)) "BIMNL" else "MNL", pchisq(2 * (m1$loglik - m0$loglik), 6, lower.tail = FALSE)))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n")
cat("\n"); print(sessionInfo())
sink()
