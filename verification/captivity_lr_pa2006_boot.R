## The reference of captivity_test(type = "lr") on pa2006 (DESIGN.md Phase 4e
## and its round-6 amendments; ?captivity_test and the article quote this
## run): one respondent carries half of the share's expected information there.
## A parametric bootstrap from the logit, boot_lr_test(m0, imnl(f | 1)) with
## B = 199 and set.seed(1). Reported: the share of replicate statistics at
## zero (the mixture: 0.5), their 0.90 / 0.95 / 0.99 quantiles (the mixture:
## 1.64 / 2.71 / 5.41), the rejection rates at those critical values (0.10 /
## 0.05 / 0.01), and the observed statistic's bootstrap p-value. The
## outcomes are drawn in the master, so the result does not depend on the
## number of workers. First run as a probe on a development build of round 6
## (probes/2026-11-06/r6_pa2006_chibar_bootstrap); here from round 7 so that
## it reruns on every frozen build.
## Run from the package root:  Rscript verification/captivity_lr_pa2006_boot.R
## Env: CT_WORKERS (default 6; PSOCK workers load polychoice from R_LIBS).
suppressPackageStartupMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."
log <- file.path(out_dir, "captivity_lr_pa2006_boot.log"); rds <- file.path(out_dir, "captivity_lr_pa2006_boot.rds")
cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
W <- as.integer(Sys.getenv("CT_WORKERS", "6"))
data(pa2006)
f <- vote ~ canvass + phone + democrat + republican + canvass:democrat + canvass:republican + phone:democrat + phone:republican + household + age + female + voted2004 + district156
m0 <- mnl(f, data = pa2006)
mc <- suppressMessages(imnl(vote ~ canvass + phone + democrat + republican + canvass:democrat + canvass:republican + phone:democrat + phone:republican +
                              household + age + female + voted2004 + district156 | 1, data = pa2006))
set.seed(1); t0 <- proc.time()[["elapsed"]]
b <- boot_lr_test(m0, mc, B = 199, cores = W)
capture.output(print(b), file = log, append = TRUE); print(b)
s <- attr(b, "boot_statistics")["lr", ]
crit <- stats::qchisq(c(0.8, 0.9, 0.98), 1)
say("%d replicates (%.1f min, %d workers): share of LR = 0 (below 1e-6) %.3f (the mixture: 0.5)", length(s), (proc.time()[["elapsed"]] - t0) / 60, W, mean(s < 1e-6))
say("quantiles 0.90 / 0.95 / 0.99: %.2f / %.2f / %.2f (the mixture: %.2f / %.2f / %.2f)", stats::quantile(s, 0.9), stats::quantile(s, 0.95), stats::quantile(s, 0.99), crit[1], crit[2], crit[3])
say("rejection at the mixture's critical values: %.3f / %.3f / %.3f (nominal 0.10 / 0.05 / 0.01)", mean(s > crit[1]), mean(s > crit[2]), mean(s > crit[3]))
say("observed LR %.3f; asymptotic (mixture) p %.3g; bootstrap p %.4f; replicates with a full refit below the restricted one: %d",
    b$statistic[1], b$p.value[1], b$p.boot[1], attr(b, "n_negative"))
saveRDS(b, rds)
say("polychoice %s from %s", as.character(packageVersion("polychoice")), find.package("polychoice"))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n", file = log, append = TRUE)
