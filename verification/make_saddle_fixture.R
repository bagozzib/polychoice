## Writes tests/testthat/fixtures/liou_saddle_theta.rds: the point at which the
## multistart of mnprobit() on Liou's Table 1 formula (`sanctions`) stops at a
## saddle (log-likelihood -1159.4773 on the development builds), with the
## saddle step switched off (.saddle at its cap). test-phase4d.R starts the
## fit there, so that the saddle check and the step run on every platform
## whatever path the optimizer takes (round 7, C-N3).
## Run from the package root:  Rscript verification/make_saddle_fixture.R [library]
a <- commandArgs(TRUE); if (length(a)) .libPaths(c(a[1], .libPaths()))
suppressPackageStartupMessages(library(polychoice))
data(sanctions)
f <- outcome ~ nonviolent * opposition + violent * opposition + growth + inflation + leg_election + exec_election + gdppc +
  trade + milex + democracy + human_rights + coldwar + peaceyears + I(peaceyears^2) + I(peaceyears^3)
mp <- suppressWarnings(suppressMessages(mnprobit(f, data = sanctions)))
ms <- suppressWarnings(suppressMessages(polychoice:::.pc_fit_design(mp$design, se = "analytic", .saddle = 3L)))
cat(sprintf("with the step: %.6f; stopped at the saddle: %.6f (saddle flag %s)\n", mp$loglik, ms$loglik, ms$saddle))
stopifnot(isTRUE(ms$saddle), abs(ms$loglik - (-1159.4773)) < 1e-3)
saveRDS(list(theta = unname(ms$theta), loglik = ms$loglik, made = format(Sys.time()), version = as.character(packageVersion("polychoice"))),
        "tests/testthat/fixtures/liou_saddle_theta.rds")
cat("saved\n")
