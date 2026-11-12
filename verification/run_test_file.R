## Run one testthat file (NOT_CRAN set) against the installed package, printing the start and
## end of every test as it happens (so a hang is attributable), then the
## failure messages and per-test timings.
args <- commandArgs(TRUE)
f <- args[1]
pkg <- if (nzchar(Sys.getenv("PC_PKG"))) Sys.getenv("PC_PKG") else getwd()   # run from the package root
options(width = 150)
library(testthat)
Sys.setenv(NOT_CRAN = "true")
suppressPackageStartupMessages(library(polychoice))
cat(sprintf("##### build %s; polychoice %s from %s; %s\n", Sys.getenv("PC_BUILD", "unrecorded"), packageVersion("polychoice"), find.package("polychoice"), R.version.string))
cat(sprintf("\n##### %s  started %s\n", basename(f), format(Sys.time(), "%H:%M:%S")))
rep <- MultiReporter$new(list(LocationReporter$new(), ListReporter$new()))
res <- testthat::test_file(file.path(pkg, "tests", "testthat", f), reporter = rep,
                           package = "polychoice", load_package = "installed")
cat("\n##### failure / error / warning messages\n")
for (t in res) for (e in t$results)
  if (inherits(e, c("expectation_failure", "expectation_error", "expectation_warning")))
    cat(sprintf("[%s] %s\n  %s\n", class(e)[1], t$test, conditionMessage(e)))
df <- as.data.frame(res)
cat("\n##### per-test summary\n")
print(df[, c("test", "nb", "failed", "warning", "error", "real")], row.names = FALSE)
cat(sprintf("##### %s  finished %s  total %.1f s\n", basename(f), format(Sys.time(), "%H:%M:%S"), sum(df$real)))
cat(sprintf("##### %s  total expectations %d failed %d warnings %d errors %d skipped %d\n", basename(f), sum(df$nb), sum(df$failed), sum(df$warning), sum(df$error), sum(df$skipped)))
