## Precompute the vignettes: knit each vignettes/<name>.Rmd.orig, with the
## installed package, to <name>.Rmd, whose code and output are then plain
## text, so that R CMD build and CRAN's checks render them without fitting the
## models (the check time; the author's decision of 2026-11-06). Rerun after
## any change to the package. verification/run_verification.sh knits them
## again on the frozen build and reports any difference from the committed
## files.
##   Rscript vignettes/precompute.R [output folder]   (from the package root;
##   default: vignettes/ itself)
a <- commandArgs(TRUE)
root <- getwd()
out <- normalizePath(if (length(a)) a[1] else "vignettes", mustWork = TRUE)
suppressPackageStartupMessages(library(polychoice))
src <- list.files("vignettes", pattern = "[.]Rmd[.]orig$")
if (!length(src)) stop("no vignettes/*.Rmd.orig sources")
setwd("vignettes")                                  # the chunks read their files relative to vignettes/
for (f in src) {
  o <- file.path(out, sub("[.]orig$", "", f))
  set.seed(1)                                       # each vignette sets its own seeds; this only fixes the start
  knitr::knit(f, output = o, quiet = TRUE, envir = new.env(parent = globalenv()))
  cat("knit ", f, " -> ", o, "\n", sep = "")
}
setwd(root)
cat("polychoice", as.character(packageVersion("polychoice")), "from", find.package("polychoice"), "\n")
