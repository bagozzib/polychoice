## Size of boot_lr_test() for two nested rologit() fits (DESIGN.md, section 1,
## Phases 4c and 4d): rankings drawn under the restricted model (the added
## attribute's coefficient 0) in five designs of 4 alternatives:
##   regular  50 complete rankings, 2 attributes plus the tested one;
##   poast    8 complete rankings, 3 attributes plus the tested one (the
##            regime of the research note's Poast replication: about one
##            choice per coefficient; with continuous attributes it almost
##            never separates, 0 of 300 data sets in the probe of 2026-11-06,
##            so the sparse design carries the separated case);
##   ties     50 rankings cut at depth 2, the top two tied in 30 percent of
##            them, independently of the utilities (the coarsening the
##            bootstrap reproduces);
##   asc      50 complete rankings with alternative-specific constants (0,
##            0.5, 0.25, -0.25) and 2 attributes plus the tested one;
##   binary   50 complete rankings, 2 attributes plus a tested binary attribute
##            of prevalence 0.1 (rarely separated: 50 rankings of 4 carry
##            enough comparisons);
##   sparse   10 rankings cut at depth 2, the top two tied in 30 percent of
##            them, 2 binary attributes of prevalence 0.2 plus a tested binary
##            attribute of prevalence 0.1 (round 6: separation and zero
##            comparisons common).
## Rejection rates at the 5 and 10 percent levels of the chi-square p-value of
## lr_test() (which compares suprema under separation) and of the bootstrap
## p-value, with the Monte Carlo standard error; the shares of data sets
## separated and with a zero comparison (one the design cannot tell apart),
## the rejection rates among the separated data sets, and the mean share of
## separated replicates. With B replicates
## the bootstrap test rejects at level a exactly floor(a (B + 1)) / (B + 1)
## under the null.
## Run from the package root:  Rscript verification/rologit_boot_size_mc.R
## Env: RBS_REPS (default 300), RBS_B (default 99).
suppressPackageStartupMessages(library(polychoice))
out_dir <- if (dir.exists("verification")) "verification" else "."
log <- file.path(out_dir, "rologit_boot_size_mc.log"); rds <- file.path(out_dir, "rologit_boot_size_mc.rds")
cat("", file = log)
say <- function(...) { msg <- sprintf(...); cat(msg, "\n", file = log, append = TRUE); cat(msg, "\n") }
R <- as.integer(Sys.getenv("RBS_REPS", "300")); B <- as.integer(Sys.getenv("RBS_B", "99"))
designs <- list(regular = list(n = 50L, K = 2L, beta = c(0.8, -0.6)), poast = list(n = 8L, K = 3L, beta = c(0.8, -0.6, 0.4)),
                ties = list(n = 50L, K = 2L, beta = c(0.8, -0.6), ties = TRUE),
                asc = list(n = 50L, K = 2L, beta = c(0.8, -0.6), asc = c(0, 0.5, 0.25, -0.25)),
                binary = list(n = 50L, K = 2L, beta = c(0.8, -0.6), binary = TRUE),
                sparse = list(n = 10L, K = 2L, beta = c(0.8, -0.6), binary = TRUE, all_binary = TRUE, ties = TRUE))
J <- 4L
draw <- function(ds) {
  n <- ds$n; K <- ds$K
  d <- data.frame(id = rep(seq_len(n), each = J), alt = factor(rep(seq_len(J), n)))
  X <- matrix(rnorm(n * J * (K + 1L)), n * J); colnames(X) <- paste0("w", seq_len(K + 1L))
  if (isTRUE(ds$all_binary)) X[, seq_len(K)] <- rbinom(n * J * K, 1L, 0.2)
  if (isTRUE(ds$binary)) X[, K + 1L] <- rbinom(n * J, 1L, 0.1)
  d <- cbind(d, X)
  u <- drop(X[, seq_len(K), drop = FALSE] %*% ds$beta) - log(-log(runif(n * J)))    # the tested attribute has coefficient 0
  if (!is.null(ds$asc)) u <- u + ds$asc[as.integer(d$alt)]
  d$rank <- ave(-u, d$id, FUN = rank)
  if (isTRUE(ds$ties)) {
    ## depth 2, and the top two tied in 30 percent of the rankings
    d$rank[d$rank > 2] <- 0
    tie <- rep(runif(n) < 0.3, each = J); d$rank[tie & d$rank == 2] <- 1
  }
  d
}
all <- list(); t0 <- Sys.time()
for (dn in names(designs)) {
  ds <- designs[[dn]]; rows <- vector("list", R)
  cst <- if (is.null(ds$asc)) "0" else character(0)        # alternative-specific constants only in the asc design
  f0 <- stats::reformulate(c(cst, paste0("w", seq_len(ds$K))), "rank"); f1 <- stats::reformulate(c(cst, paste0("w", seq_len(ds$K + 1L))), "rank")
  for (r in seq_len(R)) {
    set.seed(60000 + 1000 * match(dn, names(designs)) + r)
    d <- draw(ds)
    b <- tryCatch({
      m0 <- suppressWarnings(rologit(f0, data = d, id = "id", alt = "alt"))
      m1 <- suppressWarnings(rologit(f1, data = d, id = "id", alt = "alt"))
      s1 <- separation(m1)
      list(sep = s1$separated, n_zero = s1$n_zero, lr = lr_test(m0, m1), bt = suppressWarnings(suppressMessages(boot_lr_test(m0, m1, B = B))))
    }, error = function(e) { say("%s rep %d failed: %s", dn, r, conditionMessage(e)); NULL })
    if (is.null(b)) next
    rows[[r]] <- data.frame(design = dn, rep = r, separated = b$sep, n_zero = b$n_zero, lr = b$lr$statistic, p_chisq = b$lr$p.value, p_boot = b$bt$p.boot,
                            B_used = attr(b$bt, "B"), sep_share = attr(b$bt, "n_separated") / attr(b$bt, "B"),
                            n_negative = attr(b$bt, "n_negative"), n_failed = attr(b$bt, "n_failed"))
    if (r %% 25L == 0L) saveRDS(do.call(rbind, c(all, rows)), rds)
  }
  x <- do.call(rbind, rows); all[[dn]] <- x
  mcse <- function(p) sqrt(p * (1 - p) / nrow(x))
  say("%s: %d rankings of %d alternatives, %d + 1 attributes; %d of %d replications; data separated in %.3f; with a zero comparison in %.3f; mean share of separated replicates %.3f",
      dn, ds$n, J, ds$K, nrow(x), R, mean(x$separated), mean(x$n_zero > 0L), mean(x$sep_share))
  for (a in c(0.05, 0.10)) {
    ex <- mean(floor(a * (x$B_used + 1)) / (x$B_used + 1))
    say("  level %.2f: chi-square rejects %.3f (MC se %.3f); bootstrap rejects %.3f (exact level %.3f, MC se %.3f); among the %d separated: chi-square %.3f, bootstrap %.3f",
        a, mean(x$p_chisq <= a), mcse(a), mean(x$p_boot <= a), ex, mcse(ex),
        sum(x$separated), if (any(x$separated)) mean(x$p_chisq[x$separated] <= a) else NA_real_, if (any(x$separated)) mean(x$p_boot[x$separated] <= a) else NA_real_)
  }
  say("  replicates below the restricted fit: %d; failed: %d", sum(x$n_negative), sum(x$n_failed))
}
saveRDS(do.call(rbind, all), rds)
say("total %.0f s; polychoice %s from %s", as.numeric(Sys.time() - t0, units = "secs"), as.character(packageVersion("polychoice")), find.package("polychoice"))
cat("build:", Sys.getenv("PC_BUILD", "unrecorded (not run by run_verification.sh)"), "\n", file = log, append = TRUE)
capture.output(sessionInfo(), file = log, append = TRUE)
