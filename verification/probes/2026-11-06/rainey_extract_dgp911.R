## The covariates of the design behind Rainey's (2024) Table 2, from his
## replication archive (Harvard Dataverse, doi:10.7910/DVN/6EYRJG, CC0 1.0;
## rainey-2023-pa.zip, file id 7265498): output/sims-info.rds, the design
## matrix of power_fn_id 911 (n = 1000; s = 1 for five cases; three controls).
## The file holds the same X at every beta_s of that design (checked below).
## Writes verification/data/rainey2024_dgp911.csv for firth_mc.R's F6 cells.
## Reads the archive's data files only; runs none of its code.
##   Rscript rainey_extract_dgp911.R <archive zip> <archive output folder> <out csv>
a <- commandArgs(TRUE); zip <- a[1]; od <- a[2]; out <- a[3]
si <- readRDS(file.path(od, "scenario-info.rds"))
ids <- si$scenario_id[si$power_fn_id == 911]
info <- readRDS(file.path(od, "sims-info.rds"))
Xs <- lapply(ids, function(k) as.matrix(info$design_matrix[[which(info$scenario_id == k)]]))
stopifnot(length(Xs) == 23L, all(vapply(Xs, identical, logical(1), Xs[[1]])))   # the same X at every beta_s
X <- Xs[[1]]
stopifnot(identical(dim(X), c(1000L, 5L)), all(X[, 1] == 1), sum(X[, 2]) == 5, all(X[, 2] %in% 0:1))
d <- data.frame(s = as.integer(X[, 2]), z1 = X[, 3], z2 = X[, 4], z3 = X[, 5])
md5 <- unname(tools::md5sum(c(zip, file.path(od, "sims-info.rds"))))
hdr <- c("# The covariates of design 911 of Rainey (2024), Hypothesis tests under separation,",
         "# Political Analysis 32(2), 172-185, doi:10.1017/pan.2023.28: the design behind his Table 2.",
         "# Source: his replication archive, Harvard Dataverse doi:10.7910/DVN/6EYRJG (CC0 1.0),",
         sprintf("# rainey-2023-pa.zip (md5 %s), output/sims-info.rds (md5 %s),", md5[1], md5[2]),
         "# the design_matrix of power_fn_id 911 (identical at all 23 values of beta_s).",
         "# Columns: s (the binary covariate), z1-z3 (the controls, coefficient 1). Constant -4.1.",
         "# Written by verification/probes/2026-11-06/rainey_extract_dgp911.R.")
writeLines(hdr, out)
suppressWarnings(utils::write.table(d, out, sep = ",", row.names = FALSE, append = TRUE, quote = FALSE))   # the header line after the comments
back <- utils::read.csv(out, comment.char = "#")
stopifnot(isTRUE(all.equal(as.matrix(back), as.matrix(d), check.attributes = FALSE, tolerance = 0)))
cat(sprintf("wrote %s: %d rows; sum s %d; control sd %s; correlations %s\n", out, nrow(d), sum(d$s),
            paste(round(apply(d[, 2:4], 2, stats::sd), 3), collapse = " "),
            paste(round(stats::cor(d[, 2:4])[upper.tri(diag(3))], 3), collapse = " ")))
cat("md5:", md5, "\n")
