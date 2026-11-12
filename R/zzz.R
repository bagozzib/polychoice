## Optional texreg integration (screenreg / texreg / htmlreg). texreg is in
## Suggests; its extract() method is registered on load only when texreg is
## installed, wrapped so a registration hiccup can never block loading. The S3
## class vectors are registered with setOldClass() so S4 dispatch on the
## texreg generic reaches the shared "choicemodel" method from every estimator
## class. modelsummary users are served by the broom methods in broom.R.

.extract_choicemodel <- function(model, include.fe = FALSE, ...) {
  s <- summary(model); ct <- s$coefficients
  keep <- rep(TRUE, nrow(ct))
  if (!include.fe) keep <- keep & s$block != "fixed effect"
  nm <- rownames(ct)[keep]; bl <- s$block[keep]
  nm <- ifelse(bl == "inflation", paste0("Inflation: ", sub("^infl_", "", nm)), nm)
  gof <- c(-2 * model$loglik + 2 * model$df, -2 * model$loglik + log(model$nobs_weighted) * model$df, model$loglik, model$n)
  gofn <- c("AIC", "BIC", "Log Likelihood", "Num. cases")
  gofd <- c(TRUE, TRUE, TRUE, FALSE)
  if (model$inflated) {
    gof <- c(gof, stats::weighted.mean(model$y %in% .pc_kc(model), model$weights))
    gofn <- c(gofn, "Share in inflated alt."); gofd <- c(gofd, TRUE)
  }
  gof <- c(gof, as.numeric(isTRUE(model$converged)), as.numeric(isTRUE(model$boundary) || isTRUE(model$ill_conditioned)))
  gofn <- c(gofn, "Converged", "Boundary/ill-conditioned"); gofd <- c(gofd, FALSE, FALSE)
  ## a penalized fit: its log-likelihood and criteria are at the penalized estimates
  if (.pc_penalized(model)) { gof <- c(gof, 1); gofn <- c(gofn, paste0("Penalized (", if (isTRUE(model$firth)) "Firth" else "prior", ")")); gofd <- c(gofd, FALSE) }
  texreg::createTexreg(coef.names = nm, coef = unname(ct[keep, 1L]), se = unname(ct[keep, 2L]),
                       pvalues = unname(ct[keep, 4L]), gof.names = gofn, gof = gof, gof.decimal = gofd)
}

.onLoad <- function(libname, pkgname) {
  tryCatch({
    for (cl in c("mnl", "imnl", "nlogit", "hev", "mnprobit", "remnl", "rologit")) methods::setOldClass(c(cl, "choicemodel"))
  }, error = function(e) invisible())
  if (!requireNamespace("texreg", quietly = TRUE)) return(invisible())
  tryCatch(methods::setMethod(texreg::extract, "choicemodel", .extract_choicemodel), error = function(e) invisible())
}
