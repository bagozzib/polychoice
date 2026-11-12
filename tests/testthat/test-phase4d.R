## the round-5 rules of DESIGN.md, section 1, Phase 4d not covered in
## test-suprema.R and test-prior.R: the boundary step one parameter at a time,
## the saddle probe's point, the [sep] marks, texreg on rankings, the score
## tests' separation warning, the probit's limit

test_that("the boundary step tries the parameters one at a time when moving them together lowers the log-likelihood", {
  ## two structural log-parameters: the first has an interior optimum at 4,
  ## the second runs off (the objective falls towards +infinity in it)
  D <- list(idx = list(str = 1:2), blocks = c("scale", "scale"), par_names = c("scale:a", "scale:b"))
  bd <- function(th) pmin(pmax(th, -8), 8)
  fn <- list(nll = function(th) { th <- bd(th); (th[1] - 4)^2 + exp(-th[2]) },
             grad = function(th) { th <- bd(th); c(2 * (th[1] - 4), -exp(-th[2])) }, bound = bd)
  op <- list(par = c(4, 5), value = exp(-5), converged = FALSE, max_grad = exp(-5))
  r <- .pc_boundary_step(fn, D, op, c(1, 1), 1000L, 1e-10)
  expect_equal(r$par, c(4, 8), tolerance = 1e-6)                     # the second at the bound, the first kept at its optimum
  expect_lt(r$value, op$value)
  expect_identical(r$boundary_fixed, 2L)
})

test_that(".pc_saddle() returns the probe point where the objective fell most", {
  nll <- function(th) th[1]^2 - th[2]^2                              # a saddle at 0, falling along the second axis
  H <- diag(c(2, -2)); dn <- sqrt(abs(diag(H)))
  s <- .pc_saddle(nll, c(0, 0), c(1, 1), H / tcrossprod(dn), diag(H), c(TRUE, TRUE), 0, 1e-10, 1e-8)
  expect_identical(s$type, "saddle")
  expect_equal(abs(s$point), c(0, 1 / sqrt(2)), tolerance = 1e-12)  # the largest probe step, one normalized unit
  expect_lt(nll(s$point), 0)
})

test_that("summary() of a penalized fit marks the separated coefficients", {
  data(defection)
  m3 <- suppressMessages(mnl(outcome ~ economic * indoctrination + age_joined + birth_year + school_age + male, data = defection, ref = "captured", firth = TRUE))
  out <- capture.output(print(summary(m3)))
  ## a marked row (the footnote starts with "[sep]:"), and each penalty's own
  ## wording (round 7, C-N1 / A-3)
  expect_true(any(grepl("switched:economic \\[sep\\]", out)))
  expect_true(any(grepl("can be far off \\(Heinze and Schemper 2002\\)", out)))
  m3c <- suppressWarnings(suppressMessages(mnl(outcome ~ economic * indoctrination + age_joined + birth_year + school_age + male, data = defection,
                                               ref = "captured", prior = coef_prior())))
  outc <- capture.output(print(summary(m3c)))
  expect_true(any(grepl("switched:economic \\[sep\\]", outc)))
  expect_true(any(grepl("under a prior can have little power", outc)))
  ## the coefficient names themselves are unchanged
  expect_false(any(grepl("\\[sep\\]", rownames(summary(m3)$coefficients))))
})

test_that("simulate(ties = \"observed\") keeps each case's rank groups, as boot_lr_test() draws them", {
  set.seed(3)
  dr <- suppressMessages(rrologit(120, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6), d = c(0.2, 0.4)), depth = 3, ties = 0.3))
  r <- suppressMessages(rologit(rank ~ 1 | x1, data = dr, id = "id", alt = "alt"))
  s <- simulate(r, nsim = 2, seed = 4, ties = "observed")
  R0 <- r$design$str$rank; A <- r$design$avail == 1L
  grp <- function(R) lapply(seq_len(nrow(R)), function(i) { x <- R[i, A[i, ]]; tabulate(x[x > 0L]) })
  expect_identical(grp(unname(s$sim_1)), grp(R0))
  expect_identical(grp(unname(s$sim_2)), grp(R0))
  su <- simulate(r, nsim = 1, seed = 4)                                        # the default: untied, same number of positions
  expect_false(any(apply(su$sim_1, 1L, function(x) anyDuplicated(x[x > 0L]) > 0L)))
})

test_that("the score tests warn on separated fits; mnprobit() refuses more than 21 alternatives; texreg reads rankings", {
  data(defection)
  m <- suppressWarnings(suppressMessages(mnl(outcome ~ economic * indoctrination + age_joined + birth_year + school_age + male, data = defection, ref = "captured")))
  ## the warning comes first (here the singular information then stops the test)
  w <- character(0)
  withCallingHandlers(try(captivity_test(m), silent = TRUE), warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_true(any(grepl("diverging coefficients", w)))
  set.seed(1)
  d <- data.frame(y = factor(rep(letters[1:22], 3)), x = rnorm(66))
  expect_error(mnprobit(y ~ x, data = d), "at most 21 alternatives")
  skip_if_not_installed("texreg")
  set.seed(2)
  dr <- rrologit(100, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6)), gamma = 0.7)
  r <- rologit(rank ~ w1 | x1, data = dr, id = "id", alt = "alt")
  expect_output(print(texreg::screenreg(r)), "w1")
})

test_that("the saddle step, end to end: Liou's probit on sanctions leaves the saddle the search stops at (round 6, C-1)", {
  skip_on_cran()
  data(sanctions)
  f <- outcome ~ nonviolent * opposition + violent * opposition + growth + inflation + leg_election + exec_election + gdppc +
    trade + milex + democracy + human_rights + coldwar + peaceyears + I(peaceyears^2) + I(peaceyears^3)
  ## the multistart stops at a saddle (-1159.4773 in the stress run of build
  ## 014e91568fbd). The fit starts at that point (the fixture, written by
  ## verification/make_saddle_fixture.R), so that the saddle check and the step
  ## run whatever path the optimizer takes on a platform (round 7, C-N3); the
  ## step along the direction of ascent reaches the maximum
  fx <- readRDS(test_path("fixtures", "liou_saddle_theta.rds"))
  expect_message(mp <- suppressWarnings(mnprobit(f, data = sanctions, start = fx$theta)), "stopped at a saddle")
  expect_equal(mp$loglik, -1151.006, tolerance = 1e-3 / 1151)
  expect_false(mp$saddle)
  ## the maximum the step reached is in start_logliks (local review 2026-11-06, finding 15)
  expect_equal(max(mp$start_logliks), mp$loglik, tolerance = 1e-10)
  ## the record of the search's ends holds both points, the reported one first
  expect_equal(mp$search$loglik, sort(mp$start_logliks, decreasing = TRUE), tolerance = 1e-10)
  expect_identical(mp$search$reported, c(TRUE, FALSE))
})
