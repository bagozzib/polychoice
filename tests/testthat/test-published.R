## Regression tests against PUBLISHED results: Bagozzi and Marchetti (2017, PSRM),
## Tables 1/2 (log-likelihoods) and A.3/A.4 (coefficients) on the bundled data.

test_that("pa2006 (Replication I): MNL and BIMNL reproduce the published fits", {
  data(pa2006)
  f <- vote ~ canvass + phone + democrat + republican + canvass:democrat + canvass:republican +
    phone:democrat + phone:republican + household + age + female + voted2004 + district156
  m0 <- mnl(f, data = pa2006)
  expect_equal(m0$n, 1998L)
  expect_equal(m0$loglik, -2110.59, tolerance = 1e-5)                       # published -2110.59
  ## Table A.3, MNL: R/O vs abstain canvass 1.222 (0.455); NR vs abstain age 0.019 (0.006); R/O constant -4.677 (0.534)
  s0 <- summary(m0)$coefficients
  ## absolute agreement with the printed values (the article's Table rep: the MNL within 0.004, round 7, H7-2)
  expect_lte(max(abs(unname(s0["republican_other:canvass", 1:2]) - c(1.222, 0.455))), 0.004)
  expect_lte(max(abs(unname(s0["nonresponse:age", 1:2]) - c(0.019, 0.006))), 0.004)
  expect_lte(max(abs(unname(s0["republican_other:(Intercept)", 1:2]) - c(-4.677, 0.534))), 0.004)
  skip_on_cran()
  m1 <- imnl(vote ~ canvass + phone + democrat + republican + canvass:democrat + canvass:republican +
               phone:democrat + phone:republican + household + age + female + voted2004 + district156 |
               household + age + female + voted2004 + district156, data = pa2006)
  expect_true(m1$converged); expect_false(m1$boundary)
  expect_equal(m1$loglik, -2069.97, tolerance = 1e-5)                       # published -2069.97
  expect_gte(m1$loglik, -2069.9711)                                          # at least the published optimum (to its precision)
  expect_equal(m1$n_modes, 2L)                                               # the likelihood has a second, lower mode (-2070.99)
  s1 <- summary(m1)$coefficients
  ## Table A.3, BIMNL: R/O canvass 1.566 (0.558); inflation voted2004 2.410 (0.389), district156 0.641 (0.216), age -0.036 (0.016)
  ## the BIMNL within 0.002 of the printed values (Table rep)
  expect_lte(max(abs(unname(s1["republican_other:canvass", 1:2]) - c(1.566, 0.558))), 0.002)
  expect_lte(max(abs(unname(s1["infl_voted2004", 1:2]) - c(2.410, 0.389))), 0.002)
  expect_lte(max(abs(unname(s1["infl_district156", 1:2]) - c(0.641, 0.216))), 0.002)
  expect_lte(max(abs(unname(s1["infl_age", 1:2]) - c(-0.036, 0.016))), 0.002)
  ## all starts reach the same optimum (within the optimizer tolerance) except possibly the lowest shift
  expect_equal(max(m1$start_logliks), m1$loglik, tolerance = 1e-8)
  cm <- compare_models(MNL = m0, BIMNL = m1)
  expect_lt(cm$AIC[2], cm$AIC[1]); expect_lt(cm$BIC[2], cm$BIC[1])          # Table 1: all fit statistics favor BIMNL
  v <- vuong(m1, m0)
  expect_equal(v$statistic[v$correction == "raw"], 4.48, tolerance = 1e-3)   # Table 1: Vuong favours BIMNL
  expect_gt(classification(m1)$pre, classification(m0)$pre)                   # PRE 0.316 vs 0.299
})

test_that("pres2004 (Replication II): MNL and BIMNL reproduce the published fits", {
  data(pres2004)
  f <- vote ~ educ + female + black + hispanic + pid7 + bush_values + battleground + south + gmb +
    agegroup + evangelical + evangelical:gmb + catholic + catholic:gmb + other_religion + secular + secular:gmb + ideology + mobilization
  m0 <- mnl(f, data = pres2004)
  expect_equal(m0$n, 1341L)
  expect_equal(m0$loglik, -553.50, tolerance = 1e-5)                         # published -553.50
  s0 <- summary(m0)$coefficients
  expect_lte(max(abs(unname(s0["kerry:pid7", 1:2]) - c(0.344, 0.075))), 0.001)      # Table A.4 MNL
  expect_lte(max(abs(unname(s0["bush:pid7", 1:2]) - c(-0.386, 0.076))), 0.001)
  skip_on_cran()
  m1 <- imnl(vote ~ educ + female + black + hispanic + pid7 + bush_values + battleground + south + gmb +
               agegroup + evangelical + evangelical:gmb + catholic + catholic:gmb + other_religion + secular + secular:gmb + ideology + mobilization |
               educ + agegroup + mobilization + secular, data = pres2004)
  expect_true(m1$converged); expect_false(m1$boundary)
  expect_equal(m1$loglik, -545.27, tolerance = 2e-5)                         # published -545.27
  s1 <- summary(m1)$coefficients
  ## Table A.4, BIMNL inflation stage: constant -4.938 (1.985), educ 1.507 (0.528), agegroup 0.743 (0.437), mobilization 1.139 (0.313), secular -1.295 (0.789)
  ## the inflation-stage coefficients within 0.001 of the printed values and
  ## their standard errors within 0.0015 (Table rep; round 7, H7-2)
  pub <- rbind(`infl_(Intercept)` = c(-4.938, 1.985), infl_educ = c(1.507, 0.528), infl_agegroup = c(0.743, 0.437),
               infl_mobilization = c(1.139, 0.313), infl_secular = c(-1.295, 0.789))
  got <- s1[rownames(pub), 1:2]
  expect_lte(max(abs(got[, 1] - pub[, 1])), 0.001); expect_lte(max(abs(got[, 2] - pub[, 2])), 0.0015)
  expect_lte(max(abs(unname(s1["kerry:ideology", 1:2]) - c(-0.349, 0.233))), 0.001)
  expect_lte(max(abs(unname(s1["bush:pid7", 1:2]) - c(-0.390, 0.092))), 0.001)
  ## Table 2: mixed support -- AIC favors BIMNL, BIC favors MNL
  cm <- compare_models(MNL = m0, BIMNL = m1)
  expect_lt(cm$AIC[2], cm$AIC[1]); expect_gt(cm$BIC[2], cm$BIC[1])
  ## Table 2's Vuong and PRE recomputed from the fitted probabilities
  ## (verification/table2_check.R): the Vuong statistic is 1.680 towards BIMNL
  ## (two-sided p = 0.093) and the two models
  ## classify the same share of cases correctly (PRE 0.7045 for both, a tie;
  ## the case-level predictions differ in a few cases); LR 16.44 on 5 df
  v <- vuong(m1, m0)
  expect_equal(v$statistic[v$correction == "raw"], 1.680, tolerance = 1e-3)
  c0 <- classification(m0); c1 <- classification(m1)
  expect_equal(c1$pre, c0$pre); expect_equal(c0$pre, 0.7045, tolerance = 1e-3)
  expect_equal(2 * (m1$loglik - m0$loglik), 16.44, tolerance = 1e-3)
})
