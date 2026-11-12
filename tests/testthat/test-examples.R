## The example datasets from published replication archives reproduce the
## published specifications (data-raw/make_examples.R builds them).

test_that("fluvax reproduces Liu et al. (2024), Table 3, Model 1", {
  data(fluvax)
  m <- mnl(vaccinate ~ condition + age + female + educ + income + race + party + knowledge, data = fluvax, se = "robust")
  expect_equal(m$n, 1334L)
  nm <- c("yes:conditionus_china", "yes:conditionus_japan", "yes:partyrepublican", "dont_know:partyrepublican", "yes:knowledge")
  ## printed to three decimals: within half a unit of the last digit
  expect_true(all(abs(coef(m)[nm] - c(-0.385, 0.347, -1.321, -0.678, 0.527)) <= 5e-4 + 1e-6))
  expect_true(all(abs(sqrt(diag(vcov(m)))[nm] - c(0.171, 0.174, 0.167, 0.179, 0.243)) <= 5e-4 + 1e-6))
})

test_that("sanctions reproduces Liou (2024), Table 1", {
  data(sanctions)
  ## the cubic in raw years has one-SD effects above 15 but |z| > 4: not a
  ## separation (stress finding F14)
  expect_no_warning(m <- mnl(outcome ~ nonviolent * opposition + violent * opposition + growth + inflation + leg_election + exec_election + gdppc +
                               trade + milex + democracy + human_rights + coldwar + peaceyears + I(peaceyears^2) + I(peaceyears^3),
                             data = sanctions, cluster = sanctions$ccode, se = "cluster"))
  expect_false(m$boundary)
  expect_true(m$converged)
  expect_equal(m$loglik, -1161.5576, tolerance = 1e-4 / 1161)
  expect_equal(unname(coef(m)["sanctions:nonviolent:opposition"]), 8.146855, tolerance = 1e-4)
  ## Stata's clustered standard errors carry the factor sqrt(G / (G - 1)), G = 137
  expect_equal(unname(sqrt(diag(vcov(m)))["sanctions:nonviolent:opposition"]) * sqrt(137 / 136), 2.936951, tolerance = 1e-3)
})

test_that("bes1987 reproduces Glasgow (2001), Appendix A, Table C", {
  data(bes1987)
  m <- mnl(chosen ~ defense + phillips + taxation + nationalization + redistribution + crime + welfare |
             south + midlands + north + wales + scotland + union + public + working_class + female + age + homeowner + income + education +
             inflation + taxes + unemployment, data = bes1987, id = "id", alt = "party", ref = "alliance")
  expect_equal(round(m$loglik, 2), -1477.56)
  ## printed to two decimals
  expect_true(all(abs(coef(m)[c("defense", "phillips", "taxation", "nationalization", "redistribution", "crime", "welfare")] -
                        c(-0.18, -0.11, -0.16, -0.18, -0.08, -0.10, -0.14)) <= 5e-3 + 1e-6))
})

test_that("cabinets: a conditional logit with unlabelled alternatives (Baeck et al. 2025, Table 1, Model 3)", {
  data(cabinets)
  m <- mnl(chosen ~ 0 + incumbent + incumbent_newpm + minority + min_winning + bargaining + strongest + cross_cutting + het_economic +
             het_social + pea_res + rej_res | 1, data = cabinets, id = "election", alt = "cabinet")
  expect_equal(m$n, 125L)
  expect_true(m$converged)
  ## survival::clogit(method = "exact") on the same rows: -81.0040, coefficients equal to 1e-14
  expect_equal(round(m$loglik, 4), -81.0040)
  ## the comparison itself (round 7, H7-7): clogit(method = "exact") is
  ## coxph() with unit times and exact ties, called here directly because
  ## clogit() evaluates its coxph() call where survival may not be attached;
  ## Surv() and strata() are found in the formula's environment
  skip_if_not_installed("survival")
  dd <- m$data; dd$.time <- 1
  fc <- Surv(.time, chosen) ~ incumbent + incumbent_newpm + minority + min_winning + bargaining + strongest + cross_cutting + het_economic +
    het_social + pea_res + rej_res + strata(election)
  environment(fc) <- list2env(list(Surv = survival::Surv, strata = survival::strata), parent = environment())
  cl <- survival::coxph(fc, data = dd, method = "exact")
  expect_equal(m$loglik, cl$loglik[2L], tolerance = 1e-8)
  expect_equal(unname(coef(m)), unname(coef(cl)[names(coef(m))]), tolerance = 1e-6)
})

test_that("bes2010 and afghan_leaders fit as conditional and rank-ordered logits", {
  data(bes2010); data(afghan_leaders)
  m <- mnl(chosen ~ feeling + incumbent + distance + deprivation | 1, data = bes2010, id = "id", alt = "party", ref = "labour",
           cluster = "constituency", se = "cluster")
  expect_equal(m$n, 674L)
  expect_equal(round(m$loglik, 4), -328.4719)
  ma <- rologit(rank ~ 1 | treatment + female + age, data = afghan_leaders, id = "id", alt = "attribute", ref = "educated")
  expect_equal(ma$n, 2479L)
  expect_true(ma$converged)
})
