# Changelog

## polychoice 0.1.0

Initial CRAN release. Features:

- Estimators
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  (multinomial / conditional logit),
  [`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md)
  (baseline-inflated multinomial logit of Bagozzi and Marchetti 2017;
  with several inflated alternatives the parametrized logit captivity
  model of Swait and Ben-Akiva 1987),
  [`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md)
  (two-level nested logit),
  [`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md)
  (heteroskedastic extreme value),
  [`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md)
  (multinomial probit),
  [`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)
  (random alternative-specific intercepts, adaptive Gauss-Hermite
  quadrature),
  [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
  (rank-ordered logit for complete, partial, and tied rankings, with the
  exact tie likelihood of Allison and Christakis 1994).
- Firth’s penalized likelihood for
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) and
  [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
  (`firth = TRUE`): finite estimates under separation, with tied
  rankings through the information of the exact tie probability; profile
  penalized-likelihood intervals
  ([`confint()`](https://rdrr.io/r/stats/confint.html)) and penalized
  likelihood-ratio tests
  ([`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md),
  [`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
  [`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
  [`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md)).
  The penalized likelihood can have more than one local maximum, so a
  penalized fit also runs from the maximum-likelihood estimate and keeps
  the higher.
- Priors for
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) and
  [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
  (`prior = coef_prior(...)`): normal, Cauchy or t priors centred at
  zero on coefficients scaled to their columns (Gelman et al. 2008),
  fitted as the posterior mode, with the same penalized machinery as
  `firth = TRUE`; fits under several priors show how much of a separated
  estimate comes from the prior.
- [`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md):
  an exact check of separation for
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) and
  [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
  data (a linear program on the comparisons the choices or rankings
  report, ties included, through lpSolveAPI), which their
  maximum-likelihood fits use to name diverging coefficients. Under
  separation the maximum-likelihood forms of
  [`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md),
  [`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
  [`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
  and
  [`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
  compare the suprema of the log-likelihoods, not where the optimizer
  stopped.
- [`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
  also calibrates the comparison of two nested
  [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
  fits, drawing rankings with each case’s observed rank groups.
- Wide and long data layouts, varying choice sets, frequency weights,
  offsets, unit fixed effects with the never-chosen rule and the
  split-panel jackknife (`bias_correct = "jackknife"`), analytic /
  outer-product / robust / cluster / bootstrap standard errors.
- Quantities of interest:
  [`predict()`](https://rdrr.io/r/stats/predict.html) (probabilities,
  classes, utilities, outcome-stage and regime probabilities, posterior
  regime probabilities),
  [`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md),
  [`rrr()`](https://bagozzib.github.io/polychoice/reference/rrr.md).
- Fit and comparison:
  [`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md),
  [`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md),
  [`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
  [`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md),
  [`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md),
  [`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md);
  broom, texreg, and modelsummary support.
- IIA and specification tests:
  [`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
  [`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
  [`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
  [`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md),
  [`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md),
  [`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
  [`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md),
  [`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md),
  [`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
  (the parametric bootstrap of the tests whose null is the fitted
  model);
  [`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md)
  (stability of preferences across ranking stages);
  [`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md)
  (score test of the multinomial logit against captive choosers, the
  dogit of Gaudry and Dagenais 1979, and with `type = "lr"` the
  likelihood ratio against a constant captive share, which
  [`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
  also takes).
- Simulators
  [`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md),
  [`rnlogit()`](https://bagozzib.github.io/polychoice/reference/rnlogit.md),
  [`rhev()`](https://bagozzib.github.io/polychoice/reference/rhev.md),
  [`rmnp()`](https://bagozzib.github.io/polychoice/reference/rmnp.md),
  [`rrologit()`](https://bagozzib.github.io/polychoice/reference/rrologit.md);
  [`mundlak()`](https://bagozzib.github.io/polychoice/reference/mundlak.md).
- Data sets `pa2006` and `pres2004` (Bagozzi and Marchetti 2017), and
  `sanctions`, `fluvax`, `bes1987`, `defection`, `bes2010`, `cabinets`,
  `afghan_leaders` from published replication archives, with the
  vignette *Seven published analyses*.
