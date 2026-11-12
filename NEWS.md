# polychoice 0.1.0

Initial CRAN release. Features:

* Estimators `mnl()` (multinomial / conditional logit), `imnl()`
  (baseline-inflated multinomial logit of Bagozzi and Marchetti 2017; with
  several inflated alternatives the parametrized logit captivity model of
  Swait and Ben-Akiva 1987),
  `nlogit()` (two-level nested logit), `hev()` (heteroskedastic extreme
  value), `mnprobit()` (multinomial probit), `remnl()` (random
  alternative-specific intercepts, adaptive Gauss-Hermite quadrature),
  `rologit()` (rank-ordered logit for complete, partial, and tied rankings,
  with the exact tie likelihood of Allison and Christakis 1994).
* Firth's penalized likelihood for `mnl()` and `rologit()` (`firth = TRUE`):
  finite estimates under separation, with tied rankings through the
  information of the exact tie probability; profile penalized-likelihood
  intervals (`confint()`) and penalized likelihood-ratio tests
  (`term_test()`, `combine_test()`, `lr_test()`, `rank_test()`). The
  penalized likelihood can have more than one local maximum, so a penalized
  fit also runs from the maximum-likelihood estimate and keeps the higher.
* Priors for `mnl()` and `rologit()` (`prior = coef_prior(...)`): normal,
  Cauchy or t priors centred at zero on coefficients scaled to their columns
  (Gelman et al. 2008), fitted as the posterior mode, with the same penalized
  machinery as `firth = TRUE`; fits under several priors show how much of a
  separated estimate comes from the prior.
* `separation()`: an exact check of separation for `mnl()` and `rologit()`
  data (a linear program on the comparisons the choices or rankings report,
  ties included, through lpSolveAPI), which their maximum-likelihood fits use
  to name diverging coefficients. Under separation the maximum-likelihood
  forms of `term_test()`, `combine_test()`, `lr_test()` and `boot_lr_test()`
  compare the suprema of the log-likelihoods, not where the optimizer
  stopped.
* `boot_lr_test()` also calibrates the comparison of two nested `rologit()`
  fits, drawing rankings with each case's observed rank groups.
* Wide and long data layouts, varying choice sets, frequency weights, offsets,
  unit fixed effects with the never-chosen rule and the split-panel jackknife
  (`bias_correct = "jackknife"`), analytic / outer-product / robust / cluster /
  bootstrap standard errors.
* Quantities of interest: `predict()` (probabilities, classes, utilities,
  outcome-stage and regime probabilities, posterior regime probabilities),
  `first_difference()`, `rrr()`.
* Fit and comparison: `fit_stats()`, `classification()`, `lr_test()`,
  `vuong()`, `boot_lr_test()`, `compare_models()`; broom, texreg, and
  modelsummary support.
* IIA and specification tests: `hausman_mcfadden()`, `small_hsiao()`,
  `fp_tests()`, `iia_tests()`, `nested_lm_test()`, `combine_test()`,
  `term_test()`, `hl_test()`, `boot_pvalue()` (the parametric bootstrap of
  the tests whose null is the fitted model); `rank_test()` (stability of
  preferences across ranking stages); `captivity_test()` (score test of the
  logit against captive choosers, the dogit of Gaudry and Dagenais 1979, and
  with `type = "lr"` the likelihood ratio against a constant captive share,
  which `lr_test()` also takes).
* Simulators `rmnl()`, `rnlogit()`, `rhev()`, `rmnp()`, `rrologit()`; `mundlak()`.
* Data sets `pa2006` and `pres2004` (Bagozzi and Marchetti 2017), and
  `sanctions`, `fluvax`, `bes1987`, `defection`, `bes2010`, `cabinets`,
  `afghan_leaders` from published replication archives, with the vignette
  *Seven published analyses*.
