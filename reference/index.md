# Package index

## Package overview

- [`polychoice`](https://bagozzib.github.io/polychoice/reference/polychoice-package.md)
  [`polychoice-package`](https://bagozzib.github.io/polychoice/reference/polychoice-package.md)
  : polychoice: Multinomial Choice, Inflated Multinomial Logit, and
  Tests of IIA

## Estimators

One formula interface, one class, one set of methods.

- [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) :
  Multinomial and conditional logit regression
- [`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md) :
  Inflated multinomial logit regression
- [`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md)
  : Nested logit regression
- [`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md) :
  Heteroskedastic extreme value regression
- [`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md)
  : Multinomial probit regression
- [`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)
  : Multinomial logit with random alternative-specific intercepts
- [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
  : Rank-ordered (exploded) logit regression

## Separation and priors

- [`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md)
  : Separation in multinomial, conditional and rank-ordered logit data
- [`coef_prior()`](https://bagozzib.github.io/polychoice/reference/coef_prior.md)
  : Prior distributions for the coefficients of a penalized fit

## Inference and summaries

- [`summary(`*`<choicemodel>`*`)`](https://bagozzib.github.io/polychoice/reference/summary.choicemodel.md)
  : Summarize a choice model
- [`confint(`*`<choicemodel>`*`)`](https://bagozzib.github.io/polychoice/reference/confint.choicemodel.md)
  : Confidence intervals for a choice model
- [`update(`*`<choicemodel>`*`)`](https://bagozzib.github.io/polychoice/reference/update.choicemodel.md)
  : Update and refit a choice model
- [`residuals(`*`<choicemodel>`*`)`](https://bagozzib.github.io/polychoice/reference/residuals.choicemodel.md)
  : Residuals of a choice model

## Quantities of interest

- [`predict(`*`<choicemodel>`*`)`](https://bagozzib.github.io/polychoice/reference/predict.choicemodel.md)
  : Predictions from a choice model
- [`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md)
  : First differences in predicted probabilities
- [`rrr()`](https://bagozzib.github.io/polychoice/reference/rrr.md) :
  Relative risk ratios

## Fit and model comparison

- [`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md)
  : Fit statistics for a choice model
- [`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md)
  : Classification table of a choice model
- [`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md)
  : Compare several choice models side by side
- [`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
  : Likelihood-ratio test for nested choice models
- [`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)
  : Vuong test for non-nested or overlapping choice models
- [`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
  : Parametric-bootstrap likelihood-ratio (and Vuong) test of nested
  fits

## Tests of IIA

- [`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)
  : The battery of IIA tests
- [`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md)
  : Hausman-McFadden test of IIA
- [`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md)
  : Small-Hsiao test of IIA
- [`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md)
  : Fok-Paap tests of multinomial logit misspecification
- [`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md)
  : Nested-logit score test of IIA

## Further specification tests

- [`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md)
  : Tests for captive (inflated) choosers: score and likelihood ratio
- [`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md)
  : Test whether two alternatives can be combined
- [`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)
  : Joint tests of each term across equations
- [`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md)
  : Hosmer-Lemeshow-type goodness-of-fit test for multinomial models
- [`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md)
  : Stability of preferences across ranking stages
- [`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
  : Parametric-bootstrap p-values for a specification test

## Panel tools

- [`mundlak()`](https://bagozzib.github.io/polychoice/reference/mundlak.md)
  : Mundlak (correlated random effects) device

## Simulation

- [`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md) :
  Simulate multinomial-logit (optionally inflated) choice data
- [`rnlogit()`](https://bagozzib.github.io/polychoice/reference/rnlogit.md)
  : Simulate nested-logit choice data
- [`rhev()`](https://bagozzib.github.io/polychoice/reference/rhev.md) :
  Simulate heteroskedastic extreme value choice data
- [`rmnp()`](https://bagozzib.github.io/polychoice/reference/rmnp.md) :
  Simulate multinomial probit choice data
- [`rrologit()`](https://bagozzib.github.io/polychoice/reference/rrologit.md)
  : Simulate rank-ordered choice data
- [`simulate(`*`<choicemodel>`*`)`](https://bagozzib.github.io/polychoice/reference/simulate.choicemodel.md)
  : Simulate outcomes from a fitted choice model

## broom methods

- [`tidy(`*`<choicemodel>`*`)`](https://bagozzib.github.io/polychoice/reference/tidy.choicemodel.md)
  : Tidy a choice model (broom method)
- [`glance(`*`<choicemodel>`*`)`](https://bagozzib.github.io/polychoice/reference/glance.choicemodel.md)
  : Glance at a choice model (broom method)
- [`augment(`*`<choicemodel>`*`)`](https://bagozzib.github.io/polychoice/reference/augment.choicemodel.md)
  : Augment the estimation cases with fitted classes and probabilities
  (broom method)

## Data

- [`pa2006`](https://bagozzib.github.io/polychoice/reference/pa2006.md)
  : 2006 Pennsylvania state House vote choice (Arceneaux and Kolodny
  2009)
- [`pres2004`](https://bagozzib.github.io/polychoice/reference/pres2004.md)
  : 2004 presidential vote choice (Campbell and Monson 2008)
- [`sanctions`](https://bagozzib.github.io/polychoice/reference/sanctions.md)
  : Economic sanctions, military force and domestic unrest (Liou 2024)
- [`fluvax`](https://bagozzib.github.io/polychoice/reference/fluvax.md)
  : Intended flu vaccination and international collaboration (Liu et al.
  2024)
- [`bes1987`](https://bagozzib.github.io/polychoice/reference/bes1987.md)
  : Vote choice and issue distances in the 1987 British general election
  (Glasgow 2001)
- [`defection`](https://bagozzib.github.io/polychoice/reference/defection.md)
  : Who leaves insurgent groups: Colombian ex-combatants (Oppenheim et
  al. 2015; Cook, Niehaus and Zuhlke 2018)
- [`bes2010`](https://bagozzib.github.io/polychoice/reference/bes2010.md)
  : Candidate-voter distance and vote choice in England, 2010 (Arzheimer
  and Evans 2012)
- [`cabinets`](https://bagozzib.github.io/polychoice/reference/cabinets.md)
  : Government formation in the German states, 1990-2023 (Baeck, Debus
  and Imre 2026)
- [`afghan_leaders`](https://bagozzib.github.io/polychoice/reference/afghan_leaders.md)
  : Rankings of leader attributes in northern Afghanistan (Bhatia and
  Monroe 2024)
