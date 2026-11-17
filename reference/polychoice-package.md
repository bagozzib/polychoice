# polychoice: Multinomial Choice, Inflated Multinomial Logit, and Tests of IIA

Estimation, inference, quantities of interest, and specification tests
for unordered polytomous (multinomial) choice outcomes, behind one
interface and one class of fitted objects (`"choicemodel"`).

## Estimators

[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) fits
the multinomial / conditional logit with case-specific covariates
(alternative-specific coefficients), alternative-specific covariates
(generic coefficients), alternative-specific constants, varying choice
sets, frequency weights, offsets, unit fixed effects (with the
split-panel jackknife `bias_correct = "jackknife"`), and analytic,
outer-product, robust, cluster-robust, or bootstrap standard errors.
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md) fits
the baseline-inflated multinomial logit of Bagozzi and Marchetti (2017)
for any inflated alternative, with a logit or probit split equation, and
with several inflated alternatives the parametrized logit captivity
model of Swait and Ben-Akiva (1987), which generalizes the dogit of
Gaudry and Dagenais (1979).
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md)
fits the two-level nested logit,
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md) the
heteroskedastic extreme value model,
[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md)
the multinomial probit with a free covariance, and
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)
the multinomial logit with random alternative-specific intercepts per
unit; these accept a wide layout (one row per case) and a long layout
(one row per case-alternative pair; see the `id` and `alt` arguments).
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
fits the rank-ordered (exploded) logit to complete, partial, or tied
rankings in the long layout.

## Separation

[`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md)
checks exactly, by a linear program, whether the data of an
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) or
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
fit are separated (no finite maximum-likelihood estimate). The
maximum-likelihood likelihood-ratio tests then compare the suprema of
the log-likelihoods. `firth = TRUE` (Firth's penalty) and
`prior = coef_prior(...)` (normal, Cauchy or t priors) give finite
penalized estimates, with profile intervals and penalized
likelihood-ratio tests.

## Quantities of interest and comparison

[`predict.choicemodel()`](https://bagozzib.github.io/polychoice/reference/predict.choicemodel.md),
[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md),
[`rrr()`](https://bagozzib.github.io/polychoice/reference/rrr.md),
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md),
[`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md),
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md),
[`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md),
[`mundlak()`](https://bagozzib.github.io/polychoice/reference/mundlak.md),
the simulators
[`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md),
[`rnlogit()`](https://bagozzib.github.io/polychoice/reference/rnlogit.md),
[`rhev()`](https://bagozzib.github.io/polychoice/reference/rhev.md),
[`rmnp()`](https://bagozzib.github.io/polychoice/reference/rmnp.md),
[`rrologit()`](https://bagozzib.github.io/polychoice/reference/rrologit.md)
and
[`simulate.choicemodel()`](https://bagozzib.github.io/polychoice/reference/simulate.choicemodel.md),
and broom / texreg / modelsummary support.

## Tests of the independence of irrelevant alternatives and of specification

[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md)
(Fok and Paap 2025), the battery
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md),
the nested-logit score test
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md)
(McFadden 1987),
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md)
(Cramer and Ridder 1991),
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md),
the goodness-of-fit test
[`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md)
(Fagerland, Hosmer and Bofin 2008), the stability test of rank-ordered
fits
[`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md),
the score and likelihood-ratio tests for captive choosers
[`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md),
and the parametric-bootstrap calibration
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md).

## Data

[pa2006](https://bagozzib.github.io/polychoice/reference/pa2006.md) and
[pres2004](https://bagozzib.github.io/polychoice/reference/pres2004.md),
the two applications of Bagozzi and Marchetti (2017), whose published
fits are regression tests of the package, and seven data sets from
published replication archives:
[sanctions](https://bagozzib.github.io/polychoice/reference/sanctions.md),
[fluvax](https://bagozzib.github.io/polychoice/reference/fluvax.md),
[bes1987](https://bagozzib.github.io/polychoice/reference/bes1987.md),
[defection](https://bagozzib.github.io/polychoice/reference/defection.md),
[bes2010](https://bagozzib.github.io/polychoice/reference/bes2010.md),
[cabinets](https://bagozzib.github.io/polychoice/reference/cabinets.md),
and
[afghan_leaders](https://bagozzib.github.io/polychoice/reference/afghan_leaders.md)
(see the vignette "Seven published analyses").

## See also

Useful links:

- <https://github.com/bagozzib/polychoice>

- Report bugs at <https://github.com/bagozzib/polychoice/issues>

## Author

Benjamin E. Bagozzi <bagozzib@udel.edu>
