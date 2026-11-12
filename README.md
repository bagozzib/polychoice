# polychoice

Multinomial choice models, the inflated multinomial logit, and tests of the
independence of irrelevant alternatives (IIA), for R.

**Status: development version, not yet released.** See `DESIGN.md` in the repository for the
scope, the parameterisation, the fixed definitions of every test, and the
verification protocol.

## What it does

* `mnl()`: multinomial / conditional logit with case-specific covariates
  (alternative-specific coefficients), alternative-specific covariates (generic
  coefficients), alternative-specific constants, varying choice sets,
  frequency weights, offsets, unit fixed effects, and analytic, outer-product (BHHH), robust,
  cluster-robust, or bootstrap standard errors. Wide (one row per case) and
  long (one row per case-alternative pair) layouts.
* `imnl()`: the baseline-inflated multinomial logit of Bagozzi and Marchetti
  (2017, PSRM), for any inflated alternative, with a logit or probit split
  equation; with several inflated alternatives, the parametrized logit
  captivity model of Swait and Ben-Akiva (1987), which generalizes the dogit
  of Gaudry and Dagenais (1979); and
  `captivity_test()`, score and likelihood-ratio tests of the multinomial logit against
  captive choosers.
* `nlogit()`: two-level nested logit (RU2 / RUM-consistent form, as mlogit
  and Stata), free or shared inclusive-value parameters.
* `hev()`: heteroskedastic extreme value model (alternative-specific scales),
  with quadrature checked against integrate() for scale ratios up to 1.6e5 and a boundary flag when a
  scale runs to zero.
* `mnprobit()`: multinomial probit with a free covariance of the utility
  differences (exact bivariate normal for three alternatives, GHK with Halton
  draws beyond).
* `remnl()`: multinomial logit with random alternative-specific intercepts per
  unit (adaptive Gauss-Hermite quadrature; at two alternatives it reproduces
  `lme4::glmer`).
* `rologit()`: rank-ordered (exploded) logit for complete, partial, and tied
  rankings (exact tie likelihood; reproduces mlogit's rank-ordered logit), with
  `rank_test()` for the stability of preferences across ranking stages.
* Separation:
  - `separation()` is an exact check by linear programming (ties as groups).
  - The maximum-likelihood likelihood-ratio tests compare the suprema of the
    log-likelihoods under separation. `boot_lr_test()` calibrates them,
    rank-ordered fits included.
  - `firth = TRUE` (Firth's penalty) and `prior = coef_prior(...)` (normal,
    Cauchy or t priors, Gelman et al.'s default scaling) give finite
    estimates in `mnl()` and `rologit()`, tied rankings included. They come
    with profile penalized-likelihood intervals (`confint()`) and penalized
    likelihood-ratio tests (`term_test()`, `combine_test()`, `lr_test()`,
    `rank_test()`).
* Fixed-effects panels: unit dummies with the never-chosen rule, and
  `bias_correct = "jackknife"` (split-panel jackknife of the slopes with a
  time-homogeneity gate).
* Quantities of interest: `predict()` (observed-alternative, outcome-stage, and
  regime probabilities), `first_difference()` (delta-method or simulation
  intervals), `rrr()`.
* Fit and comparison: `fit_stats()`, `classification()`, `lr_test()`,
  `vuong()`, `boot_lr_test()` (parametric-bootstrap LR and Vuong tests for
  restrictions on the boundary, such as the multinomial logit against the inflated multinomial logit),
  `compare_models()`; broom, texreg, and modelsummary support.
* IIA / misspecification tests: `hausman_mcfadden()` (variance-difference and
  joint/suest forms, base drop handled by exact re-basing), `small_hsiao()`
  (the correct sqrt(1/2)-weighted form, repeatable over splits), `fp_tests()`
  (Fok and Paap 2025: composite-likelihood Hausman and GMM overidentification),
  and the battery `iia_tests()`.
* Further specification tests: `nested_lm_test()` (McFadden 1987 score test
  of IIA against a nested logit), `combine_test()` (Cramer-Ridder pooling of
  alternatives), `term_test()` (joint Wald / LR tests per term), `hl_test()`
  (Fagerland-Hosmer-Bofin goodness of fit), and `boot_pvalue()` (parametric
  bootstrap calibration of the tests whose null is the fitted model).
* Simulators: `rmnl()`, `rnlogit()`, `rhev()`, `rmnp()`, `rrologit()`.
* Panel tools: `mundlak()`.
* Data: `pa2006` and `pres2004`, the two applications of Bagozzi and
  Marchetti (2017); and seven datasets from published replication archives,
  `sanctions` (Liou 2024), `fluvax` (Liu et al. 2024), `bes1987` (Glasgow
  2001), `defection` (Oppenheim et al. 2015, as re-analysed by Cook, Niehaus
  and Zuhlke 2018), `bes2010` (Arzheimer and Evans 2012), `cabinets` (Baeck,
  Debus and Imre 2026) and `afghan_leaders` (Bhatia and Monroe 2024). The
  published fits of the first four are regression tests.

```r
library(polychoice)
data(pres2004)
f <- vote ~ educ + female + black + hispanic + pid7 + bush_values + battleground + south + gmb +
  agegroup + evangelical + evangelical:gmb + catholic + catholic:gmb + other_religion +
  secular + secular:gmb + ideology + mobilization
m0 <- mnl(f, data = pres2004)
m1 <- imnl(vote ~ educ + female + black + hispanic + pid7 + bush_values + battleground + south + gmb +
             agegroup + evangelical + evangelical:gmb + catholic + catholic:gmb + other_religion +
             secular + secular:gmb + ideology + mobilization | educ + agegroup + mobilization + secular,
           data = pres2004)
compare_models(MNL = m0, BIMNL = m1)
vuong(m1, m0)
boot_lr_test(m0, m1, B = 99, cores = 4)   # calibrated (the multinomial logit is the boundary case of the inflated model); slow:
                                          # every replicate refits the inflated model's multistart
iia_tests(m0, seed = 1)
```

## Installation (development)

```r
# from the package directory
install.packages(".", repos = NULL, type = "source")
```

Requires a C++ compiler (Rtools on Windows).

## Citation and licence

`citation("polychoice")` gives the reference for the package and for the
methods it implements. The package is free software under the GNU General
Public License, version 3 (GPL-3).
