# Multinomial and conditional logit regression

Fits the multinomial logit (case-specific covariates,
alternative-specific coefficients), the conditional logit
(alternative-specific covariates, generic coefficients), or their
combination by maximum likelihood, with alternative-specific constants,
varying choice sets, frequency weights, offsets, unit fixed effects, and
analytic, robust, cluster-robust, or bootstrap standard errors. This is
the baseline for
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md) and
the model the tests of IIA
([`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md))
take.

## Usage

``` r
mnl(
  formula,
  data,
  id = NULL,
  alt = NULL,
  ref = NULL,
  weights = NULL,
  offset = NULL,
  fe = NULL,
  bias_correct = c("none", "jackknife"),
  se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"),
  cluster = NULL,
  nboot = 200L,
  cores = 1L,
  start = NULL,
  maxit = 1000L,
  reltol = 1e-10,
  firth = FALSE,
  prior = NULL
)
```

## Arguments

- formula:

  A formula; see the section "Two data layouts" of `mnl()`.

- data:

  A data frame.

- id, alt:

  For the long layout, the names of the case and alternative columns of
  `data`. Leave both `NULL` for the wide layout.

- ref:

  The base (reference) alternative: a level of the response. Default:
  the first level.

- weights:

  Optional frequency weights: a column name in `data` or a numeric
  vector with one entry per row (constant within a case in the long
  layout).

- offset:

  Optional offset on the utility scale: in the long layout a column name
  or a vector with one entry per row; in the wide layout an \\n \times
  J\\ matrix (columns in level order).

- fe:

  Optional column name in `data` identifying units that receive
  alternative-specific fixed effects.

- bias_correct:

  With `fe`: `"none"` (default) or `"jackknife"`, the split-panel
  jackknife of Dhaene and Jochmans (2015) for the slopes; see the
  fixed-effects section.

- se:

  Standard errors: `"analytic"` (inverse observed information, default),
  `"opg"` (inverse outer product of the scores, the BHHH form that
  `mlogit` reports for the nested logit and HEV), `"robust"`,
  `"cluster"` (needs `cluster`), `"bootstrap"` (nonparametric over
  cases, or over clusters when `cluster` is given, with `nboot` refits
  from the full-sample estimate; the replicates are kept in
  `object$boot` and [`confint()`](https://rdrr.io/r/stats/confint.html)
  can report percentile intervals), or `"none"`. With `fe`, the
  bootstrap resamples units (a unit drawn twice counts as two units,
  each with its own fixed effects; a case bootstrap would give units
  random numbers of cases and a second incidental-parameters bias) and
  covers the slopes only: the alternative-specific constants and the
  unit dummies refer to reference units that change between replicates,
  so their bootstrap standard errors are `NA`, and
  [`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md)
  refuses intervals that need them (use analytic, robust, or cluster
  standard errors for those). A replicate whose fixed effects cannot be
  built (an alternative left unconnected to the base) counts as failed,
  the standard errors are then conditional on resamples in which they
  can be built, and the fit stops when fewer than half the replicates
  are usable. A `cluster` given with `fe` must nest the units (whole
  units are drawn; one that cuts across units is refused). The
  replicates are not jackknife-corrected.

- cluster:

  Optional cluster identifier (a column name in `data` or a vector
  aligned to its rows, constant within a case); supplying it selects
  `se = "cluster"` unless `se` is given explicitly.

- nboot:

  Number of bootstrap refits for `se = "bootstrap"` (default 200).

- cores:

  Number of parallel workers for the bootstrap (a PSOCK cluster; default
  1). [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  makes the resamples reproducible; serial and parallel runs use
  different random streams.

- start:

  Optional starting values: the full parameter vector in the order of
  [`coef()`](https://rdrr.io/r/stats/coef.html) but on the INTERNAL
  scale, which differs from
  [`coef()`](https://rdrr.io/r/stats/coef.html) for the structural
  parameters (the logarithms of the inclusive values, scales, and
  standard deviations; the free Cholesky elements with a log diagonal
  for the multinomial probit covariance). `object$theta` of an earlier
  fit is on that scale. Rarely needed for `mnl()`, whose log-likelihood
  is globally concave; the other estimators start from the `mnl()`
  estimates (and
  [`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md)
  from a multistart grid), so `start` serves to probe another basin of a
  non-concave likelihood.

- maxit, reltol:

  Controls for the BFGS stage passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) (the plain
  `mnl()` is solved by damped Newton steps first, with BFGS as the
  fallback); the Newton polish that follows runs up to 25 further
  iterations with the analytic Hessian of the logit family or central
  differences of the analytic gradient otherwise.

- firth:

  If `TRUE`, Firth's penalized likelihood (see the section "Penalized
  (Firth) estimation"); default `FALSE`, maximum likelihood. Weights are
  read as frequency weights and must total more than 2.

- prior:

  Optional: priors for the coefficients, made by
  [`coef_prior()`](https://bagozzib.github.io/polychoice/reference/coef_prior.md)
  (normal, Cauchy or t, centred at zero, on coefficients scaled to their
  columns); the fit is the posterior mode, a penalized fit like
  `firth = TRUE` (profile intervals, penalized likelihood-ratio tests,
  the separated coefficients in `separation`). Not with `firth = TRUE`
  or unit fixed effects.

## Value

An object of class `c("mnl", "choicemodel")`: a list with
`coefficients`, `vcov`, `loglik`, `fitted.values` (the \\n \times J\\
matrix of choice probabilities, rows named by case), `converged`,
`boundary`, `ill_conditioned`, `singular`, `n_modes`, `search`,
`start_logliks`, the estimation `data`, and the design used by the
methods. See
[`summary.choicemodel()`](https://bagozzib.github.io/polychoice/reference/summary.choicemodel.md),
[`predict.choicemodel()`](https://bagozzib.github.io/polychoice/reference/predict.choicemodel.md),
[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md),
[`rrr()`](https://bagozzib.github.io/polychoice/reference/rrr.md),
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md).

## Details

With case \\i\\, alternatives \\j = 1, \dots, J\\ (the levels of the
response), base alternative \\b\\, and availability \\a\_{ij}\\,
\$\$V\_{ij} = o\_{ij} + x_i'\beta_j + w\_{ij}'\gamma, \quad \beta_b = 0,
\qquad P(y_i = j) = \frac{a\_{ij} e^{V\_{ij}}}{\sum_k a\_{ik}
e^{V\_{ik}}}.\$\$ \\x_i\\ are the case-specific covariates, \\w\_{ij}\\
the alternative-specific covariates, \\o\_{ij}\\ an offset: McFadden's
(1974) conditional logit, with the multinomial logit as the case \\w =
0\\ (Train 2009, chapter 3). The alternative-specific constants are the
coefficients of an intercept in \\x_i\\: the formula's intercept in the
wide layout, the first part's in the long layout (section "Two data
layouts"). The coefficients are reported as `"<alternative>:<term>"` for
\\\beta_j\\ and `"<term>"` for \\\gamma\\.

## Two data layouts

**Wide** (one row per case): `data` is a plain data frame, the response
is the chosen alternative (a factor, or integer codes; a character
vector is refused because its alphabetical order is rarely the intended
one), and the formula is `y ~ x1 + x2`. All alternatives are available
to every case.

**Long** (one row per case-alternative pair): give `id` (the case
column) and `alt` (the alternative column; its levels define the
alternatives). The response is the 0/1 (or logical) choice indicator,
exactly one per case. The formula has up to two parts separated by `|`:
`choice ~ w1 + w2 | x1 + x2`, alternative-specific covariates first
(generic coefficients; the intercept of this part switches the
alternative-specific constants on, `0 +` switches them off),
case-specific covariates second (they must be constant within a case).
Alternatives without a row for a case are unavailable to it. An
alternative-specific coefficient on an alternative-specific variable is
the interaction with the alternative factor, `w:alt`.

The constants follow the first part only:

- `choice ~ w`, `choice ~ w | 1` and `choice ~ w | x` have them;

- `choice ~ 0 + w`, `choice ~ w - 1` and `choice ~ 0 + w | x` do not.

The second part's intercept has no effect. `| 1` is an empty
case-specific part (constants only, when the first part has them), and
`| 0` does not drop the constants, unlike in mlogit (a message says so).
The wide-layout formula `y ~ x1 + x2` is the second part alone, its
intercept giving the constants (`y ~ 0 + x1` drops them).

## Missing values, weights, standard errors

Rows with a missing value in any model variable (or in `weights`,
`offset`, `cluster`, `fe`) are dropped before fitting; in the long
layout a missing value on any row of a case drops the whole case (a
message gives the count), and cases left with a single available
alternative are dropped too, as are cases of weight zero. `object$data`
holds the retained rows. `weights` are frequency weights throughout: a
fit with weight 2 equals a fit on duplicated cases in the estimates, the
model-based standard errors, the robust and cluster sandwiches (meat
\\\sum_i w_i s_i s_i'\\), the information criteria (the weight total is
\\n\\), the Vuong and likelihood-ratio statistics, the IIA tests (except
the Small-Hsiao test, whose half-samples balance cases rather than
weight), the averaged quantities of interest, and the bootstraps
(`se = "bootstrap"` draws multinomial replicate weights with the weight
total as size, over units for
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md);
with a `cluster` the clusters are drawn and carry their rows' weights;
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
and
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
expand the data to one row per case and need integer weights).
Non-integer weights are accepted as frequency weights (a message says
so). Survey (probability) weights are NOT supported: the point estimates
would be the pseudo-maximum-likelihood ones, but the sandwich needs
\\\sum_i w_i^2 s_i s_i'\\ and the information criteria, tests, and
bootstraps have no design-based meaning (the cluster sandwich with the
primary sampling units as clusters, \\\sum_g (\sum\_{i \in g} w_i
s_i)(\sum\_{i \in g} w_i s_i)'\\, has the with-replacement linearization
form, but strata and finite-population corrections are not handled).
`se = "robust"` is the sandwich over cases, `se = "cluster"` over
clusters (no small-sample degrees-of-freedom adjustment, as in the `iop`
and `underdisp` packages; Stata's `vce(cluster)` for maximum-likelihood
commands multiplies by \\G/(G-1)\\).
[`nobs()`](https://rdrr.io/r/stats/nobs.html) returns the number of
cases; [`logLik()`](https://rdrr.io/r/stats/logLik.html) carries the
weight total as its `nobs` attribute, which is what
[`BIC()`](https://rdrr.io/r/stats/AIC.html) uses.

## Unit fixed effects

`fe = "<unit column>"` (the name of the column identifying the units,
such as `fe = "country"`) adds an alternative-specific constant per
unit: one dummy per unit in every non-base equation. A unit that never
chooses alternative \\j\\ has that constant at \\-\infty\\; the limiting
likelihood treats \\j\\ as unavailable to that unit, which is what the
fit does (with a message), and cases left with a single alternative drop
out (a captive choice says nothing about preferences and does not count
as a choice; as dropping cases can change what a unit chose, the rule
and the drop are repeated until nothing changes). A unit that chose a
single alternative therefore drops out entirely. A unit that never
chooses the base likewise has the base unavailable; its constants are
then identified only up to a common shift, so one of its dummies is
omitted. The common constant of alternative \\j\\ refers to a reference
unit: the first unit (in level order) that chose both \\j\\ and the
base, or else one that links \\j\\ to an alternative already anchored;
an alternative whose choosers are not connected to the base in this way
has no identified constant, and the fit stops and names it. Without
alternative-specific constants (a formula with `- 1`) there is nothing
to anchor, and every unit keeps all its dummies (one fewer for a unit
that never chose the base). The dummies are dense columns, so a fit with
thousands of units and many cases needs memory for an \\nJ \times\\
(units \\\times\\ alternatives) matrix. Maximum-likelihood unit
constants carry the incidental-parameters bias when units have few cases
(the fit warns below a median of 10). For short panels the alternatives
are the jackknife correction below and the correlated-random-effects
model:
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md) on
the formula and data that
[`mundlak()`](https://bagozzib.github.io/polychoice/reference/mundlak.md)
returns (the pooled `mnl()` fit of the same formula gives
population-averaged slopes, attenuated towards zero).

`bias_correct = "jackknife"` removes the leading \\1/T\\ term of that
bias from the slopes by the split-panel jackknife (Dhaene and Jochmans
2015, Review of Economic Studies): the model is refit on the first and
second halves of every unit's cases (rows are split in the order
supplied, which should be temporal order; a unit with an odd number of
cases contributes both possible splits, averaged) and the corrected
slope is `2 * full - mean(halves)`. Alternative-specific constants and
the unit dummies are left at their maximum-likelihood values, as are the
log-likelihood,
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html),
[`residuals()`](https://rdrr.io/r/stats/residuals.html), and the
standard errors; [`coef()`](https://rdrr.io/r/stats/coef.html),
[`predict()`](https://rdrr.io/r/stats/predict.html),
[`rrr()`](https://bagozzib.github.io/polychoice/reference/rrr.md),
[`confint()`](https://rdrr.io/r/stats/confint.html), and
[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md)
use the corrected slopes together with the ML constants. The ML slopes
are kept in `$uncorrected` and the half-panel estimates in `$jackknife`.
The identity requires the two halves to estimate the same parameter
(time homogeneity), so the correction is REFUSED, with a warning naming
the check and the ML fit returned, when the half-panel slopes disagree
beyond 2.5 times the mean disagreement of five placebo splits, random
halves of every unit's cases that ignore the time order (a regime
change), or when the full fit's residuals trend with the within-unit
time order for some alternative (an omitted trend, at the 1 percent
level with a Bonferroni correction over alternatives: put a time control
in the formula, a trend rather than period dummies, which are not
identified within the half panels). The placebo splits have the time
split's structure (all units, halves of the same lengths), so under time
homogeneity the ratio of the two disagreements is centred near 1 at
every panel length, while a mid-panel regime change separates the time
halves and is balanced in every random split; the splits come from fixed
seeds and do not disturb the caller's random-number stream. The constant
2.5 is the smallest of 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, and 4 that kept
the refusal rate of time-homogeneous panels at or below 5 percent at T =
6, 10, and 20 in a calibration Monte Carlo (150 units, J = 3), under a
rule written down before that run, and the choice is checked on fresh
seeds (`verification/jackknife_gate_mc.R` and `jackknife_gate_holdout.R`
in the repository, results in the verification ledger). In those
simulations a regime change was refused in 92 percent of the panels at T
= 6 and in all of them at T = 10 and 20, and an omitted trend in all
panels at T = 10 and 20 but only about half at T = 6. A refusal leaves
`$bias_correct` at `"none"` and the reason in `$jackknife_refused`,
which [`summary()`](https://rdrr.io/r/base/summary.html) prints, so that
a refused correction cannot pass for a corrected one. It warns when the
median panel length is below 6 or the implied correction is large. The
correction removes the leading \\1/T\\ term only: in the same
simulations it over-corrects at T = 6 and 10 (the corrected slopes keep
a fraction of the maximum-likelihood bias, with the opposite sign), yet
it lowers the root mean squared error and brings the coverage of Wald
intervals (with the maximum-likelihood standard errors) closer to
nominal at every T studied.

## Fit-quality flags

Every fit carries `converged` (`FALSE` when the optimizer stopped short
or at a saddle, where the log-likelihood rises along a direction in
which the Hessian curves upward; a near-zero curvature along which it
does not rise, as on the flat ridge of an unidentified model, is not a
saddle. At a saddle the fit restarts from the point where the
log-likelihood rose most along that direction, up to three times, with a
message when that reaches a higher maximum; the saddle warning stays
when it does not), `diverging` (the names of the diverging
coefficients), `phi_boundary` (an inflated fit's split equation at one
regime) and `saddle`, `boundary` (a structural parameter at zero or
infinity, a split equation that puts every case in one regime, or
diverging coefficients: separation, or a unit constant at infinity; in
`mnl()` and
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
fits without unit fixed effects separation is decided exactly by the
linear program of
[`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md),
elsewhere a coefficient counts as diverging when its standardized size
exceeds 15 and its Wald \|z\| is below 1.96), `ill_conditioned` (a
singular or nearly singular information matrix; `rcond` holds its
reciprocal condition number, taken on the matrix scaled to a unit
diagonal, and `singular` is `TRUE` when that number is below 1e-13 or
the matrix could not be inverted, and a generalized inverse took its
place), and `n_modes`, the number of distinct local maxima the search
found. Every end of the search (one for a single start; for the
multistart fits every polished start) is classified by the rules of
these flags: a point with a diverging coefficient, a structural
parameter at zero or infinity, or a split equation at one regime is a
boundary, not a maximum, and so are a saddle and an end that did not
converge. Two ends are the same point when their log-likelihoods differ
by less than 0.01 and either their estimates differ by less than a tenth
of a standard error or every fitted probability by less than 1e-4: ends
at other estimates with the same probabilities lie on a ridge where the
data do not identify the estimates, and `detail` says how many ends
joined the point that way. Maxima whose information matrix is
ill-conditioned are one point when their log-likelihoods differ by less
than 0.01, whatever their estimates (a nearly flat direction lets the
search stop anywhere along it), and `detail` counts those ends too; two
different maxima of that kind at one height would be counted once, and
only the highest end's estimates are kept. `search` lists the distinct
points (`loglik`, `type`: `"maximum"`, `"boundary"`, `"saddle"` or
`"not converged"`, `reached`: the number of ends there, `reported`,
`detail`), and the rows of `search_theta` hold their estimates on the
internal scale of `start`; `start_logliks` and `screen_logliks` hold the
log-likelihoods of every polished and every screened start. `n_modes` is
0 when no end is a local maximum (the reported fit then carries one of
the flags above). Each raises a warning or message at fit time;
[`summary.choicemodel()`](https://bagozzib.github.io/polychoice/reference/summary.choicemodel.md)
prints them and
[`glance.choicemodel()`](https://bagozzib.github.io/polychoice/reference/glance.choicemodel.md)
returns them. What they mean for the results: after `converged = FALSE`
neither the estimates nor the standard errors should be used (refit from
other `start` values, or simplify the model); with `boundary = TRUE` the
fit is at a limit (it keeps `converged = TRUE` when the optimizer met
its criterion there): the estimates of the parameters at the boundary
and their standard errors are not meaningful, the others usually are,
and the restricted model is the one to compare with (the search holds
log scales, log inclusive values, the multinomial probit's log Cholesky
diagonals and log standard deviations within +/- 8; where the
log-likelihood still rises at that bound the warning says so, and the
reported log-likelihood is then below the model's supremum, so a
comparison with the restricted model understates it);
`ill_conditioned = TRUE` (a reciprocal condition number below 1e-9, or a
generalized inverse in the covariance) means some combination of
parameters is barely identified and its standard errors are unreliable,
and with `singular = TRUE` some parameters are not identified at all:
their estimates and standard errors (from the generalized inverse) are
not meaningful, and the summary's footer says so, since their rows in
the table look like the others. The information of a model whose
parameters are not identified is singular at the maximum itself, and
close to it its condition number is of the order of the gradient. A fit
whose information is ill-conditioned is therefore brought to a
stationary point by a few Newton steps before the flag is set, so that
the flag does not depend on how closely the search converged;
`n_modes > 1` means the reported fit is the highest of several maxima
found. The reported fit is the highest end of the search, maximum or
not: where it is at a boundary and a lower end is a local maximum, a
message says so, and that maximum's row of `search_theta` refits it as
`start`. An alternative that is never chosen (where another is
available) is separation: in `mnl()` and
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
fits without unit fixed effects the fit runs and reports it
(`firth = TRUE` or a prior then gives finite estimates), while
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md),
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)
and fixed-effects fits refuse such data, and
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md)
does when the alternative has its own constant or slopes (with generic
coefficients only it refuses a nest of which nothing is chosen). Data
that overlap can still give very large, imprecise estimates (near
separation), which raise no flag: their Wald statistics understate the
evidence against zero (Hauck and Donner 1977), so the likelihood-ratio
tests of
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)
(`type = "lr"`) and
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
are the ones to use there, or `firth = TRUE`. A fixed-effects fit also
carries `bias_correct` and, when a requested correction was refused, the
reason in `jackknife_refused`.

## Penalized (Firth) estimation

Under separation (a covariate, or a combination of covariates, that
predicts some choices perfectly; Cook, Niehaus and Zuhlke 2018 for the
multinomial logit) the maximum-likelihood estimates do not exist, and
the reported ones are where the optimizer stopped. `firth = TRUE`
maximizes the log-likelihood penalized by half the log-determinant of
the Fisher information (Jeffreys' prior), \\\ell^\*(\theta) =
\ell(\theta) + \tfrac12 \log \|I(\theta)\|\\, which for the multinomial
and conditional logit is Firth's (1993) bias-reducing estimator (Bull,
Mak and Greenwood 2002; Kosmidis and Firth 2011): its estimates are
finite under separation and have smaller small-sample bias. The penalty
is invariant to the choice of base alternative. It shrinks the fitted
probabilities towards equal probabilities, so their averages no longer
equal the observed shares (Puhr et al. 2017).

A penalized fit reports `loglik`, the ordinary log-likelihood at the
penalized estimates, and `loglik_penalized`, the maximized objective;
the information criteria and pseudo-R-squared measures built on `loglik`
are not those of maximum likelihood. The standard errors are from the
inverse information at the penalized estimates (the robust and cluster
sandwiches use the per-case scores plus their share of the penalty's
gradient; the bootstrap refits with the penalty). The exact separation
check of
[`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md)
runs on the data: the names of the coefficients along the separating
direction are in `separation`, printed and announced with a message,
since their penalized estimates are then determined largely by the
penalty. The penalized likelihood is not always concave and can have
more than one local maximum (with a high-leverage covariate, even in a
binary logit), so without `start` the fit runs from zero and from the
maximum-likelihood estimate (halved towards zero where the penalty is
not finite) and keeps the higher maximum, with a warning when the two
runs end at different maxima (`start_logliks` holds both values,
`n_modes` counts the maxima, and a penalized fit keeps no `search`).
Recommended inference: the profile penalized-likelihood intervals of
[`confint()`](https://rdrr.io/r/stats/confint.html) (its default for a
penalized fit) and the penalized likelihood-ratio tests of
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md),
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
and
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
(two penalized fits), which maximize the full model's penalized
likelihood under the restriction (Heinze and Schemper 2001, 2002). The
tests derived for maximum-likelihood estimates refuse a penalized fit:
[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md),
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md),
[`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md),
and
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md).
`firth = TRUE` is not available with unit fixed effects, and the other
estimators do not offer it (in
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md) and
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)
the penalty tends to minus infinity where the model reduces to the
multinomial logit, so it would favour the larger model).

Under separation the penalty decides how large the separated estimates
are (Rainey 2016). Jeffreys' prior is one default, justified by bias
reduction; `prior = coef_prior(...)` replaces it with independent
normal, Cauchy or t priors centred at zero on coefficients scaled to
their columns (Gelman et al.'s 2008 default is
[`coef_prior()`](https://bagozzib.github.io/polychoice/reference/coef_prior.md)),
and fits under several priors show how much of a separated estimate
comes from the prior. A fit with a prior is a penalized fit in every
respect above (two starts, `separation`, profile intervals, penalized
likelihood ratios, the refusals); its covariance is the inverse
curvature of the log posterior at the mode. For maximum-likelihood fits
of separated data,
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md),
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
and
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
compare the suprema of the log-likelihoods (likelihood-ratio tests
remain valid under separation where Wald tests do not; Rainey 2024).

## References

Bull, S.B., Mak, C. and Greenwood, C.M.T. (2002). A modified score
function estimator for multinomial logistic regression in small samples.
*Computational Statistics & Data Analysis*, 39, 57-74.
[doi:10.1016/S0167-9473(01)00048-2](https://doi.org/10.1016/S0167-9473%2801%2900048-2)

Cook, S.J., Niehaus, J. and Zuhlke, S. (2018). A warning on separation
in multinomial logistic models. *Research & Politics*, 5(2).
[doi:10.1177/2053168018769510](https://doi.org/10.1177/2053168018769510)

Dhaene, G. and Jochmans, K. (2015). Split-panel jackknife estimation of
fixed-effect models. *Review of Economic Studies*, 82, 991-1030.
[doi:10.1093/restud/rdv007](https://doi.org/10.1093/restud/rdv007)

Firth, D. (1993). Bias reduction of maximum likelihood estimates.
*Biometrika*, 80, 27-38.
[doi:10.1093/biomet/80.1.27](https://doi.org/10.1093/biomet/80.1.27)

Gelman, A., Jakulin, A., Pittau, M.G. and Su, Y.-S. (2008). A weakly
informative default prior distribution for logistic and other regression
models. *Annals of Applied Statistics*, 2, 1360-1383.
[doi:10.1214/08-AOAS191](https://doi.org/10.1214/08-AOAS191)

Hauck, W.W. and Donner, A. (1977). Wald's test as applied to hypotheses
in logit analysis. *Journal of the American Statistical Association*,
72, 851-853.
[doi:10.1080/01621459.1977.10479969](https://doi.org/10.1080/01621459.1977.10479969)

Heinze, G. and Schemper, M. (2001). A solution to the problem of
monotone likelihood in Cox regression. *Biometrics*, 57, 114-119.
[doi:10.1111/j.0006-341X.2001.00114.x](https://doi.org/10.1111/j.0006-341X.2001.00114.x)

Heinze, G. and Schemper, M. (2002). A solution to the problem of
separation in logistic regression. *Statistics in Medicine*, 21,
2409-2419. [doi:10.1002/sim.1047](https://doi.org/10.1002/sim.1047)

Kosmidis, I. and Firth, D. (2011). Multinomial logit bias reduction via
the Poisson log-linear model. *Biometrika*, 98, 755-759.
[doi:10.1093/biomet/asr026](https://doi.org/10.1093/biomet/asr026)

McFadden, D. (1974). Conditional logit analysis of qualitative choice
behavior. In P. Zarembka (ed.), *Frontiers in Econometrics*, 105-142.
New York: Academic Press.

Puhr, R., Heinze, G., Nold, M., Lusa, L. and Geroldinger, A. (2017).
Firth's logistic regression with rare events: accurate effect estimates
AND predictions? *Statistics in Medicine*, 36, 2302-2317.
[doi:10.1002/sim.7273](https://doi.org/10.1002/sim.7273)

Rainey, C. (2016). Dealing with separation in logistic regression
models. *Political Analysis*, 24, 339-355.
[doi:10.1093/pan/mpw014](https://doi.org/10.1093/pan/mpw014)

Rainey, C. (2024). Hypothesis tests under separation. *Political
Analysis*, 32, 172-185.
[doi:10.1017/pan.2023.28](https://doi.org/10.1017/pan.2023.28)

Train, K.E. (2009). *Discrete Choice Methods with Simulation*, 2nd ed.
Cambridge University Press.

## See also

[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md),
[`mundlak()`](https://bagozzib.github.io/polychoice/reference/mundlak.md)

Other estimators:
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md),
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md),
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)

## Examples

``` r
set.seed(1)
d <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
m <- mnl(y ~ x1 + x2, data = d)
summary(m)
#> Multinomial logit 
#> Call: mnl(formula = y ~ x1 + x2, data = d)
#> Alternatives: a, b, c  (base: a)
#> 
#> Outcome equations (alternative:term):
#>               Estimate Std. Error z value Pr(>|z|)    
#> b:(Intercept)   0.4719     0.1236   3.817 0.000135 ***
#> b:x1           -0.8316     0.1487  -5.592 2.25e-08 ***
#> b:x2            0.7032     0.1314   5.352 8.72e-08 ***
#> c:(Intercept)  -0.3910     0.1625  -2.406 0.016124 *  
#> c:x1            0.8982     0.1622   5.538 3.05e-08 ***
#> c:x2           -0.6328     0.1432  -4.419 9.92e-06 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Standard errors: analytic. Log-likelihood -413.012 on 6 parameters, 500 cases. AIC 838.02, BIC 863.31.
head(predict(m))
#>           a         b          c
#> 1 0.2371501 0.6758323 0.08701760
#> 2 0.3247436 0.3626646 0.31259189
#> 3 0.3254302 0.4548467 0.21972314
#> 4 0.2356573 0.1010445 0.66329826
#> 5 0.2542398 0.6223149 0.12344533
#> 6 0.0921925 0.8969239 0.01088358
rrr(m)
#>           term   block   rrr lower upper      z  p.value
#>  b:(Intercept) outcome 1.603 1.258 2.043  3.817 0.000135
#>           b:x1 outcome 0.435 0.325 0.583 -5.592 2.25e-08
#>           b:x2 outcome 2.020 1.562 2.614  5.352 8.72e-08
#>  c:(Intercept) outcome 0.676 0.492 0.930 -2.406 0.016124
#>           c:x1 outcome 2.455 1.787 3.374  5.538 3.05e-08
#>           c:x2 outcome 0.531 0.401 0.703 -4.419 9.92e-06
## Firth's penalized likelihood, with profile intervals
mf <- mnl(y ~ x1 + x2, data = d, firth = TRUE)
confint(mf)
#>                    2.5 %      97.5 %
#> b:(Intercept)  0.2284361  0.71111992
#> b:x1          -1.1192064 -0.53902765
#> b:x2           0.4440831  0.95665219
#> c:(Intercept) -0.7070561 -0.07340702
#> c:x1           0.5781281  1.21027365
#> c:x2          -0.9091239 -0.35080670
```
