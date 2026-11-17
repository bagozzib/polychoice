# Fok-Paap tests of multinomial logit misspecification

The two tests of Fok and Paap (2025), built on the property that under
the multinomial logit the choice between any two alternatives \\j, m\\,
given that one of them is chosen, is a binary logit in the utility
difference \\V_j - V_m\\. Pairs are the \\J - 1\\ base pairs or all
\\J(J-1)/2\\ pairs.

## Usage

``` r
fp_tests(
  object,
  pairs = c("all", "base"),
  tests = c("cl", "gmm"),
  gmm = c("two-step", "iterated", "projected"),
  vcov = c("joint", "difference"),
  maxit = 50L
)
```

## Arguments

- object:

  A fitted
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  model with at least three alternatives. A fit with diverging
  coefficients (separation) draws a warning. Along a diverging direction
  neither estimator has an estimate, so \\Q_1\\ leaves those directions
  out: the composite-likelihood estimator keeps the fit's values along
  them, the two estimators are compared on the other directions, and the
  degrees of freedom are reduced by their number
  (`attr(, "details")$q1_left_out` names the coefficients that move
  along them and `$q1_left_out_rank` counts the directions; a separated
  reference cell is one direction through an intercept and a dummy,
  whose sum stays in the comparison). The statistic then does not depend
  on where the optimizer stopped. \\Q_2\\ is not computed (`NA`): its
  moment conditions have no finite solution along a diverging direction,
  and the moments at the point where the optimizer stopped depend on
  that point. In
  [`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
  a replicate that separates has fewer degrees of freedom than the
  observed statistic it is compared with.

- pairs:

  `"all"` (default) or `"base"`.

- tests:

  Which tests to compute: `"cl"`, `"gmm"`, or both (default).

- gmm:

  `"two-step"` (default), `"iterated"`, or `"projected"`; see Details.

- vcov:

  The covariance in \\Q_1\\: `"joint"` (default) or `"difference"`; see
  Details.

- maxit:

  Maximum number of weight-matrix updates for `gmm = "iterated"`.

## Value

A data frame of class `"iia_test"` with one row per test: `test`,
`statistic`, `df`, `p.value`, `psd`, and `pairs`. The attribute
`"details"` holds the CL and GMM estimates and the reciprocal condition
number of the moment covariance; on a separated fit also `q1_left_out`
and `q1_left_out_rank`, and `theta_cl` then equals the fit along the
directions left out, where `vcov_cl` has no variance.

## Details

**Composite-likelihood Hausman test (Q1).** The pairwise composite
likelihood estimator \\\hat\theta\_{CL}\\ maximizes the sum over pairs
and cases of the binary logit log-likelihoods; its covariance is the
Godambe sandwich. Under correct specification the MLE and the CL
estimator are both consistent and only the MLE is efficient, so \\Q_1 =
(\hat\theta\_{ML} - \hat\theta\_{CL})'(V\_{CL} -
V\_{ML})^{-1}(\hat\theta\_{ML} - \hat\theta\_{CL})\\ is chi-square on
\\\dim\theta\\ degrees of freedom (F&P eq. 18). Under the multinomial
logit the composite-likelihood estimator is nearly efficient (its
variance is within a few percent of the MLE's in most directions), so
\\V\_{CL} - V\_{ML}\\ is a small difference of two separately estimated
matrices, and the sampling noise of the sandwich makes it indefinite in
nearly every sample: in the package's simulations (`verification/` in
the development repository) it was never positive semi-definite, 44
percent of the statistics were negative, and the rejection rates were
0.033, 0.042, and 0.049 at the 1, 5, and 10 percent levels. The default
`vcov = "joint"` therefore estimates the covariance of
\\\hat\theta\_{ML} - \hat\theta\_{CL}\\ directly, from the stacked
influence functions \\V\_{ML} s^{ML}\_i - H\_{CL}^{-1} s^{CL}\_i\\ of
the two estimators (summed by cluster when the fit has a `cluster`), as
the joint form of
[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md)
does; it is positive semi-definite by construction and its rejection
rates in the same simulations were 0.015, 0.062, and 0.112.
`vcov = "difference"` is the formula of F&P eq. 18, for comparison with
published values; `psd` reports whether its covariance difference was
positive semi-definite, and the p-value refers `max(statistic, 0)` to
the chi-square (a negative statistic is no evidence against the model).

**GMM overidentification test (Q2).** With moment conditions \\E\[(1\\y
= j\\ - \Lambda(V_j - V_m))\\ 1\\y \in \\j, m\\\\\\ \zeta\] = 0\\ for
every pair, instruments \\\zeta = (x, w_j - w_m)\\, the parameters are
overidentified and Hansen's J statistic at the efficient GMM estimate is
chi-square on (number of moments \\-\\ number of parameters) degrees of
freedom (F&P eq. 19). `gmm = "two-step"` (default) is Hansen's two-step
estimator: the weight matrix is the inverse empirical covariance of the
moments at the MLE (F&P eq. 15), the parameters are re-estimated under
it, and the statistic is evaluated at that estimate with its own weight
matrix. `gmm = "iterated"` repeats the update until the estimate stops
moving, as Fok and Paap do; with many moments relative to the sample the
iteration can fail to settle (a warning is issued and the last iterate
is used), which is why it is not the default. `gmm = "projected"`
evaluates the moments at the MLE and projects out the estimated
directions (Newey's conditional-moment / C form), which has the same
asymptotic distribution without a GMM optimization. All three have the
same degrees of freedom. With case-specific covariates only, the
overidentifying restrictions carry information only through the
covariates (Fok and Paap 2025, Section 2.2): with few or weak covariates
the moment covariance is near singular and the test has little power;
its reciprocal condition number is reported. Where the moments are
linearly dependent (one binary covariate, say) the moment covariance is
singular: the statistic is then computed on its generalized inverse, its
degrees of freedom are the rank of the moment covariance minus the
number of parameters, and a warning says so.

**Units of the covariates.** Where they are computed, neither statistic
depends on the units of the covariates, on their centering, or on a
factor's reference level. The instruments of \\Q_2\\ are put on a
canonical scale before the moments are formed, and the contrast of
\\Q_1\\ is whitened by the covariance of the maximum-likelihood
estimator, so that the directions a generalized inverse drops are the
same in any parametrization.

**When the tests are not computed.** Both tests compare estimators at
the fit. Where the fit's information matrix is ill-conditioned on the
directions compared (reciprocal condition number of its normalized form
below 1e-9, the rule of the fit's `ill_conditioned` flag), the estimates
are determined to few digits along some direction and the statistics
would depend on where the optimizer stopped: both are `NA`, with a
warning. That condition number does not depend on the covariates' units,
but a covariate far from its origin (a calendar year) is nearly
collinear with the constant and can bring it below the bound; centering
it helps. \\Q_1\\ is also `NA` when the composite-likelihood estimator
does not exist (its sensitivity matrix is singular).

## References

Fok, D. and Paap, R. (2025). New misspecification tests for multinomial
logit models. *Journal of Choice Modelling*, 54, 100531.
[doi:10.1016/j.jocm.2024.100531](https://doi.org/10.1016/j.jocm.2024.100531)

Varin, C., Reid, N. and Firth, D. (2011). An overview of composite
likelihood methods. *Statistica Sinica*, 21, 5-42.

Hansen, L.P. (1982). Large sample properties of generalized method of
moments estimators. *Econometrica*, 50, 1029-1054.
[doi:10.2307/1912775](https://doi.org/10.2307/1912775)

## See also

[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)

Other IIA tests:
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md),
[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md),
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md)

## Examples

``` r
set.seed(1)
d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
m <- mnl(y ~ x1 + x2, data = d)
fp_tests(m)
#> Fok-Paap (2025) misspecification tests, all pairs 
#>                                            test statistic df p.value  psd pairs
#>                 Fok-Paap CL Hausman (Q1, joint)     8.593  9   0.476 TRUE   all
#>  Fok-Paap GMM overidentification (Q2, two-step)    18.101  9   0.034   NA   all
fp_tests(m, gmm = "projected")
#> Fok-Paap (2025) misspecification tests, all pairs 
#>                                             test statistic df p.value  psd
#>                  Fok-Paap CL Hausman (Q1, joint)     8.593  9  0.4757 TRUE
#>  Fok-Paap GMM overidentification (Q2, projected)    15.346  9  0.0819   NA
#>  pairs
#>    all
#>    all
fp_tests(m, gmm = "iterated")
#> Warning: The iterated GMM estimator did not settle within 50 weight-matrix updates; the last iterate is used. Consider gmm = "two-step" or "projected".
#> Fok-Paap (2025) misspecification tests, all pairs 
#>                                            test statistic df p.value  psd pairs
#>                 Fok-Paap CL Hausman (Q1, joint)     8.593  9  0.4757 TRUE   all
#>  Fok-Paap GMM overidentification (Q2, iterated)    19.204  9  0.0235   NA   all
```
