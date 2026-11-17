# Heteroskedastic extreme value regression

Fits the heteroskedastic extreme value (HEV) model of Bhat (1995):
independent Gumbel errors whose scale \\\theta_j\\ differs across
alternatives (the base alternative's scale is fixed at 1), so that the
substitution pattern is no longer proportional and IIA is relaxed
through the scales rather than through nests.

## Usage

``` r
hev(
  formula,
  data,
  id = NULL,
  alt = NULL,
  ref = NULL,
  weights = NULL,
  offset = NULL,
  fe = NULL,
  se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"),
  cluster = NULL,
  nboot = 200L,
  cores = 1L,
  start = NULL,
  maxit = 1000L,
  reltol = 1e-10
)
```

## Arguments

- formula:

  A formula; see the section "Two data layouts" of
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md).

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
  fit is on that scale. Rarely needed for
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
  whose log-likelihood is globally concave; the other estimators start
  from the
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  estimates (and
  [`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md)
  from a multistart grid), so `start` serves to probe another basin of a
  non-concave likelihood.

- maxit, reltol:

  Controls for the BFGS stage passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) (the plain
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) is
  solved by damped Newton steps first, with BFGS as the fallback); the
  Newton polish that follows runs up to 25 further iterations with the
  analytic Hessian of the logit family or central differences of the
  analytic gradient otherwise.

## Value

An object of class `c("hev", "choicemodel")`; as
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md), with
the scale parameters named `"scale:<alternative>"`.

## Details

The choice probability is a one-dimensional integral over the standard
Gumbel variable \\w = \epsilon_j/\theta_j\\, \$\$P(y_i = j) = \int
\exp\Big(-\sum\_{k \ne j} e^{-(V\_{ij} - V\_{ik})/\theta_k -
(\theta_j/\theta_k) w}\Big) e^{-w - e^{-w}}\\ dw,\$\$ evaluated on a
window around the mode of the (log-concave) integrand: by the
trapezoidal rule when every scale ratio \\\theta_j/\theta_k\\ is at most
20 (the integrand is analytic with exponential tails, so the rule
converges spectrally), and by adaptive Gauss-Kronrod quadrature beyond
that ratio, where a competitor with a much smaller scale makes the
integrand a near step (the initial panels carry breakpoints around each
step, so the log-likelihood and its gradient are smooth in the
parameters); probabilities agree with
[`integrate()`](https://rdrr.io/r/stats/integrate.html) to within 1e-9
(the test tolerance) for scale ratios up to 1.6e5. With every \\\theta_j
= 1\\ the model is the multinomial logit (reproduced exactly);
`lr_test(mnl_fit, hev_fit)` tests homoskedasticity. The scales are
estimated on the log scale and reported on the natural scale with
delta-method standard errors; the z statistics in the summary test
\\\theta = 1\\, and [`confint()`](https://rdrr.io/r/stats/confint.html)
gives intervals on the log scale transformed back. A scale that runs to
zero (the alternative's utility becomes deterministic relative to the
others) is flagged as a boundary fit with the information matrix
singular; the likelihood is then usually flat in that direction and the
standard errors are not interpretable (the travel-mode data of the
vignette are an example). Where it still rises at the search's bound (a
log scale of +/- 8), the warning says so, and the reported
log-likelihood is below the model's supremum. With two alternatives the
single scale is identified only through the functional form and a
warning says so. Everything else is as in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md).

## References

Bhat, C.R. (1995). A heteroscedastic extreme value model of intercity
travel mode choice. *Transportation Research Part B*, 29, 471-483.
[doi:10.1016/0191-2615(95)00015-6](https://doi.org/10.1016/0191-2615%2895%2900015-6)

## See also

[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)

Other estimators:
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md),
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md),
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)

## Examples

``` r
set.seed(4)
d <- rhev(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
          scale = c(1, 0.6, 1.5), gamma = 1)   # w1 identifies the scales
m <- hev(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
summary(m)
#> Heteroskedastic extreme value model 
#> Call: hev(formula = chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
#> Alternatives: a, b, c  (base: a)
#> 
#> Outcome equations (alternative:term):
#>               Estimate Std. Error z value Pr(>|z|)   
#> b:(Intercept)  0.27424    0.36105   0.760  0.44753   
#> b:x1          -1.02366    0.36905  -2.774  0.00554 **
#> b:x2           1.16335    0.36183   3.215  0.00130 **
#> c:(Intercept)  0.02232    0.34317   0.065  0.94813   
#> c:x1           0.61771    0.28267   2.185  0.02887 * 
#> c:x2          -0.13493    0.22185  -0.608  0.54305   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Generic coefficients (alternative-specific covariates):
#>    Estimate Std. Error z value Pr(>|z|)    
#> w1   1.0205     0.2409   4.236 2.27e-05 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Gumbel scale parameters (base alternative = 1; z tests are against 1; the joint test is lr_test() against the mnl() fit):
#>         Estimate Std. Error z value Pr(>|z|)
#> scale:b   0.9384     0.4648  -0.132    0.895
#> scale:c   1.0640     0.3525   0.182    0.856
#> 
#> Standard errors: analytic. Log-likelihood -217.972 on 9 parameters, 300 cases. AIC 453.94, BIC 487.28.
```
