# Nested logit regression

Fits the two-level nested logit (McFadden 1978) in its
utility-maximization- consistent form, the one mlogit and Stata's
`nlogit` (RUM-consistent) use: alternatives are partitioned into nests,
and each nest \\m\\ has an inclusive-value parameter \\\lambda_m\\ that
governs the correlation of the errors within the nest (\\\lambda_m =
1\\: no correlation, the multinomial logit; values in \\(0, 1)\\ are
consistent with utility maximization for all covariate values).

## Usage

``` r
nlogit(
  formula,
  data,
  nests,
  equal_lambda = FALSE,
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

- nests:

  A named list of character vectors assigning every alternative
  (response level) to exactly one nest, e.g.
  `list(fast = c("air", "train"), slow = c("bus", "car"))`.

- equal_lambda:

  If `TRUE`, one inclusive-value parameter is shared by all
  (non-singleton) nests.

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

An object of class `c("nlogit", "choicemodel")`; as
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md), with
the inclusive-value parameters named `"lambda:<nest>"` (or `"lambda"`).
`predict(type = "prob_conditional")` gives the within-nest
probabilities.

## Details

With \\V\_{ij}\\ as in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) and
\\I\_{im} = \log\sum\_{k \in B_m} e^{V\_{ik}/\lambda_m}\\ over the
alternatives of nest \\m\\ available to case \\i\\, \$\$P(y_i = j) =
\frac{e^{V\_{ij}/\lambda_m}}{\sum\_{k \in B_m} e^{V\_{ik}/\lambda_m}}
\cdot \frac{e^{\lambda_m I\_{im}}}{\sum_l e^{\lambda_l I\_{il}}}, \qquad
j \in B_m.\$\$ The inclusive-value parameters are estimated on the log
scale and reported on the natural scale with delta-method standard
errors; the z statistics of the summary test \\\lambda = 1\\ (the
multinomial logit within that nest). Test IIA jointly with
`lr_test(mnl_fit, nlogit_fit)` (the restriction \\\lambda = 1\\ is
interior). A singleton nest's parameter does not enter the likelihood
and is fixed at 1 (with a message). Everything else – layouts, weights,
offsets, fixed effects, standard errors – is as in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md).
Starting values come from the multinomial logit fit; the fit is run from
\\\lambda = 1\\ and \\\lambda = 0.7\\ and the better optimum kept.

## References

McFadden, D. (1978). Modelling the choice of residential location. In A.
Karlqvist et al. (eds.), *Spatial Interaction Theory and Planning
Models*, 75-96. North-Holland.

Train, K.E. (2009). *Discrete Choice Methods with Simulation*, 2nd ed.,
chapter 4.

## See also

[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)

Other estimators:
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md),
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md),
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)

## Examples

``` r
set.seed(3)
d <- rnlogit(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)),
             nests = list(n1 = c("a", "b"), n2 = c("c", "d")), lambda = c(0.6, 0.8))
m <- nlogit(y ~ x1 + x2, data = d, nests = list(n1 = c("a", "b"), n2 = c("c", "d")))
#> Inclusive-value parameter(s) above 1 (lambda:n2): outside the range consistent with random-utility maximization for all covariate values; the nesting may be wrong or unnecessary (lr_test() against the mnl() fit).
summary(m)
#> Nested logit (2 nests) 
#> Call: nlogit(formula = y ~ x1 + x2, data = d, nests = list(n1 = c("a", 
#>     "b"), n2 = c("c", "d")))
#> Alternatives: a, b, c, d  (base: a)
#> 
#> Outcome equations (alternative:term):
#>               Estimate Std. Error z value Pr(>|z|)    
#> b:(Intercept)   0.3304     0.1450   2.279 0.022667 *  
#> b:x1           -0.7676     0.3098  -2.478 0.013228 *  
#> b:x2            0.6039     0.2532   2.385 0.017060 *  
#> c:(Intercept)  -1.0516     0.6864  -1.532 0.125518    
#> c:x1            0.7067     0.1873   3.774 0.000161 ***
#> c:x2           -0.9434     0.3282  -2.874 0.004047 ** 
#> d:(Intercept)  -0.3608     0.3945  -0.915 0.360305    
#> d:x1            0.6553     0.2038   3.216 0.001301 ** 
#> d:x2            0.2190     0.3173   0.690 0.489955    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Inclusive-value parameters lambda (z tests are against 1 = multinomial logit; the joint test is lr_test() against the mnl() fit):
#>           Estimate Std. Error z value Pr(>|z|)    
#> lambda:n1   0.4146     0.1751  -3.343 0.000829 ***
#> lambda:n2   1.3286     0.6937   0.474 0.635727    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Standard errors: analytic. Log-likelihood -1554.187 on 11 parameters, 1500 cases. AIC 3130.37, BIC 3188.82.
lr_test(mnl(y ~ x1 + x2, data = d), m)
#> Likelihood-ratio test
#>   LR = 10.5298 on 2 df, p = 0.00517
```
