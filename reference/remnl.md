# Multinomial logit with random alternative-specific intercepts

Fits the multinomial (or conditional) logit of
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) with a
normally distributed random intercept per unit in every non-base
equation: a mixed logit whose only random coefficients are the
constants, for panels and other grouped data where the cases of a unit
share unobserved tastes for the alternatives.

## Usage

``` r
remnl(
  formula,
  data,
  unit,
  nodes = 12L,
  id = NULL,
  alt = NULL,
  ref = NULL,
  weights = NULL,
  offset = NULL,
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

- unit:

  Column name in `data` identifying the units whose cases share the
  random intercepts.

- nodes:

  Gauss-Hermite points per dimension (default 12). Because the rule is
  adaptive its error does not grow with the number of cases per unit:
  against `lme4::glmer(nAGQ = 25)` at \\J = 2\\ (150 units, 6 to 50
  cases each, sd 1 or 2) 12 points put the log-likelihood within 1.2e-3
  and the standard deviation within 0.01 percent, 8 points within 0.02
  and 0.15 percent (`verification/remnl_adaptive_check.R` in the
  development repository). Raise it for large standard deviations with
  short panels. The analytic gradient is the quadrature of the exact
  gradient, which differs from the derivative of the quadrature
  likelihood by the quadrature error; with few nodes (6 or fewer), or
  large standard deviations and short panels, that difference can keep
  the optimizer from meeting its criterion (`converged = FALSE` with a
  small gradient), and more nodes cure it. The population-averaged
  probabilities of
  [`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
  [`predict()`](https://rdrr.io/r/stats/predict.html), an integral over
  the prior that no data narrow, use a finer rule: Gauss-Hermite with
  \\\max(\mathrm{nodes}, 10\sigma^2)\\ points per dimension (rounded up
  to a power of two) while that is at most 128, and beyond that the
  trapezoidal rule with a step of at most \\1/\sigma\\ standard
  deviations, whose size grows only linearly in \\\sigma\\; both are
  accurate to about \\10^{-7}\\ or better. The rule has at most 50000
  points in all, or \\\mathrm{nodes}^{J-1}\\ if that is larger, and at
  the cap the rule with the smaller error bound is used. The cap binds
  at \\J = 3\\ for \\\sigma\\ above about 7 (the error bound passes
  \\10^{-6}\\ near 18), at \\J = 4\\ above about 1.8 (bound above
  \\10^{-6}\\ from 2.7), and at \\J \ge 5\\ above about 0.9 (from about
  1.6); the fit warns, with the bound and the number of points, whenever
  the bound exceeds \\10^{-6}\\.

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

- se:

  Standard errors: `"analytic"` (default), `"opg"`, `"robust"` (over
  units, the observation of this likelihood), `"cluster"` (over clusters
  that nest the units), `"bootstrap"` (units resampled, a unit drawn
  twice counting as two units with their own intercepts; with frequency
  weights, multinomial unit counts), or `"none"`.

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

An object of class `c("remnl", "choicemodel")`.

## Details

\$\$V\_{ij} = o\_{ij} + x_i^\top \beta_j + w\_{ij}^\top \gamma +
u\_{g(i) j}, \quad u_g \sim N(0, \mathrm{diag}(\sigma_j^2)),\$\$ with
\\u\_{g b} = 0\\ for the base alternative and independence across units
and alternatives. With three or more alternatives the base defines the
model, not only the labels. Intercepts independent relative to one base
are correlated relative to another (\\u_j - u_c\\ and \\u_k - u_c\\
share \\u_c\\), which the diagonal model cannot represent. Fits with
different `ref` are therefore different models, with different
log-likelihoods and different tests against the multinomial logit. In
the example below the log-likelihood is -488.49 with base `a` and
-490.50 with base `b` (a fit whose standard deviation for `a` runs to
zero, with a boundary warning), and the likelihood ratio against
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) is
4.53 and 0.50. Choose the base before testing, on substantive grounds:
the alternative whose utility is least expected to vary across units (an
outside option, abstention). Fits with different bases have the same
number of parameters and can be compared descriptively by their
log-likelihoods, but a base chosen that way makes the test against the
multinomial logit a maximum over the bases, for which the chi-square
reference does not hold. The likelihood of unit \\g\\ integrates the
product of the MNL probabilities of its cases over \\u_g\\ by adaptive
product Gauss-Hermite quadrature (Liu and Pierce 1994; Pinheiro and
Bates 1995): the nodes are centred at the mode of each unit's integrand
and scaled by its curvature, with `nodes` points per dimension (\\J -
1\\ dimensions, so \\\mathrm{nodes}^{J-1}\\ points per unit; a warning
above 20000 points). The standard deviations are reported as
`"sd:<alt>"` (estimated on the log scale, delta-method standard errors).
The observation is the unit: `object$lli` has one entry per unit, robust
and outer-product standard errors are over units, `cluster` must nest
units, frequency weights must be constant within a unit, and
`se = "bootstrap"` resamples units. Fitted and predicted probabilities
are population-averaged (integrated over the random intercepts).
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
against the
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) fit of
the same covariates tests \\\sigma = 0\\, a boundary restriction whose
chi-square p-value is conservative;
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)
against an MNL aggregates its contributions over units. At \\J = 2\\ the
model is the binomial generalized linear mixed model with a random
intercept, and the fit matches
[`lme4::glmer()`](https://rdrr.io/pkg/lme4/man/glmer.html) with adaptive
quadrature. Everything else is as in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
except that there is no `fe` or `bias_correct` (the unit effects are the
random intercepts) and that
[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws one
intercept vector per unit and then the unit's cases given it, while
[`fitted()`](https://rdrr.io/r/stats/fitted.values.html) and
[`predict()`](https://rdrr.io/r/stats/predict.html) are the
population-averaged probabilities. The IIA tests and
[`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md)
refuse these fits (the cases of a unit are dependent).
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md)
and the information criteria use the number of cases as \\n\\.

## References

Liu, Q. and Pierce, D.A. (1994). A note on Gauss-Hermite quadrature.
*Biometrika*, 81(3), 624-629.
[doi:10.1093/biomet/81.3.624](https://doi.org/10.1093/biomet/81.3.624)

Pinheiro, J.C. and Bates, D.M. (1995). Approximations to the
log-likelihood function in the nonlinear mixed-effects model. *Journal
of Computational and Graphical Statistics*, 4(1), 12-35.
[doi:10.1080/10618600.1995.10474663](https://doi.org/10.1080/10618600.1995.10474663)

## See also

[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)

Other estimators:
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md),
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)

## Examples

``` r
set.seed(1)
G <- 60; T <- 5; unit <- rep(seq_len(G), each = T)
x1 <- rnorm(G * T); u <- rnorm(G, 0, 1)
V <- cbind(0, 0.3 + 0.8 * x1 + u[unit], -0.2 - 0.6 * x1)
P <- exp(V) / rowSums(exp(V))
y <- apply(P, 1, function(p) sample.int(3, 1, prob = p))
d <- data.frame(y = factor(c("a", "b", "c")[y]), x1, unit = factor(unit))
m <- remnl(y ~ x1, data = d, unit = "unit", nodes = 8)
summary(m)
#> Multinomial logit with random alternative-specific intercepts (60 units, 8-point Gauss-Hermite per dimension) 
#> Call: remnl(formula = y ~ x1, data = d, unit = "unit", nodes = 8)
#> Alternatives: a, b, c  (base: a)
#> 
#> Outcome equations (alternative:term):
#>               Estimate Std. Error z value Pr(>|z|)    
#> b:(Intercept)  0.56630    0.16204   3.495 0.000474 ***
#> b:x1           0.73601    0.18854   3.904 9.47e-05 ***
#> c:(Intercept) -0.07182    0.19165  -0.375 0.707834    
#> c:x1          -0.65526    0.20311  -3.226 0.001255 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Random alternative-specific intercepts: standard deviations (no z test: sd = 0 is on the boundary; use lr_test() against the mnl() fit):
#>      Estimate Std. Error z value Pr(>|z|)
#> sd:b   0.4065     0.3232      NA       NA
#> sd:c   0.3018     0.4706      NA       NA
#> 
#> Standard errors: analytic. Log-likelihood -282.432 on 6 parameters, 300 cases. AIC 576.86, BIC 599.09.
lr_test(mnl(y ~ x1, data = d), m)
#> The restriction sd = 0 is on the boundary of the parameter space: the chi-square p-value is conservative (the null distribution is a mixture of chi-squares with fewer degrees of freedom); boot_lr_test() gives a calibrated p-value.
#> Likelihood-ratio test
#>   LR = 0.8967 on 2 df, p = 0.6387
```
