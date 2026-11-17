# Multinomial probit regression

Fits the multinomial probit (Hausman and Wise 1978; Daganzo 1979):
utilities \\U\_{ij} = V\_{ij} + \epsilon\_{ij}\\ with multivariate
normal errors whose covariance is estimated, so that any pattern of
correlation and heteroskedasticity across alternatives, and hence any
substitution pattern, is admitted.

## Usage

``` r
mnprobit(
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
  reltol = 1e-10,
  draws = 200L
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

- draws:

  Number of GHK draws per case for four or more alternatives (default
  200). The simulator's Halton sequences use one prime per dimension of
  the differenced utilities, the first 20, so the model takes at most 21
  alternatives.

## Value

An object of class `c("mnprobit", "choicemodel")`.

## Details

Only utility differences are identified. The model is parameterized
through the covariance of the differences from the base alternative,
\\\eta_k = \epsilon_k - \epsilon_b\\, \\\tilde\Sigma = LL'\\ with \\L\\
lower triangular, its first element fixed at 1 (the scale normalization)
and the remaining \\J(J-1)/2 - 1\\ elements free (Train 2009, section
5.2; the same normalization as mlogit). The reported structural
parameters are the elements of \\\tilde\Sigma\\ named `"var:<alt>"` and
`"cov:<alt>,<alt>"`, with delta-method standard errors.

Probabilities: with three alternatives the exact bivariate normal (the
Drezner-Genz algorithm, absolute error about 1e-15); with more, the GHK
simulator with `draws` Halton draws per case, held fixed across
evaluations (the simulated log-likelihood is smooth away from a singular
covariance; the simulation error of a probability has standard deviation
of order \\1/\sqrt{R}\\ and the bias of its logarithm is of order
\\1/R\\, both vanishing with the number of draws; GHK probabilities of a
case need not sum exactly to one, and the draws of a case are tied to
its row position, so predictions for the same case in another row use
other draws). Gradients are analytic throughout, including the
covariance parameters (forward-mode differentiation through the Cholesky
factor and the GHK recursion); the Hessian is obtained by central
differences of the analytic gradient. The multinomial probit does not
nest the multinomial logit, so compare them with
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md) or
the information criteria. Multinomial probit coefficients are on a
smaller scale than multinomial logit coefficients (the error-difference
variance is normalized to 1 rather than \\\pi^2/3\\). Starting values
come from the multinomial logit fit (scaled) with the covariance implied
by iid errors, and from independent differences. Everything else is as
in [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
except that the returned object has class `"mnprobit"` (not `"mnp"`,
which the CRAN package MNP uses).

Identification. With case-specific covariates only (alternative-specific
coefficients, the wide layout) the covariance parameters are identified
solely through the functional form and are fragile (Keane 1992): the
likelihood is nearly flat along directions where an alternative's
coefficients and its error variance shrink together, fits may stop at a
near-singular covariance with a convergence warning, and the covariance
estimates are not to be interpreted. Alternative-specific covariates
with generic coefficients (the long layout) identify the covariance; a
fit with a variance running to zero is flagged as a boundary fit. Near a
singular covariance the GHK simulated likelihood is not smooth: with the
draws held fixed, one draw's contribution can jump over a small change
in the parameters, so a boundary fit with four or more alternatives may
stop with a large gradient, and the boundary warning says that the
convergence check is unreliable there.

## References

Daganzo, C. (1979). *Multinomial Probit: The Theory and Its Application
to Demand Forecasting*. New York: Academic Press.

Hausman, J.A. and Wise, D.A. (1978). A conditional probit model for
qualitative choice. *Econometrica*, 46, 403-426.
[doi:10.2307/1913909](https://doi.org/10.2307/1913909)

Train, K.E. (2009). *Discrete Choice Methods with Simulation*, 2nd ed.,
chapter 5.

Keane, M.P. (1992). A note on identification in the multinomial probit
model. *Journal of Business and Economic Statistics*, 10, 193-200.
[doi:10.1080/07350015.1992.10509898](https://doi.org/10.1080/07350015.1992.10509898)

## See also

[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)

Other estimators:
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md),
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)

## Examples

``` r
set.seed(5)
d <- rmnp(1500, beta = cbind(b = c(0.3, -0.6, 0.5), c = c(-0.2, 0.4, -0.3)),
          sigma = matrix(c(1, 0.6, 0.6, 1.5), 2), gamma = 1)
m <- mnprobit(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
summary(m)
#> Multinomial probit (exact bivariate normal) 
#> Call: mnprobit(formula = chosen ~ w1 | x1 + x2, data = d, id = "id", 
#>     alt = "alt")
#> Alternatives: a, b, c  (base: a)
#> 
#> Outcome equations (alternative:term):
#>               Estimate Std. Error z value Pr(>|z|)    
#> b:(Intercept)  0.27346    0.06630   4.124 3.72e-05 ***
#> b:x1          -0.71473    0.06997 -10.215  < 2e-16 ***
#> b:x2           0.67302    0.06256  10.757  < 2e-16 ***
#> c:(Intercept) -0.42931    0.10330  -4.156 3.24e-05 ***
#> c:x1           0.66639    0.09850   6.765 1.33e-11 ***
#> c:x2          -0.28412    0.07951  -3.573 0.000352 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Generic coefficients (alternative-specific covariates):
#>    Estimate Std. Error z value Pr(>|z|)    
#> w1   1.2074     0.0749   16.12   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Covariance of the utility differences from the base (first variance fixed at 1; z tests are against the iid values, variances 1 and covariances 0.5):
#>         Estimate Std. Error z value Pr(>|z|)   
#> cov:c,b   0.4836     0.1832   -0.09  0.92845   
#> var:c     2.4069     0.4467    3.15  0.00163 **
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Standard errors: analytic. Log-likelihood -856.215 on 9 parameters, 1500 cases. AIC 1730.43, BIC 1778.25.
```
