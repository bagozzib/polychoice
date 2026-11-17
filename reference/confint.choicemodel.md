# Confidence intervals for a choice model

Wald intervals from the fit's covariance, percentile intervals from the
bootstrap replicates when the model was fit with `se = "bootstrap"`, or,
for a penalized fit (`firth = TRUE` or a `prior` in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) and
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)),
profile penalized-likelihood intervals. The Wald intervals of the
positive structural parameters (inclusive values, scales, standard
deviations, variances) are computed on the log scale and transformed
back, so they stay positive. With unit fixed effects and the bootstrap,
the alternative-specific constants and the unit effects have no
replicates and their intervals are `NA`; percentile intervals are
refused for a jackknife-corrected fit (the replicates are not
corrected).

## Usage

``` r
# S3 method for class 'choicemodel'
confint(
  object,
  parm,
  level = 0.95,
  type = c("wald", "percentile", "profile"),
  ...
)
```

## Arguments

- object:

  A `"choicemodel"` fit.

- parm:

  Parameters (names or positions); default all.

- level:

  Confidence level.

- type:

  `"wald"`, `"percentile"` (bootstrap fits only), or `"profile"`
  (penalized fits only). Default: `"profile"` for a penalized fit,
  `"wald"` otherwise.

- ...:

  Unused.

## Value

A two-column matrix.

## Details

Profile intervals (Heinze and Schemper 2002) hold the values \\c\\ of a
coefficient for which twice the drop of the penalized log-likelihood,
maximized over the other parameters with the coefficient held at \\c\\,
is at most the chi-square quantile on one degree of freedom (the full
model's penalty throughout). Under separation they are asymmetric and
wider than the Wald intervals on the side towards which the
maximum-likelihood estimate runs off. A limit not reached within 50
standard errors is `NA`, with a warning. With few cases per coefficient,
and under separation, the penalized likelihood with one coefficient held
fixed can have several local maxima. Each limit is therefore checked by
maximizing afresh at it, from the estimate and from zero; a limit at
which a higher maximum turns up (after a second search that takes the
highest of the three starts at every point), or at which the
maximization does not converge, is `NA` with a warning, since the number
found would not be a limit of the interval.

## References

Heinze, G. and Schemper, M. (2002). A solution to the problem of
separation in logistic regression. *Statistics in Medicine*, 21,
2409-2419. [doi:10.1002/sim.1047](https://doi.org/10.1002/sim.1047)
