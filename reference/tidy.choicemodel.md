# Tidy a choice model (broom method)

Tidy a choice model (broom method)

## Usage

``` r
# S3 method for class 'choicemodel'
tidy(x, conf.int = FALSE, conf.level = 0.95, ...)
```

## Arguments

- x:

  A `"choicemodel"` fit.

- conf.int:

  If `TRUE`, add `conf.low`/`conf.high` (see
  [`confint.choicemodel()`](https://bagozzib.github.io/polychoice/reference/confint.choicemodel.md)).

- conf.level:

  Confidence level for the interval.

- ...:

  Unused.

## Value

A data frame with one row per parameter: `term`, `estimate`,
`std.error`, `statistic`, `p.value`, and `component` (`"outcome"`,
`"generic"`, `"fixed effect"`, `"inflation"`, `"nest"`, `"scale"`,
`"covariance"`, or `"random intercept"`); inflation terms keep their
`infl_` prefix so both equations' terms stay distinct in tables. The
statistic and the intervals follow
[`summary.choicemodel()`](https://bagozzib.github.io/polychoice/reference/summary.choicemodel.md)
and
[`confint.choicemodel()`](https://bagozzib.github.io/polychoice/reference/confint.choicemodel.md)
(z against 1 for inclusive values and scales, intervals on the log scale
for positive structural parameters, profile penalized-likelihood
intervals for a penalized fit).

## Multiple imputation

[`mice::pool()`](https://amices.org/mice/reference/pool.html) combines a
list of fits, one per completed data set, through these methods: the
pooled estimates and standard errors are Rubin's rules. Three points
need care.

- Fit the model on each completed data set yourself:
  `with(imp, mnl(...))` does not work, because the estimators take the
  data as an argument.

- Give `pool()` the complete-data degrees of freedom,
  `dfcom = nobs(fit) - length(coef(fit))`. Left to itself, mice 3.19
  counts one residual per case and alternative, about J times too many,
  and its reference t distributions are then too narrow.

- The pooled statistic tests each parameter against 0 on the scale
  reported. For inclusive values, scales, and covariance parameters that
  is not the test
  [`summary.choicemodel()`](https://bagozzib.github.io/polychoice/reference/summary.choicemodel.md)
  reports. The fits' flags (`converged`, `boundary`, `singular`,
  `n_modes`) are not pooled either: check them on each fit first, as the
  example does.

## See also

[`summary.choicemodel()`](https://bagozzib.github.io/polychoice/reference/summary.choicemodel.md);
[`modelsummary::modelsummary()`](https://modelsummary.com/man/modelsummary.html)
consumes these methods, and
[`texreg::screenreg()`](https://rdrr.io/pkg/texreg/man/screenreg.html)
uses the package's `extract()` method; both carry the fit's
`converged`/`boundary`/`ill_conditioned` flags.

Other broom methods:
[`augment.choicemodel()`](https://bagozzib.github.io/polychoice/reference/augment.choicemodel.md),
[`glance.choicemodel()`](https://bagozzib.github.io/polychoice/reference/glance.choicemodel.md)

## Examples

``` r
if (requireNamespace("broom", quietly = TRUE)) {
  set.seed(1)
  d <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  m <- mnl(y ~ x1 + x2, data = d)
  broom::tidy(m, conf.int = TRUE)
  broom::glance(m)
  head(broom::augment(m))
}
#>   y         x1          x2 .fitted   .prob_a   .prob_b    .prob_c
#> 1 a -0.6264538  0.07730312       b 0.2371501 0.6758323 0.08701760
#> 2 c  0.1836433 -0.29686864       b 0.3247436 0.3626646 0.31259189
#> 3 c -0.8356286 -1.18324224       b 0.3254302 0.4548467 0.21972314
#> 4 c  1.5952808  0.01129269       c 0.2356573 0.1010445 0.66329826
#> 5 b  0.3295078  0.99160104       b 0.2542398 0.6223149 0.12344533
#> 6 b -0.8204684  1.59396745       b 0.0921925 0.8969239 0.01088358
# \donttest{
## multiple imputation: one fit per completed data set, pooled by mice
if (requireNamespace("mice", quietly = TRUE)) {
  set.seed(1)
  d <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
  d$x1[sample(500, 80)] <- NA
  imp <- mice::mice(d, m = 5, printFlag = FALSE, seed = 1)
  fits <- lapply(1:5, function(i) mnl(y ~ x1 + x2, data = mice::complete(imp, i)))
  all(vapply(fits, function(f) f$converged && !f$boundary && !f$singular, logical(1)))
  summary(mice::pool(fits, dfcom = nobs(fits[[1]]) - length(coef(fits[[1]]))))
}
#>            term   estimate std.error statistic        df      p.value
#> 1 b:(Intercept)  0.4821263 0.1232712  3.911103 475.78603 1.052452e-04
#> 2          b:x1 -0.7073566 0.1777954 -3.978486  31.11162 3.859728e-04
#> 3          b:x2  0.6748004 0.1309133  5.154561 441.90770 3.841501e-07
#> 4 c:(Intercept) -0.3711618 0.1679398 -2.210088 253.25992 2.799368e-02
#> 5          c:x1  0.9470919 0.1783153  5.311331 119.41062 5.113435e-07
#> 6          c:x2 -0.6314056 0.1471864 -4.289837 330.41730 2.350910e-05
# }
```
