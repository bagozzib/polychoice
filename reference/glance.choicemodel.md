# Glance at a choice model (broom method)

Glance at a choice model (broom method)

## Usage

``` r
# S3 method for class 'choicemodel'
glance(x, ...)
```

## Arguments

- x:

  A `"choicemodel"` fit.

- ...:

  Unused.

## Value

A one-row data frame: `logLik`, `AIC`, `BIC`, `df`, `nobs`, for inflated
models `inflate` and `share_inflated`, the fit-quality flags
`converged`, `boundary`, `ill_conditioned`, and `n_modes` (see the
section "Fit-quality flags" of
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)), for
fixed-effects fits `bias_correct` and `jackknife_refused` (`TRUE` when a
requested jackknife correction was refused), and for penalized fits
(`firth = TRUE` or a prior) `penalized = TRUE`: their `logLik` is the
ordinary log-likelihood at the penalized estimates, and `AIC` and `BIC`
are not maximum-likelihood criteria.

## See also

Other broom methods:
[`augment.choicemodel()`](https://bagozzib.github.io/polychoice/reference/augment.choicemodel.md),
[`tidy.choicemodel()`](https://bagozzib.github.io/polychoice/reference/tidy.choicemodel.md)
