# Joint tests of each term across equations

For every term of the model, tests that its coefficients are zero in all
equations at once: the Wald form from the fitted covariance, or the
likelihood-ratio form by refitting without the term (Stata's
`mlogtest, wald` and `lr`). For a generic coefficient on an
alternative-specific covariate the test has one degree of freedom. The
terms are those of the outcome equations; for an inflated model the
split-equation coefficients are not included (test them by
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
between two
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md) fits
that differ in the split equation, or from the `infl_` rows of
[`summary.choicemodel()`](https://bagozzib.github.io/polychoice/reference/summary.choicemodel.md)).

## Usage

``` r
term_test(object, type = c("wald", "lr"))
```

## Arguments

- object:

  A fitted choice model.

- type:

  `"wald"` (default) or `"lr"`.

## Value

A data frame of class `"iia_test"` with one row per term: `term`,
`statistic`, `df`, `p.value`.

## Details

The likelihood-ratio form works as follows:

- The model without the term is fitted as the estimators fit a model,
  multistart and saddle step included, so that both log-likelihoods are
  maxima; a warning says when that fit did not converge or still stops
  at a saddle.

- For a maximum-likelihood fit of separated data
  ([`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
  [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
  without unit fixed effects or inflation), it compares the suprema of
  the log-likelihoods, not where the optimizer stopped (see
  [`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)).
  The method line says so.

- For a penalized fit (`firth = TRUE` or a prior), it is the penalized
  likelihood ratio: the term's coefficients are set to zero in the full
  model's penalized log-likelihood.

The Wald form of a penalized fit uses the penalized estimates and their
covariance; for a separated coefficient it can be far off (see
[`summary.choicemodel()`](https://bagozzib.github.io/polychoice/reference/summary.choicemodel.md)).

The likelihood-ratio form assumes independent cases (a fit's `cluster`
does not enter it); the Wald form uses the fit's covariance,
cluster-robust when the fit has one.

## See also

[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)

Other specification tests:
[`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md),
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
[`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md),
[`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md)

## Examples

``` r
set.seed(1)
d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
term_test(mnl(y ~ x1 + x2, data = d))
#> Joint tests of each term across equations (Wald) 
#>  term statistic df p.value
#>    x1   145.392  2  <2e-16
#>    x2   107.417  2  <2e-16
```
