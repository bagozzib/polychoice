# Fit statistics for a choice model

Log-likelihood, information criteria, and the pseudo-R-squared measures
in common use for multinomial models, all relative to the same null: the
model with alternative-specific constants only (equal choice
probabilities over the available alternatives when the model has no
constants). Inflated and uninflated fits share this null, so their
McFadden measures are comparable.

## Usage

``` r
fit_stats(object)
```

## Arguments

- object:

  A `"choicemodel"` fit.

## Value

A one-row data frame of class `"fit_stats"` with `loglik`,
`loglik_null`, `df`, `n` (cases), `AIC`, `BIC` (with the weight total as
\\n\\), `mcfadden` (\\1 - \ell/\ell_0\\), `mcfadden_adj` (\\1 - (\ell -
K)/\ell_0\\), `cox_snell`, `nagelkerke`, `count` (the weighted share of
cases whose modal predicted alternative is the chosen one), and
`count_adj` (the proportional reduction in error relative to always
predicting the most frequent alternative, Long 1997).

## Details

The null fit keeps the fit's availability pattern, including the
alternatives that unit fixed effects switched off for never-chosen
pairs. For a random-intercept fit `n` is the number of cases (not
units), and so is the \\n\\ of the BIC. For a penalized fit
(`firth = TRUE` or a prior) the log-likelihood is the ordinary one at
the penalized estimates. The information criteria and R-squared measures
built on it are therefore not maximum-likelihood ones (the print says
so). The null is the maximum-likelihood constants-only fit, and
`loglik_penalized` holds the maximized objective.

## References

Long, J.S. (1997). *Regression Models for Categorical and Limited
Dependent Variables*. Sage.

## See also

[`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md),
[`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)

Other model comparison:
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md),
[`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md),
[`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)

## Examples

``` r
set.seed(1)
d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
fit_stats(mnl(y ~ x1 + x2, data = d))
#> Log-likelihood -682.248 (null -835.167), 6 parameters, 800 cases; AIC 1376.497, BIC 1404.605
#> McFadden R2 0.183 (adjusted 0.176); Cox-Snell 0.318; Nagelkerke 0.363
#> Count R2 0.618 (adjusted: 0.248 proportional reduction in error over the modal alternative)
```
