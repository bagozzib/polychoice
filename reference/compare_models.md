# Compare several choice models side by side

Compare several choice models side by side

## Usage

``` r
compare_models(...)
```

## Arguments

- ...:

  Named or unnamed `"choicemodel"` fits.

## Value

A data frame of class `"compare_models"` with one row per model: the
model type, inflated alternative, cases, parameters, log-likelihood,
AIC, BIC, McFadden R-squared, and the `converged`, `boundary`, and
`ill_conditioned` flags and `n_modes`; with a penalized fit among them a
`penalized` column, and the print notes that their log-likelihoods and
criteria are at the penalized estimates, not maximum-likelihood ones.

## See also

[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md)

Other model comparison:
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md),
[`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md),
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)

## Examples

``` r
set.seed(2)
d <- rmnl(1000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
          delta = c(0.5, 1.2), inflate = "a")
compare_models(mnl = mnl(y ~ x1 + x2, data = d),
               imnl = imnl(y ~ x1 + x2 | z1, data = d, inflate = "a"))
#>  model                                                              type
#>    mnl                                                 Multinomial logit
#>   imnl Inflated multinomial logit (inflated alternative: a, logit split)
#>  inflate    n df   loglik      AIC      BIC mcfadden converged boundary
#>          1000  6 -857.498 1726.996 1756.443    0.117      TRUE    FALSE
#>        a 1000  8 -792.033 1600.065 1639.327    0.185      TRUE    FALSE
#>  ill_conditioned n_modes
#>            FALSE       1
#>            FALSE       1
```
