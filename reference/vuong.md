# Vuong test for non-nested or overlapping choice models

Vuong's (1989) test comparing two fits to the same cases; positive
values favour the first model. It is the comparison the inflation
literature uses between an inflated model and its uninflated counterpart
(Bagozzi and Marchetti 2017; Greene 1994), where the restriction \\\phi
= 1\\ sits at the boundary so that the likelihood-ratio statistic is not
chi-square. Its standard normal reference does not hold there: when the
uninflated model is true the two models are nested at the boundary, the
per-case differences degenerate, and the statistic is not N(0, 1)
(Wilson 2015, on the same use for zero-inflated counts; a size study
under the MNL null is part of the package's verification,
`verification/vuong_size_mc.R` in the development repository). The
calibrated alternative is the parametric bootstrap of the
likelihood-ratio and Vuong statistics under the restricted model,
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md).
Between strictly non-nested models the normal reference holds; for
overlapping ones (sharing some covariates) Vuong's procedure first tests
whether the per-case differences have zero variance (the models are then
equivalent on these data), which this function reports only as an error
when the differences are identical. With frequency weights the statistic
equals the one computed on duplicated cases. Reported with the raw
numerator and with the AIC and BIC corrections for the difference in
parameter counts. When the two fits are practically the same model on
the data (one at the other's boundary, such as an inflated fit with
\\\phi\\ at 1), the per-case differences are round-off, every statistic
divides by their spread, and the corrected ones run to minus millions:
when that spread totals less than 0.01 (\\\sqrt{n}\\\hat\omega \<
0.01\\) the statistics are `NA` and a note says the models cannot be
told apart. When either model is a
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)
fit, whose likelihood is per unit, the contributions of both are summed
by unit and the statistic (and the \\n\\ of its BIC correction) is over
units, while
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md)
counts cases.

## Usage

``` r
vuong(m1, m2)
```

## Arguments

- m1, m2:

  Two `"choicemodel"` fits to the same cases.

## Value

A data frame of class `"vuong"` with one row per correction:
`statistic`, `p_m1_better` (one-sided, m1 preferred), `p_m2_better`, and
`p_two_sided`.

## Details

The statistics assume independent cases: a fit's `cluster` does not
enter them.

## References

Bagozzi, B.E. and Marchetti, K. (2017). Distinguishing occasional
abstention from routine indifference in models of vote choice.
*Political Science Research and Methods*, 5(2), 277-294.
[doi:10.1017/psrm.2015.42](https://doi.org/10.1017/psrm.2015.42)

Greene, W.H. (1994). Accounting for excess zeros and sample selection in
Poisson and negative binomial regression models. Working Paper EC-94-10,
Stern School of Business, New York University.

Vuong, Q.H. (1989). Likelihood ratio tests for model selection and
non-nested hypotheses. *Econometrica*, 57, 307-333.
[doi:10.2307/1912557](https://doi.org/10.2307/1912557)

Wilson, P. (2015). The misuse of the Vuong test for non-nested models to
test for zero-inflation. *Economics Letters*, 127, 51-53.
[doi:10.1016/j.econlet.2014.12.029](https://doi.org/10.1016/j.econlet.2014.12.029)

## See also

[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
[`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md)

Other model comparison:
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md),
[`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md),
[`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md),
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)

## Examples

``` r
set.seed(2)
d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
          delta = c(0.5, 1.2), inflate = "a")
vuong(imnl(y ~ x1 + x2 | z1, data = d, inflate = "a"), mnl(y ~ x1 + x2, data = d))
#> Vuong test:  m1 = imnl(y ~ x1 + x2 | z1, data = d, inflate = "a") [ Inflated multinomial logit (inflated alternative: a, logit split) ]
#>              m2 = mnl(y ~ x1 + x2, data = d) [ Multinomial logit ]
#> positive statistics favor m1
#> 
#>  correction statistic p_m1_better p_m2_better p_two_sided
#>         raw     5.492           0           1           0
#>         AIC     5.294           0           1           0
#>         BIC     4.830           0           1           0
#> Note: The inflated multinomial logit and the multinomial logit are nested at the boundary phi = 1, where the normal reference does not hold for testing the uninflated model; boot_lr_test() calibrates the statistic. 
```
