# Relative risk ratios

Exponentiated coefficients of a choice model: for an
alternative-specific coefficient \\\beta\_{jm}\\, \\e^{\beta\_{jm}}\\ is
the factor by which a unit increase in covariate \\m\\ multiplies the
odds of alternative \\j\\ relative to the base (Stata's `mlogit, rrr`);
for a generic coefficient it multiplies the odds of any alternative
relative to any other per unit of the alternative-specific covariate;
for the split equation it multiplies the odds of the choice regime
(logit link only).

## Usage

``` r
rrr(object, level = 0.95)
```

## Arguments

- object:

  A `"choicemodel"` fit.

- level:

  Confidence level for the (Wald, exponentiated) interval.

## Value

A data frame with `term`, `block`, `rrr`, `lower`, `upper`, and the z
statistic and p-value of the underlying coefficient.

## Details

The odds reading holds for the multinomial logit, the inflated
multinomial logit (whose outcome stage is a multinomial logit within the
choice regime), and the random-intercept multinomial logit, where it is
conditional on the unit's intercepts (the population-averaged odds ratio
is closer to one). It does not hold for the nested logit (within a nest
the ratio is \\e^{\beta/\lambda}\\, across nests it is not constant),
the HEV model (the odds depend on the scales), or the multinomial
probit, so `rrr()` refuses those fits: use
[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md).

## See also

[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md),
[`summary.choicemodel()`](https://bagozzib.github.io/polychoice/reference/summary.choicemodel.md)

Other quantities of interest:
[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md)

## Examples

``` r
set.seed(1)
d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
rrr(mnl(y ~ x1 + x2, data = d))
#>           term   block   rrr lower upper      z  p.value
#>  b:(Intercept) outcome 1.752 1.443 2.127  5.671 1.42e-08
#>           b:x1 outcome 0.374 0.299 0.467 -8.649  < 2e-16
#>           b:x2 outcome 2.030 1.653 2.492  6.763 1.35e-11
#>  c:(Intercept) outcome 0.840 0.664 1.064 -1.447    0.148
#>           c:x1 outcome 1.637 1.303 2.056  4.234 2.29e-05
#>           c:x2 outcome 0.632 0.511 0.781 -4.241 2.23e-05
```
