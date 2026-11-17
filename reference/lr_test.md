# Likelihood-ratio test for nested choice models

For a restricted model that is a special case of the full one: the same
estimator (both
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) or
both [`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md)
with the same inflated alternative and link) on the same cases, with the
restricted parameters a subset of the full ones (checked by name; two
nested logits must share their nests, and a common inclusive value,
`equal_lambda = TRUE`, nests the free ones); or an
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) fit
against the
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md), or
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)
fit with the same covariates, which tests the inclusive-value parameters
against 1 (the nested-logit test of IIA), the scale parameters against 1
(homoskedastic errors), or the random-intercept standard deviations
against 0.

## Usage

``` r
lr_test(restricted, full)
```

## Arguments

- restricted, full:

  Two `"choicemodel"` fits on the same cases.

## Value

A data frame of class `"lr_test"` with the statistic, degrees of
freedom, and p-value; the attribute `"suprema"` is `TRUE` when a
supremum replaced a separated fit's value. For the constant-share
inflated multinomial logit, `df` is q + 1 (the share and the q further
restrictions), the attribute `"mixture"` holds q, and `p.value` is `NA`,
with a warning, when the inflated fit is degenerate (outcome
coefficients diverging, or a saddle).

## Details

Separation. When the data of an
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) or
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
fit (without unit fixed effects or inflation) are separated, the
maximum-likelihood estimates do not exist and the fit stops somewhere
along the separating direction; the statistic then compares the suprema
of the two log-likelihoods, which
[`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md)
makes computable: 0 under complete separation, and under quasi-complete
separation the value reached by moving along the separating direction
and re-maximizing. Likelihood-ratio tests stay valid under separation
where Wald tests do not (Rainey 2024), but with few cases or rankings
per coefficient the chi-square reference can be far off;
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
calibrates it.
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)
and
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md)
use the same suprema.

Penalized fits. Two penalized fits (`firth = TRUE`, or the same prior)
are compared by the penalized likelihood ratio: the restricted maximum
is taken over the full model's penalized log-likelihood (Heinze and
Schemper 2001). The priors count as the same when they share the family,
degrees of freedom, autoscaling and constants' scale, and the scales of
the coefficients both fits have. A scale named for a coefficient only
the full model has is ignored by the restricted fit, with a message. A
penalized fit against an unpenalized one, or two different penalties,
are refused.

Boundaries. Two of these restrictions sit on the boundary of the
parameter space and the chi-square reference does not hold exactly. For
the random-intercept model (\\\sigma = 0\\) the null distribution is a
mixture of chi-squares with fewer degrees of freedom, so the reported
p-value is conservative and the test can be used as such (a message says
so). For the inflated model against the multinomial logit (\\\phi = 1\\)
the split-equation slopes are not identified under the null, so no
simple correction exists: the function refuses that pair;
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
calibrates it by simulation (the normal reference of
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md) is
not justified there either). The exception is an inflated model with one
inflated alternative and an intercept-only split (`y ~ x | 1`, the dogit
with a constant captive share): its captive odds at zero are a regular
boundary point, and the reference is the 50:50 mixture of chi-square(q)
and chi-square(q + 1), q the further restrictions (Self and Liang 1987),
as in `captivity_test(type = "lr")`. The multinomial probit does not
nest the multinomial logit: use
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md).
Two fits of the same family other than the multinomial logit (two HEV,
nested-logit, multinomial probit, or random-intercept fits) are accepted
when the restricted parameters are a subset of the full ones and the
fits use the same quadrature nodes or simulation draws. The function
also refuses fits with different choice sets, and a fixed-effects fit
against one without fixed effects; two fixed-effects fits are compared
with a message when the panels are short, since the
incidental-parameters bias shifts the statistic.

Clustered data. The statistic assumes independent cases: a fit's
`cluster` does not enter it.

## References

Rainey, C. (2024). Hypothesis tests under separation. *Political
Analysis*, 32, 172-185.
[doi:10.1017/pan.2023.28](https://doi.org/10.1017/pan.2023.28)

Self, S.G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*, 82,
605-610.
[doi:10.1080/01621459.1987.10478472](https://doi.org/10.1080/01621459.1987.10478472)

## See also

[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md),
[`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md)

Other model comparison:
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md),
[`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md),
[`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md),
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)

## Examples

``` r
set.seed(1)
d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
lr_test(mnl(y ~ x1, data = d), mnl(y ~ x1 + x2, data = d))
#> Likelihood-ratio test
#>   LR = 135.9978 on 2 df, p = < 2.2e-16
```
