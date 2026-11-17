# Nested-logit score test of IIA

The Lagrange-multiplier (score) test of the multinomial logit against
the nested logit with the given nests, evaluated at the
multinomial-logit estimates where every inclusive-value parameter equals
1 (McFadden 1987; Hausman and McFadden 1984, section 4). Only the
restricted (MNL) fit is needed: the score of the nested-logit
log-likelihood with respect to the inclusive-value parameters at
\\\lambda = 1\\ is combined with the information matrix (by default its
outer-product form, McFadden's regression form, positive semidefinite by
construction; or the differenced Hessian with
`information = "hessian"`). The test is directional, towards the stated
nesting structure, and the sign of each nest's score shows whether the
data pull its \\\lambda\\ below 1 (correlated errors within the nest) or
above.

## Usage

``` r
nested_lm_test(
  object,
  nests,
  equal_lambda = FALSE,
  information = c("opg", "hessian")
)
```

## Arguments

- object:

  A fitted
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  model. A fit with diverging coefficients (separation) draws a warning:
  the statistic is then evaluated where the optimizer stopped, at
  estimates that do not exist.

- nests:

  A named list assigning every alternative to a nest, as in
  [`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md).

- equal_lambda:

  If `TRUE`, one shared inclusive-value parameter (one degree of
  freedom).

- information:

  `"opg"` (default; the outer product of the scores, the form of
  McFadden's regression-based test, positive semidefinite by
  construction) or `"hessian"` (the observed information at the
  restricted estimate, which is not guaranteed positive definite there;
  a negative statistic is reported with a warning and no p-value).

## Value

A data frame of class `"iia_test"` with the statistic, degrees of
freedom, and p-value; the attribute `"scores"` holds the per-nest
scores.

## Details

The statistic assumes independent cases: a fit's `cluster` does not
enter it.

## References

Hausman, J.A. and McFadden, D. (1984). Specification tests for the
multinomial logit model. *Econometrica*, 52, 1219-1240.
[doi:10.2307/1910997](https://doi.org/10.2307/1910997)

McFadden, D. (1987). Regression-based specification tests for the
multinomial logit model. *Journal of Econometrics*, 34, 63-82.
[doi:10.1016/0304-4076(87)90067-4](https://doi.org/10.1016/0304-4076%2887%2990067-4)

## See also

[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md),
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)

Other IIA tests:
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md)

## Examples

``` r
set.seed(2)
d <- rnlogit(2000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)),
             nests = list(n1 = c("a", "b"), n2 = c("c", "d")), lambda = c(0.5, 0.5))
m <- mnl(y ~ x1 + x2, data = d)
nested_lm_test(m, nests = list(n1 = c("a", "b"), n2 = c("c", "d")))
#> Lagrange-multiplier test of IIA against a nested logit 
#>                        test statistic df p.value
#>  Nested-logit LM test (opg)    19.441  2   6e-05
```
