# The battery of IIA tests

Runs
[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md)
(both covariance forms),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
and
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md)
on one fitted
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) model
and collects the results in one table. See the individual tests for
definitions and caveats; in particular, with case-specific covariates
only, every test draws its power from the covariates' nonlinear
information, and the tests are not powerful against every departure from
IIA.

## Usage

``` r
iia_tests(
  object,
  drop = NULL,
  splits = 10L,
  seed = NULL,
  pairs = c("all", "base")
)
```

## Arguments

- object:

  A fitted
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  model with at least three alternatives. A fit with diverging
  coefficients (separation) draws a warning: the contrasts then compare
  estimates that do not exist and can reject spuriously. The Fok-Paap
  \\Q_1\\ row leaves the diverging directions out of its comparison,
  with its degrees of freedom reduced, and the printed table says so
  (see
  [`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md)).

- drop:

  Alternatives for the Hausman-McFadden and Small-Hsiao tests; default
  all.

- splits:

  Number of Small-Hsiao half-splits (default 10; the battery reports the
  median statistic over splits with the chi-square p-value of that
  median, a descriptive figure).

- seed:

  Optional seed for the Small-Hsiao splits.

- pairs:

  Pairs for the Fok-Paap tests (`"all"` or `"base"`).

## Value

A data frame of class `"iia_battery"` with `test`, `alternative`
(dropped; blank for the pairwise tests), `statistic`, `df`, `p.value`,
`psd` (whether the Hausman-McFadden covariance difference is positive
semidefinite; `NA` where it does not apply), and `share_separated` (for
the Small-Hsiao rows, the share of half-sample splits with a separated
fit; `NA` otherwise).

## See also

[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md)

Other IIA tests:
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md)

## Examples

``` r
set.seed(1)
d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
iia_tests(mnl(y ~ x1 + x2, data = d), seed = 1)
#> Tests of IIA / misspecification of the multinomial logit 
#>                                            test alternative statistic df
#>                   Hausman-McFadden (difference)           a     3.749  6
#>                   Hausman-McFadden (difference)           b     4.574  6
#>                   Hausman-McFadden (difference)           c    10.549  6
#>                   Hausman-McFadden (difference)           d     4.144  6
#>                        Hausman-McFadden (joint)           a     6.888  6
#>                        Hausman-McFadden (joint)           b     6.262  6
#>                        Hausman-McFadden (joint)           c     5.407  6
#>                        Hausman-McFadden (joint)           d     5.755  6
#>                Small-Hsiao (median over splits)           a     3.266  6
#>                Small-Hsiao (median over splits)           b     2.823  6
#>                Small-Hsiao (median over splits)           c     2.810  6
#>                Small-Hsiao (median over splits)           d     3.193  6
#>                 Fok-Paap CL Hausman (Q1, joint)                 8.593  9
#>  Fok-Paap GMM overidentification (Q2, two-step)                18.101  9
#>  p.value   psd share_separated
#>    0.711  TRUE              NA
#>    0.600  TRUE              NA
#>    0.103 FALSE              NA
#>    0.657 FALSE              NA
#>    0.331  TRUE              NA
#>    0.394  TRUE              NA
#>    0.493  TRUE              NA
#>    0.451  TRUE              NA
#>    0.775    NA               0
#>    0.831    NA               0
#>    0.832    NA               0
#>    0.784    NA               0
#>    0.476  TRUE              NA
#>    0.034    NA              NA
#> psd = FALSE: the covariance difference was not positive semi-definite for that row. A negative statistic is no evidence against the model (p-value 1); a positive one from a variance-difference form is unreliable (see the test's help).
#> Small-Hsiao (median over splits): descriptive, not a test at its level (it rejected in none of 300 samples in the package's size study); for the test, run one split with its seed fixed in advance (small_hsiao(m, seed = ), or iia_tests(m, splits = 1, seed = )).
```
