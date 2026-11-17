# Hosmer-Lemeshow-type goodness-of-fit test for multinomial models

The test of Fagerland, Hosmer and Bofin (2008): cases are sorted by the
fitted probability of the base alternative and split into `groups`
groups of equal size; the weighted observed and expected counts of every
alternative in every group give a Pearson statistic that is
approximately chi-square on \\(g - 2)(J - 1)\\ degrees of freedom under
a correct model. Small expected counts (below 5) are reported; the test
is a calibration check, not a test of IIA.

## Usage

``` r
hl_test(object, groups = 10L)
```

## Arguments

- object:

  A fitted choice model with independent single choices: not a
  [`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)
  or
  [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
  fit, and not a penalized fit (`firth = TRUE` or a prior). The cases
  are taken as independent; a fit's `cluster` does not enter the
  statistic.

- groups:

  Number of groups (default 10).

## Value

A list of class `"hl_test"` with `statistic`, `df`, `p.value`,
`observed` and `expected` (group by alternative), and `n_small`
(expected counts below 5).

## References

Fagerland, M.W., Hosmer, D.W. and Bofin, A.M. (2008). Multinomial
goodness-of-fit tests for logistic regression models. *Statistics in
Medicine*, 27, 4238-4253.
[doi:10.1002/sim.3202](https://doi.org/10.1002/sim.3202)

## See also

Other specification tests:
[`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md),
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
[`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md),
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)

## Examples

``` r
set.seed(1)
d <- rmnl(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
hl_test(mnl(y ~ x1 + x2, data = d))
#> Fagerland-Hosmer-Bofin goodness-of-fit test (10 groups by the base alternative's fitted probability)
#>   statistic = 10.332 on 16 df, p = 0.849
#>   1 expected count(s) below 5.
```
