# Simulate multinomial probit choice data

Draws a data set from the multinomial probit of
[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md):
utilities from `beta` (and `gamma`) as in
[`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md) plus
multivariate normal errors whose differences from the base alternative
have covariance `sigma` (\\(J-1) \times (J-1)\\, first element 1 to
match the estimator's normalization). With case-specific covariates only
the covariance is fragilely identified (Keane 1992);
alternative-specific covariates (`gamma`) identify it.

## Usage

``` r
rmnp(n, beta, sigma, gamma = NULL, x = NULL, levels = NULL)
```

## Arguments

- n:

  Number of cases.

- beta:

  A \\p \times (J-1)\\ matrix: rows are the intercept and the \\p - 1\\
  covariates, columns the non-base alternatives (column names become
  their labels; the base is `"a"` unless `levels` is given).

- sigma:

  Covariance of the error differences from the base, a \\(J-1) \times
  (J-1)\\ matrix.

- gamma:

  Optional generic coefficients on alternative-specific covariates `w1`,
  `w2`, ... (standard normal draws for every case and alternative). When
  given, the data are returned in the long layout: one row per case and
  alternative with columns `id`, `alt`, `chosen` (0/1), the `w`'s, the
  `x`'s (and, for inflated draws of
  [`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md),
  the `z`'s and `s`), ready for `chosen ~ w1 | x1 + x2` with
  `id = "id", alt = "alt"`.

- x:

  Optional data frame of case-specific covariates (`n` rows) to use
  instead of standard normal draws.

- levels:

  Optional labels of all \\J\\ alternatives (base first).

## Value

A data frame as in
[`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md); the
truth in the attribute `"truth"`.

## See also

[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md),
[`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md)

## Examples

``` r
set.seed(1)
d <- rmnp(500, beta = cbind(b = c(0.3, -0.6, 0.5), c = c(-0.2, 0.4, -0.3)),
          sigma = matrix(c(1, 0.5, 0.5, 1), 2), gamma = 1)
head(d)
#>   id alt chosen         w1         x1          x2
#> 1  1   a      1  1.1349651 -0.6264538  0.07730312
#> 2  1   b      0  0.8500435 -0.6264538  0.07730312
#> 3  1   c      0 -0.8861496 -0.6264538  0.07730312
#> 4  2   a      1  1.1119318  0.1836433 -0.29686864
#> 5  2   b      0 -0.9253130  0.1836433 -0.29686864
#> 6  2   c      0 -1.9222549  0.1836433 -0.29686864
```
