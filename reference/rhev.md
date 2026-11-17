# Simulate heteroskedastic extreme value choice data

Draws a data set from the HEV model of
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md):
utilities from `beta` (and `gamma`) as in
[`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md) plus
independent Gumbel errors with scales `scale` (one per level, in level
order; the base alternative's scale should be 1 to match the estimator's
normalization).

## Usage

``` r
rhev(n, beta, scale, gamma = NULL, x = NULL, levels = NULL)
```

## Arguments

- n:

  Number of cases.

- beta:

  A \\p \times (J-1)\\ matrix: rows are the intercept and the \\p - 1\\
  covariates, columns the non-base alternatives (column names become
  their labels; the base is `"a"` unless `levels` is given).

- scale:

  Gumbel scale parameters, one per level.

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

[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md)

## Examples

``` r
set.seed(1)
d <- rhev(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)), scale = c(1, 0.6, 1.5))
table(d$y)
#> 
#>   a   b   c 
#> 124 194 182 
```
