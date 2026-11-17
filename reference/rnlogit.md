# Simulate nested-logit choice data

Draws a data set from the nested logit of
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md):
utilities from `beta` (and `gamma`) as in
[`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md),
nests as in
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
and inclusive-value parameters `lambda` (one per nest, in the order of
`nests`; a singleton nest's value is irrelevant). With case-specific
covariates only, `lambda` is identified just through the curvature of
the inclusive values and recovered imprecisely; alternative-specific
covariates (`gamma`) identify it sharply.

## Usage

``` r
rnlogit(n, beta, nests, lambda, gamma = NULL, x = NULL, levels = NULL)
```

## Arguments

- n:

  Number of cases.

- beta:

  A \\p \times (J-1)\\ matrix: rows are the intercept and the \\p - 1\\
  covariates, columns the non-base alternatives (column names become
  their labels; the base is `"a"` unless `levels` is given).

- nests:

  A named list assigning every level to a nest.

- lambda:

  Inclusive-value parameters, one per nest.

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

[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md)

## Examples

``` r
set.seed(1)
d <- rnlogit(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
             nests = list(n1 = c("a", "b"), n2 = "c"), lambda = c(0.5, 1))
table(d$y)
#> 
#>   a   b   c 
#>  92 245 163 
```
