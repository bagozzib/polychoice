# Simulate rank-ordered choice data

Draws rankings from the rank-ordered logit of
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md):
utilities from `beta` (and `gamma`) as in
[`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md) plus
independent standard Gumbel errors; each case ranks its alternatives by
utility. `depth` keeps only the first positions (the rest unranked, rank
0), and `ties` merges each ranked position into the one above it with
that probability, independently of the order within the merged group,
which is the coarsening under which the exact tie likelihood of
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
is the likelihood of the reported ranking.

## Usage

``` r
rrologit(
  n,
  beta,
  gamma = NULL,
  x = NULL,
  levels = NULL,
  depth = NULL,
  ties = 0
)
```

## Arguments

- n:

  Number of cases.

- beta:

  A \\p \times (J-1)\\ matrix: rows are the intercept and the \\p - 1\\
  covariates, columns the non-base alternatives (column names become
  their labels; the base is `"a"` unless `levels` is given).

- gamma:

  Optional generic coefficients on alternative-specific covariates `w1`,
  `w2`, ... (standard normal draws for every case and alternative). The
  data are in the long layout with or without them (see Value).

- x:

  Optional data frame of case-specific covariates (`n` rows) to use
  instead of standard normal draws.

- levels:

  Optional labels of all \\J\\ alternatives (base first).

- depth:

  Optional number of ranked positions (default: every alternative is
  ranked).

- ties:

  Probability that a ranked position is tied with the position above it
  (default 0).

## Value

A data frame in the long layout, one row per case and alternative: `id`,
`alt`, `rank` (1 = most preferred; equal values are ties; 0 = unranked),
the `w`'s and the `x`'s, ready for `rank ~ w1 | x1 + x2` with
`id = "id", alt = "alt"`. The true parameters are in the attribute
`"truth"`.

## See also

[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md),
[`rmnl()`](https://bagozzib.github.io/polychoice/reference/rmnl.md)

## Examples

``` r
set.seed(1)
d <- rrologit(300, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6), d = c(0.2, 0.4)), depth = 2)
head(d, 8)
#>   id alt rank         x1
#> 1  1   a    2 -0.6264538
#> 2  1   b    1 -0.6264538
#> 3  1   c    0 -0.6264538
#> 4  1   d    0 -0.6264538
#> 5  2   a    1  0.1836433
#> 6  2   b    0  0.1836433
#> 7  2   c    0  0.1836433
#> 8  2   d    2  0.1836433
```
