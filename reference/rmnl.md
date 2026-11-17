# Simulate multinomial-logit (optionally inflated) choice data

Draws a data set from the model of
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) /
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md):
standard normal case-specific covariates `x1`, `x2`, ... (or those
supplied in `x`), alternative-specific coefficients `beta` with the
first level as the base, optionally alternative-specific covariates
`w1`, `w2`, ... with generic coefficients `gamma`, and, when `delta` is
given, a split equation on `z1`, `z2`, ... (or `z`) that puts cases in
the captive regime, where they always choose the inflated alternative.

## Usage

``` r
rmnl(
  n,
  beta,
  delta = NULL,
  inflate = NULL,
  link = c("logit", "probit"),
  gamma = NULL,
  x = NULL,
  z = NULL,
  levels = NULL
)
```

## Arguments

- n:

  Number of cases.

- beta:

  A \\p \times (J-1)\\ matrix: rows are the intercept and the \\p - 1\\
  covariates, columns the non-base alternatives (column names become
  their labels; the base is `"a"` unless `levels` is given).

- delta:

  Optional split-equation coefficients (intercept first) on `z`; `NULL`
  for a plain MNL. With several inflated alternatives, a matrix with one
  column per inflated alternative: the coefficients of \\\log(\pi\_{i0}
  / \pi\_{ic})\\, the log-odds of the choice regime against captivity to
  that alternative (see
  [`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md)).

- inflate:

  The inflated alternative (a label); default the base. Several labels
  draw several captive alternatives (logit split only); the output then
  also has the column `regime` (`"choice"` or the captive alternative).

- link:

  Link of the split equation.

- gamma:

  Optional generic coefficients on alternative-specific covariates `w1`,
  `w2`, ... (standard normal draws for every case and alternative). When
  given, the data are returned in the long layout: one row per case and
  alternative with columns `id`, `alt`, `chosen` (0/1), the `w`'s, the
  `x`'s (and, for inflated draws of `rmnl()`, the `z`'s and `s`), ready
  for `chosen ~ w1 | x1 + x2` with `id = "id", alt = "alt"`.

- x:

  Optional data frame of case-specific covariates (`n` rows) to use
  instead of standard normal draws.

- z:

  Optional data frame of split-equation covariates (`n` rows) to use
  instead of standard normal draws.

- levels:

  Optional labels of all \\J\\ alternatives (base first).

## Value

A data frame: in the wide layout `y` (factor) and the covariates, plus
the latent regime `s` (1 = choice regime) for inflated draws; in the
long layout as described under `gamma`. The true parameters are stored
in the attribute `"truth"`.

## See also

[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
[`simulate.choicemodel()`](https://bagozzib.github.io/polychoice/reference/simulate.choicemodel.md)

## Examples

``` r
set.seed(1)
d <- rmnl(1000, beta = cbind(b = c(0.2, -1, 0.5), c = c(-0.5, 0.8, -0.3)))
table(d$y)
#> 
#>   a   b   c 
#> 286 445 269 
di <- rmnl(1000, beta = cbind(b = c(0.2, -1, 0.5), c = c(-0.5, 0.8, -0.3)),
           delta = c(0.5, 1), inflate = "a")
table(di$y, di$s)
#>    
#>       0   1
#>   a 390 175
#>   b   0 270
#>   c   0 165
dl <- rmnl(500, beta = cbind(b = c(0.2, -1, 0.5), c = c(-0.5, 0.8, -0.3)), gamma = 0.8)
head(dl)
#>   id alt chosen          w1         x1         x2
#> 1  1   a      0 -0.79480342  0.8871888 -1.3254177
#> 2  1   b      0  0.26370340  0.8871888 -1.3254177
#> 3  1   c      1  0.99010104  0.8871888 -1.3254177
#> 4  2   a      1  0.60608460 -0.3776280  0.9519797
#> 5  2   b      0 -0.82945185 -0.3776280  0.9519797
#> 6  2   c      0 -0.06643542 -0.3776280  0.9519797
```
