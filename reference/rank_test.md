# Stability of preferences across ranking stages

Tests the assumption of the rank-ordered logit that the relative
preference for two alternatives is the same at every stage of the
ranking: the model is refit with a second set of coefficients for the
stages after position `after` (every coefficient, including the
constants), and the hypothesis that the two sets are equal is tested by
the likelihood ratio, as Hausman and Ruud (1987, section 3) do, and by
the Wald statistic, as Allison and Christakis (1994, section 5.4) do;
both are chi-square with as many degrees of freedom as the model has
parameters. A coefficient that one block of stages does not identify (an
alternative that no ranking orders there, such as one always ranked
within the first `after` positions) is left out of both statistics and
their degrees of freedom, with a warning; a penalized fit
(`firth = TRUE` or a prior) is then refused, since its fit with separate
early and late coefficients has no finite penalty. Chapman and Staelin's
(1982) pooling test is the sequential version that compares the first
\\E\\ stages with stage \\E + 1\\. Within a tie group the stages take
the positions of each order in the sum over orders.

## Usage

``` r
rank_test(object, after = 1L)
```

## Arguments

- object:

  A
  [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
  fit.

- after:

  The last position of the early stages (default 1: the first choice
  against the rest of the ranking).

## Value

A data frame of class `"iia_test"` with the LR and Wald rows
(`statistic`, `df`, `p.value`); the attribute `"coefficients"` holds the
early and late estimates and their difference with its standard error.
The Wald statistic uses the fit's covariance type (analytic for a
bootstrap fit).

## Details

A rejection says that later rankings follow different coefficients. That
can mean careless ranking further down (the coefficients then shrink
towards zero in the late stages, the reading of Chapman and Staelin and
of Hausman and Ruud) or a change in valuation; Allison and Christakis
also note that apparent differences between stages can arise from
heterogeneity across respondents rather than from a change within
respondents. Fitting with `depth = after` uses only the early stages.

The likelihood-ratio statistic assumes independent rankings (a fit's
`cluster` does not enter it); the Wald statistic uses the fit's
covariance, cluster-robust when the fit has one. Their chi-square
reference is asymptotic; with few rankings per coefficient it can be far
off, and no parametric bootstrap of this test is offered. The split
design (separate coefficients after the first stage) is outside
[`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md)'s
scope: when a split fit's data are separated, its maximum-likelihood
statistic uses the value where its optimizer stopped, not the supremum
(the penalized form is unaffected).

## References

Allison, P.D. and Christakis, N.A. (1994). Logit models for sets of
ranked items. *Sociological Methodology*, 24, 199-228.
[doi:10.2307/270983](https://doi.org/10.2307/270983)

Chapman, R.G. and Staelin, R. (1982). Exploiting rank ordered choice set
data within the stochastic utility model. *Journal of Marketing
Research*, 19, 288-301.
[doi:10.1177/002224378201900302](https://doi.org/10.1177/002224378201900302)

Hausman, J.A. and Ruud, P.A. (1987). Specifying and testing econometric
models for rank-ordered data. *Journal of Econometrics*, 34, 83-104.
[doi:10.1016/0304-4076(87)90068-6](https://doi.org/10.1016/0304-4076%2887%2990068-6)

## See also

[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)

Other specification tests:
[`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md),
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
[`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md),
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)

## Examples

``` r
set.seed(5)
d <- rrologit(800, beta = cbind(b = c(0.5, -1), c = c(-0.3, 0.6), d = c(0.2, 0.4)))
rank_test(rologit(rank ~ 1 | x1, data = d, id = "id", alt = "alt"))
#> Stability of preferences across ranking stages (position 1 against later; 800 ranking(s) with later stages) 
#>              test statistic df p.value
#>  Likelihood ratio     7.771  6   0.255
#>              Wald     7.756  6   0.257
```
