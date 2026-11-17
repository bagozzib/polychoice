# First differences in predicted probabilities

The change in the predicted probability of every alternative when a
covariate moves from `from` to `to`, holding the other covariates at a
profile (the estimation cases' weighted means and modes; the default in
the wide layout), at a supplied one-case `newdata`, or averaged over the
estimation cases (`average = TRUE`, the average first difference /
average marginal effect; the default in the long layout). Intervals by
the delta method (numerical Jacobian of the quantity with respect to the
parameters) or by simulation from the asymptotic normal distribution of
the estimates (the practice of Bagozzi and Marchetti 2017; this is not a
bootstrap). A fit with an offset carries it along: the estimation
offsets for `average = TRUE`, their weighted mean at a profile, and
`offset` must be given with `newdata`. For a fixed-effects fit the
profile holds the unit constants at the modal unit's values, an
arbitrary choice in a balanced panel: `average = TRUE` (the effect
averaged over the estimation cases and their units) is the usual
quantity there.

## Usage

``` r
first_difference(
  object,
  var,
  from,
  to,
  newdata = NULL,
  average = NULL,
  alternative = NULL,
  type = c("prob", "prob_outcome", "regime"),
  stage = c("both", "outcome", "inflation"),
  ci = c("delta", "sim", "none"),
  level = 0.95,
  R = 1000,
  offset = NULL
)
```

## Arguments

- object:

  A `"choicemodel"` fit.

- var:

  Name of the covariate to change.

- from, to:

  Its two values (a factor level for a factor).

- newdata:

  Optional data frame giving the profile (one case: one row in the wide
  layout, that case's rows in the long layout).

- average:

  If `TRUE`, average the first differences over the estimation cases
  (weighted by the frequency weights). `NULL` (default) means `FALSE` in
  the wide layout (a profile) and `TRUE` in the long layout, where there
  is no natural one-case profile.

- alternative:

  Long layout: the alternative whose rows receive the change (for an
  alternative-specific covariate); `NULL` changes all rows.

- type:

  `"prob"` (default), `"prob_outcome"`, or `"regime"` (inflated models).

- stage:

  `"both"` (default), `"outcome"`, or `"inflation"`: which equation(s)
  see the change (inflated models).

- ci:

  `"delta"` (default), `"sim"` (simulation from the asymptotic normal
  distribution of the estimates), or `"none"`.

- level:

  Confidence level.

- R:

  Number of simulation draws for `ci = "sim"`.

- offset:

  Offset for `newdata` when the fit has one (an n x J matrix in the wide
  layout, a vector aligned to the rows in the long layout).

## Value

A data frame of class `"pc_fd"` with one row per alternative (for
`type = "regime"`, one row per regime: one with a single inflated
alternative, the choice regime and each captive regime with several):
`from` and `to` are the predicted probabilities at the two covariate
values, `diff` their difference, and `lower`, `upper` the interval for
`diff`.

## Details

For inflated models, `type = "prob"` moves the probabilities of the
observed alternatives (regime mixture included) and
`type = "prob_outcome"` the outcome-stage probabilities, that is, the
effects among cases in the choice regime, which is what Bagozzi and
Marchetti (2017) report when comparing BIMNL with MNL effects;
`type = "regime"` reports the change in \\P(s = 1)\\ (with several
inflated alternatives, the probabilities of the choice regime and of
captivity to each). `stage` restricts the covariate change to one
equation when the covariate enters both: with `stage = "outcome"` the
covariate moves in the outcome equations while the split equation keeps
it at `from`, and conversely for `"inflation"`. The two stage effects
need not add up to the total (`"both"`): the mixture is nonlinear, so
the remainder is an interaction of the two changes.

In the long layout the covariate is changed on every row (a
case-specific covariate) or only on the rows of one alternative
(`alternative =`, for an alternative-specific covariate), and
`average = TRUE` is the default since there is no natural one-case
profile.

## References

Bagozzi, B.E. and Marchetti, K. (2017). Distinguishing occasional
abstention from routine indifference in models of vote choice.
*Political Science Research and Methods*, 5(2), 277-294.
[doi:10.1017/psrm.2015.42](https://doi.org/10.1017/psrm.2015.42)

## See also

[`predict.choicemodel()`](https://bagozzib.github.io/polychoice/reference/predict.choicemodel.md),
[`rrr()`](https://bagozzib.github.io/polychoice/reference/rrr.md)

Other quantities of interest:
[`rrr()`](https://bagozzib.github.io/polychoice/reference/rrr.md)

## Examples

``` r
set.seed(1)
d <- rmnl(1000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
m <- mnl(y ~ x1 + x2, data = d)
first_difference(m, "x1", 0, 1)
#> First differences: x1 from 0 to 1 at the covariate profile (delta intervals at 95%)
#>  alternative  from    to   diff  lower  upper
#>            a 0.271 0.324  0.053  0.018  0.089
#>            b 0.471 0.200 -0.271 -0.303 -0.239
#>            c 0.258 0.475  0.218  0.176  0.259
first_difference(m, "x1", 0, 1, average = TRUE, ci = "sim", R = 500)
#> First differences: x1 from 0 to 1, averaged over the estimation cases (simulation (R=500) intervals at 95%)
#>  alternative  from    to   diff  lower  upper
#>            a 0.254 0.299  0.045  0.012  0.081
#>            b 0.474 0.238 -0.236 -0.264 -0.208
#>            c 0.271 0.463  0.191  0.153  0.232
```
