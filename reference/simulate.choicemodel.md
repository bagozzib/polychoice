# Simulate outcomes from a fitted choice model

Draws `nsim` outcome vectors from the fitted model (holding the
covariates fixed): the parametric bootstrap of the model. For every
family but the random-intercept one the draw of a case is from its
fitted probabilities of the observed alternatives (the GHK probabilities
of a multinomial probit case are renormalized to one); for
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)
one intercept vector is drawn per unit and the unit's cases are then
drawn given it, so the within-unit dependence of the model is
reproduced. Frequency weights are left as they are: one draw per row, so
for a weighted fit the draws are not the parametric bootstrap of the
weighted data (a row of weight w would need w draws;
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
expands the rows for that reason). For nominal outcomes there is no
meaningful residual ordering, so these draws serve test calibration and
predictive checks rather than DHARMa-style quantile residuals.

## Usage

``` r
# S3 method for class 'choicemodel'
simulate(object, nsim = 1, seed = NULL, ..., ties = c("untied", "observed"))
```

## Arguments

- object:

  A `"choicemodel"` fit.

- nsim:

  Number of simulated outcome vectors.

- seed:

  Optional seed (as in
  [`stats::simulate()`](https://rdrr.io/r/stats/simulate.html); the
  caller's random-number stream is restored afterwards).

- ...:

  Unused.

- ties:

  For a
  [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
  fit: `"untied"` (default), every case ranking as many positions as it
  ranked in the data, without ties; or `"observed"`, each case's
  observed pattern of rank groups (the same group sizes in the same
  order, a tie staying a tie at its positions), as
  [`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
  draws rankings.

## Value

A data frame with `nsim` factor columns (`sim_1`, ...), one row per
estimation case. For a
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
fit, a list of `nsim` integer matrices (cases by alternatives) of
simulated ranks (0 for the unranked and unavailable alternatives).
