# Predictions from a choice model

Predictions from a choice model

## Usage

``` r
# S3 method for class 'choicemodel'
predict(
  object,
  newdata = NULL,
  type = c("prob", "class", "link", "prob_outcome", "regime", "prob_conditional",
    "posterior", "ranking"),
  offset = NULL,
  ...
)
```

## Arguments

- object:

  A `"choicemodel"` fit.

- newdata:

  Optional data frame (same layout as the fit: for the long layout it
  must contain the `id` and `alt` columns, one row per available
  alternative; the choice indicator is not needed). It must hold every
  covariate of every equation (for an inflated fit the split-equation
  variables too) and, for a fixed-effects fit, the unit column (only
  units seen in the fit). Factor levels must be among those of the fit.
  With `NULL`, the estimation cases.

- type:

  `"prob"` (default): the \\n \times J\\ matrix of probabilities of the
  observed alternatives, \\P^\*\\ (for a random-intercept fit the
  population-averaged probabilities); `"class"`: the modal alternative
  (a factor); `"link"`: the \\n \times J\\ matrix of utilities
  \\V\_{ij}\\ (\\-\infty\\ where unavailable); for nested-logit fits
  `"prob_conditional"`: the probabilities of the alternatives within
  their nests, \\P(j \mid \text{nest})\\; for inflated models also
  `"prob_outcome"`: the outcome-stage probabilities \\P\_{ij}\\ (the
  choice probabilities of a case in the choice regime, which is what
  Bagozzi and Marchetti (2017) hold the split equation at when comparing
  BIMNL and MNL effects), `"regime"`: \\P(s_i = 1 \mid z_i)\\, the
  probability of not being an inflated (always-\\k\\) case (with several
  inflated alternatives the matrix of the regime probabilities: the
  choice regime, then captivity to each inflated alternative), and
  `"posterior"` (estimation cases only): the probability of the choice
  regime given the observed choice, \\\phi_i P\_{ik} / P^\*\_{ik}\\ for
  a case that chose \\k\\ and 1 otherwise, which separates routine from
  occasional choosers of \\k\\ (with several inflated alternatives the
  matrix of the regime posteriors); for rank-ordered logit fits also
  `"ranking"` (estimation cases only): the probability of each observed
  ranking. For a rank-ordered fit `"prob"` and `"class"` describe first
  choices.

- offset:

  Offset for `newdata` when the fit has one (an \\n \times J\\ matrix in
  the wide layout, a vector aligned to the rows in the long layout; a
  long-layout fit whose offset was a column name finds it in `newdata`).
  A fit with an offset refuses `newdata` without one.

- ...:

  Unused.

## Value

A matrix (rows named by case), a factor, or a numeric vector.
