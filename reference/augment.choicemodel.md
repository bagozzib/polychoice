# Augment the estimation cases with fitted classes and probabilities (broom method)

Augment the estimation cases with fitted classes and probabilities
(broom method)

## Usage

``` r
# S3 method for class 'choicemodel'
augment(x, ...)
```

## Arguments

- x:

  A `"choicemodel"` fit.

- ...:

  Unused.

## Value

The case-level estimation data (one row per case) with `.fitted` (the
modal alternative), `.prob_<level>` columns, and for inflated models
`.regime` (fitted \\P(s = 1)\\).

## See also

Other broom methods:
[`glance.choicemodel()`](https://bagozzib.github.io/polychoice/reference/glance.choicemodel.md),
[`tidy.choicemodel()`](https://bagozzib.github.io/polychoice/reference/tidy.choicemodel.md)
