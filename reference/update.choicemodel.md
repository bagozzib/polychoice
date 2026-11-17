# Update and refit a choice model

Works on the stored call as
[`stats::update()`](https://rdrr.io/r/stats/update.html) does, with
multi-part formulas handled by
[`Formula::update.Formula()`](https://rdrr.io/pkg/Formula/man/Formula.html):
`update(fit, . ~ . + x3)` changes the first part,
`update(fit, . ~ . | . + z2)` the second.

## Usage

``` r
# S3 method for class 'choicemodel'
update(object, formula., ..., evaluate = TRUE)
```

## Arguments

- object:

  A `"choicemodel"` fit.

- formula.:

  Changes to the formula (optional).

- ...:

  Other arguments of the original call to replace.

- evaluate:

  If `FALSE`, return the updated call instead of refitting.

## Value

The refitted model (or the call).
