# Classification table of a choice model

Cross-tabulates the chosen alternative against the modal predicted
alternative (weighted by the frequency weights), with the share
correctly classified and the proportional reduction in error relative to
always predicting the most frequent alternative. Modal classification
understates a model's value for minority alternatives: a case with
probability 0.4 on a rare alternative is classified to the common one
even though 0.4 may be far above its base rate; see the
probability-based measures in
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md).

## Usage

``` r
classification(object, newdata = NULL, weights = NULL, offset = NULL)
```

## Arguments

- object:

  A `"choicemodel"` fit.

- newdata:

  Optional data frame (same layout as the fit) with the response;
  default the estimation cases. A hold-out set may lack some alternative
  or covariate level; only the response and the predictors are read from
  it.

- weights:

  Optional frequency weights for `newdata` (a column name or a vector
  aligned to its rows); the estimation cases use the fit's weights.

- offset:

  Offset for `newdata` when the fit has one (see
  [`predict.choicemodel()`](https://bagozzib.github.io/polychoice/reference/predict.choicemodel.md)).

## Value

An object of class `"pc_classification"`: a list with `table` (observed
in rows, predicted in columns), `correct`, and `pre`.

## See also

[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md)

Other model comparison:
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md),
[`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md),
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)

## Examples

``` r
set.seed(1)
d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
classification(mnl(y ~ x1 + x2, data = d))
#>         predicted
#> observed   a   b   c
#>        a   0  98  89
#>        b   0 347  46
#>        c   0  73 147
#> 
#> Correctly classified: 0.618; proportional reduction in error over the modal alternative: 0.248
```
