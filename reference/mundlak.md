# Mundlak (correlated random effects) device

Augments a data frame with the unit-level means of the numeric
case-specific covariates in the formula (every part of a multipart
formula; in the long layout only variables that are constant within a
case are averaged, over the unit's cases), and returns the augmented
formula and data. A covariate constant within every unit is not averaged
(it is its own unit mean; its coefficient is the between-unit effect),
with a message. Fitting
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md) on
the result implements the Mundlak / correlated-random-effects
specification: the coefficients on the original covariates are the
within-unit effects while those on the means capture (and test) the
correlation between the covariates and the unit effect, without the
incidental-parameters bias of unit dummies in short panels. A pooled
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) fit of
the augmented formula with `cluster = "unit"` is the population-averaged
version: its slopes are attenuated towards zero by the unit-effect
variance (in the vignette's example by about a quarter), as in any
nonlinear model, and are consistent for average effects, not for the
within-unit coefficients.

## Usage

``` r
mundlak(formula, data, unit, id = NULL, weights = NULL, suffix = "_mean")
```

## Arguments

- formula:

  A model formula, possibly with `|` parts.

- data:

  A data frame.

- unit:

  Column name identifying the unit.

- id:

  Long layout: the case column (so that each case counts once in the
  unit mean); `NULL` for the wide layout.

- weights:

  Optional frequency weights (a column name or a vector aligned to the
  rows of `data`); the unit means are then weighted means, as the fit
  will be.

- suffix:

  Suffix of the added mean columns (default `"_mean"`).

## Value

A list with `formula` (the augmented formula, keeping the `|`
structure), `data`, and `added` (the names of the mean columns).

## References

Mundlak, Y. (1978). On the pooling of time series and cross section
data. *Econometrica*, 46, 69-85.
[doi:10.2307/1913646](https://doi.org/10.2307/1913646)

## See also

[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md); the
`fe` and `cluster` arguments of the estimators.

## Examples

``` r
## 40 units of 15 cases; each unit's effect on alternative b is
## correlated with its mean of x1
set.seed(1)
unit <- rep(1:40, each = 15); u <- rnorm(40)
x1 <- rnorm(600) + u[unit]; x2 <- rnorm(600)
V <- cbind(0, 0.5 - x1 + 0.8 * x2 + u[unit], -0.3 + 0.6 * x1 - 0.5 * x2)
y <- apply(exp(V) / rowSums(exp(V)), 1, function(p) sample(c("a", "b", "c"), 1, prob = p))
d <- data.frame(y = factor(y), x1, x2, unit)
m <- mundlak(y ~ x1 + x2, d, unit = "unit")
m$formula
#> y ~ x1 + x2 + x1_mean + x2_mean
#> <environment: 0x5de62a803248>
fit <- mnl(m$formula, m$data, cluster = "unit")
## the unit mean of x1 picks up the unit effect in b's equation; without
## the means x1's slope absorbs it. This pooled fit estimates average
## effects: its x1 slope is close to the within-unit -1 of the simulation
## only because the unit mean takes up nearly all of the unit effect here
## (15 cases per unit). remnl() on m$formula estimates the within-unit
## coefficients.
coef(fit)[grep("x1", names(coef(fit)))]
#>       b:x1  b:x1_mean       c:x1  c:x1_mean 
#> -1.0542526  0.7851798  0.5203237 -0.1141086 
coef(mnl(y ~ x1 + x2, d, cluster = "unit"))[c("b:x1", "c:x1")]
#>       b:x1       c:x1 
#> -0.6107817  0.4668951 
```
