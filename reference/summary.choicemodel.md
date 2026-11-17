# Summarize a choice model

Summarize a choice model

## Usage

``` r
# S3 method for class 'choicemodel'
summary(object, ...)
```

## Arguments

- object:

  A `"choicemodel"` fit.

- ...:

  Unused.

## Value

An object of class `"summary.choicemodel"` with the coefficient table
(`coefficients`: estimate, standard error, z, p), the parameter `block`
of every row (`"outcome"`, `"generic"`, `"fixed effect"`, `"inflation"`,
`"nest"`, `"scale"`, `"covariance"`, `"random intercept"`), and the fit
statistics. The z statistics of the inclusive-value and scale blocks
test the value 1 (the multinomial logit); those of the multinomial
probit covariance block test the values under independent, identically
distributed errors (variances 1, covariances 0.5 under the normalization
of the first variance to 1); the random-intercept standard deviations
get no z test (their null is on the boundary); a parameter whose
variance is negative (the Hessian is not positive definite) has standard
error `NA`. The printed table of a penalized fit (`firth = TRUE` or a
prior) marks `[sep]` the coefficients that move along a separating
direction of the data
([`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md)).
Their Wald tests rest on a quadratic approximation of the penalized
log-likelihood that fails along that direction and can be far off
(Heinze and Schemper 2002); under a prior they can have little power
(Rainey 2024). Test them by `term_test(type = "lr")` or
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md).
A maximum-likelihood fit names its diverging coefficients in a warning
at fit time instead.

## See also

[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md)
for the pseudo-R-squared measures.
