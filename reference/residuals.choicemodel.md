# Residuals of a choice model

Residuals of a choice model

## Usage

``` r
# S3 method for class 'choicemodel'
residuals(object, type = c("response", "pearson"), ...)
```

## Arguments

- object:

  A `"choicemodel"` fit.

- type:

  `"response"` (the \\n \times J\\ matrix of indicator minus fitted
  probability) or `"pearson"` (response residuals divided by
  \\\sqrt{P(1-P)}\\).

- ...:

  Unused.

## Value

An \\n \times J\\ matrix.
