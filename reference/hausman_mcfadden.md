# Hausman-McFadden test of IIA

For each dropped alternative, compares the multinomial logit estimated
on the cases that did not choose it (restricted choice set) with the
full-sample estimate of the same parameters: under IIA both are
consistent and the restricted one is less efficient, so \\(\hat\beta_r -
\hat\beta_f)'(V_r - V_f)^{-1}(\hat\beta_r - \hat\beta_f)\\ is
asymptotically chi-square with as many degrees of freedom as restricted
parameters (Hausman and McFadden 1984).

## Usage

``` r
hausman_mcfadden(object, drop = NULL, vcov = c("difference", "joint"))
```

## Arguments

- object:

  A fitted
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  model with at least three alternatives. A fit with diverging
  coefficients (separation) draws a warning: the contrasts then compare
  estimates that do not exist and can reject spuriously.

- drop:

  Alternatives to drop, one at a time (labels or indices); default all
  of them.

- vcov:

  `"difference"` (default) or `"joint"`; see Details.

## Value

A data frame of class `"iia_test"` with one row per dropped alternative:
`alternative`, `statistic`, `df`, `p.value`, `psd` (was the covariance
of the difference positive semi-definite), and `n_restricted`.

## Details

Dropping the base alternative is handled exactly: the full model is
re-expressed with the first remaining alternative as the base through
the linear map \\\beta_j^{(c)} = \beta_j - \beta_c\\, and the restricted
model is estimated with that base, so the compared parameters coincide
(the misalignment that produces spurious rejections for the base-drop
row in some implementations does not arise).

`vcov = "difference"` is the classic form, with both covariances the
inverse observed information (the form the Hausman lemma requires; the
fit's own robust or cluster covariance is not used here). The difference
of covariance matrices need not be positive semi-definite in finite
samples; the result then reports `psd = FALSE`, a negative statistic is
shown as computed, and its p-value is 1 (Hausman and McFadden 1984, p.
1226, and Long and Freese 2014 read a negative statistic as no evidence
against IIA). A positive statistic with `psd = FALSE` is unreliable: in
the package's size study the difference form rejected well above the
nominal level. It stays the default because it is the form the
literature and Stata's `hausman` report, so that published values can be
reproduced; the joint form is the one to rely on. `vcov = "joint"`
instead estimates the covariance of the difference from the stacked
per-case scores of the two fits (the seemingly-unrelated-estimation
device of Weesie 1999; Stata's `mlogtest, suest`): \\V_r + V_f - C -
C'\\ with the sandwich covariances and cross-covariance, positive
semi-definite by construction and valid with heteroskedasticity and,
when the fit has a `cluster`, clustering.

## References

Hausman, J.A. and McFadden, D. (1984). Specification tests for the
multinomial logit model. *Econometrica*, 52, 1219-1240.
[doi:10.2307/1910997](https://doi.org/10.2307/1910997)

Weesie, J. (1999). Seemingly unrelated estimation and the
cluster-adjusted sandwich estimator. *Stata Technical Bulletin*, 52,
34-47.

Long, J.S. and Freese, J. (2014). *Regression Models for Categorical
Dependent Variables Using Stata*, 3rd ed. Stata Press.

## See also

[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)

Other IIA tests:
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md),
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md)

## Examples

``` r
set.seed(1)
d <- rmnl(1200, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
m <- mnl(y ~ x1 + x2, data = d)
hausman_mcfadden(m)
#> Hausman-McFadden test of IIA (variance difference) 
#>  alternative statistic df p.value   psd n_restricted
#>            a     7.173  6   0.305 FALSE          971
#>            b     0.027  6   1.000 FALSE          731
#>            c     1.059  6   0.983  TRUE          952
#>            d   -66.907  6   1.000 FALSE          946
#> Note: psd = FALSE: the covariance difference is not positive semi-definite. A negative statistic is no evidence against IIA (p-value 1); a positive one is unreliable (the difference form over-rejects); use vcov = "joint". 
hausman_mcfadden(m, drop = "b", vcov = "joint")
#> Hausman-McFadden test of IIA (joint (suest) covariance) 
#>  alternative statistic df p.value  psd n_restricted
#>            b     3.354  6   0.763 TRUE          731
```
