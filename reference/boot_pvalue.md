# Parametric-bootstrap p-values for a specification test

Calibrates any test in the package (or a user function returning a data
frame with a `statistic` column) by the parametric bootstrap: `B`
outcome vectors are drawn from the fitted model with the covariates held
fixed, the model is refit on each (with the fit's own covariance type; a
bootstrap covariance is replaced by the analytic one inside the
replicates), the test is recomputed, and the p-value of every row is
\\(1 + k) / (B' + 1)\\, with \\B'\\ the replicates that succeeded and
\\k\\ the number of their statistics at least as large as the observed
one. This is the reference distribution that remains valid when the
asymptotic chi-square is unreliable (small samples, rare alternatives,
Hausman statistics with a non-positive-definite covariance difference).

## Usage

``` r
boot_pvalue(object, test, B = 199L, cores = 1L, ...)
```

## Arguments

- object:

  A fitted
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  model.

- test:

  A test function taking the fitted
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  model as its first argument, whose null hypothesis is that model, and
  returning a data frame with a `statistic` column (e.g.
  [`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
  [`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
  [`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
  [`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md),
  [`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md);
  not
  [`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md),
  which returns a list, and not
  [`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)
  or
  [`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
  whose null is a restriction of the model).

- B:

  Number of bootstrap replicates (default 199).

- cores:

  Number of parallel workers.

- ...:

  Further arguments to `test`.

## Value

The test's data frame with the column `p.boot` added (and `p.value`
kept), the attribute `"boot_statistics"` holding the \\B\\ replicate
statistics per row, `"B"`, `"n_failed"` (replicates whose refit or test
failed, dropped), `"n_flagged"` (replicates whose refit raised a
warning; they are kept, and the count is printed with the method line),
`"flagged"` (which replicates, aligned with the columns of
`"boot_statistics"`), `"flag_reasons"` (a table of the first warning of
each flagged refit), and, when some but not all replicates are flagged,
`"p.boot_unflagged"`, the bootstrap p-values over the unflagged
replicates. When the observed fit raised no such warning, those are the
replicates that resemble the data: a separated simulated sample can
produce an extreme statistic that says nothing about the data at hand.

## Details

Frequency weights (integers) are expanded to one row per case, for the
simulated data and for the observed statistic alike, so that both follow
one definition; a test that refits from the data rather than from the
fitted design cannot be calibrated with weights (the package's own tests
all work from the fitted design; a user-written test that refits with
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) from
`object$data` is one that refits from the data), and a test that calls
[`update()`](https://rdrr.io/r/stats/update.html) on the fit
re-evaluates the original call on the original data, so it must work
from the fitted object (its `data` component carries the simulated
response). Replicates whose refit warns (a separated simulated sample,
an ill-conditioned information matrix) are kept and counted
(`n_flagged`): a few are expected with rare alternatives; when they are
many, the statistic itself is unstable on data like these, and the
p-value should be read with that in mind. The attribute
`p.boot_unflagged` leaves the flagged replicates out. It is a p-value
conditional on unflagged samples, whose size has not been studied, and
is no substitute for `p.boot`. A replicate in which an alternative is
never chosen, and whose maximum-likelihood estimate therefore does not
exist, is dropped and counted in `n_failed`, so the bootstrap
distribution is conditional on the estimate existing. For
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
without fixed effects,
[`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md)
decides that (an alternative without a constant can be never chosen and
still have an estimate); for the other models an alternative with
coefficients of its own must be chosen
([`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
keeps such replicates where the separation check applies). For a
penalized fit, whose estimates exist without every alternative chosen,
such replicates are kept and refit under the same penalty.

The replicates are drawn from the fitted model, so the test's null
hypothesis must be that model: the IIA tests,
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md),
[`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md),
and user-written tests of that kind.
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)
and
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md)
test a restriction of the fitted model, which the replicates would not
satisfy, and are refused; for a term, fit the model without it and use
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md).
The package's tests derived for maximum-likelihood estimates refuse a
penalized fit (`firth = TRUE` or a prior), so they cannot be calibrated
on one here.

Every replicate's outcomes are drawn in the session before any refit
runs, so a result under
[`set.seed()`](https://rdrr.io/r/base/Random.html) does not depend on
`cores` (for a test with random draws of its own, such as
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
when its `seed` is fixed). The replicates draw cases independently: with
clustered data the bootstrap distribution ignores the clustering.

## See also

[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md)

Other IIA tests:
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md),
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md),
[`small_hsiao()`](https://bagozzib.github.io/polychoice/reference/small_hsiao.md)

## Examples

``` r
# \donttest{
set.seed(3)
d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
m <- mnl(y ~ x1 + x2, data = d)
boot_pvalue(m, hausman_mcfadden, B = 99, drop = "b")
#> Hausman-McFadden test of IIA (variance difference) -- parametric bootstrap, 99 replicates 
#>  alternative statistic df p.value   psd n_restricted p.boot
#>            b     1.537  6   0.957 FALSE          483   0.39
#> Note: psd = FALSE: the covariance difference is not positive semi-definite. A negative statistic is no evidence against IIA (p-value 1); a positive one is unreliable (the difference form over-rejects); use vcov = "joint". 
# }
```
