# Test whether two alternatives can be combined

Tests the hypothesis that two alternatives are indistinguishable with
respect to the covariates: all coefficients other than the constants are
equal in the two equations (when one of the two is the base, that the
other's coefficients are zero). The Wald form uses the fitted
covariance; the likelihood-ratio form (Cramer and Ridder 1991) refits
the multinomial logit with the two alternatives pooled and adds the
binomial constant-only log-likelihood of the split within the pooled
alternative: \\\ell_R = \ell\_{pooled} + n_j \log(n_j/n\_{jk}) + n_k
\log(n_k/n\_{jk})\\. Both are chi-square on the number of non-constant
terms (Stata's `mlogtest, combine` and `lrcomb`): the restriction leaves
the two constants free, through the pooled equation's constant and the
split share. An alternative that is never chosen contributes \\0 \log 0
= 0\\.

## Usage

``` r
combine_test(object, alternatives = NULL, type = c("wald", "lr"))
```

## Arguments

- object:

  A fitted
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  model (for the LR form: case-specific covariates only, in either
  layout, and every alternative available; the Wald form works for any
  multinomial-logit fit).

- alternatives:

  Two alternatives (labels). `NULL` (default) tests every pair.

- type:

  `"wald"` (default) or `"lr"`.

## Value

A data frame of class `"iia_test"` with one row per pair: `pair`,
`statistic`, `df`, `p.value`.

## Details

The likelihood-ratio form under separation and penalties:

- For a maximum-likelihood fit of separated data it compares the suprema
  of the log-likelihoods (see
  [`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)).

- For a penalized fit (`firth = TRUE` or a prior) it is the penalized
  likelihood ratio, the restriction imposed on the full model's
  penalized log-likelihood (no pooled refit).

The likelihood-ratio form assumes independent cases (a fit's `cluster`
does not enter it); the Wald form uses the fit's covariance,
cluster-robust when the fit has one.

## References

Cramer, J.S. and Ridder, G. (1991). Pooling states in the multinomial
logit model. *Journal of Econometrics*, 47, 267-272.
[doi:10.1016/0304-4076(91)90102-J](https://doi.org/10.1016/0304-4076%2891%2990102-J)

Long, J.S. and Freese, J. (2014). *Regression Models for Categorical
Dependent Variables Using Stata*, 3rd ed. Stata Press.

## See also

[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)

Other specification tests:
[`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md),
[`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md),
[`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md),
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)

## Examples

``` r
set.seed(1)
## alternatives b and c share their slopes
d <- rmnl(1000, beta = cbind(b = c(0.5, -1, 0.8), c = c(0.3, -1, 0.8)))
m <- mnl(y ~ x1 + x2, data = d)
combine_test(m)
#> Test of combining alternatives (Wald) 
#>   pair statistic df p.value
#>  a & b   137.974  2  <2e-16
#>  a & c   110.967  2  <2e-16
#>  b & c     1.839  2   0.399
combine_test(m, c("b", "c"), type = "lr")
#> Test of combining alternatives (likelihood ratio, Cramer-Ridder) 
#>   pair statistic df p.value
#>  b & c     1.845  2   0.397
```
