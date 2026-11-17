# Small-Hsiao test of IIA

The likelihood-ratio form of the IIA test with a random half-split of
the cases (Small and Hsiao 1985). The full model is estimated on both
halves, \\\hat\theta_1\\ and \\\hat\theta_2\\, and combined as
\\\hat\theta\_{12} = \sqrt{1/2}\\\hat\theta_1 + (1 -
\sqrt{1/2})\\\hat\theta_2\\; the restricted model (alternative \\d\\
dropped) is estimated on the second half, \\\hat\theta_r\\, and the
statistic is \\-2\[\ell_2^r(\hat\theta\_{12}) -
\ell_2^r(\hat\theta_r)\]\\ with the second-half restricted
log-likelihood \\\ell_2^r\\, chi-square on as many degrees of freedom as
restricted parameters under IIA. The weights \\\sqrt{1/2}\\ and \\1 -
\sqrt{1/2}\\ are what make the statistic chi-square; the variant that
drops them and also estimates the restricted model on both halves (as in
some earlier implementations) is twice a chi-square and rejects far too
often.

## Usage

``` r
small_hsiao(object, drop = NULL, splits = 1L, seed = NULL)
```

## Arguments

- object:

  A fitted
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  model with at least three alternatives. A fit with diverging
  coefficients (separation) draws a warning: the contrasts then compare
  estimates that do not exist and can reject spuriously.

- drop:

  Alternatives to drop, one at a time; default all.

- splits:

  Number of independent half-splits (default 1, the test; more describe
  how the statistic varies across splits, see Details).

- seed:

  Optional seed for the splits; the caller's random-number stream is
  restored afterwards.

## Value

A data frame of class `"iia_test"` with one row per split and dropped
alternative: `split`, `alternative`, `statistic`, `df`, `p.value`, and
`separated` (a fit behind that row has diverging coefficients in an
equation the statistic uses: the quasi-separation of a half-sample that
a rare alternative with binary covariates invites, after which the
statistic may reflect the separation rather than IIA; a divergence
confined to the dropped alternative's own equation does not mark the
row); with `splits > 1` the attribute `"summary"` holds per alternative
the median statistic over the splits, `median_p` (the chi-square p-value
of that median, as
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)
reports it; descriptive, since the median of split statistics is not
chi-square, and with an even number of splits not the median of the
split p-values), the share of splits rejecting at 5 percent, and the
share of separated splits.
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
calibrates the median (see the IIA vignette).

## Details

The split is stratified by chosen alternative. Results can depend
heavily on the random split: one split may give a large statistic that
the next does not. The test is one split chosen before the data are
seen: fix `seed` in advance (`small_hsiao(m, seed = )`) and report that
split. In the package's size study one split held its level (rejection
0.040 to 0.047 at 0.05), while the median over ten splits rejected in
none of 300 samples, so the median is not a test at its level.
`splits > 1` repeats the statistic on independent splits and reports
every one together with the median and the share of splits rejecting at
the 5 percent level, a description of how much the statistic varies
(Long and Freese 2014 recommend inspecting this variability;
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)
reports the median over ten splits);
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
can calibrate the median. When the pre-specified split's half-sample
fits are separated (`separated`), its statistic reflects the separation
rather than IIA: report the row as uninformative, or calibrate it with
`boot_pvalue(m, small_hsiao, seed = )` and the same seed. A seed chosen
after seeing a result is not a test. Dropping the base alternative is
handled by re-basing, as in
[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md).
The test is not frequency-weight exact: the halves are balanced in
cases, not in weight. The cases are taken as independent: a fit's
`cluster` does not enter the statistic. Not for fixed-effects fits.

## References

Long, J.S. and Freese, J. (2014). *Regression Models for Categorical
Dependent Variables Using Stata*, 3rd ed. College Station, TX: Stata
Press.

Small, K.A. and Hsiao, C. (1985). Multinomial logit specification tests.
*International Economic Review*, 26, 619-627.
[doi:10.2307/2526707](https://doi.org/10.2307/2526707)

## See also

[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)

Other IIA tests:
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md),
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md),
[`hausman_mcfadden()`](https://bagozzib.github.io/polychoice/reference/hausman_mcfadden.md),
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md),
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md)

## Examples

``` r
set.seed(1)
d <- rmnl(1200, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.1, 0.3, 0.9)))
m <- mnl(y ~ x1 + x2, data = d)
small_hsiao(m, seed = 1)
#> Small-Hsiao test of IIA (stratified half-split, sqrt(1/2) weighting) 
#>  split alternative statistic df p.value separated
#>      1           a     6.703  6   0.349     FALSE
#>      1           b     3.603  6   0.730     FALSE
#>      1           c     1.479  6   0.961     FALSE
#>      1           d     5.479  6   0.484     FALSE
small_hsiao(m, drop = "c", splits = 5, seed = 1)
#> Small-Hsiao test of IIA (stratified half-split, sqrt(1/2) weighting) 
#>  split alternative statistic df p.value separated
#>      1           c     1.479  6   0.961     FALSE
#>      2           c     6.583  6   0.361     FALSE
#>      3           c     1.779  6   0.939     FALSE
#>      4           c     4.139  6   0.658     FALSE
#>      5           c     6.484  6   0.371     FALSE
#> 
#> Across splits:
#>  alternative median_statistic df median_p share_reject_05 share_separated
#>            c            4.139  6    0.658               0               0
```
