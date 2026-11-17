# Testing IIA and specification

The multinomial logit’s independence of irrelevant alternatives (IIA) is
a testable restriction, and the tests of it have a troubled history: the
most cited one, in its classic form, often yields a covariance
difference that is not positive semidefinite and a negative statistic,
and implementations of the others differ in details that change the
answer (which half of the sample is used, which weights, which
covariance). `polychoice` implements the existing tests with their
definitions fixed in the package’s design document, checks them by
simulation, and reports what the simulation found. This vignette walks
through the battery on the Pennsylvania vote data and then through the
further specification tests.

``` r

library(polychoice)
data(pa2006)
m <- mnl(vote ~ canvass + phone + democrat + republican + age + female + voted2004,
         data = pa2006)
```

## The IIA battery

### Hausman and McFadden (1984)

Under IIA, dropping an alternative and refitting on the remaining cases
leaves the common coefficients consistent, so the full-sample and
restricted estimates should agree. The classic statistic uses the
*difference* of the two covariance matrices,
$`q = (\hat\beta_r - \hat\beta_f)'(V_r - V_f)^{-1}(\hat\beta_r - \hat\beta_f)`$,
which need not be positive semidefinite in finite samples; a negative
statistic is reported with a p-value of 1 and a flag. The *joint* form
estimates the covariance of $`\hat\beta_r - \hat\beta_f`$ from the
stacked scores of the two fits (a seemingly-unrelated-estimation
sandwich), which is always positive semidefinite. Dropping the base
alternative re-bases the model exactly, so every alternative can be
dropped.

``` r

hausman_mcfadden(m)
#> Hausman-McFadden test of IIA (variance difference) 
#>       alternative statistic df p.value   psd n_restricted
#>           abstain    -1.027 16       1 FALSE         1267
#>          democrat     1.586 16       1 FALSE         1419
#>  republican_other   -56.438 16       1 FALSE         1435
#>       nonresponse    -1.394 16       1 FALSE         1873
#> Note: psd = FALSE: the covariance difference is not positive semi-definite. A negative statistic is no evidence against IIA (p-value 1); a positive one is unreliable (the difference form over-rejects); use vcov = "joint".
hausman_mcfadden(m, vcov = "joint")
#> Hausman-McFadden test of IIA (joint (suest) covariance) 
#>       alternative statistic df p.value  psd n_restricted
#>           abstain    26.165 16  0.0517 TRUE         1267
#>          democrat    14.386 16  0.5700 TRUE         1419
#>  republican_other    11.353 16  0.7872 TRUE         1435
#>       nonresponse    11.398 16  0.7843 TRUE         1873
```

### Small and Hsiao (1985)

The sample is split in half. The full model is fit on both halves and
the two estimates are combined with the weights $`\sqrt{1/2}`$ and
$`1 - \sqrt{1/2}`$; the model without the dropped alternative is fit on
the second half; the statistic is minus twice the difference between the
second-half restricted log-likelihood evaluated at the combined
estimates and at its own maximum. It depends on the split, so `splits`
repeats it over several random splits and `seed` makes the result
reproducible.

``` r

small_hsiao(m, splits = 10, seed = 1)
#> Small-Hsiao test of IIA (stratified half-split, sqrt(1/2) weighting) 
#>  split      alternative statistic df  p.value separated
#>      1          abstain    61.525 16 2.89e-07      TRUE
#>      1         democrat    59.897 16 5.45e-07      TRUE
#>      1 republican_other    59.622 16 6.06e-07      TRUE
#>      1      nonresponse    19.684 16   0.2348     FALSE
#>      2          abstain    27.239 16   0.0389     FALSE
#>      2         democrat    18.066 16   0.3200     FALSE
#>      2 republican_other    18.881 16   0.2749     FALSE
#>      2      nonresponse    28.923 16   0.0245     FALSE
#>      3          abstain     9.334 16   0.8990     FALSE
#>      3         democrat    12.213 16   0.7292     FALSE
#>      3 republican_other     5.763 16   0.9905     FALSE
#>      3      nonresponse    10.432 16   0.8431     FALSE
#>      4          abstain    52.419 16 9.37e-06      TRUE
#>      4         democrat    51.690 16 1.23e-05      TRUE
#>      4 republican_other    56.511 16 2.00e-06      TRUE
#>      4      nonresponse    14.950 16   0.5283     FALSE
#>      5          abstain    52.492 16 9.12e-06      TRUE
#>      5         democrat    52.499 16 9.10e-06      TRUE
#>      5 republican_other    52.305 16 9.78e-06      TRUE
#>      5      nonresponse     7.147 16   0.9703     FALSE
#>      6          abstain    18.735 16   0.2826      TRUE
#>      6         democrat    20.176 16   0.2124      TRUE
#>      6 republican_other    19.087 16   0.2642      TRUE
#>      6      nonresponse    11.524 16   0.7760     FALSE
#>      7          abstain    12.771 16   0.6894     FALSE
#>      7         democrat    18.122 16   0.3168     FALSE
#>      7 republican_other    13.246 16   0.6547     FALSE
#>      7      nonresponse    11.352 16   0.7872     FALSE
#>      8          abstain    14.254 16   0.5798     FALSE
#>      8         democrat    11.058 16   0.8059     FALSE
#>      8 republican_other    12.012 16   0.7432     FALSE
#>      8      nonresponse    16.070 16   0.4481     FALSE
#>      9          abstain     9.231 16   0.9036     FALSE
#>      9         democrat    17.040 16   0.3830     FALSE
#>      9 republican_other    15.078 16   0.5189     FALSE
#>      9      nonresponse    13.059 16   0.6684     FALSE
#>     10          abstain    62.620 16 1.88e-07      TRUE
#>     10         democrat    62.588 16 1.91e-07      TRUE
#>     10 republican_other    57.798 16 1.22e-06      TRUE
#>     10      nonresponse    25.730 16   0.0579     FALSE
#> 
#> Across splits:
#>       alternative median_statistic df median_p share_reject_05 share_separated
#>           abstain           22.987 16    0.114             0.5             0.5
#>          democrat           19.149 16    0.261             0.4             0.5
#>  republican_other           18.984 16    0.269             0.4             0.5
#>       nonresponse           14.004 16    0.598             0.1             0.0
#> separated = TRUE: a fit behind that row has diverging coefficients (quasi-separation in a half-sample); its statistic may reflect the separation rather than IIA.
```

The split matters, and on these data for a reason that has nothing to do
with IIA. Only two of the 125 nonresponders did not vote in 2004. In
four of the ten splits (1, 4, 5, and 10) both fall in the second half,
so in the first half the nonresponse equation separates on `voted2004`
(its coefficient diverges); the combined estimate then gives those two
respondents a near-zero probability in the second half, and the
statistics of the three drops that keep nonresponse exceed 50 on 16
degrees of freedom, while the nonresponse drop, which removes that
equation, stays below 26. The `separated` column marks every row whose
statistic uses a diverging half-sample fit: the three drops that keep
nonresponse in those four splits and in split 6 (where both fall in the
first half; there the separated half carries the smaller weight and the
evaluation half has neither case, so those statistics are small), but
not the nonresponse drops, which remove the separated equation. In the
five splits without separation only two of the 20 statistics reject at
the 5 percent level, both in split 2. A single split can therefore
suggest a decisive rejection or none, and which one is luck: on these
data a split drawn with `small_hsiao(m, seed = s)` is separated for 12
of the seeds 1 to 20. The test is the split fixed before the data are
seen (`small_hsiao(m, seed = )`): report it, and when it is separated,
as split 1 is here, report it as uninformative or calibrate it by the
parametric bootstrap with the same seed. A seed chosen after seeing a
result is not a test. The splits and their summary (the median, the
share of splits rejecting, and the share separated) describe how much
the statistic varies; when many splits separate, the test has little to
say about IIA on data like these. The summary’s `median_p` is the
chi-square p-value of the median statistic, a descriptive figure: the
median of several split statistics is not itself chi-square distributed
(in the package’s size study a single split had nominal size, while the
battery’s median rule never rejected under the null: it is conservative
by construction). When the decision matters, calibrate the median by the
parametric bootstrap: give
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
a function that returns the median statistic as a one-row data frame,
with a fixed seed so that every replicate draws its splits from the same
random stream (they are stratified on each replicate’s simulated
outcome; slow: ten splits refit the model thirty times per replicate, so
this is not run here).

``` r

sh_median <- function(object, drop) {
  s <- attr(small_hsiao(object, drop = drop, splits = 10, seed = 1), "summary")
  data.frame(test = "Small-Hsiao median", statistic = s$median_statistic)
}
boot_pvalue(m, sh_median, B = 99, drop = "abstain")
```

### Fok and Paap (2025)

Fok and Paap (2025) take the pairwise conditional likelihoods of the
multinomial logit and build two tests: $`Q_1`$, a Hausman contrast
between the maximum-likelihood and the composite-likelihood estimates
(the latter with a Godambe sandwich covariance), and $`Q_2`$, a GMM
overidentification test on the stacked pairwise scores. $`Q_2`$ is
computed in the two-step form by default (`gmm = "iterated"` and the
projected C-statistic form are available); with the pairs that involve
the base alternative only and case-specific covariates alone there is no
overidentification and only $`Q_1`$ is reported. Under the multinomial
logit the composite-likelihood estimator is nearly as efficient as the
maximum-likelihood one, so the difference of their two separately
estimated covariance matrices, which Fok and Paap’s $`Q_1`$ uses, is a
small difference of noisy matrices and is indefinite in nearly every
sample (the package’s simulations found 44 percent of these statistics
negative).
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md)
therefore estimates the covariance of the contrast directly from the two
estimators’ stacked influence functions by default, as the joint
Hausman-McFadden form does; `vcov = "difference"` gives Fok and Paap’s
formula.

``` r

fp_tests(m)
#> Fok-Paap (2025) misspecification tests, all pairs 
#>                                            test statistic df p.value  psd pairs
#>                 Fok-Paap CL Hausman (Q1, joint)    42.692 24 0.01078 TRUE   all
#>  Fok-Paap GMM overidentification (Q2, two-step)    49.820 24 0.00149   NA   all
```

### The battery

[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)
runs all of the above with one call, with ten Small-Hsiao splits
reported by their median (the row’s p-value is the chi-square p-value of
the median statistic, a descriptive figure).

``` r

iia_tests(m, seed = 1)
#> Tests of IIA / misspecification of the multinomial logit 
#>                                            test      alternative statistic df
#>                   Hausman-McFadden (difference)          abstain    -1.027 16
#>                   Hausman-McFadden (difference)         democrat     1.586 16
#>                   Hausman-McFadden (difference) republican_other   -56.438 16
#>                   Hausman-McFadden (difference)      nonresponse    -1.394 16
#>                        Hausman-McFadden (joint)          abstain    26.165 16
#>                        Hausman-McFadden (joint)         democrat    14.386 16
#>                        Hausman-McFadden (joint) republican_other    11.353 16
#>                        Hausman-McFadden (joint)      nonresponse    11.398 16
#>                Small-Hsiao (median over splits)          abstain    22.987 16
#>                Small-Hsiao (median over splits)         democrat    19.149 16
#>                Small-Hsiao (median over splits) republican_other    18.984 16
#>                Small-Hsiao (median over splits)      nonresponse    14.004 16
#>                 Fok-Paap CL Hausman (Q1, joint)                     42.692 24
#>  Fok-Paap GMM overidentification (Q2, two-step)                     49.820 24
#>  p.value   psd share_separated
#>  1.00000 FALSE              NA
#>  1.00000 FALSE              NA
#>  1.00000 FALSE              NA
#>  1.00000 FALSE              NA
#>  0.05174  TRUE              NA
#>  0.57001  TRUE              NA
#>  0.78718  TRUE              NA
#>  0.78430  TRUE              NA
#>  0.11408    NA             0.5
#>  0.26101    NA             0.5
#>  0.26950    NA             0.5
#>  0.59838    NA             0.0
#>  0.01078  TRUE              NA
#>  0.00149    NA              NA
#> psd = FALSE: the covariance difference was not positive semi-definite for that row. A negative statistic is no evidence against the model (p-value 1); a positive one from a variance-difference form is unreliable (see the test's help).
#> share_separated > 0: some Small-Hsiao splits had diverging half-sample fits; see small_hsiao(..., splits = ) for the split-level rows.
#> Small-Hsiao (median over splits): descriptive, not a test at its level (it rejected in none of 300 samples in the package's size study); for the test, run one split with its seed fixed in advance (small_hsiao(m, seed = ), or iia_tests(m, splits = 1, seed = )).
```

How to read this table for these data: the joint Hausman-McFadden rows
do not reject at the 5 percent level (the abstention row only narrowly),
and both Fok-Paap tests do. The Small-Hsiao medians are descriptive (the
note under the table), and in three of their rows half the splits are
separated. The difference-form Hausman-McFadden rows rest on covariance
differences that are not positive semidefinite (`psd = FALSE`); one rule
covers every such row: a negative statistic is no evidence against IIA
(its p-value is 1), and a positive one is unreliable, because the
difference form over-rejects (next section), so those rows are
descriptive. The tests look at different departures, so disagreement is
informative rather than contradictory: the Fok-Paap tests use all
pairwise comparisons at once, while the drop tests look at one
alternative at a time (and the Small-Hsiao medians here rest partly on
separated splits). A reasonable report states the battery, notes which
rows had an indefinite covariance difference, and, when the decision
matters, calibrates the decisive statistic by the parametric bootstrap
(next section) or fits the alternative model directly (nested logit,
HEV, random intercepts; the *Beyond IIA* vignette) and compares by
likelihood ratio.

### What the size study found

The package’s verification ledger records a 300-replication study of the
battery on data generated from the multinomial logit (three
alternatives, four covariates, ten parameters, 1000 cases; Monte Carlo
standard error about 0.013), with these rejection rates at the nominal 5
percent level: the Hausman-McFadden difference form was not positive
semidefinite in 79 percent of replications and rejected 8 percent of the
time for a non-base drop and 16 percent for the base drop; its joint
form 2 and 5 percent; Small-Hsiao (one split) 4 and 5 percent, while the
battery’s median over ten splits never rejected; $`Q_1`$ in its joint
form 6 percent; $`Q_2`$ 6 percent. Fok and Paap’s variance-difference
form of $`Q_1`$ rejected 4 percent of the time, but a separate
1,000-replication check of the same design showed that apparent
calibration to be accidental: its covariance difference was never
positive semidefinite, 44 percent of its statistics were negative, and
it rejected 3.3, 4.2, and 4.9 percent of the time at the 1, 5, and 10
percent levels, while the joint form that
[`fp_tests()`](https://bagozzib.github.io/polychoice/reference/fp_tests.md)
uses by default rejected 1.5, 6.2, and 11.2 percent. An indefinite
covariance difference can give a statistic of either sign: every form
refers $`\max(\text{statistic}, 0)`$ to the chi-square, so a negative
statistic has a p-value of 1. The practical recommendations that follow
are:

- read the difference form of Hausman-McFadden as descriptive and rely
  on its joint form, a single Small-Hsiao split with its seed fixed in
  advance (`small_hsiao(m, seed = )`; the median over several splits,
  which
  [`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)
  reports, is descriptive: in the size study it rejected in none of 300
  samples), $`Q_1`$ in its joint form, and $`Q_2`$;
- calibrate any of them with the parametric bootstrap when the decision
  matters (below);
- remember that these are tests of IIA against unspecified alternatives;
  the score test against a nested logit (below) and the likelihood-ratio
  tests against the fitted nested logit, HEV, or random-intercept model
  (vignette *Beyond IIA*) are the sharper instruments when the
  alternative is known.

## Parametric bootstrap calibration

[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
takes a test function whose null hypothesis is the fitted model (the IIA
tests,
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md),
[`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md))
and returns a data frame with a `statistic` column, with its arguments,
simulates `B` data sets from the fitted model, refits, recomputes the
statistic, and reports the bootstrap p-value $`(1 + k) / (B' + 1)`$,
with $`B'`$ the replicates that succeeded and $`k`$ the number of their
statistics at least as large as the observed one, next to the asymptotic
p-value. Replicates whose refit fails or loses an alternative are
dropped and counted (`n_failed`); replicates whose refit raises a
warning (diverging coefficients, a boundary, an ill-conditioned
information matrix) are kept and counted (`n_flagged`, with the warnings
tabulated in the attribute `flag_reasons`); the method line reports each
count that is not zero. It is the remedy for a test whose finite-sample
distribution is in doubt. Here it calibrates the two Fok-Paap tests, the
rows of the battery that reject:

``` r

set.seed(2)
bf <- boot_pvalue(m, fp_tests, B = 399)
bf
#> Fok-Paap (2025) misspecification tests, all pairs -- parametric bootstrap, 399 replicates, 47 with a warning in the refit (attribute "flag_reasons"), up to 47 with an NA statistic left out of a row (attribute "n_na") 
#>                                            test statistic df p.value  psd pairs
#>                 Fok-Paap CL Hausman (Q1, joint)    42.692 24 0.01078 TRUE   all
#>  Fok-Paap GMM overidentification (Q2, two-step)    49.820 24 0.00149   NA   all
#>  p.boot
#>  0.0075
#>  0.0170
attr(bf, "flag_reasons")
#> The data are quasi-completely separated (see separation()) 
#>                                                         47
attr(bf, "p.boot_unflagged")
#> [1] 0.008498584 0.016997167
```

Calibrated, both Fok-Paap tests still reject. Two of the 399 replicate
statistics reach the observed $`Q_1`$, so its `p.boot` is 0.0075. The 47
flagged replicates are simulated samples in which the nonresponse
equation separates on `voted2004`, the same mechanism as in the
Small-Hsiao splits above: its intercept falls and its `voted2004`
coefficient rises without bound, their sum staying finite. Along that
one direction neither estimator has an estimate, so $`Q_1`$ leaves it
out and compares the two estimators on the other 23, and $`Q_2`$ is not
computed in those samples
([`?fp_tests`](https://bagozzib.github.io/polychoice/reference/fp_tests.md)).
Five of the other 352 replicate statistics reach the observed $`Q_2`$,
and its `p.boot` is 0.017. The observed data do not separate, and over
the unflagged replicates the two p-values are 0.0085 and 0.017.

## The score test against a nested logit

McFadden’s (1987) Lagrange-multiplier test asks whether the
inclusive-value parameters of a specified nest structure differ from 1
(the multinomial logit) without fitting the nested model;
[`nested_lm_test()`](https://bagozzib.github.io/polychoice/reference/nested_lm_test.md)
computes it from the multinomial-logit scores (the outer-product
information by default, the form of McFadden’s regression test). Here
the nests group the two forms of non-participation against the two forms
of voting.

``` r

nests <- list(out = c("abstain", "nonresponse"), vote = c("democrat", "republican_other"))
nested_lm_test(m, nests)
#> Lagrange-multiplier test of IIA against a nested logit 
#>                        test statistic df p.value
#>  Nested-logit LM test (opg)     2.299  2   0.317
```

## Further specification tests

[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md)
asks whether two alternatives can be pooled (Cramer and Ridder 1991):
the Wald form tests the equality of their slope coefficients (the
constants stay free), the likelihood-ratio form compares the pooled fit
with a binomial split term.
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)
tests every covariate jointly across equations (Wald or LR).
[`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md)
is the Hosmer-Lemeshow-type goodness-of-fit test for multinomial models
of Fagerland, Hosmer, and Bofin (2008), with the observed and expected
counts by decile of risk.

``` r

combine_test(m, c("democrat", "republican_other"))
#> Test of combining alternatives (Wald) 
#>                         pair statistic df p.value
#>  democrat & republican_other    318.28  7  <2e-16
term_test(m)
#> Joint tests of each term across equations (Wald) 
#>        term statistic df  p.value
#>     canvass     3.469  3  0.32486
#>       phone     0.557  3  0.90616
#>    democrat    98.148  3  < 2e-16
#>  republican   134.614  3  < 2e-16
#>         age    70.603  3 3.17e-15
#>      female    17.997  3  0.00044
#>   voted2004   129.610  3  < 2e-16
hl_test(m)
#> Fagerland-Hosmer-Bofin goodness-of-fit test (10 groups by the base alternative's fitted probability)
#>   statistic = 69.950 on 24 df, p = 2.22e-06
#>   1 expected count(s) below 5.
```

On these data the Hosmer-Lemeshow-type test rejects clearly: the fitted
probabilities of the multinomial logit are miscalibrated in some deciles
of risk, a sign of misspecification that the IIA tests above do not
locate.

The IIA tests, the score test,
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
and
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)
return data frames of class `"iia_test"` with a `statistic` column, so
their results can be bound and filtered.
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
calibrates those whose null is the fitted model;
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)
and
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md)
test a restriction of it and are refused (for a term, fit the model
without it and use
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)).
[`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md)
returns a list with the observed and expected tables, and
[`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)
a data frame of class `"iia_battery"`.

## References

Cramer, J.S. and Ridder, G. (1991). Pooling states in the multinomial
logit model. *Journal of Econometrics*, 47, 267-272.

Fagerland, M.W., Hosmer, D.W. and Bofin, A.M. (2008). Multinomial
goodness-of-fit tests for logistic regression models. *Statistics in
Medicine*, 27, 4238-4253.

Fok, D. and Paap, R. (2025). New misspecification tests for multinomial
logit models. *Journal of Choice Modelling*, 54, 100531.
<doi:10.1016/j.jocm.2024.100531>.

Hausman, J.A. and McFadden, D. (1984). Specification tests for the
multinomial logit model. *Econometrica*, 52, 1219-1240.

McFadden, D. (1987). Regression-based specification tests for the
multinomial logit model. *Journal of Econometrics*, 34, 63-82.

Small, K.A. and Hsiao, C. (1985). Multinomial logit specification tests.
*International Economic Review*, 26, 619-627.
