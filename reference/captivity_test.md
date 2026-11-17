# Tests for captive (inflated) choosers: score and likelihood ratio

Tests the multinomial logit against the dogit model of Gaudry and
Dagenais (1979), in which a share of the decision makers is captive to
an alternative: \\P^\*\_{ij} = (P\_{ij} + \theta_j) / (1 + \sum_m
\theta_m)\\, \\\theta_j \ge 0\\, with the sum over the case's choice
set. The multinomial logit is \\\theta = 0\\. The joint statistic is the
Lagrange multiplier statistic of Tse (1987, eqs. 4-5), with the expected
information, extended here to choice sets that vary across cases and to
alternative-specific covariates; the one-sided z statistic for each
alternative is this package's addition. The score form needs only the
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) fit;
the likelihood-ratio form (`type = "lr"`, below) refits. Both need at
least three alternatives. The score of case \\i\\ in \\\theta_j\\ at
\\\theta = 0\\ is \\1\\y_i = j\\/P\_{ij} - 1\\: it is positive when the
cases that chose \\j\\ are ones the multinomial logit gave a low
probability to, which is what captive choosers produce. The expected
information at the multinomial logit estimates, with the coefficients
partialled out, gives a one-sided z statistic for each alternative
(captivity cannot be negative) and a joint chi-square statistic over the
alternatives tested. A rejection for the base alternative of a
vote-choice model is evidence of routine choosers of it, the case the
inflated multinomial logit
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md)
models.

## Usage

``` r
captivity_test(object, alternatives = NULL, type = c("score", "lr"))
```

## Arguments

- object:

  A fitted
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
  model. A fit with diverging coefficients (separation) draws a warning:
  the statistic is then evaluated where the optimizer stopped, at
  estimates that do not exist.

- alternatives:

  The alternatives to test (labels); default all.

- type:

  `"score"` (default; computed from the multinomial logit fit alone) or
  `"lr"` (the likelihood ratio against the constant-share inflated
  multinomial logit, one refit per alternative).

## Value

A data frame of class `"iia_test"`. With `type = "score"`, one row per
alternative (`statistic` the one-sided z, `p.value` its upper-tail
p-value; `NA` for an alternative no case chose, which carries no
evidence on captivity to it and is left out of the joint test) and a
last row for the joint test (`statistic` the chi-square, computed on the
correlation scale of the scores, `df`, `p.value`). With `type = "lr"`,
one row per alternative (`statistic` the likelihood ratio, `df` 1,
`p.value` from the mixture, `share` the fitted captive share), and the
attribute `fits`, one row per alternative:

- `alternative`, `loglik` (the inflated fit's), `captive_share`;

- `converged`, `boundary`, `n_modes` (the search's distinct maxima);

- `top_share` (above);

- `warnings` (the inflated fit's, joined), `note`;

- `status`:

  - `"logit"`: the fit reached the multinomial logit; statistic 0,
    p-value 1, and share 0, or NA when the share is not identified;

  - `"interior"`: the mixture p-value;

  - `"degenerate"`: diverging outcome coefficients or a saddle; the
    statistic is kept and the p-value is NA;

  - `"failed"`: the fit failed or was refused; statistic and p-value NA.
    The attribute `note` collects the fits' notes.

## Details

The constant captivity parameters make the multinomial logit an ordinary
boundary point of the dogit (\\\theta = 0\\), so the score test is
regular there, unlike the comparison with
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
whose split coefficients are not identified at \\\phi = 1\\; it can
detect covariate-dependent captivity through its average (but see the
loss of power below). Tse refers the joint statistic to the chi-square
with \\J\\ degrees of freedom, reads the boundary null as a test for
arbitrarily small \\\theta\\, and found its size at 0.05 acceptable in
his simulations (three alternatives, 50 and 300 cases). When the
multinomial logit's fitted probabilities of the tested alternatives do
not vary across cases (no case-specific and no alternative-specific
covariates), captivity is not identified and the statistic is undefined.

The weights \\1/P\_{ij}\\ make the statistic heavy-tailed when some
fitted probabilities are small: a single case that chose an alternative
the multinomial logit gave a probability of 0.001 adds about 1000 to its
score. In the package's simulations under the multinomial logit (three
and four alternatives, 500 and 2,000 cases, smallest fitted
probabilities around 0.0005) the one-sided tests rejected in 2.8 to 6.3
percent of samples at the 0.05 level and the joint test in 5.7 to 7.3
percent, but at the 0.01 level they rejected in 1.3 to 3.7 percent: for
small significance levels calibrate by simulation with
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md).
Against a baseline captive share of about 15 percent the test of that
alternative rejected in 47 percent of samples of 500 cases and 94
percent of samples of 2,000. The expected information is used: the outer
product of the scores, which would also allow clustering, sized far
worse in the same simulations (one-sided rejection at 0.05 of 0.3 to 2.3
percent, joint 16 to 23 percent).

The same weights cost power. A case whose fitted probability of an
alternative is small and that did not choose it adds -1 to the score but
\\1/P\_{ij} - 1\\ to its variance. On
[pa2006](https://bagozzib.github.io/polychoice/reference/pa2006.md)
(Replication I's multinomial logit) the one-sided z of abstention is
0.50, while
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md) with
a constant split, the test's own alternative, raises the log-likelihood
by 10.5 on one parameter, at a captive share of 0.20. The score is 37.3
(a numerical derivative of the log-likelihood agrees), but its expected
information (before the coefficients are partialled out) is 10,186, half
of it from one respondent whose fitted probability of abstaining is
0.0002 and who did not abstain; and along the captive share the
log-likelihood is convex at zero, so the gain lies away from the null,
where a local test does not look. A small statistic is therefore no
evidence against captivity. (On
[fluvax](https://bagozzib.github.io/polychoice/reference/fluvax.md) the
two agree: z = 0.63 for "don't know", and the constant-share fit gains
0.41.)

`type = "lr"` measures that gain. For each alternative it fits the
inflated multinomial logit with that alternative inflated and an
intercept-only split equation, which is the dogit with a constant
captive share of it. The fit uses the multinomial logit's own cases,
weights, offsets and choice sets,
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md)'s
search, and a start at the multinomial logit itself. It reports the
likelihood ratio against the multinomial logit and the captive share.
The multinomial logit is a captive share of zero, a boundary point of
the share's space. The point is a regular one when the other parameters
are identified there and the share's score is not a combination of
theirs, that is, when the fitted probabilities of the alternative vary
in a way its constants and slopes cannot absorb. The reference is then
the 50:50 mixture of chi-square(0) and chi-square(1) (Self and Liang
1987), whose p-value is half the chi-square(1) tail (and 1 for a
statistic of 0).

- On `pa2006` it gives 20.95 for abstention (p = 2.4e-6, a captive share
  of 0.196).

- There, where one respondent carries half the information, a parametric
  bootstrap from the multinomial logit (199 replicates) put 0.61 of the
  statistics at zero (the mixture puts half there) and exceeded the
  mixture's 0.10, 0.05 and 0.01 critical values in 0.101, 0.040 and
  0.010 of the replicates: the upper tail agreed with the mixture within
  Monte Carlo error, and the reference was, if anything, conservative.

- Where few cases dominate the information (the attribute `fits` reports
  `top_share`, the largest row's share; with frequency weights a row of
  weight w stands for w cases),
  [`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md),
  which accepts this pair, calibrates it.

- A saturated outcome design (a single factor, say) absorbs a constant
  share into the alternative's constants: the ratio is then 0 and the
  share is reported as not identified.

- A fit at which outcome coefficients diverge, or that stops at a
  saddle, keeps its statistic with a p-value of NA and a note, since the
  reference does not hold there. The fits' warnings are in
  `attr(, "fits")`.

The form refits the inflated model once per alternative (about 20
seconds each on `pa2006`) and has no joint row: with several captive
shares on the boundary the reference is a chi-square mixture whose
weights depend on the information. Each row is a separate test, so with
every alternative tested (the default) one can reject by chance. The
comparison with an inflated model whose split equation has covariates is
not of this kind;
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
calibrates it. Like the other likelihood-ratio and score statistics of
the package, both forms assume independent cases and ignore a fit's
`cluster`.

## References

Gaudry, M.J.I. and Dagenais, M.G. (1979). The dogit model.
*Transportation Research Part B*, 13, 105-111.
[doi:10.1016/0191-2615(79)90028-6](https://doi.org/10.1016/0191-2615%2879%2990028-6)

Self, S.G. and Liang, K.-Y. (1987). Asymptotic properties of maximum
likelihood estimators and likelihood ratio tests under nonstandard
conditions. *Journal of the American Statistical Association*, 82,
605-610.
[doi:10.1080/01621459.1987.10478472](https://doi.org/10.1080/01621459.1987.10478472)

Tse, Y.K. (1987). A diagnostic test for the multinomial logit model.
*Journal of Business and Economic Statistics*, 5, 283-286.
[doi:10.1080/07350015.1987.10509587](https://doi.org/10.1080/07350015.1987.10509587)

## See also

[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md),
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)

Other specification tests:
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
[`hl_test()`](https://bagozzib.github.io/polychoice/reference/hl_test.md),
[`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md),
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md)

## Examples

``` r
set.seed(6)
d <- rmnl(2000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
          delta = c(1.5, 0), inflate = "a")
m <- mnl(y ~ x1 + x2, data = d)
captivity_test(m)
#> Score test of the multinomial logit against captive choosers (dogit; one-sided z per alternative, joint chi-square) 
#>        test statistic df  p.value
#>  captive: a     4.512 NA 3.22e-06
#>  captive: b    -2.395 NA 0.991691
#>  captive: c    -1.056 NA 0.854540
#>       joint    20.951  3 0.000108
# \donttest{
## the likelihood ratio against a constant captive share of a (one refit)
captivity_test(m, alternatives = "a", type = "lr")
#> Likelihood-ratio test of the multinomial logit against a constant captive share of each alternative (the inflated multinomial logit with an intercept-only split; p-value from the 50:50 mixture of chi-square(0) and chi-square(1)) 
#>        test statistic df  p.value share
#>  captive: a    22.032  1 1.34e-06 0.212
# }
```
