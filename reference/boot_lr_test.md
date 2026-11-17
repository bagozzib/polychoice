# Parametric-bootstrap likelihood-ratio (and Vuong) test of nested fits

Calibrates the comparison of a restricted fit with a fit that nests it
by simulation under the restricted model: `B` outcome vectors are drawn
from the restricted fit, both models are refit on each, and the
bootstrap distribution of the likelihood-ratio statistic (and, for an
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) versus
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md)
pair, of the raw Vuong statistic) gives the p-values. This is the
calibrated comparison for restrictions on the boundary of the parameter
space, where neither the chi-square reference of
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
nor the normal reference of
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)
holds: the inflated model against the multinomial logit (\\\phi = 1\\,
where the split coefficients are not identified), the random-intercept
model against the multinomial logit (\\\sigma = 0\\). It applies to the
pairs
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
accepts (a subset of the parameters within one family, with the nest
structure shared for the nested logit, or the multinomial logit with the
same covariates inside a nested logit, an HEV, or a random-intercept
fit), to the multinomial logit inside the inflated multinomial logit,
and to an inflated multinomial logit inside one with more inflated
alternatives; and to two nested
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
fits of the same rankings, where few rankings per coefficient make the
chi-square reference of
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
unreliable. Penalized fits (`firth = TRUE` or a prior) are refused.

## Usage

``` r
boot_lr_test(restricted, full, B = 199L, cores = 1L)
```

## Arguments

- restricted, full:

  Two `"choicemodel"` fits to the same cases, the second nesting the
  first.

- B:

  Number of bootstrap replicates (default 199).

- cores:

  Number of parallel workers.

## Value

A data frame of class `"iia_test"` with one row per statistic (`LR`, and
`Vuong` for an mnl / imnl pair): the observed statistic, `p.value`, and
`p.boot`. For the LR, `p.value` is the chi-square p-value with `df`
degrees of freedom for an interior restriction, the (conservative)
chi-square p-value for \\\sigma = 0\\, the 50:50 chi-square mixture's
for the mnl / imnl pair with one inflated alternative and an
intercept-only split (with
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)'s
rule: `NA` when the inflated fit is degenerate, its outcome coefficients
diverging or at a saddle, as the method line then says), and `NA` for
any other mnl / imnl pair, which has no valid asymptotic reference; for
the Vuong statistic `p.value` is `NA`, and `p.boot` is one-sided in
favour of the full model. The attribute `"boot_statistics"` holds the
replicates; `n_failed` (dropped), `n_nonconverged` (replicates in which
either refit did not converge; kept), and `n_negative` (replicates whose
full refit fell short of the restricted one by more than the embedding's
own gap at a boundary restriction; set to 0 and kept) are attributes
too, and the method line reports each that is not zero. A few
non-converged replicates are common for the inflated model on data
simulated without inflation, where the split equation runs towards
\\\phi = 1\\ and the optimizer stops short of a maximum that is not
attained; `n_negative` guards against an optimizer that ends above its
own start and should be 0. The outcomes are drawn independently case by
case (unit by unit for a random-intercept restricted fit), so the
calibration assumes independent cases; with clustered data the bootstrap
p-value is not calibrated. The split equation must contain an intercept
named `(Intercept)`; a constant column under another name is not
recognized.

## Details

Each replicate refits both models the way the original fits were
obtained (the engine's own starting values and, for the inflated model,
its multistart), and in addition from the original estimates and, for
the full model, from the restricted replicate estimates embedded in the
full model at the restriction (zero for the added coefficients,
\\\lambda = 1\\, scales 1, a standard deviation near zero, or the split
equation at \\\phi = 1 - 10^{-8}\\). The full refit therefore cannot end
below the restricted one, and the replicate statistic is maximized as
thoroughly as the observed one, which is computed the same way. A
replicate takes about as long as the two original fits; `cores` runs
replicates in parallel. The outcomes, and the seeds of the ranking
draws, come from the session's generator before the refits run, so a
result under [`set.seed()`](https://rdrr.io/r/base/Random.html) does not
depend on `cores`.

The outcomes are drawn from the restricted fit's probabilities, case by
case, or for a random-intercept restricted fit unit by unit (one
intercept vector per unit, then the unit's cases). Rankings are drawn as
the order of the restricted fit's utilities plus independent standard
Gumbel errors, cut to each case's observed pattern: the same rank
groups, with the same sizes in the same order (a tie stays a tie at its
positions), and the alternatives below them unranked; the truncation and
the ties are thus treated as independent of the utilities.

Under separation (in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) and
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
fits without unit fixed effects or inflation) the statistic compares the
suprema of the two log-likelihoods, not the values where the optimizer
stopped (see
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md));
replicates whose refit is separated are kept and counted, since
separation can be common under the null when the rankings or choices are
few or the attributes discrete (a replicate in which an alternative is
never chosen, or never ranked above another, is such a replicate). When
the restricted fit is itself separated, the outcomes are drawn from the
point at which its supremum is approached, where the separated
comparisons are near-deterministic (the restricted model's limit), not
from where its optimizer stopped. Outside that scope a replicate in
which an alternative with coefficients of its own is never chosen fails.
Replicates whose refit fails are dropped and counted; the p-value is
\\(1 + k) / (B' + 1)\\ with \\B'\\ the successful replicates and \\k\\
the number of replicate statistics at least as large as the observed
one. Frequency weights are honoured by expanding the design to one row
per case (integer weights only; not with random intercepts, where a
weighted unit stands for several units with their own intercepts).
Fixed-effects fits are refused (the unit constants would be re-estimated
on simulated outcomes that may make them unidentified).

## See also

[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md),
[`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)

Other model comparison:
[`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md),
[`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md),
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)

## Examples

``` r
# \donttest{
set.seed(2)
d <- rmnl(400, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
m0 <- mnl(y ~ x1 + x2, data = d)
m1 <- imnl(y ~ x1 + x2 | 1, data = d)           # a constant captive share; none in the data
#> Warning: The fit is at a boundary: the split equation puts every case in one regime (the inflated model collapses onto the MNL); compare with mnl() by boot_lr_test(), or by lr_test() when the split has only an intercept.
boot_lr_test(m0, m1, B = 19)
#> Parametric bootstrap under the restricted model, 19 replicates; the LR's asymptotic p-value is the 50:50 chi-square mixture (a constant captive share) 
#>                                                   test statistic df p.value
#>   LR (parametric bootstrap under the restricted model)      0.00  1       1
#>  Vuong, raw (full vs restricted; parametric bootstrap)     -0.21 NA      NA
#>  p.boot
#>     1.0
#>     0.9
# }
```
