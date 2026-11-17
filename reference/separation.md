# Separation in multinomial, conditional and rank-ordered logit data

Checks whether the data of an
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) or
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
fit are separated, so that the maximum-likelihood estimates do not
exist. Each observation reports that some alternatives are preferred to
others: the choice to every other available alternative, or, in a
ranking, each alternative of a rank group to each of the next (ties as
groups, the unranked alternatives last). A finite maximum-likelihood
estimate exists if and only if no direction of the coefficients weakly
favours every one of these comparisons and strictly favours at least one
(overlap; Albert and Anderson 1984 for logistic regression, Speckman,
Lee and Sun 2009 for choice models, Sablica, Hornik and Rusch 2026 for
rankings). Along such a direction the log-likelihood rises towards its
supremum and the estimates run off. Overlap is first sought by a
certificate: positive weights on the comparisons that balance them
exactly (Stiemke's lemma), built from the fitted probabilities at the
fit's estimates, where the gradient of the log-likelihood is such a
weighted sum (not for tied rankings). Otherwise the check solves a
linear program (lp_solve, through lpSolveAPI). Both are exact up to
floating-point tolerances.

## Usage

``` r
separation(object)
```

## Arguments

- object:

  A fit from
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) or
  [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md),
  without unit fixed effects or inflation (maximum likelihood or
  penalized: separation is a property of the data and the model, not of
  the fit).

## Value

An object of class `polychoice_separation`: `separated` (logical);
`type` (`"overlap"`, `"quasi-complete"`, or `"complete"`: every
comparison separated); `terms`, the coefficients that move along some
separating direction (a term can have a zero coordinate in the reported
`direction` and still move along another); `direction`, a separating
direction (unit length; along it every comparison some separating
direction can favour is favoured); `cases`, the identifiers of the cases
with a separated comparison; `n_comparisons` (all the comparisons the
data report), `n_separated`, `n_zero` (the comparisons with identical
design rows, up to rounding); and `method`, `"certificate"` or
`"linear program"`.

## Details

The maximum-likelihood fits of
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) and
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
(without unit fixed effects or inflation) run this check and name the
diverging coefficients in a warning; penalized fits (`firth = TRUE` or a
prior) store them in `separation`, and
[`summary()`](https://rdrr.io/r/base/summary.html) marks them `[sep]`.
The other estimators judge divergence by the size of a coefficient and
its Wald statistic.

A comparison between two alternatives with the same design row in a case
(possible without alternative-specific constants), or with rows that
differ only by rounding (each attribute's difference within 64 machine
epsilons of the larger of its two values and the attribute's typical
size, the median of its nonzero absolute values; such as attributes
computed as 0.1 + 0.2 and 0.3, centred or not), is favoured by no
direction. It is counted (`n_zero`), and with one the data are at most
quasi-completely separated, since that choice's probability stays below
one along every direction.

## References

Albert, A. and Anderson, J.A. (1984). On the existence of maximum
likelihood estimates in logistic regression models. *Biometrika*, 71,
1-10. [doi:10.1093/biomet/71.1.1](https://doi.org/10.1093/biomet/71.1.1)

Sablica, L., Hornik, K. and Rusch, T. (2026). The geometry of existence
and uniqueness of maximum likelihood estimation in categorical response
models. arXiv:2610.08100.

Speckman, P.L., Lee, J. and Sun, D. (2009). Existence of the MLE and
propriety of posteriors for a general multinomial choice model.
*Statistica Sinica*, 19, 731-748.

## Examples

``` r
data(defection)
## Cook, Niehaus and Zuhlke's Model 3: no captured combatant joined for
## economic reasons without being indoctrinated (an empty cell)
m3 <- mnl(outcome ~ economic * indoctrination + age_joined + birth_year + school_age + male,
          data = defection, ref = "captured", firth = TRUE)
#> The maximum-likelihood estimates of demobilized:economic, demobilized:economic:indoctrination, switched:economic, switched:economic:indoctrination diverge (separation); their penalized estimates are determined largely by the penalty. Profile penalized-likelihood intervals (confint()) and penalized likelihood-ratio tests are the recommended inference.
separation(m3)
#> Quasi-complete separation: 15 of the 1164 comparisons the data report are separated, in 15 case(s);
#> the maximum-likelihood estimates do not exist.
#> Coefficients moving along some separating direction (the reported direction's coordinates):
#>                demobilized:economic demobilized:economic:indoctrination 
#>                                 0.5                                -0.5 
#>                   switched:economic    switched:economic:indoctrination 
#>                                 0.5                                -0.5 
```
