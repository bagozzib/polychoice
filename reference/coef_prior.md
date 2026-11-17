# Prior distributions for the coefficients of a penalized fit

Specifies independent priors centred at zero for the coefficients of
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) and
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
(`prior = coef_prior(...)`): the fit maximizes the log-likelihood plus
the log prior density, the posterior mode. Under separation the
maximum-likelihood estimates do not exist, and the prior (like Firth's
penalty, `firth = TRUE`, which is Jeffreys' prior) decides how large the
separated estimates are (Rainey 2016). Jeffreys' prior is a default
justified by bias reduction; it depends on the design and is not
independent across coefficients, and can move an estimate against the
direction of the separation (Beiser-McGrath 2022). Independent priors
make the assumption explicit, and fitting under several (a skeptical, an
informative and an enthusiastic scale, Firth's penalty, and Gelman et
al.'s default) shows how much of a separated estimate comes from the
prior.

## Usage

``` r
coef_prior(
  family = c("cauchy", "normal", "t"),
  scale = 2.5,
  df = NULL,
  intercept_scale = 10,
  autoscale = TRUE
)
```

## Arguments

- family:

  `"cauchy"` (default), `"normal"`, or `"t"` (Student t with `df`
  degrees of freedom).

- scale:

  The prior scale of the slopes (default 2.5). A named vector sets the
  scale of the coefficients it names (by their names in
  [`coef()`](https://rdrr.io/r/stats/coef.html)), and an unnamed
  element, if any, the default for the others: a separated coefficient
  can be given a prior of its own (Rainey's partial prior). With
  `autoscale`, a named scale is divided by its column's scale like the
  default. A name the model does not have is ignored with a message, so
  the restricted model of a penalized likelihood-ratio test
  ([`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md))
  can take the full model's prior.

- df:

  Degrees of freedom of the t prior (`family = "t"` only).

- intercept_scale:

  The prior scale of the constants (default 10).

- autoscale:

  If `TRUE` (default), the slopes' scales are relative to their columns
  (see Details); if `FALSE`, they are on the coefficients as they are.
  The column scales use the frequency weights, so a fit whose weights
  total at most 2 (normalized or survey weights, which are not frequency
  weights) is refused under `autoscale = TRUE`: rescale the weights to
  counts, or set `autoscale = FALSE`. A scale named for a coefficient
  the model lacks is ignored with a message, which suggests a near match
  when the name looks like a misspelling.

## Value

An object of class `"coef_prior"`.

## Details

The defaults are Gelman et al.'s (2008): Cauchy priors with scale 2.5
for the slopes and 10 for the constants, on coefficients scaled to the
input. With `autoscale = TRUE` the scale of coefficient \\k\\ is divided
by \\c_k\\: twice the standard deviation of its column when the column
takes more than two values, the difference of its two values when it
takes two (a dummy), as `arm::bayesglm` does. The column scales are
taken over the estimation cases for a case-specific covariate and over
the available alternatives for an alternative-specific (generic) one
(the spread across all available alternatives, between cases as well as
within them; Gelman et al.'s inputs are case-level, and the conditional
logit's within-case spread would be the closer analogue), with the
frequency weights (a weighted fit and the fit to the rows repeated by
their weights get the same prior). The constants (the intercept of the
case-specific part, which gives the alternative-specific constants) take
`intercept_scale`, unscaled. The covariates are not centred, unlike
Gelman et al.'s inputs, so the constants' prior concerns the log odds
where every covariate is zero: centre a covariate whose zero lies far
from its values (a calendar year) for that prior to mean what it says.
With scale 2.5 a normal prior puts 95 percent of its mass on log odds
within about 5 per two standard deviations of a covariate; the default
Cauchy, with its heavier tails, puts 30 percent beyond 5. The prior's
normalizing constants are left out of the penalized log-likelihood (they
do not depend on the coefficients).

The fit reports the log-likelihood at the posterior mode (`loglik`), the
maximized objective (`loglik_penalized`), and the inverse of the
curvature of the log posterior at the mode as the covariance (the
model's information plus the prior's curvature). `arm::bayesglm` differs
for t and Cauchy priors: its approximate EM algorithm reports a related,
always positive definite matrix at its own fixed point, and the two
agree for the normal prior. For a t or Cauchy prior far out in its tail
the exact curvature can fail to be positive definite (a separated
coefficient); the prior's negative curvature is then left out, with a
warning. Profile penalized-likelihood intervals
([`confint()`](https://rdrr.io/r/stats/confint.html)) and penalized
likelihood-ratio tests
([`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md),
[`combine_test()`](https://bagozzib.github.io/polychoice/reference/combine_test.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
between fits with the same prior,
[`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md))
are the recommended inference under separation, as for Firth's penalty.

## References

Beiser-McGrath, L.F. (2022). Separation and rare events. *Political
Science Research and Methods*, 10, 428-437.
[doi:10.1017/psrm.2020.46](https://doi.org/10.1017/psrm.2020.46)

Gelman, A., Jakulin, A., Pittau, M.G. and Su, Y.-S. (2008). A weakly
informative default prior distribution for logistic and other regression
models. *Annals of Applied Statistics*, 2, 1360-1383.
[doi:10.1214/08-AOAS191](https://doi.org/10.1214/08-AOAS191)

Rainey, C. (2016). Dealing with separation in logistic regression
models. *Political Analysis*, 24, 339-355.
[doi:10.1093/pan/mpw014](https://doi.org/10.1093/pan/mpw014)

## See also

[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md),
[`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md)

## Examples

``` r
data(defection)
## the birth year centred, so that the constants' prior concerns typical respondents
defection$birth_c <- defection$birth_year - round(mean(defection$birth_year, na.rm = TRUE))
f <- outcome ~ economic * indoctrination + age_joined + birth_c + school_age + male
mc <- mnl(f, data = defection, ref = "captured", prior = coef_prior())
#> The maximum-likelihood estimates of demobilized:economic, demobilized:economic:indoctrination, switched:economic, switched:economic:indoctrination diverge (separation); their penalized estimates are determined largely by the penalty. Profile penalized-likelihood intervals (confint()) and penalized likelihood-ratio tests are the recommended inference.
mn <- mnl(f, data = defection, ref = "captured", prior = coef_prior("normal", scale = 4.5))
#> The maximum-likelihood estimates of demobilized:economic, demobilized:economic:indoctrination, switched:economic, switched:economic:indoctrination diverge (separation); their penalized estimates are determined largely by the penalty. Profile penalized-likelihood intervals (confint()) and penalized likelihood-ratio tests are the recommended inference.
cbind(cauchy = coef(mc), normal = coef(mn))
#>                                          cauchy       normal
#> demobilized:(Intercept)              0.17715211  0.100825330
#> demobilized:economic                 0.88614176  1.624505407
#> demobilized:indoctrination           0.05477713  0.035840505
#> demobilized:age_joined               0.08275838  0.085627138
#> demobilized:birth_c                 -0.00365362 -0.003613069
#> demobilized:school_age               0.01850121  0.020474830
#> demobilized:male                     0.72002003  0.755107194
#> demobilized:economic:indoctrination -0.56163059 -1.272167890
#> switched:(Intercept)                -1.49336669 -1.571074016
#> switched:economic                    0.67597856  1.289617685
#> switched:indoctrination             -1.94387003 -2.086935508
#> switched:age_joined                  0.02568705  0.025911060
#> switched:birth_c                    -0.08576122 -0.088346355
#> switched:school_age                  0.07800518  0.083534119
#> switched:male                        1.09097791  1.206152446
#> switched:economic:indoctrination     0.74712803  0.248615776
```
