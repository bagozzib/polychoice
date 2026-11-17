# Inflated multinomial logit regression

Fits the baseline-inflated multinomial logit of Bagozzi and Marchetti
(2017) and its generalization to any inflated alternative: a mixture of
a degenerate distribution at the inflated alternative \\k\\ and the
multinomial / conditional logit of
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md), with
a binary split equation deciding which regime a case is in. The
canonical application is vote choice with abstention as the base:
routine nonvoters (the captive regime, always abstain) are mixed with
occasional nonvoters and voters (the choice regime, the multinomial
process), and the split equation separates them.

## Usage

``` r
imnl(
  formula,
  data,
  inflate = NULL,
  id = NULL,
  alt = NULL,
  ref = NULL,
  link = c("logit", "probit"),
  weights = NULL,
  offset = NULL,
  fe = NULL,
  se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"),
  cluster = NULL,
  nboot = 200L,
  cores = 1L,
  start = NULL,
  maxit = 1000L,
  reltol = 1e-10,
  crosscheck = FALSE
)
```

## Arguments

- formula:

  A formula; see the section "Two data layouts" of
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md).

- data:

  A data frame.

- inflate:

  The inflated alternative: a level of the response. Default: `ref`, the
  base alternative (baseline inflation). Several levels give the model
  with a captive share for each (see Details).

- id, alt:

  For the long layout, the names of the case and alternative columns of
  `data`. Leave both `NULL` for the wide layout.

- ref:

  The base (reference) alternative: a level of the response. Default:
  the first level.

- link:

  Link of the split equation: `"logit"` (default; the published model)
  or `"probit"`.

- weights:

  Optional frequency weights: a column name in `data` or a numeric
  vector with one entry per row (constant within a case in the long
  layout).

- offset:

  Optional offset on the utility scale: in the long layout a column name
  or a vector with one entry per row; in the wide layout an \\n \times
  J\\ matrix (columns in level order).

- fe:

  Not available for the inflated model (see Details); an error if given.

- se:

  Standard errors: `"analytic"` (inverse observed information, default),
  `"opg"`, `"robust"`, `"cluster"` (needs `cluster`), `"bootstrap"`
  (nonparametric over cases, or over clusters when `cluster` is given,
  with `nboot` refits from the full-sample estimate), or `"none"`.

- cluster:

  Optional cluster identifier (a column name in `data` or a vector
  aligned to its rows, constant within a case); supplying it selects
  `se = "cluster"` unless `se` is given explicitly.

- nboot:

  Number of bootstrap refits for `se = "bootstrap"` (default 200).

- cores:

  Number of parallel workers for the bootstrap (a PSOCK cluster; default
  1). [`set.seed()`](https://rdrr.io/r/base/Random.html) before the call
  makes the resamples reproducible; serial and parallel runs use
  different random streams.

- start:

  Optional starting values, which replace the multistart: the full
  parameter vector in the order of
  [`coef()`](https://rdrr.io/r/stats/coef.html), the outcome
  coefficients and then the split coefficients (`infl_`), on the same
  scale as [`coef()`](https://rdrr.io/r/stats/coef.html) (`object$theta`
  of an earlier fit).

- maxit, reltol:

  Controls for the BFGS stage passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) (the plain
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) is
  solved by damped Newton steps first, with BFGS as the fallback); the
  Newton polish that follows runs up to 25 further iterations with the
  analytic Hessian of the logit family or central differences of the
  analytic gradient otherwise.

- crosscheck:

  With one inflated alternative, also search under the other split link
  and polish this link's fit from its rescaled estimates (see Details);
  about doubles the fitting time. `FALSE` (default) skips it. Ignored
  with several inflated alternatives (which take the logit link only)
  and when `start` is given. The parametric-bootstrap refits of
  [`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
  and
  [`boot_pvalue()`](https://bagozzib.github.io/polychoice/reference/boot_pvalue.md)
  search as the fit did.

## Value

An object of class `c("imnl", "choicemodel")`; as
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) plus
`phi` (the fitted \\P(s_i = 1)\\), `prob_outcome` (the \\n \times J\\
outcome-stage probabilities \\P\_{ij}\\), `inflate`, and `link`; with
several inflated alternatives also `regime`, the matrix of the regime
probabilities.

## Details

With \\s_i \in \\0, 1\\\\, \\P(s_i = 1 \mid z_i) = \phi_i =
F(z_i'\delta)\\ (\\F\\ logistic by default, or standard normal), and the
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
probabilities \\P\_{ij}\\, \$\$P(y_i = k) = (1 - \phi_i) + \phi_i
P\_{ik}, \qquad P(y_i = j) = \phi_i P\_{ij}, \\ j \ne k.\$\$ In the
discrete-choice literature, a model with a logit split equation is the
parametrized logit captivity model of Swait and Ben-Akiva (1987, eq. 3)
with captivity to \\k\\ only. That model generalizes the dogit of Gaudry
and Dagenais (1979), which adds to the multinomial logit a share of
decision makers captive to each alternative, by making the captivity
odds \\\theta\_{ik} = (1 - \phi_i)/\phi_i\\ a function of covariates;
here \\\theta\_{ik} = e^{-z_i'\delta}\\, so the split coefficients are
the negatives of Swait and Ben-Akiva's captivity coefficients.

Several inflated alternatives (`inflate = c("abstain", "A")`) give the
parametrized logit captivity model with captive shares for each of them
(Swait and Ben-Akiva 1987), which suits, say, routine nonvoters together
with partisans captive to their party. The regimes (choosing, or captive
to one of the inflated alternatives \\c \in C\\) follow a multinomial
logit in \\z\\: with \\\eta\_{ic} = z_i'\delta_c\\, \$\$\pi\_{i0} =
\frac{1}{1 + \sum\_{c} e^{-\eta\_{ic}}}, \qquad \pi\_{ic} =
\frac{e^{-\eta\_{ic}}}{1 + \sum\_{c} e^{-\eta\_{ic}}}, \qquad P(y_i = j)
= \pi\_{i0} P\_{ij} + 1\\j \in C\\\\ \pi\_{ij},\$\$ so \\\eta\_{ic}\\ is
the log-odds of choosing against captivity to \\c\\: a positive
coefficient makes captivity to \\c\\ less likely, the reading of the
single split equation, which is the case \\\|C\| = 1\\. Every inflated
alternative gets its own coefficients on the same \\z\\ (named
`infl_<alternative>:<term>`; Swait and Ben-Akiva allow each alternative
its own captivity covariates, while here every inflated alternative's
equation has the same \\z\\); an inflated alternative outside a case's
choice set is not one of its regimes. Several inflated alternatives take
the logit split only. `predict(type = "regime")` then returns the matrix
of the regime probabilities and `type = "posterior"` their posteriors
given the observed choice; `phi` is the probability of the choice
regime. A captive share that collapses to zero for every case is flagged
`boundary`; the model with fewer inflated alternatives is nested at that
boundary, so
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
compares the two (as it compares the inflated model with
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)). Each
captive share is identified mainly through the mixture's form, as the
single one is; an exclusion restriction helps. In the package's
simulation (four alternatives, two of them captive, 200 samples per
design), with a split covariate outside the outcome equations the 95
percent intervals covered 0.92-0.975 at n = 5,000; at n = 2,000 two
cells missed, the constant of `c`'s split equation (coverage 0.915) and
its slope (a bias of -3.1 Monte Carlo standard errors). With the split
on an outcome covariate instead, 3 of 200 fits at n = 5,000 had
diverging choice coefficients (flagged `boundary`) and the intervals of
that covariate's outcome coefficients covered 0.89-0.90; at n = 2,000,
17 of 200 fits were at a boundary and those intervals covered 0.80-0.81.

The MNL is the boundary case \\\phi \to 1\\, where the split
coefficients are not identified, so neither the chi-square reference of
a likelihood-ratio test nor the normal reference of
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)
holds for the comparison with
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md):
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
calibrates both statistics by a parametric bootstrap under the MNL (the
information criteria are a descriptive complement). The exception is a
single inflated alternative with an intercept-only split (`y ~ x | 1`, a
constant captive share), whose comparison with the MNL has a regular
boundary:
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
and `captivity_test(type = "lr")` refer it to the 50:50 chi-square
mixture. The split-equation coefficients are reported with the prefix
`infl_`. Identification rests on the mixture structure, as for the
zero-inflated count and inflated ordered models: an exclusion
restriction (a variable in `z` not in `x`, or vice versa) is advisable
but not required (Bagozzi and Marchetti 2017). The split equation models
\\P(s_i = 1)\\, the probability of being in the CHOICE regime: a
positive coefficient means the covariate makes a case less likely to be
a routine chooser of the inflated alternative.
`predict(type = "regime")` returns \\\phi_i\\ (and 1 for a case whose
choice set lacks \\k\\: the captive regime is impossible there, so such
cases count as in the choice regime; the boundary check and
[`summary()`](https://rdrr.io/r/base/summary.html)'s mean of \\\phi_i\\
leave them out), `type = "posterior"` the regime probability given the
observed choice (\\\phi_i P\_{ik} / P^\*\_{ik}\\ for a case that chose
\\k\\, 1 otherwise), which classifies the observed choosers of \\k\\
into routine and occasional ones.

The formula has the split equation as its **last** part: in the wide
layout `y ~ x1 + x2 | z1 + z2`, in the long layout
`choice ~ w1 | x1 + x2 | z1 + z2` (write `| 1 |` for an empty
case-specific part). Everything else – layouts, weights, offsets,
standard errors – is as in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md). Unit
fixed effects are not available: with a constant per unit in the outcome
equations the level of \\\phi\\ is identified only through the
functional form of the mixture, units that always chose the inflated
alternative (the routine choosers the model is designed for) drop out,
units with few such choices have them all put in the captive regime and
their constants run to infinity, and units that never chose it pull
\\\phi\\ towards 1. In the package's simulations the split coefficients
came out 25 to 42 percent too large at 12 to 20 cases per unit (and two
to three times too large at 4), with a third or more of the unit
constants at infinity, so `imnl()` refuses `fe`; for panels, fit without
it with `cluster` set to the unit, or use
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) with
`fe` for the outcome equations.

The inflated likelihood can have several local maxima (it has two on
each of `pa2006` and `pres2004`, and many where the split equation can
separate a few cases into the inflated regime). Estimation therefore
runs a grid of fourteen starts (the
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
coefficients, or an MNL that down-weights the inflated-alternative
cases, combined with a binary fit of \\1\\y \ne k\\\\ on \\z\\ at full,
halved, doubled, tripled, or zero slopes and several shifts of the split
intercept), screens them by quasi-Newton, and polishes every screened
start to convergence with the analytic Hessian. The incumbent is then
challenged with its split slopes halved, enlarged by half, and doubled,
again from each new incumbent (at most three rounds). With one inflated
alternative, `crosscheck = TRUE` adds before the challengers the same
search under the other split link, whose estimates, the split
coefficients rescaled by 1.6 (the ratio of the logistic and normal
densities at zero), are polished under this link as one more candidate:
the two links' surfaces differ, and on the vote-choice data of Glasgow
(2001) the logit split's own search ends 4.9 log-likelihood units below
the maximum it reaches from the probit split's estimates. It about
doubles the fitting time, so it is not the default; use it when the
message reports several maxima. The highest maximum found is reported;
it need not be the global one, so probe other basins with `start` too.
The challengers start from the best polished fit only, so where the
split coefficients diverge (the likelihood then has many flat plateaus,
each a different hard classification of a few cases into the regimes)
the plateau reported depends on the search path; such a fit is flagged
`boundary`, and its split estimates are not maxima to interpret. The
log-likelihoods of the polished starts (converged or not), of the
cross-check and of the challengers that converged or improved, and the
maximum reached by a saddle step, are in `object$start_logliks`, those
of all screened starts in `object$screen_logliks`, the distinct points
where these searches ended, each classified as a local maximum, a
boundary (diverging coefficients or a split equation at one regime), a
saddle, or not converged, in `object$search` (with their estimates in
`object$search_theta`), the number of distinct local maxima among them
in `object$n_modes` (see the section "Fit-quality flags" of
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)), and
with `crosscheck = TRUE` the cross-check in `object$crosscheck` (the
other link's maximum, the maximum this link reached from it, and whether
that improved on this link's own search); a message reports
multimodality. A fit whose split equation puts every case in one regime
has collapsed onto the MNL; it is flagged `boundary` with a warning.

The likelihood can also rise without bound along a direction of the
outcome coefficients, to a supremum at infinity: separation within the
mixture. The outcome stage can then predict some cases of the inflated
alternative with certainty while the split equation puts its other cases
in the captive regime, which overlapping data under the multinomial
logit do not prevent
([`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md)
checks the multinomial logit and the rank-ordered logit only). Strictly,
the maximum-likelihood estimate then does not exist, and the reported
fit is the highest end the search reaches. On `pres2004` a search from
200 starts finds two interior maxima, the published fit (-545.279) and
-550.259, and from two starts points where the Kerry and Bush
coefficients diverge together, the higher at -542.851: 21 abstainers are
predicted with certainty by the outcome stage and the other 84 are
captive. The multistart does not reach that region; a fit started there
is flagged (diverging coefficients, `boundary`, a singular information
matrix). On `pa2006` the same search finds the multistart's two maxima
(-2069.971 and -2070.990) and no higher point.

## References

Bagozzi, B.E. and Marchetti, K. (2017). Distinguishing occasional
abstention from routine indifference in models of vote choice.
*Political Science Research and Methods*, 5(2), 277-294.
[doi:10.1017/psrm.2015.42](https://doi.org/10.1017/psrm.2015.42)

Gaudry, M.J.I. and Dagenais, M.G. (1979). The dogit model.
*Transportation Research Part B*, 13, 105-111.
[doi:10.1016/0191-2615(79)90028-6](https://doi.org/10.1016/0191-2615%2879%2990028-6)

Swait, J. and Ben-Akiva, M. (1987). Empirical test of a constrained
choice discrete model: mode choice in Sao Paulo, Brazil. *Transportation
Research Part B*, 21, 103-115.
[doi:10.1016/0191-2615(87)90010-5](https://doi.org/10.1016/0191-2615%2887%2990010-5)

## See also

[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`predict.choicemodel()`](https://bagozzib.github.io/polychoice/reference/predict.choicemodel.md),
[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md),
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)

Other estimators:
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md),
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md),
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)

## Examples

``` r
set.seed(2)
d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
          delta = c(0.5, 1.2), inflate = "a")
m <- imnl(y ~ x1 + x2 | z1, data = d, inflate = "a")
summary(m)
#> Inflated multinomial logit (inflated alternative: a, logit split) 
#> Call: imnl(formula = y ~ x1 + x2 | z1, data = d, inflate = "a")
#> Alternatives: a, b, c  (base: a)
#> 
#> Outcome equations (alternative:term):
#>               Estimate Std. Error z value Pr(>|z|)    
#> b:(Intercept)   0.1103     0.2484   0.444 0.657001    
#> b:x1           -0.9214     0.1858  -4.959 7.10e-07 ***
#> b:x2            0.8229     0.1575   5.224 1.75e-07 ***
#> c:(Intercept)  -0.8921     0.2630  -3.392 0.000693 ***
#> c:x1            0.7548     0.1765   4.277 1.90e-05 ***
#> c:x2           -0.5546     0.1701  -3.259 0.001117 ** 
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Split equation (P(choice regime) = F(z'delta)):
#>                  Estimate Std. Error z value Pr(>|z|)    
#> infl_(Intercept)   0.7931     0.3404   2.330   0.0198 *  
#> infl_z1            1.4749     0.2994   4.927 8.37e-07 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Standard errors: analytic. Log-likelihood -625.233 on 8 parameters, 800 cases. AIC 1266.47, BIC 1303.94.
#> Inflated alternative 'a': observed share 0.559; mean fitted P(choice regime) 0.660.
## the inflated model nests the multinomial logit at a boundary (phi = 1): compare the
## two by information criteria here, or by boot_lr_test(m0, m)
m0 <- mnl(y ~ x1 + x2, data = d)
compare_models(m0, m)
#>  model                                                              type
#>     m0                                                 Multinomial logit
#>      m Inflated multinomial logit (inflated alternative: a, logit split)
#>  inflate   n df   loglik      AIC      BIC mcfadden converged boundary
#>          800  6 -680.707 1373.414 1401.522    0.121      TRUE    FALSE
#>        a 800  8 -625.233 1266.467 1303.944    0.193      TRUE    FALSE
#>  ill_conditioned n_modes
#>            FALSE       1
#>            FALSE       1
```
