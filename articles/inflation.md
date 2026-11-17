# The inflated multinomial logit

## Two kinds of abstainers

In a model of vote choice the base alternative, abstention, holds two
kinds of people: *routine* nonvoters, who would not have voted whatever
the campaign did, and *occasional* nonvoters, who weighed the candidates
and stayed home this time. A standard multinomial logit treats them as
one group and, when the mixture is real, misstates how the covariates
move people between abstention and the candidates. The baseline-inflated
multinomial logit (BIMNL) of Bagozzi and Marchetti (2017) adds a split
equation that sorts cases into a regime where the base alternative is
observed with certainty and a regime where the multinomial choice
operates:

``` math
P(s_i = 1 \mid z_i) = \phi_i = F(z_i'\delta), \qquad
P(y_i = b) = (1 - \phi_i) + \phi_i P_{ib}, \qquad
P(y_i = j) = \phi_i P_{ij} \ (j \ne b),
```

with $`F`$ logistic (the default) or standard normal and $`P_{ij}`$ the
multinomial logit probabilities. The same structure fits any alternative
that mixes two populations: `inflate` names it (default: the base).

In the discrete-choice literature, with the logistic split equation this
is the parametrized logit captivity model of Swait and Ben-Akiva (1987)
with captivity to the inflated alternative only. Their model generalizes
the dogit of Gaudry and Dagenais (1979), which adds to the multinomial
logit a share of decision makers who are captive to an alternative, by
letting the captivity odds depend on covariates. Here the captives of
the inflated alternative are the routine nonvoters, their share is
$`1 - \phi_i`$, and their odds are
$`\theta_i = (1 - \phi_i)/\phi_i = e^{-z_i'\delta}`$.

## Fitting the model

The formula has two parts: the outcome equation and, after `|`, the
split equation.
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md) fits
by maximum likelihood with a multistart because the mixture likelihood
can have several local maxima (see below).

``` r

library(polychoice)
data(pa2006)
f <- vote ~ canvass + phone + democrat + republican + canvass:democrat + canvass:republican +
  phone:democrat + phone:republican + household + age + female + voted2004 + district156
m0 <- mnl(f, data = pa2006)
m1 <- imnl(vote ~ canvass + phone + democrat + republican + canvass:democrat +
             canvass:republican + phone:democrat + phone:republican + household + age +
             female + voted2004 + district156 |
             household + age + female + voted2004 + district156, data = pa2006)
summary(m1)
#> Inflated multinomial logit (inflated alternative: abstain, logit split) 
#> Call: imnl(formula = vote ~ canvass + phone + democrat + republican + 
#>     canvass:democrat + canvass:republican + phone:democrat + 
#>     phone:republican + household + age + female + voted2004 + 
#>     district156 | household + age + female + voted2004 + district156, 
#>     data = pa2006)
#> Alternatives: abstain, democrat, republican_other, nonresponse  (base: abstain)
#> 
#> Outcome equations (alternative:term):
#>                                     Estimate Std. Error z value Pr(>|z|)    
#> democrat:(Intercept)                -5.37153    0.73917  -7.267 3.68e-13 ***
#> democrat:canvass                     0.67773    0.51288   1.321 0.186368    
#> democrat:phone                       1.04057    0.50273   2.070 0.038467 *  
#> democrat:democrat                    1.73797    0.58696   2.961 0.003067 ** 
#> democrat:republican                  0.09916    0.59988   0.165 0.868712    
#> democrat:household                   0.17930    0.11474   1.563 0.118134    
#> democrat:age                         0.09262    0.01580   5.864 4.53e-09 ***
#> democrat:female                     -0.77904    0.21264  -3.664 0.000249 ***
#> democrat:voted2004                   0.88497    0.61392   1.442 0.149443    
#> democrat:district156                -0.08333    0.22233  -0.375 0.707800    
#> democrat:canvass:democrat           -0.34538    0.69987  -0.493 0.621664    
#> democrat:canvass:republican         -1.34547    0.73304  -1.835 0.066437 .  
#> democrat:phone:democrat             -1.00262    0.65287  -1.536 0.124607    
#> democrat:phone:republican           -1.66623    0.69048  -2.413 0.015815 *  
#> republican_other:(Intercept)        -6.01545    0.77389  -7.773 7.66e-15 ***
#> republican_other:canvass             1.56428    0.55710   2.808 0.004987 ** 
#> republican_other:phone               1.31136    0.55440   2.365 0.018012 *  
#> republican_other:democrat            0.66029    0.70107   0.942 0.346283    
#> republican_other:republican          1.99997    0.59882   3.340 0.000838 ***
#> republican_other:household           0.20776    0.11376   1.826 0.067816 .  
#> republican_other:age                 0.08320    0.01549   5.372 7.78e-08 ***
#> republican_other:female             -0.43432    0.20126  -2.158 0.030930 *  
#> republican_other:voted2004           1.20780    0.62310   1.938 0.052578 .  
#> republican_other:district156        -0.19960    0.22044  -0.905 0.365238    
#> republican_other:canvass:democrat   -1.46272    0.82586  -1.771 0.076537 .  
#> republican_other:canvass:republican -1.52719    0.69218  -2.206 0.027359 *  
#> republican_other:phone:democrat     -1.08746    0.77149  -1.410 0.158670    
#> republican_other:phone:republican   -1.54665    0.67476  -2.292 0.021897 *  
#> nonresponse:(Intercept)             -6.45231    1.06892  -6.036 1.58e-09 ***
#> nonresponse:canvass                  0.37835    0.56656   0.668 0.504266    
#> nonresponse:phone                    0.45892    0.54428   0.843 0.399133    
#> nonresponse:democrat                -1.06053    0.91882  -1.154 0.248406    
#> nonresponse:republican              -0.94327    0.80406  -1.173 0.240743    
#> nonresponse:household                0.09617    0.15462   0.622 0.533945    
#> nonresponse:age                      0.07678    0.01609   4.772 1.83e-06 ***
#> nonresponse:female                   0.13555    0.27677   0.490 0.624299    
#> nonresponse:voted2004                2.15289    0.91502   2.353 0.018631 *  
#> nonresponse:district156             -0.46994    0.27063  -1.736 0.082479 .  
#> nonresponse:canvass:democrat         0.16538    1.08378   0.153 0.878718    
#> nonresponse:canvass:republican      -0.06202    0.94506  -0.066 0.947675    
#> nonresponse:phone:democrat           0.37673    0.99603   0.378 0.705258    
#> nonresponse:phone:republican        -0.18810    0.89386  -0.210 0.833331    
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Split equation (P(choice regime) = F(z'delta)):
#>                  Estimate Std. Error z value Pr(>|z|)    
#> infl_(Intercept)  1.09988    1.35398   0.812  0.41660    
#> infl_household    0.08104    0.12321   0.658  0.51068    
#> infl_age         -0.03697    0.01651  -2.240  0.02509 *  
#> infl_female       0.20719    0.23931   0.866  0.38660    
#> infl_voted2004    2.41018    0.39064   6.170 6.84e-10 ***
#> infl_district156  0.64112    0.21653   2.961  0.00307 ** 
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Standard errors: analytic. Log-likelihood -2069.971 on 48 parameters, 1998 cases. AIC 4235.94, BIC 4504.74.
#> Inflated alternative 'abstain': observed share 0.366; mean fitted P(choice regime) 0.828.
#> Note: 2 distinct local maxima among the ends of the search from 14 starts; the highest is reported (see object$search).
```

The summary shows the outcome equations (`alternative:term`), the split
equation (`infl_term`, the coefficients of $`z_i'\delta`$ on the scale
of the link), the observed share of the inflated alternative, and the
mean fitted probability of the choice regime. The sign convention of the
split equation: $`\delta`$ models $`P(s_i = 1)`$, the probability of
being in the *choice* regime, so the negative `infl_age` printed above
means that older respondents are less likely to be in the choice regime,
that is, more likely to be routine abstainers. The estimates agree with
the published Table A.3 of Bagozzi and Marchetti (2017) to within 0.002
in the coefficients and standard errors that the package’s tests compare
(for example 1.564 (0.557) for the Republican/other canvassing
coefficient against the published 1.566 (0.558); the package’s tests pin
a selection of them at tolerances of one to ten percent), and the
published log-likelihoods of both models to their two printed decimals.

## Is the split equation needed?

The plain multinomial logit is the special case $`\phi_i = 1`$, a
boundary of the parameter space, so the likelihood-ratio test does not
have its usual distribution (and the split-equation slopes are not
identified under the null, so no simple correction exists; a split with
an intercept only is the exception, below under “Tests for captive
choosers”). The literature compares the two models with the Vuong (1989)
test, the information criteria, and the classification tables:
[`compare_models()`](https://bagozzib.github.io/polychoice/reference/compare_models.md)
lines up the information criteria,
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)
and
[`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md)
give the other two.

``` r

compare_models(MNL = m0, BIMNL = m1)
#>  model                                                                    type
#>    MNL                                                       Multinomial logit
#>  BIMNL Inflated multinomial logit (inflated alternative: abstain, logit split)
#>  inflate    n df    loglik      AIC      BIC mcfadden converged boundary
#>          1998 42 -2110.593 4305.187 4540.382    0.160      TRUE    FALSE
#>  abstain 1998 48 -2069.971 4235.942 4504.737    0.176      TRUE    FALSE
#>  ill_conditioned n_modes
#>            FALSE       1
#>            FALSE       2
vuong(m1, m0)
#> Vuong test:  m1 = m1 [ Inflated multinomial logit (inflated alternative: abstain, logit split) ]
#>              m2 = m0 [ Multinomial logit ]
#> positive statistics favor m1
#> 
#>  correction statistic p_m1_better p_m2_better p_two_sided
#>         raw     4.479       0.000       1.000       0.000
#>         AIC     3.817       0.000       1.000       0.000
#>         BIC     1.965       0.025       0.975       0.049
#> Note: The inflated multinomial logit and the multinomial logit are nested at the boundary phi = 1, where the normal reference does not hold for testing the uninflated model; boot_lr_test() calibrates the statistic.
```

A positive Vuong statistic favours the first model; the three rows apply
no correction, the AIC correction, and the BIC correction for the number
of parameters. The normal reference of the Vuong test does not hold for
models nested at a boundary (Wilson 2015), so the calibrated comparison
is the parametric bootstrap under the multinomial logit,
`boot_lr_test(m0, m1)`, which refits both models on data simulated from
the MNL and reports bootstrap p-values for the likelihood-ratio and
Vuong statistics. It is slow: every replicate refits the inflated model
with its whole search (about 40 seconds on these data), so the default
`B = 199` takes more than two hours on one core, and `cores` spreads the
replicates over workers; it is not run here.

## Quantities of interest by stage

Predictions can be made for the observed outcome (`type = "prob"`), for
the multinomial stage alone (`type = "prob_outcome"`, the probabilities
among cases in the choice regime), and for the regime itself
(`type = "regime"`, the probability of being in the choice regime). The
labels `infl_` and `"inflation"` (the `block` of
[`summary()`](https://rdrr.io/r/base/summary.html), the `component` of
`tidy()`, and `stage = "inflation"` in
[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md))
all mean the split equation.

``` r

head(predict(m1, type = "regime"), 3)
#>         1         2         3 
#> 0.9751583 0.9475315 0.8309676
head(predict(m1, type = "prob_outcome"), 3)
#>      abstain   democrat republican_other nonresponse
#> 1 0.75254430 0.02545574        0.1955426  0.02645737
#> 2 0.34433867 0.08146521        0.5134851  0.06071103
#> 3 0.03832782 0.16109000        0.7232341  0.07734805
```

[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md)
has the same `type` switch (the change in the observed probabilities, in
the outcome-stage probabilities, or in the regime probability) and a
`stage` argument that says through which equations the covariate is
allowed to move: both (default), the outcome equation only, or the split
equation only. A covariate that appears in both equations, such as `age`
here, can thus be decomposed.

``` r

first_difference(m1, "canvass", from = 0, to = 1)
#> First differences: canvass from 0 to 1 at the covariate profile (delta intervals at 95%)
#>       alternative  from    to   diff  lower upper
#>           abstain 0.296 0.261 -0.035 -0.096 0.025
#>          democrat 0.335 0.273 -0.062 -0.147 0.024
#>  republican_other 0.304 0.389  0.085 -0.009 0.180
#>       nonresponse 0.066 0.077  0.011 -0.042 0.065
first_difference(m1, "age", from = 30, to = 60, type = "regime")
#> First differences: age from 30 to 60 at the covariate profile (delta intervals at 95%)
#> Quantity: P(choice regime) 
#>       alternative  from    to   diff  lower  upper
#>  P(choice regime) 0.937 0.832 -0.106 -0.129 -0.083
first_difference(m1, "age", from = 30, to = 60, stage = "outcome")
#> First differences: age from 30 to 60 at the covariate profile (delta intervals at 95%; outcome stage only)
#>       alternative  from    to   diff  lower  upper
#>           abstain 0.628 0.158 -0.470 -0.536 -0.403
#>          democrat 0.143 0.388  0.245  0.203  0.288
#>  republican_other 0.185 0.378  0.193  0.147  0.239
#>       nonresponse 0.045 0.076  0.031  0.011  0.051
```

## Local maxima

Mixture likelihoods are not globally concave. On `pa2006` and on
`pres2004` a search from 200 starts finds two interior maxima of the
BIMNL likelihood. The fit screens fourteen starting values and polishes
every one of them. It then challenges the incumbent with rescaled split
slopes. Each point where these searches end is classified by the rules
of the fit-quality flags: a local maximum, a boundary (diverging
coefficients, or a split equation that puts every case in one regime), a
saddle, or not converged. `m1$search` lists the distinct points, and the
number of distinct local maxima is reported, with a message when there
is more than one; the highest found need not be the global maximum. Then
`crosscheck = TRUE` runs the same search under the other split link and
polishes the fit under its own link from the rescaled estimates (which
about doubles the time), and `start` probes other basins (a row of
`m1$search_theta` refits that point).

``` r

m1$n_modes
#> [1] 2
m1$search
#>      loglik    type reached reported detail
#> 1 -2069.971 maximum       5     TRUE       
#> 2 -2070.990 maximum      12    FALSE
```

A fit that reports more than one mode is still the highest maximum
found;
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md) says
so in a message (hidden in this vignette’s output), so that no analysis
rests on a lower maximum unknowingly. Changing the split-equation
covariates can change the number and location of the maxima.

## The second application: the 2004 presidential vote

`pres2004` holds the Campbell and Monson (2008) respondents of the 2004
Election Panel Study (Replication II of Bagozzi and Marchetti 2017).

``` r

data(pres2004)
f2 <- vote ~ educ + female + black + hispanic + pid7 + bush_values + battleground + south + gmb +
  agegroup + evangelical + evangelical:gmb + catholic + catholic:gmb + other_religion +
  secular + secular:gmb + ideology + mobilization
p0 <- mnl(f2, data = pres2004)
p1 <- imnl(vote ~ educ + female + black + hispanic + pid7 + bush_values + battleground + south +
             gmb + agegroup + evangelical + evangelical:gmb + catholic + catholic:gmb +
             other_religion + secular + secular:gmb + ideology + mobilization |
             educ + agegroup + mobilization + secular, data = pres2004)
compare_models(MNL = p0, BIMNL = p1)
#>  model                                                                    type
#>    MNL                                                       Multinomial logit
#>  BIMNL Inflated multinomial logit (inflated alternative: abstain, logit split)
#>  inflate    n df   loglik      AIC      BIC mcfadden converged boundary
#>          1341 40 -553.500 1187.000 1395.046    0.548      TRUE    FALSE
#>  abstain 1341 45 -545.279 1180.557 1414.610    0.555      TRUE    FALSE
#>  ill_conditioned n_modes
#>            FALSE       1
#>            FALSE       2
vuong(p1, p0)
#> Vuong test:  m1 = p1 [ Inflated multinomial logit (inflated alternative: abstain, logit split) ]
#>              m2 = p0 [ Multinomial logit ]
#> positive statistics favor m1
#> 
#>  correction statistic p_m1_better p_m2_better p_two_sided
#>         raw     1.680       0.047       0.953       0.093
#>         AIC     0.658       0.255       0.745       0.510
#>         BIC    -1.998       0.977       0.023       0.046
#> Note: The inflated multinomial logit and the multinomial logit are nested at the boundary phi = 1, where the normal reference does not hold for testing the uninflated model; boot_lr_test() calibrates the statistic.
```

The log-likelihoods of the published Table A.4 reproduce: the inflated
model’s -545.279 agrees with the published -545.27 when truncated to its
printed digits, and a search from 200 starts finds no interior maximum
above it. So do its split-equation coefficients and the outcome-stage
coefficients the package’s tests check. Two of the search’s starts end
higher, at points where the Kerry and Bush coefficients diverge together
(the higher at -542.851): the outcome stage predicts 21 abstainers with
certainty, and the split equation puts the other 84 in the captive
regime. There the likelihood rises towards a supremum at infinity, a
separation within the mixture that the multinomial logit’s data do not
show (see
[`?imnl`](https://bagozzib.github.io/polychoice/reference/imnl.md)). The
multistart does not reach that region, and the published fit is the
highest interior maximum. The raw Vuong statistic computed here is 1.68
(positive values favour the first model, the inflated one; two-sided p =
0.09): it leans towards the inflated model without being decisive. Its
BIC-corrected version, which charges the inflated model for its extra
parameters, is -2.00 and leans the other way, and the information
criteria disagree with each other, as in the published analysis, which
reports mixed support for inflation in these data. Which row to read:
the raw statistic is Vuong’s, the AIC and BIC rows are his corrections
for the difference in parameters; but since the multinomial logit is the
boundary case of the inflated model, none of them has a valid normal
reference here, and `boot_lr_test(p0, p1)` gives the calibrated
p-values.

## Tests for captive choosers

Whether the multinomial logit leaves room for routine choosers can be
asked before any inflated model is fit.
[`captivity_test()`](https://bagozzib.github.io/polychoice/reference/captivity_test.md)
is the score test of the multinomial logit against the dogit with a
constant captive share for each alternative, computed from the
multinomial logit fit alone: its joint statistic is the Lagrange
multiplier test of Tse (1987), and the one-sided statistic for each
alternative is added here. The score of a case for alternative $`j`$ is
$`1\{y_i = j\}/P_{ij} - 1`$: it grows when the cases that chose $`j`$
are ones the multinomial logit thought unlikely to, which is what
captive choosers produce. Each alternative gets a one-sided z statistic
(captive shares cannot be negative), and a joint test follows.

``` r

captivity_test(m0)
#> Score test of the multinomial logit against captive choosers (dogit; one-sided z per alternative, joint chi-square) 
#>                       test statistic df p.value
#>           captive: abstain     0.503 NA  0.3075
#>          captive: democrat    -1.376 NA  0.9156
#>  captive: republican_other     1.282 NA  0.0999
#>       captive: nonresponse    -0.038 NA  0.5153
#>                      joint     3.859  4  0.4254
```

Because the multinomial logit is an ordinary boundary point of the
dogit’s captivity parameters, the score test is regular there, unlike
the comparison with the inflated fit. Its weights $`1/P_{ij}`$ make it
heavy-tailed when some fitted probabilities are small;
`boot_pvalue(m0, captivity_test)` calibrates it by simulation (not run
here). The same weights cost it power, and here it does not reject
although the inflated multinomial logit with a constant captive share of
abstention, its own alternative, raises the log-likelihood by 10.5 on
one parameter. The score of abstention is 37.3, but its expected
information (before the coefficients are partialled out) is 10,186, half
of it from one respondent whose fitted probability of abstaining is
0.0002 and who did not abstain, and along the captive share the
log-likelihood is convex at zero, so the gain lies away from the null,
where a local test does not look. A small statistic is no evidence
against captivity. The likelihood-ratio form measures the gain itself:
it fits the inflated multinomial logit with an intercept-only split for
each requested alternative and refers the likelihood ratio to the 50:50
mixture of chi-square(0) and chi-square(1), since a captive share of
zero is a regular boundary point when the alternative’s fitted
probabilities vary in a way its constants and slopes cannot absorb
(otherwise the share is not identified, and the test says so; one refit
per alternative):

``` r

captivity_test(m0, alternatives = "abstain", type = "lr")
#> Likelihood-ratio test of the multinomial logit against a constant captive share of each alternative (the inflated multinomial logit with an intercept-only split; p-value from the 50:50 mixture of chi-square(0) and chi-square(1)) 
#>              test statistic df  p.value share
#>  captive: abstain    20.951  1 2.36e-06 0.196
```

A split equation with covariates, as in Replication I, is not of this
kind (its slopes are not identified at the boundary), and
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md)
calibrates that comparison.

## Several inflated alternatives

`inflate` can name several alternatives. A case is then either in the
choice regime, where the multinomial logit operates, or captive to one
of the inflated alternatives, with the regime probabilities given by a
multinomial logit in the split covariates: the parametrized logit
captivity model with a captive share for each inflated alternative. In
vote choice this holds routine nonvoters and party loyalists in one
model. Each inflated alternative gets its own split coefficients,
`infl_<alternative>:<term>`, the log-odds of choosing against captivity
to that alternative (positive: less captivity).

``` r

set.seed(1)
dc <- rmnl(2000, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)),
           delta = cbind(a = c(1, 0.8), c = c(1.8, -0.7)), inflate = c("a", "c"))
mc <- imnl(y ~ x1 + x2 | z1, data = dc, inflate = c("a", "c"))
summary(mc)$coefficients[grep("^infl_", names(coef(mc))), ]
#>                      Estimate Std. Error   z value     Pr(>|z|)
#> infl_a:(Intercept)  0.7305562 0.18670464  3.912898 9.119505e-05
#> infl_a:z1           0.7072399 0.09649689  7.329147 2.316215e-13
#> infl_c:(Intercept)  1.9608265 0.23896982  8.205331 2.299558e-16
#> infl_c:z1          -0.8044669 0.13505839 -5.956438 2.577953e-09
head(predict(mc, type = "regime"), 3)
#>      choice         a          c
#> 1 0.4083044 0.5747278 0.01696775
#> 2 0.6484533 0.2001526 0.15139415
#> 3 0.3824451 0.6036018 0.01395317
```

The split coefficients lie within 1.5 standard errors of the values used
to simulate them (1 and 0.8 for `a`, 1.8 and -0.7 for `c`). The
likelihood can also peak near the limit in which the choice regime never
picks `a`, with every choice of `a` ascribed to captivity. In a scan of
20 samples of this design one did (seed 9):

- the outcome constants ran to about 10, with standard errors above 13;
- the split constant of `a` fell to 0.34, with a standard error of
  0.065, against its true value of 1.

That is the data’s own maximum: a refit from the true values reaches it,
and one from the outcome constants raised by 5 returns to it. The fit
raises no flag. Large outcome constants with larger standard errors are
the sign. At such a fit the split coefficients’ small standard errors
are the misleading part, and the bootstrap does not repair them: 200
refits give the split constant of `a` a standard error of 0.084, against
the analytic 0.065, for an estimate 0.34 whose true value is 1, although
they put the outcome constants’ standard errors above 300. Compare
instead with the model that has fewer inflated alternatives, by
[`boot_lr_test()`](https://bagozzib.github.io/polychoice/reference/boot_lr_test.md).

The model with fewer inflated alternatives is the boundary at which a
captive share is zero, so
`boot_lr_test(imnl(y ~ x1 + x2 | z1, data = dc, inflate = "a"), mc)`
compares the two by simulation, as it compares the inflated model with
the multinomial logit.

The split covariate here, `z1`, is not in the outcome equations. In the
simulation reported in
[`?imnl`](https://bagozzib.github.io/polychoice/reference/imnl.md),
splitting on an outcome covariate instead left the intervals of that
covariate’s outcome coefficients at 0.89-0.90 coverage with 5,000 cases
and 0.80-0.81 with 2,000.

## Links and boundaries

`link = "probit"` replaces the logistic split equation by a standard
normal one (with one inflated alternative). When the split equation puts
every case in one regime the inflated model collapses onto the
multinomial logit, and with several inflated alternatives a captive
share can collapse to zero; the fit then warns that it is at a boundary,
and the standard errors of the split coefficients are not to be
interpreted.

## References

Bagozzi, B.E. and Marchetti, K. (2017). Distinguishing occasional
abstention from routine indifference in models of vote choice.
*Political Science Research and Methods*, 5(2), 277-294.

Campbell, D.E. and Monson, J.Q. (2008). The religion card: Gay marriage
and the 2004 presidential election. *Public Opinion Quarterly*, 72,
399-419.

Gaudry, M.J.I. and Dagenais, M.G. (1979). The dogit model.
*Transportation Research Part B*, 13, 105-111.

Swait, J. and Ben-Akiva, M. (1987). Empirical test of a constrained
choice discrete model: Mode choice in Sao Paulo, Brazil. *Transportation
Research Part B*, 21, 103-115.

Tse, Y.K. (1987). A diagnostic test for the multinomial logit model.
*Journal of Business and Economic Statistics*, 5, 283-286.

Vuong, Q.H. (1989). Likelihood ratio tests for model selection and
non-nested hypotheses. *Econometrica*, 57, 307-333.

Wilson, P. (2015). The misuse of the Vuong test for non-nested models to
test for zero-inflation. *Economics Letters*, 127, 51-53.
