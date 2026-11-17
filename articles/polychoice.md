# Getting started with polychoice

## The problem

Many outcomes in the social sciences are unordered choices among several
alternatives: a vote for one of several candidates or for none, a travel
mode, a type of political action, a treaty design. The workhorse model
is the multinomial logit (Theil 1969; McFadden 1974), which comes in two
notations that `polychoice` treats as one model:

- the **multinomial logit** of the wide layout, one row per case, with
  case-specific covariates whose coefficients differ by alternative
  (`y ~ x1 + x2`);
- the **conditional logit** of the long layout, one row per case and
  alternative, with alternative-specific covariates whose coefficients
  are common to all alternatives (`chosen ~ w | x`, with `id` and
  `alt`).

Both are the same choice probability

``` math
P(y_i = j) = \frac{a_{ij} \exp(V_{ij})}{\sum_k a_{ik} \exp(V_{ik})}, \qquad
V_{ij} = o_{ij} + x_i'\beta_j + w_{ij}'\gamma,
```

with $`a_{ij}`$ indicating whether alternative $`j`$ is available to
case $`i`$, $`\beta_b = 0`$ for the base alternative, and an offset
$`o_{ij}`$.

The model’s convenience comes with the independence of irrelevant
alternatives (IIA): the odds between any two alternatives do not depend
on the others. `polychoice` provides the estimators that relax IIA
(nested logit, heteroskedastic extreme value, multinomial probit, random
intercepts), the baseline-inflated multinomial logit of Bagozzi and
Marchetti (2017) for a base alternative that mixes two kinds of cases,
and the full battery of tests of IIA and of specification, all behind
one interface. This vignette covers the interface; the companion
vignettes cover the inflated model, the tests, and the models beyond
IIA.

## One interface, seven estimators

| Function | Model |
|----|----|
| [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) | multinomial / conditional logit |
| [`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md) | inflated multinomial logit (one or several inflated alternatives) |
| [`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md) | two-level nested logit |
| [`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md) | heteroskedastic extreme value |
| [`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md) | multinomial probit with a free covariance |
| [`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md) | multinomial logit with random alternative-specific intercepts |
| [`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md) | rank-ordered (exploded) logit for rankings |

Every estimator takes the same arguments for the data layout (`id`,
`alt`), the base alternative (`ref`), frequency `weights`, an `offset`,
unit fixed effects (`fe`; not
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
whose split equation they leave identified only through the functional
form, nor
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md),
whose unit effects are the random intercepts, nor
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)),
and standard errors (`se`: inverse observed information, outer product
of the scores, robust, cluster-robust, or bootstrap), and returns an
object of class `"choicemodel"` with the usual methods (`summary`,
`coef`, `vcov`, `predict`, `confint`, `logLik`, `update`, `simulate`)
plus `broom` tidiers and `texreg` / `modelsummary` support.

## A first fit: the 2006 Pennsylvania state House vote

`pa2006` holds the survey of Arceneaux and Kolodny (2009) after a field
experiment in which an interest group endorsed the Democratic candidate
by canvassing or by telephone. The outcome has four alternatives, with
abstention as the base.

``` r

library(polychoice)
data(pa2006)
table(pa2006$vote)
#> 
#>          abstain         democrat republican_other      nonresponse 
#>              731              579              563              125
m <- mnl(vote ~ canvass + phone + democrat + republican + age + female + voted2004,
         data = pa2006)
summary(m)
#> Multinomial logit 
#> Call: mnl(formula = vote ~ canvass + phone + democrat + republican + 
#>     age + female + voted2004, data = pa2006)
#> Alternatives: abstain, democrat, republican_other, nonresponse  (base: abstain)
#> 
#> Outcome equations (alternative:term):
#>                                Estimate Std. Error z value Pr(>|z|)    
#> democrat:(Intercept)         -2.9127340  0.3292465  -8.847  < 2e-16 ***
#> democrat:canvass             -0.0092838  0.2063913  -0.045 0.964122    
#> democrat:phone               -0.0428430  0.1907265  -0.225 0.822267    
#> democrat:democrat             0.7685832  0.1453229   5.289 1.23e-07 ***
#> democrat:republican          -1.1743575  0.1743220  -6.737 1.62e-11 ***
#> democrat:age                  0.0269982  0.0035721   7.558 4.09e-14 ***
#> democrat:female              -0.3931622  0.1296498  -3.032 0.002425 ** 
#> democrat:voted2004            1.7253797  0.2219251   7.775 7.57e-15 ***
#> republican_other:(Intercept) -3.5103280  0.3469268 -10.118  < 2e-16 ***
#> republican_other:canvass      0.3462358  0.2064414   1.677 0.093511 .  
#> republican_other:phone        0.1096400  0.1925853   0.569 0.569149    
#> republican_other:democrat    -0.6899875  0.1860326  -3.709 0.000208 ***
#> republican_other:republican   0.6707407  0.1441201   4.654 3.25e-06 ***
#> republican_other:age          0.0221353  0.0034814   6.358 2.04e-10 ***
#> republican_other:female      -0.1423174  0.1273234  -1.118 0.263668    
#> republican_other:voted2004    2.0453916  0.2436497   8.395  < 2e-16 ***
#> nonresponse:(Intercept)      -5.1867961  0.8219736  -6.310 2.79e-10 ***
#> nonresponse:canvass           0.0006947  0.3420094   0.002 0.998379    
#> nonresponse:phone             0.0490322  0.3137972   0.156 0.875832    
#> nonresponse:democrat         -1.0333245  0.2579053  -4.007 6.16e-05 ***
#> nonresponse:republican       -1.0603861  0.2360428  -4.492 7.04e-06 ***
#> nonresponse:age               0.0186464  0.0057459   3.245 0.001174 ** 
#> nonresponse:female            0.4688209  0.2287007   2.050 0.040371 *  
#> nonresponse:voted2004         3.0050713  0.7209875   4.168 3.07e-05 ***
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Standard errors: analytic. Log-likelihood -2125.695 on 24 parameters, 1998 cases. AIC 4299.39, BIC 4433.79.
```

The coefficient names are `alternative:term`: each non-base alternative
has its own equation against abstention.
[`rrr()`](https://bagozzib.github.io/polychoice/reference/rrr.md)
reports relative risk ratios (exponentiated coefficients) with
confidence intervals; [`summary()`](https://rdrr.io/r/base/summary.html)
also shows the information criteria and the standard-error type.

``` r

head(rrr(m), 6)
#>                  term   block   rrr lower upper      z  p.value
#>  democrat:(Intercept) outcome 0.054 0.028 0.104 -8.847  < 2e-16
#>      democrat:canvass outcome 0.991 0.661 1.485 -0.045    0.964
#>        democrat:phone outcome 0.958 0.659 1.392 -0.225    0.822
#>     democrat:democrat outcome 2.157 1.622 2.867  5.289 1.23e-07
#>   democrat:republican outcome 0.309 0.220 0.435 -6.737 1.62e-11
#>          democrat:age outcome 1.027 1.020 1.035  7.558 4.09e-14
```

### Predicted probabilities and first differences

[`predict()`](https://rdrr.io/r/stats/predict.html) returns the
probability of every alternative for the estimation rows or for
`newdata`; `type = "class"` picks the most likely alternative.

``` r

head(predict(m), 3)
#>     abstain   democrat republican_other nonresponse
#> 1 0.5149576 0.05840060        0.3747743  0.05186757
#> 2 0.3988234 0.07973611        0.4620162  0.05942424
#> 3 0.2399238 0.11691714        0.5770150  0.06614410
nd <- data.frame(canvass = c(0, 1), phone = 0, democrat = 1, republican = 0, age = 45,
                 female = 1, voted2004 = 1)
predict(m, newdata = nd)
#>     abstain  democrat republican_other nonresponse
#> 1 0.3428170 0.5129521       0.09332101  0.05090984
#> 2 0.3315748 0.4915459       0.12760476  0.04927454
```

[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md)
is the quantity of interest most analyses report: the change in every
alternative’s probability when one covariate moves between two values.
By default in the wide layout the other covariates are held at a profile
of their weighted means and modes; `newdata` gives another one-case
profile; `average = TRUE` (the default in the long layout) averages the
first difference over the estimation cases, the average effect many
applied reports prefer to the effect at a profile. The columns `from`
and `to` are the two predicted probabilities. Delta-method intervals are
the default; `ci = "sim"` draws from the asymptotic normal distribution
of the estimates instead.

``` r

first_difference(m, "canvass", from = 0, to = 1)
#> First differences: canvass from 0 to 1 at the covariate profile (delta intervals at 95%)
#>       alternative  from    to   diff  lower upper
#>           abstain 0.423 0.385 -0.038 -0.118 0.041
#>          democrat 0.271 0.244 -0.027 -0.096 0.042
#>  republican_other 0.248 0.318  0.071 -0.006 0.147
#>       nonresponse 0.058 0.053 -0.005 -0.038 0.028
first_difference(m, "canvass", from = 0, to = 1, average = TRUE)
#> First differences: canvass from 0 to 1, averaged over the estimation cases (delta intervals at 95%)
#>       alternative  from    to   diff  lower upper
#>           abstain 0.375 0.344 -0.031 -0.096 0.033
#>          democrat 0.296 0.275 -0.022 -0.081 0.037
#>  republican_other 0.264 0.324  0.059 -0.004 0.123
#>       nonresponse 0.064 0.058 -0.006 -0.041 0.028
```

### Fit and classification

[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md)
gives the log-likelihood, information criteria, and the pseudo-$`R^2`$
measures in common use (all relative to the null model with
alternative-specific constants only);
[`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md)
cross-tabulates the observed and predicted alternatives.

``` r

fit_stats(m)
#> Log-likelihood -2125.695 (null -2511.718), 24 parameters, 1998 cases; AIC 4299.390, BIC 4433.788
#> McFadden R2 0.154 (adjusted 0.144); Cox-Snell 0.321; Nagelkerke 0.349
#> Count R2 0.553 (adjusted: 0.295 proportional reduction in error over the modal alternative)
classification(m)
#>                   predicted
#> observed           abstain democrat republican_other nonresponse
#>   abstain              365      194              172           0
#>   democrat             101      402               76           0
#>   republican_other     109      116              338           0
#>   nonresponse           41       54               30           0
#> 
#> Correctly classified: 0.553; proportional reduction in error over the modal alternative: 0.295
```

## The long layout: alternative-specific covariates

When covariates vary across alternatives (travel time, price, distance
to a polling place), the data come as one row per case and alternative
and the formula has two parts: `chosen ~ w | x`, where `w` are the
alternative-specific covariates with a common coefficient and `x` the
case-specific ones with alternative-specific coefficients. The intercept
of the first part switches the alternative-specific constants on (the
default); `chosen ~ 0 + w | x` switches them off. The classic example is
the Australian travel mode data distributed with the `AER` package.

``` r

data("TravelMode", package = "AER")
TM <- TravelMode
TM$chosen <- as.integer(TM$choice == "yes")
head(TM, 4)
#>   individual  mode choice wait vcost travel gcost income size chosen
#> 1          1   air     no   69    59    100    70     35    1      0
#> 2          1 train     no   34    31    372    71     35    1      0
#> 3          1   bus     no   35    25    417    70     35    1      0
#> 4          1   car    yes    0    10    180    30     35    1      1
mc <- mnl(chosen ~ wait + vcost + travel + gcost | income + size, data = TM,
          id = "individual", alt = "mode")
summary(mc)
#> Conditional logit 
#> Call: mnl(formula = chosen ~ wait + vcost + travel + gcost | income + 
#>     size, data = TM, id = "individual", alt = "mode")
#> Alternatives: air, train, bus, car  (base: air)
#> 
#> Outcome equations (alternative:term):
#>                    Estimate Std. Error z value Pr(>|z|)    
#> train:(Intercept)  0.421795   1.073173   0.393   0.6943    
#> train:income      -0.067576   0.016455  -4.107 4.01e-05 ***
#> train:size         0.693524   0.378587   1.832   0.0670 .  
#> bus:(Intercept)   -0.570171   1.210541  -0.471   0.6376    
#> bus:income        -0.028019   0.017093  -1.639   0.1012    
#> bus:size           0.290801   0.465637   0.625   0.5323    
#> car:(Intercept)   -5.286500   1.202630  -4.396 1.10e-05 ***
#> car:income        -0.008078   0.013419  -0.602   0.5472    
#> car:size           0.530701   0.321034   1.653   0.0983 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Generic coefficients (alternative-specific covariates):
#>         Estimate Std. Error z value Pr(>|z|)    
#> wait   -0.102555   0.011376  -9.015  < 2e-16 ***
#> vcost  -0.053353   0.025282  -2.110  0.03484 *  
#> travel -0.010250   0.003418  -2.999  0.00271 ** 
#> gcost   0.046426   0.024900   1.865  0.06225 .  
#> ---
#> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
#> 
#> Standard errors: analytic. Log-likelihood -170.688 on 13 parameters, 210 cases. AIC 367.38, BIC 410.89.
```

A covariate in the first part can be given an alternative-specific
coefficient with `w:alt` (for example `vcost:mode`, a cost coefficient
per mode), so no third formula part is needed. Alternatives absent from
a case’s rows are treated as unavailable to that case. This is also how
a stated consideration set enters: when a survey asks which parties a
respondent would consider voting for, keeping only those parties’ rows
makes [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
the model of the choice within the considered set (which parties are
considered is then a separate set of binary models; Oscarsson and Rosema
2019 recommend such direct measures). A consideration set inferred from
the choices alone, without such a question, is weakly identified in the
designs political scientists have: the package’s verification found the
standard errors of the consideration parameters as large as the
parameters at 2,000 cases with respondent characteristics alone,
recovery failing when the same voter-party distances drive both stages,
and, even with a distance that shapes consideration only, intervals that
covered 0.88-0.91 instead of 0.95 in a realistic design. The package
does not fit it.

## Rankings

When respondents rank the alternatives instead of choosing one (a survey
that asks voters to order the parties, the ballots of a ranked-choice
election),
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md)
fits the rank-ordered logit (Beggs, Cardell and Hausman 1981), also
called the exploded logit after Chapman and Staelin’s (1982) “explosion”
of a ranking into independent choices. A ranking is read as a sequence
of multinomial-logit choices: the first-ranked alternative from the
whole choice set, the second from those left, and so on. Partial
rankings (only the first few positions) and tied ranks are handled
exactly (Allison and Christakis 1994). The data are in the long layout
with the rank as the response, smaller preferred and 0 for an
alternative left unranked;
[`rrologit()`](https://bagozzib.github.io/polychoice/reference/rrologit.md)
simulates such data.

``` r

set.seed(8)
dk <- rrologit(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)),
               gamma = 0.7, depth = 3)
head(dk, 4)
#>   id alt rank         w1          x1       x2
#> 1  1   a    2  0.9824732 -0.08458607 0.492746
#> 2  1   b    1 -0.4606749 -0.08458607 0.492746
#> 3  1   c    3  0.7402164 -0.08458607 0.492746
#> 4  1   d    0 -0.4267571 -0.08458607 0.492746
mro <- rologit(rank ~ w1 | x1 + x2, data = dk, id = "id", alt = "alt")
mfc <- rologit(rank ~ w1 | x1 + x2, data = dk, id = "id", alt = "alt", depth = 1)   # first choices only
cbind(rankings = sqrt(diag(vcov(mro))), first_choices = sqrt(diag(vcov(mfc))))[c(1:3, 10), ]
#>                 rankings first_choices
#> b:(Intercept) 0.07348825    0.11778461
#> b:x1          0.08144944    0.12364952
#> b:x2          0.08078852    0.12496424
#> w1            0.03324162    0.05316619
```

The rankings estimate the same coefficients as the first choices with
smaller standard errors. The model assumes that preferences are the same
at every stage of the ranking;
[`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md)
refits with separate coefficients for the stages after the first and
tests the difference (Hausman and Ruud 1987, section 3; Allison and
Christakis 1994, section 5.4):

``` r

rank_test(mro, after = 1)
#> Stability of preferences across ranking stages (position 1 against later; 800 ranking(s) with later stages) 
#>              test statistic df p.value
#>  Likelihood ratio    17.721 10  0.0599
#>              Wald    17.659 10  0.0610
```

## Standard errors, weights, and panels

`se = "robust"` gives the sandwich over cases;
`cluster = "<cluster column>"` the cluster-robust version (and selects
`se = "cluster"`); `se = "bootstrap"` resamples cases, or clusters when
`cluster` is given (units with unit fixed effects or random intercepts),
with `nboot` refits. `weights` are frequency weights: a case with weight
2 counts as two identical cases in the estimates, the model-based
standard errors, the information criteria, and the averaged quantities
of interest. Survey (probability) weights are not supported:
`se = "robust"` is the frequency-weight sandwich, not a design-based
variance (see
[`?mnl`](https://bagozzib.github.io/polychoice/reference/mnl.md)).

``` r

mr <- mnl(vote ~ canvass + phone + democrat + republican + age + female + voted2004,
          data = pa2006, se = "robust")
cbind(analytic = sqrt(diag(vcov(m)))[1:4], robust = sqrt(diag(vcov(mr)))[1:4])
#>                       analytic    robust
#> democrat:(Intercept) 0.3292465 0.3296794
#> democrat:canvass     0.2063913 0.2053533
#> democrat:phone       0.1907265 0.1901162
#> democrat:democrat    0.1453229 0.1450103
```

For panels, `fe = "<unit column>"` (the name of the column that
identifies the units) adds alternative-specific constants per unit (with
the never-chosen rule and, optionally, `bias_correct = "jackknife"`),
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)
integrates out random alternative-specific intercepts, and
[`mundlak()`](https://bagozzib.github.io/polychoice/reference/mundlak.md)
adds unit means of the covariates for a correlated random-effects
specification. These are covered in the vignette *Beyond IIA*.

## Reporting

[`broom::tidy()`](https://generics.r-lib.org/reference/tidy.html) and
[`broom::glance()`](https://generics.r-lib.org/reference/glance.html)
work on every fit, and `texreg` and `modelsummary` tables can be built
directly from the fitted objects.

``` r

head(broom::tidy(m, conf.int = TRUE), 4)
#>                   term     estimate std.error   statistic      p.value
#> 1 democrat:(Intercept) -2.912733973 0.3292465 -8.84666761 9.017146e-19
#> 2     democrat:canvass -0.009283797 0.2063913 -0.04498153 9.641220e-01
#> 3       democrat:phone -0.042843030 0.1907265 -0.22463074 8.222665e-01
#> 4    democrat:democrat  0.768583198 0.1453229  5.28879633 1.231239e-07
#>   component   conf.low  conf.high
#> 1   outcome -3.5580452 -2.2674228
#> 2   outcome -0.4138033  0.3952357
#> 3   outcome -0.4166601  0.3309740
#> 4   outcome  0.4837556  1.0534108
broom::glance(m)
#>      logLik     AIC      BIC df nobs converged boundary ill_conditioned
#> 1 -2125.695 4299.39 4433.788 24 1998      TRUE    FALSE           FALSE
#>   singular n_modes
#> 1    FALSE       1
```

With multiply imputed data, fit the model on each completed data set and
pass the list of fits to
[`mice::pool()`](https://amices.org/mice/reference/pool.html), which
combines them by Rubin’s rules through the same methods. The estimators
take the data as an argument, so `with(imp, mnl(...))` does not work.
The fits’ flags are not pooled, so check each fit first. And give
`pool()` the complete-data degrees of freedom, cases minus parameters:
left to itself, mice 3.19 counts one residual per case and alternative.
The pooled statistic tests each parameter against zero
([`?tidy.choicemodel`](https://bagozzib.github.io/polychoice/reference/tidy.choicemodel.md)).

``` r

set.seed(1)
dm <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
dm$x1[sample(500, 80)] <- NA
imp <- mice::mice(dm, m = 5, printFlag = FALSE, seed = 1)
fits <- lapply(1:5, function(i) mnl(y ~ x1 + x2, data = mice::complete(imp, i)))
all(vapply(fits, function(f) f$converged && !f$boundary && !f$singular, logical(1)))
#> [1] TRUE
pooled <- mice::pool(fits, dfcom = nobs(fits[[1]]) - length(coef(fits[[1]])))
summary(pooled)[, c("term", "estimate", "std.error", "statistic", "df", "p.value")]
#>            term   estimate std.error statistic        df      p.value
#> 1 b:(Intercept)  0.4821263 0.1232712  3.911103 475.78603 1.052452e-04
#> 2          b:x1 -0.7073566 0.1777954 -3.978486  31.11162 3.859728e-04
#> 3          b:x2  0.6748004 0.1309133  5.154561 441.90770 3.841501e-07
#> 4 c:(Intercept) -0.3711618 0.1679398 -2.210088 253.25992 2.799368e-02
#> 5          c:x1  0.9470919 0.1783153  5.311331 119.41062 5.113435e-07
#> 6          c:x2 -0.6314056 0.1471864 -4.289837 330.41730 2.350910e-05
```

## Where next

- *The inflated multinomial logit* fits the model of Bagozzi and
  Marchetti (2017) and reproduces the published applications.
- *Testing IIA and specification* covers
  [`iia_tests()`](https://bagozzib.github.io/polychoice/reference/iia_tests.md)
  and the further tests, with the size study that motivates the
  recommendations.
- *Beyond IIA* covers the nested logit, HEV, multinomial probit, random
  intercepts, and the fixed-effects tools.
- *Seven published analyses* works through the example datasets drawn
  from replication archives (`sanctions`, `fluvax`, `bes1987`,
  `defection`, `bes2010`, `cabinets`, `afghan_leaders`), four of them
  reproducing their published models, one of them with Firth’s penalized
  likelihood under separation.

## References

Allison, P.D. and Christakis, N.A. (1994). Logit models for sets of
ranked items. *Sociological Methodology*, 24, 199-228.

Arceneaux, K. and Kolodny, R. (2009). Educating the least informed:
Group endorsements in a grassroots campaign. *American Journal of
Political Science*, 53, 755-770.

Bagozzi, B.E. and Marchetti, K. (2017). Distinguishing occasional
abstention from routine indifference in models of vote choice.
*Political Science Research and Methods*, 5(2), 277-294.

Beggs, S., Cardell, S. and Hausman, J. (1981). Assessing the potential
demand for electric cars. *Journal of Econometrics*, 17, 1-19.

Chapman, R.G. and Staelin, R. (1982). Exploiting rank ordered choice set
data within the stochastic utility model. *Journal of Marketing
Research*, 19, 288-301.

Hausman, J.A. and Ruud, P.A. (1987). Specifying and testing econometric
models for rank-ordered data. *Journal of Econometrics*, 34, 83-104.

McFadden, D. (1974). Conditional logit analysis of qualitative choice
behavior. In P. Zarembka (ed.), *Frontiers in Econometrics*, 105-142.
Academic Press.

Oscarsson, H. and Rosema, M. (2019). Consideration set models of
electoral choice: Theory, method, and application. *Electoral Studies*,
57, 256-262.

Theil, H. (1969). A multinomial extension of the linear logit model.
*International Economic Review*, 10, 251-259.
