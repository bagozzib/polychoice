# Rank-ordered (exploded) logit regression

Fits the rank-ordered logit (Beggs, Cardell and Hausman 1981, who call
it the ordered logit), also called the exploded logit after Chapman and
Staelin's (1982) "explosion" of a ranking into independent choices
(Allison and Christakis 1994): the model for data in which each case
ranks some or all of its available alternatives, such as survey rankings
of parties or candidates, or the ballots of a ranked-choice election.
With the utilities of
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md), a
ranking is a sequence of multinomial-logit choices: the first-ranked
alternative is chosen from the whole choice set, the second from those
left, and so on. Rankings carry more information than first choices, so
the same coefficients are estimated more precisely than from the first
choices alone.

## Usage

``` r
rologit(
  formula,
  data,
  id,
  alt,
  ref = NULL,
  weights = NULL,
  offset = NULL,
  depth = NULL,
  se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"),
  cluster = NULL,
  nboot = 200L,
  cores = 1L,
  start = NULL,
  maxit = 1000L,
  reltol = 1e-10,
  firth = FALSE,
  prior = NULL
)
```

## Arguments

- formula:

  `rank ~ w | x` (see Details); `rank ~ w` without case-specific
  covariates.

- data:

  A data frame.

- id, alt:

  Names of the case and alternative columns of `data`.

- ref:

  The base alternative: a level of the alternative column `alt`
  (default: its first level), whose constant and case-specific
  coefficients are zero.

- weights:

  Optional frequency weights: a column name in `data` or a numeric
  vector with one entry per row (constant within a case in the long
  layout).

- offset:

  Optional offset on the utility scale: a column name of `data` or a
  vector with one entry per row.

- depth:

  Optional: use only the first `depth` positions of every ranking;
  alternatives ranked below become unranked (a tie group that starts
  within the first `depth` positions is kept whole). `depth = 1` is the
  multinomial logit of the first choices.

- se:

  Standard errors: `"analytic"` (inverse observed information, default),
  `"opg"`, `"robust"`, `"cluster"` (needs `cluster`), `"bootstrap"`
  (nonparametric over cases, or over clusters when `cluster` is given,
  with `nboot` refits), or `"none"`.

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

  Optional starting values: the full parameter vector in the order of
  [`coef()`](https://rdrr.io/r/stats/coef.html) (`object$theta` of an
  earlier fit); rarely needed, since the log-likelihood is globally
  concave.

- maxit, reltol:

  Controls for the BFGS stage passed to
  [`stats::optim()`](https://rdrr.io/r/stats/optim.html) (the plain
  [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md) is
  solved by damped Newton steps first, with BFGS as the fallback); the
  Newton polish that follows runs up to 25 further iterations with the
  analytic Hessian of the logit family or central differences of the
  analytic gradient otherwise.

- firth:

  If `TRUE`, Firth's penalized likelihood (see Details); default
  `FALSE`, maximum likelihood. Weights are read as frequency weights and
  must total more than 2.

- prior:

  Optional: priors for the coefficients, made by
  [`coef_prior()`](https://bagozzib.github.io/polychoice/reference/coef_prior.md);
  the fit is the posterior mode, a penalized fit like `firth = TRUE`.
  Not with `firth = TRUE`.

## Value

An object of class `c("rologit", "choicemodel")`, with the components of
an [`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
fit; `design$str$rank` holds the rankings used (dense ranks by case, 0
for unranked or unavailable alternatives).

## Details

Case \\i\\ ranks the alternatives \\r\_{i1}, r\_{i2}, \dots\\ of its
choice set \\A_i\\; with \\R\_{il}\\ the alternatives not ranked before
position \\l\\, \$\$\log L_i = \sum_l \Big\[ V\_{i r\_{il}} - \log
\sum\_{m \in R\_{il}} e^{V\_{im}} \Big\],\$\$ with \\V\_{ij}\\ as in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md). The
log-likelihood is globally concave. Partial rankings (the first few
positions ranked, the rest left unranked) are handled exactly: the
unranked alternatives are tied below every ranked one and contribute
only through the denominators. Ties among ranked alternatives are
handled by the exact likelihood of Allison and Christakis (1994, eq. 8):
the tied alternatives are taken to have an unknown strict order, and the
likelihood sums over every order (at most 7 tied alternatives above
others; Breslow's approximation, the default of proportional-hazards
software, attenuates the estimates and is not used).

The model assumes that the relative preference for two alternatives does
not depend on which alternatives have already been ranked: independence
of irrelevant alternatives at every stage. Respondents may rank the
first positions more carefully than the later ones: Chapman and Staelin
(1982) and Hausman and Ruud (1987) argue, and find in their
applications, that later ranks can be noisier, Hausman and Ruud adding
that the pattern depends on the survey.
[`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md)
tests whether the coefficients differ between early and late stages, and
`depth` limits the fit to the first positions (Chapman and Staelin's
explosion depth, Hausman and Ruud's top-P estimator; with only the first
position, the model is the multinomial logit of the first choices).

The data are in the long layout: one row per case and available
alternative, with `id` and `alt` naming the case and alternative
columns, as in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md). The
response is the rank: smaller numbers are preferred, equal numbers are
ties (the numbers need only be in the right order), and 0 marks an
alternative left unranked.

- An unranked alternative is in the choice set, below every ranked one:
  it enters the denominators of every stage (the respondent saw it and
  did not rank it).

- An alternative without a row for the case is unavailable and enters
  none. Drop a row only where the alternative was not on offer.

- Without a penalty, ties are exact up to 7 alternatives tied above
  others. A larger tie group above other alternatives is refused, and a
  tie at the bottom, below no other alternative, is unlimited.

The formula is `rank ~ w | x`: alternative-specific covariates with
generic coefficients, then case-specific covariates with
alternative-specific coefficients. Cases in which no alternative is
ranked above another carry no information and are dropped with a
message. Weights, offsets, clusters, and standard errors are as in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md); a
case (one ranking) is the observation. Unit fixed effects are not
available.

Fitted values, [`predict()`](https://rdrr.io/r/stats/predict.html),
[`first_difference()`](https://bagozzib.github.io/polychoice/reference/first_difference.md),
and [`rrr()`](https://bagozzib.github.io/polychoice/reference/rrr.md)
describe first choices: the multinomial-logit probabilities over the
choice set. `predict(type = "ranking")` gives the probability of each
observed ranking. `y` holds each case's first choice (with a tie at the
top, the first of the tied alternatives in level order), which is what
[`classification()`](https://bagozzib.github.io/polychoice/reference/classification.md),
the count measures of
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md),
and [`residuals()`](https://rdrr.io/r/stats/residuals.html) use; the
log-likelihood measures of
[`fit_stats()`](https://bagozzib.github.io/polychoice/reference/fit_stats.md)
are those of the rankings, against the rank-ordered logit with
alternative-specific constants only.
[`simulate()`](https://rdrr.io/r/stats/simulate.html) draws rankings.
The tests of the independence of irrelevant alternatives and the
goodness-of-fit test apply to single-choice fits and refuse a
rank-ordered one;
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
and
[`vuong()`](https://bagozzib.github.io/polychoice/reference/vuong.md)
compare two rank-ordered fits of the same rankings.

Separation occurs in rankings as in single choices: when a combination
of the covariates ranks every observed choice first among the
alternatives still unranked at its position, the maximum-likelihood
estimates run off (the monotone likelihood of the Cox model, of which
this is the version stratified by case; Bryson and Johnson 1981).
Ranking fewer positions (`depth`, or partial rankings) leaves fewer such
constraints and makes it more likely. The maximum-likelihood fit decides
it exactly (the linear program of
[`separation()`](https://bagozzib.github.io/polychoice/reference/separation.md),
ties as groups) and names the diverging coefficients in a warning; an
alternative never ranked above another available one is the simplest
case (its constant runs to minus infinity). `firth = TRUE` penalizes the
log-likelihood by half the log-determinant of its information, the sum
over the stages of each ranking of the information of the stage's choice
given the alternatives left: Heinze and Schemper's (2001) penalty for
the Cox model, with one stratum per case. A tie group of \\d\\
alternatives enters with the information of the exact tie probability
(the probability that the tied set is the top \\d\\ of the alternatives
left), the expected information of that outcome over the possible sets,
which reduces to the conditional-logit information of an untied
position. The estimates are finite under separation; the rest is as for
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)
(section "Penalized (Firth) estimation"): the log-likelihood at the
penalized estimates, standard errors from the inverse information, the
separated coefficients in `separation`, the run from a second start,
profile intervals from
[`confint()`](https://rdrr.io/r/stats/confint.html), and penalized
likelihood-ratio tests from
[`term_test()`](https://bagozzib.github.io/polychoice/reference/term_test.md),
[`lr_test()`](https://bagozzib.github.io/polychoice/reference/lr_test.md)
and
[`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md).
Under the penalty or a prior a tie group is refused when more than 62
alternatives remain at its positions, when it has more than \\10^5\\
subsets of them to enumerate, or (Firth's penalty) when the second
derivatives of its subsets would need more than 200 MB; without a
penalty only ties of more than 7 alternatives ranked above others are
refused.

## References

Allison, P.D. and Christakis, N.A. (1994). Logit models for sets of
ranked items. *Sociological Methodology*, 24, 199-228.
[doi:10.2307/270983](https://doi.org/10.2307/270983)

Beggs, S., Cardell, S. and Hausman, J. (1981). Assessing the potential
demand for electric cars. *Journal of Econometrics*, 17, 1-19.
[doi:10.1016/0304-4076(81)90056-7](https://doi.org/10.1016/0304-4076%2881%2990056-7)

Bryson, M.C. and Johnson, M.E. (1981). The incidence of monotone
likelihood in the Cox model. *Technometrics*, 23, 381-383.
[doi:10.1080/00401706.1981.10487683](https://doi.org/10.1080/00401706.1981.10487683)

Chapman, R.G. and Staelin, R. (1982). Exploiting rank ordered choice set
data within the stochastic utility model. *Journal of Marketing
Research*, 19, 288-301.
[doi:10.1177/002224378201900302](https://doi.org/10.1177/002224378201900302)

Hausman, J.A. and Ruud, P.A. (1987). Specifying and testing econometric
models for rank-ordered data. *Journal of Econometrics*, 34, 83-104.
[doi:10.1016/0304-4076(87)90068-6](https://doi.org/10.1016/0304-4076%2887%2990068-6)

Heinze, G. and Schemper, M. (2001). A solution to the problem of
monotone likelihood in Cox regression. *Biometrics*, 57, 114-119.
[doi:10.1111/j.0006-341X.2001.00114.x](https://doi.org/10.1111/j.0006-341X.2001.00114.x)

## See also

[`rank_test()`](https://bagozzib.github.io/polychoice/reference/rank_test.md),
[`rrologit()`](https://bagozzib.github.io/polychoice/reference/rrologit.md),
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md)

Other estimators:
[`hev()`](https://bagozzib.github.io/polychoice/reference/hev.md),
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md),
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
[`mnprobit()`](https://bagozzib.github.io/polychoice/reference/mnprobit.md),
[`nlogit()`](https://bagozzib.github.io/polychoice/reference/nlogit.md),
[`remnl()`](https://bagozzib.github.io/polychoice/reference/remnl.md)

## Examples

``` r
set.seed(4)
d <- rrologit(600, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)),
              gamma = 0.7, depth = 3)
head(d)
#>   id alt rank          w1         x1         x2
#> 1  1   a    2  0.03349124  0.2167549 -2.1532042
#> 2  1   b    0 -0.50828566  0.2167549 -2.1532042
#> 3  1   c    1 -0.28867304  0.2167549 -2.1532042
#> 4  1   d    3  0.85654205  0.2167549 -2.1532042
#> 5  2   a    3 -1.45849674 -0.5424926 -0.5727423
#> 6  2   b    1  1.65373350 -0.5424926 -0.5727423
m <- rologit(rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
summary(m)
#> Rank-ordered logit 
#> Call: rologit(formula = rank ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
#> Alternatives: a, b, c, d  (base: a)
#> 
#> Outcome equations (alternative:term):
#>               Estimate Std. Error z value Pr(>|z|)    
#> b:(Intercept)  0.44660    0.08554   5.221 1.78e-07 ***
#> b:x1          -1.03207    0.10025 -10.295  < 2e-16 ***
#> b:x2           0.84399    0.09302   9.073  < 2e-16 ***
#> c:(Intercept) -0.44848    0.08447  -5.309 1.10e-07 ***
#> c:x1           0.64317    0.08868   7.253 4.08e-13 ***
#> c:x2          -0.73272    0.08889  -8.243  < 2e-16 ***
#> d:(Intercept)  0.22291    0.08008   2.784 0.005377 ** 
#> d:x1           0.38290    0.08490   4.510 6.49e-06 ***
#> d:x2           0.32191    0.08343   3.858 0.000114 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Generic coefficients (alternative-specific covariates):
#>    Estimate Std. Error z value Pr(>|z|)    
#> w1  0.77039    0.04178   18.44   <2e-16 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Standard errors: analytic. Log-likelihood -1429.932 on 10 parameters, 600 cases. AIC 2879.86, BIC 2923.83.
rank_test(m, after = 1)
#> Stability of preferences across ranking stages (position 1 against later; 600 ranking(s) with later stages) 
#>              test statistic df p.value
#>  Likelihood ratio     3.915 10   0.951
#>              Wald     3.891 10   0.952
```
