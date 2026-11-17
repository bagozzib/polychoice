# Government formation in the German states, 1990-2023 (Baeck, Debus and Imre 2026)

Every potential cabinet (every combination of the parties in parliament)
at 125 government formations after German state elections, with the
cabinet that formed (Baeck, Debus and Imre 2026). One row per formation
and potential cabinet: the long layout of a conditional logit whose
alternatives are unlabelled and differ from one formation to the next (3
to 63 potential cabinets), so the model has generic coefficients and no
alternative-specific constants.

## Usage

``` r
data(cabinets)
```

## Format

A data frame with 3,352 rows and 18 variables:

- election:

  Formation opportunity (state and election year).

- state:

  German state.

- year:

  Election year.

- cabinet:

  Factor, the potential cabinet's parties (see above).

- chosen:

  1 for the cabinet that formed.

- incumbent:

  1 if the cabinet is the incumbent coalition and the prime minister did
  not change.

- incumbent_newpm:

  1 if the cabinet is the incumbent coalition with a replaced prime
  minister.

- minority:

  1 for a minority cabinet.

- min_winning:

  1 for a minimal winning coalition.

- bargaining:

  1 for the bargaining-proposition cabinet (a minimal winning coalition
  with the fewest parties).

- strongest:

  1 if the cabinet includes the largest party.

- cross_cutting:

  1 for a cross-cutting coalition.

- het_economic, het_social:

  Policy heterogeneity within the cabinet on the economic and the
  sociocultural dimension.

- pea:

  1 if the parties had a pre-electoral alliance.

- rej:

  1 if one of the parties had ruled out governing with another before
  the election.

- pea_res, rej_res:

  First-stage residuals of `pea` and `rej` (two-stage residual
  inclusion).

## Source

Baeck, H., Debus, M. and Imre, M. (2026). Incumbency advantages, prime
minister replacements and government formation. *Party Politics*, 32(2),
183-196.
[doi:10.1177/13540688251344561](https://doi.org/10.1177/13540688251344561)

## Provenance and terms

From the replication archive of Baeck, Debus and Imre (2026) at the
Harvard Dataverse
([doi:10.7910/DVN/CJIVXC](https://doi.org/10.7910/DVN/CJIVXC) , version
1.0, released under CC0 1.0): the data frame `data` of
`baeck_et_al_2025_data.RDATA` (md5 e046b31d8e83dea2e5f01ed498f0cf32).
The residuals `pea_res` and `rej_res` of the authors' two-stage residual
inclusion are rebuilt as their script builds them (logit regressions of
`pea` and `rej` on the other covariates over the potential coalitions of
two or more parties, predicted for every row). Two formations are left
out: Thuringia 2019, at which two cabinets are coded as formed, and
Hesse 2023, at which none is. The cabinet labels are the party slots of
the deposit (`cdu`, `spd`, `fdp`, `gru`, `pds`, `rep`, `dvu`, `afb`,
`sch`, `sta`, `npd`, `ssw`, `pir`, `afd`), and in a few states a slot
holds another party (in Bavaria the CSU is coded `cdu` and the Free
Voters `npd`); the labels identify cabinets within a formation only. The
rebuilding script is `data-raw/make_examples.R` of the development
repository; the conditional logit of Table 1, Model 3 on these rows
agrees with
[`survival::clogit()`](https://rdrr.io/pkg/survival/man/clogit.html)
(log-likelihood -81.0040, which the package's tests check).

## See also

Other datasets:
[`afghan_leaders`](https://bagozzib.github.io/polychoice/reference/afghan_leaders.md),
[`bes1987`](https://bagozzib.github.io/polychoice/reference/bes1987.md),
[`bes2010`](https://bagozzib.github.io/polychoice/reference/bes2010.md),
[`defection`](https://bagozzib.github.io/polychoice/reference/defection.md),
[`fluvax`](https://bagozzib.github.io/polychoice/reference/fluvax.md),
[`pa2006`](https://bagozzib.github.io/polychoice/reference/pa2006.md),
[`pres2004`](https://bagozzib.github.io/polychoice/reference/pres2004.md),
[`sanctions`](https://bagozzib.github.io/polychoice/reference/sanctions.md)

## Examples

``` r
data(cabinets)
table(table(cabinets$election))          # potential cabinets per formation
#> 
#>  3  7  8 15 31 63 
#>  1 21  1 38 44 20 
m <- mnl(chosen ~ 0 + incumbent + incumbent_newpm + minority + min_winning + bargaining +
           strongest + cross_cutting + het_economic + het_social + pea_res + rej_res | 1,
         data = cabinets, id = "election", alt = "cabinet")
summary(m)
#> Conditional logit 
#> Call: mnl(formula = chosen ~ 0 + incumbent + incumbent_newpm + minority + 
#>     min_winning + bargaining + strongest + cross_cutting + het_economic + 
#>     het_social + pea_res + rej_res | 1, data = cabinets, id = "election", 
#>     alt = "cabinet")
#> Alternatives: afb, afd, cdu, cdu+afb, cdu+afd, cdu+dvu, cdu+dvu+ssw, cdu+fdp, ... (296 alternatives)  (base: afb)
#> 
#> Generic coefficients (alternative-specific covariates):
#>                 Estimate Std. Error z value Pr(>|z|)    
#> incumbent         3.7724     0.6666   5.660 1.52e-08 ***
#> incumbent_newpm   2.5570     0.8423   3.036 0.002399 ** 
#> minority         -0.3939     0.9740  -0.404 0.685865    
#> min_winning       3.9542     0.8111   4.875 1.09e-06 ***
#> bargaining        2.0078     0.7941   2.528 0.011463 *  
#> strongest         2.1467     0.6396   3.356 0.000790 ***
#> cross_cutting    -0.7981     0.3386  -2.357 0.018416 *  
#> het_economic     -0.7999     0.1739  -4.599 4.25e-06 ***
#> het_social       -0.5306     0.1795  -2.955 0.003125 ** 
#> pea_res           3.8530     0.7348   5.244 1.57e-07 ***
#> rej_res          -6.4980     1.6873  -3.851 0.000118 ***
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Standard errors: analytic. Log-likelihood -81.004 on 11 parameters, 125 cases. AIC 184.01, BIC 215.12.
```
