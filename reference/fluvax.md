# Intended flu vaccination and international collaboration (Liu et al. 2024)

A survey experiment with 1,334 US adults (an opt-in online panel):
respondents were asked whether they would take a flu vaccine developed
by the United States alone, through a US-China collaboration, or through
a US-Japan collaboration, the version randomly assigned (Liu et al.
2024). The answer is no, don't know, or yes; "don't know" may mix
respondents undecided between the other answers with respondents who
give it whatever the vaccine, a candidate inflated alternative for
[`imnl()`](https://bagozzib.github.io/polychoice/reference/imnl.md).

## Usage

``` r
data(fluvax)
```

## Format

A data frame with 1,334 rows and 9 variables:

- vaccinate:

  Factor: `no`, `dont_know`, `yes` (would take the vaccine this year).

- condition:

  Factor, the randomized origin of the vaccine: `us_alone`, `us_china`,
  `us_japan`.

- party:

  Factor: `democrat`, `independent`, `republican` (leaners counted with
  their party), `other` (another party or unclassified).

- age:

  Age in years (81 = older than 80).

- female:

  1 if female (sex assigned at birth).

- educ:

  Education: 1 (8th grade or less), 2 (some high school), 3 (high
  school), 4 (some college), 5 (associate's degree), 7 (bachelor's
  degree), 8 (graduate degree); the survey's codes, used as a score in
  the published model.

- income:

  Household income in 16 brackets, 1 (less than 10,000 dollars) to 16
  (150,000 dollars or more).

- race:

  Factor: `white`, `black` (each alone), `other` (any other or several).

- knowledge:

  Share of four political-knowledge questions answered correctly.

## Source

Liu, R., Nagao, H., Hatungimana, W., Zhang, J.J. and Kennedy, J.J.
(2024). The politics of flu vaccines: International collaboration and
political partisanship. *Japanese Journal of Political Science*, 25,
61-74.
[doi:10.1017/S1468109924000021](https://doi.org/10.1017/S1468109924000021)

## Provenance and terms

From the replication archive of Liu et al. (2024) at the Harvard
Dataverse ([doi:10.7910/DVN/BYAFJR](https://doi.org/10.7910/DVN/BYAFJR)
, version 1.0, released under CC0 1.0): the file `Data.dta` (md5
40cad902c1e17f2de59b6dac0aacc5da), the raw survey items, from which
every variable is rebuilt as the deposit's `DoFile.do` builds it for
Table 3, Model 1 (the codes are those of the deposit's codebook). The
rows are the 1,334 respondents of that model. The rebuilding script is
`data-raw/make_examples.R` of the development repository; with
`se = "robust"` the coefficients and standard errors of Table 3, Model 1
reproduce to the printed three decimals.

## See also

Other datasets:
[`afghan_leaders`](https://bagozzib.github.io/polychoice/reference/afghan_leaders.md),
[`bes1987`](https://bagozzib.github.io/polychoice/reference/bes1987.md),
[`bes2010`](https://bagozzib.github.io/polychoice/reference/bes2010.md),
[`cabinets`](https://bagozzib.github.io/polychoice/reference/cabinets.md),
[`defection`](https://bagozzib.github.io/polychoice/reference/defection.md),
[`pa2006`](https://bagozzib.github.io/polychoice/reference/pa2006.md),
[`pres2004`](https://bagozzib.github.io/polychoice/reference/pres2004.md),
[`sanctions`](https://bagozzib.github.io/polychoice/reference/sanctions.md)

## Examples

``` r
data(fluvax)
table(fluvax$vaccinate, fluvax$condition)
#>            
#>             us_alone us_china us_japan
#>   no             130      154      108
#>   dont_know      103      111      101
#>   yes            216      174      237
m <- mnl(vaccinate ~ condition + age + female + educ + income + race + party + knowledge,
         data = fluvax, se = "robust")
summary(m)
#> Multinomial logit 
#> Call: mnl(formula = vaccinate ~ condition + age + female + educ + income + 
#>     race + party + knowledge, data = fluvax, se = "robust")
#> Alternatives: no, dont_know, yes  (base: no)
#> 
#> Outcome equations (alternative:term):
#>                              Estimate Std. Error z value Pr(>|z|)    
#> dont_know:(Intercept)       -0.193930   0.394889  -0.491 0.623357    
#> dont_know:conditionus_china -0.120920   0.183972  -0.657 0.511007    
#> dont_know:conditionus_japan  0.155645   0.193430   0.805 0.421015    
#> dont_know:age                0.010226   0.004990   2.049 0.040455 *  
#> dont_know:female             0.004734   0.175022   0.027 0.978423    
#> dont_know:educ              -0.034025   0.047230  -0.720 0.471268    
#> dont_know:income             0.023732   0.021535   1.102 0.270462    
#> dont_know:raceblack         -0.023611   0.248868  -0.095 0.924415    
#> dont_know:raceother          0.369651   0.210084   1.760 0.078485 .  
#> dont_know:partyindependent  -0.657632   0.226830  -2.899 0.003741 ** 
#> dont_know:partyrepublican   -0.677829   0.179401  -3.778 0.000158 ***
#> dont_know:partyother        -0.234531   0.445106  -0.527 0.598255    
#> dont_know:knowledge         -0.344941   0.277127  -1.245 0.213240    
#> yes:(Intercept)             -0.610506   0.366698  -1.665 0.095937 .  
#> yes:conditionus_china       -0.384848   0.170969  -2.251 0.024386 *  
#> yes:conditionus_japan        0.346517   0.174194   1.989 0.046673 *  
#> yes:age                      0.009736   0.004584   2.124 0.033666 *  
#> yes:female                  -0.361242   0.160148  -2.256 0.024091 *  
#> yes:educ                     0.194579   0.043071   4.518 6.25e-06 ***
#> yes:income                   0.030262   0.018958   1.596 0.110430    
#> yes:raceblack               -0.357783   0.232926  -1.536 0.124529    
#> yes:raceother                0.226817   0.195169   1.162 0.245173    
#> yes:partyindependent        -1.277758   0.208181  -6.138 8.37e-10 ***
#> yes:partyrepublican         -1.321323   0.166804  -7.921 2.35e-15 ***
#> yes:partyother              -1.329099   0.454544  -2.924 0.003455 ** 
#> yes:knowledge                0.526741   0.242727   2.170 0.030000 *  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Standard errors: robust. Log-likelihood -1281.466 on 26 parameters, 1334 cases. AIC 2614.93, BIC 2750.03.
first_difference(m, "condition", from = "us_alone", to = "us_china")
#> First differences: condition from us_alone to us_china at the covariate profile (delta intervals at 95%)
#>  alternative  from    to   diff  lower  upper
#>           no 0.186 0.238  0.052  0.001  0.102
#>    dont_know 0.207 0.234  0.027 -0.028  0.083
#>          yes 0.607 0.528 -0.079 -0.151 -0.007
# \donttest{
## "don't know" as a mixture of the undecided and of routine non-answerers
mi <- imnl(vaccinate ~ condition + age + female + educ + income + race + party + knowledge |
             age + educ + knowledge, data = fluvax, inflate = "dont_know")
#> The inflated likelihood has 2 distinct local maxima among the ends of the search from 14 starts (log-likelihoods about -1276.9, -1279.6); the highest is reported. The search also ended at 1 point that is not a maximum (a boundary, a saddle, or no convergence). See object$search.
compare_models(MNL = m, inflated = mi)
#>     model
#>       MNL
#>  inflated
#>                                                                       type
#>                                                          Multinomial logit
#>  Inflated multinomial logit (inflated alternative: dont_know, logit split)
#>    inflate    n df    loglik      AIC      BIC mcfadden converged boundary
#>            1334 26 -1281.466 2614.931 2750.026    0.090      TRUE    FALSE
#>  dont_know 1334 30 -1276.870 2613.740 2769.618    0.093      TRUE    FALSE
#>  ill_conditioned n_modes
#>            FALSE       1
#>            FALSE       2
# }
```
