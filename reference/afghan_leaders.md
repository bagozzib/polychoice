# Rankings of leader attributes in northern Afghanistan (Bhatia and Monroe 2024)

Respondents of a household survey experiment conducted in 2016-2017 in
Balkh, Kunduz and Sar-e Pol provinces ranked six attributes of political
leaders from most important (1) to least important (6), after hearing
one of four randomly assigned vignettes: none (control), on bribery, on
nepotism, or on insecurity (Bhatia and Monroe 2024). One row per
respondent and attribute (the long layout of
[`rologit()`](https://bagozzib.github.io/polychoice/reference/rologit.md));
2,479 respondents with a complete ranking.

## Usage

``` r
data(afghan_leaders)
```

## Format

A data frame with 14,874 rows (2,479 respondents times 6 attributes) and
9 variables:

- id:

  Respondent identifier.

- attribute:

  Factor: `religious` (strong religious values), `educated` (well
  educated), `mujahideen` (experience as a mujahideen fighter),
  `government` (many years of experience in government), `military`
  (experience in the Afghan military), `honest` (personally honest and
  incorruptible).

- rank:

  1 (most important) to 6 (least important).

- treatment:

  Factor, the randomized vignette: `control`, `bribery`, `nepotism`,
  `insecurity`.

- female:

  1 if female (recorded by the enumerator).

- age:

  Age in years.

- education:

  Factor: `none`, `primary`, `secondary`, `vocational` (post-secondary
  vocational training), `some_university`, `university`, `madrassa`.

- ethnicity:

  Factor: `pashtun`, `tajik`, `uzbek`, `turkmen`, `hazara`, `baloch`,
  `other`.

- province:

  Factor: `Balkh`, `Kunduz`, `Sar-e Pol`.

## Source

Bhatia, J. and Monroe, S.L. (2024). Insecurity and support for female
leadership in conflict states: Evidence from Afghanistan. *British
Journal of Political Science*, 54, 771-792.
[doi:10.1017/S000712342300056X](https://doi.org/10.1017/S000712342300056X)

## Provenance and terms

From the replication archive of Bhatia and Monroe (2024) at the Harvard
Dataverse ([doi:10.7910/DVN/OIV11K](https://doi.org/10.7910/DVN/OIV11K)
, version 1.0, released under CC0 1.0): the file
`Afghan_Data_All_Clean.csv` (md5 8c9d3d11d81efdad4ded7770a1dd6e77), with
the codes of the deposit's codebook. Kept: the four arms of the vignette
experiment (`Treatment` 1-4), the six attribute rankings
`Att_Strong_Rel`, `Att_Well_Ed`, `Att_Muj`, `Att_Gov_Exp`,
`Att_Afg_Mil`, `Att_Hon` where they form a complete ranking without
ties, and five respondent characteristics; the respondent's village,
district and place of birth are not carried over. The rebuilding script
is `data-raw/make_examples.R` of the development repository. The paper
analyses the conjoint experiment of the same survey; these rankings are
not among its models. In the deposit 1 is the most important, as the
codebook defines the question; the paper's code reverses the scale of
the importance items it uses (6 the most important).

## See also

Other datasets:
[`bes1987`](https://bagozzib.github.io/polychoice/reference/bes1987.md),
[`bes2010`](https://bagozzib.github.io/polychoice/reference/bes2010.md),
[`cabinets`](https://bagozzib.github.io/polychoice/reference/cabinets.md),
[`defection`](https://bagozzib.github.io/polychoice/reference/defection.md),
[`fluvax`](https://bagozzib.github.io/polychoice/reference/fluvax.md),
[`pa2006`](https://bagozzib.github.io/polychoice/reference/pa2006.md),
[`pres2004`](https://bagozzib.github.io/polychoice/reference/pres2004.md),
[`sanctions`](https://bagozzib.github.io/polychoice/reference/sanctions.md)

## Examples

``` r
data(afghan_leaders)
## mean rank of each attribute by vignette
round(tapply(afghan_leaders$rank, afghan_leaders[c("attribute", "treatment")], mean), 2)
#>             treatment
#> attribute    control bribery nepotism insecurity
#>   religious     3.06    3.24     3.14       3.11
#>   educated      2.26    2.14     2.23       2.22
#>   mujahideen    4.80    4.95     4.89       4.82
#>   government    4.15    4.10     4.06       4.15
#>   military      4.51    4.43     4.54       4.53
#>   honest        2.22    2.14     2.14       2.17
m <- rologit(rank ~ 1 | treatment + female + age, data = afghan_leaders, id = "id",
             alt = "attribute", ref = "educated")
summary(m)
#> Rank-ordered logit 
#> Call: rologit(formula = rank ~ 1 | treatment + female + age, data = afghan_leaders, 
#>     id = "id", alt = "attribute", ref = "educated")
#> Alternatives: religious, educated, mujahideen, government, military, honest  (base: educated)
#> 
#> Outcome equations (alternative:term):
#>                                  Estimate Std. Error z value Pr(>|z|)    
#> religious:(Intercept)          -0.4133878  0.1136817  -3.636 0.000277 ***
#> religious:treatmentbribery     -0.2170247  0.1035505  -2.096 0.036097 *  
#> religious:treatmentnepotism    -0.0699909  0.1030982  -0.679 0.497216    
#> religious:treatmentinsecurity  -0.0287493  0.0992497  -0.290 0.772072    
#> religious:female               -0.6916090  0.0726527  -9.519  < 2e-16 ***
#> religious:age                   0.0006388  0.0022385   0.285 0.775374    
#> mujahideen:(Intercept)         -1.8133037  0.1282583 -14.138  < 2e-16 ***
#> mujahideen:treatmentbribery    -0.1798406  0.1173997  -1.532 0.125556    
#> mujahideen:treatmentnepotism   -0.0158410  0.1158680  -0.137 0.891256    
#> mujahideen:treatmentinsecurity -0.0336242  0.1126929  -0.298 0.765420    
#> mujahideen:female              -1.2002348  0.0822687 -14.589  < 2e-16 ***
#> mujahideen:age                  0.0090992  0.0025084   3.628 0.000286 ***
#> government:(Intercept)         -1.1633188  0.1171398  -9.931  < 2e-16 ***
#> government:treatmentbribery    -0.0801560  0.1059816  -0.756 0.449457    
#> government:treatmentnepotism    0.0413401  0.1053577   0.392 0.694778    
#> government:treatmentinsecurity -0.0732843  0.1027200  -0.713 0.475575    
#> government:female              -0.6951866  0.0746466  -9.313  < 2e-16 ***
#> government:age                  0.0028092  0.0022989   1.222 0.221716    
#> military:(Intercept)           -1.5135619  0.1194192 -12.674  < 2e-16 ***
#> military:treatmentbribery      -0.0460144  0.1082345  -0.425 0.670737    
#> military:treatmentnepotism      0.0117029  0.1075044   0.109 0.913314    
#> military:treatmentinsecurity   -0.0643341  0.1050230  -0.613 0.540160    
#> military:female                -0.7103013  0.0761995  -9.322  < 2e-16 ***
#> military:age                    0.0056097  0.0023317   2.406 0.016134 *  
#> honest:(Intercept)              0.2698808  0.1121366   2.407 0.016097 *  
#> honest:treatmentbribery         0.0072947  0.1007816   0.072 0.942299    
#> honest:treatmentnepotism        0.0973385  0.1002896   0.971 0.331760    
#> honest:treatmentinsecurity      0.0133651  0.0974709   0.137 0.890937    
#> honest:female                  -0.7525936  0.0711787 -10.573  < 2e-16 ***
#> honest:age                      0.0040639  0.0021823   1.862 0.062573 .  
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> 
#> Standard errors: analytic. Log-likelihood -13593.565 on 30 parameters, 2479 cases. AIC 27247.13, BIC 27421.60.
## do the first choices follow the same preferences as the later ranks?
rank_test(m, after = 1)
#> Stability of preferences across ranking stages (position 1 against later; 2479 ranking(s) with later stages) 
#>              test statistic df p.value
#>  Likelihood ratio   401.151 30  <2e-16
#>              Wald   291.255 30  <2e-16
```
