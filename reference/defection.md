# Who leaves insurgent groups: Colombian ex-combatants (Oppenheim et al. 2015; Cook, Niehaus and Zuhlke 2018)

Ex-combatants of Colombian left-wing guerrilla groups surveyed by the
Fundacion Ideas para la Paz, with how they left the group (captured,
individually demobilized, or switched to the paramilitaries), why they
had joined, whether the group indoctrinated them or took food from
peasants by force, and their background (Oppenheim et al. 2015). Rare
outcomes crossed with rare binary covariates leave empty cells (no
captured combatant who joined for economic reasons without
indoctrination, none who joined for ideological reasons in a group that
abused peasants): the maximum-likelihood estimates of the interactions
in Oppenheim et al.'s multinomial logits do not exist, which Cook,
Niehaus and Zuhlke (2018) show, re-estimating the models with Firth's
penalized likelihood (in
[`mnl()`](https://bagozzib.github.io/polychoice/reference/mnl.md),
`firth = TRUE`).

## Usage

``` r
data(defection)
```

## Format

A data frame with 582 rows and 11 variables:

- id:

  Case number, 1 to 582 (the survey's own identifiers are not
  distributed).

- outcome:

  Factor: `captured` (49), `demobilized` (individually demobilized,
  506), `switched` (to the paramilitaries, 27).

- economic:

  1 if the respondent joined for economic reasons.

- ideological:

  1 if the respondent joined for ideological reasons.

- indoctrination:

  1 if the group gave political indoctrination.

- abuse:

  1 if the group took food from peasants by force.

- age_joined:

  Age when the respondent joined the first armed group.

- birth_year:

  Year of birth.

- school_age:

  Age when the respondent left school (0: never attended); Cook et al.'s
  "Education".

- male:

  1 if male.

- besieged:

  Degree to which the group was besieged by the Colombian armed forces,
  0 to 4 (not in the models above).

## Source

Oppenheim, B., Steele, A., Vargas, J.F. and Weintraub, M. (2015). True
believers, deserters, and traitors: Who leaves insurgent groups and why.
*Journal of Conflict Resolution*, 59(5), 794-823.
[doi:10.1177/0022002715576750](https://doi.org/10.1177/0022002715576750)
. Cook, S.J., Niehaus, J. and Zuhlke, S. (2018). A warning on separation
in multinomial logistic models. *Research & Politics*, 5(2).
[doi:10.1177/2053168018769510](https://doi.org/10.1177/2053168018769510)

## Provenance and terms

From the replication archive of Cook, Niehaus and Zuhlke (2018) at the
Harvard Dataverse
([doi:10.7910/DVN/GOYSFK](https://doi.org/10.7910/DVN/GOYSFK) , version
1, released under CC0 1.0): the file `replication data.dta` (md5
c08c38bcb3369ae38940b0a0e1723244), the data of Oppenheim et al. (2015),
restricted as Cook et al.'s replication script restricts them (a
department of residence recorded) to the complete cases of their models.
Original names: `outcome` from `DV_cap_inddemob_switcher`, `economic`
from `join_ec_need`, `ideological` from `join_ideo`, `indoctrination`
from `p409dummy`, `abuse` from `p402g`, `age_joined` from `p104a`,
`birth_year` from `yrbirth`, `school_age` from `p203`, `male` from
`p101`, `besieged` from `p410_recoded`. The rebuilding script is
`data-raw/make_examples.R` of the development repository. The penalized
estimates and standard errors of Cook et al.'s Models 3 to 5 (their
Table 3, computed with brglm2) reproduce.

## See also

Other datasets:
[`afghan_leaders`](https://bagozzib.github.io/polychoice/reference/afghan_leaders.md),
[`bes1987`](https://bagozzib.github.io/polychoice/reference/bes1987.md),
[`bes2010`](https://bagozzib.github.io/polychoice/reference/bes2010.md),
[`cabinets`](https://bagozzib.github.io/polychoice/reference/cabinets.md),
[`fluvax`](https://bagozzib.github.io/polychoice/reference/fluvax.md),
[`pa2006`](https://bagozzib.github.io/polychoice/reference/pa2006.md),
[`pres2004`](https://bagozzib.github.io/polychoice/reference/pres2004.md),
[`sanctions`](https://bagozzib.github.io/polychoice/reference/sanctions.md)

## Examples

``` r
data(defection)
with(defection, table(outcome, economic, indoctrination))   # an empty cell
#> , , indoctrination = 0
#> 
#>              economic
#> outcome         0   1
#>   captured      4   0
#>   demobilized  35  13
#>   switched      7   2
#> 
#> , , indoctrination = 1
#> 
#>              economic
#> outcome         0   1
#>   captured     40   5
#>   demobilized 383  75
#>   switched     11   7
#> 
## Cook et al.'s Model 3: maximum likelihood runs off, the penalized fit does not
m3 <- mnl(outcome ~ economic * indoctrination + age_joined + birth_year + school_age + male,
          data = defection, ref = "captured", firth = TRUE)
#> The maximum-likelihood estimates of demobilized:economic, demobilized:economic:indoctrination, switched:economic, switched:economic:indoctrination diverge (separation); their penalized estimates are determined largely by the penalty. Profile penalized-likelihood intervals (confint()) and penalized likelihood-ratio tests are the recommended inference.
coef(m3)[c("demobilized:economic", "demobilized:economic:indoctrination")]   # 0.903, -0.632
#>                demobilized:economic demobilized:economic:indoctrination 
#>                           0.9031805                          -0.6321131 
m3$separation
#> [1] "demobilized:economic"                "demobilized:economic:indoctrination"
#> [3] "switched:economic"                   "switched:economic:indoctrination"   
```
