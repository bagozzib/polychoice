# Economic sanctions, military force and domestic unrest (Liou 2024)

Country-years 1975-2005 of 137 countries with the government's choice
between no coercion of another state, imposing economic sanctions, and
initiating military force, and the domestic challenges that might push
it to divert attention abroad: mass nonviolent or violent campaigns and
the opposition's seat share (Liou 2024). The rows are the 2,985 complete
cases of the paper's Table 1.

## Usage

``` r
data(sanctions)
```

## Format

A data frame with 2,985 rows and 17 variables (covariates are lagged one
year unless stated):

- outcome:

  Factor: `none` (no sanctions or militarized conflict), `sanctions`
  (imposes economic sanctions, TIES 4.0), `force` (initiates military
  force, a MID of hostility level 4 or 5; a year with both is coded
  `force`).

- ccode:

  Correlates of War country code.

- year:

  Year.

- nonviolent, violent:

  1 if a nonviolent or a violent mass campaign was under way (NAVCO
  2.1).

- opposition:

  Opposition parties' share of the legislative seats (0 to 1; DPI).

- growth:

  GDP growth, percent.

- inflation:

  Inflation (log).

- leg_election, exec_election:

  1 if a legislative or an executive election is held in the following
  year.

- gdppc:

  GDP per capita (log).

- trade:

  Trade openness, percent of GDP (log).

- milex:

  Military expenditure (log).

- democracy:

  V-Dem electoral democracy index.

- human_rights:

  Latent human rights protection score (Fariss).

- coldwar:

  1 during the Cold War.

- peaceyears:

  Years since the last sanction or force initiation.

## Source

Liou, R.Y.-L. (2024). Mass mobilization, elite competition, and
diversionary use of sanctions. *Foreign Policy Analysis*, 20(4),
orae025. [doi:10.1093/fpa/orae025](https://doi.org/10.1093/fpa/orae025)

## Provenance and terms

From the replication archive of Liou (2024) at the Harvard Dataverse
([doi:10.7910/DVN/9QQYPB](https://doi.org/10.7910/DVN/9QQYPB) , version
1.0, released under CC0 1.0): the file `DiversionarySanctions_Main.dta`
(md5 5230d7cda77acc7165ff4e459bbe0508), restricted as
`Analysis_Diversionary Sanctions.do` restricts it for Table 1. Original
names: `outcome` from `diversion`, `nonviolent` from
`lbinonvio_campaign`, `violent` from `lbivio_campaign`, `opposition`
from `loppseatshare`, `growth` from `lgdpgr`, `inflation` from
`llninflation`, `leg_election` from `Flegelec`, `exec_election` from
`Fexelec`, `gdppc` from `llngdppc`, `trade` from `ltradeopen`, `milex`
from `llnmilex`, `democracy` from `lvdem_polyarchy`, `human_rights` from
`lagtheta_mean`, `coldwar`, and `peaceyears` from `yrssince_divert` (the
archive's square and cube of it are the squares and cubes of the bundled
column). The rebuilding script is `data-raw/make_examples.R` of the
development repository. The published log-likelihood (-1161.5576) and
coefficients reproduce; the clustered standard errors are the published
ones divided by the factor \\\sqrt{G/(G-1)}\\ (\\G = 137\\) of Stata's
small-sample correction.

## See also

Other datasets:
[`afghan_leaders`](https://bagozzib.github.io/polychoice/reference/afghan_leaders.md),
[`bes1987`](https://bagozzib.github.io/polychoice/reference/bes1987.md),
[`bes2010`](https://bagozzib.github.io/polychoice/reference/bes2010.md),
[`cabinets`](https://bagozzib.github.io/polychoice/reference/cabinets.md),
[`defection`](https://bagozzib.github.io/polychoice/reference/defection.md),
[`fluvax`](https://bagozzib.github.io/polychoice/reference/fluvax.md),
[`pa2006`](https://bagozzib.github.io/polychoice/reference/pa2006.md),
[`pres2004`](https://bagozzib.github.io/polychoice/reference/pres2004.md)

## Examples

``` r
data(sanctions)
table(sanctions$outcome)
#> 
#>      none sanctions     force 
#>      2400       169       416 
m <- mnl(outcome ~ nonviolent * opposition + violent * opposition + growth + inflation +
           leg_election + exec_election + gdppc + trade + milex + democracy + human_rights +
           coldwar + peaceyears + I(peaceyears^2) + I(peaceyears^3),
         data = sanctions, cluster = sanctions$ccode, se = "cluster")
m$loglik                                # Table 1: -1161.5576
#> [1] -1161.558
coef(m)["sanctions:nonviolent:opposition"]
#> sanctions:nonviolent:opposition 
#>                        8.146855 
```
