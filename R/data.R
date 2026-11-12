#' 2006 Pennsylvania state House vote choice (Arceneaux and Kolodny 2009)
#'
#' Survey respondents in two Pennsylvania state House districts after a field
#' experiment in which a liberal interest group endorsed the Democratic
#' candidate by door-to-door canvassing or by telephone (Arceneaux and Kolodny
#' 2009). The outcome is the reported vote: abstained, Democratic candidate,
#' Republican or other candidate, or nonresponse. Abstention is the base
#' alternative of the published analyses and mixes routine nonvoters with
#' occasional nonvoters, the baseline-inflated outcome of Bagozzi and
#' Marchetti (2017, Replication I), whose Table A.3 these data reproduce. Rows
#' are the 1,998 complete cases of that specification.
#'
#' @section Provenance and terms:
#' Taken from the replication archive of Bagozzi and Marchetti (2017),
#' deposited at the Harvard Dataverse (\doi{10.7910/DVN/3JR2YL}, version 1.1,
#' released under CC0 1.0): the file `Replication I/Least_Informed.dta` (md5
#' 0a0964fc9183e1506e893695000675a2) of `PSRM Replication Files.zip`, which
#' redistributes the Arceneaux and Kolodny (2009) replication data; the
#' rebuilding script
#' `data-raw/make_data.R` of the development repository keeps the analysis
#' variables of the published specification and the complete cases the
#' published table uses. Original variable names: `vote` from
#' `votechoice_dk` (1 abstain, 2 democrat, 0 republican_other, 3 nonresponse),
#' `canvass` from `canvass_all`, `phone` from `phone_only`, `democrat` from
#' `sid_dem`, `republican` from `sid_rep`, `household` from `hhcount`, `age`,
#' `female`, `voted2004` from `vg04`, `district156` from `d156`. The published
#' interaction terms are products of the bundled columns.
#' @format A data frame with 1998 rows and 10 variables:
#' \describe{
#'   \item{vote}{Factor: `abstain`, `democrat`, `republican_other`, `nonresponse`.}
#'   \item{canvass}{1 if the respondent received the canvassing treatment.}
#'   \item{phone}{1 if the respondent received the telephone treatment.}
#'   \item{democrat, republican}{Party identification indicators (independents are the omitted group).}
#'   \item{household}{Number of household members.}
#'   \item{age}{Age in years.}
#'   \item{female}{1 if female.}
#'   \item{voted2004}{1 if the respondent reported voting in 2004.}
#'   \item{district156}{1 if in District 156 (open seat), 0 if in District 161.}
#' }
#' @source Arceneaux, K. and Kolodny, R. (2009). Educating the least informed:
#'   Group endorsements in a grassroots campaign. \emph{American Journal of
#'   Political Science}, 53, 755-770. As analyzed in Bagozzi, B.E. and
#'   Marchetti, K. (2017). Distinguishing occasional abstention from routine
#'   indifference in models of vote choice. \emph{Political Science Research
#'   and Methods}, 5(2), 277-294.
#' @docType data
#' @name pa2006
#' @usage data(pa2006)
#' @keywords datasets
#' @family datasets
#' @examples
#' data(pa2006)
#' table(pa2006$vote)
#' f <- vote ~ canvass + phone + democrat + republican + canvass:democrat + canvass:republican +
#'   phone:democrat + phone:republican + household + age + female + voted2004 + district156
#' m0 <- mnl(f, data = pa2006)
#' m0$loglik                                       # published -2110.59
#' ## the published BIMNL (Table A.3) is fitted in vignette("inflation")
NULL

#' 2004 presidential vote choice (Campbell and Monson 2008)
#'
#' Respondents of the 2004 Election Panel Study used by Campbell and Monson
#' (2008) to study the effect of state gay-marriage-ban ballot initiatives on
#' the 2004 presidential vote. The outcome is abstained, voted for Kerry, or
#' voted for Bush; abstention is the base alternative and, as the mixture of
#' routine and occasional nonvoters, the inflated alternative of Bagozzi and
#' Marchetti (2017, Replication II), whose Table A.4 these data reproduce. Rows
#' are the 1,341 complete cases of that specification.
#'
#' @section Provenance and terms:
#' Taken from the replication archive of Bagozzi and Marchetti (2017),
#' deposited at the Harvard Dataverse (\doi{10.7910/DVN/3JR2YL}, version 1.1,
#' released under CC0 1.0): the file `Replication II/replicationdata.dta` (md5
#' 29ee7e0d84871d898fb19e1450b5c9a0) of `PSRM Replication Files.zip`, which
#' redistributes the Campbell and Monson (2008) replication data; the
#' rebuilding script
#' `data-raw/make_data.R` of the development repository keeps the analysis
#' variables of the published specification and the complete cases the
#' published table uses. Original variable names: `vote` from `vote_turn` (0
#' abstain, 1 kerry, 2 bush), `educ`, `female`, `black`, `hispanic`, `pid7`
#' from `party7`, `bush_values` from `w3gbvalu2`, `battleground` from
#' `presbatt`, `south`, `gmb` from `gayban2`, `agegroup` from `agegroup2`,
#' `evangelical` from `bornagn_w`, `catholic`, `other_religion` from
#' `other_rel`, `secular`, `ideology` from `ideo`, `mobilization` from
#' `w3mobidx`. The published interaction terms are products of the bundled
#' columns. The log-likelihoods (the inflated model's -545.279 agrees with the
#' printed -545.27 when truncated) and the split-equation coefficients of
#' Table A.4 reproduce (and the outcome-stage coefficients the package's tests
#' check); the raw Vuong statistic of the inflated against the plain model,
#' `vuong(imnl_fit, mnl_fit)`, is 1.68 on these data (positive favours the
#' inflated model; two-sided p = 0.09).
#' @format A data frame with 1341 rows and 17 variables:
#' \describe{
#'   \item{vote}{Factor: `abstain`, `kerry`, `bush`.}
#'   \item{educ}{Education, 1 (elementary) to 9 (advanced graduate work).}
#'   \item{female, black, hispanic}{Indicators.}
#'   \item{pid7}{Party identification, 1 (strong Republican) to 7 (strong Democrat).}
#'   \item{bush_values}{"Shares my values" describes Bush: 1 (not well at all) to 4 (extremely well).}
#'   \item{battleground}{1 if in a presidential battleground state.}
#'   \item{south}{1 if in a Southern state.}
#'   \item{gmb}{1 if the state had a gay-marriage-ban initiative on the 2004 ballot.}
#'   \item{agegroup}{Age cohort: 1 (18-29), 2 (30-54), 3 (55 and over).}
#'   \item{evangelical}{1 if white born-again / evangelical Christian.}
#'   \item{catholic, other_religion, secular}{Religious affiliation indicators.}
#'   \item{ideology}{1 liberal, 2 moderate, 3 conservative.}
#'   \item{mobilization}{Additive index (0-5) of campaign contacts.}
#' }
#' @source Campbell, D.E. and Monson, J.Q. (2008). The religion card: Gay
#'   marriage and the 2004 presidential election. \emph{Public Opinion
#'   Quarterly}, 72, 399-419. As analyzed in Bagozzi, B.E. and Marchetti, K.
#'   (2017). Distinguishing occasional abstention from routine indifference in
#'   models of vote choice. \emph{Political Science Research and Methods},
#'   5(2), 277-294.
#' @docType data
#' @name pres2004
#' @usage data(pres2004)
#' @keywords datasets
#' @family datasets
#' @examples
#' data(pres2004)
#' f <- vote ~ educ + female + black + hispanic + pid7 + bush_values + battleground + south + gmb +
#'   agegroup + evangelical + evangelical:gmb + catholic + catholic:gmb + other_religion +
#'   secular + secular:gmb + ideology + mobilization
#' m0 <- mnl(f, data = pres2004)
#' \donttest{
#' m1 <- imnl(vote ~ educ + female + black + hispanic + pid7 + bush_values + battleground +
#'   south + gmb + agegroup + evangelical + evangelical:gmb + catholic + catholic:gmb +
#'   other_religion + secular + secular:gmb + ideology + mobilization |
#'   educ + agegroup + mobilization + secular, data = pres2004)
#' vuong(m1, m0)
#' }
NULL

#' Rankings of leader attributes in northern Afghanistan (Bhatia and Monroe 2024)
#'
#' Respondents of a household survey experiment conducted in 2016-2017 in
#' Balkh, Kunduz and Sar-e Pol provinces ranked six attributes of political
#' leaders from most important (1) to least important (6), after hearing one
#' of four randomly assigned vignettes: none (control), on bribery, on
#' nepotism, or on insecurity (Bhatia and Monroe 2024). One row per respondent
#' and attribute (the long layout of [rologit()]); 2,479 respondents with a
#' complete ranking.
#'
#' @section Provenance and terms:
#' From the replication archive of Bhatia and Monroe (2024) at the Harvard
#' Dataverse (\doi{10.7910/DVN/OIV11K}, version 1.0, released under CC0 1.0):
#' the file `Afghan_Data_All_Clean.csv` (md5 8c9d3d11d81efdad4ded7770a1dd6e77),
#' with the codes of the deposit's codebook. Kept: the four arms of the
#' vignette experiment (`Treatment` 1-4), the six attribute rankings
#' `Att_Strong_Rel`, `Att_Well_Ed`, `Att_Muj`, `Att_Gov_Exp`, `Att_Afg_Mil`,
#' `Att_Hon` where they form a complete ranking without ties, and five
#' respondent characteristics; the respondent's village, district and place of
#' birth are not carried over. The rebuilding script is
#' `data-raw/make_examples.R` of the development repository. The paper
#' analyses the conjoint experiment of the same survey; these rankings are
#' not among its models. In the deposit 1 is the most important, as the
#' codebook defines the question; the paper's code reverses the scale of the
#' importance items it uses (6 the most important).
#' @format A data frame with 14,874 rows (2,479 respondents times 6 attributes) and 9 variables:
#' \describe{
#'   \item{id}{Respondent identifier.}
#'   \item{attribute}{Factor: `religious` (strong religious values), `educated` (well educated), `mujahideen` (experience as a mujahideen fighter),
#'     `government` (many years of experience in government), `military` (experience in the Afghan military), `honest` (personally honest and incorruptible).}
#'   \item{rank}{1 (most important) to 6 (least important).}
#'   \item{treatment}{Factor, the randomized vignette: `control`, `bribery`, `nepotism`, `insecurity`.}
#'   \item{female}{1 if female (recorded by the enumerator).}
#'   \item{age}{Age in years.}
#'   \item{education}{Factor: `none`, `primary`, `secondary`, `vocational` (post-secondary vocational training), `some_university`, `university`, `madrassa`.}
#'   \item{ethnicity}{Factor: `pashtun`, `tajik`, `uzbek`, `turkmen`, `hazara`, `baloch`, `other`.}
#'   \item{province}{Factor: `Balkh`, `Kunduz`, `Sar-e Pol`.}
#' }
#' @source Bhatia, J. and Monroe, S.L. (2024). Insecurity and support for
#'   female leadership in conflict states: Evidence from Afghanistan.
#'   \emph{British Journal of Political Science}, 54, 771-792.
#'   \doi{10.1017/S000712342300056X}
#' @docType data
#' @name afghan_leaders
#' @usage data(afghan_leaders)
#' @keywords datasets
#' @family datasets
#' @examples
#' data(afghan_leaders)
#' ## mean rank of each attribute by vignette
#' round(tapply(afghan_leaders$rank, afghan_leaders[c("attribute", "treatment")], mean), 2)
#' m <- rologit(rank ~ 1 | treatment + female + age, data = afghan_leaders, id = "id",
#'              alt = "attribute", ref = "educated")
#' summary(m)
#' ## do the first choices follow the same preferences as the later ranks?
#' rank_test(m, after = 1)
NULL

#' Intended flu vaccination and international collaboration (Liu et al. 2024)
#'
#' A survey experiment with 1,334 US adults (an opt-in online panel):
#' respondents were asked whether they would take a flu vaccine developed by
#' the United States alone, through a US-China collaboration, or through a
#' US-Japan collaboration, the version randomly assigned (Liu et al. 2024).
#' The answer is no, don't know, or yes; "don't know" may mix respondents
#' undecided between the other answers with respondents who give it whatever
#' the vaccine, a candidate inflated alternative for [imnl()].
#'
#' @section Provenance and terms:
#' From the replication archive of Liu et al. (2024) at the Harvard Dataverse
#' (\doi{10.7910/DVN/BYAFJR}, version 1.0, released under CC0 1.0): the file
#' `Data.dta` (md5 40cad902c1e17f2de59b6dac0aacc5da), the raw survey items,
#' from which every variable is rebuilt as the deposit's `DoFile.do` builds it
#' for Table 3, Model 1 (the codes are those of the deposit's codebook). The
#' rows are the 1,334 respondents of that model. The rebuilding script is
#' `data-raw/make_examples.R` of the development repository; with
#' `se = "robust"` the coefficients and standard errors of Table 3, Model 1
#' reproduce to the printed three decimals.
#' @format A data frame with 1,334 rows and 9 variables:
#' \describe{
#'   \item{vaccinate}{Factor: `no`, `dont_know`, `yes` (would take the vaccine this year).}
#'   \item{condition}{Factor, the randomized origin of the vaccine: `us_alone`, `us_china`, `us_japan`.}
#'   \item{party}{Factor: `democrat`, `independent`, `republican` (leaners counted with their party), `other` (another party or unclassified).}
#'   \item{age}{Age in years (81 = older than 80).}
#'   \item{female}{1 if female (sex assigned at birth).}
#'   \item{educ}{Education: 1 (8th grade or less), 2 (some high school), 3 (high school), 4 (some college), 5 (associate's degree),
#'     7 (bachelor's degree), 8 (graduate degree); the survey's codes, used as a score in the published model.}
#'   \item{income}{Household income in 16 brackets, 1 (less than 10,000 dollars) to 16 (150,000 dollars or more).}
#'   \item{race}{Factor: `white`, `black` (each alone), `other` (any other or several).}
#'   \item{knowledge}{Share of four political-knowledge questions answered correctly.}
#' }
#' @source Liu, R., Nagao, H., Hatungimana, W., Zhang, J.J. and Kennedy, J.J.
#'   (2024). The politics of flu vaccines: International collaboration and
#'   political partisanship. \emph{Japanese Journal of Political Science}, 25,
#'   61-74. \doi{10.1017/S1468109924000021}
#' @docType data
#' @name fluvax
#' @usage data(fluvax)
#' @keywords datasets
#' @family datasets
#' @examples
#' data(fluvax)
#' table(fluvax$vaccinate, fluvax$condition)
#' m <- mnl(vaccinate ~ condition + age + female + educ + income + race + party + knowledge,
#'          data = fluvax, se = "robust")
#' summary(m)
#' first_difference(m, "condition", from = "us_alone", to = "us_china")
#' \donttest{
#' ## "don't know" as a mixture of the undecided and of routine non-answerers
#' mi <- imnl(vaccinate ~ condition + age + female + educ + income + race + party + knowledge |
#'              age + educ + knowledge, data = fluvax, inflate = "dont_know")
#' compare_models(MNL = m, inflated = mi)
#' }
NULL

#' Government formation in the German states, 1990-2023 (Baeck, Debus and Imre 2026)
#'
#' Every potential cabinet (every combination of the parties in parliament)
#' at 125 government formations after German state elections, with the
#' cabinet that formed (Baeck, Debus and Imre 2026). One row per formation and
#' potential cabinet: the long layout of a conditional logit whose
#' alternatives are unlabelled and differ from one formation to the next (3
#' to 63 potential cabinets), so the model has generic coefficients and no
#' alternative-specific constants.
#'
#' @section Provenance and terms:
#' From the replication archive of Baeck, Debus and Imre (2026) at the Harvard
#' Dataverse (\doi{10.7910/DVN/CJIVXC}, version 1.0, released under CC0 1.0):
#' the data frame `data` of `baeck_et_al_2025_data.RDATA` (md5
#' e046b31d8e83dea2e5f01ed498f0cf32). The residuals `pea_res` and `rej_res`
#' of the authors' two-stage residual inclusion are rebuilt as their script
#' builds them (logit regressions of `pea` and `rej` on the other covariates
#' over the potential coalitions of two or more parties, predicted for every
#' row). Two formations are left out: Thuringia 2019, at which two cabinets are
#' coded as formed, and Hesse 2023, at which none is. The cabinet labels are
#' the party slots of the deposit (`cdu`, `spd`, `fdp`, `gru`, `pds`, `rep`,
#' `dvu`, `afb`, `sch`, `sta`, `npd`, `ssw`, `pir`, `afd`), and in a few states
#' a slot holds another party (in Bavaria the CSU is coded `cdu` and the Free
#' Voters `npd`); the labels identify cabinets within a formation only. The
#' rebuilding script is `data-raw/make_examples.R` of the development
#' repository; the conditional logit of Table 1, Model 3 on these rows agrees
#' with `survival::clogit()` (log-likelihood -81.0040, which the package's
#' tests check).
#' @format A data frame with 3,352 rows and 18 variables:
#' \describe{
#'   \item{election}{Formation opportunity (state and election year).}
#'   \item{state}{German state.}
#'   \item{year}{Election year.}
#'   \item{cabinet}{Factor, the potential cabinet's parties (see above).}
#'   \item{chosen}{1 for the cabinet that formed.}
#'   \item{incumbent}{1 if the cabinet is the incumbent coalition and the prime minister did not change.}
#'   \item{incumbent_newpm}{1 if the cabinet is the incumbent coalition with a replaced prime minister.}
#'   \item{minority}{1 for a minority cabinet.}
#'   \item{min_winning}{1 for a minimal winning coalition.}
#'   \item{bargaining}{1 for the bargaining-proposition cabinet (a minimal winning coalition with the fewest parties).}
#'   \item{strongest}{1 if the cabinet includes the largest party.}
#'   \item{cross_cutting}{1 for a cross-cutting coalition.}
#'   \item{het_economic, het_social}{Policy heterogeneity within the cabinet on the economic and the sociocultural dimension.}
#'   \item{pea}{1 if the parties had a pre-electoral alliance.}
#'   \item{rej}{1 if one of the parties had ruled out governing with another before the election.}
#'   \item{pea_res, rej_res}{First-stage residuals of `pea` and `rej` (two-stage residual inclusion).}
#' }
#' @source Baeck, H., Debus, M. and Imre, M. (2026). Incumbency advantages,
#'   prime minister replacements and government formation. \emph{Party
#'   Politics}, 32(2), 183-196. \doi{10.1177/13540688251344561}
#' @docType data
#' @name cabinets
#' @usage data(cabinets)
#' @keywords datasets
#' @family datasets
#' @examples
#' data(cabinets)
#' table(table(cabinets$election))          # potential cabinets per formation
#' m <- mnl(chosen ~ 0 + incumbent + incumbent_newpm + minority + min_winning + bargaining +
#'            strongest + cross_cutting + het_economic + het_social + pea_res + rej_res | 1,
#'          data = cabinets, id = "election", alt = "cabinet")
#' summary(m)
NULL

#' Economic sanctions, military force and domestic unrest (Liou 2024)
#'
#' Country-years 1975-2005 of 137 countries with the government's choice
#' between no coercion of another state, imposing economic sanctions, and
#' initiating military force, and the domestic challenges that might push it
#' to divert attention abroad: mass nonviolent or violent campaigns and the
#' opposition's seat share (Liou 2024). The rows are the 2,985 complete cases
#' of the paper's Table 1.
#'
#' @section Provenance and terms:
#' From the replication archive of Liou (2024) at the Harvard Dataverse
#' (\doi{10.7910/DVN/9QQYPB}, version 1.0, released under CC0 1.0): the file
#' `DiversionarySanctions_Main.dta` (md5 5230d7cda77acc7165ff4e459bbe0508),
#' restricted as `Analysis_Diversionary Sanctions.do` restricts it for Table 1.
#' Original names: `outcome` from `diversion`, `nonviolent` from
#' `lbinonvio_campaign`, `violent` from `lbivio_campaign`, `opposition` from
#' `loppseatshare`, `growth` from `lgdpgr`, `inflation` from `llninflation`,
#' `leg_election` from `Flegelec`, `exec_election` from `Fexelec`, `gdppc` from
#' `llngdppc`, `trade` from `ltradeopen`, `milex` from `llnmilex`, `democracy`
#' from `lvdem_polyarchy`, `human_rights` from `lagtheta_mean`, `coldwar`, and
#' `peaceyears` from `yrssince_divert` (the archive's square and cube of it
#' are the squares and cubes of the bundled column). The rebuilding script is
#' `data-raw/make_examples.R` of the development repository. The published
#' log-likelihood (-1161.5576) and coefficients reproduce; the clustered
#' standard errors are the published ones divided by the factor
#' \eqn{\sqrt{G/(G-1)}} (\eqn{G = 137}) of Stata's small-sample correction.
#' @format A data frame with 2,985 rows and 17 variables (covariates are lagged one year unless stated):
#' \describe{
#'   \item{outcome}{Factor: `none` (no sanctions or militarized conflict), `sanctions` (imposes economic sanctions, TIES 4.0),
#'     `force` (initiates military force, a MID of hostility level 4 or 5; a year with both is coded `force`).}
#'   \item{ccode}{Correlates of War country code.}
#'   \item{year}{Year.}
#'   \item{nonviolent, violent}{1 if a nonviolent or a violent mass campaign was under way (NAVCO 2.1).}
#'   \item{opposition}{Opposition parties' share of the legislative seats (0 to 1; DPI).}
#'   \item{growth}{GDP growth, percent.}
#'   \item{inflation}{Inflation (log).}
#'   \item{leg_election, exec_election}{1 if a legislative or an executive election is held in the following year.}
#'   \item{gdppc}{GDP per capita (log).}
#'   \item{trade}{Trade openness, percent of GDP (log).}
#'   \item{milex}{Military expenditure (log).}
#'   \item{democracy}{V-Dem electoral democracy index.}
#'   \item{human_rights}{Latent human rights protection score (Fariss).}
#'   \item{coldwar}{1 during the Cold War.}
#'   \item{peaceyears}{Years since the last sanction or force initiation.}
#' }
#' @source Liou, R.Y.-L. (2024). Mass mobilization, elite competition, and
#'   diversionary use of sanctions. \emph{Foreign Policy Analysis}, 20(4),
#'   orae025. \doi{10.1093/fpa/orae025}
#' @docType data
#' @name sanctions
#' @usage data(sanctions)
#' @keywords datasets
#' @family datasets
#' @examples
#' data(sanctions)
#' table(sanctions$outcome)
#' m <- mnl(outcome ~ nonviolent * opposition + violent * opposition + growth + inflation +
#'            leg_election + exec_election + gdppc + trade + milex + democracy + human_rights +
#'            coldwar + peaceyears + I(peaceyears^2) + I(peaceyears^3),
#'          data = sanctions, cluster = sanctions$ccode, se = "cluster")
#' m$loglik                                # Table 1: -1161.5576
#' coef(m)["sanctions:nonviolent:opposition"]
NULL

#' Vote choice and issue distances in the 1987 British general election (Glasgow 2001)
#'
#' Respondents of the 1987 British General Election Study who voted
#' Conservative, Labour or Alliance, with their distance to each party on
#' seven issues and their social and economic characteristics (Glasgow 2001,
#' following Alvarez, Nagler and Bowler 2000). One row per respondent and party
#' (the long layout); the classic example of the IIA debate in political
#' science and of the mixed logit as its remedy.
#'
#' @section Provenance and terms:
#' From the replication archive of Glasgow (2001) at the Harvard Dataverse
#' (\doi{10.7910/DVN/PCPVWB}, version 1.0, released under CC0 1.0): the file
#' `glasgow.data.asc` of `glasgow.web.zip` (md5
#' ab0c27430d5b281eed1859ebb88b64bc), its 39 columns read in the order of the
#' archive's GAUSS program. The survey is the British General Election Study,
#' 1987: Cross-Section Survey (Heath, Jowell and Curtice; UK Data Service SN
#' 2568, \doi{10.5255/UKDA-SN-2568-1}), which the UK Data Service distributes
#' under the Creative Commons Attribution 4.0 licence; use of these data should
#' cite it. The rebuilding script is `data-raw/make_examples.R` of the
#' development repository. The multinomial logit of the paper's Appendix A,
#' Table C (log-likelihood -1477.56) reproduces.
#' @format A data frame with 6,393 rows (2,131 respondents times 3 parties) and 26 variables:
#' \describe{
#'   \item{id}{Respondent identifier.}
#'   \item{party}{Factor: `conservative`, `labour`, `alliance`.}
#'   \item{chosen}{1 for the party voted for.}
#'   \item{defense, phillips, taxation, nationalization, redistribution, crime, welfare}{Distance between the respondent and the party on the issue:
#'     the absolute difference between the respondent's self-placement and the party's mean placement by all respondents, on 11-point scales
#'     (defense: nuclear weapons; phillips: inflation versus unemployment; crime: civil rights versus fighting crime).}
#'   \item{south, midlands, north, wales, scotland}{Region of residence.}
#'   \item{union}{1 if a member of a trade union or staff association.}
#'   \item{public}{1 if employed in the public sector.}
#'   \item{working_class}{1 if in a working-class (manual, personal service or agricultural) occupation.}
#'   \item{female}{1 if female.}
#'   \item{age}{Age in years divided by 10.}
#'   \item{homeowner}{1 if the household owns its residence.}
#'   \item{income}{Family income, 1 to 12.}
#'   \item{education}{Age at the end of continuous full-time education, divided by 10.}
#'   \item{inflation, taxes, unemployment}{Perceived change since the 1983 election in prices, taxes and unemployment, 1 (increased a lot) to 5 (fallen a lot).}
#' }
#' @source Glasgow, G. (2001). Mixed logit models for multiparty elections.
#'   \emph{Political Analysis}, 9, 116-136.
#'   Data: Heath, A., Jowell, R. and Curtice, J.K. British General Election
#'   Study, 1987: Cross-Section Survey. UK Data Service, SN 2568.
#'   \doi{10.5255/UKDA-SN-2568-1}
#' @docType data
#' @name bes1987
#' @usage data(bes1987)
#' @keywords datasets
#' @family datasets
#' @examples
#' data(bes1987)
#' f <- chosen ~ defense + phillips + taxation + nationalization + redistribution + crime + welfare |
#'   south + midlands + north + wales + scotland + union + public + working_class + female + age +
#'   homeowner + income + education + inflation + taxes + unemployment
#' m <- mnl(f, data = bes1987, id = "id", alt = "party", ref = "alliance")
#' m$loglik                                # Appendix A, Table C: -1477.56
#' \donttest{
#' iia_tests(m, seed = 1)
#' }
NULL

#' Candidate-voter distance and vote choice in England, 2010 (Arzheimer and Evans 2012)
#'
#' Respondents of the 2010 British Election Study face-to-face survey in
#' English constituencies, with the driving distance from their home to each
#' party's candidate, whether the candidate was the sitting MP, the difference
#' in local deprivation between the respondent's and the candidate's ward, and
#' their feelings about the parties before and after the election (Arzheimer
#' and Evans 2012). One row per respondent and party (Labour, Conservative,
#' Liberal Democrat; a party whose candidate's distance is unknown is left
#' out, so the choice sets vary).
#'
#' @section Provenance and terms:
#' From the replication archive of Arzheimer and Evans (2012) at the Harvard
#' Dataverse (\code{hdl:1902.1/17940}, version 1.0, released under CC0 1.0):
#' the file `replication.dta` (md5 d82e986aef3f1389c6ee3babb0a15566). The rows
#' are those of the archive's fourth conditional logit (`replication.do`, line
#' 52), dropped as Stata's `clogit` drops them; the respondents' ward
#' identifiers and the candidates' postcodes are not carried over. The
#' survey is the British Election Study, 2010: Face-to-Face Survey (Whiteley
#' and Sanders; UK Data Service SN 7529, \doi{10.5255/UKDA-SN-7529-1}), which
#' the UK Data Service distributes under the Creative Commons Attribution 4.0
#' licence; use of these data should cite it. The rebuilding script is
#' `data-raw/make_examples.R` of the development repository.
#' @format A data frame with 1,788 rows (674 respondents) and 9 variables:
#' \describe{
#'   \item{id}{Respondent identifier.}
#'   \item{party}{Factor: `labour`, `conservative`, `libdem`.}
#'   \item{chosen}{1 for the party voted for.}
#'   \item{feeling}{Feeling about the party before the election, 0 to 10.}
#'   \item{incumbent}{1 if the party's candidate was the sitting MP.}
#'   \item{distance}{Driving distance between the respondent's home and the candidate's.}
#'   \item{deprivation}{Absolute difference in the Index of Multiple Deprivation 2010 between the respondent's and the candidate's ward.}
#'   \item{rating}{Feeling about the party after the election, 0 to 10 (`NA` if not asked or not answered).}
#'   \item{constituency}{Constituency name.}
#' }
#' @source Arzheimer, K. and Evans, J. (2012). Geolocation and voting:
#'   Candidate-voter distance effects on party choice in the 2010 UK general
#'   election in England. \emph{Political Geography}, 31, 301-310.
#'   \doi{10.1016/j.polgeo.2012.03.004}. Data: Whiteley, P.F. and Sanders, D.
#'   British Election Study, 2010: Face-to-Face Survey. UK Data Service, SN
#'   7529. \doi{10.5255/UKDA-SN-7529-1}
#' @docType data
#' @name bes2010
#' @usage data(bes2010)
#' @keywords datasets
#' @family datasets
#' @examples
#' data(bes2010)
#' table(table(bes2010$id))                # parties in each respondent's choice set
#' m <- mnl(chosen ~ feeling + incumbent + distance + deprivation | 1, data = bes2010, id = "id",
#'          alt = "party", ref = "labour", cluster = "constituency", se = "cluster")
#' summary(m)
#' ## the post-election ratings as rankings (equal ratings tie)
#' r <- bes2010[!is.na(bes2010$rating), ]
#' r$rank <- stats::ave(-r$rating, r$id, FUN = function(v) match(v, sort(unique(v))))
#' keep <- stats::ave(r$rank, r$id, FUN = function(v) length(v) >= 2 & length(unique(v)) >= 2) == 1
#' mr <- rologit(rank ~ incumbent + distance + deprivation | 1, data = r[keep, ], id = "id",
#'               alt = "party", ref = "labour")
#' summary(mr)
NULL

#' Who leaves insurgent groups: Colombian ex-combatants (Oppenheim et al. 2015; Cook, Niehaus and Zuhlke 2018)
#'
#' Ex-combatants of Colombian left-wing guerrilla groups surveyed by the
#' Fundacion Ideas para la Paz, with how they left the group (captured,
#' individually demobilized, or switched to the paramilitaries), why they had
#' joined, whether the group indoctrinated them or took food from peasants by
#' force, and their background (Oppenheim et al. 2015). Rare outcomes crossed
#' with rare binary covariates leave empty cells (no captured combatant who
#' joined for economic reasons without indoctrination, none who joined for
#' ideological reasons in a group that abused peasants): the
#' maximum-likelihood estimates of the interactions in Oppenheim et al.'s
#' multinomial logits do not exist, which Cook, Niehaus and Zuhlke (2018)
#' show, re-estimating the models with Firth's penalized likelihood (in
#' [mnl()], `firth = TRUE`).
#'
#' @section Provenance and terms:
#' From the replication archive of Cook, Niehaus and Zuhlke (2018) at the
#' Harvard Dataverse (\doi{10.7910/DVN/GOYSFK}, version 1, released under CC0
#' 1.0): the file `replication data.dta` (md5
#' c08c38bcb3369ae38940b0a0e1723244), the data of Oppenheim et al. (2015),
#' restricted as Cook et al.'s replication script restricts them (a department
#' of residence recorded) to the complete cases of their models. Original
#' names: `outcome` from `DV_cap_inddemob_switcher`,
#' `economic` from `join_ec_need`, `ideological` from `join_ideo`,
#' `indoctrination` from `p409dummy`, `abuse` from `p402g`, `age_joined` from
#' `p104a`, `birth_year` from `yrbirth`, `school_age` from `p203`, `male` from
#' `p101`, `besieged` from `p410_recoded`. The rebuilding script is
#' `data-raw/make_examples.R` of the development repository. The penalized
#' estimates and standard errors of Cook et al.'s Models 3 to 5 (their Table 3,
#' computed with brglm2) reproduce.
#' @format A data frame with 582 rows and 11 variables:
#' \describe{
#'   \item{id}{Case number, 1 to 582 (the survey's own identifiers are not distributed).}
#'   \item{outcome}{Factor: `captured` (49), `demobilized` (individually demobilized, 506), `switched` (to the paramilitaries, 27).}
#'   \item{economic}{1 if the respondent joined for economic reasons.}
#'   \item{ideological}{1 if the respondent joined for ideological reasons.}
#'   \item{indoctrination}{1 if the group gave political indoctrination.}
#'   \item{abuse}{1 if the group took food from peasants by force.}
#'   \item{age_joined}{Age when the respondent joined the first armed group.}
#'   \item{birth_year}{Year of birth.}
#'   \item{school_age}{Age when the respondent left school (0: never attended); Cook et al.'s "Education".}
#'   \item{male}{1 if male.}
#'   \item{besieged}{Degree to which the group was besieged by the Colombian armed forces, 0 to 4 (not in the models above).}
#' }
#' @source Oppenheim, B., Steele, A., Vargas, J.F. and Weintraub, M. (2015).
#'   True believers, deserters, and traitors: Who leaves insurgent groups and
#'   why. \emph{Journal of Conflict Resolution}, 59(5), 794-823.
#'   \doi{10.1177/0022002715576750}. Cook, S.J., Niehaus, J. and Zuhlke, S.
#'   (2018). A warning on separation in multinomial logistic models.
#'   \emph{Research & Politics}, 5(2). \doi{10.1177/2053168018769510}
#' @docType data
#' @name defection
#' @usage data(defection)
#' @keywords datasets
#' @family datasets
#' @examples
#' data(defection)
#' with(defection, table(outcome, economic, indoctrination))   # an empty cell
#' ## Cook et al.'s Model 3: maximum likelihood runs off, the penalized fit does not
#' m3 <- mnl(outcome ~ economic * indoctrination + age_joined + birth_year + school_age + male,
#'           data = defection, ref = "captured", firth = TRUE)
#' coef(m3)[c("demobilized:economic", "demobilized:economic:indoctrination")]   # 0.903, -0.632
#' m3$separation
NULL
