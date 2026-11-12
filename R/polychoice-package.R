#' polychoice: Multinomial Choice, Inflated Multinomial Logit, and Tests of IIA
#'
#' Estimation, inference, quantities of interest, and specification tests for
#' unordered polytomous (multinomial) choice outcomes, behind one interface and
#' one class of fitted objects (`"choicemodel"`).
#'
#' @section Estimators:
#' [mnl()] fits the multinomial / conditional logit with case-specific
#' covariates (alternative-specific coefficients), alternative-specific
#' covariates (generic coefficients), alternative-specific constants, varying
#' choice sets, frequency weights, offsets, unit fixed effects (with the
#' split-panel jackknife `bias_correct = "jackknife"`), and analytic,
#' outer-product, robust, cluster-robust, or bootstrap standard errors.
#' [imnl()] fits the baseline-inflated multinomial logit of Bagozzi and
#' Marchetti (2017) for any inflated alternative, with a logit or probit split
#' equation, and with several inflated alternatives the parametrized logit
#' captivity model of Swait and Ben-Akiva (1987), which generalizes the dogit
#' of Gaudry and Dagenais (1979). [nlogit()] fits the two-level nested logit, [hev()] the
#' heteroskedastic extreme value model, [mnprobit()] the multinomial probit
#' with a free covariance, and [remnl()] the multinomial logit with random
#' alternative-specific intercepts per unit; these accept a wide layout (one row
#' per case) and a long layout (one row per case-alternative pair; see the `id`
#' and `alt` arguments). [rologit()] fits the rank-ordered (exploded) logit to
#' complete, partial, or tied rankings in the long layout.
#'
#' @section Separation:
#' [separation()] checks exactly, by a linear program, whether the data of an
#' [mnl()] or [rologit()] fit are separated (no finite maximum-likelihood
#' estimate). The maximum-likelihood likelihood-ratio tests then compare the
#' suprema of the log-likelihoods. `firth = TRUE` (Firth's penalty) and
#' `prior = coef_prior(...)` (normal, Cauchy or t priors) give finite
#' penalized estimates, with profile intervals and penalized likelihood-ratio
#' tests.
#'
#' @section Quantities of interest and comparison:
#' [predict.choicemodel()], [first_difference()], [rrr()], [fit_stats()],
#' [classification()], [lr_test()], [vuong()], [boot_lr_test()],
#' [compare_models()], [mundlak()], the simulators [rmnl()], [rnlogit()],
#' [rhev()], [rmnp()], [rrologit()] and [simulate.choicemodel()], and broom / texreg /
#' modelsummary support.
#'
#' @section Tests of the independence of irrelevant alternatives and of specification:
#' [hausman_mcfadden()], [small_hsiao()], [fp_tests()] (Fok and Paap 2025), the
#' battery [iia_tests()], the nested-logit score test [nested_lm_test()]
#' (McFadden 1987), [combine_test()] (Cramer and Ridder 1991), [term_test()],
#' the goodness-of-fit test [hl_test()] (Fagerland, Hosmer and Bofin 2008), the
#' stability test of rank-ordered fits [rank_test()], the score and
#' likelihood-ratio tests for
#' captive choosers [captivity_test()], and the parametric-bootstrap
#' calibration [boot_pvalue()].
#'
#' @author Benjamin E. Bagozzi \email{bagozzib@udel.edu}
#'
#' @section Data:
#' [pa2006] and [pres2004], the two applications of Bagozzi and Marchetti
#' (2017), whose published fits are regression tests of the package, and seven
#' data sets from published replication archives: [sanctions], [fluvax],
#' [bes1987], [defection], [bes2010], [cabinets], and [afghan_leaders] (see the
#' vignette "Seven published analyses").
#'
#' @useDynLib polychoice, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom stats model.frame model.matrix model.response terms delete.response reformulate
#' @importFrom stats glm.fit binomial optim na.omit na.pass complete.cases update as.formula
#' @importFrom stats nobs predict confint fitted residuals coef logLik AIC BIC vcov simulate
#' @importFrom stats plogis qlogis dlogis pnorm qnorm dnorm pchisq qchisq printCoefmat setNames
#' @importFrom stats quantile sd var weighted.mean rbinom runif rmultinom ave median
#' @importFrom Formula Formula
"_PACKAGE"
