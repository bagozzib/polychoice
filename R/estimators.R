## ---------------------------------------------------------------------------
## The exported estimators: mnl, imnl, nlogit, hev, mnprobit, remnl. Each checks its
## arguments and hands formula, data, and family to the engine (.pc_fit).
## ---------------------------------------------------------------------------

#' Multinomial and conditional logit regression
#'
#' Fits the multinomial logit (case-specific covariates, alternative-specific
#' coefficients), the conditional logit (alternative-specific covariates, generic
#' coefficients), or their combination by maximum likelihood, with
#' alternative-specific constants, varying choice sets, frequency weights,
#' offsets, unit fixed effects, and analytic, robust, cluster-robust, or
#' bootstrap standard errors. This is the baseline for [imnl()] and the model
#' the tests of IIA ([hausman_mcfadden()], [small_hsiao()], [fp_tests()]) take.
#'
#' @details
#' With case \eqn{i}, alternatives \eqn{j = 1, \dots, J} (the levels of the
#' response), base alternative \eqn{b}, and availability \eqn{a_{ij}},
#' \deqn{V_{ij} = o_{ij} + x_i'\beta_j + w_{ij}'\gamma, \quad \beta_b = 0, \qquad
#'       P(y_i = j) = \frac{a_{ij} e^{V_{ij}}}{\sum_k a_{ik} e^{V_{ik}}}.}
#' \eqn{x_i} are the case-specific covariates, \eqn{w_{ij}} the
#' alternative-specific covariates, \eqn{o_{ij}} an offset: McFadden's (1974)
#' conditional logit, with the multinomial logit as the case
#' \eqn{w = 0} (Train 2009, chapter 3). The
#' alternative-specific constants are the coefficients of an intercept in
#' \eqn{x_i}: the formula's intercept in the wide layout, the first part's in
#' the long layout (section "Two data layouts"). The coefficients are reported as
#' `"<alternative>:<term>"` for \eqn{\beta_j} and `"<term>"` for \eqn{\gamma}.
#'
#' @section Two data layouts:
#' **Wide** (one row per case): `data` is a plain data frame, the response is
#' the chosen alternative (a factor, or integer codes; a character vector is
#' refused because its alphabetical order is rarely the intended one), and the
#' formula is `y ~ x1 + x2`. All alternatives are available to every case.
#'
#' **Long** (one row per case-alternative pair): give `id` (the case column)
#' and `alt` (the alternative column; its levels define the alternatives). The
#' response is the 0/1 (or logical) choice indicator, exactly one per case. The
#' formula has up to two parts separated by `|`: `choice ~ w1 + w2 | x1 + x2`,
#' alternative-specific covariates first (generic coefficients; the intercept
#' of this part switches the alternative-specific constants on, `0 +` switches
#' them off), case-specific covariates second (they must be constant within a
#' case). Alternatives without a row for a case are unavailable to it. An
#' alternative-specific coefficient on an alternative-specific variable is the
#' interaction with the alternative factor, `w:alt`.
#'
#' The constants follow the first part only:
#' * `choice ~ w`, `choice ~ w | 1` and `choice ~ w | x` have them;
#' * `choice ~ 0 + w`, `choice ~ w - 1` and `choice ~ 0 + w | x` do not.
#'
#' The second part's intercept has no effect. `| 1` is an empty case-specific
#' part (constants only, when the first part has them), and `| 0` does not
#' drop the constants, unlike in \pkg{mlogit} (a message says so). The
#' wide-layout formula `y ~ x1 + x2` is the second part alone, its intercept
#' giving the constants (`y ~ 0 + x1` drops them).
#'
#' @section Missing values, weights, standard errors:
#' Rows with a missing value in any model variable (or in `weights`, `offset`,
#' `cluster`, `fe`) are dropped before fitting; in the long layout a missing
#' value on any row of a case drops the whole case (a message gives the count),
#' and cases left with a single available alternative are dropped too, as
#' are cases of weight zero. `object$data` holds the retained
#' rows. `weights` are frequency weights throughout: a fit with weight 2 equals
#' a fit on duplicated cases in the estimates, the model-based standard errors,
#' the robust and cluster sandwiches (meat \eqn{\sum_i w_i s_i s_i'}), the
#' information criteria (the weight total is \eqn{n}), the Vuong and
#' likelihood-ratio statistics, the IIA tests (except the Small-Hsiao test,
#' whose half-samples balance cases rather than weight), the averaged quantities of
#' interest, and the bootstraps (`se = "bootstrap"` draws multinomial replicate
#' weights with the weight total as size, over units for [remnl()]; with a
#' `cluster` the clusters are drawn and carry their rows' weights;
#' [boot_pvalue()] and [boot_lr_test()] expand the data to one row per case and
#' need integer weights). Non-integer weights are accepted as frequency weights
#' (a message says so). Survey
#' (probability) weights are NOT supported: the point estimates would be the
#' pseudo-maximum-likelihood ones, but the sandwich needs \eqn{\sum_i w_i^2 s_i
#' s_i'} and the information criteria, tests, and bootstraps have no
#' design-based meaning (the cluster sandwich with the primary sampling units
#' as clusters, \eqn{\sum_g (\sum_{i \in g} w_i s_i)(\sum_{i \in g} w_i s_i)'},
#' has the with-replacement linearization form, but strata and
#' finite-population corrections are not handled). `se = "robust"` is the sandwich over cases,
#' `se = "cluster"` over clusters (no small-sample degrees-of-freedom
#' adjustment, as in the `iop` and `underdisp` packages; Stata's
#' `vce(cluster)` for maximum-likelihood commands multiplies by
#' \eqn{G/(G-1)}). `nobs()` returns the number of cases; `logLik()` carries the
#' weight total as its `nobs` attribute, which is what `BIC()` uses.
#'
#' @section Unit fixed effects:
#' `fe = "<unit column>"` (the name of the column identifying the units, such as
#' `fe = "country"`) adds an alternative-specific constant per unit: one dummy per
#' unit in every non-base equation. A unit that never chooses alternative
#' \eqn{j} has that constant at \eqn{-\infty}; the limiting likelihood treats
#' \eqn{j} as unavailable to that unit, which is what the fit does (with a
#' message), and cases left with a single alternative drop out (a captive
#' choice says nothing about preferences and does not count as a choice; as
#' dropping cases can change what a unit chose, the rule and the drop are
#' repeated until nothing changes). A unit that chose a single alternative
#' therefore drops out entirely. A unit that
#' never chooses the base likewise has the base unavailable; its constants are
#' then identified only up to a common shift, so one of its dummies is
#' omitted. The common constant of alternative \eqn{j} refers to a reference
#' unit: the first unit (in level order) that chose both \eqn{j} and the base,
#' or else one that links \eqn{j} to an alternative already anchored; an
#' alternative whose choosers are not connected to the base in this way has no
#' identified constant, and the fit stops and names it. Without
#' alternative-specific constants (a formula with `- 1`) there is nothing to
#' anchor, and every unit keeps all its dummies (one fewer for a unit that never
#' chose the base). The dummies are dense columns, so a fit with thousands of
#' units and many cases needs memory for an \eqn{nJ \times} (units \eqn{\times}
#' alternatives) matrix. Maximum-likelihood unit
#' constants carry the incidental-parameters bias when units have few cases
#' (the fit warns below a median of 10). For short panels the alternatives are
#' the jackknife correction below and the correlated-random-effects model:
#' [remnl()] on the formula and data that [mundlak()] returns (the pooled
#' [mnl()] fit of the same formula gives population-averaged slopes, attenuated
#' towards zero).
#'
#' `bias_correct = "jackknife"` removes the leading \eqn{1/T} term of that bias
#' from the slopes by the split-panel jackknife (Dhaene and Jochmans 2015,
#' Review of Economic Studies): the model is refit on the first and second
#' halves of every unit's cases (rows are split in the order supplied, which
#' should be temporal order; a unit with an odd number of cases contributes
#' both possible splits, averaged) and the corrected slope is `2 * full -
#' mean(halves)`. Alternative-specific constants and the unit dummies are left
#' at their maximum-likelihood values, as are the log-likelihood, `fitted()`,
#' `residuals()`, and the standard errors; `coef()`, `predict()`, `rrr()`,
#' `confint()`, and `first_difference()` use the corrected slopes together with
#' the ML constants. The ML slopes are kept in `$uncorrected` and the
#' half-panel estimates in `$jackknife`. The identity requires the two halves
#' to estimate the same parameter (time homogeneity), so the correction is
#' REFUSED, with a warning naming the check and the ML fit returned, when the
#' half-panel slopes disagree beyond 2.5 times the mean disagreement of five
#' placebo splits, random halves of every unit's cases that ignore the time
#' order (a regime change), or when the full fit's residuals trend with the
#' within-unit time order for some alternative (an omitted trend, at the 1
#' percent level with a Bonferroni correction over alternatives: put a time
#' control in the formula, a trend rather than period dummies, which are not
#' identified within the half panels). The placebo splits have the time
#' split's structure (all units, halves of the same lengths), so under time
#' homogeneity the ratio of the two disagreements is centred near 1 at every
#' panel length, while a mid-panel regime change separates the time halves and
#' is balanced in every random split; the splits come from fixed seeds and do
#' not disturb the caller's random-number stream. The constant 2.5 is the
#' smallest of 1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, and 4 that kept the refusal
#' rate of time-homogeneous panels at or below 5 percent at T = 6, 10, and 20
#' in a calibration Monte Carlo (150 units, J = 3), under a rule written down
#' before that run, and the choice is checked on fresh seeds
#' (`verification/jackknife_gate_mc.R` and `jackknife_gate_holdout.R` in the
#' repository, results in the verification ledger). In those simulations a
#' regime change was refused in 92 percent of the panels at T = 6 and in all
#' of them at T = 10 and 20, and an omitted trend in all panels at T = 10 and
#' 20 but only about half at T = 6. A refusal leaves `$bias_correct` at
#' `"none"` and the reason in `$jackknife_refused`, which `summary()` prints,
#' so that a refused correction cannot pass for a corrected one. It warns when
#' the median panel length is below 6 or the implied correction is large. The
#' correction removes the leading \eqn{1/T} term only: in the same simulations
#' it over-corrects at T = 6 and 10 (the corrected slopes keep a fraction of
#' the maximum-likelihood bias, with the opposite sign), yet it lowers the
#' root mean squared error and brings the coverage of Wald intervals (with the
#' maximum-likelihood standard errors) closer to nominal at every T studied.
#'
#' @section Fit-quality flags:
#' Every fit carries `converged` (`FALSE` when the optimizer stopped short or
#' at a saddle, where the log-likelihood rises along a direction in which the
#' Hessian curves upward; a near-zero curvature along which it does not rise,
#' as on the flat ridge of an unidentified model, is not a saddle. At a saddle
#' the fit restarts from the point where the log-likelihood rose most along
#' that direction, up to three times, with a message when that reaches a
#' higher maximum; the saddle warning stays when it does not),
#' `diverging` (the names of the diverging coefficients), `phi_boundary` (an
#' inflated fit's split equation at one regime) and `saddle`,
#' `boundary` (a structural parameter at zero or infinity, a split equation
#' that puts every case in one regime, or diverging coefficients: separation,
#' or a unit constant at infinity; in [mnl()] and [rologit()] fits without
#' unit fixed effects separation is decided exactly by the linear program of
#' [separation()], elsewhere a coefficient counts as diverging when its
#' standardized size exceeds 15 and its Wald |z| is below 1.96), `ill_conditioned` (a singular or nearly
#' singular information matrix; `rcond` holds its reciprocal condition
#' number, taken on the matrix scaled to a unit diagonal, and `singular` is
#' `TRUE` when that number is below 1e-13 or the matrix could not be inverted, and
#' a generalized inverse took its place), and `n_modes`, the number of distinct local maxima the search
#' found. Every end of the search (one for a single start; for the multistart
#' fits every polished start) is classified by the rules of these flags: a
#' point with a diverging coefficient, a structural parameter at zero or
#' infinity, or a split equation at one regime is a boundary, not a maximum,
#' and so are a saddle and an end that did not converge. Two ends are the
#' same point when their log-likelihoods differ by less than 0.01 and either
#' their estimates differ by less than a tenth of a standard error or every
#' fitted probability by less than 1e-4: ends at other estimates with the
#' same probabilities lie on a ridge where the data do not identify the
#' estimates, and `detail` says how many ends joined the point that way. Maxima
#' whose information matrix is ill-conditioned are one point when their
#' log-likelihoods differ by less than 0.01, whatever their estimates (a nearly
#' flat direction lets the search stop anywhere along it), and `detail` counts
#' those ends too; two different maxima of that kind at one height would be
#' counted once, and only the highest end's estimates are kept. `search` lists the
#' distinct points (`loglik`, `type`: `"maximum"`, `"boundary"`, `"saddle"`
#' or `"not converged"`, `reached`: the number of ends there, `reported`,
#' `detail`), and the rows of `search_theta` hold their estimates on the
#' internal scale of `start`; `start_logliks` and `screen_logliks` hold the
#' log-likelihoods of every polished and every screened start. `n_modes` is 0
#' when no end is a local maximum (the reported fit then carries one of the
#' flags above). Each raises a warning or message at fit time;
#' [summary.choicemodel()] prints them and [glance.choicemodel()] returns
#' them. What they mean for the results: after `converged = FALSE` neither the
#' estimates nor the standard errors should be used (refit from other
#' `start` values, or simplify the model); with `boundary = TRUE` the fit is
#' at a limit (it keeps `converged = TRUE` when the optimizer met its criterion
#' there): the estimates of the parameters at the boundary and their standard
#' errors are not meaningful, the others usually are, and the restricted model
#' is the one to compare with (the search holds log scales, log inclusive
#' values, the multinomial probit's log Cholesky diagonals and log standard deviations
#' within +/- 8; where the log-likelihood still rises at that bound the
#' warning says so, and the reported log-likelihood is then below the model's
#' supremum, so a comparison with the restricted model understates it);
#' `ill_conditioned = TRUE` (a reciprocal condition
#' number below 1e-9, or a generalized inverse in the covariance) means some
#' combination of parameters is barely identified and its standard errors are
#' unreliable, and with `singular = TRUE` some parameters are not identified
#' at all: their estimates and standard errors (from the generalized inverse)
#' are not meaningful, and the summary's footer says so, since their rows in
#' the table look like the others. The information of a model whose
#' parameters are not identified is singular at the maximum itself, and close
#' to it its condition number is of the order of the gradient. A fit whose
#' information is ill-conditioned is therefore brought to a stationary point
#' by a few Newton steps before the flag is set, so that the flag does not
#' depend on how closely the search converged;
#' `n_modes > 1` means the reported fit is the highest of several
#' maxima found. The reported fit is the highest end of the search, maximum or
#' not: where it is at a boundary and a lower end is a local maximum, a message
#' says so, and that maximum's row of `search_theta` refits it as `start`. An
#' alternative that is never chosen (where another is
#' available) is separation: in [mnl()] and [rologit()] fits without unit fixed
#' effects the fit runs and reports it (`firth = TRUE` or a prior then gives
#' finite estimates), while [imnl()], [hev()], [mnprobit()],
#' [remnl()] and fixed-effects fits refuse such data, and [nlogit()] does when
#' the alternative has its own constant or slopes (with generic coefficients
#' only it refuses a nest of which nothing is chosen). Data that overlap can
#' still give very large, imprecise
#' estimates (near separation), which raise no flag: their Wald statistics
#' understate the evidence against zero (Hauck and Donner 1977), so the
#' likelihood-ratio tests of [term_test()] (`type = "lr"`) and [lr_test()] are
#' the ones to use there, or `firth = TRUE`. A fixed-effects fit also carries `bias_correct` and, when a requested
#' correction was refused, the reason in `jackknife_refused`.
#'
#' @section Penalized (Firth) estimation:
#' Under separation (a covariate, or a combination of covariates, that
#' predicts some choices perfectly; Cook, Niehaus and Zuhlke 2018 for the
#' multinomial logit) the maximum-likelihood estimates do not exist, and the
#' reported ones are where the optimizer stopped. `firth = TRUE` maximizes the
#' log-likelihood penalized by half the log-determinant of the Fisher
#' information (Jeffreys' prior),
#' \eqn{\ell^*(\theta) = \ell(\theta) + \tfrac12 \log |I(\theta)|}, which for the
#' multinomial and conditional logit is Firth's (1993) bias-reducing estimator
#' (Bull, Mak and Greenwood 2002; Kosmidis and Firth 2011): its estimates are
#' finite under separation and have smaller small-sample bias. The penalty is
#' invariant to the choice of base alternative. It shrinks the fitted
#' probabilities towards equal probabilities, so their averages no longer
#' equal the observed shares (Puhr et al. 2017).
#'
#' A penalized fit reports `loglik`, the ordinary log-likelihood at the
#' penalized estimates, and `loglik_penalized`, the maximized objective; the
#' information criteria and pseudo-R-squared measures built on `loglik` are
#' not those of maximum likelihood. The standard errors are from the inverse
#' information at the penalized estimates (the robust and cluster sandwiches
#' use the per-case scores plus their share of the penalty's gradient; the
#' bootstrap refits with the penalty). The exact separation check of
#' [separation()] runs on the data: the names of the coefficients along the
#' separating direction are in `separation`, printed and announced with a
#' message, since their penalized estimates are then determined largely by the
#' penalty. The penalized likelihood is not always concave and can have more
#' than one local maximum (with a high-leverage covariate, even in a binary
#' logit), so without `start` the fit runs from zero and from the
#' maximum-likelihood estimate (halved towards zero where the penalty is not
#' finite) and keeps the higher maximum, with a warning when the two runs end
#' at different maxima (`start_logliks` holds both values, `n_modes` counts
#' the maxima, and a penalized fit keeps no `search`). Recommended inference: the profile penalized-likelihood
#' intervals of [confint()] (its default for a penalized fit) and the
#' penalized likelihood-ratio tests of [term_test()], [combine_test()], and
#' [lr_test()] (two penalized fits), which maximize the full model's
#' penalized likelihood under the restriction (Heinze and Schemper 2001,
#' 2002). The tests derived for maximum-likelihood estimates refuse a
#' penalized fit: [hausman_mcfadden()], [small_hsiao()], [fp_tests()],
#' [iia_tests()], [nested_lm_test()], [captivity_test()], [vuong()], and
#' [boot_lr_test()]. `firth = TRUE` is not available with unit fixed effects,
#' and the other estimators do not offer it (in [imnl()] and [remnl()] the
#' penalty tends to minus infinity where the model reduces to the
#' multinomial logit, so it would favour the larger model).
#'
#' Under separation the penalty decides how large the separated estimates are
#' (Rainey 2016). Jeffreys' prior is one default, justified by bias reduction;
#' `prior = coef_prior(...)` replaces it with independent normal, Cauchy or t
#' priors centred at zero on coefficients scaled to their columns (Gelman et
#' al.'s 2008 default is `coef_prior()`), and fits under several priors show
#' how much of a separated estimate comes from the prior. A fit with a prior is
#' a penalized fit in every respect above (two starts, `separation`, profile
#' intervals, penalized likelihood ratios, the refusals); its covariance is the
#' inverse curvature of the log posterior at the mode. For maximum-likelihood
#' fits of separated data, [term_test()], [combine_test()], [lr_test()] and
#' [boot_lr_test()] compare the suprema of the log-likelihoods (likelihood-ratio
#' tests remain valid under separation where Wald tests do not; Rainey 2024).
#'
#' @param formula A formula; see the section "Two data layouts" of [mnl()].
#' @param data A data frame.
#' @param ref The base (reference) alternative: a level of the response.
#'   Default: the first level.
#' @param id,alt For the long layout, the names of the case and alternative
#'   columns of `data`. Leave both `NULL` for the wide layout.
#' @param weights Optional frequency weights: a column name in `data` or a
#'   numeric vector with one entry per row (constant within a case in the long
#'   layout).
#' @param offset Optional offset on the utility scale: in the long layout a
#'   column name or a vector with one entry per row; in the wide layout an
#'   \eqn{n \times J} matrix (columns in level order).
#' @param fe Optional column name in `data` identifying units that receive
#'   alternative-specific fixed effects.
#' @param bias_correct With `fe`: `"none"` (default) or `"jackknife"`, the
#'   split-panel jackknife of Dhaene and Jochmans (2015) for the slopes; see
#'   the fixed-effects section.
#' @param se Standard errors: `"analytic"` (inverse observed information,
#'   default), `"opg"` (inverse outer product of the scores, the BHHH form that
#'   `mlogit` reports for the nested logit and HEV), `"robust"`, `"cluster"`
#'   (needs `cluster`), `"bootstrap"`
#'   (nonparametric over cases, or over clusters when `cluster` is given, with
#'   `nboot` refits from the full-sample estimate; the replicates are kept in
#'   `object$boot` and [confint()] can report percentile intervals), or
#'   `"none"`. With `fe`, the bootstrap resamples units (a unit drawn twice
#'   counts as two units, each with its own fixed effects; a case bootstrap
#'   would give units random numbers of cases and a second incidental-parameters
#'   bias) and covers the slopes only: the alternative-specific constants and
#'   the unit dummies refer to reference units that change between replicates,
#'   so their bootstrap standard errors are `NA`, and [first_difference()]
#'   refuses intervals that need them (use analytic, robust, or cluster
#'   standard errors for those). A replicate whose fixed effects cannot be
#'   built (an alternative left unconnected to the base) counts as failed, the
#'   standard errors are then conditional on resamples in which they can be
#'   built, and the fit stops when fewer than half the replicates are usable. A
#'   `cluster` given with `fe` must nest the units (whole units are drawn; one
#'   that cuts across units is refused). The replicates are not
#'   jackknife-corrected.
#' @param cluster Optional cluster identifier (a column name in `data` or a
#'   vector aligned to its rows, constant within a case); supplying it selects
#'   `se = "cluster"` unless `se` is given explicitly.
#' @param nboot Number of bootstrap refits for `se = "bootstrap"` (default 200).
#' @param cores Number of parallel workers for the bootstrap (a PSOCK cluster;
#'   default 1). `set.seed()` before the call makes the resamples reproducible;
#'   serial and parallel runs use different random streams.
#' @param start Optional starting values: the full parameter vector in the order
#'   of `coef()` but on the INTERNAL scale, which differs from `coef()` for the
#'   structural parameters (the logarithms of the inclusive values, scales, and
#'   standard deviations; the free Cholesky elements with a log diagonal for
#'   the multinomial probit covariance). `object$theta` of an earlier fit is on that scale.
#'   Rarely needed for `mnl()`, whose log-likelihood is globally concave; the
#'   other estimators start from the `mnl()` estimates (and `imnl()` from a
#'   multistart grid), so `start` serves to probe another basin of a
#'   non-concave likelihood.
#' @param maxit,reltol Controls for the BFGS stage passed to [stats::optim()]
#'   (the plain `mnl()` is solved by damped Newton steps first, with BFGS as the
#'   fallback); the Newton polish that follows runs up to 25 further iterations
#'   with the analytic Hessian of the logit family or central differences of the
#'   analytic gradient otherwise.
#' @param firth If `TRUE`, Firth's penalized likelihood (see the section
#'   "Penalized (Firth) estimation"); default `FALSE`, maximum likelihood.
#'   Weights are read as frequency weights and must total more than 2.
#' @param prior Optional: priors for the coefficients, made by [coef_prior()]
#'   (normal, Cauchy or t, centred at zero, on coefficients scaled to their
#'   columns); the fit is the posterior mode, a penalized fit like
#'   `firth = TRUE` (profile intervals, penalized likelihood-ratio tests, the
#'   separated coefficients in `separation`). Not with `firth = TRUE` or unit
#'   fixed effects.
#' @return An object of class `c("mnl", "choicemodel")`: a list with
#'   `coefficients`, `vcov`, `loglik`, `fitted.values` (the \eqn{n \times J}
#'   matrix of choice probabilities, rows named by case), `converged`,
#'   `boundary`, `ill_conditioned`, `singular`, `n_modes`, `search`, `start_logliks`, the
#'   estimation `data`, and the design used by the methods. See [summary.choicemodel()],
#'   [predict.choicemodel()], [first_difference()], [rrr()], [fit_stats()].
#' @references
#' Bull, S.B., Mak, C. and Greenwood, C.M.T. (2002). A modified score function
#' estimator for multinomial logistic regression in small samples.
#' \emph{Computational Statistics & Data Analysis}, 39, 57-74.
#' \doi{10.1016/S0167-9473(01)00048-2}
#'
#' Cook, S.J., Niehaus, J. and Zuhlke, S. (2018). A warning on separation in
#' multinomial logistic models. \emph{Research & Politics}, 5(2).
#' \doi{10.1177/2053168018769510}
#'
#' Dhaene, G. and Jochmans, K. (2015). Split-panel jackknife estimation of
#' fixed-effect models. \emph{Review of Economic Studies}, 82, 991-1030.
#' \doi{10.1093/restud/rdv007}
#'
#' Firth, D. (1993). Bias reduction of maximum likelihood estimates.
#' \emph{Biometrika}, 80, 27-38. \doi{10.1093/biomet/80.1.27}
#'
#' Gelman, A., Jakulin, A., Pittau, M.G. and Su, Y.-S. (2008). A weakly
#' informative default prior distribution for logistic and other regression
#' models. \emph{Annals of Applied Statistics}, 2, 1360-1383.
#' \doi{10.1214/08-AOAS191}
#'
#' Hauck, W.W. and Donner, A. (1977). Wald's test as applied to hypotheses in
#' logit analysis. \emph{Journal of the American Statistical Association},
#' 72, 851-853. \doi{10.1080/01621459.1977.10479969}
#'
#' Heinze, G. and Schemper, M. (2001). A solution to the problem of monotone
#' likelihood in Cox regression. \emph{Biometrics}, 57, 114-119.
#' \doi{10.1111/j.0006-341X.2001.00114.x}
#'
#' Heinze, G. and Schemper, M. (2002). A solution to the problem of separation
#' in logistic regression. \emph{Statistics in Medicine}, 21, 2409-2419.
#' \doi{10.1002/sim.1047}
#'
#' Kosmidis, I. and Firth, D. (2011). Multinomial logit bias reduction via the
#' Poisson log-linear model. \emph{Biometrika}, 98, 755-759.
#' \doi{10.1093/biomet/asr026}
#'
#' McFadden, D. (1974). Conditional logit analysis of qualitative choice
#' behavior. In P. Zarembka (ed.), \emph{Frontiers in Econometrics}, 105-142.
#' New York: Academic Press.
#'
#' Puhr, R., Heinze, G., Nold, M., Lusa, L. and Geroldinger, A. (2017).
#' Firth's logistic regression with rare events: accurate effect estimates
#' AND predictions? \emph{Statistics in Medicine}, 36, 2302-2317.
#' \doi{10.1002/sim.7273}
#'
#' Rainey, C. (2016). Dealing with separation in logistic regression models.
#' \emph{Political Analysis}, 24, 339-355. \doi{10.1093/pan/mpw014}
#'
#' Rainey, C. (2024). Hypothesis tests under separation. \emph{Political
#' Analysis}, 32, 172-185. \doi{10.1017/pan.2023.28}
#'
#' Train, K.E. (2009). \emph{Discrete Choice Methods with Simulation}, 2nd ed.
#' Cambridge University Press.
#' @examples
#' set.seed(1)
#' d <- rmnl(500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)))
#' m <- mnl(y ~ x1 + x2, data = d)
#' summary(m)
#' head(predict(m))
#' rrr(m)
#' ## Firth's penalized likelihood, with profile intervals
#' mf <- mnl(y ~ x1 + x2, data = d, firth = TRUE)
#' confint(mf)
#' @seealso [imnl()], [hausman_mcfadden()], [small_hsiao()], [fp_tests()],
#'   [iia_tests()], [vuong()], [mundlak()]
#' @family estimators
#' @export
mnl <- function(formula, data, id = NULL, alt = NULL, ref = NULL, weights = NULL, offset = NULL, fe = NULL,
                bias_correct = c("none", "jackknife"),
                se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"), cluster = NULL, nboot = 200L,
                cores = 1L, start = NULL, maxit = 1000L, reltol = 1e-10, firth = FALSE, prior = NULL) {
  se_missing <- missing(se); se <- match.arg(se); bias_correct <- match.arg(bias_correct)
  .pc_fit(formula, data, inflated = FALSE, ref = ref, id = id, alt = alt, weights = weights, offset = offset,
          cluster = cluster, fe = fe, se = se, nboot = nboot, cores = cores, start = start, maxit = maxit,
          reltol = reltol, call = match.call(), se_missing = se_missing, bias_correct = bias_correct, firth = firth, prior = prior)
}

#' Inflated multinomial logit regression
#'
#' Fits the baseline-inflated multinomial logit of Bagozzi and Marchetti (2017)
#' and its generalization to any inflated alternative: a mixture of a
#' degenerate distribution at the inflated alternative \eqn{k} and the
#' multinomial / conditional logit of [mnl()], with a binary split equation
#' deciding which regime a case is in. The canonical application is vote
#' choice with abstention as the base: routine nonvoters (the captive regime, always
#' abstain) are mixed with occasional nonvoters and voters (the choice regime, the
#' multinomial process), and the split equation separates them.
#'
#' @details
#' With \eqn{s_i \in \{0, 1\}}, \eqn{P(s_i = 1 \mid z_i) = \phi_i = F(z_i'\delta)}
#' (\eqn{F} logistic by default, or standard normal), and the [mnl()]
#' probabilities \eqn{P_{ij}},
#' \deqn{P(y_i = k) = (1 - \phi_i) + \phi_i P_{ik}, \qquad P(y_i = j) = \phi_i P_{ij}, \; j \ne k.}
#' In the discrete-choice literature, a model with a logit split equation is
#' the parametrized logit captivity model of Swait and Ben-Akiva (1987, eq. 3)
#' with captivity to \eqn{k} only. That model generalizes the dogit of Gaudry
#' and Dagenais (1979), which adds to the multinomial logit a share of decision
#' makers captive to each alternative, by making the captivity odds
#' \eqn{\theta_{ik} = (1 - \phi_i)/\phi_i} a function of covariates; here
#' \eqn{\theta_{ik} = e^{-z_i'\delta}}, so the split coefficients are the
#' negatives of Swait and Ben-Akiva's captivity coefficients.
#'
#' Several inflated alternatives (`inflate = c("abstain", "A")`) give the
#' parametrized logit captivity model with captive shares for each of them
#' (Swait and Ben-Akiva 1987), which suits, say, routine nonvoters together
#' with partisans captive to their party. The regimes
#' (choosing, or captive to one of the inflated alternatives \eqn{c \in C})
#' follow a multinomial logit in \eqn{z}: with
#' \eqn{\eta_{ic} = z_i'\delta_c},
#' \deqn{\pi_{i0} = \frac{1}{1 + \sum_{c} e^{-\eta_{ic}}}, \qquad \pi_{ic} = \frac{e^{-\eta_{ic}}}{1 + \sum_{c} e^{-\eta_{ic}}},
#'   \qquad P(y_i = j) = \pi_{i0} P_{ij} + 1\{j \in C\}\, \pi_{ij},}
#' so \eqn{\eta_{ic}} is the log-odds of choosing against captivity to
#' \eqn{c}: a positive coefficient makes captivity to \eqn{c} less likely, the
#' reading of the single split equation, which is the case \eqn{|C| = 1}.
#' Every inflated alternative gets its own coefficients on the same \eqn{z}
#' (named `infl_<alternative>:<term>`; Swait and Ben-Akiva allow each
#' alternative its own captivity covariates, while here every inflated
#' alternative's equation has the same \eqn{z}); an inflated alternative outside a
#' case's choice set is not one of its regimes. Several inflated alternatives
#' take the logit split only. `predict(type = "regime")` then returns the
#' matrix of the regime probabilities and `type = "posterior"` their
#' posteriors given the observed choice; `phi` is the probability of the
#' choice regime. A captive share that collapses to zero for every case is
#' flagged `boundary`; the model with fewer inflated alternatives is nested
#' at that boundary, so [boot_lr_test()] compares the two (as it compares
#' the inflated model with [mnl()]). Each captive share is identified mainly
#' through the mixture's form, as the single one is; an exclusion
#' restriction helps. In the package's simulation (four
#' alternatives, two of them captive, 200 samples per design), with a split
#' covariate outside the outcome equations the 95 percent intervals covered
#' 0.92-0.975 at n = 5,000; at n = 2,000 two cells missed, the constant of
#' `c`'s split equation (coverage 0.915) and its slope (a bias of -3.1 Monte
#' Carlo standard errors). With the split on an outcome covariate instead,
#' 3 of 200 fits at n = 5,000 had diverging choice coefficients (flagged
#' `boundary`) and the intervals of that covariate's outcome coefficients
#' covered 0.89-0.90; at n = 2,000, 17 of 200 fits were at a boundary and
#' those intervals covered 0.80-0.81.
#'
#' The MNL is the boundary case \eqn{\phi \to 1}, where the split
#' coefficients are not identified, so neither the chi-square reference of a
#' likelihood-ratio test nor the normal reference of [vuong()] holds for the
#' comparison with [mnl()]: [boot_lr_test()] calibrates both statistics by a
#' parametric bootstrap under the MNL (the information criteria are a
#' descriptive complement). The exception is a single inflated alternative
#' with an intercept-only split (`y ~ x | 1`, a constant captive share), whose
#' comparison with the MNL has a regular boundary: [lr_test()] and
#' `captivity_test(type = "lr")` refer it to the 50:50 chi-square mixture. The split-equation coefficients are reported
#' with the prefix `infl_`. Identification rests on the mixture structure, as
#' for the zero-inflated count and inflated ordered models: an exclusion
#' restriction (a variable in `z` not in `x`, or vice versa) is advisable but
#' not required (Bagozzi and Marchetti 2017). The split equation models
#' \eqn{P(s_i = 1)}, the probability of being in the CHOICE regime: a positive
#' coefficient means the covariate makes a case less likely to be a routine
#' chooser of the inflated alternative. `predict(type = "regime")` returns
#' \eqn{\phi_i} (and 1 for a case whose choice set lacks \eqn{k}: the captive regime is
#' impossible there, so such cases count as in the choice regime; the boundary
#' check and `summary()`'s mean of \eqn{\phi_i} leave them out), `type =
#' "posterior"` the regime probability given the observed
#' choice (\eqn{\phi_i P_{ik} / P^*_{ik}} for a case that chose \eqn{k}, 1
#' otherwise), which classifies the observed choosers of \eqn{k} into routine
#' and occasional ones.
#'
#' The formula has the split equation as its **last** part: in the wide
#' layout `y ~ x1 + x2 | z1 + z2`, in the long layout
#' `choice ~ w1 | x1 + x2 | z1 + z2` (write `| 1 |` for an empty case-specific
#' part). Everything else -- layouts, weights, offsets, standard errors -- is
#' as in [mnl()]. Unit fixed effects are not available: with a constant per
#' unit in the outcome equations the level of \eqn{\phi} is identified only
#' through the functional form of the mixture, units that always chose the
#' inflated alternative (the routine choosers the model is designed for) drop
#' out, units with few such choices have them all put in the captive regime and their
#' constants run to infinity, and units that never chose it pull \eqn{\phi}
#' towards 1. In the package's simulations the split coefficients came out
#' 25 to 42 percent too large at 12 to 20 cases per unit (and two to three
#' times too large at 4), with a third or more of the unit constants at
#' infinity, so `imnl()` refuses `fe`; for panels, fit without it with
#' `cluster` set to the unit, or use [mnl()] with `fe` for the outcome
#' equations.
#'
#' The inflated likelihood can have several local maxima (it has two on each
#' of `pa2006` and `pres2004`, and many where the split equation can separate a
#' few cases into the inflated regime). Estimation therefore runs a grid of
#' fourteen starts (the [mnl()] coefficients, or an MNL that down-weights the
#' inflated-alternative cases, combined with a binary fit of
#' \eqn{1\{y \ne k\}} on \eqn{z} at full, halved, doubled, tripled, or zero
#' slopes and several shifts of the split intercept), screens them by
#' quasi-Newton, and polishes every screened start to convergence with the
#' analytic Hessian. The incumbent is then challenged with its split slopes
#' halved, enlarged by half, and doubled, again from each new incumbent (at
#' most three rounds). With one inflated alternative, `crosscheck = TRUE`
#' adds before the challengers the same search under the other split link,
#' whose estimates, the split coefficients rescaled by 1.6 (the ratio of the
#' logistic and normal densities at zero), are polished under this link as one
#' more candidate: the two links' surfaces differ, and on the vote-choice data
#' of Glasgow (2001) the logit split's own search ends 4.9 log-likelihood units below
#' the maximum it reaches from the probit split's estimates. It about doubles the
#' fitting time, so it is not the default; use it when the message reports
#' several maxima. The highest maximum found is reported; it need not be the
#' global one, so probe other basins with `start` too. The challengers start
#' from the best polished fit only, so where the split coefficients diverge
#' (the likelihood then has many flat plateaus, each a different hard
#' classification of a few cases into the regimes) the plateau reported
#' depends on the search path; such a fit is flagged `boundary`, and its split
#' estimates are not maxima to interpret. The
#' log-likelihoods of the polished starts (converged or not), of the
#' cross-check and of the challengers that converged or improved, and the
#' maximum reached by a saddle step, are in
#' `object$start_logliks`, those of all screened starts in
#' `object$screen_logliks`, the distinct points where these searches ended,
#' each classified as a local maximum, a boundary (diverging coefficients or a
#' split equation at one regime), a saddle, or not converged, in
#' `object$search` (with their estimates in `object$search_theta`), the number
#' of distinct local maxima among them in `object$n_modes` (see the section
#' "Fit-quality flags" of [mnl()]), and with `crosscheck = TRUE` the
#' cross-check in `object$crosscheck` (the other link's maximum, the maximum
#' this link reached from it, and whether that improved on this link's own
#' search); a message
#' reports multimodality. A fit whose split equation puts every case in one
#' regime has collapsed onto the MNL; it is flagged `boundary` with a warning.
#'
#' The likelihood can also rise without bound along a direction of the
#' outcome coefficients, to a supremum at infinity: separation within the
#' mixture. The outcome stage can then predict some cases of the inflated
#' alternative with certainty while the split equation puts its other cases
#' in the captive regime, which overlapping data under the multinomial logit do not
#' prevent ([separation()] checks the multinomial logit and the rank-ordered logit only).
#' Strictly, the maximum-likelihood estimate then does not exist, and the
#' reported fit is the highest end the search reaches. On `pres2004` a search from
#' 200 starts finds two interior maxima, the published fit (-545.279) and
#' -550.259, and from two starts points where the Kerry and Bush coefficients
#' diverge together, the higher at -542.851: 21 abstainers are predicted
#' with certainty by the outcome stage and the other 84 are captive. The
#' multistart does not reach that region; a fit started there is flagged
#' (diverging coefficients, `boundary`, a singular information matrix). On
#' `pa2006` the same search finds the multistart's two maxima (-2069.971 and
#' -2070.990) and no higher point.
#'
#' @inheritParams mnl
#' @param fe Not available for the inflated model (see Details); an error if
#'   given.
#' @param inflate The inflated alternative: a level of the response. Default:
#'   `ref`, the base alternative (baseline inflation). Several levels give the
#'   model with a captive share for each (see Details).
#' @param link Link of the split equation: `"logit"` (default; the published
#'   model) or `"probit"`.
#' @param crosscheck With one inflated alternative, also search under the other
#'   split link and polish this link's fit from its rescaled estimates (see
#'   Details); about doubles the fitting time. `FALSE` (default) skips it. Ignored with
#'   several inflated alternatives (which take the logit link only) and when
#'   `start` is given. The parametric-bootstrap refits of [boot_lr_test()] and
#'   [boot_pvalue()] search as the fit did.
#' @param se Standard errors: `"analytic"` (inverse observed information,
#'   default), `"opg"`, `"robust"`, `"cluster"` (needs `cluster`),
#'   `"bootstrap"` (nonparametric over cases, or over clusters when `cluster`
#'   is given, with `nboot` refits from the full-sample estimate), or `"none"`.
#' @param start Optional starting values, which replace the multistart: the
#'   full parameter vector in the order of `coef()`, the outcome coefficients
#'   and then the split coefficients (`infl_`), on the same scale as `coef()`
#'   (`object$theta` of an earlier fit).
#' @return An object of class `c("imnl", "choicemodel")`; as [mnl()] plus
#'   `phi` (the fitted \eqn{P(s_i = 1)}), `prob_outcome` (the \eqn{n \times J}
#'   outcome-stage probabilities \eqn{P_{ij}}), `inflate`, and `link`; with
#'   several inflated alternatives also `regime`, the matrix of the regime
#'   probabilities.
#' @references
#' Bagozzi, B.E. and Marchetti, K. (2017). Distinguishing occasional abstention
#' from routine indifference in models of vote choice. \emph{Political Science
#' Research and Methods}, 5(2), 277-294. \doi{10.1017/psrm.2015.42}
#'
#' Gaudry, M.J.I. and Dagenais, M.G. (1979). The dogit model.
#' \emph{Transportation Research Part B}, 13, 105-111.
#' \doi{10.1016/0191-2615(79)90028-6}
#'
#' Swait, J. and Ben-Akiva, M. (1987). Empirical test of a constrained choice
#' discrete model: mode choice in Sao Paulo, Brazil. \emph{Transportation
#' Research Part B}, 21, 103-115. \doi{10.1016/0191-2615(87)90010-5}
#' @examples
#' set.seed(2)
#' d <- rmnl(800, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
#'           delta = c(0.5, 1.2), inflate = "a")
#' m <- imnl(y ~ x1 + x2 | z1, data = d, inflate = "a")
#' summary(m)
#' ## the inflated model nests the multinomial logit at a boundary (phi = 1): compare the
#' ## two by information criteria here, or by boot_lr_test(m0, m)
#' m0 <- mnl(y ~ x1 + x2, data = d)
#' compare_models(m0, m)
#' @seealso [mnl()], [predict.choicemodel()], [first_difference()], [vuong()]
#' @family estimators
#' @export
imnl <- function(formula, data, inflate = NULL, id = NULL, alt = NULL, ref = NULL, link = c("logit", "probit"),
                 weights = NULL, offset = NULL, fe = NULL,
                 se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"), cluster = NULL, nboot = 200L,
                 cores = 1L, start = NULL, maxit = 1000L, reltol = 1e-10, crosscheck = FALSE) {
  if (!is.null(fe))
    stop("imnl() does not take unit fixed effects: with a constant per unit in the outcome equations the split equation is identified only ",
         "through the functional form, units that always choose the inflated alternative drop out, and the split coefficients come out ",
         "heavily biased towards phi = 1 at any usual panel length (see ?imnl). Fit without 'fe' and with 'cluster' set to the unit, ",
         "or use mnl(fe = ) for the outcome equations.", call. = FALSE)
  se_missing <- missing(se); se <- match.arg(se); link <- match.arg(link)
  .pc_fit(formula, data, inflated = TRUE, inflate = inflate, ref = ref, link = link, id = id, alt = alt,
          weights = weights, offset = offset, cluster = cluster, fe = fe, se = se, nboot = nboot, cores = cores,
          start = start, maxit = maxit, reltol = reltol, call = match.call(), se_missing = se_missing,
          crosscheck = crosscheck)
}

#' Nested logit regression
#'
#' Fits the two-level nested logit (McFadden 1978) in its utility-maximization-
#' consistent form, the one \pkg{mlogit} and Stata's `nlogit` (RUM-consistent)
#' use: alternatives are partitioned into nests, and each nest \eqn{m} has an
#' inclusive-value parameter \eqn{\lambda_m} that governs the correlation of
#' the errors within the nest (\eqn{\lambda_m = 1}: no correlation, the
#' multinomial logit; values in \eqn{(0, 1)} are consistent with utility
#' maximization for all covariate values).
#'
#' @details
#' With \eqn{V_{ij}} as in [mnl()] and \eqn{I_{im} = \log\sum_{k \in B_m} e^{V_{ik}/\lambda_m}}
#' over the alternatives of nest \eqn{m} available to case \eqn{i},
#' \deqn{P(y_i = j) = \frac{e^{V_{ij}/\lambda_m}}{\sum_{k \in B_m} e^{V_{ik}/\lambda_m}}
#'   \cdot \frac{e^{\lambda_m I_{im}}}{\sum_l e^{\lambda_l I_{il}}}, \qquad j \in B_m.}
#' The inclusive-value parameters are estimated on the log scale and reported
#' on the natural scale with delta-method standard errors; the z statistics of
#' the summary test \eqn{\lambda = 1} (the multinomial logit within that nest).
#' Test IIA jointly with `lr_test(mnl_fit, nlogit_fit)` (the restriction
#' \eqn{\lambda = 1} is interior). A singleton nest's parameter does not enter
#' the likelihood and is fixed at 1 (with a message). Everything else --
#' layouts, weights, offsets, fixed effects, standard errors -- is as in
#' [mnl()]. Starting values come from the multinomial logit fit; the fit is
#' run from \eqn{\lambda = 1} and \eqn{\lambda = 0.7} and the better optimum
#' kept.
#'
#' @inheritParams mnl
#' @param nests A named list of character vectors assigning every alternative
#'   (response level) to exactly one nest, e.g.
#'   `list(fast = c("air", "train"), slow = c("bus", "car"))`.
#' @param equal_lambda If `TRUE`, one inclusive-value parameter is shared by
#'   all (non-singleton) nests.
#' @return An object of class `c("nlogit", "choicemodel")`; as [mnl()], with
#'   the inclusive-value parameters named `"lambda:<nest>"` (or `"lambda"`).
#'   `predict(type = "prob_conditional")` gives the within-nest probabilities.
#' @references
#' McFadden, D. (1978). Modelling the choice of residential location. In
#' A. Karlqvist et al. (eds.), \emph{Spatial Interaction Theory and Planning
#' Models}, 75-96. North-Holland.
#'
#' Train, K.E. (2009). \emph{Discrete Choice Methods with Simulation}, 2nd ed.,
#' chapter 4.
#' @examples
#' set.seed(3)
#' d <- rnlogit(1500, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5), d = c(0.2, 0.4, 0.3)),
#'              nests = list(n1 = c("a", "b"), n2 = c("c", "d")), lambda = c(0.6, 0.8))
#' m <- nlogit(y ~ x1 + x2, data = d, nests = list(n1 = c("a", "b"), n2 = c("c", "d")))
#' summary(m)
#' lr_test(mnl(y ~ x1 + x2, data = d), m)
#' @seealso [mnl()], [hev()], [lr_test()]
#' @family estimators
#' @export
nlogit <- function(formula, data, nests, equal_lambda = FALSE, id = NULL, alt = NULL, ref = NULL, weights = NULL, offset = NULL,
                   fe = NULL, se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"), cluster = NULL, nboot = 200L,
                   cores = 1L, start = NULL, maxit = 1000L, reltol = 1e-10) {
  se_missing <- missing(se); se <- match.arg(se)
  if (missing(nests)) stop("'nests' is required: a named list assigning every alternative to a nest.", call. = FALSE)
  .pc_fit(formula, data, inflated = FALSE, ref = ref, id = id, alt = alt, weights = weights, offset = offset,
          cluster = cluster, fe = fe, se = se, nboot = nboot, cores = cores, start = start, maxit = maxit,
          reltol = reltol, call = match.call(), se_missing = se_missing,
          family = "nlogit", str = list(nests = nests, equal = isTRUE(equal_lambda)))
}

#' Heteroskedastic extreme value regression
#'
#' Fits the heteroskedastic extreme value (HEV) model of Bhat (1995):
#' independent Gumbel errors whose scale \eqn{\theta_j} differs across
#' alternatives (the base alternative's scale is fixed at 1), so that the
#' substitution pattern is no longer proportional and IIA is relaxed through
#' the scales rather than through nests.
#'
#' @details
#' The choice probability is a one-dimensional integral over the standard
#' Gumbel variable \eqn{w = \epsilon_j/\theta_j},
#' \deqn{P(y_i = j) = \int \exp\Big(-\sum_{k \ne j} e^{-(V_{ij} - V_{ik})/\theta_k - (\theta_j/\theta_k) w}\Big) e^{-w - e^{-w}}\, dw,}
#' evaluated on a window around the mode of the (log-concave) integrand: by the
#' trapezoidal rule when every scale ratio \eqn{\theta_j/\theta_k} is at most
#' 20 (the integrand is analytic with exponential tails, so the rule converges
#' spectrally), and by adaptive Gauss-Kronrod quadrature beyond that ratio,
#' where a competitor with a much smaller scale makes the integrand a near
#' step (the initial panels carry breakpoints around each step, so the
#' log-likelihood and its gradient are smooth in the parameters); probabilities agree with `integrate()` to within 1e-9 (the test tolerance) for
#' scale ratios up to 1.6e5. With every \eqn{\theta_j = 1} the model is the multinomial logit
#' (reproduced exactly); `lr_test(mnl_fit, hev_fit)` tests homoskedasticity.
#' The scales are estimated on the log scale and reported on the natural scale
#' with delta-method standard errors; the z statistics in the summary test
#' \eqn{\theta = 1}, and [confint()] gives intervals on the log scale
#' transformed back. A scale that runs to zero (the alternative's utility
#' becomes deterministic relative to the others) is flagged as a boundary fit
#' with the information matrix singular; the likelihood is then usually flat in
#' that direction and the standard errors are not interpretable (the
#' travel-mode data of the vignette are an example). Where it still rises at
#' the search's bound (a log scale of +/- 8), the warning says so, and the
#' reported log-likelihood is below the model's supremum. With two alternatives the single
#' scale is identified only through the functional form and a warning says so.
#' Everything else is as in [mnl()].
#'
#' @inheritParams mnl
#' @return An object of class `c("hev", "choicemodel")`; as [mnl()], with the
#'   scale parameters named `"scale:<alternative>"`.
#' @references
#' Bhat, C.R. (1995). A heteroscedastic extreme value model of intercity
#' travel mode choice. \emph{Transportation Research Part B}, 29, 471-483.
#' \doi{10.1016/0191-2615(95)00015-6}
#' @examples
#' set.seed(4)
#' d <- rhev(300, beta = cbind(b = c(0.5, -1, 0.8), c = c(-0.3, 0.6, -0.5)),
#'           scale = c(1, 0.6, 1.5), gamma = 1)   # w1 identifies the scales
#' m <- hev(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
#' summary(m)
#' @seealso [mnl()], [nlogit()], [lr_test()]
#' @family estimators
#' @export
hev <- function(formula, data, id = NULL, alt = NULL, ref = NULL, weights = NULL, offset = NULL, fe = NULL,
                se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"), cluster = NULL, nboot = 200L,
                cores = 1L, start = NULL, maxit = 1000L, reltol = 1e-10) {
  se_missing <- missing(se); se <- match.arg(se)
  .pc_fit(formula, data, inflated = FALSE, ref = ref, id = id, alt = alt, weights = weights, offset = offset,
          cluster = cluster, fe = fe, se = se, nboot = nboot, cores = cores, start = start, maxit = maxit,
          reltol = reltol, call = match.call(), se_missing = se_missing, family = "hev", str = list())
}

#' Multinomial probit regression
#'
#' Fits the multinomial probit (Hausman and Wise 1978; Daganzo 1979):
#' utilities \eqn{U_{ij} = V_{ij} + \epsilon_{ij}} with multivariate normal
#' errors whose covariance is estimated, so that any pattern of correlation
#' and heteroskedasticity across alternatives, and hence any substitution
#' pattern, is admitted.
#'
#' @details
#' Only utility differences are identified. The model is parameterized through
#' the covariance of the differences from the base alternative,
#' \eqn{\eta_k = \epsilon_k - \epsilon_b}, \eqn{\tilde\Sigma = LL'} with
#' \eqn{L} lower triangular, its first element fixed at 1 (the scale
#' normalization) and the remaining \eqn{J(J-1)/2 - 1} elements free (Train
#' 2009, section 5.2; the same normalization as \pkg{mlogit}). The reported
#' structural parameters are the elements of \eqn{\tilde\Sigma} named
#' `"var:<alt>"` and `"cov:<alt>,<alt>"`, with delta-method standard errors.
#'
#' Probabilities: with three alternatives the exact bivariate normal (the
#' Drezner-Genz algorithm, absolute error about 1e-15); with more, the GHK
#' simulator with `draws` Halton draws per case, held fixed across evaluations
#' (the simulated log-likelihood is smooth away from a singular covariance;
#' the simulation error of a
#' probability has standard deviation of order \eqn{1/\sqrt{R}} and the bias of
#' its logarithm is of order \eqn{1/R}, both vanishing with the number of
#' draws; GHK probabilities of a case need not sum exactly to one, and the
#' draws of a case are tied to its row position, so predictions for the same
#' case in another row use other draws). Gradients are analytic throughout,
#' including the covariance parameters (forward-mode differentiation through
#' the Cholesky factor and the GHK recursion); the Hessian is obtained by
#' central differences of the analytic gradient. The multinomial probit does not nest the
#' multinomial logit, so compare them with [vuong()] or the information criteria. Multinomial probit
#' coefficients are on a smaller scale than multinomial logit coefficients (the
#' error-difference variance is normalized to 1 rather than \eqn{\pi^2/3}).
#' Starting values come from the multinomial logit fit (scaled) with the
#' covariance implied by iid errors, and from independent differences.
#' Everything else is as in [mnl()], except that the returned object has class
#' `"mnprobit"` (not `"mnp"`, which the CRAN package MNP uses).
#'
#' Identification. With case-specific covariates only (alternative-specific
#' coefficients, the wide layout) the covariance parameters are identified
#' solely through the functional form and are fragile (Keane 1992): the
#' likelihood is nearly flat along directions where an alternative's
#' coefficients and its error variance shrink together, fits may stop at a
#' near-singular covariance with a convergence warning, and the covariance
#' estimates are not to be interpreted. Alternative-specific covariates with
#' generic coefficients (the long layout) identify the covariance; a fit with
#' a variance running to zero is flagged as a boundary fit. Near a singular
#' covariance the GHK simulated likelihood is not smooth: with the draws held
#' fixed, one draw's contribution can jump over a small change in the
#' parameters, so a boundary fit with four or more alternatives may stop with a
#' large gradient, and the boundary warning says that the convergence check is
#' unreliable there.
#'
#' @inheritParams mnl
#' @param draws Number of GHK draws per case for four or more alternatives
#'   (default 200). The simulator's Halton sequences use one prime per
#'   dimension of the differenced utilities, the first 20, so the model takes
#'   at most 21 alternatives.
#' @return An object of class `c("mnprobit", "choicemodel")`.
#' @references
#' Daganzo, C. (1979). \emph{Multinomial Probit: The Theory and Its Application
#' to Demand Forecasting}. New York: Academic Press.
#'
#' Hausman, J.A. and Wise, D.A. (1978). A conditional probit model for
#' qualitative choice. \emph{Econometrica}, 46, 403-426. \doi{10.2307/1913909}
#'
#' Train, K.E. (2009). \emph{Discrete Choice Methods with Simulation}, 2nd ed.,
#' chapter 5.
#'
#' Keane, M.P. (1992). A note on identification in the multinomial probit
#' model. \emph{Journal of Business and Economic Statistics}, 10, 193-200.
#' \doi{10.1080/07350015.1992.10509898}
#' @examples
#' set.seed(5)
#' d <- rmnp(1500, beta = cbind(b = c(0.3, -0.6, 0.5), c = c(-0.2, 0.4, -0.3)),
#'           sigma = matrix(c(1, 0.6, 0.6, 1.5), 2), gamma = 1)
#' m <- mnprobit(chosen ~ w1 | x1 + x2, data = d, id = "id", alt = "alt")
#' summary(m)
#' @seealso [mnl()], [nlogit()], [hev()], [vuong()]
#' @family estimators
#' @export
mnprobit <- function(formula, data, id = NULL, alt = NULL, ref = NULL, weights = NULL, offset = NULL, fe = NULL,
                     se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"), cluster = NULL, nboot = 200L,
                     cores = 1L, start = NULL, maxit = 1000L, reltol = 1e-10, draws = 200L) {
  se_missing <- missing(se); se <- match.arg(se)
  .pc_fit(formula, data, inflated = FALSE, ref = ref, id = id, alt = alt, weights = weights, offset = offset,
          cluster = cluster, fe = fe, se = se, nboot = nboot, cores = cores, start = start, maxit = maxit,
          reltol = reltol, call = match.call(), se_missing = se_missing, family = "mnp", str = list(draws = draws))
}

#' Multinomial logit with random alternative-specific intercepts
#'
#' Fits the multinomial (or conditional) logit of [mnl()] with a normally
#' distributed random intercept per unit in every non-base equation: a mixed
#' logit whose only random coefficients are the constants, for panels and other
#' grouped data where the cases of a unit share unobserved tastes for the
#' alternatives.
#'
#' @details
#' \deqn{V_{ij} = o_{ij} + x_i^\top \beta_j + w_{ij}^\top \gamma + u_{g(i) j}, \quad u_g \sim N(0, \mathrm{diag}(\sigma_j^2)),}
#' with \eqn{u_{g b} = 0} for the base alternative and independence across
#' units and alternatives. With three or more alternatives the base defines
#' the model, not only the labels. Intercepts independent relative to one base
#' are correlated relative to another (\eqn{u_j - u_c} and \eqn{u_k - u_c}
#' share \eqn{u_c}), which the diagonal model cannot represent. Fits with
#' different `ref` are therefore different models, with different
#' log-likelihoods and different tests against the multinomial logit. In the example
#' below the log-likelihood is -488.49 with base `a` and -490.50 with base
#' `b` (a fit whose standard deviation for `a` runs to zero, with a boundary
#' warning), and the likelihood ratio against [mnl()] is 4.53 and 0.50. Choose the
#' base before testing, on substantive
#' grounds: the alternative whose utility is least expected to vary across
#' units (an outside option, abstention). Fits with different bases have the
#' same number of parameters and can be compared descriptively by their
#' log-likelihoods, but a base chosen that way makes the test against the
#' multinomial logit a maximum over the bases, for which the chi-square reference does
#' not hold.
#' The likelihood of unit \eqn{g} integrates the
#' product of the MNL probabilities of its cases over \eqn{u_g} by adaptive
#' product Gauss-Hermite quadrature (Liu and Pierce 1994; Pinheiro and Bates
#' 1995): the nodes are centred at the mode of each unit's integrand and scaled
#' by its curvature, with `nodes` points per dimension (\eqn{J - 1}
#' dimensions, so \eqn{\mathrm{nodes}^{J-1}} points per unit; a warning above
#' 20000 points). The standard deviations are reported as `"sd:<alt>"`
#' (estimated on the log scale, delta-method standard errors). The observation
#' is the unit: `object$lli` has one entry per unit, robust and outer-product
#' standard errors are over units, `cluster` must nest units, frequency weights
#' must be constant within a unit, and `se = "bootstrap"` resamples units.
#' Fitted and predicted probabilities are population-averaged (integrated over
#' the random intercepts). [lr_test()] against the [mnl()] fit of the same
#' covariates tests \eqn{\sigma = 0}, a boundary restriction whose chi-square
#' p-value is conservative; [vuong()] against an MNL aggregates its
#' contributions over units. At \eqn{J = 2} the model is the binomial
#' generalized linear mixed model with a random intercept, and the fit matches
#' `lme4::glmer()` with adaptive quadrature. Everything else is as in [mnl()],
#' except that there is no `fe` or `bias_correct` (the unit effects are the
#' random intercepts) and that `simulate()` draws one intercept vector per unit
#' and then the unit's cases given it, while `fitted()` and `predict()` are the
#' population-averaged probabilities. The IIA tests and [hl_test()] refuse these
#' fits (the cases of a unit are dependent). [fit_stats()] and the information
#' criteria use the number of cases as \eqn{n}.
#'
#' @inheritParams mnl
#' @param se Standard errors: `"analytic"` (default), `"opg"`, `"robust"`
#'   (over units, the observation of this likelihood), `"cluster"` (over
#'   clusters that nest the units), `"bootstrap"` (units resampled, a unit drawn
#'   twice counting as two units with their own intercepts; with frequency
#'   weights, multinomial unit counts), or `"none"`.
#' @param unit Column name in `data` identifying the units whose cases share
#'   the random intercepts.
#' @param nodes Gauss-Hermite points per dimension (default 12). Because the rule is
#'   adaptive its error does not grow with the number of cases per unit:
#'   against `lme4::glmer(nAGQ = 25)` at \eqn{J = 2} (150 units, 6 to 50 cases
#'   each, sd 1 or 2) 12 points put the log-likelihood within 1.2e-3 and the
#'   standard deviation within 0.01 percent, 8 points within 0.02 and 0.15
#'   percent (`verification/remnl_adaptive_check.R` in the development
#'   repository). Raise it for large standard deviations with short panels. The
#'   analytic gradient is the quadrature of the exact gradient, which differs
#'   from the derivative of the quadrature likelihood by the quadrature error;
#'   with few nodes (6 or fewer), or large standard deviations and short
#'   panels, that difference can keep the optimizer from meeting its criterion
#'   (`converged = FALSE` with a small gradient), and more nodes cure it.
#'   The population-averaged probabilities of `fitted()` and `predict()`, an
#'   integral over the prior that no data narrow, use a finer rule: Gauss-Hermite
#'   with \eqn{\max(\mathrm{nodes}, 10\sigma^2)} points per dimension (rounded up
#'   to a power of two) while that is at most 128, and beyond that the
#'   trapezoidal rule with a step of at most \eqn{1/\sigma} standard deviations,
#'   whose size grows only linearly in \eqn{\sigma}; both are accurate to about
#'   \eqn{10^{-7}} or better. The rule has at most 50000 points in all, or
#'   \eqn{\mathrm{nodes}^{J-1}} if that is larger, and at the cap the rule with
#'   the smaller error bound is used. The cap binds at \eqn{J = 3} for
#'   \eqn{\sigma} above about 7 (the error bound passes \eqn{10^{-6}} near 18),
#'   at \eqn{J = 4} above about 1.8 (bound above \eqn{10^{-6}} from 2.7), and
#'   at \eqn{J \ge 5} above about 0.9 (from about 1.6); the fit warns, with the
#'   bound and the number of points, whenever the bound exceeds \eqn{10^{-6}}.
#' @return An object of class `c("remnl", "choicemodel")`.
#' @references
#' Liu, Q. and Pierce, D.A. (1994). A note on Gauss-Hermite quadrature.
#' \emph{Biometrika}, 81(3), 624-629. \doi{10.1093/biomet/81.3.624}
#'
#' Pinheiro, J.C. and Bates, D.M. (1995). Approximations to the log-likelihood
#' function in the nonlinear mixed-effects model. \emph{Journal of
#' Computational and Graphical Statistics}, 4(1), 12-35.
#' \doi{10.1080/10618600.1995.10474663}
#' @examples
#' set.seed(1)
#' G <- 60; T <- 5; unit <- rep(seq_len(G), each = T)
#' x1 <- rnorm(G * T); u <- rnorm(G, 0, 1)
#' V <- cbind(0, 0.3 + 0.8 * x1 + u[unit], -0.2 - 0.6 * x1)
#' P <- exp(V) / rowSums(exp(V))
#' y <- apply(P, 1, function(p) sample.int(3, 1, prob = p))
#' d <- data.frame(y = factor(c("a", "b", "c")[y]), x1, unit = factor(unit))
#' m <- remnl(y ~ x1, data = d, unit = "unit", nodes = 8)
#' summary(m)
#' lr_test(mnl(y ~ x1, data = d), m)
#' @seealso [mnl()], [nlogit()], [hev()], [mnprobit()], [lr_test()]
#' @family estimators
#' @export
remnl <- function(formula, data, unit, nodes = 12L, id = NULL, alt = NULL, ref = NULL, weights = NULL, offset = NULL,
                  se = c("analytic", "opg", "robust", "cluster", "bootstrap", "none"), cluster = NULL, nboot = 200L,
                  cores = 1L, start = NULL, maxit = 1000L, reltol = 1e-10) {
  se_missing <- missing(se); se <- match.arg(se)
  if (missing(unit)) stop("'unit' is required: the column of 'data' that identifies the units sharing random intercepts.", call. = FALSE)
  .pc_fit(formula, data, inflated = FALSE, ref = ref, id = id, alt = alt, weights = weights, offset = offset,
          cluster = cluster, fe = NULL, se = se, nboot = nboot, cores = cores, start = start, maxit = maxit,
          reltol = reltol, call = match.call(), se_missing = se_missing, family = "remnl", str = list(unit = unit, nodes = nodes))
}
