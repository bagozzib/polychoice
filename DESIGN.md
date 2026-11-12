# polychoice: design document

Working title: **polychoice** (polytomous choice). `polytomous` was a CRAN package
(Arppe, archived 2015-10-29) and the CRAN Repository Policy forbids a name that
matches any current or past CRAN package, so it cannot be reused. `polychoice`,
`polychotomous`, `nomchoice`, and `multichoice` were free on CRAN on
2026-10-13. Renaming later is one folder rename plus one `sed` over
DESCRIPTION, `R/polychoice-package.R`, `tests/testthat.R`, README, and this
file: the name is deliberately kept out of every function and class name
(internal prefix `.pc_`, C++ prefix `pc_`, S3 base class `"choicemodel"`).

Started 2026-10-13. Sibling of the `underdisp` (count) and `iop` (ordered)
packages; same author, same engineering conventions, same review protocol.
Package and git repository: `Underdispersion/polychoice/polychoice/`, with this
document and `verification/` inside it (both excluded from the built package).

"The ledger" below is the author's working verification ledger, a dated
record of every check, review round and stress run. The public repository's
`verification/ledger.md` summarizes the entries that the vignettes and help
pages cite and records the release build's frozen run; the rest of the
working ledger, the stress battery and the one-off probes are kept by the
author.

## 1. Scope

### Phase 1 (this build): the MNL core, the inflated MNL, the existing IIA battery
* `mnl()`: multinomial / conditional logit by maximum likelihood with
  case-specific covariates (alternative-specific coefficients), alternative-
  specific covariates (generic coefficients), alternative-specific constants,
  varying choice sets, frequency weights, offsets, unit fixed effects (dummies),
  analytic / robust / cluster-robust / bootstrap covariances.
* `imnl()`: the baseline-inflated multinomial logit of Bagozzi and Marchetti
  (2017, PSRM), generalised to any inflated category and either layout, with a
  logit or probit split equation; a multistart from the MNL fit.
* One engine, one object class, one set of methods: `print`, `summary`,
  `coef`, `vcov`, `logLik`, `nobs`, `fitted`, `residuals`, `predict`,
  `confint`, `update`, `simulate`, broom (`tidy`/`glance`/`augment`), texreg.
* Quantities of interest: `predict()` for observed, outcome-stage, and regime
  probabilities; `first_difference()` (delta-method or simulation intervals,
  with the `stage` switch of `iop`); `rrr()` relative risk ratios.
* Model fit and comparison: `fit_stats()` (McFadden, Cox-Snell, Nagelkerke,
  count and adjusted count R-squared, AIC, BIC), `classification()`,
  `lr_test()`, `vuong()` (with AIC/BIC corrections), `compare_models()`.
* IIA / misspecification tests, existing tests only:
  * `hausman_mcfadden()` (Hausman and McFadden 1984), with the classic
    variance-difference form and the seemingly-unrelated-estimation joint
    covariance (Weesie 1999; Stata `mlogtest, suest`), dropping any alternative
    including the base (handled by an exact re-basing of the full model).
  * `small_hsiao()` (Small and Hsiao 1985) in the correct form: full model on
    both halves, the sqrt(1/2) weighting, restricted model on the second half;
    optionally repeated over several random splits.
  * `fp_tests()`: Fok and Paap (2025, J. Choice Modelling) composite-likelihood
    Hausman test (Q1; default covariance the joint influence-function estimate,
    F&P's variance difference as an option, section 5) and GMM
    overidentification test (Q2; default Hansen's two-step efficient GMM, the
    iterated form of F&P and the projected form at the MLE as options), base
    pairs or all pairs.
  * `iia_tests()`: the battery in one table.
* Panel tools: `mundlak()`.
* Data: `pres2004` (Campbell and Monson 2008 as analysed by Bagozzi and
  Marchetti 2017, Replication II) and `pa2006` (Arceneaux and Kolodny 2009,
  Replication I), so the published fits are regression tests.

### Phase 2 (in progress from 2026-10-13, second session)
Done:
* Family abstraction: every model shares the design, the beta / gamma blocks,
  the optimizer, covariances, bootstrap, and methods; a family adds a
  structural block (internal scale, natural-scale reporting by the delta
  method) and its own C++ kernel; families without an analytic Hessian use
  central differences of the analytic gradient in the standardized
  coordinates.
* `nlogit()`: two-level nested logit in the RU2 / RUM-consistent form of
  mlogit and Stata; free or shared inclusive-value parameters; singleton nests
  fixed at 1; exact agreement with mlogit on TravelMode. The nest-level
  term lambda_m I_m is split as Vmax_m + lambda_m log n_m + O(1) and the nests
  are compared part by part: the summed form rounded log P to exactly 0 at
  lambda = e^40 (with equal lambdas and nest sizes the O(lambda) parts cancel
  and the O(1) parts decide), a false optimum the equal-lambda fit found.
  mlogit reports BHHH (outer-product) standard errors for the nested logit and
  HEV; `se = "opg"` reproduces them, the default remains the observed
  information.
* `hev()`: heteroskedastic extreme value. Bhat's Gauss-Laguerre scheme in
  u = e^{-w} was implemented first and found to converge only algebraically
  when scale ratios differ from 1 (u^r has a singular derivative at 0 for
  r < 1): 40 vs 80 nodes differed by 0.5 in a 200-case log-likelihood.
  Replaced by quadrature in the Gumbel variable w. The log-integrand is
  strictly concave, so per case and alternative its mode is found by
  bisection and the window is where the log-integrand has fallen 36 below
  its maximum (at most 200 units of w either side). When every scale ratio
  theta_j / theta_k is at most 20 the trapezoidal rule with step 0.35 / r_max
  converges spectrally (analytic integrand, exponential tails) and is used;
  a competitor with a much smaller scale makes exp(-h_k) a near step of width
  1 / r_k, which a uniform rule needs O(r_k) nodes to resolve, so beyond
  ratio 20 the window is integrated by adaptive Gauss-Kronrod 7-15 panels
  (bisect the panel with the largest error until the integrand and two
  derivative companions are settled to 1e-10; at most 2048 panels), O(log r)
  panels per step, starting from panels with breakpoints at
  w*_k + {-4, -1, 0, 1, 4, 12, 40} / r_k around each step w*_k = -a_k / r_k
  of a competitor with r_k > 1. Checked against R's integrate() to the test tolerance of 1e-9 for
  scale ratios up to 1.6e5 (test-families.R). Three earlier rules failed: a window anchored at the root of
  S(w) = 1 scaled with 1 / theta_j and asked for 1e11 nodes when the
  optimizer pushed a log-scale towards -20 (a 22 GB process), and a bounded
  trapezoid grid with the ratio capped at 60 under-resolved the step and let
  the TravelMode fit run to a spurious optimum at a car scale of 0.004, and
  the adaptive rule started from uniform 5-unit panels could leave a step in
  a node-free sliver at the end of a panel after a bisection, where both
  rules missed it and the error estimate stayed small: the integral came out
  short by up to 1e-5 at some parameter values and not at others, so the
  log-likelihood and its gradient jumped and the finite-difference Hessian
  of the TravelMode boundary fit showed spurious negative curvature (found by
  the saddle check, 2026-10-14). The error companion of the scale derivative
  uses sqrt(1 + e^2), not |e| (whose kink sat at the peak of the spike).
  log P is returned by log-sum-exp, finite where P underflows. mlogit's HEV
  is not an oracle: on TravelMode its reported log-likelihood at its own
  solution is 3.8 above the exact value (ledger), so the kernel is checked
  against R's integrate() and the boundary optimum instead. A log scale (or
  log inclusive value) beyond +/- 7 raises the boundary flag.
* `mnprobit()`: multinomial probit; covariance of the differences from the
  base through a Cholesky factor with the (1,1) element fixed at 1; exact
  bivariate normal for three alternatives (ported Drezner-Genz routine from
  iop), GHK with Halton draws beyond, forward-mode derivatives through the
  Cholesky differential and the GHK recursion (gradient to 3e-9 of numDeriv
  with fixed draws; eight times faster than the finite-difference version).
  Near a singular covariance the GHK simulated likelihood is not smooth: a
  truncation point is divided by a Cholesky diagonal near zero, so one
  draw's term switches between about 0 and 1 over a small change in the
  parameters (stress finding F24). The boundary warning of a fit with four
  or more alternatives says so, and the probit's warning points to
  `vuong()`, since `lr_test()` refuses the probit against the logit.
* Further existing tests: `combine_test()` (Wald and Cramer-Ridder LR),
  `term_test()` (Wald / LR per term), `hl_test()` (Fagerland-Hosmer-Bofin),
  `nested_lm_test()` (score test of IIA against a nested logit at
  lambda = 1), and `boot_pvalue()` (parametric-bootstrap calibration of a
  test whose null is the fitted model; round 7 refuses `term_test()` and
  `combine_test()`, whose null is a restriction of it). `lr_test()` accepts
  mnl vs nlogit / hev.
* `se = "opg"` (BHHH) for every estimator; simulators `rmnl()`, `rnlogit()`,
  `rhev()`, `rmnp()` with `gamma` (alternative-specific covariates, long-layout
  output), because the structural parameters of the nested logit, HEV, and
  probit are weakly identified from case-specific covariates alone (ledger;
  Keane 1992 for the probit).
* Split-panel jackknife for fixed-effects slopes (`mnl(fe = , bias_correct =
  "jackknife")`, design below). The gate was redesigned and recalibrated on
  2026-10-14: random within-unit placebo splits, kappa = 2.5, alpha = 0.01
  (selection rule, selection run, and holdout in the ledger); the jackknife
  over-corrects at T = 6 and 10 but lowers the RMSE and improves Wald
  coverage at every T studied (ledger).
* Random alternative-specific intercepts: `remnl(formula, data, unit, nodes)`,
  family `"remnl"`, kernel `src/remnl_ll.cpp`, implemented as designed below
  (the observation is the unit: `.pc_obs()` supplies unit weights and
  clusters to the covariance code, the bootstrap resamples units and relabels
  repeated draws, `vuong()` aggregates an MNL's contributions over units,
  `lr_test()` admits mnl within remnl with a boundary note, the IIA tests
  refuse non-mnl families). Verification (glmer oracle at J = 2, numDeriv,
  MNL reduction, recovery) in tests/testthat/test-re.R; results in the ledger.
* Four vignettes (getting started; the inflated model; testing IIA and
  specification; beyond IIA), built with knitr/rmarkdown.
Still to do:
* The new IIA tests: closed (author's decision 2026-10-18; section 5.1). Done since: the gate recalibration
  (2026-10-14) and three rounds of the simulated-referee review (nine
  perspectives, as for iop), recorded in the ledger.

Designs for the two remaining estimators (written before the code):
* Split-panel jackknife (`bias_correct = "jackknife"` with `fe`), mirroring
  `underdisp::cpb_fe()`: the unit dummies carry an incidental-parameters bias
  of order 1/T in the slope coefficients; the model is refit on the first and
  second halves of every unit's cases (rows split in the order supplied,
  which should be temporal order) and the corrected slopes are
  2 * full - mean(halves) (Dhaene and Jochmans 2015). The half designs are
  rebuilt from the full design (`.pc_half_design()`): the fixed-effect columns
  are regenerated on the subset so never-chosen pairs and empty dummies are
  handled exactly as in a fresh fit. Validity gate, refusing the correction
  (returning the ML fit with a warning naming the check) when: the two
  half-panel slope vectors disagree beyond kappa = 2.5 (calibrated, ledger)
  times the mean disagreement of five random splits of every unit's cases into
  halves of the time split's sizes (placebo noise, on the time split's footing
  under time homogeneity; asymmetric heterogeneity; the first design used five
  deterministic unit splits with kappa 2.5, then 4),
  or the full fit's residuals y_ij - P_ij trend with the within-unit time
  order for some alternative (Bonferroni over alternatives at 0.01; symmetric
  heterogeneity, the cure is a time control in the formula). Warns when the
  median panel length is below 6. Corrected values replace `coefficients`
  (slopes only; the unit constants are nuisance), the ML values are kept in
  `$uncorrected`, standard errors are the ML ones, `logLik` refers to the ML
  fit, and the Monte Carlo of the bias and of the gate is a verification item
  (size of the refusal on stationary panels, refusal on regime changes).
* Random alternative-specific intercepts (`re = "unit"`, family `"remnl"`):
  V_ij = ... + u_gj for unit g, u_g ~ N(0, diag(sigma_j^2)) over the non-base
  alternatives, independent across units and alternatives (a correlated
  version can follow through a Cholesky block). The likelihood is per unit:
  L_g = int prod_{i in g} P_i(y_i | u) phi(u) du by adaptive product
  Gauss-Hermite (nodes centred at the mode of each unit's log-integrand,
  found by damped Newton, and scaled by the Cholesky factor of its curvature;
  `nodes`, default 12 per dimension; a warning above 20000 points), with the
  log-sum-exp over nodes and the gradient as the mixture-weighted average of
  the per-node derivatives with the nodes held fixed. The first version used
  a fixed rule centred at zero; it aliased once the integrand was narrower
  than the node spacing (long panels, large sigma) and was replaced. The engine's observation becomes the unit for
  this family: `lli` and `scores` are per unit, frequency weights must be
  constant within a unit, `cluster` must nest units, and `vuong()` against an
  MNL sums the MNL contributions within units. `predict()` returns the
  population-averaged probabilities (integrated over u). Oracles: at J = 2
  the model is a binomial GLMM with a random intercept, so
  `lme4::glmer(nAGQ = 25)` must be matched to 1e-6 in the log-likelihood and
  estimates; reduction to the MNL as sigma -> 0; derivatives against
  numDeriv; recovery; mlogit's simulated mixed logit (`rpar`, `panel = TRUE`)
  as a loose check at J = 3.

Out of scope: mixed logit with random coefficients (apollo, mixl, gmnl,
logitr), Bayesian estimation (MNP, bayesm), ordered alternatives (iop).
Ranked alternatives moved into scope on 2026-10-19 (Phase 3).

### Phase 3 (from 2026-10-19, author's decision): captivity, rankings, consideration sets

Status (2026-10-19): A0 done; B (rologit) implemented; A (several inflated
alternatives, captivity_test()) implemented; their Monte Carlo studies wait
for a frozen build. C stopped at Gate C1, which failed in every design
(ledger); the latent consideration-set model is not implemented, and the
author decides on any further work. Gate C1 covered case-specific covariates
only; a design with a voter-party (alternative-specific) covariate was not
tested (see C below).

Designs written before any code. Order of work: (A0) the citation of the
captivity literature for `imnl()`; (B) the rank-ordered logit, the cheapest
and lowest-risk; (A) several captive alternatives and the captivity score
test; (C) the consideration-set model, which passes two gates before any
package code is written. Every item ends with its Monte Carlo study in
`verification/`, run on the frozen build with the rest.

**A0. What `imnl()` is in the discrete-choice literature.** The dogit model
(Gaudry and Dagenais 1979) adds to the MNL a captive share of the population
for each alternative:

    P*_ij = (P_ij + theta_j) / (1 + sum_{m in A_i} theta_m),   theta_j >= 0

(P_ij the MNL probabilities over the choice set A_i). With captivity to k only,
P*_ik = phi P_ik + (1 - phi), P*_ij = phi P_ij, phi = 1 / (1 + theta_k): the
inflated multinomial logit, and theta_k = exp(-z' delta) gives phi = Lambda(z' delta),
the logit split equation of `imnl()` exactly (derived here; Bagozzi and
Marchetti's model with a logit split is the dogit with covariate-dependent
captivity to one alternative). Swait and Ben-Akiva (1987a, b) made the
captivity parameters functions of covariates (the "parametrized logit
captivity" model) as a case of probabilistic choice-set generation.
Read 2026-10-19 (the author supplied the papers; ledger): the PLC is their
eq. (3) in 1987b, P_i = e^{DX_i} / (1 + sum_l e^{DX_l}) + [1 / (1 + sum_l
e^{DX_l})] P_{i|C}, captivity odds exp(DX_i) with a constant for every
alternative and covariates that may differ by alternative; so `imnl()` with a
logit split is the PLC with captivity to k only, with delta = -D, and
several inflated alternatives give the PLC with the same z in every captivity
equation. The dogit (their eq. 1) has constant theta_i >= 0 summed over all
alternatives; 1987a derives it from independent random constraints (eq. 24).
The dogit was proposed as an IIA relaxation: P*_ij / P*_im depends on the
other alternatives' utilities. Tse (1987) gives a Lagrange multiplier test of
the MNL against the dogit with constant theta (eqs. 4-5, expected
information, chi-square on J df, the boundary null read as a test for
arbitrarily small theta); the joint statistic of `captivity_test()` is his,
extended to varying choice sets and alternative-specific covariates, and the
one-sided z per alternative is added. The documentation cites both papers as
such. Bagozzi and Marchetti (2017) cite neither line of work.

**B. Rank-ordered logit: `rologit()`.** The exploded logit (Beggs, Cardell
and Hausman 1981; Chapman and Staelin 1982; Allison and Christakis 1994; the
Plackett-Luce model with covariates). Case i ranks L_i of its available
alternatives, r_i1 (most preferred) ... r_iL_i; the rest are unranked and
tied below the ranked ones (a partial ranking). With R_il = A_i minus
{r_i1, ..., r_i,l-1}:

    log L_i = sum_{l=1}^{L_i} [ V_{i r_il} - log sum_{m in R_il} exp(V_im) ]

with V_ij as in the MNL (stage l is an MNL choice of r_il from R_il; a stage
with one alternative left contributes 0). Score: sum over stages of the MNL
score of that stage, d log L_i / dV_ij = sum_l 1{j in R_il} (1{j = r_il} -
P_ij^(l)); Hessian: minus the sum over stages of the stage covariance of
dV/dtheta. Globally concave, as the MNL.
* Ties (Allison and Christakis 1994, section 4 and eq. 8, read 2026-10-19):
  tied alternatives are taken to have an unknown strict order, so a tie group
  T of d alternatives, ranked above the set B of the alternatives below it,
  contributes the sum over the d! orders p of T of
  prod_{s=1}^{d} exp(V_{p_s}) / (sum_{t >= s} exp(V_{p_t}) + sum_{m in B} exp(V_m)).
  The unranked alternatives form the bottom group, whose term is 1 (no
  enumeration). This is the exact likelihood; Breslow's approximation, the
  default of proportional-hazards software, attenuates the estimates
  (their Table 1) and is not offered. Middle tie groups larger than 7 are
  refused (5040 orders). Score and Hessian of a tie group: with l_p the
  log-probability of order p, g_p and H_p its derivatives in V, and
  omega_p = exp(l_p) / sum exp(l_q), the group's derivatives are sum omega_p g_p
  and sum omega_p (H_p + g_p g_p') - gbar gbar'.
* Interface: long layout only, `rologit(rank ~ w | x, data, id, alt, ref,
  weights, offset, cluster, se, nboot, cores, start, maxit, reltol, depth)`.
  The response is the rank: smaller is preferred, any increasing numbers
  (equal numbers are ties, as Allison and Christakis code them), 0 =
  unranked (tied at the bottom); a missing rank is a missing value (the case
  drops, as everywhere). `depth` keeps only the first `depth` positions (the
  alternatives ranked below become unranked; a tie group straddling the cut
  stays whole), Chapman and Staelin's explosion depth and Hausman and Ruud's
  top-P estimator: both argue, and find in their applications, that later
  ranks can be noisier (Hausman and Ruud: the pattern depends on the survey). No
  `fe` (refused) in this version.
* Stability of preferences across ranking stages (the LR test of Hausman
  and Ruud 1987, section 3, and the Wald form of Allison and Christakis,
  section 5.4; Chapman and Staelin's 1982 pooling test is the sequential
  version, stages 1..E against E + 1; papers read 2026-10-19):
  `rank_test(fit, after = s)` refits with separate coefficients for the
  stages after position s (V_ij^(l) uses theta_late for l > s; within a tie
  group the position is the order's own) and reports the LR and Wald tests
  of equal coefficients. Allison and Christakis note (their footnote 13)
  that differences between stages can also come from heterogeneity across
  respondents; the documentation says so.
* Implementation: a case-level kernel `src/rologit_ll.cpp` (lli, gradient,
  scores, analytic Hessian, probabilities; per case the derivatives are
  formed in V (a J-vector and a J x J matrix) and mapped to theta through the
  sparse design rows, as in `pc_hess_cpp`), on the MNL design with the ranks
  as an n x J integer matrix in the family's structural slot; the observation
  is the case (the ranking), so weights, clusters, bootstrap, `vuong()` and
  the covariances need no change. `y` holds the top choice. `predict()`
  returns first-choice probabilities (the MNL over the choice set) and
  `type = "ranking"` the probability of the observed ranking; `classification()`
  and `fit_stats()` are about first choices except the log-likelihood-based
  measures (ranking likelihood; null model: the rank-ordered logit with
  alternative-specific constants only, or equal probabilities at every stage
  without constants). `simulate()` draws rankings (Gumbel utilities), and a
  simulator `rrologit()` produces long-layout data with a rank column. The IIA
  battery, `hl_test()` and `combine_test()` refuse a rologit fit;
  `lr_test()` and `vuong()` compare two rologit fits only (a rologit and an
  MNL fit have different responses).
* Oracles: (1) `depth = 1` equals `mnl()` on the first choices exactly
  (log-likelihood and estimates) when there is no tie at the top; (2) without
  ties, the exploded data (one pseudo-case per stage, availability R_il)
  fitted by `mnl()` give identical estimates and log-likelihood; (3) a tie
  group's term equals the brute-force sum over its orders computed in R, and
  a case whose alternatives are all tied contributes 0; (4) J = 2 is the
  binary logit; (5) gradient and Hessian against numDeriv, with and without
  ties; (6) mlogit's rank-ordered logit where its current version supports
  it; (7) recovery and coverage Monte Carlo (`rologit_mc.R`: complete,
  partial, and tied rankings, J = 4, 5; criteria as recovery_coverage_mc:
  coverage 0.92-0.98, bias within 3 MC standard errors) and the size of
  `rank_test()` under the model.

**A. Several captive alternatives (the dogit with covariates) in `imnl()`.**
`inflate` may name several alternatives C (|C| >= 2; logit split only, a
probit split is refused with several). Regimes: choose (0) or captive to c in
C. With eta_ic = z_i' delta_c,

    pi_i0 = 1 / D_i,   pi_ic = exp(-eta_ic) / D_i,   D_i = 1 + sum_{c in C, a_ic = 1} exp(-eta_ic)
    P*_ij = pi_i0 P_ij + 1{j in C} pi_ij

so eta_ic = log(pi_i0 / pi_ic): a positive coefficient makes captivity to c
less likely relative to choosing, the reading of `imnl()`'s split equation,
and |C| = 1 is `imnl()` with the logit link exactly (theta_ic = exp(-eta_ic)).
A captive alternative not in a case's choice set is not a regime for that
case (the dogit sums over the choice set; `imnl()`'s convention for k).
Derivatives (derived here; checked against numDeriv and, at |C| = 1, against
the `imnl()` kernel): with the posteriors tau_0 = pi_0 P_y / P*_y and
rho_c = 1{y = c} pi_c / P*_y,

    dl/dV_j          = tau_0 (1{j = y} - P_j)
    dl/deta_c        = pi_c - rho_c
    d2l/dV_j dV_m    = tau_0 (1 - tau_0)(1{m = y} - P_m)(1{j = y} - P_j) - tau_0 P_j (1{j = m} - P_m)
    d2l/dV_j deta_c  = tau_0 rho_c (1{j = y} - P_j)
    d2l/deta_c deta_d = pi_c (pi_d - 1{c = d}) - rho_c (rho_d - 1{c = d})

computed on the log scale (log pi by log-sum-exp over the case's regimes,
log P*_y as a two-term log-sum-exp).
* Kernel `src/captive_ll.cpp`; `imnl()` with one inflated alternative keeps
  the existing kernel, so its results do not change. Parameters `infl_<alt>:<term>`
  (one block per captive alternative, the same z). Starts as for `imnl()`, with
  a binary logit of 1{y != c} on z for each c and common shifts of the
  intercepts; challengers as for `imnl()`. Boundary: the existing rules plus
  a captive share below 1e-6 for every case where c is available (that
  alternative's captivity has collapsed; the fit is then the model without
  it). `predict(type = "regime")` returns the n x (1 + |C|) regime
  probabilities, `type = "posterior"` their posteriors given the choice;
  `rmnl()` takes several `inflate` with a matrix `delta`. Comparisons: MNL
  against the captivity model, and the single- against the multi-captive
  model, are boundary restrictions, calibrated by `boot_lr_test()`.
* Captivity score test `captivity_test(mnl_fit, alternatives)`: the score
  test of the MNL against the dogit with constant captivity parameters at
  theta = 0, which is a boundary point in theta but a regular one (the
  captivity parameterization makes the MNL a finite point, unlike delta at
  infinity in `imnl()`). Per-case score s_ij = a_ij (1{y_i = j} / P_ij - 1);
  information (expected, at the MNL estimate) I_tt[j, m] = sum_i w_i a_ij a_im
  (1{j = m} / P_ij - 1), I_tb[j, .] = sum_i w_i a_ij (D_ij - Dbar_i)', I_bb the
  MNL information, with D_ij = dV_ij / dbeta and Dbar_i = sum_l P_il D_il;
  efficient information I_tt - I_tb I_bb^{-1} I_bt. Reported: a one-sided z
  for each alternative (captivity cannot be negative) and the joint LM
  chi-square. The weights 1 / P_ij make the statistic heavy-tailed when some
  P_ij are small, so the size Monte Carlo decides whether the normal
  reference is usable; `boot_pvalue()` calibrates it in any case.
  Probe (2026-10-19, before the Monte Carlo; scratchpad): the score equals
  the numerical derivative of the dogit log-likelihood; under the logit
  (J = 3, n = 1500, 2000 draws) the z of an alternative whose fitted
  probabilities reach 3e-4 has skewness 17, IQR-based sd 0.68, 99th
  percentile 3.07, while its one-sided rejection at 0.05 is 0.042 (0.053 and
  0.064 for the other two). So a second form, `information = "opg"` (the
  outer product of the per-case efficient scores, summed by cluster when the
  fit has one), is offered, and the rule for the default is fixed now: the
  expected form stays the default unless its one-sided size at 0.05 leaves
  0.05 +/- 2 Monte Carlo standard errors in some design of
  captivity_test_mc.R while the opg form's does not. The documentation
  reports the size at 0.01 as well.
* Oracles and verification: the kernel at |C| = 1 against `imnl()`'s to 1e-12
  (lli, gradient, Hessian); numDeriv; the dogit closed form with an
  intercept-only split; recovery MC (`captive_mc.R`: J = 4, two captive
  alternatives, n = 2000 and 5000, with and without an exclusion
  restriction; criteria as recovery_coverage_mc, and the shares of boundary
  fits reported); size and power MC of `captivity_test()` (`captivity_test_mc.R`:
  size under the MNL at n = 500, 2000; power against a single-captive imnl).

**C. Consideration sets (latent): gated.** Manski's (1977) two-stage model
with independent consideration (the independent-availability model; Swait
and Ben-Akiva 1987a; Ben-Akiva and Boccara 1995): alternative j in A_i is
considered with probability q_ij = Lambda(z_i' alpha_j), independently, and
the case chooses by the MNL among the considered set C:

    P*(y_i = j) = sum_{C subset A_i, j in C} Q_i(C) exp(V_ij) / sum_{m in C} exp(V_im)
    Q_i(C) = prod_{m in C} q_im prod_{m in A_i \ C} (1 - q_im) / (1 - prod_{m in A_i} (1 - q_im))

or with alternatives the user declares always considered (abstention, a
status quo), which removes the conditioning on a non-empty set. The MNL is
the boundary q = 1. Political-science uses: Paap et al. (2005), Oscarsson and
Rosema (2019) and that Electoral Studies symposium. The identification
results in the economics literature (Abaluck and Adams-Prassl 2021;
Barseghyan et al. 2021; the survey of Crawford, Griffith and Iaria 2021) lean
on alternative-specific characteristics; with case-specific covariates only
(the usual political-science design) the consideration and utility equations
are separated only by functional form and exclusion restrictions. Hence two
gates before any package code.
Corrections after reading the papers (2026-10-19; ledger): Paap et al. use
stated consideration sets in estimation (a multivariate probit for
consideration, a multinomial probit with a diagonal covariance for the choice,
the same respondent covariates in both), but political science also
estimates latent consideration sets: Oscarsson and Rosema (2019, p. 260)
name Steenbergen and Hangartner (2008), De Vries et al. (2009), and Moral
and Zhirnov (2018, AJPS) for that strand (not read here), while recommending
direct measures. Of the economics papers, Abaluck and Adams-Prassl identify
from asymmetric responses to an attribute that varies across alternatives
(prices), Barseghyan et al. give partial identification (point identification
needs a restricted choice-set distribution, variation in x, and exclusion
restrictions), and Crawford, Griffith and Iaria's sufficient sets need panel
data or groups facing the same choice set; none of them analyses a design with
respondent covariates only, so the sentence above is this design's inference,
not their statement. Gate C1 covered case-specific covariates only; designs
with a voter-party covariate (an alternative-specific w_ij such as
ideological distance), where the economics results apply in spirit, were not
tested.
* Gate C1 (identification, no fitting): the exact expected information of
  the model (summing over the J outcomes per case) at designs J = 3, 4, 5,
  n = 2000, slope sd 0.5 and 1, consideration probabilities averaging 0.6-0.9,
  (a) z = x, (b) one variable excluded from each side, (c) one
  always-considered alternative; report the smallest eigenvalue of the
  normalized information and the standard error of every parameter relative
  to its true value. Pass: in design (b) or (c), every consideration
  parameter's standard error below half its true magnitude at n = 2000.
* Gate C2 (recovery, a prototype outside the package): 100 replications at
  the passing design; pass: median bias of every parameter within 10 percent
  (or 0.05 absolute), Wald coverage 0.90-0.98, at most 5 percent of fits at a
  boundary. If C1 fails in every design, the model is not implemented, and
  the documentation shows instead how a stated consideration set (a survey
  question) enters `mnl()` as the choice set; the author decides.
* If both pass: kernel `src/consider_ll.cpp` (enumeration over the 2^(|A|-1)
  sets containing y, |A| at most 10; the gradient as the posterior-weighted
  complete-data score; Hessian by differences of the gradient), function
  `consider()`, simulator, recovery MC.
* Gate C1' (fixed 2026-10-19, after C1 failed, at the author's go-ahead;
  before any computation): the same exact expected information, now on
  designs with a voter-party covariate, calibrated on real data. Datasets:
  those of the replication stress test (verification/stress/) that have an
  alternative-specific covariate varying across voters (a voter-party
  distance): Alvarez and Nagler (1995), Glasgow (2001), Arzheimer and Evans
  (2012), and Steiner and Hillen (2021), whichever have their data. In each,
  the cases, the covariates and the choice sets are the paper's estimation
  sample; every alternative-specific covariate is standardized to sd 1 over
  the available rows. Utility: the paper's specification, with the true
  values at the conditional / multinomial logit estimates on that sample.
  Consideration: logit q_ij = a_j + alpha' v_ij with an intercept for every
  alternative not always considered (set so that the mean of q over the
  cases runs from 0.6 to 0.9 across the alternatives, in level order) and
  common slopes alpha of magnitude 0.5 and of 1 per sd, negative (more
  distance, less consideration). Designs: (d) v = the utility's
  alternative-specific covariates (the same distances in both equations;
  the asymmetry identification of Abaluck and Adams-Prassl 2021); (e) v = an
  alternative-specific covariate excluded from the utility (a second
  distance in the data, moved from the utility to consideration; skipped
  where the data have only one); (f) (e) with the paper's base alternative
  always considered. Information over the paper's cases, rescaled to
  n = 2000 (the information is a sum over cases) and also reported at the
  paper's n. Pass (as C1): in some design of some dataset, every
  consideration parameter's standard error below half its true magnitude at
  n = 2000. A pass sends that design to Gate C2; the documentation would
  then state the designs that pass and those that do not. A failure
  everywhere closes the latent model.
* Gate C2 on the C1' passes (fixed 2026-10-19 after C1' passed in 29 of 30
  cells, before any C2 run): three designs, each the weakest passing cell of
  its kind: Alvarez-Nagler (d) alpha -0.5; Glasgow (d) alpha -0.5; Arzheimer
  (e) with the driving distance to the candidate in consideration only
  (drdistp), alpha -0.5. 100 replications each at n = 2000: cases drawn with
  replacement from the paper's cases, outcomes simulated from the model at
  the C1' truth. Fit: maximum likelihood (BFGS, numerical gradient) from two
  starts, the truth and the conditional logit with every q at 0.8 (alpha = 0),
  the higher kept; a fit is at the boundary when a consideration intercept
  exceeds 10 or a slope 10 in magnitude; standard errors from the inverse
  numerical Hessian. Criteria as C2 above, for every parameter; the Monte
  Carlo error of each criterion is reported beside it.
* Gate C2 result (2026-10-20, ledger): passes in the design with a
  consideration-only covariate (Arzheimer (e)), fails in both designs where
  the same distances drive both equations (Alvarez-Nagler (d): 14 percent
  boundary fits, consideration coverage 0.80-0.91; Glasgow (d): 12 percent
  boundary fits). The author chose (2026-10-20) to implement the model
  restricted to designs with consideration-only covariates.
* Implementation (fixed 2026-10-20, before any code), within what passed:
  - `consider(choice ~ w | x | v, data, id, alt, ref = NULL, always = NULL,
    weights, offset, cluster, se, nboot, cores, start, maxit, reltol)`, long
    layout only (`| 1 |` for an empty case-specific part). Utility as in
    `mnl()`. Consideration: q_ij = Lambda(a_j + v_ij' alpha), an intercept
    a_j for every alternative not named in `always` (always considered,
    q = 1) and generic slopes alpha on the columns of v. Every variable of v
    must vary across the alternatives within cases (a case-specific
    consideration covariate is identified only through the functional form,
    Gate C1) and must not appear in the utility parts (the designs that
    shared them failed Gate C2); both are refused with the reason. The
    considered set is conditioned to be non-empty when no always-considered
    alternative is available. No `fe` (refused).
  - Parameters: the utility blocks of `mnl()`, then the structural block
    "consideration": `consider:<alternative>` (the intercepts, in level
    order) and `consider:<column>` (the slopes); scales 1 for the intercepts
    and the column standard deviations over the available rows for the slopes.
  - Kernel `src/consider_ll.cpp`: log P*_iy by enumeration over the subsets
    of A_i that contain y (lli, gradient, scores), all non-empty subsets for
    the probabilities; |A_i| at most 10 (larger choice sets refused). The
    gradient as in consider_gate2.R (checked there against numDeriv):
    omega(C) = Q(C) pi_y(C) / P*_y, dl/dV_m = sum_C omega 1{m in C}(1{m = y}
    - pi_m(C)), dl/deta_m = sum_C omega (1{m in C} - q_m) - q_m P0 / N (the
    last term only without an available always-considered alternative); the
    Hessian by differences of the gradient (the engine's default).
  - Starts: the conditional logit on the utility with every free q at 0.6,
    0.8 and 0.95 (intercepts at their logits) and zero slopes, each screened
    and polished (the engine's multistart).
  - Boundary: the separation rule of the fit, plus an alternative considered
    with q > 1 - 1e-6 by every case where it is available (the model without
    that alternative's consideration stage), flagged with a warning.
  - `predict(type = "prob")` returns P*, `type = "consideration"` the n x J
    matrix of q (1 for always-considered alternatives, NA where
    unavailable); `simulate()` draws from P*; simulator `rconsider()`.
  - Oracles (tests): with every intercept large (q near 1) the probabilities
    and log-likelihood equal the conditional logit's; the probabilities equal
    the independent R enumeration of Gate C1' (`pstar`) and sum to one; the
    gradient against numDeriv, with and without an always-considered
    alternative and with varying choice sets.
  - Recovery Monte Carlo `verification/consider_mc.R`, fixed now: designs
    calibrated on real data as in Gate C2, with the package's own starts (no
    start at the truth): (A) Arzheimer (e), drdistp in consideration,
    alpha -0.5; (B) Glasgow (e) with disdef moved to consideration (its
    weakest (e) cell in C1'), alpha -0.5; (C) Glasgow (f), the same with the
    base (Alliance) always considered. 200 replications each at n = 2000.
    Criteria as Gate C2 (median bias within 10 percent or 0.05, coverage
    0.90-0.98, at most 5 percent boundary fits), with their Monte Carlo
    errors; the time per fit reported. A design that fails is reported in
    the documentation as such, and the author decides on the scope.
* Recovery result (2026-10-20, ledger): (A) passes; (B) fails (consideration
  coverage 0.877-0.913, two utility biases; at n = 5000 one consideration
  interval still covers 0.88); (C) fails on one utility bias. The standard
  errors are the inverse observed information (checked against numDeriv);
  the shortfall is weak identification near q = 1. **The author dropped
  the latent model (2026-10-20)**: `consider()`, `rconsider()` and
  src/consider_ll.cpp are removed (their last state is commit 5d41a66), the
  recovery study is archived in verification/archive/consider_mc.*, and the
  documentation shows the stated consideration set entering `mnl()` as the
  choice set.

### Phase 4 (from 2026-11-05, author's decision): Firth's penalty (`firth = TRUE`)

Definitions fixed 2026-11-05, before any package code. Derivations: the
author's working note (outside the repository); every formula below was
checked numerically first (verification/probes/2026-11-05: firth_proto.R,
firth_proto_se.R, coxphf_window.R, with their outputs). Motivation:
separation (Cook, Niehaus and Zuhlke 2018 for the multinomial logit); the
package flags it (`boundary`, the IIA tests' separation warnings) but had no
remedy.

**Scope.** `mnl()` and `rologit()` gain `firth = FALSE`. The other
estimators get no such argument: in `imnl()` and `remnl()` the restricted
model sits where the information is singular (phi -> 1, sigma -> 0), so the
Jeffreys penalty tends to -Inf there and would bias every comparison with
the multinomial logit towards the larger model; `nlogit()`, `hev()` and
`mnprobit()` are not canonical (the penalty is not Firth's adjusted score
there, and the GHK likelihood is not smooth). `firth = TRUE` with `fe` is
refused (and so the jackknife).

**The penalized log-likelihood.** l*(theta) = l(theta) + (1/2) log|I(theta)|,
maximized; with D_R the design rows of the alternatives in R,
W_R = diag(P_R) - P_R P_R':
* `mnl()` (wide and long, varying choice sets, offsets, frequency weights):
  I = sum_i w_i D_i' W_i D_i, the Fisher information (observed = expected).
  The model is a canonical exponential family, so this is Firth's (1993)
  adjusted score exactly (Bull, Mak and Greenwood 2002; Kosmidis and Firth
  2011).
* `rologit()`: I = sum_i w_i sum_G D_{R_G}' Omega_G D_{R_G} over every group
  G of a ranking (an untied position is a group of one) with an alternative
  below it, R_G the group and the alternatives below it, d_G its size:
  Omega_G = W_{R_G} for d_G = 1; for a tie group Omega_G =
  sum_{S subset R_G, |S| = d_G} pi_S g_S g_S', where pi_S is the exact tie
  probability (Allison and Christakis 1994, eq. 8) of S as the top-d set of
  R_G and g_S = d log pi_S / dV_{R_G}. Untied rankings give Heinze and
  Schemper's (2001) penalty for the Cox model stratified by case, each stage's
  information given its remaining set (Heinze and Dunkler 2008, p. 6459); the
  tie term is the expected information of the tie group's outcome given R_G
  and d_G, positive semidefinite and independent of which set occurred (the
  observed information of the tie term is not used: it depends on the set and
  need not be positive semidefinite). Only the stages the fit uses (`depth`)
  enter; under `rank_test()` the split model's information.
* Gradient (closed forms, derived and checked against numDeriv): per untied
  stage (1/2) w_i D' W a with a = h - 2 H P, H = D I^{-1} D', h = diag(H);
  per tie group (1/2) w_i D_R' sum_S pi_S [(g_S' H g_S) g_S + 2 G_S H g_S],
  G_S = d^2 log pi_S / dV dV'. Two passes per evaluation: I at the current
  point, its inverse, then the per-case terms (C++).
* Tie groups: pi_S, g_S and G_S by the recursion f(U) = sum_{j in U}
  f(U \ j) lambda_j / (Lambda_R - lambda(U \ j)) over the subsets of R_G of
  size at most d_G (positive terms only; not inclusion-exclusion, which
  cancels). A group with more than 10^5 such subsets is refused under
  `firth = TRUE`, naming the case.
* Optimizer: BFGS with the analytic penalized gradient from the usual start,
  then the Newton polish with the Hessian by central differences of the
  analytic penalized gradient (no Newton-first: concavity of l* is not known);
  convergence by `.pc_gtol()` on the penalized objective; the saddle check on
  its Hessian. Amended 2026-11-05 during the first tests, before any Monte
  Carlo: on a completely separated ranking design BFGS's first steps followed
  the separated direction far out and stalled on a flat stretch (34.7 where
  the maximum, confirmed by coxphf and by the gradient, is at 4.6), so a
  penalized fit now starts with Fisher scoring (the damped Newton iteration
  with the information I as its curvature, up to 200 steps, the polish's step
  cap and step-halving; as Bull, Mak and Greenwood 2002, logistf and coxphf),
  and BFGS and the polish continue from its endpoint only when it stops
  short. The restricted and profile maxima use the same scoring.

**What the fit reports.** `coefficients` = argmax l*; `loglik` = l at those
estimates (the ordinary log-likelihood, below the maximum-likelihood value);
`loglik_penalized` = l*; `firth = TRUE`. Covariances: `"analytic"` =
I(theta_hat)^{-1}, the model's information at the penalized estimates (as
brglm2, logistf and coxphf; Kosmidis and Firth 2011, p. 758); `"opg"`,
`"robust"`, `"cluster"` from the per-case adjusted scores (the score plus the
case's share of the penalty gradient; they sum to zero at the estimate) with
that bread; `"bootstrap"` refits every replicate with the penalty. The
separation diagnostic runs on the maximum-likelihood fit of the same design
(computed alongside; both log-likelihoods are concave, the tie term by
Prekopa's theorem): the names of its diverging coefficients are stored
(`separation`) and printed, saying that their penalized estimates are
determined largely by the penalty. Amended 2026-11-05 in the first tests:
a coefficient is reported as separated when the rule flags its
maximum-likelihood estimate and that estimate lies more than 3 penalized
standard errors from the penalized one (on Cook, Niehaus and Zuhlke's
Model 3 the rule also flagged the constant of an uncentered year of birth,
171.8 against 165.5 with a standard error of 77, which is not separated).

**Inference with a penalized fit.**
* `confint()`: `type = "profile"` (the existing `type` argument gains the
  value; written `method` in the first version of this section), the default
  for a penalized fit: the values c of theta_r with
  2 [l*(theta_hat) - max_{theta_r = c} l*(theta)] <= qchisq(level, 1), the
  full model's penalty in both terms (Heinze and Schemper 2002); found by
  root finding outward from the Wald limits, NA with a warning when no root
  lies within 50 Wald standard errors. Maximum-likelihood fits keep Wald
  intervals.
* `summary()`: Wald z from I(theta_hat)^{-1}, with a line naming the penalty
  and pointing to profile intervals. Rule fixed now: if the Monte Carlo below
  finds Wald coverage outside [0.92, 0.98] in a design where profile coverage
  is inside, `summary()` reports penalized likelihood-ratio p-values instead.
* Penalized likelihood ratio (Heinze and Schemper 2001; `coxphftest`):
  2 [l*(theta_hat) - max_{theta in Theta_0} l*(theta)], the restricted
  maximum taken over the full model's penalized log-likelihood (two
  separately penalized fits do not give it: their penalties differ). Used by
  the LR forms of `term_test()` (the term's coefficients zero) and
  `combine_test()` (Cramer and Ridder's restriction: all coefficients but the
  constants equal for the two alternatives), by `lr_test()` for two
  penalized fits nested by name (the restricted parameters zero; a penalized
  and an unpenalized fit are refused), and by `rank_test()` (early equal to
  late); the Wald forms use I(theta_hat)^{-1}.
* Refused for penalized fits, with the reason: `hausman_mcfadden()`,
  `small_hsiao()`, `fp_tests()`, `iia_tests()`, `nested_lm_test()`,
  `captivity_test()`, `vuong()`, `boot_lr_test()`.
* Unchanged: `predict()`, `fitted()`, `residuals()`, `first_difference()`,
  `rrr()`, `classification()`, `hl_test()`, `simulate()`, `update()`,
  `coef()`, `vcov()`, `nobs()`; `boot_pvalue()` refits every replicate with
  the penalty. The documentation says that the penalty shrinks fitted
  probabilities towards equality (Puhr et al. 2017).
* `fit_stats()`: the log-likelihood at the penalized estimates and the
  penalized log-likelihood; AIC and BIC from the former, labelled as not
  maximum-likelihood information criteria. `logLik()` the former, with
  attribute `penalized`; `glance()` adds `penalized = TRUE`;
  `compare_models()` marks penalized fits.

**Oracles (tests; the reference packages in Suggests, each test
`skip_if_not_installed`).**
1. `mnl()`, case-specific covariates, with and without weights:
   `brglm2::brmultinom()` (estimates to its tolerance; standard errors
   I^{-1}).
2. J = 2: `logistf::logistf()` on the differenced data (case-specific and
   alternative-specific covariates): estimates, profile intervals, penalized
   LR p-values.
3. `rologit()` without ties (complete, partial, `depth = 1`, varying choice
   sets) and `mnl()` in the long layout: `coxphf::coxphf()` with each case on
   its own disjoint time window, Surv(start, stop, event) (reproduces a
   stratified Cox fit exactly; verified 2026-11-05): estimates, variance,
   profile intervals (`pl = TRUE`), `coxphftest()`; iteration settings
   raised for separated designs.
4. numDeriv: the penalized gradient (mnl wide and long with weights, offsets
   and varying choice sets; rologit with partial rankings and ties of two and
   three).
5. Tie information: the C++ Omega against an independent R enumeration over
   subsets and orders; Omega = W at d = 1; the subset recursion against the
   sum over orders.
6. Invariance and identities: fitted probabilities and l* identical across
   `ref`; weight duplication; `firth = FALSE` leaves every existing fit
   bit-identical (the existing tests unchanged; the stress battery's mnl and
   rologit fits rerun and compared bitwise).
7. Separation: finite estimates on a constructed separated MNL and ranking
   design, and on Cook, Niehaus and Zuhlke's data (`defection`, added at the
   author's decision of 2026-11-05 after a licence check, CC0 1.0): their
   penalized Models 3 to 5 (brglm2, from their replication log) reproduce.

**Monte Carlo `verification/firth_mc.R`** (fixed now; run on the final
frozen build with the rest), 1000 replications per cell:
* F1, multinomial logit, J = 3, wide, n = 50 and 200 (after Bull, Mak and
  Greenwood 2002, Table 6, models 9 and 11): x1 = 1{z1 > 0}, x2 = z2,
  (z1, z2) bivariate normal with correlation 0.752 (cor(x1, x2) = 0.6), and
  a null covariate x3 ~ N(0, 1); coefficients (intercept, x1, x2, x3) for
  alternatives 2 and 3: model 9 (-1.4, 1, 1, 0) and (-1.4, 1, -1, 0);
  model 11 (-1.4, 2, 1, 0) and (0, 2, 1, 0).
* F2, rank-ordered logit, J = 4, n = 100 and 300 cases, complete rankings
  and the same data cut at depth 1 and 2: constants (0.5, 0.25, -0.25, base
  0), a normal attribute (gamma 0.8), a binary attribute with prevalence 0.1
  (gamma 1.5) and a null normal attribute (gamma 0).
* F3, F2 at n = 100 with ties: in each ranking, independently of the
  utilities, positions 2 and 3 are merged into a tie with probability 0.3
  and positions 1 to 3 with probability 0.1 (a coarsening that does not
  depend on the order, the exact likelihood's assumption).
* Reported per cell: the share of samples whose maximum-likelihood fit is
  flagged separated; for the penalized fit, convergence, mean bias, RMSE and
  median absolute deviation over all samples; maximum likelihood's over the
  samples it converges in; coverage of the 95 percent profile and Wald
  intervals for every coefficient; the size at 0.05 of the penalized LR test
  (`term_test(test = "lr")`) of the null coefficient; the time per fit.
* Criteria: (a) every penalized fit converges; (b) profile coverage within
  [0.93, 0.97] for every coefficient at the larger n of F1 and F2 and in F3;
  (c) the penalized LR test's size at 0.05 within [0.035, 0.065] in every
  cell; (d) at the larger n, every penalized estimate's mean bias within 3
  Monte Carlo standard errors or 5 percent of its true value. A cell that
  fails is reported in the documentation as such, and the author decides on
  the scope. The `summary()` rule above is applied from (b)'s Wald column.

**Documentation.** `?mnl` and `?rologit` sections on penalized estimation;
NEWS; a vignette section with a separated example; the article (a paragraph
and the software table). Datasets from published separated analyses only at
the author's decision, after a licence check.

### Phase 4b (2026-11-05 evening, author's decision): an exact separation check, two starts for penalized fits, and the rank-ordered refusal

Written before the code. Motivation (research-note checks of the same day,
`Research Notes/rologit_firth/checks/`, outside the repository): on 240
simulated ranking data sets the engine's separation rule (`.pc_diverge`:
|theta scale| > 15, or a 0/1 column's |theta| > 15, with Wald |z| < 1.96)
missed one quasi-complete separation (a binary attribute's coefficient
stopped at 11.4 with a standard error of 2.1e4, fit reported converged) and
flagged one overlapping data set whose estimates are large and finite (23 to
30); a linear program decided both correctly and agreed with
detectseparation on 549 binary data sets. The penalized likelihood can have
two local maxima (binary logit examples, research notes section 9). And
`rologit()` refused data in which an alternative is never ranked above
another, with or without `firth = TRUE` and with or without alternative
constants, although that is separation (or nothing at all), not
non-identification.

**Separation check.**
* Scope: `mnl()` and `rologit()` fits without unit fixed effects and without
  inflation (the conditional, multinomial and rank-ordered logit), fitted
  with diagnostics (`check = TRUE`), and the IIA tests' refits of such
  designs (amended below); not `rank_test()`'s split designs. Other families
  and fixed effects keep the rule `.pc_diverge` (their scale and nesting
  parameters and unit constants are outside this criterion).
* Criterion: a finite maximum-likelihood estimate exists iff the structure
  vectors overlap (Albert and Anderson 1984; Speckman, Lee and Sun 2009 for
  choice models; Sablica, Hornik and Rusch 2026, Theorem 1, for rankings,
  partial and tied ones included: a tie group is the event that each tied
  alternative beats each alternative below it).
* Structure vectors, in theta coordinates (D_ij: the case-specific block of
  non-base alternative m equals X_i when j is that alternative, zero
  otherwise; then the generic part W_ij), for every case with positive
  weight: `mnl()`, D_iy - D_il for every other available l; `rologit()`,
  D_ia - D_ib for a in a rank group and b in the next one (ties as groups,
  the available unranked alternatives last; transitivity leaves the cone
  unchanged). Offsets and positive weights do not enter. Zero rows and
  duplicates dropped; columns scaled by their largest absolute value.
* Linear program (lpSolveAPI): maximize sum_i tau_i subject to
  s_i'v - tau_i >= 0, 0 <= tau_i <= 1, v free. The cone of separating
  directions is closed under addition and positive scaling, so one v reaches
  tau_i = 1 on every comparison that some separating direction makes
  strictly positive: the optimum counts those comparisons. Overlap iff the
  optimum is below 0.5; separated comparisons are those with tau_i > 0.5
  (reported by case); complete separation iff all comparisons are separated;
  the reported direction is v normalized (it separates every separable
  comparison, so it lies in the relative interior of the cone of separating
  directions), and the separated terms are its coordinates with |v_k| above
  1e-8 times the largest.
  Amended 2026-11-05 at the first full test run and timing, before any Monte
  Carlo: (a) that program has a variable per comparison, and lp_solve took
  2.1 s at 6,000 comparisons and 88 s at 30,000 (multinomial logit, J = 4,
  K = 18-21; `lp_bench.R`, kept with the research-note checks), and it
  returned "infeasible" (status 2, impossible for a program that v = 0,
  tau = 0 satisfies) on five fits of the test suite (the cabinets
  conditional logit, TravelMode, a design with covariates scaled by 1e-6,
  a long-layout bootstrap, a rank_test() fit). The program is now in the K
  coefficients: maximize the sum of s_i'v over the comparisons not yet found
  separable, subject to s_i'v >= 0 for all and -1 <= v_k <= 1; the
  comparisons it makes positive (above 1e-9) are separable, and it is solved
  again over the rest until it adds none. A separating direction positive on
  a remaining comparison would make the objective positive, so every
  separable comparison is found, and the sum of the solutions, positive on
  all of them, is the reported direction (0.05 s and 0.1-0.7 s on the same
  designs, the same separated comparisons; no solver failure in the suite).
  (b) The separated terms are the coefficients not orthogonal to the span of
  the cone of separating directions, which is the null space of the
  comparisons no separating direction favours (singular values at most 1e-9
  of the largest): every coefficient that moves along some separating
  direction, which the coordinates of one solution need not show (a vertex
  can have zero coordinates). (c) Duplicate comparisons are found by
  sorting, exactly. (d) Applied to the stress battery's 63 in-scope designs
  (`verification/separation_stress_check.R`), the program took up to 27 s
  (101,249 cases, 202,498 comparisons) and 59 s (777,688 comparisons),
  against at most 6 s for the fits themselves, so overlap is first sought by
  a certificate: by Stiemke's lemma the comparisons overlap iff some u > 0
  has S'u = 0, and at any theta the gradient of the log-likelihood is S'u
  for positive weights u from the fitted probabilities (in `mnl()`
  w_i p_il on the choice over l; in rankings without ties the stage
  probabilities summed along the chain of consecutive comparisons; not with
  ties, whose gradient is not such a combination). At a maximum-likelihood
  estimate S'u is near zero; the least-norm correction that makes it vanish
  (twice) certifies overlap when every corrected weight stays above 1e-12 of
  the largest and the remaining S'u is at most 1e-9 of sum |S|'u. Otherwise
  the program decides. The fits pass their estimates (a penalized fit the
  maximum-likelihood ones of its second start); `separation()` the fit's.
  (e) A re-solve after the objective changes starts from the previous basis;
  where lp_solve fails from there (status 3 on Kardasheva 2013, where a
  fresh model solves under six settings) the program is rebuilt once.
* Use: (i) in-scope maximum-likelihood fits: the diverging coefficients are
  the separated terms (replacing `.pc_diverge` there; an overlapping data set
  has none, however large its estimates); the warning names the separation
  type and terms and points to `firth = TRUE` and `separation()`; (ii)
  in-scope penalized fits: `separation` holds the separated terms (replacing
  the maximum-likelihood rule with the 3-standard-error distance); (iii) an
  exported `separation(object)` returns the full result (overlap, type,
  terms, direction, separated cases, counts) for any in-scope fit. Errors in
  the solver fall back to the rule with a warning.

**Two starts for a penalized fit.** Without a user `start`, a penalized fit
runs from zero (the default) and from the endpoint of the maximum-likelihood
fit of the same design (Newton first; when the penalized objective is not
finite there, as along a separating direction far out, the point is halved
towards zero until it is, at most 60 times), and keeps the higher maximum. If
both runs converge to points that differ (largest |theta scale| difference
above 1e-3) with penalized log-likelihoods that differ by more than
1e-6 (1 + |l*|), a warning reports both maxima. Reason: in the research-note
examples the zero start finds the lower maximum in one design (separated)
and the maximum-likelihood start the lower one in the other (not separated).
`rank_test()`'s penalized fit with separate early and late coefficients
(started at the fit's estimates) runs from the split design's
maximum-likelihood estimate as well, under the same rule; the restricted
maxima already run from two starts (`.pc_firth_restricted`).

**`rologit()`'s refusal.** The stop for an alternative never ranked above
another available alternative is removed. With alternative-specific
coefficients it is separation (the alternative's constant and case-specific
coefficients run to minus infinity), which the check reports for maximum
likelihood and the penalty handles; without them it is no restriction. The
identification check (`.pc_identification_check`) and the singular
information warning remain. `mnl()`'s check (an alternative chosen only where
no other is available) is unchanged.

Amended the same evening, at the author's decisions after the first tests
(before any Monte Carlo):
* `mnl()` without unit fixed effects or inflation no longer refuses an
  alternative never chosen (or chosen only where no other is available, or
  only by cases of weight zero): the same separation as in rankings, reported
  by the check (its constant and slopes; the other coefficients tend to the
  fit without that alternative) and finite under `firth = TRUE`. `imnl()`,
  `nlogit()`, `hev()`, `mnprobit()`, `remnl()` and fixed-effects fits keep
  the refusal: there the alternative's own parameters (split equation,
  nests, scales, variances, unit constants) would not be identified at all.
* Bootstrap replicates (`se = "bootstrap"`, the bootstrap tests' simulated
  responses): a resample in which an alternative with coefficients of its
  own is never chosen (in a ranking, never ranked above another available
  alternative, when the model has alternative-specific coefficients) is a
  failed maximum-likelihood replicate, as before; a penalized replicate is
  kept (finite).
* The IIA tests judge their half-sample and restricted refits of `mnl()`
  by the same exact check (`.pc_diverging`), so one criterion holds wherever
  an in-scope design is fitted.
* Large but finite estimates raise no flag; `?mnl` says that their Wald
  statistics understate the evidence (Hauck and Donner 1977) and points to
  the likelihood-ratio tests and `firth = TRUE`.

Amended 2026-11-05 (night), at the first replication for the research note
(Poast 2007, outside the repository), before any Monte Carlo: the penalty's
log-determinant and the inverse in its gradient are those of the
equilibrated information I_s = I / (d d'), d = sqrt(diag(I)):
log|I| = log|I_s| + 2 sum log d and I^{-1} = I_s^{-1} / (d d'). The same
quantities; in theta units I carries the covariates' scales (Poast's bid
covariates run from 0.3 to 1.7e8: reciprocal condition number 1e-18, against
0.02 for I_s), solve() refused it, and the generalized inverse in its place
truncated it, so the penalty's gradient was wrong (by up to 1.7 in the
standardized coordinates) and every penalized fit of the archive's models
stopped unconverged. The covariance was not affected: it inverts I in the
standardized coordinates. Test: a covariate multiplied by 1e9 in `rologit()`
and `mnl()` leaves the penalized fit equivariant.

**Tests.** The check against detectseparation (binary data through `mnl()`
and through `rologit()` with two alternatives); the two data sets the rule
got wrong, kept as fixed data; the log-likelihood nondecreasing along the
reported direction; never-ranked-above data (maximum likelihood warns, the
penalized fit is finite, the fit without constants is unaffected); the two
two-maximum designs (penalized maxima -1.922670 and -0.157299 found, with
the warning); the existing tests updated where the separation report changes.
The time of the check on the stress battery's in-scope data sets is measured
and reported (not a test).

### Phase 4c (2026-11-06, author's decision): likelihood-ratio statistics from suprema, the parametric bootstrap for rankings, and priors

Written before the code. Motivation: the research note's analysis of Poast's
(2007) unreported models (`Research Notes/rologit_firth/replications/`,
outside the repository).
1. `term_test(type = "lr")` on a completely separated maximum-likelihood
   fit reported LR = 2.50 for a term whose removal leaves the data
   completely separated. Both suprema are 0; the restricted refit had
   stopped at -1.25.
2. With about one to three choices per coefficient, the chi-square
   reference of the likelihood ratio was far off. Bootstrap medians were 13
   to 30 on 10 or 11 degrees of freedom; one bootstrap p-value of 0.145 had
   a chi-square p-value of 0.0003. `boot_lr_test()` refused rank-ordered
   fits.
3. Under separation the penalty decides the size of the separated
   estimates (Rainey 2016). Jeffreys' prior is one default; it is not
   independent across coefficients and can move an estimate against the
   direction of the separation (Beiser-McGrath 2022). Rainey recommends
   showing the estimates under several priors.

**Suprema.** The maximum-likelihood forms of `term_test()`, `lr_test()`,
`combine_test()` and `boot_lr_test()` compare suprema of the log-likelihood,
not the values where the optimizer stopped. For an in-scope fit (the scope of
the separation check):
* **Overlap:** the fit's maximum (unchanged).
* **Complete separation:** 0. Along the reported direction v every reported
  comparison is favoured, so every choice probability, and a tie group's
  exact probability, tends to one.
* **Quasi-complete separation.** Along v the separated comparisons' terms
  vanish, while the others do not change, since v is orthogonal to every
  comparison it does not favour. The supremum is approached as follows:
  1. From the fit's point, move along v until every separated comparison's
     utility difference is at least 40.
  2. Re-maximize from there by the engine's optimizer.
  3. Repeat the move and the re-maximization until the log-likelihood rises
     by less than 1e-10, at most five times.

  The value is below the supremum by at most the number of separated
  comparisons times e^-40 plus the optimizer's tolerance.
* **Out-of-scope designs** (unit fixed effects, inflation, the other
  families, `rank_test()`'s split designs): the fit's value, as before.
  There the rule `.pc_diverge` judges divergence and no direction is
  available.

A result whose statistic used a supremum in place of a fit's value says so in
its method line.

**The parametric bootstrap for rank-ordered fits.** `boot_lr_test()` accepts
two nested `rologit()` maximum-likelihood fits of the same rankings.
* **Draws.** Each replicate draws, for every case, the restricted fit's
  utilities plus independent standard Gumbel errors over the available
  alternatives, and orders them. Each case then gets its observed pattern:
  the same rank groups, with the same sizes in the same order (a tie stays a
  tie at its positions), and the alternatives below them unranked. This
  treats the coarsening (truncation, ties) as independent of the utilities.
* **Statistic.** 2 (sup_full - sup_restricted), with the suprema above. The
  full refit also starts from the restricted replicate's estimates embedded
  in it, so it cannot end below the restricted one.
* **Separated replicates are kept.** The statistic is defined through the
  suprema, and in small rankings separation is common under the null: in
  Poast's subsamples, 121 to 199 of 199 replicates. They are counted in the
  method line.
* **Penalized fits** are still refused. Their bootstrap needs penalized
  refits of both models on every replicate and is not part of this phase.

**Priors (`prior =` in `mnl()` and `rologit()`).**
* **Constructor:**
  `coef_prior(family = c("normal", "cauchy", "t"), scale = 2.5, df = NULL,
  intercept_scale = 10, autoscale = TRUE)`, location 0. The defaults are
  Gelman et al.'s (2008) Cauchy(0, 2.5) for slopes and scale 10 for
  constants.
* **Scope and penalty.** It is not allowed with `firth = TRUE`, and has the
  same scope as `firth`: `mnl()` and `rologit()`, no unit fixed effects.
  The objective is the log-likelihood plus the log prior density of each
  coefficient, up to constants:
  - normal: -(theta/s)^2/2;
  - t(df): -(df + 1)/2 log(1 + (theta/s)^2/df);
  - Cauchy: t(1).
* **Scales.** In theta units the scale is s_k = scale_k / c_k.
  - With `autoscale`, c_k is twice the standard deviation of coefficient
    k's column for a column with more than two values. For a two-valued
    column it is the difference of the two values. Both are unweighted,
    over the estimation rows: the cases for a case-specific covariate, the
    available alternative rows for a generic one. Otherwise c_k = 1.
  - The constants (the intercept of the case-specific part, which gives the
    alternative-specific constants) take `intercept_scale`, with c_k = 1.
  - `scale` may be a named vector by coefficient name. A name set there
    overrides the default scale for that coefficient, so a separated
    coefficient can carry its own prior (Rainey's partial prior). An
    unknown name is an error.
* **Estimation:** as a penalized fit.
  - Two starts: zero, and the maximum-likelihood endpoint.
  - Fisher scoring with I(theta) + max(H_k, 0) as the curvature, where H_k
    is the negative second derivative of coefficient k's log prior.
* **Output.**
  - `loglik` is the log-likelihood at the estimates; `loglik_penalized`
    is the maximized objective.
  - The covariance is (I + H)^{-1} at the estimates: the curvature of the
    log posterior at its mode, as `arm::bayesglm` reports. Firth's fits
    report I^{-1}, since Jeffreys' prior is a bias correction, not
    information.
  - Where I + H is not positive definite (a t prior far out on a
    separated coefficient), the covariance is (I + max(H, 0))^{-1} and a
    warning points to profile intervals.
* **Inference.** As for Firth's fits, with the prior's penalty:
  - profile penalized-likelihood intervals (`confint()`, the default);
  - penalized likelihood ratios in `term_test()`, `lr_test()` (two fits
    with the same prior specification) and `rank_test()`.
  - The functions that refuse a penalized fit refuse a prior fit the same
    way.
* **Separation.** It is reported in `separation`, as for Firth's fits.

**Tests.**
* **Suprema:**
  - a completely separated ranking in which `term_test()` gives 0 for a
    term whose removal leaves the data completely separated (the Poast
    pattern), where the stopping values gave a positive statistic;
  - a quasi-complete multinomial design whose supremum equals the
    maximum-likelihood fit of the data with the separated cases removed by
    hand;
  - `lr_test()` between separated fits.
* **Bootstrap:**
  - `boot_lr_test()` on two `rologit()` fits returns p-values in (0, 1];
  - the drawn rankings keep each case's groups;
  - a seed reproduces the result;
  - penalized fits are refused.
* **Priors:**
  - the normal and t modes equal a direct maximization of the
    log-likelihood plus log prior written in the test, for `mnl()` and for
    a tied `rologit()` design;
  - the autoscale constants, and named scales;
  - `firth` with `prior` refused;
  - finite estimates under separation, with profile intervals and penalized
    likelihood ratios available;
  - a very large scale reproduces maximum likelihood on overlapping data.
* **Size check.** A size check of the rank-ordered bootstrap is a
  verification script, run at the frozen run.

### Phase 4d (2026-11-06, author's decisions on the round-5 review and the stress rerun)

Written with the code, before any test of it ran (internal_reviews/round5/
triage.md; the stress comparison in the ledger, 2026-11-06).
* **Zero comparisons.** A comparison whose structure vector is zero (two
  available alternatives with the same design row in a case: no
  alternative-specific constants and equal attributes) is favoured by no
  direction.
  - It is left out of the linear program and the certificate, as before,
    but counted (`n_zero`; `n_comparisons` now counts every comparison).
  - With one, the data are at most quasi-completely separated: that choice's
    term stays below zero along every direction. "Complete" requires every
    comparison separated and none zero.
  - The supremum then comes from the push and refit, which leaves zero rows
    alone. Phase 4c's 0 for a design with zero rows was wrong; the round-5
    reviewers' examples have suprema log(1/2), and the LR of a term was 0
    where it is 2 log 2.
* **The separation check's scaling.**
  - After the columns, every row of the comparisons is divided by its
    largest absolute entry. This changes neither the cone nor which
    comparisons are separated, and it makes the program's 1e-9 tolerance
    relative to each row.
  - The overlap certificate projects the weights by a pivoted QR (the
    residual on the columns), not by the normal equations, which square the
    condition number.
* **Suprema for rankings.** The push makes every separated consecutive
  comparison's utility difference at least 40. The re-maximization and the
  rounds take up the comparisons a stage also contains (the e^-40 bound is
  exact for the multinomial logit only).
* **The saddle step.**
  - When the saddle check of the fit diagnostics finds a saddle, the search
    restarts from the probe point where the objective fell most along the
    eigenvector of the most negative eigenvalue, then takes the boundary
    step.
  - When that refit ends higher (by more than 1e-8 relative), the fit is
    redone from its maximum: at most three times, the multistart's record
    kept (the saddle among its ends, classified as a saddle, or as a
    boundary where coefficients diverge there; Phase 4f), with a message.
  - Otherwise the saddle warning stands, saying that the step did not end
    higher.
  - The check and the step run before the bootstrap and before any warning,
    on the information already computed for the diagnostics. Internal refits
    (`check = FALSE`) skip both.
  - Case: Liou (2024) Table 1 `mnprobit()` stopped at a saddle 8.47 below
    the maximum that a step of 0.3 normalized units reaches.
* **The boundary step** tries the structural log-parameters beyond 3 first
  together, as before. When that lowers the log-likelihood, it tries them
  one at a time, largest first, each kept with those already kept when the
  log-likelihood does not fall.
* **The bootstrap.**
  - Replicate screening treats prior fits as penalized (they keep resamples
    with a never-chosen alternative, as Firth's fits do).
  - Weighted fits' replicates (multinomial weights) pass the same checks as
    unweighted resamples.
  - `boot_lr_test()`:
    - Within the separation check's scope it keeps single-choice replicates
      with a never-chosen alternative (separation, whose supremum enters),
      as for rankings. Out of scope they fail, and the error names the
      commonest cause.
    - When the restricted fit is separated, the outcomes are drawn from the
      point its supremum is approached at (the moved point under complete
      separation, the moved and re-maximized point under quasi-complete), not
      from where its optimizer stopped. The method line says so.
    - Each ranking replicate's seed is set under the session's generator
      kinds, so serial and parallel runs agree.
* **Priors.**
  - The column scales are frequency-weighted over the rows of positive
    weight, so a weighted fit and its rows repeated by their weights get the
    same prior.
  - A scale named for a coefficient the model lacks is ignored with a
    message (it was an error), so the restricted model of a penalized LR can
    take the full model's prior.
  - `lr_test()` compares the penalties on the coefficients the fits share:
    the same family, degrees of freedom, autoscaling, constants' scale, and
    scales.
  - The documentation says:
    - that the covariates are not centred;
    - that the covariance is the inverse of the exact curvature, where
      `arm::bayesglm`'s approximate EM reports a related matrix for t and
      Cauchy priors;
    - what the Cauchy's tail allows.
* **Profile intervals.**
  - A non-finite deviance (an information underflowing far out on a
    separated coefficient) lies beyond the limit, not a failed bracket.
  - Each limit is checked by solving the restricted maximum there afresh,
    from the estimate and from zero. The warning (the interval may be too
    narrow) is raised when that maximum does not converge or its deviance
    is more than 1e-3 from the quantile. The points of the search itself
    start from the previous point's solution, and one of them stopping
    short is not reported: on `defection` such a point, next to a limit,
    stopped at a scaled gradient of 0.013 (tolerance 0.002), while the
    limits solved afresh reproduce the quantile to six decimals, as plain
    `optim()` does.
  - `rank_test()`'s penalized restricted fit warns likewise.
* **summary() of a penalized fit** marks the separated coefficients' rows
  `[sep]` and says that their Wald tests can be far off (Heinze and
  Schemper 2002); the test to use is the likelihood ratio. Under a prior the
  note says that Wald tests can have little power (Rainey 2024). (The
  wording changed in rounds 6 and 7 and after the local review.)
* **Smaller rules.**
  - `combine_test(type = "lr")` takes 0 log 0 = 0 (a never-chosen
    alternative); the Wald warning is for the Wald form only.
  - `captivity_test()` and `nested_lm_test()` warn on separated fits.
  - `mnprobit()` refuses more than 21 alternatives up front.
  - Firth's tie kernel computes the subsets' second derivatives only for the
    gradient, caps their memory at 200 MB, checks for interrupts, and names
    the case by its position in the estimation sample.
  - A restricted model with no parameters is evaluated, not optimized.
* **Monte Carlo.**
  - `firth_mc.R` gains F4, a design following that in Rainey (2024, section 3.1), at J = 2 and J = 3 (beta_s
    in {0, +/-2}), with criterion (e): at beta_s = 0 the penalized and the
    maximum-likelihood ratio (suprema) reject in at most 0.08 of samples.
  - `rologit_boot_size_mc.R` gains designs with ties and truncation,
    alternative-specific constants and a rare binary attribute; B = 99.
* **Not changed.**
  - A coordinate exactly at the search's bound stays free (round 5, I6).
    The reviewer's `>=` would hold a coordinate at the bound even where its
    gradient points inward, to save one Newton iteration.
  - The quasi-complete re-maximization in the complement of the separating
    direction (I4) waits for the run that would show a cost.
  - `remnl()`'s dependence on the base alternative is documented, not
    reparametrized.

### Phase 4e (2026-11-06, author's decision): the likelihood-ratio test of captivity

Why. On `pa2006` the score test `captivity_test()` gives z = 0.50 for
abstention while the inflated multinomial logit with a constant captive share of
abstention, the score test's own alternative, gains 10.5 in log-likelihood
on one parameter (ledger, 2026-11-06). Its expected information is
dominated by cases with small fitted probabilities, and along the captive
share the log-likelihood is convex at zero. The likelihood-ratio statistic
measures the gain itself.

Definitions, fixed before the code:
* **The model.** For an alternative k, the inflated multinomial logit with k inflated
  and an intercept-only split equation: P(y = k) = (1 - phi) + phi P_k,
  P(y = j) = phi P_j, phi constant. In the captive odds
  theta = (1 - phi) / phi >= 0 it is the dogit with captivity to k alone,
  and the logit is theta = 0. That is a boundary point of the parameter
  space, but a regular one: the model is smooth in theta at 0, and every
  other parameter is identified there. The covariate-dependent split, whose
  slopes are not identified at phi = 1, is not.
* **The reference distribution.** With q further restrictions in the
  interior (outcome coefficients of the inflated model that the logit
  lacks), 2 (l_full - l_restricted) is asymptotically the 50:50 mixture of
  chi-square(q) and chi-square(q + 1) (Self and Liang 1987;
  chi-square(0) is the point mass at 0). The p-value is
  1/2 P(chi2_q >= LR) + 1/2 P(chi2_{q+1} >= LR). With q = 0 it is
  1/2 P(chi2_1 >= LR) for LR > 0, and 1 for LR = 0. The link of the split
  is immaterial: with a constant split, logit and probit parameterize the
  same phi.
* **A statistic of zero.** The statistic is floored at 0. A value below
  1e-6 is the inflated fit reaching the logit itself, at its boundary
  (phi -> 1), and is reported as 0. A restricted log-likelihood above the
  full one by more than 1e-6 warns, as `lr_test()` does.
* **`lr_test(mnl_fit, imnl_fit)`** accepts the pair when the inflated fit
  has one inflated alternative and an intercept-only split
  (`y ~ x | 1`), the same cases, weights and choice sets, and outcome
  covariates that include the logit's. Several inflated alternatives
  (several boundary parameters, a chi-bar-square whose weights depend on
  the information) and any split covariate stay refused, pointing to
  `boot_lr_test()`.
* **`captivity_test(object, type = "lr")`.** For each requested
  alternative, the constant-share inflated multinomial logit is fitted to the logit
  fit's own design, with the split block (one intercept) added: the same
  cases, weights, offsets and choice sets, by `imnl()`'s search
  (multistart and challengers). The table holds one row per alternative:
  the statistic, df 1, and the mixture p-value. There is no joint row (the
  joint alternative has several boundary parameters, as above). The
  attribute `fits` holds the log-likelihoods, captive shares and
  convergence flags.
  - An alternative whose inflated fit fails, or that the inflated model
    refuses (never chosen; two alternatives), gets NA, with the reason.
  - Unit fixed effects are refused, as in `imnl()`. Penalized fits are
    refused, as in every test of `captivity_test()`.
  - The score form stays the default: it needs no refits.
* **Monte Carlo (`captivity_test_mc.R`).** The LR form joins the size
  designs (under the logit) and the power designs, at the 0.01, 0.05 and
  0.10 levels. The documentation reports what it finds.
* **Tests.** The statistic against the difference of an `imnl(y ~ x | 1)`
  fit and the logit, on the wide and the long layout and with weights. The
  mixture p-value against its formula. LR = 0 at a boundary fit gives
  p = 1. `lr_test()` gives the same statistic and p-value, and refuses a
  split with covariates or several inflated alternatives.

### Round 6 amendments (2026-11-06, from the round-6 review)

Phase 4e:
* **Identification.** "Every other parameter is identified" is not enough.
  The share's score must also not be a combination of the other
  parameters' scores: the fitted probabilities of k must vary in a way its
  constants and slopes cannot absorb. Under a saturated outcome design the
  information is singular. The share is then reported as NA, with LR 0,
  p = 1, status "logit" and a note.
* **Start at the boundary.** Each alternative's search also starts from the
  logit embedded at the boundary (the logit's estimates, split intercept
  30). When that start ends above `imnl()`'s result by more than the
  tolerance, the fit is redone from it.
* **Tolerance.** max(1e-6, 2e-8 max(1, |l|)) replaces 1e-6, in
  `captivity_test()` and `lr_test()`.
* **Status.** Each fit gets one of four statuses:
  - "logit": the inflated fit at phi -> 1 or LR below the tolerance; LR 0,
    p 1, and no warning;
  - "interior": the mixture p-value;
  - "degenerate": LR above the tolerance with diverging outcome
    coefficients or a saddle (the pres2004 mechanism, ledger 2026-11-06).
    p is NA, with a note pointing to `boot_lr_test()`; `lr_test()` warns
    and gives NA;
  - "failed": NA, with the reason.
* **Diagnostics.** `fits` keeps the inflated fits' warnings. `top_share`
  is the largest single case's share of the share's information at the
  logit, sum over cases of w (1/P_k - 1) over the available cases. A large
  value points to `boot_lr_test()`. On pa2006 it is 0.53, and the
  bootstrap there (ledger) found the mixture adequate.
* **`boot_lr_test()`** with the constant-split pair reports the mixture
  p-value as its asymptotic one.
* **Monte Carlo.** The script logs the raw statistic before the floor, the
  refits below the logit, the degenerate fits, the warnings and
  `top_share`. It adds two size designs with small fitted probabilities
  (coefficients tripled; a rare alternative).

Elsewhere:
* **Separation.** A comparison whose largest entry is at most 1e-10 after
  column scaling counts as zero, before the rows are scaled, so rounding
  cannot turn a zero comparison into a separated one. The certificate's
  positivity is tested in the rows' own units.
* **Priors.** `autoscale` with a weight total of at most 2 (normalized or
  survey weights) is refused. A misspelled coefficient name gets the near
  matches.
* **Engine.**
  - The greedy boundary step skips a candidate with |parameter| <= 3.
  - The saddle step runs only when the fit is checked.
  - Fits carry `diverging`, `phi_boundary` and `saddle`.
* **`summary()`'s `[sep]` note.** It has one wording per penalty (Phase 4d
  above).

### Round 7 amendments (2026-11-06, from the final review round)

* **One rule for the captive statistic** (`.pc_captive_rule()`, shared by
  `lr_test()` and `boot_lr_test()`). It replaces round 6's rule in
  `lr_test()`, which set the statistic to 0 whenever the share collapsed,
  for any q. With q further restrictions that is wrong: a collapsed fit with
  q > 0 is the logit with the further coefficients, and its ratio is theirs
  (probe r7_lr_test_collapse_q1: 0 instead of 54-121 in 8 of 15 data sets).
  Now:
  - a collapse counts as 0 only when q = 0;
  - a statistic within the tolerance is 0;
  - degenerate means diverging *outcome* coefficients or a saddle; the split
    intercept diverges at a collapse and does not count.
  `boot_lr_test()` applies the same rule to its asymptotic p-value (NA for a
  degenerate fit, noted in the method line), with the flags of the point it
  uses.
* **The LR form's identification** is decided by the score form's rule, on
  the share's information after the logit's coefficients are partialled out
  (`.pc_captivity_info()`), not by the inflated fit's conditioning.
* **The LR form's fits.**
  - The embedded start is used when the search fails too, and when it ends
    above the search by more than the tolerance in the ratio's units.
  - Each fit keeps its own warnings.
  - `n_modes` is the search's count.
* **`boot_pvalue()`** draws from the fitted model, so it calibrates only
  tests whose null is that model. It refuses `term_test()` and
  `combine_test()`, whose null is a restriction of it (probe
  r7_boot_pvalue_term_test: p.boot 0.55 for a Wald statistic of 185).
* **Zero comparisons.** A generic difference within 64 machine epsilons of
  its operands is set to 0 when the comparisons are built, and a comparison
  is zero when all its entries are. This replaces round 6's 1e-10 of the
  column's largest entry, which zeroed genuine differences in wide columns
  and reported false separation (probe r7_separation_edges).
* **The certificate** also requires positivity in the scaled units.
* **`hl_test()`** refuses penalized fits.
* **`se = "bootstrap"`** refuses frequency weights totalling at most 2.
* **The prior's near-match hint** is not given for a term that is a
  variable of the data the model leaves out.
* **Monte Carlo and verification.**
  - `firth_mc.R` gains F5 (40 cases with s = 1), since F4's five cases carry
    almost no information.
  - `captivity_lr_pa2006_boot.R` joins the frozen run.

### Local review amendments (2026-11-06, a max-effort code review run locally)
The review's report and reproductions are in
`internal_reviews/local_review_2026-11-06/`; each fix has a test in
`tests/testthat/test-localreview.R` (findings 1 and 15 on the Liou fixture in
`test-phase4d.R`).
* **Likelihood ratios with internal refits.** `term_test(type = "lr")` fits
  the model without the term as the estimators fit a model, multistart and
  saddle step included (`.pc_refit_checked()`), and warns when that fit does
  not converge or still stops at a saddle. On the Liou probit the old refit
  of the model without `human_rights` stopped at a saddle (LR 18.90,
  p = 7.9e-5); with the step it reaches -1152.021 (LR 2.03, p = 0.36).
* **A saddle step's maximum** joins `start_logliks`.
* **Zero comparisons** (`.pc_sep_struct()`): a difference within 64
  machine epsilons of the larger of its operands and the column's typical
  size (the median of its nonzero absolute values) is zero, so a tie up to
  rounding stays a tie after centring, and an outlier does not zero genuine
  differences.
* **Never-chosen alternatives.**
  - An alternative never chosen, below the chosen one in separated
    comparisons only, has probability zero at the supremum
    (`.pc_absent_alts()`): its case-specific coefficients do not count in the
    degrees of freedom of `lr_test()` and `term_test(type = "lr")`; the
    output says so.
  - `nested_lm_test()` takes it out (its coefficients and any nest it
    leaves with one alternative).
  - `captivity_test()`'s score form does not test an alternative no case
    chose; the joint statistic is computed on the correlation scale with one
    eigenvalue rule for its rank and its inverse.
  - In the separation check's scope, a resample (`se = "bootstrap"`,
    `boot_pvalue()`) with an unchosen alternative is dropped only when the
    linear program finds it separated.
* **`rank_test()`** leaves out coefficients that a block of stages does not
  identify, with a warning and fewer degrees of freedom; a penalized fit is
  then refused.
* **The captive rule** (`.pc_captive_rule()`, now used by `lr_test()`,
  `boot_lr_test()` and `captivity_test(type = "lr")`): a collapsed share
  silences only the "restricted fit is higher" warning; a statistic above
  the tolerance stands. `boot_lr_test()` takes the full fit's flags from the
  engine and reports the statistic its p-values use; its replicates use the
  same tolerance.
* **Smaller rules.**
  - `combine_test()` refuses the same alternative twice.
  - `boot_pvalue()` leaves an NA replicate statistic out of its row only
    (attribute `n_na`).
  - A failed separation check in a supremum warns.
  - Profile intervals: the root to 1e-6 of a standard error.
  - `firth = TRUE` refuses weights totalling at most 2; the bootstrap's
    weight check runs before the fit and covers `remnl()`.
  - `hl_test()`'s degrees of freedom count the groups and alternatives with
    expected counts.
  - Prior `df` is compared numerically.
  - Firth's untied stage guards its 64-bit mask beyond 64 alternatives.
  - `.pc_best_fit()` skips its further starts for a concave log-likelihood
    that converged.
* **`summary()`'s `[sep]` note** of a Firth fit cites Heinze and Schemper
  2002 (the author's decision).
* **F6** in `firth_mc.R` (the author's decision): the covariates of design
  911 of Rainey's (2024) replication archive, with criterion (f) on the ML
  LR rate.

### Phase 4f (2026-11-10, author's decision): the record of the search's ends

Why. `n_modes` counted the distinct log-likelihoods, rounded to 0.1, of the
converged polished fits. An end at which coefficients diverge can meet the
gradient criterion (the gradient vanishes along the diverging direction), so
it counted as a maximum: on `fluvax` the count was 3, of which one end has
`infl_knowledge` diverging and the split at one regime for 87 percent of the
cases (verification/probes/2026-11-10/nmodes_classify.R and .txt). Rounding could also
merge two maxima less than 0.1 apart.

* **The record.** Every end of the search (the polished multistart fits, the
  challengers and the cross-check that are kept, before a saddle step the
  original search's, and the reported fit) enters `.pc_search_record()`.
  Ends are taken in decreasing order of log-likelihood. Two ends are the
  same point when their log-likelihoods differ by less than 0.01 and either
  every estimate differs by less than a tenth of the first's
  inverse-information standard error (by less than 0.01 in the standardized
  coordinates where it has none) or every fitted probability of every case
  by less than 1e-4 (the estimate test is tried against every point first).
  The second test is the data's own: two ends with the same probabilities
  have the same likelihood whatever their estimates, so they lie on a ridge
  of non-identification, and the point's `detail` says how many ends joined
  it that way. The end the reported fit came from belongs to the reported
  point even when the boundary step has moved the fit.
* **Amendment (2026-11-10, author's decision, before release).** As first
  implemented (build 28b870609bea) the record grouped by the estimates
  alone. The stress rerun of that build found the saturated designs of
  DiGiuseppe et al. (2024) (one binary covariate) reported as up to eleven
  "distinct local maxima" at one log-likelihood: the search ends at
  different points of a curved ridge, with estimates several units apart
  and fitted probabilities equal to 1e-8 (the old rounded count said 1).
  A test on the straight line between two ends would not do, since the
  log-likelihood falls along the chord of a curved ridge (by up to 25 at
  the midpoint there). The tolerance 1e-4 is four orders of magnitude above
  the agreement of the ends on the ridge and a tenth of the last decimal
  that `summary()` prints for a probability. The probability test can only
  merge points the estimate test kept apart; on a simulated saturated
  design (one binary covariate, J = 3, two inflated alternatives, six seeds
  at n = 600 and 1500) the count is 1 where the estimate test gave 3 to 11
  (verification/probes/2026-11-10/ridge_probe.R, ridge_sim.R and their
  outputs).
* **The classification**, at each point's highest end, by the rules of the
  reported fit's flags (`.pc_end_type()`), in this order:
  - `"boundary"`: a diverging coefficient (the exact check of
    `separation()` where it applies, the size-and-Wald rule otherwise), a
    structural log-parameter beyond 7, or a split equation at one regime;
  - `"saddle"`: `.pc_saddle()` on the normalized Hessian over the
    parameters with curvature;
  - `"not converged"`;
  - `"maximum"` otherwise (an ill-conditioned information is noted in
    `detail`).
  The reported point takes its type from the fit's own flags.
* **The count.** `n_modes` is the number of points of type `"maximum"`; 0
  when there is none (a separated logit, a fit at a boundary with no
  interior maximum found). A penalized fit keeps its count of the two runs'
  maxima and no record; refits with `check = FALSE` keep neither (`n_modes`
  NA).
* **What is reported** (the author's decision of 2026-11-10): the highest end
  of the search, maximum or not. An inflated fit's message names the maxima
  when there are several, with the number of other ends, and when the
  reported fit is not a maximum but a lower end is, it says so and points to
  `search_theta`, whose row refits that maximum as `start`.
* **What does not change.** The search, the reported fit, its estimates,
  log-likelihood, standard errors, flags and warnings, and `start_logliks`
  are as before; on twelve fits (three bundled data sets, six replications
  of the captive design without an exclusion restriction, three logits)
  they are identical to build 9446ac9410f9 bit for bit. The classification
  costs one Hessian per point other than the reported one (0.1 to 0.3 s on
  the bundled data, against fits of 8 to 20 s).

## 2. Model and parameterisation

Case i = 1..n, alternatives j = 1..J (the levels of the response), base b,
availability a_ij in {0, 1} (all ones in the wide layout).

    V_ij = o_ij + x_i' beta_j + w_ij' gamma,      beta_b = 0
    P_ij = a_ij exp(V_ij) / sum_k a_ik exp(V_ik)

x_i: case-specific covariates (intercept included: the alternative-specific
constants), beta_j alternative-specific. w_ij: alternative-specific covariates,
gamma generic. o_ij: offset. This is McFadden's conditional logit with
Theil's multinomial logit as the special case w = 0, which is how Fok and Paap
(2025, eq. 2) write it.

Inflation (imnl): s_i in {0, 1}, P(s_i = 1 | z_i) = phi_i = F(z_i' delta), F
logistic (default; Bagozzi and Marchetti 2017, eq. 5) or standard normal.
Inflated category k (default: the base):

    P*(y_i = k) = (1 - phi_i) + phi_i P_ik
    P*(y_i = j) = phi_i P_ij,   j != k

Log-likelihood sum_i w_i log P*(y_i), w_i frequency weights. The MNL is the
boundary case phi -> 1 (delta_0 -> +Inf), so MNL vs imnl is a Vuong / AIC /
BIC comparison, not a standard LR test (as in the ZiOP literature).

Parameter vector theta (internal and reported order):
1. beta, alternative-major: for each non-base alternative in level order, the
   p case-specific coefficients, named `"<alt>:<term>"`;
2. gamma, named `"<term>"`;
3. delta (imnl only), named `"infl_<term>"`.

Per-case score (C++):

    c_i = 1                       if y_i != k
    c_i = phi P_ik / P*_ik        if y_i == k
    r_ij = c_i (1{j = y_i} - P_ij)                           dlogP*/dV_ij
    u_i  = f(eta_i)/phi_i             if y_i != k            dlogP*/deta_i
    u_i  = -f(eta_i)(1 - P_ik)/P*_ik  if y_i == k
    s_i(beta_m) = r_im x_i;  s_i(gamma) = sum_j r_ij w_ij;  s_i(delta) = u_i z_i

Gradient of the weighted log-likelihood: sum_i w_i s_i. Hessian: central
finite differences of the analytic gradient in the column-scaled
parameterisation (the `iop` recipe); robust covariance
H^{-1} (sum_i w_i s_i s_i') H^{-1}; cluster: sum over clusters of the summed
weighted scores. All derivatives checked against numDeriv in the tests.

## 3. Interfaces

Two data layouts, one formula grammar: the parts of the right-hand side are
separated by `|` (package Formula) and **the last part is the inflation
equation for `imnl()`**.

Wide layout (one row per case, the response is the chosen alternative):

    mnl(vote ~ educ + age + female, data = d)
    imnl(vote ~ educ + age + female | educ + age, data = d, inflate = "abstain")

Long layout (one row per case-alternative pair; `id` and `alt` name the case
and alternative columns; the response is the 0/1 choice indicator):

    mnl(chosen ~ cost + time | income, data = long, id = "case", alt = "mode")
    imnl(chosen ~ cost + time | income | income, data = long, id = "case", alt = "mode")

Part 1 of the long layout: alternative-specific covariates with generic
coefficients; its intercept means alternative-specific constants. Part 2:
case-specific covariates (constant within case; checked), alternative-specific
coefficients. Alternatives absent from a case's rows are unavailable to it.
Alternative-specific coefficients on an alternative-specific variable are the
interaction `w:alt`, so no third covariate block is needed.

Common arguments (names follow `iop`): `ref` (base alternative; default the
first level), `inflate` (imnl; default `ref`), `link` for the split equation
(`"logit"`, `"probit"`), `weights`, `offset`, `fe`, `se = c("analytic",
"opg", "robust", "cluster", "bootstrap", "none")`, `cluster`, `nboot`, `cores`,
`start`, `maxit`, `reltol`.

Semantics carried over from `iop`/`underdisp`: `weights` are frequency weights
(a fit with weight 2 equals a fit on duplicated rows, in estimates, standard
errors, information criteria, Vuong, averaged QoIs); survey (probability)
weights are not supported (the robust meat is the frequency-weight form
sum w s s', see mnl.Rd). Rows with a missing value in any model variable
are dropped before fitting; `object$data` holds the estimation rows. A
character response is refused (alphabetical ordering hazard). `cluster` alone
selects `se = "cluster"`.

Fixed effects (`fe = "unit"`): one dummy per unit in every non-base equation.
A unit that never chooses alternative j has that dummy at -Inf; the limiting
likelihood is the one with j unavailable to that unit, so those rows get
a_ij = 0 and the dummy is omitted (message). A unit that never chooses the
base likewise has the base unavailable; its constants are then identified only
up to a common shift, so one of its dummies is omitted. The common constant of
j refers to a reference unit: the first that chose both j and the base, else
(breadth first) one that links j to an already anchored alternative, in which
case that unit omits the dummy of the alternative it anchors through; an
alternative not connected to the base through the units' choices is refused
by name. Captive cases (one alternative available in the data) do not count
as choices, cases of weight zero are dropped first, and cases left with a
single alternative drop out; dropping cases can change what a unit chose, so
the rules and the drop are repeated to a fixed point (in the fit and in every
half-panel or bootstrap design), a unit that chose a single alternative drops
out entirely and links nothing, and an alternative chosen only by such units
(or only in captive cases) is refused with that cause. Without
alternative-specific constants nothing is anchored: every unit keeps all its
dummies, one fewer for a unit that never chose the base. imnl() refuses
fixed effects (author's decision 2026-10-18, ledger): with a constant per
unit in the outcome equations the level of phi is identified only through
the mixture's form, the routine choosers of k drop out, units with few k
choices put them all in regime 0 (their constants diverge), and the split
coefficients came out 25-42 percent too large at T = 12-20 (imnl_fe_mc,
build 1245fc27f8e3; archived; the split intercept at T = 12 by 42 percent,
corrected from "25-40" in round 5). The bootstrap of a fixed-effects fit resamples units (a unit drawn
twice is two units) and covers the slopes; a replicate whose fixed effects
cannot be built counts as failed. Short panels: the incidental-parameters
warning, the split-panel jackknife, and remnl() on the Mundlak-augmented
formula.

## 4. Numerical strategy

* Starting values: MNL from zero (the log-likelihood is globally concave given
  full column rank and no separation); imnl from fourteen starts (the MNL fit
  or an MNL down-weighting the inflated cases, with a logit of 1{y != k} on z
  for delta at full, halved, doubled, tripled or zero slopes and several
  shifts of the split intercept), screened by a short quasi-Newton run and
  every one polished (2026-10-20, stress finding F12: a cap of five polished
  starts let added starts lower the maximum); with one inflated alternative
  the same search under the other split link, its estimates rescaled by 1.6
  and polished under this link (`crosscheck = TRUE`; not the default, at the
  author's decision of 2026-10-20, since it doubles the time again); then
  challengers rescaling delta by 0.5, 1.5 and 2, repeated from a new
  incumbent while one improves (at most three rounds). The polished,
  cross-check and challenger log-likelihoods are stored (`start_logliks`),
  the screened ones too (`screen_logliks`), and the cross-check
  (`crosscheck`); the distinct points where they end, each classified, are
  in `search` and `search_theta`, and `n_modes` counts the local maxima
  among them (Phase 4f). The challengers run from the best polished fit only, so
  the search is not monotone in its starts where the split coefficients
  diverge: the likelihood then has many flat plateaus (each a hard
  classification of a few cases into the regimes), and which one the
  challengers reach depends on the incumbent (stress finding F19: two fits,
  0.14 and 1.24 units). Such fits are flagged `boundary`; challenging the
  best three maxima would remove the dependence at 40 to 100 percent more
  fitting time, which the author declined (2026-10-21).
* Optimiser: damped Newton first for the plain MNL (analytic Hessian), BFGS
  (`optim`) with the analytic gradient on column-standardised parameters
  otherwise, then a Newton polish (analytic Hessian for the logit family,
  finite-difference Hessian of the analytic gradient otherwise), a step
  accepted when the objective does not increase (within 1e-14 relative), the
  step capped at 2 standardized units, convergence declared when the largest
  scaled gradient is below 1e-5 max(1, |ll|) (`.pc_gtol()`, the one rule of
  the polish, the separation rule and the search finish below). The search
  holds the structural log-parameters within +/- 8, the value the boundary
  step moves them to (`.pc_str_bound`): a point beyond is evaluated at the
  bound with a zero gradient in that coordinate, whose Hessian row and
  column are those of a fixed parameter (stress finding F22). Past the
  boundary flag the likelihood is usually flat; where it still rises at the
  bound (the gradient points out of it by more than the tolerance), the
  boundary warning names the parameter and says that the reported
  log-likelihood is below the model's supremum, so a comparison with the
  restricted model understates it (on Palik et al. 2025 by 0.67 or more; the
  author chose this over a wider bound, a measured shortfall or no bound,
  2026-11-04). Under hev()
  and mnprobit() BFGS and the Newton polish run in coordinates that measure
  each non-base alternative's coefficients in units of the standard
  deviation of its utility difference from the base, relative to that at the
  iid start (`.pc_search_map()`; sqrt((1 + theta_j^2) / 2) and sqrt(Sigma_jj)):
  where an alternative's scale runs to infinity its coefficients grow in
  proportion, a ridge curved in the coefficients and straight in these
  coordinates (stress finding F23). Convergence is judged in the model's
  coordinates; where the search's point misses the tolerance there, up to 25
  Newton iterations in the model's coordinates finish the polish (Newton in
  the model's coordinates alone is ill-conditioned along the ridge, and the
  search's gradient, which weighs the coefficients by the utility's standard
  deviation, has another tolerance).
* Time limits: R raises a limit set by the caller (`setTimeLimit()`, or a
  timeout built on it) as an error at the next interrupt check and clears it
  as it fires; every error handler of the package (a failed start, a
  singular solve, a failed replicate) signals that error again, so the limit
  stops the fit rather than one attempt, after which nothing would (stress
  finding F21). In compiled code with interrupt checks (the remnl(), hev()
  and mnprobit() likelihoods, every 100 cases in the HEV and probit loops)
  the limit arrives as an interrupt, as does a user's interrupt.
* Diagnostics stored on the object, computed from the information matrix
  whatever the standard errors are: `converged` (false also at a saddle:
  a negative eigenvalue of the diagonally normalized Hessian over the
  parameters not at a boundary, beyond ten times the asymmetry of a
  finite-difference Hessian over the same parameters, along whose
  eigenvector the log-likelihood rises by more than the optimizer's relative
  tolerance at a step of 1e-3 to 1 normalized units; on the curved flat ridge
  of an unidentified model the eigenvalue can be slightly negative while the
  log-likelihood falls both ways, stress finding F18; a negative eigenvalue
  whose predicted gain the smallest step would show, with the log-likelihood
  curving down there instead, is reported as an inaccurate Hessian),
  `boundary` (a coefficient whose one-SD effect exceeds 15, or a 0/1 column
  whose coefficient exceeds 15, and whose Wald |z| is below 1.96, a log scale /
  inclusive value / sd beyond +/- 7, or the inflated model collapsing onto the
  MNL), `ill_conditioned` (reciprocal condition number of the normalized
  Hessian below 1e-9, or a generalised inverse used), the rank check of every
  design block, and the scale-free identification check of the choice
  design. `singular` (2026-11-10, the author's decision) records that the
  generalised inverse was used; `summary()` then says in its footer that some
  parameters are not identified and that their estimates and standard errors
  are not meaningful, since their rows in the table look like the others (the
  stress battery of build 9446ac9410f9 has 208 such fits, 18 of them not at a
  boundary; a simulated saturated design printed a z of 127 for a split
  intercept the data do not identify).
* Log-sum-exp with the row maximum; log P* for the inflated category through
  a two-term log-sum-exp; `plogis`/`pnorm` with `log.p` for the split equation;
  1 - P_ik computed as the sum of the other probabilities.
* Scale-free checks: the standardization floor is relative to the column's
  magnitude (a covariate measured in tiny units is still scaled), and the
  identification check compares each alternative-specific column's within-case
  variation with its own centred sum of squares (a shifted column, w + 1e9, is
  identified exactly as w is).
* Probit probabilities below 1e-300 are held at that floor with a zero
  gradient (the floored log-likelihood is flat), for the bivariate-normal and
  the GHK paths alike. The remnl mode search runs up to 100 damped Newton
  steps. The HEV adaptive integration counts the integrals that stop at its
  panel or depth limit, and the fit warns when any does at the estimate.
* remnl population-averaged probabilities (no data narrow the integrand):
  Gauss-Hermite with max(nodes, 7 sigma^2) points per dimension on the 2^k
  ladder up to 128, then the trapezoidal rule on [-9, 9] standard deviations
  with a power-of-two number of points at least 18 sigma (step at most
  1 / sigma: error about exp(-2 pi^2 / (sigma h)), below 1e-8), so the rule
  grows linearly in sigma and needs no eigenproblem; Gauss-Hermite rules are
  memoised. At most 5e4 points in all (or nodes^(J-1)); the fit warns with
  the error bound when that cap binds (from J = 4 at large sigma).

## 5. IIA tests: definitions fixed in advance

Hausman-McFadden. Restricted sample: cases with y != d (alternative d
dropped), restricted MNL with the same base (or, when d is the base, base c =
the first remaining level; the full model is re-expressed with base c by the
exact linear map beta_j^(c) = beta_j - beta_c, V^(c) = A V A'). Statistic
(beta_r - beta_f)' (V_r - V_f)^{-1} (beta_r - beta_f) on dim(beta_r) df;
`vcov = "difference"` inverts the difference (Moore-Penrose with a warning when
it is not positive semi-definite; negative statistics reported as such with
the Hausman-McFadden/Long reading "support for the null", p-value 1; a
positive statistic with an indefinite difference is unreliable, since the form
over-rejects); kept as the default because it is the form the literature and
Stata's `hausman` report; `vcov = "joint"`
uses Var(beta_r - beta_f) = V_r + V_f - C - C' with the cross term from the
stacked per-case scores (suest), positive semi-definite by construction and
valid with robust/cluster covariances.

Small-Hsiao. Stratified random half-split (by chosen alternative). Full MNL
on half 1 and on half 2: theta_1, theta_2; theta_12 = sqrt(1/2) theta_1 +
(1 - sqrt(1/2)) theta_2. Restricted MNL on half 2 without alternative d:
theta_r. Statistic -2 [l_2^r(theta_12 restricted to the remaining parameters)
- l_2^r(theta_r)] on dim(theta_r) df. The prior project's lineage (restricted
model on both halves, no weighting) is 2 chi-square and is not implemented.

Fok-Paap. Pairs {j, m}: base pairs (J - 1) or all pairs J(J-1)/2. Pairwise
binary logit P(y = j | y in {j, m}) = Lambda(V_j - V_m). CL estimator maximises
the sum over pairs and cases in the pair of the binary log-likelihood;
Godambe covariance (H^{-1} J H^{-1})/n with J from the per-case summed pair
scores. Q1 = (theta_ML - theta_CL)' V_d^{-1} (theta_ML - theta_CL), df =
dim(theta) (or the rank of V_d when a generalized inverse is needed). Default
V_d (`vcov = "joint"`): the outer product of the stacked influence functions
psi_i = V_ML s_i^ML - H_CL^{-1} s_i^CL (summed by cluster with a cluster),
positive semi-definite by construction; `vcov = "difference"`: V_CL - V_ML,
F&P eq. 18, indefinite in nearly every sample because the CL estimator is
nearly efficient under the MNL (ledger, q1_covariance_check.R). GMM moments g_i = stack over pairs of
(1{y_i = j} - Lambda(V_ij - V_im)) 1{y_i in {j, m}} zeta_ijm with instruments
zeta = (x_i, w_ij - w_im); efficient GMM from the MLE with the empirical
weight matrix (F&P eq. 15): `gmm = "two-step"` (default; one weight update,
Hansen's standard two-step estimator), `gmm = "iterated"` (F&P's practice;
warns when the iteration does not settle, which happens with many moments
relative to the sample, e.g. pa2006 with 42 parameters and 84 moments);
Q2 = n g_bar' W g_bar at the GMM estimate, df = #moments - dim(theta). Option
`gmm = "projected"`: the C-statistic at the MLE,
n g_bar'[W - W D (D' W D)^{-1} D' W] g_bar, same df.

### 5.1 New tests: design and Gate 1 (written 2026-10-18, before any code)

Status: designed, not pursued as a paper (author's decision 2026-10-18, after
the prior-work finding below); an implementation in the package is not
scheduled. The search for a new IIA test is closed (author's decision the
same day): on one choice set with case-specific covariates IIA is not
identified apart from the logit's functional form (any positive choice
probabilities are a logit in utilities log(P_j / P_b)), and a local power
envelope (ledger, 2026-10-18; not checked by simulation, by decision) puts the
existing tests within a factor of two to four of the best test that knows the
direction, Small-Hsiao excepted, while that best test itself needs thousands
of cases for moderate violations. Tests that use choice-set variation (the
two-data-set idea of McFadden, Tye and Train 1977) or second choices (the
rank-consistency tests) exist and fit few political-science data sets.

Source and authority: the author's validation plan for the link test
(2026-09-02, not in this repository) fixes the gates and their pass/fail
rules; nothing here changes a rule. This section fixes the definitions for
polychoice and records Gate 1 (theory and invariance). Amendments to the
plan are marked as such, with the reason, before any result they could
affect.

Common setup. An `mnl()` fit with base b, full choice sets, no `fe`; e_ij =
V_ij - V_ib for j != b, the fitted utility differences from the base
(`predict(type = "link")` minus its base column; they include the generic
terms w_ij' gamma in the long layout). X is the case-specific design (the
regressors that carry alternative-specific coefficients).

Link test (L2, the definition polychoice implements). Q_i = all products
e_im e_il, m <= l, m, l != b: (J - 1) J / 2 case-specific columns. The
augmented model adds Q with alternative-specific coefficients gamma_j to every
non-base equation; H0: gamma = 0. LR = 2 (l_aug - l_MNL). Degrees of freedom
(J - 1) r, where r is the rank of Q after projecting out the column space of X
(weighted by the frequency weights; columns scaled by their norms, singular
values above 1e-7 of the largest). The augmented model adds an orthonormal
basis of that r-dimensional span, not a subset of the products: the
likelihood, and so LR and LM, does not depend on the basis, and a subset is
ill-conditioned when two equations' slope vectors are nearly collinear (the
first probe, which kept a pivoted subset, drew divergence warnings in a design
with a cosine of -0.995 between two slope vectors; the basis version does
not).
Score form (L3): the Rao statistic of the same restriction at the MNL fit,
LM = s_g' [I_gg - I_gt I_tt^{-1} I_tg]^{-1} s_g, with the MNL information,
which does not depend on y (canonical link), so the observed and expected
forms coincide; robust and cluster variants replace the bracket by the outer
product of the per-case (per-cluster) scores of the gamma block after
partialling out the theta block. J = 2 is refused (IIA is vacuous there).
Varying choice sets are refused in the first version (e_ij is undefined for
an unavailable j; a definition is a design decision for later).
The squares-only variant (L1, Q = {e_im^2}) is not base-invariant (below) and
is not exported; verification scripts keep it only to reproduce the
September 2026 runs.

Prior work (literature check 2026-10-18; it corrects the September 2026
review's "not proposed for IIA before"). Murteira and Ramalho (2016,
Econometric Reviews 35(4), 515-552, doi:10.1080/07474938.2013.806849; read in
the March 2012 preprint, section 4.1.1) augment multinomial-logit-type mean
models with "powers and cross-products of all the elements of X beta" in every
argument, a RESET-type test that is in substance L2 (their worked versions add
own squares only), as LM tests by artificial regression under the
multinomial quasi-likelihood, which is the MNL likelihood for 0/1 choices,
and report its power against nested-logit and random-parameter designs; they
frame it as a general functional-form test, not as a test of IIA, and their
simulations are fractional responses. Murteira (2016, Comm. Statist. Theory
Methods 45(24), 7367-7375, doi:10.1080/03610926.2014.980518) generalises the
goodness-of-link test to multivariate models with multinomial choice as the
illustration (abstract only; to be read). McFadden (1987) sets out
auxiliary-regression tests for "functional misspecification" of the MNL
(abstract only; to be read). McFadden and Train (2000, JAE 15(5), 447-470)
give the LM test for mixing through artificial variables in the fitted
probabilities (Gate 1b's random-constant direction), a different test. Not
found in the literature checked: the test analysed as a test of IIA, the base
dependence of the squares-only form against the invariance of the full
quadratic, the rank df, and a maximum-over-nests bootstrap version of
McFadden's nested-logit test (in the check's reading of the 1981 working
paper of Hausman and McFadden, the tests over trees are dependent and no
prescription is given; Murteira and Ramalho's NESTEDg is a joint chi-square
over partitions). Any novelty statement is limited to
these and waits for McFadden (1987) and Murteira (2016) to be read in full.

What the link test is. With case-specific covariates only, a violation of IIA
shows only as a departure of P(y | z) from the multinomial-logit form in z;
the test rejects for any such departure, an omitted nonlinearity in z
included. It is a functional-form test with power against IIA violations, as
the Fok-Paap moments are (Gate 3d measures how much each test reacts to an
omitted x^2).

Max-over-nests test (L4). For every admissible nest S (2 <= |S| <= J - 1, the
other alternatives singletons; 2^J - J - 2 nests: 3 at J = 3, 10 at J = 4, 25
at J = 5) the McFadden (1987) score statistic LM_S of `nested_lm_test()`;
statistic max_S LM_S, calibrated by the parametric bootstrap from the fitted
MNL (the `boot_pvalue()` machinery); the maximising nest is reported as a
direction, not as a test of that nest.

Gate 1a, base invariance (derived). A change of base from b to c maps the
index differences linearly and invertibly: e'_j = e_j - e_c (j != b, c),
e'_b = -e_c; the MNL fit is equivariant under the same map (the exact
re-basing already used for Hausman-McFadden). Homogeneous quadratics in e'
are homogeneous quadratics in e and conversely, so Q and Q' span the same
columns; the augmented models are the same model, and LR, LM, and df are
identical. L1 is not invariant: e'_j^2 = e_j^2 - 2 e_j e_c + e_c^2 needs the
cross product, which L1 lacks (already at J = 3). Confirmed numerically
(verification/probes/2026-10-18/link_gate1a.R, n = 3000, every base): L2's LR
agrees across the four bases to 1e-12 at J = 4 with three covariates (12.689,
df 18), while L1's moves from 2.52 to 7.25 (df 9). L1 and L2 coincide when the
squares already span the quadratic forms (J - 1 >= k(k + 1)/2: one or two
covariates at J = 4, binary covariates).

Gate 1a', the degrees of freedom (derived; amendment to the plan's
definitions, written before Gate 2; the probe above finds the derived ranks:
r = 6, 1, 3 at J = 4 with three, one, two continuous covariates, r = 1 at J = 3
with two binary covariates). The plan's df (J - 1)^2 for L1 and
(J - 1)^2 J / 2 for L2 count columns, not directions. With e_j = a_j + x' b_j,
e_m e_l = constant + linear in x + (x' b_m)(x' b_l): beyond the span of X the
products add only quadratic forms in the continuous covariates, at most
k(k + 1)/2 of them for k continuous covariates, fewer with binary covariates
(x^2 = x) or few equations. The df is therefore (J - 1) times a numerical
rank. Check against the September 2026 null runs (300 replications each,
the plan's diagnostics of 2026-09-02, not in this repository, L1): at
k_x = 1 the mean statistic is 2.00 (J = 3) and 2.97 (J = 4), the rank df 2 and
3; rejection at 5 percent is 0.050 and 0.047 at the rank df against the 0.007
and 0.000 reported at the column-count df 4 and 9. At k_x = 4 and 8 the two
counts agree and so do the rates (0.050-0.077). The "undersized link test at
one covariate" of September was the df count.

Gate 1b, local directions (derived; numerical R^2 to follow). Each departure's
score at the MNL point has the form sum_j (y_ij - P_ij) c_ij(z), so its
overlap with the link test is the overlap of c with Q after the MNL score is
partialled out (residuals are orthogonal to constants and to X at the MLE).
* Nested logit, nest S, lambda -> 1: d log P_j / d lambda = -log P_{j|S} -
  P_S H_S for j in S and -P_S H_S otherwise (H_S the within-nest entropy), so
  c_ij = -1{j in S} log P_{ij|S}: McFadden's (1987) artificial variable. For
  S = {m, l}, -log P_{m|S} = log 2 - d/2 + d^2/8 + O(d^4) with d = e_m - e_l;
  the constant and linear parts are annihilated by the MNL score, and the
  quadratic part, (y_m - P_m + y_l - P_l) d^2 / 8, lies in L2's span (d^2 =
  e_m^2 - 2 e_m e_l + e_l^2; e_b = 0 when the base is in S). L2 contains the
  nested-logit direction to second order in the within-nest index difference;
  L1 contains it only when the base is in the nest. The higher-order terms
  lie outside the span, which is where L4, the exact direction, gains.
* Random alternative-specific intercepts (case level, covariance sigma^2 S0,
  sigma^2 -> 0): the score in sigma^2 is (1/2) sum_jk S0_jk [(y_j - P_j)(y_k -
  P_k) - (1{j = k} P_j - P_j P_k)], which with one-hot y equals sum_j (y_j -
  P_j) c_j with c_j = S0_jj / 2 - sum_k S0_jk P_k: the McFadden-Train (2000)
  artificial variable for random constants, linear in the fitted
  probabilities. P is close to affine in z when effects are weak, so this
  direction is close to the MNL score: random-constant heterogeneity with
  case-specific covariates is weakly identified whatever the test. It is the
  direction of the retired Moran moments (residual times fitted probability),
  whose 10-13 percent power in September is consistent with this.
* Heteroskedastic extreme value, scales -> 1: the score has the same form with
  c a nonlinear function of the index differences (a Gumbel integral); its
  overlap with Q is computed numerically from the `hev()` kernel at unit
  scales (pending).
* Multinomial probit correlation: the probit does not contain the MNL, so
  there is no score at the MNL point (amendment to the plan's 1b, written
  before any result: the probit enters only Gate 3's power runs).

Gate 1c, generated regressors (derived). The gamma score at the fit is
s_g = sum_i (y_i - P_i(theta_hat)) (x) Q_i(theta_hat). The derivative with
respect to the theta inside Q is sum_i (y_i - P_i) (x) dQ_i/dtheta, whose mean
is zero under H0 (the residual is mean-zero given z), so replacing Q(theta) by
Q(theta_hat) changes s_g / sqrt(n) by o_p(1); the theta inside P is handled by
the partialling in the LM statistic. LR and LM are therefore asymptotically
chi-square on the rank df under H0, the condition of Pagan (1984): the
generated regressors' coefficients are zero under the null. Gate 1c's
large-sample check (N = 50,000, mean statistic / df and the 95th percentile)
verifies this numerically.

Gate 1d, calibration. Every bootstrap in polychoice is the parametric
bootstrap from the fitted null (y regenerated from the fitted MNL with X
fixed), so the observed statistic is never exchangeable with the bootstrap
draws under the alternative. The plan's shuffle check (y permuted against x,
under which the constants-only MNL is true) and the strong-violation check
(p-value to zero as N grows) are run with Gate 2.

Verification on implementation (in addition to the plan's Gates 2-6): the
augmented-model gradient against numDeriv; LR against LM at large n; the base
invariance as a regression test; the rank df with one continuous covariate,
with binary covariates only, and at k(k + 1)/2 < (J - 1) J / 2.

## 6. Verification protocol (the standard set by iop and underdisp)

1. Oracle agreement: `mnl()` equals `nnet::multinom()` (coefficients,
   log-likelihood, covariance) and `mlogit::mlogit()` (conditional logit with
   alternative-specific covariates, varying choice sets) to 1e-6 on simulated
   and bundled data.
2. Published results: Bagozzi and Marchetti (2017) Tables A.3/A.4 and the
   log-likelihoods (MNL -2110.59 / BIMNL -2069.97 on pa2006; MNL -553.50 /
   BIMNL -545.27 on pres2004) reproduce from the bundled data.
3. Derivatives: analytic gradient and scores against numDeriv on every
   layout x inflation x link combination; weight-duplication identity; offset
   identity; re-basing invariance (same log-likelihood and probabilities for
   every `ref`).
4. Parameter recovery on `rmnl()` data (bias within Monte Carlo error at
   n = 5000); coverage of analytic and robust intervals.
5. Tests: HM invariance to the choice of base; HM/SH/Q1/Q2 null rejection
   rates within [0.03, 0.08] at nominal 0.05 in a 300-replication check on a
   J = 3, k_x = 4 design (a band fixed before the check was run), with the
   per-replication statistics saved under `verification/`.
6. Adversarial review rounds (simulated referee panel: applied user,
   econometrician, numerical analyst, reproducibility reviewer) before any
   release, with every finding logged in `verification/ledger.md`.
7. `R CMD check --as-cran` clean on Windows (local), GitHub Actions matrix,
   and the sanitizer workflow from `underdisp` for the C++ (the last two wait
   for the GitHub repository, which the author has deferred; they are
   release requirements).

## 7. Decisions log

* 2026-10-13 Name: polychoice (polytomous unavailable, CRAN policy).
* 2026-10-13 Grammar: multipart formula, last part = inflation; long layout
  via `id`/`alt` arguments on a plain data frame (no dfidx dependency).
* 2026-10-13 Split equation: logit default (the published BIMNL), probit
  option; inflated category any level, default the base.
* 2026-10-13 Hessian by finite differences of the analytic gradient (as iop),
  not analytic, to keep one code path for both models.
* 2026-10-13 IIA tests take a fitted `mnl` object (not formula + data), so a
  test is always on the model the user actually reports.
* 2026-10-13 Hessian: analytic in C++ after all (the finite-difference Hessian
  was off by up to 15 percent on TravelMode under an inverted column scaling);
  one code path still, since the inflated model's second derivatives are
  closed-form.
* 2026-10-13 Standardization: ts = theta * scale, with the within-case
  standard deviation for alternative-specific columns; boundary flag on
  |theta * scale| > 15 (effect of a one-SD change).
* 2026-10-13 Multistart for imnl(): twelve starts (ten on a grid of beta /
  delta starts plus two challengers of the incumbent) because the inflated
  likelihood is bimodal on both bundled data sets; `n_modes` reported.
* 2026-10-13 Hausman-McFadden default stays the classic variance-difference
  form (it is what the literature reports), with the psd flag and a pointer to
  the joint form; in the 60-replication null check the difference form
  over-rejects (0.12), the joint form is conservative (0.02), Small-Hsiao
  (0.03) and Q2 (0.05) are at nominal size.
* 2026-10-13 Fok-Paap with base pairs and case-specific covariates only has no
  overidentifying restrictions; `fp_tests()` then reports Q1 only, with a
  warning (an error when only the GMM test was requested).
* 2026-10-14 Fok-Paap Q1 default covariance: the joint (stacked influence
  function) estimate of Var(theta_ML - theta_CL), positive semi-definite by
  construction; F&P's variance difference V_CL - V_ML stays available as
  `vcov = "difference"`. The CL estimator is within a few percent of efficient
  under the MNL, so the difference of the two separately estimated matrices is
  indefinite in nearly every sample (44 percent negative statistics, size
  0.033 / 0.042 / 0.049 at 1 / 5 / 10 percent); the joint form had
  0.015 / 0.062 / 0.112 (`verification/q1_covariance_check.R`, 1000
  replications).
* 2026-10-14 Small-Hsiao rows carry `separated` (a fit of the split has
  diverging coefficients); the summary's `median_p` is the chi-square p-value
  of the median statistic, as in the battery.
* 2026-10-14 boot_lr_test(): the Vuong row's normal p-value is NA for the
  mnl / imnl pair (no valid normal reference at the boundary), as the LR's
  chi-square one is; vuong() prints a note for that pair.
* 2026-11-05 Firth's penalty (author's decision, before the release review
  rounds): `firth = TRUE` for mnl() and rologit() only, ties included
  through the expected information of the tie group (the author asked for
  ties in this version); validated against brglm2, logistf and coxphf (the
  last through disjoint time windows, since it has no strata); definitions,
  refusals and the Monte Carlo in section 1, Phase 4.
* 2026-11-05 (author's decision, before the review rounds): an exact
  linear-programming separation check for mnl() and rologit() (lpSolveAPI in
  Imports; exported `separation()`), two starts for penalized fits, and the
  removal of rologit()'s never-ranked-above refusal; section 1, Phase 4b.
* 2026-11-06 (author's decision, during the stress battery's heavy part and
  before the review rounds):
  - the maximum-likelihood likelihood-ratio forms compare suprema under
    separation (fixing term_test()'s statistic from where the optimizer
    stopped);
  - boot_lr_test() for two rologit() fits;
  - priors for mnl() and rologit() (exported `coef_prior()`).

  Section 1, Phase 4c.
* 2026-11-06 (author, on the round-5 review and the stress rerun: "Agreed
  with all decisions/suggestions"):
  - the saddle step;
  - the zero-comparison correction of the separation type and suprema;
  - the boundary step one parameter at a time;
  - bootstrap draws from a separated restricted model's limit;
  - single-choice separated replicates kept;
  - weighted prior scales;
  - named prior scales of absent coefficients ignored;
  - `remnl()`'s base dependence documented;
  - Rainey's design in `firth_mc.R`;
  - ties, constants and a binary attribute in the ranking-bootstrap size
    check;
  - the article to gain a separation-and-priors subsection and the
    captivity-test interpretation;
  - its reproducibility through a tagged public commit (changed on
    2026-11-06, below: supplementary material).

  Section 1, Phase 4d.
* 2026-11-06 (author, on the article's audit and the replications):
  - the likelihood-ratio test of the logit against a constant captive share
    (section 1, Phase 4e);
  - the article's reproducibility: the scripts and run logs behind its
    simulation tables as supplementary material, for now (not a tagged
    public commit);
  - the review rounds 6 and 7 first, the stress battery's rerun after them;
  - terminology: the inflated multinomial logit is never called "the
    inflated logit" for short, in the documentation, the article or
    messages. Inflated binary logits exist, and "inflated logit" usually
    means one of them.

## Candidates for later releases (author, 2026-11-06 and 2026-11-07; not for this release)

A running list of candidates, not commitments. Each is reviewed before any
work starts, and some will not be taken up.

### Binary relatives of the inflated multinomial logit

* **The zero-inflated binary logit**: `imnl()` with two alternatives, the
  inflated category 0 and a logit split,

      P(y = 1) = Lambda(z'gamma) Lambda(x'beta).

  One likelihood under several names:
  - the logistic regression model with a cure fraction (Diop, Diop and
    Dupuy 2011, Electronic Journal of Statistics 5, 460-483,
    doi:10.1214/11-EJS616), the primary reference;
  - the one-trial case of the zero-inflated binomial of Hall (2000,
    Biometrics 56(4), 1030-1039), a model for bounded counts with random
    effects;
  - the split-population logit of Beger, DeMeritt, Hwang and Moore (2011,
    SSRN 1773594, unpublished), with the splitting coefficients' signs
    reversed (their eqs. 8-9 and 14-15, pp. 36-37);
  - partial observability with independent parts in its logit form (the
    probit form is Abowd and Farber 1982, Industrial and Labor Relations
    Review 35(3), 354-367);
  - the two-path "and" Boolean logit (Braumoeller 2003, Political Analysis
    11(3), 209-233, doi:10.1093/pan/mpg012).

  Inflating category 1 is the same model for 1 - y. With the same
  covariates in both parts the likelihood is symmetric in (beta, gamma), so
  the model is identified only up to exchanging the parts; Tomo, Eguchi and
  Yoneoka (2026, arXiv:2604.20322) prove this for the logit-logit model,
  give sufficient conditions for the existence of estimates (double
  separation) and a relabeling heuristic, and do no inference after
  relabeling. Diop et al. require a continuous covariate in one part only.
  With an intercept-only split the model is not identified. So an
  implementation refuses an intercept-only split and, by default, a shared
  design; checks for double separation; and starts from asymmetric points
  that include the exchanged pair.

  A probit split, Phi(z'delta) Lambda(x'beta), is the mirror image of the
  zero-inflated logit probit model of Pho (2024, Communications in
  Statistics - Theory and Methods 53(18), 6580-6599,
  doi:10.1080/03610926.2023.2248325), whose eq. (1) is
  Lambda(w'alpha) Phi(x'mu): a logit split and a probit outcome. Whether a
  difference in links identifies the parts under a shared design is open
  (Tomo et al. cover logit-logit only). Software: VGAM `zibinomial()`,
  gamlss.dist `ZIBI`, brms, GLMMadaptive and glmmTMB fit the likelihood;
  none has the identification checks.
* **The partial-observability bivariate probit** (Poirier 1980, Journal of
  Econometrics 12(2), 209-217, doi:10.1016/0304-4076(80)90007-X): y = y1 y2
  with correlated probit parts,

      P(y = 1) = Phi_2(z'gamma, x'beta; rho).

  rho = 0 gives Abowd and Farber's (1982) model, the probit analogue of the
  zero-inflated binary logit. Identification needs a covariate excluded from
  one equation (Poirier 1980). In political science: Przeworski and Vreeland
  (2002, Political Analysis 10(2), 101-112), on IMF agreements, and Nieman's
  strategic probit with partial observability (2015, Political Analysis
  23(3), 429-448). Rainey and Jackson (2018, Political Science Research and
  Methods 6(2), 381-391, doi:10.1017/psrm.2017.3) show by simulation that
  small misspecifications of functional form bias the estimates badly when
  a covariate enters both equations. Software: Stata `biprobit, partial`;
  GJRM (`Model = "BPO"`, and `"BPO0"` for rho = 0); BiProbitPartial was
  archived from CRAN in 2023.
* **Two-sided misclassification of the outcome** (Hausman, Abrevaya and
  Scott-Morton 1998, Journal of Econometrics 87(2), 239-269,
  doi:10.1016/S0304-4076(98)00015-3). With false-positive and false-negative
  rates alpha_0 and alpha_1,

      P(y = 1) = alpha_0 + (1 - alpha_0 - alpha_1) F(x'beta),

  constant rates with alpha_0 + alpha_1 < 1, identified through the
  functional form of F. With alpha_0 = 0 and alpha_1 = 1 - Lambda(z'gamma)
  it is the zero-inflated binary logit above. The multinomial version
  replaces the rates by a misreport matrix (Poterba and Summers 1995, Review
  of Economics and Statistics 77(2), 207-216, with rates from reinterview
  data; Molinari 2008, Journal of Econometrics 144(1), 81-117, partial
  identification). `imnl()` is already the case in which a true j is
  reported as k with probability 1 - phi(z); a free matrix is identified
  only by functional form (compare Gate C1 of Phase 3), so validation data
  or known rates would be needed. In political science: over-reporting a
  vote for the winner (Atkeson 1999, Political Behavior 21(3), 197-215);
  Katz and Katz (2010, American Journal of Political Science 54(3),
  815-835) correct binary turnout misreports with auxiliary information.
  Software: logistic4p (binary, constant rates).

### Multinomial models

* **A three-level nested logit**, extending nlogit()'s two levels (Train
  2009, ch. 4; Stata's nlogit fits more than two levels).
* **The latent-class rank-ordered logit** of Fok, Paap and van Dijk (2012,
  Journal of Applied Econometrics 27, 831-846, doi:10.1002/jae.1223). Latent
  segments differ in how many ranks a respondent reports reliably; it nests
  rologit() and allows a likelihood-ratio test of the segments.
* **The multinomial logit with a covariate-dependent scale** (heterogeneous
  choice): the utilities are scaled by a function of respondent covariates,
  so that choices are more or less predictable from one respondent to the
  next. Swait and Louviere (1993, Journal of Marketing Research 30(3),
  305-314) on the scale parameter; DeShazo and Fermo (2002, Journal of
  Environmental Economics and Management 44(1), 123-143) parametrize it.
  Political science has the tradition of Alvarez and Brehm (1995, American
  Journal of Political Science 39(4), 1055-1082) for binary and ordered
  choices, and Mauerer and Tutz (2023, Statistical Methods & Applications
  32(1), 129-148, doi:10.1007/s10260-022-00642-5) apply a heterogeneous
  multinomial logit to vote choice (GLES 2017; replication code, no
  package). It is the respondent-level counterpart of hev(). The scale
  equation has no constant, and scale and coefficient heterogeneity are
  confounded when a covariate enters both.
* **The fractional multinomial logit for shares**: the MNL quasi-likelihood
  with shares in place of 0/1 choices (Papke and Wooldridge 1996, Journal of
  Applied Econometrics 11(6), 619-632, for one fraction; Mullahy 2015,
  Journal of Econometric Methods 4(1), 71-100, for several), e.g. party vote
  shares by district, the data that Katz and King (1999, American Political
  Science Review 93(1), 15-32) and Tomz, Tucker and Wittenberg (2002,
  Political Analysis 10(1), 66-83) model through log-ratios. Martins and
  Murteira (2026, "Measuring electoral persistence in multiparty systems: a
  multivariate fractional regression approach", CeBER Working Paper
  2026-05) fit the dogit, that is imnl(), to vote shares, so a fractional
  imnl() comes with it. Only sandwich standard errors are valid, and the
  likelihood-based tools (LR tests, vuong(), AIC, the IIA tests) refuse
  these fits.
* **The conditional fixed-effects multinomial logit for panels** (Chamberlain
  1980, Review of Economic Studies 47(1), 225-238, doi:10.2307/2297110;
  Stata `femlogit`, Pforr 2014, Stata Journal 14(4), 847-862, with an
  application to British election panel data). It conditions the unit
  effects out, avoiding the incidental-parameters bias of unit dummies in
  short panels. The likelihood conditions on each unit's choice counts and
  so needs a recursion over orderings; it yields odds ratios but not
  probabilities, and units that never change and time-invariant covariates
  drop out. No R implementation found (2026-11-07).

### Features

* **Counterfactual choice sets, elasticities and diversion ratios, with
  intervals.** The question of what would have happened had a party not run
  (Alvarez and Nagler 1995, American Journal of Political Science 39(3),
  714-744; Lacy and Burden 1999, American Journal of Political Science
  43(1), 233-255), with simulation intervals (King, Tomz and Wittenberg
  2000, American Journal of Political Science 44(2), 347-361), under every
  family, including those that relax IIA (nlogit, hev, mnprobit, imnl).
  predict() on long data already drops an alternative whose rows are left
  out; intervals, an availability argument for wide data, and elasticities
  are missing. For imnl() the counterfactual has to say what respondents
  captive to a removed alternative do.
* **Survey designs and multiple imputation.** Design-based standard errors
  (Lumley 2004, Journal of Statistical Software 9(8)) for weighted election
  studies: polychoice's weights are frequency weights, and non-integer
  survey weights are treated as frequency weights. svyVGAM covers the
  multinomial logit with case-specific covariates only. Probability weights
  change what the log-likelihood means, so LR tests, vuong() and AIC refuse
  them or use design-based versions; WESML for choice-based samples (Manski
  and Lerman 1977, Econometrica 45(8), 1977-1988) comes with them. Multiple
  imputation (Amelia, Honaker, King and Blackwell 2011, Journal of
  Statistical Software 45(7); mice): coefficients pooled by Rubin's rules
  through tidy(), and pooled quantities of interest.

## 8. References

Abaluck, J. and Adams-Prassl, A. (2021). What do consumers consider before they
choose? Identification from asymmetric demand responses. Quarterly Journal of
Economics, 136, 1611-1663. doi:10.1093/qje/qjab008
Albert, A. and Anderson, J.A. (1984). On the existence of maximum likelihood
estimates in logistic regression models. Biometrika, 71, 1-10.
doi:10.1093/biomet/71.1.1
Allison, P.D. and Christakis, N.A. (1994). Logit models for sets of ranked
items. Sociological Methodology, 24, 199-228. doi:10.2307/270983
Bagozzi, B.E. and Marchetti, K. (2017). Distinguishing occasional abstention
from routine indifference in models of vote choice. Political Science Research
and Methods, 5(2), 277-294.
Barseghyan, L., Coughlin, M., Molinari, F. and Teitelbaum, J.C. (2021).
Heterogeneous choice sets and preferences. Econometrica, 89, 2015-2048.
doi:10.3982/ECTA17448
Beggs, S., Cardell, S. and Hausman, J. (1981). Assessing the potential demand
for electric cars. Journal of Econometrics, 17, 1-19.
doi:10.1016/0304-4076(81)90056-7
Ben-Akiva, M. and Boccara, B. (1995). Discrete choice models with latent choice
sets. International Journal of Research in Marketing, 12, 9-24.
doi:10.1016/0167-8116(95)00002-J
Bhat, C.R. (1995). A heteroscedastic extreme value model of intercity travel
mode choice. Transportation Research B, 29, 471-483.
Bull, S.B., Mak, C. and Greenwood, C.M.T. (2002). A modified score function
estimator for multinomial logistic regression in small samples. Computational
Statistics & Data Analysis, 39, 57-74. doi:10.1016/S0167-9473(01)00048-2
Chapman, R.G. and Staelin, R. (1982). Exploiting rank ordered choice set data
within the stochastic utility model. Journal of Marketing Research, 19,
288-301. doi:10.1177/002224378201900302
Cheng, S. and Long, J.S. (2007). Testing for IIA in the multinomial logit
model. Sociological Methods & Research, 35, 583-600.
Cook, S.J., Niehaus, J. and Zuhlke, S. (2018). A warning on separation in
multinomial logistic models. Research & Politics, 5(2).
doi:10.1177/2053168018769510
Crawford, I., Griffith, R. and Iaria, A. (2021). A survey of preference
estimation with unobserved choice set heterogeneity. Journal of Econometrics,
222, 4-43. doi:10.1016/j.jeconom.2020.07.024
Firth, D. (1993). Bias reduction of maximum likelihood estimates. Biometrika,
80, 27-38. doi:10.1093/biomet/80.1.27
Fok, D. and Paap, R. (2025). New misspecification tests for multinomial logit
models. Journal of Choice Modelling, 54, 100531.
Fry, T.R.L. and Harris, M.N. (1996). A Monte Carlo study of tests for the
independence of irrelevant alternatives property. Transportation Research B,
30, 19-30.
Gaudry, M.J.I. and Dagenais, M.G. (1979). The dogit model. Transportation
Research B, 13, 105-111. doi:10.1016/0191-2615(79)90028-6
Hausman, J.A. and McFadden, D. (1984). Specification tests for the multinomial
logit model. Econometrica, 52, 1219-1240.
Hausman, J.A. and Ruud, P.A. (1987). Specifying and testing econometric models
for rank-ordered data. Journal of Econometrics, 34, 83-104.
doi:10.1016/0304-4076(87)90068-6
Heinze, G. and Dunkler, D. (2008). Avoiding infinite estimates of
time-dependent effects in small-sample survival studies. Statistics in
Medicine, 27, 6455-6469. doi:10.1002/sim.3418
Heinze, G. and Schemper, M. (2001). A solution to the problem of monotone
likelihood in Cox regression. Biometrics, 57, 114-119.
Heinze, G. and Schemper, M. (2002). A solution to the problem of separation in
logistic regression. Statistics in Medicine, 21, 2409-2419.
Kosmidis, I. and Firth, D. (2011). Multinomial logit bias reduction via the
Poisson log-linear model. Biometrika, 98, 755-759. doi:10.1093/biomet/asr026
Long, J.S. and Freese, J. (2014). Regression Models for Categorical Dependent
Variables Using Stata, 3rd ed. Stata Press.
Manski, C.F. (1977). The structure of random utility models. Theory and
Decision, 8, 229-254. doi:10.1007/BF00133443
McFadden, D. (1974). Conditional logit analysis of qualitative choice
behavior. In Zarembka (ed.), Frontiers in Econometrics.
McFadden, D. (1987). Regression-based specification tests for the multinomial
logit model. Journal of Econometrics, 34, 63-82.
Oscarsson, H. and Rosema, M. (2019). Consideration set models of electoral
choice: Theory, method, and application. Electoral Studies, 57, 256-262.
doi:10.1016/j.electstud.2018.08.003
Paap, R., van Nierop, E., van Heerde, H.J., Wedel, M., Franses, P.H. and
Alsem, K.J. (2005). Consideration sets, intentions and the inclusion of
"don't know" in a two-stage model for voter choice. International Journal of
Forecasting, 21, 53-71. doi:10.1016/j.ijforecast.2004.02.004
Puhr, R., Heinze, G., Nold, M., Lusa, L. and Geroldinger, A. (2017). Firth's
logistic regression with rare events: accurate effect estimates AND
predictions? Statistics in Medicine, 36, 2302-2317. doi:10.1002/sim.7273
Sablica, L., Hornik, K. and Rusch, T. (2026). The geometry of existence and
uniqueness of maximum likelihood estimation in categorical response models.
arXiv:2610.08100.
Small, K.A. and Hsiao, C. (1985). Multinomial logit specification tests.
International Economic Review, 26, 619-627.
Speckman, P.L., Lee, J. and Sun, D. (2009). Existence of the MLE and
propriety of posteriors for a general multinomial choice model. Statistica
Sinica, 19, 731-748.
Stiemke, E. (1915). Über positive Lösungen homogener linearer Gleichungen.
Mathematische Annalen, 76, 340-342.
Swait, J. and Ben-Akiva, M. (1987a). Incorporating random constraints in
discrete models of choice set generation. Transportation Research B, 21,
91-102. doi:10.1016/0191-2615(87)90009-9
Swait, J. and Ben-Akiva, M. (1987b). Empirical test of a constrained choice
discrete model: mode choice in Sao Paulo, Brazil. Transportation Research B,
21, 103-115. doi:10.1016/0191-2615(87)90010-5
Train, K.E. (2009). Discrete Choice Methods with Simulation, 2nd ed. CUP.
Tse, Y.K. (1987). A diagnostic test for the multinomial logit model. Journal
of Business and Economic Statistics, 5, 283-286.
doi:10.1080/07350015.1987.10509587
Varin, C., Reid, N. and Firth, D. (2011). An overview of composite likelihood
methods. Statistica Sinica, 21, 5-42.
Vuong, Q.H. (1989). Likelihood ratio tests for model selection and non-nested
hypotheses. Econometrica, 57, 307-333.
Weesie, J. (1999). Seemingly unrelated estimation and the cluster-adjusted
sandwich estimator. Stata Technical Bulletin, 52, 34-47.
