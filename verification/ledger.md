# polychoice verification ledger: the results the documentation cites

The vignettes and help pages cite a verification ledger for four results. This file summarizes the
entries of the author's working ledger on those four, with their dates, and then records the frozen run
of the release build, which is the source of the numbers the documentation quotes. The working ledger
itself (every check, review round and stress run since 2026-10-13) is kept by the author, as are the
one-off probe scripts named below; their results are given here. The permanent checks are the package
tests under `tests/testthat/` and the scripts in `verification/`, which
`verification/run_verification.sh` runs on a frozen copy of the package and records in
`verification/runs/<build id>/`.

## 1. The HEV fit of the travel-mode data and mlogit's value (2026-10-13, 10-14)

* `hev(chosen ~ wait + gcost | income + size)` on AER::TravelMode converges to a log-likelihood of
  -169.969869 with the car scale going to zero. The profile in the log car scale, the other parameters
  at their estimates, is -170.127 at 0.0719, -169.979 at 0.018, -169.970 at 0.0025 and flat below 3e-4:
  the maximum is at the boundary, and the boundary flag fires.
* mlogit 2.0.0 (`heterosc = TRUE`) reports -172.4258 at scales 0.83021 / 0.56441 / 0.07192. The
  log-likelihood at those estimates is -176.1917, both by the package's kernel and by R's
  `integrate()` (rel.tol 1e-12, case by case): the reported value is 3.77 above the log-likelihood at
  mlogit's own solution, the error of its Gauss-Laguerre quadrature at a scale ratio of 14.
* The package's kernel (src/hev_ll.cpp) uses the trapezoid rule where every scale ratio is at most 20
  and adaptive Gauss-Kronrod 7-15 panels beyond, with breakpoints around each competitor's step. Against
  `integrate()` on the test design (n = 200, J = 4), with one log scale at -3 to -12 and +4, +8 (ratios
  20 to 1.6e5), every probability agrees to 7e-12 and rows sum to one within 1.1e-11. At the travel-mode
  boundary fit the analytic gradient agrees with numDeriv to 6e-9 in standardized units.
* tests/testthat/test-families.R therefore checks the kernel against `integrate()` at mlogit's point,
  the boundary optimum and the profile, not against mlogit's numbers.

## 2. remnl() against lme4::glmer() (2026-10-13, 10-14)

* At J = 2 the model is the binomial GLMM with a random intercept. On 150 units of 6 cases (sd 1.1),
  `glmer(nAGQ = 25)` gives -518.027156; `remnl()` gives -518.027160 with 20 nodes and -518.027156 with
  40 (a difference of 5e-9).
* verification/remnl_adaptive_check.R (run at every frozen build): J = 2, 150 units, T in {6, 20, 50},
  sd in {1, 2}. With the default 12 adaptive nodes the log-likelihood is within 1.2e-3 of glmer's
  (worst: T = 6, sd 2) and the standard deviation within 0.006 percent in every cell; with 40 nodes
  within 2e-6.
* Recovery at J = 3 (400 units of 5 cases, sd 0.8 / 1.2): all eight parameters within 1.4 standard
  errors of the truth.

## 3. The gate of the split-panel jackknife (2026-10-14)

* **Design.** The gate compares the time split's distance between half-panel slopes with the noise of
  five placebo splits (random halves of every unit's cases, from fixed seeds), and tests for a residual
  trend at alpha = 0.01. Simulations: 150 units, J = 3, T in {6, 10, 20}; 200 stationary, 40
  regime-change and 40 omitted-trend panels per T.
* **Selection rule, written at 13:15 before any simulation of this design was run.** kappa is the
  smallest value on {1.25, 1.5, 1.75, 2, 2.5, 3, 3.5, 4} whose refusal rate on stationary panels,
  counting every reason the gate refuses, is at most 0.05 at every T (verification/jackknife_gate_mc.R,
  seed base 5000). verification/jackknife_gate_holdout.R (seed base 9000) confirms it when its
  stationary refusal rate is at most 0.081 (0.05 plus two Monte Carlo standard errors) at every T.
* **Selection run (16:27).** Stationary refusal at T = 6 / 10 / 20:

  | kappa | 1.25 | 1.5 | 1.75 | 2 | 2.5 | 3 | 3.5 | 4 |
  |---|---|---|---|---|---|---|---|---|
  | T = 6 | 0.285 | 0.200 | 0.120 | 0.085 | 0.025 | 0.015 | 0.010 | 0.010 |
  | T = 10 | 0.285 | 0.135 | 0.090 | 0.055 | 0.030 | 0.030 | 0.020 | 0.020 |
  | T = 20 | 0.210 | 0.130 | 0.055 | 0.040 | 0.010 | 0.000 | 0.000 | 0.000 |

  By the rule, kappa = 2.5. At that value a regime change is refused in 0.925 / 1.000 / 1.000 of the
  panels and an omitted trend in 0.525 / 1.000 / 1.000.
* **Both estimators on the accepted stationary panels** (b:x1, b:x2, c:x1, c:x2; truth 0.8, -0.5, -0.6,
  0.7; coverage of 95 percent Wald intervals with the ML standard errors):

  | T | | signed bias | RMSE | coverage |
  |---|---|---|---|---|
  | 6 | ML | +0.240 -0.137 -0.189 +0.230 | 0.302 0.209 0.258 0.278 | 0.63 0.77 0.74 0.68 |
  | 6 | jackknife | -0.117 +0.094 +0.073 -0.089 | 0.234 0.190 0.217 0.210 | 0.78 0.83 0.81 0.81 |
  | 10 | ML | +0.109 -0.072 -0.120 +0.132 | 0.153 0.121 0.163 0.174 | 0.76 0.86 0.76 0.74 |
  | 10 | jackknife | -0.061 +0.045 +0.030 -0.035 | 0.114 0.105 0.119 0.120 | 0.88 0.94 0.93 0.88 |
  | 20 | ML | +0.051 -0.036 -0.049 +0.057 | 0.081 0.071 0.083 0.088 | 0.88 0.91 0.88 0.87 |
  | 20 | jackknife | -0.011 +0.004 +0.010 -0.008 | 0.062 0.059 0.063 0.066 | 0.94 0.97 0.95 0.95 |

  The jackknife over-corrects at T = 6 and 10.
* An earlier design of the gate (deterministic placebo splits, kappa = 4 by the same rule) was replaced
  by this one the same day, before release. The selection run's own log and results are
  verification/jackknife_gate_selection_run.log and .rds; in the log, the path of the temporary library
  the run used is shortened.

## 4. The release build's frozen run, with the size study of the IIA tests

Frozen run of polychoice 0.1.0, 2026-11-15 to 11-16, build f46ebb007be6, run folder
verification/runs/f46ebb007be6. The build id is the fingerprint of the package's files
(fingerprint.txt in the run folder); status.txt has one header line per segment of the run.

* **Tests and check.** 1,454 expectations, 0 failed, warnings, errors or skipped. R CMD check
  --as-cran: 0 errors, 0 warnings, 3 NOTEs (a new submission; pandoc and V8 missing locally). The
  precomputed vignettes knitted again on the build reproduce the shipped ones.
* **Size study of the IIA tests** (size_check_iia.R, q1_covariance_check.R; the multinomial logit, J = 3, four
  covariates, n = 1,000, 300 replications). Rejection at 0.05 / 0.10: Hausman-McFadden difference form
  0.077 / 0.163, joint form 0.020 / 0.047; Small-Hsiao with one split 0.040 / 0.047, ten-split median
  0; Q2 0.063 / 0.063 / 0.060. Q1 at 0.05: joint 0.057, difference 0.043; at 1,000 replications and
  0.01 / 0.05 / 0.10, joint 0.015 / 0.062 / 0.112 and difference 0.033 / 0.042 / 0.049. The
  Hausman-McFadden difference form is not positive semi-definite in 0.79 of the replications, and
  0.436 of Q1's difference statistics are negative.
* **Jackknife gate.** jackknife_gate_mc.R reproduces the selection run (kappa 2.5; stationary refusal
  0.025 / 0.030 / 0.010, regime change 0.925 / 1 / 1, omitted trend 0.525 / 1 / 1). On the holdout's
  fresh seeds the stationary refusal is at most 0.025.
* **remnl_adaptive_check.** Within 1.9e-06 of glmer's log-likelihood at 40 nodes in every design and
  within 1.2e-03 at 12. table2_check: Replication II's Vuong statistic -1.680, two-sided p 0.093.
* **Recovery and coverage.** mnl coverage 0.945-0.965 and imnl 0.920-0.970, biases within 2.1 se.
  rologit coverage 0.920-0.983. mnprobit: at n = 1500 coverage 0.920-0.965; at n = 3000 0.905-0.970
  (c:x2 below the band) with b:x1 biased by -0.010 (-3.4 se). The two cells share their seeds; on 200
  fresh seeds at n = 3000 (a one-off probe on an earlier build with the same estimator) every bias is
  within 1.42 Monte Carlo se and coverage is 0.940-0.985.
* **Captivity.** captivity_test_mc, score form: size 0.028-0.063 one-sided and 0.057-0.073 joint at
  0.05, power 0.466 / 0.942. Likelihood-ratio form: size 0.036-0.064 in the four regular designs,
  0.017-0.027 with near-zero fitted probabilities, 0.041-0.061 with an alternative chosen by about 3
  percent; power 0.498 / 0.958. captivity_lr_pa2006_boot: rejection at the mixture's 10, 5, 1 percent
  critical values 0.101 / 0.040 / 0.010; observed LR 20.951, bootstrap p 0.005.
* **captive_mc.** With the exclusion restriction: every criterion met at n = 5000 (coverage
  0.920-0.975); two misses at n = 2000 (coverage 0.915, a bias of -3.1 se). Without one: at n = 2000,
  17 of 200 fits at a boundary, the x1 outcome slopes of the others cover 0.798-0.814 with bias +0.31
  to +0.33, and the search reports several maxima in 87 fits; at n = 5000, 3 at a boundary, coverage
  0.893-0.898, bias +0.08 to +0.10, several maxima in 55.
* **Bootstrap tests.** boot_lr_size_mc: LR 0.047 / 0.117 and Vuong 0.040 / 0.090 at 0.05 / 0.10.
  rologit_boot_size_mc: 0.033-0.057 at 0.05 and 0.077-0.123 at 0.10 in the six designs; in the sparse
  design the chi-square reference rejects 0.088 / 0.180 and the bootstrap 0.044 / 0.102.
  vuong_size_mc: the raw statistic rejects 0.003 under the null; the LR against the 50:50 mixture
  0.082.
* **firth_mc.** Every pre-registered criterion holds: profile coverage 0.940-0.968, penalized LR size
  0.039-0.060, penalized bias at most 0.11 in absolute value. Rare covariate with 5 cases: Wald size
  0.084 / 0.060, penalized LR 0.064 / 0.055 (J = 2 / 3); with 40 cases, 0.057-0.081 and 0.040-0.047,
  power 0.86-0.94 and 0.80-0.91 at +2.
