# Release checklist (polychoice 0.1.0)

Written 2026-10-18. The release is deferred by the author (GitHub and CRAN);
nothing below starts until the author says so. Steps are in order; a step that
changes any file under R/, src/, tests/, vignettes/, man/, data/, inst/,
DESCRIPTION or NAMESPACE after the final frozen run means a new frozen run of
every Monte Carlo (standing rule) before the tarball is built.

## Done (2026-10-18)

* Name: `polychoice` is not on CRAN, not in the CRAN archive (`polytomous`
  is archived and cannot be reused), and `bagozzib/polychoice` does not exist
  on GitHub (checked 2026-10-18).
* DESCRIPTION: title of 64 characters since 2026-11-06 (the limit is 65); every
  method in the Description has a DOI, each resolved on Crossref to the right
  paper (Bagozzi and Marchetti 2017, Hausman and McFadden 1984, Small and
  Hsiao 1985, Fok and Paap 2025, McFadden 1987, Cramer and Ridder 1991,
  Fagerland, Hosmer and Bofin 2008, Vuong 1989); 'C++' quoted.
* inst/CITATION title matches the new title.
* cran-comments.md drafted with [FILL] markers.
* Words the CRAN spell check will flag in the Description (hunspell en_US):
  surnames plus GMM, heteroskedastic, IIA, overidentification (listed in
  cran-comments.md).

## Done (2026-11-05)

* The author confirmed that polychoice goes to CRAN before the R Journal
  article is submitted.
* GitHub Actions workflows written in `.github/workflows/` (not yet pushed;
  `.github` is .Rbuildignored): `R-CMD-check.yaml` (the r-lib matrix with
  NOT_CRAN set, so the full suite runs there), `sanitizer.yaml` (clang-UBSAN
  and clang-ASAN r-hub containers; every vignette rendered directly with
  tools::buildVignettes(), all output grepped for sanitizer lines, as for
  underdisp), and the generic `rhub.yaml`.
* Check-time budget, measured from the suite log of build 014e91568fbd
  (verification/runs/014e91568fbd/tests.log, a quiet machine): the tests
  that run on CRAN (no skip_on_cran()) take 127 s of the suite's 621 s.
  Four of them take 82 s: the HEV kernel derivative and integrate() sweep
  (test-families.R, 36 s; the five scale regimes and two numDeriv Hessians),
  boot_lr_test() nesting (test-review2.R, 25 s), imnl() with two inflated
  alternatives (test-captive.R, 14 s), and the weighted fixed-effects
  bootstrap (test-review2.R, 7 s). The October checks under load: examples
  30 s, \donttest examples 82-155 s, tests 154-188 s, vignettes 179-219 s
  (beyond.Rmd 109-130 s). To do on a quiet machine once the stress rerun has
  finished: R CMD check --as-cran with per-vignette and per-example times,
  then trim to a total well under ten minutes (candidates: the four tests
  above, by moving their slow parts under skip_on_cran() while a fast core
  stays on CRAN; beyond.Rmd's HEV, probit and jackknife chunks).

## Done (2026-11-10)

* The record of the search's ends and the count of local maxima (DESIGN.md
  Phase 4f, the author's decision of 2026-11-10), with the documentation
  fixes that followed the frozen run of 9446ac9410f9 (the ledger's
  reconciliation entries of 2026-11-09).
* Step 4 below: version 0.1.0 in DESCRIPTION, NEWS.md "# polychoice 0.1.0"
  with "Initial CRAN release."
* So the final frozen run (step 5) and the stress rerun are on the release
  commit itself (the standing rule above).
* The stress rerun of that commit (build 28b870609bea) found ends on a ridge
  of non-identification counted as separate maxima; the author stopped its
  runs for the grouping amendment of Phase 4f and the `singular` flag (the
  ledger's entry of 2026-11-10), so the stress rerun and the final frozen run
  are on the commit that follows.

## Done (2026-11-12)

* The stress rerun and the final frozen run (step 5 below) on the release
  commit 8a8de05, build beb30c4880b2 (the ledger's entry of 2026-11-12):
  every saved result of the frozen run identical to 9446ac9410f9's apart
  from the several-maxima counts; the article's pres2004 bootstrap and
  timings rerun on the build.
* The clean-up for the public repository (the author's decisions of
  2026-11-10 and 11-11): files outside the built package only, so the build
  id is unchanged.

## Release sequence (the author's order, 2026-11-05)

1. The stress rerun on 014e91568fbd: compare with d2c679895bf9, report to
   the author, address what it raises, commit its records.
1b. Firth's penalty for mnl() and rologit() (author's decision 2026-11-05,
   before the reviews): definitions fixed in DESIGN.md Phase 4 (commit
   d23a7a3); implementation, oracle tests, documentation and the article
   paragraph; its Monte Carlo runs with the rest at step 3. Done 2026-11-05:
   0588634 (code, tests, the `defection` data, docs; full suite on build
   82c388c1537e) and d427a86 (verification/firth_mc.R, in the runner).
2. Three adversarial review rounds of the package, its documentation, and
   the R Journal article (internal_reviews/round5 to round7, outside the
   repository): ten simulated referees in five reports per round (the nine
   roles of rounds 1-4, updated, plus an R Journal referee); fixes after
   each round; each later round first verifies the previous round's fixes.
   Then check in with the author before step 3.
3. The final frozen run (tests, R CMD check --as-cran with the check-time
   trimming, every Monte Carlo) on the release build.
4. GitHub (with the author's go-ahead and decision on what the public
   repository holds), win-builder, cran-comments; the author submits.
5. After CRAN acceptance: the article's [FILL] numbers, rjtools checks,
   the motivation letter; the author submits to the R Journal. Before the
   [FILL] numbers: re-knit `bagozzi.Rmd` against the release build, diff
   every chunk's output against the previous render, and re-read each
   sentence that quotes an output that changed; run
   `scripts/pres2004_boot_lr.R` on the release build; name the release
   build in the captions of the size and Monte Carlo tables; assemble the
   supplementary material (the `verification/` scripts and run logs behind
   Tables size and mc, the article's scripts; author's decision of
   2026-11-06, listed in `Software Article/rjournal/supplementary/README.md`).

## To do, in order

1. Author decisions: the name (polychoice), the title, and the go-ahead for a
   public GitHub repository.
2. GitHub: create `bagozzib/polychoice`, push `main` (the verification/ and
   DESIGN.md stay in the repository; both are .Rbuildignored). The internal
   notes (`polychoice/internal_notes.md`, `internal_reviews/`) stay outside the
   repository. The URLs in DESCRIPTION, inst/CITATION and
   R/polychoice-package.R then resolve.
3. GitHub Actions: R-CMD-check matrix (ubuntu release/devel/oldrel-1, macOS
   release, windows release) and the sanitizer workflow copied from
   `underdisp` (.github/workflows/sanitizer.yaml, with the control-branch
   trick); read the `.Rout` files for `runtime error:` lines, since a
   sanitizer `Status: OK` hides them (underdisp 0.1.2 lesson).
4. Version 0.0.0.9000 -> 0.1.0 in DESCRIPTION; NEWS.md header
   "# polychoice 0.1.0" with "Initial CRAN release."
5. Final frozen run (verification/run_verification.sh) on the release commit:
   tests, R CMD check --as-cran, every MC; reconcile every quoted number.
6. win-builder R-devel and R-release on the final tarball (plain FTP or the
   web form); check the total check time (the CRAN incoming limit is ten
   minutes; underdisp 0.1.0 was auto-rejected at 13). Render the vignettes
   with tools::buildVignettes() once directly, since R CMD build/check discard
   the vignette child's output on success.
7. Fill every [FILL] in cran-comments.md from those results; re-read the
   whole file against the tarball.
8. Submit through the CRAN web form (the author); then the confirmation
   e-mail link; then wait (pretest, then the newbies queue).
9. On acceptance: tag v0.1.0 on the submitted commit, GitHub release, README
   badges and install.packages() first, CRAN URL into DESCRIPTION URL and
   inst/CITATION at the next version.

## For 0.1.1

The author's decision of 2026-11-15: the items listed here on 2026-11-14 and
11-15 are done before the first release (DESIGN.md Phase 4g; ledger entry of
2026-11-15). Done: model names in full; Q2 and Q1 of fp_tests() in other
units; Q1 on a separated fit; the GMM step's silent exit; a `singular` flag
that does not depend on where the search stopped (a polish to stationarity,
not the rule on the scores); Wald tests on a singular fit; profile limits
that do not verify; the second start of a penalized fit; the saddle step's
second start (koyama2025_t4walter, mnprobit: -492.6971 again).

Traced, with no change to the package:

* mehrez2024_t4, mnprobit: alone on the machine on build 8dce1271fc30 it
  takes 1,572 s and ends where build beb30c4880b2 did (-422.3993, largest
  scaled gradient 23.4, not converged); it reached the 1,800 s cap only
  beside the other streams (probes/2026-11-15/). The stress chain's step
  cap is 3,600 s from the next run on, which also covers jost2025_a25's
  inflated fits (750 to 1,400 s alone).

Open:

* First for 0.1.1 (the author, 2026-11-16): tests on a clustered covariance
  with fewer clusters than coefficients. Q1's joint covariance then has the
  rank of the clusters, the statistic is computed on that column space with
  its rank as the degrees of freedom, and nothing warns: eifert2010_t4 (22
  clusters, 84 coefficients) 10.44 on 21 degrees of freedom, goldring2023_all5
  (16, 32) 90.65 on 14, kardasheva2013_t5m5 (19, 57) 28.33 on 18; on build
  8dce1271fc30 150.03, 189.24 and 22.55, changing with the covariates' units
  (probes/2026-11-16/q1_units_7.txt, q1_units_8.txt). Q2 on eifert2010_t4 is
  0.65 on 126 degrees of freedom. One rule for Q1, Q2, the joint
  Hausman-McFadden test and the Wald tests of term_test() and combine_test()
  (not yet examined in that case).
* fp_tests(), Q2 with gmm = "projected": in two small simulations of the
  adversarial review of 2026-11-15 it rejected 0.8 and 0 percent at the 5
  percent level (two-step: 5.2 and 11.7 percent); probes/2026-11-15/review/
  r4_mc.R and r5_mc2.R. The size study of the frozen run covers the default.
* A fit of an unidentified model that stops with a reciprocal condition
  number above 1e-6 is not polished and not flagged (DESIGN.md Phase 4g,
  item 4).
* `start_logliks` is not updated after the polish to stationarity (a
  difference of the order of 1e-10).
* The gate of fp_tests() on an ill-conditioned information depends on a
  covariate's origin (a calendar year); the help page says so.
* Robust standard errors of diverging coefficients: on one separated
  rank-ordered fit with tied rankings (the ledger's entry of 2026-11-15) they
  change from 1.60 to 2.01 with a change of the
  estimates in the eighth digit (analytic: 4,079). The fit warns of the
  separation; the sandwich there could be reported as NA.
