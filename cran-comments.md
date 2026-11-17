## First submission notes

This is the first submission of `polychoice` (0.1.0).

## Test environments

* local Windows 11, R 4.6.1: R CMD check --as-cran, 0 errors, 0 warnings,
  3 notes (the one below, and two that are properties of the local machine:
  pandoc and the 'V8' package are not installed).
* win-builder R-devel (2026-11-16 r90914 ucrt): 1 note (below); installation
  49 s, check 275 s (examples 26 s, tests 99 s).
* win-builder R-release (R 4.6.2): 1 note (below); installation 49 s, check
  239 s (examples 27 s, tests 100 s).
* GitHub Actions (r-lib/actions): ubuntu-latest (R release, devel,
  oldrel-1), macOS-latest (release), windows-latest (release), all passing
  with the full test suite (NOT_CRAN = "true").
* clang-ASAN and clang-UBSAN (the R-hub containers, on GitHub Actions) on
  the 'C++' likelihoods, with every vignette rendered: no sanitizer
  diagnostics in the vignette, build or check output.

## R CMD check results

0 errors | 0 warnings | 1 note

* CRAN incoming feasibility: "New submission". The words flagged as possibly
  misspelled in the Description are author surnames (Bagozzi, Beggs, Bofin,
  Cardell, Christakis, Cramer, Dagenais, Fagerland, Firth, Fok, Gaudry,
  Gelman, Hausman, Hosmer, Hsiao, Jakulin, Liang, Marchetti, Paap, Pittau,
  Rainey, Ridder, Su, Swait, Tse, Vuong) and standard terms of the field
  (dogit, GMM, heteroskedastic, IIA, overidentification, suprema).

## Notes for the reviewer

* The likelihoods of the seven estimators and their gradients are compiled
  'C++' (Rcpp), built with R's default C++ standard (there is no
  src/Makevars); no system requirements.
* The examples run in about 27 s on win-builder (about 90 s locally with
  \donttest{}); the examples in \donttest{} are the parametric-bootstrap
  calibrations (boot_lr_test, boot_pvalue), the likelihood-ratio form of
  captivity_test, the published-application fits on `pa2006` and `pres2004`,
  the inflated fit on `fluvax`, and the IIA battery on `bes1987`.
  The CRAN-visible test suite runs in about 100 s; the full suite
  (NOT_CRAN = "true", with the published-application replications and the
  bootstrap calibrations) runs locally and on GitHub Actions. The five
  vignettes are precomputed, so their code is not run at check time.
* Nine data sets are bundled from Harvard Dataverse replication archives,
  all released under CC0 1.0: the two applications of Bagozzi and Marchetti
  2017 (doi:10.7910/DVN/3JR2YL; data-raw/make_data.R) and seven from other
  published studies (data-raw/make_examples.R); both scripts check the md5
  sums of the source files. Two of the seven are extracts of British Election
  Study surveys that the UK Data Service distributes under CC BY 4.0, which
  their help pages cite. All strings are ASCII (the state name
  Baden-Wuerttemberg in `cabinets` is transliterated). The installed data
  directory is under 300 KB.
