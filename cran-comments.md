## First submission notes

This is the first submission of `polychoice` (0.1.0).

<!-- DRAFT written 2026-10-18 while the release is deferred. Every line marked
     [FILL] must be replaced by a result from the final tarball before the
     package is submitted; the rest must be re-read against that tarball.
     verification/release_checklist.md lists the steps. -->

## Test environments

* local Windows 11, R 4.6.1 (R CMD check --as-cran, development build
  0e818674026d of 2026-11-06, git 3ff0bdf): 0 errors, 0 warnings, 4 notes
  (the three below, and the `remnl` example at 5.1 s).
  [FILL: the same check on the final 0.1.0 tarball]
* [FILL: GitHub Actions (r-lib/actions): ubuntu-latest (R release, devel,
  oldrel-1), macOS-latest (release), windows-latest (release)]
* [FILL: win-builder R-devel: version, status, examples / tests / vignette
  times]
* [FILL: win-builder R-release]
* [FILL: sanitizer workflow (gcc ASAN/UBSAN) on the C++ likelihoods, as for
  `underdisp`]

## R CMD check results

[FILL: final counts] (development build 0e818674026d: 0 errors | 0 warnings | 4 notes)

* CRAN incoming feasibility: "New submission". The words flagged as possibly
  misspelled in the Description are author surnames (Albert, Allison,
  Anderson, Bagozzi, Beggs, Ben-Akiva, Bofin, Cardell, Christakis, Cramer,
  Dagenais, Fagerland, Firth, Fok, Gaudry, Gelman, Hausman, Hosmer, Hsiao,
  Jakulin, Liang, Marchetti, McFadden, Paap, Pittau, Rainey, Ridder, Self,
  Su, Swait, Tse, Vuong) and standard terms of the field (dogit, GMM,
  heteroskedastic, IIA, overidentification) [FILL: the release's list].
  On the development build this note also listed the development version
  number and the GitHub URLs, which did not exist yet; both are gone in the
  release [FILL: confirm].
* "Files 'README.md' or 'NEWS.md' cannot be checked without 'pandoc' being
  installed" and "Skipping checking math rendering: package 'V8' unavailable"
  are properties of the local machine [FILL: confirm absent on win-builder].

## Notes for the reviewer

* The likelihoods of the seven estimators and their gradients are compiled
  'C++' (Rcpp), built with R's default C++ standard (there is no
  src/Makevars); no system requirements.
* Examples run in [FILL] s in total ([FILL] s with \donttest{}); the examples in
  \donttest{} are the parametric-bootstrap calibrations (boot_lr_test,
  boot_pvalue), the likelihood-ratio form of captivity_test, the
  published-application fits on `pa2006` and `pres2004`, the inflated fit on
  `fluvax`, and the IIA battery on `bes1987`.
  The CRAN-visible test suite runs in about [FILL] s (209 s on the
  development build 0e818674026d, before the release's trimming; the
  vignettes took 247 s and the \donttest{} examples 157 s there); the full suite
  (NOT_CRAN = "true", with the published-application replications and the
  bootstrap calibrations) runs locally. The five vignettes build in about
  [FILL] s. [FILL: the final tarball's times; the total must stay well below
  ten minutes on win-builder.]
* Nine data sets are bundled from Harvard Dataverse replication archives,
  all released under CC0 1.0: the two applications of Bagozzi and Marchetti
  2017 (doi:10.7910/DVN/3JR2YL; data-raw/make_data.R) and seven from other
  published studies (data-raw/make_examples.R); both scripts check the md5
  sums of the source files. Two of the seven are extracts of British Election
  Study surveys that the UK Data Service distributes under CC BY 4.0, which
  their help pages cite. All strings are ASCII (the state name
  Baden-Wuerttemberg in `cabinets` is transliterated). [FILL: installed data
  size.]
