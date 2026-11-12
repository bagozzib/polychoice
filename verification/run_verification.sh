#!/bin/bash
# Reruns the verification on one frozen copy of the package.
#
#   bash verification/run_verification.sh [--tests] [--check] [script.R ...]   (from the package root)
#
#   --tests   every tests/testthat file with NOT_CRAN set (verification/run_test_file.R)
#   --check   R CMD build (with vignettes) and R CMD check --as-cran of the frozen copy
#   script.R  verification scripts by file name; with no arguments at all,
#             the tests, the check, and every script below are run
#
# The frozen copy is git HEAD when the package is a git checkout (uncommitted
# changes are left out, with a warning recorded in status.txt) and the working
# tree otherwise. The md5 sums of the package's own files (R/, src/, data/,
# inst/, man/, tests/, vignettes/, DESCRIPTION, NAMESPACE), taken before
# anything is compiled, give the build id PC_BUILD, which every script prints
# in its log; the verification scripts are fingerprinted separately
# (scripts_fingerprint.txt), so a change to them does not change the build id.
# The copy is installed into a private library. Everything a run writes goes to
# verification/runs/<PC_BUILD>/: the fingerprints, the install log, status.txt
# (git commit, whether the tree had uncommitted changes, the replication
# settings passed through the environment, the pandoc version, every step's
# exit status), the test log (rewritten by each --tests run), the build and
# check logs, and the scripts' logs, results, and stderr. The scripts run from
# the frozen copy; their committed outputs there are removed before each runs,
# so a script that fails early cannot leave an old file to be copied.
# R CMD build needs pandoc for the vignettes: RSTUDIO_PANDOC is set to the
# RStudio copy when it is not set and that copy exists, and the check step stops
# when rmarkdown cannot find pandoc.
# Replication settings are environment variables (JK_REPS_NULL, JK_REPS_ALT,
# JK_SEED_BASE, VUONG_REPS, COV_REPS_MNL, COV_REPS_IMNL, IFE_REPS, BLR_REPS,
# BLR_B, Q1_REPS, MNP_REPS, RO_REPS, CT_REPS_SIZE, CT_REPS_POWER, CAP_REPS, FIRTH_REPS, FIRTH_WORKERS, RBS_REPS, RBS_B); the defaults are those the ledger reports.
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
cd "$ROOT" || exit 1
to_win() { if command -v cygpath > /dev/null 2>&1; then cygpath -m "$1"; else echo "$1"; fi; }
DO_TESTS=0; DO_CHECK=0; SCRIPTS=()
for a in "$@"; do case "$a" in --tests) DO_TESTS=1 ;; --check) DO_CHECK=1 ;; *) SCRIPTS+=("$a") ;; esac; done
if [ $# -eq 0 ]; then
  DO_TESTS=1; DO_CHECK=1
  SCRIPTS=(remnl_adaptive_check.R table2_check.R size_check_iia.R q1_covariance_check.R jackknife_gate_mc.R jackknife_gate_holdout.R vuong_size_mc.R
           recovery_coverage_mc.R mnprobit_recovery_mc.R boot_lr_size_mc.R rologit_mc.R captivity_test_mc.R captivity_lr_pa2006_boot.R captive_mc.R firth_mc.R rologit_boot_size_mc.R)
fi
WORK="$(mktemp -d)"; SRC="$WORK/polychoice"; mkdir -p "$SRC" "$WORK/lib"
DIRTY="no"
if git -C "$ROOT" rev-parse --git-dir > /dev/null 2>&1; then
  COMMIT="$(git -C "$ROOT" rev-parse HEAD)"
  if [ -n "$(git -C "$ROOT" status --porcelain --untracked-files=no)" ]; then
    DIRTY="yes (uncommitted changes are NOT in the frozen copy)"; echo "warning: uncommitted changes are not in the frozen copy (git HEAD $COMMIT)"
  fi
  git -C "$ROOT" archive HEAD | tar -x -C "$SRC"
else
  COMMIT="none (not a git checkout)"; cp -r "$ROOT/." "$SRC/"
fi
rm -rf "$SRC/verification/runs"; rm -f "$SRC"/src/*.o "$SRC"/src/*.dll "$SRC"/src/*.so
( cd "$SRC" && find R src data inst man tests vignettes DESCRIPTION NAMESPACE -type f 2>/dev/null | sort | xargs md5sum ) > "$WORK/fingerprint.txt"
( cd "$SRC" && find verification -type f | sort | xargs md5sum ) > "$WORK/scripts_fingerprint.txt"
PC_BUILD="$(md5sum "$WORK/fingerprint.txt" | cut -c1-12)"
R CMD INSTALL --no-multiarch -l "$(to_win "$WORK/lib")" "$(to_win "$SRC")" > "$WORK/install.log" 2>&1 || { echo "install failed: $WORK/install.log"; exit 1; }
OUT="$ROOT/verification/runs/$PC_BUILD"; mkdir -p "$OUT"
cp "$WORK/fingerprint.txt" "$WORK/scripts_fingerprint.txt" "$WORK/install.log" "$OUT/"
export PC_BUILD PC_COMMIT="$COMMIT" R_LIBS="$(to_win "$WORK/lib")" PC_PKG="$(to_win "$SRC")"
if [ -z "${RSTUDIO_PANDOC:-}" ] && [ -d "/c/Program Files/RStudio/resources/app/bin/quarto/bin/tools" ]; then
  export RSTUDIO_PANDOC="C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools"
fi
PANDOC="$(Rscript -e 'cat(if (rmarkdown::pandoc_available()) as.character(rmarkdown::pandoc_version()) else "none")' 2>/dev/null)"
STATUS="$OUT/status.txt"
{ echo "build $PC_BUILD  git $COMMIT  uncommitted changes: $DIRTY  started $(date '+%F %T')  R_LIBS=$R_LIBS  pandoc $PANDOC"
  env | grep -E '^(JK_|VUONG_|COV_|IFE_|BLR_|Q1_|MNP_|RO_|CT_|CAP_|FIRTH_|RBS_)' | sed 's/^/  env /'; } >> "$STATUS"
cd "$SRC" || exit 1
if [ $DO_TESTS -eq 1 ]; then
  : > "$OUT/tests.log"
  for f in tests/testthat/test-*.R; do
    Rscript verification/run_test_file.R "$(basename "$f")" >> "$OUT/tests.log" 2>&1
    echo "  tests $(basename "$f") exit $?" >> "$STATUS"
  done
  grep "total expectations" "$OUT/tests.log" | awk '{n += $5; f += $7; w += $9; e += $11; s += $13} END {printf "  tests: %d expectations, %d failed, %d warnings, %d errors, %d skipped\n", n, f, w, e, s}' >> "$STATUS"
fi
if [ $DO_CHECK -eq 1 ] && [ "$PANDOC" = "none" ]; then
  echo "  R CMD build/check NOT RUN: pandoc not available (set RSTUDIO_PANDOC)" >> "$STATUS"
elif [ $DO_CHECK -eq 1 ]; then
  ## the precomputed vignettes (vignettes/precompute.R) knitted again on this
  ## build and compared with the committed ones: a difference means they are
  ## stale (the author's decision of 2026-11-06 to precompute them)
  if ls vignettes/*.Rmd.orig > /dev/null 2>&1; then
    mkdir -p "$WORK/vig"
    Rscript vignettes/precompute.R "$(to_win "$WORK/vig")" > "$OUT/vignettes_precompute.log" 2>&1
    vdiff=0; for v in vignettes/*.Rmd.orig; do b="$(basename "${v%.orig}")"; diff -q --strip-trailing-cr "vignettes/$b" "$WORK/vig/$b" > /dev/null 2>&1 || { vdiff=1; diff --strip-trailing-cr "vignettes/$b" "$WORK/vig/$b" > "$OUT/vignette_diff_$b.txt" 2>&1; }; done
    echo "  precomputed vignettes $( [ $vdiff -eq 0 ] && echo 'current' || echo 'STALE (vignette_diff_*.txt)' )" >> "$STATUS"
  fi
  mkdir -p "$WORK/check"
  ( cd "$WORK/check" && R CMD build "$(to_win "$SRC")" > build.log 2>&1 ); bcode=$?
  echo "  R CMD build exit $bcode" >> "$STATUS"
  if [ $bcode -eq 0 ]; then
    ( cd "$WORK/check" && R CMD check --as-cran --timings polychoice_*.tar.gz > check.log 2>&1 )
    echo "  R CMD check exit $?" >> "$STATUS"
  fi
  cp "$WORK/check/build.log" "$OUT/r_cmd_build.log" 2>/dev/null
  cp "$WORK/check/polychoice.Rcheck/00check.log" "$OUT/r_cmd_check_00check.log" 2>/dev/null
  cp "$WORK/check/polychoice.Rcheck/tests/testthat.Rout" "$OUT/r_cmd_check_testthat.Rout" 2>/dev/null
  cp "$WORK/check/polychoice.Rcheck/polychoice-Ex.timings" "$OUT/r_cmd_check_Ex.timings" 2>/dev/null
fi
for s in "${SCRIPTS[@]}"; do
  name="${s%.R}"; start="$(date '+%F %T')"
  rm -f verification/"$name"*.log verification/"$name"*.rds verification/"$name"*.csv
  Rscript "verification/$s" > /dev/null 2> "$OUT/$name.err"
  code=$?
  echo "  $s exit $code  $start -> $(date '+%F %T')" >> "$STATUS"
  for f in verification/"$name"*.log verification/"$name"*.rds verification/"$name"*.csv; do [ -f "$f" ] && cp "$f" "$OUT/"; done
  [ -s "$OUT/$name.err" ] || rm -f "$OUT/$name.err"
done
echo "  finished $(date '+%F %T')" >> "$STATUS"
echo "done: $OUT"
