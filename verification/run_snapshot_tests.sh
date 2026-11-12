#!/bin/bash
# The working-tree snapshot variant of run_verification.sh --tests (round 5, H7):
# the run records of builds 014e91568fbd, 1535b95732ac, 69c432153ae2,
# 82c388c1537e, 8d8adbd2a40f, bf4c913715e2 and c32e2c338c0f came from it (their
# status.txt says "working-tree snapshot"). The release build runs
# run_verification.sh on git HEAD.
# The full polychoice test suite (every tests/testthat file, NOT_CRAN set) on a
# snapshot of the working tree, as verification/run_verification.sh --tests
# does for git HEAD: the snapshot is `git stash create` (a commit object of the
# working tree, refs untouched), exported with git archive, fingerprinted the
# same way, installed into a private library. Logs go to
# verification/runs/<build>/ (status.txt, tests.log, fingerprint, install log).
#   bash verification/run_snapshot_tests.sh <work dir outside the repository>
set -u
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
RB="/c/Program Files/R/R-4.6.1/bin"
WORK="$1"; rm -rf "$WORK"; SRC="$WORK/polychoice"; mkdir -p "$SRC" "$WORK/lib"
to_win() { cygpath -m "$1"; }
HEAD="$(git -C "$ROOT" rev-parse HEAD)"
SNAP="$(git -C "$ROOT" stash create)"; [ -n "$SNAP" ] || SNAP="$HEAD"
git -C "$ROOT" archive "$SNAP" | tar -x -C "$SRC"
rm -rf "$SRC/verification/runs"; rm -f "$SRC"/src/*.o "$SRC"/src/*.dll "$SRC"/src/*.so
( cd "$SRC" && find R src data inst man tests vignettes DESCRIPTION NAMESPACE -type f 2>/dev/null | sort | xargs md5sum ) > "$WORK/fingerprint.txt"
( cd "$SRC" && find verification -type f | sort | xargs md5sum ) > "$WORK/scripts_fingerprint.txt"
PC_BUILD="$(md5sum "$WORK/fingerprint.txt" | cut -c1-12)"
echo "$PC_BUILD" > "$WORK/build_id.txt"
"$RB/R.exe" CMD INSTALL --no-multiarch -l "$(to_win "$WORK/lib")" "$(to_win "$SRC")" > "$WORK/install.log" 2>&1 || { echo "install failed: $WORK/install.log" > "$WORK/FAILED"; exit 1; }
OUT="$ROOT/verification/runs/$PC_BUILD"; mkdir -p "$OUT"
cp "$WORK/fingerprint.txt" "$WORK/scripts_fingerprint.txt" "$WORK/install.log" "$OUT/"
export PC_BUILD PC_COMMIT="$HEAD" R_LIBS="$(to_win "$WORK/lib")" PC_PKG="$(to_win "$SRC")"
STATUS="$OUT/status.txt"
echo "build $PC_BUILD  git $HEAD  working-tree snapshot $SNAP (git stash create; uncommitted changes ARE in this copy)  started $(date '+%F %T')  R_LIBS=$R_LIBS" >> "$STATUS"
cd "$SRC" || exit 1
: > "$OUT/tests.log"
for f in tests/testthat/test-*.R; do
  "$RB/Rscript.exe" verification/run_test_file.R "$(basename "$f")" >> "$OUT/tests.log" 2>&1
  echo "  tests $(basename "$f") exit $?" >> "$STATUS"
done
grep "total expectations" "$OUT/tests.log" | awk '{n += $5; f += $7; w += $9; e += $11; s += $13} END {printf "  tests: %d expectations, %d failed, %d warnings, %d errors, %d skipped\n", n, f, w, e, s}' >> "$STATUS"
echo "finished $(date '+%F %T')" >> "$STATUS"
echo done > "$WORK/DONE"
