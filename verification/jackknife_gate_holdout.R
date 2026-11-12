## The fresh-seed holdout of the jackknife gate calibration (ledger, gate
## selection rule): the calibration Monte Carlo with seed base 9000 (the
## selection run used 5000), written to jackknife_gate_holdout.log / .rds.
## Run from the package root:  Rscript verification/jackknife_gate_holdout.R
Sys.setenv(JK_SEED_BASE = Sys.getenv("JK_SEED_BASE_HOLDOUT", "9000"), JK_TAG = "jackknife_gate_holdout")
source(file.path(if (dir.exists("verification")) "verification" else ".", "jackknife_gate_mc.R"))
