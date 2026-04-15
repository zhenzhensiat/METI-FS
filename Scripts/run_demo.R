#!/usr/bin/env Rscript
# ==============================================================================
# run_demo.R — METI-FS Quick-Start Demo
#
# This script demonstrates the complete METI-FS pipeline execution.
# Before running, ensure:
#   1. PROJECT_DIR points to a directory containing data_raw/
#   2. data_raw/ contains: {prefix}_counts.csv, {prefix}_tpm.csv, {prefix}_sample_info.csv
#   3. All R dependencies are installed (see README.md)
#
# Usage:
#   PROJECT_DIR <- "/path/to/your/project"
#   source("scripts/run_demo.R")
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  METI-FS Pipeline — Quick Start Demo\n")
cat("================================================================\n\n")

# ---- Check PROJECT_DIR ----
if (!exists("PROJECT_DIR") || !dir.exists(PROJECT_DIR)) {
  stop("Please set PROJECT_DIR to a valid directory before running this demo.\n",
       "  Example: PROJECT_DIR <- '/path/to/your/project'\n",
       "  The directory should contain a data_raw/ subdirectory with input files.")
}

cat(sprintf("  Project: %s\n", PROJECT_DIR))
cat(sprintf("  Started: %s\n\n", format(Sys.time())))

# ---- Locate pipeline scripts ----
# Assumes this script is in scripts/ and R/ is a sibling directory
script_dir <- dirname(sys.frame(1)$ofile %||% ".")
r_dir <- file.path(dirname(script_dir), "R")
if (!dir.exists(r_dir)) {
  r_dir <- "R"  # Fallback: assume working directory is repo root
}

run_step <- function(script_name, description) {
  cat(sprintf("  [%s] %s ... ", format(Sys.time(), "%H:%M:%S"), description))
  t0 <- Sys.time()
  source(file.path(r_dir, script_name), local = FALSE)
  elapsed <- round(difftime(Sys.time(), t0, units = "secs"), 1)
  cat(sprintf("done (%.1fs)\n", elapsed))
}

# ---- Run pipeline ----
t_start <- Sys.time()

run_step("00_setup.R",            "Step 0:  Environment setup")
run_step("01_data_import.R",      "Step 1:  Data import")
run_step("02_preprocessing.R",    "Step 2:  Preprocessing")
run_step("03_normalization_QC.R", "Step 3:  Normalization + QC")
run_step("04_DEG_analysis.R",     "Step 4:  DEG analysis")
run_step("06_maSigPro_trends.R",  "Step 6:  Temporal trends")
run_step("08_WGCNA.R",            "Step 8:  WGCNA network")
run_step("09A_candidate_pool.R",  "Step 9A: Candidate pool")
run_step("09C_ML_stability_selection.R", "Step 9C: ML stability")
run_step("09D_gap_union_selection.R",    "Step 9D: Gap-union")
run_step("09F_PPI_hub_selection.R",      "Step 9F: PPI hubs")
run_step("10_integration.R",      "Step 10: Integration")

total_min <- round(as.numeric(difftime(Sys.time(), t_start, units = "mins")), 1)

cat("\n")
cat("================================================================\n")
cat(sprintf("  COMPLETE — Total runtime: %.1f minutes\n", total_min))
cat(sprintf("  Results saved to: %s/data/\n", PROJECT_DIR))
cat("================================================================\n")
