#!/usr/bin/env Rscript
# ==============================================================================
# S01b_benchmark_timing.R — Test
#
# ： →→pipeline→ ，
# ，150Time，。
#
# Usage：
#   setwd(file.path(METHODS_BASE, "scripts"))
#   source("S01b_benchmark_timing.R")
#
# Output：
# + Time
#   BENCH_DIR/timing_test_results.txt
# ==============================================================================

source("S_config.R")
source("S01_simulation_engine.R")

cat("\n")
cat("============================================================\n")
cat("  METI-FS Benchmark Timing Test\n")
cat("  Running 1 scenario (medium_medium_medium) end-to-end\n")
cat("============================================================\n\n")

timing <- list()
test_dir <- file.path(RUN_DIR, "_timing_test")

# ---- Step 1: Data ----
cat("[TIMING] Step 1: Generate simulation data...\n")
t1 <- system.time({
  sim <- generate_simulation(
    n_genes = 13000,
    n_timepoints = 4,
    time_values = c(4, 7, 14, 21),
    n_reps_ind = 4,
    n_reps_ctrl = 3,
    n_true_temporal = 20,
    n_true_main = 50,
    n_true_timeonly = 100,
    snr = "medium",
    dispersion_source = "deseq2",
    correlation_structure = "independent",
    output_dir = file.path(test_dir, "simulation"),
    seed = 9999,
    verbose = TRUE
  )
})
timing$step1_simulation <- t1["elapsed"]
cat(sprintf("  → Step 1 done: %.1f seconds\n\n", t1["elapsed"]))

# ---- Step 2: pipelineInputFormat (S02) ----
cat("[TIMING] Step 2: Adapt to pipeline format (S02)...\n")
t2 <- system.time({
 # S02: Createdata_raw/Directory
  pipeline_dir <- file.path(test_dir, "pipeline")
  data_raw_dir <- file.path(pipeline_dir, "data_raw")
  dir.create(data_raw_dir, recursive = TRUE, showWarnings = FALSE)
  
 # Savecounts (S04Format)
  saveRDS(sim$counts, file.path(data_raw_dir, "counts_matrix.rds"))
  saveRDS(sim$tpm, file.path(data_raw_dir, "tpm_matrix.rds"))
  saveRDS(sim$sample_info, file.path(data_raw_dir, "sample_info.rds"))
  
 # ground truthdata/
  data_dir <- file.path(pipeline_dir, "data")
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  saveRDS(sim$ground_truth, file.path(data_dir, "ground_truth.rds"))
})
timing$step2_adapt <- t2["elapsed"]
cat(sprintf("  → Step 2 done: %.1f seconds\n\n", t2["elapsed"]))

# ---- Step 3: PipelineRun (S04) — ----
cat("[TIMING] Step 3: Full pipeline run (S04)...\n")
cat("  *** This is the bottleneck step ***\n")
cat("  If this takes >60 minutes, consider reducing bootstrap to 50.\n\n")

# CheckS04Available
s04_file <- file.path(METHODS_SCRIPTS, "S04_run_pipeline_wrapper.R")
if (!file.exists(s04_file)) {
  s04_file <- "S04_run_pipeline_wrapper.R"
}

if (file.exists(s04_file)) {
  t3 <- system.time({
    tryCatch({
 # S04SetupWorking directory
 # TestsourceRun
      cat("  Attempting to source S04 and run pipeline...\n")
      cat("  (If S04 requires interactive setup, this may need manual execution)\n")
      
 # : S04Run，Time
      cat("\n  *** MANUAL TIMING NEEDED ***\n")
      cat("  Please run the following commands manually and time them:\n")
      cat("  ─────────────────────────────────────────────\n")
      cat(sprintf("  setwd('%s')\n", METHODS_SCRIPTS))
      cat(sprintf("  source('S04_run_pipeline_wrapper.R')\n"))
      cat(sprintf("  # Set run_dir to: %s\n", pipeline_dir))
      cat("  # Record the total elapsed time\n")
      cat("  ─────────────────────────────────────────────\n")
      
    }, error = function(e) {
      cat(sprintf("  S04 could not auto-run: %s\n", e$message))
    })
  })
 timing$step3_pipeline <- NA # 
} else {
  cat("  S04 not found at expected path. Manual timing required.\n")
  timing$step3_pipeline <- NA
}
cat("\n")

# ---- Step 4: Result (S05) — ----
cat("[TIMING] Step 4: Result collection (S05) — skipped (needs pipeline output)\n")
timing$step4_collect <- NA

# ---- ----
cat("\n")
cat("============================================================\n")
cat("  TIMING RESULTS\n")
cat("============================================================\n")
cat(sprintf("  Step 1 (Simulation generation):  %.1f sec\n", timing$step1_simulation))
cat(sprintf("  Step 2 (Format adaptation):      %.1f sec\n", timing$step2_adapt))
cat(sprintf("  Step 3 (Pipeline execution):     %s\n", 
            ifelse(is.na(timing$step3_pipeline), "MANUAL TIMING NEEDED", 
                   sprintf("%.1f sec", timing$step3_pipeline))))
cat(sprintf("  Step 4 (Result collection):      %s\n",
            ifelse(is.na(timing$step4_collect), "~10-30 sec (estimated)",
                   sprintf("%.1f sec", timing$step4_collect))))

cat("\n------------------------------------------------------------\n")
cat("  PROJECTIONS (fill in Step 3 time after manual run):\n")
cat("------------------------------------------------------------\n")

# 
MANUAL_STEP3_MINUTES <- NA # ← S04

if (!is.na(MANUAL_STEP3_MINUTES)) {
  per_scenario_min <- timing$step1_simulation / 60 + 
                      timing$step2_adapt / 60 + 
                      MANUAL_STEP3_MINUTES + 
 0.5 # Step 4
  
  total_150 <- per_scenario_min * 150
  
  cat(sprintf("  Per scenario:   ~%.0f minutes\n", per_scenario_min))
  cat(sprintf("  150 scenarios:  ~%.0f minutes = %.1f hours = %.1f days\n",
              total_150, total_150 / 60, total_150 / 60 / 24))
  cat("\n")
  
  if (total_150 / 60 > 100) {
    cat("  ⚠️  WARNING: >100 hours estimated.\n")
    cat("  Consider:\n")
    cat("    - Reducing bootstrap from 100 to 50 in S04 wrapper\n")
    cat("    - Removing 'dense' marker scenarios (saves 45 runs)\n")
    cat("    - Running overnight / weekend\n")
  } else if (total_150 / 60 > 48) {
    cat("  ℹ️  48-100 hours: manageable with overnight runs over a few days.\n")
  } else {
    cat("  ✅  <48 hours: very manageable. Can finish in 2-3 overnight runs.\n")
  }
} else {
  cat("  Step 3 time not yet measured.\n")
  cat("  After running S04 manually on the timing test data, update:\n")
  cat("    MANUAL_STEP3_MINUTES <- XX  (line ~130 of this script)\n")
  cat("  Then re-run this script to get projections.\n")
  cat("\n")
  cat("  Quick estimates based on typical hardware:\n")
  cat("    If Step 3 ≈ 20 min → 150 runs ≈ 53 hours (2.2 days)\n")
  cat("    If Step 3 ≈ 40 min → 150 runs ≈ 103 hours (4.3 days)\n")
  cat("    If Step 3 ≈ 60 min → 150 runs ≈ 153 hours (6.4 days)\n")
}

cat("============================================================\n")

# ---- SavetimingResult ----
timing_file <- file.path(BENCH_DIR, "timing_test_results.txt")
sink(timing_file)
cat("METI-FS Benchmark Timing Test Results\n")
cat(sprintf("Date: %s\n", Sys.time()))
cat(sprintf("Step 1 (Simulation): %.1f sec\n", timing$step1_simulation))
cat(sprintf("Step 2 (Adaptation): %.1f sec\n", timing$step2_adapt))
cat(sprintf("Step 3 (Pipeline):   %s\n", 
            ifelse(is.na(timing$step3_pipeline), "PENDING", 
                   sprintf("%.1f sec", timing$step3_pipeline))))
cat(sprintf("Test directory: %s\n", test_dir))
sink()
cat(sprintf("\nTiming results saved: %s\n", timing_file))

# ---- Cleanup ----
cat(sprintf("\nTest files created in: %s\n", test_dir))
cat("You can delete this directory after timing is complete:\n")
cat(sprintf("  unlink('%s', recursive = TRUE)\n", test_dir))
