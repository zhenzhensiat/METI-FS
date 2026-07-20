#!/usr/bin/env Rscript
# ==============================================================================
# run_rf_weighted_demo.R -- METI-FS A9_RF_weighted ensemble vignette
#
# This script demonstrates the weighted Random Forest ensemble option, which
# provides an alternative to the default binary RF exclusion in Layer 5.
# Instead of discarding all RF-selected genes when RF stability is low, the
# weighted ensemble retains RF genes proportionally to their Nogueira stability
# index relative to LASSO and SVM-RFE.
#
# When to use binary exclusion (default):
#   - WGCNA-enriched candidate pools with strong co-expression structure
#   - RF Nogueira index consistently below 0.5 across bootstrap iterations
#   - When minimizing false positives is the primary objective
#
# When to use weighted ensemble (A9_RF_weighted):
#   - Datasets with weak or uncertain co-expression enrichment
#   - When WGCNA layer is bypassed (skip_wgcna = TRUE)
#   - When RF may capture nonlinear interactions missed by linear methods
#   - Exploratory analyses where sensitivity is prioritized over specificity
#
# Usage:
#   source("scripts/run_rf_weighted_demo.R")
#
# Dependencies: METI-FS R/ and scripts/ directories must be on the search path
#   or in the working directory.
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  METI-FS A9_RF_weighted Ensemble Vignette\n")
cat("  Weighted Random Forest integration in bootstrap stability selection\n")
cat("================================================================\n\n")

# ---- Locate scripts ----
file_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(sub("--file=", "", file_arg[1]))
} else {
  script_dir <- "."
}
r_dir <- file.path(dirname(script_dir), "R")
if (!dir.exists(r_dir)) r_dir <- "R"

# Source ablation framework (contains apply_rf_weighted_ensemble)
ablation_path <- file.path(script_dir, "S06_ablation_study.R")
if (!file.exists(ablation_path)) {
  stop("S06_ablation_study.R not found. Ensure this script is in the scripts/ directory.")
}
source(ablation_path)

# Source config for paths
config_path <- file.path(script_dir, "S_config.R")
if (file.exists(config_path)) source(config_path)

# ---- Step 1: Generate test data ----
cat("Step 1: Generating test dataset with weak co-expression structure\n\n")

# Use a weaker correlation structure: this is where weighted RF may be preferred
# because RF stability is less degraded by co-expression enrichment
test_dir <- tempfile("rf_weighted_demo_")
dir.create(test_dir, recursive = TRUE)

if (exists("generate_simulation")) {
  sim <- generate_simulation(
    n_genes = 5000,
    n_timepoints = 4,
    time_values = c(0, 6, 24, 48),
    n_reps_ind = 5,
    n_reps_ctrl = 5,
    n_true_temporal = 20,
    n_true_main = 50,
    n_true_timeonly = 50,
    snr = "medium",
    dispersion_source = "deseq2",
    correlation_structure = "independent",
    output_dir = file.path(test_dir, "simulation"),
    seed = 42,
    verbose = FALSE
  )
} else {
  stop("Simulation engine not available. Source S01_simulation_engine.R first.")
}

# ---- Step 2: Run pipeline (binary exclusion, default) ----
cat("\nStep 2: Running pipeline with default binary RF exclusion\n\n")

pipeline_dir <- file.path(test_dir, "pipeline_default")
data_raw_dir <- file.path(pipeline_dir, "data_raw")
data_dir     <- file.path(pipeline_dir, "data")
fig_dir      <- file.path(pipeline_dir, "Figure")
for (d in c(data_raw_dir, data_dir, fig_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

prefix <- "RFDemo"
write.csv(as.data.frame(sim$counts),
          file.path(data_raw_dir, paste0(prefix, "_counts.csv")), row.names = TRUE)
write.csv(as.data.frame(sim$tpm),
          file.path(data_raw_dir, paste0(prefix, "_tpm.csv")), row.names = TRUE)
write.csv(sim$sample_info,
          file.path(data_raw_dir, paste0(prefix, "_sample_info.csv")), row.names = FALSE)
saveRDS(sim$ground_truth, file.path(data_dir, "ground_truth.rds"))

# Run pipeline with default binary RF exclusion
old_project <- if (exists("PROJECT_DIR", envir = .GlobalEnv)) {
  get("PROJECT_DIR", envir = .GlobalEnv)
} else NULL
assign("PROJECT_DIR", pipeline_dir, envir = .GlobalEnv)

source(file.path(r_dir, "00_setup.R"), local = FALSE)
source(file.path(r_dir, "01_data_import.R"), local = FALSE)
source(file.path(r_dir, "02_preprocessing.R"), local = FALSE)
source(file.path(r_dir, "03_normalization_QC.R"), local = FALSE)
source(file.path(r_dir, "04_DEG_analysis.R"), local = FALSE)
source(file.path(r_dir, "06_maSigPro_trends.R"), local = FALSE)
source(file.path(r_dir, "08_WGCNA.R"), local = FALSE)
source(file.path(r_dir, "09A_candidate_pool.R"), local = FALSE)
source(file.path(r_dir, "09C_ML_stability_selection.R"), local = FALSE)
source(file.path(r_dir, "09D_gap_union_selection.R"), local = FALSE)

# Capture binary exclusion results
binary_ml_genes <- get("final_ml_genes", envir = .GlobalEnv)

# Capture Nogueira-like stability (frequency-based proxy)
stab_data <- get("stab_result", envir = .GlobalEnv)
if (!is.null(stab_data)) {
  binary_rf_excluded <- TRUE  # Default: RF is hard-excluded
  cat(sprintf("\n  Binary exclusion: RF excluded from union (default)\n"))
  cat(sprintf("  ML genes (LASSO + SVM only): %d\n", length(binary_ml_genes)))
}

# Save upstream for weighted ensemble step
upstream_saved <- list(
  stab = stab_data,
  gap_union = list(final_gene_ids = binary_ml_genes)
)

# ---- Step 3: Apply A9_RF_weighted ensemble ----
cat("\nStep 3: Applying A9_RF_weighted ensemble\n\n")

weighted_result <- apply_rf_weighted_ensemble(upstream_saved)

cat(sprintf("  Nogueira-based RF weight: %.3f\n", weighted_result$rf_weight))
cat(sprintf("  RF genes retained (weight >= %.0f threshold): %d\n",
            ifelse(weighted_result$rf_weight >= 0.5, 0.5, 0),
            length(weighted_result$rf_genes)))
cat(sprintf("  LASSO + SVM genes: %d\n",
            length(setdiff(weighted_result$gene_ids, weighted_result$rf_genes))))
cat(sprintf("  Total ML genes (weighted ensemble): %d\n",
            length(weighted_result$gene_ids)))

# ---- Step 4: Compare results ----
cat("\nStep 4: Comparison\n\n")

only_binary <- setdiff(binary_ml_genes, weighted_result$gene_ids)
only_weighted <- setdiff(weighted_result$gene_ids, binary_ml_genes)
shared <- intersect(binary_ml_genes, weighted_result$gene_ids)

cat(sprintf("  Genes selected by BOTH methods:  %d\n", length(shared)))
cat(sprintf("  Genes EXCLUSIVE to binary exclusion: %d\n", length(only_binary)))
cat(sprintf("  Genes EXCLUSIVE to weighted ensemble: %d\n", length(only_weighted)))

if (length(only_weighted) > 0) {
  cat("\n  Genes uniquely added by weighted RF ensemble:\n")
  for (g in only_weighted) {
    # Check ground truth if available
    is_true <- if (!is.null(sim$ground_truth)) {
      g %in% sim$ground_truth$TRUE_TEMPORAL
    } else NA
    cat(sprintf("    %s  (true temporal: %s)\n", g, is_true))
  }
}

# ---- Step 5: Decision guidance ----
cat("\n================================================================\n")
cat("  DECISION GUIDANCE\n")
cat("================================================================\n\n")

cat("When to prefer binary exclusion (default):\n")
cat("  - WGCNA-enriched candidate pools where co-expression is strong\n")
cat("  - RF Nogueira index < 0.5, indicating below-moderate stability\n")
cat("  - Goal: minimize false positives, maximize selection reliability\n\n")

cat("When to prefer weighted ensemble (A9_RF_weighted):\n")
cat("  - Weak or absent co-expression structure (skip_wgcna = TRUE)\n")
cat("  - RF Nogueira index >= 0.5, indicating adequate stability\n")
cat("  - Goal: capture nonlinear interactions, maximize sensitivity\n")
cat("  - Exploratory analyses where biological validation will follow\n\n")

cat("How to use in the METI-FS R package:\n")
cat("  1. Run pipeline normally to obtain upstream results\n")
cat("  2. Examine RF Nogueira index in the stability diagnostics\n")
cat("  3. If RF index >= 0.5 or co-expression is weak, call:\n")
cat("     result <- apply_rf_weighted_ensemble(upstream)\n")
cat("  4. The function returns weighted gene IDs plus diagnostic weight\n\n")

# Cleanup
if (!is.null(old_project)) {
  assign("PROJECT_DIR", old_project, envir = .GlobalEnv)
}
unlink(test_dir, recursive = TRUE)

cat(sprintf("Demo complete. Temporary files cleaned up.\n"))
cat("For integration into automated pipeline, set rf_mode = 'weighted' in 00_setup.R\n")
cat("See scripts/S06_ablation_study.R lines 465-540 for full implementation.\n")
