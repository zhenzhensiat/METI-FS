#!/usr/bin/env Rscript
# ==============================================================================
# S_config.R — METI-FS Shared Configuration (v3)
#
# Extended co-expression structures, temporal depth extension, DESeq2-style
# dispersion engine, and unified benchmarking framework.
#
# All paths and global parameters defined here. Other S-series scripts
# source this file to inherit configuration.
#
# Usage:
#   Set METI_FS_PROJECT environment variable to your project root, or run
#   scripts from within the project directory.
#   Sys.setenv(METI_FS_PROJECT = "/path/to/METI-FS")
#   source("scripts/S_config.R")
#
# Directory structure (created automatically):
#   <project>/
#     R/                    <- Core pipeline scripts
#     scripts/              <- Simulation, benchmark, analysis scripts
#     simulations/          <- Synthetic data output
#     geo_datasets/         <- Downloaded GEO public datasets
#     pipeline_runs/        <- Pipeline execution results per dataset
#     benchmark_results/    <- Aggregated benchmark/ablation results
#     figures/              <- Manuscript figures
#     tables/               <- Manuscript tables
#     search_results/       <- GEO search output
# ==============================================================================

# ---- 1. Project root ----

METHODS_BASE <- Sys.getenv("METI_FS_PROJECT", unset = getwd())

# ---- 2. Core pipeline scripts ----

PIPELINE_SCRIPTS <- file.path(METHODS_BASE, "R")

# ---- 3. Output subdirectories ----

METHODS_SCRIPTS <- file.path(METHODS_BASE, "scripts")
SIM_DIR         <- file.path(METHODS_BASE, "simulations")
GEO_DIR         <- file.path(METHODS_BASE, "geo_datasets")
GEO_DOWNLOAD    <- file.path(GEO_DIR, "downloads")
GEO_METADATA    <- file.path(GEO_DIR, "metadata")
RUN_DIR         <- file.path(METHODS_BASE, "pipeline_runs")
BENCH_DIR       <- file.path(METHODS_BASE, "benchmark_results")
FIG_DIR_METHODS <- file.path(METHODS_BASE, "figures")
TAB_DIR_METHODS <- file.path(METHODS_BASE, "tables")
SEARCH_DIR      <- file.path(METHODS_BASE, "search_results")

# ---- 4. Create all directories ----

all_dirs <- c(METHODS_SCRIPTS, SIM_DIR, GEO_DIR, GEO_DOWNLOAD, GEO_METADATA,
              RUN_DIR, BENCH_DIR, FIG_DIR_METHODS, TAB_DIR_METHODS, SEARCH_DIR)

for (d in all_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- 5. Simulation benchmark parameters ----

BENCHMARK_PARAMS <- list(
  # Signal-to-noise ratio levels
  snr_levels = c("low", "medium", "high"),

  # Sample size configurations (n_reps_induced, n_reps_control)
  sample_configs = list(
    small  = list(n_reps_ind = 2, n_reps_ctrl = 2),   # 16 samples
    medium = list(n_reps_ind = 4, n_reps_ctrl = 3),   # 28 samples
    large  = list(n_reps_ind = 5, n_reps_ctrl = 5)    # 40 samples
  ),

  # True marker gene counts per category
  marker_configs = list(
    sparse = list(n_true_temporal = 10, n_true_main = 30,  n_true_timeonly = 60),
    medium = list(n_true_temporal = 20, n_true_main = 50,  n_true_timeonly = 100),
    dense  = list(n_true_temporal = 50, n_true_main = 100, n_true_timeonly = 200)
  ),

  # Number of replicate simulations per scenario
  n_repeats = 5,

  # Number of genes (post-filterByExpr)
  n_genes = 13000,

  # Timepoint design (default: MSC-like)
  n_timepoints = 4,
  time_values = c(4, 7, 14, 21),

  # ---- DESeq2-style dispersion parameters ----
  # alpha(mu) = disp_intercept + disp_slope / mu + log-normal noise
  # Ref: powsimR (Vieth et al. 2017); compcodeR (Soneson 2014)
  dispersion_source = "deseq2",
  disp_intercept = 0.05,
  disp_slope = 4.0,
  disp_noise_sd = 0.5,

  # ---- Co-expression structure (v3) ----
  # Five architectures spanning biological complexity
  # block4_rho0.4:  4 modules, within-module rho = 0.4 (default)
  # block4_rho0.2:  4 modules, within-module rho = 0.2 (weak correlation)
  # block4_rho0.6:  4 modules, within-module rho = 0.6 (strong correlation)
  # hierarchical_8: 8 sub-modules in 2 super-modules (rho_super=0.3, rho_sub=0.5)
  # overlapping:    6 modules, 20% genes belong to >=2 modules
  include_correlation = TRUE,
  correlation_structures = c("block4_rho0.4", "block4_rho0.2", "block4_rho0.6",
                              "hierarchical_8", "overlapping"),
  n_modules = 4,       # Backward-compatible
  rho_within = 0.4,    # Backward-compatible

  # ---- Temporal depth extension (v3) ----
  include_temporal_extension = TRUE,
  temporal_extension_tps = c(6, 8),
  temporal_extension_n_repeats = 3
)

# ---- 6. GEO dataset registry ----

GEO_REGISTRY <- list(
  GSE197067 = list(
    gse = "GSE197067",
    prefix = "Tcell",
    domain = "immune_activation",
    n_samples = 48,
    n_timepoints = 6,
    time_labels = c("0h", "6h", "12h", "24h", "48h", "72h"),
    design = "4 healthy donors x activated vs non-activated x 6 timepoints",
    note = "Pan T-cells anti-CD3/CD28 activation"
  ),
  GSE303975 = list(
    gse = "GSE303975",
    prefix = "PCa",
    domain = "cancer_drug_response",
    n_samples = 18,
    n_timepoints = 3,
    time_labels = c("8h", "24h", "72h"),
    design = "LNCaP, Combination vs DMSO x 3 timepoints x 3 reps",
    note = "Prostate cancer combination therapy"
  ),
  GSE307424 = list(
    gse = "GSE307424",
    prefix = "Lung",
    domain = "cancer_drug_response",
    n_samples = 18,
    n_timepoints = 3,
    time_labels = c("6h", "48h", "72h"),
    design = "NCI-H1693, PRT3789 vs DMSO x 3 timepoints x 3 reps",
    note = "Lung cancer SMARCA2 degrader"
  )
)

# ---- 7. Logging utility ----

methods_log <- function(step, msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] %s: %s", timestamp, step, msg)
  cat(log_line, "\n")
  log_file <- file.path(METHODS_BASE, "methods_analysis_log.txt")
  cat(log_line, "\n", file = log_file, append = TRUE)
}

# ---- 8. Print configuration summary ----

cat("============================================================\n")
cat("  METI-FS Methods Paper — Configuration Loaded (v3)\n")
cat("============================================================\n")
cat("  Project root:      ", METHODS_BASE, "\n")
cat("  Pipeline scripts:  ", PIPELINE_SCRIPTS, "\n")
cat("  Simulations:       ", SIM_DIR, "\n")
cat("  GEO datasets:      ", GEO_DIR, "\n")
cat("  Pipeline runs:     ", RUN_DIR, "\n")
cat("  Benchmark results: ", BENCH_DIR, "\n")
cat("  Figures:           ", FIG_DIR_METHODS, "\n")
cat("------------------------------------------------------------\n")
cat("  Simulation params (v3):\n")
cat("    Dispersion:     ", BENCHMARK_PARAMS$dispersion_source, "\n")
cat("    Structures:     ", paste(BENCHMARK_PARAMS$correlation_structures, collapse=", "), "\n")
cat("============================================================\n")
