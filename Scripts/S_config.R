#!/usr/bin/env Rscript
# ==============================================================================
# S_config.R — METI-FS Methods paperConfigure (v2 — BIB)
#
# v2 :
# - Simulation parameters correlation_structure dispersion_source
# - : 27 core + 3 correlation = 30 unique × 5 repeats = 150 runs
#
# PathParameters。
# S source("S_config.R") 。
#
# Directory：
# METHODS_BASE/ <- Methods paperDirectory
# METHODS_BASE/Scripts/ <- S
# file.path(METHODS_BASE, "simulations")/ <- DataOutput
# METHODS_BASE/geo_datasets/ <- GEOData
# file.path(METHODS_BASE, "pipeline_runs")/ <- DatapipelineRunResult
# file.path(METHODS_BASE, "benchmark")_results/ <- benchmark/ablationResult
# METHODS_BASE/figures/ <- Methods paper
# METHODS_BASE/tables/ <- Methods paper
#
# file.path(METHODS_BASE, "..")/Scripts/ <- pipeline（，）
# file.path(METHODS_BASE, "..")/{Lineage_A,...}/ <- Data（）
#
# ：
# - Methods paperOutput METHODS_BASE/ 
# - pipeline（file.path(METHODS_BASE, "..")/Scripts/）source
# - Data directoryWriteFile
# ==============================================================================

# ---- 1. Path ----

METHODS_BASE    <- file.path(dirname(sys.frame(1)$ofile %||% "."), "..")
METHODS_SCRIPTS <- file.path(METHODS_BASE, "Scripts")
PIPELINE_SCRIPTS <- file.path(METHODS_BASE, "R")

# ---- 2. OutputDirectory ----

SIM_DIR         <- file.path(METHODS_BASE, "simulations")
GEO_DIR         <- file.path(METHODS_BASE, "geo_datasets")
GEO_DOWNLOAD    <- file.path(GEO_DIR, "downloads")
GEO_METADATA    <- file.path(GEO_DIR, "metadata")
RUN_DIR         <- file.path(METHODS_BASE, "pipeline_runs")
BENCH_DIR       <- file.path(METHODS_BASE, "benchmark_results")
FIG_DIR_METHODS <- file.path(METHODS_BASE, "figures")
TAB_DIR_METHODS <- file.path(METHODS_BASE, "tables")
SEARCH_DIR      <- file.path(METHODS_BASE, "search_results")

# ---- 3. CreateDirectory ----

all_dirs <- c(METHODS_SCRIPTS, SIM_DIR, GEO_DIR, GEO_DOWNLOAD, GEO_METADATA,
              RUN_DIR, BENCH_DIR, FIG_DIR_METHODS, TAB_DIR_METHODS, SEARCH_DIR)

for (d in all_dirs) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- 4. benchmarkParameters ----
# ，Parameters

BENCHMARK_PARAMS <- list(
 # 
  snr_levels = c("low", "medium", "high"),
  
 # SampleConfigure (n_reps_induced, n_reps_control)
  sample_configs = list(
    small  = list(n_reps_ind = 2, n_reps_ctrl = 2),   # 16 samples
 medium = list(n_reps_ind = 4, n_reps_ctrl = 3), # 28 samples (matches original study)
    large  = list(n_reps_ind = 5, n_reps_ctrl = 5)    # 40 samples
  ),
  
 # true markerNumberConfigure
  marker_configs = list(
    sparse = list(n_true_temporal = 10, n_true_main = 30,  n_true_timeonly = 60),
    medium = list(n_true_temporal = 20, n_true_main = 50,  n_true_timeonly = 100),
    dense  = list(n_true_temporal = 50, n_true_main = 100, n_true_timeonly = 200)
  ),
  
 # ReplicateCount
  n_repeats = 5,
  
 # Gene（filterByExpr）
  n_genes = 13000,
  
 # Time（MSCExperiment）
  n_timepoints = 4,
  time_values = c(4, 7, 14, 21),
  
 # ---- [v2] DispersionParameters ----
  # DESeq2-style: α(μ) = α_intercept + α_slope / μ + log-normal noise
 # bulk RNA-seq:
 # α_intercept = 0.05 (ExpressionGenedispersion)
 # α_slope = 4.0 (ExpressionGenedispersion)
 # noise_sd = 0.5 (gene-to-gene)
 # : powsimR (Vieth et al. 2017, Bioinformatics 33(21):3486)
  #        compcodeR (Soneson 2014, Bioinformatics 30(18):2670)
  dispersion_source = "deseq2",
  disp_intercept = 0.05,
  disp_slope = 4.0,
  disp_noise_sd = 0.5,
  
  # ---- [v2] Correlation structure ----
 # : medium SNR × medium sample × medium marker Testblock correlation
 # 3 SNR × 5 repeats = 15 Run
 # : Hédou et al. 2024 (Stabl) correlated features (R≈0.5) Test
 # Langfelder & Horvath 2008 (WGCNA) ModuleExpression
  include_correlation = TRUE,
 n_modules = 4, # TRUE_TEMPORALGene4ExpressionModule
 rho_within = 0.4 # Module (, StablR≈0.5Test)
)

# ---- 5. GEOData ----
# Data

GEO_REGISTRY <- list(
  GSE197067 = list(
    gse = "GSE197067",
    prefix = "Tcell",
    domain = "immune_activation",
    n_samples = 48,
    n_timepoints = 6,
    time_labels = c("0h", "6h", "12h", "24h", "48h", "72h"),
    design = "4 healthy donors x activated vs non-activated x 6 timepoints",
    note = "Pan T-cells anti-CD3/CD28 activation. Best dataset."
  ),
  GSE303975 = list(
    gse = "GSE303975",
    prefix = "PCa",
    domain = "cancer_drug_response",
    n_samples = 18,
    n_timepoints = 3,
    time_labels = c("8h", "24h", "72h"),
    design = "LNCaP, Combination vs DMSO x 3 timepoints x 3 reps",
    note = "Prostate cancer combination therapy. Use Combination vs DMSO only."
  ),
  GSE307424 = list(
    gse = "GSE307424",
    prefix = "Lung",
    domain = "cancer_drug_response",
    n_samples = 18,
    n_timepoints = 3,
    time_labels = c("6h", "48h", "72h"),
    design = "NCI-H1693, PRT3789 vs DMSO x 3 timepoints x 3 reps",
    note = "Lung cancer SMARCA2 degrader. Smallest dataset."
  )
)

# ---- 6. Function ----

methods_log <- function(step, msg) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] %s: %s", timestamp, step, msg)
  cat(log_line, "\n")
  log_file <- file.path(METHODS_BASE, "methods_analysis_log.txt")
  cat(log_line, "\n", file = log_file, append = TRUE)
}

# ---- 7. Configure ----

cat("============================================================\n")
cat("  METI-FS Methods Paper — Configuration Loaded (v2)\n")
cat("============================================================\n")
cat("  Methods base:      ", METHODS_BASE, "\n")
cat("  Methods scripts:   ", METHODS_SCRIPTS, "\n")
cat("  Pipeline scripts:  ", PIPELINE_SCRIPTS, "\n")
cat("  Simulations:       ", SIM_DIR, "\n")
cat("  GEO datasets:      ", GEO_DIR, "\n")
cat("  Pipeline runs:     ", RUN_DIR, "\n")
cat("  Benchmark results: ", BENCH_DIR, "\n")
cat("  Figures:           ", FIG_DIR_METHODS, "\n")
cat("------------------------------------------------------------\n")
cat("  Simulation params (v2):\n")
cat("    Dispersion:     ", BENCHMARK_PARAMS$dispersion_source, "\n")
cat("    Correlation:    ", ifelse(BENCHMARK_PARAMS$include_correlation, 
                                   "independent + block", "independent only"), "\n")
cat("    Core scenarios:  27 × 5 = 135\n")
if (BENCHMARK_PARAMS$include_correlation) {
  cat("    Corr scenarios:  3 × 5 = 15\n")
  cat("    Total:           150 runs\n")
} else {
  cat("    Total:           135 runs\n")
}
cat("============================================================\n")
