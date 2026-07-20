#!/usr/bin/env Rscript
# METI-FS Computational Runtime Benchmark
# Writes timing data directly to CSV as each step completes.
# No in-memory result compilation; avoids all data-structure issues.

# ---- Locate project root ----
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
if (length(file_arg) > 0) {
  script_dir <- dirname(sub("--file=", "", file_arg[1]))
} else {
  script_dir <- "."
}
proj_root <- normalizePath(file.path(script_dir, ".."), winslash = "/")
stopifnot(dir.exists(file.path(proj_root, "R")))
setwd(proj_root)

# ---- Output file ----
out_dir <- "C:/Temp/metifs_bench"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
csv_file <- file.path(out_dir, "timing_benchmark_results.csv")
cat("config,step,description,elapsed_sec,mem_mb\n", file = csv_file)

# ---- Params ----
BOOTSTRAP_N <- 100
N_GENES     <- 13000
N_TP        <- 4L
TP_LABELS   <- c("4d", "7d", "14d", "21d")
TP_VALUES   <- c(4, 7, 14, 21)

TEST_CONFIGS <- list(
  small  = list(n_ind = 3, n_ctrl = 2, label = "n=20"),
  medium = list(n_ind = 4, n_ctrl = 3, label = "n=28"),
  large  = list(n_ind = 5, n_ctrl = 5, label = "n=40")
)

PIPELINE_STEPS <- list(
  list(file = "01_data_import.R",      desc = "01_data_import"),
  list(file = "02_preprocessing.R",    desc = "02_preprocessing"),
  list(file = "03_normalization_QC.R", desc = "03_norm_QC"),
  list(file = "04_DEG_analysis.R",     desc = "04_DEG_L1"),
  list(file = "06_maSigPro_trends.R",  desc = "06_maSigPro_L2"),
  list(file = "08_WGCNA.R",            desc = "08_WGCNA_L3"),
  list(file = "09A_candidate_pool.R",  desc = "09A_candidate_pool_L4"),
  list(file = "09C_ML_stability_selection.R", desc = "09C_ML_stability_L5a"),
  list(file = "09D_gap_union_selection.R",    desc = "09D_gap_union_L5b")
)

get_mem <- function() sum(gc(reset = TRUE)[, "(Mb)"])

# ---- Generate synthetic data ----
suppressPackageStartupMessages(library(org.Hs.eg.db))
all_ensembl <- grep("^ENSG", keys(org.Hs.eg.db, keytype = "ENSEMBL"), value = TRUE)

generate_data <- function(n_genes, n_ind, n_ctrl, n_tp, seed = 42) {
  set.seed(seed)
  n_samp <- (n_ind + n_ctrl) * n_tp
  use_genes <- sample(all_ensembl, n_genes)
  base_mean <- rgamma(n_genes, shape = 2, scale = 500)
  dispersion <- pmax(0.001, (0.05 + 4.0 / base_mean) * exp(rnorm(n_genes, 0, 0.3)))
  lib_size <- rnorm(n_samp, 5e6, 5e5); lib_size <- pmax(lib_size, 1e6)
  counts <- matrix(0L, n_genes, n_samp)
  for (i in seq_len(n_samp))
    counts[, i] <- rnbinom(n_genes, mu = base_mean * lib_size[i] / mean(lib_size), size = 1 / dispersion)
  n_signal <- round(n_genes * 0.10)
  signal_idx <- sample(n_genes, n_signal)
  half <- n_signal / 2
  fc_up <- c(1.0, 1.5, 2.0, 3.0)
  for (tp in seq_len(n_tp)) {
    sc <- (tp - 1L) * (n_ind + n_ctrl) + 1L; ec <- sc + n_ind - 1L
    for (col in sc:ec) {
      counts[signal_idx[1:half], col] <- round(counts[signal_idx[1:half], col] * fc_up[tp])
      counts[signal_idx[(half + 1):n_signal], col] <- round(counts[signal_idx[(half + 1):n_signal], col] / fc_up[tp])
    }
  }
  counts <- pmax(counts, 0L)
  rownames(counts) <- use_genes
  tpm <- counts
  for (i in seq_len(n_samp)) tpm[, i] <- counts[, i] / sum(counts[, i]) * 1e6
  rownames(tpm) <- use_genes
  list(counts = counts, tpm = tpm)
}

# ---- Run one config ----
run_one_config <- function(cfg, cfg_name) {

  n_ind <- cfg$n_ind; n_ctrl <- cfg$n_ctrl
  n_samp <- (n_ind + n_ctrl) * N_TP
  prefix <- "TimingTest"

  cat(sprintf("\n===== %s: %d genes, %d samples =====\n\n", cfg$label, N_GENES, n_samp))

  out_dir <- file.path("C:/Temp/metifs_bench", cfg_name)
  unlink(out_dir, recursive = TRUE)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # A: Generate data
  cat("[A] Generating data ... ")
  tA <- system.time({ sim <- generate_data(N_GENES, n_ind, n_ctrl, N_TP, seed = 42) })
  cat(sprintf("%.0f sec\n", tA["elapsed"]))
  append_line(cfg_name, "A_generate", "Data generation", tA["elapsed"], get_mem())

  # B: Prepare pipeline input
  pipeline_dir <- file.path(out_dir, prefix)
  data_raw <- file.path(pipeline_dir, "data_raw")
  for (d in c(data_raw, file.path(pipeline_dir, "data"), file.path(pipeline_dir, "Figure")))
    dir.create(d, recursive = TRUE, showWarnings = FALSE)

  sample_names <- character(0)
  for (tp in TP_LABELS) {
    for (r in seq_len(n_ind)) sample_names <- c(sample_names, paste0(prefix, tp, r))
    for (r in seq_len(n_ctrl)) sample_names <- c(sample_names, paste0(prefix, tp, "C", r))
  }
  colnames(sim$counts) <- sample_names; colnames(sim$tpm) <- sample_names

  write.table(as.data.frame(sim$counts),
              file.path(data_raw, paste0(prefix, "_all_counts_with_order.tsv")),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)
  write.table(as.data.frame(sim$tpm),
              file.path(data_raw, paste0(prefix, "_all_tpm.tsv")),
              sep = "\t", quote = FALSE, row.names = TRUE, col.names = NA)

  # C: Source 00_setup then each pipeline step
  old_project <- if (exists("PROJECT_DIR", envir = .GlobalEnv)) get("PROJECT_DIR", envir = .GlobalEnv) else NULL
  old_wd <- setwd(proj_root)
  assign("PROJECT_DIR", pipeline_dir, envir = .GlobalEnv)

  cat("[C] 00_setup.R ... ")
  t0 <- system.time({ source(file.path(proj_root, "R", "00_setup.R"), local = FALSE) })
  cat(sprintf("%.0f sec\n", t0["elapsed"]))

  for (step in PIPELINE_STEPS) {
    step_key <- step$file
    cat(sprintf("  [%s] %s ... ", format(Sys.time(), "%H:%M:%S"), step$desc))
    gc(reset = TRUE)
    t1 <- system.time({
      ok <- tryCatch({
        source(file.path(proj_root, "R", step$file), local = FALSE)
        TRUE
      }, error = function(e) { cat(sprintf("SKIP: %s\n", e$message)); FALSE })
    })
    if (ok) cat(sprintf("%.0f sec\n", t1["elapsed"]))
    append_line(cfg_name, gsub(".R$", "", step_key), step$desc, t1["elapsed"], get_mem())
  }

  setwd(old_wd)
  if (!is.null(old_project)) assign("PROJECT_DIR", old_project, envir = .GlobalEnv)
  unlink(out_dir, recursive = TRUE)
}

# ---- Direct-to-CSV helper ----
append_line <- function(cfg, step, desc, elapsed, mem) {
  cat(sprintf("%s,%s,%s,%.1f,%.0f\n", cfg, step, desc, elapsed, mem),
      file = csv_file, append = TRUE)
}

# ---- Main ----
cat("\n============================================================\n")
cat(sprintf("  METI-FS Runtime Benchmark  |  %s\n", format(Sys.time())))
cat(sprintf("  R %s.%s  |  %d cores  |  B=%d\n",
            R.version$major, R.version$minor, parallel::detectCores(), BOOTSTRAP_N))
cat("============================================================\n")

for (cfg_name in names(TEST_CONFIGS)) {
  cfg <- TEST_CONFIGS[[cfg_name]]
  tryCatch(run_one_config(cfg, cfg_name),
           error = function(e) cat(sprintf("FATAL %s: %s\n", cfg_name, e$message)))
  gc()
}

cat(sprintf("\nDone. Results: %s\n", csv_file))
