#!/usr/bin/env Rscript
# ==============================================================================
# S05_benchmark_collector.R — Pipeline run result collection and performance evaluation
#
# Functions:
#   1. Extract screening layer data from completed pipeline run directories
#   2. Simulation mode: compare with ground truth, calculate precision/recall/FDR/F1
#   3. Real data mode: collect funnel layer gene counts, stability structure, gap quality
#   4. Unified output format for S06 (ablation study) and S07 (summary figures)
#
# Input:
#   - Pipeline run directory (with .rds and .csv result files under data/ subdirectory)
#   - For simulation data also requires data/ground_truth.rds
#
# Output:
#   BENCH_DIR/benchmark_results_{run_id}.rds   — Complete metrics for single run
#   BENCH_DIR/benchmark_master.csv              — Summary table for all runs (append mode)
#
# Usage:
#   source("S05_benchmark_collector.R")
#   # Single run directory:
#   res <- collect_run(file.path(METHODS_BASE, "pipeline_runs", "GEO_GSE307424_Lung"))
#   # Batch collect simulation data:
#   collect_all_simulations(file.path(METHODS_BASE, "simulations/benchmark"))
#   # Batch collect GEO data:
#   collect_all_geo()
#
# Literature basis:
#   - pipeComp (Germain et al. 2020, Genome Biology 21:227):
#     multi-level evaluation metrics for pipeline benchmarking
#   - Stabl (Hédou et al. 2024, Nature Biotechnology 42:1581-1593):
#     FDR/stability evaluation, precision/recall on synthetic data
#   - Nogueira, Sechidis & Brown 2018, JMLR 18(174):1-54:
#     stability index (using stabm R package implementation, Bommert 2021)
#
# v2 changes:
#   [2.2] Added in collect_stability() pairwise Jaccard index
#         Ref: Real & Vargas (1996) Syst. Biol. 45(3):380-391
#   [2.3] Added in build_summary_row() jaccard_mean/sd columns
# ==============================================================================

# ---- 0. Load configuration ----
if (file.exists("S_config.R")) {
  source("S_config.R")
} else if (file.exists(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))) {
  source(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))
}

# S01 evaluate_selection function (needed for simulation data evaluation）
if (file.exists("S01_simulation_engine.R")) {
  source("S01_simulation_engine.R")
} else if (file.exists(file.path(METHODS_SCRIPTS, "S01_simulation_engine.R"))) {
  source(file.path(METHODS_SCRIPTS, "S01_simulation_engine.R"))
}


# ==============================================================================
# Core function: collect_run()
# ==============================================================================

#' Collect all performance metrics from a single pipeline run directory
#'
#' @param run_dir Pipeline run directory (containing data/ subdirectory)
#' @param run_id  Run identifier (for summary table)
#' @param mode    "auto" | "simulation" | "real"
#'                autoauto-detect whether ground_truth.rds exists
#' @return list: metrics (data.frame), funnel, stability, gap_details
collect_run <- function(run_dir, run_id = basename(run_dir), mode = "auto") {

  data_dir <- file.path(run_dir, "data")
  if (!dir.exists(data_dir)) {
    stop(sprintf("[S05] data/ directory not found: %s", data_dir))
  }

  methods_log("S05_COLLECT", sprintf("Collecting from: %s (id=%s)", run_dir, run_id))

  # ---- 0. Determine mode ----
  gt_file <- file.path(data_dir, "ground_truth.rds")
  has_gt <- file.exists(gt_file)
  if (mode == "auto") mode <- ifelse(has_gt, "simulation", "real")
  methods_log("S05_COLLECT", sprintf("Mode: %s (ground_truth: %s)", mode, has_gt))

  # ---- 1. Screening funnel layer data ----
  funnel <- collect_funnel(data_dir)

  # ---- 2. ML stability structure ----
  stability <- collect_stability(data_dir)

  # ---- 3. Gap analysis quality ----
  gap <- collect_gap_quality(data_dir)

  # ---- 4. PPI hub information ----
  ppi <- collect_ppi(data_dir)

  # ---- 5. Final candidate genes ----
  final <- collect_final_candidates(data_dir)

  # ---- 6. Performance metrics (simulation data only) ----
  perf <- NULL
  if (mode == "simulation" && has_gt) {
    perf <- evaluate_against_ground_truth(data_dir)
  }

  # ---- 7. Assemble metrics summary row (one row, for appending to master CSV) ----
  summary_row <- build_summary_row(
    run_id = run_id, mode = mode,
    funnel = funnel, stability = stability,
    gap = gap, ppi = ppi, final = final, perf = perf
  )

  # ---- 8. Package output ----
  result <- list(
    run_id    = run_id,
    run_dir   = run_dir,
    mode      = mode,
    summary   = summary_row,
    funnel    = funnel,
    stability = stability,
    gap       = gap,
    ppi       = ppi,
    final     = final,
    perf      = perf
  )

  # Save single run results
  out_file <- file.path(BENCH_DIR, sprintf("benchmark_results_%s.rds", run_id))
  saveRDS(result, out_file)
  methods_log("S05_COLLECT", sprintf("Saved: %s", basename(out_file)))

  # Append to master CSV
  append_to_master(summary_row)

  return(result)
}


# ==============================================================================
# Sub-collection functions
# ==============================================================================

#' Collect gene counts at each funnel layer
collect_funnel <- function(data_dir) {

  result <- list()

  # screening_funnel_data.csv（10_integration output）
  funnel_file <- file.path(data_dir, "screening_funnel_data.csv")
  if (file.exists(funnel_file)) {
    df <- read.csv(funnel_file, stringsAsFactors = FALSE)
    result$funnel_df <- df
    result$n_steps <- nrow(df)
    for (i in seq_len(nrow(df))) {
      key <- gsub("[^a-zA-Z0-9]", "_", tolower(df$step[i]))
      result[[key]] <- df$n[i]
    }
  }

  # Candidate pool details (09A output)
  pool_file <- file.path(data_dir, "candidate_pool.rds")
  if (file.exists(pool_file)) {
    pool <- readRDS(pool_file)
    result$n_masigpro_raw         <- pool$n_masigpro_raw
    result$n_masigpro_interaction <- pool$n_masigpro_interaction
    result$n_wgcna_key_genes      <- length(pool$wgcna_key_genes)
    result$n_pool_pre_effect      <- pool$n_pool_pre_effect
    result$n_pool_post_effect     <- pool$n_pool_post_effect
    result$n_candidate_pool       <- length(pool$candidate_pool)
    result$effect_filter_applied  <- pool$effect_filter_applied
  }

  return(result)
}


#' Collect ML stability selection structure
collect_stability <- function(data_dir) {

  result <- list()

  stab_file <- file.path(data_dir, "ml_stability_selection.rds")
  if (!file.exists(stab_file)) return(result)

  stab <- readRDS(stab_file)

  # Frequency distribution summary for each algorithm
  for (algo in c("lasso", "rf", "svm")) {
    freq_vec <- stab[[paste0(algo, "_freq")]]
    if (is.null(freq_vec)) next

    n_nonzero <- sum(freq_vec > 0)
    n_stable  <- sum(freq_vec >= stab$params$freq_threshold)
    max_freq  <- ifelse(length(freq_vec) > 0, max(freq_vec), 0)
    # Gini coefficient of freq distribution (concentration measure)
    freq_sorted <- sort(freq_vec[freq_vec > 0], decreasing = TRUE)
    gini <- if (length(freq_sorted) > 1) {
      n <- length(freq_sorted)
      sum((2 * seq_len(n) - n - 1) * freq_sorted) / (n * sum(freq_sorted))
    } else { 1.0 }

    result[[paste0(algo, "_n_nonzero")]] <- n_nonzero
    result[[paste0(algo, "_n_stable")]]  <- n_stable
    result[[paste0(algo, "_max_freq")]]  <- round(max_freq, 4)
    result[[paste0(algo, "_gini")]]      <- round(gini, 4)
  }

  # Selection matrices (for Nogueira metric calculation)
  result$has_selection_matrices <- all(c(
    "lasso_selection_matrix", "rf_selection_matrix", "svm_selection_matrix"
  ) %in% names(stab))

  if (result$has_selection_matrices) {
    result$lasso_sel_mat <- stab$lasso_selection_matrix
    result$rf_sel_mat    <- stab$rf_selection_matrix
    result$svm_sel_mat   <- stab$svm_selection_matrix
    result$n_bootstrap   <- stab$params$n_bootstrap

    # [v2] Pairwise Jaccard index for each algorithm's selection matrix
    # Ref: Real & Vargas (1996); Nogueira et al. (2018) recommend comparing
    #        adjusted (Nogueira) and unadjusted (Jaccard) stability metrics
    for (algo in c("lasso", "rf", "svm")) {
      mat_name <- paste0(algo, "_sel_mat")
      if (!is.null(result[[mat_name]])) {
        ji <- compute_pairwise_jaccard(result[[mat_name]])
        result[[paste0(algo, "_jaccard_mean")]] <- ji$mean
        result[[paste0(algo, "_jaccard_sd")]]   <- ji$sd
      }
    }
  }

  return(result)
}


#' Collect gap thresholding quality metrics
collect_gap_quality <- function(data_dir) {

  result <- list()

  gap_file <- file.path(data_dir, "ML_09D_gap_analysis.csv")
  if (!file.exists(gap_file)) return(result)

  gap_df <- read.csv(gap_file, stringsAsFactors = FALSE)
  result$gap_df <- gap_df

  # 09D complete results
  gap_rds <- file.path(data_dir, "ml_gap_union.rds")
  if (file.exists(gap_rds)) {
    gap_data <- readRDS(gap_rds)
    result$n_ml_final <- length(gap_data$final_gene_ids)
    result$final_gene_ids <- gap_data$final_gene_ids
    if (!is.null(gap_data$final_genes)) {
      result$final_genes_df <- gap_data$final_genes
    }
  }

  # Gap quality for each algorithm
  for (algo in unique(gap_df$algorithm)) {
    sub <- gap_df[gap_df$algorithm == algo, ]
    freqs <- sub$freq
    gaps  <- sub$gap_below[!is.na(sub$gap_below)]

    max_gap    <- ifelse(length(gaps) > 0, max(gaps), 0)
    max_freq   <- max(freqs)
    n_selected <- sum(freqs > 0.5)  # rough count

    # Gap-to-signal ratio: max_gap / max_freq
    gap_signal_ratio <- ifelse(max_freq > 0, max_gap / max_freq, 0)

    result[[paste0(algo, "_max_gap")]]   <- round(max_gap, 4)
    result[[paste0(algo, "_max_freq")]]  <- round(max_freq, 4)
    result[[paste0(algo, "_gap_ratio")]] <- round(gap_signal_ratio, 4)
    result[[paste0(algo, "_n_scanned")]] <- nrow(sub)
  }

  return(result)
}


#' Collect PPI hub information
collect_ppi <- function(data_dir) {

  result <- list(n_ppi_hubs = 0, hub_symbols = character(0))

  hub_file <- file.path(data_dir, "PPI_09F_hub_genes.csv")
  if (!file.exists(hub_file)) return(result)

  hub_df <- read.csv(hub_file, stringsAsFactors = FALSE)
  result$n_ppi_hubs   <- nrow(hub_df)
  result$hub_symbols  <- hub_df$symbol
  result$hub_df       <- hub_df

  # PPI complete results
  ppi_rds <- file.path(data_dir, "ppi_hub_selection.rds")
  if (file.exists(ppi_rds)) {
    ppi_data <- readRDS(ppi_rds)
    if (!is.null(ppi_data$network_stats)) {
      result$ppi_n_nodes <- ppi_data$network_stats$n_nodes
      result$ppi_n_edges <- ppi_data$network_stats$n_edges
      result$ppi_density <- ppi_data$network_stats$density
    }
  }

  return(result)
}


#' Collect final candidate genes
collect_final_candidates <- function(data_dir) {

  result <- list(n_final = 0)

  final_file <- file.path(data_dir, "Final_candidate_genes.csv")
  if (!file.exists(final_file)) return(result)

  df <- read.csv(final_file, stringsAsFactors = FALSE)
  result$n_final     <- nrow(df)
  result$final_df    <- df
  result$gene_ids    <- df$ensembl_id
  result$symbols     <- df$symbol

  # ML vs PPI source statistics
  result$n_from_ml   <- sum(df$in_ML == TRUE | df$in_ML == "TRUE", na.rm = TRUE)
  result$n_from_ppi  <- sum(df$is_WGCNA_hub == TRUE | df$is_WGCNA_hub == "TRUE", na.rm = TRUE)
  result$n_ml_and_ppi <- sum(
    (df$in_ML == TRUE | df$in_ML == "TRUE") &
    (df$is_WGCNA_hub == TRUE | df$is_WGCNA_hub == "TRUE"),
    na.rm = TRUE
  )
  result$n_ml_only  <- result$n_from_ml - result$n_ml_and_ppi
  result$n_ppi_only <- result$n_from_ppi - result$n_ml_and_ppi

  return(result)
}


# ==============================================================================
# Ground truth evaluation (simulation data only)
# ==============================================================================

#' Compare with simulation ground truth, calculate per-layer and final performance
evaluate_against_ground_truth <- function(data_dir) {

  gt <- readRDS(file.path(data_dir, "ground_truth.rds"))
  result <- list()

  # --- Final candidates vs TRUE_TEMPORAL ---
  final_file <- file.path(data_dir, "Final_candidate_genes.csv")
  if (file.exists(final_file)) {
    final_df <- read.csv(final_file, stringsAsFactors = FALSE)
    final_ids <- final_df$ensembl_id
    result$final <- evaluate_selection(final_ids, gt, "TRUE_TEMPORAL")
    result$final_any <- evaluate_selection(final_ids, gt, "TRUE_TEMPORAL")
    # Also compute a relaxed version: selected genes matching any signal (including TRUE_MAIN, TRUE_TIME)
    result$final$precision_any_signal <- sum(final_ids %in%
      gt$gene_id[gt$class != "NULL"]) / max(length(final_ids), 1)
  }

  # --- ML markers (09D) vs TRUE_TEMPORAL ---
  ml_rds <- file.path(data_dir, "ml_gap_union.rds")
  if (file.exists(ml_rds)) {
    ml_data <- readRDS(ml_rds)
    ml_ids <- ml_data$final_gene_ids
    result$ml <- evaluate_selection(ml_ids, gt, "TRUE_TEMPORAL")
  }

  # --- PPI hubs (09F) vs TRUE_TEMPORAL ---
  ppi_file <- file.path(data_dir, "PPI_09F_hub_genes.csv")
  if (file.exists(ppi_file)) {
    ppi_df <- read.csv(ppi_file, stringsAsFactors = FALSE)
    # PPI hubs use symbol, ground truth uses gene_id → need matching
    # In simulation data gene_id = symbol, so direct matching works
    ppi_ids <- ppi_df$ensembl_id
    if (all(is.na(ppi_ids)) || length(ppi_ids) == 0) {
      ppi_ids <- ppi_df$symbol  # fallback
    }
    result$ppi <- evaluate_selection(ppi_ids, gt, "TRUE_TEMPORAL")
  }

  # --- Candidate pool (09A) vs TRUE_TEMPORAL ---
  pool_file <- file.path(data_dir, "candidate_pool.rds")
  if (file.exists(pool_file)) {
    pool <- readRDS(pool_file)
    pool_ids <- pool$candidate_pool
    result$pool <- evaluate_selection(pool_ids, gt, "TRUE_TEMPORAL")
  }

  # --- maSigPro filtered vs TRUE_TEMPORAL ---
  masigpro_file <- file.path(data_dir, "masigpro_results.rds")
  if (file.exists(masigpro_file)) {
    masigpro <- readRDS(masigpro_file)
    if (!is.null(masigpro$gene_clusters)) {
      msig_ids <- unique(masigpro$gene_clusters$ensembl_id)
      result$masigpro <- evaluate_selection(msig_ids, gt, "TRUE_TEMPORAL")
    }
  }

  # --- Summarize per-layer recall decay (funnel tracking) ---
  layer_names <- c("masigpro", "pool", "ml", "ppi", "final")
  recall_track <- data.frame(
    layer = layer_names,
    recall = sapply(layer_names, function(ln) {
      if (!is.null(result[[ln]])) result[[ln]]$recall else NA
    }),
    precision = sapply(layer_names, function(ln) {
      if (!is.null(result[[ln]])) result[[ln]]$precision else NA
    }),
    n_selected = sapply(layer_names, function(ln) {
      if (!is.null(result[[ln]])) result[[ln]]$n_selected else NA
    }),
    stringsAsFactors = FALSE
  )
  result$recall_track <- recall_track

  return(result)
}


# ==============================================================================
# Summary row construction
# ==============================================================================

#' Build one summary row (for master CSV)
build_summary_row <- function(run_id, mode, funnel, stability,
                               gap, ppi, final, perf) {

  row <- data.frame(
    run_id              = run_id,
    mode                = mode,
    timestamp           = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    # Funnel
    n_transcriptome     = funnel$n_candidate_pool %||% NA,
    n_masigpro_raw      = funnel$n_masigpro_raw %||% NA,
    n_masigpro_interact = funnel$n_masigpro_interaction %||% NA,
    n_wgcna_key         = funnel$n_wgcna_key_genes %||% NA,
    n_pool_pre_effect   = funnel$n_pool_pre_effect %||% NA,
    n_pool_post_effect  = funnel$n_pool_post_effect %||% NA,
    n_candidate_pool    = funnel$n_candidate_pool %||% NA,
    effect_filter       = funnel$effect_filter_applied %||% NA,
    # Stability
    lasso_n_nonzero     = stability$lasso_n_nonzero %||% NA,
    lasso_n_stable      = stability$lasso_n_stable %||% NA,
    lasso_max_freq      = stability$lasso_max_freq %||% NA,
    lasso_gini          = stability$lasso_gini %||% NA,
    rf_n_nonzero        = stability$rf_n_nonzero %||% NA,
    rf_n_stable         = stability$rf_n_stable %||% NA,
    rf_max_freq         = stability$rf_max_freq %||% NA,
    rf_gini             = stability$rf_gini %||% NA,
    svm_n_nonzero       = stability$svm_n_nonzero %||% NA,
    svm_n_stable        = stability$svm_n_stable %||% NA,
    svm_max_freq        = stability$svm_max_freq %||% NA,
    svm_gini            = stability$svm_gini %||% NA,
    # Gap
    n_ml_final          = gap$n_ml_final %||% NA,
    # PPI
    n_ppi_hubs          = ppi$n_ppi_hubs %||% 0,
    ppi_n_nodes         = ppi$ppi_n_nodes %||% NA,
    ppi_n_edges         = ppi$ppi_n_edges %||% NA,
    # Final
    n_final             = final$n_final %||% 0,
    n_from_ml           = final$n_from_ml %||% 0,
    n_from_ppi          = final$n_from_ppi %||% 0,
    n_ml_and_ppi        = final$n_ml_and_ppi %||% 0,
    # [v2] Jaccard stability
    lasso_jaccard_mean  = stability$lasso_jaccard_mean %||% NA,
    lasso_jaccard_sd    = stability$lasso_jaccard_sd %||% NA,
    rf_jaccard_mean     = stability$rf_jaccard_mean %||% NA,
    rf_jaccard_sd       = stability$rf_jaccard_sd %||% NA,
    svm_jaccard_mean    = stability$svm_jaccard_mean %||% NA,
    svm_jaccard_sd      = stability$svm_jaccard_sd %||% NA,
    stringsAsFactors = FALSE
  )

  # Simulation-data-only columns
  if (!is.null(perf) && !is.null(perf$final)) {
    row$precision    <- perf$final$precision
    row$recall       <- perf$final$recall
    row$F1           <- perf$final$F1
    row$FDR          <- perf$final$FDR
    row$TP           <- perf$final$TP
    row$FP           <- perf$final$FP
    row$FN           <- perf$final$FN
    row$precision_any <- perf$final$precision_any_signal
  }

  return(row)
}

# ==============================================================================
# [v2] Pairwise Jaccard index
# ==============================================================================

#' Calculate pairwise Jaccard index from bootstrap selection matrix
#' @param sel_matrix  binary matrix (n_bootstrap × n_genes), 1=selected
#' @return list(mean, sd, n_pairs)
compute_pairwise_jaccard <- function(sel_matrix, freq_threshold = 0) {
  if (is.null(sel_matrix) || nrow(sel_matrix) < 2) {
    return(list(mean = NA_real_, sd = NA_real_, n_pairs = 0))
  }
  n_boot <- nrow(sel_matrix)
  if (freq_threshold > 0) {
    col_freqs <- colMeans(sel_matrix)
    keep <- col_freqs >= freq_threshold
    if (sum(keep) == 0) return(list(mean = 0, sd = 0, n_pairs = 0))
    sel_matrix <- sel_matrix[, keep, drop = FALSE]
  }
  # Sample up to 500 pairs for efficiency
  if (n_boot > 32) {
    n_pairs <- 500
    pairs <- matrix(0, n_pairs, 2)
    for (p in seq_len(n_pairs)) pairs[p, ] <- sample(n_boot, 2, replace = FALSE)
  } else {
    pairs <- t(combn(n_boot, 2))
  }
  jaccards <- numeric(nrow(pairs))
  for (p in seq_len(nrow(pairs))) {
    a <- which(sel_matrix[pairs[p, 1], ] == 1)
    b <- which(sel_matrix[pairs[p, 2], ] == 1)
    intersection <- length(intersect(a, b))
    union_size   <- length(union(a, b))
    jaccards[p] <- if (union_size > 0) intersection / union_size else 0
  }
  return(list(mean = round(mean(jaccards), 4), sd = round(sd(jaccards), 4),
              n_pairs = nrow(pairs)))
}


# NULL-coalescing operator
`%||%` <- function(x, y) if (is.null(x)) y else x


# ==============================================================================
# Master CSV management
# ==============================================================================

#' Append one row to benchmark_master.csv
append_to_master <- function(row_df) {

  master_file <- file.path(BENCH_DIR, "benchmark_master.csv")

  if (file.exists(master_file)) {
    existing <- read.csv(master_file, stringsAsFactors = FALSE)
    # If run_id already exists, replace
    existing <- existing[existing$run_id != row_df$run_id, , drop = FALSE]
    # Align column names (new columns may not exist in old data)
    all_cols <- union(colnames(existing), colnames(row_df))
    for (col in setdiff(all_cols, colnames(existing))) existing[[col]] <- NA
    for (col in setdiff(all_cols, colnames(row_df))) row_df[[col]] <- NA
    combined <- rbind(existing[, all_cols], row_df[, all_cols])
  } else {
    combined <- row_df
  }

  write.csv(combined, master_file, row.names = FALSE)
  methods_log("S05_COLLECT", sprintf("Master CSV updated: %d runs (%s)",
                                      nrow(combined), master_file))
}


# ==============================================================================
# Batch collection functions
# ==============================================================================

#' Batch collect all simulation pipeline results
#'
#' @param sim_base_dir Simulation benchmark root directory
#'        Expected structure: sim_base_dir/{scenario_name}/data/{pipeline outputs}
#'        i.e., directory after S01 generates simulation → S02 adapts → S04 runs pipeline
collect_all_simulations <- function(
    sim_base_dir = file.path(RUN_DIR, "simulations"),
    pattern = "^(low|medium|high)_") {

  if (!dir.exists(sim_base_dir)) {
    methods_log("S05_BATCH", sprintf("Directory not found: %s", sim_base_dir))
    return(invisible(NULL))
  }

  # Find all scenarios with data/ subdirectory
  all_dirs <- list.dirs(sim_base_dir, recursive = FALSE, full.names = TRUE)
  run_dirs <- all_dirs[sapply(all_dirs, function(d) {
    dir.exists(file.path(d, "data")) &&
    file.exists(file.path(d, "data", "ground_truth.rds"))
  })]

  if (length(run_dirs) == 0) {
    methods_log("S05_BATCH", "No simulation runs found with data/ground_truth.rds")
    return(invisible(NULL))
  }

  methods_log("S05_BATCH", sprintf("Found %d simulation runs", length(run_dirs)))

  results <- list()
  for (rd in run_dirs) {
    run_name <- basename(rd)
    tryCatch({
      results[[run_name]] <- collect_run(rd, run_id = run_name, mode = "simulation")
    }, error = function(e) {
      methods_log("S05_BATCH", sprintf("ERROR in %s: %s", run_name, e$message))
    })
  }

  methods_log("S05_BATCH", sprintf("Collected %d/%d runs successfully",
                                    length(results), length(run_dirs)))
  return(invisible(results))
}


#' Batch collect pipeline results for all GEO datasets
collect_all_geo <- function() {

  geo_runs <- list(
    GSE307424 = file.path(RUN_DIR, "GEO_GSE307424_Lung"),
    GSE197067 = file.path(RUN_DIR, "GEO_GSE197067_Tcell")
  )

  results <- list()
  for (gse in names(geo_runs)) {
    rd <- geo_runs[[gse]]
    if (dir.exists(file.path(rd, "data"))) {
      tryCatch({
        results[[gse]] <- collect_run(rd, run_id = gse, mode = "real")
      }, error = function(e) {
        methods_log("S05_GEO", sprintf("ERROR in %s: %s", gse, e$message))
      })
    } else {
      methods_log("S05_GEO", sprintf("Not found: %s", rd))
    }
  }

  return(invisible(results))
}


# ==============================================================================
# Diagnostic print
# ==============================================================================

#' Print diagnostic summary for single run
print_run_summary <- function(result) {

  cat("\n")
  cat("================================================================\n")
  cat(sprintf("  Run: %s (%s)\n", result$run_id, result$mode))
  cat("================================================================\n")

  # Funnel
  f <- result$funnel
  cat("\n  [FUNNEL]\n")
  if (!is.null(f$funnel_df)) {
    for (i in seq_len(nrow(f$funnel_df))) {
      cat(sprintf("    %s: %d\n", f$funnel_df$step[i], f$funnel_df$n[i]))
    }
  }

  # Stability
  s <- result$stability
  cat("\n  [STABILITY]\n")
  for (algo in c("lasso", "rf", "svm")) {
    cat(sprintf("    %s: %d nonzero, %d stable, max_freq=%.2f, gini=%.2f\n",
                toupper(algo),
                s[[paste0(algo, "_n_nonzero")]] %||% 0,
                s[[paste0(algo, "_n_stable")]] %||% 0,
                s[[paste0(algo, "_max_freq")]] %||% 0,
                s[[paste0(algo, "_gini")]] %||% 0))
  }

  # Gap
  cat(sprintf("\n  [GAP-UNION] ML markers: %d\n", result$gap$n_ml_final %||% 0))

  # PPI
  cat(sprintf("  [PPI] Hubs: %d\n", result$ppi$n_ppi_hubs %||% 0))

  # Final
  cat(sprintf("  [FINAL] Candidates: %d (ML=%d, PPI=%d, overlap=%d)\n",
              result$final$n_final,
              result$final$n_from_ml,
              result$final$n_from_ppi,
              result$final$n_ml_and_ppi))

  # Performance (simulation only)
  if (!is.null(result$perf) && !is.null(result$perf$final)) {
    p <- result$perf$final
    cat(sprintf("\n  [PERFORMANCE vs ground truth]\n"))
    cat(sprintf("    Precision: %.3f\n", p$precision))
    cat(sprintf("    Recall:    %.3f\n", p$recall))
    cat(sprintf("    F1:        %.3f\n", p$F1))
    cat(sprintf("    FDR:       %.3f\n", p$FDR))
    cat(sprintf("    TP=%d, FP=%d, FN=%d\n", p$TP, p$FP, p$FN))

    if (!is.null(result$perf$recall_track)) {
      cat("\n    [Recall by layer]\n")
      rt <- result$perf$recall_track
      for (i in seq_len(nrow(rt))) {
        if (!is.na(rt$recall[i])) {
          cat(sprintf("      %s: recall=%.3f, precision=%.3f (n=%d)\n",
                      rt$layer[i], rt$recall[i], rt$precision[i], rt$n_selected[i]))
        }
      }
    }
  }

  cat("\n================================================================\n")
}


# ==============================================================================
# Demo when running directly
# ==============================================================================

if (sys.nframe() == 0) {
  cat("\n")
  cat("================================================================\n")
  cat("  S05_benchmark_collector.R — Usage\n")
  cat("================================================================\n")
  cat("\n")
  cat("  # Collect single GEO dataset result:\n")
  cat("  res <- collect_run('METHODS_BASE')\n")
  cat("  print_run_summary(res)\n")
  cat("\n")
  cat("  # Batch collect all GEO:\n")
  cat("  geo_results <- collect_all_geo()\n")
  cat("\n")
  cat("  # Batch collect all simulation data:\n")
  cat("  sim_results <- collect_all_simulations()\n")
  cat("\n")
  cat("  # View master summary table:\n")
  cat("  master <- read.csv(file.path(BENCH_DIR, 'benchmark_master.csv'))\n")
  cat("================================================================\n")
}
