#!/usr/bin/env Rscript
# ==============================================================================
# S01_simulation_engine.R — RNA-seqSimulation data generation (v2 — BIB)
#
# ：Knownground truthRNA-seq countsData，
# EvaluationMETI-FS pipelinesensitivity/specificity/FDR。
#
# v2 (v1):
# [1.1] DispersionTrue: DESeq2-style α(μ) = α₀ + α₁/μ Parameters
# : Love et al. 2014 (Genome Biology); 
#               Soneson (2014, compcodeR, Bioinformatics);
#               Vieth et al. 2017 (powsimR, Bioinformatics)
# [1.2] Expressionblock: TRUE_TEMPORALGeneExpressionModule
# : Stabl (Hédou et al. 2024, Nat Biotechnol) Testcorrelated features;
# Langfelder & Horvath 2008 (BMC Bioinformatics) WGCNAModule
# [1.1+] library size (biological + technical variation)
#
# （）：
# - NBcounts: Love et al. 2014 (DESeq2), Nueda et al. 2014 (maSigPro)
# - Expression: Spies et al. 2019 (Briefings in Bioinformatics) 
# - DispersionParameters: DESeq2Parameters α(μ) = α₀ + α₁/μ
# Soneson 2014 / powsimR (Vieth 2017):
#       bulk RNA-seq: α₀ ≈ 0.01-0.1, α₁ ≈ 1-10
# - : Stabl (Hédou et al. 2024, Nature Biotechnology) Data
# - Expression: Choleskyblock correlation
#
# GeneClassification（ground truth）：
# - Class 1: TRUE_TEMPORAL — Process×TimeSignificant (pipeline)
# - Class 2: TRUE_MAIN — ProcessSignificantTime
# - Class 3: TRUE_TIME — Time，
# - Class 4: NULL — DifferentialGene
#
# Output：
#   - counts_matrix.rds, tpm_matrix.rds, sample_info.rds
# - ground_truth.rds ( module_id )
# - simulation_params.rds ( dispersion_source, correlation_structure)
# ==============================================================================

# ---- 0. LoadConfigure ----
if (file.exists("S_config.R")) {
  source("S_config.R")
} else if (file.exists(file.path(file.path(METHODS_BASE, "scripts"), "S_config.R"))) {
  source(file.path(file.path(METHODS_BASE, "scripts"), "S_config.R"))
}

# ---- 0a. ----
suppressPackageStartupMessages({
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
  library(MASS)
})

# ==============================================================================
# Function: DESeq2-style dispersion
# ==============================================================================

#' DESeq2-style dispersionParameters
#' 
#' DESeq2mean-dispersion: α(μ) = α_intercept + α_slope / μ
#' ExpressionGenedispersion（），ExpressionGenedispersion（noisy）
#' gene-to-gene (log-normal noise)
#'
#' :
#'   Love et al. 2014, Genome Biology 15:550 (DESeq2)
#'   Soneson 2014, Bioinformatics 30(18):2670 (compcodeR benchmark)
#'   Vieth et al. 2017, Bioinformatics 33(21):3486 (powsimR)
#'
#' @param base_means ，GeneExpression
#' @param alpha_intercept dispersion (ExpressionGene), 0.01-0.1
#' @param alpha_slope dispersion (ExpressionGenedispersion), 1-10
#' @param noise_sd gene-to-gene log-normalSD, 0.3-0.8
#' @return ，Genedispersion
generate_deseq2_dispersions <- function(base_means, 
                                         alpha_intercept = 0.05, 
                                         alpha_slope = 4.0,
                                         noise_sd = 0.5) {
 # : α(μ) = α₀ + α₁/μ
  trend <- alpha_intercept + alpha_slope / base_means
  
 # gene-to-gene (DESeq2dispersionlog-normal)
  # log(α_gene) = log(α_trend) + ε, ε ~ N(0, σ²)
  n <- length(base_means)
  log_dispersions <- log(trend) + rnorm(n, mean = 0, sd = noise_sd)
  dispersions <- exp(log_dispersions)
  
 # : [0.001, 5.0]
  dispersions <- pmax(dispersions, 0.001)
  dispersions <- pmin(dispersions, 5.0)
  
  return(dispersions)
}


# ==============================================================================
# Function: Expressionblock
# ==============================================================================

#' log-spaceGeneblock correlation
#'
#' nGene × mSamplelog2-expression，
#' Module。
#'
#' : 
#' 1. Gene n_modules block
#' 2. blocklatent factor
#'   3. gene_expr = sqrt(ρ) * latent + sqrt(1-ρ) * independent
#' Module ρ，Module 0
#'
#' :
#'   Langfelder & Horvath 2008, BMC Bioinformatics 9:559 (WGCNA)
#'   Hédou et al. 2024, Nat Biotechnol 42:1581 (Stabl, R≈0.5 correlated)
#'
#' @param lfc_matrix Gene×Time log2FC ()
#' @param gene_ids GeneID
#' @param n_modules Module
#' @param rho_within Module (0-1), 0.4
#' @param n_samples Sample (latent factors)
#' @return list(lfc_matrix_corr, module_assignments)
inject_block_correlation <- function(lfc_matrix, gene_ids, 
                                      n_modules = 4, 
                                      rho_within = 0.4,
                                      n_timepoints = ncol(lfc_matrix)) {
  
  n_genes <- length(gene_ids)
  
  if (n_genes < n_modules) {
 # Gene，Module
    return(list(
      lfc_matrix = lfc_matrix,
      module_ids = rep("M1", n_genes)
    ))
  }
  
 # 1. GeneModule
  module_ids <- paste0("M", rep(1:n_modules, length.out = n_genes))
  
 # 2. Module，pattern
 # ModuleGeneTime
  for (m in 1:n_modules) {
    idx <- which(module_ids == paste0("M", m))
    if (length(idx) < 2) next
    
 # Modulelatent temporal perturbation
 # Time，ModuleGene
    latent_perturbation <- rnorm(n_timepoints, mean = 0, sd = 0.3)
    
    for (g in idx) {
      # gene_lfc = sqrt(1-ρ) * original_lfc + sqrt(ρ) * module_shared
 # ，Module
      original <- lfc_matrix[g, ]
      shared <- latent_perturbation * sign(mean(original[original != 0]) + 0.01)
      
 # : Raw，Module
      lfc_matrix[g, ] <- sqrt(1 - rho_within) * original + 
                          sqrt(rho_within) * shared * max(abs(original))
    }
  }
  
  return(list(
    lfc_matrix = lfc_matrix,
    module_ids = module_ids
  ))
}


# ==============================================================================
# Function: generate_simulation() — v2 
# ==============================================================================

generate_simulation <- function(
 # ---- ExperimentParameters ----
    n_genes       = 13000,
    n_timepoints  = 4,
    time_values   = NULL,
    n_reps_ind    = 4,
    n_reps_ctrl   = 3,
    
 # ---- Parameters ----
    n_true_temporal  = 20,
    n_true_main      = 50,
    n_true_timeonly   = 100,
    
 # ---- Parameters ----
    snr = "medium",
    lfc_range = NULL,
    
 # ---- GeneExpressionParameters ----
    base_mean_range = c(50, 5000),
 dispersion_range = c(0.1, 0.5), # dispersion_source="linear" 
    
    # ---- [v2 NEW] Dispersion source ----
 # "deseq2": DESeq2-style α(μ) = α₀ + α₁/μ + noise (Recommended)
 # "linear": v1 (Map, )
    dispersion_source = "deseq2",
    
 # DESeq2 dispersionParameters ( dispersion_source="deseq2" )
 # bulk RNA-seq: α₀=0.05, α₁=4.0, noise_sd=0.5
 # : powsimR (Vieth 2017) ENCODE/GTExParameters
    disp_intercept = 0.05,
    disp_slope     = 4.0,
    disp_noise_sd  = 0.5,
    
    # ---- [v2 NEW] Correlation structure ----
 # "independent": Gene (v1, )
 # "block": TRUE_TEMPORALGeneExpressionModule
    correlation_structure = "independent",
 n_modules = 4, # Module ( correlation_structure="block")
 rho_within = 0.4, # Module ( correlation_structure="block")
    
 # ---- Parameters ----
    temporal_patterns = c("sustained_up", "sustained_down", 
                          "early_peak", "late_onset",
                          "transient_up", "transient_down"),
    
 # ---- Output ----
    output_dir = NULL,
    seed = 42,
    verbose = TRUE
) {
  
  set.seed(seed)
  
 # ---- 1. Parameters ----
  if (is.null(time_values)) {
    time_values <- seq_len(n_timepoints)
  }
  stopifnot(length(time_values) == n_timepoints)
  
 # -> log2FC
  if (is.null(lfc_range)) {
    lfc_range <- switch(snr,
      "low"    = c(0.5, 1.5),
      "medium" = c(1.0, 3.0),
      "high"   = c(2.0, 5.0),
      stop("snr must be 'low', 'medium', or 'high'")
    )
  }
  
  n_samples_per_tp <- n_reps_ind + n_reps_ctrl
  n_total_samples <- n_samples_per_tp * n_timepoints
  
  if (verbose) {
    cat(sprintf("[SIM] Generating simulation (v2):\n"))
    cat(sprintf("  Genes: %d (true_temporal=%d, true_main=%d, true_time=%d, null=%d)\n",
                n_genes, n_true_temporal, n_true_main, n_true_timeonly,
                n_genes - n_true_temporal - n_true_main - n_true_timeonly))
    cat(sprintf("  Design: %d timepoints × (Ind=%d + Ctrl=%d) = %d samples\n",
                n_timepoints, n_reps_ind, n_reps_ctrl, n_total_samples))
    cat(sprintf("  SNR: %s (lfc_range: [%.1f, %.1f])\n", snr, lfc_range[1], lfc_range[2]))
    cat(sprintf("  Dispersion: %s", dispersion_source))
    if (dispersion_source == "deseq2") {
      cat(sprintf(" (α₀=%.3f, α₁=%.1f, σ=%.2f)", disp_intercept, disp_slope, disp_noise_sd))
    }
    cat("\n")
    cat(sprintf("  Correlation: %s", correlation_structure))
    if (correlation_structure == "block") {
      cat(sprintf(" (%d modules, ρ=%.2f)", n_modules, rho_within))
    }
    cat("\n")
    cat(sprintf("  Seed: %d\n", seed))
  }
  
 # ---- 2. SampleInfo ----
  sample_info <- data.frame(
    sample_id  = character(0),
    Treatment  = character(0),
    Time       = character(0),
    time_num   = numeric(0),
    replicate  = integer(0),
    stringsAsFactors = FALSE
  )
  
  for (t in seq_len(n_timepoints)) {
    tp_label <- paste0("T", t)
    for (r in seq_len(n_reps_ind)) {
      sample_info <- rbind(sample_info, data.frame(
        sample_id = sprintf("Ind_T%d_R%d", t, r),
        Treatment = "Induced",
        Time = tp_label,
        time_num = time_values[t],
        replicate = r,
        stringsAsFactors = FALSE
      ))
    }
    for (r in seq_len(n_reps_ctrl)) {
      sample_info <- rbind(sample_info, data.frame(
        sample_id = sprintf("Ctrl_T%d_R%d", t, r),
        Treatment = "Control",
        Time = tp_label,
        time_num = time_values[t],
        replicate = r,
        stringsAsFactors = FALSE
      ))
    }
  }
  rownames(sample_info) <- sample_info$sample_id
  
 # ---- 3. GeneParameters ----
  gene_ids <- sprintf("GENE_%05d", seq_len(n_genes))
  
 # Expression（log-normal）
  base_means <- exp(runif(n_genes, 
                          log(base_mean_range[1]), 
                          log(base_mean_range[2])))
  
 # [v2] DispersionParameters — sourceSelection
  if (dispersion_source == "deseq2") {
    dispersions <- generate_deseq2_dispersions(
      base_means, 
      alpha_intercept = disp_intercept,
      alpha_slope     = disp_slope,
      noise_sd        = disp_noise_sd
    )
    if (verbose) {
      cat(sprintf("  Dispersion range: [%.4f, %.4f] (median=%.4f)\n",
                  min(dispersions), max(dispersions), median(dispersions)))
    }
  } else {
 # v1: Map
    dispersions <- dispersion_range[2] - 
      (dispersion_range[2] - dispersion_range[1]) * 
      (log(base_means) - log(base_mean_range[1])) / 
      (log(base_mean_range[2]) - log(base_mean_range[1]))
    dispersions <- pmax(dispersions, dispersion_range[1])
  }
  
 # ---- 4. Gene ----
  n_null <- n_genes - n_true_temporal - n_true_main - n_true_timeonly
  stopifnot(n_null > 0)
  
  gene_class <- c(
    rep("TRUE_TEMPORAL", n_true_temporal),
    rep("TRUE_MAIN", n_true_main),
    rep("TRUE_TIME", n_true_timeonly),
    rep("NULL", n_null)
  )
  shuffle_idx <- sample(n_genes)
  gene_class <- gene_class[shuffle_idx]
  
 # ---- 5. Expression ----
  
  generate_temporal_profile <- function(pattern, n_tp, lfc_max) {
    t_norm <- seq(0, 1, length.out = n_tp)
    
    profile <- switch(pattern,
      "sustained_up" = {
        lfc_max * (1 - exp(-5 * t_norm))
      },
      "sustained_down" = {
        -lfc_max * (1 - exp(-5 * t_norm))
      },
      "early_peak" = {
        lfc_max * t_norm * exp(1 - t_norm) * exp(1)
      },
      "late_onset" = {
        lfc_max * pmax(0, (t_norm - 0.5) * 2)^2
      },
      "transient_up" = {
        lfc_max * dnorm(t_norm, mean = 0.4, sd = 0.15) / 
          dnorm(0.4, mean = 0.4, sd = 0.15)
      },
      "transient_down" = {
        -lfc_max * dnorm(t_norm, mean = 0.5, sd = 0.2) / 
          dnorm(0.5, mean = 0.5, sd = 0.2)
      }
    )
    return(profile)
  }
  
  gene_params <- data.frame(
    gene_id = gene_ids,
    class = gene_class,
    base_mean = base_means,
    dispersion = dispersions,
    lfc_max = NA_real_,
    pattern = NA_character_,
    stringsAsFactors = FALSE
  )
  
  lfc_temporal_matrix <- matrix(0, nrow = n_genes, ncol = n_timepoints)
  rownames(lfc_temporal_matrix) <- gene_ids
  colnames(lfc_temporal_matrix) <- paste0("T", seq_len(n_timepoints))
  
  for (i in seq_len(n_genes)) {
    cls <- gene_class[i]
    
    if (cls == "TRUE_TEMPORAL") {
      lfc_max <- runif(1, lfc_range[1], lfc_range[2])
      if (runif(1) < 0.5) lfc_max <- -lfc_max
      
      pattern <- sample(temporal_patterns, 1)
      lfc_profile <- generate_temporal_profile(
        pattern = gsub("_down$", "_up", pattern),
        n_tp = n_timepoints, 
        lfc_max = abs(lfc_max)
      )
      if (lfc_max < 0) lfc_profile <- -lfc_profile
      
      lfc_temporal_matrix[i, ] <- lfc_profile
      gene_params$lfc_max[i] <- lfc_max
      gene_params$pattern[i] <- pattern
      
    } else if (cls == "TRUE_MAIN") {
      lfc_const <- runif(1, lfc_range[1] * 0.5, lfc_range[2] * 0.5)
      if (runif(1) < 0.5) lfc_const <- -lfc_const
      
      lfc_temporal_matrix[i, ] <- rep(lfc_const, n_timepoints)
      gene_params$lfc_max[i] <- lfc_const
      gene_params$pattern[i] <- "constant"
      
    } else if (cls == "TRUE_TIME") {
      lfc_temporal_matrix[i, ] <- rep(0, n_timepoints)
      gene_params$lfc_max[i] <- 0
      time_lfc <- runif(1, 0.5, 2.0) * (1 - exp(-3 * seq(0, 1, length.out = n_timepoints)))
      if (runif(1) < 0.5) time_lfc <- -time_lfc
      gene_params$pattern[i] <- "time_shared"
      attr(gene_params, "time_effects") <- if (is.null(attr(gene_params, "time_effects"))) {
        list()
      } else {
        attr(gene_params, "time_effects")
      }
      attr(gene_params, "time_effects")[[gene_ids[i]]] <- time_lfc
      
    } else {
      lfc_temporal_matrix[i, ] <- rep(0, n_timepoints)
      gene_params$lfc_max[i] <- 0
      gene_params$pattern[i] <- "null"
    }
  }
  
 # ---- 5b. [v2] Expressionblock ----
  module_ids <- rep(NA_character_, n_genes)
  
  if (correlation_structure == "block") {
    temporal_idx <- which(gene_class == "TRUE_TEMPORAL")
    
    if (length(temporal_idx) >= n_modules) {
      block_result <- inject_block_correlation(
        lfc_matrix   = lfc_temporal_matrix[temporal_idx, , drop = FALSE],
        gene_ids     = gene_ids[temporal_idx],
        n_modules    = n_modules,
        rho_within   = rho_within,
        n_timepoints = n_timepoints
      )
      lfc_temporal_matrix[temporal_idx, ] <- block_result$lfc_matrix
      module_ids[temporal_idx] <- block_result$module_ids
      
      if (verbose) {
        cat(sprintf("  Block correlation injected: %d temporal genes → %d modules (ρ=%.2f)\n",
                    length(temporal_idx), n_modules, rho_within))
      }
    }
  }
  
 # ---- 6. counts ----
  counts_matrix <- matrix(0L, nrow = n_genes, ncol = n_total_samples)
  rownames(counts_matrix) <- gene_ids
  colnames(counts_matrix) <- sample_info$sample_id
  
  time_effects_list <- attr(gene_params, "time_effects")
  
 # [v2] Samplelibrary size (Differential)
 # log-normal: Median1.0, SampleCV20%
  lib_size_factors <- exp(rnorm(n_total_samples, mean = 0, sd = 0.2))
  
  for (i in seq_len(n_genes)) {
    base_mu <- base_means[i]
    disp <- dispersions[i]
    cls <- gene_class[i]
    
    for (j in seq_len(n_total_samples)) {
      trt <- sample_info$Treatment[j]
      tp_idx <- match(sample_info$Time[j], paste0("T", seq_len(n_timepoints)))
      
 mu <- base_mu * lib_size_factors[j] # [v2] lib size variation
      
 # Time (Class 3)
      if (cls == "TRUE_TIME" && !is.null(time_effects_list[[gene_ids[i]]])) {
        time_eff <- time_effects_list[[gene_ids[i]]]
        mu <- mu * 2^(time_eff[tp_idx])
      }
      
 # Process (Induced)
      if (trt == "Induced") {
        lfc <- lfc_temporal_matrix[i, tp_idx]
        mu <- mu * 2^(lfc)
      }
      
      mu <- max(mu, 1)
      
 # NB: mean = mu, size = 1/dispersion
      size_param <- 1 / disp
      counts_matrix[i, j] <- rnbinom(1, mu = mu, size = size_param)
    }
  }
  
 # ---- 7. TPM ----
  gene_lengths <- rep(2000, n_genes)
  rpk_matrix <- counts_matrix / (gene_lengths / 1000)
  tpm_matrix <- apply(rpk_matrix, 2, function(x) x / sum(x) * 1e6)
  
 # ---- 8. ground truth ----
  ground_truth <- data.frame(
    gene_id = gene_ids,
    class = gene_class,
    is_true_marker = gene_class == "TRUE_TEMPORAL",
    is_any_signal = gene_class != "NULL",
    base_mean = base_means,
    dispersion = dispersions,
    lfc_max = gene_params$lfc_max,
    pattern = gene_params$pattern,
 module_id = module_ids, # [v2] Module
    stringsAsFactors = FALSE
  )
  
 # ---- 9. Output ----
  simulation <- list(
    counts = counts_matrix,
    tpm = tpm_matrix,
    sample_info = sample_info,
    ground_truth = ground_truth,
    lfc_matrix = lfc_temporal_matrix,
    params = list(
      n_genes = n_genes,
      n_timepoints = n_timepoints,
      time_values = time_values,
      n_reps_ind = n_reps_ind,
      n_reps_ctrl = n_reps_ctrl,
      n_total_samples = n_total_samples,
      n_true_temporal = n_true_temporal,
      n_true_main = n_true_main,
      n_true_timeonly = n_true_timeonly,
      n_null = n_null,
      snr = snr,
      lfc_range = lfc_range,
      base_mean_range = base_mean_range,
      dispersion_source = dispersion_source,          # [v2]
      correlation_structure = correlation_structure,   # [v2]
      n_modules = if (correlation_structure == "block") n_modules else NA,
      rho_within = if (correlation_structure == "block") rho_within else NA,
      seed = seed
    )
  )
  
  if (verbose) {
    cat(sprintf("\n[SIM] Generated %d genes × %d samples\n", n_genes, n_total_samples))
    cat(sprintf("  Counts range: [%d, %d]\n", min(counts_matrix), max(counts_matrix)))
    cat(sprintf("  Median library size: %.0f\n", median(colSums(counts_matrix))))
    cat(sprintf("  True temporal markers: %d (%.1f%%)\n", 
                n_true_temporal, 100 * n_true_temporal / n_genes))
    
    true_lfc <- ground_truth$lfc_max[ground_truth$is_true_marker]
    if (length(true_lfc) > 0) {
      cat(sprintf("  True marker |lfc| range: [%.2f, %.2f]\n",
                  min(abs(true_lfc)), max(abs(true_lfc))))
    }
    
    if (n_true_temporal > 0) {
      pat_tab <- table(ground_truth$pattern[ground_truth$is_true_marker])
      cat("  Temporal patterns: ")
      cat(paste(sprintf("%s=%d", names(pat_tab), pat_tab), collapse = ", "))
      cat("\n")
    }
    
 # [v2] ModuleInfo
    if (correlation_structure == "block" && any(!is.na(module_ids))) {
      mod_tab <- table(module_ids[!is.na(module_ids)])
      cat("  Modules: ")
      cat(paste(sprintf("%s=%d", names(mod_tab), mod_tab), collapse = ", "))
      cat("\n")
    }
    
 # [v2] Dispersion
    if (dispersion_source == "deseq2") {
      cat(sprintf("  Dispersion (DESeq2-style): median=%.4f, IQR=[%.4f, %.4f]\n",
                  median(dispersions), quantile(dispersions, 0.25), quantile(dispersions, 0.75)))
    }
  }
  
 # SaveFile
  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    saveRDS(simulation$counts, file.path(output_dir, "counts_matrix.rds"))
    saveRDS(simulation$tpm, file.path(output_dir, "tpm_matrix.rds"))
    saveRDS(simulation$sample_info, file.path(output_dir, "sample_info.rds"))
    saveRDS(simulation$ground_truth, file.path(output_dir, "ground_truth.rds"))
    saveRDS(simulation$params, file.path(output_dir, "simulation_params.rds"))
    saveRDS(simulation, file.path(output_dir, "simulation_full.rds"))
    if (verbose) cat(sprintf("  [SAVED] All files to %s/\n", output_dir))
  }
  
  return(simulation)
}


# ==============================================================================
# ：EvaluationParameters
# ==============================================================================

#' benchmark (v2 )
#' 
#' Parameters: 3 SNR × 3 Sample × 3 markerDensity = 27 × 5Replicate = 135
#' [v2] : 3 SNR × block correlation × mediumSample × mediumDensity = 3 × 5Replicate = 15
#' Total: 150 (15，Calculate)

generate_benchmark_scenarios <- function(
    output_base_dir = if (exists("SIM_DIR")) file.path(SIM_DIR, "benchmark") else "simulation_benchmark",
    n_repeats = 5,
 include_correlation = TRUE, # [v2] 
    verbose = TRUE
) {
  
  snr_levels <- c("low", "medium", "high")
  
  sample_configs <- list(
    small  = list(n_reps_ind = 2, n_reps_ctrl = 2),
    medium = list(n_reps_ind = 4, n_reps_ctrl = 3),
    large  = list(n_reps_ind = 5, n_reps_ctrl = 5)
  )
  
  marker_configs <- list(
    sparse = list(n_true_temporal = 10, n_true_main = 30, n_true_timeonly = 60),
    medium = list(n_true_temporal = 20, n_true_main = 50, n_true_timeonly = 100),
    dense  = list(n_true_temporal = 50, n_true_main = 100, n_true_timeonly = 200)
  )
  
 # ---- : 27 × 5 = 135 ----
  scenarios <- expand.grid(
    snr = snr_levels,
    sample_size = names(sample_configs),
    marker_density = names(marker_configs),
    repeat_id = seq_len(n_repeats),
    stringsAsFactors = FALSE
  )
  scenarios$correlation <- "independent"
  
 # ---- [v2] : 3 SNR × 1 × 1 × 5 = 15 ----
  if (include_correlation) {
    corr_scenarios <- expand.grid(
      snr = snr_levels,
 sample_size = "medium", # mediumSample
 marker_density = "medium", # medium markerDensity
      repeat_id = seq_len(n_repeats),
      stringsAsFactors = FALSE
    )
    corr_scenarios$correlation <- "block"
    scenarios <- rbind(scenarios, corr_scenarios)
  }
  
  n_scenarios <- nrow(scenarios)
  if (verbose) {
    cat(sprintf("[BENCHMARK v2] Generating %d simulation scenarios\n", n_scenarios))
    cat(sprintf("  Core: 27 scenarios × %d repeats = %d\n", 
                n_repeats, 27 * n_repeats))
    if (include_correlation) {
      cat(sprintf("  Correlation: 3 scenarios × %d repeats = %d\n",
                  n_repeats, 3 * n_repeats))
    }
    cat(sprintf("  Output: %s/\n\n", output_base_dir))
  }
  
  dir.create(output_base_dir, recursive = TRUE, showWarnings = FALSE)
  
  scenario_log <- list()
  
  for (i in seq_len(n_scenarios)) {
    sc <- scenarios[i, ]
    
 # : correlation
    if (sc$correlation == "block") {
      sc_name <- sprintf("%s_%s_%s_block_rep%d", 
                         sc$snr, sc$sample_size, sc$marker_density, sc$repeat_id)
    } else {
      sc_name <- sprintf("%s_%s_%s_rep%d", 
                         sc$snr, sc$sample_size, sc$marker_density, sc$repeat_id)
    }
    
    sc_dir <- file.path(output_base_dir, sc_name)
    sc_sample <- sample_configs[[sc$sample_size]]
    sc_marker <- marker_configs[[sc$marker_density]]
    
    if (verbose && i %% 10 == 1) {
      cat(sprintf("  [%d/%d] %s ...\n", i, n_scenarios, sc_name))
    }
    
    sim <- generate_simulation(
      n_genes = 13000,
      n_timepoints = 4,
      time_values = c(4, 7, 14, 21),
      n_reps_ind = sc_sample$n_reps_ind,
      n_reps_ctrl = sc_sample$n_reps_ctrl,
      n_true_temporal = sc_marker$n_true_temporal,
      n_true_main = sc_marker$n_true_main,
      n_true_timeonly = sc_marker$n_true_timeonly,
      snr = sc$snr,
 dispersion_source = "deseq2", # [v2] DESeq2-style
      correlation_structure = sc$correlation,             # [v2]
      n_modules = 4,
      rho_within = 0.4,
      output_dir = sc_dir,
      seed = i * 1000 + sc$repeat_id,
      verbose = FALSE
    )
    
    scenario_log[[i]] <- data.frame(
      scenario_id = i,
      scenario_name = sc_name,
      snr = sc$snr,
      sample_size = sc$sample_size,
      marker_density = sc$marker_density,
      correlation = sc$correlation,                       # [v2]
      repeat_id = sc$repeat_id,
      n_samples = sim$params$n_total_samples,
      n_true_temporal = sim$params$n_true_temporal,
      median_lib_size = median(colSums(sim$counts)),
      median_dispersion = median(sim$ground_truth$dispersion),  # [v2]
      stringsAsFactors = FALSE
    )
  }
  
  scenario_log_df <- do.call(rbind, scenario_log)
  write.csv(scenario_log_df, file.path(output_base_dir, "scenario_log.csv"), row.names = FALSE)
  
  if (verbose) {
    cat(sprintf("\n[BENCHMARK v2 DONE] %d scenarios generated\n", n_scenarios))
    cat(sprintf("  Log saved: %s/scenario_log.csv\n", output_base_dir))
    cat(sprintf("  Dispersion: DESeq2-style (α₀=0.05, α₁=4.0)\n"))
    if (include_correlation) {
      cat(sprintf("  Correlation scenarios: block (4 modules, ρ=0.4)\n"))
    }
  }
  
  return(scenario_log_df)
}


# ==============================================================================
# EvaluationFunction: pipelineOutputground truth ()
# ==============================================================================

evaluate_selection <- function(selected_genes, ground_truth, 
                               target_class = "TRUE_TEMPORAL") {
  
  true_positives <- ground_truth$gene_id[ground_truth$class == target_class]
  all_positives <- ground_truth$gene_id[ground_truth$class != "NULL"]
  
  tp <- sum(selected_genes %in% true_positives)
  fp <- sum(!selected_genes %in% true_positives)
  fn <- sum(!true_positives %in% selected_genes)
  tn <- sum(!ground_truth$gene_id %in% c(selected_genes, true_positives))
  
  precision <- ifelse(tp + fp > 0, tp / (tp + fp), 0)
  recall    <- ifelse(tp + fn > 0, tp / (tp + fn), 0)
  f1        <- ifelse(precision + recall > 0, 
                      2 * precision * recall / (precision + recall), 0)
  fdr       <- ifelse(tp + fp > 0, fp / (tp + fp), 0)
  
  tp_any <- sum(selected_genes %in% all_positives)
  precision_any <- ifelse(length(selected_genes) > 0, 
                          tp_any / length(selected_genes), 0)
  
  return(data.frame(
    n_selected = length(selected_genes),
    n_true = length(true_positives),
    TP = tp, FP = fp, FN = fn,
    precision = round(precision, 4),
    recall = round(recall, 4),
    F1 = round(f1, 4),
    FDR = round(fdr, 4),
    precision_any_signal = round(precision_any, 4),
    stringsAsFactors = FALSE
  ))
}


# ==============================================================================
# Demo: Run
# ==============================================================================

if (sys.nframe() == 0) {
  cat("\n========== Running demo simulation (v2) ==========\n\n")
  
  # Demo 1: DESeq2 dispersion + independent
  cat("--- Demo 1: DESeq2 dispersion, independent genes ---\n")
  sim1 <- generate_simulation(
    n_genes = 5000, n_timepoints = 4, time_values = c(4, 7, 14, 21),
    n_reps_ind = 4, n_reps_ctrl = 3,
    n_true_temporal = 15, n_true_main = 40, n_true_timeonly = 80,
    snr = "medium",
    dispersion_source = "deseq2",
    correlation_structure = "independent",
    output_dir = if (exists("SIM_DIR")) file.path(SIM_DIR, "demo_v2_indep") else "demo_v2_indep",
    seed = 42
  )
  
  cat("\n--- Demo 2: DESeq2 dispersion, block correlation ---\n")
  sim2 <- generate_simulation(
    n_genes = 5000, n_timepoints = 4, time_values = c(4, 7, 14, 21),
    n_reps_ind = 4, n_reps_ctrl = 3,
    n_true_temporal = 15, n_true_main = 40, n_true_timeonly = 80,
    snr = "medium",
    dispersion_source = "deseq2",
    correlation_structure = "block",
    n_modules = 4, rho_within = 0.4,
    output_dir = if (exists("SIM_DIR")) file.path(SIM_DIR, "demo_v2_block") else "demo_v2_block",
    seed = 42
  )
  
  cat("\n--- Ground truth comparison ---\n")
  cat("Independent:\n")
  print(table(sim1$ground_truth$class))
  cat("\nBlock correlated:\n")
  print(table(sim2$ground_truth$class))
  
 # VerifyblockModule
  if (any(!is.na(sim2$ground_truth$module_id))) {
    temporal_genes <- sim2$ground_truth$gene_id[sim2$ground_truth$is_true_marker]
    tpm_temporal <- sim2$tpm[temporal_genes, ]
    log_tpm <- log2(tpm_temporal + 1)
    cor_mat <- cor(t(log_tpm))
    
    modules <- sim2$ground_truth$module_id[sim2$ground_truth$is_true_marker]
    within_cors <- c()
    between_cors <- c()
    for (ii in 1:(length(modules)-1)) {
      for (jj in (ii+1):length(modules)) {
        if (modules[ii] == modules[jj]) {
          within_cors <- c(within_cors, cor_mat[ii, jj])
        } else {
          between_cors <- c(between_cors, cor_mat[ii, jj])
        }
      }
    }
    cat(sprintf("\nBlock correlation validation:\n"))
    cat(sprintf("  Within-module mean cor: %.3f (n=%d pairs)\n", 
                mean(within_cors, na.rm=TRUE), length(within_cors)))
    cat(sprintf("  Between-module mean cor: %.3f (n=%d pairs)\n",
                mean(between_cors, na.rm=TRUE), length(between_cors)))
  }
  
  cat("\n[DEMO v2 DONE]\n")
}
