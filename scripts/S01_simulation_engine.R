#!/usr/bin/env Rscript
# ==============================================================================
# S01_simulation_engine.R — 时序RNA-seq模拟数据生成引擎 (v2 — BIB增强版)
#
# 目的：生成带已知ground truth的模拟时序RNA-seq counts数据，
#       用于评估METI-FS pipeline各组件的sensitivity/specificity/FDR。
#
# v2 变更 (相对v1):
#   [1.1] Dispersion真实化: DESeq2-style α(μ) = α₀ + α₁/μ 参数化
#         文献: Love et al. 2014 (Genome Biology); 
#               Soneson (2014, compcodeR, Bioinformatics);
#               Vieth et al. 2017 (powsimR, Bioinformatics)
#   [1.2] 共表达block结构: TRUE_TEMPORAL基因可分成共表达模块
#         文献: Stabl (Hédou et al. 2024, Nat Biotechnol) 测试correlated features;
#               Langfelder & Horvath 2008 (BMC Bioinformatics) WGCNA模块假设
#   [1.1+] 添加library size变异性 (biological + technical variation)
#
# 模拟策略（文献依据）：
#   - NB分布模拟counts: Love et al. 2014 (DESeq2), Nueda et al. 2014 (maSigPro)
#   - 时序表达模式: Spies et al. 2019 (Briefings in Bioinformatics) 的模拟框架
#   - Dispersion参数: DESeq2经典参数化 α(μ) = α₀ + α₁/μ
#     典型值来自 Soneson 2014 / powsimR (Vieth 2017):
#       bulk RNA-seq: α₀ ≈ 0.01-0.1, α₁ ≈ 1-10
#   - 信噪比控制: Stabl (Hédou et al. 2024, Nature Biotechnology) 的合成数据方案
#   - 共表达结构: Cholesky分解注入block correlation
#
# 基因分类（ground truth）：
#   - Class 1: TRUE_TEMPORAL — 处理×时间交互效应显著 (pipeline真阳性)
#   - Class 2: TRUE_MAIN    — 处理主效应显著但不随时间变化
#   - Class 3: TRUE_TIME    — 只有时间效应，两组变化相同
#   - Class 4: NULL          — 无差异基因
#
# 输出：
#   - counts_matrix.rds, tpm_matrix.rds, sample_info.rds
#   - ground_truth.rds (含 module_id 列)
#   - simulation_params.rds (含 dispersion_source, correlation_structure)
# ==============================================================================

# ---- 0. 加载配置 ----
if (file.exists("S_config.R")) {
  source("S_config.R")

  
}

# ---- 0a. 依赖 ----
suppressPackageStartupMessages({
  if (!requireNamespace("MASS", quietly = TRUE)) install.packages("MASS")
  library(MASS)
})

# ==============================================================================
# 辅助函数: DESeq2-style dispersion
# ==============================================================================

#' 生成DESeq2-style dispersion参数
#' 
#' DESeq2的mean-dispersion关系: α(μ) = α_intercept + α_slope / μ
#' 即高表达基因dispersion低（更精确），低表达基因dispersion高（更noisy）
#' 加上gene-to-gene随机变异 (log-normal noise)
#'
#' 文献:
#'   Love et al. 2014, Genome Biology 15:550 (DESeq2)
#'   Soneson 2014, Bioinformatics 30(18):2670 (compcodeR benchmark)
#'   Vieth et al. 2017, Bioinformatics 33(21):3486 (powsimR)
#'
#' @param base_means 向量，每个基因的基础表达量
#' @param alpha_intercept 渐近dispersion (高表达基因的下限), 典型值 0.01-0.1
#' @param alpha_slope     dispersion斜率 (低表达基因额外dispersion), 典型值 1-10
#' @param noise_sd        gene-to-gene log-normal噪声SD, 典型值 0.3-0.8
#' @return 向量，每个基因的dispersion值
generate_deseq2_dispersions <- function(base_means, 
                                         alpha_intercept = 0.05, 
                                         alpha_slope = 4.0,
                                         noise_sd = 0.5) {
  # 趋势线: α(μ) = α₀ + α₁/μ
  trend <- alpha_intercept + alpha_slope / base_means
  
  # 加入gene-to-gene变异 (DESeq2中dispersion在趋势线周围log-normal散布)
  # log(α_gene) = log(α_trend) + ε, ε ~ N(0, σ²)
  n <- length(base_means)
  log_dispersions <- log(trend) + rnorm(n, mean = 0, sd = noise_sd)
  dispersions <- exp(log_dispersions)
  
  # 合理范围约束: [0.001, 5.0]
  dispersions <- pmax(dispersions, 0.001)
  dispersions <- pmin(dispersions, 5.0)
  
  return(dispersions)
}


# ==============================================================================
# 辅助函数: 共表达block结构
# ==============================================================================

#' 在log-space为一组基因注入block correlation
#'
#' 对于一组n个基因 × m个样本的log2-expression矩阵，
#' 注入模块内相关结构。
#'
#' 方法: 
#'   1. 将基因分成 n_modules 个block
#'   2. 每个block共享一个latent factor
#'   3. gene_expr = sqrt(ρ) * latent + sqrt(1-ρ) * independent
#'   这产生模块内相关约 ρ，模块间相关约 0
#'
#' 文献:
#'   Langfelder & Horvath 2008, BMC Bioinformatics 9:559 (WGCNA)
#'   Hédou et al. 2024, Nat Biotechnol 42:1581 (Stabl, R≈0.5 correlated)
#'
#' @param lfc_matrix  基因×时间点 的log2FC矩阵 (已有信号)
#' @param gene_ids    基因ID向量
#' @param n_modules   模块数
#' @param rho_within  模块内相关系数 (0-1), 默认0.4
#' @param n_samples   样本数 (用于生成latent factors)
#' @return list(lfc_matrix_corr, module_assignments)
inject_block_correlation <- function(lfc_matrix, gene_ids, 
                                      n_modules = 4, 
                                      rho_within = 0.4,
                                      n_timepoints = ncol(lfc_matrix)) {
  
  n_genes <- length(gene_ids)
  
  if (n_genes < n_modules) {
    # 基因太少，不分模块
    return(list(
      lfc_matrix = lfc_matrix,
      module_ids = rep("M1", n_genes)
    ))
  }
  
  # 1. 将基因均匀分配到模块
  module_ids <- paste0("M", rep(1:n_modules, length.out = n_genes))
  
  # 2. 对每个模块，生成共享的时序噪声pattern
  #    这使得同模块基因在各时间点的行为更一致
  for (m in 1:n_modules) {
    idx <- which(module_ids == paste0("M", m))
    if (length(idx) < 2) next
    
    # 模块共享的latent temporal perturbation
    # 在每个时间点上，同模块基因有共同的微扰
    latent_perturbation <- rnorm(n_timepoints, mean = 0, sd = 0.3)
    
    for (g in idx) {
      # gene_lfc = sqrt(1-ρ) * original_lfc + sqrt(ρ) * module_shared
      # 保持信号方向，但加入模块一致性
      original <- lfc_matrix[g, ]
      shared <- latent_perturbation * sign(mean(original[original != 0]) + 0.01)
      
      # 混合: 保留原始信号主体，叠加模块共享分量
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
# [v3] 层级化共表达结构生成器
# ==============================================================================

#' 生成层级化共表达模块结构
#'
#' 模拟多层调控网络: 8个模块组织为2个super-module
#' Super-module内ρ=0.3, 子模块内ρ=0.5, 反映多层级生物调控
#'
#' 文献:
#'   Saelens et al. 2018, Nature Communications 9:4358 (模块检测方法基准)
#'   Langfelder & Horvath 2008, BMC Bioinformatics 9:559 (WGCNA)
inject_hierarchical_correlation <- function(lfc_matrix, gene_ids,
                                              n_submodules = 8,
                                              n_supermodules = 2,
                                              rho_super = 0.3,
                                              rho_sub = 0.5,
                                              n_timepoints = ncol(lfc_matrix)) {

  n_genes <- length(gene_ids)
  if (n_genes < n_submodules) {
    return(list(lfc_matrix = lfc_matrix,
                module_ids = rep("M1", n_genes),
                super_ids = rep("S1", n_genes)))
  }

  # 1. 分配基因到子模块 (均匀分布)
  sub_ids <- paste0("M", rep(1:n_submodules, length.out = n_genes))
  # 子模块 → super-module 映射
  sub_per_super <- n_submodules / n_supermodules
  super_ids <- paste0("S", ceiling(as.numeric(gsub("M", "", sub_ids)) / sub_per_super))

  # 2. 生成super-module级别的共享latent factor
  super_factors <- list()
  for (s in unique(super_ids)) {
    super_factors[[s]] <- rnorm(n_timepoints, mean = 0, sd = 0.3)
  }

  # 3. 生成子模块级别的共享latent factor
  sub_factors <- list()
  for (m in unique(sub_ids)) {
    sub_factors[[m]] <- rnorm(n_timepoints, mean = 0, sd = 0.2)
  }

  # 4. 混合: sqrt(rho_sub)*sub_factor + sqrt(rho_super)*super_factor + sqrt(1-rho_sub-rho_super)*original
  total_rho <- rho_sub + rho_super
  for (g in seq_len(n_genes)) {
    if (g > nrow(lfc_matrix)) next
    original <- lfc_matrix[g, ]
    s_id <- super_ids[g]
    m_id <- sub_ids[g]
    sub_signal <- sub_factors[[m_id]] * max(abs(original[original != 0]), 0.01)
    super_signal <- super_factors[[s_id]] * max(abs(original[original != 0]), 0.01) * 0.5
    lfc_matrix[g, ] <- sqrt(1 - total_rho) * original +
                        sqrt(rho_sub) * sub_signal +
                        sqrt(rho_super) * super_signal
  }

  return(list(
    lfc_matrix = lfc_matrix,
    module_ids = sub_ids,
    super_ids  = super_ids
  ))
}


#' 生成重叠模块共表达结构
#'
#' 模拟真实生物中基因的多功能参与: 20%基因同时属于2个模块
#' 使用混合latent factor, 反映combinatorial gene regulation
#'
#' 文献:
#'   Saelens et al. 2018, Nature Communications (overlapping module benchmark)
#'   PLOS Comp Biol 2024 (SBM analysis of transcriptome modularity)
inject_overlapping_correlation <- function(lfc_matrix, gene_ids,
                                             n_modules = 6,
                                             rho_within = 0.4,
                                             overlap_fraction = 0.2,
                                             n_timepoints = ncol(lfc_matrix)) {

  n_genes <- length(gene_ids)
  if (n_genes < n_modules) {
    return(list(lfc_matrix = lfc_matrix,
                module_ids = rep("M1", n_genes),
                is_overlap = rep(FALSE, n_genes)))
  }

  # 1. 分配基因: 80% 单模块, 20% 双模块
  n_overlap <- round(n_genes * overlap_fraction)
  n_single <- n_genes - n_overlap
  n_per_module_single <- ceiling(n_single / n_modules)

  # 2. 生成模块共享的latent factor
  module_factors <- list()
  for (m in 1:n_modules) {
    module_factors[[m]] <- rnorm(n_timepoints, mean = 0, sd = 0.4)
  }

  # 3. 分配单模块基因
  module_ids <- rep(NA_character_, n_genes)
  overlap_modules <- vector("list", n_genes)
  is_overlap <- rep(FALSE, n_genes)

  single_idx <- 1:n_single
  module_ids[single_idx] <- paste0("M", rep(1:n_modules, length.out = n_single))

  # 4. 分配双模块基因
  for (i in (n_single + 1):n_genes) {
    mods <- sample(1:n_modules, 2, replace = FALSE)
    module_ids[i] <- paste0("M", mods[1], "_M", mods[2])
    overlap_modules[[i]] <- mods
    is_overlap[i] <- TRUE
  }

  # 5. 注入相关结构
  for (g in seq_len(n_genes)) {
    if (g > nrow(lfc_matrix)) next
    original <- lfc_matrix[g, ]

    if (is_overlap[g]) {
      # 双模块基因: 混合两个模块的latent factor
      mods <- overlap_modules[[g]]
      latent1 <- module_factors[[mods[1]]]
      latent2 <- module_factors[[mods[2]]]
      mixed_latent <- (latent1 + latent2) / sqrt(2)
      lfc_matrix[g, ] <- sqrt(1 - rho_within) * original +
                          sqrt(rho_within) * mixed_latent * max(abs(original[original != 0]), 0.01)
    } else {
      # 单模块基因: 标准injection
      mod <- as.numeric(gsub("M", "", module_ids[g]))
      latent <- module_factors[[mod]]
      lfc_matrix[g, ] <- sqrt(1 - rho_within) * original +
                          sqrt(rho_within) * latent * max(abs(original[original != 0]), 0.01)
    }
  }

  return(list(
    lfc_matrix = lfc_matrix,
    module_ids = module_ids,
    is_overlap = is_overlap,
    overlap_modules = overlap_modules
  ))
}


#' 共表达结构生成分发器
#'
#' 根据structure_type参数调度到不同的结构生成器
generate_coexpression_structure <- function(lfc_matrix, gene_ids, structure_type,
                                              n_timepoints = ncol(lfc_matrix)) {
  switch(structure_type,
    "independent" = list(
      lfc_matrix = lfc_matrix,
      module_ids = rep("M0", length(gene_ids))
    ),
    "block4_rho0.4" = inject_block_correlation(lfc_matrix, gene_ids, 4, 0.4, n_timepoints),
    "block4_rho0.2" = inject_block_correlation(lfc_matrix, gene_ids, 4, 0.2, n_timepoints),
    "block4_rho0.6" = inject_block_correlation(lfc_matrix, gene_ids, 4, 0.6, n_timepoints),
    "hierarchical_8" = inject_hierarchical_correlation(lfc_matrix, gene_ids, 8, 2, 0.3, 0.5, n_timepoints),
    "overlapping"    = inject_overlapping_correlation(lfc_matrix, gene_ids, 6, 0.4, 0.2, n_timepoints),
    # 默认: 无共表达结构
    list(lfc_matrix = lfc_matrix,
         module_ids = rep("M0", length(gene_ids)))
  )
}


# ==============================================================================
# 核心函数: generate_simulation() — v3 增强版
# ==============================================================================

generate_simulation <- function(
    # ---- 实验设计参数 ----
    n_genes       = 13000,
    n_timepoints  = 4,
    time_values   = NULL,
    n_reps_ind    = 4,
    n_reps_ctrl   = 3,
    
    # ---- 信号参数 ----
    n_true_temporal  = 20,
    n_true_main      = 50,
    n_true_timeonly   = 100,
    
    # ---- 效应量参数 ----
    snr = "medium",
    lfc_range = NULL,
    
    # ---- 基因表达参数 ----
    base_mean_range = c(50, 5000),
    dispersion_range = c(0.1, 0.5),  # 仅在 dispersion_source="linear" 时使用
    
    # ---- [v2 NEW] Dispersion source ----
    # "deseq2": DESeq2-style α(μ) = α₀ + α₁/μ + noise (推荐)
    # "linear": 原v1行为 (线性映射, 向后兼容)
    dispersion_source = "deseq2",
    
    # DESeq2 dispersion参数 (仅 dispersion_source="deseq2" 时使用)
    # 典型bulk RNA-seq值: α₀=0.05, α₁=4.0, noise_sd=0.5
    # 参考: powsimR (Vieth 2017) 从ENCODE/GTEx估计的参数
    disp_intercept = 0.05,
    disp_slope     = 4.0,
    disp_noise_sd  = 0.5,
    
    # ---- [v3 NEW] Correlation structure ----
    # "independent": 基因独立 (原v1行为, 向后兼容)
    # "block4_rho0.4/0.2/0.6": TRUE_TEMPORAL基因分成共表达模块
    # "hierarchical_8": 层级化8模块 (2 super-module x 4 sub-module)
    # "overlapping": 20%基因属2+模块的重叠结构
    correlation_structure = "independent",
    n_modules     = 4,     # 模块数 (向后兼容, v3中自动根据structure_type设置)
    rho_within    = 0.4,   # 模块内相关 (向后兼容)

    # ---- 时序模式参数 ----
    temporal_patterns = c("sustained_up", "sustained_down", 
                          "early_peak", "late_onset",
                          "transient_up", "transient_down"),
    
    # ---- 输出控制 ----
    output_dir = NULL,
    seed = 42,
    verbose = TRUE
) {
  
  set.seed(seed)
  
  # ---- 1. 解析参数 ----
  if (is.null(time_values)) {
    time_values <- seq_len(n_timepoints)
  }
  stopifnot(length(time_values) == n_timepoints)
  
  # 信噪比 -> log2FC范围
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
    cat(sprintf("  [v3] Correlation: %s", correlation_structure))
    if (correlation_structure != "independent") {
      cat(sprintf(" (structure=%s)", correlation_structure))
    }
    cat("\n")
    cat(sprintf("  Seed: %d\n", seed))
  }
  
  # ---- 2. 构建样本信息表 ----
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
  
  # ---- 3. 生成基因参数 ----
  gene_ids <- sprintf("GENE_%05d", seq_len(n_genes))
  
  # 基础表达量（log-normal分布）
  base_means <- exp(runif(n_genes, 
                          log(base_mean_range[1]), 
                          log(base_mean_range[2])))
  
  # [v2] Dispersion参数 — 根据source选择生成方式
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
    # 原v1行为: 线性映射
    dispersions <- dispersion_range[2] - 
      (dispersion_range[2] - dispersion_range[1]) * 
      (log(base_means) - log(base_mean_range[1])) / 
      (log(base_mean_range[2]) - log(base_mean_range[1]))
    dispersions <- pmax(dispersions, dispersion_range[1])
  }
  
  # ---- 4. 分配基因类别 ----
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
  
  # ---- 5. 生成时序表达模式 ----
  
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
  
  # ---- 5b. [v3] 注入共表达结构 (通过分发器支持多种结构) ----
  module_ids <- rep(NA_character_, n_genes)
  super_ids  <- rep(NA_character_, n_genes)   # [v3] 层级结构上级模块
  
  if (correlation_structure != "independent") {
    temporal_idx <- which(gene_class == "TRUE_TEMPORAL")
    
    if (length(temporal_idx) >= 4) {
      struct_result <- generate_coexpression_structure(
        lfc_matrix      = lfc_temporal_matrix[temporal_idx, , drop = FALSE],
        gene_ids        = gene_ids[temporal_idx],
        structure_type  = correlation_structure,
        n_timepoints    = n_timepoints
      )
      lfc_temporal_matrix[temporal_idx, ] <- struct_result$lfc_matrix
      module_ids[temporal_idx] <- struct_result$module_ids
      if (!is.null(struct_result$super_ids)) {
        super_ids[temporal_idx] <- struct_result$super_ids
      }
      
      if (verbose) {
        n_mods <- length(unique(na.omit(module_ids[temporal_idx])))
        cat(sprintf("  [v3] Co-expression structure '%s' injected: %d temporal genes -> %d modules
",
                    correlation_structure, length(temporal_idx), n_mods))
      }
    }
  }
  
  # ---- 6. 生成counts矩阵 ----
  counts_matrix <- matrix(0L, nrow = n_genes, ncol = n_total_samples)
  rownames(counts_matrix) <- gene_ids
  colnames(counts_matrix) <- sample_info$sample_id
  
  time_effects_list <- attr(gene_params, "time_effects")
  
  # [v2] 样本特异的library size乘数 (模拟测序深度差异)
  # log-normal: 中位数1.0, 样本间CV约20%
  lib_size_factors <- exp(rnorm(n_total_samples, mean = 0, sd = 0.2))
  
  for (i in seq_len(n_genes)) {
    base_mu <- base_means[i]
    disp <- dispersions[i]
    cls <- gene_class[i]
    
    for (j in seq_len(n_total_samples)) {
      trt <- sample_info$Treatment[j]
      tp_idx <- match(sample_info$Time[j], paste0("T", seq_len(n_timepoints)))
      
      mu <- base_mu * lib_size_factors[j]  # [v2] 加入lib size variation
      
      # 时间效应 (Class 3)
      if (cls == "TRUE_TIME" && !is.null(time_effects_list[[gene_ids[i]]])) {
        time_eff <- time_effects_list[[gene_ids[i]]]
        mu <- mu * 2^(time_eff[tp_idx])
      }
      
      # 处理效应 (仅Induced组)
      if (trt == "Induced") {
        lfc <- lfc_temporal_matrix[i, tp_idx]
        mu <- mu * 2^(lfc)
      }
      
      mu <- max(mu, 1)
      
      # NB采样: mean = mu, size = 1/dispersion
      size_param <- 1 / disp
      counts_matrix[i, j] <- rnbinom(1, mu = mu, size = size_param)
    }
  }
  
  # ---- 7. 生成模拟TPM矩阵 ----
  gene_lengths <- rep(2000, n_genes)
  rpk_matrix <- counts_matrix / (gene_lengths / 1000)
  tpm_matrix <- apply(rpk_matrix, 2, function(x) x / sum(x) * 1e6)
  
  # ---- 8. 生成ground truth标签 ----
  ground_truth <- data.frame(
    gene_id = gene_ids,
    class = gene_class,
    is_true_marker = gene_class == "TRUE_TEMPORAL",
    is_any_signal = gene_class != "NULL",
    base_mean = base_means,
    dispersion = dispersions,
    lfc_max = gene_params$lfc_max,
    pattern = gene_params$pattern,
    module_id = module_ids,           # [v2] 模块归属
    stringsAsFactors = FALSE
  )
  
  # ---- 9. 汇总输出 ----
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
      correlation_structure = correlation_structure,   # [v3]
      n_modules = if (correlation_structure != "independent") length(unique(na.omit(module_ids))) else NA,
      rho_within = if (correlation_structure != "independent") rho_within else NA,
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
    
    # [v3] 模块信息
    if (correlation_structure != "independent" && any(!is.na(module_ids))) {
      mod_tab <- table(module_ids[!is.na(module_ids)])
      cat("  Modules: ")
      cat(paste(sprintf("%s=%d", names(mod_tab), mod_tab), collapse = ", "))
      cat("\n")
    }
    
    # [v2] Dispersion诊断
    if (dispersion_source == "deseq2") {
      cat(sprintf("  Dispersion (DESeq2-style): median=%.4f, IQR=[%.4f, %.4f]\n",
                  median(dispersions), quantile(dispersions, 0.25), quantile(dispersions, 0.75)))
    }
  }
  
  # 保存文件
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
# 预定义场景：覆盖文章中需要评估的参数空间
# ==============================================================================

#' 生成完整的模拟benchmark矩阵 (v2 增强版)
#' 
#' 核心参数空间: 3 SNR × 3 样本量 × 3 marker密度 = 27场景 × 5重复 = 135次
#' [v2] 附加相关场景: 3 SNR × block correlation × medium样本 × medium密度 = 3场景 × 5重复 = 15次
#' 总计: 150次模拟 (仅比原来多15次，计算开销可控)

generate_benchmark_scenarios <- function(
    output_base_dir = if (exists("SIM_DIR")) file.path(SIM_DIR, "benchmark") else "simulation_benchmark",
    n_repeats = 5,
    include_correlation = TRUE,   # [v2] 是否包含相关结构场景
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
  
  # ---- 核心场景: 27 × 5 = 135 ----
  scenarios <- expand.grid(
    snr = snr_levels,
    sample_size = names(sample_configs),
    marker_density = names(marker_configs),
    repeat_id = seq_len(n_repeats),
    stringsAsFactors = FALSE
  )
  scenarios$correlation <- "independent"
  
  # ---- [v3] 附加共表达结构场景: 3 SNR x 5 structures x medium x medium x 3 reps = 45 ----
  if (include_correlation) {
    corr_structures <- c("block4_rho0.4", "block4_rho0.2", "block4_rho0.6", 
                         "hierarchical_8", "overlapping")
    corr_scenarios <- expand.grid(
      snr = snr_levels,
      sample_size = "medium",
      marker_density = "medium",
      correlation = corr_structures,
      repeat_id = seq_len(n_repeats),
      stringsAsFactors = FALSE
    )
    scenarios <- rbind(scenarios, corr_scenarios)
  }
  

  # ---- Temporal depth extension: 3 SNR x 2 TPs x medium x medium x 3 reps = 18 ----
  if (exists("BENCHMARK_PARAMS") &&
      isTRUE(BENCHMARK_PARAMS$include_temporal_extension)) {
    tp_scenarios <- expand.grid(
      snr = snr_levels,
      sample_size = "medium",
      marker_density = "medium",
      n_timepoints = BENCHMARK_PARAMS$temporal_extension_tps,
      repeat_id = seq_len(BENCHMARK_PARAMS$temporal_extension_n_repeats),
      stringsAsFactors = FALSE
    )
    tp_scenarios$correlation <- "block4_rho0.4"
    scenarios <- rbind(scenarios, tp_scenarios)
    if (verbose) {
      cat(sprintf("  Temporal depth: %d extra scenarios (6 or 8 timepoints)
",
                  nrow(tp_scenarios)))
    }
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
    
    # 场景命名: 含correlation标记
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
      n_timepoints = if (!is.null(sc$n_timepoints)) sc$n_timepoints else 4,
      time_values = c(4, 7, 14, 21),
      n_reps_ind = sc_sample$n_reps_ind,
      n_reps_ctrl = sc_sample$n_reps_ctrl,
      n_true_temporal = sc_marker$n_true_temporal,
      n_true_main = sc_marker$n_true_main,
      n_true_timeonly = sc_marker$n_true_timeonly,
      snr = sc$snr,
      dispersion_source = "deseq2",                     # [v2] 始终用DESeq2-style
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
# 评估函数: 比较pipeline输出与ground truth (不变)
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
# [v3] 扩展评估: 支持overlapping ground truth (补丁B)
# ==============================================================================

#' 扩展评估函数: 支持overlapping ground truth
#'
#' 当基因属于多个模块时, 选中任一模块即算true positive (relaxed precision)
#' 与严格评估 (精确模块匹配) 互补
#'
#' @param selected_genes 选中的基因ID向量
#' @param ground_truth   ground truth data.frame (需含 module_id, is_overlap列)
#' @param target_class   目标基因类别
#' @return data.frame with strict + relaxed metrics
evaluate_selection_overlap <- function(selected_genes, ground_truth,
                                        target_class = "TRUE_TEMPORAL") {

  # 标准严格评估
  strict <- evaluate_selection(selected_genes, ground_truth, target_class)

  # 宽松评估: 基因属于至少1个signal模块
  true_positives <- ground_truth$gene_id[ground_truth$class == target_class]
  all_signal <- ground_truth$gene_id[ground_truth$class != "NULL"]

  tp_any <- sum(selected_genes %in% all_signal)
  precision_any <- ifelse(length(selected_genes) > 0,
                          tp_any / length(selected_genes), 0)

  # 模块级宽松: 选中的基因属于任意TRUE_TEMPORAL模块
  if ("module_id" %in% colnames(ground_truth) &&
      any(!is.na(ground_truth$module_id[ground_truth$is_true_marker]))) {

    true_modules <- unique(na.omit(ground_truth$module_id[ground_truth$is_true_marker]))
    selected_in_truth <- ground_truth[ground_truth$gene_id %in% selected_genes, ]

    # Module-aware: a gene selected from ANY signal module counts as TP
    tp_module_aware <- sum(selected_genes %in%
      ground_truth$gene_id[ground_truth$class == target_class |
                           ground_truth$module_id %in% true_modules])

    precision_module <- ifelse(length(selected_genes) > 0,
                               tp_module_aware / length(selected_genes), 0)
  } else {
    precision_module <- NA
  }

  result <- strict
  result$precision_any_signal <- round(precision_any, 4)
  result$precision_module_aware <- round(precision_module, 4)

  return(result)
}



# ==============================================================================
# Demo: 直接运行本脚本
# ==============================================================================

if (sys.nframe() == 0) {
  cat("\n========== Running demo simulation (v2) ==========\n\n")
  
  # Demo 1: DESeq2 dispersion + independent
  cat("--- Demo 1: DESeq2 dispersion, independent genes ---\n")
  sim1 <- generate_simulation(
    n_genes = 5000, n_timepoints = if (!is.null(sc$n_timepoints)) sc$n_timepoints else 4, time_values = c(4, 7, 14, 21),
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
  
  # 验证block结构确实产生了模块内相关
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
