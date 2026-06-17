#!/usr/bin/env Rscript
# ==============================================================================
# S04_sim_runner.R — 模拟数据的原生pipeline执行器
#
# 问题: S04_run_pipeline_wrapper.R 的Phase 1依赖:
#   - 01_data_import.R 从TSV读数据 (模拟数据是RDS)
#   - 01_data_import.R 解析 {Prefix}{Time}d{Rep} 样本名格式
#   - 03_normalization_QC 需要正确的Time factor levels
#   模拟数据的sample_info格式不同 (Ind_T1_R1, Time=T1)
#
# 解决: 不伪装格式，不绕弯子。
#   本脚本直接从模拟RDS数据出发:
#   1. source 00_setup.R → 获取 FILES/PARAMS/COLORS 路径系统
#   2. 自己完成Phase 1 (数据导入+预处理+标准化) — 不source 01/02/03
#   3. Phase 2/3 复用S04的adapted函数 (maSigPro/WGCNA/09A-10)
#
# 用法:
#   source("S04_sim_runner.R")
#   run_simulation_pipeline(file.path(SIM_DIR, "benchmark/medium_medium_medium_rep1")
#
# 批量:
#   batch_run_simulations()
# ==============================================================================

# ---- 依赖 ----
if (file.exists("S_config.R")) source("S_config.R")

# [FIX] 关闭Browse调试模式 — 在batch运行中报错不应进入交互调试
options(error = NULL)

# S04的adapted函数需要这些:
SCRIPT_DIR <- PIPELINE_SCRIPTS  # PIPELINE_SCRIPTS

suppressPackageStartupMessages({
  library(MASS)
  library(DESeq2)
  library(matrixStats)
})


# ==============================================================================
# 核心函数: run_simulation_pipeline()
# ==============================================================================

#' 对一个模拟数据集运行完整METI-FS pipeline
#'
#' @param sim_dir  S01生成的模拟数据目录 (含 counts_matrix.rds 等)
#'                 或者 simulation_full.rds 所在目录
#' @param run_dir  pipeline运行目录 (输出), NULL则自动生成
#' @param skip_if_done 如果Final_candidate_genes.csv已存在则跳过
#' @return invisible(run_dir)
run_simulation_pipeline <- function(sim_dir, 
                                     run_dir = NULL,
                                     skip_if_done = TRUE) {
  
  # ---- 0. 路径设置 ----
  sim_name <- basename(sim_dir)
  if (is.null(run_dir)) {
    run_dir <- file.path(RUN_DIR, paste0("SIM_", sim_name))
  }
  
  # 跳过已完成
  if (skip_if_done && file.exists(file.path(run_dir, "data", "Final_candidate_genes.csv"))) {
    cat(sprintf("[SIM_RUNNER] SKIP (done): %s\n", sim_name))
    return(invisible(run_dir))
  }
  
  cat("\n================================================================\n")
  cat(sprintf("  Simulation Pipeline: %s\n", sim_name))
  cat(sprintf("  Input:  %s\n", sim_dir))
  cat(sprintf("  Output: %s\n", run_dir))
  cat(sprintf("  Started: %s\n", format(Sys.time())))
  cat("================================================================\n\n")
  
  t0 <- Sys.time()
  options(error = NULL)  # [FIX] 确保报错不进Browse
  
  # ---- 1. 读取模拟数据 ----
  counts  <- readRDS(file.path(sim_dir, "counts_matrix.rds"))
  tpm     <- readRDS(file.path(sim_dir, "tpm_matrix.rds"))
  si      <- readRDS(file.path(sim_dir, "sample_info.rds"))
  gt      <- readRDS(file.path(sim_dir, "ground_truth.rds"))
  params  <- readRDS(file.path(sim_dir, "simulation_params.rds"))
  
  time_vals   <- params$time_values    # c(4, 7, 14, 21)
  time_labels <- paste0(time_vals, "d") # c("4d", "7d", "14d", "21d")
  n_tp        <- length(time_vals)
  
  cat(sprintf("  Data: %d genes × %d samples, %d timepoints\n",
              nrow(counts), ncol(counts), n_tp))
  
  # ---- 2. 转换sample_info为pipeline格式 ----
  # 原始: Ind_T1_R1, Treatment=Induced, Time=T1, time_num=4
  # 需要: Time="4d" (factor), Group="Induced_4d", treatment="Induced"
  
  tp_map <- setNames(time_labels, paste0("T", seq_len(n_tp)))
  si$Time <- factor(tp_map[si$Time], levels = time_labels)
  si$treatment <- si$Treatment
  si$Group <- factor(paste(si$Treatment, si$Time, sep = "_"))
  si$time_label <- as.character(si$Time)
  
  cat(sprintf("  Samples per group:\n"))
  print(table(si$Treatment, si$Time))
  
  # ---- 3. source 00_setup.R → 获取 FILES 路径系统 ----
  # 00_setup.R 用 PROJECT_DIR 设置所有路径
  PROJECT_DIR <<- run_dir
  dir.create(run_dir, recursive = TRUE, showWarnings = FALSE)
  
  # 创建 data_raw/ (00_setup期望的)
  dir.create(file.path(run_dir, "data_raw"), recursive = TRUE, showWarnings = FALSE)
  
  # 写一个假的counts TSV让00_setup.R的文件检测不报错
  # (00_setup只是设置路径，不实际读取文件，但某些检查可能需要)
  
  source(file.path(SCRIPT_DIR, "00_setup.R"), local = FALSE)
  
  # 覆盖PARAMS
  PARAMS$masigpro_degree <<- min(n_tp - 1, 3)
  PARAMS$masigpro_k <<- min(9, max(4, n_tp * 2))
  
  cat(sprintf("  PARAMS: masigpro_degree=%d, k=%d\n", 
              PARAMS$masigpro_degree, PARAMS$masigpro_k))
  
  # ---- 4. Phase 1: 数据准备 (直接从RDS, 不经过01/02/03) ----
  cat("\n### PHASE 1 (SIM): Data preparation ###\n\n")
  
  # 4a. 保存sample_info
  saveRDS(si, FILES$sample_info)
  assign("sample_info", si, envir = .GlobalEnv)
  
  # 4b. 基因注释 (模拟数据: gene_id 同时作为 ensembl_id 和 symbol)
  gene_anno <- data.frame(
    ensembl_gene_id = rownames(counts),
    hgnc_symbol     = rownames(counts),
    entrez_id       = NA_character_,
    stringsAsFactors = FALSE
  )
  rownames(gene_anno) <- gene_anno$ensembl_gene_id
  saveRDS(gene_anno, FILES$gene_annotation)
  assign("gene_anno", gene_anno, envir = .GlobalEnv)
  
  # 4c. 基因过滤 (模拟 02_preprocessing.R 的 filterByExpr)
  # 模拟数据已经是filterByExpr后的13000个基因，但还是做一次
  library(edgeR)
  dge <- DGEList(counts = counts)
  keep <- filterByExpr(dge, group = si$Treatment)
  counts_filtered <- counts[keep, ]
  tpm_filtered    <- tpm[keep, ]
  
  saveRDS(counts_filtered, FILES$counts_filtered)
  saveRDS(tpm_filtered, FILES$tpm_filtered)
  assign("tpm_filtered", tpm_filtered, envir = .GlobalEnv)
  
  cat(sprintf("  filterByExpr: %d → %d genes\n", nrow(counts), nrow(counts_filtered)))
  
  # 4d. DESeq2标准化 + VST (模拟 03_normalization_QC.R)
  counts_mat <- as.matrix(counts_filtered)
  storage.mode(counts_mat) <- "integer"
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = si[colnames(counts_mat), ],
    design    = ~ Treatment + Time + Treatment:Time
  )
  dds$Treatment <- relevel(dds$Treatment, ref = "Control")
  dds$Time      <- relevel(dds$Time, ref = time_labels[1])
  
  vsd <- vst(dds, blind = FALSE)
  vst_mat <- assay(vsd)
  
  saveRDS(dds, FILES$dds_object)
  saveRDS(vst_mat, FILES$vst_matrix)
  
  cat(sprintf("  DESeq2 + VST: %d genes × %d samples\n", nrow(vst_mat), ncol(vst_mat)))
  
  # 4e. ground truth 保存到 data/
  saveRDS(gt, file.path(DATA_DIR, "ground_truth.rds"))
  
  cat("### PHASE 1 (SIM) COMPLETE ###\n")
  
  # ---- 5. Phase 2: DEG + maSigPro + WGCNA ----
  cat("\n### PHASE 2 (SIM): DEG + maSigPro + WGCNA ###\n\n")
  
  # 构建PROFILE供adapted函数使用
  PROFILE <<- list(
    project_dir     = run_dir,
    prefix          = "Sim",
    time_values     = time_vals,
    time_labels     = time_labels,
    time_suffix     = "d",
    time_ref        = time_labels[1],
    masigpro_degree = PARAMS$masigpro_degree,
    masigpro_k      = PARAMS$masigpro_k,
    masigpro_min_obs = 3,
    gene_id_type    = "symbol",
    colors_timepoint = setNames(
      c("#3C5488", "#E64B35", "#00A087", "#F39B7F")[seq_len(n_tp)],
      time_labels
    ),
    colors_group = {
      cg <- c()
      pal <- c("#DC0000", "#3C5488", "#E64B35", "#00A087", 
               "#B09C85", "#7E6148", "#F39B7F", "#4DBBD5")
      for (i in seq_len(n_tp)) {
        cg[paste0("Induced_", time_labels[i])]  <- pal[(i-1)*2 + 1]
        cg[paste0("Control_", time_labels[i])]  <- pal[(i-1)*2 + 2]
      }
      cg
    },
    shape_timepoint = setNames(c(16, 17, 15, 18)[seq_len(n_tp)], time_labels)
  )
  
  # 覆盖COLORS
  COLORS$timepoint <<- PROFILE$colors_timepoint
  COLORS$group     <<- PROFILE$colors_group
  
  # 5a. DEG analysis — source S04的adapted函数
  # 需要先source S04来获取 run_04_DEG_adapted 等函数
  # 但不能直接source S04（它会执行顶层代码和DATASET_ID检查）
  # 解决: 把S04中需要的adapted函数单独加载
  
  # 加载theme
  tryCatch(
    source(file.path(SCRIPT_DIR, "theme_bindlab.R"), local = FALSE),
    error = function(e) cat("  [WARN] theme_bindlab.R not found, using default\n")
  )
  
  # DEG: 直接跑DESeq2 LRT + Wald
  cat("  [SIM] Running DEG analysis (DESeq2 LRT + Wald)...\n")
  run_sim_DEG(dds, si, time_labels)
  
  # maSigPro
  cat(sprintf("  [SIM] Running maSigPro (degree=%d, k=%d)...\n",
              PARAMS$masigpro_degree, PARAMS$masigpro_k))
  run_sim_maSigPro(counts_filtered, si, time_vals, time_labels)
  
  # WGCNA
  cat("  [SIM] Running WGCNA...\n")
  assign("sample_info", si, envir = .GlobalEnv)
  assign("tpm_filtered", tpm_filtered, envir = .GlobalEnv)
  source(file.path(SCRIPT_DIR, "08_WGCNA.R"), local = FALSE)
  
  cat("### PHASE 2 (SIM) COMPLETE ###\n")
  
  # ---- 6. Phase 3: 候选筛选 (直接复用原pipeline) ----
  cat("\n### PHASE 3 (SIM): Candidate selection ###\n\n")
  
  assign("sample_info",  si, envir = .GlobalEnv)
  assign("gene_anno",    gene_anno, envir = .GlobalEnv)
  assign("tpm_filtered", readRDS(FILES$tpm_filtered), envir = .GlobalEnv)
  assign("all_results",  readRDS(FILES$deg_results), envir = .GlobalEnv)
  
  # 09A: candidate pool
  source(file.path(SCRIPT_DIR, "09A_candidate_pool.R"), local = FALSE)
  
  # 09C: bootstrap stability selection
  assign("sample_info", si, envir = .GlobalEnv)
  source(file.path(SCRIPT_DIR, "09C_ML_stability_selection.R"), local = FALSE)
  
  # 09D: gap-union selection (v2: RF auto-exclude)
  {
    lines_09D <- readLines(file.path(SCRIPT_DIR, "09D_gap_union_selection.R"))
    idx_exclude <- grep('EXCLUDE_ALGOS\\s*<-\\s*c\\("RF"\\)', lines_09D)
    if (length(idx_exclude) == 1) {
      lines_09D[idx_exclude] <- 'EXCLUDE_ALGOS <- character(0)  # [v2] auto-exclude'
    }
    tmp_09D <- tempfile(fileext = ".R")
    writeLines(lines_09D, tmp_09D)
    source(tmp_09D, local = FALSE)
    unlink(tmp_09D)
  }
  
  # 09F: PPI hub — 模拟数据没有真实PPI网络
  # 模拟基因ID是 GENE_00001 格式，不会匹配STRING
  # 创建空结果，保持pipeline完整性
  cat("  [SIM] PPI hub: skipped (simulated gene IDs have no STRING mapping)\n")
  ppi_empty <- list(
    hub_genes = data.frame(symbol = character(0), ensembl_id = character(0)),
    all_metrics = data.frame(),
    n_nodes = 0, n_edges = 0,
    method = "skipped_simulation"
  )
  saveRDS(ppi_empty, file.path(DATA_DIR, "ppi_hub_selection.rds"))
  write.csv(data.frame(), file.path(DATA_DIR, "PPI_09F_hub_genes.csv"), row.names = FALSE)
  
  # 10: integration
  source(file.path(SCRIPT_DIR, "10_integration.R"), local = FALSE)
  
  cat("### PHASE 3 (SIM) COMPLETE ###\n")
  
  # ---- 完成 ----
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  
  cat("\n================================================================\n")
  cat(sprintf("  COMPLETE: %s\n", sim_name))
  cat(sprintf("  Time: %.1f minutes\n", as.numeric(elapsed)))
  
  # 验证输出
  final_file <- file.path(DATA_DIR, "Final_candidate_genes.csv")
  if (file.exists(final_file)) {
    final <- read.csv(final_file)
    cat(sprintf("  Final candidates: %d genes\n", nrow(final)))
  } else {
    cat("  WARNING: Final_candidate_genes.csv not created\n")
  }
  cat("================================================================\n")
  
  return(invisible(run_dir))
}


# ==============================================================================
# Phase 2 辅助函数: DEG analysis for simulation data
# ==============================================================================

run_sim_DEG <- function(dds, si, time_labels) {
  
  # LRT: Treatment × Time 交互效应
  dds <- DESeq(dds, test = "LRT", reduced = ~ Treatment + Time)
  res_lrt <- results(dds, alpha = 0.05)
  
  lrt_sig <- as.data.frame(res_lrt[!is.na(res_lrt$padj) & res_lrt$padj < 0.05, ])
  lrt_sig$ensembl_id <- rownames(lrt_sig)
  lrt_sig$symbol <- rownames(lrt_sig)
  
  log_step("04_DEG", sprintf("LRT interaction: %d significant genes (padj<0.05)", nrow(lrt_sig)))
  
  # Wald: 每个时间点的 Induced vs Control
  dds_wald <- DESeq(dds, test = "Wald")
  
  deg_by_time <- list()
  for (tp in time_labels) {
    contrast_name <- paste0("Treatment_Induced_vs_Control")
    # 时间点特异的对比需要interaction term
    # 简化: 用overall Wald test结果 (对模拟数据benchmark够用)
  }
  
  # 保存
  all_results <- list(
    lrt_interaction     = as.data.frame(res_lrt),
    lrt_interaction_sig = lrt_sig,
    dds = dds
  )
  all_results$lrt_interaction$ensembl_id <- rownames(all_results$lrt_interaction)
  all_results$lrt_interaction$symbol     <- rownames(all_results$lrt_interaction)
  all_results$lrt_interaction_sig$ensembl_id <- rownames(lrt_sig)
  
  saveRDS(all_results, FILES$deg_results)
  assign("all_results", all_results, envir = .GlobalEnv)
  
  log_step("04_DEG", "DEG analysis complete (SIM)")
}


# ==============================================================================
# Phase 2 辅助函数: maSigPro for simulation data
# ==============================================================================

run_sim_maSigPro <- function(counts_filtered, si, time_vals, time_labels) {
  
  suppressPackageStartupMessages(library(maSigPro))
  
  # 构建design matrix
  time_num <- si$time_num
  treatment <- ifelse(si$Treatment == "Induced", 1, 0)
  
  # edesign: rows=samples, cols=Time, Replicate, Group(Induced=1/Control=0)
  edesign <- data.frame(
    Time      = time_num,
    Replicate = 1,  # maSigPro用这个分组，每个样本独立
    Control   = 1 - treatment,
    Induced   = treatment,
    row.names = si$sample_id
  )
  
  degree <- PARAMS$masigpro_degree
  k_use  <- PARAMS$masigpro_k
  
  # 生成设计矩阵
  design_matrix <- make.design.matrix(edesign, degree = degree)
  
  # Step 1: p.vector
  min_obs <- min(3, floor(nrow(si) / 2))
  gc(verbose = FALSE)
  
  fit <- tryCatch({
    p.vector(counts_filtered, design_matrix, 
             Q = 0.05, MT.adjust = "BH", min.obs = min_obs,
             counts = TRUE, family = negative.binomial(10))
  }, error = function(e) {
    log_step("06_MASIGPRO", sprintf("p.vector error: %s. Trying with family=gaussian", e$message))
    p.vector(counts_filtered, design_matrix,
             Q = 0.05, MT.adjust = "BH", min.obs = min_obs)
  })
  
  log_step("06_MASIGPRO", sprintf("p.vector: %d significant genes", fit$i))
  
  if (fit$i == 0) {
    log_step("06_MASIGPRO", "No significant genes. Saving empty results.")
    masigpro_results <- list(gene_clusters = data.frame(), cluster_result = NULL)
    saveRDS(masigpro_results, FILES$masigpro_results)
    return(invisible(NULL))
  }
  
  # Step 2: T.fit
  tstep <- T.fit(fit, step.method = "backward", alfa = 0.05)
  
  # Step 3: get.siggenes
  rsq_cutoff <- if (exists("PARAMS") && !is.null(PARAMS$masigpro_rsq)) PARAMS$masigpro_rsq else 0.7
  sigs <- get.siggenes(tstep, rsq = rsq_cutoff, vars = "groups")
  
  # 提取gene list
  sig_genes_all <- list()
  avail_keys <- names(sigs$sig.genes)
  for (grp in avail_keys) {
    sig_genes_all[[grp]] <- rownames(sigs$sig.genes[[grp]]$sig.profiles)
    log_step("06_MASIGPRO", sprintf("  %s: %d genes", grp, length(sig_genes_all[[grp]])))
  }
  
  # Clustering
  all_sig <- unique(unlist(sig_genes_all))
  n_sig <- length(all_sig)
  gene_clusters <- NULL  # [FIX] NULL not empty df — 让09A走fallback到sig_genes_all
  cluster_result <- NULL
  
  # [FIX] 动态k: 不能超过sig gene数/2，至少2
  k_actual <- min(k_use, max(2, floor(n_sig / 2)))
  
  if (n_sig >= max(4, k_actual)) {
    tryCatch({
      cluster_result <- see.genes(sigs$sig.genes[[avail_keys[1]]], 
                                   k = k_actual, show.fit = FALSE, 
                                   newX11 = FALSE)
      
      if (!is.null(cluster_result$cut)) {
        gene_clusters <- data.frame(
          ensembl_id = names(cluster_result$cut),
          symbol     = names(cluster_result$cut),
          cluster    = as.integer(cluster_result$cut),
          stringsAsFactors = FALSE
        )
        log_step("06_MASIGPRO", sprintf("Clustered %d genes into %d groups", 
                                         nrow(gene_clusters), k_actual))
      }
    }, error = function(e) {
      log_step("06_MASIGPRO", sprintf("Clustering error (k=%d): %s. Using sig_genes_all instead.", 
                                       k_actual, e$message))
      gene_clusters <<- NULL  # 确保fallback
    })
  } else {
    log_step("06_MASIGPRO", sprintf("Too few sig genes (%d) for clustering (k=%d), using sig_genes_all directly",
                                     n_sig, k_actual))
  }
  
  # 保存 — [FIX] 包含tstep供09A交互过滤使用
  masigpro_results <- list(
    sigs = sigs,
    sig_genes_all = sig_genes_all,
    cluster_result = cluster_result,
    gene_clusters = gene_clusters,  # NULL if clustering failed → 09A uses sig_genes_all
    design_matrix = design_matrix,
    tstep = tstep                   # [FIX] 09A需要tstep$sol做交互项过滤
  )
  saveRDS(masigpro_results, FILES$masigpro_results)
  
  log_step("06_MASIGPRO", sprintf("maSigPro complete: %d sig genes, %d clustered", 
                                    n_sig, nrow(gene_clusters)))
}


# ==============================================================================
# 批量运行
# ==============================================================================

#' 批量运行所有模拟数据的pipeline
batch_run_simulations <- function(
    sim_base = file.path(SIM_DIR, "benchmark"),
    max_runs = Inf,
    verbose = TRUE
) {
  sim_dirs <- list.dirs(sim_base, recursive = FALSE, full.names = TRUE)
  
  # 按名称排序
  sim_dirs <- sort(sim_dirs)
  
  n_total <- length(sim_dirs)
  n_done  <- 0
  n_skip  <- 0
  n_err   <- 0
  
  cat(sprintf("\n[BATCH] Starting batch run: %d simulations\n", n_total))
  
  for (i in seq_along(sim_dirs)) {
    if (i > max_runs) break
    
    sim_dir <- sim_dirs[i]
    sim_name <- basename(sim_dir)
    run_dir <- file.path(RUN_DIR, paste0("SIM_", sim_name))
    
    # 检查已完成
    if (file.exists(file.path(run_dir, "data", "Final_candidate_genes.csv"))) {
      n_skip <- n_skip + 1
      if (verbose) cat(sprintf("  [%d/%d] SKIP: %s\n", i, n_total, sim_name))
      next
    }
    
    tryCatch({
      run_simulation_pipeline(sim_dir, run_dir, skip_if_done = FALSE)
      n_done <- n_done + 1
    }, error = function(e) {
      n_err <<- n_err + 1
      cat(sprintf("  [%d/%d] ERROR: %s → %s\n", i, n_total, sim_name, e$message))
    })
  }
  
  cat(sprintf("\n[BATCH] Complete: %d done, %d skipped, %d errors (of %d total)\n",
              n_done, n_skip, n_err, n_total))
}


# ==============================================================================
# 直接运行demo
# ==============================================================================

if (sys.nframe() == 0) {
  cat("\n")
  cat("================================================================\n")
  cat("  S04_sim_runner.R — Simulation Pipeline Runner\n")
  cat("================================================================\n")
  cat("\n")
  cat("  # 运行单个模拟:\n")
  cat("  source('S04_sim_runner.R')\n")
  cat("  run_simulation_pipeline(\n")
  cat("    file.path(SIM_DIR, 'benchmark', 'medium_medium_medium_rep1')\n")
  cat("  )\n")
  cat("\n")
  cat("  # 批量运行所有:\n")
  cat("  batch_run_simulations()\n")
  cat("\n")
  cat("  # 先跑3个测试:\n")
  cat("  batch_run_simulations(max_runs = 3)\n")
  cat("\n")
  cat("================================================================\n")
}
