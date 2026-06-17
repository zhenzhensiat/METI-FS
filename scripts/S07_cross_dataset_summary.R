#!/usr/bin/env Rscript
# ==============================================================================
# S07_cross_dataset_summary.R — v2 (修复+整合版)
#
# 变更记录 (v2 vs v1):
#   - FIX: Nogueira selection matrix转置 (genes×boots → boots×genes)
#   - NEW: PART 1b — batch_nogueira() 批量计算132个模拟场景
#   - NEW: PART 3b — build_ablation_paper_table() 合并A0-A8消融表
#   - FIX: run_id解析兼容 SIM_ 前缀和 block_corr 场景
#   - FIX: ablation config名称匹配 (A0_FULL vs A0_full)
#   - NEW: Fig.5 — Precision vs Recall + precision_any 高亮
#   - IMPROVED: generate_all_v2() 一键运行全部流程
#
# 运行方式:
#   See S_config.R for path configuration
#   source("S07_cross_dataset_summary.R")
#
#   # 一键全部 (推荐):
#   generate_all_v2()
#
#   # 或分步:
#   batch_nogueira()                    # ~10-20分钟
#   build_ablation_paper_table()        # <1分钟
#   plot_simulation_benchmark()         # Fig.2
#   plot_ablation_results_v2()          # Fig.3
#   plot_stability_nogueira()           # Fig.4
#   plot_precision_recall_highlight()   # Fig.5
#   generate_table1()                   # Table 1
#
# 文献依据:
#   - Nogueira, Sechidis & Brown (2018) JMLR 18(174):1-54: stability index
#   - stabm R包: Bommert et al. (2021) JOSS 6(59):3010
#   - pipeComp: Germain et al. (2020) Genome Biology 21:227
#   - Spooner et al. (2023) BMC Bioinformatics 24:9
#   - Stabl: Hédou et al. (2024) Nature Biotechnology 42:1581-1593
# ==============================================================================

# ---- 0. 依赖 ----
if (file.exists("S_config.R")) {
  source("S_config.R")

  
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(cowplot)
  library(RColorBrewer)
})

# 主题
tryCatch({
  source(file.path(PIPELINE_SCRIPTS, "theme_bindlab.R"))
}, error = function(e) {
  theme_bindlab <- function(base_size = 12) {
    theme_classic(base_size = base_size) %+replace%
      theme(
        plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
        axis.title = element_text(face = "bold"),
        legend.key = element_blank()
      )
  }
  assign("theme_bindlab", theme_bindlab, envir = .GlobalEnv)
})

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || all(is.na(x))) y else x

# 配色
METHODS_COLORS <- list(
  snr       = c("low" = "#3C5488", "medium" = "#00A087", "high" = "#E64B35"),
  sample    = c("small" = "#F39B7F", "medium" = "#4DBBD5", "large" = "#3C5488"),
  ablation  = c(
    "Full pipeline"         = "#00A087",
    "-maSigPro interaction" = "#E64B35",
    "-WGCNA"                = "#3C5488",
    "-Effect size"          = "#F39B7F",
    "-Gap-union"            = "#B09C85",
    "-PPI hub"              = "#7E6148",
    "-Bootstrap stability"  = "#DC9FB4",
    "Pool cap = 1000"       = "#91D1C2",
    "Pool cap = 3000"       = "#8491B4"
  ),
  algo      = c("LASSO" = "#E64B35", "RF" = "#4DBBD5", "SVM" = "#00A087")
)

# 图片保存
save_methods_fig <- function(plot_obj, filename, width = 8, height = 6, dpi = 300) {
  dir.create(FIG_DIR_METHODS, recursive = TRUE, showWarnings = FALSE)
  for (ext in c("pdf", "png")) {
    fpath <- file.path(FIG_DIR_METHODS, paste0(filename, ".", ext))
    if (ext == "pdf") {
      ggsave(fpath, plot = plot_obj, width = width, height = height,
             device = if (.Platform$OS.type == "windows") cairo_pdf else "pdf")
    } else {
      ggsave(fpath, plot = plot_obj, width = width, height = height,
             dpi = dpi, device = "png")
    }
  }
  methods_log("S07_FIG", sprintf("Saved: %s (.pdf + .png)", filename))
}


# ==============================================================================
# 辅助: run_id解析
# ==============================================================================

#' 解析 run_id 字符串为场景参数
#' 兼容格式:
#'   "high_large_dense_rep1"            → snr=high, sample=large, density=dense
#'   "SIM_high_large_dense_rep1"        → 同上 (去掉SIM_前缀)
#'   "medium_medium_medium_block_rep1"  → block_corr=TRUE
parse_run_id <- function(rid) {
  # 去掉SIM_前缀
  rid_clean <- sub("^SIM_", "", rid)
  parts <- strsplit(rid_clean, "_")[[1]]
  
  is_block <- "block" %in% parts
  # 去掉 "block" 和 "repN" 后取前3个
  core <- parts[!parts %in% c("block") & !grepl("^rep\\d+$", parts)]
  
  data.frame(
    snr            = if (length(core) >= 1) core[1] else NA_character_,
    sample_size    = if (length(core) >= 2) core[2] else NA_character_,
    marker_density = if (length(core) >= 3) core[3] else NA_character_,
    block_corr     = is_block,
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# PART 1: Nogueira稳定性指标 (修复版)
# ==============================================================================

# stabm包安装
install_if_missing_stabm <- function() {
  if (!requireNamespace("stabm", quietly = TRUE)) {
    install.packages("stabm", repos = "https://cloud.r-project.org")
  }
}

#' 计算Nogueira stability index
#'
#' @param sel_matrix M×p binary matrix (rows=bootstraps, cols=features)
#'        ★ 调用前必须已转置为正确方向!
#' @return numeric stability index
nogueira_stability <- function(sel_matrix) {
  if (is.null(sel_matrix) || nrow(sel_matrix) < 2) return(NA_real_)
  
  install_if_missing_stabm()
  
  M <- nrow(sel_matrix)
  p <- ncol(sel_matrix)
  
  k_bar <- mean(rowSums(sel_matrix))
  if (k_bar == 0 || k_bar == p) return(0)
  
  feature_lists <- lapply(seq_len(M), function(i) which(sel_matrix[i, ] == 1))
  non_empty <- feature_lists[vapply(feature_lists, length, integer(1)) > 0]
  if (length(non_empty) < 2) return(NA_real_)
  
  phi <- stabm::stabilityNogueira(features = non_empty, p = p)
  return(round(phi, 6))
}


#' 计算单个运行目录的Nogueira指标 (修复版: 添加转置)
#'
#' @param run_dir pipeline运行目录 (SIM_* 或 GEO_*)
#' @return data.frame with method, nogueira, k_bar, p, M
compute_nogueira_for_run <- function(run_dir) {
  
  data_dir <- file.path(run_dir, "data")
  stab_file <- file.path(data_dir, "ml_stability_selection.rds")
  if (!file.exists(stab_file)) return(NULL)
  
  stab <- tryCatch(readRDS(stab_file), error = function(e) NULL)
  if (is.null(stab)) return(NULL)
  
  results <- list()
  
  for (algo in c("lasso", "rf", "svm")) {
    mat_name <- paste0(algo, "_selection_matrix")
    if (!mat_name %in% names(stab)) next
    
    mat_raw <- stab[[mat_name]]
    if (is.null(mat_raw) || length(dim(mat_raw)) != 2) next
    
    # ════════════════════════════════════════════════════════════════
    # ★ 关键修复: 转置 ★
    # rds中存储: genes × bootstraps (行=基因名, 列=V1..V100)
    # 函数需要: bootstraps × genes (行=每轮, 列=每个基因)
    # ════════════════════════════════════════════════════════════════
    mat <- t(mat_raw)
    
    if (nrow(mat) < 2 || ncol(mat) < 1) next
    
    phi <- tryCatch(nogueira_stability(mat), error = function(e) NA_real_)
    k_bar <- mean(rowSums(mat))
    
    results[[algo]] <- data.frame(
      method   = toupper(algo),
      nogueira = phi,
      k_bar    = round(k_bar, 2),
      p        = ncol(mat),   # 基因数
      M        = nrow(mat),   # bootstrap次数
      stringsAsFactors = FALSE
    )
  }
  
  # Gap-union联合稳定性 (同样需要转置)
  gap_file <- file.path(data_dir, "ml_gap_union.rds")
  if (file.exists(gap_file)) {
    gap_data <- tryCatch(readRDS(gap_file), error = function(e) NULL)
    if (!is.null(gap_data) && !is.null(gap_data$final_gene_ids)) {
      final_ids <- gap_data$final_gene_ids
      
      # 从各算法转置后的矩阵中提取final genes列
      union_rows <- list()
      for (algo in c("lasso", "svm")) {
        mat_name <- paste0(algo, "_selection_matrix")
        if (!mat_name %in% names(stab)) next
        mat <- t(stab[[mat_name]])  # ★ 转置后: boots×genes, colnames=基因名
        common <- intersect(final_ids, colnames(mat))
        if (length(common) > 0) {
          union_rows[[algo]] <- mat[, common, drop = FALSE]
        }
      }
      
      if (length(union_rows) > 0) {
        M <- nrow(union_rows[[1]])
        all_final_genes <- unique(unlist(lapply(union_rows, colnames)))
        combined <- matrix(0L, nrow = M, ncol = length(all_final_genes))
        colnames(combined) <- all_final_genes
        for (algo_mat in union_rows) {
          for (g in colnames(algo_mat)) {
            combined[, g] <- pmax(combined[, g], algo_mat[, g])
          }
        }
        
        phi_union <- tryCatch(nogueira_stability(combined), error = function(e) NA_real_)
        results[["gap_union"]] <- data.frame(
          method   = "Gap-Union",
          nogueira = phi_union,
          k_bar    = round(mean(rowSums(combined)), 2),
          p        = ncol(combined),
          M        = M,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  
  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}


# ==============================================================================
# PART 1b: 批量Nogueira计算 (132个模拟场景)
# ==============================================================================

#' 批量计算所有SIM_*场景的Nogueira稳定性
#'
#' @return data.frame with run_id, method, nogueira, k_bar, p, M, snr, ...
#'         同时保存CSV到 BENCH_DIR/
batch_nogueira <- function() {
  
  cat("\n============================================================\n")
  cat("  PART 1b: Batch Nogueira Stability\n")
  cat("============================================================\n\n")
  
  # 发现SIM_*目录
  sim_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
  sim_dirs <- sim_dirs[grepl("^SIM_", basename(sim_dirs))]
  cat(sprintf("Found %d SIM_* directories\n", length(sim_dirs)))
  
  if (length(sim_dirs) == 0) {
    methods_log("S07_NOG", "No SIM_* directories found")
    return(invisible(NULL))
  }
  
  has_stab <- vapply(sim_dirs, function(d) {
    file.exists(file.path(d, "data", "ml_stability_selection.rds"))
  }, logical(1))
  sim_dirs <- sim_dirs[has_stab]
  cat(sprintf("  %d have ml_stability_selection.rds\n\n", length(sim_dirs)))
  
  all_results <- list()
  t0 <- Sys.time()
  n_done <- 0; n_err <- 0
  
  for (i in seq_along(sim_dirs)) {
    rd <- sim_dirs[i]
    run_id <- basename(rd)
    
    if (i %% 20 == 1) {
      elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
      cat(sprintf("[%d/%d] %s (%.1f min)\n", i, length(sim_dirs), run_id, elapsed))
    }
    
    res <- tryCatch(compute_nogueira_for_run(rd), error = function(e) {
      cat(sprintf("  [ERROR] %s: %s\n", run_id, e$message))
      n_err <<- n_err + 1
      NULL
    })
    
    if (!is.null(res)) {
      res$run_id <- sub("^SIM_", "", run_id)  # 与benchmark_master一致
      all_results[[run_id]] <- res
      n_done <- n_done + 1
    }
  }
  
  elapsed_total <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat(sprintf("\n[DONE] %d/%d computed, %d errors. %.1f min\n",
              n_done, length(sim_dirs), n_err, elapsed_total))
  
  if (length(all_results) == 0) return(invisible(NULL))
  
  nogueira_all <- do.call(rbind, all_results)
  rownames(nogueira_all) <- NULL
  
  # 解析场景参数
  parsed <- do.call(rbind, lapply(nogueira_all$run_id, parse_run_id))
  nogueira_all <- cbind(nogueira_all, parsed)
  
  # 保存完整结果
  out_all <- file.path(BENCH_DIR, "nogueira_stability_all.csv")
  write.csv(nogueira_all, out_all, row.names = FALSE)
  cat(sprintf("Saved: %s (%d rows)\n", out_all, nrow(nogueira_all)))
  
  # 按算法汇总
  summary_df <- nogueira_all %>%
    group_by(method) %>%
    summarise(
      n       = n(),
      mean    = round(mean(nogueira, na.rm = TRUE), 4),
      sd      = round(sd(nogueira, na.rm = TRUE), 4),
      median  = round(median(nogueira, na.rm = TRUE), 4),
      k_bar   = round(mean(k_bar, na.rm = TRUE), 1),
      .groups = "drop"
    )
  write.csv(summary_df, file.path(BENCH_DIR, "nogueira_summary.csv"), row.names = FALSE)
  
  # 按SNR×算法交叉
  summary_snr <- nogueira_all %>%
    filter(!is.na(snr)) %>%
    group_by(method, snr) %>%
    summarise(n = n(), mean = round(mean(nogueira, na.rm = TRUE), 4),
              sd = round(sd(nogueira, na.rm = TRUE), 4), .groups = "drop")
  write.csv(summary_snr, file.path(BENCH_DIR, "nogueira_by_snr.csv"), row.names = FALSE)
  
  cat("\n=== Nogueira Summary ===\n")
  print(as.data.frame(summary_df), row.names = FALSE)
  
  # 与Jaccard对比
  master_file <- file.path(BENCH_DIR, "benchmark_master.csv")
  if (file.exists(master_file)) {
    master <- read.csv(master_file, stringsAsFactors = FALSE)
    jaccard_cols <- grep("jaccard_mean$", names(master), value = TRUE)
    if (length(jaccard_cols) > 0) {
      cat("\n=== Jaccard vs Nogueira ===\n")
      for (jc in jaccard_cols) {
        algo <- toupper(gsub("_jaccard_mean", "", jc))
        jac_mean <- mean(master[[jc]], na.rm = TRUE)
        nog_sub <- nogueira_all %>% filter(method == algo)
        nog_mean <- mean(nog_sub$nogueira, na.rm = TRUE)
        cat(sprintf("  %s: Jaccard=%.4f  Nogueira=%.4f\n", algo, jac_mean, nog_mean))
      }
    }
  }
  
  methods_log("S07_NOG", sprintf("Batch Nogueira: %d runs, saved to %s", n_done, out_all))
  return(invisible(nogueira_all))
}


#' 诊断单个rds文件的矩阵方向 (debug用)
diagnose_stab_rds <- function(run_dir) {
  stab_file <- file.path(run_dir, "data", "ml_stability_selection.rds")
  if (!file.exists(stab_file)) { cat("File not found\n"); return(invisible(NULL)) }
  stab <- readRDS(stab_file)
  cat(sprintf("=== %s ===\n", basename(run_dir)))
  cat("Elements:", paste(names(stab), collapse=", "), "\n\n")
  for (algo in c("lasso", "rf", "svm")) {
    mat_name <- paste0(algo, "_selection_matrix")
    if (mat_name %in% names(stab)) {
      mat <- stab[[mat_name]]
      cat(sprintf("  %s: %d × %d (raw: genes×boots)\n", mat_name, nrow(mat), ncol(mat)))
      cat(sprintf("    → t(): %d boots × %d genes\n", ncol(mat), nrow(mat)))
      cat(sprintf("    rownames[1:3]: %s\n", paste(head(rownames(mat), 3), collapse=", ")))
      cat(sprintf("    colnames[1:3]: %s\n", paste(head(colnames(mat), 3), collapse=", ")))
    }
  }
}


# ==============================================================================
# PART 2: 模拟Benchmark可视化
# ==============================================================================

#' Fig.2: F1 by SNR × sample_size × density
plot_simulation_benchmark <- function(
    master_file = file.path(BENCH_DIR, "benchmark_master.csv")) {
  
  if (!file.exists(master_file)) {
    methods_log("S07_FIG", "benchmark_master.csv not found, skip Fig.2")
    return(invisible(NULL))
  }
  
  df <- read.csv(master_file, stringsAsFactors = FALSE)
  sim_df <- df[df$mode == "simulation" & !is.na(df$precision), ]
  if (nrow(sim_df) == 0) {
    methods_log("S07_FIG", "No simulation results"); return(invisible(NULL))
  }
  
  # 解析场景参数
  parsed <- do.call(rbind, lapply(sim_df$run_id, parse_run_id))
  sim_df <- cbind(sim_df, parsed)
  
  sim_df$snr            <- factor(sim_df$snr, levels = c("low", "medium", "high"))
  sim_df$sample_size    <- factor(sim_df$sample_size, levels = c("small", "medium", "large"))
  sim_df$marker_density <- factor(sim_df$marker_density, levels = c("sparse", "medium", "dense"))
  
  # 排除block场景 (单独分析)
  sim_ind <- sim_df[!sim_df$block_corr, ]
  
  # 汇总
  agg <- sim_ind %>%
    filter(!is.na(snr), !is.na(sample_size), !is.na(marker_density)) %>%
    group_by(snr, sample_size, marker_density) %>%
    summarise(
      F1_mean  = mean(F1, na.rm = TRUE),
      F1_sd    = sd(F1, na.rm = TRUE),
      Prec_mean = mean(precision, na.rm = TRUE),
      Rec_mean = mean(recall, na.rm = TRUE),
      n_runs   = n(),
      .groups  = "drop"
    )
  
  # Panel A: F1 bar
  p_f1 <- ggplot(agg, aes(x = snr, y = F1_mean, fill = sample_size)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, F1_mean - F1_sd),
                      ymax = pmin(1, F1_mean + F1_sd)),
                  position = position_dodge(0.8), width = 0.2, linewidth = 0.4) +
    facet_wrap(~ marker_density, labeller = labeller(
      marker_density = c("sparse" = "Sparse (10)", "medium" = "Medium (20)", "dense" = "Dense (50)")
    )) +
    scale_fill_manual(values = METHODS_COLORS$sample,
                      labels = c("small" = "n=16", "medium" = "n=28", "large" = "n=40")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "METI-FS Performance Across Simulation Scenarios",
         subtitle = sprintf("F1 score (mean ± SD, %d scenarios)", nrow(sim_ind)),
         x = "Signal-to-noise ratio", y = "F1 Score", fill = "Sample size") +
    theme_bindlab() +
    theme(legend.position = "bottom")
  
  save_methods_fig(p_f1, "Fig02A_F1_by_scenario", width = 10, height = 5)
  
  # Panel B: Precision vs Recall scatter
  p_pr <- ggplot(sim_ind, aes(x = recall, y = precision, color = snr, shape = sample_size)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_manual(values = METHODS_COLORS$snr) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Precision-Recall Trade-off",
         x = "Recall", y = "Precision",
         color = "SNR", shape = "Sample size") +
    theme_bindlab()
  
  save_methods_fig(p_pr, "Fig02B_precision_recall", width = 7, height = 6)
  
  # 保存Table 1
  tab1 <- agg %>%
    mutate(
      F1   = sprintf("%.3f ± %.3f", F1_mean, F1_sd),
      Prec = sprintf("%.3f", Prec_mean),
      Rec  = sprintf("%.3f", Rec_mean)
    ) %>%
    select(snr, sample_size, marker_density, n_runs, F1, Prec, Rec)
  
  dir.create(TAB_DIR_METHODS, recursive = TRUE, showWarnings = FALSE)
  write.csv(tab1, file.path(TAB_DIR_METHODS, "Table01_benchmark_summary.csv"),
            row.names = FALSE)
  methods_log("S07_TABLE", sprintf("Table 1: %d rows", nrow(tab1)))
  
  return(invisible(list(f1_plot = p_f1, pr_plot = p_pr, agg = agg)))
}


# ==============================================================================
# PART 3: 消融实验可视化 (v2: 兼容A1/A2 full + A3-A8 fast)
# ==============================================================================

#' Fig.3: 消融ΔF1 barplot
#'
#' 优先使用 ablation_paper_table.csv (build_ablation_paper_table产出),
#' 回退到 ablation_all_simulations.csv (A3-A8 only)
plot_ablation_results_v2 <- function() {
  
  # 优先用合并后的论文表
  paper_file <- file.path(BENCH_DIR, "ablation_paper_table.csv")
  if (file.exists(paper_file)) {
    agg <- read.csv(paper_file, stringsAsFactors = FALSE)
    if ("delta_F1" %in% names(agg) && "Label" %in% names(agg)) {
      ablation_data <- agg %>%
        filter(config != "A0_full" & config != "A0_FULL") %>%
        mutate(
          label     = Label,
          dF1_mean  = delta_F1,
          dF1_sd    = ifelse(is.na(F1_sd), 0, F1_sd)
        ) %>%
        arrange(dF1_mean) %>%
        mutate(label = factor(label, levels = label))
      
      return(.plot_ablation_bars(ablation_data))
    }
  }
  
  # 回退: ablation_all_simulations.csv
  abl_file <- file.path(BENCH_DIR, "ablation_all_simulations.csv")
  if (!file.exists(abl_file)) {
    methods_log("S07_FIG", "No ablation data found, skip Fig.3")
    return(invisible(NULL))
  }
  
  df <- read.csv(abl_file, stringsAsFactors = FALSE)
  if (!"F1" %in% names(df)) return(invisible(NULL))
  
  # config名称可能是 A0_FULL 或 A0_full — 统一为大写
  baseline_name <- if ("A0_FULL" %in% df$config) "A0_FULL" else "A0_full"
  
  df_delta <- df %>%
    group_by(run_id) %>%
    mutate(
      F1_base  = F1[config == baseline_name][1],
      delta_F1 = F1 - F1_base
    ) %>%
    ungroup() %>%
    filter(config != baseline_name)
  
  agg <- df_delta %>%
    group_by(config, label) %>%
    summarise(
      dF1_mean = mean(delta_F1, na.rm = TRUE),
      dF1_sd   = sd(delta_F1, na.rm = TRUE),
      n = n(), .groups = "drop"
    ) %>%
    arrange(dF1_mean) %>%
    mutate(label = factor(label, levels = label))
  
  .plot_ablation_bars(agg)
}

# 内部辅助: 画消融柱状图
.plot_ablation_bars <- function(agg) {
  p_abl <- ggplot(agg, aes(x = label, y = dF1_mean, fill = label)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = dF1_mean - dF1_sd, ymax = dF1_mean + dF1_sd),
                  width = 0.2, linewidth = 0.4) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = METHODS_COLORS$ablation, guide = "none") +
    coord_flip() +
    labs(title = "Component Ablation: Impact on F1 Score",
         subtitle = "ΔF1 relative to full pipeline (mean ± SD)",
         x = NULL, y = "ΔF1 Score") +
    theme_bindlab()
  
  save_methods_fig(p_abl, "Fig03_ablation_barplot", width = 8, height = 5)
  
  # Fig.3B: SNR热图 (如果有ablation_by_snr.csv)
  snr_file <- file.path(BENCH_DIR, "ablation_by_snr.csv")
  if (file.exists(snr_file)) {
    by_snr <- read.csv(snr_file, stringsAsFactors = FALSE)
    baseline_name2 <- if ("A0_FULL" %in% by_snr$config) "A0_FULL" else "A0_full"
    
    # 计算ΔF1 per SNR
    base_snr <- by_snr %>% filter(config == baseline_name2) %>% select(snr, F1_base = F1)
    heat_data <- by_snr %>%
      filter(config != baseline_name2) %>%
      left_join(base_snr, by = "snr") %>%
      mutate(dF1 = F1 - F1_base)
    
    # 需要label列
    if (!"label" %in% names(heat_data)) {
      heat_data$label <- heat_data$config
    }
    
    if (nrow(heat_data) > 0) {
      p_heat <- ggplot(heat_data, aes(x = snr, y = label, fill = dF1)) +
        geom_tile(color = "white", linewidth = 0.5) +
        geom_text(aes(label = sprintf("%+.2f", dF1)), size = 3.5) +
        scale_fill_gradient2(low = "#E64B35", mid = "white", high = "#00A087",
                             midpoint = 0, name = "ΔF1") +
        labs(title = "Ablation Impact by SNR", x = "SNR", y = NULL) +
        theme_bindlab() +
        theme(panel.grid = element_blank())
      
      save_methods_fig(p_heat, "Fig03B_ablation_by_snr", width = 7, height = 5)
    }
  }
  
  return(invisible(agg))
}


# ==============================================================================
# PART 3b: 合并消融结果为论文Table (A0-A8)
# ==============================================================================

#' 合并 ablation_A1_A2_full.csv 和 ablation_all_simulations.csv
#' 输出 ablation_paper_table.csv (9行)
build_ablation_paper_table <- function() {
  
  cat("\n============================================================\n")
  cat("  PART 3b: Ablation Paper Table\n")
  cat("============================================================\n\n")
  
  f_a12 <- file.path(BENCH_DIR, "ablation_A1_A2_full.csv")
  f_all <- file.path(BENCH_DIR, "ablation_all_simulations.csv")
  
  if (!file.exists(f_a12) && !file.exists(f_all)) {
    methods_log("S07_ABL", "No ablation CSV found")
    return(invisible(NULL))
  }
  
  # A1/A2完整消融
  if (file.exists(f_a12)) {
    ab_a12 <- read.csv(f_a12, stringsAsFactors = FALSE)
    cat(sprintf("  A1/A2 full: %d rows, configs: %s\n",
                nrow(ab_a12), paste(unique(ab_a12$config), collapse=", ")))
  } else {
    ab_a12 <- NULL
  }
  
  # A0-A8快速消融
  if (file.exists(f_all)) {
    ab_all <- read.csv(f_all, stringsAsFactors = FALSE)
    cat(sprintf("  All sims:   %d rows, configs: %s\n",
                nrow(ab_all), paste(unique(ab_all$config), collapse=", ")))
  } else {
    ab_all <- NULL
  }
  
  # 标签映射
  config_labels <- c(
    "A0_FULL"            = "Full pipeline (A0)",
    "A1_no_maSigPro"     = "-maSigPro interaction (A1)",
    "A2_no_WGCNA"        = "-WGCNA (A2)",
    "A3_no_EffectSize"   = "-Effect size filter (A3)",
    "A4_no_GapUnion"     = "-Gap-union threshold (A4)",
    "A5_no_PPI"          = "-PPI hub (A5)",
    "A6_no_Bootstrap"    = "-Bootstrap stability (A6)",
    "A7_pool_1000"       = "Pool cap=1000 (A7)",
    "A8_pool_3000"       = "Pool cap=3000 (A8)"
  )
  
  # 合并策略: A0/A1/A2 from ab_a12 (完整重跑, 更准确); A3-A8 from ab_all
  a012_configs <- c("A0_FULL", "A1_no_maSigPro", "A2_no_WGCNA")
  a38_configs  <- setdiff(names(config_labels), a012_configs)
  
  parts <- list()
  if (!is.null(ab_a12)) {
    parts$a012 <- ab_a12 %>% filter(config %in% a012_configs)
  }
  if (!is.null(ab_all)) {
    # 如果ab_a12不存在, 也从ab_all取A0/A1/A2
    if (is.null(ab_a12)) {
      parts$a012 <- ab_all %>% filter(config %in% a012_configs)
    }
    parts$a38 <- ab_all %>% filter(config %in% a38_configs)
  }
  
  combined <- bind_rows(parts)
  if (nrow(combined) == 0) {
    methods_log("S07_ABL", "No data after merge"); return(invisible(NULL))
  }
  
  # 确保n_final列存在
  if (!"n_final" %in% names(combined)) combined$n_final <- NA_integer_
  
  # 汇总
  summary_all <- combined %>%
    group_by(config) %>%
    summarise(
      n_scenarios    = n(),
      Precision_mean = round(mean(precision, na.rm = TRUE), 3),
      Precision_sd   = round(sd(precision, na.rm = TRUE), 3),
      Recall_mean    = round(mean(recall, na.rm = TRUE), 3),
      Recall_sd      = round(sd(recall, na.rm = TRUE), 3),
      F1_mean        = round(mean(F1, na.rm = TRUE), 3),
      F1_sd          = round(sd(F1, na.rm = TRUE), 3),
      FDR_mean       = round(mean(FDR, na.rm = TRUE), 3),
      n_final_mean   = round(mean(n_final, na.rm = TRUE), 1),
      .groups = "drop"
    )
  
  # ΔF1
  f1_a0 <- summary_all$F1_mean[summary_all$config == "A0_FULL"]
  if (length(f1_a0) == 0) f1_a0 <- NA
  summary_all$delta_F1 <- round(summary_all$F1_mean - f1_a0, 3)
  
  # 标签
  summary_all$Label <- config_labels[summary_all$config]
  summary_all$Label[is.na(summary_all$Label)] <- summary_all$config[is.na(summary_all$Label)]
  
  # 排序
  summary_all <- summary_all %>%
    mutate(config = factor(config, levels = names(config_labels))) %>%
    arrange(config) %>%
    mutate(config = as.character(config))
  
  paper_table <- summary_all %>%
    select(Label, config, n_scenarios,
           Precision_mean, Precision_sd,
           Recall_mean, Recall_sd,
           F1_mean, F1_sd, delta_F1,
           FDR_mean, n_final_mean)
  
  out_file <- file.path(BENCH_DIR, "ablation_paper_table.csv")
  write.csv(paper_table, out_file, row.names = FALSE)
  cat(sprintf("\nSaved: %s (%d rows)\n", out_file, nrow(paper_table)))
  
  # 按SNR分组
  parsed_snr <- do.call(rbind, lapply(combined$run_id, parse_run_id))
  combined$snr_parsed <- parsed_snr$snr
  
  by_snr <- combined %>%
    filter(!is.na(snr_parsed)) %>%
    group_by(config, snr = snr_parsed) %>%
    summarise(
      n  = n(),
      F1 = round(mean(F1, na.rm = TRUE), 3),
      F1_sd = round(sd(F1, na.rm = TRUE), 3),
      Prec = round(mean(precision, na.rm = TRUE), 3),
      .groups = "drop"
    )
  
  snr_file <- file.path(BENCH_DIR, "ablation_by_snr.csv")
  write.csv(by_snr, snr_file, row.names = FALSE)
  cat(sprintf("Saved: %s\n", snr_file))
  
  # 打印
  cat("\n=== Ablation Paper Table ===\n")
  print(as.data.frame(paper_table %>% select(Label, delta_F1, F1_mean, Precision_mean, Recall_mean)),
        row.names = FALSE)
  
  methods_log("S07_ABL", sprintf("Ablation table: %d configs", nrow(paper_table)))
  return(invisible(paper_table))
}


# ==============================================================================
# PART 4: 稳定性可视化 (Nogueira from batch results)
# ==============================================================================

#' Fig.4: Nogueira stability + Jaccard 对比
plot_stability_nogueira <- function() {
  
  nog_file <- file.path(BENCH_DIR, "nogueira_stability_all.csv")
  if (!file.exists(nog_file)) {
    methods_log("S07_FIG", "nogueira_stability_all.csv not found, skip Fig.4")
    return(invisible(NULL))
  }
  
  nog_df <- read.csv(nog_file, stringsAsFactors = FALSE)
  
  # 只取 LASSO/RF/SVM (排除Gap-Union, 单独处理)
  algo_df <- nog_df %>% filter(method %in% c("LASSO", "RF", "SVM"))
  
  # Panel A: Nogueira分布 (boxplot by algorithm)
  p_nog <- ggplot(algo_df, aes(x = method, y = nogueira, fill = method)) +
    geom_boxplot(width = 0.6, outlier.size = 1, outlier.alpha = 0.5) +
    scale_fill_manual(values = METHODS_COLORS$algo, guide = "none") +
    scale_y_continuous(limits = c(-0.1, 1), breaks = seq(0, 1, 0.2)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = "Feature Selection Stability (Nogueira Index)",
         subtitle = sprintf("Across %d simulation scenarios",
                            length(unique(algo_df$run_id))),
         x = "Algorithm", y = "Nogueira Stability Index") +
    theme_bindlab()
  
  # Panel B: Jaccard (from benchmark_master)
  master_file <- file.path(BENCH_DIR, "benchmark_master.csv")
  p_combined <- p_nog  # default
  
  if (file.exists(master_file)) {
    master <- read.csv(master_file, stringsAsFactors = FALSE)
    jac_long <- master %>%
      select(run_id,
             LASSO = lasso_jaccard_mean,
             RF    = rf_jaccard_mean,
             SVM   = svm_jaccard_mean) %>%
      pivot_longer(-run_id, names_to = "method", values_to = "jaccard") %>%
      filter(!is.na(jaccard))
    
    if (nrow(jac_long) > 0) {
      p_jac <- ggplot(jac_long, aes(x = method, y = jaccard, fill = method)) +
        geom_boxplot(width = 0.6, outlier.size = 1, outlier.alpha = 0.5) +
        scale_fill_manual(values = METHODS_COLORS$algo, guide = "none") +
        scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
        labs(title = "Pairwise Jaccard Index",
             subtitle = "Between bootstrap iterations",
             x = "Algorithm", y = "Jaccard Index") +
        theme_bindlab()
      
      p_combined <- plot_grid(p_nog, p_jac, nrow = 1, labels = c("A", "B"))
    }
  }
  
  save_methods_fig(p_combined, "Fig04_stability_comparison", width = 10, height = 5)
  
  # Nogueira by SNR (if data available)
  if ("snr" %in% names(algo_df) && sum(!is.na(algo_df$snr)) > 0) {
    algo_df$snr <- factor(algo_df$snr, levels = c("low", "medium", "high"))
    p_snr <- ggplot(algo_df %>% filter(!is.na(snr)),
                    aes(x = snr, y = nogueira, fill = method)) +
      geom_boxplot(width = 0.7, position = position_dodge(0.8),
                   outlier.size = 1, outlier.alpha = 0.5) +
      scale_fill_manual(values = METHODS_COLORS$algo) +
      scale_y_continuous(limits = c(-0.1, 1)) +
      labs(title = "Stability by SNR Level",
           x = "SNR", y = "Nogueira Index", fill = "Algorithm") +
      theme_bindlab() +
      theme(legend.position = "bottom")
    
    save_methods_fig(p_snr, "Fig04B_stability_by_snr", width = 8, height = 5)
  }
  
  methods_log("S07_FIG", "Fig.4 stability saved")
  return(invisible(algo_df))
}


# ==============================================================================
# PART 5: precision_any高亮 + 典型条件性能
# ==============================================================================

#' Fig.5: Precision_any高亮 + 典型条件汇总
plot_precision_recall_highlight <- function() {
  
  master_file <- file.path(BENCH_DIR, "benchmark_master.csv")
  if (!file.exists(master_file)) {
    methods_log("S07_FIG", "No master file, skip Fig.5"); return(invisible(NULL))
  }
  
  df <- read.csv(master_file, stringsAsFactors = FALSE)
  sim_df <- df[df$mode == "simulation" & !is.na(df$precision), ]
  if (nrow(sim_df) == 0) return(invisible(NULL))
  
  parsed <- do.call(rbind, lapply(sim_df$run_id, parse_run_id))
  sim_df <- cbind(sim_df, parsed)
  
  # Panel A: Precision (strict) vs Precision_any scatter
  if ("precision_any" %in% names(sim_df)) {
    sim_df$snr <- factor(sim_df$snr, levels = c("low", "medium", "high"))
    
    p_pany <- ggplot(sim_df, aes(x = precision, y = precision_any, color = snr)) +
      geom_point(size = 2, alpha = 0.7) +
      scale_color_manual(values = METHODS_COLORS$snr) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
      scale_x_continuous(limits = c(0, 1)) +
      scale_y_continuous(limits = c(0, 1.05)) +
      annotate("text", x = 0.1, y = 0.95,
               label = sprintf("Mean precision_any = %.3f", mean(sim_df$precision_any, na.rm = TRUE)),
               hjust = 0, size = 4, fontface = "bold", color = "#00A087") +
      labs(title = "Strict vs Relaxed Precision",
           subtitle = "precision_any: selected gene is ANY true signal (not just TRUE_TEMPORAL)",
           x = "Precision (strict: TRUE_TEMPORAL only)",
           y = "Precision_any (any TRUE_* category)",
           color = "SNR") +
      theme_bindlab()
    
    save_methods_fig(p_pany, "Fig05A_precision_any", width = 7, height = 6)
  }
  
  # Panel B: 典型条件 (medium SNR × medium/large sample) 性能分布
  typical <- sim_df %>%
    filter(snr == "medium",
           sample_size %in% c("medium", "large"),
           !block_corr)
  
  if (nrow(typical) > 0) {
    typical_long <- typical %>%
      select(run_id, Precision = precision, Recall = recall, F1) %>%
      pivot_longer(-run_id, names_to = "Metric", values_to = "Value")
    
    p_typical <- ggplot(typical_long, aes(x = Metric, y = Value, fill = Metric)) +
      geom_boxplot(width = 0.5) +
      geom_jitter(width = 0.15, size = 1.5, alpha = 0.5) +
      scale_fill_brewer(palette = "Set2", guide = "none") +
      scale_y_continuous(limits = c(0, 1)) +
      annotate("text", x = 3.3, y = 0.95,
               label = sprintf("n = %d scenarios", nrow(typical)),
               hjust = 1, size = 3.5, color = "grey40") +
      labs(title = "Performance Under Typical Experimental Conditions",
           subtitle = "Medium SNR, medium/large sample size (closest to real experiments)",
           x = NULL, y = "Score") +
      theme_bindlab()
    
    save_methods_fig(p_typical, "Fig05B_typical_conditions", width = 6, height = 5)
  }
  
  methods_log("S07_FIG", "Fig.5 saved")
  return(invisible(sim_df))
}


# ==============================================================================
# PART 6: 完整Table生成
# ==============================================================================

#' Table 1: 在 plot_simulation_benchmark() 中已生成
#' Table 2: ablation_paper_table.csv 已在 build_ablation_paper_table() 中生成

#' Table 3: 跨数据集一致性 (GEO数据, 有数据时才运行)
generate_table3 <- function(run_dirs = NULL) {
  
  if (is.null(run_dirs)) {
    run_dirs <- list.dirs(RUN_DIR, recursive = FALSE, full.names = TRUE)
    run_dirs <- run_dirs[grepl("^GEO_", basename(run_dirs))]
  }
  
  if (length(run_dirs) == 0) {
    methods_log("S07_TABLE", "No GEO run dirs, skip Table 3")
    return(invisible(NULL))
  }
  
  rows <- list()
  for (rd in run_dirs) {
    run_id <- basename(rd)
    data_dir <- file.path(rd, "data")
    row <- data.frame(dataset = run_id, stringsAsFactors = FALSE)
    
    pool_f <- file.path(data_dir, "candidate_pool.rds")
    if (file.exists(pool_f)) {
      pool <- readRDS(pool_f)
      row$n_pool <- length(pool$candidate_pool)
    }
    
    stab_f <- file.path(data_dir, "ml_stability_selection.rds")
    if (file.exists(stab_f)) {
      stab <- readRDS(stab_f)
      row$lasso_nonzero <- sum(stab$lasso_freq > 0)
      row$rf_nonzero    <- sum(stab$rf_freq > 0)
      row$svm_nonzero   <- sum(stab$svm_freq > 0)
    }
    
    ng <- compute_nogueira_for_run(rd)
    if (!is.null(ng)) {
      for (j in seq_len(nrow(ng))) {
        row[[paste0("nogueira_", tolower(ng$method[j]))]] <- ng$nogueira[j]
      }
    }
    
    fin_f <- file.path(data_dir, "Final_candidate_genes.csv")
    if (file.exists(fin_f)) row$n_final <- nrow(read.csv(fin_f))
    
    rows[[run_id]] <- row
  }
  
  if (length(rows) == 0) return(invisible(NULL))
  
  tab3 <- do.call(rbind, rows)
  write.csv(tab3, file.path(TAB_DIR_METHODS, "Table03_cross_dataset_consistency.csv"),
            row.names = FALSE)
  methods_log("S07_TABLE", sprintf("Table 3: %d datasets", nrow(tab3)))
  return(tab3)
}


# ==============================================================================
# PART 7: 一键生成 (v2)
# ==============================================================================

#' 完整运行: Nogueira计算 → 消融合并 → 出图 → 出表
#'
#' @param skip_nogueira 如果nogueira_stability_all.csv已存在, 跳过重算
generate_all_v2 <- function(skip_nogueira = FALSE) {
  
  t0 <- Sys.time()
  cat("\n================================================================\n")
  cat("  S07 v2: Generate All Figures and Tables\n")
  cat(sprintf("  Started: %s\n", format(t0)))
  cat("================================================================\n\n")
  
  dir.create(FIG_DIR_METHODS, recursive = TRUE, showWarnings = FALSE)
  dir.create(TAB_DIR_METHODS, recursive = TRUE, showWarnings = FALSE)
  
  # ---- Step 1: Nogueira batch ----
  nog_file <- file.path(BENCH_DIR, "nogueira_stability_all.csv")
  if (skip_nogueira && file.exists(nog_file)) {
    methods_log("S07_MAIN", "Nogueira: SKIPPED (file exists, skip_nogueira=TRUE)")
  } else {
    methods_log("S07_MAIN", "Step 1: Nogueira batch computation...")
    batch_nogueira()
  }
  
  # ---- Step 2: Ablation paper table ----
  methods_log("S07_MAIN", "Step 2: Ablation paper table...")
  build_ablation_paper_table()
  
  # ---- Step 3: Figures ----
  methods_log("S07_MAIN", "Step 3: Figures...")
  
  methods_log("S07_MAIN", "  Fig.2: Simulation benchmark...")
  tryCatch(plot_simulation_benchmark(), error = function(e)
    methods_log("S07_ERR", paste("Fig.2:", e$message)))
  
  methods_log("S07_MAIN", "  Fig.3: Ablation results...")
  tryCatch(plot_ablation_results_v2(), error = function(e)
    methods_log("S07_ERR", paste("Fig.3:", e$message)))
  
  methods_log("S07_MAIN", "  Fig.4: Stability comparison...")
  tryCatch(plot_stability_nogueira(), error = function(e)
    methods_log("S07_ERR", paste("Fig.4:", e$message)))
  
  methods_log("S07_MAIN", "  Fig.5: Precision_any + typical conditions...")
  tryCatch(plot_precision_recall_highlight(), error = function(e)
    methods_log("S07_ERR", paste("Fig.5:", e$message)))
  
  # ---- Step 4: GEO Table (if available) ----
  methods_log("S07_MAIN", "Step 4: Cross-dataset table (if GEO data exists)...")
  tryCatch(generate_table3(), error = function(e)
    methods_log("S07_ERR", paste("Table 3:", e$message)))
  
  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  cat("\n================================================================\n")
  cat(sprintf("  ALL DONE. Total time: %.1f minutes\n", elapsed))
  cat(sprintf("  Figures: %s\n", FIG_DIR_METHODS))
  cat(sprintf("  Tables:  %s\n", TAB_DIR_METHODS))
  cat(sprintf("  Data:    %s\n", BENCH_DIR))
  cat("================================================================\n")
  
  methods_log("S07_MAIN", sprintf("generate_all_v2 complete (%.1f min)", elapsed))
}


# ==============================================================================
# 直接运行提示
# ==============================================================================

if (sys.nframe() == 0) {
  cat("\n")
  cat("================================================================\n")
  cat("  S07_cross_dataset_summary.R v2\n")
  cat("================================================================\n")
  cat("\n")
  cat("  # 一键运行全部 (推荐):\n")
  cat("  generate_all_v2()\n")
  cat("\n")
  cat("  # 跳过Nogueira (已有结果时):\n")
  cat("  generate_all_v2(skip_nogueira = TRUE)\n")
  cat("\n")
  cat("  # 分步:\n")
  cat("  batch_nogueira()                    # Nogueira (~10-20 min)\n")
  cat("  build_ablation_paper_table()        # 消融汇总 (<1 min)\n")
  cat("  plot_simulation_benchmark()         # Fig.2\n")
  cat("  plot_ablation_results_v2()          # Fig.3\n")
  cat("  plot_stability_nogueira()           # Fig.4\n")
  cat("  plot_precision_recall_highlight()   # Fig.5\n")
  cat("  generate_table3()                   # Table 3 (GEO)\n")
  cat("\n")
  cat("  # 诊断:\n")
  cat("  diagnose_stab_rds(file.path(RUN_DIR, 'SIM_xxx'))  # check matrix\n")
  cat("================================================================\n")
}
