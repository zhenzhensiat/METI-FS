#!/usr/bin/env Rscript
# ==============================================================================
# S06_ablation_study.R — Pipeline组件消融实验
#
# 目的：
#   系统量化METI-FS pipeline每个组件的贡献。
#   逐层移除组件，观察最终marker set的变化：
#     - 模拟数据: precision/recall/FDR的变化
#     - 真实数据: marker数量/稳定性/重叠度的变化
#
# 消融方案（7种配置）：
#   A0: FULL         — 完整pipeline（基准线）
#   A1: -maSigPro    — 跳过maSigPro交互过滤，候选池 = WGCNA key genes
#   A2: -WGCNA       — 跳过WGCNA，候选池 = maSigPro filtered genes
#   A3: -EffectSize  — 跳过效应量过滤，直接用maSigPro∩WGCNA
#   A4: -GapUnion    — 用固定阈值(freq≥0.7, ≥2 algos)替代gap-union
#   A5: -PPI         — 只用ML markers，不加PPI hub union
#   A6: -Bootstrap   — 用单次ML运行替代100次bootstrap stability
#
# 实现策略：
#   不修改原pipeline脚本。每种消融在09A之后介入：
#   - A1/A2/A3: 修改candidate_pool.rds后重跑09C→09D→09F→10
#   - A4: 替换09D的选择逻辑
#   - A5: 修改10_integration的union逻辑
#   - A6: 运行单次ML替代bootstrap
#
# 输入：
#   - 已完成的FULL pipeline运行（data/目录下有所有中间文件）
#   - 模拟数据时：data/ground_truth.rds
#
# 输出：
#   BENCH_DIR/ablation_{run_id}.rds     — 完整消融结果
#   BENCH_DIR/ablation_{run_id}.csv     — 消融汇总表
#
# 文献依据：
#   - pipeComp (Germain et al. 2020, Genome Biology 21:227):
#     pipeline benchmark框架, LOCO (Leave-One-Component-Out) 消融范式
#     "enabling the exploration of combinations of parameters and of the
#      robustness of methods to various changes in other parts of a pipeline"
#   - Mangul et al. 2019, Nature Biotechnology 37:1127-1133:
#     benchmarking best practices, "systematic evaluation of bioinformatics
#     software should include ablation studies"
#   - Spooner et al. 2023, BMC Bioinformatics 24:9:
#     data-driven thresholding vs fixed threshold comparison (A4消融的对照依据)
#   - Nogueira, Sechidis & Brown 2018, JMLR 18(174):1-54:
#     stability as evaluation metric for ablation (引用年份修正: 2018非2017)
# ==============================================================================

# ---- 0. 依赖 ----
if (file.exists("S_config.R")) {
  source("S_config.R")

  
}
if (file.exists("S01_simulation_engine.R")) {
  source("S01_simulation_engine.R")
} else if (file.exists(file.path(METHODS_SCRIPTS, "S01_simulation_engine.R"))) {
  source(file.path(METHODS_SCRIPTS, "S01_simulation_engine.R"))
}

suppressPackageStartupMessages({
  library(glmnet)
  library(randomForest)
  library(caret)
})


# ==============================================================================
# 消融配置注册表
# ==============================================================================

ABLATION_CONFIGS <- list(
  A0_FULL = list(
    label       = "Full pipeline",
    description = "Complete METI-FS (baseline)",
    skip        = character(0)
  ),
  A1_no_maSigPro = list(
    label       = "-maSigPro interaction",
    description = "Skip maSigPro interaction filter; pool = WGCNA key genes only",
    skip        = "masigpro_interaction"
  ),
  A2_no_WGCNA = list(
    label       = "-WGCNA",
    description = "Skip WGCNA; pool = maSigPro filtered genes only",
    skip        = "wgcna"
  ),
  A3_no_EffectSize = list(
    label       = "-Effect size",
    description = "Skip lfcThreshold=1 test; pool = maSigPro ∩ WGCNA without effect filter",
    skip        = "effect_size"
  ),
  A4_no_GapUnion = list(
    label       = "-Gap-union",
    description = "Replace gap thresholding with fixed threshold (freq>=0.7, >=2 algos)",
    skip        = "gap_union"
  ),
  A5_no_PPI = list(
    label       = "-PPI hub",
    description = "Final = ML markers only, no PPI hub union",
    skip        = "ppi_union"
  ),
  A6_no_Bootstrap = list(
    label       = "-Bootstrap stability",
    description = "Single ML run instead of 100× bootstrap",
    skip        = "bootstrap"
  ),
  # [v2] Pool size sensitivity — Mangul et al. 2019 Nat Biotechnol
  A7_pool_1000 = list(
    label       = "Pool cap = 1000",
    description = "Truncate candidate pool to top 1000 by padj, then re-run ML",
    skip        = "pool_resize",
    pool_cap    = 1000
  ),
  A8_pool_3000 = list(
    label       = "Pool cap = 3000",
    description = "Expand candidate pool cap to 3000, then re-run ML",
    skip        = "pool_resize",
    pool_cap    = 3000
  ),
  # RF weighted ensemble (based on Nogueira index ratios)
  A9_RF_weighted = list(
    label       = "RF weighted ensemble",
    description = "Weighted RF inclusion by Nogueira ratio instead of binary exclusion",
    skip        = "rf_binary_exclusion"
  )
)


# ==============================================================================
# 核心函数: run_ablation()
# ==============================================================================

#' 在一个已完成的pipeline运行上执行全套消融实验
#'
#' @param run_dir 已完成FULL pipeline的运行目录
#' @param run_id  运行标识符
#' @param configs 要执行的消融配置（默认全部7种）
#' @return list of ablation results
run_ablation <- function(run_dir, run_id = basename(run_dir),
                          configs = names(ABLATION_CONFIGS)) {

  data_dir <- file.path(run_dir, "data")
  gt_file  <- file.path(data_dir, "ground_truth.rds")
  has_gt   <- file.exists(gt_file)

  methods_log("S06_ABLATION", sprintf(
    "=== Ablation study: %s (%d configs, ground_truth=%s) ===",
    run_id, length(configs), has_gt))

  # 加载共享上游数据
  upstream <- load_upstream_data(data_dir)
  gt <- if (has_gt) readRDS(gt_file) else NULL

  results <- list()

  for (cfg_name in configs) {
    cfg <- ABLATION_CONFIGS[[cfg_name]]
    methods_log("S06_ABLATION", sprintf("--- %s: %s ---", cfg_name, cfg$label))

    tryCatch({
      res <- run_single_ablation(
        cfg_name  = cfg_name,
        cfg       = cfg,
        upstream  = upstream,
        data_dir  = data_dir,
        gt        = gt
      )
      results[[cfg_name]] <- res
      methods_log("S06_ABLATION", sprintf(
        "  → %d final markers%s",
        res$n_final,
        if (!is.null(res$perf)) sprintf(
          " (P=%.3f, R=%.3f, F1=%.3f)",
          res$perf$precision, res$perf$recall, res$perf$F1
        ) else ""
      ))
    }, error = function(e) {
      methods_log("S06_ABLATION", sprintf("  ERROR: %s", e$message))
      results[[cfg_name]] <<- list(error = e$message, n_final = NA)
    })
  }

  # 汇总表
  summary_df <- build_ablation_summary(results, has_gt)

  # 保存
  output <- list(
    run_id   = run_id,
    run_dir  = run_dir,
    has_gt   = has_gt,
    results  = results,
    summary  = summary_df
  )

  out_rds <- file.path(BENCH_DIR, sprintf("ablation_%s.rds", run_id))
  out_csv <- file.path(BENCH_DIR, sprintf("ablation_%s.csv", run_id))
  saveRDS(output, out_rds)
  write.csv(summary_df, out_csv, row.names = FALSE)

  methods_log("S06_ABLATION", sprintf("Saved: %s, %s", basename(out_rds), basename(out_csv)))

  # 打印汇总
  print_ablation_summary(summary_df, has_gt)

  return(output)
}


# ==============================================================================
# 上游数据加载
# ==============================================================================

load_upstream_data <- function(data_dir) {

  up <- list()

  # 09A候选池（含各层中间数据）
  pool_file <- file.path(data_dir, "candidate_pool.rds")
  if (file.exists(pool_file)) {
    pool <- readRDS(pool_file)
    up$candidate_pool       <- pool$candidate_pool
    up$masigpro_ensembl     <- pool$masigpro_ensembl_filtered
    up$wgcna_key_genes      <- pool$wgcna_key_genes
    up$n_pool_pre_effect    <- pool$n_pool_pre_effect
    up$n_pool_post_effect   <- pool$n_pool_post_effect
    up$effect_filter_applied <- pool$effect_filter_applied
    # maSigPro全集（交互过滤前）
    up$masigpro_all         <- pool$masigpro_ensembl_all %||%
                                pool$masigpro_ensembl_filtered
  }

  # 09C stability selection matrices
  stab_file <- file.path(data_dir, "ml_stability_selection.rds")
  if (file.exists(stab_file)) {
    stab <- readRDS(stab_file)
    up$stab <- stab
  }

  # 09D gap-union结果（A0 baseline）
  gap_file <- file.path(data_dir, "ml_gap_union.rds")
  if (file.exists(gap_file)) {
    up$gap_union <- readRDS(gap_file)
  }

  # 09F PPI hub结果
  ppi_file <- file.path(data_dir, "ppi_hub_selection.rds")
  if (file.exists(ppi_file)) {
    up$ppi <- readRDS(ppi_file)
  }

  # TPM + gene_anno（某些消融需要重跑ML）
  tpm_file <- file.path(data_dir, "tpm_filtered.rds")
  if (file.exists(tpm_file)) up$tpm <- readRDS(tpm_file)

  anno_file <- file.path(data_dir, "gene_annotation.rds")
  if (file.exists(anno_file)) up$gene_anno <- readRDS(anno_file)

  sample_file <- file.path(data_dir, "sample_info.rds")
  if (file.exists(sample_file)) up$sample_info <- readRDS(sample_file)

  # DEG results（用于效应量过滤消融）
  deg_file <- file.path(data_dir, "deg_results.rds")
  if (file.exists(deg_file)) up$deg <- readRDS(deg_file)

  return(up)
}


# ==============================================================================
# 单次消融执行
# ==============================================================================

run_single_ablation <- function(cfg_name, cfg, upstream, data_dir, gt) {

  skip <- cfg$skip

  # ====== Step 1: 构建候选池 ======
  pool <- build_ablated_pool(skip, upstream, cfg = cfg, data_dir = data_dir)

  # ====== Step 2: ML feature selection ======
  ml_result <- run_ablated_ml(skip, pool, upstream, data_dir)

  # ====== Step 3: PPI hub ======
  ppi_result <- get_ablated_ppi(skip, upstream)

  # ====== Step 4: Final union ======
  final_ids <- build_ablated_final(skip, ml_result, ppi_result)

  # ====== Step 5: Evaluate ======
  perf <- NULL
  if (!is.null(gt)) {
    perf <- evaluate_selection(final_ids, gt, "TRUE_TEMPORAL")
  }

  return(list(
    config    = cfg_name,
    label     = cfg$label,
    n_pool    = length(pool),
    n_ml      = length(ml_result$gene_ids),
    n_ppi     = length(ppi_result$gene_ids),
    n_final   = length(final_ids),
    final_ids = final_ids,
    ml_ids    = ml_result$gene_ids,
    ppi_ids   = ppi_result$gene_ids,
    perf      = perf
  ))
}


# ==============================================================================
# 消融子函数：候选池构建
# ==============================================================================

build_ablated_pool <- function(skip, upstream, cfg = NULL, data_dir = NULL) {

  masigpro_genes <- upstream$masigpro_ensembl
  wgcna_genes    <- upstream$wgcna_key_genes
  full_pool      <- upstream$candidate_pool

  if ("masigpro_interaction" %in% skip) {
    # A1: 不用maSigPro过滤，只用WGCNA key module genes
    pool <- wgcna_genes
  } else if ("wgcna" %in% skip) {
    # A2: 不用WGCNA，只用maSigPro filtered genes
    pool <- masigpro_genes
  } else if ("effect_size" %in% skip) {
    # A3: 用maSigPro∩WGCNA但跳过效应量过滤
    pool <- intersect(masigpro_genes, wgcna_genes)
  } else if ("pool_resize" %in% skip && !is.null(cfg) && !is.null(cfg$pool_cap)) {
    # [v2] A7/A8: 改变候选池大小上限
    pool_cap <- cfg$pool_cap
    pool <- full_pool  # 从完整候选池开始
    
    if (length(pool) > pool_cap) {
      # 尝试按LRT padj截断
      deg_file <- file.path(data_dir, "deg_results.rds")
      if (!is.null(data_dir) && file.exists(deg_file)) {
        all_results <- readRDS(deg_file)
        lrt_full <- all_results$lrt_interaction_sig
        id_col <- if ("ensembl_id" %in% colnames(lrt_full)) "ensembl_id" else "ensembl_gene_id"
        top_ids <- lrt_full[lrt_full[[id_col]] %in% pool, ]
        top_ids <- top_ids[order(top_ids$padj), ]
        pool <- head(top_ids[[id_col]], pool_cap)
      } else {
        # Fallback: 随机截断 (保留信息)
        pool <- sample(pool, pool_cap)
      }
    }
    # 如果pool本身 < pool_cap，不做改变
  } else {
    # A0/A4/A5/A6: 使用完整候选池
    pool <- full_pool
  }

  pool <- pool[!is.na(pool)]
  return(pool)
}


# ==============================================================================
# 消融子函数：ML特征选择
# ==============================================================================

run_ablated_ml <- function(skip, pool, upstream, data_dir) {

  # A6: 单次ML（不bootstrap）→ 用单次运行的结果
  if ("bootstrap" %in% skip) {
    return(run_single_ml(pool, upstream))
  }

  # 如果候选池与原始不同（A1/A2/A3），需要重跑bootstrap
  # 但重跑很贵（每次30-60分钟），所以用近似：
  # 从原始stability矩阵中提取子集，重新计算频率
  if (!identical(sort(pool), sort(upstream$candidate_pool))) {
    return(recompute_stability_subset(pool, upstream, skip))
  }

  # A4: 候选池不变，但改变09D的选择逻辑
  if ("gap_union" %in% skip) {
    return(apply_fixed_threshold(upstream))
  }

  # A9: RF weighted ensemble — 基于Nogueira指数比率加权保留RF基因
  if ("rf_binary_exclusion" %in% skip) {
    return(apply_rf_weighted_ensemble(upstream))
  }

  # A0/A5: 使用原始ML结果
  ml_ids <- upstream$gap_union$final_gene_ids
  return(list(gene_ids = ml_ids))
}


#' A6消融：单次ML运行（无bootstrap）
run_single_ml <- function(pool, upstream) {

  tpm <- upstream$tpm
  si  <- upstream$sample_info

  pool <- pool[pool %in% rownames(tpm)]
  if (length(pool) < 10) return(list(gene_ids = character(0)))

  tpm_sub <- t(log2(tpm[pool, ] + 1))
  labels  <- factor(si$Treatment, levels = c("Control", "Induced"))

  selected <- character(0)

  # LASSO (single run)
  tryCatch({
    cv_fit <- cv.glmnet(tpm_sub, labels, family = "binomial",
                         alpha = 1, nfolds = min(10, nrow(tpm_sub)))
    coefs <- coef(cv_fit, s = "lambda.min")[-1, ]
    lasso_genes <- names(coefs)[coefs != 0]
    selected <- union(selected, lasso_genes)
  }, error = function(e) {})

  # SVM-RFE (single run)
  tryCatch({
    ctrl <- rfeControl(functions = caretFuncs, method = "cv",
                       number = min(5, nrow(tpm_sub)))
    rfe_result <- rfe(tpm_sub, labels, sizes = c(5, 10, 20, 50),
                      rfeControl = ctrl, method = "svmRadial")
    svm_genes <- predictors(rfe_result)
    selected <- union(selected, svm_genes)
  }, error = function(e) {})

  return(list(gene_ids = selected))
}


#' 候选池变更时：从原始selection矩阵中提取子集重算频率
recompute_stability_subset <- function(pool, upstream, skip) {

  stab <- upstream$stab
  if (is.null(stab) || !all(c("lasso_selection_matrix", "svm_selection_matrix") %in% names(stab))) {
    # 无selection矩阵，退回到原始结果中筛选
    orig_ids <- upstream$gap_union$final_gene_ids
    return(list(gene_ids = intersect(orig_ids, pool)))
  }

  # 取出各算法的selection matrix列（基因）中属于新pool的子集
  freq_threshold <- stab$params$freq_threshold  # 0.7

  algo_stable <- list()
  for (algo in c("lasso", "rf", "svm")) {
    mat <- stab[[paste0(algo, "_selection_matrix")]]
    if (is.null(mat)) next
    # mat: bootstrap × genes (binary)
    pool_genes <- intersect(colnames(mat), pool)
    if (length(pool_genes) == 0) {
      algo_stable[[algo]] <- character(0)
      next
    }
    freqs <- colMeans(mat[, pool_genes, drop = FALSE])
    algo_stable[[algo]] <- names(freqs)[freqs >= freq_threshold]
  }

  # 简化union: ≥1个算法stable即入选
  all_stable <- unique(unlist(algo_stable))
  return(list(gene_ids = all_stable))
}


#' A9消融: RF加权集成 — 以Nogueira指数比率加权保留RF基因
#'
#' 替代二元排除: RF基因按 RF_Nogueira/max(LASSO_Nogueira, SVM_Nogueira) 加权
#' 权重 >= 0.5 的基因保留, 加入最终ML gene set
apply_rf_weighted_ensemble <- function(upstream) {

  stab <- upstream$stab
  if (is.null(stab)) return(list(gene_ids = character(0)))

  # 获取各算法频率向量
  freq_lasso <- stab[["lasso_freq"]]
  freq_svm   <- stab[["svm_freq"]]
  freq_rf    <- stab[["rf_freq"]]

  # 获取各算法选择矩阵用于计算Nogueira指数 (简化: 用频率向量的SD估计稳定性)
  # 完整Nogueira需selection matrix; 这里用bootstrap frequency的方差近似
  compute_approx_nogueira <- function(freq_vec) {
    if (is.null(freq_vec) || length(freq_vec) < 2) return(0)
    # 近似: 基于频率的稳定性 proxy (越高频越稳定, 但需要chance correction)
    # 简化版: mean frequency of top-quartile genes
    top_genes <- freq_vec[freq_vec >= quantile(freq_vec, 0.75, na.rm = TRUE)]
    if (length(top_genes) < 2) return(0)
    mean(top_genes, na.rm = TRUE)
  }

  # 尝试从stab获取完整selection matrices计算Nogueira
  nogueira_lasso <- NA; nogueira_svm <- NA; nogueira_rf <- NA

  for (algo in c("lasso", "svm", "rf")) {
    mat <- stab[[paste0(algo, "_selection_matrix")]]
    if (!is.null(mat) && nrow(mat) >= 2 && ncol(mat) >= 2) {
      # 使用stabm包计算Nogueira指数 (如果可用)
      nogueira_val <- tryCatch({
        if (requireNamespace("stabm", quietly = TRUE)) {
          stabm::stabilityNogueira(mat)
        } else {
          compute_approx_nogueira(stab[[paste0(algo, "_freq")]])
        }
      }, error = function(e) compute_approx_nogueira(stab[[paste0(algo, "_freq")]]))
      assign(paste0("nogueira_", algo), nogueira_val)
    } else {
      assign(paste0("nogueira_", algo),
             compute_approx_nogueira(stab[[paste0(algo, "_freq")]]))
    }
  }

  # 计算RF权重
  max_ref <- max(nogueira_lasso, nogueira_svm, na.rm = TRUE)
  if (is.na(max_ref) || max_ref == 0) {
    rf_weight <- 0
  } else {
    rf_weight <- nogueira_rf / max_ref
  }

  # 获取各算法选中的基因
  lasso_genes <- if (!is.null(freq_lasso)) names(freq_lasso)[freq_lasso >= 0.20] else character(0)
  svm_genes   <- if (!is.null(freq_svm))   names(freq_svm)[freq_svm >= 0.20]   else character(0)

  # RF基因: 按权重过滤
  if (rf_weight >= 0.5 && !is.null(freq_rf)) {
    # 保留高频RF基因, 按权重比例降低阈值
    rf_threshold <- 0.20 + (1 - rf_weight) * 0.3  # 权重越低阈值越高
    rf_genes <- names(freq_rf)[freq_rf >= rf_threshold]
  } else {
    rf_genes <- character(0)
  }

  # Union: LASSO + SVM + weighted RF
  all_genes <- unique(c(lasso_genes, svm_genes, rf_genes))

  return(list(
    gene_ids = all_genes,
    rf_weight = rf_weight,
    rf_genes = rf_genes
  ))
}

#' A4消融：固定阈值替代gap-union
apply_fixed_threshold <- function(upstream) {

  stab <- upstream$stab
  if (is.null(stab)) return(list(gene_ids = character(0)))

  freq_threshold <- 0.7
  min_algos <- 2

  stable_per_algo <- list()
  for (algo in c("lasso", "rf", "svm")) {
    freq_vec <- stab[[paste0(algo, "_freq")]]
    if (is.null(freq_vec)) next
    stable_per_algo[[algo]] <- names(freq_vec)[freq_vec >= freq_threshold]
  }

  # 计算每个基因被多少算法选中
  all_genes <- unique(unlist(stable_per_algo))
  if (length(all_genes) == 0) return(list(gene_ids = character(0)))

  n_algos <- sapply(all_genes, function(g) {
    sum(sapply(stable_per_algo, function(s) g %in% s))
  })

  selected <- names(n_algos)[n_algos >= min_algos]
  return(list(gene_ids = selected))
}


# ==============================================================================
# 消融子函数：PPI hub
# ==============================================================================

get_ablated_ppi <- function(skip, upstream) {

  if ("ppi_union" %in% skip) {
    return(list(gene_ids = character(0)))
  }

  if (!is.null(upstream$ppi) && !is.null(upstream$ppi$hub_genes)) {
    hub_df <- upstream$ppi$hub_genes
    ids <- if ("ensembl_id" %in% colnames(hub_df)) hub_df$ensembl_id else hub_df$symbol
    return(list(gene_ids = ids))
  }

  return(list(gene_ids = character(0)))
}


# ==============================================================================
# 消融子函数：最终union
# ==============================================================================

build_ablated_final <- function(skip, ml_result, ppi_result) {

  ml_ids  <- ml_result$gene_ids
  ppi_ids <- ppi_result$gene_ids

  final <- union(ml_ids, ppi_ids)
  final <- final[!is.na(final) & final != ""]
  return(final)
}


# ==============================================================================
# 汇总与输出
# ==============================================================================

build_ablation_summary <- function(results, has_gt) {

  rows <- list()
  for (cfg_name in names(results)) {
    res <- results[[cfg_name]]
    if (is.null(res$n_final) || is.na(res$n_final)) next

    row <- data.frame(
      config   = cfg_name,
      label    = res$label %||% cfg_name,
      n_pool   = res$n_pool %||% NA,
      n_ml     = res$n_ml %||% NA,
      n_ppi    = res$n_ppi %||% NA,
      n_final  = res$n_final,
      stringsAsFactors = FALSE
    )

    if (has_gt && !is.null(res$perf)) {
      row$precision <- res$perf$precision
      row$recall    <- res$perf$recall
      row$F1        <- res$perf$F1
      row$FDR       <- res$perf$FDR
      row$TP        <- res$perf$TP
      row$FP        <- res$perf$FP
      row$FN        <- res$perf$FN
    }

    rows[[cfg_name]] <- row
  }

  do.call(rbind, rows)
}


print_ablation_summary <- function(summary_df, has_gt) {

  cat("\n")
  cat("================================================================\n")
  cat("  ABLATION STUDY RESULTS\n")
  cat("================================================================\n\n")

  if (has_gt) {
    cat(sprintf("  %-25s %5s %5s %5s %7s %7s %7s %7s\n",
                "Config", "Pool", "ML", "Final", "Prec", "Recall", "F1", "FDR"))
    cat(paste(rep("-", 78), collapse = ""), "\n")
    for (i in seq_len(nrow(summary_df))) {
      r <- summary_df[i, ]
      cat(sprintf("  %-25s %5d %5d %5d %7.3f %7.3f %7.3f %7.3f\n",
                  r$label, r$n_pool, r$n_ml, r$n_final,
                  r$precision, r$recall, r$F1, r$FDR))
    }
  } else {
    cat(sprintf("  %-25s %5s %5s %5s %5s\n",
                "Config", "Pool", "ML", "PPI", "Final"))
    cat(paste(rep("-", 52), collapse = ""), "\n")
    for (i in seq_len(nrow(summary_df))) {
      r <- summary_df[i, ]
      cat(sprintf("  %-25s %5d %5d %5d %5d\n",
                  r$label, r$n_pool, r$n_ml, r$n_ppi, r$n_final))
    }
  }

  # Delta vs baseline
  if ("A0_FULL" %in% summary_df$config && has_gt) {
    cat("\n  [Delta vs Full pipeline]\n")
    base <- summary_df[summary_df$config == "A0_FULL", ]
    for (i in seq_len(nrow(summary_df))) {
      r <- summary_df[i, ]
      if (r$config == "A0_FULL") next
      dp <- r$precision - base$precision
      dr <- r$recall - base$recall
      df1 <- r$F1 - base$F1
      cat(sprintf("    %-25s ΔPrec=%+.3f  ΔRecall=%+.3f  ΔF1=%+.3f  Δn=%+d\n",
                  r$label, dp, dr, df1, r$n_final - base$n_final))
    }
  }

  cat("\n================================================================\n")
}


# ==============================================================================
# 批量消融（模拟数据全场景）
# ==============================================================================

#' 对所有模拟场景执行消融实验
ablate_all_simulations <- function(
    sim_run_dir = file.path(RUN_DIR, "simulations"),
    configs = names(ABLATION_CONFIGS)) {

  all_dirs <- list.dirs(sim_run_dir, recursive = FALSE, full.names = TRUE)
  run_dirs <- all_dirs[sapply(all_dirs, function(d) {
    file.exists(file.path(d, "data", "ground_truth.rds")) &&
    file.exists(file.path(d, "data", "ml_gap_union.rds"))
  })]

  methods_log("S06_BATCH", sprintf("Found %d completed simulation runs", length(run_dirs)))

  all_results <- list()
  for (rd in run_dirs) {
    run_name <- basename(rd)
    tryCatch({
      all_results[[run_name]] <- run_ablation(rd, run_id = run_name, configs = configs)
    }, error = function(e) {
      methods_log("S06_BATCH", sprintf("ERROR in %s: %s", run_name, e$message))
    })
  }

  # 合并所有消融汇总
  all_summaries <- do.call(rbind, lapply(names(all_results), function(nm) {
    df <- all_results[[nm]]$summary
    df$run_id <- nm
    # 从run_id解析场景参数
    parts <- strsplit(nm, "_")[[1]]
    if (length(parts) >= 4) {
      df$snr            <- parts[1]
      df$sample_size    <- parts[2]
      df$marker_density <- parts[3]
      df$repeat_id      <- as.integer(gsub("rep", "", parts[4]))
    }
    df
  }))

  out_file <- file.path(BENCH_DIR, "ablation_all_simulations.csv")
  write.csv(all_summaries, out_file, row.names = FALSE)
  methods_log("S06_BATCH", sprintf("Saved: %s (%d rows)", basename(out_file), nrow(all_summaries)))

  return(invisible(all_results))
}


# ==============================================================================
# 直接运行时的使用说明
# ==============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

if (sys.nframe() == 0) {
  cat("\n")
  cat("================================================================\n")
  cat("  S06_ablation_study.R — Usage\n")
  cat("================================================================\n")
  cat("\n")
  cat("  # 单个数据集消融:\n")
  cat("  res <- run_ablation(file.path(RUN_DIR, 'GEO_GSE307424_Lung'))\n")
  cat("\n")
  cat("  # 批量模拟数据消融:\n")
  cat("  ablate_all_simulations()\n")
  cat("\n")
  cat("  # 查看消融配置:\n")
  cat("  str(ABLATION_CONFIGS)\n")
  cat("================================================================\n")
}
