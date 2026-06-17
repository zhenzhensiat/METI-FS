#!/usr/bin/env Rscript
# ==============================================================================
# S04_run_pipeline_wrapper.R — GEO数据集pipeline适配运行器

# ---- 关闭报错后进入Browse调试模式 ----
options(error = NULL)
#
# 目的：在不修改原pipeline脚本的前提下，将外部GEO数据集跑通pipeline。
#       Preserves time labels, factor levels from the core pipeline.
#       颜色方案、maSigPro参数等针对MSC 4d/7d/14d/21d设计硬编码，
#       本wrapper通过"source后覆盖"策略解决。
#
# 适配方式：
#   Level 1: source 00_setup.R 后覆盖 PARAMS / COLORS
#   Level 2: source 01_data_import.R 后修正 sample_info Time/Group factors
#   Level 3: 每个下游脚本 source 后修正其硬编码的 timepoints 变量
#   Level 4: 基因ID适配（symbol vs Ensembl）
#
# 用法：
#   1. 在本脚本中选择DATASET_ID（见下方注册表）
#   2. 直接 source("S04_run_pipeline_wrapper.R") 或逐phase执行
#   3. 不改原pipeline任何一行代码
#
# 关键原则：
#   - 原pipeline脚本只读引用
#   - 所有适配在本文件中完成
#   - 方法学文章不使用MSC数据
# ==============================================================================

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  STEP 0: 选择数据集 + 加载配置                                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ---- 选择要跑的数据集 ----
# 方式1：在source之前设置 DATASET_ID <- "GSE197067"
# 方式2：直接修改下面的默认值
if (!exists("DATASET_ID")) {
  DATASET_ID <- "GSE307424"   # 默认：肺癌SMARCA2，18样本，3时间点
}

# ---- 数据集注册表（每个数据集的适配参数） ----
DATASET_PROFILES <- list(

  GSE307424 = list(
    project_dir   = file.path(RUN_DIR, "GEO_GSE307424_Lung"),
    prefix        = "Lung",
    # 时间配置
    time_values   = c(6, 48, 72),              # 数值型时间（小时）
    time_labels   = c("6h", "48h", "72h"),     # 显示标签
    time_suffix   = "d",                       # 样本名中的后缀（Lung6d1 → 保持"d"）
    time_ref      = "6h",                      # DESeq2参考水平（最早时间点）
    # maSigPro
    masigpro_degree  = 2,   # 3个时间点 → degree = nTP - 1 = 2
    masigpro_k       = 6,   # 聚类数减少（基因少）
    masigpro_min_obs = 3,   # p.vector()默认值; User's Guide示例用20但需≤n_samples
    # 基因ID
    gene_id_type  = "symbol",  # 非Ensembl！
    # WGCNA
    wgcna_note    = "18 samples → FAQ power table fallback likely",
    # 可视化
    colors_timepoint = c("6h"  = "#3C5488",
                         "48h" = "#F39B7F",
                         "72h" = "#E64B35"),
    colors_group = c("Induced_6h"  = "#DC0000", "Control_6h"  = "#3C5488",
                     "Induced_48h" = "#E64B35", "Control_48h" = "#00A087",
                     "Induced_72h" = "#B09C85", "Control_72h" = "#7E6148"),
    shape_timepoint = c("6h" = 16, "48h" = 17, "72h" = 15)
  ),

  GSE197067 = list(
    project_dir   = file.path(RUN_DIR, "GEO_GSE197067_Tcell"),
    prefix        = "Tcell",
    # 时间配置（排除0h：0h只有Control组，无Induced，交互项无法估计）
    time_values   = c(6, 12, 24, 48, 72),       # 排除0h
    time_labels   = c("6h", "12h", "24h", "48h", "72h"),
    time_suffix   = "d",
    time_ref      = "6h",                        # 最早时间点作参考
    # 排除0h的样本名模式（data_raw中Tcell0dC1-4需要在导入时去掉）
    exclude_samples_pattern = "^Tcell0d",        # 匹配Tcell0dC1等
    # maSigPro
    masigpro_degree  = 4,   # 5个时间点 → degree = nTP - 1 = 4
    masigpro_k       = 9,
    masigpro_min_obs = 20,  # 40样本，沿用User's Guide推荐值
    # 基因ID
    gene_id_type  = "ensembl",
    # WGCNA
    wgcna_note    = "40 samples → adequate for WGCNA",
    # 可视化
    colors_timepoint = c("6h"  = "#3C5488",
                         "12h" = "#00A087",
                         "24h" = "#F39B7F",
                         "48h" = "#E64B35",
                         "72h" = "#B09C85"),
    colors_group = c("Induced_6h"   = "#DC0000", "Control_6h"   = "#3C5488",
                     "Induced_12h"  = "#F39B7F", "Control_12h"  = "#4DBBD5",
                     "Induced_24h"  = "#E64B35", "Control_24h"  = "#00A087",
                     "Induced_48h"  = "#B09C85", "Control_48h"  = "#7E6148",
                     "Induced_72h"  = "#DC9FB4", "Control_72h"  = "#7570B3"),
    shape_timepoint = c("6h" = 16, "12h" = 17, "24h" = 15, "48h" = 18, "72h" = 8)
  ),

  GSE236646 = list(
    project_dir   = file.path(RUN_DIR, "GEO_GSE236646_NPC"),
    prefix        = "NPC",
    time_values   = c(3, 5, 7),
    time_labels   = c("3d", "5d", "7d"),
    time_suffix   = "d",
    time_ref      = "3d",
    masigpro_degree  = 2,
    masigpro_k       = 6,
    masigpro_min_obs = 3,
    gene_id_type  = "entrez",
    wgcna_note    = "17 samples -> power table fallback likely",
    colors_timepoint = c("3d" = "#3C5488", "5d" = "#F39B7F", "7d" = "#E64B35"),
    colors_group = c("Induced_3d"  = "#DC0000", "Control_3d"  = "#3C5488",
                     "Induced_5d"  = "#E64B35", "Control_5d"  = "#00A087",
                     "Induced_7d"  = "#B09C85", "Control_7d"  = "#7E6148"),
    shape_timepoint = c("3d" = 16, "5d" = 17, "7d" = 15)
  ),

  GSE150411 = list(
    project_dir   = file.path(RUN_DIR, "GEO_GSE150411_Chon"),
    prefix        = "Chon",
    time_values   = c(3, 6, 18),
    time_labels   = c("3h", "6h", "18h"),
    time_suffix   = "d",
    time_ref      = "3h",
    masigpro_degree  = 2,
    masigpro_k       = 6,
    masigpro_min_obs = 3,
    gene_id_type  = "entrez",
    wgcna_note    = "18 samples -> power table fallback likely",
    colors_timepoint = c("3h" = "#3C5488", "6h" = "#F39B7F", "18h" = "#E64B35"),
    colors_group = c("Induced_3h"  = "#DC0000", "Control_3h"  = "#3C5488",
                     "Induced_6h"  = "#E64B35", "Control_6h"  = "#00A087",
                     "Induced_18h" = "#B09C85", "Control_18h" = "#7E6148"),
    shape_timepoint = c("3h" = 16, "6h" = 17, "18h" = 15)
  )
)

# ---- 加载当前数据集的profile ----
PROFILE <- DATASET_PROFILES[[DATASET_ID]]
if (is.null(PROFILE)) stop("Unknown DATASET_ID: ", DATASET_ID)

cat("================================================================\n")
cat("  METI-FS Pipeline Wrapper\n")
cat(sprintf("  Dataset: %s\n", DATASET_ID))
cat(sprintf("  Project: %s\n", PROFILE$project_dir))
cat(sprintf("  Times:   %s\n", paste(PROFILE$time_labels, collapse = ", ")))
cat(sprintf("  Degree:  %d\n", PROFILE$masigpro_degree))
cat(sprintf("  GeneID:  %s\n", PROFILE$gene_id_type))
cat("================================================================\n")


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  HELPER: 通用修正函数                                                      ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

#' 修正sample_info的Time和Group因子水平
#' 在每次source pipeline脚本后调用（因为部分脚本会重新从files加载sample_info）
fix_sample_info_factors <- function(si) {
  actual_time_labels <- PROFILE$time_labels   # c("6h","48h","72h")
  actual_time_values <- PROFILE$time_values   # c(6, 48, 72)
  
  # 01_data_import.R 把 Lung6d1 解析为 time_num=6, time_label="6d"
  # 但实际时间单位可能是小时(h)不是天(d)
  # 用 time_num → PROFILE$time_labels 建映射表
  time_map <- setNames(actual_time_labels, as.character(actual_time_values))
  # e.g. "6" → "6h", "48" → "48h", "72" → "72h"
  
  # 重建 time_label
  si$time_label <- time_map[as.character(si$time_num)]
  
  # 修正 Time factor
  si$Time <- factor(si$time_label, levels = actual_time_labels)
  
  # 修正 Group factor
  group_levels_ordered <- c()
  for (tl in actual_time_labels) {
    group_levels_ordered <- c(group_levels_ordered,
                              paste0("Control_", tl),
                              paste0("Induced_", tl))
  }
  si$Group <- factor(paste(si$treatment, si$time_label, sep = "_"),
                     levels = group_levels_ordered)
  
  return(si)
}

#' 覆盖COLORS中的timepoint和group映射
fix_colors <- function() {
  COLORS$timepoint <<- PROFILE$colors_timepoint
  COLORS$group     <<- PROFILE$colors_group
}

#' 覆盖PARAMS中的maSigPro参数
fix_params <- function() {
  PARAMS$masigpro_degree <<- PROFILE$masigpro_degree
  PARAMS$masigpro_k      <<- PROFILE$masigpro_k
}

#' 修正prefix问题：在data_raw下创建00_setup.R期望的文件别名
#' 
#' 根因：00_setup.R用basename(PROJECT_DIR)作prefix → "GEO_GSE307424_Lung"
#'        但S03b创建的文件用"Lung"前缀
#'        而且每个pipeline脚本开头都re-source 00_setup.R，覆盖无效
#' 
#' 解决：在data_raw下创建从expected名到actual名的副本（file.copy）
#'        这样无论00_setup.R被source多少次，FILES路径都能解析
#'        只需调用一次（run_phase1开头），之后文件持久存在
fix_prefix <- function() {
  expected_prefix <- basename(PROFILE$project_dir)  # "GEO_GSE307424_Lung"
  actual_prefix   <- PROFILE$prefix                  # "Lung"
  
  if (expected_prefix == actual_prefix) {
    cat("  [PREFIX] No alias needed, prefixes match\n")
    return(invisible(TRUE))
  }
  
  raw_dir <- file.path(PROFILE$project_dir, "data_raw")
  
  # 需要创建别名的4个文件
  suffixes <- c("_all_counts_with_order.tsv",
                "_all_tpm.tsv",
                "_metadata.tsv",
                "_gene_annotation.tsv")
  
  for (sfx in suffixes) {
    actual_file   <- file.path(raw_dir, paste0(actual_prefix, sfx))
    expected_file <- file.path(raw_dir, paste0(expected_prefix, sfx))
    
    if (file.exists(expected_file)) {
      # 别名已存在，跳过
      next
    }
    
    if (!file.exists(actual_file)) {
      cat(sprintf("  [WARN] Source file missing: %s\n", basename(actual_file)))
      next
    }
    
    # 创建副本（Windows不可靠地支持symlink，用file.copy更安全）
    file.copy(actual_file, expected_file, overwrite = FALSE)
    cat(sprintf("  [ALIAS] %s → %s\n", basename(actual_file), basename(expected_file)))
  }
  
  # 同步更新lineage显示变量（这些不会被re-source覆盖，因为00_setup用的是diff_type）
  # 但diff_type本身会被每次re-source重置，所以我们接受它为"GEO_GSE307424_Lung"
  # 只要FILES路径能找到文件就行
  
  cat(sprintf("  [PREFIX] File aliases created: '%s_*' → '%s_*'\n",
              actual_prefix, expected_prefix))
  invisible(TRUE)
}

#' 获取当前数据集的时间点标签向量（替代硬编码的c("4d","7d","14d","21d")）
get_timepoints <- function() {
  PROFILE$time_labels
}

#' 获取当前数据集的时间数值向量
get_time_values <- function() {
  PROFILE$time_values
}

#' 获取当前数据集的Induced组目标groups
get_target_groups <- function() {
  paste0("Induced_", PROFILE$time_labels)
}

#' 确保全局环境变量已设置（支持单独调用任意phase）
ensure_env <- function() {
  options(error = NULL)  # 防止Browse调试模式
  PROJECT_DIR <<- PROFILE$project_dir
  SCRIPT_DIR  <<- PIPELINE_SCRIPTS
  source(file.path(SCRIPT_DIR, "00_setup.R"), local = FALSE)
  fix_params()
  fix_prefix()
  fix_colors()
}

#' 自动将 geo_datasets/{GSE}/data_raw/ 下的文件复制到 pipeline_runs/{project}/data_raw/
#' 如果 pipeline_runs 下已有文件则跳过。
#' 
#' S03b 将数据准备到 geo_datasets/ 下，但 pipeline (00_setup.R) 读取 pipeline_runs/ 下的文件。
#' 本函数弥合这个路径差异，对所有数据集通用——新增数据集时无需手动复制。
auto_stage_data <- function() {
  dst_dir <- file.path(PROFILE$project_dir, "data_raw")
  
  # 如果目标目录已有文件，跳过
  if (dir.exists(dst_dir) && length(list.files(dst_dir)) > 0) {
    cat(sprintf("  [STAGE] data_raw already populated (%d files), skipping\n",
                length(list.files(dst_dir))))
    return(invisible(TRUE))
  }
  
  # 推断 geo_datasets 源路径:
  #   pipeline_runs 和 geo_datasets 都在 METHODS_BASE 下
  #   PROFILE$project_dir = file.path(RUN_DIR, "GEO_XXX_YYY")
  #   → METHODS_BASE = dirname(dirname(PROFILE$project_dir))
  #   → geo source  = METHODS_BASE/geo_datasets/DATASET_ID/data_raw
  methods_base <- dirname(dirname(PROFILE$project_dir))  # 上两级 = methods/
  src_dir <- file.path(methods_base, "geo_datasets", DATASET_ID, "data_raw")
  
  if (!dir.exists(src_dir)) {
    cat(sprintf("  [STAGE] Source not found: %s\n", src_dir))
    cat("  [STAGE] Please run S03b first to prepare the dataset.\n")
    stop("Cannot find prepared data for ", DATASET_ID)
  }
  
  src_files <- list.files(src_dir, full.names = TRUE)
  if (length(src_files) == 0) {
    stop(sprintf("Source directory is empty: %s", src_dir))
  }
  
  dir.create(dst_dir, recursive = TRUE, showWarnings = FALSE)
  file.copy(src_files, dst_dir, overwrite = FALSE)
  
  cat(sprintf("  [STAGE] Copied %d files: %s -> %s\n",
              length(src_files), src_dir, dst_dir))
  invisible(TRUE)
}

#' 统一样本名时间后缀为 "d"（pipeline约定格式）
#' 
#' 01_data_import.R 用硬编码regex解析 {prefix}{number}d{rep} 格式。
#' S03b 对小时级数据集生成的样本名用 "h" 后缀 (如 Chon3h1, Chon18hC2)，
#' 导致解析失败。本函数在 01_data_import 之前自动检测并修正。
#' 
#' 只修改 data_raw 下的 counts/tpm/metadata 文件中的列名/Sample列。
#' 已经用 "d" 后缀的数据集（GSE236646, GSE307424等）自动跳过。
normalize_sample_suffix <- function() {
  raw_dir   <- file.path(PROFILE$project_dir, "data_raw")
  actual_pfx <- PROFILE$prefix
  
  # 读取任一counts文件的表头以检测后缀
  for (pfx in c(basename(PROFILE$project_dir), actual_pfx)) {
    cf <- file.path(raw_dir, paste0(pfx, "_all_counts_with_order.tsv"))
    if (file.exists(cf)) { counts_file <- cf; break }
  }
  if (!exists("counts_file") || !file.exists(counts_file)) {
    cat("  [SUFFIX] No counts file found, skipping suffix check\n")
    return(invisible(FALSE))
  }
  
  header <- names(read.delim(counts_file, nrows = 1, check.names = FALSE))
  sample_cols <- setdiff(header, "gene_id")
  
  # 检测: {prefix}{digits}{单字母} 中的那个字母
  pat <- paste0("^", actual_pfx, "\\d+([a-zA-Z])")
  matches <- regmatches(sample_cols, regexec(pat, sample_cols))
  detected <- unique(sapply(matches, function(x) if (length(x) >= 2) x[2] else NA))
  detected <- detected[!is.na(detected)]
  
  if (length(detected) == 0 || all(detected == "d")) {
    cat("  [SUFFIX] Sample names already use 'd' suffix, OK\n")
    return(invisible(TRUE))
  }
  
  old_sfx <- detected[1]
  cat(sprintf("  [SUFFIX] Detected suffix '%s' in sample names, normalizing to 'd'...\n", old_sfx))
  
  # gsub 替换模式: {prefix}{digits}{old_suffix} → {prefix}{digits}d
  gsub_pat <- paste0("^(", actual_pfx, "\\d+)", old_sfx)
  gsub_rep <- "\\1d"
  
  # 处理所有 data_raw 下的文件
  file_suffixes <- c("_all_counts_with_order.tsv", "_all_tpm.tsv", "_metadata.tsv")
  prefixes_to_try <- unique(c(basename(PROFILE$project_dir), actual_pfx))
  n_renamed <- 0
  
  for (fs in file_suffixes) {
    for (pfx in prefixes_to_try) {
      f <- file.path(raw_dir, paste0(pfx, fs))
      if (!file.exists(f)) next
      
      dat <- read.delim(f, check.names = FALSE, stringsAsFactors = FALSE)
      
      if (fs == "_metadata.tsv") {
        # metadata: 修改 Sample 列
        if ("Sample" %in% names(dat)) {
          dat$Sample <- gsub(gsub_pat, gsub_rep, dat$Sample)
        }
      } else {
        # counts/tpm: 修改列名
        names(dat) <- gsub(gsub_pat, gsub_rep, names(dat))
      }
      
      write.table(dat, f, sep = "\t", quote = FALSE, row.names = FALSE)
      n_renamed <- n_renamed + 1
    }
  }
  
  cat(sprintf("  [SUFFIX] Normalized %d files: '%s' -> 'd'\n", n_renamed, old_sfx))
  invisible(TRUE)
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 1: 00_setup.R + 01_data_import.R + 02 + 03                         ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_phase1 <- function() {
  cat("\n### PHASE 1: Data import + preprocessing + normalization ###\n\n")
  
  # ---- 0. 自动复制 geo_datasets → pipeline_runs ----
  auto_stage_data()
  
  # ---- 1a. source 00_setup.R ----
  PROJECT_DIR <<- PROFILE$project_dir
  SCRIPT_DIR  <<- PIPELINE_SCRIPTS
  source(file.path(SCRIPT_DIR, "00_setup.R"), local = FALSE)
  
  # 立即覆盖PARAMS和COLORS
  fix_params()
  fix_prefix()
  fix_colors()
  
  cat(sprintf("  [OVERRIDE] masigpro_degree = %d (was 3)\n", PARAMS$masigpro_degree))
  cat(sprintf("  [OVERRIDE] masigpro_k = %d (was 9)\n", PARAMS$masigpro_k))
  cat(sprintf("  [OVERRIDE] COLORS$timepoint: %s\n",
              paste(names(COLORS$timepoint), collapse = ", ")))
  
  # ---- 1a-bis. 统一样本名后缀为 "d" ----
  # 01_data_import.R 硬编码解析 "{prefix}{number}d{rep}" 格式
  # S03b 生成的小时级数据集用 "h" 后缀 (Chon3h1) → 解析失败
  # 本步骤在 01_data_import 之前自动修正
  normalize_sample_suffix()
  
  # ---- 1b. source 01_data_import.R ----
  # 这个脚本会解析样本名并设置Time/Group factor levels
  # 问题：L77硬编码 levels = c("4d","7d","14d","21d")
  # 策略：source后立即修正
  source(file.path(SCRIPT_DIR, "01_data_import.R"), local = FALSE)
  
  # ---- 排除不完整时间点的样本（如GSE197067的0h只有Control） ----
  if (!is.null(PROFILE$exclude_samples_pattern)) {
    excl_pat <- PROFILE$exclude_samples_pattern
    excl_idx <- grepl(excl_pat, colnames(counts_raw))
    n_excl <- sum(excl_idx)
    if (n_excl > 0) {
      excl_names <- colnames(counts_raw)[excl_idx]
      cat(sprintf("\n  [EXCLUDE] Removing %d samples matching '%s': %s\n",
                  n_excl, excl_pat, paste(excl_names, collapse = ", ")))
      
      # 从counts和tpm矩阵中删除
      counts_raw <<- counts_raw[, !excl_idx, drop = FALSE]
      tpm_raw    <<- tpm_raw[, !excl_idx, drop = FALSE]
      
      # 从sample_info中删除
      sample_info <<- sample_info[!excl_idx, , drop = FALSE]
      
      cat(sprintf("  [EXCLUDE] Remaining: %d samples\n", ncol(counts_raw)))
    }
  }
  
  # 修正Time和Group因子水平
  sample_info <<- fix_sample_info_factors(sample_info)
  
  # 验证修正结果
  cat("\n  [FIX] sample_info$Time levels after correction:\n")
  print(table(sample_info$Time))
  cat("  [FIX] sample_info$Group levels after correction:\n")
  print(table(sample_info$Group))
  
  # 检查是否有NA
  if (any(is.na(sample_info$Time))) {
    stop("FATAL: sample_info$Time contains NA after fix! Check time_labels in PROFILE.")
  }
  
  # 重新保存修正后的sample_info
  save_data(sample_info, FILES$sample_info)
  
  # 如果有样本被排除，必须把修剪后的counts/tpm写回磁盘
  # 否则02_preprocessing从FILES$counts_raw读到的还是原始列数
  if (!is.null(PROFILE$exclude_samples_pattern)) {
    cat("  [WRITE] Overwriting raw data files with excluded samples removed...\n")
    counts_df <- data.frame(gene_id = rownames(counts_raw), counts_raw, check.names = FALSE)
    write.table(counts_df, FILES$counts_raw, sep = "\t", row.names = FALSE, quote = FALSE)
    tpm_df <- data.frame(gene_id = rownames(tpm_raw), tpm_raw, check.names = FALSE)
    write.table(tpm_df, FILES$tpm_raw, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  [WRITE] %s: %d genes × %d samples\n",
                basename(FILES$counts_raw), nrow(counts_raw), ncol(counts_raw)))
  }
  
  # ---- 1c. 基因ID适配 ----
  # 如果数据集用gene symbol（非Ensembl），需要修正gene_annotation
  if (PROFILE$gene_id_type == "symbol") {
    cat("\n  [ADAPT] Gene ID type = symbol (not Ensembl)\n")
    # gene_anno已由01_data_import.R创建
    # 但它试图用ensembl_gene_id去查org.Hs.eg.db，symbol当作Ensembl会失败
    # 修正：用SYMBOL作为keytype重新映射
    
    gene_anno_fixed <- data.frame(
      ensembl_gene_id = rownames(counts_raw),  # 实际是gene symbol
      hgnc_symbol     = rownames(counts_raw),  # 同样是gene symbol
      stringsAsFactors = FALSE
    )
    
    # 尝试获取Entrez ID（用SYMBOL keytype）
    tryCatch({
      sym2entrez <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = unique(gene_anno_fixed$hgnc_symbol),
        columns = c("ENTREZID"),
        keytype = "SYMBOL"
      )
      sym2entrez <- sym2entrez[!duplicated(sym2entrez$SYMBOL), ]
      gene_anno_fixed$entrez_id <- sym2entrez$ENTREZID[
        match(gene_anno_fixed$hgnc_symbol, sym2entrez$SYMBOL)]
      
      n_mapped <- sum(!is.na(gene_anno_fixed$entrez_id))
      cat(sprintf("  [ADAPT] Symbol->Entrez mapped: %d/%d (%.1f%%)\n",
                  n_mapped, nrow(gene_anno_fixed),
                  100 * n_mapped / nrow(gene_anno_fixed)))
    }, error = function(e) {
      cat(sprintf("  [WARN] Entrez mapping failed: %s\n", e$message))
      gene_anno_fixed$entrez_id <<- NA
    })
    
    gene_anno <<- gene_anno_fixed
    save_data(gene_anno, FILES$gene_annotation)
    cat("  [ADAPT] gene_annotation overwritten with symbol-based mapping\n")
  }
  
  # ---- 1c-bis. Entrez → Ensembl ID mapping (NCBI counts) ----
  else if (PROFILE$gene_id_type == "entrez") {
    cat("\n  [ADAPT] Gene ID type = entrez (NCBI counts)\n")
    cat("  [ADAPT] Performing Entrez -> Ensembl ID mapping...\n")
    
    entrez_ids <- rownames(counts_raw)
    n_total <- length(entrez_ids)
    cat(sprintf("  [ADAPT] Input: %d genes with Entrez IDs\n", n_total))
    
    tryCatch({
      mapping <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys    = entrez_ids,
        keytype = "ENTREZID",
        columns = c("ENSEMBL", "SYMBOL")
      )
      # 每个Entrez ID只保留第一个Ensembl映射
      mapping <- mapping[!duplicated(mapping$ENTREZID), ]
      
      n_mapped   <- sum(!is.na(mapping$ENSEMBL))
      n_unmapped <- sum(is.na(mapping$ENSEMBL))
      cat(sprintf("  [ADAPT] Mapped: %d/%d (%.1f%%), Unmapped: %d\n",
                  n_mapped, n_total, 100 * n_mapped / n_total, n_unmapped))
      
      # 只保留有Ensembl映射的基因
      mapping_valid <- mapping[!is.na(mapping$ENSEMBL), ]
      
      # 处理反向重复: 多个Entrez -> 同一Ensembl
      ens_dup <- duplicated(mapping_valid$ENSEMBL)
      if (any(ens_dup)) {
        cat(sprintf("  [ADAPT] %d duplicate Ensembl IDs removed\n", sum(ens_dup)))
        mapping_valid <- mapping_valid[!ens_dup, ]
      }
      
      keep_entrez <- mapping_valid$ENTREZID
      new_ensembl <- mapping_valid$ENSEMBL
      new_symbol  <- mapping_valid$SYMBOL
      
      # 更新 counts_raw 和 tpm_raw 的行名
      counts_raw <<- counts_raw[keep_entrez, , drop = FALSE]
      rownames(counts_raw) <<- new_ensembl
      
      tpm_raw <<- tpm_raw[keep_entrez, , drop = FALSE]
      rownames(tpm_raw) <<- new_ensembl
      
      cat(sprintf("  [ADAPT] counts_raw: %d -> %d genes (after mapping)\n",
                  n_total, nrow(counts_raw)))
      
      # 重写磁盘上的 counts/tpm 文件
      cat("  [ADAPT] Overwriting counts/tpm files with Ensembl IDs...\n")
      counts_df <- data.frame(gene_id = rownames(counts_raw), counts_raw, check.names = FALSE)
      write.table(counts_df, FILES$counts_raw, sep = "\t", row.names = FALSE, quote = FALSE)
      tpm_df <- data.frame(gene_id = rownames(tpm_raw), tpm_raw, check.names = FALSE)
      write.table(tpm_df, FILES$tpm_raw, sep = "\t", row.names = FALSE, quote = FALSE)
      
      # 构建 gene_annotation
      gene_anno_fixed <- data.frame(
        ensembl_gene_id = new_ensembl,
        hgnc_symbol     = ifelse(is.na(new_symbol), new_ensembl, new_symbol),
        entrez_id       = keep_entrez,
        stringsAsFactors = FALSE
      )
      
      gene_anno <<- gene_anno_fixed
      save_data(gene_anno, FILES$gene_annotation)
      
      # 更新 data_raw 下的 annotation 文件 (fix_prefix 别名)
      for (pfx in c(basename(PROFILE$project_dir), PROFILE$prefix)) {
        anno_path <- file.path(PROFILE$project_dir, "data_raw",
                                paste0(pfx, "_gene_annotation.tsv"))
        if (file.exists(anno_path)) {
          write.table(gene_anno_fixed, anno_path,
                      sep = "\t", row.names = FALSE, quote = FALSE)
        }
      }
      
      cat(sprintf("  [ADAPT] Gene annotation: %d genes (Ensembl + Symbol + Entrez)\n",
                  nrow(gene_anno_fixed)))
      cat("  [ADAPT] Entrez -> Ensembl mapping COMPLETE\n")
      
    }, error = function(e) {
      stop(sprintf("[FATAL] Entrez -> Ensembl mapping failed: %s\n", e$message))
    })
  }
  
  # ---- 1d. source 02_preprocessing.R ----
  source(file.path(SCRIPT_DIR, "02_preprocessing.R"), local = FALSE)
  
  # ---- 1e. 03_normalization_QC (ADAPTED) ----
  # 问题：L38 relevel(ref="4d") + L81/L110 shape硬编码
  # 不能source原脚本（relevel在source中途执行会报错），写adapted版
  cat("  [ADAPT] Running adapted normalization + QC (03)...\n")
  run_03_normalization_adapted()
  
  cat("\n### PHASE 1 COMPLETE ###\n")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  ADAPTED SCRIPT: 03_normalization_QC                                       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_03_normalization_adapted <- function() {
  ensure_env()
  
  log_step("03_NORM", "Starting normalization and QC (ADAPTED)...")
  
  counts_filtered <- readRDS(FILES$counts_filtered)
  sample_info     <- readRDS(FILES$sample_info)
  sample_info     <- fix_sample_info_factors(sample_info)
  
  # 1. DESeqDataSet
  counts_mat <- as.matrix(counts_filtered)
  storage.mode(counts_mat) <- "integer"
  stopifnot(all(colnames(counts_mat) == rownames(sample_info)))
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = sample_info,
    design    = ~ Treatment + Time + Treatment:Time
  )
  dds$Treatment <- relevel(dds$Treatment, ref = "Control")
  dds$Time      <- relevel(dds$Time, ref = PROFILE$time_ref)  # ★ ADAPTED ★
  
  log_step("03_NORM", sprintf("DESeqDataSet: %d genes × %d samples, Time ref = %s",
                               nrow(dds), ncol(dds), levels(dds$Time)[1]))
  
  # 2. VST
  vsd <- vst(dds, blind = FALSE)
  vst_mat <- assay(vsd)
  log_step("03_NORM", "VST transformation complete")
  
  save_data(dds, FILES$dds_object)
  save_data(vst_mat, FILES$vst_matrix)
  
  # 3. PCA
  log_step("03_NORM", "PCA analysis...")
  ntop <- 500
  rv <- matrixStats::rowVars(vst_mat)
  select_genes <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca_data <- prcomp(t(vst_mat[select_genes, ]), center = TRUE, scale. = FALSE)
  
  pca_df <- as.data.frame(pca_data$x[, 1:min(5, ncol(pca_data$x))])
  pca_df$sample    <- rownames(pca_df)
  pca_df$Treatment <- sample_info$Treatment
  pca_df$Time      <- sample_info$Time
  pca_df$Group     <- sample_info$Group
  pct_var <- round(100 * (pca_data$sdev^2 / sum(pca_data$sdev^2)), 1)
  
  # Shape mapping (ADAPTED)
  shape_vals <- PROFILE$shape_timepoint
  
  p_pca_main <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Treatment, shape = Time), size = 3.5, stroke = 0.8) +
    scale_color_treatment() +
    scale_shape_manual(values = shape_vals) +  # ★ ADAPTED ★
    stat_ellipse(aes(group = Treatment, color = Treatment),
                 type = "t", level = 0.95, linetype = 2, linewidth = 0.5) +
    labs(title = "PCA of Gene Expression",
         subtitle = sprintf("Top %d variable genes, VST-transformed", ntop),
         x = sprintf("PC1 (%s%% variance)", pct_var[1]),
         y = sprintf("PC2 (%s%% variance)", pct_var[2])) +
    theme_bindlab_minimal()
  save_pub_fig(p_pca_main, "PCA_Treatment_Time", "02_PCA_Clustering", width = 8, height = 6)
  
  p_pca_group <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Group), size = 3.5) +
    scale_color_manual(values = COLORS$group) +
    ggrepel::geom_text_repel(aes(label = sample), size = 2.2, max.overlaps = 20) +
    labs(title = "PCA — All Groups",
         x = sprintf("PC1 (%s%% variance)", pct_var[1]),
         y = sprintf("PC2 (%s%% variance)", pct_var[2])) +
    theme_bindlab_minimal()
  save_pub_fig(p_pca_group, "PCA_AllGroups_labeled", "02_PCA_Clustering", width = 10, height = 7)
  
  # Scree plot
  scree_df <- data.frame(PC = paste0("PC", 1:min(10, length(pct_var))),
                         Variance = pct_var[1:min(10, length(pct_var))])
  scree_df$PC <- factor(scree_df$PC, levels = scree_df$PC)
  p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "#3C5488", width = 0.6) +
    geom_line(aes(group = 1), color = "#E64B35", linewidth = 0.8) +
    geom_point(color = "#E64B35", size = 2) +
    geom_text(aes(label = paste0(Variance, "%")), vjust = -0.5, size = 3) +
    labs(title = "PCA Scree Plot", x = NULL, y = "Variance Explained (%)") +
    theme_bindlab() +
    expand_limits(y = max(scree_df$Variance) * 1.15)
  save_pub_fig(p_scree, "PCA_scree_plot", "02_PCA_Clustering", width = 8, height = 5)
  
  # 4. Sample correlation heatmap
  log_step("03_NORM", "Sample correlation heatmap...")
  cor_mat <- cor(vst_mat, method = "spearman")
  anno_col <- data.frame(
    Treatment = sample_info$Treatment,
    Time      = sample_info$Time,
    row.names = rownames(sample_info)
  )
  anno_colors <- list(
    Treatment = COLORS$treatment,
    Time      = COLORS$timepoint
  )
  
  save_heatmap_fig(
    draw_func = function() {
      print(pheatmap::pheatmap(cor_mat,
               clustering_distance_rows = as.dist(1 - cor_mat),
               clustering_distance_cols = as.dist(1 - cor_mat),
               clustering_method = "complete",
               color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
               breaks = seq(min(cor_mat), 1, length.out = 101),
               annotation_col = anno_col,
               annotation_colors = anno_colors,
               show_rownames = TRUE, show_colnames = TRUE,
               fontsize = 8,
               main = "Sample Correlation (Spearman)"))
    },
    filename = "Correlation_heatmap_spearman",
    subdir = "02_PCA_Clustering", width = 10, height = 9
  )
  
  # 5. Hierarchical clustering
  log_step("03_NORM", "Hierarchical clustering...")
  dist_mat <- as.dist(1 - cor_mat)
  hc <- hclust(dist_mat, method = "complete")
  dend <- as.dendrogram(hc)
  labels_order <- labels(dend)
  label_colors <- ifelse(sample_info[labels_order, "Treatment"] == "Induced",
                         COLORS$treatment["Induced"],
                         COLORS$treatment["Control"])
  dend <- dend %>%
    dendextend::set("labels_cex", 0.7) %>%
    dendextend::set("labels_col", label_colors)
  
  save_heatmap_fig(
    draw_func = function() {
      par(mar = c(8, 4, 3, 1))
      plot(dend, main = "Sample Hierarchical Clustering",
           ylab = "1 - Spearman Correlation", xlab = "")
    },
    filename = "Sample_clustering_dendrogram",
    subdir = "02_PCA_Clustering", width = 12, height = 6
  )
  
  # 6. Distance heatmap
  sample_dist <- dist(t(vst_mat))
  sample_dist_mat <- as.matrix(sample_dist)
  save_heatmap_fig(
    draw_func = function() {
      print(pheatmap::pheatmap(sample_dist_mat,
               clustering_distance_rows = sample_dist,
               clustering_distance_cols = sample_dist,
               color = colorRampPalette(c("#00A087", "white", "#E64B35"))(100),
               annotation_col = anno_col,
               annotation_colors = anno_colors,
               show_rownames = TRUE, show_colnames = FALSE,
               fontsize = 7,
               main = "Sample Distance (Euclidean on VST)"))
    },
    filename = "Sample_distance_heatmap",
    subdir = "02_PCA_Clustering", width = 10, height = 9
  )
  
  log_step("03_NORM", "Step 03 COMPLETE (ADAPTED)")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  ADAPTED SCRIPT: 06_maSigPro_trends                                       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_06_maSigPro_adapted <- function() {
  ensure_env()
  
  log_step("06_MASIGPRO", sprintf(
    "Starting maSigPro (ADAPTED: degree=%d, k=%d)...",
    PARAMS$masigpro_degree, PARAMS$masigpro_k))
  
  counts_filtered <- readRDS(FILES$counts_filtered)
  sample_info     <- readRDS(FILES$sample_info)
  sample_info     <- fix_sample_info_factors(sample_info)
  gene_anno       <- readRDS(FILES$gene_annotation)
  
  # 1. edesign
  edesign <- data.frame(
    Time      = sample_info$time_num,
    Replicate = sample_info$replicate,
    Control   = as.integer(sample_info$Treatment == "Control"),
    Induced   = as.integer(sample_info$Treatment == "Induced"),
    row.names = rownames(sample_info)
  )
  edesign <- edesign[colnames(counts_filtered), ]
  
  log_step("06_MASIGPRO", sprintf("Design: %d samples, Time: %s",
    nrow(edesign), paste(sort(unique(edesign$Time)), collapse = ", ")))
  
  # 2. maSigPro analysis — degree ADAPTED
  data_mat <- as.matrix(counts_filtered)
  design_matrix <- make.design.matrix(edesign, degree = PARAMS$masigpro_degree)
  
  min_obs_use <- PROFILE$masigpro_min_obs
  log_step("06_MASIGPRO", sprintf("Step 1: p.vector (min.obs=%d)...", min_obs_use))
  fit <- p.vector(data_mat, design_matrix,
                  Q = PARAMS$masigpro_alfa,
                  counts = TRUE,
                  min.obs = min_obs_use)
  log_step("06_MASIGPRO", sprintf("Step 1: %d significant genes", fit$i))
  
  log_step("06_MASIGPRO", "Step 2: T.fit (stepwise regression)...")
  tstep <- T.fit(fit, step.method = "backward", alfa = PARAMS$masigpro_alfa)
  
  log_step("06_MASIGPRO", sprintf("Step 3: Filtering R² >= %.2f...", PARAMS$masigpro_rsq))
  sigs <- get.siggenes(tstep, rsq = PARAMS$masigpro_rsq, vars = "groups")
  sig_genes_all <- sigs$summary
  
  for (grp in names(sig_genes_all)) {
    log_step("06_MASIGPRO", sprintf("  %s: %d genes", grp, length(sig_genes_all[[grp]])))
  }
  
  # 3. Clustering
  cluster_result <- NULL
  tryCatch({
    sig_genes_for_cluster <- NULL
    if (!is.null(sigs$sig.genes)) {
      avail_keys <- names(sigs$sig.genes)
      log_step("06_MASIGPRO", sprintf("sig.genes keys: %s", paste(avail_keys, collapse = ", ")))
      if ("InducedvsControl" %in% avail_keys) {
        sig_genes_for_cluster <- sigs$sig.genes$InducedvsControl
      } else if (length(avail_keys) > 0) {
        sig_genes_for_cluster <- sigs$sig.genes[[avail_keys[1]]]
      }
    }
    
    if (!is.null(sig_genes_for_cluster)) {
      # 动态调整k：如果基因太少，减小k
      n_sig <- if (is.list(sig_genes_for_cluster)) 
        nrow(sig_genes_for_cluster$sig.profiles) else nrow(sig_genes_for_cluster)
      k_use <- min(PARAMS$masigpro_k, max(2, floor(n_sig / 5)))
      log_step("06_MASIGPRO", sprintf("Clustering %d genes into %d groups...", n_sig, k_use))
      
      pdf(file.path(FIG_DIR, "06_maSigPro", "maSigPro_cluster_profiles.pdf"),
          width = 14, height = 10)
      cluster_result <- see.genes(sig_genes_for_cluster,
                                  edesign = edesign,
                                  groups.vector = design_matrix$groups.vector,
                                  show.fit = TRUE,
                                  dis = design_matrix$dis,
                                  cluster.method = "hclust",
                                  cluster.data = 1,
                                  k = k_use,
                                  newX11 = FALSE)
      dev.off()
      log_step("06_MASIGPRO", "Cluster profiles plotted")
    } else {
      log_step("06_MASIGPRO", "WARNING: No sig genes for clustering")
    }
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    log_step("06_MASIGPRO", sprintf("WARNING: Clustering error: %s", e$message))
  })
  
  # 4. Extract clusters
  gene_clusters <- NULL
  if (!is.null(cluster_result)) {
    gene_clusters <- data.frame(
      ensembl_id = names(cluster_result$cut),
      cluster    = cluster_result$cut,
      stringsAsFactors = FALSE
    )
    gene_clusters$symbol <- gene_anno$hgnc_symbol[
      match(gene_clusters$ensembl_id, gene_anno$ensembl_gene_id)]
    log_step("06_MASIGPRO", "Genes per cluster:")
    print(table(gene_clusters$cluster))
    write.csv(gene_clusters, file.path(DATA_DIR, "maSigPro_gene_clusters.csv"),
              row.names = FALSE)
  }
  
  # 5. Custom trend plots (ADAPTED: use actual time values for axis)
  if (!is.null(cluster_result)) {
    tpm_filtered <- readRDS(FILES$tpm_filtered)
    time_vals <- sort(unique(sample_info$time_num))
    time_labs <- PROFILE$time_labels
    
    for (cl in sort(unique(gene_clusters$cluster))) {
      cl_genes <- gene_clusters$ensembl_id[gene_clusters$cluster == cl]
      if (length(cl_genes) < 2) next
      
      tpm_cl_log <- log2(tpm_filtered[cl_genes, , drop = FALSE] + 1)
      
      plot_data <- tpm_cl_log %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
        dplyr::left_join(
          sample_info %>% tibble::rownames_to_column("sample") %>%
            dplyr::select(sample, Treatment, Time, time_num),
          by = "sample") %>%
        dplyr::group_by(Treatment, Time, time_num) %>%
        dplyr::summarise(
          mean_expr = mean(expr), se_expr = sd(expr) / sqrt(dplyr::n()),
          .groups = "drop")
      
      p_trend <- ggplot(plot_data, aes(x = time_num, y = mean_expr,
                                       color = Treatment, group = Treatment)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2.5) +
        geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr),
                      width = max(time_vals) * 0.03, linewidth = 0.5) +
        scale_color_treatment() +
        scale_x_continuous(breaks = time_vals, labels = time_labs) +  # ★ ADAPTED ★
        labs(title = sprintf("Cluster %d (%d genes)", cl, length(cl_genes)),
             x = "Time", y = "Mean log2(TPM + 1)") +
        theme_bindlab()
      
      save_pub_fig(p_trend, sprintf("maSigPro_Cluster%d_trend", cl),
                   "06_maSigPro", width = 6, height = 4.5)
    }
  }
  
  # 6. Save
  masigpro_results <- list(
    fit            = fit,
    tstep          = tstep,
    sigs           = sigs,
    sig_genes_all  = sig_genes_all,
    cluster_result = cluster_result,
    gene_clusters  = gene_clusters,
    design_matrix  = design_matrix
  )
  save_data(masigpro_results, FILES$masigpro_results)
  
  log_step("06_MASIGPRO", "Step 06 COMPLETE (ADAPTED)")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 2: DEG + maSigPro + WGCNA                                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_phase2 <- function() {
  cat("\n### PHASE 2: DEG analysis + time series + WGCNA ###\n\n")
  ensure_env()
  
  # ---- 2a. 04_DEG_analysis.R ----
  # 问题：L65 timepoints <- c("4d","7d","14d","21d")
  #       L137 target_groups <- c("Induced_4d", ...)
  # 策略：source后让timepoints和target_groups指向正确值
  #       但04内部用了局部变量，source时已经执行完了
  #       关键问题：L65的timepoints决定了Wald检验做哪些对比
  #       如果Group factor已修正，c("group","Induced_6h","Control_6h")是合法contrast
  #       但04_DEG硬编码了for loop里的timepoints!
  #
  # 最安全的方案：临时覆盖04_DEG的核心逻辑
  # 但由于"不改原脚本"原则，我们需要另一个策略：
  #   → 04_DEG.R 在执行时，timepoints变量在其内部定义
  #   → 我们无法在source前注入局部变量
  #   → 需要写一个替代版本的04_DEG for GEO数据

  cat("  [ADAPT] Running adapted DEG analysis (04_DEG)...\n")
  run_04_DEG_adapted()
  
  # ---- 2b. 05_DEG_visualization.R ----
  # 同理，hardcoded timepoints
  cat("  [ADAPT] Running adapted DEG visualization (05)...\n")
  run_05_DEGvis_adapted()
  
  # ---- 2c. 06_maSigPro_trends.R ----
  # 06 re-source 00_setup.R → degree被重置为3 → 3时间点+degree=3溢出
  # 需要adapted版本
  cat(sprintf("  [ADAPT] Running adapted maSigPro (degree=%d, k=%d)...\n",
      PROFILE$masigpro_degree, PROFILE$masigpro_k))
  run_06_maSigPro_adapted()
  
  # ---- 2d. 08_WGCNA.R ----
  # WGCNA从VST矩阵和sample_info读取，没有硬编码时间
  # group_levels在L67-68硬编码了8组，但它用的是Group列的实际值
  # 等一下——L67实际上是：
  #   group_levels <- c("Control_4d", ... "Induced_21d")  ← 硬编码！
  # 需要处理
  cat("  [ADAPT] Running WGCNA...\n")
  run_08_WGCNA_adapted()
  
  cat("\n### PHASE 2 COMPLETE ###\n")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 3: 候选筛选（09A → 09C → 09D → 09F → 10）                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_phase3 <- function() {
  cat("\n### PHASE 3: Candidate selection pipeline ###\n\n")
  ensure_env()
  
  # 预加载全局变量（原pipeline在RUN_GUIDE里从01一路source下来，这些变量自然存在）
  # 用assign确保放在.GlobalEnv中（<<-在嵌套环境中可能指向错误位置）
  assign("sample_info",  fix_sample_info_factors(readRDS(FILES$sample_info)), envir = .GlobalEnv)
  assign("gene_anno",    readRDS(FILES$gene_annotation), envir = .GlobalEnv)
  assign("tpm_filtered", readRDS(FILES$tpm_filtered), envir = .GlobalEnv)
  assign("all_results",  readRDS(FILES$deg_results), envir = .GlobalEnv)
  
  # 09A: candidate pool
  source(file.path(SCRIPT_DIR, "09A_candidate_pool.R"), local = FALSE)
  
  # 重新注入sample_info（确保在.GlobalEnv中）
  assign("sample_info", fix_sample_info_factors(readRDS(FILES$sample_info)), envir = .GlobalEnv)
  
  # 09C: bootstrap stability selection
  source(file.path(SCRIPT_DIR, "09C_ML_stability_selection.R"), local = FALSE)
  
  # 09D: gap-union selection
  # [v2] 不预排除任何算法 — RF通过 MIN_FREQ_SIGNAL=0.20 自动排除
  # 文献: Strobl et al. 2007 BMC Bioinf; Nicodemus et al. 2010 BMC Bioinf
  #   RF在p>>n场景下VIM不稳定，但应由数据驱动排除，而非预设
  {
    lines_09D <- readLines(file.path(SCRIPT_DIR, "09D_gap_union_selection.R"))
    idx_exclude <- grep('EXCLUDE_ALGOS\\s*<-\\s*c\\("RF"\\)', lines_09D)
    if (length(idx_exclude) == 1) {
      lines_09D[idx_exclude] <- 'EXCLUDE_ALGOS <- character(0)  # [v2] auto-exclude by MIN_FREQ_SIGNAL'
      cat("  [v2 PATCH] 09D: RF pre-exclusion removed → auto-exclude by MIN_FREQ_SIGNAL\n")
    }
    tmp_09D <- tempfile(pattern = "09D_v2_", fileext = ".R")
    writeLines(lines_09D, tmp_09D)
    source(tmp_09D, local = FALSE)
    unlink(tmp_09D)
  }
  
  # 09F: PPI hub — 需要STRING本地数据库
  # 检查stringdb_cache是否存在
  string_cache <- file.path(DATA_DIR, "stringdb_cache")
  links_file <- file.path(string_cache, "9606.protein.links.v12.0.txt.gz")
  
  if (file.exists(links_file)) {
    cat("  [OK] STRING local database found\n")
    source(file.path(SCRIPT_DIR, "09F_PPI_hub_selection.R"), local = FALSE)
  } else {
    cat("  [SKIP] STRING local database not found at:\n")
    cat(sprintf("    %s\n", string_cache))
    cat("  09F_PPI_hub_selection skipped. Copy STRING files to run.\n")
    cat("  Required files:\n")
    cat("    9606.protein.links.v12.0.txt.gz\n")
    cat("    9606.protein.info.v12.0.txt.gz\n")
  }
  
  # 10_integration.R — 含hardcoded timepoints on L124
  cat("  [ADAPT] Running integration (10)...\n")
  run_10_integration_adapted()
  
  cat("\n### PHASE 3 COMPLETE ###\n")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  ADAPTED SCRIPT FUNCTIONS                                                  ║
# ║  这些函数替代原pipeline中含硬编码时间引用的脚本                              ║
# ║  逻辑完全相同，只是将 c("4d","7d","14d","21d") 替换为 get_timepoints()       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ------------------------------------------------------------------
# 04_DEG_analysis adapted
# ------------------------------------------------------------------
run_04_DEG_adapted <- function() {
  ensure_env()
  
  log_step("04_DEG", "Starting DEG analysis (ADAPTED for GEO)...")
  
  dds         <- readRDS(FILES$dds_object)
  sample_info <- readRDS(FILES$sample_info)
  gene_anno   <- readRDS(FILES$gene_annotation)
  
  # --- 1a. LRT: Treatment:Time interaction ---
  log_step("04_DEG", "LRT for Treatment:Time interaction...")
  dds_lrt_interaction <- DESeq(dds, test = "LRT",
                               reduced = ~ Treatment + Time)
  res_lrt_interaction <- results(dds_lrt_interaction, alpha = PARAMS$padj_cutoff)
  n_sig_lrt <- sum(res_lrt_interaction$padj < PARAMS$padj_cutoff, na.rm = TRUE)
  log_step("04_DEG", sprintf("LRT Interaction: %d genes with padj < %.2f",
                              n_sig_lrt, PARAMS$padj_cutoff))
  
  # --- 1b. LRT: Treatment main effect ---
  log_step("04_DEG", "LRT for Treatment effect...")
  dds_lrt_treatment <- DESeq(dds, test = "LRT",
                              reduced = ~ Time)
  res_lrt_treatment <- results(dds_lrt_treatment, alpha = PARAMS$padj_cutoff)
  n_sig_trt <- sum(res_lrt_treatment$padj < PARAMS$padj_cutoff, na.rm = TRUE)
  log_step("04_DEG", sprintf("LRT Treatment: %d genes with padj < %.2f",
                              n_sig_trt, PARAMS$padj_cutoff))
  
  # --- 2. Wald test per timepoint (ADAPTED) ---
  log_step("04_DEG", "Wald tests for each timepoint...")
  dds_group <- dds
  dds_group$group <- factor(paste(dds_group$Treatment, dds_group$Time, sep = "_"))
  design(dds_group) <- ~ group
  dds_group <- DESeq(dds_group)
  
  # ★ 关键适配点：使用实际时间标签 ★
  timepoints <- get_timepoints()
  # 排除只有Control组的时间点（如GSE197067的0h）
  if (!is.null(PROFILE$has_0h_control_only) && PROFILE$has_0h_control_only) {
    timepoints <- setdiff(timepoints, "0h")
    cat("  [ADAPT] Excluding 0h from Wald contrasts (Control-only)\n")
  }
  
  deg_by_time <- list()
  deg_by_time_lfc <- list()
  
  for (tp in timepoints) {
    contrast_name <- paste0("Induced_", tp, "_vs_Control_", tp)
    contrast_vec <- c("group", paste0("Induced_", tp), paste0("Control_", tp))
    
    res <- results(dds_group, contrast = contrast_vec, alpha = PARAMS$padj_cutoff)
    
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("ensembl_id") %>%
      dplyr::left_join(gene_anno, by = c("ensembl_id" = "ensembl_gene_id")) %>%
      dplyr::arrange(padj) %>%
      dplyr::mutate(
        timepoint = tp,
        regulation = dplyr::case_when(
          padj < PARAMS$padj_cutoff & log2FoldChange > PARAMS$lfc_cutoff ~ "Up",
          padj < PARAMS$padj_cutoff & log2FoldChange < -PARAMS$lfc_cutoff ~ "Down",
          TRUE ~ "NS"
        )
      )
    deg_by_time[[tp]] <- res_df
    
    # Effect size test (lfcThreshold=1)
    res_lfc <- results(dds_group, contrast = contrast_vec,
                       alpha = PARAMS$padj_cutoff,
                       lfcThreshold = PARAMS$lfc_cutoff)
    res_lfc_df <- as.data.frame(res_lfc) %>%
      tibble::rownames_to_column("ensembl_id") %>%
      dplyr::select(ensembl_id, padj_lfc = padj)
    deg_by_time_lfc[[tp]] <- res_lfc_df
    
    n_up <- sum(res_df$regulation == "Up", na.rm = TRUE)
    n_down <- sum(res_df$regulation == "Down", na.rm = TRUE)
    log_step("04_DEG", sprintf("  %s: Up=%d, Down=%d", tp, n_up, n_down))
  }
  
  # --- 3. 8-group model: stage-specific contrasts (ADAPTED) ---
  log_step("04_DEG", "8-Group Model: Stage-Specific Planned Contrasts")
  dds_8g <- dds
  dds_8g$group <- factor(paste(dds_8g$Treatment, dds_8g$Time, sep = "_"))
  design(dds_8g) <- ~ 0 + group
  dds_8g <- DESeq(dds_8g)
  
  group_names <- levels(dds_8g$group)
  target_groups <- get_target_groups()
  # 排除只有Control的时间点
  if (!is.null(PROFILE$has_0h_control_only) && PROFILE$has_0h_control_only) {
    target_groups <- setdiff(target_groups, "Induced_0h")
  }
  
  n_groups <- length(group_names)
  stage_contrasts <- list()
  stage_deg <- list()
  
  for (target in target_groups) {
    contrast_vec <- rep(0, n_groups)
    names(contrast_vec) <- paste0("group", group_names)
    contrast_vec[paste0("group", target)] <- 1
    other_groups <- setdiff(group_names, target)
    for (og in other_groups) {
      contrast_vec[paste0("group", og)] <- -1 / (n_groups - 1)
    }
    
    res_stage <- results(dds_8g, contrast = contrast_vec, alpha = PARAMS$padj_cutoff)
    res_df <- as.data.frame(res_stage) %>%
      tibble::rownames_to_column("ensembl_id") %>%
      dplyr::left_join(gene_anno, by = c("ensembl_id" = "ensembl_gene_id")) %>%
      dplyr::arrange(padj) %>%
      dplyr::mutate(
        target_group = target,
        regulation = dplyr::case_when(
          padj < PARAMS$padj_cutoff & log2FoldChange > PARAMS$lfc_cutoff ~ "Up",
          padj < PARAMS$padj_cutoff & log2FoldChange < -PARAMS$lfc_cutoff ~ "Down",
          TRUE ~ "NS"
        )
      )
    stage_contrasts[[target]] <- contrast_vec
    stage_deg[[target]] <- res_df
    
    n_up <- sum(res_df$regulation == "Up", na.rm = TRUE)
    n_down <- sum(res_df$regulation == "Down", na.rm = TRUE)
    log_step("04_DEG", sprintf("  %s vs rest: Up=%d, Down=%d", target, n_up, n_down))
  }
  
  # Construct derived objects that downstream scripts expect
  lrt_interaction_df <- as.data.frame(res_lrt_interaction) %>%
    tibble::rownames_to_column("ensembl_id") %>%
    dplyr::left_join(gene_anno, by = c("ensembl_id" = "ensembl_gene_id")) %>%
    dplyr::arrange(padj)
  
  lrt_interaction_sig <- lrt_interaction_df %>%
    dplyr::filter(padj < PARAMS$lrt_padj)
  
  deg_all_timepoints <- dplyr::bind_rows(deg_by_time)
  
  log_step("04_DEG", sprintf("LRT interaction sig genes: %d", nrow(lrt_interaction_sig)))
  
  # Save — field names must match what 09A/09C/10 expect
  all_results <- list(
    lrt_interaction     = lrt_interaction_df,
    lrt_interaction_sig = lrt_interaction_sig,   # ★ 09C expects this ★
    lrt_treatment       = as.data.frame(res_lrt_treatment) %>%
                            tibble::rownames_to_column("ensembl_id"),
    wald_by_time        = deg_by_time,           # ★ original field name ★
    wald_by_time_lfc    = deg_by_time_lfc,       # ★ 09A expects this ★
    wald_all            = deg_all_timepoints,
    stage_deg           = stage_deg,
    stage_contrasts     = stage_contrasts,
    timepoints          = timepoints,
    dds_group           = dds_group
  )
  
  save_data(all_results, FILES$deg_results)
  
  # Save LRT separately (some scripts load it independently)
  save_data(res_lrt_interaction, FILES$lrt_results)
  
  # Save DEG summary table
  deg_all_timepoints <- dplyr::bind_rows(deg_by_time)
  write.csv(deg_all_timepoints,
            file.path(DATA_DIR, "DEG_all_timepoints.csv"),
            row.names = FALSE)
  
  log_step("04_DEG", "Step 04 COMPLETE (ADAPTED)")
}


# ------------------------------------------------------------------
# 05_DEG_visualization adapted (simplified)
# ------------------------------------------------------------------
run_05_DEGvis_adapted <- function() {
  ensure_env()
  
  log_step("05_VIS", "Starting DEG visualization (ADAPTED)...")
  
  all_results <- readRDS(FILES$deg_results)
  sample_info <- readRDS(FILES$sample_info)
  gene_anno   <- readRDS(FILES$gene_annotation)
  
  timepoints <- all_results$timepoints  # 使用04保存的实际时间点
  
  # DEG count bar chart
  deg_summary <- data.frame()
  for (tp in timepoints) {
    res <- all_results$deg_by_time[[tp]]
    if (is.null(res)) next
    n_up   <- sum(res$regulation == "Up", na.rm = TRUE)
    n_down <- sum(res$regulation == "Down", na.rm = TRUE)
    deg_summary <- rbind(deg_summary, data.frame(
      Timepoint = tp,
      Direction = c("Up", "Down"),
      Count = c(n_up, -n_down),
      stringsAsFactors = FALSE
    ))
  }
  
  if (nrow(deg_summary) > 0) {
    deg_summary$Timepoint <- factor(deg_summary$Timepoint, levels = timepoints)
    
    p_bar <- ggplot(deg_summary, aes(x = Timepoint, y = Count, fill = Direction)) +
      geom_bar(stat = "identity", width = 0.6, position = "identity") +
      scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#3C5488")) +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(title = "DEG Counts by Timepoint",
           subtitle = sprintf("|log2FC| > %.1f, padj < %.2f",
                              PARAMS$lfc_cutoff, PARAMS$padj_cutoff),
           x = "Timepoint", y = "Number of DEGs (Down / Up)") +
      theme_bindlab()
    
    save_pub_fig(p_bar, "DEG_barplot", "05_DEG_vis")
  }
  
  # Volcano plots per timepoint
  for (tp in timepoints) {
    res <- all_results$deg_by_time[[tp]]
    if (is.null(res)) next
    
    p_vol <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
      geom_point(size = 0.8, alpha = 0.6) +
      scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#3C5488", "NS" = "grey70")) +
      geom_vline(xintercept = c(-PARAMS$lfc_cutoff, PARAMS$lfc_cutoff),
                 linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(PARAMS$padj_cutoff),
                 linetype = "dashed", color = "grey40") +
      labs(title = sprintf("Volcano Plot — %s", tp),
           x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
      theme_bindlab()
    
    save_pub_fig(p_vol, sprintf("Volcano_%s", gsub("[^a-zA-Z0-9]", "", tp)),
                 "05_DEG_vis", width = 7, height = 6)
  }
  
  log_step("05_VIS", "Step 05 COMPLETE (ADAPTED)")
}


# ------------------------------------------------------------------
# 08_WGCNA adapted
# ------------------------------------------------------------------
run_08_WGCNA_adapted <- function() {
  ensure_env()
  
  log_step("08_WGCNA", "Starting WGCNA (ADAPTED)...")
  
  # 注入正确的sample_info到全局环境
  si <- fix_sample_info_factors(readRDS(FILES$sample_info))
  assign("sample_info", si, envir = .GlobalEnv)
  
  # 直接source原08_WGCNA.R（已合并修复：power R²>=0.80 + dynamic group_levels）
  source(file.path(SCRIPT_DIR, "08_WGCNA.R"), local = FALSE)
  
  log_step("08_WGCNA", "Step 08 COMPLETE (ADAPTED)")
}


# ------------------------------------------------------------------
# 10_integration adapted
# ------------------------------------------------------------------
run_10_integration_adapted <- function() {
  # 10_integration.R L124: for (tp in c("4d","7d","14d","21d"))
  # 这段是用deg_by_time结果做注释，如果key不存在会返回NULL
  # 实际上deg_results现在保存了adapted的timepoints
  # 但10_integration硬编码了for loop → 结果丢失该注释
  # 影响：不致命，只是integration表里的per-timepoint LFC注释列为空
  
  # 策略：直接source，接受per-timepoint注释可能丢失
  # 核心筛选逻辑（ML + PPI union）不依赖这些注释
  
  cat("  [INFO] Running original 10_integration.R...\n")
  cat("  [INFO] Per-timepoint LFC annotation may be incomplete (hardcoded MSC timepoints)\n")
  source(file.path(SCRIPT_DIR, "10_integration.R"), local = FALSE)
  
  log_step("10_INT", "Step 10 COMPLETE (via wrapper)")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  主执行入口                                                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_all <- function() {
  cat("\n================================================================\n")
  cat(sprintf("  METI-FS Full Pipeline Run: %s\n", DATASET_ID))
  cat(sprintf("  Started: %s\n", format(Sys.time())))
  cat("================================================================\n")
  
  t0 <- Sys.time()
  
  run_phase1()
  run_phase2()
  run_phase3()
  
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  
  cat("\n================================================================\n")
  cat(sprintf("  COMPLETE: %s\n", DATASET_ID))
  cat(sprintf("  Total time: %.1f minutes\n", as.numeric(elapsed)))
  cat(sprintf("  Output: %s\n", PROFILE$project_dir))
  cat("================================================================\n")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  STRING数据库准备（辅助函数）                                               ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

#' 从MSC项目复制STRING本地数据库到当前GEO项目
#' 或者创建符号链接
setup_string_cache <- function(source_project = METHODS_BASE) {
  src_cache <- file.path(source_project, "data", "stringdb_cache")
  dst_cache <- file.path(PROFILE$project_dir, "data", "stringdb_cache")
  
  if (!dir.exists(src_cache)) {
    cat("  [ERROR] Source STRING cache not found: ", src_cache, "\n")
    return(invisible(FALSE))
  }
  
  if (dir.exists(dst_cache)) {
    cat("  [OK] Destination STRING cache already exists\n")
    return(invisible(TRUE))
  }
  
  dir.create(dst_cache, recursive = TRUE, showWarnings = FALSE)
  
  # 复制文件
  files_to_copy <- list.files(src_cache, full.names = TRUE)
  for (f in files_to_copy) {
    file.copy(f, file.path(dst_cache, basename(f)), overwrite = FALSE)
    cat(sprintf("  [COPIED] %s\n", basename(f)))
  }
  
  cat(sprintf("  [DONE] STRING cache prepared: %s\n", dst_cache))
  invisible(TRUE)
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  使用说明（直接运行本脚本时显示）                                            ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

if (sys.nframe() == 0) {
  cat("\n")
  cat("================================================================\n")
  cat("  S04_run_pipeline_wrapper.R — Usage\n")
  cat("================================================================\n")
  cat("\n")
  cat("  # 逐phase运行（推荐，可以检查中间结果）：\n")
  cat("  DATASET_ID <- 'GSE307424'    # 选择数据集\n")
  cat("  source('S04_run_pipeline_wrapper.R')  # 加载配置\n")
  cat("  run_phase1()                 # 数据导入+预处理+标准化\n")
  cat("  run_phase2()                 # DEG+maSigPro+WGCNA\n")
  cat("  run_phase3()                 # 候选筛选+整合\n")
  cat("\n")
  cat("  # 一键运行：\n")
  cat("  run_all()\n")
  cat("\n")
  cat("  # STRING数据库准备（09F需要）：\n")
  cat("  setup_string_cache()         # 从AdipoMSC复制\n")
  cat("\n")
  cat("  # 查看当前数据集profile：\n")
  cat("  str(PROFILE)\n")
  cat("\n")
  cat("================================================================\n")
}
