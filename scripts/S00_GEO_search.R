#!/usr/bin/env Rscript
# ==============================================================================
# S00_GEO_search_v2.R — 增强版GEO时序RNA-seq数据集搜索
#
# 目的: 为METI-FS方法学文章寻找额外的外部验证数据集
#
# Pipeline硬性要求:
#   1. Homo sapiens bulk RNA-seq (非scRNA-seq, 非microarray)
#   2. Treatment vs Control 两组设计
#   3. ≥3个时间点
#   4. 每组每时间点 ≥2个生物学重复
#   5. 总样本量 ≥18 (WGCNA要求)
#   6. Raw counts可获取 (NCBI预计算 或 作者提供)
#
# 已有数据集 (不重复搜索):
#   - GSE197067 (T细胞激活, 40样本, 5TP) ✓
#   - GSE307424 (肺癌SMARCA2, 18样本, 3TP) ✓
#   - GSE303975 (前列腺癌, 18样本, 3TP) — 已注册但未跑
#
# 搜索策略: "宽进严出"
#   Step 1: 用多组关键词在GEO搜索, 收集候选GSE
#   Step 2: 解析每个GSE的样本数和元数据
#   Step 3: 输出候选列表, 标注结构适合度, 供人工审核
#
# 用法:
#   Usage: source("scripts/S_config.R") first
#   source("S00_GEO_search_v2.R")
#   # Results saved to search_results/
# ==============================================================================

# ---- 0. 依赖 ----
for (pkg in c("rentrez", "xml2", "dplyr", "stringr")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
library(rentrez)
library(xml2)
library(dplyr)
library(stringr)

# ---- 1. 配置 ----
# 如果S_config.R存在就加载, 否则直接定义输出路径
if (file.exists("S_config.R")) {
  source("S_config.R")
  OUT_DIR <- SEARCH_DIR
} else {
  OUT_DIR <- "."
}
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 已有数据集 — 跳过
KNOWN_GSE <- c("GSE197067", "GSE307424", "GSE303975")

cat("============================================================\n")
cat("  METI-FS GEO Dataset Search v2\n")
cat("  Requirements: Human, bulk RNA-seq, treatment vs control,\n")
cat("                ≥3 timepoints, ≥2 bio reps/group/TP, ≥18 samples\n")
cat("  Output: ", OUT_DIR, "\n")
cat("============================================================\n\n")

# ---- 2. 搜索查询组 ----
# 使用多组关键词覆盖不同领域
queries <- list(
  # 药物/处理响应时序
  drug_timecourse = paste0(
    '"Homo sapiens"[Organism] AND "Expression profiling by high throughput sequencing"[DataSet Type] ',
    'AND ("time course"[Description] OR "time series"[Description] OR "time point"[Description]) ',
    'AND ("drug"[Description] OR "treatment"[Description] OR "treated"[Description]) ',
    'AND ("control"[Description] OR "untreated"[Description] OR "vehicle"[Description] OR "DMSO"[Description])'
  ),
  
  # 刺激响应 (细胞因子/生长因子)
  stimulus_response = paste0(
    '"Homo sapiens"[Organism] AND "Expression profiling by high throughput sequencing"[DataSet Type] ',
    'AND ("time course"[Description] OR "time series"[Description]) ',
    'AND ("stimulat"[Description] OR "cytokine"[Description] OR "growth factor"[Description]) ',
    'AND ("unstimulated"[Description] OR "control"[Description])'
  ),
  
  # 分化实验
  differentiation = paste0(
    '"Homo sapiens"[Organism] AND "Expression profiling by high throughput sequencing"[DataSet Type] ',
    'AND ("time course"[Description] OR "time series"[Description] OR "time point"[Description]) ',
    'AND ("differentiation"[Description] OR "differentiated"[Description]) ',
    'AND ("undifferentiated"[Description] OR "control"[Description])'
  ),
  
  # 感染模型
  infection = paste0(
    '"Homo sapiens"[Organism] AND "Expression profiling by high throughput sequencing"[DataSet Type] ',
    'AND ("time course"[Description] OR "time series"[Description]) ',
    'AND ("infect"[Description] OR "virus"[Description] OR "pathogen"[Description]) ',
    'AND ("mock"[Description] OR "uninfected"[Description] OR "control"[Description])'
  ),
  
  # 宽泛: 任何 time course + treatment
  broad_tc = paste0(
    '"Homo sapiens"[Organism] AND "Expression profiling by high throughput sequencing"[DataSet Type] ',
    'AND "time course"[Description] ',
    'AND ("replicate"[Description] OR "biological replicate"[Description] OR "triplicate"[Description])'
  )
)

# ---- 3. 执行搜索 ----
all_gse_ids <- character(0)

for (qname in names(queries)) {
  cat(sprintf("[SEARCH] Query: %s ...\n", qname))
  
  tryCatch({
    # 搜索GEO DataSets
    search_result <- entrez_search(
      db = "gds",
      term = queries[[qname]],
      retmax = 200,
      use_history = FALSE
    )
    
    cat(sprintf("  Found: %d results (showing up to 200)\n", search_result$count))
    
    if (search_result$count > 0 && length(search_result$ids) > 0) {
      all_gse_ids <- c(all_gse_ids, search_result$ids)
    }
    
    Sys.sleep(0.5)  # NCBI rate limit
    
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}

all_gse_ids <- unique(all_gse_ids)
cat(sprintf("\n[TOTAL] %d unique GDS/GSE IDs from all queries\n\n", length(all_gse_ids)))

# ---- 4. 获取每个数据集的详细信息 ----
cat("[PARSE] Fetching details for each dataset...\n")

results <- list()
batch_size <- 20

for (batch_start in seq(1, length(all_gse_ids), by = batch_size)) {
  batch_end <- min(batch_start + batch_size - 1, length(all_gse_ids))
  batch_ids <- all_gse_ids[batch_start:batch_end]
  
  cat(sprintf("  Batch %d-%d / %d\n", batch_start, batch_end, length(all_gse_ids)))
  
  tryCatch({
    summaries <- entrez_summary(db = "gds", id = batch_ids)
    
    # entrez_summary可能返回单个或列表
    if (inherits(summaries, "esummary")) {
      summaries <- list(summaries)
    }
    
    for (s in summaries) {
      tryCatch({
        gse_acc  <- ifelse(!is.null(s$accession), s$accession, NA)
        title    <- ifelse(!is.null(s$title), s$title, NA)
        summary  <- ifelse(!is.null(s$summary), s$summary, "")
        gpl      <- ifelse(!is.null(s$gpl), s$gpl, NA)
        n_samples <- ifelse(!is.null(s$n_samples), as.integer(s$n_samples), NA)
        taxon    <- ifelse(!is.null(s$taxon), s$taxon, NA)
        gdstype  <- ifelse(!is.null(s$gdstype), s$gdstype, NA)
        pdat     <- ifelse(!is.null(s$pdat), s$pdat, NA)
        
        # 跳过非GSE (GDS记录)
        if (is.na(gse_acc)) next
        
        # 跳过已知数据集
        if (gse_acc %in% KNOWN_GSE) next
        
        # 跳过非人类
        if (!is.na(taxon) && !grepl("Homo sapiens", taxon, ignore.case = TRUE)) next
        
        # 跳过样本太少
        if (!is.na(n_samples) && n_samples < 18) next
        
        # 启发式: 检查描述中是否包含时间和处理关键词
        desc_lower <- tolower(paste(title, summary))
        
        has_timewords <- grepl("time.?course|time.?series|time.?point|\\d+\\s*h(our|r)?|\\d+\\s*d(ay)?|\\d+\\s*min", desc_lower)
        has_treatment <- grepl("treat|drug|stimulat|infect|induc|differenti|expos|challenged|agonist|inhibitor|compound", desc_lower)
        has_control   <- grepl("control|untreated|vehicle|dmso|mock|uninfect|unstimulat|undifferent|baseline", desc_lower)
        has_replicate <- grepl("replic|triplic|duplic|\\bn\\s*=\\s*[2-9]|\\bn\\s*=\\s*\\d{2}", desc_lower)
        
        # 排除scRNA-seq
        is_sc <- grepl("single.?cell|scRNA|10x|chromium|droplet|smart-?seq", desc_lower)
        
        # 计算适合度分数
        score <- sum(c(has_timewords, has_treatment, has_control, has_replicate, !is_sc))
        
        results[[length(results) + 1]] <- data.frame(
          GSE = gse_acc,
          Title = substr(title, 1, 120),
          n_samples = n_samples,
          Score = score,
          TimeWords = has_timewords,
          Treatment = has_treatment,
          Control = has_control,
          Replicate = has_replicate,
          SingleCell = is_sc,
          PubDate = pdat,
          Summary_short = substr(summary, 1, 300),
          stringsAsFactors = FALSE
        )
        
      }, error = function(e) NULL)
    }
    
    Sys.sleep(0.4)  # rate limit
    
  }, error = function(e) {
    cat(sprintf("  Batch error: %s\n", e$message))
  })
}

# ---- 5. 汇总与排序 ----
if (length(results) == 0) {
  cat("\n[WARN] No results found. Check network connection and queries.\n")
} else {
  df <- bind_rows(results) %>%
    distinct(GSE, .keep_all = TRUE) %>%
    filter(!SingleCell) %>%        # 排除single-cell
    filter(n_samples >= 18) %>%    # 样本量下限
    filter(Score >= 3) %>%         # 至少3/5适合度
    arrange(desc(Score), desc(n_samples))
  
  cat(sprintf("\n[RESULTS] %d candidate datasets after filtering\n", nrow(df)))
  cat(sprintf("  Score 5 (perfect): %d\n", sum(df$Score == 5)))
  cat(sprintf("  Score 4 (good):    %d\n", sum(df$Score == 4)))
  cat(sprintf("  Score 3 (maybe):   %d\n", sum(df$Score == 3)))
  
  # 输出完整结果
  out_file <- file.path(OUT_DIR, sprintf("GEO_candidates_v2_%s.csv", format(Sys.Date(), "%Y%m%d")))
  write.csv(df, out_file, row.names = FALSE)
  cat(sprintf("\n[SAVED] Full results: %s\n", out_file))
  
  # 输出top candidates的详细信息
  cat("\n============================================================\n")
  cat("  TOP CANDIDATES (Score >= 4, n_samples >= 18)\n")
  cat("============================================================\n\n")
  
  top <- df %>% filter(Score >= 4)
  if (nrow(top) > 0) {
    for (i in seq_len(min(30, nrow(top)))) {
      r <- top[i, ]
      cat(sprintf("--- %s (n=%d, Score=%d, %s) ---\n", r$GSE, r$n_samples, r$Score, r$PubDate))
      cat(sprintf("    %s\n", r$Title))
      cat(sprintf("    Time:%s Treat:%s Ctrl:%s Rep:%s\n",
                  ifelse(r$TimeWords, "✓", "✗"),
                  ifelse(r$Treatment, "✓", "✗"),
                  ifelse(r$Control, "✓", "✗"),
                  ifelse(r$Replicate, "✓", "✗")))
      cat(sprintf("    %s\n\n", r$Summary_short))
    }
  } else {
    cat("  No Score>=4 candidates. Check Score>=3 in output CSV.\n")
  }
}

# ---- 6. 人工审核指南 ----
cat("\n============================================================\n")
cat("  MANUAL REVIEW GUIDE\n")
cat("============================================================\n")
cat("
对每个候选GSE, 需要人工确认以下内容:

1. 到 https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSExxxxxx 查看
2. 确认是 bulk RNA-seq (不是 scRNA-seq / microarray / ATAC-seq 等)
3. 确认有 treatment vs control 两组 (不是 dose-response / multi-arm)
4. 数清楚:
   - 时间点数 (≥3)
   - 每组每时间点的重复数 (≥2, 最好≥3)
   - 总样本数 (≥18)
5. 确认 raw counts 可获取:
   - 方法1: NCBI预计算counts (大多数近年human RNA-seq都有)
     检查: https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSExxxxxx
   - 方法2: 作者在supplementary files中提供了count矩阵
6. 如果是multi-arm (多种药物), 选取其中一个arm + vehicle/control

适合的数据集格式示例 (可直接加入S_config.R):

  GSExxxxxx = list(
    gse = 'GSExxxxxx',
    prefix = 'ShortName',
    domain = 'cancer_drug_response / immune_activation / differentiation / ...',
    n_samples = XX,
    n_timepoints = X,
    time_labels = c('Xh', 'Xh', ...),
    design = 'CellLine, DrugA vs DMSO x 3 timepoints x 3 reps',
    note = '...'
  )

请将审核通过的GSE告知Claude, 我来写S03b适配脚本并运行pipeline。
")

cat("\n[DONE] Search complete.\n")
