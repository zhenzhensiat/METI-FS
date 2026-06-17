# ==============================================================================
# S17_ppi_analysis.R — PPI Hub Validation using Local STRING v12.0
#
# 目的:
#   对缺PPI输出的GEO数据集(GSE197067, GSE236646, GSE150411),
#   使用本地STRING v12.0缓存运行PPI hub分析(与原pipeline Layer 6完全一致),
#   然后比较PPI-hub基因 vs ML-only基因的文献支持率。
#
# 输入:
#   - 各GEO数据集的 Final_candidate_genes*.csv
#   - Local STRING v12.0 cache (file.path(RUN_DIR, "GEO_GSE307424_Lung", "data", "stringdb_cache"))
#   - 手稿Table 5文献证据 (hardcoded)
#
# 输出:
#   - Console: 每个数据集的分步状态 + 跨数据集汇总表
#
# 文献依据:
#   - STRINGdb: Szklarczyk et al. 2019, Nucleic Acids Res
#   - 5 centrality metrics: 与原pipeline 09_PPI_hub.R 完全一致
# ==============================================================================
source("S_config.R")

suppressPackageStartupMessages({
  library(STRINGdb)
  library(igraph)
})

# ---- 1. 日志函数 (与原pipeline风格一致) ----
log_e6 <- function(msg) {
  ts <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] E6_PPI: %s\n", ts, msg))
}

# ---- 2. 文献证据 (手稿Table 5) ----
LIT <- list(
  GSE197067 = c(FKBP4="Direct", JARID2="Direct", TRAF3IP3="Direct",
                CYPA="Indirect", SMIM5="None", TMEM51="None",
                HRCT1="None", NXPH3="None"),
  GSE307424 = c(ADAMTS18="Direct", ADGRF4="Direct", GJB3="Direct",
                ITGA5="Direct", MRGPRX3="Indirect", KCNK3="Indirect",
                RNF152="Indirect", FAM83A="None"),
  GSE236646 = c(LAMP2="Direct", THBD="Direct", JAK1="Direct",
                STAT2="Direct", IFIT1="Indirect", IRF9="Indirect",
                MX1="Indirect", ADAR="None", ZNFX1="None"),
  GSE150411 = c(CD38="Direct", NKX3_2="Direct", NOD2="Direct",
                RIPK2="Direct", MMP13="Direct", COL2A1="Direct",
                ACAN="Direct", SOX9="Direct", COMP="Direct",
                ADAMTS4="Direct", ADAMTS5="Direct", IL6="Direct",
                TNF="Direct", PTGS2="Direct", MMP3="Indirect",
                MMP1="Indirect", TIMP1="Indirect", FN1="Indirect",
                VCAN="Indirect", CHAD="Indirect", OGN="Indirect",
                ASPN="Indirect", COL1A1="None", COL3A1="None",
                COL5A1="None", COL6A1="None", COL11A1="None",
                COL12A1="None", HAPLN1="None", LUM="None")
)

# ---- 3. 待处理数据集 ----
TARGETS <- list(
  GSE197067 = file.path(RUN_DIR, "GEO_GSE197067_Tcell"),
  GSE236646 = file.path(RUN_DIR, "GEO_GSE236646_NPC"),
  GSE150411 = file.path(RUN_DIR, "GEO_GSE150411_Chon")
)

STRING_CACHE <- file.path(RUN_DIR, "GEO_GSE307424_Lung/data/stringdb_cache")

# ---- 4. 初始化 STRINGdb (本地缓存, 与09_PPI_hub.R完全一致) ----
log_e6("==========================================")
log_e6("STEP 1: Initialize STRINGdb v12.0 (local cache)")
log_e6(sprintf("  Cache dir: %s", STRING_CACHE))
cached_files <- list.files(STRING_CACHE, pattern = "\\.gz$")
log_e6(sprintf("  Cached files: %d (%s)", length(cached_files), paste(cached_files, collapse=", ")))

string_db <- STRINGdb$new(
  version = "12.0",
  species = 9606,
  score_threshold = 400,
  input_directory = STRING_CACHE
)
log_e6("  STRINGdb initialized successfully")

# ---- 5. 逐数据集处理 ----
ppi_support <- c(Direct=0L, Indirect=0L, None=0L)
ml_support  <- c(Direct=0L, Indirect=0L, None=0L)
ppi_total <- 0L; ml_total <- 0L

for (ds_name in names(TARGETS)) {
  log_e6("==========================================")
  log_e6(sprintf("STEP 2.%s: Process %s", ds_name, ds_name))

  project_dir <- TARGETS[[ds_name]]
  data_dir <- file.path(project_dir, "data")

  # 2a. 加载最终候选基因
  final_files <- list.files(data_dir, pattern = "^Final_candidate", full.names = TRUE)
  if (length(final_files) == 0) {
    log_e6("  SKIP: no Final_candidate_genes*.csv found")
    next
  }
  fin <- read.csv(final_files[1])
  gene_col <- if ("symbol" %in% colnames(fin)) "symbol" else colnames(fin)[1]
  all_genes <- unique(as.character(fin[[gene_col]]))
  all_genes <- all_genes[!is.na(all_genes) & all_genes != ""]
  log_e6(sprintf("  Final genes: %d (col: %s)", length(all_genes), gene_col))

  # 2b. 映射到STRING
  log_e6("  Mapping genes to STRING IDs...")
  mapped <- tryCatch({
    string_db$map(data.frame(gene = all_genes, stringsAsFactors = FALSE),
                  "gene", removeUnmappedRows = TRUE)
  }, error = function(e) NULL)

  if (is.null(mapped) || nrow(mapped) == 0) {
    log_e6("  FAILED: STRING mapping returned 0 genes")
    next
  }
  log_e6(sprintf("  Mapped: %d/%d genes (%.1f%%)",
                 nrow(mapped), length(all_genes),
                 100 * nrow(mapped) / length(all_genes)))

  # 2c. 提取PPI子网络
  log_e6(sprintf("  Extracting PPI subnetwork (%d nodes)...", nrow(mapped)))
  subgraph <- tryCatch({
    string_db$get_subnetwork(mapped$STRING_id)
  }, error = function(e) NULL)

  if (is.null(subgraph)) {
    log_e6("  FAILED: could not extract subnetwork")
    next
  }
  n_nodes <- vcount(subgraph)
  n_edges <- ecount(subgraph)
  log_e6(sprintf("  PPI network: %d nodes, %d edges", n_nodes, n_edges))

  # 2d. 计算5个中心性指标
  log_e6("  Computing 5 centrality metrics...")
  g <- subgraph
  deg <- degree(g)
  bet <- betweenness(g, normalized = TRUE)
  pr  <- page_rank(g)$vector
  eig <- eigen_centrality(g)$vector
  hub <- hub_score(g)$vector
  log_e6(sprintf("  Centrality computed for %d nodes", n_nodes))

  # 2e. 识别hub基因 (top N, >=3 metrics)
  n_top <- min(length(all_genes), 5)
  top_deg <- names(sort(deg, decreasing = TRUE))[1:n_top]
  top_bet <- names(sort(bet, decreasing = TRUE))[1:n_top]
  top_pr  <- names(sort(pr,  decreasing = TRUE))[1:n_top]
  top_eig <- names(sort(eig, decreasing = TRUE))[1:n_top]
  top_hub <- names(sort(hub, decreasing = TRUE))[1:n_top]

  all_tops <- c(top_deg, top_bet, top_pr, top_eig, top_hub)
  hub_counts <- table(all_tops)
  hub_ids <- names(hub_counts[hub_counts >= 3])
  hub_symbols <- mapped$gene[match(hub_ids, mapped$STRING_id)]
  hub_symbols <- hub_symbols[!is.na(hub_symbols)]

  log_e6(sprintf("  PPI hub genes (>=3 metrics, top-%d): %d", n_top, length(hub_symbols)))
  if (length(hub_symbols) > 0) {
    cat("    ")
    cat(paste(hub_symbols, collapse = ", "))
    cat("\n")
  }

  # 2f. 与文献证据交叉对比
  ml_genes <- setdiff(toupper(all_genes), toupper(hub_symbols))
  ev <- LIT[[ds_name]]; ev_names <- toupper(names(ev))

  matched_ppi <- 0; matched_ml <- 0
  for (g in toupper(hub_symbols)) {
    idx <- which(ev_names == g)
    if (length(idx) > 0 && g %in% toupper(all_genes)) {
      ppi_support[ev[idx[1]]] <- ppi_support[ev[idx[1]]] + 1L
      ppi_total <- ppi_total + 1L; matched_ppi <- matched_ppi + 1
    }
  }
  for (g in ml_genes) {
    idx <- which(ev_names == g)
    if (length(idx) > 0) {
      ml_support[ev[idx[1]]] <- ml_support[ev[idx[1]]] + 1L
      ml_total <- ml_total + 1L; matched_ml <- matched_ml + 1
    }
  }
  log_e6(sprintf("  Literature matched: PPI=%d, ML=%d", matched_ppi, matched_ml))
}

# ---- 6. 加入GSE307424 (已有PPI输出) ----
log_e6("==========================================")
log_e6("STEP 3: Add GSE307424 (existing PPI output)")

fin424 <- read.csv(file.path(RUN_DIR, "GEO_GSE307424_Lung/data/Final_candidate_genesGSE307424.csv"))
ppi424 <- read.csv(file.path(RUN_DIR, "GEO_GSE307424_Lung/data/PPI_09F_hub_genesGSE307424.csv"))
ppi_genes_424 <- toupper(as.character(ppi424$symbol))
ml_genes_424 <- setdiff(toupper(fin424$symbol), ppi_genes_424)
ev424 <- LIT$GSE307424; ev424_names <- toupper(names(ev424))

for (g in ppi_genes_424) {
  idx <- which(ev424_names == g)
  if (length(idx) > 0) {
    ppi_support[ev424[idx[1]]] <- ppi_support[ev424[idx[1]]] + 1L
    ppi_total <- ppi_total + 1L
  }
}
for (g in ml_genes_424) {
  idx <- which(ev424_names == g)
  if (length(idx) > 0) {
    ml_support[ev424[idx[1]]] <- ml_support[ev424[idx[1]]] + 1L
    ml_total <- ml_total + 1L
  }
}
log_e6(sprintf("  PPI hub: %d, ML-only: %d", length(ppi_genes_424), length(ml_genes_424)))

# ---- 7. 跨数据集汇总 ----
log_e6("==========================================")
log_e6("STEP 4: Cross-Dataset Summary")
cat("\n")
cat(sprintf("  %-10s %8s %8s %8s %8s\n", "Group", "Direct", "Indirect", "None", "Total"))
cat(sprintf("  %-10s %8d %8d %8d %8d\n",
            "PPI-hub", ppi_support[1], ppi_support[2], ppi_support[3], ppi_total))
cat(sprintf("  %-10s %8d %8d %8d %8d\n\n",
            "ML-only", ml_support[1], ml_support[2], ml_support[3], ml_total))

if (ppi_total > 0) {
  ppi_any <- ppi_support["Direct"] + ppi_support["Indirect"]
  cat(sprintf("  PPI-hub support rate: %d/%d = %.1f%%\n",
              ppi_any, ppi_total, 100 * ppi_any / ppi_total))
}
if (ml_total > 0) {
  ml_any <- ml_support["Direct"] + ml_support["Indirect"]
  cat(sprintf("  ML-only support rate: %d/%d = %.1f%%\n",
              ml_any, ml_total, 100 * ml_any / ml_total))
}

if (ppi_total > 1 && ml_total > 1) {
  tbl <- matrix(c(ppi_any, ppi_total - ppi_any,
                   ml_any,  ml_total  - ml_any), nrow = 2)
  ft <- fisher.test(tbl)
  cat(sprintf("  Fisher exact test p = %.4f\n", ft$p.value))
} else {
  cat("  (insufficient data for Fisher test)\n")
}

log_e6("E6 complete")
