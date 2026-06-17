# ==============================================================================
# S16_wgcna_enrichment.R — WGCNA Module GO/KEGG Enrichment Analysis
#
# Validates that WGCNA modules at n=16-18 capture biologically relevant functions.
# ==============================================================================
source("S_config.R")

suppressPackageStartupMessages({
  library(WGCNA)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
allowWGCNAThreads(2)

# ---- 1. 日志 ----
log_e5 <- function(msg) {
  cat(sprintf("[%s] E5_ENRICH: %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg))
}

# ---- 2. 数据集配置 ----
DSETS <- list(
  GSE236646 = list(
    label = "NPC viral infection (n=16)",
    vst = file.path(RUN_DIR, "GEO_GSE236646_NPC/data/vst_matrix.rds"),
    si  = file.path(RUN_DIR, "GEO_GSE236646_NPC/data/sample_info.rds"),
    context = c("HSV", "herpes", "viral", "infection", "antiviral", "immune")
  ),
  GSE307424 = list(
    label = "SMARCA2 degrader (n=18)",
    vst = file.path(RUN_DIR, "GEO_GSE307424_Lung/data/vst_matrix.rds"),
    si  = file.path(RUN_DIR, "GEO_GSE307424_Lung/data/sample_info.rds"),
    context = c("lung", "cancer", "NSCLC", "SMARCA2", "apoptosis", "cell cycle")
  )
)

all_results <- list()

for (ds_name in names(DSETS)) {
  ds <- DSETS[[ds_name]]
  log_e5(sprintf("=========================================="))
  log_e5(sprintf("STEP: %s — %s", ds_name, ds$label))

  # 2a. 加载数据
  vst <- readRDS(ds$vst)
  si  <- readRDS(ds$si)
  n <- ncol(vst)

  # MAD top50% (与08_WGCNA.R一致)
  gene_mad <- apply(vst, 1, mad)
  cutoff <- sort(gene_mad, decreasing = TRUE)[round(nrow(vst) * 0.5)]
  genes_use <- names(gene_mad)[gene_mad >= cutoff]
  datExpr <- t(vst[genes_use, ])
  log_e5(sprintf("  Input: %d genes x %d samples (MAD top50%%)", ncol(datExpr), n))

  # 2b. pick soft threshold
  sft <- pickSoftThreshold(datExpr, powerVector = 1:20,
                           networkType = "signed hybrid", verbose = 0)
  r2 <- max(-sign(sft$fitIndices[,3]) * sft$fitIndices[,2], na.rm = TRUE)
  if (r2 >= 0.8) {
    power <- sft$fitIndices$Power[which.max(-sign(sft$fitIndices[,3]) * sft$fitIndices[,2])]
  } else {
    power <- if (n < 20) 9 else 8
  }
  log_e5(sprintf("  Power: %d (R²=%.3f, %s)", power, r2, if(r2>=0.8) "auto" else "FAQ"))

  # 2c. Build network + detect modules (fixed seed for reproducibility)
  set.seed(42)
  adj <- adjacency(datExpr, power = power, type = "signed hybrid")
  TOM <- TOMsimilarity(adj)
  tree <- hclust(as.dist(1 - TOM), method = "average")
  mods <- cutreeDynamic(dendro = tree, distM = 1 - TOM, deepSplit = 2,
                         minClusterSize = 30, pamRespectsDendro = FALSE)
  names(mods) <- genes_use
  cols <- labels2colors(mods)
  n_mods <- length(unique(cols))
  log_e5(sprintf("  Modules detected: %d", n_mods))

  # 2d. Module-trait correlation
  MEs <- moduleEigengenes(datExpr, cols)$eigengenes
  trt <- as.numeric(si$Treatment == "Induced")
  pv <- corPvalueStudent(cor(MEs, trt, use = "p"), n)
  sig_mods <- names(which(pv[,1] < 0.05))
  log_e5(sprintf("  Condition-significant modules: %d/%d", length(sig_mods), n_mods))

  if (length(sig_mods) == 0) {
    log_e5("  No significant modules — skipping enrichment")
    next
  }

  # 2e. Extract genes using mods (numeric) not cols (color labels)
  me_to_mod <- setNames(match(gsub("^ME", "", colnames(MEs)), unique(cols)), colnames(MEs))
  sig_genes_all <- unique(unlist(lapply(sig_mods, function(sm) {
    mod_num <- me_to_mod[sm]
    names(mods)[mods == mod_num]
  })))
  sig_genes_all <- sig_genes_all[!is.na(sig_genes_all)]
  log_e5(sprintf("  Genes in sig modules: %d", length(sig_genes_all)))

  # 2f. GO BP enrichment (auto-detect ID type)
  log_e5("  Running GO BP enrichment...")
  # Detect gene ID type
  key_type <- if (any(grepl("^ENSG", sig_genes_all[1:min(10, length(sig_genes_all))]))) "ENSEMBL" else "SYMBOL"
  log_e5(sprintf("  Gene ID type: %s", key_type))
  ego <- tryCatch({
    enrichGO(gene = sig_genes_all, OrgDb = org.Hs.eg.db,
             keyType = key_type, ont = "BP",
             pAdjustMethod = "BH", pvalueCutoff = 0.01, qvalueCutoff = 0.05)
  }, error = function(e) NULL)

  if (!is.null(ego) && nrow(ego@result) > 0) {
    n_go <- nrow(ego@result)
    log_e5(sprintf("  GO BP terms enriched: %d", n_go))

    # Show top terms and check context relevance
    top_terms <- head(ego@result[order(ego@result$pvalue), ], 10)
    cat(sprintf("  Top 10 GO terms:\n"))
    for (i in 1:nrow(top_terms)) {
      cat(sprintf("    %s (p=%.1e)\n", top_terms$Description[i], top_terms$pvalue[i]))
    }

    # Count context-relevant terms
    context_hits <- sum(sapply(ds$context, function(kw) {
      sum(grepl(kw, ego@result$Description, ignore.case = TRUE))
    }))
    log_e5(sprintf("  Context-relevant terms: %d", context_hits))
  } else {
    log_e5("  No significant GO BP enrichment found")
    n_go <- 0; context_hits <- 0
  }

  all_results[[ds_name]] <- data.frame(
    dataset = ds_name, label = ds$label, n_samples = n,
    n_genes_MAD = ncol(datExpr), power = power, R2 = r2,
    n_modules = n_mods, n_sig_modules = length(sig_mods),
    n_sig_genes = length(sig_genes_all),
    n_GO_BP = n_go, n_context_GO = context_hits,
    stringsAsFactors = FALSE
  )
}

# ---- 3. Summary ----
log_e5("==========================================")
log_e5("SUMMARY")
df <- do.call(rbind, all_results)
print(df)

write.csv(df, file.path(TAB_DIR_METHODS, "Supp_WGCNA_module_enrichment.csv"), row.names = FALSE)
log_e5(sprintf("Saved: tables/Supp_WGCNA_module_enrichment.csv"))
log_e5("E5 complete")
