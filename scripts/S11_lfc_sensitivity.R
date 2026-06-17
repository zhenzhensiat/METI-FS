#!/usr/bin/env Rscript
# ==============================================================================
# S11_lfc_sensitivity.R — lfcThreshold Sensitivity Analysis
#
# 文献: Love, Huber & Anders 2014, Genome Biology 15:550
#   "If lfcThreshold is specified, the results are for Wald tests."
#   "DESeq2 offers tests for composite null hypotheses |beta| <= theta."
#
# 注意: dds_object.rds 是 DESeq() 前的原始状态 (无size factors/dispersions)
#   需完整运行 DESeq(dds, test="Wald") 获取 coefficient names 用于 lfcThreshold
# ==============================================================================

if (file.exists("S_config.R")) source("S_config.R")
suppressPackageStartupMessages({library(DESeq2); library(ggplot2)})

LFC_DATASETS <- list(
  GSE197067 = list(file = file.path(RUN_DIR, "GEO_GSE197067_Tcell/data/dds_object.rds"),
                   label = "T cell activation", n = 40, n_tp = 5),
  GSE307424 = list(file = file.path(RUN_DIR, "GEO_GSE307424_Lung/data/dds_object.rds"),
                   label = "SMARCA2 degrader", n = 18, n_tp = 3),
  GSE236646 = list(file = file.path(RUN_DIR, "GEO_GSE236646_NPC/data/dds_object.rds"),
                   label = "NPC viral infection", n = 16, n_tp = 3),
  GSE150411 = list(file = file.path(RUN_DIR, "GEO_GSE150411_Chon/data/dds_object.rds"),
                   label = "Chondrocyte inflammation", n = 18, n_tp = 3)
)

run_lfc_sensitivity <- function(thresholds = c(0, 0.5, 1.0, 1.5), verbose = TRUE) {
  all_rows <- list()

  for (ds_name in names(LFC_DATASETS)) {
    ds <- LFC_DATASETS[[ds_name]]
    if (!file.exists(ds$file)) { cat(sprintf("[SKIP] %s\n", ds_name)); next }

    cat(sprintf("\n%s (%s): fitting DESeq2 Wald...\n", ds_name, ds$label))
    dds <- readRDS(ds$file)
    dds <- DESeq(dds, test = "Wald")  # full pipeline: sizeFactors → dispersions → Wald

    for (thr in thresholds) {
      res <- if (thr == 0) {
        results(dds, contrast = c("Treatment", "Induced", "Control"), alpha = 0.05)
      } else {
        results(dds, contrast = c("Treatment", "Induced", "Control"),
                lfcThreshold = thr, altHypothesis = "greaterAbs", alpha = 0.05)
      }
      n_pass <- sum(res$padj < 0.05, na.rm = TRUE)

      all_rows[[length(all_rows) + 1]] <- data.frame(
        dataset = ds_name, label = ds$label,
        n_samples = ds$n, n_tp = ds$n_tp,
        lfc_threshold = thr, n_pass = n_pass,
        stringsAsFactors = FALSE
      )
      if (verbose) cat(sprintf("  thr=%.1f: %d genes\n", thr, n_pass))
    }
  }

  results_df <- do.call(rbind, all_rows)

  dir.create(TAB_DIR_METHODS, recursive = TRUE, showWarnings = FALSE)
  write.csv(results_df, file.path(TAB_DIR_METHODS, "Supp_lfcThreshold_sensitivity.csv"),
            row.names = FALSE)

  p <- ggplot(results_df, aes(x = factor(lfc_threshold), y = n_pass,
                               group = dataset, color = dataset)) +
    geom_line(linewidth = 1) + geom_point(size = 3) +
    facet_wrap(~ label, scales = "free_y") +
    labs(title = "Sensitivity of Gene Counts to lfcThreshold",
         subtitle = "DESeq2 Wald test: |log2FC| > threshold [Love et al. 2014, Genome Biology]",
         x = "lfcThreshold", y = "Genes passing") +
    theme_classic(12) + theme(legend.position = "bottom")

  dir.create(FIG_DIR_METHODS, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(FIG_DIR_METHODS, "Supp_lfcThreshold_sensitivity.pdf"), p, width = 10, height = 6)
  ggsave(file.path(FIG_DIR_METHODS, "Supp_lfcThreshold_sensitivity.png"), p, width = 10, height = 6, dpi = 300)

  return(invisible(results_df))
}

if (sys.nframe() == 0) {
  cat("\n========== S11: lfcThreshold Sensitivity ==========\n")
  results <- run_lfc_sensitivity()
  cat("\n"); print(results); cat("\n[DONE]\n")
}
