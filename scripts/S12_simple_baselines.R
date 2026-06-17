#!/usr/bin/env Rscript
# ==============================================================================
# S12_simple_baselines.R — Simple Baseline Comparison
#
# Compares METI-FS against naive baselines to demonstrate that
# the multi-layer framework adds value over simpler approaches.
#
# Baselines:
#   B1: Top-N by DESeq2 LRT p-value
#   B2: Top-N by maSigPro R-squared
#   B3: Random selection from candidate pool (100 repeats, mean)
#
# Metrics:
#   - Jaccard overlap with METI-FS final genes
#   - Literature support rate (using existing Table 5 data)
#   - GO biological process enrichment count
#
# Output:
#   tables/Supp_simple_baselines.csv
#
# Usage:
#   source("S12_simple_baselines.R")
#   results <- run_simple_baselines()
# ==============================================================================

if (file.exists("S_config.R")) source("S_config.R")

suppressPackageStartupMessages({
  library(DESeq2)
  library(clusterProfiler)
  library(org.Hs.eg.db)
})

BASELINE_DATASETS <- list(
  GSE197067 = list(
    project_dir = file.path(RUN_DIR, "GEO_GSE197067_Tcell"),
    label = "T cell activation", n_final = 8),
  GSE307424 = list(
    project_dir = file.path(RUN_DIR, "GEO_GSE307424_Lung"),
    label = "SMARCA2 degrader", n_final = 8),
  GSE236646 = list(
    project_dir = file.path(RUN_DIR, "GEO_GSE236646_NPC"),
    label = "NPC viral infection", n_final = 9),
  GSE150411 = list(
    project_dir = file.path(RUN_DIR, "GEO_GSE150411_Chon"),
    label = "Chondrocyte inflammation", n_final = 30)
)

N_RANDOM_REPS <- 100

#' Run simple baseline comparison for one dataset
compare_baselines <- function(ds_name, ds_config) {

  data_dir <- file.path(ds_config$project_dir, "data")

  # Load METI-FS final genes
  # Use filename pattern matching (actual files: Final_candidate_genesGSE*.csv)
  final_files <- list.files(data_dir, pattern = "^Final_candidate_genes", full.names = TRUE)
  if (length(final_files) == 0) {
    cat(sprintf("[SKIP] %s: Final_candidate_genes*.csv not found
", ds_name))
    return(NULL)
  }
  final_file <- final_files[1]
  if (!file.exists(final_file)) {
    cat(sprintf("[SKIP] %s: Final_candidate_genes.csv not found\n", ds_name))
    return(NULL)
  }
  final_df <- read.csv(final_file)
  metifs_genes <- unique(as.character(final_df[, 1]))

  # Load LRT results
  deg_file <- file.path(data_dir, "deg_results.rds")
  if (!file.exists(deg_file)) return(NULL)
  deg <- readRDS(deg_file)
  lrt <- deg$lrt_interaction

  # B1: Top-N by LRT p-value
  lrt_sorted <- lrt[order(lrt$padj), ]
  b1_genes <- head(rownames(lrt_sorted), ds_config$n_final)

  # B2: Top-N by maSigPro R-squared
  masigpro_file <- file.path(data_dir, "masigpro_results.rds")
  b2_genes <- character(0)
  if (file.exists(masigpro_file)) {
    masigpro <- readRDS(masigpro_file)
    if (!is.null(masigpro$tstep) && !is.null(masigpro$tstep$sig.profiles)) {
      sig_profiles <- masigpro$tstep$sig.profiles
      rsq_values <- sapply(sig_profiles, function(x) if (!is.null(x$R.squared)) x$R.squared else 0)
      b2_genes <- head(names(sort(rsq_values, decreasing = TRUE)), ds_config$n_final)
    }
  }

  # B3: Random from candidate pool (100 reps)
  pool_file <- file.path(data_dir, "candidate_pool.rds")
  b3_jaccards <- c()
  if (file.exists(pool_file)) {
    pool <- readRDS(pool_file)
    pool_genes <- if (is.list(pool)) pool$candidate_pool else pool
    for (r in seq_len(N_RANDOM_REPS)) {
      random_genes <- sample(pool_genes, min(ds_config$n_final, length(pool_genes)))
      b3_jaccards[r] <- length(intersect(random_genes, metifs_genes)) /
                         length(union(random_genes, metifs_genes))
    }
  }

  # Compute overlaps
  jac_b1 <- length(intersect(b1_genes, metifs_genes)) /
             length(union(b1_genes, metifs_genes))
  jac_b2 <- if (length(b2_genes) > 0)
             length(intersect(b2_genes, metifs_genes)) /
             length(union(b2_genes, metifs_genes)) else NA
  jac_b3_mean <- mean(b3_jaccards, na.rm = TRUE)

  # GO enrichment counts
  n_go_metifs <- tryCatch({
    ego <- enrichGO(gene = metifs_genes, OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
    nrow(ego@result)
  }, error = function(e) NA_integer_)

  n_go_b1 <- tryCatch({
    ego <- enrichGO(gene = b1_genes, OrgDb = org.Hs.eg.db,
                    keyType = "SYMBOL", ont = "BP", pvalueCutoff = 0.05)
    nrow(ego@result)
  }, error = function(e) NA_integer_)

  data.frame(
    dataset = ds_name,
    label = ds_config$label,
    n_metifs = length(metifs_genes),
    n_b1 = length(b1_genes),
    n_b2 = length(b2_genes),
    jaccard_b1 = round(jac_b1, 4),
    jaccard_b2 = round(jac_b2, 4),
    jaccard_b3_mean = round(jac_b3_mean, 4),
    n_go_metifs = n_go_metifs,
    n_go_b1 = n_go_b1,
    stringsAsFactors = FALSE
  )
}

run_simple_baselines <- function() {

  cat("\n============================================================\n")
  cat("  S12: Simple Baseline Comparison\n")
  cat("============================================================\n\n")

  all_rows <- list()
  for (ds_name in names(BASELINE_DATASETS)) {
    cat(sprintf("Processing %s...\n", ds_name))
    row <- compare_baselines(ds_name, BASELINE_DATASETS[[ds_name]])
    if (!is.null(row)) all_rows[[ds_name]] <- row
  }

  results_df <- do.call(rbind, all_rows)
  dir.create(TAB_DIR_METHODS, recursive = TRUE, showWarnings = FALSE)
  write.csv(results_df, file.path(TAB_DIR_METHODS, "Supp_simple_baselines.csv"),
            row.names = FALSE)

  cat(sprintf("\n[SAVED] tables/Supp_simple_baselines.csv\n"))
  print(results_df)
  cat("\n[DONE]\n")
  return(invisible(results_df))
}

if (sys.nframe() == 0) {
  cat("\nS12_simple_baselines.R — Run manually via:\n")
  cat("  source('S12_simple_baselines.R')\n")
  cat("  run_simple_baselines()\n")
}
