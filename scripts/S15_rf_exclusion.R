# S15_rf_exclusion.R — RF Auto-Exclusion Status from GEO pipeline bootstrap results
# Data structure: ml_stability_selection.rds contains:
#   lasso_freq, svm_freq, rf_freq: named numeric vectors (gene -> frequency)
#   lasso_selection_matrix, svm_selection_matrix, rf_selection_matrix: genes x bootstraps
source("S_config.R")

geo_labels <- c(
  GSE197067_Tcell = "T cell activation",
  GSE307424_Lung = "SMARCA2 degrader",
  GSE236646_NPC = "NPC viral infection",
  GSE150411_Chon = "Chondrocyte inflammation"
)

compute_stability <- function(freq_vec) {
  # Mean frequency of top-10% most selected genes (proxy for Nogueira)
  if (is.null(freq_vec) || all(freq_vec == 0)) return(0)
  top_cut <- quantile(freq_vec[freq_vec > 0], 0.9, na.rm = TRUE)
  top_genes <- freq_vec[freq_vec >= top_cut]
  if (length(top_genes) < 2) return(mean(freq_vec[freq_vec > 0]))
  mean(top_genes)
}

cat("=== RF Auto-Exclusion Status per GEO Dataset ===\n\n")
cat(sprintf("%-25s %8s %8s %8s %10s\n", "Dataset", "LASSO", "SVM", "RF", "Excluded?"))

for (nm in names(geo_labels)) {
  stab_file <- sprintf(file.path(RUN_DIR, "GEO_%s/data/ml_stability_selection.rds"), nm)
  if (!file.exists(stab_file)) {
    cat(sprintf("%-25s (no stability data)\n", nm))
    next
  }
  stab <- readRDS(stab_file)

  n_lasso <- compute_stability(stab$lasso_freq)
  n_svm   <- compute_stability(stab$svm_freq)
  n_rf    <- compute_stability(stab$rf_freq)

  excluded <- (n_rf < 0.5 && n_rf < min(n_lasso, n_svm))

  cat(sprintf("%-25s %8.3f %8.3f %8.3f %10s\n",
              geo_labels[[nm]], n_lasso, n_svm, n_rf,
              if (excluded) "YES" else "NO"))
}
cat("\n=== T6 DONE ===\n")
