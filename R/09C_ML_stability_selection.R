#!/usr/bin/env Rscript
# ==============================================================================
# 09C_ML_stability_selection.R — Bootstrap Stability Selection ()
#
# and 09B_ML_feature_selection.RRun， using 。
# ： types MLN times bootstrap，Statistics Gene in N times in 
# by in （stability score）， to ""。
#
# Literature basis:
#   - Meinshausen & Bühlmann (2010) JRSS-B: Stability Selection
#   - Abeel et al. (2011) Pattern Recog Letters: ensemble SVM-RFE
#   - Park et al. (2016) BMC Genomics: ensemble L1-SVM for biomarker
#   - Pes (2019) Neural Computing Appl: ensemble FS stability analysis
#
# and 09B ：
# 09B: types times → (≥2/3)
# 09C: types B times bootstrap → Statisticsselection frequency → 
#        Cross-algorithm aggregation → stability scoreFilter
# ==============================================================================

source(file.path(SCRIPT_DIR, "00_setup.R"))
source(file.path(SCRIPT_DIR, "theme_bindlab.R"))

library(glmnet)
library(randomForest)
library(caret)

# ── 09C using ──
PARAMS_09C <- list(
  n_bootstrap      = 100,     # Bootstrap resampling iterations
  bootstrap_ratio  = 0.8,     # 80% subsample per iteration (no replacement)
 # ↑ Meinshausen & Bühlmann (2010) subsample without replacement
 # ↑ 0.5-0.9 can , 0.8 is using in 
  
 # within Threshold
  lasso_lambda     = "lambda.min",  # Or "lambda.1se" (more conservative)
  rf_ntree         = 500,           # Not too many trees needed per bootstrap
  rf_top_n         = 50,            # RF takes top 50 by MDG each time
  svm_opt_only     = TRUE,          # SVM-RFE takes only optimal features
  
  # Stability SelectionThreshold
  # Meinshausen & Bühlmann (2010): π_thr ∈ [0.6, 0.9]
 # using within frequency ≥ 0.7 Define "Feature"
  freq_threshold   = 0.7,           # Within single algorithm: selected in 70% of bootstraps
  
  # Cross-algorithm aggregation
 # Final: by 2 types is "Feature"
  min_algos        = 2,
  
 # 
  max_top_report   = 50             # Report at most top 50 genes
)

# ── above DependenciesCheck ──
fig_subdir <- "09C_ML_Stability"
dir.create(file.path(FIG_DIR, fig_subdir), showWarnings = FALSE, recursive = TRUE)

pool_file <- file.path(DATA_DIR, "candidate_pool.rds")
if (!file.exists(pool_file)) {
  stop("candidate_pool.rds not found — run 09A_candidate_pool.R first")
}
pool_result    <- readRDS(pool_file)
ml_candidates  <- pool_result$candidate_pool
all_results    <- readRDS(FILES$deg_results)
tpm_filtered   <- readRDS(FILES$tpm_filtered)
gene_anno      <- readRDS(FILES$gene_annotation)

log_step("09C_ML_Stability", sprintf("Starting Bootstrap Stability Selection (B=%d, ratio=%.1f)",
                                      PARAMS_09C$n_bootstrap, PARAMS_09C$bootstrap_ratio))

# ==============================================================================
# 1. ( and 09B)
# ==============================================================================
ml_candidates <- ml_candidates[ml_candidates %in% rownames(tpm_filtered)]

# (>2000)，LRT padj top 2000
if (length(ml_candidates) > 2000) {
  lrt_full <- all_results$lrt_interaction_sig
  id_col <- if ("ensembl_id" %in% colnames(lrt_full)) "ensembl_id" else "ensembl_gene_id"
  top_ids <- lrt_full %>%
    dplyr::filter(.data[[id_col]] %in% ml_candidates) %>%
    dplyr::arrange(padj) %>%
    head(2000) %>%
    dplyr::pull(.data[[id_col]])
  ml_candidates <- top_ids
  log_step("09C_ML_Stability", sprintf("Truncated to top 2000 by padj"))
}

X <- t(log2(tpm_filtered[ml_candidates, ] + 1))
y <- factor(ifelse(sample_info[rownames(X), "Treatment"] == "Induced", "Induced", "Control"),
            levels = c("Control", "Induced"))

n_samples <- nrow(X)
n_features <- ncol(X)
n_sub <- floor(n_samples * PARAMS_09C$bootstrap_ratio)

log_step("09C_ML_Stability", sprintf("Data: %d samples × %d features, subsample size=%d",
                                      n_samples, n_features, n_sub))

# ==============================================================================
# 2. Bootstrap Stability Selection — LASSO
# ==============================================================================
log_step("09C_ML_Stability", sprintf("=== LASSO Bootstrap (B=%d) ===", PARAMS_09C$n_bootstrap))

lasso_selection_matrix <- matrix(0L, nrow = n_features, ncol = PARAMS_09C$n_bootstrap,
                                  dimnames = list(colnames(X), NULL))

set.seed(42)
lasso_errors <- 0
for (b in seq_len(PARAMS_09C$n_bootstrap)) {
  # Subsample without replacement (Meinshausen & Bühlmann 2010)
  idx <- sample(n_samples, n_sub, replace = FALSE)
  X_b <- X[idx, , drop = FALSE]
  y_b <- y[idx]
  
 # Skip type Sample
  if (length(unique(y_b)) < 2) next
  
  tryCatch({
    # LOOCV on subsample
    cv_b <- cv.glmnet(X_b, y_b, family = "binomial", alpha = 1,
                       nfolds = min(length(y_b), 10), type.measure = "class")
    coef_b <- coef(cv_b, s = PARAMS_09C$lasso_lambda)
    nonzero <- rownames(coef_b)[coef_b[, 1] != 0]
    nonzero <- setdiff(nonzero, "(Intercept)")
    
    if (length(nonzero) > 0) {
      lasso_selection_matrix[nonzero, b] <- 1L
    }
  }, error = function(e) {
    lasso_errors <<- lasso_errors + 1
  })
  
  if (b %% 20 == 0) {
    log_step("09C_ML_Stability", sprintf("  LASSO bootstrap %d/%d done", b, PARAMS_09C$n_bootstrap))
  }
}

lasso_freq <- rowMeans(lasso_selection_matrix)
lasso_stable <- names(lasso_freq[lasso_freq >= PARAMS_09C$freq_threshold])

log_step("09C_ML_Stability", sprintf("LASSO: %d/%d bootstraps successful, %d errors",
                                      PARAMS_09C$n_bootstrap - lasso_errors, PARAMS_09C$n_bootstrap, lasso_errors))
log_step("09C_ML_Stability", sprintf("LASSO: %d genes with freq ≥ %.1f (stable), max freq=%.2f",
                                      length(lasso_stable), PARAMS_09C$freq_threshold, max(lasso_freq)))

# ==============================================================================
# 3. Bootstrap Stability Selection — Random Forest
# ==============================================================================
log_step("09C_ML_Stability", sprintf("=== RF Bootstrap (B=%d) ===", PARAMS_09C$n_bootstrap))

rf_selection_matrix <- matrix(0L, nrow = n_features, ncol = PARAMS_09C$n_bootstrap,
                               dimnames = list(colnames(X), NULL))

set.seed(42)
rf_errors <- 0
for (b in seq_len(PARAMS_09C$n_bootstrap)) {
  idx <- sample(n_samples, n_sub, replace = FALSE)
  X_b <- X[idx, , drop = FALSE]
  y_b <- y[idx]
  
  if (length(unique(y_b)) < 2) next
  
  tryCatch({
    rf_b <- randomForest(X_b, y_b, ntree = PARAMS_09C$rf_ntree, importance = TRUE)
    imp_b <- importance(rf_b)
    mdg <- imp_b[, "MeanDecreaseGini"]
    
 # top N by MDG
    top_genes <- names(sort(mdg, decreasing = TRUE))[1:min(PARAMS_09C$rf_top_n, length(mdg))]
 # outside Filter：MDA > 0 ( need to )
    mda <- imp_b[top_genes, "MeanDecreaseAccuracy"]
    top_genes <- top_genes[mda > 0]
    
    if (length(top_genes) > 0) {
      rf_selection_matrix[top_genes, b] <- 1L
    }
  }, error = function(e) {
    rf_errors <<- rf_errors + 1
  })
  
  if (b %% 20 == 0) {
    log_step("09C_ML_Stability", sprintf("  RF bootstrap %d/%d done", b, PARAMS_09C$n_bootstrap))
  }
}

rf_freq <- rowMeans(rf_selection_matrix)
rf_stable <- names(rf_freq[rf_freq >= PARAMS_09C$freq_threshold])

log_step("09C_ML_Stability", sprintf("RF: %d/%d successful, %d errors",
                                      PARAMS_09C$n_bootstrap - rf_errors, PARAMS_09C$n_bootstrap, rf_errors))
log_step("09C_ML_Stability", sprintf("RF: %d genes with freq ≥ %.1f (stable), max freq=%.2f",
                                      length(rf_stable), PARAMS_09C$freq_threshold, max(rf_freq)))

# ==============================================================================
# 4. Bootstrap Stability Selection — SVM-RFE
# ==============================================================================
log_step("09C_ML_Stability", sprintf("=== SVM-RFE Bootstrap (B=%d) ===", PARAMS_09C$n_bootstrap))

# Note:SVM-RFELASSO/RF。Feature>1000Sample, times can can need to 。
# ： using LASSO to top 500SVM-RFE

svm_prefilter <- 500  # SVM-RFE pre-filter: take top genes from LASSO+RF union before each bootstrap
# not is "DEGWGCNA" types —— is pipeline within Calculate,
# (500)FinalSelection

svm_selection_matrix <- matrix(0L, nrow = n_features, ncol = PARAMS_09C$n_bootstrap,
                                dimnames = list(colnames(X), NULL))

set.seed(42)
svm_errors <- 0

# SVM-RFE：LASSO freq + RF freq top genes
combined_freq <- (lasso_freq + rf_freq) / 2
svm_candidate_genes <- names(sort(combined_freq, decreasing = TRUE))[1:min(svm_prefilter, n_features)]

log_step("09C_ML_Stability", sprintf("SVM-RFE: pre-filtered to %d genes (by LASSO+RF freq)",
                                      length(svm_candidate_genes)))

for (b in seq_len(PARAMS_09C$n_bootstrap)) {
  idx <- sample(n_samples, n_sub, replace = FALSE)
  X_b <- X[idx, svm_candidate_genes, drop = FALSE]
  y_b <- y[idx]
  
  if (length(unique(y_b)) < 2) next
  
  tryCatch({
    ctrl_b <- rfeControl(functions = caretFuncs, method = "cv", number = 5, verbose = FALSE)
 # using 5-fold CVLOOCV to bootstrap
    rfe_b <- rfe(
      x = X_b, y = y_b,
      sizes = c(5, 10, 20, 30, 50),
      rfeControl = ctrl_b,
      method = "svmRadial"
    )
    selected_b <- predictors(rfe_b)
    
    if (length(selected_b) > 0) {
      svm_selection_matrix[selected_b, b] <- 1L
    }
  }, error = function(e) {
    svm_errors <<- svm_errors + 1
  })
  
  if (b %% 10 == 0) {
    log_step("09C_ML_Stability", sprintf("  SVM-RFE bootstrap %d/%d done", b, PARAMS_09C$n_bootstrap))
  }
}

svm_freq <- rowMeans(svm_selection_matrix)
svm_stable <- names(svm_freq[svm_freq >= PARAMS_09C$freq_threshold])

log_step("09C_ML_Stability", sprintf("SVM-RFE: %d/%d successful, %d errors",
                                      PARAMS_09C$n_bootstrap - svm_errors, PARAMS_09C$n_bootstrap, svm_errors))
log_step("09C_ML_Stability", sprintf("SVM-RFE: %d genes with freq ≥ %.1f (stable), max freq=%.2f",
                                      length(svm_stable), PARAMS_09C$freq_threshold, max(svm_freq)))

# ==============================================================================
# 5. Cross-algorithm aggregation — Stability Score
# ==============================================================================
log_step("09C_ML_Stability", "=== Cross-algorithm aggregation ===")

stability_df <- data.frame(
  ensembl_id   = colnames(X),
  lasso_freq   = lasso_freq,
  rf_freq      = rf_freq,
  svm_freq     = svm_freq[colnames(X)],   # svm_freq is only a subset of svm_candidate_genes
  stringsAsFactors = FALSE
)
# SVM freq not yet and Gene is 0
stability_df$svm_freq[is.na(stability_df$svm_freq)] <- 0

# types is to ""Threshold
stability_df$lasso_stable <- stability_df$lasso_freq >= PARAMS_09C$freq_threshold
stability_df$rf_stable    <- stability_df$rf_freq    >= PARAMS_09C$freq_threshold
stability_df$svm_stable   <- stability_df$svm_freq   >= PARAMS_09C$freq_threshold

# by types is ""
stability_df$n_stable_algos <- as.integer(stability_df$lasso_stable) +
                                as.integer(stability_df$rf_stable) +
                                as.integer(stability_df$svm_stable)

# Comprehensivestability score = types freq ( types etc.)
stability_df$stability_score <- (stability_df$lasso_freq + 
                                  stability_df$rf_freq + 
                                  stability_df$svm_freq) / 3

# Addgene symbol
stability_df <- stability_df %>%
  dplyr::left_join(gene_anno %>% dplyr::select(ensembl_gene_id, hgnc_symbol),
                   by = c("ensembl_id" = "ensembl_gene_id")) %>%
  dplyr::rename(symbol = hgnc_symbol) %>%
  dplyr::arrange(dplyr::desc(stability_score))

# FinalGene： by 2 types is 
final_stable <- stability_df %>%
  dplyr::filter(n_stable_algos >= PARAMS_09C$min_algos)

log_step("09C_ML_Stability", sprintf(
  "Final stable genes: %d (≥%d algos with freq≥%.1f)",
  nrow(final_stable), PARAMS_09C$min_algos, PARAMS_09C$freq_threshold))
log_step("09C_ML_Stability", sprintf(
  "  3-algo stable: %d, 2-algo stable: %d, 1-algo: %d, 0-algo: %d",
  sum(stability_df$n_stable_algos == 3),
  sum(stability_df$n_stable_algos == 2),
  sum(stability_df$n_stable_algos == 1),
  sum(stability_df$n_stable_algos == 0)))

# ==============================================================================
# 6. and 09B
# ==============================================================================
log_step("09C_ML_Stability", "=== Comparison with 09B (single-run consensus) ===")

ml_09b_file <- FILES$ml_results
if (file.exists(ml_09b_file)) {
  ml_09b <- readRDS(ml_09b_file)
  
  genes_09b_double <- ml_09b$ml_double_genes    # 09B double-algorithm consensus
  genes_09c_stable <- final_stable$ensembl_id    # 09C stable genes
  
  overlap <- intersect(genes_09b_double, genes_09c_stable)
  only_09b <- setdiff(genes_09b_double, genes_09c_stable)
  only_09c <- setdiff(genes_09c_stable, genes_09b_double)
  
  log_step("09C_ML_Stability", sprintf(
    "09B ≥2-consensus: %d genes | 09C stable: %d genes | overlap: %d | only-09B: %d | only-09C: %d",
    length(genes_09b_double), length(genes_09c_stable),
    length(overlap), length(only_09b), length(only_09c)))
  
 # 09BGene in 09C in stability score
  if (length(genes_09b_double) > 0) {
    comparison_df <- stability_df %>%
      dplyr::filter(ensembl_id %in% genes_09b_double) %>%
      dplyr::select(ensembl_id, symbol, lasso_freq, rf_freq, svm_freq, stability_score, n_stable_algos)
    log_step("09C_ML_Stability", "09B genes' stability scores in 09C:")
    print(comparison_df)
  }
} else {
  log_step("09C_ML_Stability", "09B results not found, skipping comparison")
}

# ==============================================================================
# 7. Visualization
# ==============================================================================

# --- 7a. Stability frequency heatmap (top genes) ---
top_n_plot <- min(PARAMS_09C$max_top_report, nrow(stability_df %>% dplyr::filter(stability_score > 0)))
plot_df <- stability_df %>% head(top_n_plot)

if (nrow(plot_df) > 0) {
  row_labels <- ifelse(is.na(plot_df$symbol), substr(plot_df$ensembl_id, 1, 15), plot_df$symbol)
  
  heat_mat <- as.matrix(plot_df[, c("lasso_freq", "rf_freq", "svm_freq")])
  rownames(heat_mat) <- make.unique(row_labels)
  colnames(heat_mat) <- c("LASSO", "Random Forest", "SVM-RFE")
  
 # is 09C stable
  row_anno <- data.frame(
    Stable = factor(ifelse(plot_df$n_stable_algos >= PARAMS_09C$min_algos, "Yes", "No"),
                    levels = c("Yes", "No"))
  )
  rownames(row_anno) <- rownames(heat_mat)
  anno_colors <- list(Stable = c(Yes = "#E64B35", No = "grey85"))
  
  save_heatmap_fig(
    draw_func = function() {
      print(pheatmap(heat_mat,
               color = colorRampPalette(c("grey98", "#FFF2CC", "#F6B73C", "#E64B35"))(50),
               breaks = seq(0, 1, length.out = 51),
               cluster_cols = FALSE, cluster_rows = FALSE,
               fontsize_row = 7, fontsize_col = 10,
               annotation_row = row_anno,
               annotation_colors = anno_colors,
               border_color = "white",
               display_numbers = TRUE,
               number_format = "%.2f",
               number_color = "black",
               main = sprintf("Bootstrap Stability Selection (B=%d, threshold=%.1f)",
                              PARAMS_09C$n_bootstrap, PARAMS_09C$freq_threshold)))
    },
    filename = "09C_stability_heatmap_top",
    subdir = fig_subdir,
    width = 8, height = max(5, nrow(plot_df) * 0.3 + 2)
  )
}

# --- 7b. Stability score distribution ---
p_dist <- ggplot(stability_df %>% dplyr::filter(stability_score > 0),
                 aes(x = stability_score)) +
  geom_histogram(bins = 50, fill = "#3C5488", color = "white", linewidth = 0.2) +
  geom_vline(xintercept = PARAMS_09C$freq_threshold / 3 * 2,
             linetype = 2, color = "#E64B35", linewidth = 0.8) +
  annotate("text", x = PARAMS_09C$freq_threshold / 3 * 2 + 0.02, y = Inf,
           label = sprintf("≥2 algos stable\n(n=%d)", nrow(final_stable)),
           vjust = 1.5, hjust = 0, color = "#E64B35", size = 3.5) +
  labs(title = sprintf("Stability Score Distribution (B=%d)", PARAMS_09C$n_bootstrap),
       subtitle = sprintf("%d genes with score > 0 (out of %d)", 
                          sum(stability_df$stability_score > 0), nrow(stability_df)),
       x = "Stability Score (mean freq across 3 algorithms)",
       y = "Number of genes") +
  theme_bindlab()
save_pub_fig(p_dist, "09C_stability_score_distribution", fig_subdir, width = 8, height = 5)

# --- 7c. Per-algorithm frequency comparison (scatter) ---
if (nrow(final_stable) > 0) {
  p_scatter <- ggplot(stability_df %>% dplyr::filter(lasso_freq > 0 | rf_freq > 0),
                      aes(x = lasso_freq, y = rf_freq)) +
    geom_point(aes(color = factor(n_stable_algos)), alpha = 0.5, size = 1.5) +
    geom_hline(yintercept = PARAMS_09C$freq_threshold, linetype = 2, color = "grey60") +
    geom_vline(xintercept = PARAMS_09C$freq_threshold, linetype = 2, color = "grey60") +
    scale_color_manual(values = c("0" = "grey80", "1" = "#F0B27A", "2" = "#E67E22", "3" = "#E64B35"),
                       name = "N stable\nalgos") +
    ggrepel::geom_text_repel(
      data = final_stable %>% head(15),
      aes(x = lasso_freq, y = rf_freq, label = ifelse(is.na(symbol), ensembl_id, symbol)),
      size = 2.5, max.overlaps = 20, color = "#E64B35") +
    labs(title = "LASSO vs RF Selection Frequency",
         subtitle = sprintf("Dashed lines = stability threshold (%.1f)", PARAMS_09C$freq_threshold),
         x = "LASSO selection frequency", y = "RF selection frequency") +
    theme_bindlab()
  save_pub_fig(p_scatter, "09C_LASSO_vs_RF_frequency", fig_subdir, width = 8, height = 7)
}

# --- 7d. Venn-style: 3Gene overlap ---
stable_list <- list(
  LASSO = lasso_stable,
  RF    = rf_stable,
  SVM   = svm_stable
)
tryCatch({
  p_venn <- ggVennDiagram::ggVennDiagram(stable_list, label_alpha = 0) +
    scale_fill_gradient(low = "white", high = "#E64B35") +
    labs(title = sprintf("Stable Genes per Algorithm (freq ≥ %.1f, B=%d)",
                         PARAMS_09C$freq_threshold, PARAMS_09C$n_bootstrap)) +
    theme(legend.position = "none")
  save_pub_fig(p_venn, "09C_Venn_stable_genes", fig_subdir, width = 7, height = 6)
}, error = function(e) {
  log_step("09C_ML_Stability", "ggVennDiagram not available, skipping Venn plot")
})

# ==============================================================================
# 8. Save
# ==============================================================================
output_09c <- list(
 # 
  params          = PARAMS_09C,
  
 # types 
  lasso_freq      = lasso_freq,
  rf_freq         = rf_freq,
  svm_freq        = svm_freq,
  
 # types Gene
  lasso_stable    = lasso_stable,
  rf_stable       = rf_stable,
  svm_stable      = svm_stable,
  
 # 
  stability_df    = stability_df,
  
 # FinalGene
  final_stable    = final_stable,
  final_genes     = final_stable$ensembl_id,
  
 # SelectionMatrix ( using after )
  lasso_selection_matrix = lasso_selection_matrix,
  rf_selection_matrix    = rf_selection_matrix,
  svm_selection_matrix   = svm_selection_matrix
)

saveRDS(output_09c, file.path(DATA_DIR, "ml_stability_selection.rds"))
write.csv(stability_df, file.path(DATA_DIR, "ML_stability_scores.csv"), row.names = FALSE)
write.csv(final_stable, file.path(DATA_DIR, "ML_stability_final_genes.csv"), row.names = FALSE)

log_step("09C_ML_Stability", sprintf(
  "Step 09C COMPLETE — %d final stable genes saved", nrow(final_stable)))

# ==============================================================================
# 9. need to 
# ==============================================================================
cat("\n")
cat("==============================================================\n")
cat("  09C Bootstrap Stability Selection — SUMMARY\n")
cat("==============================================================\n")
cat(sprintf("  Input: %d genes × %d samples\n", n_features, n_samples))
cat(sprintf("  Bootstrap: B=%d, subsample ratio=%.1f (n_sub=%d)\n",
            PARAMS_09C$n_bootstrap, PARAMS_09C$bootstrap_ratio, n_sub))
cat(sprintf("  Stability threshold: freq ≥ %.1f per algorithm\n", PARAMS_09C$freq_threshold))
cat(sprintf("  Cross-algorithm: ≥ %d algorithms\n", PARAMS_09C$min_algos))
cat("--------------------------------------------------------------\n")
cat(sprintf("  LASSO stable genes:   %d (max freq=%.2f)\n", length(lasso_stable), max(lasso_freq)))
cat(sprintf("  RF stable genes:      %d (max freq=%.2f)\n", length(rf_stable), max(rf_freq)))
cat(sprintf("  SVM-RFE stable genes: %d (max freq=%.2f)\n", length(svm_stable), max(svm_freq)))
cat("--------------------------------------------------------------\n")
cat(sprintf("  FINAL stable (≥%d algos): %d genes\n", PARAMS_09C$min_algos, nrow(final_stable)))
cat(sprintf("    3-algo: %d | 2-algo: %d\n",
            sum(final_stable$n_stable_algos == 3),
            sum(final_stable$n_stable_algos == 2)))
cat("--------------------------------------------------------------\n")
if (nrow(final_stable) > 0) {
  cat("  Top genes by stability score:\n")
  top_show <- head(final_stable, 10)
  for (i in seq_len(nrow(top_show))) {
    cat(sprintf("    %2d. %-12s score=%.3f  L=%.2f R=%.2f S=%.2f  (%d algos)\n",
                i,
                ifelse(is.na(top_show$symbol[i]), top_show$ensembl_id[i], top_show$symbol[i]),
                top_show$stability_score[i],
                top_show$lasso_freq[i],
                top_show$rf_freq[i],
                top_show$svm_freq[i],
                top_show$n_stable_algos[i]))
  }
}
cat("==============================================================\n")
