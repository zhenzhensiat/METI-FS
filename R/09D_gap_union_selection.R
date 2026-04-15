#!/usr/bin/env Rscript
# =============================================================================
# 09D_gap_union_selection.R — Algorithm-internal gap thresholding + union
#
# Replaces: 09D_stability_ranking_BFE.R (deprecated)
#
# Method (literature):
#   1. Bootstrap Stability Selection (Meinshausen & Bühlmann 2010, JRSS-B)
#      → 09C produces per-gene stability frequencies for each algorithm
#   2. Within each algorithm, identify natural frequency gaps as data-driven
#      thresholds (Spooner et al. 2023, BMC Bioinformatics)
#   3. Union across algorithms (not intersection) to preserve complementary
#      information (Feltes et al. 2022, Knowledge-Based Systems; PMC6420823)
#   4. Algorithms with no discernible signal are excluded
#
# Requires: source("00_setup.R") already executed
#
# Input:  DATA_DIR/ml_stability_selection.rds  (09C)
#         DATA_DIR/ML_stability_scores.csv     (09C)
# Output: DATA_DIR/ml_gap_union.rds
#         DATA_DIR/ML_09D_final_genes.csv
#         DATA_DIR/ML_09D_gap_analysis.csv
#         DATA_DIR/ML_09D_selection_summary.txt
#         FIG_DIR/09D_GapUnion/*.pdf|png
# =============================================================================

STEP_NAME  <- "09D_GapUnion"
fig_subdir <- "09D_GapUnion"
dir.create(file.path(FIG_DIR, fig_subdir), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Upstream dependency check
# ---------------------------------------------------------------------------

stab_file  <- file.path(DATA_DIR, "ml_stability_selection.rds")
score_file <- file.path(DATA_DIR, "ML_stability_scores.csv")

check_upstream(STEP_NAME,
  upstream_files = c(
    "09C_stability_rds" = stab_file,
    "09C_stability_csv" = score_file
  ),
  output_files = c(
    file.path(DATA_DIR, "ml_gap_union.rds"),
    file.path(DATA_DIR, "ML_09D_final_genes.csv")
  )
)

log_step(STEP_NAME, "Starting gap-union selection")

# ---------------------------------------------------------------------------
# 2. Load data
# ---------------------------------------------------------------------------

output_09c <- readRDS(stab_file)
scores     <- read.csv(score_file, stringsAsFactors = FALSE)

log_step(STEP_NAME, sprintf("Loaded %d genes from 09C stability scores", nrow(scores)))

required_cols <- c("ensembl_id", "symbol", "lasso_freq", "rf_freq", "svm_freq")
missing_cols  <- setdiff(required_cols, colnames(scores))
if (length(missing_cols) > 0) {
  stop(sprintf("[%s] Missing columns: %s", STEP_NAME, paste(missing_cols, collapse = ", ")))
}

# ---------------------------------------------------------------------------
# 3. Parameters
# ---------------------------------------------------------------------------
# These are method-level parameters, not in PARAMS (which covers pipeline-wide
# settings). They are specific to the gap-union algorithm.

# MIN_FREQ_SIGNAL: below this, an algorithm is declared failed (no signal).
# Meinshausen (2010) suggests freq in [0.6, 0.9] for strong selection, but in
# p>>n settings, we use 0.20 as the floor for *any* detectable signal.
MIN_FREQ_SIGNAL <- 0.20

# MIN_GAP_RATIO: biggest gap must be >= this × max_freq to count as real.
MIN_GAP_RATIO <- 0.15

# TOP_N_SCAN: scan top N nonzero genes for gaps.
TOP_N_SCAN <- 30

# MIN_GENES_PER_ALGO: minimum genes from any algorithm with signal.
MIN_GENES_PER_ALGO <- 1

# EXCLUDE_ALGOS: algorithms to exclude from gap analysis and union.
# RF excluded based on diagnostic evidence from both lineages:
#   - Lineage_A: RF max_freq=0.14 (below MIN_FREQ_SIGNAL, auto-excluded)
#   - Lineage_B:   RF max_freq=0.86 but max_gap=0.07, 30 genes selected (uniform spray)
#   - Root cause: p>>n (2868:28 or 1514:28), RF cannot distinguish signal from noise
#   - Literature: PMC4387916 (RF high-dim non-informative splits),
#                 BMC Bioinformatics 2016 (VIM instability in high-dim small-sample)
# To re-enable RF: change to EXCLUDE_ALGOS <- character(0)
EXCLUDE_ALGOS <- c("RF")

log_step(STEP_NAME, sprintf(
  "Params: MIN_FREQ_SIGNAL=%.2f, MIN_GAP_RATIO=%.2f, TOP_N_SCAN=%d",
  MIN_FREQ_SIGNAL, MIN_GAP_RATIO, TOP_N_SCAN))
if (length(EXCLUDE_ALGOS) > 0) {
  log_step(STEP_NAME, sprintf("Excluded algorithms: %s", paste(EXCLUDE_ALGOS, collapse = ", ")))
}

# ---------------------------------------------------------------------------
# 4. Core: gap analysis for one algorithm
# ---------------------------------------------------------------------------

analyze_algorithm_gaps <- function(freq_vec, gene_symbols, ensembl_ids,
                                   algo_name) {
  
  ord <- order(freq_vec, decreasing = TRUE)
  freq_sorted <- freq_vec[ord]
  sym_sorted  <- gene_symbols[ord]
  ens_sorted  <- ensembl_ids[ord]
  max_freq    <- freq_sorted[1]
  
  result <- list(
    algo_name    = algo_name,
    max_freq     = max_freq,
    has_signal   = FALSE,
    status       = "",
    selected     = data.frame(
      rank = integer(), symbol = character(), ensembl_id = character(),
      freq = numeric(), gap_below = numeric(), stringsAsFactors = FALSE
    ),
    gap_details  = data.frame(),
    cutoff_freq  = NA_real_,
    biggest_gap  = NA_real_,
    gap_position = NA_integer_
  )
  
  # Check 1: max frequency above signal threshold?
  if (max_freq < MIN_FREQ_SIGNAL) {
    result$status <- sprintf(
      "FAILED: max_freq=%.2f < %.2f (uniform/noise).", max_freq, MIN_FREQ_SIGNAL)
    log_step(STEP_NAME, sprintf("  [%s] %s", algo_name, result$status))
    return(result)
  }
  
  # Top N nonzero
  n_scan   <- min(TOP_N_SCAN, length(freq_sorted))
  top_freq <- freq_sorted[1:n_scan]
  top_sym  <- sym_sorted[1:n_scan]
  top_ens  <- ens_sorted[1:n_scan]
  
  nz_mask  <- top_freq > 0
  if (sum(nz_mask) < 2) {
    if (sum(nz_mask) == 1) {
      result$has_signal  <- TRUE
      result$status      <- sprintf("SIGNAL: single gene with freq=%.2f", max_freq)
      result$selected    <- data.frame(
        rank = 1L, symbol = top_sym[1], ensembl_id = top_ens[1],
        freq = top_freq[1], gap_below = NA_real_, stringsAsFactors = FALSE)
      result$cutoff_freq <- top_freq[1]
    } else {
      result$status <- "FAILED: 0 nonzero genes."
    }
    log_step(STEP_NAME, sprintf("  [%s] %s", algo_name, result$status))
    return(result)
  }
  
  nz_freq <- top_freq[nz_mask]
  nz_sym  <- top_sym[nz_mask]
  nz_ens  <- top_ens[nz_mask]
  
  # Consecutive gaps
  gaps <- abs(diff(nz_freq))
  
  gap_df <- data.frame(
    rank       = seq_along(nz_freq),
    symbol     = nz_sym,
    ensembl_id = nz_ens,
    freq       = nz_freq,
    gap_below  = c(gaps, NA_real_),
    stringsAsFactors = FALSE
  )
  result$gap_details <- gap_df
  
  # Find biggest gap
  biggest_idx <- which.max(gaps)
  biggest_val <- gaps[biggest_idx]
  gap_ratio   <- biggest_val / max_freq
  
  result$biggest_gap  <- biggest_val
  result$gap_position <- biggest_idx
  
  # Check 2: gap large enough?
  if (gap_ratio < MIN_GAP_RATIO) {
    fb <- nz_freq >= MIN_FREQ_SIGNAL
    if (sum(fb) >= 1) {
      result$has_signal  <- TRUE
      result$selected    <- gap_df[fb, , drop = FALSE]
      result$cutoff_freq <- min(nz_freq[fb])
      result$status <- sprintf(
        "WEAK_GAP: gap=%.3f (ratio=%.3f<%.3f). Fallback freq>=%.2f -> %d genes.",
        biggest_val, gap_ratio, MIN_GAP_RATIO, MIN_FREQ_SIGNAL, sum(fb))
    } else {
      result$status <- sprintf(
        "FAILED: gap=%.3f (ratio=%.3f), no genes above %.2f.",
        biggest_val, gap_ratio, MIN_FREQ_SIGNAL)
    }
    log_step(STEP_NAME, sprintf("  [%s] %s", algo_name, result$status))
    return(result)
  }
  
  # Select genes above the gap
  n_selected <- max(biggest_idx, MIN_GENES_PER_ALGO)
  result$has_signal  <- TRUE
  result$selected    <- gap_df[1:n_selected, , drop = FALSE]
  result$cutoff_freq <- nz_freq[n_selected]
  result$status <- sprintf(
    "SIGNAL: gap=%.3f at rank %d->%d (%.2f->%.2f). Selected %d genes (freq>=%.2f).",
    biggest_val, biggest_idx, biggest_idx + 1,
    nz_freq[biggest_idx], nz_freq[biggest_idx + 1],
    n_selected, result$cutoff_freq)
  log_step(STEP_NAME, sprintf("  [%s] %s", algo_name, result$status))
  
  return(result)
}

# ---------------------------------------------------------------------------
# 5. Run gap analysis
# ---------------------------------------------------------------------------

log_step(STEP_NAME, "=== Per-algorithm gap analysis ===")

algo_results <- list()

algo_results[["SVM-RFE"]] <- analyze_algorithm_gaps(
  scores$svm_freq, scores$symbol, scores$ensembl_id, "SVM-RFE")

algo_results[["LASSO"]] <- analyze_algorithm_gaps(
  scores$lasso_freq, scores$symbol, scores$ensembl_id, "LASSO")

if ("RF" %in% EXCLUDE_ALGOS) {
  log_step(STEP_NAME, "  [RF] EXCLUDED by EXCLUDE_ALGOS parameter")
  # Still run diagnostic for logging purposes
  rf_freqs <- sort(scores$rf_freq, decreasing = TRUE)
  rf_nonzero <- rf_freqs[rf_freqs > 0]
  log_step(STEP_NAME, "=== RF quality diagnostic (excluded, for reference) ===")
  log_step(STEP_NAME, sprintf("  RF max_freq = %.2f", rf_freqs[1]))
  log_step(STEP_NAME, sprintf("  RF nonzero genes = %d (of %d)", length(rf_nonzero), nrow(scores)))
  log_step(STEP_NAME, sprintf("  RF genes freq>=0.50 = %d", sum(rf_freqs >= 0.50)))
  # Compute max gap for reference
  top_nz <- head(rf_nonzero, TOP_N_SCAN)
  if (length(top_nz) >= 2) {
    rf_gaps <- abs(diff(top_nz))
    log_step(STEP_NAME, sprintf("  RF max_gap = %.3f (at position %d)", max(rf_gaps), which.max(rf_gaps)))
  }
  log_step(STEP_NAME, "=== End RF diagnostic ===")
  algo_results[["RF"]] <- list(algo_name = "RF", has_signal = FALSE,
    status = "EXCLUDED: removed from analysis (p>>n degradation, see EXCLUDE_ALGOS)",
    max_freq = rf_freqs[1], selected = data.frame(), gap_details = data.frame(),
    cutoff_freq = NA_real_, biggest_gap = NA_real_, gap_position = NA_integer_)
} else {
  algo_results[["RF"]] <- analyze_algorithm_gaps(
    scores$rf_freq, scores$symbol, scores$ensembl_id, "RF")

  # ---------------------------------------------------------------------------
  # 5b. RF quality diagnostic (for manual review)
  #     RF is known to degrade in p>>n settings (PMC4387916, BMC Bioinf 2016).
  #     Output diagnostic metrics to help decide if RF should be trusted.
  # ---------------------------------------------------------------------------

rf_res <- algo_results[["RF"]]
rf_freqs <- sort(scores$rf_freq, decreasing = TRUE)
rf_nonzero <- rf_freqs[rf_freqs > 0]
n_rf_nonzero <- length(rf_nonzero)
n_rf_above_50 <- sum(rf_freqs >= 0.50)
n_rf_above_20 <- sum(rf_freqs >= 0.20)

# Meinshausen 2010 PFER diagnostic: E(V) <= q^2 / ((2*pi_thr - 1) * p)
# q = avg number of features selected per bootstrap (approximate from nonzero count)
# p = candidate pool size
p_pool <- nrow(scores)
# Estimate q: in our bootstrap, RF selects features via importance ranking
# A rough proxy: number of genes with freq > 0 indicates how many distinct
# genes RF ever selects across 100 bootstraps
q_est <- min(n_rf_nonzero, p_pool)  # upper bound

log_step(STEP_NAME, "=== RF quality diagnostic ===")
log_step(STEP_NAME, sprintf("  Candidate pool p = %d", p_pool))
log_step(STEP_NAME, sprintf("  RF max_freq = %.2f", rf_res$max_freq))
log_step(STEP_NAME, sprintf("  RF max_gap = %.3f (at position %s)",
  ifelse(is.na(rf_res$biggest_gap), 0, rf_res$biggest_gap),
  ifelse(is.na(rf_res$gap_position), "NA", as.character(rf_res$gap_position))))
log_step(STEP_NAME, sprintf("  RF nonzero genes = %d (of %d)", n_rf_nonzero, p_pool))
log_step(STEP_NAME, sprintf("  RF genes freq>=0.50 = %d", n_rf_above_50))
log_step(STEP_NAME, sprintf("  RF genes freq>=0.20 = %d", n_rf_above_20))
if (rf_res$has_signal) {
  log_step(STEP_NAME, sprintf("  RF selected %d genes via gap thresholding", nrow(rf_res$selected)))
  if (nrow(rf_res$selected) > 15) {
    log_step(STEP_NAME, sprintf(
      "  WARNING: RF selected %d genes (>15). In p>>n=%d:28 (ratio %.0f:1),",
      nrow(rf_res$selected), p_pool, p_pool/28))
    log_step(STEP_NAME, "  this may indicate RF cannot distinguish signal from noise (PMC4387916).")
    log_step(STEP_NAME, "  Consider excluding RF from union if SVM/LASSO have clear gaps.")
    # Meinshausen PFER estimate at pi_thr = min selected freq
    pi_thr_est <- min(rf_res$selected$freq)
    if (pi_thr_est > 0.5) {
      pfer_est <- q_est^2 / ((2 * pi_thr_est - 1) * p_pool)
      log_step(STEP_NAME, sprintf(
        "  Meinshausen PFER estimate: E(V) <= q²/((2π-1)p) = %d²/((2×%.2f-1)×%d) = %.1f",
        q_est, pi_thr_est, p_pool, pfer_est))
      log_step(STEP_NAME, sprintf(
        "  (PFER>>1 means many false positives expected; PFER<1 means well-controlled)"))
    }
  }
} else {
  log_step(STEP_NAME, sprintf("  RF status: %s", rf_res$status))
}
log_step(STEP_NAME, "=== End RF diagnostic ===")
} # end else (RF not excluded)

# ---------------------------------------------------------------------------
# 6. Cross-algorithm union
# ---------------------------------------------------------------------------

log_step(STEP_NAME, "=== Cross-algorithm union ===")

algo_source_list <- list()

for (aname in names(algo_results)) {
  res <- algo_results[[aname]]
  if (res$has_signal && nrow(res$selected) > 0) {
    sel <- res$selected[, c("ensembl_id", "symbol", "freq"), drop = FALSE]
    sel$source_algo <- aname
    algo_source_list[[aname]] <- sel
    log_step(STEP_NAME, sprintf("  [%s] -> %d genes: %s", aname, nrow(sel),
      paste(sprintf("%s(%.2f)", sel$symbol, sel$freq), collapse = ", ")))
  } else {
    log_step(STEP_NAME, sprintf("  [%s] excluded: %s", aname, res$status))
  }
}

if (length(algo_source_list) == 0) {
  log_step(STEP_NAME, "ERROR: No algorithm produced signal.")
  stop("No algorithm produced signal. Check 09C results.")
}

all_selected <- do.call(rbind, algo_source_list)
rownames(all_selected) <- NULL

# ---------------------------------------------------------------------------
# 6b. Determine BFE core vs gap extension tier per algorithm
# ---------------------------------------------------------------------------
# BFE core: freq >= BFE_CORE_THRESHOLD (very high stability, algorithm-agnostic)
# Gap extension: above gap but below BFE core threshold
BFE_CORE_THRESHOLD <- 0.90  # freq >= 0.90 = BFE core (near-perfect bootstrap)

for (aname in names(algo_source_list)) {
  df <- algo_source_list[[aname]]
  df$selection_tier <- ifelse(df$freq >= BFE_CORE_THRESHOLD, "BFE_core", "gap_extension")
  algo_source_list[[aname]] <- df
  n_core <- sum(df$selection_tier == "BFE_core")
  n_ext  <- sum(df$selection_tier == "gap_extension")
  log_step(STEP_NAME, sprintf("  [%s] Tier: %d BFE_core + %d gap_extension",
    aname, n_core, n_ext))
}

# Aggregate: one row per gene, note all contributing algorithms
unique_ids <- unique(all_selected$ensembl_id)

# Rebuild all_selected with tier info
all_selected_tiered <- do.call(rbind, algo_source_list)
rownames(all_selected_tiered) <- NULL

final_rows <- lapply(unique_ids, function(eid) {
  g <- all_selected_tiered[all_selected_tiered$ensembl_id == eid, , drop = FALSE]
  # A gene is BFE_core if ANY contributing algorithm classifies it as core
  best_tier <- if (any(g$selection_tier == "BFE_core")) "BFE_core" else "gap_extension"
  data.frame(
    ensembl_id     = eid,
    symbol         = g$symbol[1],
    source_algos   = paste(sort(unique(g$source_algo)), collapse = " + "),
    n_algos        = length(unique(g$source_algo)),
    max_freq       = max(g$freq),
    freq_detail    = paste(sprintf("%s=%.2f", g$source_algo, g$freq), collapse = "; "),
    selection_tier = best_tier,
    stringsAsFactors = FALSE
  )
})
final_genes_df <- do.call(rbind, final_rows)
final_genes_df <- final_genes_df[order(-final_genes_df$n_algos, -final_genes_df$max_freq), ]
rownames(final_genes_df) <- NULL

# Merge full stability scores
final_genes_df <- merge(
  final_genes_df,
  scores[, c("ensembl_id", "lasso_freq", "rf_freq", "svm_freq", "stability_score")],
  by = "ensembl_id", all.x = TRUE
)
# Re-sort after merge
final_genes_df <- final_genes_df[order(-final_genes_df$n_algos, -final_genes_df$max_freq), ]

log_step(STEP_NAME, sprintf("RESULT: %d genes from %d algorithm(s)",
  nrow(final_genes_df), length(algo_source_list)))
log_step(STEP_NAME, sprintf("Genes: %s",
  paste(final_genes_df$symbol, collapse = ", ")))

# ---------------------------------------------------------------------------
# 7. Save outputs
# ---------------------------------------------------------------------------

# 7a. Full result object
output_09d <- list(
  params = list(
    MIN_FREQ_SIGNAL    = MIN_FREQ_SIGNAL,
    MIN_GAP_RATIO      = MIN_GAP_RATIO,
    TOP_N_SCAN         = TOP_N_SCAN,
    MIN_GENES_PER_ALGO = MIN_GENES_PER_ALGO
  ),
  algo_results   = algo_results,
  final_genes    = final_genes_df,
  final_gene_ids = final_genes_df$ensembl_id,
  method = "gap_thresholding_union",
  literature = c(
    "Meinshausen & Buhlmann 2010 JRSS-B (stability selection)",
    "Spooner et al. 2023 BMC Bioinformatics (data-driven thresholding)",
    "Feltes et al. 2022 Knowledge-Based Systems (hybrid EFS, union)",
    "PMC6420823 (multi-algorithm union for complementary information)"
  )
)

saveRDS(output_09d, file.path(DATA_DIR, "ml_gap_union.rds"))
log_step(STEP_NAME, "Saved RDS -> ml_gap_union.rds")

# 7b. CSV
write.csv(final_genes_df,
          file.path(DATA_DIR, "ML_09D_final_genes.csv"), row.names = FALSE)
log_step(STEP_NAME, "Saved CSV -> ML_09D_final_genes.csv")

# 7c. Gap analysis details
gap_all_list <- lapply(names(algo_results), function(aname) {
  gd <- algo_results[[aname]]$gap_details
  if (is.data.frame(gd) && nrow(gd) > 0) {
    gd$algorithm <- aname
    return(gd)
  }
  return(NULL)
})
gap_all <- do.call(rbind, Filter(Negate(is.null), gap_all_list))

if (!is.null(gap_all) && nrow(gap_all) > 0) {
  write.csv(gap_all, file.path(DATA_DIR, "ML_09D_gap_analysis.csv"), row.names = FALSE)
  log_step(STEP_NAME, "Saved gap analysis -> ML_09D_gap_analysis.csv")
}

# 7d. Summary report
summary_file <- file.path(DATA_DIR, "ML_09D_selection_summary.txt")
sink(summary_file)
cat("=============================================================\n")
cat(sprintf("09D Gap-Union Selection Summary -- %s\n", PARAMS$diff_type))
cat(sprintf("Date: %s\n", Sys.time()))
cat("=============================================================\n\n")

cat("METHOD:\n")
cat("  Bootstrap Stability Selection (Meinshausen & Buhlmann, 2010)\n")
cat("  -> Algorithm-internal gap thresholding (Spooner et al., 2023)\n")
cat("  -> Cross-algorithm union (Feltes et al., 2022; PMC6420823)\n\n")

cat(sprintf("PARAMETERS:\n"))
cat(sprintf("  MIN_FREQ_SIGNAL = %.2f\n", MIN_FREQ_SIGNAL))
cat(sprintf("  MIN_GAP_RATIO   = %.2f\n", MIN_GAP_RATIO))
cat(sprintf("  TOP_N_SCAN      = %d\n\n", TOP_N_SCAN))

cat("PER-ALGORITHM RESULTS:\n")
cat("-------------------------------------------------------------\n")
for (aname in names(algo_results)) {
  res <- algo_results[[aname]]
  cat(sprintf("\n[%s]\n", aname))
  cat(sprintf("  Max frequency: %.2f\n", res$max_freq))
  cat(sprintf("  Status: %s\n", res$status))
  if (res$has_signal && nrow(res$selected) > 0) {
    cat(sprintf("  Selected %d genes:\n", nrow(res$selected)))
    for (i in seq_len(nrow(res$selected))) {
      cat(sprintf("    %d. %s (freq=%.2f)\n", i,
                  res$selected$symbol[i], res$selected$freq[i]))
    }
  }
}

cat("\n\nFINAL UNION:\n")
cat("-------------------------------------------------------------\n")
cat(sprintf("Total: %d genes\n\n", nrow(final_genes_df)))
for (i in seq_len(nrow(final_genes_df))) {
  r <- final_genes_df[i, ]
  cat(sprintf("  %d. %-12s  source=%-15s  tier=%-14s  LASSO=%.2f RF=%.2f SVM=%.2f  score=%.4f\n",
              i, r$symbol, r$source_algos, r$selection_tier,
              r$lasso_freq, r$rf_freq, r$svm_freq, r$stability_score))
}
sink()
log_step(STEP_NAME, "Saved summary -> ML_09D_selection_summary.txt")


# ---------------------------------------------------------------------------
# 8. Diagnostic plots
# ---------------------------------------------------------------------------

log_step(STEP_NAME, "=== Diagnostic plots ===")

suppressPackageStartupMessages(library(ggplot2))

algo_freq_cols <- c("SVM-RFE" = "svm_freq", "LASSO" = "lasso_freq", "RF" = "rf_freq")
algo_colors    <- c("SVM-RFE" = "#2E86AB", "LASSO" = "#27AE60", "RF" = "#888780")

for (aname in names(algo_freq_cols)) {
  fcol <- algo_freq_cols[[aname]]
  res  <- algo_results[[aname]]
  acol <- algo_colors[[aname]]
  
  df <- data.frame(symbol = scores$symbol, freq = scores[[fcol]],
                   stringsAsFactors = FALSE)
  df <- df[df$freq > 0, ]
  df <- df[order(-df$freq), ]
  df <- head(df, 25)
  if (nrow(df) == 0) next
  df$rank <- seq_len(nrow(df))
  
  sel_syms <- if (res$has_signal && nrow(res$selected) > 0) {
    res$selected$symbol
  } else {
    character(0)
  }
  df$selected <- ifelse(df$symbol %in% sel_syms, "Selected", "Not selected")
  
  p <- ggplot(df, aes(x = reorder(symbol, -rank), y = freq, fill = selected)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Selected" = acol, "Not selected" = "#CCCCCC"),
                      name = "") +
    coord_flip() +
    labs(title = sprintf("%s stability frequency (top %d)", aname, nrow(df)),
         subtitle = res$status,
         x = NULL, y = "Bootstrap frequency (100 iterations)") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom",
          panel.grid.major.y = element_blank(),
          plot.subtitle = element_text(size = 9, color = "grey50"))
  
  if (!is.na(res$cutoff_freq)) {
    p <- p + geom_hline(yintercept = res$cutoff_freq - 0.01,
                        linetype = "dashed", color = "#E74C3C", linewidth = 0.6)
  }
  
  fname <- sprintf("09D_%s_freq", gsub("-", "", aname))
  save_pub_fig(p, fname, fig_subdir, width = 8, height = 6)
}

# Final overview
if (nrow(final_genes_df) > 0) {
  p_final <- ggplot(final_genes_df,
    aes(x = reorder(symbol, stability_score),
        y = stability_score, fill = source_algos)) +
    geom_col(width = 0.7) +
    scale_fill_manual(
      values = c("SVM-RFE" = "#2E86AB", "LASSO" = "#27AE60", "RF" = "#888780",
                 "LASSO + SVM-RFE" = "#9B59B6", "RF + SVM-RFE" = "#1ABC9C",
                 "LASSO + RF" = "#F39C12", "LASSO + RF + SVM-RFE" = "#E74C3C"),
      name = "Source") +
    coord_flip() +
    labs(title = "Final ML biomarkers (gap-union)",
         x = NULL, y = "Stability score") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "bottom")
  
  save_pub_fig(p_final, "09D_final_overview", fig_subdir, width = 8, height = 5)
}


# ---------------------------------------------------------------------------
# 9. Done
# ---------------------------------------------------------------------------

log_step(STEP_NAME, "========================================")
log_step(STEP_NAME, sprintf("COMPLETE: %d ML markers -> %s",
  nrow(final_genes_df), paste(final_genes_df$symbol, collapse = ", ")))
log_step(STEP_NAME, sprintf("  BFE_core: %d, gap_extension: %d",
  sum(final_genes_df$selection_tier == "BFE_core"),
  sum(final_genes_df$selection_tier == "gap_extension")))
for (aname in names(algo_results)) {
  res <- algo_results[[aname]]
  n_s <- if (res$has_signal) nrow(res$selected) else 0
  log_step(STEP_NAME, sprintf("  %-8s: %d genes [%s]", aname, n_s,
    ifelse(res$has_signal, "SIGNAL", "FAILED")))
}
log_step(STEP_NAME, "========================================")
