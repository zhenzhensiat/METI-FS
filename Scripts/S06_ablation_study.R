#!/usr/bin/env Rscript
# ==============================================================================
# S06_ablation_study.R — Pipeline component ablation study
#
# Purpose:
#   Systematically quantify the contribution of each pipeline component.
#   Remove components layer by layer, observe changes in final marker set:
#     - simulation data: precision/recall/FDRchanges
#     - real data: marker count/stability/overlap changes
#
# Ablation configurations (7 configs):
#   A0: FULL         — Full pipeline (baseline)
#   A1: -maSigPro    — Skip maSigPro interaction filtering, pool = WGCNA key genes
#   A2: -WGCNA       — Skip WGCNA, pool = maSigPro filtered genes
#   A3: -EffectSize  — Skip effect size filtering, usemaSigPro∩WGCNA
#   A4: -GapUnion    — Use fixed threshold (freq>=0.7, >=2 algos) instead of gap-union
#   A5: -PPI         — Only ML markers, no PPI hub union
#   A6: -Bootstrap   — Use single ML run instead of 100x bootstrap stability
#
# Implementation strategy:
#   No modification to original pipeline scripts. Each ablation intervenes after 09A:
#   - A1/A2/A3: Modify candidate_pool.rds then re-run09C→09D→09F→10
#   - A4: Replace 09D selection logic
#   - A5: Modify 10_integration union logic
#   - A6: Run single ML instead of bootstrap
#
# Input:
#   - Completed FULL pipeline run (all intermediate files under data/ directory)
#   - For simulation data：data/ground_truth.rds
#
# Output:
#   BENCH_DIR/ablation_{run_id}.rds     — Complete ablation results
#   BENCH_DIR/ablation_{run_id}.csv     — Ablation summary table
#
# Literature basis:
#   - pipeComp (Germain et al. 2020, Genome Biology 21:227):
#     pipeline benchmark framework, LOCO (Leave-One-Component-Out) ablation paradigm
#     "enabling the exploration of combinations of parameters and of the
#      robustness of methods to various changes in other parts of a pipeline"
#   - Mangul et al. 2019, Nature Biotechnology 37:1127-1133:
#     benchmarking best practices, "systematic evaluation of bioinformatics
#     software should include ablation studies"
#   - Spooner et al. 2023, BMC Bioinformatics 24:9:
#     data-driven thresholding vs fixed threshold comparison (A4control basis for ablation)
#   - Nogueira, Sechidis & Brown 2018, JMLR 18(174):1-54:
#     stability as evaluation metric for ablation (citation year correction: 2018 not 2017)
# ==============================================================================

# ---- 0. Dependencies ----
if (file.exists("S_config.R")) {
  source("S_config.R")
} else if (file.exists(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))) {
  source(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))
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
# Ablation configuration registry
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
  )
)


# ==============================================================================
# Core function: run_ablation()
# ==============================================================================

#' Execute complete ablation study on a completed pipeline run
#'
#' @param run_dir Completed FULL pipeline run directory
#' @param run_id  Run identifier
#' @param configs Ablation configs to execute (default: all 7)
#' @return list of ablation results
run_ablation <- function(run_dir, run_id = basename(run_dir),
                          configs = names(ABLATION_CONFIGS)) {

  data_dir <- file.path(run_dir, "data")
  gt_file  <- file.path(data_dir, "ground_truth.rds")
  has_gt   <- file.exists(gt_file)

  methods_log("S06_ABLATION", sprintf(
    "=== Ablation study: %s (%d configs, ground_truth=%s) ===",
    run_id, length(configs), has_gt))

  # Load shared upstream data
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

  # Summary table
  summary_df <- build_ablation_summary(results, has_gt)

  # Save
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

  # Print summary
  print_ablation_summary(summary_df, has_gt)

  return(output)
}


# ==============================================================================
# Upstream data loading
# ==============================================================================

load_upstream_data <- function(data_dir) {

  up <- list()

  # 09A candidate pool (with layer intermediate data)
  pool_file <- file.path(data_dir, "candidate_pool.rds")
  if (file.exists(pool_file)) {
    pool <- readRDS(pool_file)
    up$candidate_pool       <- pool$candidate_pool
    up$masigpro_ensembl     <- pool$masigpro_ensembl_filtered
    up$wgcna_key_genes      <- pool$wgcna_key_genes
    up$n_pool_pre_effect    <- pool$n_pool_pre_effect
    up$n_pool_post_effect   <- pool$n_pool_post_effect
    up$effect_filter_applied <- pool$effect_filter_applied
    # maSigPro full set (before interaction filtering)
    up$masigpro_all         <- pool$masigpro_ensembl_all %||%
                                pool$masigpro_ensembl_filtered
  }

  # 09C stability selection matrices
  stab_file <- file.path(data_dir, "ml_stability_selection.rds")
  if (file.exists(stab_file)) {
    stab <- readRDS(stab_file)
    up$stab <- stab
  }

  # 09D gap-union results (A0 baseline)
  gap_file <- file.path(data_dir, "ml_gap_union.rds")
  if (file.exists(gap_file)) {
    up$gap_union <- readRDS(gap_file)
  }

  # 09F PPI hub results
  ppi_file <- file.path(data_dir, "ppi_hub_selection.rds")
  if (file.exists(ppi_file)) {
    up$ppi <- readRDS(ppi_file)
  }

  # TPM + gene_anno (some ablations need to re-run ML)
  tpm_file <- file.path(data_dir, "tpm_filtered.rds")
  if (file.exists(tpm_file)) up$tpm <- readRDS(tpm_file)

  anno_file <- file.path(data_dir, "gene_annotation.rds")
  if (file.exists(anno_file)) up$gene_anno <- readRDS(anno_file)

  sample_file <- file.path(data_dir, "sample_info.rds")
  if (file.exists(sample_file)) up$sample_info <- readRDS(sample_file)

  # DEG results (for effect size filtering ablation)
  deg_file <- file.path(data_dir, "deg_results.rds")
  if (file.exists(deg_file)) up$deg <- readRDS(deg_file)

  return(up)
}


# ==============================================================================
# Single ablation execution
# ==============================================================================

run_single_ablation <- function(cfg_name, cfg, upstream, data_dir, gt) {

  skip <- cfg$skip

  # ====== Step 1: Build candidate pool ======
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
# Ablation sub-function: candidate pool construction
# ==============================================================================

build_ablated_pool <- function(skip, upstream, cfg = NULL, data_dir = NULL) {

  masigpro_genes <- upstream$masigpro_ensembl
  wgcna_genes    <- upstream$wgcna_key_genes
  full_pool      <- upstream$candidate_pool

  if ("masigpro_interaction" %in% skip) {
    # A1: No maSigPro filtering, only WGCNA key module genes
    pool <- wgcna_genes
  } else if ("wgcna" %in% skip) {
    # A2: No WGCNA, only maSigPro filtered genes
    pool <- masigpro_genes
  } else if ("effect_size" %in% skip) {
    # A3: Use maSigPro∩WGCNA but skip effect size filtering
    pool <- intersect(masigpro_genes, wgcna_genes)
  } else if ("pool_resize" %in% skip && !is.null(cfg) && !is.null(cfg$pool_cap)) {
    # [v2] A7/A8: Change candidate pool size cap
    pool_cap <- cfg$pool_cap
    pool <- full_pool  # Start from full candidate pool
    
    if (length(pool) > pool_cap) {
      # Try truncating by LRT padj
      deg_file <- file.path(data_dir, "deg_results.rds")
      if (!is.null(data_dir) && file.exists(deg_file)) {
        all_results <- readRDS(deg_file)
        lrt_full <- all_results$lrt_interaction_sig
        id_col <- if ("ensembl_id" %in% colnames(lrt_full)) "ensembl_id" else "ensembl_gene_id"
        top_ids <- lrt_full[lrt_full[[id_col]] %in% pool, ]
        top_ids <- top_ids[order(top_ids$padj), ]
        pool <- head(top_ids[[id_col]], pool_cap)
      } else {
        # Fallback: Random truncation (preserving information)
        pool <- sample(pool, pool_cap)
      }
    }
    # If pool itself < pool_cap，no change
  } else {
    # A0/A4/A5/A6: Use full candidate pool
    pool <- full_pool
  }

  pool <- pool[!is.na(pool)]
  return(pool)
}


# ==============================================================================
# Ablation sub-function: ML feature selection
# ==============================================================================

run_ablated_ml <- function(skip, pool, upstream, data_dir) {

  # A6: Single ML (no bootstrap) → use single run result
  if ("bootstrap" %in% skip) {
    return(run_single_ml(pool, upstream))
  }

  # If candidate pool differs from original (A1/A2/A3), need to re-run bootstrap
  # But re-running is expensive (30-60 min each), so use approximation:
  # Extract subset from original stability matrix, recalculate frequencies
  if (!identical(sort(pool), sort(upstream$candidate_pool))) {
    return(recompute_stability_subset(pool, upstream, skip))
  }

  # A4: Pool unchanged, but change 09D selection logic
  if ("gap_union" %in% skip) {
    return(apply_fixed_threshold(upstream))
  }

  # A0/A5: Use original ML results
  ml_ids <- upstream$gap_union$final_gene_ids
  return(list(gene_ids = ml_ids))
}


#' A6ablation: single ML run (no bootstrap)
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


#' When pool changes: extract subset from original selection matrix, recalculate frequencies
recompute_stability_subset <- function(pool, upstream, skip) {

  stab <- upstream$stab
  if (is.null(stab) || !all(c("lasso_selection_matrix", "svm_selection_matrix") %in% names(stab))) {
    # No selection matrix, fall back to filtering original results
    orig_ids <- upstream$gap_union$final_gene_ids
    return(list(gene_ids = intersect(orig_ids, pool)))
  }

  # Extract columns (genes) belonging to new pool from each algorithm's selection matrix
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

  # Simplified union: ≥1algorithms stable to be included
  all_stable <- unique(unlist(algo_stable))
  return(list(gene_ids = all_stable))
}


#' A4ablation: fixed threshold instead of gap-union
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

  # Calculate how many algorithms selected each gene
  all_genes <- unique(unlist(stable_per_algo))
  if (length(all_genes) == 0) return(list(gene_ids = character(0)))

  n_algos <- sapply(all_genes, function(g) {
    sum(sapply(stable_per_algo, function(s) g %in% s))
  })

  selected <- names(n_algos)[n_algos >= min_algos]
  return(list(gene_ids = selected))
}


# ==============================================================================
# Ablation sub-function: PPI hub
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
# Ablation sub-function: final union
# ==============================================================================

build_ablated_final <- function(skip, ml_result, ppi_result) {

  ml_ids  <- ml_result$gene_ids
  ppi_ids <- ppi_result$gene_ids

  final <- union(ml_ids, ppi_ids)
  final <- final[!is.na(final) & final != ""]
  return(final)
}


# ==============================================================================
# Summary and output
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
# Batch ablation (all simulation scenarios)
# ==============================================================================

#' Execute ablation study on all simulation scenarios
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

  # Merge all ablation summaries
  all_summaries <- do.call(rbind, lapply(names(all_results), function(nm) {
    df <- all_results[[nm]]$summary
    df$run_id <- nm
    # Parse scenario parameters from run_id
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
# Usage instructions when running directly
# ==============================================================================

`%||%` <- function(x, y) if (is.null(x)) y else x

if (sys.nframe() == 0) {
  cat("\n")
  cat("================================================================\n")
  cat("  S06_ablation_study.R — Usage\n")
  cat("================================================================\n")
  cat("\n")
  cat("  # Single dataset ablation:\n")
  cat("  res <- run_ablation('METHODS_BASE')\n")
  cat("\n")
  cat("  # Batch simulation data ablation:\n")
  cat("  ablate_all_simulations()\n")
  cat("\n")
  cat("  # View ablation configuration:\n")
  cat("  str(ABLATION_CONFIGS)\n")
  cat("================================================================\n")
}
