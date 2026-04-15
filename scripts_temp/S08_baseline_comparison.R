#!/usr/bin/env Rscript
# ==============================================================================
# S08_baseline_comparison.R  v3 — Baseline Method Comparison for METI-FS Paper
#
# Purpose:
#   Compare METI-FS against 4 established feature selection methods on 4 GEO
#   datasets. All alternative methods receive the FULL transcriptome as input;
#   METI-FS uses its multi-evidence candidate pool. This end-to-end design
#   follows the benchmark strategy of Stabl (Hedou et al. 2024 Nat Biotech).
#
# Methods:
#   B1  Single LASSO            Tibshirani 1996, JRSSB
#   B2  Single Elastic Net      Zou & Hastie 2005, JRSSB
#   B3  Stability Selection     Meinshausen & Buhlmann 2010, JRSSB (q-based)
#   B4  Boruta                  Kursa & Rudnicki 2010, JSS
#   B5  METI-FS                 This paper
#
# Key design choices:
#   - B1-B4 operate on full transcriptome (p = 14,000–22,000)
#   - B5 re-runs bootstrap LASSO + gap-union on candidate pool (p = 400–6,200)
#   - Nogueira stability computed with SAME p (full transcriptome) for all
#     methods to ensure comparability
#   - Additionally reports METI-FS actual pipeline n_final (LASSO + SVM union)
#
# Stability Selection note (B3):
#   With n = 17–40 and p > 14,000, standard SS may select 0 genes. This is
#   expected: n/2 subsamples can fit at most ~n/2 LASSO features, and
#   reaching pi_thr = 0.6 across 100 subsamples is extremely unlikely in
#   this regime. We use q-based path selection (targeting ~q features per
#   subsample) to give SS the best chance, but the fundamental small-n
#   limitation persists. This result motivates METI-FS's approach of
#   reducing p before applying stability-based ML selection.
#
# Usage:
#   source("S08_baseline_comparison.R")
#   quick_test()                        # ~2 min
#   results <- run_all_comparisons()    # ~30-60 min (parallel)
#
# Output:  BENCH_DIR/baseline/baseline_comparison_summary.csv
#          BENCH_DIR/baseline/baseline_raw_selections.csv
#          BENCH_DIR/baseline/baseline_paper_table.csv
# ==============================================================================

# ---- 0. Config + Dependencies ----

if (file.exists("S_config.R")) {
  source("S_config.R")
} else if (file.exists(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))) {
  source(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))
}

for (pkg in c("glmnet", "Boruta", "randomForest", "parallel")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[S08] Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(glmnet); library(Boruta); library(randomForest); library(parallel)
})

BASELINE_DIR <- file.path(BENCH_DIR, "baseline")
dir.create(BASELINE_DIR, recursive = TRUE, showWarnings = FALSE)


# ==============================================================================
# 1. Dataset registry + parameters
# ==============================================================================

BASELINE_DATASETS <- list(
  GSE197067 = list(
    project_dir = file.path(METHODS_BASE, "pipeline_runs", "GEO_GSE197067_Tcell"),
    label = "T cell activation", n_samples = 40, n_tp = 5),
  GSE307424 = list(
    project_dir = file.path(METHODS_BASE, "pipeline_runs", "GEO_GSE307424_Lung"),
    label = "SMARCA2 degradation", n_samples = 18, n_tp = 3),
  GSE236646 = list(
    project_dir = file.path(METHODS_BASE, "pipeline_runs/GEO_GSE236646_HSV"),
    label = "HSV-1 infection", n_samples = 17, n_tp = 3),
  GSE150411 = list(
    project_dir = file.path(METHODS_BASE, "pipeline_runs/GEO_GSE150411_Chondro"),
    label = "Chondrocyte FN-f", n_samples = 18, n_tp = 3)
)

# ---- Tuning parameters ----
N_SEEDS          <- 10
SEEDS            <- 42:(42 + N_SEEDS - 1)
N_CORES          <- max(1, detectCores() - 3)   # ~75% utilisation
SS_B             <- 100          # stability selection: number of subsamples
SS_CUTOFF        <- 0.6          # M&B threshold (standard: 0.6–0.9)
SS_Q             <- 20           # target features per subsample (q-based path)
BORUTA_MAX_RUNS  <- 300
BORUTA_NTREE     <- 200
METIFS_B         <- 100          # METI-FS bootstrap iterations
METIFS_SUBSAMPLE <- 0.8          # subsample ratio (matches 09C: 80% without replacement)
METIFS_MIN_FREQ  <- 0.20         # gap-union MIN_FREQ_SIGNAL (matches 09D)

cat(sprintf("[S08] Cores: %d / %d | Output: %s\n",
            N_CORES, detectCores(), BASELINE_DIR))


# ==============================================================================
# 2. Data loading
# ==============================================================================

#' Load full transcriptome TPM matrix + Treatment labels
#' @return list(X, y, gene_names, n, p)
load_transcriptome <- function(project_dir) {

  data_dir <- file.path(project_dir, "data")
  tpm <- readRDS(file.path(data_dir, "tpm_filtered.rds"))   # genes x samples
  si  <- readRDS(file.path(data_dir, "sample_info.rds"))

  # Robust sample matching
  shared <- intersect(colnames(tpm), rownames(si))
  if (!length(shared)) shared <- intersect(colnames(tpm), si$sample_id)
  if (!length(shared) && ncol(tpm) == nrow(si)) shared <- colnames(tpm)
  stopifnot(length(shared) > 0)
  tpm <- tpm[, shared, drop = FALSE]

  trt_col <- intersect(c("Treatment", "treatment"), colnames(si))[1]
  trt <- if (length(shared) == nrow(si)) si[[trt_col]] else
    si[[trt_col]][match(shared, si$sample_id)]
  y <- ifelse(trt == "Induced", 1L, 0L)

  X <- t(as.matrix(tpm))
  X <- log2(X + 1)
  X <- scale(X)
  X[is.na(X)] <- 0

  list(X = X, y = y, gene_names = colnames(X), n = nrow(X), p = ncol(X))
}


#' Load candidate pool gene IDs for METI-FS re-run
load_candidate_pool <- function(project_dir) {
  hub_file <- list.files(file.path(project_dir, "data"),
                         "^hub_genes_final", full.names = TRUE)[1]
  if (is.na(hub_file)) stop("hub_genes_final not found")
  hub <- read.csv(hub_file, stringsAsFactors = FALSE)
  pool <- hub$ensembl_id[hub$in_candidate_pool == TRUE]
  pool <- pool[!is.na(pool) & pool != ""]
  if (!length(pool)) {
    pool <- hub$symbol[hub$in_candidate_pool == TRUE]
    pool <- pool[!is.na(pool) & pool != ""]
  }
  pool
}


#' Load actual METI-FS pipeline output (final candidate genes)
load_metifs_actual <- function(project_dir) {
  f <- list.files(file.path(project_dir, "data"),
                  "^Final_candidate_genes", full.names = TRUE)[1]
  if (is.na(f)) return(character(0))
  df <- read.csv(f, stringsAsFactors = FALSE)
  id <- if ("ensembl_id" %in% names(df)) df$ensembl_id else df$symbol
  id[!is.na(id) & id != ""]
}


# ==============================================================================
# 3. Method implementations
# ==============================================================================

safe_nf <- function(n) max(3L, min(5L, as.integer(n / 4)))


# ---- B1: Single LASSO ----

run_single_lasso <- function(X, y, seed) {
  set.seed(seed)
  tryCatch({
    fit <- cv.glmnet(X, y, family = "binomial", alpha = 1,
                     nfolds = safe_nf(nrow(X)))
    cf <- as.matrix(coef(fit, s = "lambda.min"))
    setdiff(rownames(cf)[abs(cf[, 1]) > 0], "(Intercept)")
  }, error = function(e) character(0))
}


# ---- B2: Single Elastic Net ----

run_single_en <- function(X, y, seed) {
  set.seed(seed)
  tryCatch({
    fit <- cv.glmnet(X, y, family = "binomial", alpha = 0.5,
                     nfolds = safe_nf(nrow(X)))
    cf <- as.matrix(coef(fit, s = "lambda.min"))
    setdiff(rownames(cf)[abs(cf[, 1]) > 0], "(Intercept)")
  }, error = function(e) character(0))
}


# ---- B3: Stability Selection with LASSO (q-based path) ----
#
# Standard SS with cv.glmnet lambda.min fails when n < 20 because each
# n/2-subsample can fit at most ~n/2 features, making it nearly impossible
# for any gene to reach the 0.6 frequency threshold across 100 subsamples.
#
# q-based approach: use the full glmnet regularization path and select the
# lambda yielding approximately q non-zero features per subsample.
# This follows M&B's original formulation where q controls per-subsample
# model complexity.

run_stability_selection <- function(X, y, seed,
                                    B = SS_B, pi_thr = SS_CUTOFF,
                                    q = SS_Q) {
  set.seed(seed)
  n <- nrow(X); p <- ncol(X); half <- floor(n / 2)
  q_eff <- min(q, half - 1)

  freq <- setNames(numeric(p), colnames(X))
  ok <- 0L

  for (b in seq_len(B)) {
    idx <- sample(n, half)
    Xb <- X[idx, , drop = FALSE]; yb <- y[idx]
    if (length(unique(yb)) < 2) next

    tryCatch({
      fit <- glmnet(Xb, yb, family = "binomial", alpha = 1)
      target_idx <- which.min(abs(fit$df - q_eff))
      cf <- as.matrix(coef(fit, s = fit$lambda[target_idx]))
      nz <- setdiff(rownames(cf)[abs(cf[, 1]) > 0], "(Intercept)")
      nz <- nz[nz %in% names(freq)]
      freq[nz] <- freq[nz] + 1
      ok <- ok + 1L
    }, error = function(e) NULL)
  }

  if (ok == 0) return(character(0))
  freq <- freq / ok
  names(freq[freq >= pi_thr])
}


# ---- B4: Boruta ----

run_boruta <- function(X, y, seed,
                       maxRuns = BORUTA_MAX_RUNS, ntree = BORUTA_NTREE) {
  set.seed(seed)
  df <- as.data.frame(X)
  df$.Y <- factor(y, 0:1, c("C", "I"))
  tryCatch({
    bor <- Boruta(.Y ~ ., data = df, maxRuns = maxRuns,
                  num.trees = ntree, doTrace = 0)
    getSelectedAttributes(bor, withTentative = TRUE)
  }, error = function(e) character(0))
}


# ---- B5: METI-FS bootstrap LASSO + gap-union on candidate pool ----
#
# Mirrors 09C + 09D pipeline logic:
#   - 80% subsample without replacement (Meinshausen & Buhlmann 2010)
#   - LASSO with lambda.min
#   - Gap-union data-driven threshold (MIN_FREQ_SIGNAL = 0.20)
#
# This re-run uses LASSO only (the dominant algorithm in 3/4 datasets).
# The actual METI-FS pipeline additionally runs SVM-RFE and RF, with RF
# auto-exclusion. The actual pipeline output (n_actual) is reported
# separately from the re-run stability metrics.

run_metifs_ml <- function(X_pool, y, seed,
                          B = METIFS_B,
                          subsample_ratio = METIFS_SUBSAMPLE,
                          min_freq = METIFS_MIN_FREQ) {
  set.seed(seed)
  n <- nrow(X_pool); p <- ncol(X_pool)
  n_sub <- floor(n * subsample_ratio)

  freq <- setNames(numeric(p), colnames(X_pool))
  ok <- 0L

  for (b in seq_len(B)) {
    # 80% subsample WITHOUT replacement (matches 09C)
    idx <- sample(n, n_sub, replace = FALSE)
    Xb <- X_pool[idx, , drop = FALSE]; yb <- y[idx]
    if (length(unique(yb)) < 2) next

    tryCatch({
      fit <- cv.glmnet(Xb, yb, family = "binomial", alpha = 1,
                       nfolds = safe_nf(n_sub))
      cf <- as.matrix(coef(fit, s = "lambda.min"))
      nz <- setdiff(rownames(cf)[abs(cf[, 1]) > 0], "(Intercept)")
      nz <- nz[nz %in% names(freq)]
      freq[nz] <- freq[nz] + 1
      ok <- ok + 1L
    }, error = function(e) NULL)
  }

  if (ok == 0) return(character(0))
  freq <- freq / ok

  # ---- Gap-union threshold (mirrors 09D) ----
  above <- sort(freq[freq >= min_freq], decreasing = TRUE)
  if (length(above) <= 1) return(names(above))

  gaps <- -diff(above)
  if (all(gaps == 0)) return(names(above))   # all identical freq

  cut_idx <- which.max(gaps)
  thr <- above[cut_idx + 1]
  names(freq[freq >= thr])
}


# ==============================================================================
# 4. Stability metrics
# ==============================================================================

#' Nogueira stability index (Nogueira, Sechidis & Brown 2018, JMLR 18(174))
#'
#' @param sel_list  list of character vectors (selected genes per seed)
#' @param p         total feature universe size — MUST be the same for all
#'                  methods to ensure comparability
compute_nogueira <- function(sel_list, p) {
  M <- length(sel_list); if (M < 2) return(NA_real_)
  all_g <- unique(unlist(sel_list))
  if (!length(all_g)) return(NA_real_)
  if (p < length(all_g)) p <- length(all_g)

  freq <- setNames(numeric(length(all_g)), all_g)
  k_vec <- integer(M)
  for (m in seq_len(M)) {
    s <- sel_list[[m]]; k_vec[m] <- length(s)
    for (g in s) if (g %in% names(freq)) freq[g] <- freq[g] + 1
  }
  phat <- freq / M; k_bar <- mean(k_vec)
  if (k_bar == 0 || k_bar == p) return(NA_real_)

  s2 <- sum(phat * (1 - phat)) / p
  denom <- (k_bar / p) * (1 - k_bar / p)
  if (denom == 0) return(NA_real_)
  1 - (M / (M - 1)) * s2 / denom
}


#' Mean pairwise Jaccard index
compute_mean_jaccard <- function(sel_list) {
  M <- length(sel_list); if (M < 2) return(NA_real_)
  jvec <- numeric(M * (M - 1) / 2); k <- 0L
  for (i in 1:(M - 1)) for (j in (i + 1):M) {
    k <- k + 1L
    u <- length(union(sel_list[[i]], sel_list[[j]]))
    jvec[k] <- if (u == 0) 0 else length(intersect(sel_list[[i]], sel_list[[j]])) / u
  }
  mean(jvec)
}


# ==============================================================================
# 5. Parallel runner
# ==============================================================================

#' Run a single method across seeds, using parallel for slow methods
run_method_parallel <- function(method_name, X, y, seeds,
                                X_pool = NULL, n_cores = N_CORES) {

  fast <- method_name %in% c("SingleLASSO", "SingleEN")
  use_par <- !fast && n_cores > 1 && length(seeds) > 1

  run_one <- function(seed) {
    switch(method_name,
      SingleLASSO  = run_single_lasso(X, y, seed),
      SingleEN     = run_single_en(X, y, seed),
      StabilitySel = run_stability_selection(X, y, seed),
      Boruta       = run_boruta(X, y, seed),
      METIFS       = if (!is.null(X_pool)) run_metifs_ml(X_pool, y, seed)
                     else character(0))
  }

  if (!use_par) return(lapply(seeds, run_one))

  nc <- min(n_cores, length(seeds))
  cl <- makeCluster(nc)
  on.exit(stopCluster(cl), add = TRUE)

  clusterExport(cl, c("X", "y", "X_pool",
                       "run_single_lasso", "run_single_en",
                       "run_stability_selection", "run_boruta",
                       "run_metifs_ml", "safe_nf",
                       "SS_B", "SS_CUTOFF", "SS_Q",
                       "BORUTA_MAX_RUNS", "BORUTA_NTREE",
                       "METIFS_B", "METIFS_SUBSAMPLE", "METIFS_MIN_FREQ"),
                envir = environment())

  clusterEvalQ(cl, suppressPackageStartupMessages({
    library(glmnet); library(Boruta); library(randomForest)
  }))

  parLapply(cl, seeds, run_one)
}


# ==============================================================================
# 6. Single-dataset runner
# ==============================================================================

run_comparison_single <- function(gse_id, seeds = SEEDS) {

  ds <- BASELINE_DATASETS[[gse_id]]
  stopifnot(!is.null(ds))

  cat("\n================================================================\n")
  cat(sprintf("  %s (%s) | %d seeds | %d cores\n",
              gse_id, ds$label, length(seeds), N_CORES))
  cat("================================================================\n")

  # ---- Load ----
  cat("[LOAD] ")
  dat <- load_transcriptome(ds$project_dir)
  p_universe <- dat$p   # SAME p for all Nogueira calculations
  cat(sprintf("%d samples x %d genes (%d/%d Ind/Ctrl)\n",
              dat$n, dat$p, sum(dat$y), sum(!dat$y)))

  pool_genes <- tryCatch(load_candidate_pool(ds$project_dir),
                         error = function(e) character(0))
  pool_in <- intersect(pool_genes, dat$gene_names)
  X_pool <- if (length(pool_in) > 10) dat$X[, pool_in, drop = FALSE] else NULL

  metifs_actual <- load_metifs_actual(ds$project_dir)
  cat(sprintf("       Pool: %d | METI-FS actual: %d genes\n",
              length(pool_in), length(metifs_actual)))

  # ---- Run methods ----
  methods <- c("SingleLASSO", "SingleEN", "StabilitySel", "Boruta", "METIFS")
  results <- list()

  for (meth in methods) {
    t0 <- Sys.time()
    cat(sprintf("[%-13s] ", meth))
    sels <- run_method_parallel(meth, dat$X, dat$y, seeds, X_pool)
    el <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    nv <- sapply(sels, length)
    cat(sprintf("n = %s  (%.0fs, %.1fs/seed)\n",
                paste(nv, collapse = ","), el, el / length(seeds)))
    results[[meth]] <- sels
  }

  # ---- Compute metrics (ALL use p_universe for Nogueira) ----
  cat("\n[SUMMARY]\n")
  summ <- data.frame(
    dataset = character(), method = character(),
    n_mean = numeric(), n_sd = numeric(),
    n_min = integer(), n_max = integer(), n_cv = numeric(),
    nogueira = numeric(), jaccard = numeric(),
    n_actual = integer(),     # actual pipeline output (METI-FS only)
    p_pool = integer(),       # feature space size used by method
    stringsAsFactors = FALSE
  )

  for (meth in methods) {
    nv  <- sapply(results[[meth]], length)
    nog <- compute_nogueira(results[[meth]], p_universe)  # consistent p
    jac <- compute_mean_jaccard(results[[meth]])

    n_act <- if (meth == "METIFS") length(metifs_actual) else NA_integer_
    p_used <- if (meth == "METIFS") length(pool_in) else p_universe

    summ <- rbind(summ, data.frame(
      dataset = gse_id, method = meth,
      n_mean = mean(nv), n_sd = sd(nv),
      n_min = min(nv), n_max = max(nv),
      n_cv = ifelse(mean(nv) > 0, sd(nv) / mean(nv), NA),
      nogueira = nog, jaccard = jac,
      n_actual = n_act, p_pool = p_used,
      stringsAsFactors = FALSE))

    actual_note <- if (!is.na(n_act)) sprintf(" [actual=%d]", n_act) else ""
    cat(sprintf("  %-14s n=%5.1f +/- %4.1f  Nogueira=%6.3f  Jaccard=%5.3f%s\n",
                meth, mean(nv), sd(nv), nog, jac, actual_note))
  }

  # ---- Raw long table ----
  raw <- do.call(rbind, lapply(methods, function(meth)
    do.call(rbind, lapply(seq_along(seeds), function(i) {
      g <- results[[meth]][[i]]
      if (!length(g)) return(NULL)
      data.frame(dataset = gse_id, method = meth,
                 seed = seeds[i], gene = g, stringsAsFactors = FALSE)
    }))))

  list(summary = summ, raw = raw, selections = results)
}


# ==============================================================================
# 7. Full run
# ==============================================================================

run_all_comparisons <- function(seeds = SEEDS) {

  cat("\n########################################################\n")
  cat(sprintf("  METI-FS Baseline Comparison\n"))
  cat(sprintf("  %d datasets x 5 methods x %d seeds | %d cores\n",
              length(BASELINE_DATASETS), length(seeds), N_CORES))
  cat("########################################################\n")

  all_summ <- all_raw <- NULL
  t0 <- Sys.time()

  for (gse in names(BASELINE_DATASETS)) {
    res <- run_comparison_single(gse, seeds)
    all_summ <- rbind(all_summ, res$summary)
    all_raw  <- rbind(all_raw, res$raw)
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  # ---- Save ----
  write.csv(all_summ, file.path(BASELINE_DIR, "baseline_comparison_summary.csv"),
            row.names = FALSE)
  write.csv(all_raw, file.path(BASELINE_DIR, "baseline_raw_selections.csv"),
            row.names = FALSE)

  # Paper table: for METI-FS, show actual pipeline n alongside re-run stability
  pt <- all_summ
  pt$n_display   <- ifelse(!is.na(pt$n_actual),
                           sprintf("%.0f +/- %.0f [%d]", pt$n_mean, pt$n_sd, pt$n_actual),
                           sprintf("%.0f +/- %.0f", pt$n_mean, pt$n_sd))
  pt$nog_display <- sprintf("%.3f", pt$nogueira)
  pt$jac_display <- sprintf("%.3f", pt$jaccard)
  write.csv(pt, file.path(BASELINE_DIR, "baseline_paper_table.csv"),
            row.names = FALSE)

  # ---- Cross-dataset summary ----
  cat("\n\n============================================================\n")
  cat("  CROSS-DATASET SUMMARY (mean across 4 datasets)\n")
  cat("============================================================\n\n")
  cat(sprintf("  %-14s %10s  %8s  %8s  %s\n",
              "Method", "n (mean)", "Nogueira", "Jaccard", "Note"))
  cat("  ", strrep("-", 62), "\n")

  for (meth in unique(all_summ$method)) {
    sub <- all_summ[all_summ$method == meth, ]
    note <- ""
    if (meth == "METIFS") {
      actual_vals <- sub$n_actual[!is.na(sub$n_actual)]
      if (length(actual_vals)) note <- sprintf("actual=%s",
                                                paste(actual_vals, collapse = "/"))
    }
    cat(sprintf("  %-14s %5.1f +/- %3.0f  %8.3f  %8.3f  %s\n",
                meth,
                mean(sub$n_mean), mean(sub$n_sd),
                mean(sub$nogueira, na.rm = TRUE),
                mean(sub$jaccard, na.rm = TRUE),
                note))
  }

  cat(sprintf("\nTotal: %.1f min | Output: %s\n", elapsed, BASELINE_DIR))
  invisible(list(summary = all_summ, raw = all_raw))
}


# ==============================================================================
# 8. Quick test
# ==============================================================================

quick_test <- function() {
  cat("\n=== QUICK TEST (GSE307424, 3 seeds) ===\n")
  res <- run_comparison_single("GSE307424", seeds = 42:44)
  cat("\n=== Test complete ===\n")
  print(res$summary[, c("method", "n_mean", "n_sd", "nogueira", "jaccard", "n_actual")])
  invisible(res)
}


# ==============================================================================
# Entry
# ==============================================================================

if (sys.nframe() == 0) {
  cat("\n================================================================\n")
  cat("  S08_baseline_comparison.R v3\n")
  cat(sprintf("  Cores: %d | Output: %s\n", N_CORES, BASELINE_DIR))
  cat("================================================================\n\n")
  cat("  quick_test()               ~2 min\n")
  cat("  run_all_comparisons()      ~30-60 min\n\n")
  cat("  Methods: LASSO | EN | StabilitySel | Boruta | METI-FS\n")
  cat("  Nogueira p: full transcriptome (consistent across methods)\n")
  cat("================================================================\n")
}
