#!/usr/bin/env Rscript
# ==============================================================================
# S08_baseline_comparison.R — Two-Layer Comparison Framework for METI-FS
#
# Purpose:
#   Evaluate METI-FS contributions through two complementary comparisons:
#
#   Layer 2 — Frontend pipeline comparison
#     3 candidate-pool strategies × 2 ML backends × 4 datasets × 50 seeds
#     Shows that maSigPro∩WGCNA (METI-FS front-end) produces a better
#     candidate pool than DEG-only or DEG∩WGCNA (the current common practice,
#     cf. Cai et al. 2022, BMC Biology).
#
#   Layer 3 — Fair ML comparison
#     5 ML methods × 1 pool (METI-FS candidate pool) × 4 datasets × 50 seeds
#     Shows that bootstrap LASSO + gap-union outperforms single-run alternatives
#     when all methods receive the same input.
#
# Layer 2 pipeline definitions:
#   P1  DEG only         LRT interaction significant genes
#   P2  DEG ∩ WGCNA      LRT genes ∩ WGCNA significant-module genes
#   P3  METI-FS pool     maSigPro ∩ WGCNA candidate pool
#
# Layer 2 ML backends:
#   LASSO   (sparse selection philosophy)
#   Boruta  (all-relevant selection philosophy)
#
# Layer 3 methods (all receive P3 as input):
#   M1  Single LASSO
#   M2  Single Elastic Net
#   M3  Stability Selection (M&B 2010)
#   M4  Boruta
#   M5  METI-FS bootstrap LASSO + gap-union
#
# Note: Layer 2 P3-LASSO ≡ Layer 3 M1 (natural bridge between layers)
#
# Stability metrics:
#   Layer 2: Jaccard (primary, no p needed), Nogueira with each pool's own p
#   Layer 3: Jaccard + Nogueira with unified p (= METI-FS pool size)
#
# References:
#   Nogueira, Sechidis & Brown 2018, JMLR 18(174)
#   Meinshausen & Bühlmann 2010, JRSSB
#   Hédou et al. 2024, Nature Biotechnology (Stabl: 50-seed benchmark)
#   Cai et al. 2022, BMC Biology (DEG+WGCNA vs WGCNA+DEG)
#   Bommert et al. 2020, CSDA (fair FS benchmarking)
#
# Usage:
#   source("S08_baseline_comparison.R")
#   quick_test_v5()                     # ~3 min, GSE307424 only, 3 seeds
#   results <- run_all_v5()             # ~2-3 hours (parallel)
#
# Output:
#   BENCH_DIR/baseline_v5/layer2_summary.csv
#   BENCH_DIR/baseline_v5/layer3_summary.csv
#   BENCH_DIR/baseline_v5/layer2_raw.csv
#   BENCH_DIR/baseline_v5/layer3_raw.csv
#   BENCH_DIR/baseline_v5/layer2_paper_table.csv
#   BENCH_DIR/baseline_v5/layer3_paper_table.csv
# ==============================================================================


# ---- 0. Config + Dependencies ------------------------------------------------

if (file.exists("S_config.R")) {
  source("S_config.R")

  
}

for (pkg in c("glmnet", "Boruta", "randomForest", "parallel")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("[S08v5] Installing %s...\n", pkg))
    install.packages(pkg, repos = "https://cran.r-project.org")
  }
}

suppressPackageStartupMessages({
  library(glmnet); library(Boruta); library(randomForest); library(parallel)
})

OUT_DIR <- file.path(BENCH_DIR, "baseline_v5")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)


# ---- 1. Dataset registry + parameters ---------------------------------------

DATASETS_V5 <- list(
  GSE197067 = list(
    project_dir = file.path(RUN_DIR, "GEO_GSE197067_Tcell"),
    label       = "T cell activation",
    n_samples   = 40,
    n_tp        = 5),
  GSE307424 = list(
    project_dir = file.path(RUN_DIR, "GEO_GSE307424_Lung"),
    label       = "SMARCA2 degradation",
    n_samples   = 18,
    n_tp        = 3),
  GSE236646 = list(
    project_dir = file.path(RUN_DIR, "GEO_GSE236646_NPC"),
    label       = "NPC infection",
    n_samples   = 16,
    n_tp        = 3),
  GSE150411 = list(
    project_dir = file.path(RUN_DIR, "GEO_GSE150411_Chon"),
    label       = "Chondrocyte FN-f",
    n_samples   = 18,
    n_tp        = 3)
)

# ---- Tuning parameters ----
N_SEEDS          <- 50           # Hédou 2024 Nat Biotech standard
SEEDS            <- 42:(42 + N_SEEDS - 1)
N_CORES          <- max(1, detectCores() - 3)   # i5-10500: 12 threads → 9 workers
SS_B             <- 100          # Stability Selection: subsamples
SS_CUTOFF        <- 0.6          # M&B threshold
SS_Q             <- 20           # target features per subsample
BORUTA_MAX_RUNS  <- 300
BORUTA_NTREE     <- 200
METIFS_B         <- 100          # METI-FS bootstrap iterations
METIFS_SUBSAMPLE <- 0.8          # 80% subsample without replacement (09C)
METIFS_MIN_FREQ  <- 0.20         # gap-union MIN_FREQ_SIGNAL (09D)

cat(sprintf("[S08v5] Cores: %d / %d | Seeds: %d | Output: %s\n",
            N_CORES, detectCores(), N_SEEDS, OUT_DIR))


# ==============================================================================
# 2. Data loading — three gene pools per dataset
# ==============================================================================

#' Detect gene ID type used in a vector of names
#' @return "ensembl" or "symbol"
detect_id_type <- function(ids) {
  if (any(grepl("^ENSG", ids[1:min(5, length(ids))]))) "ensembl" else "symbol"
}


#' Extract WGCNA significant-module genes (all genes in sig_modules, not just hubs)
#'
#' @param wgcna_signed  The $signed component of wgcna_results.rds
#' @return Character vector of gene IDs (ensembl or symbol, matching datExpr)
extract_wgcna_key_genes <- function(wgcna_signed) {
  sig_mods   <- wgcna_signed$sig_modules              # e.g. "ME7","ME1",...
  mod_nums   <- as.numeric(gsub("ME", "", sig_mods))   # e.g. 7, 1, ...
  mod_labels <- wgcna_signed$moduleLabels              # named numeric vector
  # genes whose module label is in sig_mods
  key_genes  <- names(mod_labels)[mod_labels %in% mod_nums]
  key_genes
}


#' Load all data needed for v5 comparison
#'
#' Returns list with:
#'   X          — scaled log2(TPM+1) matrix, samples × genes
#'   y          — binary treatment label (1=Induced, 0=Control)
#'   pools      — list(P1, P2, P3) of gene ID character vectors
#'   pool_sizes — named integer vector
#'   gene_id_type — "ensembl" or "symbol"
#'   metifs_actual_n — actual pipeline output gene count
load_dataset_v5 <- function(project_dir) {

  data_dir <- file.path(project_dir, "data")
  cat("  [LOAD] ")

  # ---- 2a. Expression matrix & treatment labels (from v4, proven to work) ----
  tpm <- readRDS(file.path(data_dir, "tpm_filtered.rds"))   # genes × samples
  si  <- readRDS(file.path(data_dir, "sample_info.rds"))

  shared <- intersect(colnames(tpm), rownames(si))
  if (!length(shared)) shared <- intersect(colnames(tpm), si$sample_id)
  if (!length(shared) && ncol(tpm) == nrow(si)) shared <- colnames(tpm)
  stopifnot(length(shared) > 0)
  tpm <- tpm[, shared, drop = FALSE]

  trt_col <- intersect(c("Treatment", "treatment"), colnames(si))[1]
  trt <- if (length(shared) == nrow(si)) si[[trt_col]] else
    si[[trt_col]][match(shared, si$sample_id)]
  y <- ifelse(trt == "Induced", 1L, 0L)

  X <- t(as.matrix(tpm))          # samples × genes
  X <- log2(X + 1)
  X <- scale(X)
  X[is.na(X)] <- 0
  all_genes <- colnames(X)
  gene_id_type <- detect_id_type(all_genes)

  cat(sprintf("%d samples × %d genes [%s] (%d/%d Ind/Ctrl)\n",
              nrow(X), ncol(X), gene_id_type, sum(y), sum(!y)))

  # ---- 2b. WGCNA key module genes ----
  wgcna <- readRDS(file.path(data_dir, "wgcna_results.rds"))
  ws    <- wgcna$signed
  wgcna_key <- extract_wgcna_key_genes(ws)
  wgcna_id_type <- detect_id_type(wgcna_key)

  # If WGCNA uses different ID type from expression matrix, convert via hub_genes table
  if (wgcna_id_type != gene_id_type) {
    hg <- ws$hub_genes
    if (gene_id_type == "ensembl" && "symbol" %in% names(hg)) {
      # WGCNA uses symbols, expression uses ensembl → map via hub_genes
      map <- setNames(hg$ensembl_id, hg$symbol)
      wgcna_key <- unname(map[wgcna_key])
      wgcna_key <- wgcna_key[!is.na(wgcna_key)]
    } else if (gene_id_type == "symbol" && "ensembl_id" %in% names(hg)) {
      map <- setNames(hg$symbol, hg$ensembl_id)
      wgcna_key <- unname(map[wgcna_key])
      wgcna_key <- wgcna_key[!is.na(wgcna_key)]
    }
  }
  cat(sprintf("  [WGCNA] %d sig-module genes (ID type: %s)\n",
              length(wgcna_key), wgcna_id_type))

  # ---- 2c. LRT interaction genes (significant only: padj < 0.05) ----
  lrt <- readRDS(file.path(data_dir, "lrt_interaction_results.rds"))
  n_lrt_total <- nrow(lrt)
  if ("padj" %in% names(lrt)) {
    lrt <- lrt[!is.na(lrt$padj) & lrt$padj < 0.05, ]
  } else {
    warning("No padj column in LRT results — using all genes (no filter)")
  }
  lrt_genes <- rownames(lrt)
  lrt_id_type <- detect_id_type(lrt_genes)

  # Convert LRT gene IDs if needed
  if (lrt_id_type != gene_id_type) {
    # Try to use hub_genes_final for mapping
    hgf_file <- list.files(data_dir, "^hub_genes_final", full.names = TRUE)[1]
    if (!is.na(hgf_file)) {
      hgf <- read.csv(hgf_file, stringsAsFactors = FALSE)
      if (gene_id_type == "symbol" && "symbol" %in% names(hgf)) {
        map <- setNames(hgf$symbol, hgf$ensembl_id)
        lrt_genes <- unname(map[lrt_genes])
        lrt_genes <- lrt_genes[!is.na(lrt_genes) & lrt_genes != ""]
      }
    }
  }
  cat(sprintf("  [LRT]   %d / %d interaction genes pass padj<0.05 (ID type: %s)\n",
              length(lrt_genes), n_lrt_total, lrt_id_type))

  # ---- 2d. METI-FS candidate pool ----
  # Try candidate_pool.rds first, then hub_genes_final.csv fallback
  pool_file <- file.path(data_dir, "candidate_pool.rds")
  if (file.exists(pool_file)) {
    cp <- readRDS(pool_file)
    if (is.data.frame(cp) && "ensembl_id" %in% names(cp)) {
      pool_genes <- cp$ensembl_id
    } else if (is.list(cp) && "candidate_pool" %in% names(cp)) {
      pool_genes <- cp$candidate_pool
      if (is.data.frame(pool_genes)) pool_genes <- pool_genes$ensembl_id
    } else {
      pool_genes <- character(0)
    }
  } else {
    pool_genes <- character(0)
  }

  # Fallback: extract from hub_genes_final
  if (!length(pool_genes)) {
    hgf_file <- list.files(data_dir, "^hub_genes_final", full.names = TRUE)[1]
    if (is.na(hgf_file)) {
      # Recursive search (v4 fix for GSE236646)
      hgf_file <- list.files(project_dir, "^hub_genes_final.*\\.csv$",
                              recursive = TRUE, full.names = TRUE)[1]
    }
    if (!is.na(hgf_file)) {
      hgf <- read.csv(hgf_file, stringsAsFactors = FALSE)
      pool_genes <- hgf$ensembl_id[hgf$in_candidate_pool == TRUE]
      pool_genes <- pool_genes[!is.na(pool_genes) & pool_genes != ""]
      # If ensembl_id column has symbols (GSE307424 quirk), try symbol column
      if (!length(pool_genes) || detect_id_type(pool_genes) != gene_id_type) {
        pool_genes_alt <- hgf$symbol[hgf$in_candidate_pool == TRUE]
        pool_genes_alt <- pool_genes_alt[!is.na(pool_genes_alt) & pool_genes_alt != ""]
        if (length(pool_genes_alt) > length(pool_genes)) pool_genes <- pool_genes_alt
      }
    }
  }

  # Convert pool gene IDs if needed
  pool_id_type <- if (length(pool_genes)) detect_id_type(pool_genes) else gene_id_type
  if (pool_id_type != gene_id_type && length(pool_genes)) {
    hgf <- if (exists("hgf")) hgf else {
      f <- list.files(data_dir, "^hub_genes_final", full.names = TRUE)[1]
      if (!is.na(f)) read.csv(f, stringsAsFactors = FALSE) else NULL
    }
    if (!is.null(hgf)) {
      if (gene_id_type == "symbol") {
        map <- setNames(hgf$symbol, hgf$ensembl_id)
      } else {
        map <- setNames(hgf$ensembl_id, hgf$symbol)
      }
      pool_genes <- unname(map[pool_genes])
      pool_genes <- pool_genes[!is.na(pool_genes) & pool_genes != ""]
    }
  }

  # ---- 2e. Construct three pools (intersect with expression matrix) ----
  P1 <- intersect(lrt_genes, all_genes)
  P2 <- intersect(P1, wgcna_key)
  P3 <- intersect(pool_genes, all_genes)

  cat(sprintf("  [POOLS] P1(DEG)=%d  P2(DEG∩WGCNA)=%d  P3(METI-FS)=%d\n",
              length(P1), length(P2), length(P3)))

  # Safety checks
  if (length(P1) < 50) warning(sprintf("P1 very small (%d) — check ID matching", length(P1)))
  if (length(P2) < 50) warning(sprintf("P2 very small (%d) — check ID matching", length(P2)))
  if (length(P3) < 10) warning(sprintf("P3 very small (%d) — check ID matching", length(P3)))
  # P2 must be subset of P1; P3 should largely overlap with P2
  stopifnot(all(P2 %in% P1))

  # ---- 2f. Actual METI-FS output ----
  final_file <- list.files(data_dir, "^Final_candidate_genes", full.names = TRUE)[1]
  if (is.na(final_file)) {
    final_file <- list.files(project_dir, "^Final_candidate_genes.*\\.csv$",
                              recursive = TRUE, full.names = TRUE)[1]
  }
  metifs_actual_n <- 0L
  if (!is.na(final_file)) {
    fc <- read.csv(final_file, stringsAsFactors = FALSE)
    metifs_actual_n <- nrow(fc)
  }

  list(
    X              = X,
    y              = y,
    all_genes      = all_genes,
    pools          = list(P1 = P1, P2 = P2, P3 = P3),
    pool_sizes     = c(P1 = length(P1), P2 = length(P2), P3 = length(P3)),
    gene_id_type   = gene_id_type,
    metifs_actual_n = metifs_actual_n
  )
}


# ==============================================================================
# 3. ML method implementations (reused from v4, proven correct)
# ==============================================================================

safe_nf <- function(n) max(3L, min(5L, as.integer(n / 4)))


# ---- Single LASSO ----
run_single_lasso <- function(X, y, seed) {
  set.seed(seed)
  tryCatch({
    fit <- cv.glmnet(X, y, family = "binomial", alpha = 1,
                     nfolds = safe_nf(nrow(X)))
    cf <- as.matrix(coef(fit, s = "lambda.min"))
    setdiff(rownames(cf)[abs(cf[, 1]) > 0], "(Intercept)")
  }, error = function(e) character(0))
}


# ---- Single Elastic Net (alpha = 0.5) ----
run_single_en <- function(X, y, seed) {
  set.seed(seed)
  tryCatch({
    fit <- cv.glmnet(X, y, family = "binomial", alpha = 0.5,
                     nfolds = safe_nf(nrow(X)))
    cf <- as.matrix(coef(fit, s = "lambda.min"))
    setdiff(rownames(cf)[abs(cf[, 1]) > 0], "(Intercept)")
  }, error = function(e) character(0))
}


# ---- Stability Selection (M&B 2010, q-based path) ----
run_stability_selection <- function(X, y, seed,
                                    B = SS_B, pi_thr = SS_CUTOFF, q = SS_Q) {
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


# ---- Boruta ----
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


# ---- METI-FS bootstrap LASSO + gap-union ----
run_metifs_ml <- function(X, y, seed,
                          B = METIFS_B,
                          subsample_ratio = METIFS_SUBSAMPLE,
                          min_freq = METIFS_MIN_FREQ) {
  set.seed(seed)
  n <- nrow(X); p <- ncol(X)
  n_sub <- floor(n * subsample_ratio)

  freq <- setNames(numeric(p), colnames(X))
  ok <- 0L

  for (b in seq_len(B)) {
    idx <- sample(n, n_sub, replace = FALSE)
    Xb <- X[idx, , drop = FALSE]; yb <- y[idx]
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

  # Gap-union threshold (mirrors 09D)
  above <- sort(freq[freq >= min_freq], decreasing = TRUE)
  if (length(above) <= 1) return(names(above))

  gaps <- -diff(above)
  if (all(gaps == 0)) return(names(above))

  cut_idx <- which.max(gaps)
  thr <- above[cut_idx + 1]
  names(freq[freq >= thr])
}


# ==============================================================================
# 4. Stability metrics
# ==============================================================================

#' Nogueira stability index (Nogueira, Sechidis & Brown 2018, JMLR)
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
# 5. Parallel infrastructure — persistent cluster per dataset
# ==============================================================================

#' Dispatch a single ML method on a given expression subset
run_one_method <- function(method_name, X, y, seed) {
  switch(method_name,
    SingleLASSO   = run_single_lasso(X, y, seed),
    SingleEN      = run_single_en(X, y, seed),
    StabilitySel  = run_stability_selection(X, y, seed),
    Boruta        = run_boruta(X, y, seed),
    METIFS_ML     = run_metifs_ml(X, y, seed),
    stop("Unknown method: ", method_name)
  )
}


#' Create a reusable cluster with all ML functions pre-loaded
#'
#' Cluster is created ONCE per dataset and reused across all methods.
#' Only the expression matrix (X) and labels (y) need updating when
#' the gene pool changes.
create_persistent_cluster <- function(n_cores = N_CORES) {
  nc <- min(n_cores, N_SEEDS)
  cat(sprintf("  [CLUSTER] Creating %d workers... ", nc))
  cl <- makeCluster(nc)

  # Export all ML functions and tuning parameters (one-time cost)
  clusterExport(cl, c("run_one_method",
                       "run_single_lasso", "run_single_en",
                       "run_stability_selection", "run_boruta",
                       "run_metifs_ml", "safe_nf",
                       "SS_B", "SS_CUTOFF", "SS_Q",
                       "BORUTA_MAX_RUNS", "BORUTA_NTREE",
                       "METIFS_B", "METIFS_SUBSAMPLE", "METIFS_MIN_FREQ"),
                envir = globalenv())

  clusterEvalQ(cl, suppressPackageStartupMessages({
    library(glmnet); library(Boruta); library(randomForest)
  }))

  cat("done\n")
  cl
}


#' Update expression matrix and labels on existing cluster
update_cluster_data <- function(cl, X, y) {
  clusterExport(cl, c("X", "y"), envir = environment())
}


#' Run method across seeds on a pre-initialised cluster
#'
#' @param cl        Persistent cluster (NULL = run serial)
#' @param method_name  ML method name
#' @param seeds     Seed vector
run_on_cluster <- function(cl, method_name, seeds) {
  run_fn <- function(seed) run_one_method(method_name, X, y, seed)

  if (is.null(cl)) {
    # Serial fallback
    return(lapply(seeds, function(s) {
      run_one_method(method_name, X, y, s)
    }))
  }

  # Ensure method_name is visible on workers
  mn <- method_name
  clusterExport(cl, "mn", envir = environment())
  parLapply(cl, seeds, function(seed) {
    run_one_method(mn, X, y, seed)
  })
}


# ==============================================================================
# 6. Layer 2 — Frontend pipeline comparison (persistent cluster)
# ==============================================================================

#' Run Layer 2 for a single dataset
#'
#' Creates ONE cluster, reuses for all 6 configs (3 pools × 2 ML backends).
#' @return list(summary, raw)
run_layer2_single <- function(gse_id, dat, seeds = SEEDS) {

  cat(sprintf("\n  --- Layer 2: %s ---\n", gse_id))

  # Create persistent cluster once for this dataset
  cl <- tryCatch(create_persistent_cluster(), error = function(e) {
    cat("  [CLUSTER] Failed, falling back to serial\n"); NULL
  })
  on.exit(if (!is.null(cl)) stopCluster(cl), add = TRUE)

  pipelines   <- c("P1_DEG", "P2_DEG_WGCNA", "P3_METIFS")
  ml_backends <- c("SingleLASSO", "Boruta")

  summ <- NULL
  raw  <- NULL

  for (pipe in pipelines) {
    pool_key <- c(P1_DEG = "P1", P2_DEG_WGCNA = "P2", P3_METIFS = "P3")[[pipe]]
    pool_ids <- dat$pools[[pool_key]]
    p_pool   <- length(pool_ids)

    if (p_pool < 10) {
      cat(sprintf("    [%s] SKIP — pool too small (%d)\n", pipe, p_pool))
      next
    }

    # Subset expression matrix for this pool
    X_sub <- dat$X[, pool_ids, drop = FALSE]

    for (ml in ml_backends) {
      config_label <- sprintf("%s × %s", pipe, ml)
      t0 <- Sys.time()
      cat(sprintf("    [%-30s] p=%5d  ", config_label, p_pool))

      # Update data on cluster workers (cheap: only matrix + label)
      X <- X_sub; y <- dat$y
      if (!is.null(cl)) update_cluster_data(cl, X, y)

      sels <- run_on_cluster(cl, ml, seeds)
      el   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      nv   <- sapply(sels, length)

      nog <- compute_nogueira(sels, p_pool)
      jac <- compute_mean_jaccard(sels)

      cat(sprintf("n=%5.1f±%4.1f  Jac=%.3f  Nog=%.3f  (%.0fs)\n",
                  mean(nv), sd(nv), jac, nog, el))

      summ <- rbind(summ, data.frame(
        dataset    = gse_id,
        pipeline   = pipe,
        ml_backend = ml,
        pool_size  = p_pool,
        n_mean     = mean(nv),
        n_sd       = sd(nv),
        n_min      = min(nv),
        n_max      = max(nv),
        jaccard    = jac,
        nogueira   = nog,
        stringsAsFactors = FALSE
      ))

      # Raw selections
      for (i in seq_along(seeds)) {
        g <- sels[[i]]
        if (length(g)) {
          raw <- rbind(raw, data.frame(
            dataset    = gse_id,
            pipeline   = pipe,
            ml_backend = ml,
            seed       = seeds[i],
            gene       = g,
            stringsAsFactors = FALSE
          ))
        }
      }
    }
  }

  list(summary = summ, raw = raw)
}


# ==============================================================================
# 7. Layer 3 — Fair ML comparison (persistent cluster, all methods on P3)
# ==============================================================================

#' Run Layer 3 for a single dataset
run_layer3_single <- function(gse_id, dat, seeds = SEEDS) {

  cat(sprintf("\n  --- Layer 3: %s ---\n", gse_id))

  pool_ids <- dat$pools$P3
  p_pool   <- length(pool_ids)

  if (p_pool < 10) {
    cat(sprintf("    SKIP — P3 pool too small (%d)\n", p_pool))
    return(list(summary = NULL, raw = NULL))
  }

  X_pool <- dat$X[, pool_ids, drop = FALSE]

  # Create persistent cluster — one cluster for all 5 methods
  cl <- tryCatch(create_persistent_cluster(), error = function(e) {
    cat("  [CLUSTER] Failed, falling back to serial\n"); NULL
  })
  on.exit(if (!is.null(cl)) stopCluster(cl), add = TRUE)

  # Export X and y once (same for all methods in Layer 3)
  X <- X_pool; y <- dat$y
  if (!is.null(cl)) update_cluster_data(cl, X, y)

  methods <- c("SingleLASSO", "SingleEN", "StabilitySel", "Boruta", "METIFS_ML")
  summ <- NULL
  raw  <- NULL

  for (meth in methods) {
    t0 <- Sys.time()
    cat(sprintf("    [%-14s] p=%5d  ", meth, p_pool))

    sels <- run_on_cluster(cl, meth, seeds)
    el   <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
    nv   <- sapply(sels, length)

    nog <- compute_nogueira(sels, p_pool)
    jac <- compute_mean_jaccard(sels)

    actual_note <- if (meth == "METIFS_ML")
      sprintf("  [actual=%d]", dat$metifs_actual_n) else ""

    cat(sprintf("n=%5.1f±%4.1f  Jac=%.3f  Nog=%.3f  (%.0fs)%s\n",
                mean(nv), sd(nv), jac, nog, el, actual_note))

    summ <- rbind(summ, data.frame(
      dataset   = gse_id,
      method    = meth,
      pool_size = p_pool,
      n_mean    = mean(nv),
      n_sd      = sd(nv),
      n_min     = min(nv),
      n_max     = max(nv),
      jaccard   = jac,
      nogueira  = nog,
      n_actual  = if (meth == "METIFS_ML") dat$metifs_actual_n else NA_integer_,
      stringsAsFactors = FALSE
    ))

    for (i in seq_along(seeds)) {
      g <- sels[[i]]
      if (length(g)) {
        raw <- rbind(raw, data.frame(
          dataset = gse_id,
          method  = meth,
          seed    = seeds[i],
          gene    = g,
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  list(summary = summ, raw = raw)
}


# ==============================================================================
# 8. Full orchestrator
# ==============================================================================

run_all_v5 <- function(seeds = SEEDS) {

  cat("\n########################################################\n")
  cat("  S08 v5 — Two-Layer Comparison Framework\n")
  cat(sprintf("  %d datasets | %d seeds | %d cores\n",
              length(DATASETS_V5), length(seeds), N_CORES))
  cat("  Layer 2: 3 pipelines × 2 ML backends\n")
  cat("  Layer 3: 5 ML methods on METI-FS pool\n")
  cat("########################################################\n")

  l2_summ <- l2_raw <- NULL
  l3_summ <- l3_raw <- NULL
  t0 <- Sys.time()

  for (gse in names(DATASETS_V5)) {
    ds <- DATASETS_V5[[gse]]
    cat(sprintf("\n================================================================\n"))
    cat(sprintf("  %s (%s)\n", gse, ds$label))
    cat(sprintf("================================================================\n"))

    dat <- load_dataset_v5(ds$project_dir)

    # Layer 2
    res2 <- run_layer2_single(gse, dat, seeds)
    l2_summ <- rbind(l2_summ, res2$summary)
    l2_raw  <- rbind(l2_raw,  res2$raw)

    # Layer 3
    res3 <- run_layer3_single(gse, dat, seeds)
    l3_summ <- rbind(l3_summ, res3$summary)
    l3_raw  <- rbind(l3_raw,  res3$raw)
  }

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))

  # ---- Save results ----
  write.csv(l2_summ, file.path(OUT_DIR, "layer2_summary.csv"), row.names = FALSE)
  write.csv(l3_summ, file.path(OUT_DIR, "layer3_summary.csv"), row.names = FALSE)

  if (!is.null(l2_raw))
    write.csv(l2_raw, file.path(OUT_DIR, "layer2_raw.csv"), row.names = FALSE)
  if (!is.null(l3_raw))
    write.csv(l3_raw, file.path(OUT_DIR, "layer3_raw.csv"), row.names = FALSE)

  # ---- Paper tables ----
  if (!is.null(l2_summ)) {
    pt2 <- l2_summ
    pt2$n_display   <- sprintf("%.1f ± %.1f", pt2$n_mean, pt2$n_sd)
    pt2$jac_display <- sprintf("%.3f", pt2$jaccard)
    pt2$nog_display <- sprintf("%.3f", pt2$nogueira)
    write.csv(pt2, file.path(OUT_DIR, "layer2_paper_table.csv"), row.names = FALSE)
  }

  if (!is.null(l3_summ)) {
    pt3 <- l3_summ
    pt3$n_display <- ifelse(!is.na(pt3$n_actual),
                            sprintf("%.1f ± %.1f [%d]", pt3$n_mean, pt3$n_sd, pt3$n_actual),
                            sprintf("%.1f ± %.1f", pt3$n_mean, pt3$n_sd))
    pt3$jac_display <- sprintf("%.3f", pt3$jaccard)
    pt3$nog_display <- sprintf("%.3f", pt3$nogueira)
    write.csv(pt3, file.path(OUT_DIR, "layer3_paper_table.csv"), row.names = FALSE)
  }

  # ---- Console summary ----
  cat("\n\n############################################################\n")
  cat("  LAYER 2 SUMMARY (Frontend Pipeline Comparison)\n")
  cat("############################################################\n\n")
  if (!is.null(l2_summ)) {
    cat(sprintf("  %-26s  %6s  %8s  %8s  %8s\n",
                "Config", "p_pool", "n_sel", "Jaccard", "Nogueira"))
    cat("  ", strrep("-", 68), "\n")
    for (i in seq_len(nrow(l2_summ))) {
      r <- l2_summ[i, ]
      cat(sprintf("  %-10s %-8s %-7s %6d  %5.1f±%3.0f  %8.3f  %8.3f\n",
                  r$dataset, r$pipeline, r$ml_backend,
                  r$pool_size, r$n_mean, r$n_sd, r$jaccard, r$nogueira))
    }

    # Cross-dataset average per config
    cat("\n  --- Cross-dataset averages ---\n")
    for (pipe in unique(l2_summ$pipeline)) {
      for (ml in unique(l2_summ$ml_backend)) {
        sub <- l2_summ[l2_summ$pipeline == pipe & l2_summ$ml_backend == ml, ]
        if (nrow(sub) == 0) next
        cat(sprintf("  %-26s  p̄=%5.0f  n̄=%5.1f  Jac=%.3f  Nog=%.3f\n",
                    sprintf("%s × %s", pipe, ml),
                    mean(sub$pool_size), mean(sub$n_mean),
                    mean(sub$jaccard, na.rm = TRUE),
                    mean(sub$nogueira, na.rm = TRUE)))
      }
    }
  }

  cat("\n\n############################################################\n")
  cat("  LAYER 3 SUMMARY (Fair ML Comparison on METI-FS Pool)\n")
  cat("############################################################\n\n")
  if (!is.null(l3_summ)) {
    cat(sprintf("  %-14s  %8s  %8s  %8s\n",
                "Method", "n_sel", "Jaccard", "Nogueira"))
    cat("  ", strrep("-", 50), "\n")
    for (meth in unique(l3_summ$method)) {
      sub <- l3_summ[l3_summ$method == meth, ]
      note <- ""
      if (meth == "METIFS_ML") {
        vals <- sub$n_actual[!is.na(sub$n_actual)]
        if (length(vals)) note <- sprintf("  actual=%s", paste(vals, collapse="/"))
      }
      cat(sprintf("  %-14s  %5.1f±%3.0f  %8.3f  %8.3f%s\n",
                  meth,
                  mean(sub$n_mean), mean(sub$n_sd),
                  mean(sub$jaccard, na.rm = TRUE),
                  mean(sub$nogueira, na.rm = TRUE),
                  note))
    }
  }

  cat(sprintf("\nTotal: %.1f min | Output: %s\n", elapsed, OUT_DIR))
  invisible(list(layer2 = list(summary = l2_summ, raw = l2_raw),
                 layer3 = list(summary = l3_summ, raw = l3_raw)))
}


# ==============================================================================
# 9. Quick test
# ==============================================================================

quick_test_v5 <- function() {
  cat("\n=== QUICK TEST v5 (GSE307424, 3 seeds) ===\n")
  ds  <- DATASETS_V5$GSE307424
  dat <- load_dataset_v5(ds$project_dir)

  test_seeds <- 42:44

  res2 <- run_layer2_single("GSE307424", dat, test_seeds)
  res3 <- run_layer3_single("GSE307424", dat, test_seeds)

  cat("\n--- Layer 2 ---\n")
  if (!is.null(res2$summary))
    print(res2$summary[, c("pipeline", "ml_backend", "pool_size",
                            "n_mean", "jaccard", "nogueira")])

  cat("\n--- Layer 3 ---\n")
  if (!is.null(res3$summary))
    print(res3$summary[, c("method", "pool_size", "n_mean",
                            "jaccard", "nogueira", "n_actual")])

  cat("\n=== Quick test complete ===\n")
  invisible(list(layer2 = res2, layer3 = res3))
}


# ==============================================================================
# Entry point
# ==============================================================================

if (sys.nframe() == 0) {
  cat("\n================================================================\n")
  cat("  S08_baseline_comparison.R\n")
  cat("  Two-Layer Comparison Framework for METI-FS\n")
  cat(sprintf("  Cores: %d | Seeds: %d | Output: %s\n",
              N_CORES, N_SEEDS, OUT_DIR))
  cat("================================================================\n\n")
  cat("  quick_test_v5()          ~3 min  (GSE307424, 3 seeds)\n")
  cat("  run_all_v5()             ~2-3 h  (4 datasets, 50 seeds)\n")
  cat("\n")
  cat("  Layer 2: P1(DEG) / P2(DEG∩WGCNA) / P3(METI-FS) × LASSO+Boruta\n")
  cat("  Layer 3: LASSO / EN / SS / Boruta / METI-FS ML on P3 pool\n")
  cat("================================================================\n")
}
