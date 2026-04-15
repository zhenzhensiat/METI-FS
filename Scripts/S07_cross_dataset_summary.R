#!/usr/bin/env Rscript
# ==============================================================================
# S07_cross_dataset_summary.R — Cross-dataset summary, Nogueira stability, methods paper figures
#
# Functions:
#   1. Nogueira stability metric calculation（Nogueira et al. 2017, JMLR）
#   2. Cross-dataset pipeline behavior consistency statistics
#   3. Simulation benchmark visualization(precision/recall vs SNR/sample size/density)
#   4. Ablation study visualization（ΔF1 barplot）
#   5. Real data application summary（Funnel comparison, algorithm complementarity）
#   6. Generate methods paper'sFigure 2-5 and Table 1-3
#
# Input:
#   BENCH_DIR/benchmark_master.csv       — S05 summary
#   BENCH_DIR/benchmark_results_*.rds    — S05 per-run details
#   BENCH_DIR/ablation_*.csv             — S06 ablation summary
#   BENCH_DIR/ablation_all_simulations.csv — S06 batch ablation
#
# Output:
#   FIG_DIR_METHODS/Fig02_simulation_benchmark.pdf
#   FIG_DIR_METHODS/Fig03_ablation_heatmap.pdf
#   FIG_DIR_METHODS/Fig04_stability_comparison.pdf
#   FIG_DIR_METHODS/Fig05_real_data_summary.pdf
#   TAB_DIR_METHODS/Table01_benchmark_summary.csv
#   TAB_DIR_METHODS/Table02_ablation_summary.csv
#   TAB_DIR_METHODS/Table03_cross_dataset_consistency.csv
#
# Literature basis:
#   - Nogueira, Sechidis & Brown 2018, JMLR 18(174):1-54: stability index
#     (submitted 9/2017, published 4/2018; using stabm R package standard implementation)
#   - stabm R package: Bommert et al. 2021, JOSS 6(59):3010
#   - Meinshausen & Bühlmann 2010, JRSS-B: PFER bound
#   - pipeComp: Germain et al. 2020, Genome Biology 21:227
#     (pipeline benchmark framework, LOCO ablation paradigm)
#   - Spooner et al. 2023, BMC Bioinformatics 24:9:
#     "biggest gap"method is intuitive but may lead to extreme subsets
#     We mitigate this through MIN_GAP_RATIO and MIN_FREQ_SIGNAL parameter constraints
#   - Stabl: Hédou et al. 2024, Nature Biotechnology 42:1581-1593
#     noise injection + data-driven θ threshold (different from gap-union path)
#     Discussionqualitative comparison in Discussion, no head-to-head experiment
# ==============================================================================

# ---- 0. Dependencies ----
if (file.exists("S_config.R")) {
  source("S_config.R")
} else if (file.exists(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))) {
  source(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))
}

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(cowplot)
  library(RColorBrewer)
})

# Theme loading
tryCatch({
  source(file.path(PIPELINE_SCRIPTS, "theme_bindlab.R"))
}, error = function(e) {
  # fallback: simplified theme
  theme_bindlab <- function(base_size = 12) {
    theme_classic(base_size = base_size) %+replace%
      theme(
        plot.title = element_text(size = base_size + 2, face = "bold", hjust = 0),
        axis.title = element_text(face = "bold"),
        legend.key = element_blank()
      )
  }
  assign("theme_bindlab", theme_bindlab, envir = .GlobalEnv)
})

# Methods paper color scheme
METHODS_COLORS <- list(
  snr       = c("low" = "#3C5488", "medium" = "#00A087", "high" = "#E64B35"),
  sample    = c("small" = "#F39B7F", "medium" = "#4DBBD5", "large" = "#3C5488"),
  ablation  = c("Full pipeline" = "#00A087",
                "-maSigPro interaction" = "#E64B35",
                "-WGCNA" = "#3C5488",
                "-Effect size" = "#F39B7F",
                "-Gap-union" = "#B09C85",
                "-PPI hub" = "#7E6148",
                "-Bootstrap stability" = "#DC9FB4"),
  algo      = c("LASSO" = "#E64B35", "RF" = "#4DBBD5", "SVM-RFE" = "#00A087"),
  dataset   = c("GSE307424" = "#3C5488", "GSE197067" = "#E64B35", "Reference" = "#00A087")
)

# Figure saving utility
save_methods_fig <- function(plot_obj, filename,
                              width = 8, height = 6, dpi = 300) {
  for (ext in c("pdf", "png")) {
    fpath <- file.path(FIG_DIR_METHODS, paste0(filename, ".", ext))
    if (ext == "pdf") {
      ggsave(fpath, plot = plot_obj, width = width, height = height,
             device = if (.Platform$OS.type == "windows") cairo_pdf else "pdf")
    } else {
      ggsave(fpath, plot = plot_obj, width = width, height = height,
             dpi = dpi, device = "png")
    }
  }
  methods_log("S07_FIG", sprintf("Saved: %s (.pdf + .png)", filename))
}


# ==============================================================================
# PART 1: NogueiraStability metrics
# ==============================================================================

# Using stabm package standard implementation（Bommert et al. 2021, JOSS 6(59):3010）
# Based on Nogueira, Sechidis & Brown (2018) JMLR 18(174):1-54
# Formula: φ = 1 - (1/p × Σ_f s²_f) / (k̄/p × (1 - k̄/p))
#   where s²_f = M/(M-1) × p̂_f(1 - p̂_f) is the unbiased variance estimate of the f-th feature selection frequency
#   k̄ = average number of selected features, p = total number of features, M = number of bootstrap iterations
# Properties（Nogueira et al. 20185 necessary properties proven by）:
#   1. Fully defined: defined for any M feature subset combinations
#   2. Strict monotonicity: lower selection variance means higher stability
#   3. Bounds: range [-1, 1], lower bound asymptotically approaches 0 as M -> infinity
#   4. Maximum: reaches maximum of 1 when all subsets are identical
#   5. Correction for chance: expected value is 0 for random selection

install_if_missing_stabm <- function() {
  if (!requireNamespace("stabm", quietly = TRUE)) {
    install.packages("stabm", repos = "https://cloud.r-project.org")
  }
}

#' CalculateNogueira stability index
#'
#' Wrapper for stabm::stabilityNogueira(), accepts M×p binary selection matrix
#' and converts to stabm required list-of-feature-indices format.
#'
#' @param sel_matrix M×p binary matrix (rows=bootstrap, cols=features)
#' @return numeric stability index (via stabm)
nogueira_stability <- function(sel_matrix) {

  if (is.null(sel_matrix) || nrow(sel_matrix) < 2) return(NA_real_)

  install_if_missing_stabm()

  M <- nrow(sel_matrix)
  p <- ncol(sel_matrix)

  # Degenerate case: all bootstraps selected no genes, or selected all
  k_bar <- mean(rowSums(sel_matrix))
  if (k_bar == 0 || k_bar == p) return(0)

  # stabmrequires list of feature index vectors
  feature_lists <- lapply(seq_len(M), function(i) {
    which(sel_matrix[i, ] == 1)
  })

  # Filter out empty sets (stabm does not accept empty sets)
  non_empty <- feature_lists[sapply(feature_lists, length) > 0]
  if (length(non_empty) < 2) return(NA_real_)

  phi <- stabm::stabilityNogueira(features = non_empty, p = p)
  return(round(phi, 6))
}


#' Calculate Nogueira metric for all algorithms and gap-union
#'
#' @param run_dir Pipeline run directory
#' @return data.frame with method, nogueira_index, k_bar, p
compute_nogueira_for_run <- function(run_dir) {

  data_dir <- file.path(run_dir, "data")
  stab_file <- file.path(data_dir, "ml_stability_selection.rds")
  if (!file.exists(stab_file)) return(NULL)

  stab <- readRDS(stab_file)
  results <- list()

  for (algo in c("lasso", "rf", "svm")) {
    mat_name <- paste0(algo, "_selection_matrix")
    if (!mat_name %in% names(stab)) next

    mat <- stab[[mat_name]]
    phi <- nogueira_stability(mat)
    k_bar <- mean(rowSums(mat))

    results[[algo]] <- data.frame(
      method   = toupper(algo),
      nogueira = phi,
      k_bar    = round(k_bar, 2),
      p        = ncol(mat),
      M        = nrow(mat),
      stringsAsFactors = FALSE
    )
  }

  # Gap-union result stability: construct "pseudo selection matrix" from 09D BFE_core gene frequencies
  # This is an approximation: the true gap-union applies threshold after each bootstrap
  gap_file <- file.path(data_dir, "ml_gap_union.rds")
  if (file.exists(gap_file)) {
    gap_data <- readRDS(gap_file)
    if (!is.null(gap_data$final_genes)) {
      # Extract final genes columns from each algorithm's selection matrix, calculate joint stability
      final_ids <- gap_data$final_gene_ids
      union_rows <- list()
      for (algo in c("lasso", "svm")) {
        mat_name <- paste0(algo, "_selection_matrix")
        if (!mat_name %in% names(stab)) next
        mat <- stab[[mat_name]]
        common <- intersect(final_ids, colnames(mat))
        if (length(common) > 0) {
          union_rows[[algo]] <- mat[, common, drop = FALSE]
        }
      }
      if (length(union_rows) > 0) {
        # For each bootstrap, mark as 1 if the gene was selected by any algorithm
        M <- nrow(stab$lasso_selection_matrix)
        all_final_genes <- unique(unlist(lapply(union_rows, colnames)))
        combined <- matrix(0L, nrow = M, ncol = length(all_final_genes))
        colnames(combined) <- all_final_genes
        for (algo_mat in union_rows) {
          for (g in colnames(algo_mat)) {
            combined[, g] <- pmax(combined[, g], algo_mat[, g])
          }
        }
        phi_union <- nogueira_stability(combined)
        results[["gap_union"]] <- data.frame(
          method   = "Gap-Union",
          nogueira = phi_union,
          k_bar    = round(mean(rowSums(combined)), 2),
          p        = ncol(combined),
          M        = M,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(results) == 0) return(NULL)
  do.call(rbind, results)
}


# ==============================================================================
# PART 2: Simulation benchmark visualization
# ==============================================================================

#' Fig.2: Precision/Recall/F1 vs SNR × sample_size × marker_density
plot_simulation_benchmark <- function(master_file = file.path(BENCH_DIR, "benchmark_master.csv")) {

  if (!file.exists(master_file)) {
    methods_log("S07_FIG", "benchmark_master.csv not found, skipping Fig.2")
    return(invisible(NULL))
  }

  df <- read.csv(master_file, stringsAsFactors = FALSE)
  sim_df <- df[df$mode == "simulation" & !is.na(df$precision), ]

  if (nrow(sim_df) == 0) {
    methods_log("S07_FIG", "No simulation results with performance data")
    return(invisible(NULL))
  }

  # Parse scenario parameters from run_id
  parts <- strsplit(sim_df$run_id, "_")
  sim_df$snr     <- sapply(parts, `[`, 1)
  sim_df$sample  <- sapply(parts, `[`, 2)
  sim_df$density <- sapply(parts, `[`, 3)

  sim_df$snr     <- factor(sim_df$snr, levels = c("low", "medium", "high"))
  sim_df$sample  <- factor(sim_df$sample, levels = c("small", "medium", "large"))
  sim_df$density <- factor(sim_df$density, levels = c("sparse", "medium", "dense"))

  # Summary: mean ± sd per scenario (across 5 repeats)
  agg <- sim_df %>%
    group_by(snr, sample, density) %>%
    summarise(
      precision_mean = mean(precision, na.rm = TRUE),
      precision_sd   = sd(precision, na.rm = TRUE),
      recall_mean    = mean(recall, na.rm = TRUE),
      recall_sd      = sd(recall, na.rm = TRUE),
      F1_mean        = mean(F1, na.rm = TRUE),
      F1_sd          = sd(F1, na.rm = TRUE),
      FDR_mean       = mean(FDR, na.rm = TRUE),
      FDR_sd         = sd(FDR, na.rm = TRUE),
      n_final_mean   = mean(n_final, na.rm = TRUE),
      n_runs         = n(),
      .groups = "drop"
    )

  # Panel A: F1 by SNR × sample_size (faceted by density)
  p_f1 <- ggplot(agg, aes(x = snr, y = F1_mean, fill = sample)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_errorbar(aes(ymin = pmax(0, F1_mean - F1_sd),
                      ymax = pmin(1, F1_mean + F1_sd)),
                  position = position_dodge(0.8), width = 0.2, linewidth = 0.4) +
    facet_wrap(~ density, labeller = labeller(
      density = c("sparse" = "Sparse markers", "medium" = "Medium", "dense" = "Dense")
    )) +
    scale_fill_manual(values = METHODS_COLORS$sample,
                      labels = c("small" = "16 samples", "medium" = "28 samples",
                                 "large" = "40 samples")) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "METI-FS Performance Across Simulation Scenarios",
         subtitle = "F1 score (mean ± SD across 5 repeats)",
         x = "Signal-to-noise ratio", y = "F1 Score", fill = "Sample size") +
    theme_bindlab() +
    theme(legend.position = "bottom")

  save_methods_fig(p_f1, "Fig02_simulation_benchmark", width = 10, height = 5)

  # Panel B: Precision vs Recall scatter
  p_pr <- ggplot(sim_df, aes(x = recall, y = precision, color = snr, shape = sample)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_manual(values = METHODS_COLORS$snr) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
    scale_x_continuous(limits = c(0, 1)) +
    scale_y_continuous(limits = c(0, 1)) +
    labs(title = "Precision-Recall Trade-off",
         x = "Recall", y = "Precision",
         color = "SNR", shape = "Sample size") +
    theme_bindlab()

  save_methods_fig(p_pr, "Fig02B_precision_recall_scatter", width = 7, height = 6)

  return(invisible(list(f1_plot = p_f1, pr_plot = p_pr, agg = agg)))
}


# ==============================================================================
# PART 3: Ablation study visualization
# ==============================================================================

#' Fig.3: Ablation studyΔF1 barplot
plot_ablation_results <- function(
    ablation_file = file.path(BENCH_DIR, "ablation_all_simulations.csv")) {

  if (!file.exists(ablation_file)) {
    methods_log("S07_FIG", "ablation_all_simulations.csv not found, skipping Fig.3")
    return(invisible(NULL))
  }

  df <- read.csv(ablation_file, stringsAsFactors = FALSE)
  if (!"F1" %in% colnames(df)) {
    methods_log("S07_FIG", "No F1 column in ablation data (real data mode?)")
    return(invisible(NULL))
  }

  # Calculate ΔF1 for each scenario (relative to A0_FULL)
  df_delta <- df %>%
    group_by(run_id) %>%
    mutate(
      F1_baseline = F1[config == "A0_FULL"],
      delta_F1 = F1 - F1_baseline
    ) %>%
    ungroup() %>%
    filter(config != "A0_FULL")

  # Summary: mean ΔF1 per ablation config
  agg <- df_delta %>%
    group_by(config, label) %>%
    summarise(
      delta_F1_mean = mean(delta_F1, na.rm = TRUE),
      delta_F1_sd   = sd(delta_F1, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    arrange(delta_F1_mean)

  agg$label <- factor(agg$label, levels = agg$label)

  p_abl <- ggplot(agg, aes(x = label, y = delta_F1_mean, fill = label)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_errorbar(aes(ymin = delta_F1_mean - delta_F1_sd,
                      ymax = delta_F1_mean + delta_F1_sd),
                  width = 0.2, linewidth = 0.4) +
    geom_hline(yintercept = 0, linewidth = 0.5) +
    scale_fill_manual(values = METHODS_COLORS$ablation, guide = "none") +
    coord_flip() +
    labs(title = "Component Ablation: Impact on F1 Score",
         subtitle = "ΔF1 relative to full pipeline (mean ± SD across scenarios)",
         x = NULL, y = "ΔF1 Score") +
    theme_bindlab()

  save_methods_fig(p_abl, "Fig03_ablation_barplot", width = 8, height = 5)

  # Fig.3B: Ablation heatmap (SNR × ablation)
  if (all(c("snr", "sample_size") %in% colnames(df_delta))) {
    heatmap_data <- df_delta %>%
      group_by(label, snr) %>%
      summarise(delta_F1_mean = mean(delta_F1, na.rm = TRUE), .groups = "drop")

    p_heat <- ggplot(heatmap_data, aes(x = snr, y = label, fill = delta_F1_mean)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = sprintf("%.2f", delta_F1_mean)), size = 3.5) +
      scale_fill_gradient2(low = "#E64B35", mid = "white", high = "#00A087",
                           midpoint = 0, name = "ΔF1") +
      labs(title = "Ablation Impact by Signal-to-Noise Ratio",
           x = "SNR", y = NULL) +
      theme_bindlab() +
      theme(panel.grid = element_blank(),
            axis.text.y = element_text(size = 10))

    save_methods_fig(p_heat, "Fig03B_ablation_heatmap", width = 7, height = 5)
  }

  return(invisible(agg))
}


# ==============================================================================
# PART 4: Stability comparison visualization
# ==============================================================================

#' Fig.4: Gap-union vs fixed thresholdstability comparison + Nogueira metric
plot_stability_comparison <- function(run_dirs = NULL) {

  if (is.null(run_dirs)) {
    # Auto-discover completed runs
    run_dirs <- c(
      file.path(RUN_DIR, "GEO_GSE307424_Lung"),
      file.path(RUN_DIR, "GEO_GSE197067_Tcell")
    )
    run_dirs <- run_dirs[dir.exists(run_dirs)]
  }

  if (length(run_dirs) == 0) {
    methods_log("S07_FIG", "No run directories found, skipping Fig.4")
    return(invisible(NULL))
  }

  # Collect Nogueira metrics
  all_nogueira <- list()
  for (rd in run_dirs) {
    run_id <- basename(rd)
    ng <- compute_nogueira_for_run(rd)
    if (!is.null(ng)) {
      ng$run_id <- run_id
      all_nogueira[[run_id]] <- ng
    }
  }

  if (length(all_nogueira) == 0) {
    methods_log("S07_FIG", "No Nogueira data computed")
    return(invisible(NULL))
  }

  ng_df <- do.call(rbind, all_nogueira)

  # Save table
  write.csv(ng_df, file.path(TAB_DIR_METHODS, "nogueira_stability_index.csv"),
            row.names = FALSE)

  # Fig.4: Nogueira stability by method and dataset
  p_nog <- ggplot(ng_df, aes(x = method, y = nogueira, fill = run_id)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%.3f", nogueira)),
              position = position_dodge(0.8), vjust = -0.3, size = 3) +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(title = "Feature Selection Stability (Nogueira Index)",
         subtitle = "Higher = more stable selection across bootstrap iterations",
         x = "Method", y = "Nogueira Stability Index", fill = "Dataset") +
    theme_bindlab() +
    theme(legend.position = "bottom")

  save_methods_fig(p_nog, "Fig04_nogueira_stability", width = 8, height = 5.5)

  return(invisible(ng_df))
}


# ==============================================================================
# PART 5: Real data application summary
# ==============================================================================

#' Fig.5: Cross-dataset pipeline behavior consistency
plot_real_data_summary <- function() {

  master_file <- file.path(BENCH_DIR, "benchmark_master.csv")
  if (!file.exists(master_file)) {
    methods_log("S07_FIG", "benchmark_master.csv not found")
    return(invisible(NULL))
  }

  df <- read.csv(master_file, stringsAsFactors = FALSE)
  real_df <- df[df$mode == "real", ]

  if (nrow(real_df) < 2) {
    methods_log("S07_FIG", "Need ≥2 real datasets for comparison")
    return(invisible(NULL))
  }

  # Panel A: Algorithm frequency distribution comparison
  freq_data <- real_df %>%
    select(run_id,
           LASSO = lasso_n_nonzero, RF = rf_n_nonzero, `SVM-RFE` = svm_n_nonzero) %>%
    pivot_longer(-run_id, names_to = "Algorithm", values_to = "n_nonzero")

  p_freq <- ggplot(freq_data, aes(x = Algorithm, y = n_nonzero, fill = run_id)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Algorithm Signal Concentration",
         subtitle = "Number of genes with nonzero selection frequency",
         x = NULL, y = "Genes with freq > 0", fill = "Dataset") +
    theme_bindlab() +
    theme(legend.position = "bottom")

  # Panel B: Funnel comparison
  funnel_data <- real_df %>%
    select(run_id, n_candidate_pool, n_ml_final, n_ppi_hubs, n_final) %>%
    pivot_longer(-run_id, names_to = "Step", values_to = "Count") %>%
    mutate(Step = factor(Step,
                         levels = c("n_candidate_pool", "n_ml_final", "n_ppi_hubs", "n_final"),
                         labels = c("Candidate pool", "ML markers", "PPI hubs", "Final union")))

  p_funnel <- ggplot(funnel_data, aes(x = Step, y = Count, fill = run_id)) +
    geom_bar(stat = "identity", position = position_dodge(0.8), width = 0.7) +
    geom_text(aes(label = Count),
              position = position_dodge(0.8), vjust = -0.3, size = 3) +
    scale_fill_brewer(palette = "Set2") +
    labs(title = "Pipeline Funnel Across Datasets",
         x = NULL, y = "Number of genes", fill = "Dataset") +
    theme_bindlab() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(angle = 20, hjust = 1))

  # Combine panels
  p_combined <- plot_grid(p_freq, p_funnel, nrow = 1, labels = c("A", "B"),
                           rel_widths = c(1, 1.2))

  save_methods_fig(p_combined, "Fig05_real_data_summary", width = 12, height = 5.5)

  return(invisible(real_df))
}


# ==============================================================================
# PART 6: Summary tables
# ==============================================================================

#' Table 1: Simulation benchmark summary
generate_table1 <- function() {

  master_file <- file.path(BENCH_DIR, "benchmark_master.csv")
  if (!file.exists(master_file)) return(NULL)

  df <- read.csv(master_file, stringsAsFactors = FALSE)
  sim_df <- df[df$mode == "simulation" & !is.na(df$precision), ]
  if (nrow(sim_df) == 0) return(NULL)

  parts <- strsplit(sim_df$run_id, "_")
  sim_df$snr     <- sapply(parts, `[`, 1)
  sim_df$sample  <- sapply(parts, `[`, 2)
  sim_df$density <- sapply(parts, `[`, 3)

  tab1 <- sim_df %>%
    group_by(snr, sample, density) %>%
    summarise(
      n_runs     = n(),
      precision  = sprintf("%.3f ± %.3f", mean(precision), sd(precision)),
      recall     = sprintf("%.3f ± %.3f", mean(recall), sd(recall)),
      F1         = sprintf("%.3f ± %.3f", mean(F1), sd(F1)),
      FDR        = sprintf("%.3f ± %.3f", mean(FDR), sd(FDR)),
      n_selected = sprintf("%.1f ± %.1f", mean(n_final), sd(n_final)),
      .groups    = "drop"
    ) %>%
    arrange(factor(snr, levels = c("low", "medium", "high")),
            factor(sample, levels = c("small", "medium", "large")))

  write.csv(tab1, file.path(TAB_DIR_METHODS, "Table01_benchmark_summary.csv"),
            row.names = FALSE)
  methods_log("S07_TABLE", sprintf("Table 1: %d rows", nrow(tab1)))
  return(tab1)
}


#' Table 2: Ablation study summary
generate_table2 <- function() {

  abl_file <- file.path(BENCH_DIR, "ablation_all_simulations.csv")
  if (!file.exists(abl_file)) return(NULL)

  df <- read.csv(abl_file, stringsAsFactors = FALSE)
  if (!"F1" %in% colnames(df)) return(NULL)

  df_delta <- df %>%
    group_by(run_id) %>%
    mutate(F1_base = F1[config == "A0_FULL"],
           delta_F1 = F1 - F1_base) %>%
    ungroup()

  tab2 <- df_delta %>%
    group_by(config, label) %>%
    summarise(
      n_final_mean  = sprintf("%.1f ± %.1f", mean(n_final), sd(n_final)),
      F1_mean       = sprintf("%.3f ± %.3f", mean(F1), sd(F1)),
      delta_F1      = sprintf("%+.3f ± %.3f", mean(delta_F1), sd(delta_F1)),
      precision     = sprintf("%.3f ± %.3f", mean(precision), sd(precision)),
      recall        = sprintf("%.3f ± %.3f", mean(recall), sd(recall)),
      .groups       = "drop"
    ) %>%
    arrange(match(config, names(ABLATION_CONFIGS)))

  write.csv(tab2, file.path(TAB_DIR_METHODS, "Table02_ablation_summary.csv"),
            row.names = FALSE)
  methods_log("S07_TABLE", sprintf("Table 2: %d rows", nrow(tab2)))
  return(tab2)
}


#' Table 3: Cross-dataset consistency
generate_table3 <- function(run_dirs = NULL) {

  if (is.null(run_dirs)) {
    run_dirs <- c(
      file.path(RUN_DIR, "GEO_GSE307424_Lung"),
      file.path(RUN_DIR, "GEO_GSE197067_Tcell")
    )
    run_dirs <- run_dirs[dir.exists(run_dirs)]
  }

  rows <- list()
  for (rd in run_dirs) {
    run_id <- basename(rd)
    data_dir <- file.path(rd, "data")

    row <- data.frame(dataset = run_id, stringsAsFactors = FALSE)

    # Candidate pool
    pool_f <- file.path(data_dir, "candidate_pool.rds")
    if (file.exists(pool_f)) {
      pool <- readRDS(pool_f)
      row$n_pool <- length(pool$candidate_pool)
    }

    # Stability scores
    stab_f <- file.path(data_dir, "ml_stability_selection.rds")
    if (file.exists(stab_f)) {
      stab <- readRDS(stab_f)
      row$lasso_nonzero <- sum(stab$lasso_freq > 0)
      row$rf_nonzero    <- sum(stab$rf_freq > 0)
      row$svm_nonzero   <- sum(stab$svm_freq > 0)
      row$rf_max_freq   <- round(max(stab$rf_freq), 3)
    }

    # ML final
    ml_f <- file.path(data_dir, "ml_gap_union.rds")
    if (file.exists(ml_f)) {
      ml <- readRDS(ml_f)
      row$n_ml <- length(ml$final_gene_ids)
    }

    # PPI
    ppi_f <- file.path(data_dir, "PPI_09F_hub_genes.csv")
    if (file.exists(ppi_f)) {
      ppi <- read.csv(ppi_f, stringsAsFactors = FALSE)
      row$n_ppi <- nrow(ppi)
    }

    # Final
    fin_f <- file.path(data_dir, "Final_candidate_genes.csv")
    if (file.exists(fin_f)) {
      fin <- read.csv(fin_f, stringsAsFactors = FALSE)
      row$n_final <- nrow(fin)
    }

    # Nogueira
    ng <- compute_nogueira_for_run(rd)
    if (!is.null(ng)) {
      for (i in seq_len(nrow(ng))) {
        col_name <- paste0("nogueira_", tolower(ng$method[i]))
        row[[col_name]] <- ng$nogueira[i]
      }
    }

    rows[[run_id]] <- row
  }

  tab3 <- do.call(rbind, rows)
  write.csv(tab3, file.path(TAB_DIR_METHODS, "Table03_cross_dataset_consistency.csv"),
            row.names = FALSE)
  methods_log("S07_TABLE", sprintf("Table 3: %d datasets", nrow(tab3)))
  return(tab3)
}


# ==============================================================================
# PART 7: Generate all figures at once
# ==============================================================================

generate_all <- function(run_dirs = NULL) {

  methods_log("S07_MAIN", "=== Generating all figures and tables ===")

  dir.create(FIG_DIR_METHODS, recursive = TRUE, showWarnings = FALSE)
  dir.create(TAB_DIR_METHODS, recursive = TRUE, showWarnings = FALSE)

  # Tables
  methods_log("S07_MAIN", "--- Tables ---")
  tab1 <- generate_table1()
  tab2 <- generate_table2()
  tab3 <- generate_table3(run_dirs)

  # Figures
  methods_log("S07_MAIN", "--- Figures ---")
  plot_simulation_benchmark()
  plot_ablation_results()
  plot_stability_comparison(run_dirs)
  plot_real_data_summary()

  methods_log("S07_MAIN", "=== All outputs generated ===")
  methods_log("S07_MAIN", sprintf("  Figures: %s", FIG_DIR_METHODS))
  methods_log("S07_MAIN", sprintf("  Tables:  %s", TAB_DIR_METHODS))
}


# ==============================================================================
# Direct execution
# ==============================================================================

if (sys.nframe() == 0) {
  cat("\n")
  cat("================================================================\n")
  cat("  S07_cross_dataset_summary.R — Usage\n")
  cat("================================================================\n")
  cat("\n")
  cat("  # Generate all figures:\n")
  cat("  generate_all()\n")
  cat("\n")
  cat("  # Calculate Nogueira stability individually:\n")
  cat("  ng <- compute_nogueira_for_run('path/to/run')\n")
  cat("\n")
  cat("  # Generate individual figure:\n")
  cat("  plot_simulation_benchmark()\n")
  cat("  plot_ablation_results()\n")
  cat("  plot_stability_comparison()\n")
  cat("  plot_real_data_summary()\n")
  cat("\n")
  cat("  # Generate individual table:\n")
  cat("  generate_table1()  # Simulation benchmark\n")
  cat("  generate_table2()  # Ablation study\n")
  cat("  generate_table3()  # Cross-dataset consistency\n")
  cat("================================================================\n")
}
