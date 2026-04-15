#!/usr/bin/env Rscript
# ==============================================================================
# 10_integration.R — Multi-evidence integration table + funnel
# Updated: 2026-03-23 (09D gap-union version)
#
# Changes from 2026-03-19 version:
#   - ML source: ml_gap_union.rds (09D) replaces ml_feature_selection.rds (09B)
#   - PPI is post-hoc annotation only (09E), NOT a parallel selection line
#   - No ML+PPI union — final candidates = ML markers from 09D
#   - Funnel reflects: txome -> maSigPro -> WGCNA -> effect size -> ML gap-union
# ==============================================================================

source(file.path(SCRIPT_DIR, "00_setup.R"))
source(file.path(SCRIPT_DIR, "theme_bindlab.R"))

log_step("10_INT", "Starting integration (2026-03-23: gap-union version)")

# ---------------------------------------------------------------------------
# 0. Upstream check
# ---------------------------------------------------------------------------

ml_gap_file <- file.path(DATA_DIR, "ml_gap_union.rds")
pool_file   <- file.path(DATA_DIR, "candidate_pool.rds")

# 09E PPI is optional (needs internet)
ppi_interaction_file <- file.path(DATA_DIR, "ML_09E_ppi_interactions.csv")
ppi_enrichment_file  <- file.path(DATA_DIR, "ML_09E_enrichment.csv")
ppi_neighbor_file    <- file.path(DATA_DIR, "ML_09E_ppi_neighbor_summary.csv")

check_upstream("10_INT",
  upstream_files = c(
    "09D_gap_union" = ml_gap_file,
    "09A_pool"      = pool_file
  ),
  output_files = c(
    FILES$hub_genes,
    file.path(DATA_DIR, "Final_candidate_genes.csv"),
    file.path(DATA_DIR, "screening_funnel_data.csv")
  )
)

# ---------------------------------------------------------------------------
# 1. Load all upstream results
# ---------------------------------------------------------------------------

all_results      <- readRDS(FILES$deg_results)
wgcna_results    <- readRDS(FILES$wgcna_results)
masigpro_results <- tryCatch(readRDS(FILES$masigpro_results), error = function(e) NULL)
gene_anno        <- readRDS(FILES$gene_annotation)
tpm_filtered     <- readRDS(FILES$tpm_filtered)
kinetic_df       <- tryCatch(readRDS(FILES$kinetic_results), error = function(e) NULL)

# 09D: ML gap-union results (new)
output_09d       <- readRDS(ml_gap_file)
ml_final_genes   <- output_09d$final_genes      # data.frame with symbol, ensembl_id, source_algos...
ml_gene_ids      <- output_09d$final_gene_ids    # character vector of ensembl_ids

log_step("10_INT", sprintf("09D ML markers: %d genes (%s)",
  length(ml_gene_ids), paste(ml_final_genes$symbol, collapse = ", ")))

# 09A: candidate pool
pool_result <- readRDS(pool_file)
candidate_pool                <- pool_result$candidate_pool
n_masigpro_before_interaction <- pool_result$n_masigpro_raw
n_masigpro_after_interaction  <- pool_result$n_masigpro_interaction
n_pool_before_effect          <- pool_result$n_pool_pre_effect
n_pool_after_effect           <- pool_result$n_pool_post_effect
effect_filter_applied         <- pool_result$effect_filter_applied
masigpro_ensembl              <- pool_result$masigpro_ensembl_filtered
wgcna_key_module_genes        <- pool_result$wgcna_key_genes

log_step("10_INT", sprintf("Candidate pool from 09A: %d genes", length(candidate_pool)))

# 09E PPI annotation (optional, post-hoc)
ppi_neighbors <- NULL
ppi_enrichment <- NULL
has_ppi <- FALSE
if (file.exists(ppi_interaction_file)) {
  has_ppi <- TRUE
  log_step("10_INT", "09E PPI annotation found (post-hoc)")
  if (file.exists(ppi_neighbor_file)) {
    ppi_neighbors <- read.csv(ppi_neighbor_file, stringsAsFactors = FALSE)
  }
  if (file.exists(ppi_enrichment_file)) {
    ppi_enrichment <- read.csv(ppi_enrichment_file, stringsAsFactors = FALSE)
  }
} else {
  log_step("10_INT", "09E PPI not found (optional, run 09E if needed)")
}

# LRT sig genes
lrt_sig    <- all_results$lrt_interaction_sig
id_col     <- if ("ensembl_id" %in% colnames(lrt_sig)) "ensembl_id" else "ensembl_gene_id"
lrt_sig_ids <- unique(lrt_sig[[id_col]])

# WGCNA hub (annotation only)
wgcna_hub_symbols <- character(0)
tryCatch({
  wgcna_hub_raw <- dplyr::bind_rows(
    wgcna_results$signed$hub_genes,
    wgcna_results$unsigned$hub_genes
  )
  if (nrow(wgcna_hub_raw) > 0 && "is_hub" %in% colnames(wgcna_hub_raw)) {
    wgcna_hub_symbols <- unique(wgcna_hub_raw$symbol[wgcna_hub_raw$is_hub == TRUE & !is.na(wgcna_hub_raw$symbol)])
  }
}, error = function(e) log_step("10_INT", sprintf("WGCNA hub: %s", e$message)))

# ---------------------------------------------------------------------------
# 2. Build integration table (all genes)
# ---------------------------------------------------------------------------

log_step("10_INT", "Building integrated evidence table...")

all_genes <- data.frame(
  ensembl_id = rownames(tpm_filtered),
  stringsAsFactors = FALSE
)

# Gene annotation
all_genes$symbol    <- gene_anno$hgnc_symbol[match(all_genes$ensembl_id, gene_anno$ensembl_gene_id)]
all_genes$entrez_id <- gene_anno$entrez_id[match(all_genes$ensembl_id, gene_anno$ensembl_gene_id)]

# LFC and padj per timepoint
for (tp in PARAMS$time_labels) {
  wald_res <- tryCatch(all_results$wald_by_time[[tp]], error = function(e) NULL)
  if (!is.null(wald_res)) {
    wald_id_col <- if ("ensembl_id" %in% colnames(wald_res)) "ensembl_id" else "ensembl_gene_id"
    idx <- match(all_genes$ensembl_id, wald_res[[wald_id_col]])
    all_genes[[paste0("LFC_", tp)]]  <- wald_res$log2FoldChange[idx]
    all_genes[[paste0("padj_", tp)]] <- wald_res$padj[idx]
  }
}

# Kinetic class (annotation only)
if (!is.null(kinetic_df) && "kinetic_class" %in% colnames(kinetic_df)) {
  kc_col <- if ("ensembl_id" %in% colnames(kinetic_df)) "ensembl_id" else "ensembl_gene_id"
  all_genes$kinetic_class <- kinetic_df$kinetic_class[match(all_genes$ensembl_id, kinetic_df[[kc_col]])]
} else {
  all_genes$kinetic_class <- NA
}

# Evidence flags
all_genes <- all_genes %>%
  dplyr::mutate(
    lrt_sig           = ensembl_id %in% lrt_sig_ids,
    in_maSigPro       = ensembl_id %in% masigpro_ensembl,
    in_WGCNA_key      = ensembl_id %in% wgcna_key_module_genes,
    in_candidate_pool = ensembl_id %in% candidate_pool,
    is_WGCNA_hub      = symbol %in% wgcna_hub_symbols,
    # ML: from 09D gap-union (final selection)
    in_ML             = ensembl_id %in% ml_gene_ids,
    # PPI: post-hoc annotation flag (NOT selection)
    has_PPI_annotation = has_ppi & in_ML  # only annotated for ML genes
  )

# Add 09D source algorithm info for ML genes
ml_source_map <- setNames(ml_final_genes$source_algos, ml_final_genes$ensembl_id)
ml_freq_map   <- setNames(ml_final_genes$stability_score, ml_final_genes$ensembl_id)
all_genes$ml_source_algos   <- ml_source_map[all_genes$ensembl_id]
all_genes$ml_stability_score <- ml_freq_map[all_genes$ensembl_id]

# Add PPI neighbor info if available
if (!is.null(ppi_neighbors) && nrow(ppi_neighbors) > 0) {
  shared_partners <- ppi_neighbors$partner[ppi_neighbors$n_connections >= 2]
  all_genes$ppi_shared_neighbor <- all_genes$symbol %in% shared_partners
} else {
  all_genes$ppi_shared_neighbor <- FALSE
}

# ---------------------------------------------------------------------------
# 2b. qPCR feasibility annotation (SEQC/MAQC-III standard)
#     Ct ≤ 35 = detectable (SEQC Nat Biotechnol 2014)
#     log2(TPM) correlates linearly with Ct (Everaert 2017 Sci Rep)
#     TPM > 32 (log2 > 5) = reliably quantifiable (Vasiliu PMC7026138)
# ---------------------------------------------------------------------------

log_step("10_INT", "Adding qPCR feasibility annotation...")

sample_info <- readRDS(FILES$sample_info)
idx_ind <- which(sample_info$Treatment == "Induced")
idx_ctl <- which(sample_info$Treatment == "Control")

all_genes$Induced_mean_TPM <- rowMeans(tpm_filtered[all_genes$ensembl_id, idx_ind, drop=FALSE])
all_genes$Control_mean_TPM <- rowMeans(tpm_filtered[all_genes$ensembl_id, idx_ctl, drop=FALSE])
all_genes$Induced_max_TPM  <- apply(tpm_filtered[all_genes$ensembl_id, idx_ind, drop=FALSE], 1, max)
all_genes$FC_Ind_vs_Ctrl   <- (all_genes$Induced_mean_TPM + 0.1) / (all_genes$Control_mean_TPM + 0.1)
all_genes$log2FC_TPM       <- log2(all_genes$FC_Ind_vs_Ctrl)

# qPCR feasibility rating based on Induced group expression
# qPCR feasibility rating based on the higher-expressed group
# For upregulated genes, Induced is high; for downregulated, Control is high
# Use max(Induced, Control) as the detectable expression level
all_genes$detection_TPM <- pmax(all_genes$Induced_mean_TPM, all_genes$Control_mean_TPM)
all_genes$qPCR_feasibility <- dplyr::case_when(
  all_genes$detection_TPM >= 32 ~ "easy (Ct~25-28)",
  all_genes$detection_TPM >= 5  ~ "moderate (Ct~28-33)",
  all_genes$detection_TPM >= 1  ~ "challenging (Ct~33-35)",
  TRUE                          ~ "below detection (Ct>35)"
)

log_step("10_INT", "Evidence counts:")
log_step("10_INT", sprintf("  Filtered transcriptome: %d", nrow(all_genes)))
log_step("10_INT", sprintf("  maSigPro (R2>=%.1f): %d", PARAMS$masigpro_rsq, n_masigpro_before_interaction))
log_step("10_INT", sprintf("  maSigPro (interaction P<0.05): %d", n_masigpro_after_interaction))
log_step("10_INT", sprintf("  WGCNA key modules: %d", length(wgcna_key_module_genes)))
log_step("10_INT", sprintf("  Candidate pool: %d", length(candidate_pool)))
log_step("10_INT", sprintf("  ML markers (09D gap-union): %d", sum(all_genes$in_ML)))

# ---------------------------------------------------------------------------
# 3. Save outputs
# ---------------------------------------------------------------------------

# Full table
write.csv(all_genes, FILES$hub_genes, row.names = FALSE)
log_step("10_INT", sprintf("Full table: %s (%d genes)", FILES$hub_genes, nrow(all_genes)))

# Final candidates = ML markers
candidates <- all_genes %>% dplyr::filter(in_ML)
write.csv(candidates, file.path(DATA_DIR, "Final_candidate_genes.csv"), row.names = FALSE)
log_step("10_INT", sprintf("Final candidates: %d", nrow(candidates)))

# qPCR feasibility summary for candidates
log_step("10_INT", "--- qPCR feasibility for ML candidates ---")
for (i in seq_len(nrow(candidates))) {
  r <- candidates[i, ]
  log_step("10_INT", sprintf(
    "  %-10s Induced=%.1f Control=%.1f FC=%.1f log2FC=%.2f  [%s]",
    r$symbol, r$Induced_mean_TPM, r$Control_mean_TPM,
    r$FC_Ind_vs_Ctrl, r$log2FC_TPM, r$qPCR_feasibility))
}
log_step("10_INT", sprintf("  qPCR easy (detection TPM>=32): %d/%d",
  sum(candidates$detection_TPM >= 32), nrow(candidates)))
log_step("10_INT", sprintf("  qPCR moderate (detection TPM 5-32): %d/%d",
  sum(candidates$detection_TPM >= 5 & candidates$detection_TPM < 32), nrow(candidates)))
log_step("10_INT", sprintf("  qPCR challenging (detection TPM <5): %d/%d",
  sum(candidates$detection_TPM < 5), nrow(candidates)))

# ML detail
write.csv(candidates, file.path(DATA_DIR, "ML_feature_consensus.csv"), row.names = FALSE)

# ---------------------------------------------------------------------------
# 4. Funnel data
# ---------------------------------------------------------------------------
# Pipeline: txome -> maSigPro(R2) -> maSigPro(interaction) -> ∩WGCNA
#           -> effect size -> candidate pool -> ML | PPI parallel -> union

# Load PPI hub count
ppi_hub_file <- file.path(DATA_DIR, "ppi_hub_selection.rds")
n_ppi_hubs <- 0
if (file.exists(ppi_hub_file)) {
  ppi_result <- readRDS(ppi_hub_file)
  n_ppi_hubs <- nrow(ppi_result$hub_genes)
}

n_ml <- sum(all_genes$in_ML)

# ML and PPI union (remove duplicates)
ml_gene_syms <- ml_final_genes$symbol
ppi_gene_syms <- if (n_ppi_hubs > 0) ppi_result$hub_genes$symbol else character(0)
final_union_syms <- unique(c(ml_gene_syms, ppi_gene_syms))
n_final_union <- length(final_union_syms)
n_overlap <- length(intersect(ml_gene_syms, ppi_gene_syms))

funnel_steps <- c("Filtered transcriptome",
                  sprintf("maSigPro (R\u00b2\u2265%.1f)", PARAMS$masigpro_rsq),
                  "maSigPro (interaction P<0.05)",
                  "maSigPro(filtered) \u2229 WGCNA")
funnel_ns <- c(nrow(all_genes),
               n_masigpro_before_interaction,
               n_masigpro_after_interaction,
               n_pool_before_effect)

if (effect_filter_applied) {
  funnel_steps <- c(funnel_steps, "Effect size (|LFC|>1 test)")
  funnel_ns    <- c(funnel_ns, n_pool_after_effect)
}

# Only add "Candidate pool" if it differs from effect size result
if (!effect_filter_applied || length(candidate_pool) != n_pool_after_effect) {
  funnel_steps <- c(funnel_steps, "Candidate pool")
  funnel_ns    <- c(funnel_ns, length(candidate_pool))
}

funnel_steps <- c(funnel_steps,
                  sprintf("ML markers (gap-union): %d", n_ml),
                  sprintf("PPI hubs (5-algo consensus): %d", n_ppi_hubs),
                  sprintf("Final (ML \u222a PPI, overlap=%d): %d", n_overlap, n_final_union))
funnel_ns    <- c(funnel_ns, n_ml, n_ppi_hubs, n_final_union)

funnel_data <- data.frame(step = funnel_steps, n = funnel_ns)
write.csv(funnel_data, file.path(DATA_DIR, "screening_funnel_data.csv"), row.names = FALSE)

log_step("10_INT", "Screening funnel:")
for (i in seq_len(nrow(funnel_data))) {
  log_step("10_INT", sprintf("  %s: n = %d", funnel_data$step[i], funnel_data$n[i]))
}

# ---------------------------------------------------------------------------
# 5. Funnel plot
# ---------------------------------------------------------------------------

tryCatch({
  library(ggplot2)
  
  funnel_data$step <- factor(funnel_data$step, levels = rev(funnel_data$step))
  
  step_colors <- sapply(as.character(funnel_data$step), function(s) {
    if (grepl("Filtered transcriptome", s))      "#2E8B57"
    else if (grepl("maSigPro.*R", s))            "#FF8C00"
    else if (grepl("interaction", s))            "#FFA500"
    else if (grepl("WGCNA", s))                  "#FF6347"
    else if (grepl("Effect size", s))            "#CD5C5C"
    else if (grepl("Candidate pool", s))         "#9B59B6"
    else if (grepl("ML markers", s))             "#27AE60"
    else if (grepl("PPI hubs", s))               "#2E86AB"
    else if (grepl("Final", s))                  "#E74C3C"
    else                                         "#999999"
  })
  names(step_colors) <- as.character(funnel_data$step)
  
  p_funnel <- ggplot(funnel_data, aes(x = n, y = step)) +
    geom_col(aes(fill = step), show.legend = FALSE) +
    geom_text(aes(label = sprintf("n = %s", format(n, big.mark = ","))),
              hjust = -0.1, size = 3.8, fontface = "bold") +
    scale_fill_manual(values = step_colors) +
    labs(title = sprintf("Screening Funnel \u2014 %s", PARAMS$diff_type),
         subtitle = "maSigPro \u2229 WGCNA \u2192 effect size \u2192 ML(bootstrap+gap) | PPI(5-algo consensus) \u2192 union",
         x = "Number of genes", y = NULL) +
    theme_bindlab_box(base_size = 11) +
    theme(axis.text.y = element_text(size = 11),
          panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3)) +
    xlim(0, max(funnel_data$n) * 1.15)
  
  save_pub_fig(p_funnel, "Fig01_screening_funnel", "10_Integration",
               width = 12, height = 7)
  log_step("10_INT", "Funnel plot saved")
}, error = function(e) log_step("10_INT", sprintf("Funnel plot error: %s", e$message)))

# ---------------------------------------------------------------------------
# 5b. qPCR feasibility plot for final candidates
# ---------------------------------------------------------------------------

tryCatch({
  library(ggplot2)
  
  if (nrow(candidates) >= 2) {
    # Prepare data: Induced vs Control TPM for each candidate
    cand_plot <- candidates[, c("symbol", "Induced_mean_TPM", "Control_mean_TPM",
                                 "log2FC_TPM", "qPCR_feasibility")]
    cand_plot <- cand_plot[order(-cand_plot$Induced_mean_TPM), ]
    cand_plot$symbol <- factor(cand_plot$symbol, levels = rev(cand_plot$symbol))
    
    # Long format for grouped bar
    cand_long <- tidyr::pivot_longer(
      cand_plot[, c("symbol", "Induced_mean_TPM", "Control_mean_TPM")],
      cols = c("Induced_mean_TPM", "Control_mean_TPM"),
      names_to = "group", values_to = "TPM"
    )
    cand_long$group <- gsub("_mean_TPM", "", cand_long$group)
    cand_long$group <- factor(cand_long$group, levels = c("Induced", "Control"))
    
    # Plot: grouped bar of Induced vs Control TPM
    p_tpm <- ggplot(cand_long, aes(x = symbol, y = TPM, fill = group)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.6) +
      geom_text(aes(label = sprintf("%.1f", TPM)),
                position = position_dodge(width = 0.7), hjust = -0.1, size = 2.8) +
      geom_hline(yintercept = 32, linetype = "dashed", color = "#E64B35", linewidth = 0.5) +
      geom_hline(yintercept = 5, linetype = "dotted", color = "#F39B7F", linewidth = 0.5) +
      annotate("text", x = 0.6, y = 34, label = "TPM=32 (reliable qPCR)",
               hjust = 0, color = "#E64B35", size = 2.5, fontface = "italic") +
      annotate("text", x = 0.6, y = 6.5, label = "TPM=5 (detection limit)",
               hjust = 0, color = "#F39B7F", size = 2.5, fontface = "italic") +
      coord_flip() +
      scale_fill_manual(values = c("Induced" = "#E64B35", "Control" = "#4DBBD5"),
                        name = "Group", labels = COLORS$treatment_labels) +
      scale_y_continuous(trans = "log1p",
                         breaks = c(0, 1, 5, 10, 32, 100, 500, 2000),
                         labels = c("0","1","5","10","32","100","500","2000")) +
      labs(title = sprintf("qPCR Feasibility: ML Candidate Biomarkers — %s", PARAMS$diff_type),
           subtitle = "Induced vs Control mean TPM (log scale). Dashed line = reliable detection threshold.",
           x = NULL, y = "Mean TPM (log scale)") +
      theme_bindlab(base_size = 10) +
      theme(legend.position = "bottom",
            axis.text.y = element_text(face = "italic"))
    
    save_pub_fig(p_tpm, "Fig_qPCR_feasibility_candidates", "10_Integration",
                 width = 10, height = max(5, nrow(candidates) * 0.5 + 2))
    log_step("10_INT", "qPCR feasibility plot saved")
    
    # Plot 2: log2FC with qPCR feasibility color-coded
    cand_fc <- candidates[, c("symbol", "log2FC_TPM", "Induced_mean_TPM",
                               "Control_mean_TPM", "qPCR_feasibility")]
    cand_fc <- cand_fc[order(-abs(cand_fc$log2FC_TPM)), ]
    cand_fc$symbol <- factor(cand_fc$symbol, levels = rev(cand_fc$symbol))
    
    # For qPCR feasibility: upregulated genes check Induced TPM,
    # downregulated genes check Control TPM (that's where they're expressed)
    cand_fc$detection_TPM <- ifelse(cand_fc$log2FC_TPM >= 0,
      cand_fc$Induced_mean_TPM, cand_fc$Control_mean_TPM)
    cand_fc$feasibility <- dplyr::case_when(
      cand_fc$detection_TPM >= 32 ~ "Easy (TPM\u226532)",
      cand_fc$detection_TPM >= 5  ~ "Moderate (TPM 5-32)",
      TRUE                        ~ "Challenging (TPM<5)"
    )
    cand_fc$feasibility <- factor(cand_fc$feasibility,
      levels = c("Easy (TPM\u226532)", "Moderate (TPM 5-32)", "Challenging (TPM<5)"))
    
    # Build label: include both Induced and Control TPM
    cand_fc$label_text <- sprintf("FC=%.1f  Ind=%.0f  Ctrl=%.0f",
      2^abs(cand_fc$log2FC_TPM), cand_fc$Induced_mean_TPM, cand_fc$Control_mean_TPM)
    
    # X-axis range: find symmetric limits for labels
    fc_max <- max(abs(cand_fc$log2FC_TPM)) * 1.6
    
    # Feasibility color mapping — always show all 3 levels via drop=FALSE
    feas_all_colors <- c("Easy (TPM\u226532)" = "#00A087",
                         "Moderate (TPM 5-32)" = "#F39B7F",
                         "Challenging (TPM<5)" = "#E64B35")
    
    p_fc <- ggplot(cand_fc, aes(x = symbol, y = log2FC_TPM, fill = feasibility)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = label_text),
                hjust = ifelse(cand_fc$log2FC_TPM >= 0, -0.03, 1.03),
                size = 2.6, color = "grey30") +
      geom_hline(yintercept = 0, linewidth = 0.3, color = "grey50") +
      # Direction annotation
      annotate("text", x = nrow(cand_fc) + 0.5, y = fc_max * 0.5,
               label = "\u2191 Upregulated in Induced",
               size = 3, color = "#B22222", fontface = "bold") +
      annotate("text", x = nrow(cand_fc) + 0.5, y = -fc_max * 0.5,
               label = "\u2193 Downregulated in Induced",
               size = 3, color = "#2E5090", fontface = "bold") +
      coord_flip(clip = "off") +
      scale_fill_manual(values = feas_all_colors, drop = FALSE,
                        name = "qPCR feasibility") +
      scale_y_continuous(limits = c(-fc_max, fc_max)) +
      labs(title = sprintf("Marker Discriminative Power & qPCR Feasibility — %s", PARAMS$diff_type),
           subtitle = "Feasibility: \u2191 upregulated \u2192 check Induced TPM; \u2193 downregulated \u2192 check Control TPM",
           x = NULL, y = "log2(Induced/Control) TPM fold change") +
      theme_bindlab(base_size = 10) +
      theme(axis.text.y = element_text(face = "italic"),
            legend.position = "bottom",
            legend.box = "horizontal",
            legend.text = element_text(size = 8.5),
            legend.title = element_text(size = 9, face = "bold"),
            plot.subtitle = element_text(size = 9, color = "grey30"),
            plot.margin = margin(t = 10, r = 10, b = 10, l = 10)) +
      guides(fill = guide_legend(nrow = 1))
    
    save_pub_fig(p_fc, sprintf("%s_Fig_marker_FC_feasibility", PARAMS$diff_type),
                 "10_Integration",
                 width = 11, height = max(5.5, nrow(candidates) * 0.6 + 3))
    log_step("10_INT", "Marker FC + feasibility plot saved")
  }
}, error = function(e) log_step("10_INT", sprintf("qPCR feasibility plot error: %s", e$message)))

# ---------------------------------------------------------------------------
# 6. Console summary
# ---------------------------------------------------------------------------

cat("\n================================================================\n")
cat(sprintf("  %s Analysis \u2014 INTEGRATION COMPLETE\n", PARAMS$diff_type))
cat("================================================================\n")
cat(sprintf("  Filtered transcriptome: %d\n", nrow(all_genes)))
cat(sprintf("  maSigPro (R\u00b2\u2265%.1f): %d\n", PARAMS$masigpro_rsq, n_masigpro_before_interaction))
cat(sprintf("  maSigPro (interaction P<0.05): %d\n", n_masigpro_after_interaction))
cat(sprintf("  maSigPro(filtered) \u2229 WGCNA: %d\n", n_pool_before_effect))
if (effect_filter_applied) {
  cat(sprintf("  Effect size (|LFC|>1 test): %d\n", n_pool_after_effect))
}
cat(sprintf("  Candidate pool: %d\n", length(candidate_pool)))
cat(sprintf("  ML markers (09D gap-union): %d\n", sum(all_genes$in_ML)))
cat(sprintf("  PPI annotation: %s\n", ifelse(has_ppi, "available (post-hoc)", "not run")))
cat("================================================================\n")

# ---------------------------------------------------------------------------
# 7. Auto-call 10C
# ---------------------------------------------------------------------------

tryCatch({
  log_step("10_INT", "Running classic marker tracking (10C)...")
  source(file.path(SCRIPT_DIR, "10C_classic_gene_tracker.R"))
}, error = function(e) {
  log_step("10_INT", sprintf("10C skipped: %s", e$message))
})

log_step("10_INT", "Step 10 COMPLETE (gap-union version)")
