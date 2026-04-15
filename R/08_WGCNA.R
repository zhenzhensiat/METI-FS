#!/usr/bin/env Rscript
# ==============================================================================
# 08_WGCNA.R — ExpressionNetwork analysis ( 2026-03-19)
#
# within ：
# 1. 9 group unsigned → 2 times (unsigned+signed)
# 2. PowerSelection is WGCNA FAQDecision tree（auto + ConnectionCheck + FAQ）
# 3. Adipo_time → diff_time using 
# 4. GeneMAD top50%（WGCNAStandard）
#
# Literature basis:
#   - Langfelder & Horvath (2008) BMC Bioinformatics 9:559
# - WGCNA FAQ #6: Sample sizevspower
# - FAQ：unsigned 20-30Sample → power=8; signed → power=16
# ==============================================================================

source(file.path(SCRIPT_DIR, "00_setup.R"))
source(file.path(SCRIPT_DIR, "theme_bindlab.R"))

log_step("08_WGCNA", "Starting WGCNA (corrected: FAQ power decision tree)...")

# Load data
vst_mat     <- readRDS(FILES$vst_matrix)
sample_info <- readRDS(FILES$sample_info)
gene_anno   <- readRDS(FILES$gene_annotation)

n_samples <- ncol(vst_mat)
log_step("08_WGCNA", sprintf("VST matrix: %d genes × %d samples", nrow(vst_mat), n_samples))

# ===========================================================================
# 1. MADFilter：top50%Gene（WGCNAStandard）
# ===========================================================================

gene_mad <- apply(vst_mat, 1, mad)
gene_mad_sorted <- sort(gene_mad, decreasing = TRUE)
n_total_genes <- length(gene_mad_sorted)
n_keep <- round(n_total_genes * 0.50)
genes_use <- names(gene_mad_sorted)[1:n_keep]
mad_cutoff <- gene_mad_sorted[n_keep]

log_step("08_WGCNA", sprintf("MAD filtering: %d total → top50%% = %d genes (MAD >= %.4f)",
                              n_total_genes, n_keep, mad_cutoff))

# MADFigure
p_mad <- ggplot(data.frame(MAD = gene_mad), aes(x = MAD)) +
  geom_histogram(bins = 80, fill = "grey70", color = "white") +
  geom_vline(xintercept = mad_cutoff, linetype = 2, color = "#E64B35", linewidth = 1) +
  annotate("text", x = mad_cutoff, y = Inf, vjust = 2, hjust = -0.1,
           label = sprintf("top50%%\nn=%d", n_keep), color = "#E64B35", size = 3.5) +
  labs(title = "Gene MAD Distribution",
       subtitle = sprintf("WGCNA input: top 50%% by MAD (%d genes)", n_keep),
       x = "Median Absolute Deviation", y = "Count") +
  theme_bindlab()
save_pub_fig(p_mad, "WGCNA_MAD_filtering", "08_WGCNA", width = 8, height = 5)

# ===========================================================================
# 2. TraitMatrix（ using ， not ）
# ===========================================================================

group_levels <- c("Control_4d", "Control_7d", "Control_14d", "Control_21d",
                  "Induced_4d", "Induced_7d", "Induced_14d", "Induced_21d")

traits_group <- as.data.frame(sapply(group_levels, function(g) {
  as.numeric(sample_info$Group == g)
}))
rownames(traits_group) <- rownames(sample_info)
colnames(traits_group) <- group_levels

traits_continuous <- data.frame(
  diff_indicator = as.numeric(sample_info$Treatment == "Induced"),
  time_numeric    = sample_info$time_num,
  diff_time      = ifelse(sample_info$Treatment == "Induced", sample_info$time_num, 0),
  row.names = rownames(sample_info)
)

traits_mat_full <- cbind(traits_group, traits_continuous)

# ===========================================================================
# 3. WGCNA FAQ PowerDecision tree
# : WGCNA FAQ #6 + Langfelder Biostars 
#    Logic：
#      Step1: pickSoftThreshold auto (R²≥0.85)
# Step2: auto power > max_reasonable， or connectivity < min_conn → FAQ
# Step3: FAQ (Sample size→power)
# ===========================================================================

faq_power_decision <- function(sft, network_type, n_samples, run_label) {
  
  sft_df <- sft$fitIndices
  signed_r2 <- -sign(sft_df[, 3]) * sft_df[, 2]
  mean_conn <- sft_df[, 5]
  powers_vec <- sft_df[, 1]
  
 # FAQ power
  max_reasonable <- if (network_type == "unsigned" || network_type == "signed hybrid") 15 else 30
  
 # FAQ（Sample size → Defaultpower）
  faq_default <- if (network_type == "unsigned" || network_type == "signed hybrid") {
    if (n_samples < 20) 9
    else if (n_samples <= 30) 8
    else if (n_samples <= 40) 7
    else 6
  } else {  # signed
    if (n_samples < 20) 18
    else if (n_samples <= 30) 16
    else if (n_samples <= 40) 14
    else 12
  }
  
  min_conn_threshold <- 10  # Minimum acceptable connectivity
  
 # Step 1: auto
  auto_power <- sft$powerEstimate
  
  if (!is.na(auto_power) && auto_power <= max_reasonable) {
 # CheckConnection
    idx <- which(powers_vec == auto_power)
    conn_at_power <- if (length(idx) == 1) mean_conn[idx] else NA
    r2_at_power <- if (length(idx) == 1) signed_r2[idx] else NA
    
    if (!is.na(conn_at_power) && conn_at_power >= min_conn_threshold) {
      log_step("08_WGCNA", sprintf("[%s] Auto pick: power=%d, R²=%.3f, MeanConn=%.1f ✓",
                                    run_label, auto_power, r2_at_power, conn_at_power))
      return(list(power = auto_power, source = "auto",
                  r2 = r2_at_power, conn = conn_at_power))
    } else {
      log_step("08_WGCNA", sprintf("[%s] Auto power=%d has low connectivity (%.1f < %d), falling back to FAQ",
                                    run_label, auto_power, conn_at_power, min_conn_threshold))
    }
  } else {
    if (is.na(auto_power)) {
      log_step("08_WGCNA", sprintf("[%s] Auto pick returned NA (R² never reached 0.85)", run_label))
    } else {
      log_step("08_WGCNA", sprintf("[%s] Auto power=%d exceeds max reasonable (%d), falling back to FAQ",
                                    run_label, auto_power, max_reasonable))
    }
  }
  
 # Step 2: to FAQ
  idx_faq <- which(powers_vec == faq_default)
  r2_faq <- if (length(idx_faq) == 1) signed_r2[idx_faq] else NA
  conn_faq <- if (length(idx_faq) == 1) mean_conn[idx_faq] else NA
  
  log_step("08_WGCNA", sprintf("[%s] FAQ fallback: power=%d (n_samples=%d, %s), R²=%.3f, MeanConn=%.1f",
                                run_label, faq_default, n_samples, network_type,
                                ifelse(is.na(r2_faq), 0, r2_faq),
                                ifelse(is.na(conn_faq), 0, conn_faq)))
  
  if (!is.na(conn_faq) && conn_faq < min_conn_threshold) {
    log_step("08_WGCNA", sprintf("[%s] WARNING: FAQ power=%d still has low connectivity (%.1f). Consider adjusting gene filtering.",
                                  run_label, faq_default, conn_faq))
  }
  
  return(list(power = faq_default, source = "FAQ_table",
              r2 = ifelse(is.na(r2_faq), 0, r2_faq),
              conn = ifelse(is.na(conn_faq), 0, conn_faq)))
}

# ===========================================================================
# 4. Network construction（ using FAQDecision tree）
# ===========================================================================

build_wgcna <- function(datExpr, network_type, tom_type,
                         traits_mat, gene_anno_df, run_label,
                         min_mod_size = 25, deep_split = 2, merge_cut = 0.25) {
  
  log_step("08_WGCNA", sprintf(">>> [%s] Building %s network, %d genes <<<",
                                run_label, network_type, ncol(datExpr)))
  
  # --- PowerSelection ---
  powers_vec <- c(1:20, seq(22, 30, by = 2))
  sft <- pickSoftThreshold(datExpr,
                            powerVector = powers_vec,
                            networkType = network_type,
                            corFnc = "bicor",
                            verbose = 3)
  
  sft_df <- sft$fitIndices
  signed_r2 <- -sign(sft_df[, 3]) * sft_df[, 2]
  mean_conn <- sft_df[, 5]
  
  # FAQDecision tree
  power_decision <- faq_power_decision(sft, network_type, nrow(datExpr), run_label)
  power_sel <- power_decision$power
  
 # power 
  power_idx <- which(powers_vec == power_sel)
  sel_r2 <- if (length(power_idx) == 1) signed_r2[power_idx] else power_decision$r2
  sel_conn <- if (length(power_idx) == 1) mean_conn[power_idx] else power_decision$conn
  
 # --- ThresholdFigure ---
  tryCatch({
    p_sft <- ggplot(data.frame(Power = sft_df[, 1], R2 = signed_r2, Conn = mean_conn)) +
      geom_point(aes(Power, R2), color = "#E64B35", size = 2) +
      geom_line(aes(Power, R2), color = "#E64B35", alpha = 0.5) +
      geom_hline(yintercept = 0.85, linetype = 2, color = "grey50") +
      geom_vline(xintercept = power_sel, linetype = 2, color = "#3C5488", linewidth = 1) +
      annotate("text", x = power_sel, y = 0.3,
               label = sprintf("power=%d\n(%s)\nR²=%.3f\nConn=%.1f",
                                power_sel, power_decision$source, sel_r2, sel_conn),
               hjust = -0.1, size = 3, color = "#3C5488") +
      labs(title = sprintf("Scale-free Fit [%s]", run_label),
           x = "Soft Threshold (power)", y = "Scale Free Topology R²") +
      theme_bindlab()
    save_pub_fig(p_sft, sprintf("WGCNA_soft_threshold_%s", run_label),
                 "08_WGCNA", width = 7, height = 5)
  }, error = function(e) log_step("08_WGCNA", sprintf("SFT plot error: %s", e$message)))
  
  # --- Network construction ---
  net <- blockwiseModules(
    datExpr,
    power            = power_sel,
    networkType      = network_type,
    TOMType          = tom_type,
    corType          = "bicor",
    maxPOutliers     = 0.1,
    minModuleSize    = min_mod_size,
    deepSplit        = deep_split,
    mergeCutHeight   = merge_cut,
    numericLabels    = TRUE,
    pamRespectsDendro = FALSE,
    saveTOMs         = FALSE,
    verbose          = 3
  )
  
  modLabels <- net$colors
  modColors <- labels2colors(modLabels)
  names(modColors) <- colnames(datExpr)
  MEs <- net$MEs
  
  n_mods <- length(unique(modColors)) - ifelse("grey" %in% modColors, 1, 0)
  n_grey <- sum(modColors == "grey")
  pct_grey <- 100 * n_grey / ncol(datExpr)
  
  log_step("08_WGCNA", sprintf("[%s] Modules: %d (excl. grey), grey=%d (%.1f%%)",
                                run_label, n_mods, n_grey, pct_grey))
  
 # --- Figure ---
  tryCatch({
    save_heatmap_fig(function() {
      plotDendroAndColors(net$dendrograms[[1]], modColors[net$blockGenes[[1]]],
                          "Module colors", dendroLabels = FALSE,
                          main = sprintf("Gene Dendrogram [%s]", run_label))
    },
    filename = sprintf("WGCNA_dendrogram_%s", run_label),
    subdir = "08_WGCNA", width = 12, height = 6)
  }, error = function(e) log_step("08_WGCNA", sprintf("Dendro plot error: %s", e$message)))
  
 # --- Module-Trait ---
  moduleTraitCor <- WGCNA::cor(MEs, traits_mat[rownames(datExpr), ], use = "p", method = "pearson")
  moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  tryCatch({
    textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                        signif(moduleTraitPvalue, 1), ")", sep = "")
    dim(textMatrix) <- dim(moduleTraitCor)
    save_heatmap_fig(function() {
      labeledHeatmap(Matrix = moduleTraitCor,
                     xLabels = colnames(moduleTraitCor),
                     yLabels = rownames(moduleTraitCor),
                     ySymbols = rownames(moduleTraitCor),
                     colorLabels = FALSE,
                     colors = blueWhiteRed(50),
                     textMatrix = textMatrix, setStdMargins = FALSE,
                     cex.text = 0.5, zlim = c(-1, 1),
                     main = sprintf("Module-Trait Relationships [%s]", run_label))
    },
    filename = sprintf("WGCNA_module_trait_heatmap_%s", run_label),
    subdir = "08_WGCNA", width = 12, height = max(6, n_mods * 0.35))
  }, error = function(e) log_step("08_WGCNA", sprintf("Trait heatmap error: %s", e$message)))
  
 # --- R²（， not is ） ---
 # : WGCNA FAQ #6, PMC12457846
 # etc. level after ， not Threshold
  if (sel_r2 >= 0.80) {
    log_step("08_WGCNA", sprintf("[%s] R² DIAGNOSTIC: %.3f — GOOD (≥0.80)", run_label, sel_r2))
  } else if (sel_r2 >= 0.65) {
    log_step("08_WGCNA", sprintf("[%s] R² DIAGNOSTIC: %.3f — MARGINAL (0.65-0.80). Check dendrogram for distinct branches.", run_label, sel_r2))
  } else {
    log_step("08_WGCNA", sprintf("[%s] R² DIAGNOSTIC: %.3f — LOW (<0.65). Strong biological driver likely (FAQ #6). Check dendrogram and sample clustering tree carefully.", run_label, sel_r2))
    log_step("08_WGCNA", sprintf("[%s] Per Langfelder: low R² does not invalidate modules if branches are distinct. Power selected via FAQ table.", run_label))
  }

  # --- Significant modulesFilter ---
 # Key modulesDefine： and diff_time/diff_indicator |cor|>PARAMS$wgcna_sig_cor P<PARAMS$wgcna_sig_p（grey）
  key_traits <- intersect(c("diff_time", "diff_indicator"), colnames(traits_mat))
  sig_modules <- character(0)
  n_sig <- 0
  
  if (length(key_traits) > 0) {
    min_pval_per_ME <- apply(moduleTraitPvalue[, key_traits, drop = FALSE], 1, min)
    max_cor_per_ME  <- apply(abs(moduleTraitCor[, key_traits, drop = FALSE]), 1, max)
 # ANDLogic：pSignificant Correlation（00_setup.R）
 # greyModule：blockwiseModules using (ME0)，moduleEigengenes using (MEgrey)
    sig_mask <- (min_pval_per_ME < PARAMS$wgcna_sig_p & max_cor_per_ME > PARAMS$wgcna_sig_cor) &
                !grepl("^ME0$|^MEgrey$", rownames(moduleTraitCor))
    sig_modules <- rownames(moduleTraitCor)[sig_mask]
    n_sig <- length(sig_modules)
  }
  
  log_step("08_WGCNA", sprintf("[%s] Significant modules (p<%.2f AND |cor|>%.2f): %d — %s",
                                run_label, PARAMS$wgcna_sig_p, PARAMS$wgcna_sig_cor, n_sig,
                                if (n_sig > 0) paste(sig_modules, collapse = ", ") else "none"))
  
 # R²sig_modules， outside Info
  if (sel_r2 < 0.80 && n_sig <= 3) {
    log_step("08_WGCNA", sprintf(
      "[%s] NOTE: R²=%.3f with only %d sig modules. After pipeline completes, check:",
      run_label, sel_r2, n_sig))
    log_step("08_WGCNA", sprintf(
      "[%s]   1. Dendrogram: do modules show distinct branches?", run_label))
    log_step("08_WGCNA", sprintf(
      "[%s]   2. Module-trait heatmap: are there near-sig modules (|cor|=0.4-0.5)?", run_label))
    log_step("08_WGCNA", sprintf(
      "[%s]   3. Candidate pool size in 09A: if too small (<1000), consider relaxing sig_cor", run_label))
    
 # columns"" Module，
    if (length(key_traits) > 0) {
      near_sig_mask <- (min_pval_per_ME < 0.10 & max_cor_per_ME > 0.35) &
                       !sig_mask &
                       !grepl("^ME0$|^MEgrey$", rownames(moduleTraitCor))
      near_sig_mods <- rownames(moduleTraitCor)[near_sig_mask]
      if (length(near_sig_mods) > 0) {
        for (ns_mod in near_sig_mods) {
          ns_cor <- max(abs(moduleTraitCor[ns_mod, key_traits]))
          ns_p   <- min(moduleTraitPvalue[ns_mod, key_traits])
          log_step("08_WGCNA", sprintf(
            "[%s]   Near-sig module: %s (|cor|=%.3f, p=%.4f) — did not pass current threshold",
            run_label, ns_mod, ns_cor, ns_p))
        }
      }
    }
  }
  
 # --- HubGene（GS + MM） ---
  hub_genes_df <- data.frame()
  n_hub_total <- 0
  
  if (n_sig > 0 && "diff_time" %in% colnames(traits_mat)) {
    gs_trait <- traits_mat[rownames(datExpr), "diff_time", drop = FALSE]
    GS <- as.data.frame(bicor(datExpr, gs_trait, use = "p"))
    colnames(GS) <- "GS_diff_time"
    GS$absGS <- abs(GS$GS_diff_time)
    
    for (mod_me in sig_modules) {
      mod_label <- as.integer(gsub("^ME", "", mod_me))
      mod_color <- labels2colors(mod_label)
      mod_genes <- colnames(datExpr)[modColors == mod_color]
      if (length(mod_genes) < 5) next
      
      kME <- bicor(datExpr[, mod_genes], MEs[, mod_me, drop = FALSE], use = "p")
      
      mod_df <- data.frame(
        ensembl_id = mod_genes,
        module = mod_color,
        kME = abs(as.numeric(kME)),
        GS = GS$absGS[match(mod_genes, rownames(GS))],
        stringsAsFactors = FALSE
      )
      mod_df$symbol <- gene_anno_df$hgnc_symbol[match(mod_df$ensembl_id, gene_anno_df$ensembl_gene_id)]
      
 # HubDefine：|kME| > top10%Threshold |GS| > median
      kme_top10_cutoff <- quantile(mod_df$kME, 0.90, na.rm = TRUE)
      gs_median <- median(mod_df$GS, na.rm = TRUE)
      mod_df$is_hub <- mod_df$kME > kme_top10_cutoff & mod_df$GS > gs_median
      
      hub_genes_df <- dplyr::bind_rows(hub_genes_df, mod_df)
      n_hub_total <- n_hub_total + sum(mod_df$is_hub)
      
 # MM vs GSFigure
      tryCatch({
        p_mmgs <- ggplot(mod_df, aes(kME, GS)) +
          geom_point(alpha = 0.4, size = 1.5, color = mod_color) +
          geom_point(data = mod_df[mod_df$is_hub, ], color = "red", size = 2.5) +
          labs(title = sprintf("[%s] Module %s: MM vs GS", run_label, mod_color),
               x = "|Module Membership (kME)|",
               y = "|Gene Significance (diff_time)|") +
          theme_bindlab()
        save_pub_fig(p_mmgs, sprintf("WGCNA_MMGS_%s_%s", run_label, mod_color),
                     "08_WGCNA", width = 6, height = 5)
      }, error = function(e) NULL)
    }
  }
  
 # --- EigengeneFigure ---
  for (mod_me in sig_modules) {
    tryCatch({
      me_df <- data.frame(
        ME = MEs[, mod_me],
        Treatment = sample_info$Treatment,
        Time = sample_info$Time,
        Group = sample_info$Group,
        row.names = rownames(MEs)
      )
      p_me <- ggplot(me_df, aes(Time, ME, color = Treatment, group = Treatment)) +
        stat_summary(fun = mean, geom = "line", linewidth = 1) +
        stat_summary(fun = mean, geom = "point", size = 3) +
        stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
        scale_color_manual(values = COLORS$treatment) +
        labs(title = sprintf("[%s] Eigengene: %s", run_label, mod_me),
             x = "Time", y = "Module Eigengene") +
        theme_bindlab()
      save_pub_fig(p_me, sprintf("WGCNA_eigengene_%s_%s", run_label, sub("ME", "", mod_me)),
                   "08_WGCNA", width = 6, height = 4)
    }, error = function(e) NULL)
  }
  
 # --- diff_time ---
  diff_max_cor <- 0
  diff_best_pval <- 1
  if ("diff_time" %in% colnames(moduleTraitCor)) {
    diff_cors <- abs(moduleTraitCor[, "diff_time"])
    diff_cors <- diff_cors[names(diff_cors) != "MEgrey"]
    if (length(diff_cors) > 0) {
      diff_max_cor <- max(diff_cors)
      diff_best_pval <- min(moduleTraitPvalue[names(diff_cors), "diff_time"])
    }
  }
  
  summary_metrics <- data.frame(
    run_label       = run_label,
    network_type    = network_type,
    n_genes         = ncol(datExpr),
    power           = power_sel,
    power_source    = power_decision$source,
    R2_at_power     = sel_r2,
    mean_conn       = sel_conn,
    n_modules       = n_mods,
    pct_grey        = pct_grey,
    n_sig_modules   = n_sig,
    n_hub_genes     = n_hub_total,
    diff_max_cor   = diff_max_cor,
    diff_best_pval = diff_best_pval,
    stringsAsFactors = FALSE
  )
  
  list(
    net              = net,
    moduleLabels     = modLabels,
    moduleColors     = modColors,
    MEs              = MEs,
    power_selected   = power_sel,
    power_source     = power_decision$source,
    sft              = sft,
    moduleTraitCor   = moduleTraitCor,
    moduleTraitPvalue = moduleTraitPvalue,
    hub_genes        = hub_genes_df,
    sig_modules      = sig_modules,
    traits           = traits_mat,
    summary          = summary_metrics
  )
}

# ===========================================================================
# 5. Run：1 times unsigned（）+ 1 times signed（Reference）
# ===========================================================================

datExpr_main <- t(vst_mat[genes_use, ])
gsg <- goodSamplesGenes(datExpr_main, verbose = 0)
if (!gsg$allOK) {
  if (sum(!gsg$goodGenes) > 0) datExpr_main <- datExpr_main[, gsg$goodGenes]
  if (sum(!gsg$goodSamples) > 0) datExpr_main <- datExpr_main[gsg$goodSamples, ]
}
log_step("08_WGCNA", sprintf("datExpr: %d samples × %d genes", nrow(datExpr_main), ncol(datExpr_main)))

# --- 5a. Unsigned（， using below ML/PPI） ---
log_step("08_WGCNA", "=== Run 1/2: Unsigned network (primary) ===")
res_unsigned <- tryCatch(
  build_wgcna(
    datExpr      = datExpr_main,
    network_type = "unsigned",
    tom_type     = "unsigned",
    traits_mat   = traits_mat_full,
    gene_anno_df = gene_anno,
    run_label    = "unsigned",
    min_mod_size = 25,
    deep_split   = 2,
    merge_cut    = 0.25
  ),
  error = function(e) {
    log_step("08_WGCNA", sprintf("[unsigned] FAILED: %s", e$message))
    NULL
  }
)

# --- 5b. Signed hybrid（Reference， using ） ---
log_step("08_WGCNA", "=== Run 2/2: Signed hybrid network (reference) ===")
res_signed <- tryCatch(
  build_wgcna(
    datExpr      = datExpr_main,
    network_type = "signed hybrid",
    tom_type     = "signed",
    traits_mat   = traits_mat_full,
    gene_anno_df = gene_anno,
    run_label    = "signed",
    min_mod_size = 25,
    deep_split   = 2,
    merge_cut    = 0.25
  ),
  error = function(e) {
    log_step("08_WGCNA", sprintf("[signed] FAILED: %s", e$message))
    NULL
  }
)

# ===========================================================================
# 6. Save results（ below 09/09B/10 not ）
# ===========================================================================

wgcna_results <- list(
 # （ below using ）
  unsigned = res_unsigned,
  signed   = res_signed,
  datExpr  = datExpr_main
)

save_data(wgcna_results, FILES$wgcna_results)

# ===========================================================================
# 7. ExportHubGene
# ===========================================================================

all_hubs <- dplyr::bind_rows(
  if (!is.null(res_signed) && !is.null(res_signed$hub_genes) && nrow(res_signed$hub_genes) > 0)
    res_signed$hub_genes %>% dplyr::mutate(network = "signed") else NULL,
  if (!is.null(res_unsigned) && !is.null(res_unsigned$hub_genes) && nrow(res_unsigned$hub_genes) > 0)
    res_unsigned$hub_genes %>% dplyr::mutate(network = "unsigned") else NULL
)

if (nrow(all_hubs) > 0 && "is_hub" %in% colnames(all_hubs)) {
  hub_export <- all_hubs %>%
    dplyr::filter(is_hub) %>%
    dplyr::arrange(network, module, dplyr::desc(kME))
  
  if (nrow(hub_export) > 0) {
    write.csv(hub_export, file.path(DATA_DIR, "WGCNA_hub_genes.csv"), row.names = FALSE)
    log_step("08_WGCNA", sprintf("Hub genes exported: signed=%d, unsigned=%d",
                                  sum(hub_export$network == "signed"),
                                  sum(hub_export$network == "unsigned")))
  } else {
    # Fallback: top5% kME
    log_step("08_WGCNA", "WARNING: No strict hubs, exporting top5% kME as emergency hubs.")
    hub_export <- all_hubs %>%
      dplyr::group_by(network, module) %>%
      dplyr::mutate(kME_rank = dplyr::percent_rank(kME)) %>%
      dplyr::filter(kME_rank >= 0.95) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(is_hub = TRUE)
    write.csv(hub_export, file.path(DATA_DIR, "WGCNA_hub_genes.csv"), row.names = FALSE)
  }
} else {
  log_step("08_WGCNA", "WARNING: No hub genes found")
  write.csv(data.frame(), file.path(DATA_DIR, "WGCNA_hub_genes.csv"), row.names = FALSE)
}

# ===========================================================================
# 8. 
# ===========================================================================

all_summaries <- dplyr::bind_rows(
  if (!is.null(res_unsigned)) res_unsigned$summary else NULL,
  if (!is.null(res_signed)) res_signed$summary else NULL
)
write.csv(all_summaries, file.path(DATA_DIR, "WGCNA_summary.csv"), row.names = FALSE)

cat("\n================================================================\n")
cat(sprintf("  WGCNA — %s — COMPLETE\n", PARAMS$diff_type))
cat("================================================================\n")
if (!is.null(res_unsigned)) {
  cat(sprintf("  Unsigned: power=%d (%s), R²=%.3f, Conn=%.1f, %d modules, %.1f%% grey, %d hubs\n",
              res_unsigned$power_selected, res_unsigned$power_source,
              res_unsigned$summary$R2_at_power, res_unsigned$summary$mean_conn,
              res_unsigned$summary$n_modules, res_unsigned$summary$pct_grey,
              res_unsigned$summary$n_hub_genes))
}
if (!is.null(res_signed)) {
  cat(sprintf("  Signed:   power=%d (%s), R²=%.3f, Conn=%.1f, %d modules, %.1f%% grey, %d hubs\n",
              res_signed$power_selected, res_signed$power_source,
              res_signed$summary$R2_at_power, res_signed$summary$mean_conn,
              res_signed$summary$n_modules, res_signed$summary$pct_grey,
              res_signed$summary$n_hub_genes))
}
cat("================================================================\n")

log_step("08_WGCNA", "Step 08 COMPLETE")
