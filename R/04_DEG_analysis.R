#!/usr/bin/env Rscript
# ==============================================================================
# 04_DEG_analysis.R — DifferentialExpression
# ：
# 1. LRT Treatment×Time should （）
# 2. LRT Treatment should （ has Timepoint should ）
# 3. Timepoint Wald Induced vs Control（pair-wise）
# ： using raw counts， DESeq2 within Normalization
# ==============================================================================

source(file.path(SCRIPT_DIR, "00_setup.R"))
source(file.path(SCRIPT_DIR, "theme_bindlab.R"))

log_step("04_DEG", "Starting differential expression analysis...")

# Load
dds         <- readRDS(FILES$dds_object)
sample_info <- readRDS(FILES$sample_info)
gene_anno   <- readRDS(FILES$gene_annotation)

# ============================================================================
# 1. DESeq2 （）
# ============================================================================
# : ~ Treatment + Time + Treatment:Time
# Reference: Control, 4d

# --- 1a. LRT: Treatment×Time should ---
# full model:  ~ Treatment + Time + Treatment:Time
# reduced:     ~ Treatment + Time
# Testing issue: is ?
log_step("04_DEG", "LRT for Treatment:Time interaction...")

dds_lrt_interaction <- DESeq(dds, test = "LRT", 
                              reduced = ~ Treatment + Time)
res_lrt_interaction <- results(dds_lrt_interaction, alpha = PARAMS$padj_cutoff)

n_sig_interaction <- sum(res_lrt_interaction$padj < PARAMS$padj_cutoff, na.rm = TRUE)
log_step("04_DEG", sprintf("LRT Interaction: %d genes with padj < %.2f",
                            n_sig_interaction, PARAMS$padj_cutoff))

# --- 1b. LRT: Treatment should ---
# reduced: ~ Time
# Testing issue: in should after ， is has Significant?
log_step("04_DEG", "LRT for Treatment effect (controlling for Time)...")

dds_lrt_treatment <- DESeq(dds, test = "LRT",
                            reduced = ~ Time)
res_lrt_treatment <- results(dds_lrt_treatment, alpha = PARAMS$padj_cutoff)

n_sig_treatment <- sum(res_lrt_treatment$padj < PARAMS$padj_cutoff, na.rm = TRUE)
log_step("04_DEG", sprintf("LRT Treatment: %d genes with padj < %.2f",
                            n_sig_treatment, PARAMS$padj_cutoff))

# ============================================================================
# 2. Wald：Timepoint Induced vs Control
# ============================================================================
log_step("04_DEG", "Wald tests for each timepoint...")

# using group Timepoint Comparison
dds_group <- dds
dds_group$group <- factor(paste(dds_group$Treatment, dds_group$Time, sep = "_"))
design(dds_group) <- ~ group
dds_group <- DESeq(dds_group)

timepoints <- PARAMS$time_labels
deg_by_time <- list()
deg_by_time_lfc <- list()  # Effect size test result (lfcThreshold=1）

for (tp in timepoints) {
  contrast_name <- paste0("Induced_", tp, "_vs_Control_", tp)
  contrast_vec <- c("group", paste0("Induced_", tp), paste0("Control_", tp))
  
 # --- StandardWald（H₀: LFC=0）： using LFC and Visualization ---
  res <- results(dds_group, 
                 contrast = contrast_vec,
                 alpha = PARAMS$padj_cutoff)
  
  # AddGeneAnnotation
  res_df <- as.data.frame(res) %>%
    rownames_to_column("ensembl_id") %>%
    left_join(gene_anno, by = c("ensembl_id" = "ensembl_gene_id")) %>%
    arrange(padj) %>%
    mutate(
      timepoint = tp,
      regulation = case_when(
        padj < PARAMS$padj_cutoff & log2FoldChange > PARAMS$lfc_cutoff ~ "Up",
        padj < PARAMS$padj_cutoff & log2FoldChange < -PARAMS$lfc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  deg_by_time[[tp]] <- res_df
  
 # --- Effect size test（H₀: |LFC|≤1）： using below Effect size filtering ---
 # : Love et al. (2014) Genome Biology, "Specifying minimum effect size"
 # lfcThreshold=1 " is 2"
  res_lfc <- results(dds_group,
                     contrast = contrast_vec,
                     alpha = PARAMS$padj_cutoff,
                     lfcThreshold = PARAMS$lfc_cutoff)
  
  res_lfc_df <- as.data.frame(res_lfc) %>%
    rownames_to_column("ensembl_id") %>%
    dplyr::select(ensembl_id, padj_lfc = padj)
  
  deg_by_time_lfc[[tp]] <- res_lfc_df
  
  n_up   <- sum(res_df$regulation == "Up", na.rm = TRUE)
  n_down <- sum(res_df$regulation == "Down", na.rm = TRUE)
  n_lfc_sig <- sum(res_lfc_df$padj_lfc < PARAMS$padj_cutoff, na.rm = TRUE)
  log_step("04_DEG", sprintf("  %s: Up=%d, Down=%d (|log2FC|>%.1f, padj<%.2f) | lfcThreshold test sig: %d",
                              tp, n_up, n_down, PARAMS$lfc_cutoff, PARAMS$padj_cutoff, n_lfc_sig))
}

# ============================================================================
# 3. ★★ 8 group — SpecificityPlanned Contrasts
# ============================================================================
# : ~ 0 + group（ no ， group = group Expression）
# : Gene in Adipo_Xd7 group Significant/?
# Reference: Ekiz (2024) DESeq2 design guide; DESeq2 vignette §3.3

log_step("04_DEG", "=== 8-Group Model: Stage-Specific Planned Contrasts ===")

# no 
dds_8g <- dds
dds_8g$group <- factor(paste(dds_8g$Treatment, dds_8g$Time, sep = "_"))
design(dds_8g) <- ~ 0 + group
dds_8g <- DESeq(dds_8g)

log_step("04_DEG", sprintf("Model coefficients: %s", 
                            paste(resultsNames(dds_8g), collapse = ", ")))

# " group vs 7 group " contrasts
# no ~ 0 + group:
# Adipo4 vs rest = groupInduced_4d - (sum of other 7 groups)/7
group_names <- levels(dds_8g$group)  # 8 levels
target_groups <- c("Induced_4d", "Induced_7d", "Induced_14d", "Induced_21d")

stage_contrasts <- list()
stage_deg <- list()

for (target in target_groups) {
 # contrast： group =1，7 group =-1/7
  contrast_vec <- rep(0, length(group_names))
  names(contrast_vec) <- paste0("group", group_names)
  contrast_vec[paste0("group", target)] <- 1
  other_groups <- setdiff(group_names, target)
  for (og in other_groups) {
    contrast_vec[paste0("group", og)] <- -1/7
  }
  
 # 
  res_stage <- results(dds_8g, contrast = contrast_vec, alpha = PARAMS$padj_cutoff)
  
  res_df <- as.data.frame(res_stage) %>%
    tibble::rownames_to_column("ensembl_id") %>%
    dplyr::left_join(gene_anno, by = c("ensembl_id" = "ensembl_gene_id")) %>%
    dplyr::arrange(padj) %>%
    dplyr::mutate(
      target_group = target,
      regulation = dplyr::case_when(
        padj < PARAMS$padj_cutoff & log2FoldChange > PARAMS$lfc_cutoff ~ "Up",
        padj < PARAMS$padj_cutoff & log2FoldChange < -PARAMS$lfc_cutoff ~ "Down",
        TRUE ~ "NS"
      )
    )
  
  stage_contrasts[[target]] <- contrast_vec
  stage_deg[[target]] <- res_df
  
  n_up   <- sum(res_df$regulation == "Up", na.rm = TRUE)
  n_down <- sum(res_df$regulation == "Down", na.rm = TRUE)
  log_step("04_DEG", sprintf("  %s vs rest(7): Up=%d, Down=%d", target, n_up, n_down))
}

# "Adipo_Xd vs Ctrl_Xd" （ already in above Calculate）

# ============================================================================
# 4. 
# ============================================================================

# LRT should Gene（ using WGCNA）
lrt_interaction_df <- as.data.frame(res_lrt_interaction) %>%
  rownames_to_column("ensembl_id") %>%
  left_join(gene_anno, by = c("ensembl_id" = "ensembl_gene_id")) %>%
  arrange(padj)

lrt_interaction_sig <- lrt_interaction_df %>%
  filter(padj < PARAMS$lrt_padj)

log_step("04_DEG", sprintf("LRT interaction significant genes for WGCNA: %d", 
                            nrow(lrt_interaction_sig)))

# LRT should Gene
lrt_treatment_df <- as.data.frame(res_lrt_treatment) %>%
  rownames_to_column("ensembl_id") %>%
  left_join(gene_anno, by = c("ensembl_id" = "ensembl_gene_id")) %>%
  arrange(padj)

# Merge has Wald
deg_all_timepoints <- bind_rows(deg_by_time)

# DEGStatistics
deg_summary <- deg_all_timepoints %>%
  filter(regulation != "NS") %>%
  group_by(timepoint, regulation) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = regulation, values_from = n, values_fill = 0)

log_step("04_DEG", "DEG Summary by timepoint:")
print(deg_summary)

# ============================================================================
# 5. Save results
# ============================================================================
all_results <- list(
  lrt_interaction    = lrt_interaction_df,
  lrt_interaction_sig = lrt_interaction_sig,
  lrt_treatment      = lrt_treatment_df,
  wald_by_time       = deg_by_time,
  wald_by_time_lfc   = deg_by_time_lfc,  # Effect size test result (lfcThreshold=1）
  wald_all           = deg_all_timepoints,
  deg_summary        = deg_summary,
  dds_group          = dds_group
)

# 8 group （Save，）
group_deg_results <- list(
  dds_8g          = dds_8g,
  stage_deg       = stage_deg,
  stage_contrasts = stage_contrasts
)

save_data(all_results, FILES$deg_results)
save_data(group_deg_results, FILES$group_deg_results)
save_data(lrt_interaction_df, FILES$lrt_results)

# Export is CSV（ after Reference/）
write.csv(lrt_interaction_sig, 
          file.path(DATA_DIR, "LRT_interaction_significant_genes.csv"),
          row.names = FALSE)

for (tp in timepoints) {
  sig_genes <- deg_by_time[[tp]] %>% dplyr::filter(regulation != "NS")
  write.csv(sig_genes,
            file.path(DATA_DIR, sprintf("DEG_%s_Induced_vs_Control.csv", tp)),
            row.names = FALSE)
}

# Export8 group DEG
for (target in target_groups) {
  sig_stage <- stage_deg[[target]] %>% dplyr::filter(regulation != "NS")
  write.csv(sig_stage,
            file.path(DATA_DIR, sprintf("DEG_%s_vs_rest7.csv", target)),
            row.names = FALSE)
}

log_step("04_DEG", "Step 04 COMPLETE")
