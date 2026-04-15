#!/usr/bin/env Rscript
# ==============================================================================
# 09A_candidate_pool.R — Candidate gene pool construction（）
#
# can ：
# maSigPro ∩ WGCNA ， times ：
# 1. maSigPro item PFilter（C）
# 2. Effect size filtering（DESeq2 lfcThreshold=1）
# candidate_pool.rds09B_ML and 09_PPI
#
# Literature basis:
# - maSigPro item Filter：Conesa et al. (2006) Bioinformatics; BioBam
#   - Effect size filtering：Love et al. (2014) Genome Biology, "Specifying minimum effect size"
# - WGCNARun not Filter：Ribas et al. (2022) BMC Biology
#
# Run：04_DEG → 06_maSigPro → 08_WGCNA → [09A] → 09B_ML → 09_PPI → 10
# ==============================================================================

source(file.path(SCRIPT_DIR, "00_setup.R"))

log_step("09A_POOL", "=== Building candidate gene pool ===")

# above DependenciesCheck：deg_results and wgcna_resultscandidate_poolUpdate
check_upstream("09A_POOL",
  upstream_files = c(
    "04_DEG"      = FILES$deg_results,
    "06_maSigPro" = FILES$masigpro_results,
    "08_WGCNA"    = FILES$wgcna_results
  ),
  output_files = c(
    file.path(DATA_DIR, "candidate_pool.rds")
  )
)

# ============================================================================
# 1. Load above 
# ============================================================================
all_results <- readRDS(FILES$deg_results)
wgcna_results <- readRDS(FILES$wgcna_results)
masigpro_results <- tryCatch(
  readRDS(FILES$masigpro_results),
  error = function(e) {
    log_step("09A_POOL", sprintf("ERROR: Cannot load maSigPro results: %s", e$message))
    stop("maSigPro results required for candidate pool construction")
  }
)

# ============================================================================
# 2. maSigPro Gene（R²≥0.7）
# ============================================================================
masigpro_ensembl <- character(0)
if (!is.null(masigpro_results$gene_clusters)) {
  masigpro_ensembl <- unique(masigpro_results$gene_clusters$ensembl_id[
    !is.na(masigpro_results$gene_clusters$ensembl_id)])
} else if (!is.null(masigpro_results$sig_genes_all)) {
  masigpro_ensembl <- unique(unlist(masigpro_results$sig_genes_all))
}
n_masigpro_raw <- length(masigpro_ensembl)
log_step("09A_POOL", sprintf("maSigPro (R\u00b2\u2265%.1f): %d genes", PARAMS$masigpro_rsq, n_masigpro_raw))

# ============================================================================
# 3. C：maSigPro item PFilter
# Time×Group item Significant(P<0.05) Gene
# " group should not Significant" Gene
# ============================================================================
n_masigpro_interaction <- n_masigpro_raw

if (!is.null(masigpro_results$tstep) && !is.null(masigpro_results$tstep$sol)) {
  sol <- masigpro_results$tstep$sol
  all_p_cols <- colnames(sol)[grep("^p\\.valor_", colnames(sol))]
  
 # item columns："x"（ item ）"Group" or "Induced"
  interaction_p_cols <- all_p_cols[grepl("x", all_p_cols, ignore.case = FALSE) & 
                                   grepl("Group|Induced", all_p_cols, ignore.case = TRUE)]
  
 # ： to "x" 
  if (length(interaction_p_cols) == 0) {
    interaction_p_cols <- all_p_cols[grepl("Group|Induced", all_p_cols, ignore.case = TRUE)]
  }
  
  log_step("09A_POOL", sprintf("  Interaction P-value columns: %s", 
                                paste(interaction_p_cols, collapse = ", ")))
  
  if (length(interaction_p_cols) > 0) {
    interaction_p_mat <- sol[, interaction_p_cols, drop = FALSE]
    min_interaction_p <- apply(interaction_p_mat, 1, function(x) {
      valid <- x[!is.na(x)]
      if (length(valid) == 0) return(1)
      min(valid)
    })
    
    interaction_sig_genes <- rownames(sol)[min_interaction_p < 0.05]
    
 # sol can can is ensembl or symbol， need need to 
    if (length(intersect(interaction_sig_genes, masigpro_ensembl)) > 0) {
      masigpro_ensembl <- intersect(masigpro_ensembl, interaction_sig_genes)
    } else if (!is.null(masigpro_results$gene_clusters)) {
      gc <- masigpro_results$gene_clusters
      sig_ens <- gc$ensembl_id[gc$symbol %in% interaction_sig_genes & !is.na(gc$ensembl_id)]
      if (length(sig_ens) > 0) {
        masigpro_ensembl <- intersect(masigpro_ensembl, sig_ens)
      } else {
        log_step("09A_POOL", "WARNING: Could not map interaction sig genes, skipping filter")
      }
    }
    
    n_masigpro_interaction <- length(masigpro_ensembl)
    log_step("09A_POOL", sprintf("  Interaction filter: %d → %d (removed %d)",
                                  n_masigpro_raw, n_masigpro_interaction,
                                  n_masigpro_raw - n_masigpro_interaction))
  }
} else {
  log_step("09A_POOL", "WARNING: tstep$sol not available, skipping interaction filter")
}

# ============================================================================
# 4. WGCNA Key modulesGene
# ============================================================================
wgcna_key_genes <- character(0)
tryCatch({
  res <- NULL
  if (!is.null(wgcna_results$unsigned)) res <- wgcna_results$unsigned
  else if (!is.null(wgcna_results$signed)) res <- wgcna_results$signed
  
  if (!is.null(res)) {
    mod_colors <- res$moduleColors
    gene_ids <- names(mod_colors)
    if (!is.null(res$sig_modules)) {
 # greyModule（ME0 or MEgrey）
      sig_mods_clean <- res$sig_modules[!grepl("^ME0$|^MEgrey$", res$sig_modules)]
      log_step("09A_POOL", sprintf("  sig_modules from 08: %d (after grey exclusion: %d)",
                                    length(res$sig_modules), length(sig_mods_clean)))
      sig_labels <- as.integer(gsub("^ME", "", sig_mods_clean))
      sig_colors <- WGCNA::labels2colors(sig_labels)
      wgcna_key_genes <- gene_ids[mod_colors %in% sig_colors]
    }
  }
}, error = function(e) log_step("09A_POOL", sprintf("WGCNA: %s", e$message)))

log_step("09A_POOL", sprintf("WGCNA key module genes: %d", length(wgcna_key_genes)))

# ============================================================================
# 5. maSigPro(filtered) ∩ WGCNA
# ============================================================================
candidate_pool <- intersect(masigpro_ensembl, wgcna_key_genes)
n_pool_pre_effect <- length(candidate_pool)
log_step("09A_POOL", sprintf("maSigPro(filtered) ∩ WGCNA: %d", n_pool_pre_effect))

# ============================================================================
# 6. Effect size filtering： Timepoint lfcThreshold=1 Significant
# : Love et al. (2014) Genome Biology
# H₀: |LFC| ≤ 1 → padj_lfc < 0.05 has |LFC| > 1
# ============================================================================
n_pool_post_effect <- n_pool_pre_effect
effect_filter_applied <- FALSE

if (!is.null(all_results$wald_by_time_lfc)) {
  lfc_sig_genes <- character(0)
  for (tp in names(all_results$wald_by_time_lfc)) {
    lfc_df <- all_results$wald_by_time_lfc[[tp]]
    sig_ids <- lfc_df$ensembl_id[!is.na(lfc_df$padj_lfc) & lfc_df$padj_lfc < 0.05]
    lfc_sig_genes <- union(lfc_sig_genes, sig_ids)
  }
  
  log_step("09A_POOL", sprintf("  Genes with |LFC|>1 at any timepoint: %d (of %d in transcriptome)",
                                length(lfc_sig_genes), nrow(all_results$lrt_interaction)))
  
  candidate_pool <- intersect(candidate_pool, lfc_sig_genes)
  n_pool_post_effect <- length(candidate_pool)
  effect_filter_applied <- TRUE
  
  log_step("09A_POOL", sprintf("  Effect size filter: %d → %d (removed %d)",
                                n_pool_pre_effect, n_pool_post_effect,
                                n_pool_pre_effect - n_pool_post_effect))
} else {
  log_step("09A_POOL", "WARNING: wald_by_time_lfc not found — run 04_DEG with lfcThreshold first")
  log_step("09A_POOL", "Proceeding WITHOUT effect size filter")
}

# ============================================================================
# 7. Save candidate pool
# ============================================================================
pool_result <- list(
  candidate_pool         = candidate_pool,
  n_masigpro_raw         = n_masigpro_raw,
  n_masigpro_interaction = n_masigpro_interaction,
  n_wgcna_key            = length(wgcna_key_genes),
  n_pool_pre_effect      = n_pool_pre_effect,
  n_pool_post_effect     = n_pool_post_effect,
  effect_filter_applied  = effect_filter_applied,
 # SaveIntermediateGene10_integration using 
  masigpro_ensembl_filtered = masigpro_ensembl,
  wgcna_key_genes        = wgcna_key_genes
)

pool_file <- file.path(DATA_DIR, "candidate_pool.rds")
saveRDS(pool_result, pool_file)

# ============================================================================
# 
# ============================================================================
cat("\n================================================================\n")
cat(sprintf("  %s — CANDIDATE POOL CONSTRUCTED\n", PARAMS$diff_type))
cat("================================================================\n")
cat(sprintf("  maSigPro (R\u00b2\u2265%.1f):              %d\n", PARAMS$masigpro_rsq, n_masigpro_raw))
cat(sprintf("  maSigPro (interaction P<0.05):   %d\n", n_masigpro_interaction))
cat(sprintf("  WGCNA key modules:               %d\n", length(wgcna_key_genes)))
cat(sprintf("  maSigPro(filtered) ∩ WGCNA:      %d\n", n_pool_pre_effect))
if (effect_filter_applied) {
  cat(sprintf("  Effect size (|LFC|>1 test):      %d\n", n_pool_post_effect))
}
cat(sprintf("  *** FINAL CANDIDATE POOL:        %d ***\n", length(candidate_pool)))
cat("================================================================\n")

log_step("09A_POOL", sprintf("Step 09A COMPLETE — candidate pool: %d genes saved to %s",
                              length(candidate_pool), pool_file))
