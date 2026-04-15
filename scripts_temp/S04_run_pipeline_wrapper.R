#!/usr/bin/env Rscript
# ==============================================================================
# S04_run_pipeline_wrapper.R — GEO dataset pipeline adaptation wrapper

# ---- Disable Browse debug mode on error ----
options(error = NULL)
#
# Purpose: Run external GEO datasets through the pipeline without modifying original scripts.
#       The original pipeline scripts contain hardcoded timepoint labels, factor levels,
#       color schemes, maSigPro parameters, etc.
#       This wrapper solves this via a "source-then-override" strategy.
#
# Adaptation strategy:
#   Level 1: source 00_setup.R then override PARAMS / COLORS
#   Level 2: source 01_data_import.R then fix sample_info Time/Group factors
#   Level 3: Fix hardcoded timepoint variables after sourcing each downstream script
#   Level 4: Gene IDadaptation (symbol vs Ensembl)
#
# Usage:
#   1. Select DATASET_ID (see registry below)
#   2. source("S04_run_pipeline_wrapper.R") or run per phase
#   3. No modification to original pipeline code
#
# Key principles:
#   - Original pipeline scripts are read-only
#   - All adaptation is done in this file
#   - Methods paper uses independent validation datasets
# ==============================================================================

# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  STEP 0: STEP 0: Select dataset + load config                                             ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ---- Select dataset to run ----
# Option 1: Set before sourcing DATASET_ID <- "GSE197067"
# Option 2: Modify the default below
if (!exists("DATASET_ID")) {
  DATASET_ID <- "GSE307424"   # Default: Lung SMARCA2, 18 samples, 3 timepoints
}

# ---- Dataset registry (adaptation parameters per dataset) ----
DATASET_PROFILES <- list(

  GSE307424 = list(
    project_dir   = file.path(METHODS_BASE, "pipeline_runs", "GEO_GSE307424_Lung"),
    prefix        = "Lung",
    # Time configuration
    time_values   = c(6, 48, 72),              # Numeric time values (hours)
    time_labels   = c("6h", "48h", "72h"),     # Display labels
    time_suffix   = "d",                       # Suffix in sample names（Lung6d1 → keep "d"）
    time_ref      = "6h",                      # DESeq2 reference level (earliest timepoint)
    # maSigPro
    masigpro_degree  = 2,   # 3timepoints → degree = nTP - 1 = 2
    masigpro_k       = 6,   # Fewer clusters (fewer genes)
    masigpro_min_obs = 3,   # p.vector()default value; User's Guide example uses 20 but requires <= n_samples
    # Gene ID
    gene_id_type  = "symbol",  # Non-Ensembl!
    # WGCNA
    wgcna_note    = "18 samples → FAQ power table fallback likely",
    # Visualization
    colors_timepoint = c("6h"  = "#3C5488",
                         "48h" = "#F39B7F",
                         "72h" = "#E64B35"),
    colors_group = c("Induced_6h"  = "#DC0000", "Control_6h"  = "#3C5488",
                     "Induced_48h" = "#E64B35", "Control_48h" = "#00A087",
                     "Induced_72h" = "#B09C85", "Control_72h" = "#7E6148"),
    shape_timepoint = c("6h" = 16, "48h" = 17, "72h" = 15)
  ),

  GSE197067 = list(
    project_dir   = file.path(METHODS_BASE, "pipeline_runs", "GEO_GSE197067_Tcell"),
    prefix        = "Tcell",
    # Time configuration（Exclude 0h: 0h has only Control, no Induced, interaction term not estimable）
    time_values   = c(6, 12, 24, 48, 72),       # Exclude 0h
    time_labels   = c("6h", "12h", "24h", "48h", "72h"),
    time_suffix   = "d",
    time_ref      = "6h",                        # Earliest timepoint as reference
    # Sample name pattern to exclude 0h（data_rawTcell0dC1-4 samples need to be removed during import）
    exclude_samples_pattern = "^Tcell0d",        # Matches Tcell0dC1 etc.
    # maSigPro
    masigpro_degree  = 4,   # 5timepoints → degree = nTP - 1 = 4
    masigpro_k       = 9,
    masigpro_min_obs = 20,  # 40samples, following User's Guide recommendation
    # Gene ID
    gene_id_type  = "ensembl",
    # WGCNA
    wgcna_note    = "40 samples → adequate for WGCNA",
    # Visualization
    colors_timepoint = c("6h"  = "#3C5488",
                         "12h" = "#00A087",
                         "24h" = "#F39B7F",
                         "48h" = "#E64B35",
                         "72h" = "#B09C85"),
    colors_group = c("Induced_6h"   = "#DC0000", "Control_6h"   = "#3C5488",
                     "Induced_12h"  = "#F39B7F", "Control_12h"  = "#4DBBD5",
                     "Induced_24h"  = "#E64B35", "Control_24h"  = "#00A087",
                     "Induced_48h"  = "#B09C85", "Control_48h"  = "#7E6148",
                     "Induced_72h"  = "#DC9FB4", "Control_72h"  = "#7570B3"),
    shape_timepoint = c("6h" = 16, "12h" = 17, "24h" = 15, "48h" = 18, "72h" = 8)
  )
)

# ---- Load current dataset's profile ----
PROFILE <- DATASET_PROFILES[[DATASET_ID]]
if (is.null(PROFILE)) stop("Unknown DATASET_ID: ", DATASET_ID)

cat("================================================================\n")
cat("  METI-FS Pipeline Wrapper\n")
cat(sprintf("  Dataset: %s\n", DATASET_ID))
cat(sprintf("  Project: %s\n", PROFILE$project_dir))
cat(sprintf("  Times:   %s\n", paste(PROFILE$time_labels, collapse = ", ")))
cat(sprintf("  Degree:  %d\n", PROFILE$masigpro_degree))
cat(sprintf("  GeneID:  %s\n", PROFILE$gene_id_type))
cat("================================================================\n")


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  HELPER: Common fix functions                                                      ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

#' Fix Time and Group factor levels in sample_info
#' Call after each pipeline script source (some scripts reload sample_info from files)
fix_sample_info_factors <- function(si) {
  actual_time_labels <- PROFILE$time_labels   # c("6h","48h","72h")
  actual_time_values <- PROFILE$time_values   # c(6, 48, 72)
  
  # 01_data_import.R parses Lung6d1 as time_num=6, time_label="6d"
  # but actual time unit may be hours (h) not days (d)
  # Use time_num -> PROFILE$time_labels to build mapping table
  time_map <- setNames(actual_time_labels, as.character(actual_time_values))
  # e.g. "6" → "6h", "48" → "48h", "72" → "72h"
  
  # Rebuild time_label
  si$time_label <- time_map[as.character(si$time_num)]
  
  # Fix Time factor
  si$Time <- factor(si$time_label, levels = actual_time_labels)
  
  # Fix Group factor
  group_levels_ordered <- c()
  for (tl in actual_time_labels) {
    group_levels_ordered <- c(group_levels_ordered,
                              paste0("Control_", tl),
                              paste0("Induced_", tl))
  }
  si$Group <- factor(paste(si$treatment, si$time_label, sep = "_"),
                     levels = group_levels_ordered)
  
  return(si)
}

#' Override COLORS timepoint and group mapping
fix_colors <- function() {
  COLORS$timepoint <<- PROFILE$colors_timepoint
  COLORS$group     <<- PROFILE$colors_group
}

#' Override maSigPro parameters in PARAMS
fix_params <- function() {
  PARAMS$masigpro_degree <<- PROFILE$masigpro_degree
  PARAMS$masigpro_k      <<- PROFILE$masigpro_k
}

#' Fix prefix issue: create file aliases expected by 00_setup.R under data_raw
#' 
#' Root cause：00_setup.Ruses basename(PROJECT_DIR) as prefix → "GEO_GSE307424_Lung"
#'        but S03b-created files use "Lung" prefix
#'        and each pipeline script re-sources at the beginning 00_setup.R，making overrides ineffective
#' 
#' Solution: create copies from expected names to actual names under data_raw (file.copy)
#'        This way FILES paths resolve regardless of how many times 00_setup.R is sourced
#'        Only need to call once (at run_phase1 start), files persist afterwards
fix_prefix <- function() {
  expected_prefix <- basename(PROFILE$project_dir)  # "GEO_GSE307424_Lung"
  actual_prefix   <- PROFILE$prefix                  # "Lung"
  
  if (expected_prefix == actual_prefix) {
    cat("  [PREFIX] No alias needed, prefixes match\n")
    return(invisible(TRUE))
  }
  
  raw_dir <- file.path(PROFILE$project_dir, "data_raw")
  
  # 4 files that need aliases
  suffixes <- c("_all_counts_with_order.tsv",
                "_all_tpm.tsv",
                "_metadata.tsv",
                "_gene_annotation.tsv")
  
  for (sfx in suffixes) {
    actual_file   <- file.path(raw_dir, paste0(actual_prefix, sfx))
    expected_file <- file.path(raw_dir, paste0(expected_prefix, sfx))
    
    if (file.exists(expected_file)) {
      # Alias already exists, skip
      next
    }
    
    if (!file.exists(actual_file)) {
      cat(sprintf("  [WARN] Source file missing: %s\n", basename(actual_file)))
      next
    }
    
    # Create copy (Windows symlink unreliable, file.copy is safer)
    file.copy(actual_file, expected_file, overwrite = FALSE)
    cat(sprintf("  [ALIAS] %s → %s\n", basename(actual_file), basename(expected_file)))
  }
  
  # Sync lineage display variables (these won't be overridden by re-source, as 00_setup uses diff_type)
  # But diff_type itself resets on each re-source, so we accept it as"GEO_GSE307424_Lung"
  # As long as FILES paths can find the files
  
  cat(sprintf("  [PREFIX] File aliases created: '%s_*' → '%s_*'\n",
              actual_prefix, expected_prefix))
  invisible(TRUE)
}

#' Get timepoint labels for current datasetvector (replaces hardcodedc(<original_timepoints>)）
get_timepoints <- function() {
  PROFILE$time_labels
}

#' Get numeric time vector for current dataset
get_time_values <- function() {
  PROFILE$time_values
}

#' Get target groups for current dataset's Induced arm
get_target_groups <- function() {
  paste0("Induced_", PROFILE$time_labels)
}

#' Ensure global environment variables are set (supports calling any phase independently)
ensure_env <- function() {
  options(error = NULL)  # Prevent Browse debug mode
  PROJECT_DIR <<- PROFILE$project_dir
  SCRIPT_DIR  <<- PIPELINE_SCRIPTS
  source(file.path(SCRIPT_DIR, "00_setup.R"), local = FALSE)
  fix_params()
  fix_prefix()
  fix_colors()
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 1: 00_setup.R + 01_data_import.R + 02 + 03                         ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_phase1 <- function() {
  cat("\n### PHASE 1: Data import + preprocessing + normalization ###\n\n")
  
  # ---- 1a. source 00_setup.R ----
  PROJECT_DIR <<- PROFILE$project_dir
  SCRIPT_DIR  <<- PIPELINE_SCRIPTS
  source(file.path(SCRIPT_DIR, "00_setup.R"), local = FALSE)
  
  # Immediately override PARAMS and COLORS
  fix_params()
  fix_prefix()
  fix_colors()
  
  cat(sprintf("  [OVERRIDE] masigpro_degree = %d (was 3)\n", PARAMS$masigpro_degree))
  cat(sprintf("  [OVERRIDE] masigpro_k = %d (was 9)\n", PARAMS$masigpro_k))
  cat(sprintf("  [OVERRIDE] COLORS$timepoint: %s\n",
              paste(names(COLORS$timepoint), collapse = ", ")))
  
  # ---- 1b. source 01_data_import.R ----
  # This script parses sample names and sets Time/Group factor levels
  # Issue: L77 hardcodes levels = c(<original_timepoints>)
  # Strategy: fix immediately after source
  source(file.path(SCRIPT_DIR, "01_data_import.R"), local = FALSE)
  
  # ---- Exclude samples from incomplete timepoints (e.g., GSE197067 0h has only Control) ----
  if (!is.null(PROFILE$exclude_samples_pattern)) {
    excl_pat <- PROFILE$exclude_samples_pattern
    excl_idx <- grepl(excl_pat, colnames(counts_raw))
    n_excl <- sum(excl_idx)
    if (n_excl > 0) {
      excl_names <- colnames(counts_raw)[excl_idx]
      cat(sprintf("\n  [EXCLUDE] Removing %d samples matching '%s': %s\n",
                  n_excl, excl_pat, paste(excl_names, collapse = ", ")))
      
      # Remove from counts and tpm matrices
      counts_raw <<- counts_raw[, !excl_idx, drop = FALSE]
      tpm_raw    <<- tpm_raw[, !excl_idx, drop = FALSE]
      
      # Remove from sample_info
      sample_info <<- sample_info[!excl_idx, , drop = FALSE]
      
      cat(sprintf("  [EXCLUDE] Remaining: %d samples\n", ncol(counts_raw)))
    }
  }
  
  # Fix Time and Group factor levels
  sample_info <<- fix_sample_info_factors(sample_info)
  
  # Verify fix results
  cat("\n  [FIX] sample_info$Time levels after correction:\n")
  print(table(sample_info$Time))
  cat("  [FIX] sample_info$Group levels after correction:\n")
  print(table(sample_info$Group))
  
  # Check for NAs
  if (any(is.na(sample_info$Time))) {
    stop("FATAL: sample_info$Time contains NA after fix! Check time_labels in PROFILE.")
  }
  
  # Re-save corrected sample_info
  save_data(sample_info, FILES$sample_info)
  
  # If samples were excluded, must write trimmed counts/tpm back to disk
  # Otherwise 02_preprocessing reads original column count from FILES$counts_raw
  if (!is.null(PROFILE$exclude_samples_pattern)) {
    cat("  [WRITE] Overwriting raw data files with excluded samples removed...\n")
    counts_df <- data.frame(gene_id = rownames(counts_raw), counts_raw, check.names = FALSE)
    write.table(counts_df, FILES$counts_raw, sep = "\t", row.names = FALSE, quote = FALSE)
    tpm_df <- data.frame(gene_id = rownames(tpm_raw), tpm_raw, check.names = FALSE)
    write.table(tpm_df, FILES$tpm_raw, sep = "\t", row.names = FALSE, quote = FALSE)
    cat(sprintf("  [WRITE] %s: %d genes × %d samples\n",
                basename(FILES$counts_raw), nrow(counts_raw), ncol(counts_raw)))
  }
  
  # ---- 1c. Gene ID adaptation ----
  # If dataset uses gene symbol (not Ensembl), fix gene_annotation
  if (PROFILE$gene_id_type == "symbol") {
    cat("\n  [ADAPT] Gene ID type = symbol (not Ensembl)\n")
    # gene_anno already created by 01_data_import.R
    # But it tries to query org.Hs.eg.db with ensembl_gene_id; treating symbols as Ensembl will fail
    # Fix: remap using SYMBOL as keytype
    
    gene_anno_fixed <- data.frame(
      ensembl_gene_id = rownames(counts_raw),  # actually gene symbol
      hgnc_symbol     = rownames(counts_raw),  # also gene symbol
      stringsAsFactors = FALSE
    )
    
    # Try to get Entrez ID (using SYMBOL keytype)
    tryCatch({
      sym2entrez <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = unique(gene_anno_fixed$hgnc_symbol),
        columns = c("ENTREZID"),
        keytype = "SYMBOL"
      )
      sym2entrez <- sym2entrez[!duplicated(sym2entrez$SYMBOL), ]
      gene_anno_fixed$entrez_id <- sym2entrez$ENTREZID[
        match(gene_anno_fixed$hgnc_symbol, sym2entrez$SYMBOL)]
      
      n_mapped <- sum(!is.na(gene_anno_fixed$entrez_id))
      cat(sprintf("  [ADAPT] Symbol->Entrez mapped: %d/%d (%.1f%%)\n",
                  n_mapped, nrow(gene_anno_fixed),
                  100 * n_mapped / nrow(gene_anno_fixed)))
    }, error = function(e) {
      cat(sprintf("  [WARN] Entrez mapping failed: %s\n", e$message))
      gene_anno_fixed$entrez_id <<- NA
    })
    
    gene_anno <<- gene_anno_fixed
    save_data(gene_anno, FILES$gene_annotation)
    cat("  [ADAPT] gene_annotation overwritten with symbol-based mapping\n")
  }
  
  # ---- 1d. source 02_preprocessing.R ----
  source(file.path(SCRIPT_DIR, "02_preprocessing.R"), local = FALSE)
  
  # ---- 1e. 03_normalization_QC (ADAPTED) ----
  # Issue: L38 relevel(ref="4d") + L81/L110 shape hardcoded
  # Cannot source original script (relevel during source throws error), write adapted version
  cat("  [ADAPT] Running adapted normalization + QC (03)...\n")
  run_03_normalization_adapted()
  
  cat("\n### PHASE 1 COMPLETE ###\n")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  ADAPTED SCRIPT: 03_normalization_QC                                       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_03_normalization_adapted <- function() {
  ensure_env()
  
  log_step("03_NORM", "Starting normalization and QC (ADAPTED)...")
  
  counts_filtered <- readRDS(FILES$counts_filtered)
  sample_info     <- readRDS(FILES$sample_info)
  sample_info     <- fix_sample_info_factors(sample_info)
  
  # 1. DESeqDataSet
  counts_mat <- as.matrix(counts_filtered)
  storage.mode(counts_mat) <- "integer"
  stopifnot(all(colnames(counts_mat) == rownames(sample_info)))
  
  dds <- DESeqDataSetFromMatrix(
    countData = counts_mat,
    colData   = sample_info,
    design    = ~ Treatment + Time + Treatment:Time
  )
  dds$Treatment <- relevel(dds$Treatment, ref = "Control")
  dds$Time      <- relevel(dds$Time, ref = PROFILE$time_ref)  # ★ ADAPTED ★
  
  log_step("03_NORM", sprintf("DESeqDataSet: %d genes × %d samples, Time ref = %s",
                               nrow(dds), ncol(dds), levels(dds$Time)[1]))
  
  # 2. VST
  vsd <- vst(dds, blind = FALSE)
  vst_mat <- assay(vsd)
  log_step("03_NORM", "VST transformation complete")
  
  save_data(dds, FILES$dds_object)
  save_data(vst_mat, FILES$vst_matrix)
  
  # 3. PCA
  log_step("03_NORM", "PCA analysis...")
  ntop <- 500
  rv <- matrixStats::rowVars(vst_mat)
  select_genes <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca_data <- prcomp(t(vst_mat[select_genes, ]), center = TRUE, scale. = FALSE)
  
  pca_df <- as.data.frame(pca_data$x[, 1:min(5, ncol(pca_data$x))])
  pca_df$sample    <- rownames(pca_df)
  pca_df$Treatment <- sample_info$Treatment
  pca_df$Time      <- sample_info$Time
  pca_df$Group     <- sample_info$Group
  pct_var <- round(100 * (pca_data$sdev^2 / sum(pca_data$sdev^2)), 1)
  
  # Shape mapping (ADAPTED)
  shape_vals <- PROFILE$shape_timepoint
  
  p_pca_main <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Treatment, shape = Time), size = 3.5, stroke = 0.8) +
    scale_color_treatment() +
    scale_shape_manual(values = shape_vals) +  # ★ ADAPTED ★
    stat_ellipse(aes(group = Treatment, color = Treatment),
                 type = "t", level = 0.95, linetype = 2, linewidth = 0.5) +
    labs(title = "PCA of Gene Expression",
         subtitle = sprintf("Top %d variable genes, VST-transformed", ntop),
         x = sprintf("PC1 (%s%% variance)", pct_var[1]),
         y = sprintf("PC2 (%s%% variance)", pct_var[2])) +
    theme_bindlab_minimal()
  save_pub_fig(p_pca_main, "PCA_Treatment_Time", "02_PCA_Clustering", width = 8, height = 6)
  
  p_pca_group <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(color = Group), size = 3.5) +
    scale_color_manual(values = COLORS$group) +
    ggrepel::geom_text_repel(aes(label = sample), size = 2.2, max.overlaps = 20) +
    labs(title = "PCA — All Groups",
         x = sprintf("PC1 (%s%% variance)", pct_var[1]),
         y = sprintf("PC2 (%s%% variance)", pct_var[2])) +
    theme_bindlab_minimal()
  save_pub_fig(p_pca_group, "PCA_AllGroups_labeled", "02_PCA_Clustering", width = 10, height = 7)
  
  # Scree plot
  scree_df <- data.frame(PC = paste0("PC", 1:min(10, length(pct_var))),
                         Variance = pct_var[1:min(10, length(pct_var))])
  scree_df$PC <- factor(scree_df$PC, levels = scree_df$PC)
  p_scree <- ggplot(scree_df, aes(x = PC, y = Variance)) +
    geom_bar(stat = "identity", fill = "#3C5488", width = 0.6) +
    geom_line(aes(group = 1), color = "#E64B35", linewidth = 0.8) +
    geom_point(color = "#E64B35", size = 2) +
    geom_text(aes(label = paste0(Variance, "%")), vjust = -0.5, size = 3) +
    labs(title = "PCA Scree Plot", x = NULL, y = "Variance Explained (%)") +
    theme_bindlab() +
    expand_limits(y = max(scree_df$Variance) * 1.15)
  save_pub_fig(p_scree, "PCA_scree_plot", "02_PCA_Clustering", width = 8, height = 5)
  
  # 4. Sample correlation heatmap
  log_step("03_NORM", "Sample correlation heatmap...")
  cor_mat <- cor(vst_mat, method = "spearman")
  anno_col <- data.frame(
    Treatment = sample_info$Treatment,
    Time      = sample_info$Time,
    row.names = rownames(sample_info)
  )
  anno_colors <- list(
    Treatment = COLORS$treatment,
    Time      = COLORS$timepoint
  )
  
  save_heatmap_fig(
    draw_func = function() {
      print(pheatmap::pheatmap(cor_mat,
               clustering_distance_rows = as.dist(1 - cor_mat),
               clustering_distance_cols = as.dist(1 - cor_mat),
               clustering_method = "complete",
               color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
               breaks = seq(min(cor_mat), 1, length.out = 101),
               annotation_col = anno_col,
               annotation_colors = anno_colors,
               show_rownames = TRUE, show_colnames = TRUE,
               fontsize = 8,
               main = "Sample Correlation (Spearman)"))
    },
    filename = "Correlation_heatmap_spearman",
    subdir = "02_PCA_Clustering", width = 10, height = 9
  )
  
  # 5. Hierarchical clustering
  log_step("03_NORM", "Hierarchical clustering...")
  dist_mat <- as.dist(1 - cor_mat)
  hc <- hclust(dist_mat, method = "complete")
  dend <- as.dendrogram(hc)
  labels_order <- labels(dend)
  label_colors <- ifelse(sample_info[labels_order, "Treatment"] == "Induced",
                         COLORS$treatment["Induced"],
                         COLORS$treatment["Control"])
  dend <- dend %>%
    dendextend::set("labels_cex", 0.7) %>%
    dendextend::set("labels_col", label_colors)
  
  save_heatmap_fig(
    draw_func = function() {
      par(mar = c(8, 4, 3, 1))
      plot(dend, main = "Sample Hierarchical Clustering",
           ylab = "1 - Spearman Correlation", xlab = "")
    },
    filename = "Sample_clustering_dendrogram",
    subdir = "02_PCA_Clustering", width = 12, height = 6
  )
  
  # 6. Distance heatmap
  sample_dist <- dist(t(vst_mat))
  sample_dist_mat <- as.matrix(sample_dist)
  save_heatmap_fig(
    draw_func = function() {
      print(pheatmap::pheatmap(sample_dist_mat,
               clustering_distance_rows = sample_dist,
               clustering_distance_cols = sample_dist,
               color = colorRampPalette(c("#00A087", "white", "#E64B35"))(100),
               annotation_col = anno_col,
               annotation_colors = anno_colors,
               show_rownames = TRUE, show_colnames = FALSE,
               fontsize = 7,
               main = "Sample Distance (Euclidean on VST)"))
    },
    filename = "Sample_distance_heatmap",
    subdir = "02_PCA_Clustering", width = 10, height = 9
  )
  
  log_step("03_NORM", "Step 03 COMPLETE (ADAPTED)")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  ADAPTED SCRIPT: 06_maSigPro_trends                                       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_06_maSigPro_adapted <- function() {
  ensure_env()
  
  log_step("06_MASIGPRO", sprintf(
    "Starting maSigPro (ADAPTED: degree=%d, k=%d)...",
    PARAMS$masigpro_degree, PARAMS$masigpro_k))
  
  counts_filtered <- readRDS(FILES$counts_filtered)
  sample_info     <- readRDS(FILES$sample_info)
  sample_info     <- fix_sample_info_factors(sample_info)
  gene_anno       <- readRDS(FILES$gene_annotation)
  
  # 1. edesign
  edesign <- data.frame(
    Time      = sample_info$time_num,
    Replicate = sample_info$replicate,
    Control   = as.integer(sample_info$Treatment == "Control"),
    Induced   = as.integer(sample_info$Treatment == "Induced"),
    row.names = rownames(sample_info)
  )
  edesign <- edesign[colnames(counts_filtered), ]
  
  log_step("06_MASIGPRO", sprintf("Design: %d samples, Time: %s",
    nrow(edesign), paste(sort(unique(edesign$Time)), collapse = ", ")))
  
  # 2. maSigPro analysis — degree ADAPTED
  data_mat <- as.matrix(counts_filtered)
  design_matrix <- make.design.matrix(edesign, degree = PARAMS$masigpro_degree)
  
  min_obs_use <- PROFILE$masigpro_min_obs
  log_step("06_MASIGPRO", sprintf("Step 1: p.vector (min.obs=%d)...", min_obs_use))
  fit <- p.vector(data_mat, design_matrix,
                  Q = PARAMS$masigpro_alfa,
                  counts = TRUE,
                  min.obs = min_obs_use)
  log_step("06_MASIGPRO", sprintf("Step 1: %d significant genes", fit$i))
  
  log_step("06_MASIGPRO", "Step 2: T.fit (stepwise regression)...")
  tstep <- T.fit(fit, step.method = "backward", alfa = PARAMS$masigpro_alfa)
  
  log_step("06_MASIGPRO", sprintf("Step 3: Filtering R² >= %.2f...", PARAMS$masigpro_rsq))
  sigs <- get.siggenes(tstep, rsq = PARAMS$masigpro_rsq, vars = "groups")
  sig_genes_all <- sigs$summary
  
  for (grp in names(sig_genes_all)) {
    log_step("06_MASIGPRO", sprintf("  %s: %d genes", grp, length(sig_genes_all[[grp]])))
  }
  
  # 3. Clustering
  cluster_result <- NULL
  tryCatch({
    sig_genes_for_cluster <- NULL
    if (!is.null(sigs$sig.genes)) {
      avail_keys <- names(sigs$sig.genes)
      log_step("06_MASIGPRO", sprintf("sig.genes keys: %s", paste(avail_keys, collapse = ", ")))
      if ("InducedvsControl" %in% avail_keys) {
        sig_genes_for_cluster <- sigs$sig.genes$InducedvsControl
      } else if (length(avail_keys) > 0) {
        sig_genes_for_cluster <- sigs$sig.genes[[avail_keys[1]]]
      }
    }
    
    if (!is.null(sig_genes_for_cluster)) {
      # Dynamically adjust k: reduce k if too few genes
      n_sig <- if (is.list(sig_genes_for_cluster)) 
        nrow(sig_genes_for_cluster$sig.profiles) else nrow(sig_genes_for_cluster)
      k_use <- min(PARAMS$masigpro_k, max(2, floor(n_sig / 5)))
      log_step("06_MASIGPRO", sprintf("Clustering %d genes into %d groups...", n_sig, k_use))
      
      pdf(file.path(FIG_DIR, "06_maSigPro", "maSigPro_cluster_profiles.pdf"),
          width = 14, height = 10)
      cluster_result <- see.genes(sig_genes_for_cluster,
                                  edesign = edesign,
                                  groups.vector = design_matrix$groups.vector,
                                  show.fit = TRUE,
                                  dis = design_matrix$dis,
                                  cluster.method = "hclust",
                                  cluster.data = 1,
                                  k = k_use,
                                  newX11 = FALSE)
      dev.off()
      log_step("06_MASIGPRO", "Cluster profiles plotted")
    } else {
      log_step("06_MASIGPRO", "WARNING: No sig genes for clustering")
    }
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    log_step("06_MASIGPRO", sprintf("WARNING: Clustering error: %s", e$message))
  })
  
  # 4. Extract clusters
  gene_clusters <- NULL
  if (!is.null(cluster_result)) {
    gene_clusters <- data.frame(
      ensembl_id = names(cluster_result$cut),
      cluster    = cluster_result$cut,
      stringsAsFactors = FALSE
    )
    gene_clusters$symbol <- gene_anno$hgnc_symbol[
      match(gene_clusters$ensembl_id, gene_anno$ensembl_gene_id)]
    log_step("06_MASIGPRO", "Genes per cluster:")
    print(table(gene_clusters$cluster))
    write.csv(gene_clusters, file.path(DATA_DIR, "maSigPro_gene_clusters.csv"),
              row.names = FALSE)
  }
  
  # 5. Custom trend plots (ADAPTED: use actual time values for axis)
  if (!is.null(cluster_result)) {
    tpm_filtered <- readRDS(FILES$tpm_filtered)
    time_vals <- sort(unique(sample_info$time_num))
    time_labs <- PROFILE$time_labels
    
    for (cl in sort(unique(gene_clusters$cluster))) {
      cl_genes <- gene_clusters$ensembl_id[gene_clusters$cluster == cl]
      if (length(cl_genes) < 2) next
      
      tpm_cl_log <- log2(tpm_filtered[cl_genes, , drop = FALSE] + 1)
      
      plot_data <- tpm_cl_log %>%
        as.data.frame() %>%
        tibble::rownames_to_column("gene") %>%
        tidyr::pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
        dplyr::left_join(
          sample_info %>% tibble::rownames_to_column("sample") %>%
            dplyr::select(sample, Treatment, Time, time_num),
          by = "sample") %>%
        dplyr::group_by(Treatment, Time, time_num) %>%
        dplyr::summarise(
          mean_expr = mean(expr), se_expr = sd(expr) / sqrt(dplyr::n()),
          .groups = "drop")
      
      p_trend <- ggplot(plot_data, aes(x = time_num, y = mean_expr,
                                       color = Treatment, group = Treatment)) +
        geom_line(linewidth = 1) +
        geom_point(size = 2.5) +
        geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr),
                      width = max(time_vals) * 0.03, linewidth = 0.5) +
        scale_color_treatment() +
        scale_x_continuous(breaks = time_vals, labels = time_labs) +  # ★ ADAPTED ★
        labs(title = sprintf("Cluster %d (%d genes)", cl, length(cl_genes)),
             x = "Time", y = "Mean log2(TPM + 1)") +
        theme_bindlab()
      
      save_pub_fig(p_trend, sprintf("maSigPro_Cluster%d_trend", cl),
                   "06_maSigPro", width = 6, height = 4.5)
    }
  }
  
  # 6. Save
  masigpro_results <- list(
    fit            = fit,
    tstep          = tstep,
    sigs           = sigs,
    sig_genes_all  = sig_genes_all,
    cluster_result = cluster_result,
    gene_clusters  = gene_clusters,
    design_matrix  = design_matrix
  )
  save_data(masigpro_results, FILES$masigpro_results)
  
  log_step("06_MASIGPRO", "Step 06 COMPLETE (ADAPTED)")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 2: DEG + maSigPro + WGCNA                                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_phase2 <- function() {
  cat("\n### PHASE 2: DEG analysis + time series + WGCNA ###\n\n")
  ensure_env()
  
  # ---- 2a. 04_DEG_analysis.R ----
  # Issue: L65 timepoints <- c(<original_timepoints>)
  #       L137 target_groups <- c("Induced_4d", ...)
  # Strategy: after source, point timepoints and target_groups to correct values
  #       But 04 uses local variables, already executed when source completes
  #       Key issue: L65 timepoints determines which Wald test contrasts to make
  #       If Group factor is already fixed，c("group","Induced_6h","Control_6h")is a valid contrast
  #       But 04_DEG hardcodes timepoints in the for loop!
  #
  # Safest approach: temporarily override 04_DEG core logic
  # But due to the "no modification" principle, we need another strategy:
  #   → 04_DEG.R at execution, timepoints variable is defined internally
  #   → We cannot inject local variables before source
  #   → Need to write an adapted version of 04_DEG for GEO data

  cat("  [ADAPT] Running adapted DEG analysis (04_DEG)...\n")
  run_04_DEG_adapted()
  
  # ---- 2b. 05_DEG_visualization.R ----
  # Similarly, hardcoded timepoints
  cat("  [ADAPT] Running adapted DEG visualization (05)...\n")
  run_05_DEGvis_adapted()
  
  # ---- 2c. 06_maSigPro_trends.R ----
  # 06 re-sources 00_setup.R → degreeis reset to3 → 3timepoints+degree=3overflow
  # Need adapted version
  cat(sprintf("  [ADAPT] Running adapted maSigPro (degree=%d, k=%d)...\n",
      PROFILE$masigpro_degree, PROFILE$masigpro_k))
  run_06_maSigPro_adapted()
  
  # ---- 2d. 08_WGCNA.R ----
  # WGCNAReads from VST matrix and sample_info, no hardcoded timepoints
  # group_levels at L67-68 hardcodes 8 groups, but it uses actual values from the Group column
  # Wait — L67 actually is:
  #   group_levels <- c("Control_4d", ... "Induced_21d")  ← hardcoded!
  # Needs handling
  cat("  [ADAPT] Running WGCNA...\n")
  run_08_WGCNA_adapted()
  
  cat("\n### PHASE 2 COMPLETE ###\n")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  PHASE 3: Candidate filtering（09A → 09C → 09D → 09F → 10）                          ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_phase3 <- function() {
  cat("\n### PHASE 3: Candidate selection pipeline ###\n\n")
  ensure_env()
  
  # Pre-load global variables (original pipeline sources from 01 onward in RUN_GUIDE, these exist naturally)
  # Use assign to ensure placement in .GlobalEnv (<<- may target wrong scope in nested environments)
  assign("sample_info",  fix_sample_info_factors(readRDS(FILES$sample_info)), envir = .GlobalEnv)
  assign("gene_anno",    readRDS(FILES$gene_annotation), envir = .GlobalEnv)
  assign("tpm_filtered", readRDS(FILES$tpm_filtered), envir = .GlobalEnv)
  assign("all_results",  readRDS(FILES$deg_results), envir = .GlobalEnv)
  
  # 09A: candidate pool
  source(file.path(SCRIPT_DIR, "09A_candidate_pool.R"), local = FALSE)
  
  # Re-inject sample_info (ensure in .GlobalEnv)
  assign("sample_info", fix_sample_info_factors(readRDS(FILES$sample_info)), envir = .GlobalEnv)
  
  # 09C: bootstrap stability selection
  source(file.path(SCRIPT_DIR, "09C_ML_stability_selection.R"), local = FALSE)
  
  # 09D: gap-union selection
  # [v2] Do not pre-exclude any algorithm — RFautomatically excludes via MIN_FREQ_SIGNAL=0.20
  # Ref: Strobl et al. 2007 BMC Bioinf; Nicodemus et al. 2010 BMC Bioinf
  #   RFVIM is unstable in p>>n scenarios, but exclusion should be data-driven, not preset
  {
    lines_09D <- readLines(file.path(SCRIPT_DIR, "09D_gap_union_selection.R"))
    idx_exclude <- grep('EXCLUDE_ALGOS\\s*<-\\s*c\\("RF"\\)', lines_09D)
    if (length(idx_exclude) == 1) {
      lines_09D[idx_exclude] <- 'EXCLUDE_ALGOS <- character(0)  # [v2] auto-exclude by MIN_FREQ_SIGNAL'
      cat("  [v2 PATCH] 09D: RF pre-exclusion removed → auto-exclude by MIN_FREQ_SIGNAL\n")
    }
    tmp_09D <- tempfile(pattern = "09D_v2_", fileext = ".R")
    writeLines(lines_09D, tmp_09D)
    source(tmp_09D, local = FALSE)
    unlink(tmp_09D)
  }
  
  # 09F: PPI hub — needs local STRING database
  # Check if stringdb_cache exists
  string_cache <- file.path(DATA_DIR, "stringdb_cache")
  links_file <- file.path(string_cache, "9606.protein.links.v12.0.txt.gz")
  
  if (file.exists(links_file)) {
    cat("  [OK] STRING local database found\n")
    source(file.path(SCRIPT_DIR, "09F_PPI_hub_selection.R"), local = FALSE)
  } else {
    cat("  [SKIP] STRING local database not found at:\n")
    cat(sprintf("    %s\n", string_cache))
    cat("  09F_PPI_hub_selection skipped. Copy STRING files to run.\n")
    cat("  Required files:\n")
    cat("    9606.protein.links.v12.0.txt.gz\n")
    cat("    9606.protein.info.v12.0.txt.gz\n")
  }
  
  # 10_integration.R — contains hardcoded timepoints on L124
  cat("  [ADAPT] Running integration (10)...\n")
  run_10_integration_adapted()
  
  cat("\n### PHASE 3 COMPLETE ###\n")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  ADAPTED SCRIPT FUNCTIONS                                                  ║
# ║  These functions replace original pipeline scripts that contain hardcoded time references                              ║
# ║  Logic is identical, just replacing c(<original_timepoints>) with get_timepoints()       ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

# ------------------------------------------------------------------
# 04_DEG_analysis adapted
# ------------------------------------------------------------------
run_04_DEG_adapted <- function() {
  ensure_env()
  
  log_step("04_DEG", "Starting DEG analysis (ADAPTED for GEO)...")
  
  dds         <- readRDS(FILES$dds_object)
  sample_info <- readRDS(FILES$sample_info)
  gene_anno   <- readRDS(FILES$gene_annotation)
  
  # --- 1a. LRT: Treatment:Time interaction ---
  log_step("04_DEG", "LRT for Treatment:Time interaction...")
  dds_lrt_interaction <- DESeq(dds, test = "LRT",
                               reduced = ~ Treatment + Time)
  res_lrt_interaction <- results(dds_lrt_interaction, alpha = PARAMS$padj_cutoff)
  n_sig_lrt <- sum(res_lrt_interaction$padj < PARAMS$padj_cutoff, na.rm = TRUE)
  log_step("04_DEG", sprintf("LRT Interaction: %d genes with padj < %.2f",
                              n_sig_lrt, PARAMS$padj_cutoff))
  
  # --- 1b. LRT: Treatment main effect ---
  log_step("04_DEG", "LRT for Treatment effect...")
  dds_lrt_treatment <- DESeq(dds, test = "LRT",
                              reduced = ~ Time)
  res_lrt_treatment <- results(dds_lrt_treatment, alpha = PARAMS$padj_cutoff)
  n_sig_trt <- sum(res_lrt_treatment$padj < PARAMS$padj_cutoff, na.rm = TRUE)
  log_step("04_DEG", sprintf("LRT Treatment: %d genes with padj < %.2f",
                              n_sig_trt, PARAMS$padj_cutoff))
  
  # --- 2. Wald test per timepoint (ADAPTED) ---
  log_step("04_DEG", "Wald tests for each timepoint...")
  dds_group <- dds
  dds_group$group <- factor(paste(dds_group$Treatment, dds_group$Time, sep = "_"))
  design(dds_group) <- ~ group
  dds_group <- DESeq(dds_group)
  
  # ★ Key adaptation point: use actual time labels ★
  timepoints <- get_timepoints()
  # Exclude timepoints with only Control group (e.g., GSE197067 0h)
  if (!is.null(PROFILE$has_0h_control_only) && PROFILE$has_0h_control_only) {
    timepoints <- setdiff(timepoints, "0h")
    cat("  [ADAPT] Excluding 0h from Wald contrasts (Control-only)\n")
  }
  
  deg_by_time <- list()
  deg_by_time_lfc <- list()
  
  for (tp in timepoints) {
    contrast_name <- paste0("Induced_", tp, "_vs_Control_", tp)
    contrast_vec <- c("group", paste0("Induced_", tp), paste0("Control_", tp))
    
    res <- results(dds_group, contrast = contrast_vec, alpha = PARAMS$padj_cutoff)
    
    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("ensembl_id") %>%
      dplyr::left_join(gene_anno, by = c("ensembl_id" = "ensembl_gene_id")) %>%
      dplyr::arrange(padj) %>%
      dplyr::mutate(
        timepoint = tp,
        regulation = dplyr::case_when(
          padj < PARAMS$padj_cutoff & log2FoldChange > PARAMS$lfc_cutoff ~ "Up",
          padj < PARAMS$padj_cutoff & log2FoldChange < -PARAMS$lfc_cutoff ~ "Down",
          TRUE ~ "NS"
        )
      )
    deg_by_time[[tp]] <- res_df
    
    # Effect size test (lfcThreshold=1)
    res_lfc <- results(dds_group, contrast = contrast_vec,
                       alpha = PARAMS$padj_cutoff,
                       lfcThreshold = PARAMS$lfc_cutoff)
    res_lfc_df <- as.data.frame(res_lfc) %>%
      tibble::rownames_to_column("ensembl_id") %>%
      dplyr::select(ensembl_id, padj_lfc = padj)
    deg_by_time_lfc[[tp]] <- res_lfc_df
    
    n_up <- sum(res_df$regulation == "Up", na.rm = TRUE)
    n_down <- sum(res_df$regulation == "Down", na.rm = TRUE)
    log_step("04_DEG", sprintf("  %s: Up=%d, Down=%d", tp, n_up, n_down))
  }
  
  # --- 3. 8-group model: stage-specific contrasts (ADAPTED) ---
  log_step("04_DEG", "8-Group Model: Stage-Specific Planned Contrasts")
  dds_8g <- dds
  dds_8g$group <- factor(paste(dds_8g$Treatment, dds_8g$Time, sep = "_"))
  design(dds_8g) <- ~ 0 + group
  dds_8g <- DESeq(dds_8g)
  
  group_names <- levels(dds_8g$group)
  target_groups <- get_target_groups()
  # Exclude timepoints with only Control
  if (!is.null(PROFILE$has_0h_control_only) && PROFILE$has_0h_control_only) {
    target_groups <- setdiff(target_groups, "Induced_0h")
  }
  
  n_groups <- length(group_names)
  stage_contrasts <- list()
  stage_deg <- list()
  
  for (target in target_groups) {
    contrast_vec <- rep(0, n_groups)
    names(contrast_vec) <- paste0("group", group_names)
    contrast_vec[paste0("group", target)] <- 1
    other_groups <- setdiff(group_names, target)
    for (og in other_groups) {
      contrast_vec[paste0("group", og)] <- -1 / (n_groups - 1)
    }
    
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
    
    n_up <- sum(res_df$regulation == "Up", na.rm = TRUE)
    n_down <- sum(res_df$regulation == "Down", na.rm = TRUE)
    log_step("04_DEG", sprintf("  %s vs rest: Up=%d, Down=%d", target, n_up, n_down))
  }
  
  # Construct derived objects that downstream scripts expect
  lrt_interaction_df <- as.data.frame(res_lrt_interaction) %>%
    tibble::rownames_to_column("ensembl_id") %>%
    dplyr::left_join(gene_anno, by = c("ensembl_id" = "ensembl_gene_id")) %>%
    dplyr::arrange(padj)
  
  lrt_interaction_sig <- lrt_interaction_df %>%
    dplyr::filter(padj < PARAMS$lrt_padj)
  
  deg_all_timepoints <- dplyr::bind_rows(deg_by_time)
  
  log_step("04_DEG", sprintf("LRT interaction sig genes: %d", nrow(lrt_interaction_sig)))
  
  # Save — field names must match what 09A/09C/10 expect
  all_results <- list(
    lrt_interaction     = lrt_interaction_df,
    lrt_interaction_sig = lrt_interaction_sig,   # ★ 09C expects this ★
    lrt_treatment       = as.data.frame(res_lrt_treatment) %>%
                            tibble::rownames_to_column("ensembl_id"),
    wald_by_time        = deg_by_time,           # ★ original field name ★
    wald_by_time_lfc    = deg_by_time_lfc,       # ★ 09A expects this ★
    wald_all            = deg_all_timepoints,
    stage_deg           = stage_deg,
    stage_contrasts     = stage_contrasts,
    timepoints          = timepoints,
    dds_group           = dds_group
  )
  
  save_data(all_results, FILES$deg_results)
  
  # Save LRT separately (some scripts load it independently)
  save_data(res_lrt_interaction, FILES$lrt_results)
  
  # Save DEG summary table
  deg_all_timepoints <- dplyr::bind_rows(deg_by_time)
  write.csv(deg_all_timepoints,
            file.path(DATA_DIR, "DEG_all_timepoints.csv"),
            row.names = FALSE)
  
  log_step("04_DEG", "Step 04 COMPLETE (ADAPTED)")
}


# ------------------------------------------------------------------
# 05_DEG_visualization adapted (simplified)
# ------------------------------------------------------------------
run_05_DEGvis_adapted <- function() {
  ensure_env()
  
  log_step("05_VIS", "Starting DEG visualization (ADAPTED)...")
  
  all_results <- readRDS(FILES$deg_results)
  sample_info <- readRDS(FILES$sample_info)
  gene_anno   <- readRDS(FILES$gene_annotation)
  
  timepoints <- all_results$timepoints  # Use actual timepoints saved by step 04
  
  # DEG count bar chart
  deg_summary <- data.frame()
  for (tp in timepoints) {
    res <- all_results$deg_by_time[[tp]]
    if (is.null(res)) next
    n_up   <- sum(res$regulation == "Up", na.rm = TRUE)
    n_down <- sum(res$regulation == "Down", na.rm = TRUE)
    deg_summary <- rbind(deg_summary, data.frame(
      Timepoint = tp,
      Direction = c("Up", "Down"),
      Count = c(n_up, -n_down),
      stringsAsFactors = FALSE
    ))
  }
  
  if (nrow(deg_summary) > 0) {
    deg_summary$Timepoint <- factor(deg_summary$Timepoint, levels = timepoints)
    
    p_bar <- ggplot(deg_summary, aes(x = Timepoint, y = Count, fill = Direction)) +
      geom_bar(stat = "identity", width = 0.6, position = "identity") +
      scale_fill_manual(values = c("Up" = "#E64B35", "Down" = "#3C5488")) +
      geom_hline(yintercept = 0, linewidth = 0.5) +
      labs(title = "DEG Counts by Timepoint",
           subtitle = sprintf("|log2FC| > %.1f, padj < %.2f",
                              PARAMS$lfc_cutoff, PARAMS$padj_cutoff),
           x = "Timepoint", y = "Number of DEGs (Down / Up)") +
      theme_bindlab()
    
    save_pub_fig(p_bar, "DEG_barplot", "05_DEG_vis")
  }
  
  # Volcano plots per timepoint
  for (tp in timepoints) {
    res <- all_results$deg_by_time[[tp]]
    if (is.null(res)) next
    
    p_vol <- ggplot(res, aes(x = log2FoldChange, y = -log10(padj), color = regulation)) +
      geom_point(size = 0.8, alpha = 0.6) +
      scale_color_manual(values = c("Up" = "#E64B35", "Down" = "#3C5488", "NS" = "grey70")) +
      geom_vline(xintercept = c(-PARAMS$lfc_cutoff, PARAMS$lfc_cutoff),
                 linetype = "dashed", color = "grey40") +
      geom_hline(yintercept = -log10(PARAMS$padj_cutoff),
                 linetype = "dashed", color = "grey40") +
      labs(title = sprintf("Volcano Plot — %s", tp),
           x = "log2 Fold Change", y = "-log10(adjusted p-value)") +
      theme_bindlab()
    
    save_pub_fig(p_vol, sprintf("Volcano_%s", gsub("[^a-zA-Z0-9]", "", tp)),
                 "05_DEG_vis", width = 7, height = 6)
  }
  
  log_step("05_VIS", "Step 05 COMPLETE (ADAPTED)")
}


# ------------------------------------------------------------------
# 08_WGCNA adapted
# ------------------------------------------------------------------
run_08_WGCNA_adapted <- function() {
  ensure_env()
  
  log_step("08_WGCNA", "Starting WGCNA (ADAPTED)...")
  
  # Inject correct sample_info into global environment
  si <- fix_sample_info_factors(readRDS(FILES$sample_info))
  assign("sample_info", si, envir = .GlobalEnv)
  
  # Source original 08_WGCNA.R directly (with fixes: power R2>=0.80 + dynamic group_levels)
  source(file.path(SCRIPT_DIR, "08_WGCNA.R"), local = FALSE)
  
  log_step("08_WGCNA", "Step 08 COMPLETE (ADAPTED)")
}


# ------------------------------------------------------------------
# 10_integration adapted
# ------------------------------------------------------------------
run_10_integration_adapted <- function() {
  # 10_integration.R L124: for (tp in c(<original_timepoints>))
  # This section annotates with deg_by_time results; returns NULL if key doesn't exist
  # deg_results now stores adapted timepoints
  # But 10_integration hardcodes the for loop → annotation lost
  # Impact: not fatal, just per-timepoint LFC annotation columns are empty in integration table
  
  # Strategy: source directly, accept that per-timepoint annotation may be lost
  # Core filtering logic (ML + PPI union) does not depend on these annotations
  
  cat("  [INFO] Running original 10_integration.R...\n")
  cat("  [INFO] Per-timepoint LFC annotation may be incomplete (hardcoded timepoints from original analysis)\n")
  source(file.path(SCRIPT_DIR, "10_integration.R"), local = FALSE)
  
  log_step("10_INT", "Step 10 COMPLETE (via wrapper)")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Main execution entry point                                                                ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

run_all <- function() {
  cat("\n================================================================\n")
  cat(sprintf("  METI-FS Full Pipeline Run: %s\n", DATASET_ID))
  cat(sprintf("  Started: %s\n", format(Sys.time())))
  cat("================================================================\n")
  
  t0 <- Sys.time()
  
  run_phase1()
  run_phase2()
  run_phase3()
  
  elapsed <- difftime(Sys.time(), t0, units = "mins")
  
  cat("\n================================================================\n")
  cat(sprintf("  COMPLETE: %s\n", DATASET_ID))
  cat(sprintf("  Total time: %.1f minutes\n", as.numeric(elapsed)))
  cat(sprintf("  Output: %s\n", PROFILE$project_dir))
  cat("================================================================\n")
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  STRING database preparation (utility functions)                                               ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

#' Copy STRING local database to current GEO project
#' or create symbolic link
setup_string_cache <- function(source_project = file.path(METHODS_BASE, "..", "reference_project")) {
  src_cache <- file.path(source_project, "data", "stringdb_cache")
  dst_cache <- file.path(PROFILE$project_dir, "data", "stringdb_cache")
  
  if (!dir.exists(src_cache)) {
    cat("  [ERROR] Source STRING cache not found: ", src_cache, "\n")
    return(invisible(FALSE))
  }
  
  if (dir.exists(dst_cache)) {
    cat("  [OK] Destination STRING cache already exists\n")
    return(invisible(TRUE))
  }
  
  dir.create(dst_cache, recursive = TRUE, showWarnings = FALSE)
  
  # Copy files
  files_to_copy <- list.files(src_cache, full.names = TRUE)
  for (f in files_to_copy) {
    file.copy(f, file.path(dst_cache, basename(f)), overwrite = FALSE)
    cat(sprintf("  [COPIED] %s\n", basename(f)))
  }
  
  cat(sprintf("  [DONE] STRING cache prepared: %s\n", dst_cache))
  invisible(TRUE)
}


# ╔══════════════════════════════════════════════════════════════════════════════╗
# ║  Usage instructions (shown when running this script directly)                                            ║
# ╚══════════════════════════════════════════════════════════════════════════════╝

if (sys.nframe() == 0) {
  cat("\n")
  cat("================================================================\n")
  cat("  S04_run_pipeline_wrapper.R — Usage\n")
  cat("================================================================\n")
  cat("\n")
  cat("  # Run per phase (recommended, allows inspecting intermediate results):\n")
  cat("  DATASET_ID <- 'GSE307424'    # Select dataset\n")
  cat("  source('S04_run_pipeline_wrapper.R')  # Load configuration\n")
  cat("  run_phase1()                 # Data import + preprocessing + normalization\n")
  cat("  run_phase2()                 # DEG+maSigPro+WGCNA\n")
  cat("  run_phase3()                 # Candidate selection + integration\n")
  cat("\n")
  cat("  # Run all at once:\n")
  cat("  run_all()\n")
  cat("\n")
  cat("  # STRING database preparation (needed for 09F):\n")
  cat("  setup_string_cache()         # Copy from source project\n")
  cat("\n")
  cat("  # View current dataset profile:\n")
  cat("  str(PROFILE)\n")
  cat("\n")
  cat("================================================================\n")
}
