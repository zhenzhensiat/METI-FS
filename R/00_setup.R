#!/usr/bin/env Rscript
# ==============================================================================
# 00_setup.R — Environment setup, package installation, and global parameters
# Project: MSCtri-lineage differentiation group Pipeline
# Updated: 2026-03-20 (Parameter fixes final version)
#
# Change log:
#   2026-03-20: R²from 0.6 to 0.7(Nueda2014), STRINGfrom 400 to 700(high confidence),
#               WGCNA sig_modulefrom OR to AND(p<0.05&|cor|>0.5),
#               greyModule exclusion fix(ME0+MEgrey), treatment_labelsComplete,
#               fig_subdirsComplete07A/07B/07C/05_DEG_vis, 
#               save_pub_figAdd lineage_label initialization guard
# ==============================================================================

# ---- 0. Project path configuration ----
# Usage: set PROJECT_DIR in console first，then source this file
#   PROJECT_DIR <- "Lineage_A"   # Adipogenic
#   PROJECT_DIR <- "Lineage_B"      # Osteogenic
#   PROJECT_DIR <- "Lineage_C"      # Chondrogenic
#   source("R//00_setup.R")

# ---- 0a. Fix Windows temp directory issue with non-ASCII usernames ----
# WindowsDefault TMPDIR is under user AppData，If username contains non-ASCII characters，
# gseGO/gseKEGG/fgseapackages will fail when writing temp files
# "Cannot open file...No such file or directory"
if (.Platform$OS.type == "windows") {
  safe_tmp <- "C:/Temp"
  if (!dir.exists(safe_tmp)) dir.create(safe_tmp, recursive = TRUE)
  Sys.setenv(TMPDIR = safe_tmp)
  Sys.setenv(TMP = safe_tmp)
  Sys.setenv(TEMP = safe_tmp)
}

# If PROJECT_DIR not defined in console，Use default value
if (!exists("PROJECT_DIR")) {
  PROJECT_DIR <- "Lineage_A"
  message("PROJECT_DIR not set, using default: ", PROJECT_DIR)
}

# ScriptsDirectory independent of three lineage projects，under projects/
SCRIPT_DIR   <- "R/"
RAW_DIR      <- file.path(PROJECT_DIR, "data_raw")
DATA_DIR     <- file.path(PROJECT_DIR, "data")
FIG_DIR      <- file.path(PROJECT_DIR, "Figure")

# Create subdirectory
for (d in c(DATA_DIR, FIG_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Figurecreate subdirectory under，Categorized by analysis module
fig_subdirs <- c(
  "01_QC", "02_PCA_Clustering", "03_DEG", "04_Venn",
  "04B_StageSpecific", "04B_Kinetic",
  "05_DEG_vis", "05_Heatmap",
  "06_maSigPro",
  "07_Enrichment", "07A_Enrichment", "07B_Enrichment", "07C_Enrichment",
  "08_WGCNA",
  "09_PPI", "09B_ML", "09C_ML_Stability", "09D_GapUnion", "09E_PPI",
  "09F_PPI_Hub", "09G_RankAgg",
  "10_Integration",
  "PushPull", "CohenD"
)
for (sd in fig_subdirs) {
  dir.create(file.path(FIG_DIR, sd), showWarnings = FALSE, recursive = TRUE)
}

# ---- 1. and Load ----
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (bioc) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install(pkg, update = FALSE, ask = FALSE)
    } else {
      install.packages(pkg, repos = "https://cloud.r-project.org")
    }
  }
}

# CRAN packages
cran_pkgs <- c("tidyverse", "ggplot2", "pheatmap", "RColorBrewer",
               "ggrepel", "VennDiagram", "gridExtra", "cowplot",
               "scales", "circlize", "reshape2", "WGCNA",
               "flashClust", "dynamicTreeCut", "igraph",
               "ggpubr", "pROC", "corrplot", "dendextend",
               "glmnet", "randomForest", "caret",
               "futile.logger", "ragg")
for (p in cran_pkgs) install_if_missing(p)

# Bioconductor packages
# impute and preprocessCore is WGCNA Dependencies，
bioc_pkgs <- c("impute", "preprocessCore", "GO.db",
               "DESeq2", "edgeR", "clusterProfiler", "org.Hs.eg.db",
               "AnnotationDbi", "enrichplot", "pathview",
               "maSigPro", "ComplexHeatmap", "STRINGdb",
               "DOSE", "ReactomePA", "biomaRt")
for (p in bioc_pkgs) install_if_missing(p, bioc = TRUE)

# ---- 2. Load ----
suppressPackageStartupMessages({
  library(tidyverse)
  library(DESeq2)
  library(pheatmap)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(ggrepel)
  library(VennDiagram)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(DOSE)
  library(maSigPro)
  library(WGCNA)
  library(igraph)
  library(STRINGdb)
  library(cowplot)
  library(ggpubr)
  library(pROC)
 # glmnet, randomForest, caret — in 09B in need Load
 # caret::margin mask ggplot2::margin
  library(corrplot)
  library(scales)
  library(dendextend)
})

# WGCNA
allowWGCNAThreads()
options(stringsAsFactors = FALSE)

# ---- 3. ----
PARAMS <- list(
 # ── Project identifier (from PROJECT_DIR) ──
  diff_type = basename(PROJECT_DIR),

  # ── Timepoint configuration ──
  time_labels = strsplit("4d,7d,14d,21d", ",")[[1]], # Default; override per dataset

  # ── Gene filtering ──
  min_count      = 3,       # Minimum counts threshold
  min_samples    = 3,       # At least N samples above threshold

  # ── DifferentialExpression ──
  padj_cutoff    = 0.05,    # FDRThreshold
  lfc_cutoff     = 1,       # log2FCThreshold（|log2FC| > 1）
  lrt_padj       = 0.05,    # LRT interaction FDR threshold

  # ── WGCNA ──
  wgcna_network    = "unsigned",   # Main result (signed as reference)
  wgcna_corFnc     = "bicor",
  wgcna_minModSize = 30,
  wgcna_mergeCut   = 0.25,        # Module merge similarity threshold
  wgcna_sig_cor    = 0.5,         # Key modules: |cor| > This value (ANDLogic)
  wgcna_sig_p      = 0.05,        # Key modules: p < This value (ANDLogic)
 # ↑ ANDLogic (Langfelder & Horvath), grey/ME0

  # ── maSigPro ──
  masigpro_alfa   = 0.05,   # Step 1 significance
  masigpro_rsq    = 0.7,    # R-squared threshold (Nueda et al. 2014 default)
  masigpro_degree = 3,      # Polynomial degree (4 timepoints, degree=nTP-1=3)
  masigpro_k      = 9,      # Number of clusters

 # ── Enrichment ──
  enrich_pvalue  = 0.05,
  enrich_qvalue  = 0.05,

  # ── PPI ──
  string_version = "12.0",
  string_species = 9606,    # Human
  string_score   = 700,     # High confidence (STRING official definition >= 700)
 # ↑ from 400(medium) is 700(high confidence)

  # ── Machine learningFeature selection ──
  ml_nfolds       = -1,     # -1 = LOOCV
  ml_lasso_family = "binomial",
  ml_rf_ntree     = 1000,
  ml_svm_sizes    = c(5, 10, 20, 30, 50, 75, 100, 150, 200),

 # ── Specificity (B) ──
  tau_strict      = 0.85,
  tau_moderate    = 0.75,
  tpm_on          = 3,
  tpm_off         = 1,
  peak_ratio_min  = 3,
  consistency_min = 0.75,

 # ── Kinetic classification (A) ──
  kinetic_lfc_strong    = 1.0,
  kinetic_lfc_weak      = 0.5,
  kinetic_padj          = 0.05,
  kinetic_decay_delta   = 1.0,
  kinetic_peak_delta    = 0.75,
  kinetic_sustained_min_tp = 3,

 # ── Figure ──
  fig_width   = 8,
  fig_height  = 6,
  fig_dpi     = 300,
  fig_format  = "pdf"
)

# ---- 4. Define ----
.prefix <- PARAMS$diff_type

FILES <- list(
  # Raw data
  counts_raw    = file.path(RAW_DIR, paste0(.prefix, "_all_counts_with_order.tsv")),
  tpm_raw       = file.path(RAW_DIR, paste0(.prefix, "_all_tpm.tsv")),
  gene_anno_raw = file.path(RAW_DIR, paste0(.prefix, "_gene_annotation.tsv")),
  metadata_raw  = file.path(RAW_DIR, paste0(.prefix, "_metadata.tsv")),

 # Intermediate
  counts_filtered   = file.path(DATA_DIR, "counts_filtered.rds"),
  tpm_filtered      = file.path(DATA_DIR, "tpm_filtered.rds"),
  filter_log        = file.path(DATA_DIR, "gene_filter_log.csv"),
  sample_info       = file.path(DATA_DIR, "sample_info.rds"),
  dds_object        = file.path(DATA_DIR, "dds_object.rds"),
  vst_matrix        = file.path(DATA_DIR, "vst_matrix.rds"),
  deg_results       = file.path(DATA_DIR, "deg_results.rds"),
  group_deg_results = file.path(DATA_DIR, "group_deg_results.rds"),
  stage_specificity = file.path(DATA_DIR, "stage_specificity_results.rds"),
  kinetic_results   = file.path(DATA_DIR, "kinetic_classification.rds"),
  lrt_results       = file.path(DATA_DIR, "lrt_interaction_results.rds"),
  masigpro_results  = file.path(DATA_DIR, "masigpro_results.rds"),
  wgcna_results     = file.path(DATA_DIR, "wgcna_results.rds"),
  ml_results        = file.path(DATA_DIR, "ml_feature_selection.rds"),
  hub_genes         = file.path(DATA_DIR, "hub_genes_final.csv"),
  gene_annotation   = file.path(DATA_DIR, "gene_annotation.rds")
)

# ---- 5. Color schemeDefine（Nature） ----
# （ using Figurelegend and caption）
.lineage_label <- switch(PARAMS$diff_type,
  "Lineage_A" = "Adipogenic",
  "Lineage_B"   = "Osteogenic",
  "Lineage_C"   = "Chondrogenic",
  PARAMS$diff_type
)

# （ using treatmentMapping）
.lineage_short <- switch(PARAMS$diff_type,
  "Lineage_A" = "Adi",
  "Lineage_B"   = "Ost",
  "Lineage_C"   = "Cho",
  "Diff"
)

COLORS <- list(
  treatment = c("Induced" = "#E64B35", "Control" = "#4DBBD5"),

 # treatmentMapping（Figure in Induced/Control → Adi/Adi-Control etc.）
  treatment_labels = setNames(
    c(.lineage_short, paste0(.lineage_short, "-Control")),
    c("Induced", "Control")
  ),

  # Timepoint
  timepoint = c("4d" = "#3C5488", "7d" = "#00A087",
                "14d" = "#F39B7F", "21d" = "#E64B35"),

 # group （treatment x time）
  group = c("Induced_4d"  = "#DC0000", "Control_4d"  = "#3C5488",
            "Induced_7d"  = "#F39B7F", "Control_7d"  = "#4DBBD5",
            "Induced_14d" = "#E64B35", "Control_14d" = "#00A087",
            "Induced_21d" = "#B09C85", "Control_21d" = "#7E6148"),

 # above Down-regulated
  regulation = c("Up" = "#E64B35", "Down" = "#3C5488", "NS" = "grey70"),

  # Kinetic classification
  kinetic = c("Early-on_up"    = "#E64B35",
              "Late_up"        = "#3C5488",
              "Sustained_up"   = "#00A087",
              "Transient_up"   = "#F39B7F",
              "Early-on_down"  = "#DC9FB4",
              "Late_down"      = "#91D1C2",
              "Sustained_down" = "#7E6148",
              "None"           = "grey80"),

 # Figure
  heatmap_col = colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)
)

# ---- 6. ----

#' SaveFigure：PDF（）+ PNG（）
#' in Figure right below Addcaption
save_pub_fig <- function(plot_obj, filename, subdir,
                         width = PARAMS$fig_width,
                         height = PARAMS$fig_height,
                         dpi = PARAMS$fig_dpi) {

  lineage_label <- .lineage_label  # Pre-defined to prevent uninitialized variable outside ggplot branch

 # Addcaption（ggplot）
  if (inherits(plot_obj, "gg")) {
    plot_obj <- plot_obj +
      labs(caption = lineage_label) +
      theme(plot.caption = element_text(face = "italic", color = "grey50",
                                         size = 8, hjust = 1))
  }

  # PDF
  pdf_path <- file.path(FIG_DIR, subdir, paste0(filename, ".pdf"))
  if (.Platform$OS.type == "windows") {
    ggsave(pdf_path, plot = plot_obj, width = width, height = height,
           device = cairo_pdf)
  } else {
    ggsave(pdf_path, plot = plot_obj, width = width, height = height,
           device = "pdf", useDingbats = FALSE)
  }

  # PNG preview
 # : ggplot2 >= 3.5.0 already using ggsave(..., type=) 
 # type="cairo" in will FailedFigure not 
 # using ragg::agg_png ( already ) or Default png 
  png_path <- file.path(FIG_DIR, subdir, paste0(filename, ".png"))
  if (requireNamespace("ragg", quietly = TRUE)) {
    ggsave(png_path, plot = plot_obj, width = width, height = height,
           dpi = dpi, device = ragg::agg_png)
  } else {
    ggsave(png_path, plot = plot_obj, width = width, height = height,
           dpi = dpi, device = "png")
  }
  message(sprintf("[SAVED] %s (.pdf + .png) [%s]", filename, lineage_label))
}

#' SaveComplexHeatmap or pheatmap（ggplot）
save_heatmap_fig <- function(draw_func, filename, subdir,
                             width = PARAMS$fig_width,
                             height = PARAMS$fig_height) {
  pdf_path <- file.path(FIG_DIR, subdir, paste0(filename, ".pdf"))
  png_path <- file.path(FIG_DIR, subdir, paste0(filename, ".png"))

 # PDF
  pdf_ok <- FALSE
  tryCatch({
    grDevices::pdf(pdf_path, width = width, height = height)
    draw_func()
    dev.off()
    if (file.exists(pdf_path) && file.size(pdf_path) > 1000) {
      pdf_ok <- TRUE
    }
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message(sprintf("[WARNING] PDF failed: %s", e$message))
  })

  if (!pdf_ok && capabilities("cairo")) {
    tryCatch({
      grDevices::cairo_pdf(pdf_path, width = width, height = height)
      draw_func()
      dev.off()
    }, error = function(e) {
      try(dev.off(), silent = TRUE)
      message(sprintf("[WARNING] cairo_pdf also failed: %s", e$message))
    })
  }

 # PNG（ragg，）
  tryCatch({
    if (requireNamespace("ragg", quietly = TRUE)) {
      ragg::agg_png(png_path, width = width, height = height,
                    units = "in", res = PARAMS$fig_dpi)
    } else {
      grDevices::png(png_path, width = width, height = height,
                     units = "in", res = PARAMS$fig_dpi)
    }
    draw_func()
    dev.off()
  }, error = function(e) {
    try(dev.off(), silent = TRUE)
    message(sprintf("[WARNING] PNG failed for %s: %s", filename, e$message))
  })

  message(sprintf("[SAVED] %s (.pdf + .png) [%s]", filename, PARAMS$diff_type))
}

#' Ensembl ID -> Gene Symbol Transformation
ensembl_to_symbol <- function(ensembl_ids, anno_df) {
  matched <- anno_df$hgnc_symbol[match(ensembl_ids, anno_df$ensembl_gene_id)]
  matched[is.na(matched) | matched == ""] <- ensembl_ids[is.na(matched) | matched == ""]
  return(matched)
}

#' SaveIntermediate
save_data <- function(obj, filepath) {
  saveRDS(obj, filepath)
  message(sprintf("[DATA SAVED] %s (%.1f MB)",
                  basename(filepath),
                  file.size(filepath) / 1e6))
}

#' above DependenciesCheck + 
#'
#' Call:
#' 1. has above Dependencies in 
#' 2. above Update，Delete（）
#' 3. "0809A"
check_upstream <- function(step_name, upstream_files, output_files = character(0),
                           stop_on_missing = TRUE) {

 # 1. Check above in 
  for (nm in names(upstream_files)) {
    fpath <- upstream_files[[nm]]
    if (!file.exists(fpath)) {
      msg <- sprintf("Upstream file missing: %s (%s)", nm, basename(fpath))
      if (stop_on_missing) {
        log_step(step_name, paste("ERROR:", msg))
        stop(msg)
      } else {
        log_step(step_name, paste("WARNING:", msg))
      }
    }
  }

 # 2. Check： above is Update
  upstream_exist <- upstream_files[file.exists(upstream_files)]
  output_exist   <- output_files[file.exists(output_files)]

  if (length(upstream_exist) > 0 && length(output_exist) > 0) {
    latest_upstream <- max(file.mtime(upstream_exist))
    oldest_output   <- min(file.mtime(output_exist))

    if (latest_upstream > oldest_output) {
      stale_files <- output_exist[file.mtime(output_exist) < latest_upstream]
      upstream_newer <- names(upstream_exist)[file.mtime(upstream_exist) == latest_upstream][1]

      log_step(step_name, sprintf(
        "STALE OUTPUT DETECTED: %s updated at %s, but outputs are older",
        upstream_newer, format(latest_upstream, "%H:%M:%S")))

      for (sf in stale_files) {
        file.remove(sf)
        log_step(step_name, sprintf("  Deleted stale output: %s", basename(sf)))
      }

      log_step(step_name, sprintf("  Cleared %d stale outputs", length(stale_files)))
    } else {
      log_step(step_name, "Upstream check passed: all outputs up-to-date")
    }
  }

  invisible(TRUE)
}

#' 
log_step <- function(step_name, message_text) {
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  log_line <- sprintf("[%s] %s: %s", timestamp, step_name, message_text)
  cat(log_line, "\n")
  log_file <- file.path(DATA_DIR, "analysis_log.txt")
  cat(log_line, "\n", file = log_file, append = TRUE)
}

cat("========================================\n")
cat("  METI-FS Pipeline — Setup Complete\n")
cat("  Project: ", PROJECT_DIR, "\n")
cat("  Lineage: ", .lineage_label, "\n")
cat("  Timestamp:", format(Sys.time()), "\n")
cat("========================================\n")
cat("  Key params (2026-03-20 final):\n")
cat("    masigpro_rsq  =", PARAMS$masigpro_rsq, "(Nueda 2014)\n")
cat("    string_score  =", PARAMS$string_score, "(high confidence)\n")
cat("    wgcna_sig     = AND: |cor|>", PARAMS$wgcna_sig_cor,
    "& p<", PARAMS$wgcna_sig_p, "\n")
cat("========================================\n")
