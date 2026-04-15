#!/usr/bin/env Rscript
# ==============================================================================
# S02_pipeline_adapter.R — DataFormatAdapter
#
# ：RNA-seqDataConvertMETI-FS pipelineInputFormat
#
# pipeline (00_setup.R → 01_data_import.R) Input：
# data_raw/Directory：
# {Prefix}_all_counts_with_order.tsv — raw counts (genes × samples)
# {Prefix}_all_tpm.tsv — TPM
# {Prefix}_metadata.tsv — SampleData (Optional)
# {Prefix}_gene_annotation.tsv — GeneAnnotation (Optional)
# Sample：
# {Prefix}{Time}d{Rep} (Induced) / {Prefix}{Time}dC{Rep} (Control)
#
# ProcessData：
# Source A: Data (S01_simulation_engine.R Output)
# Source B: GEOData (NCBICalculatecounts + metadata)
# Source C: counts + metadata
#
# Usage：
#   source("S02_pipeline_adapter.R")
#   adapt_simulation(sim_dir, output_dir, prefix = "SimA")
#   adapt_geo_dataset(geo_accession, output_dir, prefix = "GeoA", ...)
#   adapt_custom(counts_file, metadata_file, output_dir, prefix = "CustA")
# ==============================================================================

# ---- 0. LoadConfigure ----
if (file.exists("S_config.R")) {
  source("S_config.R")
} else if (file.exists(file.path(file.path(METHODS_BASE, "scripts"), "S_config.R"))) {
  source(file.path(file.path(METHODS_BASE, "scripts"), "S_config.R"))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# ==============================================================================
# Source A: Data → pipelineFormat
# ==============================================================================

adapt_simulation <- function(sim_dir, output_dir, prefix = "Sim",
                              verbose = TRUE) {
  
  if (verbose) cat(sprintf("[ADAPT] Simulation data: %s -> %s\n", sim_dir, output_dir))
  
 # LoadData
  counts <- readRDS(file.path(sim_dir, "counts_matrix.rds"))
  tpm <- readRDS(file.path(sim_dir, "tpm_matrix.rds"))
  sample_info <- readRDS(file.path(sim_dir, "sample_info.rds"))
  
 # Samplepipeline
 # Raw: Ind_T1_R1, Ctrl_T1_R1
 # : Sim4d1, Sim4dC1 (4dT1)
  
 # ParametersTimeMap
  params <- readRDS(file.path(sim_dir, "simulation_params.rds"))
  time_map <- setNames(params$time_values, paste0("T", seq_along(params$time_values)))
  
  new_names <- sapply(seq_len(nrow(sample_info)), function(i) {
    trt <- sample_info$Treatment[i]
    tp <- sample_info$Time[i]
    rep <- sample_info$replicate[i]
    time_val <- time_map[tp]
    
    if (trt == "Induced") {
      sprintf("%s%dd%d", prefix, time_val, rep)
    } else {
      sprintf("%s%ddC%d", prefix, time_val, rep)
    }
  })
  
  colnames(counts) <- new_names
  colnames(tpm) <- new_names
  sample_info$new_sample_id <- new_names
  
 # CreateOutput directory
  raw_dir <- file.path(output_dir, "data_raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  
 # Writecounts
  counts_df <- data.frame(gene_id = rownames(counts), counts, check.names = FALSE)
  write.table(counts_df, 
              file.path(raw_dir, paste0(prefix, "_all_counts_with_order.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
 # WriteTPM
  tpm_df <- data.frame(gene_id = rownames(tpm), tpm, check.names = FALSE)
  write.table(tpm_df, 
              file.path(raw_dir, paste0(prefix, "_all_tpm.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
 # Writemetadata
  meta_df <- data.frame(
    sample_id = new_names,
    Treatment = sample_info$Treatment,
    Time = paste0(time_map[sample_info$Time], "d"),
    time_num = sample_info$time_num,
    replicate = sample_info$replicate,
    stringsAsFactors = FALSE
  )
  write.table(meta_df, 
              file.path(raw_dir, paste0(prefix, "_metadata.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
 # WriteGeneAnnotation（DataAnnotation）
  anno_df <- data.frame(
    ensembl_gene_id = rownames(counts),
 hgnc_symbol = rownames(counts), # DataGeneIDsymbol
    stringsAsFactors = FALSE
  )
  write.table(anno_df, 
              file.path(raw_dir, paste0(prefix, "_gene_annotation.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
 # ground truthOutput directory（Evaluation）
  ground_truth <- readRDS(file.path(sim_dir, "ground_truth.rds"))
  dir.create(file.path(output_dir, "data"), recursive = TRUE, showWarnings = FALSE)
  saveRDS(ground_truth, file.path(output_dir, "data", "ground_truth.rds"))
  
  if (verbose) {
    cat(sprintf("  Written: %s_all_counts_with_order.tsv (%d genes × %d samples)\n",
                prefix, nrow(counts), ncol(counts)))
    cat(sprintf("  Written: %s_all_tpm.tsv\n", prefix))
    cat(sprintf("  Written: %s_metadata.tsv\n", prefix))
    cat(sprintf("  Written: %s_gene_annotation.tsv\n", prefix))
    cat(sprintf("  Copied: ground_truth.rds\n"))
    cat(sprintf("  Pipeline command: PROJECT_DIR <- \"%s\"\n", normalizePath(output_dir)))
  }
  
  invisible(list(
    output_dir = output_dir,
    prefix = prefix,
    n_genes = nrow(counts),
    n_samples = ncol(counts)
  ))
}


# ==============================================================================
# Source B: GEOData → pipelineFormat
# ==============================================================================

#' GEOData
#' 
#' @param counts_file countsFile path (TSV/CSV)
#' @param metadata_df SampleDatadata.frame，:
#' - sample_id: RawSample（counts）
#' - Treatment: "Induced" "Control"
#' - Time: Time， "6h", "24h", "48h"
#' - time_num: Time
#' - replicate: Replicate
#' @param output_dir Output directory（pipelinePROJECT_DIR）
#' @param prefix （ "Drug", "Immune"）
#' @param gene_id_col countsGeneID ()
#' @param gene_id_type "ensembl", "entrez", "symbol"

adapt_geo_dataset <- function(counts_file, metadata_df, output_dir, 
                               prefix = "Geo",
                               gene_id_col = NULL,
                               gene_id_type = "ensembl",
                               tpm_file = NULL,
                               verbose = TRUE) {
  
  if (verbose) cat(sprintf("[ADAPT] GEO dataset: %s -> %s\n", 
                           basename(counts_file), output_dir))
  
 # Readcounts
  if (grepl("\\.csv$|\\.csv\\.gz$", counts_file, ignore.case = TRUE)) {
    counts_raw <- read.csv(counts_file, row.names = 1, check.names = FALSE)
  } else {
    counts_raw <- read.delim(counts_file, row.names = 1, check.names = FALSE)
  }
  
 # Sample
  common_samples <- intersect(colnames(counts_raw), metadata_df$sample_id)
  if (length(common_samples) == 0) {
    stop("No matching sample IDs between counts matrix and metadata!")
  }
  
  counts_raw <- counts_raw[, common_samples, drop = FALSE]
  metadata_df <- metadata_df[match(common_samples, metadata_df$sample_id), ]
  
  if (verbose) {
    cat(sprintf("  Matched %d samples (of %d in counts, %d in metadata)\n",
                length(common_samples), ncol(counts_raw), nrow(metadata_df)))
  }
  
 # Sample
 # TimeExtract
  metadata_df$time_label <- gsub("[^0-9]", "", metadata_df$Time)
  
  new_names <- sapply(seq_len(nrow(metadata_df)), function(i) {
    trt <- metadata_df$Treatment[i]
    tl <- metadata_df$time_label[i]
    rep <- metadata_df$replicate[i]
    
    if (trt == "Induced" | trt == "Treated" | trt == "Stimulated") {
      sprintf("%s%sd%d", prefix, tl, rep)
    } else {
      sprintf("%s%sdC%d", prefix, tl, rep)
    }
  })
  
 # Replicate
  if (any(duplicated(new_names))) {
    warning("Duplicate sample names detected! Adding suffix.")
    dup_idx <- which(duplicated(new_names))
    new_names[dup_idx] <- paste0(new_names[dup_idx], "_", seq_along(dup_idx))
  }
  
  colnames(counts_raw) <- new_names
  
 # counts
  counts_int <- round(as.matrix(counts_raw))
  mode(counts_int) <- "integer"
  
 # CalculateTPM（TPMFile）
  if (!is.null(tpm_file) && file.exists(tpm_file)) {
    tpm_raw <- read.delim(tpm_file, row.names = 1, check.names = FALSE)
    tpm_raw <- tpm_raw[, common_samples, drop = FALSE]
    colnames(tpm_raw) <- new_names
    tpm_mat <- as.matrix(tpm_raw)
  } else {
 # TPMCalculate（GeneLength2kb）
    rpk <- counts_int / 2  # 2kb
    tpm_mat <- apply(rpk, 2, function(x) x / sum(x) * 1e6)
  }
  
 # CreateOutput directory
  raw_dir <- file.path(output_dir, "data_raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  
 # Treatment
  metadata_df$Treatment <- ifelse(
    metadata_df$Treatment %in% c("Induced", "Treated", "Stimulated", "Infected"),
    "Induced", "Control"
  )
  
 # WriteFile
  counts_df <- data.frame(gene_id = rownames(counts_int), counts_int, check.names = FALSE)
  write.table(counts_df, 
              file.path(raw_dir, paste0(prefix, "_all_counts_with_order.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  tpm_df <- data.frame(gene_id = rownames(tpm_mat), tpm_mat, check.names = FALSE)
  write.table(tpm_df, 
              file.path(raw_dir, paste0(prefix, "_all_tpm.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  meta_out <- data.frame(
    sample_id = new_names,
    Treatment = metadata_df$Treatment,
    Time = metadata_df$Time,
    time_num = metadata_df$time_num,
    replicate = metadata_df$replicate,
    original_id = common_samples,
    stringsAsFactors = FALSE
  )
  write.table(meta_out, 
              file.path(raw_dir, paste0(prefix, "_metadata.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
 # GeneAnnotation
  anno_df <- data.frame(
    ensembl_gene_id = rownames(counts_int),
    hgnc_symbol = rownames(counts_int),
    gene_id_type = gene_id_type,
    stringsAsFactors = FALSE
  )
  write.table(anno_df, 
              file.path(raw_dir, paste0(prefix, "_gene_annotation.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  if (verbose) {
    cat(sprintf("  Written: %d genes × %d samples\n", nrow(counts_int), ncol(counts_int)))
    cat(sprintf("  Treatment groups: %s\n", 
                paste(names(table(meta_out$Treatment)), 
                      table(meta_out$Treatment), sep = "=", collapse = ", ")))
    cat(sprintf("  Time points: %s\n", 
                paste(sort(unique(meta_out$Time)), collapse = ", ")))
    cat(sprintf("  Pipeline: set PROJECT_DIR <- \"%s\"\n", normalizePath(output_dir)))
  }
  
  invisible(list(
    output_dir = output_dir,
    prefix = prefix,
    n_genes = nrow(counts_int),
    n_samples = ncol(counts_int),
    timepoints = sort(unique(meta_out$Time)),
    sample_table = meta_out
  ))
}


# ==============================================================================
# Function：NCBICalculatecounts
# ==============================================================================

#' NCBIGEODataCalculateRNA-seq counts
#' 
#' @param gse_accession GSE， "GSE123456"
#' @param download_dir Directory

download_ncbi_counts <- function(gse_accession, download_dir = ".") {
  
  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
  
 # NCBICalculatecountsURL
  base_url <- sprintf(
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=%s&format=file&file=%s",
    gse_accession, 
    paste0(gse_accession, "_raw_counts_GRCh38_NCBI.tsv.gz")
  )
  
  dest_file <- file.path(download_dir, 
                          paste0(gse_accession, "_raw_counts_GRCh38_NCBI.tsv.gz"))
  
  cat(sprintf("[DOWNLOAD] Attempting NCBI counts for %s...\n", gse_accession))
  cat(sprintf("  URL: %s\n", base_url))
  
  tryCatch({
    download.file(base_url, dest_file, mode = "wb", quiet = FALSE)
    cat(sprintf("  [OK] Saved to %s\n", dest_file))
    return(dest_file)
  }, error = function(e) {
    cat(sprintf("  [FAIL] %s\n", e$message))
    cat("  Try downloading manually from the GEO page.\n")
    return(NULL)
  })
}


# ==============================================================================
# 
# ==============================================================================

if (sys.nframe() == 0) {
  cat("================================================================\n")
  cat("  S02_pipeline_adapter.R — Usage Guide\n")
  cat("================================================================\n\n")
  
  cat("1. For SIMULATION data:\n")
  cat("   source('S01_simulation_engine.R')\n")
  cat("   sim <- generate_simulation(output_dir = 'sim_data/scenario1')\n")
  cat("   source('S02_pipeline_adapter.R')\n")
  cat("   adapt_simulation('sim_data/scenario1', 'pipeline_runs/Sim1', prefix='Sim')\n\n")
  
  cat("2. For GEO datasets:\n")
  cat("   source('S02_pipeline_adapter.R')\n")
  cat("   # Step 1: Download counts\n")
  cat("   download_ncbi_counts('GSE123456', 'geo_downloads')\n")
  cat("   # Step 2: Manually create metadata (see template below)\n")
  cat("   meta <- data.frame(\n")
  cat("     sample_id = c('GSM_1', 'GSM_2', ...),\n")
  cat("     Treatment = c('Treated', 'Control', ...),\n")
  cat("     Time = c('6h', '6h', '24h', '24h', ...),\n")
  cat("     time_num = c(6, 6, 24, 24, ...),\n")
  cat("     replicate = c(1, 1, 1, 1, ...)\n")
  cat("   )\n")
  cat("   # Step 3: Adapt\n")
  cat("   adapt_geo_dataset('geo_downloads/GSE123456_counts.tsv.gz',\n")
  cat("                     meta, 'pipeline_runs/Drug1', prefix='Drug')\n\n")
  
  cat("3. Then run pipeline:\n")
  cat("   PROJECT_DIR <- 'pipeline_runs/Sim1'  # or Drug1, etc.\n")
  cat("   source('Scripts/00_setup.R')\n")
  cat("   source('Scripts/01_data_import.R')\n")
  cat("   # ... etc.\n\n")
  
  cat("NOTE: For GEO datasets with non-Ensembl gene IDs,\n")
  cat("  01_data_import.R may need minor adjustments.\n")
  cat("  See comments in that script for gene ID mapping.\n")
}
