#!/usr/bin/env Rscript
# ==============================================================================
# S03_download_and_prepare_GEO.R — Download three public datasets and adapt to pipeline format
#
# Datasets:
#   1. GSE197067 — T-cell activation time-series (48 samples, 6 timepoints, 4 donors)
#   2. GSE303975 — Prostate cancer combination therapy (36 samples -> select 18, 3 timepoints, 3 replicates)
#   3. GSE307424 — Lung cancer SMARCA2 inhibitor (18 samples, 3 timepoints, 3 replicates)
#
# Download source:
#   NCBIpre-calculated counts matrices for all human RNA-seq（GRCh38 + NCBI Gene annotation）
#   URLFormat: https://www.ncbi.nlm.nih.gov/geo/download/?acc={GSE}&format=file&file={GSE}_raw_counts_GRCh38_NCBI.tsv.gz
#
# Output:
#   geo_datasets/downloads/          <- Raw downloaded files
#   geo_datasets/metadata/           <- Manually constructed sample metadata
#   pipeline_runs/{DatasetName}/     <- Pipeline-ready formatted data
#
# Usage:
#   setwd(file.path(METHODS_BASE, "Scripts"))
#   source("S03_download_and_prepare_GEO.R")
#
#   # Download and prepare individually:
#   download_all_datasets()       # Step 1: Download counts
#   prepare_GSE197067()           # Step 2: Adapt T-cell data
#   prepare_GSE303975()           # Step 3: Adapt prostate cancer data
#   prepare_GSE307424()           # Step 4: Adapt lung cancer data
#
#   # Or run all at once:
#   run_all()
# ==============================================================================

# ---- 0. Load configuration ----
if (file.exists("S_config.R")) {
  source("S_config.R")
} else if (file.exists(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))) {
  source(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))
} else {
  stop("S_config.R not found! Please run from METHODS_BASE")
}

suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
})

# ---- 1. Dataset registry ----
DATASETS <- list(
  GSE197067 = list(
    gse = "GSE197067",
    title = "T cell activation time course",
    domain = "immune_activation",
    prefix = "Tcell",
    description = "Pan T-cells from 4 healthy donors, anti-CD3/CD28 activated vs non-activated, 6 timepoints"
  ),
  GSE303975 = list(
    gse = "GSE303975",
    title = "Prostate cancer combination therapy",
    domain = "cancer_drug_response",
    prefix = "PCa",
    description = "LNCaP cells, Enzalutamide+Olaparib combination vs DMSO, 3 timepoints"
  ),
  GSE307424 = list(
    gse = "GSE307424",
    title = "Lung cancer SMARCA2 degrader",
    domain = "cancer_drug_response",
    prefix = "Lung",
    description = "NCI-H1693 cells, PRT3789 vs DMSO, 3 timepoints"
  )
)

# ==============================================================================
# Download functions
# ==============================================================================

#' Download NCBI pre-calculated counts for a single dataset
download_geo_counts <- function(gse, download_dir = GEO_DOWNLOAD) {
  dir.create(download_dir, recursive = TRUE, showWarnings = FALSE)
  
  # NCBI pre-calculated counts URL
  filename <- paste0(gse, "_raw_counts_GRCh38_NCBI.tsv.gz")
  url <- sprintf(
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=%s&format=file&file=%s",
    gse, filename
  )
  
  dest_file <- file.path(download_dir, filename)
  
  if (file.exists(dest_file)) {
    cat(sprintf("[SKIP] %s already exists (%.1f MB)\n", 
                filename, file.size(dest_file) / 1e6))
    return(dest_file)
  }
  
  cat(sprintf("[DOWNLOAD] %s ...\n", gse))
  cat(sprintf("  URL: %s\n", url))
  
  tryCatch({
    download.file(url, dest_file, mode = "wb", quiet = FALSE)
    cat(sprintf("  [OK] Saved: %s (%.1f MB)\n", 
                dest_file, file.size(dest_file) / 1e6))
    return(dest_file)
  }, error = function(e) {
    cat(sprintf("  [FAIL] %s\n", e$message))
    cat("\n  === Manual download guide ===\n")
    cat(sprintf("  1. Open in browser: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\n", gse))
    cat("  2. Click at the bottom of the page 'Download Data' or 'Supplementary file'\n")
    cat(sprintf("  3. Download raw counts file to: %s\n", download_dir))
    cat("  4. If no counts file, check 'NCBI-generated RNA-seq counts' link\n\n")
    return(NULL)
  })
}

#' Download all three datasets
download_all_datasets <- function() {
  cat("============================================================\n")
  cat("  Downloading NCBI-generated counts for 3 GEO datasets\n")
  cat("============================================================\n\n")
  
  results <- list()
  for (ds_name in names(DATASETS)) {
    results[[ds_name]] <- download_geo_counts(DATASETS[[ds_name]]$gse)
  }
  
  # Summary
  cat("\n--- Download Summary ---\n")
  for (ds_name in names(results)) {
    status <- if (!is.null(results[[ds_name]]) && file.exists(results[[ds_name]])) "OK" else "FAILED"
    cat(sprintf("  %s: %s\n", ds_name, status))
  }
  
  invisible(results)
}


# ==============================================================================
# General utility functions
# ==============================================================================

#' Read NCBI pre-calculated counts matrix
read_ncbi_counts <- function(gse) {
  filename <- paste0(gse, "_raw_counts_GRCh38_NCBI.tsv.gz")
  filepath <- file.path(GEO_DOWNLOAD, filename)
  
  if (!file.exists(filepath)) {
    stop(sprintf("Counts file not found: %s\nRun download_all_datasets() first.", filepath))
  }
  
  cat(sprintf("[READ] %s ...", filename))
  counts <- read.delim(gzfile(filepath), row.names = 1, check.names = FALSE)
  cat(sprintf(" %d genes × %d samples\n", nrow(counts), ncol(counts)))
  
  return(counts)
}

#' Write files in pipeline format
write_pipeline_files <- function(counts_mat, tpm_mat, metadata_df, gene_anno_df,
                                  output_dir, prefix) {
  
  raw_dir <- file.path(output_dir, "data_raw")
  dir.create(raw_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "data"), recursive = TRUE, showWarnings = FALSE)
  dir.create(file.path(output_dir, "Figure"), recursive = TRUE, showWarnings = FALSE)
  
  # counts
  counts_df <- data.frame(gene_id = rownames(counts_mat), counts_mat, check.names = FALSE)
  write.table(counts_df,
              file.path(raw_dir, paste0(prefix, "_all_counts_with_order.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # TPM
  tpm_df <- data.frame(gene_id = rownames(tpm_mat), tpm_mat, check.names = FALSE)
  write.table(tpm_df,
              file.path(raw_dir, paste0(prefix, "_all_tpm.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # metadata
  write.table(metadata_df,
              file.path(raw_dir, paste0(prefix, "_metadata.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  # gene annotation
  write.table(gene_anno_df,
              file.path(raw_dir, paste0(prefix, "_gene_annotation.tsv")),
              sep = "\t", row.names = FALSE, quote = FALSE)
  
  cat(sprintf("[WRITTEN] %s: %d genes × %d samples -> %s\n",
              prefix, nrow(counts_mat), ncol(counts_mat), raw_dir))
}

#' Calculate simplified TPM from counts (assuming uniform gene length)
compute_simple_tpm <- function(counts_mat, gene_length = 2000) {
  rpk <- counts_mat / (gene_length / 1000)
  tpm <- apply(rpk, 2, function(x) x / sum(x) * 1e6)
  return(tpm)
}

#' Build gene annotation table (from NCBI Gene ID)
build_gene_annotation <- function(gene_ids) {
  # NCBI counts use numeric GeneID, need mapping to symbol
  # Use GeneID as placeholder here; 01_data_import.R will handle mapping at runtime
  data.frame(
    ensembl_gene_id = as.character(gene_ids),
    hgnc_symbol = as.character(gene_ids),
    gene_id_type = "ncbi_gene_id",
    stringsAsFactors = FALSE
  )
}


# ==============================================================================
# GSE197067: T-cell activation time-series
# ==============================================================================
#
# Experimental design (from GEO description):
#   4 healthy donors × 6 timepoints (0, 6, 12, 24, 48, 72h)
#     × 2 conditions (activated, non-activated) = 48 samples
#
# Sample naming convention (needs verification from GEO GSM titles, built from description here):
#   Need to parse exact group information from GEO sample metadata
#   Column names in NCBI counts file are GSM IDs

prepare_GSE197067 <- function() {
  cat("\n============================================================\n")
  cat("  Preparing GSE197067: T cell activation time course\n")
  cat("============================================================\n\n")
  
  ds <- DATASETS$GSE197067
  counts_raw <- read_ncbi_counts(ds$gse)
  
  # ---- Step 1: Get sample information ----
  # Column names in NCBI counts are GSM IDs
  # We need to get experimental conditions for each GSM from GEO
  # Try GEOquery first; if unavailable, provide manual construction template
  
  gsm_ids <- colnames(counts_raw)
  cat(sprintf("  Found %d samples: %s ... %s\n", 
              length(gsm_ids), gsm_ids[1], gsm_ids[length(gsm_ids)]))
  
  # Try to get sample metadata via GEOquery
  metadata_file <- file.path(GEO_METADATA, paste0(ds$gse, "_metadata.rds"))
  
  if (file.exists(metadata_file)) {
    cat("  [CACHE] Loading cached metadata\n")
    sample_meta <- readRDS(metadata_file)
  } else {
    cat("  [FETCH] Downloading sample metadata from GEO...\n")
    cat("  (Need to install GEOquery package: BiocManager::install('GEOquery'))\n\n")
    
    if (!requireNamespace("GEOquery", quietly = TRUE)) {
      cat("  [INFO] GEOquery not installed. Installing...\n")
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager")
      BiocManager::install("GEOquery", update = FALSE, ask = FALSE)
    }
    
    tryCatch({
      gse_data <- GEOquery::getGEO(ds$gse, GSEMatrix = TRUE, getGPL = FALSE)
      
      # Extract phenoData
      if (is.list(gse_data)) {
        pdata <- Biobase::pData(gse_data[[1]])
      } else {
        pdata <- Biobase::pData(gse_data)
      }
      
      # Save raw pdata for debugging
      saveRDS(pdata, file.path(GEO_METADATA, paste0(ds$gse, "_pdata_raw.rds")))
      cat(sprintf("  Raw pdata: %d samples, %d columns\n", nrow(pdata), ncol(pdata)))
      cat("  Columns: ", paste(head(colnames(pdata), 20), collapse = ", "), "\n\n")
      
      # Parse condition information
      # GEO characteristics usually contain condition information
      # Print first few rows for debugging
      cat("  --- Sample preview (first 4) ---\n")
      preview_cols <- grep("title|characteristics|source|description", 
                           colnames(pdata), ignore.case = TRUE, value = TRUE)
      if (length(preview_cols) > 0) {
        print(pdata[1:min(4, nrow(pdata)), preview_cols])
      }
      cat("\n")
      
      # Build standardized metadata
      # Need to parse based on actual pdata column names
      # Save raw data first, then provide manual construction template
      sample_meta <- data.frame(
        gsm_id = rownames(pdata),
        title = pdata$title,
        stringsAsFactors = FALSE
      )
      
      # Try to parse conditions from title
      # GSE197067Expected title format similar to:
      #   "Activated_0h_Donor1" or "non-activated_6h_Donor2" etc.
      # Exact format depends on actual data
      
      # Parse from characteristics column
      char_cols <- grep("characteristics_ch1", colnames(pdata), value = TRUE)
      for (cc in char_cols) {
        sample_meta[[cc]] <- pdata[[cc]]
      }
      
      saveRDS(sample_meta, metadata_file)
      cat(sprintf("  [SAVED] Metadata cache: %s\n", metadata_file))
      
    }, error = function(e) {
      cat(sprintf("  [ERROR] GEOquery failed: %s\n\n", e$message))
      cat("  Please construct metadata manually, steps:\n")
      cat(sprintf("  1. Open https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=%s\n", ds$gse))
      cat("  2. Expand Samples list, record conditions for each GSM\n")
      cat(sprintf("  3. Create file: %s\n", file.path(GEO_METADATA, paste0(ds$gse, "_manual_metadata.csv"))))
      cat("  4. Format: gsm_id, Treatment, Time, time_num, donor_id\n")
      cat("  5. Re-run this function\n\n")
      
      sample_meta <- NULL
    })
  }
  
  if (is.null(sample_meta)) {
    cat("[STOP] Cannot proceed without metadata. See instructions above.\n")
    return(invisible(NULL))
  }
  
  # ---- Step 2: Parse experimental conditions ----
  # Print all metadata for user to verify parsing correctness
  cat("\n  --- Full metadata (please verify) ---\n")
  print(head(sample_meta, 10))
  cat("  ...\n\n")
  
  # Save a CSV template for manual editing
  template_file <- file.path(GEO_METADATA, paste0(ds$gse, "_metadata_template.csv"))
  
  # Build template (user needs to fill in correct Treatment/Time/donor info)
  template <- data.frame(
    gsm_id = colnames(counts_raw),
    Treatment = "TODO_Activated_or_Control",
    Time = "TODO_0h_6h_12h_24h_48h_72h",
    time_num = NA_real_,
    donor_id = "TODO_donor1_2_3_4",
    replicate = NA_integer_,
    stringsAsFactors = FALSE
  )
  
  # If GEOquery info is available, try auto-parsing
  if ("title" %in% colnames(sample_meta)) {
    cat("  [AUTO-PARSE] Attempting to parse conditions from sample titles...\n")
    
    for (i in seq_len(nrow(template))) {
      gsm <- template$gsm_id[i]
      idx <- match(gsm, sample_meta$gsm_id)
      if (!is.na(idx)) {
        ttl <- tolower(sample_meta$title[idx])
        
        # Parse Treatment
        if (grepl("activat|stimulat|anti.cd3|acd3", ttl)) {
          template$Treatment[i] <- "Induced"
        } else if (grepl("non.activ|unactivat|unstimulat|resting|control|no.activ", ttl)) {
          template$Treatment[i] <- "Control"
        }
        
        # Parse time
        time_match <- str_extract(ttl, "\\d+\\s*(h|hr|hour)")
        if (!is.na(time_match)) {
          time_val <- as.numeric(str_extract(time_match, "\\d+"))
          template$Time[i] <- paste0(time_val, "h")
          template$time_num[i] <- time_val
        } else if (grepl("(^|\\s)0\\s*(h|hr)|baseline|before|pre", ttl)) {
          template$Time[i] <- "0h"
          template$time_num[i] <- 0
        }
        
        # Parse donor
        donor_match <- str_extract(ttl, "donor\\s*\\d+|individual\\s*\\d+|subject\\s*\\d+|patient\\s*\\d+|\\bdon?\\d+|\\bind\\d+")
        if (!is.na(donor_match)) {
          template$donor_id[i] <- donor_match
        }
      }
    }
    
    # Check parsing results
    n_parsed_trt <- sum(!grepl("TODO", template$Treatment))
    n_parsed_time <- sum(!grepl("TODO", template$Time))
    cat(sprintf("  Auto-parsed: Treatment=%d/%d, Time=%d/%d\n\n",
                n_parsed_trt, nrow(template), n_parsed_time, nrow(template)))
  }
  
  write.csv(template, template_file, row.names = FALSE)
  cat(sprintf("  [SAVED] Metadata template: %s\n", template_file))
  cat("  >>> Please open this CSV file, check/correct Treatment, Time, donor_id columns <<<\n")
  cat("  >>> Save corrections, then run finalize_GSE197067() <<<\n\n")
  
  invisible(list(counts = counts_raw, template = template, sample_meta = sample_meta))
}


# ==============================================================================
# GSE303975: Prostate cancer combination therapy
# ==============================================================================
#
# Experimental design (confirmed from metadata):
#   4 treatments × 3 timepoints × 3 replicates = 36 samples
#   We only select Combination vs DMSO = 18 samples
#
# Known sample naming, e.g. "DMSO_8hr_1", "Combination_24hr_2", etc.

prepare_GSE303975 <- function() {
  cat("\n============================================================\n")
  cat("  Preparing GSE303975: Prostate cancer combination therapy\n")
  cat("  (Selecting: Combination vs DMSO only)\n")
  cat("============================================================\n\n")
  
  ds <- DATASETS$GSE303975
  counts_raw <- read_ncbi_counts(ds$gse)
  
  gsm_ids <- colnames(counts_raw)
  cat(sprintf("  Total samples in file: %d\n", length(gsm_ids)))
  
  # Get sample titles
  cat("  [FETCH] Getting sample metadata...\n")
  
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("GEOquery", update = FALSE, ask = FALSE)
  }
  
  pdata_file <- file.path(GEO_METADATA, paste0(ds$gse, "_pdata_raw.rds"))
  
  if (file.exists(pdata_file)) {
    pdata <- readRDS(pdata_file)
  } else {
    tryCatch({
      gse_data <- GEOquery::getGEO(ds$gse, GSEMatrix = TRUE, getGPL = FALSE)
      pdata <- Biobase::pData(if (is.list(gse_data)) gse_data[[1]] else gse_data)
      saveRDS(pdata, pdata_file)
    }, error = function(e) {
      cat(sprintf("  [ERROR] %s\n", e$message))
      pdata <- NULL
    })
  }
  
  # Parse groups from sample titles
  # Known format: "Combination_8hr_1", "DMSO_24hr_2", etc.
  if (!is.null(pdata)) {
    titles <- pdata$title
    names(titles) <- rownames(pdata)
  } else {
    cat("  [FALLBACK] Using GSM IDs directly (may need manual mapping)\n")
    titles <- setNames(gsm_ids, gsm_ids)
  }
  
  # Build metadata
  meta_all <- data.frame(
    gsm_id = names(titles),
    title = as.character(titles),
    stringsAsFactors = FALSE
  )
  
  # Parse title
  meta_all$drug <- str_extract(meta_all$title, "^[A-Za-z]+")
  meta_all$time_raw <- str_extract(meta_all$title, "\\d+hr")
  meta_all$time_num <- as.numeric(str_extract(meta_all$time_raw, "\\d+"))
  meta_all$rep <- as.integer(str_extract(meta_all$title, "\\d+$"))
  
  cat("\n  --- All samples parsed ---\n")
  print(table(meta_all$drug, meta_all$time_raw))
  cat("\n")
  
  # Only select Combination + DMSO
  meta_selected <- meta_all %>%
    filter(drug %in% c("Combination", "DMSO")) %>%
    mutate(
      Treatment = ifelse(drug == "DMSO", "Control", "Induced"),
      Time = paste0(time_num, "h"),
      replicate = rep
    )
  
  cat(sprintf("  Selected: %d samples (Combination vs DMSO)\n", nrow(meta_selected)))
  cat("  Design:\n")
  print(table(meta_selected$Treatment, meta_selected$Time))
  cat("\n")
  
  # Subset counts matrix
  counts_sel <- counts_raw[, meta_selected$gsm_id, drop = FALSE]
  
  # Rename samples: PCa8h1, PCa8hC1, etc.
  prefix <- ds$prefix
  new_names <- sapply(seq_len(nrow(meta_selected)), function(i) {
    if (meta_selected$Treatment[i] == "Induced") {
      sprintf("%s%dd%d", prefix, meta_selected$time_num[i], meta_selected$replicate[i])
    } else {
      sprintf("%s%ddC%d", prefix, meta_selected$time_num[i], meta_selected$replicate[i])
    }
  })
  
  colnames(counts_sel) <- new_names
  meta_selected$pipeline_name <- new_names
  
  # CalculateTPM
  tpm_sel <- compute_simple_tpm(as.matrix(counts_sel))
  
  # Gene annotation
  gene_anno <- build_gene_annotation(rownames(counts_sel))
  
  # pipeline metadata
  pipeline_meta <- data.frame(
    sample_id = new_names,
    Treatment = meta_selected$Treatment,
    Time = meta_selected$Time,
    time_num = meta_selected$time_num,
    replicate = meta_selected$replicate,
    original_gsm = meta_selected$gsm_id,
    stringsAsFactors = FALSE
  )
  
  # Write
  output_dir <- file.path(RUN_DIR, "GEO_GSE303975_PCa")
  write_pipeline_files(as.matrix(counts_sel), tpm_sel, pipeline_meta, gene_anno,
                       output_dir, prefix)
  
  # Save full metadata for reference
  write.csv(meta_all, file.path(GEO_METADATA, paste0(ds$gse, "_full_metadata.csv")),
            row.names = FALSE)
  write.csv(meta_selected, file.path(GEO_METADATA, paste0(ds$gse, "_selected_metadata.csv")),
            row.names = FALSE)
  
  cat(sprintf("\n  [DONE] Pipeline-ready data: %s\n", output_dir))
  cat(sprintf("  To run pipeline: PROJECT_DIR <- \"%s\"\n", output_dir))
  
  invisible(list(counts = counts_sel, tpm = tpm_sel, metadata = pipeline_meta))
}


# ==============================================================================
# GSE307424: Lung cancer SMARCA2 inhibitor
# ==============================================================================
#
# Experimental design (confirmed from description):
#   2 treatments (PRT3789 vs DMSO) × 3 timepoints (6h, 48h, 72h) × 3 replicates = 18

prepare_GSE307424 <- function() {
  cat("\n============================================================\n")
  cat("  Preparing GSE307424: Lung cancer SMARCA2 degrader\n")
  cat("============================================================\n\n")
  
  ds <- DATASETS$GSE307424
  counts_raw <- read_ncbi_counts(ds$gse)
  
  gsm_ids <- colnames(counts_raw)
  cat(sprintf("  Total samples: %d\n", length(gsm_ids)))
  
  # Get sample metadata
  pdata_file <- file.path(GEO_METADATA, paste0(ds$gse, "_pdata_raw.rds"))
  
  if (!requireNamespace("GEOquery", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    BiocManager::install("GEOquery", update = FALSE, ask = FALSE)
  }
  
  if (file.exists(pdata_file)) {
    pdata <- readRDS(pdata_file)
  } else {
    tryCatch({
      gse_data <- GEOquery::getGEO(ds$gse, GSEMatrix = TRUE, getGPL = FALSE)
      pdata <- Biobase::pData(if (is.list(gse_data)) gse_data[[1]] else gse_data)
      saveRDS(pdata, pdata_file)
    }, error = function(e) {
      cat(sprintf("  [ERROR] %s\n", e$message))
      pdata <- NULL
    })
  }
  
  if (!is.null(pdata)) {
    titles <- pdata$title
    names(titles) <- rownames(pdata)
    
    cat("  --- Sample titles ---\n")
    print(sort(titles))
    cat("\n")
  }
  
  # Build metadata
  # Expected title format similar to: "PRT3789_6h_rep1", "DMSO_48h_rep2"  etc.
  # Need to adjust parsing logic based on actual title format
  
  meta <- data.frame(
    gsm_id = names(titles),
    title = as.character(titles),
    stringsAsFactors = FALSE
  )
  
  # Parse based on common naming patterns
  meta$drug <- ifelse(grepl("PRT|prt|drug|treated", meta$title, ignore.case = TRUE),
                      "PRT3789",
                      ifelse(grepl("DMSO|dmso|control|vehicle", meta$title, ignore.case = TRUE),
                             "DMSO", "UNKNOWN"))
  
  meta$time_raw <- str_extract(meta$title, "\\d+\\s*(h|hr|hour)")
  meta$time_num <- as.numeric(str_extract(meta$time_raw, "\\d+"))
  meta$rep <- as.integer(str_extract(meta$title, "(rep|r|_)\\s*(\\d+)$"))
  
  # If rep not parsed, try assigning by within-group order
  if (all(is.na(meta$rep))) {
    meta <- meta %>%
      group_by(drug, time_num) %>%
      mutate(rep = row_number()) %>%
      ungroup() %>%
      as.data.frame()
  }
  
  meta$Treatment <- ifelse(meta$drug == "DMSO", "Control", "Induced")
  meta$Time <- paste0(meta$time_num, "h")
  meta$replicate <- meta$rep
  
  cat("  --- Parsed design ---\n")
  print(table(meta$Treatment, meta$Time))
  cat("\n")
  
  # Rename
  prefix <- ds$prefix
  new_names <- sapply(seq_len(nrow(meta)), function(i) {
    if (meta$Treatment[i] == "Induced") {
      sprintf("%s%dd%d", prefix, meta$time_num[i], meta$replicate[i])
    } else {
      sprintf("%s%ddC%d", prefix, meta$time_num[i], meta$replicate[i])
    }
  })
  
  counts_mat <- as.matrix(counts_raw[, meta$gsm_id, drop = FALSE])
  colnames(counts_mat) <- new_names
  
  tpm_mat <- compute_simple_tpm(counts_mat)
  gene_anno <- build_gene_annotation(rownames(counts_mat))
  
  pipeline_meta <- data.frame(
    sample_id = new_names,
    Treatment = meta$Treatment,
    Time = meta$Time,
    time_num = meta$time_num,
    replicate = meta$replicate,
    original_gsm = meta$gsm_id,
    stringsAsFactors = FALSE
  )
  
  output_dir <- file.path(RUN_DIR, "GEO_GSE307424_Lung")
  write_pipeline_files(counts_mat, tpm_mat, pipeline_meta, gene_anno,
                       output_dir, prefix)
  
  write.csv(meta, file.path(GEO_METADATA, paste0(ds$gse, "_metadata.csv")),
            row.names = FALSE)
  
  cat(sprintf("\n  [DONE] Pipeline-ready data: %s\n", output_dir))
  cat(sprintf("  To run pipeline: PROJECT_DIR <- \"%s\"\n", output_dir))
  
  invisible(list(counts = counts_mat, tpm = tpm_mat, metadata = pipeline_meta))
}


# ==============================================================================
# Finalization function: call after user confirms metadata
# ==============================================================================

#' Finalize adaptation after user manually corrects GSE197067 metadata CSV
finalize_GSE197067 <- function() {
  cat("\n============================================================\n")
  cat("  Finalizing GSE197067: T cell activation time course\n")
  cat("============================================================\n\n")
  
  ds <- DATASETS$GSE197067
  
  # Read user-corrected metadata
  template_file <- file.path(GEO_METADATA, paste0(ds$gse, "_metadata_template.csv"))
  if (!file.exists(template_file)) {
    stop(sprintf("Metadata template not found: %s\nRun prepare_GSE197067() first.", template_file))
  }
  
  meta <- read.csv(template_file, stringsAsFactors = FALSE)
  
  # Check for remaining TODOs
  n_todo <- sum(grepl("TODO", meta$Treatment)) + sum(grepl("TODO", meta$Time))
  if (n_todo > 0) {
    cat(sprintf("  [WARN] %d fields still contain 'TODO'!\n", n_todo))
    cat("  Please edit the CSV file and replace all TODO values.\n")
    cat(sprintf("  File: %s\n\n", template_file))
    
    # If auto-parse filled most, check remaining TODO count
    n_total <- nrow(meta) * 2  # Treatment + Time
    if (n_todo > n_total * 0.5) {
      cat("  [STOP] Too many unresolved fields. Cannot proceed.\n")
      return(invisible(NULL))
    } else {
      cat("  [WARN] Proceeding with partial metadata (some samples may be excluded).\n")
      meta <- meta[!grepl("TODO", meta$Treatment) & !grepl("TODO", meta$Time), ]
      cat(sprintf("  Keeping %d samples with resolved metadata.\n\n", nrow(meta)))
    }
  }
  
  # Read counts
  counts_raw <- read_ncbi_counts(ds$gse)
  
  # Filter to samples with metadata
  common <- intersect(colnames(counts_raw), meta$gsm_id)
  cat(sprintf("  Matched %d samples between counts and metadata\n", length(common)))
  
  counts_sel <- as.matrix(counts_raw[, common, drop = FALSE])
  meta <- meta[match(common, meta$gsm_id), ]
  
  # Assign replicate IDs (by donor number)
  if (all(is.na(meta$replicate) | meta$replicate == 0)) {
    donors <- unique(meta$donor_id)
    meta$replicate <- match(meta$donor_id, donors)
  }
  
  # Rename
  prefix <- ds$prefix
  new_names <- sapply(seq_len(nrow(meta)), function(i) {
    tn <- meta$time_num[i]
    rep <- meta$replicate[i]
    if (meta$Treatment[i] == "Induced") {
      sprintf("%s%dd%d", prefix, tn, rep)
    } else {
      sprintf("%s%ddC%d", prefix, tn, rep)
    }
  })
  
  # Check for duplicate names
  if (any(duplicated(new_names))) {
    cat("  [WARN] Duplicate names detected, adding suffix\n")
    dup <- duplicated(new_names) | duplicated(new_names, fromLast = TRUE)
    new_names[dup] <- paste0(new_names[dup], "_", seq_along(which(dup)))
  }
  
  colnames(counts_sel) <- new_names
  
  tpm_sel <- compute_simple_tpm(counts_sel)
  gene_anno <- build_gene_annotation(rownames(counts_sel))
  
  pipeline_meta <- data.frame(
    sample_id = new_names,
    Treatment = meta$Treatment,
    Time = meta$Time,
    time_num = meta$time_num,
    replicate = meta$replicate,
    original_gsm = meta$gsm_id,
    donor_id = meta$donor_id,
    stringsAsFactors = FALSE
  )
  
  output_dir <- file.path(RUN_DIR, "GEO_GSE197067_Tcell")
  write_pipeline_files(counts_sel, tpm_sel, pipeline_meta, gene_anno,
                       output_dir, prefix)
  
  cat(sprintf("\n  [DONE] Pipeline-ready data: %s\n", output_dir))
  cat(sprintf("  Design: %d samples\n", ncol(counts_sel)))
  print(table(pipeline_meta$Treatment, pipeline_meta$Time))
  cat(sprintf("\n  To run pipeline: PROJECT_DIR <- \"%s\"\n", output_dir))
  
  invisible(list(counts = counts_sel, tpm = tpm_sel, metadata = pipeline_meta))
}


# ==============================================================================
# Run all at once
# ==============================================================================

run_all <- function() {
  cat("################################################################\n")
  cat("#  METI-FS Methods Paper: GEO Dataset Download & Preparation  #\n")
  cat("################################################################\n\n")
  
  # Step 1: Download
  cat("=== STEP 1/4: Download counts matrices ===\n")
  download_all_datasets()
  
  # Step 2: Prepare GSE303975 (simplest, known complete naming format)
  cat("\n=== STEP 2/4: Prepare GSE303975 (Prostate cancer) ===\n")
  tryCatch(prepare_GSE303975(), error = function(e) {
    cat(sprintf("  [ERROR] %s\n", e$message))
  })
  
  # Step 3: Prepare GSE307424
  cat("\n=== STEP 3/4: Prepare GSE307424 (Lung cancer) ===\n")
  tryCatch(prepare_GSE307424(), error = function(e) {
    cat(sprintf("  [ERROR] %s\n", e$message))
  })
  
  # Step 4: Prepare GSE197067 (requires manual metadata verification)
  cat("\n=== STEP 4/4: Prepare GSE197067 (T cell) ===\n")
  cat("  NOTE: This dataset requires manual metadata verification.\n")
  cat("  Running prepare_GSE197067() to generate template...\n\n")
  tryCatch(prepare_GSE197067(), error = function(e) {
    cat(sprintf("  [ERROR] %s\n", e$message))
  })
  
  cat("\n################################################################\n")
  cat("#  Download & preparation complete!                           #\n")
  cat("#                                                             #\n")
  cat("#  Next steps:                                                #\n")
  cat("#  1. Check GSE303975 and GSE307424 outputs in pipeline_runs/ #\n")
  cat("#  2. For GSE197067: edit the metadata template CSV,          #\n")
  cat("#     then run finalize_GSE197067()                           #\n")
  cat("#  3. Run pipeline on each dataset:                           #\n")
  cat("#     PROJECT_DIR <- 'METHODS_BASE'   #\n")
  cat("#     source('R/00_setup.R')        #\n")
  cat("#     source('R/01_data_import.R')  #\n")
  cat("#     ...                                                     #\n")
  cat("################################################################\n")
}


# ==============================================================================
# If running this script directly
# ==============================================================================
if (sys.nframe() == 0) {
  run_all()
}
