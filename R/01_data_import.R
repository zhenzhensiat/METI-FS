#!/usr/bin/env Rscript
# ==============================================================================
# 01_data_import.R — Import and SampleInfo
# ==============================================================================

source(file.path(SCRIPT_DIR, "00_setup.R"))
source(file.path(SCRIPT_DIR, "theme_bindlab.R"))

log_step("01_IMPORT", "Starting data import...")

# ============================================================================
# 1. ImportRaw data
# ============================================================================

# CountsMatrix（ using ：DESeq2Differential analysis、WGCNA、maSigPro）
counts_raw <- read.delim(FILES$counts_raw, row.names = 1, check.names = FALSE)
log_step("01_IMPORT", sprintf("Counts loaded: %d genes × %d samples", 
                               nrow(counts_raw), ncol(counts_raw)))

# TPMMatrix（ using ：FigureVisualization、GeneExpression， not using Statistics）
tpm_raw <- read.delim(FILES$tpm_raw, row.names = 1, check.names = FALSE)
log_step("01_IMPORT", sprintf("TPM loaded: %d genes × %d samples", 
                               nrow(tpm_raw), ncol(tpm_raw)))

# 
stopifnot("Counts and TPM gene count mismatch" = nrow(counts_raw) == nrow(tpm_raw))
stopifnot("Counts and TPM sample mismatch" = all(colnames(counts_raw) == colnames(tpm_raw)))
stopifnot("Counts must be integers" = all(counts_raw == floor(counts_raw)))
log_step("01_IMPORT", "Data integrity check PASSED")

# ============================================================================
# 2. SampleInfo (colData)
# ============================================================================

# Sample then （ using ，）：
# Induced group: {Prefix}{Time}d{Rep} e.g., Sample4d1, Sample14d2
# Control group: {Prefix}{Time}dC{Rep} e.g., Sample4dC1, Sample7dC3
# before from countscolumns

sample_names <- colnames(counts_raw)

# Sample before ： Sample， d/dC
.sample_prefix <- sub("\\d+d(C)?\\d+$", "", sample_names[1])
log_step("01_IMPORT", sprintf("Auto-detected sample prefix: '%s'", .sample_prefix))

# then Expression
.re_ctrl    <- sprintf("^%s(\\d+)dC(\\d+)$", .sample_prefix)
.re_induced <- sprintf("^%s(\\d+)d(\\d+)$", .sample_prefix)

# Treatment and Time
parse_sample <- function(name) {
  if (grepl("C\\d+$", name)) {
    # Control group
    treatment <- "Control"
    time_str <- sub(.re_ctrl, "\\1", name)
    rep_str  <- sub(.re_ctrl, "\\2", name)
  } else {
    # Induced group
    treatment <- "Induced"
    time_str <- sub(.re_induced, "\\1", name)
    rep_str  <- sub(.re_induced, "\\2", name)
  }
  return(c(treatment = treatment, 
           time_num = as.integer(time_str),
           time_label = paste0(time_str, "d"),
           replicate = as.integer(rep_str)))
}

sample_info_list <- lapply(sample_names, parse_sample)
sample_info <- as.data.frame(do.call(rbind, sample_info_list), stringsAsFactors = FALSE)
rownames(sample_info) <- sample_names
sample_info$time_num <- as.integer(sample_info$time_num)
sample_info$replicate <- as.integer(sample_info$replicate)

# （ and Reference）
sample_info$Treatment <- factor(sample_info$treatment, levels = c("Control", "Induced"))
sample_info$Time <- factor(sample_info$time_label, levels = PARAMS$time_labels)
sample_info$Group <- factor(paste(sample_info$treatment, sample_info$time_label, sep = "_"),
                            levels = c("Control_4d", "Induced_4d",
                                       "Control_7d", "Induced_7d",
                                       "Control_14d", "Induced_14d",
                                       "Control_21d", "Induced_21d"))

# ValidationSampleInfo
log_step("01_IMPORT", "Sample information constructed:")
print(table(sample_info$Treatment, sample_info$Time))

# Save
save_data(sample_info, FILES$sample_info)

# ============================================================================
# 3. GeneAnnotation（Ensembl ID → Gene Symbol）
# ============================================================================
log_step("01_IMPORT", "Building gene annotation table...")

ensembl_ids <- rownames(counts_raw)

# using gene_annotation.tsv（Osteogenic/Chondrogenic has ）
if (file.exists(FILES$gene_anno_raw)) {
  log_step("01_IMPORT", sprintf("Loading local gene annotation: %s", FILES$gene_anno_raw))
  local_anno <- read.delim(FILES$gene_anno_raw, stringsAsFactors = FALSE)
  
  log_step("01_IMPORT", sprintf("Local annotation columns: %s", paste(colnames(local_anno), collapse = ", ")))
  log_step("01_IMPORT", sprintf("Local annotation rows: %d", nrow(local_anno)))
  
 # can columnsNormalization： has 
  # Ensembl IDcolumns
  ens_candidates <- c("gene_id", "ensembl_id", "Geneid", "GeneID", "gene", "Gene", "ENSEMBL")
  for (ec in ens_candidates) {
    if (ec %in% colnames(local_anno) && !"ensembl_gene_id" %in% colnames(local_anno)) {
      colnames(local_anno)[colnames(local_anno) == ec] <- "ensembl_gene_id"
      log_step("01_IMPORT", sprintf("Renamed column '%s' -> 'ensembl_gene_id'", ec))
      break
    }
  }
  
  # Gene Symbolcolumns
  sym_candidates <- c("gene_name", "symbol", "gene_symbol", "GeneName", "Symbol", 
                       "SYMBOL", "external_gene_name", "hgnc_symbol", "gene.name")
  for (sc in sym_candidates) {
    if (sc %in% colnames(local_anno) && !"hgnc_symbol" %in% colnames(local_anno)) {
      colnames(local_anno)[colnames(local_anno) == sc] <- "hgnc_symbol"
      log_step("01_IMPORT", sprintf("Renamed column '%s' -> 'hgnc_symbol'", sc))
      break
    }
  }
  
 # Validation need to columns in 
  if (!"ensembl_gene_id" %in% colnames(local_anno)) {
    log_step("01_IMPORT", sprintf("WARNING: Cannot find Ensembl ID column in local annotation! Available: %s",
                                   paste(colnames(local_anno), collapse = ", ")))
    log_step("01_IMPORT", "Falling back to online annotation...")
    local_anno <- NULL  # Trigger online retrieval below
  }
  
  if (!is.null(local_anno) && !"hgnc_symbol" %in% colnames(local_anno)) {
    log_step("01_IMPORT", sprintf("WARNING: Cannot find Symbol column in local annotation! Available: %s",
                                   paste(colnames(local_anno), collapse = ", ")))
    log_step("01_IMPORT", "Falling back to online annotation...")
    local_anno <- NULL
  }
}

if (file.exists(FILES$gene_anno_raw) && !is.null(local_anno)) {
 # using 
  gene_anno <- data.frame(ensembl_gene_id = ensembl_ids, stringsAsFactors = FALSE)
  gene_anno$hgnc_symbol <- local_anno$hgnc_symbol[match(gene_anno$ensembl_gene_id, local_anno$ensembl_gene_id)]
  
  n_matched <- sum(!is.na(gene_anno$hgnc_symbol))
  log_step("01_IMPORT", sprintf("Local annotation matched: %d/%d (%.1f%%)", 
                                 n_matched, length(ensembl_ids), 100*n_matched/length(ensembl_ids)))
  
 # Entrez ID
  ens2entrez <- tryCatch(
    AnnotationDbi::select(org.Hs.eg.db, keys = ensembl_ids,
                           columns = c("ENTREZID"), keytype = "ENSEMBL"),
    error = function(e) { message("Entrez mapping failed: ", e$message); NULL }
  )
  if (!is.null(ens2entrez)) {
    ens2entrez <- ens2entrez[!duplicated(ens2entrez$ENSEMBL), ]
    gene_anno$entrez_id <- ens2entrez$ENTREZID[match(gene_anno$ensembl_gene_id, ens2entrez$ENSEMBL)]
  } else {
    gene_anno$entrez_id <- NA
  }
  
} else {
 # in （Adipogenic has Logic， or columns not ）
  log_step("01_IMPORT", "Using org.Hs.eg.db online mapping...")
  ens2entrez <- AnnotationDbi::select(org.Hs.eg.db,
                                       keys = ensembl_ids,
                                       columns = c("ENTREZID", "SYMBOL"),
                                       keytype = "ENSEMBL")
  ens2entrez <- ens2entrez[!duplicated(ens2entrez$ENSEMBL), ]
  
  gene_anno <- data.frame(ensembl_gene_id = ensembl_ids, stringsAsFactors = FALSE)
  gene_anno$hgnc_symbol <- ens2entrez$SYMBOL[match(gene_anno$ensembl_gene_id, ens2entrez$ENSEMBL)]
  gene_anno$entrez_id   <- ens2entrez$ENTREZID[match(gene_anno$ensembl_gene_id, ens2entrez$ENSEMBL)]
}

# StatisticsMapping
n_mapped_symbol <- sum(!is.na(gene_anno$hgnc_symbol))
n_mapped_entrez <- sum(!is.na(gene_anno$entrez_id))
log_step("01_IMPORT", sprintf("Gene annotation: %d/%d (%.1f%%) mapped to Symbol, %d/%d (%.1f%%) to Entrez",
                               n_mapped_symbol, nrow(gene_anno), 100*n_mapped_symbol/nrow(gene_anno),
                               n_mapped_entrez, nrow(gene_anno), 100*n_mapped_entrez/nrow(gene_anno)))

save_data(gene_anno, FILES$gene_annotation)

# ============================================================================
# 4. 
# ============================================================================
log_step("01_IMPORT", "=== Data Overview ===")
log_step("01_IMPORT", sprintf("Total genes: %d", nrow(counts_raw)))
log_step("01_IMPORT", sprintf("Total samples: %d (Induced: %d, Control: %d)",
                               ncol(counts_raw),
                               sum(sample_info$Treatment == "Induced"),
                               sum(sample_info$Treatment == "Control")))
log_step("01_IMPORT", sprintf("Library sizes: %.1fM - %.1fM (median: %.1fM)",
                               min(colSums(counts_raw))/1e6,
                               max(colSums(counts_raw))/1e6,
                               median(colSums(counts_raw))/1e6))

log_step("01_IMPORT", "Step 01 COMPLETE")
