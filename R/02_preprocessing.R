#!/usr/bin/env Rscript
# ==============================================================================
# 02_preprocessing.R — GeneFilter and QC
# Key principles:
# - Filter using counts，Filter
# - After filtering counts and TPM step Subset
# - ， has by Delete GeneInfo
# ==============================================================================

source(file.path(SCRIPT_DIR, "00_setup.R"))
source(file.path(SCRIPT_DIR, "theme_bindlab.R"))

# Figure can using 
if (!is.null(dev.list())) dev.off()

log_step("02_PREPROC", "Starting preprocessing...")

# Load data
counts_raw  <- read.delim(FILES$counts_raw, row.names = 1, check.names = FALSE)
tpm_raw     <- read.delim(FILES$tpm_raw, row.names = 1, check.names = FALSE)
sample_info <- readRDS(FILES$sample_info)
gene_anno   <- readRDS(FILES$gene_annotation)

# ============================================================================
# 1. GeneFilter（ using edgeR::filterByExpr — Standard）
# ============================================================================
# filterByExpr in group Sample size，Calculate CPMThreshold
# Reference：Chen et al. (2016) F1000Research; Michael Love and DESeq2 using 

# group（8 group Grouping）
group_vec <- factor(paste(sample_info$Treatment, sample_info$Time, sep = "_"))

log_step("02_PREPROC", "Using edgeR::filterByExpr for design-aware filtering...")
log_step("02_PREPROC", sprintf("Groups: %s", paste(levels(group_vec), collapse = ", ")))
log_step("02_PREPROC", sprintf("Min group size: %d", min(table(group_vec))))

# filterByExpr: in group in has Expression Gene
keep_genes <- edgeR::filterByExpr(counts_raw, group = group_vec)

# using Comparison
gene_detect_count <- rowSums(counts_raw >= PARAMS$min_count)
gene_total_count  <- rowSums(counts_raw)
gene_max_count    <- apply(counts_raw, 1, max)
gene_mean_count   <- rowMeans(counts_raw)

# Filter
filter_log <- data.frame(
  ensembl_id      = rownames(counts_raw),
  symbol          = gene_anno$hgnc_symbol[match(rownames(counts_raw), gene_anno$ensembl_gene_id)],
  total_count     = gene_total_count,
  max_count       = gene_max_count,
  mean_count      = round(gene_mean_count, 2),
  n_samples_detected = gene_detect_count,
  kept            = keep_genes,
  filter_reason   = ifelse(keep_genes, "PASSED",
                    ifelse(gene_total_count == 0, "ALL_ZERO",
                    ifelse(gene_max_count < 1, "ALL_ZERO_OR_NEAR",
                           "LOW_EXPRESSION")))
)

# SaveFilter
write.csv(filter_log, FILES$filter_log, row.names = FALSE)

# Statistics
n_total    <- nrow(counts_raw)
n_kept     <- sum(keep_genes)
n_removed  <- n_total - n_kept
n_allzero  <- sum(filter_log$filter_reason == "ALL_ZERO" | 
                  filter_log$filter_reason == "ALL_ZERO_OR_NEAR")
n_lowexpr  <- sum(filter_log$filter_reason == "LOW_EXPRESSION")

log_step("02_PREPROC", "=== Filter Summary (filterByExpr) ===")
log_step("02_PREPROC", sprintf("Total genes:        %d", n_total))
log_step("02_PREPROC", sprintf("Kept genes:         %d (%.1f%%)", n_kept, 100*n_kept/n_total))
log_step("02_PREPROC", sprintf("Removed genes:      %d (%.1f%%)", n_removed, 100*n_removed/n_total))
log_step("02_PREPROC", sprintf("  - All/near zeros: %d", n_allzero))
log_step("02_PREPROC", sprintf("  - Low expression: %d", n_lowexpr))

# should using Filter
counts_filtered <- counts_raw[keep_genes, ]
tpm_filtered    <- tpm_raw[keep_genes, ]

# Save
save_data(counts_filtered, FILES$counts_filtered)
save_data(tpm_filtered, FILES$tpm_filtered)

# ============================================================================
# 2. QCVisualization
# ============================================================================

# --- 2a. ---
lib_sizes <- data.frame(
  sample = colnames(counts_filtered),
  lib_size = colSums(counts_filtered) / 1e6,
  Treatment = sample_info$Treatment,
  Time = sample_info$Time
)

p_libsize <- ggplot(lib_sizes, aes(x = sample, y = lib_size, fill = Treatment)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_treatment() +
  labs(title = "Library Size Distribution",
       subtitle = "After gene filtering",
       x = NULL, y = "Library Size (Millions)") +
  theme_bindlab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

save_pub_fig(p_libsize, "QC_library_size", "01_QC", width = 12, height = 5)

# --- 2b. Gene count ---
gene_detect <- data.frame(
  sample = colnames(counts_filtered),
  n_genes = colSums(counts_filtered > 0),
  Treatment = sample_info$Treatment,
  Time = sample_info$Time
)

p_ngenes <- ggplot(gene_detect, aes(x = sample, y = n_genes, fill = Treatment)) +
  geom_bar(stat = "identity", width = 0.7) +
  scale_fill_treatment() +
  labs(title = "Number of Detected Genes per Sample",
       subtitle = paste("Genes with count > 0 after filtering (", n_kept, "total genes)"),
       x = NULL, y = "Detected Genes") +
  theme_bindlab() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8))

save_pub_fig(p_ngenes, "QC_detected_genes", "01_QC", width = 12, height = 5)

# --- 2c. Counts（DensityFigure） ---
log_counts <- log2(counts_filtered + 1)
log_counts_long <- log_counts %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "sample", values_to = "log2count") %>%
  left_join(data.frame(sample = lib_sizes$sample, 
                       Treatment = lib_sizes$Treatment, 
                       Time = lib_sizes$Time,
                       stringsAsFactors = FALSE), by = "sample")

p_density <- ggplot(log_counts_long, aes(x = log2count, color = sample)) +
  geom_density(linewidth = 0.4, show.legend = FALSE) +
  facet_wrap(~Time, nrow = 1) +
  labs(title = "Gene Expression Density Distribution",
       subtitle = "log2(count + 1)",
       x = "log2(count + 1)", y = "Density") +
  theme_bindlab()

save_pub_fig(p_density, "QC_density_distribution", "01_QC", width = 14, height = 4)

# --- 2d. Before filtering after ---
filter_summary_df <- data.frame(
  Category = c("Passed", "All/Near Zeros", "Low Expression"),
  Count = c(n_kept, n_allzero, n_lowexpr),
  stringsAsFactors = FALSE
)
filter_summary_df$Category <- factor(filter_summary_df$Category, 
                                     levels = rev(filter_summary_df$Category))

p_filter <- ggplot(filter_summary_df, aes(x = Category, y = Count, fill = Category)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = Count), hjust = -0.1, size = 3.5) +
  scale_fill_manual(values = c("Passed" = "#00A087", "All/Near Zeros" = "#E64B35",
                                "Low Expression" = "#4DBBD5")) +
  coord_flip() +
  labs(title = "Gene Filtering Summary (filterByExpr)",
       subtitle = "Design-aware filtering by edgeR::filterByExpr",
       x = NULL, y = "Number of Genes") +
  theme_bindlab() +
  theme(legend.position = "none") +
  expand_limits(y = max(filter_summary_df$Count) * 1.15)

save_pub_fig(p_filter, "QC_filter_summary", "01_QC", width = 8, height = 4)

log_step("02_PREPROC", "Step 02 COMPLETE")
