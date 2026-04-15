#!/usr/bin/env Rscript
# ==============================================================================
# 03_normalization_QC.R — DESeq2Normalization、PCA、CorrelationFigure、SampleClustering
# Key principles:
# - PCA using VST after counts（TPM）
# - VST/rlog using Visualization
# - NormalizationDESeq2 within 
# ==============================================================================

source(file.path(SCRIPT_DIR, "00_setup.R"))
source(file.path(SCRIPT_DIR, "theme_bindlab.R"))

log_step("03_NORM", "Starting normalization and sample-level QC...")

# Load data
counts_filtered <- readRDS(FILES$counts_filtered)
sample_info     <- readRDS(FILES$sample_info)

# ============================================================================
# 1. DESeqDataSet
# ============================================================================
# counts is Matrix
counts_mat <- as.matrix(counts_filtered)
storage.mode(counts_mat) <- "integer"

# Sample
stopifnot(all(colnames(counts_mat) == rownames(sample_info)))

# Build DESeq2 object（：Treatment + Time + item ）
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat,
  colData   = sample_info,
  design    = ~ Treatment + Time + Treatment:Time
)

# Reference
dds$Treatment <- relevel(dds$Treatment, ref = "Control")
dds$Time      <- relevel(dds$Time, ref = "4d")

log_step("03_NORM", sprintf("DESeqDataSet created: %d genes × %d samples",
                             nrow(dds), ncol(dds)))

# ============================================================================
# 2. VST（ using Visualization，blind=FALSE is has Info）
# ============================================================================
# blind=FALSE: using Info， already 
# Note:28 Sample > 30Gene using VSTrlog（）
vsd <- vst(dds, blind = FALSE)
vst_mat <- assay(vsd)

log_step("03_NORM", "VST transformation complete (blind=FALSE)")

# Save
save_data(dds, FILES$dds_object)
save_data(vst_mat, FILES$vst_matrix)

# ============================================================================
# 3. PCA（ using VST after counts）
# ============================================================================
log_step("03_NORM", "PCA analysis on VST-transformed counts...")

# CalculatePCA（top 500 most variable genes）
ntop <- 500
rv <- rowVars(vst_mat)
select_genes <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
pca_data <- prcomp(t(vst_mat[select_genes, ]), center = TRUE, scale. = FALSE)

# PC and 
pca_df <- as.data.frame(pca_data$x[, 1:5])
pca_df$sample    <- rownames(pca_df)
pca_df$Treatment <- sample_info$Treatment
pca_df$Time      <- sample_info$Time
pca_df$Group     <- sample_info$Group

pct_var <- round(100 * (pca_data$sdev^2 / sum(pca_data$sdev^2)), 1)

# --- PCAFigure: PC1 vs PC2, Treatment，Time ---
p_pca_main <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Treatment, shape = Time), size = 3.5, stroke = 0.8) +
  scale_color_treatment() +
  scale_shape_manual(values = c("4d" = 16, "7d" = 17, "14d" = 15, "21d" = 18)) +
  stat_ellipse(aes(group = Treatment, color = Treatment), 
               type = "t", level = 0.95, linetype = 2, linewidth = 0.5) +
  labs(title = "PCA of Gene Expression",
       subtitle = sprintf("Top %d variable genes, VST-transformed counts", ntop),
       x = sprintf("PC1 (%s%% variance)", pct_var[1]),
       y = sprintf("PC2 (%s%% variance)", pct_var[2])) +
  theme_bindlab_minimal()

save_pub_fig(p_pca_main, "PCA_Treatment_Time", "02_PCA_Clustering", width = 8, height = 6)

# --- PCAFigure: PC1 vs PC2, Group ---
p_pca_group <- ggplot(pca_df, aes(x = PC1, y = PC2)) +
  geom_point(aes(color = Group), size = 3.5) +
  scale_color_manual(values = COLORS$group) +
  ggrepel::geom_text_repel(aes(label = sample), size = 2.2, max.overlaps = 20,
                            segment.color = "grey70", segment.size = 0.3) +
  labs(title = "PCA — All Groups",
       subtitle = sprintf("Top %d variable genes", ntop),
       x = sprintf("PC1 (%s%% variance)", pct_var[1]),
       y = sprintf("PC2 (%s%% variance)", pct_var[2])) +
  theme_bindlab_minimal()

save_pub_fig(p_pca_group, "PCA_AllGroups_labeled", "02_PCA_Clustering", width = 10, height = 7)

# --- PCAFigure: PC2 vs PC3 ---
p_pca_23 <- ggplot(pca_df, aes(x = PC2, y = PC3)) +
  geom_point(aes(color = Treatment, shape = Time), size = 3.5, stroke = 0.8) +
  scale_color_treatment() +
  scale_shape_manual(values = c("4d" = 16, "7d" = 17, "14d" = 15, "21d" = 18)) +
  labs(title = "PCA — PC2 vs PC3",
       x = sprintf("PC2 (%s%% variance)", pct_var[2]),
       y = sprintf("PC3 (%s%% variance)", pct_var[3])) +
  theme_bindlab_minimal()

save_pub_fig(p_pca_23, "PCA_PC2_PC3", "02_PCA_Clustering", width = 8, height = 6)

# --- Figure： ---
scree_df <- data.frame(PC = paste0("PC", 1:10), Variance = pct_var[1:10])
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

# ============================================================================
# 4. SampleCorrelationFigure
# ============================================================================
log_step("03_NORM", "Sample correlation heatmap...")

# Spearman（Expression）
cor_mat <- cor(vst_mat, method = "spearman")

# AnnotationInfo
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
    print(pheatmap(cor_mat,
             clustering_distance_rows = as.dist(1 - cor_mat),
             clustering_distance_cols = as.dist(1 - cor_mat),
             clustering_method = "complete",
             color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100),
             breaks = seq(min(cor_mat), 1, length.out = 101),
             annotation_col = anno_col,
             annotation_colors = anno_colors,
             show_rownames = TRUE,
             show_colnames = TRUE,
             fontsize = 8,
             fontsize_row = 7,
             fontsize_col = 7,
             main = "Sample Correlation Heatmap (Spearman)"))
  },
  filename = "Correlation_heatmap_spearman",
  subdir = "02_PCA_Clustering",
  width = 10, height = 9
)

# ============================================================================
# 5. Sample layer times Clustering
# ============================================================================
log_step("03_NORM", "Hierarchical clustering...")

# Based on1-Spearman 
dist_mat <- as.dist(1 - cor_mat)
hc <- hclust(dist_mat, method = "complete")

# using dendextend
dend <- as.dendrogram(hc)
labels_order <- labels(dend)
label_colors <- ifelse(sample_info[labels_order, "Treatment"] == "Induced",
                       COLORS$treatment["Induced"],
                       COLORS$treatment["Control"])
dend <- dend %>%
  set("labels_cex", 0.7) %>%
  set("labels_col", label_colors)

save_heatmap_fig(
  draw_func = function() {
    par(mar = c(8, 4, 3, 1))
    plot(dend, main = "Sample Hierarchical Clustering",
         ylab = "1 - Spearman Correlation", xlab = "")
    legend("topright", legend = c(COLORS$treatment_labels["Induced"], 
                                  COLORS$treatment_labels["Control"]),
           col = c(COLORS$treatment["Induced"], COLORS$treatment["Control"]),
           pch = 15, cex = 0.8, bty = "n")
  },
  filename = "Sample_clustering_dendrogram",
  subdir = "02_PCA_Clustering",
  width = 12, height = 6
)

# ============================================================================
# 6. SampleFigure（Figure）
# ============================================================================
sample_dist <- dist(t(vst_mat))
sample_dist_mat <- as.matrix(sample_dist)

save_heatmap_fig(
  draw_func = function() {
    print(pheatmap(sample_dist_mat,
             clustering_distance_rows = sample_dist,
             clustering_distance_cols = sample_dist,
             color = colorRampPalette(c("#00A087", "white", "#E64B35"))(100),
             annotation_col = anno_col,
             annotation_colors = anno_colors,
             show_rownames = TRUE,
             show_colnames = FALSE,
             fontsize = 7,
             main = "Sample Distance Heatmap (Euclidean on VST)"))
  },
  filename = "Sample_distance_heatmap",
  subdir = "02_PCA_Clustering",
  width = 10, height = 9
)

log_step("03_NORM", "Step 03 COMPLETE")
