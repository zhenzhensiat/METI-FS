#!/usr/bin/env Rscript
# ==============================================================================
# 06_maSigPro_trends.R — columns
# using maSigPro-GLMRNA-seq countscolumnsRegression analysis
# ：maSigPro using item GLMcounts
# ==============================================================================

source(file.path(SCRIPT_DIR, "00_setup.R"))
source(file.path(SCRIPT_DIR, "theme_bindlab.R"))

log_step("06_MASIGPRO", "Starting maSigPro time series analysis...")

# Load data
counts_filtered <- readRDS(FILES$counts_filtered)
sample_info     <- readRDS(FILES$sample_info)
gene_anno       <- readRDS(FILES$gene_annotation)

# ============================================================================
# 1. maSigProMatrix
# ============================================================================

# maSigPro need need to edesign:
# =Sample, columns=Time, Replicate, Group1, Group2, ...
# Groupcolumns is (0/1)

edesign <- data.frame(
  Time      = sample_info$time_num,
  Replicate = sample_info$replicate,
  Control   = as.integer(sample_info$Treatment == "Control"),
  Induced   = as.integer(sample_info$Treatment == "Induced"),
  row.names = rownames(sample_info)
)

# Sample and counts
edesign <- edesign[colnames(counts_filtered), ]

log_step("06_MASIGPRO", "Experimental design matrix:")
print(head(edesign))
log_step("06_MASIGPRO", sprintf("Design: %d samples, Time points: %s",
                                 nrow(edesign),
                                 paste(sort(unique(edesign$Time)), collapse = ", ")))

# ============================================================================
# 2. maSigPro
# ============================================================================

# Matrix: genes in rows, samples in columns
data_mat <- as.matrix(counts_filtered)

# step : p.vector — using item GLMFilterprofiles
# degree = 3 (4 Timepoint， using times item )
# counts = TRUE using item 
log_step("06_MASIGPRO", "Step 1: Identifying genes with non-flat profiles (p.vector)...")

design_matrix <- make.design.matrix(edesign, degree = PARAMS$masigpro_degree)

fit <- p.vector(data_mat, design_matrix,
                Q = PARAMS$masigpro_alfa,
                counts = TRUE,          # Use negative binomial GLM
                min.obs = 20)           # Minimum observations

log_step("06_MASIGPRO", sprintf("Step 1 result: %d genes with significant profiles (FDR < %.2f)",
                                 fit$i, PARAMS$masigpro_alfa))

# step : T.fit — step Regression， to Gene 
log_step("06_MASIGPRO", "Step 2: Stepwise regression (T.fit)...")

tstep <- T.fit(fit, step.method = "backward",
               alfa = PARAMS$masigpro_alfa)

# step : FilterR² Gene
log_step("06_MASIGPRO", sprintf("Step 3: Filtering by R² >= %.2f...", PARAMS$masigpro_rsq))

sigs <- get.siggenes(tstep, 
                     rsq = PARAMS$masigpro_rsq,
                     vars = "groups")

# group SignificantGene
sig_genes_all <- sigs$summary

log_step("06_MASIGPRO", "Significant genes per group:")
for (grp in names(sig_genes_all)) {
  log_step("06_MASIGPRO", sprintf("  %s: %d genes", grp, length(sig_genes_all[[grp]])))
}
log_step("06_MASIGPRO", "  (Only InducedvsControl used for clustering and downstream integration)")

# ============================================================================
# 3. Trend clustering
# ============================================================================
log_step("06_MASIGPRO", sprintf("Clustering into %d groups...", PARAMS$masigpro_k))

# using see.genesCluster visualization
# Induced group SignificantGeneClustering
cluster_result <- NULL

tryCatch({
 # using hclustClustering
 # ★ Check: sigs$sig.genes edesigncolumns
  sig_genes_for_cluster <- NULL
  if (!is.null(sigs$sig.genes)) {
    avail_keys <- names(sigs$sig.genes)
    log_step("06_MASIGPRO", sprintf("Available sig.genes keys: %s", paste(avail_keys, collapse = ", ")))
 # InducedvsControl, fallback to can using key
    if ("InducedvsControl" %in% avail_keys) {
      sig_genes_for_cluster <- sigs$sig.genes$InducedvsControl
    } else if (length(avail_keys) > 0) {
      sig_genes_for_cluster <- sigs$sig.genes[[avail_keys[1]]]
      log_step("06_MASIGPRO", sprintf("Using fallback key: %s", avail_keys[1]))
    }
  }
  
  if (!is.null(sig_genes_for_cluster) && 
      (is.list(sig_genes_for_cluster) || nrow(sig_genes_for_cluster) >= 10)) {
 # ★ : see.genes need need to edesign and groups.vector
 # sig.genes$InducedvsControl is list (sig.profiles/coefficients/t.score etc.)
 # see.genes Defaultedesign=data$edesign，sig.geneslist not edesign
 # from Rawedesign， then within time/repvect/groupsFailed
 # (" need need to TRUE/FALSE not can to using ")
    #   Reference: maSigPro User's Guide (Conesa & Nueda, 2006, Bioinformatics 22:1096)
    log_step("06_MASIGPRO", sprintf("Clustering significant genes into %d groups...", 
                                     PARAMS$masigpro_k))
    pdf(file.path(FIG_DIR, "06_maSigPro", "maSigPro_cluster_profiles.pdf"),
        width = 14, height = 10)
    cluster_result <- see.genes(sig_genes_for_cluster,
                               edesign = edesign,
                               groups.vector = design_matrix$groups.vector,
                               show.fit = TRUE,
                               dis = design_matrix$dis,
                               cluster.method = "hclust",
                               cluster.data = 1,  # Based on regression coefficient clustering
                               k = PARAMS$masigpro_k,
                               newX11 = FALSE)
  dev.off()
  
 # SavePNG
  png(file.path(FIG_DIR, "06_maSigPro", "maSigPro_cluster_profiles.png"),
      width = 14, height = 10, units = "in", res = 300)
  see.genes(sig_genes_for_cluster,
            edesign = edesign,
            groups.vector = design_matrix$groups.vector,
            show.fit = TRUE,
            dis = design_matrix$dis,
            cluster.method = "hclust",
            cluster.data = 1,
            k = PARAMS$masigpro_k,
            newX11 = FALSE)
  dev.off()
  
  log_step("06_MASIGPRO", "Cluster profiles plotted successfully")
  } else {
    log_step("06_MASIGPRO", "WARNING: No suitable sig.genes for clustering (NULL or <10 genes)")
  }
}, error = function(e) {
  log_step("06_MASIGPRO", sprintf("WARNING: Cluster visualization error: %s", e$message))
})

# ============================================================================
# 4. ClusteringInfoSave
# ============================================================================

if (!is.null(cluster_result)) {
 # Gene-cluster should 
  gene_clusters <- data.frame(
    ensembl_id = names(cluster_result$cut),
    cluster    = cluster_result$cut,
    stringsAsFactors = FALSE
  )
  gene_clusters$symbol <- gene_anno$hgnc_symbol[
    match(gene_clusters$ensembl_id, gene_anno$ensembl_gene_id)]
  
  log_step("06_MASIGPRO", "Genes per cluster:")
  print(table(gene_clusters$cluster))
  
  write.csv(gene_clusters, 
            file.path(DATA_DIR, "maSigPro_gene_clusters.csv"),
            row.names = FALSE)
}

# ============================================================================
# 5. CustomVisualization（Nature）
# ============================================================================
log_step("06_MASIGPRO", "Generating custom trend plots...")

if (!is.null(cluster_result)) {
 # using TPM（ ）
  tpm_filtered <- readRDS(FILES$tpm_filtered)
  
  for (cl in sort(unique(gene_clusters$cluster))) {
    cl_genes <- gene_clusters$ensembl_id[gene_clusters$cluster == cl]
    
    if (length(cl_genes) < 2) next
    
 # Gene TPM，log2
    tpm_cl <- tpm_filtered[cl_genes, ]
    tpm_cl_log <- log2(tpm_cl + 1)
    
 # Calculate and SE
    plot_data <- tpm_cl_log %>%
      as.data.frame() %>%
      rownames_to_column("gene") %>%
      pivot_longer(-gene, names_to = "sample", values_to = "expr") %>%
      left_join(sample_info %>% rownames_to_column("sample") %>%
                  dplyr::select(sample, Treatment, Time, time_num), by = "sample") %>%
      group_by(Treatment, Time, time_num) %>%
      summarise(
        mean_expr = mean(expr),
        se_expr   = sd(expr) / sqrt(n()),
        .groups = "drop"
      )
    
    p_trend <- ggplot(plot_data, aes(x = time_num, y = mean_expr, 
                                     color = Treatment, group = Treatment)) +
      geom_line(linewidth = 1) +
      geom_point(size = 2.5) +
      geom_errorbar(aes(ymin = mean_expr - se_expr, ymax = mean_expr + se_expr),
                    width = 0.8, linewidth = 0.5) +
      scale_color_treatment() +
      scale_x_continuous(breaks = c(4, 7, 14, 21), labels = PARAMS$time_labels) +
      labs(title = sprintf("Cluster %d (%d genes)", cl, length(cl_genes)),
           x = "Time (days)", y = "Mean log2(TPM + 1)") +
      theme_bindlab()
    
    save_pub_fig(p_trend, sprintf("maSigPro_Cluster%d_trend", cl), 
                 "06_maSigPro", width = 6, height = 4.5)
  }
  
 # group has clusterFigure
  log_step("06_MASIGPRO", "All cluster trend plots saved")
}

# ============================================================================
# 6. Save results
# ============================================================================
masigpro_results <- list(
  fit             = fit,
  tstep           = tstep,
  sigs            = sigs,
  sig_genes_all   = sig_genes_all,
  cluster_result  = cluster_result,
  gene_clusters   = if(exists("gene_clusters")) gene_clusters else NULL,
  design_matrix   = design_matrix
)

save_data(masigpro_results, FILES$masigpro_results)

log_step("06_MASIGPRO", "Step 06 COMPLETE")
