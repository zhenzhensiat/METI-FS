#!/usr/bin/env Rscript
# =============================================================================
# 09F_PPI_hub_selection.R — PPI hub parallel selection (LOCAL STRING database)
#
# Role: PARALLEL selection line alongside 09D ML gap-union.
#       Literature: ~85% of ML+PPI biomarker papers use parallel union
#       (Session 5 literature survey, 41 papers, 2026-03-23).
#
# Method:
#   1. Input = candidate pool from 09A (same as ML input)
#   2. Build PPI network from LOCAL STRING database files
#      (stringdb_cache/9606.protein.links.v12.0.txt.gz)
#   3. Map STRING protein IDs to gene symbols via protein.info
#   4. Compute 5 centrality metrics in igraph
#   5. Cascading consensus hub selection (5/5 -> 4/5 -> 3/5)
#
# Centrality metrics (5):
#   Degree, Betweenness, PageRank, Eigenvector, Hub score (HITS)
#   Note: Closeness centrality removed — undefined for disconnected graphs
#   (Boldi & Vigna 2014, Wasserman & Faust 1994)
#
# STRING data source:
#   Local files in DATA_DIR/stringdb_cache/ (downloaded once, fully reproducible)
#   - 9606.protein.links.v12.0.txt.gz   (interaction edges + combined_score)
#   - 9606.protein.info.v12.0.txt.gz    (protein_id -> preferred_name mapping)
#
# Requires: source("00_setup.R") already executed
#
# Input:  DATA_DIR/candidate_pool.rds  (09A)
# Output: DATA_DIR/ppi_hub_selection.rds
#         DATA_DIR/PPI_09F_hub_genes.csv
#         DATA_DIR/PPI_09F_all_metrics.csv
#         DATA_DIR/PPI_09F_selection_summary.txt
#         FIG_DIR/09F_PPI_Hub/*.pdf|png
# =============================================================================

STEP_NAME  <- "09F_PPI_Hub"
fig_subdir <- "09F_PPI_Hub"
dir.create(file.path(FIG_DIR, fig_subdir), showWarnings = FALSE, recursive = TRUE)

# ---------------------------------------------------------------------------
# 1. Upstream check
# ---------------------------------------------------------------------------

pool_file <- file.path(DATA_DIR, "candidate_pool.rds")

check_upstream(STEP_NAME,
  upstream_files = c("09A_pool" = pool_file),
  output_files = c(
    file.path(DATA_DIR, "ppi_hub_selection.rds"),
    file.path(DATA_DIR, "PPI_09F_hub_genes.csv")
  )
)

log_step(STEP_NAME, "Starting PPI hub selection (parallel line, LOCAL STRING database)")

# ---------------------------------------------------------------------------
# 2. Load candidate pool
# ---------------------------------------------------------------------------

pool_result    <- readRDS(pool_file)
candidate_pool <- pool_result$candidate_pool
gene_anno      <- readRDS(FILES$gene_annotation)

# Map ensembl -> symbol
ppi_genes <- unique(gene_anno$hgnc_symbol[
  match(candidate_pool, gene_anno$ensembl_gene_id)])
ppi_genes <- ppi_genes[!is.na(ppi_genes) & ppi_genes != ""]

log_step(STEP_NAME, sprintf("Candidate pool: %d ensembl IDs -> %d gene symbols",
  length(candidate_pool), length(ppi_genes)))

# ---------------------------------------------------------------------------
# 3. Parameters
# ---------------------------------------------------------------------------

STRING_SPECIES <- PARAMS$string_species  # 9606
STRING_SCORE   <- PARAMS$string_score    # 700

# Hub selection: 5 algorithms, cascading consensus
# Closeness removed: Boldi & Vigna 2014 (undefined for disconnected graphs)
N_ALGO         <- 5
TOP_N          <- 10
MIN_ALGO_FLOOR <- 3   # absolute minimum consensus

log_step(STEP_NAME, sprintf("STRING: species=%d, score>=%d (LOCAL database)",
  STRING_SPECIES, STRING_SCORE))
log_step(STEP_NAME, sprintf("Hub params: %d algos, top%d, cascading consensus 5->4->3/%d",
  N_ALGO, TOP_N, N_ALGO))

# ---------------------------------------------------------------------------
# 4. Load LOCAL STRING database
# ---------------------------------------------------------------------------

string_cache_dir <- file.path(DATA_DIR, "stringdb_cache")

links_file <- file.path(string_cache_dir,
  sprintf("%d.protein.links.v12.0.txt.gz", STRING_SPECIES))
info_file  <- file.path(string_cache_dir,
  sprintf("%d.protein.info.v12.0.txt.gz", STRING_SPECIES))

# Check files exist
if (!file.exists(links_file)) {
  stop(sprintf("[%s] STRING links file not found: %s\nDownload from https://string-db.org/cgi/download",
    STEP_NAME, links_file))
}
if (!file.exists(info_file)) {
  stop(sprintf("[%s] STRING info file not found: %s\nDownload from https://string-db.org/cgi/download",
    STEP_NAME, info_file))
}

# --- Load protein info for ID mapping ---
log_step(STEP_NAME, sprintf("Loading STRING protein info: %s", basename(info_file)))
protein_info <- read.delim(gzfile(info_file), stringsAsFactors = FALSE)
log_step(STEP_NAME, sprintf("  Loaded %d protein entries", nrow(protein_info)))

# Identify the protein ID column (naming varies across STRING versions)
id_col_name <- intersect(
  c("protein_external_id", "X.string_protein_id", "string_protein_id"),
  colnames(protein_info))
if (length(id_col_name) == 0) {
  id_col_name <- colnames(protein_info)[1]
  log_step(STEP_NAME, sprintf("  Using first column as protein ID: %s", id_col_name))
} else {
  id_col_name <- id_col_name[1]
}

# Build symbol <-> STRING ID mapping (only for candidate genes)
symbol_to_string <- protein_info[protein_info$preferred_name %in% ppi_genes,
                                  c(id_col_name, "preferred_name")]
colnames(symbol_to_string) <- c("string_id", "symbol")
symbol_to_string <- unique(symbol_to_string)

n_mapped <- length(unique(symbol_to_string$symbol))
log_step(STEP_NAME, sprintf("  Mapped %d/%d candidate symbols to STRING protein IDs",
  n_mapped, length(ppi_genes)))

if (n_mapped < 10) {
  log_step(STEP_NAME, "WARNING: Very few genes mapped. Check species ID and gene symbol format.")
}

# --- Load interaction links ---
log_step(STEP_NAME, sprintf("Loading STRING links: %s (this may take a moment...)",
  basename(links_file)))

all_links <- read.delim(gzfile(links_file), sep = " ", stringsAsFactors = FALSE)
log_step(STEP_NAME, sprintf("  Loaded %d total interactions (genome-wide)", nrow(all_links)))

# Filter: both proteins in candidate set AND score >= threshold
candidate_string_ids <- unique(symbol_to_string$string_id)

edge_df_raw <- all_links[
  all_links$protein1 %in% candidate_string_ids &
  all_links$protein2 %in% candidate_string_ids &
  all_links$combined_score >= STRING_SCORE,
]

log_step(STEP_NAME, sprintf("  Filtered: %d edges (both in candidate pool, score>=%d)",
  nrow(edge_df_raw), STRING_SCORE))

# Clean up large object
rm(all_links); gc(verbose = FALSE)

if (nrow(edge_df_raw) == 0) {
  log_step(STEP_NAME, "ERROR: No PPI interactions found above score threshold.")
  saveRDS(list(hub_genes = data.frame(), all_metrics = data.frame(),
               n_nodes = 0, n_edges = 0, method = "no_interactions"),
          file.path(DATA_DIR, "ppi_hub_selection.rds"))
  write.csv(data.frame(), file.path(DATA_DIR, "PPI_09F_hub_genes.csv"), row.names = FALSE)
  log_step(STEP_NAME, "Saved empty result. Step 09F COMPLETE (no network)")
} else {

# --- Map STRING protein IDs back to gene symbols ---
string_to_symbol <- setNames(symbol_to_string$symbol, symbol_to_string$string_id)

edge_df <- data.frame(
  gene_A = string_to_symbol[edge_df_raw$protein1],
  gene_B = string_to_symbol[edge_df_raw$protein2],
  score  = edge_df_raw$combined_score,
  stringsAsFactors = FALSE
)

# Remove unmapped / self-loops / duplicates
edge_df <- edge_df[!is.na(edge_df$gene_A) & !is.na(edge_df$gene_B), ]
edge_df <- edge_df[edge_df$gene_A != edge_df$gene_B, ]
edge_df <- unique(edge_df)

log_step(STEP_NAME, sprintf("PPI network: %d unique edges (gene-symbol level)", nrow(edge_df)))

# ---------------------------------------------------------------------------
# 5. Build igraph + compute centrality metrics
# ---------------------------------------------------------------------------

library(igraph)

g <- graph_from_data_frame(edge_df[, c("gene_A", "gene_B")], directed = FALSE)
g <- simplify(g)

log_step(STEP_NAME, sprintf("Graph: %d nodes, %d edges", vcount(g), ecount(g)))

# Network connectivity diagnostic
n_components <- components(g)$no
largest_comp <- max(components(g)$csize)
log_step(STEP_NAME, sprintf("Graph components: %d (largest: %d nodes, %.1f%%)",
  n_components, largest_comp, 100 * largest_comp / vcount(g)))
if (n_components > 1) {
  log_step(STEP_NAME, "NOTE: Network is disconnected — Closeness centrality excluded (Boldi & Vigna 2014)")
}

# Compute 5 centrality metrics (Closeness removed — undefined for disconnected graphs)
node_metrics <- data.frame(
  symbol       = V(g)$name,
  degree       = degree(g),
  betweenness  = betweenness(g, normalized = TRUE),
  pagerank     = page_rank(g)$vector,
  eigenvector  = tryCatch(eigen_centrality(g)$vector,
                           error = function(e) rep(NA_real_, vcount(g))),
  hub_score    = tryCatch(hub_score(g)$vector,
                           error = function(e) rep(NA_real_, vcount(g))),
  stringsAsFactors = FALSE
)

# Rank per metric (1 = best)
algo_cols <- c("degree", "betweenness", "pagerank", "eigenvector", "hub_score")
algo_labels <- c("Degree", "Betweenness", "PageRank", "Eigenvector", "Hub(HITS)")

for (ac in algo_cols) {
  rank_col <- paste0("rank_", ac)
  node_metrics[[rank_col]] <- rank(-node_metrics[[ac]], ties.method = "min", na.last = "keep")
}

# Mean rank
rank_cols <- paste0("rank_", algo_cols)
node_metrics$mean_rank <- rowMeans(node_metrics[, rank_cols], na.rm = TRUE)
node_metrics <- node_metrics[order(node_metrics$mean_rank), ]
rownames(node_metrics) <- NULL

# ---------------------------------------------------------------------------
# 6. Hub selection: cascading consensus (5/5 -> 4/5 -> 3/5)
# ---------------------------------------------------------------------------

log_step(STEP_NAME, sprintf("=== Hub selection: top%d x %d algos, cascading consensus 5->4->3/%d ===",
  TOP_N, N_ALGO, N_ALGO))

# Count how many algorithms put each gene in top-N
in_topN <- sapply(rank_cols, function(rc) {
  node_metrics[[rc]] <= TOP_N
})
node_metrics$n_topN <- rowSums(in_topN, na.rm = TRUE)

# Cascading: try 5/5 first, then 4/5, then 3/5
hub_df <- data.frame()
actual_min <- N_ALGO

for (try_min in seq(N_ALGO, MIN_ALGO_FLOOR, by = -1)) {
  candidates <- node_metrics[node_metrics$n_topN >= try_min, ]
  if (nrow(candidates) > 0) {
    hub_df <- candidates
    actual_min <- try_min
    log_step(STEP_NAME, sprintf("  Consensus %d/%d: %d genes found",
      try_min, N_ALGO, nrow(candidates)))
    break
  } else {
    log_step(STEP_NAME, sprintf("  Consensus %d/%d: 0 genes, relaxing...",
      try_min, N_ALGO))
  }
}

hub_df <- hub_df[order(hub_df$mean_rank), ]

log_step(STEP_NAME, sprintf("RESULT: %d PPI hub genes (consensus >= %d/%d)",
  nrow(hub_df), actual_min, N_ALGO))
if (nrow(hub_df) > 0) {
  log_step(STEP_NAME, sprintf("Genes: %s", paste(hub_df$symbol, collapse = ", ")))
}

# Map back to ensembl_id
hub_df$ensembl_id <- gene_anno$ensembl_gene_id[
  match(hub_df$symbol, gene_anno$hgnc_symbol)]

# ---------------------------------------------------------------------------
# 7. Save outputs
# ---------------------------------------------------------------------------

log_step(STEP_NAME, "=== Saving ===")

output_09f <- list(
  hub_genes   = hub_df,
  hub_gene_ids = hub_df$ensembl_id[!is.na(hub_df$ensembl_id)],
  all_metrics = node_metrics,
  n_nodes     = vcount(g),
  n_edges     = ecount(g),
  params      = list(
    STRING_SCORE       = STRING_SCORE,
    STRING_VERSION     = "v12.0",
    STRING_SOURCE      = "local",
    N_ALGO             = N_ALGO,
    TOP_N              = TOP_N,
    MIN_ALGO_FLOOR     = MIN_ALGO_FLOOR,
    actual_min         = actual_min,
    cascading_strategy = sprintf("%d->%d->%d/%d", N_ALGO, N_ALGO-1, MIN_ALGO_FLOOR, N_ALGO),
    algo_labels        = algo_labels
  ),
  method = "5algo_topN_consensus_local_STRING",
  literature = c(
    "Session 5 survey: ML+PPI parallel union in ~85% of biomarker papers",
    "Closeness removed: Boldi & Vigna 2014 Internet Math 10:222-262 (undefined for disconnected graphs)",
    "Wasserman & Faust 1994: closeness only defined for connected graphs",
    "Chin et al. 2014 BMC Syst Biol (cytoHubba): multi-algorithm consensus hub selection"
  )
)

saveRDS(output_09f, file.path(DATA_DIR, "ppi_hub_selection.rds"))
log_step(STEP_NAME, "Saved RDS -> ppi_hub_selection.rds")

write.csv(hub_df, file.path(DATA_DIR, "PPI_09F_hub_genes.csv"), row.names = FALSE)
log_step(STEP_NAME, "Saved hub genes CSV")

write.csv(node_metrics, file.path(DATA_DIR, "PPI_09F_all_metrics.csv"), row.names = FALSE)
log_step(STEP_NAME, "Saved all metrics CSV")

# Summary report
summary_file <- file.path(DATA_DIR, "PPI_09F_selection_summary.txt")
sink(summary_file)
cat("=============================================================\n")
cat(sprintf("09F PPI Hub Selection -- %s\n", PARAMS$diff_type))
cat(sprintf("Date: %s\n", Sys.time()))
cat("=============================================================\n\n")
cat("METHOD: 5-algorithm cascading consensus hub selection (parallel to ML)\n")
cat(sprintf("  STRING source: LOCAL database (v12.0)\n"))
cat("  Closeness centrality removed: undefined for disconnected graphs\n")
cat("  (Boldi & Vigna 2014; Wasserman & Faust 1994)\n")
cat(sprintf("  STRING: species=%d, score>=%d\n", STRING_SPECIES, STRING_SCORE))
cat(sprintf("  Algorithms: %s\n", paste(algo_labels, collapse = ", ")))
cat("  Cascading consensus: prefer 5/5, fallback 4/5, minimum 3/5\n")
cat(sprintf("  Achieved consensus: >= %d/%d\n\n", actual_min, N_ALGO))

cat(sprintf("NETWORK: %d nodes, %d edges\n", vcount(g), ecount(g)))
cat(sprintf("  Components: %d (largest: %d nodes, %.1f%%)\n",
  n_components, largest_comp, 100 * largest_comp / vcount(g)))
cat(sprintf("  Candidate symbols mapped to STRING: %d/%d\n", n_mapped, length(ppi_genes)))
cat(sprintf("HUB GENES: %d\n\n", nrow(hub_df)))

if (nrow(hub_df) > 0) {
  for (i in seq_len(nrow(hub_df))) {
    r <- hub_df[i, ]
    cat(sprintf("  %d. %-12s  degree=%d  n_topN=%d/%d  mean_rank=%.1f\n",
                i, r$symbol, r$degree, r$n_topN, N_ALGO, r$mean_rank))
  }
}

cat("\n\nTOP-10 PER ALGORITHM:\n")
cat("-------------------------------------------------------------\n")
for (j in seq_along(rank_cols)) {
  rc <- rank_cols[j]
  top10 <- node_metrics[node_metrics[[rc]] <= 10, ]
  top10 <- top10[order(top10[[rc]]), ]
  cat(sprintf("\n[%s] top 10:\n", algo_labels[j]))
  for (k in seq_len(min(10, nrow(top10)))) {
    cat(sprintf("  %2d. %s (rank=%d, value=%.4f)\n",
                k, top10$symbol[k], top10[[rc]][k], top10[[algo_cols[j]]][k]))
  }
}
sink()
log_step(STEP_NAME, "Saved summary -> PPI_09F_selection_summary.txt")

# ---------------------------------------------------------------------------
# 8. Diagnostic plots
# ---------------------------------------------------------------------------

log_step(STEP_NAME, "=== Diagnostic plots ===")

suppressPackageStartupMessages(library(ggplot2))

# 8a. Hub degree barplot
if (nrow(hub_df) > 0) {
  plot_df <- hub_df[order(-hub_df$degree), ]
  plot_df$symbol <- factor(plot_df$symbol, levels = rev(plot_df$symbol))

  p_hub <- ggplot(plot_df, aes(x = symbol, y = degree)) +
    geom_col(fill = "#E64B35", width = 0.6) +
    geom_text(aes(label = sprintf("n=%d/%d", n_topN, N_ALGO)), hjust = -0.1, size = 3) +
    coord_flip() +
    labs(title = sprintf("PPI Hub Genes (%d, consensus >= %d/%d)", nrow(hub_df), actual_min, N_ALGO),
         subtitle = sprintf("STRING v12.0 (local), score >= %d, top-%d per algorithm", STRING_SCORE, TOP_N),
         x = NULL, y = "Degree centrality") +
    theme_minimal(base_size = 12) +
    theme(legend.position = "none",
          panel.grid.major.y = element_blank()) +
    expand_limits(y = max(plot_df$degree) * 1.15)

  save_pub_fig(p_hub, "09F_hub_degree_barplot", fig_subdir, width = 8, height = max(4, nrow(hub_df) * 0.4 + 2))
}

# 8b. Consensus distribution
consensus_dist <- data.frame(table(node_metrics$n_topN))
colnames(consensus_dist) <- c("n_topN", "count")
consensus_dist$n_topN <- as.integer(as.character(consensus_dist$n_topN))
consensus_dist <- consensus_dist[consensus_dist$n_topN > 0, ]

if (nrow(consensus_dist) > 0) {
  p_dist <- ggplot(consensus_dist, aes(x = factor(n_topN), y = count)) +
    geom_col(fill = "#4DBBD5", width = 0.7) +
    geom_text(aes(label = count), vjust = -0.3, size = 3.5) +
    geom_vline(xintercept = which(levels(factor(consensus_dist$n_topN)) == as.character(actual_min)) - 0.5,
               linetype = "dashed", color = "#E64B35", linewidth = 0.6) +
    labs(title = "Distribution of algorithm consensus",
         subtitle = sprintf("How many algorithms' top-%d list each gene appears in", TOP_N),
         x = sprintf("Number of algorithms (out of %d)", N_ALGO),
         y = "Number of genes") +
    theme_minimal(base_size = 12)

  save_pub_fig(p_dist, "09F_consensus_distribution", fig_subdir, width = 8, height = 5)
}

# 8c. Rank heatmap for top 20
top20_metrics <- head(node_metrics, 20)
if (nrow(top20_metrics) > 0) {
  tryCatch({
    rank_mat <- as.matrix(top20_metrics[, rank_cols])
    rownames(rank_mat) <- top20_metrics$symbol
    colnames(rank_mat) <- algo_labels

    if (requireNamespace("pheatmap", quietly = TRUE)) {
      save_heatmap_fig(
        draw_func = function() {
          print(pheatmap::pheatmap(rank_mat,
            color = colorRampPalette(c("#E64B35", "white", "#3C5488"))(100),
            cluster_cols = FALSE,
            fontsize = 9,
            display_numbers = TRUE,
            number_format = "%.0f",
            main = sprintf("Hub Gene Centrality Rankings (Top 20) — %s", PARAMS$diff_type)))
        },
        filename = "09F_hub_centrality_heatmap",
        subdir = fig_subdir,
        width = 8, height = 8
      )
    }
  }, error = function(e) {
    log_step(STEP_NAME, sprintf("Heatmap error: %s", e$message))
  })
}

# ---------------------------------------------------------------------------
# 9. Done
# ---------------------------------------------------------------------------

log_step(STEP_NAME, "========================================")
log_step(STEP_NAME, sprintf("COMPLETE: %d PPI hub genes", nrow(hub_df)))
log_step(STEP_NAME, sprintf("  Network: %d nodes, %d edges (LOCAL STRING v12.0)", vcount(g), ecount(g)))
log_step(STEP_NAME, sprintf("  Consensus: >= %d/%d algos in top-%d", actual_min, N_ALGO, TOP_N))
if (nrow(hub_df) > 0) {
  log_step(STEP_NAME, sprintf("  Genes: %s", paste(hub_df$symbol, collapse = ", ")))
}
log_step(STEP_NAME, "========================================")

} # end if (edge_df_raw has rows)
