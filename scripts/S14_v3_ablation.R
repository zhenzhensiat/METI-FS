# S14_v3_ablation.R — Extended structure ablation (A1, A2) with JT trend test
# Reads existing intermediate files for fast re-evaluation across structures.
source("S_config.R")
source("S06_ablation_study.R")

v3_patterns <- c("block4_rho0.2","block4_rho0.6","hierarchical_8","overlapping","tp6","tp8")

# Find all v3 pipeline runs with ground truth + gap union
sim_dirs <- list.dirs(file.path(RUN_DIR), recursive=FALSE, full.names=TRUE)
v3_dirs <- c()
for (sd in sim_dirs) {
  nm <- basename(sd)
  is_v3 <- any(sapply(v3_patterns, function(p) grepl(p, nm, fixed=TRUE)))
  if (!is_v3) next
  gt <- file.path(sd, "data", "ground_truth.rds")
  gp <- file.path(sd, "data", "ml_gap_union.rds")
  if (file.exists(gt) && file.exists(gp)) v3_dirs <- c(v3_dirs, sd)
}

cat(sprintf("Found %d v3 pipeline runs with ground truth + gap union\n", length(v3_dirs)))

# Run A0, A1, A2 ablation on all v3 scenarios
all_rows <- list()
n_ok <- 0; n_err <- 0

for (rd in v3_dirs) {
  nm <- basename(rd)
  res <- tryCatch({
    run_ablation(rd, run_id=nm, configs=c("A0_FULL","A1_no_maSigPro","A2_no_WGCNA"))
  }, error = function(e) { cat(sprintf("  [ERROR] %s: %s\n", nm, e$message)); NULL })
  if (is.null(res)) { n_err <- n_err+1; next }

  # Extract ablation summary
  df <- res$summary
  if (is.null(df) || nrow(df) < 2) { n_err <- n_err+1; next }

  # Parse structure info from name (strip SIM_ prefix first)
  nm_clean <- sub("^SIM_", "", nm)
  parts <- strsplit(nm_clean, "_")[[1]]
  snr <- parts[1]
  # Identify structure
  struct <- "unknown"
  for (p in v3_patterns) {
    if (grepl(p, nm, fixed=TRUE)) { struct <- p; break }
  }

  for (i in seq_len(nrow(df))) {
    cfg_row <- df[i,]
    all_rows[[length(all_rows)+1]] <- data.frame(
      run_id = nm, config = cfg_row$config, label = cfg_row$label,
      snr = snr, structure = struct,
      n_pool = cfg_row$n_pool, n_final = cfg_row$n_final,
      precision = cfg_row$precision, recall = cfg_row$recall,
      F1 = cfg_row$F1, FDR = cfg_row$FDR,
      precision_any = if("precision_any_signal" %in% colnames(cfg_row)) cfg_row$precision_any_signal else NA,
      stringsAsFactors = FALSE
    )
  }
  n_ok <- n_ok+1
  if (n_ok %% 10 == 0) cat(sprintf("  %d/%d done\n", n_ok, length(v3_dirs)))
}

cat(sprintf("\n%d OK, %d errors\n", n_ok, n_err))

ab_df <- do.call(rbind, all_rows)
write.csv(ab_df, file.path(BENCH_DIR, "ablation_v3_all.csv"), row.names=FALSE)

# ---- Compute ΔF1 per scenario (paired) ----
cat("\n=== ΔF1 by structure ===\n")
baseline <- ab_df[ab_df$config == "A0_FULL", c("run_id","F1","precision","recall")]

# Compute dF1 for each config and accumulate
all_delta <- list()
for (cfg in c("A1_no_maSigPro","A2_no_WGCNA")) {
  cfg_data <- ab_df[ab_df$config == cfg,]
  m <- merge(baseline, cfg_data, by="run_id", suffixes=c(".base",".abl"))
  m$dF1 <- m$F1.abl - m$F1.base
  m$dPrec <- m$precision.abl - m$precision.base
  m$dRec <- m$recall.abl - m$recall.base
  m$config_label <- if(cfg=="A1_no_maSigPro") "-maSigPro" else "-WGCNA"

  for (i in seq_len(nrow(m))) {
    for (p in v3_patterns) {
      if (grepl(p, m$run_id[i], fixed=TRUE)) { m$structure[i] <- p; break }
    }
  }

  lbl <- if(cfg=="A1_no_maSigPro") "-maSigPro" else "-WGCNA"
  cat(sprintf("\n--- %s ---\n", lbl))
  cat(sprintf("%-18s %5s %7s %7s\n", "Structure","n","mean_dF1","sd_dF1"))
  for (st in v3_patterns) {
    sub <- m[m$structure==st,]
    if (nrow(sub)==0) next
    cat(sprintf("%-18s %5d %+7.3f %7.3f\n", st, nrow(sub), mean(sub$dF1), sd(sub$dF1)))
  }
  cat(sprintf("%-18s %5d %+7.3f %7.3f\n", "ALL", nrow(m), mean(m$dF1), sd(m$dF1)))
  all_delta[[cfg]] <- m
}

# Merge both configs into one dataframe
m_all <- rbind(all_delta[["A1_no_maSigPro"]], all_delta[["A2_no_WGCNA"]])

# ---- JT trend test for WGCNA contribution ----
cat("\n=== Jonckheere-Terpstra trend test (WGCNA dF1 vs structure complexity) ===\n")
order_map <- c("block4_rho0.2"=1,"block4_rho0.6"=2,"hierarchical_8"=3,"overlapping"=4)
m_a2 <- all_delta[["A2_no_WGCNA"]]
m_a2 <- m_a2[m_a2$structure %in% names(order_map),]
if (nrow(m_a2) > 0) {
  m_a2$order <- order_map[m_a2$structure]
  jt <- tryCatch({
    cor.test(m_a2$dF1, m_a2$order, method="kendall")
  }, error = function(e) NULL)
  if (!is.null(jt)) {
    cat(sprintf("Kendall tau=%.3f, p=%.4f\n", jt$estimate, jt$p.value))
    cat(sprintf("CONCLUSION: %s\n", if(jt$p.value < 0.05) "Significant trend" else "Not significant"))
  }
}

# Also JT for A1 (-maSigPro)
cat("\n=== JT trend test (maSigPro dF1 vs structure complexity) ===\n")
m_a1 <- all_delta[["A1_no_maSigPro"]]
m_a1 <- m_a1[m_a1$structure %in% names(order_map),]
if (nrow(m_a1) > 0) {
  m_a1$order <- order_map[m_a1$structure]
  jt1 <- tryCatch({
    cor.test(m_a1$dF1, m_a1$order, method="kendall")
  }, error = function(e) NULL)
  if (!is.null(jt1)) {
    cat(sprintf("Kendall tau=%.3f, p=%.4f\n", jt1$estimate, jt1$p.value))
    cat(sprintf("CONCLUSION: %s\n", if(jt1$p.value < 0.05) "Significant trend" else "Not significant"))
  }
}

# Save BOTH configs
write.csv(m_all, file.path(BENCH_DIR, "ablation_v3_deltaF1.csv"), row.names=FALSE)
cat(sprintf("\n[SAVED] ablation_v3_all.csv, ablation_v3_deltaF1.csv\n"))
cat("[DONE] T5\n")
