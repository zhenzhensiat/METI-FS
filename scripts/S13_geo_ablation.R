# T4: GEO Ablation A3-A9 — reuse existing intermediate files
source("S_config.R")
source("S06_ablation_study.R")

geo_dirs <- c(
  file.path(RUN_DIR, "GEO_GSE197067_Tcell"),
  file.path(RUN_DIR, "GEO_GSE307424_Lung"),
  file.path(RUN_DIR, "GEO_GSE236646_NPC"),
  file.path(RUN_DIR, "GEO_GSE150411_Chon")
)

# Run A3,A4,A5,A6,A9 (skip A0,A1,A2 already done in simulation)
configs <- c("A0_FULL", "A3_no_EffectSize", "A4_no_GapUnion",
             "A5_no_PPI", "A6_no_Bootstrap", "A9_RF_weighted")

for (gd in geo_dirs) {
  rid <- basename(gd)
  cat(sprintf("\n=== %s ===\n", rid))
  tryCatch({
    res <- run_ablation(gd, run_id = rid, configs = configs)
    cat(sprintf("[OK] %s\n", rid))
  }, error = function(e) {
    cat(sprintf("[ERROR] %s: %s\n", rid, e$message))
  })
}

cat("\n=== T4 DONE ===\n")
cat("Output: benchmark_results/ablation_GEO_*.csv\n")
