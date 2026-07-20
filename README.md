# METI-FS: Multi-Evidence Temporal Integration for Feature Selection

A computational framework for discovering sparse, reproducible biomarkers from
time-series transcriptomics data. METI-FS integrates differential expression,
temporal trend analysis, co-expression networks, machine-learning stability
selection, and protein-protein interaction evidence into a unified pipeline.

## Overview

METI-FS addresses the "curse of dimensionality" in temporal omics by combining
multiple complementary evidence layers:

1. **Temporal filtering** — DESeq2 interaction LRT + maSigPro trend analysis
2. **Network context** — WGCNA co-expression modules identify biologically coherent gene sets
3. **Candidate pool** — Intersection of temporal and network evidence with effect-size filtering
4. **Stability selection** — 100× bootstrap with LASSO, Random Forest, and SVM-RFE
5. **Gap-union thresholding** — Data-driven frequency cutoffs per algorithm
6. **PPI integration** — STRING network hub genes complement ML-selected markers
7. **Final union** — ML ∪ PPI hub markers as the final biomarker panel

## Repository Structure

```
METI-FS/
├── R/                        # Core pipeline scripts (run sequentially)
│   ├── 00_setup.R            # Environment, packages, global parameters
│   ├── theme_bindlab.R       # Custom ggplot2 theme
│   ├── 01_data_import.R      # Count matrix + sample info import
│   ├── 02_preprocessing.R    # Low-expression gene filtering
│   ├── 03_normalization_QC.R # DESeq2 normalization + PCA/QC
│   ├── 04_DEG_analysis.R     # Differential expression (LRT + Wald)
│   ├── 06_maSigPro_trends.R  # Temporal trend clustering
│   ├── 08_WGCNA.R            # Weighted co-expression network
│   ├── 09A_candidate_pool.R  # Multi-evidence candidate pool
│   ├── 09C_ML_stability_selection.R  # Bootstrap ML feature selection
│   ├── 09D_gap_union_selection.R     # Gap-based threshold + union
│   ├── 09F_PPI_hub_selection.R       # STRING PPI hub identification
│   └── 10_integration.R      # Final ML ∪ PPI integration
│
├── scripts/                   # Simulation, benchmarking, and analysis scripts
│   ├── S_config.R                     # Shared configuration (paths, parameters)
│   ├── S00_GEO_search.R               # GEO dataset search methodology
│   ├── S01_simulation_engine.R        # Synthetic data generator (v3)
│   ├── S01b_benchmark_timing.R        # Runtime profiling
│   ├── S02_pipeline_adapter.R         # Simulation to pipeline format adapter
│   ├── S03_download_and_prepare_GEO.R # GEO dataset download and preparation
│   ├── S04_run_pipeline_wrapper.R     # GEO dataset pipeline runner
│   ├── S04_sim_runner.R               # Simulation pipeline runner
│   ├── S05_benchmark_collector.R      # Result collection and metrics
│   ├── S06_ablation_study.R           # Component ablation analysis (LOCO)
│   ├── S07_cross_dataset_summary.R    # Cross-dataset summary and figures
│   ├── S08_baseline_comparison.R      # Baseline method comparison
│   ├── S09_manuscript_figures.R       # Publication figure generation
│   ├── S10_module_preservation.R      # Split-half WGCNA module preservation
│   ├── S11_lfc_sensitivity.R          # lfcThreshold sensitivity analysis
│   ├── S12_simple_baselines.R         # Simple baseline comparisons
│   ├── S13_geo_ablation.R             # GEO dataset ablation (A3-A9)
│   ├── S14_v3_ablation.R              # Extended structure ablation
│   ├── S15_rf_exclusion.R             # RF exclusion status analysis
│   ├── S16_wgcna_enrichment.R         # WGCNA module GO enrichment
│   ├── S17_ppi_analysis.R             # STRING PPI hub comparison
│   └── run_demo.R                     # Quick-start demo script
│
├── README.md
├── LICENSE
├── CITATION.cff
└── .gitignore

```

## Quick Start

```r
# 1. Set your project directory (containing data_raw/)
PROJECT_DIR <- "/path/to/your/project"

# 2. Source the setup script
source("R/00_setup.R")

# 3. Run pipeline steps sequentially
source("R/01_data_import.R")
source("R/02_preprocessing.R")
source("R/03_normalization_QC.R")
source("R/04_DEG_analysis.R")
source("R/06_maSigPro_trends.R")
source("R/08_WGCNA.R")
source("R/09A_candidate_pool.R")
source("R/09C_ML_stability_selection.R")
source("R/09D_gap_union_selection.R")
source("R/09F_PPI_hub_selection.R")
source("R/10_integration.R")
```

See `scripts/run_demo.R` for a complete walkthrough.

## Input Data Format

METI-FS expects the following files under `{PROJECT_DIR}/data_raw/`:

| File                       | Description                                                    |
| -------------------------- | -------------------------------------------------------------- |
| `{prefix}_counts.csv`      | Raw count matrix (genes × samples)                             |
| `{prefix}_tpm.csv`         | TPM expression matrix (genes × samples)                        |
| `{prefix}_sample_info.csv` | Sample metadata with columns: `sample_id`, `Time`, `Treatment` |

Sample naming convention: `{Prefix}{TimeLabel}{Rep}` for induced, `{Prefix}{TimeLabel}C{Rep}` for control.

## Key Parameters

All parameters are defined in `R/00_setup.R` under the `PARAMS` list and can be
overridden before sourcing downstream scripts:

| Parameter        | Default            | Description                                 |
| ---------------- | ------------------ | ------------------------------------------- |
| `time_labels`    | `c("4d","7d",...)` | Timepoint factor levels (user-configurable) |
| `padj_cutoff`    | 0.05               | FDR threshold for DEG                       |
| `lfc_cutoff`     | 1.0                | log2FC threshold                            |
| `masigpro_k`     | 9                  | Number of maSigPro clusters                 |
| `ml_bootstrap_n` | 100                | Bootstrap iterations for ML                 |
| `string_score`   | 700                | STRING confidence threshold                 |

## Advanced: Random Forest Handling

METI-FS provides two strategies for handling Random Forest (RF) in the bootstrap
stability selection layer (Layer 5), selectable via the `rf_mode` parameter.

### Binary Exclusion (default: `rf_mode = "exclude"`)

RF selections are excluded from the algorithm union when both conditions are met:
(a) the RF Nogueira chance-corrected stability index falls below the minimum of
LASSO and SVM-RFE indices for that dataset, and (b) the RF index is below 0.5.

This is the recommended default for WGCNA-enriched candidate pools, where
correlated predictors systematically degrade RF permutation importance
reliability [Strobl et al. 2007; Nicodemus et al. 2010]. In our validation on
four external GEO datasets, RF stability was consistently below 0.5, triggering
exclusion in all cases.

### Weighted Ensemble (`rf_mode = "weighted"`)

As an alternative to binary exclusion, RF genes are retained proportionally to
their stability relative to LASSO and SVM-RFE:

```
weight_RF = Nogueira_RF / max(Nogueira_LASSO, Nogueira_SVM)
```

Genes with weight >= 0.5 are added to the final ML gene set. This option is
appropriate when:

- Co-expression enrichment is weak or the WGCNA layer is bypassed (`skip_wgcna = TRUE`)
- RF may capture nonlinear interactions missed by linear methods
- Exploratory analyses where sensitivity is prioritized over specificity

**Usage example:**

```r
# Option A: Set in pipeline parameters (before sourcing pipeline scripts)
PARAMS$rf_mode <- "weighted"

# Option B: Apply directly to existing bootstrap results
source("scripts/S06_ablation_study.R")
weighted_result <- apply_rf_weighted_ensemble(upstream)
```

See `scripts/run_rf_weighted_demo.R` for a complete walkthrough comparing both
strategies on simulated data. The implementation is in
`scripts/S06_ablation_study.R` (lines 465-540, function `apply_rf_weighted_ensemble`).

## Computational Runtime

Benchmark measurements on synthetic data (13,000 genes, 4 timepoints, NB counts
with DESeq2-style dispersion, 12 CPU cores, R 4.5.2). All timings with B=100
bootstrap iterations per algorithm (LASSO, SVM-RFE, RF).

| Config | n | Pipeline Runtime | Peak Memory | L5a (Bootstrap) % |
|--------|---|-----------------|-------------|--------------------|
| Small (3+2 reps x 4 tp) | 20 | 38.5 min | 919 MB | 84.9% |
| Medium (4+3 reps x 4 tp) | 28 | 43.1 min | 946 MB | 85.9% |
| Large (5+5 reps x 4 tp) | 40 | 49.8 min | 986 MB | 87.2% |

Per-step breakdown (medium, n=28): data import (0.7 sec), preprocessing (4.2 sec),
DESeq2 norm+QC (6.5 sec), L1 DE analysis (63.7 sec), L2 maSigPro (116.4 sec),
L3 WGCNA (173.8 sec), L4 candidate pool (0.5 sec), **L5a bootstrap stability
(2219.6 sec, 85.9%)**, L5b gap-union (~2 sec). PPI (L6) and final integration
add ~7 sec combined.

Full per-step detail and system specifications in
`benchmark_results/timing_benchmark_results.md` (generated by
`scripts/S01c_benchmark_timing_full.R`).

### Key findings

- **Bootstrap stability selection (L5a) dominates total runtime** (85-87%),
  3 algorithms x B=100 iterations x 80% subsampling. Scales roughly linearly
  with n and B. Per-iteration cost: ~7 sec per algorithm at n=28.
- **WGCNA (L3) is second most expensive** (~174 sec at n=28), scaling with gene
  count rather than sample size. Pre-filtering low-expression genes is essential.
- **DESeq2 (L1) and maSigPro (L2) are lightweight** (~1-2 min each).
- **Peak memory <1 GB** for p=13,000, n<=40. 8 GB RAM sufficient.
- **Halving B to 50** halves L5a time (~20 min total for n=28), adequate for
  exploratory analyses with coarser frequency resolution.
  combined for typical datasets).
- **Peak memory usage** occurs during DESeq2 normalization or WGCNA adjacency
  matrix construction. For datasets with p < 20,000 and n < 50, 8 GB RAM is
  sufficient. For p > 30,000 genes or large cohort studies, 16 GB recommended.
- **Reducing B to 50** halves the bootstrap runtime at the cost of coarser
  frequency resolution (minimum detectable difference: 0.02 vs 0.01).
- **For time-limited users**: the effect size filter (Layer 4, lfcThreshold)
  reduces the candidate pool size and thus the bootstrap runtime. Setting
  lfcThreshold to a higher value further decreases runtime.

Benchmark script: `scripts/S01c_benchmark_timing_simple.R`.

## Dependencies

- R >= 4.2.0
- Bioconductor: DESeq2, maSigPro, clusterProfiler, org.Hs.eg.db
- CRAN: WGCNA, glmnet, randomForest, e1071, igraph, ggplot2, pheatmap
- Optional: STRINGdb (for PPI analysis), stabm (for Nogueira stability)

## Citation

If you use METI-FS in your research, please cite:

> Zhang Zhen et al. (2026). METI-FS: Multi-Evidence Temporal
> Integration for Feature Selection in Time-Series Transcriptomics.
> *BMC Bioinformatics*.

## License

This project is licensed under the MIT License — see [LICENSE](LICENSE) for details.