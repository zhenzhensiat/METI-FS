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
├── scripts/                   # Simulation, benchmarking, and manuscript scripts
│   ├── S_config.R             # Shared configuration (paths, constants)
│   ├── S01_simulation_engine.R       # Synthetic data generator
│   ├── S01b_benchmark_timing.R       # Runtime profiling
│   ├── S02_pipeline_adapter.R        # Simulation → pipeline format adapter
│   ├── S03_download_and_prepare_GEO.R # GEO dataset download/preparation
│   ├── S04_run_pipeline_wrapper.R    # GEO dataset pipeline runner
│   ├── S04_sim_runner.R              # Simulation pipeline runner
│   ├── S05_benchmark_collector.R     # Result collection + metrics
│   ├── S06_ablation_study.R          # Component ablation analysis
│   ├── S07_cross_dataset_summary.R   # Cross-dataset summary + figures
│   ├── S08_baseline_comparison.R     # Baseline method comparison
│   ├── S09_manuscript_figures.R      # Publication figure generation
│   └── run_demo.R                    # Quick-start demo script
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