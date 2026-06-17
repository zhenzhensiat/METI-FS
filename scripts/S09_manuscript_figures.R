#!/usr/bin/env Rscript
# ==============================================================================
# S09_manuscript_figures.R — Publication-ready figures for METI-FS
#
# Figure mapping:
#   Fig. 2  Simulation benchmark              generate_fig2()
#   Fig. 3  Algorithm stability               generate_fig3()
#   Fig. 4  Ablation study                    generate_fig4()
#   Fig. 5  Cross-dataset screening funnel    generate_fig5()
#   Fig. 6  Upstream strategy comparison      generate_fig6()
#   Fig. 7  ML method comparison              generate_fig7()
#
# Figure mapping (v3.0):
#   Fig. 1  Pipeline schematic               (manual)
#   Fig. 2  Simulation benchmark              generate_fig2()
#   Fig. 3  Algorithm stability               generate_fig3()
#   Fig. 4  Ablation study                    generate_fig4()
#   Fig. 5  Cross-dataset screening funnel    generate_fig5()
#   Fig. 6  Upstream strategy comparison      generate_fig6()
#   Fig. 7  ML method comparison              generate_fig7()
#   Table S1                                  generate_table_S1()
#
# Usage:
#   source("S09_manuscript_figures.R")
#   generate_all_manuscript_figures()
# ==============================================================================

# ---- 0. Configuration ----
if (file.exists("S_config.R")) {
  source("S_config.R")

  
} else {
  BENCH_DIR <- BENCH_DIR
}

FIG_DIR <- file.path(BENCH_DIR, "manuscript_figures")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr)
  library(cowplot); library(scales)
})

# ---- Robust file finder ----
find_bench_file <- function(filename) {
  direct <- file.path(BENCH_DIR, filename)
  if (file.exists(direct)) return(direct)
  hits <- list.files(BENCH_DIR, pattern = paste0("^", filename, "$"),
                     recursive = TRUE, full.names = TRUE)
  if (length(hits) > 0) {
    if (length(hits) > 1)
      cat(sprintf("[S09] Multiple matches for %s — using: %s\n", filename, hits[1]))
    return(hits[1])
  }
  return(NA_character_)
}

# ---- Publication theme ----
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      plot.title       = element_text(size = base_size + 1, face = "bold", hjust = 0,
                                      margin = margin(b = 6)),
      plot.subtitle    = element_text(size = base_size - 1, color = "grey40",
                                      margin = margin(b = 8)),
      axis.title       = element_text(face = "bold", size = base_size),
      axis.text        = element_text(size = base_size - 1, color = "grey20"),
      legend.title     = element_text(face = "bold", size = base_size - 1),
      legend.text      = element_text(size = base_size - 1),
      legend.key.size  = unit(0.4, "cm"),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text       = element_text(face = "bold", size = base_size - 1),
      panel.grid.major.y = element_line(color = "grey92", linewidth = 0.3),
      plot.margin      = margin(10, 10, 10, 10)
    )
}

# ---- Palettes ----
PAL_SNR     <- c("low" = "#3C5488", "medium" = "#00A087", "high" = "#E64B35")
PAL_SAMPLE  <- c("small" = "#F39B7F", "medium" = "#4DBBD5", "large" = "#3C5488")
PAL_DENSITY <- c("sparse" = "#F39B7F", "medium" = "#4DBBD5", "dense" = "#3C5488")
PAL_METHOD  <- c("Single LASSO" = "#2E86AB", "Single Elastic Net" = "#A23B72",
                 "Stability Selection" = "#F18F01", "Boruta" = "#C73E1D",
                 "METI-FS ML" = "#00A087")
PAL_PIPE    <- c("P1 (DESeq2 LRT only)" = "#E64B35",
                 "P2 (P1 \u2229 WGCNA)" = "#4DBBD5",
                 "P3 (METI-FS)" = "#00A087")
PAL_DATASET <- c("GSE197067" = "#E64B35", "GSE307424" = "#4DBBD5",
                 "GSE236646" = "#00A087", "GSE150411" = "#3C5488")
LAB_DATASET <- c("GSE197067" = "GSE197067 (T cell)",
                 "GSE307424" = "GSE307424 (Lung cancer)",
                 "GSE236646" = "GSE236646 (HSV-1)",
                 "GSE150411" = "GSE150411 (Chondrocyte)")
DS_LEVELS   <- c("GSE197067", "GSE307424", "GSE236646", "GSE150411")

# ---- Pipeline run directories (for funnel data etc.) ----
RUN_DIRS <- list(
  GSE197067 = file.path(RUN_DIR, "GEO_GSE197067_Tcell"),
  GSE307424 = file.path(RUN_DIR, "GEO_GSE307424_Lung"),
  GSE236646 = file.path(RUN_DIR, "GEO_GSE236646_NPC"),
  GSE150411 = file.path(RUN_DIR, "GEO_GSE150411_Chon")
)

save_fig <- function(p, name, width, height) {
  for (ext in c("pdf", "png")) {
    fpath <- file.path(FIG_DIR, paste0(name, ".", ext))
    dev <- if (ext == "pdf" && capabilities("cairo")) cairo_pdf else ext
    ggsave(fpath, plot = p, width = width, height = height,
           dpi = 300, bg = "white", device = dev)
  }
  cat(sprintf("[S09] Saved: %s (.pdf + .png)\n", name))
}


# ==============================================================================
# Fig. 2  Simulation benchmark — 4-panel
# ==============================================================================
generate_fig2 <- function() {
  f <- find_bench_file("benchmark_master.csv")
  if (is.na(f)) { cat("[S09] benchmark_master.csv not found — skip Fig.2\n"); return(invisible()) }
  df  <- read.csv(f, stringsAsFactors = FALSE)
  sim <- df[df$mode == "simulation" & !is.na(df$precision), ]
  parts <- strsplit(sim$run_id, "_")
  sim$snr     <- factor(sapply(parts, `[`, 1), levels = c("low","medium","high"))
  sim$sample  <- factor(sapply(parts, `[`, 2), levels = c("small","medium","large"))
  sim$density <- factor(sapply(parts, `[`, 3), levels = c("sparse","medium","dense"))

  pA <- ggplot(sim, aes(snr, F1, fill = snr)) +
    geom_boxplot(width=.6, outlier.size=.8, alpha=.85) + geom_jitter(width=.15, size=.6, alpha=.3) +
    scale_fill_manual(values=PAL_SNR, guide="none") +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
    labs(title="A", subtitle="F1 score by signal-to-noise ratio", x="SNR", y="F1") + theme_pub()
  pB <- ggplot(sim, aes(snr, precision_any, fill=snr)) +
    geom_boxplot(width=.6, outlier.size=.8, alpha=.85) + geom_jitter(width=.15, size=.6, alpha=.3) +
    scale_fill_manual(values=PAL_SNR, guide="none") +
    scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1,.2)) +
    labs(title="B", subtitle="Precision_any by SNR", x="SNR", y=expression(Precision[any])) + theme_pub()
  pC <- ggplot(sim, aes(sample, F1, fill=sample)) +
    geom_boxplot(width=.6, outlier.size=.8, alpha=.85) + geom_jitter(width=.15, size=.6, alpha=.3) +
    scale_fill_manual(values=PAL_SAMPLE, guide="none") +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
    labs(title="C", subtitle="F1 score by sample size", x="Sample size category", y="F1") + theme_pub()
  pD <- ggplot(sim, aes(density, recall, fill=density)) +
    geom_boxplot(width=.6, outlier.size=.8, alpha=.85) + geom_jitter(width=.15, size=.6, alpha=.3) +
    scale_fill_manual(values=PAL_DENSITY, guide="none") +
    scale_y_continuous(limits=c(0,1), breaks=seq(0,1,.2)) +
    labs(title="D", subtitle="Recall by marker density", x="Marker density", y="Recall") + theme_pub()

  p <- plot_grid(pA, pB, pC, pD, nrow=2, align="hv")
  save_fig(p, "Fig02_simulation_benchmark", 10, 8)
  invisible(p)
}


# ==============================================================================
# Fig. 3  Algorithm stability — violin + scatter  (was Fig. 4 in v2.x)
# ==============================================================================
generate_fig3 <- function() {
  f <- find_bench_file("nogueira_stability_all.csv")
  if (is.na(f)) { cat("[S09] nogueira_stability_all.csv not found — skip Fig.3\n"); return(invisible()) }
  nog <- read.csv(f, stringsAsFactors = FALSE)
  nog$Method <- recode(nog$method, "LASSO"="LASSO","SVM"="SVM-RFE","RF"="RF","Gap-Union"="Gap-Union")
  nog$Method <- factor(nog$Method, levels=c("SVM-RFE","LASSO","Gap-Union","RF"))

  mu <- nog %>% group_by(Method) %>% summarise(mv=mean(nogueira,na.rm=TRUE), .groups="drop")

  pA <- ggplot(nog, aes(Method, nogueira, fill=Method)) +
    geom_violin(alpha=.7, quantiles=c(.25,.5,.75), linewidth=.4) +
    geom_jitter(width=.12, size=.5, alpha=.2) +
    geom_text(data=mu, aes(Method, mv, label=sprintf("%.3f",mv)), vjust=-1.5, size=3.3, fontface="bold") +
    scale_fill_manual(values=c("SVM-RFE"="#00A087","LASSO"="#2E86AB",
                                "Gap-Union"="#E64B35","RF"="#F39B7F"), guide="none") +
    scale_y_continuous(limits=c(-.1,1.05), breaks=seq(0,1,.2)) +
    labs(title="A", subtitle="Nogueira stability index (132 simulation scenarios)",
         x=NULL, y="Nogueira stability index") + theme_pub()

  fm <- find_bench_file("benchmark_master.csv")
  if (!is.na(fm)) {
    bm <- read.csv(fm, stringsAsFactors=FALSE)
    sim <- bm[bm$mode=="simulation" & !is.na(bm$precision),]
    pr <- data.frame(run_id=rep(sim$run_id,3),
      Method=rep(c("LASSO","RF","SVM-RFE"), each=nrow(sim)),
      Jaccard=c(sim$lasso_jaccard_mean, sim$rf_jaccard_mean, sim$svm_jaccard_mean),
      stringsAsFactors=FALSE)
    nw <- nog %>% select(run_id, Method, nogueira) %>% mutate(Method=as.character(Method))
    pr <- merge(pr, nw, by=c("run_id","Method"), all.x=TRUE)
    pr <- pr[complete.cases(pr),]
    pr$Method <- factor(pr$Method, levels=c("LASSO","RF","SVM-RFE"))

    pB <- ggplot(pr, aes(Jaccard, nogueira, color=Method)) +
      geom_point(alpha=.4, size=1.2) +
      geom_abline(slope=1, intercept=0, linetype="dashed", color="grey50") +
      scale_color_manual(values=c("LASSO"="#2E86AB","RF"="#F39B7F","SVM-RFE"="#00A087")) +
      scale_x_continuous(limits=c(0,1.05), breaks=seq(0,1,.25)) +
      scale_y_continuous(limits=c(0,1.05), breaks=seq(0,1,.25)) +
      annotate("text", x=.85, y=.15, label="RF: high Jaccard\nbut low Nogueira",
               color="#F39B7F", fontface="italic", size=3) +
      labs(title="B", subtitle="Jaccard vs. Nogueira (per scenario)",
           x="Pairwise Jaccard index", y="Nogueira stability index", color="Algorithm") +
      theme_pub() +
      theme(legend.position = "top", legend.justification = "left")

    p <- plot_grid(pA, pB, nrow=1)
    save_fig(p, "Fig03_algorithm_stability", 10, 5)
    invisible(p)
  } else {
    save_fig(pA, "Fig03_algorithm_stability", 6, 5)
    invisible(pA)
  }
}


# ==============================================================================
# Fig. 4  Ablation — A1/A2 paired differences (2-panel)  (was Fig. 3 in v2.x)
# ==============================================================================
generate_fig4 <- function() {
  f <- find_bench_file("ablation_A1_A2_full.csv")
  if (is.na(f)) { cat("[S09] ablation_A1_A2_full.csv not found — skip Fig.4\n"); return(invisible()) }
  abl  <- read.csv(f, stringsAsFactors = FALSE)
  runs <- unique(abl$run_id)

  paired <- do.call(rbind, lapply(runs, function(rid) {
    a0 <- abl[abl$run_id == rid & abl$config == "A0_FULL", ]
    a1 <- abl[abl$run_id == rid & abl$config == "A1_no_maSigPro", ]
    a2 <- abl[abl$run_id == rid & abl$config == "A2_no_WGCNA", ]
    out <- NULL
    if (nrow(a0)==1 && nrow(a1)==1)
      out <- rbind(out, data.frame(run_id=rid, config="A1: -maSigPro",
        dF1=a1$F1-a0$F1, dPrec=a1$precision-a0$precision, dRec=a1$recall-a0$recall, stringsAsFactors=FALSE))
    if (nrow(a0)==1 && nrow(a2)==1)
      out <- rbind(out, data.frame(run_id=rid, config="A2: -WGCNA",
        dF1=a2$F1-a0$F1, dPrec=a2$precision-a0$precision, dRec=a2$recall-a0$recall, stringsAsFactors=FALSE))
    out
  }))

  summ <- paired %>% group_by(config) %>%
    summarise(mean_dF1=mean(dF1), sd_dF1=sd(dF1), mean_dPrec=mean(dPrec), sd_dPrec=sd(dPrec),
              mean_dRec=mean(dRec), sd_dRec=sd(dRec), n=n(), .groups="drop")

  # Panel A
  label_y <- max(paired$dF1, na.rm=TRUE) + 0.12
  pA <- ggplot(paired, aes(config, dF1, fill=config)) +
    geom_hline(yintercept=0, linewidth=.5, linetype="dashed", color="grey50") +
    geom_boxplot(width=.5, alpha=.3, outlier.shape=NA, color="grey40") +
    geom_jitter(width=.18, size=1, alpha=.35, shape=16) +
    geom_point(data=summ, aes(config, mean_dF1), size=4, shape=18, color="black") +
    geom_text(data=summ, aes(config, y=label_y,
              label=sprintf("mean = %+.3f +/- %.3f", mean_dF1, sd_dF1)),
              size=3.2, fontface="bold") +
    scale_fill_manual(values=c("A1: -maSigPro"="#E64B35","A2: -WGCNA"="#4DBBD5"), guide="none") +
    coord_cartesian(ylim=c(min(paired$dF1, na.rm=TRUE)-.05, label_y+.06)) +
    labs(title="A",
         subtitle=expression(paste("Per-scenario ", Delta, "F1 (132 paired comparisons)")),
         x=NULL, y=expression(paste(Delta, "F1 (ablated - complete)"))) +
    theme_pub() + theme(panel.grid.major.x=element_blank())

  # Panel B
  decomp <- paired %>% select(run_id, config, dPrec, dRec) %>%
    pivot_longer(c(dPrec, dRec), names_to="metric", values_to="delta") %>%
    mutate(metric = factor(metric, levels=c("dPrec","dRec")))
  ds <- decomp %>% group_by(config, metric) %>%
    summarise(m=mean(delta), s=sd(delta), .groups="drop")

  pB <- ggplot(ds, aes(config, m, fill=metric)) +
    geom_hline(yintercept=0, linewidth=.5, linetype="dashed", color="grey50") +
    geom_col(position=position_dodge(.7), width=.6, alpha=.85) +
    geom_errorbar(aes(ymin=m-s, ymax=m+s), position=position_dodge(.7), width=.2, linewidth=.4) +
    geom_text(aes(label=sprintf("%+.3f", m),
                  y=ifelse(m<0, m-s-.04, m+s+.04)),
              position=position_dodge(.7), size=2.8, fontface="bold") +
    scale_fill_manual(
      values=c("dPrec"="#3C5488","dRec"="#F39B7F"),
      labels=c("dPrec"=expression(paste(Delta,"Precision")),
                "dRec"=expression(paste(Delta,"Recall")))) +
    coord_cartesian(ylim=c(min(ds$m-ds$s)-.08, max(ds$m+ds$s)+.08)) +
    labs(title="B", subtitle="Decomposition into precision and recall changes",
         x=NULL, y=expression(paste("Mean ",Delta," (ablated - complete)")), fill=NULL) +
    theme_pub() +
    theme(legend.position = c(.22, .92),
          legend.background=element_rect(fill=alpha("white",.8), color=NA),
          legend.text.align=0, panel.grid.major.x=element_blank())

  p <- plot_grid(pA, pB, nrow=1, rel_widths=c(1, 1.1))
  save_fig(p, "Fig04_ablation_study", 10, 5)
  invisible(p)
}


# ==============================================================================
# Fig. 5  Cross-dataset screening funnel — 2-panel  (NEW in v3.0)
# ==============================================================================
generate_fig5 <- function() {

  # ---- Load funnel data for each dataset ----
  # File is saved by 10_integration.R as "screening_funnel_data.csv" in each
  # dataset's data/ directory.  Also check BENCH_DIR for renamed copies
  # (e.g., screening_funnel_dataGSE197067.csv).
  all_data <- list()
  for (gse in DS_LEVELS) {
    f <- NA_character_
    # Strategy 1: check BENCH_DIR (with GSE suffix)
    f <- find_bench_file(paste0("screening_funnel_data", gse, ".csv"))
    # Strategy 2: check pipeline run directory
    if (is.na(f) && !is.null(RUN_DIRS[[gse]])) {
      candidate <- file.path(RUN_DIRS[[gse]], "data", "screening_funnel_data.csv")
      if (file.exists(candidate)) f <- candidate
    }
    # Strategy 3: recursive search under pipeline run directory
    if (is.na(f) && !is.null(RUN_DIRS[[gse]]) && dir.exists(RUN_DIRS[[gse]])) {
      hits <- list.files(RUN_DIRS[[gse]], pattern = "^screening_funnel_data",
                         recursive = TRUE, full.names = TRUE)
      if (length(hits) > 0) f <- hits[1]
    }
    if (is.na(f)) { cat(sprintf("[S09] WARNING: funnel data not found for %s\n", gse)); next }
    df <- read.csv(f, stringsAsFactors = FALSE)
    df$dataset <- gse
    all_data[[gse]] <- df
    cat(sprintf("[S09] Fig5: loaded %s (%d steps) from %s\n", gse, nrow(df), basename(dirname(f))))
  }
  if (length(all_data) < 2) {
    cat("[S09] Need >= 2 datasets for Fig.5 — skip\n"); return(invisible())
  }

  raw <- bind_rows(all_data)

  # ---- Standardize step names ----
  raw$step_std <- case_when(
    grepl("Filtered transcriptome", raw$step) ~ "Transcriptome",
    grepl("maSigPro.*interaction",  raw$step) ~ "maSigPro \u2229 LRT",
    grepl("maSigPro.*R",            raw$step) &
      !grepl("interaction",         raw$step) ~ "maSigPro (R\u00b2)",
    grepl("WGCNA",                  raw$step) ~ "\u2229 WGCNA",
    grepl("Effect size",            raw$step) ~ "Effect size\nfilter",
    grepl("Final",                  raw$step) ~ "Final",
    TRUE                                      ~ NA_character_
  )
  df <- raw %>% filter(!is.na(step_std))

  # If both maSigPro R² and interaction exist, keep only interaction
  df <- df %>%
    group_by(dataset) %>%
    mutate(.n_masig = sum(grepl("maSigPro", step_std))) %>%
    ungroup() %>%
    filter(!(.n_masig > 1 & step_std == "maSigPro (R\u00b2)")) %>%
    select(-.n_masig)
  df$step_std[grepl("maSigPro", df$step_std)] <- "maSigPro \u2229 LRT"

  step_levels <- c("Transcriptome", "maSigPro \u2229 LRT", "\u2229 WGCNA",
                   "Effect size\nfilter", "Final")
  df$step_std <- factor(df$step_std, levels = step_levels)
  df <- df %>% filter(!is.na(step_std))
  df$dataset  <- factor(df$dataset, levels = DS_LEVELS)

  # ---- Per-step reduction ----
  df <- df %>%
    arrange(dataset, step_std) %>%
    group_by(dataset) %>%
    mutate(n_prev = lag(n),
           pct_red = ifelse(!is.na(n_prev) & n_prev > 0,
                            round(100 * (1 - n / n_prev), 1), NA)) %>%
    ungroup()

  # ---- Panel A: gene count (log10) ----
  if (!requireNamespace("ggrepel", quietly = TRUE))
    install.packages("ggrepel", repos = "https://cloud.r-project.org", quiet = TRUE)

  pA <- ggplot(df, aes(x = step_std, y = n, color = dataset, group = dataset)) +
    geom_line(linewidth = 0.9, alpha = 0.85) +
    geom_point(size = 2.8) +
    ggrepel::geom_text_repel(
      aes(label = format(n, big.mark = ",")),
      size = 2.6, fontface = "bold", show.legend = FALSE,
      direction = "y",              # only repel vertically
      nudge_y = 0.08,               # slight upward nudge (log scale)
      segment.size = 0.3,           # thin connector lines
      segment.color = "grey60",
      segment.alpha = 0.5,
      min.segment.length = 0.3,     # show connectors when displaced
      box.padding = 0.2,
      point.padding = 0.15,
      max.overlaps = 20,
      seed = 42
    ) +
    scale_y_log10(labels = comma, breaks = c(10, 100, 1000, 10000),
                  limits = c(4, 40000)) +
    scale_color_manual(values = PAL_DATASET, labels = LAB_DATASET, name = "Dataset") +
    annotation_logticks(sides = "l", size = 0.3, color = "grey60") +
    labs(title = "A", subtitle = "Gene count at each pipeline step (log scale)",
         x = NULL, y = "Number of genes") +
    theme_pub() + theme(axis.text.x = element_text(size = 9, lineheight = 0.9))

  # ---- Panel B: per-step reduction % (exclude Final) ----
  df_red <- df %>% filter(!is.na(pct_red), step_std != "Final")

  pB <- ggplot(df_red, aes(x = step_std, y = pct_red, fill = dataset)) +
    geom_col(position = position_dodge(0.75), width = 0.65, alpha = 0.85) +
    geom_text(aes(label = paste0(pct_red, "%")),
              position = position_dodge(0.75), vjust = -0.4,
              size = 2.5, fontface = "bold") +
    annotate("rect", xmin = 2.5, xmax = 3.5, ymin = -2, ymax = 82,
             fill = "#FFF3CD", alpha = 0.35) +
    scale_fill_manual(values = PAL_DATASET, labels = LAB_DATASET, name = "Dataset") +
    scale_y_continuous(limits = c(0, 85), expand = c(0, 0)) +
    labs(title = "B", subtitle = "Per-step gene reduction (%)",
         x = NULL, y = "Reduction from previous step (%)") +
    theme_pub() +
    theme(axis.text.x = element_text(size = 9, lineheight = 0.9),
          legend.position = "none")

  # ---- Combine ----
  legend_grob <- get_legend(
    pA + theme(legend.position = "bottom", legend.direction = "horizontal",
               legend.margin = margin(t = 5))
  )
  p_top <- plot_grid(pA + theme(legend.position = "none"), pB,
                     ncol = 2, rel_widths = c(1.1, 1), align = "h", axis = "tb")
  p_final <- plot_grid(p_top, legend_grob, ncol = 1, rel_heights = c(1, 0.08))

  save_fig(p_final, "Fig05_cross_dataset_funnel", 12, 6)

  # ---- Summary ----
  cat("\n[S09] Fig5 — effect size filter reduction by dataset:\n")
  es <- df_red %>% filter(step_std == "Effect size\nfilter") %>% select(dataset, pct_red)
  print(as.data.frame(es), row.names = FALSE)

  invisible(p_final)
}


# ==============================================================================
# Fig. 6  Upstream strategy comparison — Layer 2 (was Fig. 5 in v2.x)
# ==============================================================================
generate_fig6 <- function() {
  f <- find_bench_file("layer2_paper_table.csv")
  if (is.na(f)) {
    cat("[S09] layer2_paper_table.csv not found — skip Fig.6\n")
    return(invisible())
  }
  cat(sprintf("[S09] Using: %s\n", f))
  l2 <- read.csv(f, stringsAsFactors=FALSE)
  l2$Strategy <- recode(l2$pipeline, "P1_DEG"="P1 (DESeq2 LRT only)",
    "P2_DEG_WGCNA"="P2 (P1 \u2229 WGCNA)", "P3_METIFS"="P3 (METI-FS)")
  l2$Strategy <- factor(l2$Strategy, levels=names(PAL_PIPE))
  l2$Backend  <- recode(l2$ml_backend, "SingleLASSO"="LASSO backend","Boruta"="Boruta backend")
  l2$dataset  <- factor(l2$dataset, levels=DS_LEVELS)

  make_panel <- function(data, tl, st) {
    ggplot(data, aes(dataset, jaccard, fill=Strategy)) +
      geom_col(position=position_dodge(.8), width=.7, alpha=.85) +
      scale_fill_manual(values=PAL_PIPE) +
      scale_y_continuous(limits=c(0,1.1), breaks=seq(0,1,.2)) +
      geom_text(aes(label=sprintf("%.3f",jaccard)), position=position_dodge(.8), vjust=-.4, size=2.6) +
      labs(title=tl, subtitle=st, x=NULL, y="Pairwise Jaccard index", fill="Upstream strategy") +
      theme_pub() + theme(legend.position="bottom", axis.text.x=element_text(angle=20,hjust=1))
  }
  pA <- make_panel(l2[l2$Backend=="LASSO backend",], "A", "LASSO backend: selection stability by upstream strategy")
  pB <- make_panel(l2[l2$Backend=="Boruta backend",], "B", "Boruta backend: selection stability by upstream strategy")
  leg <- get_legend(pA + theme(legend.position="bottom", legend.box.margin=margin(t=10)))
  p_top <- plot_grid(pA+theme(legend.position="none"), pB+theme(legend.position="none"), nrow=1)
  p <- plot_grid(p_top, leg, nrow=2, rel_heights=c(1,.1))
  save_fig(p, "Fig06_layer2_comparison", 11, 5.5)
  invisible(p)
}


# ==============================================================================
# Fig. 7  ML method comparison — Layer 3 (was Fig. 6 in v2.x)
# ==============================================================================
generate_fig7 <- function() {
  f <- find_bench_file("layer3_paper_table.csv")
  if (is.na(f)) {
    cat("[S09] layer3_paper_table.csv not found — skip Fig.7\n")
    return(invisible())
  }
  cat(sprintf("[S09] Using: %s\n", f))
  l3 <- read.csv(f, stringsAsFactors=FALSE)
  l3$Method <- recode(l3$method, "SingleLASSO"="Single LASSO","SingleEN"="Single Elastic Net",
    "StabilitySel"="Stability Selection","Boruta"="Boruta","METIFS_ML"="METI-FS ML")
  l3$Method <- factor(l3$Method, levels=names(PAL_METHOD))
  l3$dataset <- factor(l3$dataset, levels=DS_LEVELS)

  avg <- l3 %>% group_by(Method) %>%
    summarise(mean_n=mean(n_mean,na.rm=TRUE), mean_jac=mean(jaccard,na.rm=TRUE),
              mean_nog=mean(nogueira,na.rm=TRUE), .groups="drop")

  pA <- ggplot(avg, aes(Method, mean_n, fill=Method)) +
    geom_col(width=.65, alpha=.85) +
    geom_text(aes(label=sprintf("%.1f",mean_n)), vjust=-.4, size=3, fontface="bold") +
    scale_fill_manual(values=PAL_METHOD, guide="none") +
    scale_y_continuous(expand=expansion(mult=c(0,.15))) +
    labs(title="A", subtitle="Mean number of selected features", x=NULL, y="Mean n selected") +
    theme_pub() + theme(axis.text.x=element_text(angle=30,hjust=1))

  pB <- ggplot(l3, aes(dataset, jaccard, fill=Method)) +
    geom_col(position=position_dodge(.8), width=.7, alpha=.85) +
    scale_fill_manual(values=PAL_METHOD) +
    scale_y_continuous(limits=c(0,1.15), breaks=seq(0,1,.2)) +
    labs(title="B", subtitle="Pairwise Jaccard index by dataset", x=NULL,
         y="Pairwise Jaccard index", fill="Method") +
    theme_pub() + theme(legend.position="none", axis.text.x=element_text(angle=20,hjust=1))

  pC <- ggplot(l3, aes(dataset, nogueira, fill=Method)) +
    geom_col(position=position_dodge(.8), width=.7, alpha=.85) +
    scale_fill_manual(values=PAL_METHOD) +
    scale_y_continuous(limits=c(0,1.15), breaks=seq(0,1,.2)) +
    labs(title="C", subtitle="Nogueira stability index by dataset", x=NULL,
         y="Nogueira stability index", fill="Method") +
    theme_pub() + theme(legend.position="none", axis.text.x=element_text(angle=20,hjust=1))

  leg <- get_legend(pB + theme(legend.position="bottom", legend.box.margin=margin(t=10)) +
                      guides(fill=guide_legend(nrow=1)))
  p_top <- plot_grid(pA, pB, pC, nrow=1, rel_widths=c(.9,1.1,1.1))
  p <- plot_grid(p_top, leg, nrow=2, rel_heights=c(1,.08))
  save_fig(p, "Fig07_layer3_comparison", 13, 5.5)
  invisible(p)
}


# ==============================================================================
# Table S1
# ==============================================================================
generate_table_S1 <- function() {
  f <- find_bench_file("benchmark_master.csv")
  if (is.na(f)) return(invisible())
  df <- read.csv(f, stringsAsFactors=FALSE)
  sim <- df[df$mode=="simulation",]
  parts <- strsplit(sim$run_id, "_")
  sim$snr <- sapply(parts,`[`,1); sim$sample <- sapply(parts,`[`,2)
  sim$density <- sapply(parts,`[`,3); sim$rep <- sapply(parts,`[`,4)
  tab <- sim %>%
    select(run_id,snr,sample,density,rep,n_transcriptome,n_candidate_pool,
           n_final,precision,recall,F1,precision_any,FDR) %>%
    arrange(factor(snr,c("low","medium","high")),
            factor(sample,c("small","medium","large")),
            factor(density,c("sparse","medium","dense")))
  out <- file.path(FIG_DIR, "TableS1_simulation_parameters.csv")
  write.csv(tab, out, row.names=FALSE)
  cat(sprintf("[S09] Saved: TableS1_simulation_parameters.csv (%d scenarios)\n", nrow(tab)))
}


# ==============================================================================
# Master
# ==============================================================================
generate_all_manuscript_figures <- function() {
  cat("\n============================================================\n")
  cat("  S09 v3.0: Generating manuscript figures for METI-FS\n")
  cat(sprintf("  BENCH_DIR: %s\n", BENCH_DIR))
  cat(sprintf("  FIG_DIR:   %s\n", FIG_DIR))
  cat("============================================================\n\n")
  cat("[NOTE] Fig. 1 (pipeline schematic) — create manually.\n\n")
  generate_fig2()   # Fig. 2: Simulation benchmark
  generate_fig3()   # Fig. 3: Algorithm stability
  generate_fig4()   # Fig. 4: Ablation study
  generate_fig5()   # Fig. 5: Cross-dataset funnel
  generate_fig6()   # Fig. 6: Upstream strategy (Layer 2)
  generate_fig7()   # Fig. 7: ML method comparison (Layer 3)
  generate_table_S1()
  cat("\n============================================================\n")
  cat("  Done. Check FIG_DIR for all outputs.\n")
  cat("============================================================\n")
}

if (sys.nframe() == 0) {
  cat("\n  source('S09_manuscript_figures.R'); generate_all_manuscript_figures()\n\n")
}
