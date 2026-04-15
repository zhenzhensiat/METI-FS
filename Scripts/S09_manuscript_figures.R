#!/usr/bin/env Rscript
# ==============================================================================
# S09_manuscript_figures.R  v2.3 — Publication-ready figures for METI-FS
# 2026-04-13
#
# v2.3 changes: Fig3B legend → top-left; Fig4B legend → outside top
# v2.2 changes: recursive file search; cairo_pdf; coord_cartesian; expression()
# v2.0 changes: Fig3 redesigned for A1/A2 only
#
# Usage:
#   source("S09_manuscript_figures.R")
#   generate_all_manuscript_figures()
# ==============================================================================

# ---- 0. Configuration ----
if (file.exists("S_config.R")) {
  source("S_config.R")
} else if (file.exists(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))) {
  source(file.path(file.path(METHODS_BASE, "Scripts"), "S_config.R"))
} else {
  BENCH_DIR <- file.path(METHODS_BASE, "benchmark_results")
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
DS_LEVELS   <- c("GSE197067", "GSE307424", "GSE236646", "GSE150411")

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
# Fig. 3  Ablation — A1/A2 paired differences (2-panel)
# ==============================================================================
generate_fig3 <- function() {
  f <- find_bench_file("ablation_A1_A2_full.csv")
  if (is.na(f)) { cat("[S09] ablation_A1_A2_full.csv not found — skip Fig.3\n"); return(invisible()) }
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

  # Panel B — legend at TOP-LEFT to avoid overlap
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
  save_fig(p, "Fig03_ablation_study", 10, 5)
  invisible(p)
}


# ==============================================================================
# Fig. 4  Algorithm stability — violin + scatter
# ==============================================================================
generate_fig4 <- function() {
  f <- find_bench_file("nogueira_stability_all.csv")
  if (is.na(f)) { cat("[S09] nogueira_stability_all.csv not found — skip Fig.4\n"); return(invisible()) }
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

    # Legend OUTSIDE plot (top) to avoid data overlap
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
    save_fig(p, "Fig04_algorithm_stability", 10, 5)
    invisible(p)
  } else {
    save_fig(pA, "Fig04_algorithm_stability", 6, 5)
    invisible(pA)
  }
}


# ==============================================================================
# Fig. 5  Layer 2 — upstream strategy comparison (2 panels)
# ==============================================================================
generate_fig5 <- function() {
  f <- find_bench_file("layer2_paper_table.csv")
  if (is.na(f)) {
    cat("[S09] layer2_paper_table.csv not found anywhere under BENCH_DIR — skip Fig.5\n")
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
  save_fig(p, "Fig05_layer2_comparison", 11, 5.5)
  invisible(p)
}


# ==============================================================================
# Fig. 6  Layer 3 — ML method comparison (3 panels)
# ==============================================================================
generate_fig6 <- function() {
  f <- find_bench_file("layer3_paper_table.csv")
  if (is.na(f)) {
    cat("[S09] layer3_paper_table.csv not found anywhere under BENCH_DIR — skip Fig.6\n")
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
  save_fig(p, "Fig06_layer3_comparison", 13, 5.5)
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
  cat("  S09 v2.3: Generating manuscript figures for METI-FS\n")
  cat(sprintf("  BENCH_DIR: %s\n", BENCH_DIR))
  cat(sprintf("  FIG_DIR:   %s\n", FIG_DIR))
  cat("============================================================\n\n")
  cat("[NOTE] Fig. 1 (pipeline schematic) — create manually.\n\n")
  generate_fig2(); generate_fig3(); generate_fig4()
  generate_fig5(); generate_fig6(); generate_table_S1()
  cat("\n============================================================\n")
  cat("  Done. Check FIG_DIR for all outputs.\n")
  cat("============================================================\n")
}

if (sys.nframe() == 0) {
  cat("\n  source('S09_manuscript_figures.R'); generate_all_manuscript_figures()\n\n")
}
