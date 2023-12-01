library(tidyverse)
library(glue)
library(ProjectSVR)
library(tidydr)
library(Seurat)
library(patchwork)
library(ggrastr)
setwd(here::here())
dir.create("results/09_summary_metrics")
#' 1. Run time
#' - model building
#' - reference mapping

#' 2. Model size

#' 3. Results
#' results
#' files:
#' - reference umap embeddings        <data_set>_ref_embeddings.csv
#' - reference cell labels, group     <data_set>_ref_cellmeta.csv
#' - query umap embeddings            <data_set>_query_embeddings.csv
#' - query cell labels (TRUE labels), <data_set>_query_cellmeta.csv
#'   group, donor id
#' calculation:
#' - step1: KNN(k=10) label transfer => predicted labels
#' - step2: true labels vs predicted labels => confusion matrix
#' - step3: confusion matrix => ACC/ARI
#' - step4: LISI
#' plots
#' projection plot (ref to symphony fig 4c,d)

ComputeMetrics <- function(software = "02_ProjectSVR", task = "PBMC") {
  ref.emb.file <- glue("results/{software}/{task}_ref_embeddings.csv")
  ref.cellmeta.file <- glue("results/{software}/{task}_ref_cellmeta.csv")
  query.emb.file <- glue("results/{software}/{task}_query_embeddings.csv")
  query.cellmeta.file <- glue("results/{software}/{task}_query_cellmeta.csv")
  ## Loading results
  ref.emb <- read.csv(ref.emb.file, row.names = 1, header = T)
  ref.cellmeta <- read.csv(ref.cellmeta.file, row.names = 1, header = T)
  query.emb <- read.csv(query.emb.file, row.names = 1, header = T)
  query.cellmeta <- read.csv(query.cellmeta.file, row.names = 1, header = T)
  ## LISI
  X <- rbind(ref.emb, query.emb)
  meta.data <- rbind(ref.cellmeta, query.cellmeta[, 1:2])
  lisi <- lisi::compute_lisi(X, meta.data, label_colnames = c("label", "group"))
  ## KNN label transfer (k=10)
  ref.labels <- ref.cellmeta$label
  names(ref.labels) <- rownames(ref.cellmeta)
  knn.pred <- ProjectSVR::KnnLabelTransfer(query.emb, ref.emb, ref.labels, k = 10)
  query.cellmeta$pred.label <- knn.pred$labels
  ## confusion matrix
  true.labels <- query.cellmeta$label
  pred.labels <- query.cellmeta$pred.label
  label.levels <- unique(c(true.labels, pred.labels))
  true.labels <- factor(true.labels, levels = label.levels)
  pred.labels <- factor(pred.labels, levels = label.levels)
  confusion.matrix <- as.matrix(table(true.labels, pred.labels))
  ## ACC/ARI
  acc <- sum(diag(confusion.matrix) / sum(confusion.matrix))
  ari <- mclust::adjustedRandIndex(true.labels, pred.labels)
  ## Runtime 
  runtime <- read.csv(glue("results/{software}/{task}_runtime.csv"), row.names = 1, header = T)
  ## Save results
  metrics <- data.frame(
    software = software,
    task = task,
    accuracy = acc,
    ARI = ari
  )
  metrics <- cbind(metrics, runtime)
  lisi$software <- software
  lisi$task <- task
  
  return(list(
    confusion.matrix = confusion.matrix,
    metrics = metrics,
    lisi = lisi
  ))
}

## Calculation
softwares <- c("01_Golden", "02_ProjectSVR", "03_symphony", "04_projectTILs", "05_Seurat", "06_scArches", "07_iNMF", "08_SCALEX")
tasks <- c("PBMC", "MFI", "mTCA")

for (tt in tasks) {
  for (sw in softwares) {
    res <- ComputeMetrics(software = sw, task = tt)
    saveRDS(res, glue("results/09_summary_metrics/{sw}_{tt}.rds"))
  }
}

## Metrics Summary Plots
softwares <- c("01_Golden", "02_ProjectSVR", "03_symphony", "04_projectTILs", "05_Seurat", "06_scArches", "07_iNMF", "08_SCALEX")
tasks <- c("PBMC", "MFI", "mTCA")

LoadSoftRes <- function(sw) {
  lapply(tasks, function(tt) {
    res <- readRDS(glue("results/09_summary_metrics/{sw}_{tt}.rds"))
    res$metrics
  }) %>% Reduce(rbind, .)
}
metrics <- lapply(softwares, LoadSoftRes) %>% Reduce(rbind, .)
metrics$task <- factor(metrics$task, levels = tasks)
metrics$software <- sapply(strsplit(metrics$software, split = "_"), function(xx) xx[2])
metrics[1:3, "software"] <- "De novo integration"
metrics$software <- factor(metrics$software, levels = unique(metrics$software))
metrics

write_tsv(metrics, "results/09_summary_metrics/summary_metrics.tsv")

## Projection Accuracy

p1 <- ggplot(metrics, aes(task, accuracy)) + 
  geom_bar(stat = 'identity', aes(fill = software), position = position_dodge(0.9)) + 
  scale_fill_manual(values = c("black", ggsci::pal_d3()(10))) + 
  ylab("Accuracy") + 
  theme_bw(base_size = 15) + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none")

p2 <- ggplot(metrics, aes(task, ARI)) + 
  geom_bar(stat = 'identity', aes(fill = software), position = position_dodge(0.9)) + 
  ylab("Adjusted rand index (ARI)") + 
  scale_fill_manual(values = c("black", ggsci::pal_d3()(10))) + 
  theme_bw(base_size = 15) + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"))
p1 + p2
ggsave("results/09_summary_metrics/performance-1.png", width = 12, height = 5)

## Speed

p3 <- ggplot(metrics, aes(task, model.building)) + 
  geom_bar(stat = 'identity', aes(fill = software), position = position_dodge(0.9)) + 
  ylab("Time (s)") + ggtitle("Model building") + 
  ggsci::scale_fill_d3() + 
  scale_y_log10() + 
  theme_bw(base_size = 15) + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        legend.position = "none",
        plot.title = element_text(hjust = .5, face = "bold"))

p4 <- ggplot(metrics, aes(task, reference.mapping)) + 
  geom_bar(stat = 'identity', aes(fill = software), position = position_dodge(0.9)) + 
  ylab("Time (s)") + ggtitle("Reference mapping") + 
  ggsci::scale_fill_d3() + 
  scale_y_log10() + 
  theme_bw(base_size = 15) + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"),
        plot.title = element_text(hjust = .5, face = "bold"))

p3 + p4
ggsave("results/09_summary_metrics/performance-2.png", width = 12, height = 5)

## LISI
LoadSoftRes <- function(sw) {
  lapply(tasks, function(tt) {
    res <- readRDS(glue("results/09_summary_metrics/{sw}_{tt}.rds"))
    res$lisi
  }) %>% Reduce(rbind, .)
}
lisi <- lapply(softwares, LoadSoftRes) %>% Reduce(rbind, .)
lisi$task <- factor(lisi$task, levels = tasks)
lisi$software <- sapply(strsplit(lisi$software, split = "_"), function(xx) xx[2])
lisi$software <- factor(lisi$software, levels = unique(lisi$software))

ggplot(lisi, aes(task, group)) + 
  geom_boxplot(aes(fill = software), outlier.shape = NA) + 
  ggsci::scale_fill_d3() + 
  ylab("LISI of group") + 
  theme_bw(base_size = 15) + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"))

ggplot(lisi, aes(task, label)) + 
  geom_boxplot(aes(fill = software), outlier.shape = NA) + 
  ggsci::scale_fill_d3() + 
  coord_cartesian(ylim = c(1, 4)) + 
  ylab("LISI of cell label") + 
  theme_bw(base_size = 15) + 
  theme(axis.title.x = element_blank(),
        axis.text = element_text(color = "black"))

## Projection Plots
ProjectionPlot <- function(software, task, colors) {
  ref.emb.file <- glue("results/{software}/{task}_ref_embeddings.csv")
  ref.cellmeta.file <- glue("results/{software}/{task}_ref_cellmeta.csv")
  query.emb.file <- glue("results/{software}/{task}_query_embeddings.csv")
  query.cellmeta.file <- glue("results/{software}/{task}_query_cellmeta.csv")
  ## Loading results
  ref.emb <- read.csv(ref.emb.file, row.names = 1, header = T)
  ref.cellmeta <- read.csv(ref.cellmeta.file, row.names = 1, header = T)
  query.emb <- read.csv(query.emb.file, row.names = 1, header = T)
  query.cellmeta <- read.csv(query.cellmeta.file, row.names = 1, header = T)
  
  ref.data <- cbind(ref.emb, ref.cellmeta[rownames(ref.emb), ])
  query.data <- cbind(query.emb, query.cellmeta[rownames(query.emb), ])
  
  if (!is.null(names(colors))) {
    ref.data$label <- factor(ref.data$label, levels = names(colors))
    query.data$label <- factor(query.data$label, levels = names(colors))
  }
  p1 <- ggplot(ref.data, aes(Dim_1, Dim_2, color = label)) + 
    geom_point(size = .1, alpha = .1) +
    scale_color_manual(values = colors) + 
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) + 
    theme_dr() + 
    theme(panel.grid = element_blank(),
          plot.title = element_text(hjust = .5),
          legend.title = element_blank(), 
          legend.key.height = unit (8, "pt"))

  p2 <- ggplot() + 
    geom_point(data = ref.data, aes(Dim_1, Dim_2), size = .1, color = "grey") + 
    geom_point(data = query.data, aes(Dim_1, Dim_2, color = label), size = .1) + 
    scale_color_manual(values = colors) + 
    theme_dr() + 
    theme(panel.grid = element_blank(), legend.position = "none")
  list(p1, p2)
}

softwares <- c("02_ProjectSVR", "03_symphony", "04_projectTILs", "05_Seurat", "06_scArches", "07_iNMF", "08_SCALEX")
tasks <- c("PBMC", "MFI", "mTCA")

for (sw in softwares) {
  ProjectionPlot(software = sw, task = "PBMC", colors = pals$disco_blood)
  ggsave(glue("results/09_summary_metrics/plots/{sw}_PBMC.png"), width = 14, height = 6)
  ProjectionPlot(software = sw, task = "MFI", colors = pals$hfmi)
  ggsave(glue("results/09_summary_metrics/plots/{sw}_MFI.png"), width = 14, height = 6)
  ProjectionPlot(software = sw, task = "mTCA", colors = pals$mtca)
  ggsave(glue("results/09_summary_metrics/plots/{sw}_mTCA.png"), width = 15, height = 6)
}

#### PBMC ####

p.list <- ProjectionPlot(software = "02_ProjectSVR", task = "PBMC", colors = pals$disco_blood)
p1 <- p.list[[1]] + ggtitle("ProjectSVR", "Reference")
p2 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "03_symphony", task = "PBMC", colors = pals$disco_blood)
p3 <- p.list[[1]] + ggtitle("Symphony", "Reference")
p4 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "04_projectTILs", task = "PBMC", colors = pals$disco_blood)
p5 <- p.list[[1]] + ggtitle("ProjectTILs", "Reference")
p6 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "05_Seurat", task = "PBMC", colors = pals$disco_blood)
p7 <- p.list[[1]] + ggtitle("Seurat", "Reference")
p8 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "06_scArches", task = "PBMC", colors = pals$disco_blood)
p9 <- p.list[[1]] + ggtitle("scArches", "Reference")
p10 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "07_iNMF", task = "PBMC", colors = pals$disco_blood)
p11 <- p.list[[1]] + ggtitle("iNMF", "Reference")
p12 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "08_SCALEX", task = "PBMC", colors = pals$disco_blood)
p13 <- p.list[[1]] + ggtitle("SCALEX", "Reference")
p14 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")

p1 + p2 + p3 + p4 + 
  p5 + p6 + p7 + p8 + 
  p9 + p10 + p11 + p12 + 
  p13 + p14 + 
  plot_layout(nrow = 2, byrow = F, guides = "collect")
ggsave("results/09_summary_metrics/Projection-PBMC.png", width = 21, height = 6)  

#### hMFI ####

p.list <- ProjectionPlot(software = "02_ProjectSVR", task = "MFI", colors = pals$hfmi)
p1 <- p.list[[1]] + ggtitle("ProjectSVR", "Reference")
p2 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "03_symphony", task = "MFI", colors = pals$hfmi)
p3 <- p.list[[1]] + ggtitle("Symphony", "Reference")
p4 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "04_projectTILs", task = "MFI", colors = pals$hfmi)
p5 <- p.list[[1]] + ggtitle("ProjectTILs", "Reference")
p6 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "05_Seurat", task = "MFI", colors = pals$hfmi)
p7 <- p.list[[1]] + ggtitle("Seurat", "Reference")
p8 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "06_scArches", task = "MFI", colors = pals$hfmi)
p9 <- p.list[[1]] + ggtitle("scArches", "Reference")
p10 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "07_iNMF", task = "MFI", colors = pals$hfmi)
p11 <- p.list[[1]] + ggtitle("iNMF", "Reference")
p12 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "08_SCALEX", task = "MFI", colors = pals$hfmi)
p13 <- p.list[[1]] + ggtitle("SCALEX", "Reference")
p14 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")

p1 + p2 + p3 + p4 + 
  p5 + p6 + p7 + p8 + 
  p9 + p10 + p11 + p12 + 
  p13 + p14 + 
  plot_layout(nrow = 2, byrow = F, guides = "collect")
ggsave("results/09_summary_metrics/Projection-MFI.png", width = 21, height = 6)  


#### mTCA ####

p.list <- ProjectionPlot(software = "02_ProjectSVR", task = "mTCA", colors = pals$mtca)
p1 <- p.list[[1]] + ggtitle("ProjectSVR", "Reference")
p2 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "03_symphony", task = "mTCA", colors = pals$mtca)
p3 <- p.list[[1]] + ggtitle("Symphony", "Reference")
p4 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "04_projectTILs", task = "mTCA", colors = pals$mtca)
p5 <- p.list[[1]] + ggtitle("ProjectTILs", "Reference")
p6 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "05_Seurat", task = "mTCA", colors = pals$mtca)
p7 <- p.list[[1]] + ggtitle("Seurat", "Reference")
p8 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "06_scArches", task = "mTCA", colors = pals$mtca)
p9 <- p.list[[1]] + ggtitle("scArches", "Reference")
p10 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "07_iNMF", task = "mTCA", colors = pals$mtca)
p11 <- p.list[[1]] + ggtitle("iNMF", "Reference")
p12 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")
p.list <- ProjectionPlot(software = "08_SCALEX", task = "mTCA", colors = pals$mtca)
p13 <- p.list[[1]] + ggtitle("SCALEX", "Reference")
p14 <- p.list[[2]] + ggtitle(label = "", subtitle = "Mapped query")

p1 + p2 + p3 + p4 + 
  p5 + p6 + p7 + p8 + 
  p9 + p10 + p11 + p12 + 
  p13 + p14 + 
  plot_layout(nrow = 2, byrow = F, guides = "collect")
ggsave("results/09_summary_metrics/Projection-mTCA.png", width = 21, height = 6)  




