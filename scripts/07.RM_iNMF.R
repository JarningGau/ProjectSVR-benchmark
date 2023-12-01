library(rliger)
library(Seurat)
library(tidyverse)
library(glue)
setwd(here::here())
dir.create("ref_models/07_iNMF")
dir.create("results/07_iNMF")
source("R/IO.R")
source("R/iNMF_utils.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ PBMC task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seu.ref <- LoadSeuratSlimData("data/DISCO_hPBMC_seurat.train.slim.qs")
liger.ref <- Seurat2LigerRef(seu.ref, split.by = "project_id", nfeatures = 2000)
gc()
liger.ref@clusters <- factor(seu.ref$cell_subtype[rownames(liger.ref@cell.data)])

#### 1. Model building ####
TIME0 <- Sys.time()
liger.ref <- online_iNMF(liger.ref, k = 20, miniBatch_size = 5000, max.epochs = 5)
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")
gc()

liger.ref <- quantile_norm(liger.ref)
liger.ref <- runUMAP(liger.ref)

# plotByDatasetAndCluster(liger.ref, axis.labels = c("UMAP1", "UMAP2"))

## 保存参考模型
qs::qsave(liger.ref, 'ref_models/07_iNMF/PBMC_ref_model.liger.qs')

#### 2. Reference mapping ####
liger.ref <- qs::qread("ref_models/07_iNMF/PBMC_ref_model.liger.qs")
seu.q <- LoadSeuratSlimData("data/DISCO_hPBMC_seurat.test.slim.qs")
liger.query <- Seurat2LigerQuery(seu.q, vfeatures = liger.ref@var.genes)

TIME0 <- Sys.time()
liger.ref <- online_iNMF(liger.ref, X_new = list(query = liger.query), k = 20, max.epochs = 1, projection = TRUE)
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

liger.ref <- quantile_norm(liger.ref)
liger.ref <- runUMAP(liger.ref)

#### 3. Save results ####
ref.emb <- liger.ref@tsne.coords[colnames(seu.ref), ] %>% as.data.frame()
query.emb <- liger.ref@tsne.coords[colnames(seu.q), ] %>% as.data.frame()
ref.cellmeta <- seu.ref@meta.data["cell_subtype"]
ref.cellmeta$group <- "reference"
query.cellmeta <- seu.q@meta.data["cell_subtype"]
query.cellmeta$group <- "query"

colnames(ref.emb) <- paste0("Dim_", 1:2)
colnames(query.emb) <- paste0("Dim_", 1:2)
colnames(ref.cellmeta) <- c("label", "group")
colnames(query.cellmeta) <- c("label", "group")

write.csv(ref.emb, "results/07_iNMF/PBMC_ref_embeddings.csv")
write.csv(ref.cellmeta, "results/07_iNMF/PBMC_ref_cellmeta.csv")
write.csv(query.emb, "results/07_iNMF/PBMC_query_embeddings.csv")
write.csv(query.cellmeta, "results/07_iNMF/PBMC_query_cellmeta.csv")

## Run time
runtime <- data.frame(
  model.building = runtime1,
  reference.mapping = runtime2
)
write.csv(runtime, "results/07_iNMF/PBMC_runtime.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ MFI task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seu.ref <- LoadSeuratSlimData("data/Vento2018.MFI.seurat.train.slim.qs")
liger.ref <- Seurat2LigerRef(seu.ref, split.by = "Method", nfeatures = 2000)
gc()

#### 1. Model building ####
TIME0 <- Sys.time()
liger.ref <- online_iNMF(liger.ref, k = 40, miniBatch_size = 5000, max.epochs = 5)
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")
gc()

liger.ref <- quantile_norm(liger.ref)
liger.ref <- runUMAP(liger.ref)
liger.ref@clusters <- factor(seu.ref$annotation[rownames(liger.ref@cell.data)])

plotByDatasetAndCluster(liger.ref, axis.labels = c("UMAP1", "UMAP2"))

## 保存参考模型
qs::qsave(liger.ref, 'ref_models/07_iNMF/MFI_ref_model.liger.qs')

#### 2. Reference mapping ####
liger.ref <- qs::qread("ref_models/07_iNMF/MFI_ref_model.liger.qs")
seu.q <- LoadSeuratSlimData("data/Vento2018.MFI.seurat.test.slim.qs")
liger.query <- Seurat2LigerQuery(seu.q, vfeatures = liger.ref@var.genes)

TIME0 <- Sys.time()
liger.ref <- online_iNMF(liger.ref, X_new = list(query = liger.query), k = 40, max.epochs = 1, projection = TRUE)
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

liger.ref <- quantile_norm(liger.ref)
liger.ref <- runUMAP(liger.ref)

#### 3. Save results ####
celltype.col = "annotation"
task.name = "MFI"

ref.emb <- liger.ref@tsne.coords[colnames(seu.ref), ] %>% as.data.frame()
query.emb <- liger.ref@tsne.coords[colnames(seu.q), ] %>% as.data.frame()
ref.cellmeta <- seu.ref@meta.data[celltype.col]
ref.cellmeta$group <- "reference"
query.cellmeta <- seu.q@meta.data[celltype.col]
query.cellmeta$group <- "query"

colnames(ref.emb) <- paste0("Dim_", 1:2)
colnames(query.emb) <- paste0("Dim_", 1:2)
colnames(ref.cellmeta) <- c("label", "group")
colnames(query.cellmeta) <- c("label", "group")
## Run time
runtime <- data.frame(
  model.building = runtime1,
  reference.mapping = runtime2
)

write.csv(ref.emb,        glue("results/07_iNMF/{task.name}_ref_embeddings.csv"))
write.csv(ref.cellmeta,   glue("results/07_iNMF/{task.name}_ref_cellmeta.csv"))
write.csv(query.emb,      glue("results/07_iNMF/{task.name}_query_embeddings.csv"))
write.csv(query.cellmeta, glue("results/07_iNMF/{task.name}_query_cellmeta.csv"))
write.csv(runtime,        glue("results/07_iNMF/{task.name}_runtime.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ mTCA task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seu.ref <- LoadSeuratSlimData("data/mTCA.seurat.train.slim.qs")
liger.ref <- Seurat2LigerRef(seu.ref, split.by = "GSE_ID", nfeatures = 2000)
gc()

#### 1. Model building ####
TIME0 <- Sys.time()
liger.ref <- online_iNMF(liger.ref, k = 40, miniBatch_size = 5000, max.epochs = 5)
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")
gc()

liger.ref <- quantile_norm(liger.ref)
liger.ref <- runUMAP(liger.ref)
liger.ref@clusters <- factor(seu.ref$Cell_type_symbol[rownames(liger.ref@cell.data)])

plotByDatasetAndCluster(liger.ref, axis.labels = c("UMAP1", "UMAP2"))

## 保存参考模型
qs::qsave(liger.ref, 'ref_models/07_iNMF/mTCA_ref_model.liger.qs')

#### 2. Reference mapping ####
liger.ref <- qs::qread("ref_models/07_iNMF/mTCA_ref_model.liger.qs")
seu.q <- LoadSeuratSlimData("data/mTCA.seurat.test.slim.qs")
liger.query <- Seurat2LigerQuery(seu.q, vfeatures = liger.ref@var.genes)

TIME0 <- Sys.time()
liger.ref <- online_iNMF(liger.ref, X_new = list(query = liger.query), k = 40, max.epochs = 1, projection = TRUE)
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

liger.ref <- quantile_norm(liger.ref)
liger.ref <- runUMAP(liger.ref)

#### 3. Save results ####
celltype.col = "Cell_type_symbol"
task.name = "mTCA"

ref.emb <- liger.ref@tsne.coords[colnames(seu.ref), ] %>% as.data.frame()
query.emb <- liger.ref@tsne.coords[colnames(seu.q), ] %>% as.data.frame()
ref.cellmeta <- seu.ref@meta.data[celltype.col]
ref.cellmeta$group <- "reference"
query.cellmeta <- seu.q@meta.data[celltype.col]
query.cellmeta$group <- "query"

colnames(ref.emb) <- paste0("Dim_", 1:2)
colnames(query.emb) <- paste0("Dim_", 1:2)
colnames(ref.cellmeta) <- c("label", "group")
colnames(query.cellmeta) <- c("label", "group")
## Run time
runtime <- data.frame(
  model.building = runtime1,
  reference.mapping = runtime2
)

write.csv(ref.emb,        glue("results/07_iNMF/{task.name}_ref_embeddings.csv"))
write.csv(ref.cellmeta,   glue("results/07_iNMF/{task.name}_ref_cellmeta.csv"))
write.csv(query.emb,      glue("results/07_iNMF/{task.name}_query_embeddings.csv"))
write.csv(query.cellmeta, glue("results/07_iNMF/{task.name}_query_cellmeta.csv"))
write.csv(runtime,        glue("results/07_iNMF/{task.name}_runtime.csv"))

