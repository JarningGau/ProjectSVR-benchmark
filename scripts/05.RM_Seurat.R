#### 导入相关的包 ####
library(Seurat)
library(tidyverse)
# set up future for parallelization
library(future)
library(future.apply)
plan("multicore", workers = 10)
options(future.globals.maxSize = 100 * 1024^3)
setwd(here::here())
source("R/IO.R")
dir.create("ref_models/05_Seurat", recursive = T)
dir.create("results/05_Seurat", recursive = T)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ PBMC task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# seu.ref <- LoadSeuratSlimData("data/DISCO_hPBMC_seurat.train.slim.qs")
# #### 1.Model building ####
# obj.list <- SplitObject(seu.ref, split.by = "project_id")
# for (i in 1:length(obj.list)) {
#   obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
#   obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
#                                         selection.method = "vst",
#                                         nfeatures = 2000,
#                                         verbose = FALSE)
# }
# vfeatures <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)
# 
# obj.list <- lapply(X = obj.list, FUN = function(x) {
#   x <- ScaleData(x, features = vfeatures, verbose = FALSE)
#   x <- RunPCA(x, features = vfeatures, verbose = FALSE)
#   x
# })
# 
# 
# TIME0 <- Sys.time()
# reference_dataset <- which(sapply(obj.list, ncol) > 5000)
# anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:20, reference = reference_dataset, reduction = "rpca")
# ## Return corrected normalized gene expression matrix storing in seu.ref[["integrated"]]@data
# seu.ref <- IntegrateData(anchorset = anchors, dims = 1:20)
# DefaultAssay(seu.ref) <- "integrated"
# TIME1 <- Sys.time()
# runtime1 <- difftime(TIME1, TIME0, units = "secs")
# 
# seu.ref <- ScaleData(seu.ref, verbose = FALSE)
# seu.ref <- RunPCA(seu.ref, npcs = 20, verbose = FALSE)
# seu.ref <- RunUMAP(seu.ref, reduction = "pca", dims = 1:20, return.model = TRUE, verbose = FALSE)
# 
# ## 保存Reference data
# qs::qsave(seu.ref, "ref_models/05_Seurat/PBMC_reference.seurat.qs")
# 
# #### 2.Reference mapping ####
# seu.q <- LoadSeuratSlimData("data/DISCO_hPBMC_seurat.test.slim.qs")
# 
# TIME0 <- Sys.time()
# ## step1: Find anchors between query and ref: ~ 2 min
# anchors <- FindTransferAnchors(
#   reference = seu.ref,
#   query = seu.q,
#   dims = 1:20,
#   reference.reduction = "pca",
#   verbose = T)
# 
# ## step2: map query to reference (fast!)
# seu.q <- MapQuery(
#   anchorset = anchors,
#   reference = seu.ref,
#   query = seu.q,
#   reference.reduction = "pca",
#   reduction.model = "umap")
# 
# TIME1 <- Sys.time()
# runtime2 <- difftime(TIME1, TIME0, units = "secs")
# 
# #### 3. Save results ++ ####
# ref.emb <- seu.ref[["umap"]]@cell.embeddings
# query.emb <- seu.q[["ref.umap"]]@cell.embeddings
# ref.cellmeta <- seu.ref@meta.data["cell_subtype"]
# ref.cellmeta$group <- "reference"
# query.cellmeta <- seu.q@meta.data["cell_subtype"]
# query.cellmeta$group <- "query"
# 
# colnames(ref.emb) <- paste0("Dim_", 1:2)
# colnames(query.emb) <- paste0("Dim_", 1:2)
# colnames(ref.cellmeta) <- c("label", "group")
# colnames(query.cellmeta) <- c("label", "group")
# ## Run time
# runtime <- data.frame(
#   model.building = runtime1,
#   reference.mapping = runtime2
# )
# 
# write.csv(ref.emb, "results/05_Seurat/PBMC_ref_embeddings.csv")
# write.csv(ref.cellmeta, "results/05_Seurat/PBMC_ref_cellmeta.csv")
# write.csv(query.emb, "results/05_Seurat/PBMC_query_embeddings.csv")
# write.csv(query.cellmeta, "results/05_Seurat/PBMC_query_cellmeta.csv")
# write.csv(runtime, "results/05_Seurat/PBMC_runtime.csv")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ MFI task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# seu.ref <- LoadSeuratSlimData("data/Vento2018.MFI.seurat.train.slim.qs")
# #### 1.Model building ####
# obj.list <- SplitObject(seu.ref, split.by = "Method")
# for (i in 1:length(obj.list)) {
#   obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
#   obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
#                                         selection.method = "vst",
#                                         nfeatures = 2000,
#                                         verbose = FALSE)
# }
# vfeatures <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)
# 
# obj.list <- lapply(X = obj.list, FUN = function(x) {
#   x <- ScaleData(x, features = vfeatures, verbose = FALSE)
#   x <- RunPCA(x, features = vfeatures, verbose = FALSE)
#   x
# })
# 
# 
# TIME0 <- Sys.time()
# anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30, reduction = "rpca")
# ## Return corrected normalized gene expression matrix storing in seu.ref[["integrated"]]@data
# seu.ref <- IntegrateData(anchorset = anchors, dims = 1:30)
# DefaultAssay(seu.ref) <- "integrated"
# TIME1 <- Sys.time()
# runtime1 <- difftime(TIME1, TIME0, units = "secs")
# 
# seu.ref <- ScaleData(seu.ref, verbose = FALSE)
# seu.ref <- RunPCA(seu.ref, npcs = 30, verbose = FALSE)
# seu.ref <- RunUMAP(seu.ref, reduction = "pca", dims = 1:30, return.model = TRUE, verbose = FALSE)
# 
# ## 保存Reference data
# qs::qsave(seu.ref, "ref_models/05_Seurat/MFI_reference.seurat.qs")
# 
# #### 2.Reference mapping ####
# seu.q <- LoadSeuratSlimData("data/Vento2018.MFI.seurat.test.slim.qs")
# 
# TIME0 <- Sys.time()
# ## step1: Find anchors between query and ref: ~ 2 min
# anchors <- FindTransferAnchors(
#   reference = seu.ref,
#   query = seu.q,
#   dims = 1:30,
#   reference.reduction = "pca",
#   verbose = T)
# 
# ## step2: map query to reference (fast!)
# seu.q <- MapQuery(
#   anchorset = anchors,
#   reference = seu.ref,
#   query = seu.q,
#   reference.reduction = "pca",
#   reduction.model = "umap")
# 
# TIME1 <- Sys.time()
# runtime2 <- difftime(TIME1, TIME0, units = "secs")
# 
# #### 3. Save results ++ ####
# ref.emb <- seu.ref[["umap"]]@cell.embeddings
# query.emb <- seu.q[["ref.umap"]]@cell.embeddings
# ref.cellmeta <- seu.ref@meta.data["annotation"]
# ref.cellmeta$group <- "reference"
# query.cellmeta <- seu.q@meta.data["annotation"]
# query.cellmeta$group <- "query"
# 
# colnames(ref.emb) <- paste0("Dim_", 1:2)
# colnames(query.emb) <- paste0("Dim_", 1:2)
# colnames(ref.cellmeta) <- c("label", "group")
# colnames(query.cellmeta) <- c("label", "group")
# ## Run time
# runtime <- data.frame(
#   model.building = runtime1,
#   reference.mapping = runtime2
# )
# 
# write.csv(ref.emb, "results/05_Seurat/MFI_ref_embeddings.csv")
# write.csv(ref.cellmeta, "results/05_Seurat/MFI_ref_cellmeta.csv")
# write.csv(query.emb, "results/05_Seurat/MFI_query_embeddings.csv")
# write.csv(query.cellmeta, "results/05_Seurat/MFI_query_cellmeta.csv")
# write.csv(runtime, "results/05_Seurat/MFI_runtime.csv")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ mTCA task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seu.ref <- LoadSeuratSlimData("data/mTCA.seurat.train.slim.qs")
#### 1.Model building ####
obj.list <- SplitObject(seu.ref, split.by = "GSE_ID")
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- NormalizeData(obj.list[[i]], verbose = FALSE)
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
                                        selection.method = "vst",
                                        nfeatures = 2000,
                                        verbose = FALSE)
}
vfeatures <- SelectIntegrationFeatures(object.list = obj.list, nfeatures = 2000)

obj.list <- lapply(X = obj.list, FUN = function(x) {
  x <- ScaleData(x, features = vfeatures, verbose = FALSE)
  x <- RunPCA(x, features = vfeatures, verbose = FALSE)
  x
})


TIME0 <- Sys.time()
anchors <- FindIntegrationAnchors(object.list = obj.list, dims = 1:30, reduction = "rpca")
## Return corrected normalized gene expression matrix storing in seu.ref[["integrated"]]@data
seu.ref <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(seu.ref) <- "integrated"
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")

seu.ref <- ScaleData(seu.ref, verbose = FALSE)
seu.ref <- RunPCA(seu.ref, npcs = 30, verbose = FALSE)
seu.ref <- RunUMAP(seu.ref, reduction = "pca", dims = 1:30, return.model = TRUE, verbose = FALSE)

## 保存Reference data
rm(obj.list)
gc()
qs::qsave(seu.ref, "ref_models/05_Seurat/mTCA_reference.seurat.qs")

#### 2.Reference mapping ####
seu.q <- LoadSeuratSlimData("data/mTCA.seurat.test.slim.qs")

TIME0 <- Sys.time()
## step1: Find anchors between query and ref: ~ 2 min
anchors <- FindTransferAnchors(
  reference = seu.ref,
  query = seu.q,
  dims = 1:30,
  reference.reduction = "pca",
  verbose = T)

## step2: map query to reference (fast!)
seu.q <- MapQuery(
  anchorset = anchors,
  reference = seu.ref,
  query = seu.q,
  reference.reduction = "pca",
  reduction.model = "umap")

TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

#### 3. Save results ####
ref.emb <- seu.ref[["umap"]]@cell.embeddings
query.emb <- seu.q[["ref.umap"]]@cell.embeddings
ref.cellmeta <- seu.ref@meta.data["Cell_type_symbol"]
ref.cellmeta$group <- "reference"
query.cellmeta <- seu.q@meta.data["Cell_type_symbol"]
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

write.csv(ref.emb, "results/05_Seurat/mTCA_ref_embeddings.csv")
write.csv(ref.cellmeta, "results/05_Seurat/mTCA_ref_cellmeta.csv")
write.csv(query.emb, "results/05_Seurat/mTCA_query_embeddings.csv")
write.csv(query.cellmeta, "results/05_Seurat/mTCA_query_cellmeta.csv")
write.csv(runtime, "results/05_Seurat/mTCA_runtime.csv")

