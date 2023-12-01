library(ProjecTILs)
library(STACAS)
library(Seurat)
library(tidyverse)
library(glue)
setwd(here::here())
dir.create("ref_models/04_projectTILs")
dir.create("results/04_projectTILs")
source("R/IO.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ PBMC task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# seu.ref <- LoadSeuratSlimData("data/DISCO_hPBMC_seurat.train.slim.qs")
# 
# #### 1. Model building ####
# obj.list <- SplitObject(seu.ref, split.by = "project_id")
# rm(seu.ref)
# gc()
# 
# TIME0 <- Sys.time()
# ## ~2549s
# seu.ref <-  Run.STACAS(object.list = obj.list, dims = 1:20, anchor.features = 2000)
# seu.ref <- RunUMAP(seu.ref, dims = 1:20)
# 
# ## 2293s
# seu.ref <- make.reference(
#   ref = seu.ref,
#   umap.method = 'umap',
#   recalculate.umap = TRUE,
#   atlas.name = "DISCO-hPBMC",
#   annotation.column = "cell_subtype")
# 
# TIME1 <- Sys.time()
# runtime1 <- difftime(TIME1, TIME0, units = "secs")
# 
# ## 保存Reference model
# qs::qsave(seu.ref, "ref_models/04_projectTILs/PBMC_reference.qs")
# 
# #### 2. Reference mapping ####
# seu.q <- LoadSeuratSlimData("data/DISCO_hPBMC_seurat.test.slim.qs")
# 
# ## 1119s
# TIME0 <- Sys.time()
# seu.q <- make.projection(query = seu.q, filter.cells = F, ref = seu.ref)
# TIME1 <- Sys.time()
# runtime2 <- difftime(TIME1, TIME0, units = "secs")
# 
# save.image("results/04_projectTILs/PBMC.rdata")
# load("results/04_projectTILs/PBMC.rdata")
# 
# #### 3. Save results ####
# ref.emb <- seu.ref[["umap"]]@cell.embeddings
# query.emb <- seu.q[["umap"]]@cell.embeddings
# ref.cellmeta <- seu.ref@meta.data["cell_subtype"]
# ref.cellmeta$group <- "reference"
# query.cellmeta <- seu.q@meta.data["cell_subtype"]
# query.cellmeta$group <- "query"
# 
# colnames(ref.emb) <- paste0("Dim_", 1:2)
# colnames(query.emb) <- paste0("Dim_", 1:2)
# colnames(ref.cellmeta) <- c("label", "group")
# colnames(query.cellmeta) <- c("label", "group")
# 
# write.csv(ref.emb, "results/04_projectTILs/PBMC_ref_embeddings.csv")
# write.csv(ref.cellmeta, "results/04_projectTILs/PBMC_ref_cellmeta.csv")
# write.csv(query.emb, "results/04_projectTILs/PBMC_query_embeddings.csv")
# write.csv(query.cellmeta, "results/04_projectTILs/PBMC_query_cellmeta.csv")
# 
# ## Run time
# runtime <- data.frame(
#   model.building = runtime1,
#   reference.mapping = runtime2
# )
# write.csv(runtime, "results/04_projectTILs/PBMC_runtime.csv")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ MFI task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seu.ref <- LoadSeuratSlimData("data/Vento2018.MFI.seurat.train.slim.qs")

#### 1. Model building ####
obj.list <- SplitObject(seu.ref, split.by = "orig.ident")
rm(seu.ref)
gc()

TIME0 <- Sys.time()
## ~2549s
seu.ref <- Run.STACAS(object.list = obj.list, dims = 1:30, anchor.features = 2000)
seu.ref <- RunUMAP(seu.ref, dims = 1:30)

## 2293s
seu.ref <- make.reference(
  ref = seu.ref,
  umap.method = 'umap',
  recalculate.umap = TRUE,
  atlas.name = "Vento-hMFI",
  annotation.column = "annotation")

TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")

## 保存Reference model
qs::qsave(seu.ref, "ref_models/04_projectTILs/MFI_reference.qs")
# seu.ref <- qs::qread("ref_models/04_projectTILs/MFI_reference.qs")
# DimPlot(seu.ref)

#### 2. Reference mapping ####
seu.q <- LoadSeuratSlimData("data/Vento2018.MFI.seurat.test.slim.qs")

## 1119s
TIME0 <- Sys.time()
seu.q <- make.projection(query = seu.q, filter.cells = F, ref = seu.ref)
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

save.image("results/04_projectTILs/MFI.rdata")
# load("results/04_projectTILs/MFI.rdata")

#### 3. Save results ####
celltype.col = "annotation"
task.name = "MFI"

ref.emb <- seu.ref[["umap"]]@cell.embeddings
query.emb <- seu.q[["umap"]]@cell.embeddings
ref.cellmeta <- seu.ref@meta.data[celltype.col]
ref.cellmeta$group <- "reference"
query.cellmeta <- seu.q@meta.data[celltype.col]
query.cellmeta$group <- "query"

colnames(ref.emb) <- paste0("Dim_", 1:2)
colnames(query.emb) <- paste0("Dim_", 1:2)
colnames(ref.cellmeta) <- c("label", "group")
colnames(query.cellmeta) <- c("label", "group")

write.csv(ref.emb,        glue("results/04_projectTILs/{task.name}_ref_embeddings.csv"))
write.csv(ref.cellmeta,   glue("results/04_projectTILs/{task.name}_ref_cellmeta.csv"))
write.csv(query.emb,      glue("results/04_projectTILs/{task.name}_query_embeddings.csv"))
write.csv(query.cellmeta, glue("results/04_projectTILs/{task.name}_query_cellmeta.csv"))

## Run time
runtime <- data.frame(
  model.building = runtime1,
  reference.mapping = runtime2
)
write.csv(runtime, glue("results/04_projectTILs/{task.name}_runtime.csv"))


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ mTCA task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# seu.ref <- LoadSeuratSlimData("data/mTCA.seurat.train.slim.qs")
# 
# #### 1. Model building ####
# obj.list <- SplitObject(seu.ref, split.by = "GSE_ID")
# rm(seu.ref)
# gc()
# 
# TIME0 <- Sys.time()
# ## ~2549s
# seu.ref <- Run.STACAS(object.list = obj.list, dims = 1:30, anchor.features = 2000)
# seu.ref <- RunUMAP(seu.ref, dims = 1:30)
# 
# ## 2293s
# seu.ref <- make.reference(
#   ref = seu.ref,
#   umap.method = 'umap',
#   recalculate.umap = TRUE,
#   atlas.name = "mTCA",
#   annotation.column = "Cell_type_symbol")
# 
# TIME1 <- Sys.time()
# runtime1 <- difftime(TIME1, TIME0, units = "secs")
# 
# ## 保存Reference model
# rm(obj.list)
# gc()
# qs::qsave(seu.ref, "ref_models/04_projectTILs/mTCA_reference.qs")
# save.image("results/04_projectTILs/mTCA.rdata")
# 
# DimPlot(seu.ref, reduction = "umap", group.by = "Cell_type_symbol", raster = F) +
#   scale_color_manual(values = ProjectSVR::pals$mtca)
# 
# #### 2. Reference mapping ####
# seu.q <- LoadSeuratSlimData("data/mTCA.seurat.test.slim.qs")
# 
# ## 1119s
# TIME0 <- Sys.time()
# seu.q <- make.projection(query = seu.q, filter.cells = F, ref = seu.ref)
# TIME1 <- Sys.time()
# runtime2 <- difftime(TIME1, TIME0, units = "secs")
# 
# save.image("results/04_projectTILs/mTCA.rdata")
# 
# #### 3. Save results ####
# celltype.col = "Cell_type_symbol"
# task.name = "mTCA"
# 
# ref.emb <- seu.ref[["umap"]]@cell.embeddings
# query.emb <- seu.q[["umap"]]@cell.embeddings
# ref.cellmeta <- seu.ref@meta.data[celltype.col]
# ref.cellmeta$group <- "reference"
# query.cellmeta <- seu.q@meta.data[celltype.col]
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
# write.csv(ref.emb,        glue("results/04_projectTILs/{task.name}_ref_embeddings.csv"))
# write.csv(ref.cellmeta,   glue("results/04_projectTILs/{task.name}_ref_cellmeta.csv"))
# write.csv(query.emb,      glue("results/04_projectTILs/{task.name}_query_embeddings.csv"))
# write.csv(query.cellmeta, glue("results/04_projectTILs/{task.name}_query_cellmeta.csv"))
# write.csv(runtime,        glue("results/04_projectTILs/{task.name}_runtime.csv"))


