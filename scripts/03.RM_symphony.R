library(symphony)
library(Seurat)
library(tidyverse)
library(glue)
setwd(here::here())
dir.create("ref_models/03_symphony")
dir.create("results/03_symphony")
source("R/IO.R")
source("R/symphony_utils.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ PBMC task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seu.ref <- LoadSeuratSlimData("data/DISCO_hPBMC_seurat.train.slim.qs")
seu.ref <- NormalizeData(seu.ref)

#### 1. Model building ####
obj.list <- SplitObject(seu.ref, split.by = "project_id")
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
                                        selection.method = "vst",
                                        nfeatures = 2000,
                                        verbose = FALSE)
}
vfeatures <- SelectIntegrationFeatures(obj.list, nfeatures = 2000)
seu.ref[["RNA"]]@var.features <- vfeatures

TIME0 <- Sys.time()
reference <- BuildReference.Symphony(
  seu.ref = seu.ref, npcs = 20, 
  batch.names = c("project_id"), 
  save_uwot_path = "ref_models/03_symphony/PBMC_ref_uwot")
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")

## 保存参考模型
qs::qsave(reference, 'ref_models/03_symphony/PBMC_ref_model.qs')

## 将harmony的结果保存到Seurat对象中
seu.ref[["harmony"]] <- CreateDimReducObject(
  embeddings = t(reference$Z_corr), key = "harmony_", assay = "RNA")
seu.ref[["umap"]] <- CreateDimReducObject(
  embeddings = reference$umap$embedding, key = "UMAP_", assay = "RNA")


#### 2. Reference mapping ####
seu.q <- LoadSeuratSlimData("data/DISCO_hPBMC_seurat.test.slim.qs")
reference <- qs::qread("ref_models/03_symphony/PBMC_ref_model.qs")

TIME0 <- Sys.time()
seu.q <- MapQuery.Symphony(seu.q, reference = reference, assay.q = "RNA")
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

#### 3. Save results ####
SaveResults(celltype.col = "cell_subtype", task.name = "PBMC")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ MFI task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seu.ref <- LoadSeuratSlimData("data/Vento2018.MFI.seurat.train.slim.qs")
seu.ref <- NormalizeData(seu.ref)

#### 1. Model building ####
obj.list <- SplitObject(seu.ref, split.by = "orig.ident")
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
                                        selection.method = "vst",
                                        nfeatures = 2000,
                                        verbose = FALSE)
}
vfeatures <- SelectIntegrationFeatures(obj.list, nfeatures = 2000)
seu.ref[["RNA"]]@var.features <- vfeatures

TIME0 <- Sys.time()
reference <- BuildReference.Symphony(
  seu.ref = seu.ref, npcs = 30, 
  batch.names = c("Method"), 
  save_uwot_path = "ref_models/03_symphony/MFI_ref_uwot")
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")

## 保存参考模型
qs::qsave(reference, 'ref_models/03_symphony/MFI_ref_model.qs')

## 将harmony的结果保存到Seurat对象中
seu.ref[["harmony"]] <- CreateDimReducObject(
  embeddings = t(reference$Z_corr), key = "harmony_", assay = "RNA")
seu.ref[["umap"]] <- CreateDimReducObject(
  embeddings = reference$umap$embedding, key = "UMAP_", assay = "RNA")

DimPlot(seu.ref, group.by = "annotation", label = T) + scale_color_manual(values = ProjectSVR::pals$hfmi)

#### 2. Reference mapping ####
seu.q <- LoadSeuratSlimData("data/Vento2018.MFI.seurat.test.slim.qs")
reference <- qs::qread("ref_models/03_symphony/MFI_ref_model.qs")

TIME0 <- Sys.time()
seu.q <- MapQuery.Symphony(seu.q, reference = reference, assay.q = "RNA")
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

#### 3. Save results ####
SaveResults(celltype.col = "annotation", task.name = "MFI")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ mTCA task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
seu.ref <- LoadSeuratSlimData("data/mTCA.seurat.train.slim.qs")
seu.ref <- NormalizeData(seu.ref)

#### 1. Model building ####
obj.list <- SplitObject(seu.ref, split.by = "GSE_ID")
for (i in 1:length(obj.list)) {
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
                                        selection.method = "vst",
                                        nfeatures = 2000,
                                        verbose = FALSE)
}
vfeatures <- SelectIntegrationFeatures(obj.list, nfeatures = 2000)
seu.ref[["RNA"]]@var.features <- vfeatures

TIME0 <- Sys.time()
reference <- BuildReference.Symphony(
  seu.ref = seu.ref, npcs = 30, 
  batch.names = c("GSE_ID"), 
  save_uwot_path = "ref_models/03_symphony/mTCA_ref_uwot")
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")

## 保存参考模型
qs::qsave(reference, 'ref_models/03_symphony/mTCA_ref_model.qs')

## 将harmony的结果保存到Seurat对象中
seu.ref[["harmony"]] <- CreateDimReducObject(
  embeddings = t(reference$Z_corr), key = "harmony_", assay = "RNA")
seu.ref[["umap"]] <- CreateDimReducObject(
  embeddings = reference$umap$embedding, key = "UMAP_", assay = "RNA")
DimPlot(seu.ref, group.by = "Cell_type_symbol", label = F, raster = F) + scale_color_manual(values = ProjectSVR::pals$mtca)

#### 2. Reference mapping ####
seu.q <- LoadSeuratSlimData("data/mTCA.seurat.test.slim.qs")
# reference <- qs::qread("ref_models/03_symphony/mTCA_ref_model.qs")

TIME0 <- Sys.time()
seu.q <- MapQuery.Symphony(seu.q, reference = reference, assay.q = "RNA")
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

#### 3. Save results ####
SaveResults(celltype.col = "Cell_type_symbol", task.name = "mTCA")
