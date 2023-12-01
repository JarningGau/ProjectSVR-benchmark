library(Seurat)
library(tidyverse)

seu <- qs::qread("~/dev/ProjectSVR-data/reference_atlas/DISCO_hPBMCs.seurat.slim.qs")
seu[["RNA"]]@counts <- seu[["RNA"]]@data
sceasy::convertFormat(seu, from = "seurat", to = "anndata", drop_single_values = F,
                      main_layer = "counts", outFile = "../ProjectDNN-explore/data/DISCO_hPBMC.h5ad")

write.table(table(seu$project_id, seu$cell_subtype), "tmp/DISCO_hPBMCs_cell_stat.txt", sep = "\t")

test.query <- "GSE175499"

seu.test <- subset(seu, project_id == test.query)
seu.train <- subset(seu, project_id == test.query, invert = T)

sceasy::convertFormat(seu.test, from = "seurat", to = "anndata", drop_single_values = F,
                      main_layer = "counts", outFile = "data/DISCO_hPBMC.test.counts.h5ad")
sceasy::convertFormat(seu.train, from = "seurat", to = "anndata", drop_single_values = F,
                      main_layer = "counts", outFile = "data/DISCO_hPBMC.train.counts.h5ad")


seu.test[["RNA"]]@counts <- matrix()
seu.train[["RNA"]]@counts <- matrix()

qs::qsave(seu.test, "data/DISCO_hPBMC_seurat.test.slim.qs")
qs::qsave(seu.train, "data/DISCO_hPBMC_seurat.train.slim.qs")

## save golden results
seu.query <- qs::qread("data/DISCO_hPBMC_seurat.test.slim.qs")
seu.ref <- qs::qread("data/DISCO_hPBMC_seurat.train.slim.qs")

ref.emb <- seu.ref[["umap"]]@cell.embeddings
query.emb <- seu.q[["umap"]]@cell.embeddings
ref.cellmeta <- seu.ref@meta.data["cell_subtype"]
ref.cellmeta$group <- "reference"
query.cellmeta <- seu.q@meta.data["cell_subtype"]
query.cellmeta$group <- "query"

colnames(ref.emb) <- paste0("Dim_", 1:2)
colnames(query.emb) <- paste0("Dim_", 1:2)
colnames(ref.cellmeta) <- c("label", "group")
colnames(query.cellmeta) <- c("label", "group")

dir.create("results/01_Golden")
write.csv(ref.emb, "results/01_Golden/PBMC_ref_embeddings.csv")
write.csv(ref.cellmeta, "results/01_Golden/PBMC_ref_cellmeta.csv")
write.csv(query.emb, "results/01_Golden/PBMC_query_embeddings.csv")
write.csv(query.cellmeta, "results/01_Golden/PBMC_query_cellmeta.csv")

runtime <- data.frame(
  model.building = 0,
  reference.mapping = 0
)
write.csv(runtime, "results/01_Golden/PBMC_runtime.csv")

##///////////////////////////////

seu <- qs::qread("~/dev/ProjectSVR-data/reference_atlas/Vento2018.MFI.seurat.slim.qs")
seu[["RNA"]]@counts <- seu[["RNA"]]@data

write.table(table(seu$orig.ident, seu$annotation), "tmp/Vento2018_MFI_cell_stat.txt", sep = "\t")

test.query <- c("FCA7196224", "FCA7196225", "FCA7511884")

seu.test <- subset(seu, orig.ident %in% test.query)
seu.train <- subset(seu, orig.ident %in% test.query, invert = T)

sceasy::convertFormat(seu.test, from = "seurat", to = "anndata", drop_single_values = F,
                      main_layer = "counts", outFile = "data/Vento2018.MFI.test.counts.h5ad")
sceasy::convertFormat(seu.train, from = "seurat", to = "anndata", drop_single_values = F,
                      main_layer = "counts", outFile = "data/Vento2018.MFI.train.counts.h5ad")

seu.test[["RNA"]]@counts <- matrix()
seu.train[["RNA"]]@counts <- matrix()

qs::qsave(seu.test, "data/Vento2018.MFI.seurat.test.slim.qs")
qs::qsave(seu.train, "data/Vento2018.MFI.seurat.train.slim.qs")

## save golden results
seu.q <- qs::qread("data/Vento2018.MFI.seurat.test.slim.qs")
seu.ref <- qs::qread("data/Vento2018.MFI.seurat.train.slim.qs")

ref.emb <- seu.ref[["umap"]]@cell.embeddings
query.emb <- seu.q[["umap"]]@cell.embeddings
ref.cellmeta <- seu.ref@meta.data["annotation"]
ref.cellmeta$group <- "reference"
query.cellmeta <- seu.q@meta.data["annotation"]
query.cellmeta$group <- "query"

colnames(ref.emb) <- paste0("Dim_", 1:2)
colnames(query.emb) <- paste0("Dim_", 1:2)
colnames(ref.cellmeta) <- c("label", "group")
colnames(query.cellmeta) <- c("label", "group")

write.csv(ref.emb, "results/01_Golden/MFI_ref_embeddings.csv")
write.csv(ref.cellmeta, "results/01_Golden/MFI_ref_cellmeta.csv")
write.csv(query.emb, "results/01_Golden/MFI_query_embeddings.csv")
write.csv(query.cellmeta, "results/01_Golden/MFI_query_cellmeta.csv")

runtime <- data.frame(
  model.building = 0,
  reference.mapping = 0
)
write.csv(runtime, "results/01_Golden/MFI_runtime.csv")


##///////////////////////////////

seu <- qs::qread("~/dev/ProjectSVR-data/reference_atlas/mTCA.seurat.slim.qs")
seu[["RNA"]]@counts <- seu[["RNA"]]@data

write.table(table(seu$GSM_ID, seu$Cell_type_symbol), "tmp/mTCA_cell_stat.txt", sep = "\t")

test.query <- c("GSM2928505", "do17825", "GSM3744444")

seu.test <- subset(seu, orig.ident %in% test.query)
seu.train <- subset(seu, orig.ident %in% test.query, invert = T)

sceasy::convertFormat(seu.test, from = "seurat", to = "anndata", drop_single_values = F,
                      main_layer = "counts", outFile = "data/mTCA.test.counts.h5ad")
sceasy::convertFormat(seu.train, from = "seurat", to = "anndata", drop_single_values = F,
                      main_layer = "counts", outFile = "data/mTCA.train.counts.h5ad")

seu.test[["RNA"]]@counts <- matrix()
seu.train[["RNA"]]@counts <- matrix()

qs::qsave(seu.test, "data/mTCA.seurat.test.slim.qs")
qs::qsave(seu.train, "data/mTCA.seurat.train.slim.qs")

## save golden results
seu.q <- qs::qread("data/mTCA.seurat.test.slim.qs")
seu.ref <- qs::qread("data/mTCA.seurat.train.slim.qs")

ref.emb <- seu.ref[["umap"]]@cell.embeddings
query.emb <- seu.q[["umap"]]@cell.embeddings
ref.cellmeta <- seu.ref@meta.data["Cell_type_symbol"]
ref.cellmeta$group <- "reference"
query.cellmeta <- seu.q@meta.data["Cell_type_symbol"]
query.cellmeta$group <- "query"

colnames(ref.emb) <- paste0("Dim_", 1:2)
colnames(query.emb) <- paste0("Dim_", 1:2)
colnames(ref.cellmeta) <- c("label", "group")
colnames(query.cellmeta) <- c("label", "group")

write.csv(ref.emb, "results/01_Golden/mTCA_ref_embeddings.csv")
write.csv(ref.cellmeta, "results/01_Golden/mTCA_ref_cellmeta.csv")
write.csv(query.emb, "results/01_Golden/mTCA_query_embeddings.csv")
write.csv(query.cellmeta, "results/01_Golden/mTCA_query_cellmeta.csv")

runtime <- data.frame(
  model.building = 0,
  reference.mapping = 0
)
write.csv(runtime, "results/01_Golden/mTCA_runtime.csv")


##///////////////////////////////

seu <- qs::qread("~/dev/ProjectSVR-data/reference_atlas/ZhengLiangtao.CD4.seurat.slim.qs")
seu[["RNA"]]@counts <- seu[["RNA"]]@data

write.table(table(seu$dataset, seu$meta.cluster.name), "tmp/ZhengLiangtao_CD4_cell_stat.txt", sep = "\t")

test.query <- c("BC.Elham2018.10X")

seu.test <- subset(seu, dataset %in% test.query)
seu.train <- subset(seu, dataset %in% test.query, invert = T)

seu.test[["RNA"]]@counts <- matrix()
seu.train[["RNA"]]@counts <- matrix()

qs::qsave(seu.test, "data/ZhengLiangtao.CD4.seurat.test.slim.qs")
qs::qsave(seu.train, "data/ZhengLiangtao.CD4.seurat.train.slim.qs")

##///////////////////////////////

seu <- qs::qread("~/dev/ProjectSVR-data/reference_atlas/ZhengLiangtao.CD8.seurat.slim.qs")
seu[["RNA"]]@counts <- seu[["RNA"]]@data

write.table(table(seu$dataset, seu$meta.cluster.name), "tmp/ZhengLiangtao_CD8_cell_stat.txt", sep = "\t")

test.query <- c("BC.Elham2018.10X")

seu.test <- subset(seu, dataset %in% test.query)
seu.train <- subset(seu, dataset %in% test.query, invert = T)

seu.test[["RNA"]]@counts <- matrix()
seu.train[["RNA"]]@counts <- matrix()

qs::qsave(seu.test, "data/ZhengLiangtao.CD8.seurat.test.slim.qs")
qs::qsave(seu.train, "data/ZhengLiangtao.CD8.seurat.train.slim.qs")
