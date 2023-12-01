library(ProjectSVR)
library(Seurat)
library(tidyverse)
data(ribo.genes)
setwd(here::here())
dir.create("ref_models/02_ProjectSVR")
dir.create("results/02_ProjectSVR")
source("R/IO.R")
source("R/ProjectSVR_utils.R")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ PBMC task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 1 .Model building ####
seu.ref <- LoadSeuratSlimData("data/DISCO_hPBMC_seurat.train.slim.qs")
seu.ref <- NormalizeData(seu.ref)

top.genes <- SelectFeatures(seu.ref)
bg.genes <- do.call(c, top.genes) %>% unique()

TIME0 <- Sys.time()
seu.ref <- ComputeModuleScore(seu.ref, gene.sets = top.genes, bg.genes = bg.genes,
                              method = "UCell", cores = 10)
DefaultAssay(seu.ref) <- "SignatureScore"

gss.mat <- FetchData(seu.ref, vars = rownames(seu.ref))
embeddings.df <- FetchData(seu.ref, vars = paste0("UMAP_", 1:2))
umap.model <- FitEnsembleSVM(feature.mat = gss.mat,
                             emb.mat = embeddings.df,
                             balance.cell.type = F,
                             batch.size = 16000,  # number of subsampled cells for each SVR model
                             n.models = 10,      # number of SVR models trained
                             cores = 10)

meta.data <- FetchData(seu.ref, vars = c(paste0("UMAP_", 1:2), "cell_type", "cell_subtype") )
reference <- CreateReference(umap.model = umap.model,
                             gene.sets = top.genes,
                             bg.genes = bg.genes,
                             meta.data = meta.data,
                             gss.method = "UCell",
                             colors = pals$disco_blood)
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")
## save reference model
qs::qsave(reference, "ref_models/02_ProjectSVR/PBMC_ref_model.qs")

#### 2. Reference mapping ####
seu.q <- LoadSeuratSlimData("data/DISCO_hPBMC_seurat.test.slim.qs")
TIME0 <- Sys.time()
seu.q <- ProjectSVR::MapQuery(seu.q, reference = reference, add.map.qual = T, ncores = 10)
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

#### 3. Save results ####
SaveResults(celltype.col = "cell_subtype", task.name = "PBMC")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ MFI task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 1 .Model building ####
seu.ref <- LoadSeuratSlimData("data/Vento2018.MFI.seurat.train.slim.qs")
seu.ref <- NormalizeData(seu.ref)

top.genes <- SelectFeatures(seu.ref)
bg.genes <- do.call(c, top.genes) %>% unique()

TIME0 <- Sys.time()
seu.ref <- ComputeModuleScore(seu.ref, gene.sets = top.genes, bg.genes = bg.genes, 
                              method = "UCell", cores = 10)
DefaultAssay(seu.ref) <- "SignatureScore"

gss.mat <- FetchData(seu.ref, vars = rownames(seu.ref))
embeddings.df <- FetchData(seu.ref, vars = paste0("UMAP_", 1:2))
umap.model <- FitEnsembleSVM(feature.mat = gss.mat,
                             emb.mat = embeddings.df, 
                             balance.cell.type = F,
                             batch.size = 16000,  # number of subsampled cells for each SVR model 
                             n.models = 10,      # number of SVR models trained
                             cores = 10)

meta.data <- FetchData(seu.ref, vars = c(paste0("UMAP_", 1:2), "annotation") )
reference <- CreateReference(umap.model = umap.model, 
                             gene.sets = top.genes, 
                             bg.genes = bg.genes,
                             meta.data = meta.data, 
                             gss.method = "UCell",
                             colors = pals$hfmi)
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")
## save reference model
qs::qsave(reference, "ref_models/02_ProjectSVR/MFI_ref_model.qs")

#### 2. Reference mapping ####
seu.q <- LoadSeuratSlimData("data/Vento2018.MFI.seurat.test.slim.qs")
TIME0 <- Sys.time()
seu.q <- ProjectSVR::MapQuery(seu.q, reference = reference, add.map.qual = T, ncores = 10)
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

#### 3. Save results ####
SaveResults(celltype.col = "annotation", task.name = "MFI")

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ mTCA task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 1 .Model building ####
seu.ref <- LoadSeuratSlimData("data/mTCA.seurat.train.slim.qs")

top.genes <- read.table("../ProjectSVR/vignettes/mTCA/K70/top100_genes.k_70.dt_0_2.txt", header = T)
names(top.genes) <- paste0("feature.", 1:length(top.genes))

TIME0 <- Sys.time()
seu.ref <- ComputeModuleScore(seu.ref, gene.sets = top.genes, 
                              method = "AUCell", cores = 10)
DefaultAssay(seu.ref) <- "SignatureScore"

gss.mat <- FetchData(seu.ref, vars = rownames(seu.ref))
embeddings.df <- FetchData(seu.ref, vars = paste0("UMAP_", 1:2))
umap.model <- FitEnsembleSVM(feature.mat = gss.mat,
                             emb.mat = embeddings.df, 
                             balance.cell.type = F,
                             batch.size = 16000,  # number of subsampled cells for each SVR model 
                             n.models = 10,      # number of SVR models trained
                             cores = 10)

meta.data <- FetchData(seu.ref, vars = c(paste0("UMAP_", 1:2), "Cell_type_symbol") )
reference <- CreateReference(umap.model = umap.model, 
                             gene.sets = top.genes, 
                             meta.data = meta.data, 
                             gss.method = "AUCell",
                             colors = pals$mtca)
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")
## save reference model
qs::qsave(reference, "ref_models/02_ProjectSVR/mTCA_ref_model.qs")

#### 2. Reference mapping ####
seu.q <- LoadSeuratSlimData("data/mTCA.seurat.test.slim.qs")
TIME0 <- Sys.time()
seu.q <- ProjectSVR::MapQuery(seu.q, reference = reference, add.map.qual = T, ncores = 10)
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

#### 3. Save results ####
SaveResults(celltype.col = "Cell_type_symbol", task.name = "mTCA")


#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ CD4TIT task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 1 .Model building ####
seu.ref <- LoadSeuratSlimData("data/ZhengLiangtao.CD4.seurat.train.slim.qs")

top.genes <- seu.ref@misc$markers
top.genes <- lapply(top.genes, function(xx) head(xx, 20))

TIME0 <- Sys.time()
seu.ref <- ComputeModuleScore(seu.ref, gene.sets = top.genes, 
                              method = "UCell", cores = 10)
DefaultAssay(seu.ref) <- "SignatureScore"

gss.mat <- FetchData(seu.ref, vars = rownames(seu.ref))
embeddings.df <- FetchData(seu.ref, vars = paste0("UMAP_", 1:2))
umap.model <- FitEnsembleSVM(feature.mat = gss.mat,
                             emb.mat = embeddings.df, 
                             balance.cell.type = F,
                             batch.size = 16000,  # number of subsampled cells for each SVR model 
                             n.models = 10,      # number of SVR models trained
                             cores = 10)

meta.data <- seu.ref@misc$data.refplot$meta.data
bg.genes <- rownames(seu.ref[["RNA"]])
reference <- CreateReference(umap.model = umap.model, 
                             gene.sets = top.genes, 
                             bg.genes = bg.genes,
                             meta.data = meta.data, 
                             gss.method = "UCell", 
                             colors = pals$`pan-cancer_cd4t`)
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")
## save reference model
qs::qsave(reference, "ref_models/02_ProjectSVR/CD4TIT_ref_model.qs")
# reference <- qs::qread("ref_models/02_ProjectSVR/CD4TIT_ref_model.qs")

#### 2. Reference mapping ####
seu.q <- LoadSeuratSlimData("data/ZhengLiangtao.CD4.seurat.test.slim.qs")
TIME0 <- Sys.time()
seu.q <- ProjectSVR::MapQuery(seu.q, reference = reference, add.map.qual = T, ncores = 10)
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

# p1 <- DimPlot(seu.q, reduction = "umap", group.by = "meta.cluster.name") + scale_color_manual(values = pals$`pan-cancer_cd4t`)
# p2 <- DimPlot(seu.q, reduction = "ref.umap", group.by = "meta.cluster.name") + scale_color_manual(values = pals$`pan-cancer_cd4t`)
# p1 + NoLegend() + p2

#### 3. Save results ####
task.name = "CD4TIT"

ref.emb <- seu.ref@misc$data.refplot$meta.data[, paste0("UMAP_", 1:2)]
query.emb <- seu.q[["ref.umap"]]@cell.embeddings
ref.cellmeta <- seu.ref@misc$data.refplot$meta.data["cluster.name"]
ref.cellmeta$group <- "reference"
query.cellmeta <- seu.q@meta.data["meta.cluster.name"]
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
## save
write.csv(ref.emb,        glue("results/02_ProjectSVR/{task.name}_ref_embeddings.csv"))
write.csv(ref.cellmeta,   glue("results/02_ProjectSVR/{task.name}_ref_cellmeta.csv"))
write.csv(query.emb,      glue("results/02_ProjectSVR/{task.name}_query_embeddings.csv"))
write.csv(query.cellmeta, glue("results/02_ProjectSVR/{task.name}_query_cellmeta.csv"))
write.csv(runtime,        glue("results/02_ProjectSVR/{task.name}_runtime.csv"))

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### ++ CD8TIT task ++ ####
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#### 1 .Model building ####
seu.ref <- LoadSeuratSlimData("data/ZhengLiangtao.CD8.seurat.train.slim.qs")

top.genes <- seu.ref@misc$markers
top.genes <- lapply(top.genes, function(xx) head(xx, 20))

TIME0 <- Sys.time()
seu.ref <- ComputeModuleScore(seu.ref, gene.sets = top.genes, 
                              method = "UCell", cores = 10)
DefaultAssay(seu.ref) <- "SignatureScore"

gss.mat <- FetchData(seu.ref, vars = rownames(seu.ref))
embeddings.df <- FetchData(seu.ref, vars = paste0("UMAP_", 1:2))
umap.model <- FitEnsembleSVM(feature.mat = gss.mat,
                             emb.mat = embeddings.df, 
                             balance.cell.type = F,
                             batch.size = 16000,  # number of subsampled cells for each SVR model 
                             n.models = 10,      # number of SVR models trained
                             cores = 10)

meta.data <- seu.ref@misc$data.refplot$meta.data
bg.genes <- rownames(seu.ref[["RNA"]])
reference <- CreateReference(umap.model = umap.model, 
                             gene.sets = top.genes, 
                             bg.genes = bg.genes,
                             meta.data = meta.data, 
                             gss.method = "UCell", 
                             colors = pals$`pan-cancer_cd8t`)
TIME1 <- Sys.time()
runtime1 <- difftime(TIME1, TIME0, units = "secs")
## save reference model
qs::qsave(reference, "ref_models/02_ProjectSVR/CD8TIT_ref_model.qs")

#### 2. Reference mapping ####
seu.q <- LoadSeuratSlimData("data/ZhengLiangtao.CD8.seurat.test.slim.qs")
TIME0 <- Sys.time()
seu.q <- ProjectSVR::MapQuery(seu.q, reference = reference, add.map.qual = T, ncores = 10)
TIME1 <- Sys.time()
runtime2 <- difftime(TIME1, TIME0, units = "secs")

#### 3. Save results ####
task.name = "CD8TIT"

ref.emb <- seu.ref@misc$data.refplot$meta.data[, paste0("UMAP_", 1:2)]
query.emb <- seu.q[["ref.umap"]]@cell.embeddings
ref.cellmeta <- seu.ref@misc$data.refplot$meta.data["cluster.name"]
ref.cellmeta$group <- "reference"
query.cellmeta <- seu.q@meta.data["meta.cluster.name"]
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
## save
write.csv(ref.emb,        glue("results/02_ProjectSVR/{task.name}_ref_embeddings.csv"))
write.csv(ref.cellmeta,   glue("results/02_ProjectSVR/{task.name}_ref_cellmeta.csv"))
write.csv(query.emb,      glue("results/02_ProjectSVR/{task.name}_query_embeddings.csv"))
write.csv(query.cellmeta, glue("results/02_ProjectSVR/{task.name}_query_cellmeta.csv"))
write.csv(runtime,        glue("results/02_ProjectSVR/{task.name}_runtime.csv"))

