rowSDs <- function(A) {
  apply(A, 1, stats::sd)
}

CalMeanSDs <- function(seu.ref, assay = "RNA") {
  vfeatures <- seu.ref[[assay]]@var.features
  ref_exp <- seu.ref[[assay]]@data[vfeatures, ]
  ## Calculate the mean and sd for each gene.
  vfeatures_means_sds <- tibble(symbol = vfeatures, mean = Matrix::rowMeans(ref_exp))
  vfeatures_means_sds$stddev <- rowSDs(ref_exp)
  return(vfeatures_means_sds)
}

scaleDataWithStats <- function(A, mean_vec, sd_vec) {
  nc <- ncol(A)
  U <- matrix(rep(mean_vec, each = nc), ncol = nc, byrow = T)
  V <- matrix(rep(sd_vec, each = nc), ncol = nc, byrow = T)
  (A - U) / V
}

DoScale <- function(seu.ref, assay = "RNA", vfeatures_means_sds) {
  data.scaled <- seu.ref[["RNA"]]@data[vfeatures_means_sds$symbol, ]
  data.scaled <- scaleDataWithStats(A = data.scaled,
                                    mean_vec = vfeatures_means_sds$mean,
                                    sd_vec = vfeatures_means_sds$stddev)
  seu.ref[["RNA"]]@scale.data <- as.matrix(data.scaled)
  return(seu.ref)
}

BuildReference.Symphony <- function(seu.ref, npcs = 20, batch.names = c("orig.ident"), save_uwot_path = "") {
  vfeatures_means_sds <- CalMeanSDs(seu.ref)
  seu.ref <- DoScale(seu.ref, vfeatures_means_sds = vfeatures_means_sds)
  seu.ref <- RunPCA(seu.ref, npcs = npcs, verbose = FALSE)
  
  ## Harmony (return with full harmony object)
  set.seed(0)
  Z_pca_ref <- seu.ref[["pca"]]@cell.embeddings
  ref_metadata <- seu.ref@meta.data
  ref_harmObj <- harmony::HarmonyMatrix(
    data_mat = Z_pca_ref,     ## PCA embedding matrix of cells
    meta_data = ref_metadata, ## data.frame with cell labels
    vars_use = batch.names,   ## variable to integrate out
    nclust = NULL,            ## number of clusters in Harmony model
    max.iter.harmony = 10,
    return_object = TRUE,     ## return the full Harmony model object
    do_pca = FALSE,           ## don't recompute PCs
  )
  
  ## 将harmony对象转变为symphony对象，用于参考映射
  vfeature_loadings <- seu.ref[["pca"]]@feature.loadings
  # save_uwot_path <- file.path(getwd(), "ref_models/03_symphony/MFI_ref_uwot")
  save_uwot_path <- file.path(getwd(), save_uwot_path)
  reference <- symphony::buildReferenceFromHarmonyObj(
    ref_harmObj,            # output object from HarmonyMatrix()
    ref_metadata,           # reference cell metadata
    vfeatures_means_sds,    # gene names, means, and std devs for scaling
    vfeature_loadings,      # genes x PCs matrix
    verbose = TRUE,         # verbose output
    do_umap = TRUE,         # Set to TRUE only when UMAP model was saved for reference
    save_uwot_path = save_uwot_path) # save_uwot_path should be full path.
  return(reference)
}

MapQuery.Symphony <- function(seu.q, reference, assay.q = "RNA") {
  vfeatures <- intersect(reference$vargenes$symbol, rownames(seu.q))
  query_exp <- seu.q[[assay.q]]@data[vfeatures, ]
  query_metadata <- seu.q@meta.data
  query <- mapQuery(query_exp,             # query gene expression (genes x cells)
                    query_metadata,        # query metadata (cells x attributes)
                    reference,             # Symphony reference object
                    vars = NULL,           # Query batch variables to harmonize over (NULL treats query as one batch)
                    do_normalize = FALSE,  # perform log(CP10k) normalization on query (set to FALSE if already normalized)
                    do_umap = TRUE)        # project query cells into reference UMAP
  colnames(query$umap) <- paste0("refUmap_", 1:2)
  seu.q[["ref.umap"]] <- CreateDimReducObject(embeddings = query$umap, key = "refUmap_", assay = assay.q)
  return(seu.q)
}

LabelTransfer.Symphony <- function(seu.q, reference, reduction.q = "ref.umap", ref.label.col = "cluster.name", k = 10) {
  ref.cellmeta <- reference$meta_data
  ref.emb <- reference$umap$embedding
  ref.labels <- ref.cellmeta[[ref.label.col]]
  names(ref.labels) <- rownames(ref.cellmeta)
  
  query.emb <- seu.q[[reduction.q]]@cell.embeddings
  knn.pred <- ProjectSVR::KnnLabelTransfer(query.emb = query.emb, ref.emb = ref.emb, ref.labels = ref.labels, k = k)
  seu.q$knn.pred.celltype <- knn.pred$labels
  seu.q$knn.pred.votes <- knn.pred$votes
  seu.q$knn.pred.perc <- knn.pred$perc
  ref.celltype.levels <- levels(ref.cellmeta[[ref.label.col]])
  if (!is.null(ref.celltype.levels)) {
    seu.q$knn.pred.celltype <- factor(seu.q$knn.pred.celltype, 
                                      levels = ref.celltype.levels)
  }
  return(seu.q)
}

SaveResults <- function(celltype.col, task.name) {
  ref.emb <- seu.ref[["umap"]]@cell.embeddings
  query.emb <- seu.q[["ref.umap"]]@cell.embeddings
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
  ## save
  write.csv(ref.emb,        glue("results/03_symphony/{task.name}_ref_embeddings.csv"))
  write.csv(ref.cellmeta,   glue("results/03_symphony/{task.name}_ref_cellmeta.csv"))
  write.csv(query.emb,      glue("results/03_symphony/{task.name}_query_embeddings.csv"))
  write.csv(query.cellmeta, glue("results/03_symphony/{task.name}_query_cellmeta.csv"))
  write.csv(runtime,        glue("results/03_symphony/{task.name}_runtime.csv"))
}