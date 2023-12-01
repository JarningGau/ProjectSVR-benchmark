`%notin%` <- Negate(`%in%`)

SelectFeatures <- function(seu, downsample.size = 200, top.n.genes = 25, n.cores = 10) {
  ## To accelerate the calculation, we downsampled 200 cells from each cluster
  seu.ref.ds <- subset(seu, downsample = downsample.size)
  seu.ref.ds <- NormalizeData(seu.ref.ds)
  
  ## Parallel calculation of the cell markers.
  all.markers <- mcFindAllMarkers(seu.ref.ds, do.flatten = F, n.cores = n.cores)
  
  top.genes <- lapply(all.markers, function(xx){
    yy <- subset(xx, p_val_adj < 1e-6 & avg_log2FC > log2(1.5))
    yy <- subset(yy, Gene.name.uniq %notin% ribo.genes)
    yy <- yy[!grepl("^MT-", yy$Gene.name.uniq), ]
    yy <- yy[!grepl("^mt-", yy$Gene.name.uniq), ]
    head(yy$Gene.name.uniq, top.n.genes)
  })
  names(top.genes) <- paste0("feature.", 1:length(top.genes))
  return(top.genes)
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
  write.csv(ref.emb,        glue("results/02_ProjectSVR/{task.name}_ref_embeddings.csv"))
  write.csv(ref.cellmeta,   glue("results/02_ProjectSVR/{task.name}_ref_cellmeta.csv"))
  write.csv(query.emb,      glue("results/02_ProjectSVR/{task.name}_query_embeddings.csv"))
  write.csv(query.cellmeta, glue("results/02_ProjectSVR/{task.name}_query_cellmeta.csv"))
  write.csv(runtime,        glue("results/02_ProjectSVR/{task.name}_runtime.csv"))
}