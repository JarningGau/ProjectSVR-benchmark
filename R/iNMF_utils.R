Seurat2LigerRef <- function(seu, split.by, nfeatures=2000) {
  seu <- NormalizeData(seu)
  obj.list <- SplitObject(seu, split.by = split.by)
  for (i in 1:length(obj.list)) {
    obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]],
                                          selection.method = "vst",
                                          nfeatures = 2000,
                                          verbose = FALSE)
  }
  vfeatures <- SelectIntegrationFeatures(obj.list, nfeatures = 2000)
  seu <- ScaleData(seu, features = vfeatures, split.by = split.by, do.center = FALSE)
  obj.list <- SplitObject(seu, split.by = split.by)
  liger <- createLiger(lapply(obj.list, function(xx) xx[["RNA"]]@counts), remove.missing = FALSE)
  liger@var.genes <- vfeatures
  liger@scale.data <- lapply(obj.list, function(xx) t(xx[["RNA"]]@scale.data))
  return(liger)
}


Seurat2LigerQuery <- function(seu, vfeatures) {
  seu <- NormalizeData(seu)
  seu <- ScaleData(seu, features = vfeatures, do.center = FALSE)
  liger <- createLiger(list(query = seu[["RNA"]]@counts), remove.missing = FALSE)
  liger@var.genes <- vfeatures
  liger@scale.data <- list(query = t(seu[["RNA"]]@scale.data))
  return(liger)
}
