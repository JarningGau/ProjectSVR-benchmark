LabelTransfer.projecTILs <- function(seu.q, seu.ref, reduction.ref = "umap", reduction.q = "umap", ref.label.col = "cluster.name", k = 10) {
  ref.emb <- seu.ref[[reduction.ref]]@cell.embeddings
  ref.labels <- seu.ref[[ref.label.col, drop=T]]
  
  query.emb <- seu.q[[reduction.q]]@cell.embeddings
  knn.pred <- ProjectSVR::KnnLabelTransfer(query.emb = query.emb, ref.emb = ref.emb, ref.labels = ref.labels, k = k)
  seu.q$knn.pred.celltype <- knn.pred$labels
  seu.q$knn.pred.votes <- knn.pred$votes
  seu.q$knn.pred.perc <- knn.pred$perc
  ref.celltype.levels <- levels(ref.labels)
  if (!is.null(ref.celltype.levels)) {
    seu.q$knn.pred.celltype <- factor(seu.q$knn.pred.celltype, 
                                      levels = ref.celltype.levels)
  }
  return(seu.q)
}