library(glue)

LoadSeuratSlimData <- function(qs.file) {
  seu <- qs::qread(qs.file)
  seu[["RNA"]]@counts <- seu[["RNA"]]@data
  return(seu)
}