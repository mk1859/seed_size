# PCA plot with coloring of continuous values

pca_continuous <- function(seurat_obj, column, excluded = NULL) {
  require (Seurat)
  require (ggplot2)
  require (viridis)
  
  # export data from the Seurat object
  plot <- cbind (as.data.frame (Embeddings(object = seurat_obj, reduction = "pca"))  [,1:2], seurat_obj@meta.data) 
  
  # exclude treatment if necessary
  if (!is.null(excluded )){
    plot <-  plot [-grep( excluded, rownames(plot)),]
  }
  
  g <- ggplot(plot, aes(x=PC_1, y= PC_2, color = .data[[column]])) +
            geom_point (size = 2) + 
            scale_colour_viridis() +
            theme_classic () +
            theme(axis.line=element_blank(),
                  axis.text=element_blank(),
                  axis.ticks=element_blank(),
                  axis.title=element_blank()) 
  
  return (g)
}
