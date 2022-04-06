# signature vs signature plot
sig_vs_sig <- function (seurat_obj, signature1, signature2, vis_obj = NULL, order) {
  # export data from Seurat object
  plot <- seurat_obj@meta.data
  
  if (!is.null(vis_obj)) {
    plot <- cbind (plot, vis_obj@SigScores)
  }
  
  # set order of treatments
  plot$timepoint <- factor(plot$timepoint, 
                           levels = order)
  
  # exclude some time points if necessary
  if (!is.null(excluded )){
    plot <-  plot [-grep( excluded, rownames(plot)),]
  }
  
  g <- ggplot(plot, aes(x=.data[[signature1]], y= .data[[signature2]], color = timepoint)) +
    geom_point(size = 1.5) +
    theme_classic() +
    scale_color_tableau("Green-Orange-Teal")
  
  return (g)
  
}
