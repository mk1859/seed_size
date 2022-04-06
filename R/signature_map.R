# function for signature plots with violin insets
signature_map <- function (seurat_obj, vis_obj = NULL, signature, order, column) {
  require (VISION)
  require (ggplot2)
  require (Seurat)
  
  # export data from Seurat object
  plot <- cbind (as.data.frame (Embeddings(object = seurat_obj, 
                                           reduction = "pca"))  [,1:2],seurat_obj@meta.data)
  if (!is.null(vis_obj)) {
    plot <- cbind (plot, vis_obj@SigScores)
  }
  
  # set order of treatments
  plot$timepoint <- factor(plot$timepoint, levels = order)
  
  # PCA plot
  g1 <- ggplot(plot, aes(x=PC_1, y= PC_2, color = .data[[signature]])) +
    geom_point (size = 2) + 
    scale_colour_viridis() +
    theme_classic () +
    theme(legend.title=element_blank(),
          axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks=element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) 
  
  # violin plot for treatments
  g2 <- ggplot(plot, aes(x=.data[[column]], y= .data[[signature]], fill = .data[[column]], color = .data[[column]])) +
    geom_violin () + 
    scale_fill_tableau("Green-Orange-Teal") +
    scale_color_tableau("Green-Orange-Teal") +
    theme_classic() +
    theme(legend.position = "none",
          strip.background = element_blank(),
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(),
          axis.text.x = element_blank())
  
  # combine plots
  g <- ggdraw() +
    draw_plot(g1) +
    draw_plot(g2, x = 0.75, y = .7, width = .25, height = .3)
  
  return (g)
}
