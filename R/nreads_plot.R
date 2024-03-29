# function to create a multipanel diagram, it takes order of treatments as input
# require seed naming: treatment.library(01-99).batch(01-99)

nreads_plot <- function(matrix, order, tableu) {
  require (cowplot)
  require (patchwork)
  # data frame with statistics per seed
  seed_attr <- data.frame(n_reads = colSums(matrix), # number of reads
                        n_gene = colSums(matrix > 0), # number of genes
                        # extract information about time point from seed name
                        timepoint = as.factor(gsub('.{0,6}$', '', colnames(matrix))))
  
  # set order of treatments
  seed_attr$timepoint <- factor(seed_attr$timepoint, 
                              levels = order)
  
  # plot for number of reads and number of genes
  g1 <- ggplot(seed_attr, aes(x=n_reads, y=n_gene, color = timepoint)) +
   geom_point(size = 1) + 
   scale_color_tableau(tableu) +
   theme_classic() +
   theme(legend.position = "none")
  
  # boxplot for number of genes
  g2 <- ggplot(seed_attr, aes(x=timepoint, y=n_gene, color = timepoint)) +
   geom_boxplot() + 
   scale_color_tableau(tableu) +
   theme_classic() +
   theme(axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank()) +
   theme(legend.position = "none")

  # boxplot for number of reads
  g3 <- ggplot(seed_attr, aes(x=timepoint, y=n_reads, color = timepoint)) +
   geom_boxplot() + 
   scale_color_tableau(tableu) +
   theme_classic() +
   theme(axis.line=element_blank(),
         axis.text.x=element_blank(),
         axis.text.y=element_blank(),
         axis.ticks=element_blank(),
         axis.title.x=element_blank(),
         axis.title.y=element_blank()) +
   coord_flip() +
   theme(legend.position = "none")
  
  # spacer plot 
  s <- plot_spacer()
  
  # create multipanel plot
  g <- plot_grid(g3, s, g1, g2, ncol = 2, align = "hv", 
               rel_heights = c(1,2), rel_widths = c(2,0.64))
  return (g)
}
