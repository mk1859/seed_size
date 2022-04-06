background_plot <- function(matrix, order, background, tableu = "Green-Orange-Teal") {
 # data frame with statistics per seed
  seed_attr <- data.frame(n_reads = colSums(matrix), # number of reads
                          background = background, # fraction of background
                        # extract information about time point from seed name
                        timepoint = as.factor(gsub('.{0,6}$', '', colnames(matrix))))
  
  # set order of treatments
  seed_attr$timepoint <- factor(seed_attr$timepoint, 
                              levels = order)
  
  # boxplot for  fraction of background
  g <- ggplot(seed_attr, aes(x=timepoint, y=background, color = timepoint)) +
   geom_boxplot() + 
   scale_color_tableau(tableu) +
   theme_classic() 

  return (g)
}
