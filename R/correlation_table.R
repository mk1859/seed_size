# function to calculate the correlation of gene's read count to the fraction of background reads
# it takes a matrix of counts, normalizes sequencing depth per seed using Seurat 
# then it extracts information about treatment from seed name 
# finally calculates Pearson correlation to the fraction of background reads for seeds in each treatment and all seeds combined

correlation_table <- function(matrix, background) {
  require (Seurat)
  
  # input for normalization
  seed_attr <- data.frame(background = background,                                       # fraction of background reads
                          timepoint = as.factor(gsub('.{0,6}$', '', colnames(matrix))))  # extract information about treatment from seed name
                          
  # normalize using Seurat
  norm_seeds <- NormalizeData(CreateSeuratObject(counts = matrix, meta.data = seed_attr), verbose = FALSE)
  
  # export normalized data
  norm_seeds <- as.matrix(GetAssayData(object = norm_seeds, slot = "data"))
  
  # correlation of gene's read count with background fraction for each treatment
  timepoints <- levels(seed_attr$timepoint)
  
  background_corr <- lapply (timepoints, function (x) { 
                            apply (norm_seeds [,which(seed_attr$timepoint == x)],1, 
                            function(z) cor.test(z,seed_attr$background
                            [which(seed_attr$timepoint == x)])$estimate)
  })
  
  # correlation of gene's read count with background fraction when all seeds combined
  background_corr$all <- apply (norm_seeds,1, function(x) {
                                cor.test(x,seed_attr$background)$estimate})
  
  # data frame from the list
  background_corr <- t(matrix(unlist(background_corr), nrow=length(background_corr), byrow=TRUE))
  colnames (background_corr) <- c (timepoints, "all")
  rownames(background_corr) <- rownames (matrix)
                           
  return (background_corr)
}
