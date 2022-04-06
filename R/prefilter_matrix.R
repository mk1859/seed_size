# function to prefilter genes and seeds, it uses Araport data frame with gene type information

prefilter_matrix <- function(matrix, mean_exp, n_reads) {
  
  filt_mat <- matrix [-c (grep ("__", rownames (matrix)),                                 # remove summary statistics
                            grep ("ATC", rownames (matrix)),                              # remove chloroplast genes
                            grep ("ATM", rownames (matrix))),]                            # remove mitochondrial genes
  
  # remove genes non-encoding proteins based on prepared earlier gene description file
  
  filt_mat <- filt_mat [which (rownames(filt_mat) %in% 
                                 Araport [Araport$type == "protein_coding","gene"]),]
                                 
  # remove genes below set average expression threshold
  
  filt_mat <- filt_mat [rowMeans(filt_mat) > mean_exp,]
  
  # remove seeds with less than set number of sequenced reads
  
  filt_mat <- filt_mat [,colSums(filt_mat)> n_reads]
  
  return (filt_mat)
}
