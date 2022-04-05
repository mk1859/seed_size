# R function to upload a set of matrices containing read counts (with the same set of genes) from a directory.
# It assumes each matrix have gene names as row names.
# It is possible to set if matrices contain column names.

import_counts <- function(directory, header = TRUE) {
  require (rlist)
  
  files <- list.files (directory)
  
  matrix <- lapply (files, function (x){
    read.csv (paste0(directory, x), sep = "\t", header = header, row.names =1)
  })
  
  matrix <- list.cbind(matrix)
  
  return (matrix)
}
