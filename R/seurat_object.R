# function to create Surat object

seurat_object <- function(matrix, background = NULL, include_background=TRUE) {
  require (Seurat)
  require (sctransform)
  
  # metadata for each seed
  seed_attr <- data.frame(n_reads = colSums(matrix),                          # number of reads
                  n_gene = colSums(matrix > 0),                               # number of genes
                  log10_reads = log10(colSums(matrix)),                       # log10 number of reads
                  timepoint = gsub('.{0,6}$', '', colnames(matrix)),          # treatment extracted from seed name
                  batch = gsub('.{0,3}$', '', colnames(matrix)),              # batch extracted from seed name
                  background = background)                                    # fraction of background reads
                  
  seurat_seeds <- CreateSeuratObject(counts = matrix, meta.data = seed_attr)
  
  # sctransform with optional background reads fraction as value to regress
  if (include_background==TRUE) {
    seurat_seeds <- SCTransform(object = seurat_seeds, verbose =FALSE, n_genes = NULL,
                                vars.to.regress = 'background', return.only.var.genes = TRUE)
  }
  else {
    seurat_seeds <- SCTransform(object = seurat_seeds, verbose =FALSE, n_genes = NULL,
                                return.only.var.genes = TRUE)
  }
  
  # run PCA
  seurat_seeds <- RunPCA(object = seurat_seeds, verbose = FALSE)
  
  return (seurat_seeds)
}
