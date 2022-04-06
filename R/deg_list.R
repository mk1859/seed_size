# a wrapper function to Seurat FindMarkers that allows to pairwise compare gene expression between conditions named in two vectors
# as input it requires Seurat object, a column of metadata containing names of conditions, two vectors with condition names as well as 
# padj and log2 fold change thresholds
# function generates list as output which elements contain data frames with significant DEGs

deg_list <- function(seurat_obj, vector1, vector2, column, padj = 0.05, log2FC_threshold = log2(2)) {
    require (Seurat)
    require (dplyr)
    
    # set column to be active identity in Seurat object
    seurat_obj@active.ident <- as.factor(seurat_obj [[column]] [,1])
    names (seurat_obj@active.ident) <- rownames (seurat_obj@meta.data)
    
    # DEG identification
    deg_list <- mapply (function (vector1,vector2) {
          FindMarkers(object = seurat_obj, ident.1 = vector2,
                      ident.2 = vector1,logfc.threshold = log2FC_threshold, 
                      test.use = "wilcox",only.pos =FALSE,
                      assay = "SCT", slot ="data") %>%
                    filter(.,p_val_adj < padj)},
          vector1,vector2, SIMPLIFY = FALSE)
    
    # set names to results
    names(deg_list) <- paste0(vector1, "_", vector2)
    
    return (deg_list)
}
