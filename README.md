Here we provide code for the preparation of the NGS figures for analysis of Arabidopsis seed size correlation with seed transcriptomes.

#### Load required R packages.

```
library (DESeq2)
library (tidyverse)
library (ggthemes)
library (eulerr)
library (gprofiler2)
library (rlist)
library (Seurat)
library (sctransform)
library (VISION)
library (cowplot)
library (patchwork)
library (viridis)
library (scran)
library (scater)
library (graph)
library (RBGL)
```

## Import the data

On GEO we deposited results of three different NGS experiments. 
1) 3'RNA-seq results for the analysis of Col-0 dry seed pools divided according to the size into small (S), medium (M) and large (L) seeds (4 biological replicas). Here we provide read counts as a single file.

``` R
rnaseq_size <- read.csv("/seed_size.tsv", header=T, sep = "\t")
```

2) single seed RNA-seq for *dog1-4* mutant and Col-0 in two time points of secondary dormancy induction (GSEXXXXX). 
   Each of the four treatments (genotype + time point) consist of 3 libraries. We provide a matrix of read counts for each library in which rows are genes and columns are one of 32 seeds used for library preparation.

3) single seed RNA-seq for Col-0 small and large seeds in two time points of secondary dormancy induction. 
   Seeds treatment and libraries preparation for that analysis were performed along with the experiment from the previous point.

We created the function called import_counts to upload the data and combine matrices into a single matrix.

Data for *dog1-4* and Col-0 single seed experiment.
``` R
data_dog1 <- import_counts ("/matrix/dog1/", header = TRUE)
```

Data for the single seed size experiment.
``` R
data_size <- import_counts ("/matrix/size/", header = TRUE)
```

# seed size 3'RNA-seq analysis

We sequenced mRNAs isolated from dry seeds divided according to their sizes.
We identified DEGs using DESeq2.

``` R
# metadata
coldata <- data.frame(lib = colnames(rnaseq_size), size = substr (colnames(rnaseq_size), 1,1), 
                      replica = substr (colnames(rnaseq_size), 2,2))

# DESeq2 analysis
dds <- DESeqDataSetFromMatrix(countData = rnaseq_size, colData = coldata, design = ~ size)
dds <- DESeq(dds)

# PCA plot
plot <- plotPCA(rlogTransformation(dds, blind=TRUE), intgroup=c('size', 'replica'), ntop =500)$data
plot$replica <- as.factor(plot$replica)

ggplot(plot , aes(x=PC1, y=PC2, color = size, label = replica)) +
  geom_point(size = 5) + 
  scale_color_tableau() +
  theme_classic() + 
  geom_text(size=8, hjust=-.25, vjust=-.25) + 
  guides(label = FALSE)
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/pca_rnaseq.jpeg" width=33% height=33%>

``` R
# identify affected genes in pairwise comparisons

deg_rnaseq <- list (SvL = results(dds, alpha = 0.05, contrast= c("size","S","L")),
                    SvM = results(dds, alpha = 0.05, contrast= c("size","S","M")),
                    MvL = results(dds, alpha = 0.05, contrast= c("size","M","L")))
                    
# plot overlaps between affected genes
plot <- lapply (deg_rnaseq, function(x) rownames(x [which(x$padj < 0.05),]))

plot(euler(plot), quantities = TRUE, fill = c("#0073C2FF", "#EFC000FF","#E15759"))
```   
 <img src="https://github.com/mk1859/seed_size/blob/main/images/venn_degs.jpeg" width=30% height=30%>

``` R
# select affected genes for small vs large seeds comparison
affected_genes <- rownames(deg_rnaseq$SvL [which(deg_rnaseq$SvL$padj < 0.05),])

# export normalized gene counts and filter for affected genes
norm_genes <- counts(dds, normalized = TRUE)
norm_genes <- norm_genes [which(rownames(norm_genes)%in% affected_genes),]

# scale gene expression
norm_genes <- as.data.frame(scale(t(norm_genes)))

# add seed size as column
norm_genes$condition <- as.factor(substr(rownames(norm_genes),1,1))

# calculate mean expression for replicas
norm_genes <- as.data.frame(t(apply (norm_genes [,-ncol(norm_genes)], 2, function (x) {
                               tapply (x, norm_genes$condition, mean)})))
  
# cluster genes
gene_clusters <- norm_genes %>% 
  dist(.) %>%
  hclust(., method = "complete") %>%
  cutree(., k = 2) %>%
  enframe(., name = "gene", value = "cluster")

# create data frame for plotting
norm_genes <- norm_genes %>%
  mutate (.,gene = rownames(norm_genes)) %>%
  pivot_longer(., cols = L:S, names_to = "size", values_to = "exp") %>%
  merge (.,gene_clusters, by= "gene") 
norm_genes$size <- factor(norm_genes$size, levels = c("S", "M", "L"))

ggplot(norm_genes, aes(size, exp, color = as.factor(cluster))) +
  geom_line(aes(group = gene), alpha = 0.3) +
  facet_wrap(~ cluster, nrow = 1)+ 
  theme_classic() + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_color_manual(values=c("#E15759", "#4E79A7"))
 ```
 <img src="https://github.com/mk1859/seed_size/blob/main/images/clusters_genes.jpeg" width=30% height=30%>
 
 ``` R
# GO term enrichment for affected genes
background_genes <- rownames(deg_rnaseq$SvL [which(deg_rnaseq$SvL$baseMean > 5),])

# genes upregulated in small seeds
go_plot (rownames(deg_rnaseq$SvL [which(deg_rnaseq$SvL$padj < .05 &deg_rnaseq$SvL$log2FoldChange > 0),]), background_genes) 
```  
 <img src="https://github.com/mk1859/seed_size/blob/main/images/go_small.jpeg" width=30% height=30%> 
 
``` R
# genes upregulated in large seeds
go_plot (rownames(deg_rnaseq$SvL [which(deg_rnaseq$SvL$padj < .05 &deg_rnaseq$SvL$log2FoldChange < 0),]), background_genes) 
``` 
 <img src="https://github.com/mk1859/seed_size/blob/main/images/go_large.jpeg" width=30% height=30%>  
  
``` R
# among affected genes we found a group with important function in seed biology
intresting_genes <- read.csv("/size/intresting_genes.txt",header=T, sep = "\t", dec =".")

plot <- norm_genes [which(norm_genes$gene %in% intresting_genes$gene),]

set_order <- plot$gene [order(plot$size,plot$exp)]
set_order <- set_order [!duplicated(set_order)]

plot$gene <- factor(plot$gene, level = set_order)

ggplot(plot, aes(size, gene, fill = exp)) +
  geom_tile () +
  theme_classic() +
  scale_fill_gradient2(midpoint=0, high="#4E79A7", mid="white", low="#E15759", space ="Lab")
```
 <img src="https://github.com/mk1859/seed_size/blob/main/images/heatmap.jpeg" width=30% height=30%>  

``` R
# we noticed that genes differentially expressed between small and large seeds 
# are similar to genes underlying germination competence index from our previous work (REF)  
germ_genes <- read.csv("/size/germ_genes_timecourse.txt",header=T, sep = "\t", dec =".")

# we repaeated gene expression scalling for the genes included in germination index
norm_genes <- counts(dds, normalized = TRUE)
norm_genes <- norm_genes [which(rownames(norm_genes)%in% germ_genes$gene),]
norm_genes <- as.data.frame(scale(t(norm_genes)))
norm_genes$condition <- as.factor(substr(rownames(norm_genes),1,1))
norm_genes <- as.data.frame(t(apply (norm_genes [,-ncol(norm_genes)], 2, function (x) {
  tapply (x, norm_genes$condition, mean)})))

norm_genes <- norm_genes %>%
  mutate (.,gene = rownames(norm_genes)) %>%
  pivot_longer(., cols = L:S, names_to = "size", values_to = "exp") %>%
  merge (.,germ_genes, by= "gene")

norm_genes$size <- factor(norm_genes$size, levels = c("S", "M", "L"))
norm_genes$gene.groups <- factor(norm_genes$gene.groups, levels = c("2", "1"))

# we ploted expression of these genes in small, medium and large seeds
ggplot(norm_genes, aes(size, exp, color = as.factor(size))) +
  geom_boxplot (size = 2) +
  facet_wrap(~ gene.groups, nrow = 1)+ 
  theme_classic() + 
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text.x = element_blank()) + 
  scale_color_manual(values=c("#E15759","#F28E2B", "#4E79A7"))
``` 
 <img src="https://github.com/mk1859/seed_size/blob/main/images/boxplot_germ.jpeg" width=30% height=30%> 
 
# single seed RNA-seq

Quality control of Col-0 and *dog1-4* experiment was performed earlier (REF and https://github.com/mk1859/single_seed).
Here we start with quality controls for smal/large seeds experiment.

Our library preparation protocol was designed to detect mRNAs. To filter out non-protein-coding genes we needed a reference file with information about gene types.

``` R
Araport <- read.csv ("Araport.txt", sep = "\t", header = TRUE)

head (Araport)
```
```
       gene           type
1 AT1G01010 protein_coding
2 AT1G01020 protein_coding
3 AT1G01030 protein_coding
4 AT1G01040 protein_coding
5 AT1G01050 protein_coding
6 AT1G01060 protein_coding
```
 
## Pre-filtering single seed data

Similarly to single-cell experiments, our reads count data is sparse. We needed to clean it by:
1) removing non-protein-coding genes
2) removing of genes encoded in organelles
3) removing summary lines at the last rows of the count matrix
4) filtering out genes with a low count number
5) filtering seeds with not enough reads

To do that we created the function prefilter_matrix and applied it to our single seed matrices. By default, it uses Araport data frame with columns described above.
We require the mean expression of a gene to be at least 1 read per seed for a gene to remain and at least 5,000 reads per seed for a seed to remain.

``` R
filtered_size <- prefilter_matrix (data_size, mean_exp=1, n_reads=5000)

dim (filtered_size) # genes / seeds remaining
```
```
[1] 11785   382
```
We wrote a function to plot the number of sequenced reads and identified genes per seed. We wanted to show treatments in the specified order.

``` R
order_lib <- c ("SD_small_3d","SD_small_7d24h","SD_large_3d","SD_large_7d24h")

nreads_plot (filtered_size, tableu = "Green-Orange-Teal", order = order_lib)
```
 <img src="https://github.com/mk1859/seed_size/blob/main/images/nreads.jpeg" width=30% height=30%> 

As visible on the plot above, our libraries vary in the number of identified genic reads. One source of this may be the different quality of sequenced libraries reflected by the ratio of target and off-target reads.
We counted the fraction of off-target reads (not in protein-coding genes) for each seed with the background_reads function which uses raw and pre-filtered matrices as input.

``` R
background_size <- background_reads (data_size, filtered_size)

# function to boxplot fraction of background reads
background_plot (filtered_size, order = order_lib, background = background_size)
```
 <img src="https://github.com/mk1859/seed_size/blob/main/images/boxplot_background.jpeg" width=30% height=30%> 
 
The abundance of background reads may imply that some counts attributed to genes may not reflect their expression.
A closer examination of reads' tracks in the browser showed that the distribution of background reads is not random and they tend to create hot spots laying both between genes and partially overlapping with them. In addition, the strength of genic peaks is negatively correlated with the number of background reads.
Based on these observations, we decided to remove from our analysis genes whose read count is strongly positively correlated with the number of background reads. As gene expression patterns are different between treatments, we calculated these correlations for each of them separately as well as for all seeds combined. To do this we wrote the function called correlation_table.

``` R
correlation_size <- correlation_table (filtered_size, background_size)

filtered_size <- filtered_size [-which (rowMaxs (correlation_size) > 0.3),]
nrow (filtered_size) # genes remaining
```
```
[1] 7287
```

After we obtained filtered matrices of counts, we created Seurat objects with sctransform normalization. To do this, we prepared a wrapper function that takes the count matrix and extracts information about seeds from their names.
``` R
seurat_size <- seurat_object (filtered_size, background = background_size)
```

We calculated PCA during the preparation of Seurat objects. Now, we plotted it to show treatments with the pca_discrete function.
This function exports dimension reduction and metadata from the Seurat object. It is possible to choose a colour pallet from the ggthemes package.
``` R
pca_discrete (seurat_size, "timepoint", order = order_lib)
```
 <img src="https://github.com/mk1859/seed_size/blob/main/images/pca_size.jpeg" width=40% height=40%> 
 
Some technical parameters like number of reads, number of identified genes and fraction of background reads may affect the position of seeds on the PCA plot.
To check continuous values on PCA plots we wrote another plotting function.
``` R
pca_continuous (seurat_size, column = "log10_reads")
pca_continuous (seurat_size, column = "n_gene")
pca_continuous (seurat_size, column = "background")
```
 <img src="https://github.com/mk1859/seed_size/blob/main/images/pca_size_reads.jpeg" width=30% height=30%>  <img src="https://github.com/mk1859/seed_size/blob/main/images/pca_size_genes.jpeg" width=30% height=30%>  <img src="https://github.com/mk1859/seed_size/blob/main/images/pca_size_background.jpeg" width=30% height=30%> 

After we showed that small and large single seed RNA sequencing has high quality, we can combine read counts from both small/large and Col-0/*dog1-4* experiments and performed their analysis together. First, we need to repeat all filtering steps on the combined matrix of counts and then create the Seurat object.
 ``` R
data_both <- cbind (data_size, data_dog1)
filtered_both <- prefilter_matrix (data_both, mean_exp=1, n_reads=5000)
background_both <- background_reads (data_both, filtered_both)
correlation_both <- correlation_table (filtered_both, background_both)
filtered_both <- filtered_both [-which (rowMaxs (correlation_both) > 0.3),]
seurat_both <- seurat_object (filtered_both, background = background_both)

order_lib <- c ("SD_small_3d","SD_small_7d24h","SD_large_3d","SD_large_7d24h", 
                "SD_dog1_3d","SD_dog1_7d24h","SD_Col0_3d","SD_Col0_7d24h")
                
pca_discrete (seurat_both, "timepoint", order = order_lib)
```
 <img src="https://github.com/mk1859/seed_size/blob/main/images/pca_both.jpeg" width=40% height=40%>
      
To better visualise seeds' grouping, we performed tSNE transformation.
``` R
seurat_both <- RunTSNE(object = seurat_both, dims = 1:15, verbose = FALSE, perplexity = 40, 
                       theta =0, max_iter = 10000)

# customize tSNE plot
plot <- cbind (as.data.frame (Embeddings(object = seurat_both, reduction = "tsne")), 
               seurat_both@meta.data)

plot$timepoint <- factor(plot$timepoint, levels = order_time)

ggplot(plot, aes(x=tSNE_1, y= tSNE_2, color = timepoint)) +
  geom_point (size = 2) + 
  scale_color_tableau("Green-Orange-Teal") +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
```    
<img src="https://github.com/mk1859/seed_size/blob/main/images/tsne.jpeg" width=40% height=40%>

## Gene expression patterns

To explain gene expression patterns underlying seed positions on the PCA map, we performed a few analyses.
First, we clustered seeds and identified genes differentially expressed between seeds' clusters.
``` R
seurat_both <- FindNeighbors(object = seurat_both, dims = 1:15, verbose = FALSE)

seurat_both <- FindClusters(object = seurat_both, verbose = FALSE, resolution = 0.2)

both_cluster <- FindAllMarkers(object = seurat_both, 
                                logfc.threshold = log2(1.5), test.use = "wilcox",only.pos =TRUE, 
                                assay = "SCT", slot ="data", verbose = FALSE)

both_cluster <- both_cluster [both_cluster$p_val_adj < .05,]

plot <- cbind (as.data.frame (Embeddings(object = seurat_both, reduction = "pca")) [,1:2], 
               seurat_both@meta.data)

ggplot(plot, aes(x=PC_1, y= PC_2, color = seurat_clusters)) +
  geom_point (size = 2) + 
  scale_color_tableau("Classic 10") +
  theme_classic() +
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank()) +
  theme(legend.title = element_blank())
``` 
<img src="https://github.com/mk1859/seed_size/blob/main/images/seed_clusters.jpeg" width=40% height=40%>

We checked if genes characteristic for seeds' clusters show any enrichment for GO terms
``` R
# first set background to genes that passed filtering steps
background_genes <- rownames(filtered_both)

# cluster 0
go_plot (both_cluster$gene [which(both_cluster$cluster == 0)], background_genes) 
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/cluster0.jpeg" width=30% height=30%>

``` R
# cluster 1
go_plot (both_cluster$gene [which(both_cluster$cluster == 1)], background_genes) 
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/cluster1.jpeg" width=20% height=20%>

``` R
# cluster 2
go_plot (both_cluster$gene [which(both_cluster$cluster == 2)], background_genes) 
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/cluster2.jpeg" width=10% height=10%>

The most direct way to identify genes affected by the seed size and *dog1-4* mutation is DGE analysis. 
To do that we created a wrapper function for Seurat FindMarkers.

``` R
deg_both <- deg_list (seurat_both, 
                      vector1 = c ("SD_small_3d","SD_small_7d24h", "SD_Col0_3d","SD_Col0_7d24h"), 
                      vector2 = c ("SD_large_3d","SD_large_7d24h", "SD_dog1_3d","SD_dog1_7d24h"), 
                      column = "timepoint", padj = 0.05, log2FC_threshold = log2(1.2))
                      
# plot number of affected genes using another function.
deg_plot (deg_both, direction = TRUE)
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/degs_single.jpeg" width=60% height=60%>

``` R
# we created Venn diagrams to show overlaps between identified genes.
plot <- list(large_3d = rownames(deg_both [[1]] [which(deg_both [[1]]$avg_log2FC > 0),]), 
             small_3d = rownames(deg_both [[1]] [which(deg_both [[1]]$avg_log2FC < 0),]),
             large_7d24h = rownames(deg_both [[2]] [which(deg_both [[2]]$avg_log2FC > 0),]), 
             small_7d24h = rownames(deg_both [[2]] [which(deg_both [[2]]$avg_log2FC < 0),]),
             dog_3d = rownames(deg_both [[3]] [which(deg_both [[3]]$avg_log2FC > 0),]), 
             col_3d = rownames(deg_both [[3]] [which(deg_both [[3]]$avg_log2FC < 0),]),
             dog_7d24h = rownames(deg_both [[4]] [which(deg_both [[4]]$avg_log2FC > 0),]), 
             col_7d24h = rownames(deg_both [[4]] [which(deg_both [[4]]$avg_log2FC < 0),]))

# 3d
plot(euler(plot[c(1,2,5,6)]), quantities = TRUE)
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/venn_3d.jpeg" width=30% height=30%>

``` R
# 7d+24h
plot(euler(plot[c(3,4,7,8)]), quantities = TRUE)

```
<img src="https://github.com/mk1859/seed_size/blob/main/images/venn_7d24h.jpeg" width=30% height=30%>

Finally, we also identified co-expressed gene groups using the coexpressed function. It calculates pairwise gene expression correlations and filters them. Next, it creates a graph object, looks for highly connected groups in it and outputs groups with gene numbers above the set threshold.

``` R
clusters <- coexpressed (seurat_both, threshold = 0.5, n_gene = 10)
               
lengths (clusters)
```
```
cluster_1 cluster_2 cluster_3 cluster_4 
      492       146        41        23 
```

Identified co-expressed gene groups were used to create signatures that were plotted on PCA maps.

``` R
seurat_both <- AddModuleScore(seurat_both, features = clusters, name = "cluster_")

# we use function to plot signatures
signature_map (seurat_both, signature = "cluster_1", order = order_lib, column = "timepoint")

signature_map (seurat_both, signature = "cluster_2", order = order_lib, column = "timepoint")

signature_map (seurat_both, signature = "cluster_3", order = order_lib, column = "timepoint")

signature_map (seurat_both, signature = "cluster_4", order = order_lib, column = "timepoint")

```
<img src="https://github.com/mk1859/seed_size/blob/main/images/sig_clust1.jpeg" width=20% height=20%> <img src="https://github.com/mk1859/seed_size/blob/main/images/sig_clust2.jpeg" width=20% height=20%> <img src="https://github.com/mk1859/seed_size/blob/main/images/sig_clust3.jpeg" width=20% height=20%> <img src="https://github.com/mk1859/seed_size/blob/main/images/sig_clust4.jpeg" width=20% height=20%>

We observe specific GO terms enrichment for two clusters.
``` R
# first set background to genes that passed filtering steps
background_genes <- rownames(filtered_both)

# cluster_1
go_plot (clusters$cluster_1, background_genes) 
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/go_cluster1.jpeg" width=10% height=10%>

``` R
# cluster_4
go_plot (clusters$cluster_4, background_genes) 
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/go_cluster4.jpeg" width=30% height=30%>

We observed that cluster_1 and cluster_2 expression are anticorelated
``` R
sig_vs_sig (seurat_both, "cluster_1", "cluster_2", order = order_lib)
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/clust1_vs_clust2.jpeg" width=30% height=30%>

To combine two signatures with opposing expression we use VISION package.
``` R
# effect of seed size
deg_size <- as.data.frame(deg_rnaseq$SvL) %>% filter (., padj < 0.05)
size_sign <-  -sign(deg_size$log2FoldChange)
size_sign <- setNames (size_sign, rownames(deg_size))
size_sign <- createGeneSignature (name = "size_sign", sigData = size_sign)

# cluster_1 and cluster_2 genes
germ_sign <- c(rep (1, 492),rep (-1, 146))
germ_sign <- setNames (germ_sign, c(clusters$cluster_1, clusters$cluster_2))
germ_sign <- createGeneSignature (name = "germ_sign", sigData = germ_sign)
vis <- Vision(seurat_both, signatures = list(size_sign,germ_sign), meta = seurat_both@meta.data, assay = "SCT")
vis <- analyze(vis)

# signature of seed size
signature_map (seurat_both, vis_obj = vis, signature = "size_sign", order = order_lib, column = "timepoint")
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/size_sign.jpeg" width=30% height=30%>

``` R
# signature of transcriptional germination competence
signature_map (seurat_both, vis_obj = vis, signature = "germ_sign", order = order_lib, column = "timepoint")
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/germ_sign.jpeg" width=30% height=30%>

``` R
# correlation of these two signatures
sig_vs_sig (seurat_both, "size_sign", "germ_sign", vis_obj = vis, order = order_lib)
```
<img src="https://github.com/mk1859/seed_size/blob/main/images/size_vs_germ.jpeg" width=30% height=30%>

```
sessionInfo()

R version 4.0.2 (2020-06-22)
Platform: x86_64-pc-linux-gnu (64-bit)
Running under: Ubuntu 20.04 LTS

Matrix products: default
BLAS/LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.8.so

locale:
 [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C               LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8     LC_MONETARY=en_US.UTF-8   
 [6] LC_MESSAGES=C              LC_PAPER=en_US.UTF-8       LC_NAME=C                  LC_ADDRESS=C               LC_TELEPHONE=C            
[11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] RBGL_1.66.0                 graph_1.68.0                scater_1.18.6               scran_1.18.7               
 [5] SingleCellExperiment_1.12.0 viridis_0.6.2               viridisLite_0.4.0           patchwork_1.1.1            
 [9] cowplot_1.1.1               VISION_2.1.0                sctransform_0.3.2.9006      SeuratObject_4.0.4         
[13] Seurat_4.0.6                rlist_0.4.6.2               gprofiler2_0.2.1            eulerr_6.1.1               
[17] ggthemes_4.2.4              forcats_0.5.1               stringr_1.4.0               dplyr_1.0.8                
[21] purrr_0.3.4                 readr_2.1.1                 tidyr_1.2.0                 tibble_3.1.6               
[25] ggplot2_3.3.5               tidyverse_1.3.1             DESeq2_1.30.1               SummarizedExperiment_1.20.0
[29] Biobase_2.50.0              MatrixGenerics_1.2.1        matrixStats_0.61.0          GenomicRanges_1.42.0       
[33] GenomeInfoDb_1.26.7         IRanges_2.24.1              S4Vectors_0.28.1            BiocGenerics_0.36.1        

loaded via a namespace (and not attached):
  [1] utf8_1.2.2                reticulate_1.22           tidyselect_1.1.2          RSQLite_2.2.9             AnnotationDbi_1.52.0     
  [6] htmlwidgets_1.5.4         grid_4.0.2                BiocParallel_1.24.1       Rtsne_0.15                munsell_0.5.0            
 [11] codetools_0.2-16          ica_1.0-2                 statmod_1.4.36            future_1.23.0             miniUI_0.1.1.1           
 [16] withr_2.5.0               colorspace_2.0-3          fastICA_1.2-3             rstudioapi_0.13           ROCR_1.0-11              
 [21] plumber_1.1.0             tensor_1.5                pbmcapply_1.5.0           listenv_0.8.0             labeling_0.4.2           
 [26] GenomeInfoDbData_1.2.4    polyclip_1.10-0           farver_2.1.0              bit64_4.0.5               parallelly_1.30.0        
 [31] vctrs_0.4.0               generics_0.1.2            R6_2.5.1                  ggbeeswarm_0.6.0          rsvd_1.0.5               
 [36] locfit_1.5-9.4            bitops_1.0-7              spatstat.utils_2.3-0      cachem_1.0.6              webutils_1.1             
 [41] DelayedArray_0.16.3       assertthat_0.2.1          promises_1.2.0.1          scales_1.1.1              beeswarm_0.4.0           
 [46] gtable_0.3.0              beachmat_2.6.4            globals_0.14.0            goftest_1.2-3             rlang_1.0.2              
 [51] genefilter_1.72.1         splines_4.0.2             lazyeval_0.2.2            spatstat.geom_2.3-1       broom_0.7.12             
 [56] reshape2_1.4.4            abind_1.4-5               modelr_0.1.8              backports_1.4.1           httpuv_1.6.5             
 [61] tools_4.0.2               logging_0.10-108          ellipsis_0.3.2            spatstat.core_2.3-2       RColorBrewer_1.1-3       
 [66] wordspace_0.2-6           ggridges_0.5.3            Rcpp_1.0.8.3              plyr_1.8.7                sparseMatrixStats_1.2.1  
 [71] zlibbioc_1.36.0           RCurl_1.98-1.5            rpart_4.1-15              deldir_1.0-6              pbapply_1.5-0            
 [76] zoo_1.8-9                 swagger_3.33.1            haven_2.4.3               ggrepel_0.9.1             cluster_2.1.0            
 [81] fs_1.5.2                  magrittr_2.0.3            data.table_1.14.2         scattermore_0.7           lmtest_0.9-40            
 [86] reprex_2.0.1              RANN_2.6.1                fitdistrplus_1.1-6        hms_1.1.1                 mime_0.12                
 [91] xtable_1.8-4              XML_3.99-0.8              sparsesvd_0.2             mclust_5.4.9              readxl_1.3.1             
 [96] gridExtra_2.3             compiler_4.0.2            KernSmooth_2.23-17        crayon_1.5.1              htmltools_0.5.2          
[101] mgcv_1.8-31               later_1.3.0               tzdb_0.2.0                geneplotter_1.68.0        lubridate_1.8.0          
[106] DBI_1.1.2                 dbplyr_2.1.1              MASS_7.3-51.6             Matrix_1.4-0              permute_0.9-5            
[111] cli_3.2.0                 igraph_1.3.0              pkgconfig_2.0.3           scuttle_1.0.4             plotly_4.10.0            
[116] spatstat.sparse_2.1-0     xml2_1.3.3                annotate_1.68.0           vipor_0.4.5               dqrng_0.3.0              
[121] iotools_0.3-2             XVector_0.30.0            rvest_1.0.2               digest_0.6.29             RcppAnnoy_0.0.19         
[126] vegan_2.5-7               spatstat.data_2.1-2       polylabelr_0.2.0          cellranger_1.1.0          leiden_0.3.9             
[131] edgeR_3.32.1              uwot_0.1.11               DelayedMatrixStats_1.12.3 loe_1.1                   shiny_1.7.1              
[136] lifecycle_1.0.1           nlme_3.1-148              jsonlite_1.8.0            BiocNeighbors_1.8.2       limma_3.46.0             
[141] fansi_1.0.3               pillar_1.7.0              lattice_0.20-41           fastmap_1.1.0             httr_1.4.2               
[146] survival_3.1-12           glue_1.6.2                png_0.1-7                 bluster_1.0.0             bit_4.0.4                
[151] stringi_1.7.6             blob_1.2.2                BiocSingular_1.6.0        memoise_2.0.1             irlba_2.3.5              
[156] future.apply_1.8.1 
```
