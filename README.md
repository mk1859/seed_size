Here we provide code for the preparation of the NGS figures for analysis of Arabidopsis seed size correlation with seed germination phenotypes.

#### Load required R packages.

```
library (DESeq2)
library (ggplot2)
library (ggthemes)



library (tidyverse)
library (ggthemes)
library (ggbeeswarm)
library (viridis)
library (sctransform)
library (Seurat)
library (rlist)
library (cowplot)
library (patchwork)
library (gprofiler2)
library (rrvgo)
library (scran)
library (scater)
library (graph)
library (RBGL)
library (data.table)
library (eulerr)
library (DESeq2)
library (VISION)
library (org.At.tair.db)
library (biomaRt)
library (GO.db)
library (scales)
library (matrixStats)
library (UpSetR)
```

## Import the data

On GEO we deposited results of three different NGS experiments. 
1) 3'RNA-seq results for the analysis of Col-0 dry seed pools divided according to the size into small (S), medium (M) and large (L) seeds (4 biological replicas). Here we provide read counts as single file.

``` R
rnaseq_size <- read.csv("/seed_size.tsv", header=T, sep = "\t")
```

2) single seed RNA-seq for *dog1-4* mutant and Col-0 in two time points of secondary dormancy induction (GSEXXXXX). 
   Each of the four treatments (genotype + time point) consist of 3 libraries. We provide a matrix of read counts for each library in which rows are genes and columns are one of 32 seeds used for library preparation.

3) single seed RNA-seq for Col-0 small and large seeds in two time points of secondary dormancy induction. 
   Seeds treatment and libraries preparation for that analysis were performed along the experiment from previous point.

We created the function called import_counts to upload the data and combine matrices into a single matrix.

Data for *dog1-4* and Col-0 single seed experiment.
``` R
data_dog1 <- import_counts ("/matrix/dog1/", header = TRUE)
```

Data for the single seed size experiment.
``` R
data_size <- import_counts ("/matrix/size/", header = TRUE)
```

## seed size 3'RNA-seq analysis

We sequenced mRNAs isolated from dry seeds divided accroding to their sizes.
We identified DEGs using DESeq2.

``` R
# metadata
coldata <- data.frame(lib = colnames(rnaseq_size), size = substr (colnames(rnaseq_size), 1,1), replica = substr (colnames(rnaseq_size), 2,2))

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
 <img src="https://github.com/mk1859/seed_size/blob/main/images/cluster_genes.jpeg" width=70% height=70%>
 
 ``` R
  

deg_rnaseq <- as.data.frame(results(dds, alpha = 0.05, contrast= c("genotype","dog1","Col0")))

# Volcano plot
ggplot(deg_dry_dog1 , aes(y=-log10(padj), x= log2FoldChange, color = padj < 0.05 , alpha = padj < 0.05)) +
  geom_point(size = 1) + 
  scale_color_tableau() +
  theme_classic() +
  scale_alpha_ordinal(range = c(0.1, 1))
  
  
  
res <- results(dds, alpha = 0.05, contrast= c("size","S","L"))
length(res [which(res$padj < 0.05 & res$log2FoldChange > log(1.5)),2])
nrow(res [which(res$padj < 0.05 & res$log2FoldChange < -log(1.5)),])
write.table (res,"D:/drop/Dropbox/nowe_polecenia/size/SvL.txt", sep = "\t", col.names = T, row.names = T, quote = FALSE)


res <- results(dds, alpha = 0.05, contrast= c("size","S","M"))
nrow(res [which(res$padj < 0.05 & res$log2FoldChange > log(1.5)),])
nrow(res [which(res$padj < 0.05 & res$log2FoldChange < -log(1.5)),])
write.table (res,"D:/drop/Dropbox/nowe_polecenia/size/SvM.txt", sep = "\t", col.names = T, row.names = T, quote = FALSE)

res <- results(dds, alpha = 0.05, contrast= c("size","M","L"))
nrow(res [which(res$padj < 0.05 & res$log2FoldChange > log(1.5)),])
nrow(res [which(res$padj < 0.05 & res$log2FoldChange < -log(1.5)),])
write.table (res,"D:/drop/Dropbox/nowe_polecenia/size/MvL.txt", sep = "\t", col.names = T, row.names = T, quote = FALSE)

```



