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
deg_dry_dog1 <- as.data.frame(results(dds, alpha = 0.05, contrast= c("genotype","dog1","Col0")))

# Volcano plot
ggplot(deg_dry_dog1 , aes(y=-log10(padj), x= log2FoldChange, color = padj < 0.05 , alpha = padj < 0.05)) +
  geom_point(size = 1) + 
  scale_color_tableau() +
  theme_classic() +
  scale_alpha_ordinal(range = c(0.1, 1))
```



