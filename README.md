Here we provide code for the preparation of the NGS figures for analysis of Arabidopsis seed size correlation with seed germination phenotypes.

#### Load required R packages.

```
library (DESeq2)
library (tidyverse)
library (ggthemes)
library (eulerr)
library (gprofiler2)
library (Seurat)
library (sctransform)
library (VISION)
library (cowplot)
library (patchwork)
library (viridis)

library (ggbeeswarm)

library (rlist)
library (cowplot)
library (patchwork)

library (rrvgo)
library (scran)
library (scater)
library (graph)
library (RBGL)
library (data.table)

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
plot <- gost(query = rownames(deg_rnaseq$SvL [which(deg_rnaseq$SvL$padj < .05 &deg_rnaseq$SvL$log2FoldChange > 0),]),
             organism = "athaliana", custom_bg = background_genes, user_threshold = 0.05, sources = "GO")$result

plot$log_pval <- -log10(plot$p_value)
plot$term_name <- factor(plot$term_name, levels = as.factor(plot$term_name [order(c(plot$source, plot$log_pval), decreasing = F)]))

ggplot(plot, aes(log_pval, term_name, fill = source)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() + 
  scale_fill_tableau()
```  
 <img src="https://github.com/mk1859/seed_size/blob/main/images/go_small.jpeg" width=30% height=30%> 
 
``` R
# genes upregulated in large seeds
plot <- gost(query = rownames(deg_rnaseq$SvL [which(deg_rnaseq$SvL$padj < .05 &deg_rnaseq$SvL$log2FoldChange < 0),]),
             organism = "athaliana", custom_bg = background_genes, user_threshold = 0.05, sources = "GO")$result

plot$log_pval <- -log10(plot$p_value)
plot$term_name <- factor(plot$term_name, levels = as.factor(plot$term_name [order(c(plot$source, plot$log_pval), decreasing = F)]))

ggplot(plot, aes(log_pval, term_name, fill = source)) +
  geom_bar(stat="identity", position=position_dodge()) +
  theme_classic() + 
  scale_fill_tableau()
``` 
 <img src="https://github.com/mk1859/seed_size/blob/main/images/go_large.jpeg" width=30% height=30%>  
  
``` R
# among affected genes we found a group with important function in seed biology
intresting_genes <- read.csv("D:/drop/Dropbox/nowe_polecenia/size/intresting_genes.txt",header=T, sep = "\t", dec =".")

plot <- norm_genes [which(norm_genes$gene %in% intresting_genes$gene),]

set_order <- plot$gene [order(plot$size,plot$exp)]
set_order <- set_order [!duplicated(set_order)]

plot$gene <- factor(plot$gene, level = set_order)

ggplot(plot, aes(size, gene, fill = exp)) +
  geom_tile () +
  theme_classic() +
  scale_fill_gradient2(midpoint=0, high="#4E79A7", mid="white",
                       low="#E15759", space ="Lab")
```
 <img src="https://github.com/mk1859/seed_size/blob/main/images/heatmap.jpeg" width=30% height=30%>  

``` R
# we noticed that genes differentially expressed between small and large seeds 
# are similar to genes underlying germination competence index from our previous work (REF)  
germ_genes <- read.csv("D:/drop/Dropbox/nowe_polecenia/size/germ_genes_timecourse.txt",header=T, sep = "\t", dec =".")

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
 
## single seed RNA-seq

Quality control of Col-0 and *dog1-4* experiment was performed earlier (REF and https://github.com/mk1859/single_seed).
Here we start with quality controls for small large seeds experiment.

First, we need a refernce file as our library preparation protocol was designed to detect mRNAs. To filter out non-protein-coding genes we needed a reference file with information about gene types.

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
 
# Pre-filtering single seed data

Similarly to single-cell experiments, our count data is sparse. We needed to clean it by:
1) removing of non-protein-coding genes
2) removing of genes encoded in organelles
3) removing of summary lines at last rows of the count matrix
4) filtering out genes with a low count number
5) filtering seeds with not enough reads

To do that we created function prefilter_matrix and applied it to our single seed matrices. By default it uses Araport data frame with columns described as above.
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
background_timecourse <- background_reads (data_size, filtered_size)

# function to boxplot fraction of background reads
background_plot (filtered_size, order = order_lib, background = background_size)
```
 <img src="https://github.com/mk1859/seed_size/blob/main/images/boxplot_background.jpeg" width=30% height=30%> 
 
 The abundance of background reads may imply that some counts attributed to genes may not reflect their expression.
Closer examination of read tracks in the browser showed that the distribution of background reads is not random and they tend to create hot spots laying both between genes and partially overlapping with them. In addition, the strength of genic peaks is negatively correlated with the number of background reads.
Based on these observations, we decided to remove from our analysis genes whose read count is strongly positively correlated with the number of background reads. As gene expression patterns are different between treatments, we calculated these correlations for each of them separately as well as for all seeds combined. To do that we wrote the function called correlation_table.

``` R
correlation_size <- correlation_table (filtered_tsize, background_size)

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



After we showed that small and large single seed RNA sequencing is high quality, we can combine read counts both small/large and Col-0/*dog1-4* experiments and performed their analysis together. First we need to repeat all filtering steps on combined matrix of counts.
 
 
