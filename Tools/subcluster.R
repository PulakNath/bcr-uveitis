# This R script loads the last stored annotated data object and does subcluster analysis
# All the subcluster analysis results are in the Results/subcluster folder
# Author : Vijay Nagarajan, NEI/NIH
load("/data/../TotalSeq/Results/annotation/SeuratObjectAfterAnnotation.RData")
#load("SeuratObjectAfterAnnotation.RData")

library(Seurat)
library(kableExtra)
library(ggplot2)
library(plotly)
library(tidyverse)
library(rmarkdown)
library(cowplot)

Idents(allsamples) <- "annotated_RNA_snn_res.0.25"

levels(allsamples)

nk=subset(x=allsamples, idents=c('C5_NK cells'))

## ----PCA, message=FALSE, warning=FALSE----------------------------------------
##################### PCA
# Run PCA with 50 pcs
nk <- RunPCA(object = nk, npcs=50)
# Print genes in the top 5 pcs
print(nk[["pca"]], dims = 1:5, nfeatures = 5)
# Identify number of PCs that explains majority of variations
pdf("nk_elbow.pdf")
ElbowPlot(nk, ndims=50, reduction = "pca")
dev.off()

# Visualize the genes in the first PC
pdf("nk_genes_pc1.pdf")
VizDimLoadings(nk, dims = 1, ncol = 1) + theme_minimal(base_size = 8)
dev.off()

# Set the identity column to WT vs TG
Idents(nk) <- "orig.ident"
# Visualize the cells after pca
pdf("nk_pca.pdf")
DimPlot(object = nk, reduction = "pca")
dev.off()

# Plot heatmap with cells=500 plotting cells with extreme cells on both ends of spectrum
pdf("nk_heatmap.pdf")
DimHeatmap(object = nk, dims = 1:6, cells = 200, balanced = TRUE)
dev.off()

## ----clustering, message=FALSE, warning=FALSE---------------------------------
##################### Clustering
# use the first 20 pc's
use.pcs = 1:20
# FindNeighbors with 20 pcs
nk <- FindNeighbors(nk, reduction="pca", dims = use.pcs)
# FindClusters with resolution starting at 0.25 and ending at 4, with 0.5 increment
nk <- FindClusters(object = nk)

head(nk)
# Count number of clusters at each resolution
sapply(grep("res",colnames(nk@meta.data),value = TRUE),
       function(x) length(unique(nk@meta.data[,x])))

set.seed(10)

## ----tsne, message=FALSE, warning=FALSE---------------------------------------
####################### TSNE
# Run tsne with 20 pcs
nk <- RunTSNE(object = nk,reduction.use = "pca",dims = use.pcs,do.fast = TRUE)

# change default identity
Idents(nk) <- "RNA_snn_res.0.75"
# list cell number in each cluster for WT vs TG
table(Idents(nk),nk$type)
# Visualize cluster at 0.8 resolution
pdf("nk_tsne.pdf")
DimPlot(object = nk, pt.size=0.5, reduction = "tsne", label = T)
dev.off()
# Color by HC vs BCR-UV
pdf("nk_tsne_grouped.pdf")
DimPlot(object = nk, pt.size=0.5, reduction = "tsne", split.by = "type" )
dev.off()


## ----umap, message=FALSE, warning=FALSE---------------------------------------
########################## UMAP
# Run umap with 20 pcs
nk <- RunUMAP(
  object = nk,
  reduction.use = "pca",
  n.neighbors = 90L,
  dims = use.pcs)
# Plot umap
pdf("nk_umap.pdf")
DimPlot(object = nk, pt.size=0.5, reduction = "umap", label = T)
dev.off()

# Color by HC vs BCR-UV
pdf("nk_umap_group.pdf")
DimPlot(object = nk, pt.size=0.5, reduction = "umap", split.by = "type")
dev.off()


# Visualizations
# expression of variable genes across clusters
# Custom list of genes to visualize
#top10
mygenes = c("GZMB","GNLY")
pdf("nk_ridgeplot.pdf")
RidgePlot(nk, features = mygenes)
dev.off()

pdf("nk_vlnplot.pdf")
VlnPlot(nk, features = mygenes, pt.size=0)
dev.off()

#VlnPlot(WT_TG_Filt_Scaled, features = mygenes, split.by = "orig.ident")
pdf("nk_featureplot.pdf")
FeaturePlot(nk, features = mygenes)
dev.off()

#FeaturePlot(WT_TG_Filt_Scaled, features = mygenes, split.by = "orig.ident")
pdf("nk_dotplot.pdf")
DotPlot(nk, features = mygenes)
dev.off()

pdf("nk_heatmap.pdf")
DoHeatmap(subset(nk), features = mygenes, size = 3)
dev.off()


# Visualize ADT with GEX data
# Specify key for rna features and ADT features
# Add rna_ or adt_ to any gene to visualize
rownames(nk[["ADT"]])
Key(nk[["RNA"]])
Key(nk[["ADT"]])

p1 <- FeaturePlot(nk, "adt_CD19.1") + ggtitle("CD19.1 protein")
p2 <- FeaturePlot(nk, "rna_CD19") + ggtitle("CD19 RNA")
pdf("nk-adt-gex-cd19.pdf")
p1 | p2
dev.off()

###### feature plot at individual patient level, for any given gene

allsamplesplot<-FeaturePlot(nk, features = "GZMB", pt.size=0, split.by="orig.ident", combine=F)
pdf("all-samples-plot.pdf")
p9 <- list()
for (i in seq_along(allsamplesplot)){
    #Change x and y label font size.
    p9[[i]] = allsamplesplot[[i]] +
      theme(legend.position = "none",
		axis.text = element_text(size = 8),
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            #axis.text.y.right = element_text(size = 8),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            title = element_text(size=8))
}
plot_grid(plotlist = p9, ncol = 4)
dev.off()


########################### Marker genes between groups
# Current identity levels
levels(nk)
head(Idents(nk))
# Change identity to HC vs BCR-UV
Idents(nk) <- "type"
head(Idents(nk))
levels(nk)
########################### tables
# Proportion and counts of cells in clusters and groups
prop.table(table(Idents(nk)))
table(Idents(nk))
# Find differentially expressed genes for UV vs BCR-UV
nkmarkershcbcruv = FindMarkers(nk, ident.1="BCR-UV", ident.2="HC")
head(nkmarkershcbcruv)
# Set marker rownames
nkmarkershcbcruv$gene=rownames(nkmarkershcbcruv)
# Save the list of differentially expressed genes as a tab delimited text file
write.table(nkmarkershcbcruv, file="nkmarkershcbcruv.txt", quote=F, row.names = FALSE, sep = "\t")

pdf("nk-hc-uv-marker-vln.pdf")
VlnPlot(nk,features=c("HLA-B","FOS"), pt.size = 0)
dev.off()

pdf("nk-hc-uv-marker-feature.pdf")
FeaturePlot(nk, features="HLA-B", split.by="type")
dev.off()

table(nk$type)

# Get all the expression count matrix
#head(GetAssayData(WT_TG_Filt_Scaled, slot="counts"))
# Get expression values for a single gene from all cells
head(FetchData(nk, vars="HLA-B"))


## ----differentialbetweenclusters, message=FALSE, warning=FALSE----------------
########################### tables
# Proportion and counts of cells in clusters and groups
# Change identity to 0.25 cluster resolution
Idents(nk) <- "RNA_snn_res.0.75"
prop.table(table(Idents(nk)))
table(Idents(nk), nk$type)

########################### Marker genes between clusters
# Current identity levels
levels(nk)
head(Idents(nk))

# Find differentially expressed genes for cluster 5 vs 15
nkmarkerscluster = FindMarkers(nk, ident.1=5, ident.2=15)
head(nkmarkerscluster)
write.table(nkmarkerscluster, file="nkmarkerscluster.txt", quote=F, row.names = FALSE, sep = "\t")

pdf("nk-5-vs-15-marker-vln.pdf")
VlnPlot(nk,features=c("GZMB","IL4I1"), pt.size = 0)
dev.off()

pdf("nk-5-vs-15-marker-feature.pdf")
FeaturePlot(nk, features="GZMB", split.by="type")
dev.off()

table(nk$orig.ident)

# Get all the expression count matrix
#head(GetAssayData(WT_TG_Filt_Scaled, slot="counts"))
# Get expression values for a single gene from all cells
head(FetchData(allsamples, vars="GZMB"))


## ----differentialbetweeneveryclusters, message=FALSE, warning=FALSE-----------
########################### Find all markers
nkmarkersall=FindAllMarkers(nk)
head(nkmarkersall)
#View(markersall)
# Save all markers as tab delimited file
write.table(nkmarkersall, file="nkmarkersall.txt", quote=F, row.names = FALSE, sep = "\t")

# Extract the top marker for each cluster
head(nkmarkersall)
nkmarkersall %>%
  group_by(cluster) %>%
    slice(1) %>%
      pull(gene) -> nk.best.gene.per.cluster

nk.best.gene.per.cluster

# Store plot as ggplot list
nkvlnplot1<-VlnPlot(nk,features=nk.best.gene.per.cluster, pt.size = 0, combine = FALSE)
# vlnplot1<-VlnPlot(allsamples,features=best.gene.per.cluster, pt.size = 0, same.y.lims = TRUE, combine = FALSE)

# Set properties for each plot in the list
p10 <- list()
for (i in seq_along(nkvlnplot1)){
    #Change x and y tick label font size.
    p10[[i]] = nkvlnplot1[[i]] +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.y = element_text(size = 6),
            axis.title.x = element_blank(),
            title = element_text(size=8))
}
# Plot all the plots in the list
plot_grid(plotlist = p10, ncol = 3)
# Save plot to file
nkvlnplot1=plot_grid(plotlist = p10, ncol = 3)
ggsave("nk-vlnplot-best-per-cluster.png", plot = nkvlnplot1, width = 20, height = 15, units = "cm")

# Plot features in umap
nkfplot1<-FeaturePlot(nk,features=nk.best.gene.per.cluster, pt.size = 0, combine = FALSE)
nkp2 <- list()
for (i in seq_along(nkfplot1)){
    #Change x and y label font size.
    nkp2[[i]] = nkfplot1[[i]] +
      theme(legend.position = "none",
            axis.text.x = element_text(size = 8),
            axis.text.y = element_text(size = 8),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            title = element_text(size=8))
}
plot_grid(plotlist = nkp2, ncol = 4)
nkfplot1=plot_grid(plotlist = nkp2, ncol = 4)
ggsave("nk-fplot-best-per-cluster.png", plot = nkfplot1, width = 20, height = 15, units = "cm")




## ----differentialgroupsinacluster, message=FALSE, warning=FALSE---------------
########################### tables
# Proportion and counts of cells in clusters and groups
# Change identity to 0.25 cluster resolution
Idents(nk) <- "RNA_snn_res.0.75"
prop.table(table(Idents(nk)))
table(Idents(nk))

# Current identity levels
levels(nk)
head(Idents(nk))

# Take all cells in cluster 1, and find markers that separate 'HC' and 'BCR-UV' that are in the 'orig.ident' column. subset.ident would the default identity of clusters
nkmarkersingroup <- FindMarkers(nk, ident.1 = "HC", ident.2 = "BCR-UV", group.by = 'type', subset.ident = 5)

head(nkmarkersingroup)
# Save markers
write.table(nkmarkersingroup, file="nkmarkersingroup.txt", quote=F, row.names = FALSE, sep = "\t")
table(nk$type)

# Get all the expression count matrix
#head(GetAssayData(allsamples, slot="counts"))
# Get expression values for a single gene from all cells
head(FetchData(nk, vars="HLA-B"))

##################### Proportion Test and Annotation
library("scProportionTest")

## ----proportiontest, message=FALSE, warning=FALSE-----------------------------
# Sample proportion table
# Proportion and counts of cells in individual samples
# Change identity to 0.25 cluster resolution

# Proportion test
# assign sample group
nk$sample=ifelse(grepl("HC",nk$type), "HC", ifelse(grepl("UV",nk$type), "UV", ""))

Idents(allsamples) <- "RNA_snn_res.0.75"
prop.table(table(Idents(nk), nk$sample))
table(Idents(nk), nk$sample)
# Save tables
write.table(prop.table(table(Idents(nk), nk$sample)), file="cluster-proportions.txt", quote=F, sep = "\t")
write.table(table(Idents(nk), nk$sample), file="cluster-counts.txt", quote=F, sep = "\t")


# Run permutation test
nkprop_test <- sc_utils(nk)
nkprop_test <- permutation_test(
	nkprop_test, cluster_identity = "RNA_snn_res.0.75",
	sample_1 = "HC", sample_2 = "UV",
	sample_identity = "sample",
	n_permutations = 10000
)

permutation_plot(nkprop_test)
png("proportion.png")
permutation_plot(nkprop_test)
dev.off()


save.image("SeuratObjectAfterNKClustering.RData")
sessionInfo()


## ----reference----------------------------------------------------------------
#https://ucdavis-bioinformatics-training.github.io/2021-March-Single-Cell-RNA-Seq-Analysis/data_analysis/scRNA_Workshop-PART4_fixed
