#install.packages("Seurat")
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GSVA")
# BiocManager::install("GSEABase")
# BiocManager::install("limma")
# BiocManager::install("SingleR")
# BiocManager::install("celldex")
# BiocManager::install("monocle")

getwd()

####################### 01. Preprocessing and correction #######################
# Load required packages
library(limma)
library(Seurat)
library(dplyr)
library(magrittr)
library(celldex)
library(SingleR)
library(monocle)
library(data.table)

logFCfilter = 1              # logFC filter threshold
adjPvalFilter = 0.05         # Adjusted p-value filter threshold
inputFile = "merge.txt"      # Single-cell expression matrix file
geneFile = "Genes.csv"       # Hub gene list file
setwd("D:/lu_singlecell/2.Seurat")    # Set working directory

# Read and organize input data
rt = fread(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)

# Convert matrix to Seurat object and filter cells
pbmc = CreateSeuratObject(counts = data, project = "seurat", min.cells = 10, min.features = 200, names.delim = "_")

# View slots in Seurat object
slotNames(pbmc)
# Assays
pbmc@assays
# Check cell metadata
dim(pbmc@meta.data)
View(pbmc@meta.data)

# QC: Calculate mitochondrial gene percentage for each cell
pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")

# Draw violin plot for features
pdf(file = "01.featureViolin.pdf", width = 15, height = 8)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

# QC: Quantile distribution statistics for each sample
gene.freq <- do.call("cbind", tapply(pbmc@meta.data$nFeature_RNA, pbmc@meta.data$orig.ident, quantile, probs = seq(0, 1, 0.05)))
rna.freq <- do.call("cbind", tapply(pbmc@meta.data$nCount_RNA, pbmc@meta.data$orig.ident, quantile, probs = seq(0, 1, 0.05)))
mt.freq <- do.call("cbind", tapply(pbmc@meta.data$percent.mt, pbmc@meta.data$orig.ident, quantile, probs = seq(0, 1, 0.05)))
freq.combine <- as.data.frame(cbind(gene.freq, rna.freq, mt.freq))
colnames(freq.combine) <- c(paste(colnames(gene.freq), "Gene", sep = "_"),
                            paste(colnames(rna.freq), "RNA", sep = "_"),
                            paste(colnames(mt.freq), "MT", sep = "_"))
write.table(freq.combine, file = "QC-gene_frequency.txt", quote = FALSE, sep = "\t")
rm(gene.freq, rna.freq, mt.freq)
View(freq.combine)

# QC: Scatter plots of nCount_RNA vs percent.mt and nCount_RNA vs nFeature_RNA
plot1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(file = "01.featureCor.pdf", width = 13, height = 7)
CombinePlots(plots = list(plot1, plot2), legend = "none")
dev.off()
rm(plot1, plot2)

# Cell filtering
cat("Before filter:", nrow(pbmc@meta.data), "cells\n")
pbmc <- subset(pbmc, 
               subset = 
                 nFeature_RNA > 500 & 
                 nCount_RNA > 1000 & 
                 nCount_RNA < 20000 & 
                 percent.mt < 10)
cat("After filter:", nrow(pbmc@meta.data), "cells\n")

# Scatter plots after filtering
plot1_1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2_1 <- FeatureScatter(pbmc, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
pdf(file = "01_1.featureCor.pdf", width = 13, height = 7)
CombinePlots(plots = list(plot1_1, plot2_1), legend = "none")
dev.off()
rm(plot1_1, plot2_1)

# Data normalization
pbmc <- NormalizeData(object = pbmc, normalization.method = "LogNormalize", scale.factor = 10000)

# Identify highly variable genes
pbmc <- FindVariableFeatures(object = pbmc, selection.method = "vst", nfeatures = 1500)
top10 <- head(x = VariableFeatures(object = pbmc), 10)
pdf(file = "01.featureVar.pdf", width = 10, height = 6)
plot1 <- VariableFeaturePlot(object = pbmc)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2), legend = "none")
dev.off()

####################### 02. Scaling and PCA analysis #######################
# Data scaling (may take time)
pbmc <- ScaleData(object = pbmc, do.scale = FALSE, do.center = FALSE, vars.to.regress = c("orig.ident", "percent.mt"))

# Run PCA
pbmc = RunPCA(object = pbmc, verbose = FALSE, npcs = 50, pc.genes = VariableFeatures(object = pbmc))

# Loadings for top PCs
pdf(file = "02.pcaGene.pdf", width = 10, height = 8)
VizDimLoadings(object = pbmc, dims = 1:4, reduction = "pca", nfeatures = 20)
dev.off()

# PCA plot
pdf(file = "02.PCA.pdf", width = 7.5, height = 5)
DimPlot(object = pbmc, reduction = "pca")
dev.off()

# Heatmap for top PCs
pdf(file = "02.pcaHeatmap.pdf", width = 10, height = 8)
DimHeatmap(object = pbmc, dims = 1:4, cells = 500, balanced = TRUE, nfeatures = 30, ncol = 2)
dev.off()

# JackStraw analysis
pbmc <- JackStraw(object = pbmc, num.replicate = 100, dims = 40)
pbmc <- ScoreJackStraw(object = pbmc, dims = 1:40)
pdf(file = "02.pcaJackStraw.pdf", width = 8, height = 6)
JackStrawPlot(object = pbmc, dims = 1:40)
dev.off()

# Elbow plot
pdf(file = "02.ElbowPlot.pdf", width = 5, height = 4)
ElbowPlot(object = pbmc, ndims = 40)
dev.off()

####################### 03. TSNE clustering and marker genes #######################
pcSelect = 20
pbmc <- FindNeighbors(object = pbmc, dims = 1:pcSelect)
pbmc <- FindClusters(object = pbmc, resolution = 0.5)
pbmc <- RunTSNE(object = pbmc, dims = 1:pcSelect, do.fast = TRUE)
pdf(file = "03.TSNE.pdf", width = 6.5, height = 6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)
dev.off()
write.table(pbmc$seurat_clusters, file = "03.tsneCluster.txt", quote = FALSE, sep = "\t", col.names = FALSE)

# Plot cells by sample (orig.ident)
pdf(file = "03.CellCluster-TSNEPlot_TSNE.pdf", width = 6.5, height = 6)
DimPlot(object = pbmc, group.by = "orig.ident", pt.size = 0.5, reduction = "tsne")
dev.off()

# Find marker genes for each cluster
pbmc.markers <- FindAllMarkers(object = pbmc, only.pos = FALSE, min.pct = 0.25, logfc.threshold = logFCfilter)
sig.markers = pbmc.markers[(abs(as.numeric(pbmc.markers$avg_log2FC)) > logFCfilter & pbmc.markers$p_val_adj < adjPvalFilter), ]
write.table(sig.markers, file = "03.clusterMarkers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

top5 <- pbmc.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)

# Heatmap of top markers
pdf(file = "03.tsneHeatmap.pdf", width = 15, height = 15)
DoHeatmap(object = pbmc, features = top5$gene) + NoLegend()
dev.off()

# DotPlot of top markers
pdf(file = "03.tsneDotplot.pdf", width = 15, height = 15)
DotPlot(pbmc, features = unique(top5$gene)) + RotatedAxis()
dev.off()

# Loop through each cluster to plot top 4 genes
marker.sig <- pbmc.markers %>% mutate(Ratio = round(pct.1 / pct.2, 3)) %>% filter(p_val_adj <= 0.05)

for (cluster_id in unique(marker.sig$cluster)) {
  cl4.genes <- marker.sig %>% filter(cluster == cluster_id) %>% arrange(desc(avg_log2FC))
  cl4.genes <- cl4.genes[1:min(nrow(cl4.genes), 4), "gene"]
  
  # Violin plot
  pvn <- VlnPlot(pbmc, features = cl4.genes, ncol = 2)
  pdf(paste0("./MarkerGene-VlnPlot_cluster", cluster_id, "_tsne_PC.pdf"), width = 7, height = 6)
  print(pvn)
  dev.off()
  
  # Feature plot
  pvn <- FeaturePlot(pbmc, features = cl4.genes, ncol = 2)
  pdf(paste0("./MarkerGene-FeaturePlot_cluster", cluster_id, "_tsne_PC.pdf"), width = 7, height = 6)
  print(pvn)
  dev.off()
  
  # Ridge plot
  pvn <- RidgePlot(pbmc, features = cl4.genes, ncol = 2)
  pdf(paste0("./MarkerGene-RidgePlot_cluster", cluster_id, "_tsne_PC.pdf"), width = 7, height = 6)
  print(pvn)
  dev.off()
}
rm(cl4.genes, cluster_id, pvn)

####################### 04. Cell type annotation using SingleR #######################
pbmc_for_SingleR <- GetAssayData(pbmc, slot = "data")
clusters <- pbmc@meta.data$seurat_clusters
ref = get(load("ref_Human_all.RData"))
# Alternative: ref = celldex::HumanPrimaryCellAtlasData()
singler = SingleR(test = pbmc_for_SingleR, ref = ref,
                  labels = ref$label.main, clusters = clusters,
                  method = "cluster", assay.type.test = "logcounts", assay.type.ref = "logcounts")
clusterAnn = as.data.frame(singler)
clusterAnn = cbind(id = row.names(clusterAnn), clusterAnn)
clusterAnn = clusterAnn[, c("id", "labels")]
write.table(clusterAnn, file = "04.clusterAnn.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Add annotations to each cell
cellAnn = c()
for (i in 1:length(pbmc$seurat_clusters)) {
  index = pbmc$seurat_clusters[i]
  cellAnn = c(cellAnn, clusterAnn[index, 2])
}
cellAnnOut = cbind(names(pbmc$seurat_clusters), cellAnn)
colnames(cellAnnOut) = c("id", "labels")
write.table(cellAnnOut, file = "04.cellAnn.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Export normalized expression matrix
expMatrix = as.data.table(pbmc@assays$RNA@layers$data)
rownames(expMatrix) <- rownames(pbmc[["RNA"]]@features)
colnames(expMatrix) <- rownames(pbmc[["RNA"]]@cells)
expMatrix = cbind(id = row.names(expMatrix), expMatrix)
write.table(expMatrix, file = "04.expMatirx.txt", quote = FALSE, sep = "\t", row.names = FALSE)

# Visualization with annotated clusters
newLabels = singler$labels
names(newLabels) = levels(pbmc)
pbmc = RenameIdents(pbmc, newLabels)
pdf(file = "04.TSNE.pdf", width = 7.5, height = 6)
TSNEPlot(object = pbmc, pt.size = 2, label = TRUE)
dev.off()

# Find differentially expressed genes by cell type
pbmc.markers = FindAllMarkers(object = pbmc, only.pos = FALSE, min.pct = 0.25, logfc.threshold = logFCfilter)
sig.cellMarkers = pbmc.markers[(abs(as.numeric(pbmc.markers$avg_log2FC)) > logFCfilter & pbmc.markers$p_val_adj < adjPvalFilter), ]
write.table(sig.cellMarkers, file = "04.cellMarkers.txt", sep = "\t", row.names = FALSE, quote = FALSE)

####################### 05. Visualization of hub genes #######################
# Read hub gene list
geneRT = read.csv(geneFile, header = TRUE, sep = ",", check.names = FALSE)
hubGenes = as.vector(geneRT[, 1])

# Violin plot
pdf(file = "05.hubGeneViolin.pdf", width = 20, height = 16)
VlnPlot(object = pbmc, features = hubGenes)
dev.off()

# Feature plot
pdf(file = "05.hubGeneScatter.pdf", width = 14, height = 10)
FeaturePlot(object = pbmc, features = hubGenes, cols = c("green", "red"))
dev.off()

# Dot plot
pdf(file = "05.hubGeneBubble.pdf", width = 15, height = 8)
DotPlot(object = pbmc, features = hubGenes)
dev.off()
