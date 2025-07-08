# GTEX_TCGA_dataanlysis

#1.GTEX_TCGA_dataarrange####
# Load required packages
library(dplyr)
library(limma)
library(stringi)
library(tidyverse)
library(data.table)
library(impute)

# Set working directory
setwd("G:\\20240429_ATAC_project\\1.GTEX+TCGA_DIFF\\1.dataarrange")

# Read gene list file
data0 = read.table("arrange.txt", sep = "\t", check.names = F, quote = "!")
m <- data0[,1]

# Read expression matrix file
data1 = fread("XENA-TCGAGTEx-LUAD.txt", sep = "\t", header = T, check.names = F, quote = "!")
data1 <- as.matrix(data1)

# Set row names
rownames(data1) = data1[,1]
geneExp = data1[,2:ncol(data1)]

# Convert to numeric matrix
matrixName = list(rownames(geneExp), colnames(geneExp))
geneExp = matrix(as.numeric(as.matrix(geneExp)), nrow = nrow(geneExp), dimnames = matrixName)

# Keep only proteins expressed in more than half of the samples
geneExp = geneExp[apply(geneExp, 1, function(x) sum(x > 0) >= 0.5 * ncol(geneExp)), ]

# Impute missing values
geneExp1 = impute.knn(geneExp)
geneExp = geneExp1$data

# Remove duplicated proteins (average replicates)
geneExp = avereps(geneExp)

# Check number of missing values
sum(is.na(geneExp))

# #Normalize between arrays (commented out)
# geneExp = normalizeBetweenArrays(as.matrix(geneExp))

# #Log2 transform (commented out)
geneExp = log2(geneExp + 1)

# Add a symbol row on top
geneExp = rbind(symbol = colnames(geneExp), geneExp)

# Select columns matching m
geneExpr <- geneExp[, m]
max(geneExpr)

# Save final matrix
write.table(geneExpr, "geneExpr.txt", sep = "\t", quote = F, col.names = F)




#2.GTEx_TCGA_data_Differential expression analysis####

#tip
# 1.This script is designed to identify differentially expressed genes (DEGs) 
# between normal tissues (n = 347) and tumor tissues (n = 515). The main steps include:

# 2.Data preprocessing: The expression matrix is read, lowly expressed genes are filtered out, 
# and duplicate genes are merged by averaging.
# 
# 3.Differential expression analysis: For each gene, 
# a Wilcoxon rank-sum test is performed to compare expression levels between tumor and normal groups. 
# Mean values, log fold changes (logFC), p-values, and median differences are calculated.
# 
# 4.Filtering significant DEGs: Genes with an absolute logFC greater than 1 and 
# a false discovery rate (FDR) below 0.05 are considered significantly differentially expressed. 
# Both the full results and the filtered DEG lists are saved.
# 
# 5.Heatmap preparation: A gene expression matrix of significant DEGs is prepared 
# and saved for subsequent heatmap visualization.


# Install limma if not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library("limma")

# Set working directory
setwd("G:\\20240429_ATAC_project\\1.GTEX+TCGA_DIFF\\11.diff")   

inputFile = "geneExpr.txt"           # Input file
fdrFilter = 0.05                      # FDR cutoff
logFCfilter = 1                      # logFC cutoff
conNum = 347   # 59+288             # Number of normal samples
treatNum = 515                      # Number of tumor samples

# Read input file
outTab = data.frame()
grade = c(rep(1, conNum), rep(2, treatNum))
rt = read.table(inputFile, sep = "\t", header = T, check.names = F)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[,2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)
data = data[rowMeans(data) > 0.5, ]    # Keep genes with mean expression > 0.5

# Differential expression analysis
for(i in row.names(data)){
  geneName = unlist(strsplit(i, "\\|", ))[1]
  geneName = gsub("\\/", "_", geneName)
  rt_tmp = rbind(expression = data[i, ], grade = grade)
  rt_tmp = as.matrix(t(rt_tmp))
  wilcoxTest <- wilcox.test(expression ~ grade, data = rt_tmp)
  conGeneMeans = mean(data[i, 1:conNum])
  treatGeneMeans = mean(data[i, (conNum+1):ncol(data)])
  logFC = treatGeneMeans - conGeneMeans
  pvalue = wilcoxTest$p.value
  conMed = median(data[i, 1:conNum])
  treatMed = median(data[i, (conNum+1):ncol(data)])
  diffMed = treatMed - conMed
  if(((logFC > 0) & (diffMed > 0)) | ((logFC < 0) & (diffMed < 0))){
    outTab = rbind(outTab, cbind(gene = i, conMean = conGeneMeans, treatMean = treatGeneMeans, logFC = logFC, pValue = pvalue))
  }
}
pValue = outTab[, "pValue"]
fdr = p.adjust(as.numeric(as.vector(pValue)), method = "fdr")
outTab = cbind(outTab, fdr = fdr)

# Save all differential analysis results
write.table(outTab, file = "all.xls", sep = "\t", row.names = F, quote = F)

# Filter significant differential expression results
outDiff = outTab[(abs(as.numeric(as.vector(outTab$logFC))) > logFCfilter & as.numeric(as.vector(outTab$fdr)) < fdrFilter), ]
write.table(outDiff, file = "diff.xls", sep = "\t", row.names = F, quote = F)

# Prepare file for heatmap plotting
heatmap = rbind(ID = colnames(data[as.vector(outDiff[,1]), ]), data[as.vector(outDiff[,1]), ])
write.table(heatmap, file = "diffGeneExp.txt", sep = "\t", col.names = F, quote = F)


#3.Draw heatmap####
## install.packages("pheatmap")  # Uncomment to install if not installed

# Set working directory
setwd("E:\\NEW_data\\Project\\GTEX_tcga\\12.heatmap")   

# Read gene expression data matrix
data = read.table("diffGeneExp.txt", sep = "\t", header = T, row.names = 1, check.names = F)

# Rank transformation: divide rank values into 10 bins and round up
Rank <- ceiling(rank(data, ties.method = "first") / ceiling(nrow(data) * ncol(data) / 10))
dimnames1 = list(rownames(data), colnames(data))
rt = matrix(as.numeric(as.matrix(Rank)), nrow = nrow(data), dimnames = dimnames1)

# Load pheatmap package
library(pheatmap)

# Define sample types (normal vs tumor)
Type = c(rep("con", 347), rep("treat", 526))    # Modify number of normal and tumor samples as needed
names(Type) = colnames(rt)
Type = as.data.frame(Type)

# Output heatmap as PDF
pdf("heatmap2.pdf", height = 15, width = 20)   

# Load color palette
library(RColorBrewer)
coul <- colorRampPalette(brewer.pal(8, "PiYG"))(25)

# Draw heatmap
pheatmap(rt, 
         annotation = Type, 
         show_colnames = F,
         color = colorRampPalette(c("green", "black", "red"))(50),  # Alternative color scheme (commented)
         color = coul,
         cluster_cols = F,
         fontsize = 10,
         fontsize_row = 1,
         fontsize_col = 3)

dev.off()

#4.Draw volcano plot####
library(ggplot2)          

inputFile = "diffGeneExp.txt"    # Input file
outFile = "vol.pdf"        # Output PDF file for volcano plot
logFCfilter = 1            # logFC threshold for filtering
fdrFilter = 0.05           # Adjusted p-value (FDR) threshold
setwd("D:\\biowolf\\Desktop\\bioR\\19.vol")   # Set working directory

# Read data
rt = read.table(inputFile, sep = "\t", header = T, check.names = F)

# Define significance groups
Significant = ifelse((rt$fdr < fdrFilter & abs(rt$logFC) > logFCfilter), 
                     ifelse(rt$logFC > logFCfilter, "Up", "Down"), 
                     "Not")

# Draw volcano plot
p = ggplot(rt, aes(logFC, -log10(fdr))) +
  geom_point(aes(col = Significant)) +
  scale_color_manual(values = c("green", "black", "red")) +
  labs(title = " ") +
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))
p = p + theme_bw()

# Save plot to PDF
pdf(outFile, width = 5.5, height = 4.5)
print(p)
dev.off()
