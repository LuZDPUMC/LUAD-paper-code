
#TCGA_ATAC_LUAD_dataanlysis####

#1.TCGA_ATAC_LUAD_prepareData####
# Set working directory
setwd("G:\\20240429_ATAC_project\\4.ATAC_LUAD\\03.prepareData")

# Clear environment
rm(list = ls())

# Set options
options(stringsAsFactors = F)

# Define project name
Project = "TCGA-LUAD"  # Example using LUAD, modify according to your cancer type

# Define data paths
ATACDataPath = "LUAD_log2norm.txt"              # ATAC-seq expression data file
IDmapFilePath = "TCGA_identifier_mapping.txt"   # TCGA ID mapping file

###### Data preprocessing

# Read ATAC-seq data
ATACData = read.table(ATACDataPath, header = T, sep = "\t", check.names = F)

# Read TCGA ID mapping file
TCGA_IDmap = read.table(IDmapFilePath, header = T, sep = "\t", check.names = F)

# Filter samples specific to LUAD (or your cancer type)
TCGA_IDmap = TCGA_IDmap[grep("LUAD.*", TCGA_IDmap$bam_prefix), ]

# Extract peak-related columns and rename columns with sample Case_ID
tidyDATA = ATACData[, c("seqnames", "start", "end", "name", "score",
                        gsub(pattern = "-", replacement = "_", TCGA_IDmap$bam_prefix))]
colnames(tidyDATA) = c("seqnames", "start", "end", "name", "score", TCGA_IDmap$Case_ID)

# Extract only peak coordinate columns
peakData = tidyDATA[, 1:5]

# Write peak coordinates to BED file (without column or row names)
write.table(peakData, file = paste0(Project, "-peakData.bed"), sep = "\t",
            col.names = F, row.names = F)

# Load tidyr package
library(tidyr)

# Create peak IDs by merging seqnames and start positions
peak = tidyr::unite(peakData, "ID", seqnames:start, sep = ":")

# Complete peak IDs by adding end positions
peak = tidyr::unite(peak, "ID", ID:end, sep = "-")

# Construct peak matrix (peaks as rows, samples as columns)
peakMatrix = data.frame(ID = peak$ID, tidyDATA[, TCGA_IDmap$Case_ID],
                        check.names = F)

# Save data as an R list object
ATAC_COAD_data <- list(peakData = peakData, peak = peak, peakMatrix = peakMatrix)
save(ATAC_COAD_data, file = "ATAC_LUAD_data.Rdata")




#2.TCGA_ATAC_LUAD_coverage####
# Install genome and karyoploteR packages if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("karyoploteR")

# Set working directory
setwd("E:\\NEW_data\\Project\\lu_tcga_ATACseq_project_LUAD\\61tcgaATAC\\04.coverage")

library(karyoploteR)

# Read ATAC-seq peak file (BED format)
rt = read.table("peak.bed", sep = "\t", header = F)

# Rename columns
colnames(rt) = c("chr", "start", "end", "peak", "value")

# Convert data frame to GRanges object for genomic plotting
data.points = makeGRangesFromDataFrame(rt)

# Assign values (signal intensity) to metadata
mcols(data.points) <- data.frame(y = rt[, 5])

# Create PDF to save the karyotype plot
pdf(file = "coverage2.pdf", width = 10, height = 7)

# Initialize karyotype plot using hg38 reference genome
kp <- plotKaryotype("hg38", plot.type = 1)

# Add background to data panel 1
kpDataBackground(kp, data.panel = 1)

# Plot area tracks representing ATAC-seq signal across genome
kpArea(kp, data = data.points, border = "red", ymin = 0, ymax = 100)

# Add base number ticks to chromosomes
kpAddBaseNumbers(kp, tick.dist = 10000000, minor.tick.dist = 1000000)

# Close the PDF device
dev.off()



#3.TCGA_ATAC_LUAD_features####
# Install necessary Bioconductor packages if not installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ChIPseeker", version = "3.8")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db", version = "3.8")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clusterProfiler", version = "3.8")

# Set working directory
setwd("E:\\NEW_data\\Project\\lu_tcga_ATACseq_project_LUAD\\61tcgaATAC\\05.features")

# Define input peak file
inputFile = "peak.bed"

# Load required libraries for peak annotation and enrichment analysis
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

# Load transcript database for hg38 genome
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate peaks with genomic features and gene annotation database
peakAnno <- annotatePeak(inputFile, tssRegion = c(-3000, 3000),
                         TxDb = txdb, annoDb = "org.Hs.eg.db")

# Plot pie chart of peak annotation distribution and save as PDF
pdf("pie.pdf", height = 7, width = 7)
plotAnnoPie(peakAnno)
dev.off()

# Plot bar chart of peak annotation distribution and save as PDF
pdf("bar.pdf", height = 6, width = 14)
plotAnnoBar(peakAnno)
dev.off()

# Plot venn pie chart showing annotation overlap and save as PDF
pdf("vennpie.pdf", height = 7, width = 8)
vennpie(peakAnno)
dev.off()


#4.TCGA_ATAC_LUAD_upsetplot####

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("org.Hs.eg.db", version = "3.8")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clusterProfiler", version = "3.8")

# Set working directory
setwd("E:\\NEW_data\\Project\\lu_tcga_ATACseq_project_LUAD\\61tcgaATAC\\06.upsetplot")

# Define input peak file
inputFile = "peak.bed"

# Load required libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

# Load transcript database for hg38 genome
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate peaks with genomic features and gene annotation database
peakAnno <- annotatePeak(inputFile, tssRegion = c(-3000, 3000),
                         TxDb = txdb, annoDb = "org.Hs.eg.db")

# Generate and save upset plot showing overlaps of annotated genomic features
pdf("upsetplot1.pdf", height = 10, width = 12)
upsetplot(peakAnno)
dev.off()

# Generate and save upset plot with additional venn pie visualization
pdf("upsetplot2.pdf", height = 10, width = 12)
upsetplot(peakAnno, vennpie = TRUE)
dev.off()


#5.TCGA_ATAC_LUAD_DistToTSS####
# Install required Bioconductor packages if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ChIPseeker", version = "3.8")

# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")
# BiocManager::install("org.Hs.eg.db", version = "3.8")
# BiocManager::install("clusterProfiler", version = "3.8")

# Set the working directory where input and output files are located
setwd("E:\\NEW_data\\Project\\lu_tcga_ATACseq_project_LUAD\\61tcgaATAC\\07.DistToTSS")

# Define the input peak file (BED format)
inputFile = "peak.bed"

# Load required libraries for peak annotation and visualization
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

# Load the transcript annotation database for the human genome (hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate peaks relative to known genomic features
# The transcription start site (TSS) region is defined as -3000 bp upstream to +3000 bp downstream
peakAnno <- annotatePeak(inputFile, 
                         tssRegion = c(-3000, 3000),
                         TxDb = txdb, 
                         annoDb = "org.Hs.eg.db")

# Generate and save the plot showing the distribution of peaks relative to TSS
pdf("DistToTSS.pdf", height = 6, width = 14)
plotDistToTSS(peakAnno,
              title = "Distribution of transcription factor-binding loci\nrelative to TSS")
dev.off()



#6.TCGA_ATAC_LUAD_tagHeatmap####

# Install required Bioconductor packages if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ChIPseeker", version = "3.8")

# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")
# BiocManager::install("org.Hs.eg.db", version = "3.8")
# BiocManager::install("clusterProfiler", version = "3.8")

# Set working directory
setwd("E:\\NEW_data\\Project\\lu_tcga_ATACseq_project_LUAD\\61tcgaATAC\\08.tagHeatmap")

# Define input peak file
inputFile = "peak.bed"

# Load required libraries for peak analysis
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)

# Load the transcript annotation database for the human genome (hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read the peak file (BED format)
peak <- readPeakFile(inputFile)

# Define promoter regions (3 kb upstream and 3 kb downstream of transcription start sites)
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)

# Create a tag matrix representing peak signals across promoter regions
tagMatrix <- getTagMatrix(peak, windows = promoter)

# Generate and save the tag heatmap as a PDF
pdf("tagHeatmap.pdf", height = 12, width = 12)
tagHeatmap(tagMatrix, xlim = c(-3000, 3000), color = "red")
dev.off()





#7.TCGA_ATAC_LUAD_aveProfile####
# Set working directory
setwd("E:\\NEW_data\\Project\\lu_tcga_ATACseq_project_LUAD\\61tcgaATAC\\09.aveProfile")

# Define input peak file
inputFile = "peak.bed"

# Load necessary libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)

# Load human genome annotation database (hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read peak file (BED format)
peak <- readPeakFile(inputFile)

# Define promoter regions as 3 kb upstream and 3 kb downstream of transcription start sites (TSS)
promoter <- getPromoters(TxDb = txdb, upstream = 3000, downstream = 3000)

# Create a tag matrix representing signal intensities over promoter regions
tagMatrix <- getTagMatrix(peak, windows = promoter)

# Plot average profile (smoothed curve of read coverage across promoters) and save as PDF
pdf("aveProfile.pdf", height = 6, width = 10)
plotAvgProf(tagMatrix, xlim = c(-3000, 3000),
            xlab = "Genomic Region (5'->3')", 
            ylab = "Read Count Frequency")
dev.off()

# Plot average profile with 95% confidence intervals by resampling, and save as PDF
pdf("aveProfile.conf.pdf", height = 6, width = 10)
plotAvgProf(tagMatrix, xlim = c(-3000, 3000), conf = 0.95, resample = 1000,
            xlab = "Genomic Region (5'->3')", 
            ylab = "Read Count Frequency")
dev.off()





#8.TCGA_ATAC_LUAD_GO####
# Install required Bioconductor and CRAN packages if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ChIPseeker", version = "3.8")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")
# BiocManager::install("org.Hs.eg.db", version = "3.8")
# BiocManager::install("clusterProfiler", version = "3.8")
# install.packages("ggplot2")

# Set working directory
setwd("E:\\NEW_data\\Project\\lu_tcga_ATACseq_project_LUAD\\61tcgaATAC\\10.GO")

# Define input peak file
inputFile = "peak.bed"

# Load necessary libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Load human genome annotation database (hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read peak file
peak <- readPeakFile(inputFile)

# Annotate peaks with genomic features and gene information (within ±3 kb of TSS)
peakAnno <- annotatePeak(inputFile, tssRegion = c(-3000, 3000),
                         TxDb = txdb, annoDb = "org.Hs.eg.db")

# Extract gene IDs associated with annotated peaks
gene = as.data.frame(peakAnno)$geneId

# Perform Gene Ontology (GO) enrichment analysis on peak-associated genes
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff = 0.05, 
               qvalueCutoff = 0.05,
               ont = "all",        # Use all GO categories: BP, CC, MF
               readable = TRUE)    # Convert gene IDs to gene symbols

# Save GO enrichment results to a text file
write.table(kk, file = "GO.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot bar plot of enriched GO terms and save as PDF
pdf(file = "barplot.pdf", width = 12, height = 7)
barplot(kk, drop = TRUE, showCategory = 10, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scale = 'free')
dev.off()

# Plot bubble (dot) plot of enriched GO terms and save as PDF
pdf(file = "bubble.pdf", width = 12, height = 7)
dotplot(kk, showCategory = 10, split = "ONTOLOGY") + 
  facet_grid(ONTOLOGY ~ ., scale = 'free')
dev.off()




#9.TCGA_ATAC_LUAD_KEGG####
# Install required Bioconductor and CRAN packages if not already installed
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("ChIPseeker", version = "3.8")
# BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")
# BiocManager::install("org.Hs.eg.db", version = "3.8")
# BiocManager::install("clusterProfiler", version = "3.8")
# install.packages("ggplot2")

# Set working directory
setwd("E:\\NEW_data\\Project\\lu_tcga_ATACseq_project_LUAD\\61tcgaATAC\\11.KEGG")

# Define input peak file
inputFile = "peak.bed"

# Load required libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)

# Load human genome annotation database (hg38)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Read ATAC-seq peak file
peak <- readPeakFile(inputFile)

# Annotate peaks with genomic features and gene information (within ±3 kb of TSS)
peakAnno <- annotatePeak(inputFile, tssRegion = c(-3000, 3000),
                         TxDb = txdb, annoDb = "org.Hs.eg.db")

# Perform KEGG pathway enrichment analysis using genes associated with peaks
gene = as.data.frame(peakAnno)$geneId
kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff = 0.05, qvalueCutoff = 0.05)

# Save KEGG enrichment results to a text file
write.table(kk, file = "KEGG.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot bar plot of top enriched KEGG pathways and save as PDF
pdf(file = "barplot.pdf", width = 12, height = 7)
barplot(kk, drop = TRUE, showCategory = 20)
dev.off()

# Plot bubble (dot) plot of top enriched KEGG pathways and save as PDF
pdf(file = "bubble.pdf", width = 12, height = 7)
dotplot(kk, showCategory = 20)
dev.off()



#10.TCGA_ATAC_LUAD_diffpeeks####

##10.1 data_preprocess####
# Set working directory
setwd("G:\\20240429_ATAC_project\\4.ATAC_LUAD\\12.diffpeeks\\1.expr_data\\LUAD_Raw_count\\data_preprocess")

# Read input raw count data
df <- read.table("TCGA_ATAC_LUAD_raw_counts_Step1.txt", header = TRUE,
                 sep = "\t", check.names = FALSE)

# Calculate the mean of every two columns for each row
# Create a result matrix with the same number of rows and half the number of columns
result <- matrix(NA, nrow = nrow(df), ncol = (ncol(df) / 2))

# Loop through pairs of columns and calculate row means
for (i in 1:ncol(result)) {
  result[, i] <- rowMeans(df[, c(2 * i, 2 * i + 1)])
}

# Convert result matrix to data frame
result_df <- as.data.frame(result)

# Get column names of every second column in the original data frame to assign to new columns
m <- colnames(df[, c(seq(2, ncol(df), by = 2))])

# Combine the first column (id column) with the result data frame
exprmean <- data.frame(df[, 1], result_df)

# Assign column names to the new data frame
colnames(exprmean) <- c("id", m)

# Save the processed expression mean data to a file
write.table(exprmean, file = "raw_count_exprmean.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)




##10.2 diffpeeks anlysis####

# Load DESeq2 package
library(DESeq2)

# Set working directory
setwd("G:\\20240429_ATAC_project\\4.ATAC_LUAD\\12.diffpeeks\\3.diffpeeks")

# Read raw count expression data
countData <- read.table("raw_count_exprmean.txt", header = TRUE, row.names = 1, check.names = FALSE)

# Read sample metadata
colData <- read.table("type.txt", header = TRUE)

# Standardize column names (first 12 characters to match sample IDs)
colnames(countData) <- substr(colnames(countData), 1, 12)

# Round expression values to integers (DESeq2 requires integer counts)
countData <- round(countData, 0)

# Check the first few rows of count data
head(countData)

# Filter out rows (genes/peaks) with zero counts across all samples
countData <- countData[rowSums(countData) != 0, ]

# Keep only features expressed in at least 70% of samples
countData <- countData[apply(countData, 1, function(x) sum(x > 0) >= 0.7 * ncol(countData)), ]

# Convert Type column in metadata to a factor
colData$Type <- factor(colData$Type)

# Create a new column grouping samples into "Early" (stage I+II) and "Late" (stage III+IV)
# Here, samples with Type equal to "1" are labeled as "Early", others as "Late"
# Adjust logic based on actual staging information if necessary
colData$Group <- ifelse(colData$Type %in% c("1"), "Early", "Late")

# Create DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ Type)

# Normalize data and perform differential expression analysis
dds <- DESeq(dds)

# Extract results
res <- results(dds)

# Filter differentially expressed features: p-value < 0.05 and |log2FoldChange| > 1
resSig <- res[which(res$pvalue < 0.05 & abs(res$log2FoldChange) > 1), ]

# Add a new column to indicate upregulated or downregulated
resSig[which(resSig$log2FoldChange > 0), 'up_down'] <- 'up'
resSig[which(resSig$log2FoldChange < 0), 'up_down'] <- 'down'

# Save significant results to CSV
write.csv(as.data.frame(resSig), file = "deseq_results_TNM_DIFF.csv")

# Save full results to CSV
write.csv(as.data.frame(res), file = "deseq_results_TNM.csv")




##10.3

# Uncomment the following lines to install required packages if not already installed
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("ChIPseeker", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("org.Hs.eg.db", version = "3.8")

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("clusterProfiler", version = "3.8")

# Set working directory
setwd("G:\\20240429_ATAC_project\\4.ATAC_LUAD\\13.diffpeeks_ann")

# Specify input file containing differential peak regions
inputFile = "diffPeak.txt"

# Load required libraries
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(clusterProfiler)

# Load transcript annotation database for hg38
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene

# Annotate peaks with genomic features, defining promoter regions as -3000 to +3000 relative to TSS
peakAnno <- annotatePeak(inputFile, tssRegion = c(-3000, 3000),
                         TxDb = txdb, annoDb = "org.Hs.eg.db")

# Save the peak annotation results to a tab-separated file
write.table(as.data.frame(peakAnno), file = "diffPeakAnn.xls", sep = "\t", quote = FALSE, row.names = FALSE)

