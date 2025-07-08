
#1.Venn diagram to select the intersection of Diffgenes and Diffpeak genes####
# 1. Venn Diagram: Identify intersection genes between DEGs (RNA-seq) and differential peak-associated genes (ATAC-seq)
# Load the VennDiagram package
library(VennDiagram)

# Set working directory
setwd("G:\\20240429_ATAC_project\\5.diffpeek_diffgene_veen")

# Get a list of .txt files in the directory
myfiles = dir()
myfiles = grep("txt", myfiles, value = TRUE)

# Initialize an empty list to store gene sets
mylist = list()

# Loop through each file and store gene names
for (i in 1:length(myfiles)) {
    onefile = myfiles[i]
    data = read.table(onefile, header = FALSE)
    
    # Extract file name prefix to use as the list name (e.g., "DEG" or "DiffPeek")
    hd = unlist(strsplit(onefile, "\\.|\\-"))
    mylist[[hd[1]]] = as.vector(data[, 1])
    
    # Print gene count in each file
    onel = length(unique(as.vector(data[, 1])))
    print(paste(hd[1], onel, sep = " "))
}

# Take intersection of gene sets from all files
vedata = Reduce(intersect, mylist)

# Save intersection gene list to file
write.table(file = "Venn.txt", vedata, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Draw Venn diagram
venn.plot <- venn.diagram(mylist, fill = c("red", "blue"), filename = NULL)
pdf("venn.pdf")
grid.draw(venn.plot)
dev.off()



# 2. Extract RNA-seq expression matrix for the intersected genes####

# Load necessary packages
library(dplyr)
library(limma)
library(stringi)
library(tidyverse)
library(data.table)
library(impute)

# Set working directory
setwd("G:\\20240429_ATAC_project\\6.SVM+RF+ANN\\12_Veen_acquire")

# Read intersection genes list (from Venn analysis)
data0 = read.table("veen_gene.txt", sep = "\t", check.names = FALSE, quote = "!")

# Extract gene names
m = data0[, 1]

# Read DEG expression matrix (RNA-seq)
data1 = fread("diffGeneExp.txt", sep = "\t", header = TRUE, check.names = FALSE, quote = "!")
data1 = as.data.frame(data1)

# Move gene IDs to rownames
data1 = column_to_rownames(data1, "ID")

# Subset expression data for the intersected genes
veengene = data1[m, ]

# Convert rownames back to a column for saving
veengene = rownames_to_column(veengene, "ID")

# Save the final expression matrix of intersected genes
write.table(veengene, "veengenexpr.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = TRUE)
