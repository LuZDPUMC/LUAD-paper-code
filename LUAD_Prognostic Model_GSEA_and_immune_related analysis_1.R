
#1.risk model GSEA analysis####
#if (!requireNamespace("BiocManager", quietly = TRUE))
# install.packages("BiocManager")
# BiocManager::install("limma")
# BiocManager::install("org.Hs.eg.db")
# BiocManager::install("DOSE")
# BiocManager::install("clusterProfiler")
# BiocManager::install("enrichplot")

# Load required packages
library(limma)
library(org.Hs.eg.db)
library(clusterProfiler)
library(enrichplot)
library(data.table)

expFile = "impute.txt"              # Expression data file
riskFile = "risk.train.txt"         # Risk file
gmtFile = "c2.cp.kegg_legacy.v2023.2.Hs.symbols.gmt"   # Gene set file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\20_2.GSEA")    # Set working directory

# Read expression data and preprocess
rt = fread(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)                      # Average duplicates
data = data[rowMeans(data) > 0.5, ]      # Filter low-expression genes

# Remove normal samples
group = sapply(strsplit(colnames(data), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
group = gsub("2", "1", group)
data = data[, group == 0]
data = t(data)
rownames(data) = substr(rownames(data), 1, 12)
data = t(avereps(data))

# Read risk file
Risk = read.table(riskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
colnames(Risk)[ncol(Risk)] <- "Risk"

# Ensure sample consistency
data = data[, row.names(Risk)]

# Compare high vs. low risk groups and calculate logFC
dataL = data[, row.names(Risk[Risk[, "Risk"] == "low", ])]
dataH = data[, row.names(Risk[Risk[, "Risk"] == "high", ])]
meanL = rowMeans(dataL)
meanH = rowMeans(dataH)
meanL[meanL < 0.00001] = 0.00001
meanH[meanH < 0.00001] = 0.00001
logFC = log2(meanH) - log2(meanL)
logFC = sort(logFC, decreasing = TRUE)
genes = names(logFC)

# Read gene set file
gmt = read.gmt(gmtFile)

# Perform GSEA analysis
kk = GSEA(logFC, TERM2GENE = gmt, pvalueCutoff = 1)
kkTab = as.data.frame(kk)
kkTab = kkTab[kkTab$pvalue < 0.05, ]
write.table(kkTab, file = "GSEA.result.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Plot pathways enriched in high-risk group
termNum = 5   # Number of pathways to display
kkUp = kkTab[kkTab$NES > 0, ]
if (nrow(kkUp) >= termNum) {
  showTerm = row.names(kkUp)[1:termNum]
  gseaplot = gseaplot2(kk, showTerm, base_size = 8, title = "Enriched in high risk group")
  pdf(file = "GSEA.highRisk.pdf", width = 7, height = 5.5)
  print(gseaplot)
  dev.off()
}

# Plot pathways enriched in low-risk group
termNum = 5   # Number of pathways to display
kkDown = kkTab[kkTab$NES < 0, ]
if (nrow(kkDown) >= termNum) {
  showTerm = row.names(kkDown)[1:termNum]
  gseaplot = gseaplot2(kk, showTerm, base_size = 8, title = "Enriched in low risk group")
  pdf(file = "GSEA.lowRisk.pdf", width = 7, height = 5.5)
  print(gseaplot)
  dev.off()
}





#2. Prognostic risk model ESTIMATE analysis####
# library(utils)
# rforge <- "http://r-forge.r-project.org"
# install.packages("estimate", repos = rforge, dependencies = TRUE)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")

# Load required libraries
library(limma)
library(estimate)

inputFile = "TCGA-LUAD.txt"      # Expression data file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\20_estimate")    # Set working directory

# Read expression data and preprocess
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)  # Average duplicated gene symbols

# Remove normal samples
group = sapply(strsplit(colnames(data), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
group = gsub("2", "1", group)
data = data[, group == 0]
out = data[rowMeans(data) > 0, ]  # Filter genes with low expression

# Output processed expression data
out = rbind(ID = colnames(out), out)
write.table(out, file = "uniq.symbol.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Run ESTIMATE package
filterCommonGenes(input.f = "uniq.symbol.txt",
                  output.f = "commonGenes.gct",
                  id = "GeneSymbol")

estimateScore(input.ds = "commonGenes.gct",
              output.ds = "estimateScore.gct",
              platform = "illumina")

# Output tumor microenvironment scores for each sample
scores = read.table("estimateScore.gct", skip = 2, header = TRUE)
rownames(scores) = scores[, 1]
scores = t(scores[, 3:ncol(scores)])
rownames(scores) = gsub("\\.", "\\-", rownames(scores))
out = rbind(ID = colnames(scores), scores)
write.table(out, file = "TMEscores.txt", sep = "\t", quote = FALSE, col.names = FALSE)




#3. the TME difference in low/high  prognostic model risk scores####
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")

# install.packages("ggpubr")

# Load required libraries
library(limma)
library(ggpubr)

scoreFile = "TMEscores.txt"    # Tumor Microenvironment (TME) scores file
riskFile = "risk.train.txt"    # Risk file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\21_TMEdiff")   # Set working directory

# Read TME scores and preprocess the data
rt = read.table(scoreFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
data = as.matrix(rt)

# Simplify sample names to first three parts separated by hyphens
rownames(data) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))

# Average duplicates
data = avereps(data)

# Read risk file
risk = read.table(riskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Merge data by intersecting samples
sameSample = intersect(row.names(data), row.names(risk))
data = data[sameSample, , drop = FALSE]
risk = risk[sameSample, "risk", drop = FALSE]
rt = cbind(data, risk)

# Set comparison groups
rt$risk = factor(rt$risk, levels = c("low", "high"))
group = levels(factor(rt$risk))
rt$risk = factor(rt$risk, levels = group)

# Create pairwise comparison list
comp = combn(group, 2)
my_comparisons = list()
for (i in 1:ncol(comp)) {
  my_comparisons[[i]] <- comp[, i]
}

# Loop through TME score columns and plot boxplots
for (i in colnames(rt)[1:3]) {
  boxplot = ggboxplot(rt, x = "risk", y = i, fill = "risk",
                      xlab = "",
                      ylab = i,
                      legend.title = "Risk",
                      palette = c("#264385", "#b22a2a")
  ) + 
    stat_compare_means(comparisons = my_comparisons)
  
  # Output the plot to PDF
  pdf(file = paste0(i, ".pdf"), width = 5, height = 4.5)
  print(boxplot)
  dev.off()
}





#4.Correlation analysis between risk score and immune cell infiltration####
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")

# install.packages("scales")
# install.packages("ggplot2")
# install.packages("ggtext")
# install.packages("tidyverse")
# install.packages("ggpubr")

# Load required libraries
library(limma)
library(scales)
library(ggplot2)
library(ggtext)
library(reshape2)
library(tidyverse)
library(ggpubr)

riskFile = "risk.train.txt"                    # Risk score file
immFile = "infiltration_estimation_for_tcga.csv"   # Immune cell infiltration file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\22_immuneCor")   # Set working directory

# Read risk score data
risk = read.table(riskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Read immune cell infiltration data
immune = read.csv(immFile, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
immune = as.matrix(immune)

# Simplify sample names to first three parts separated by hyphens
rownames(immune) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))

# Average duplicate samples
immune = avereps(immune)

# Get intersecting samples between risk and immune data
sameSample = intersect(row.names(risk), row.names(immune))

# Subset risk and immune data by intersecting samples
risk = risk[sameSample, "riskScore"]
immune = immune[sameSample, ]

# Correlation analysis between risk score and immune cell infiltration
x = as.numeric(risk)
x[x > quantile(x, 0.99)] = quantile(x, 0.99)   # Cap outliers at 99th percentile

outTab = data.frame()
for (i in colnames(immune)) {
  y = as.numeric(immune[, i])
  if (sd(y) < 0.001) { next }                  # Skip if near zero variance
  corT = cor.test(x, y, method = "spearman")   # Spearman correlation test
  cor = corT$estimate
  pvalue = corT$p.value
  if (pvalue < 0.05) {
    outTab = rbind(outTab, cbind(immune = i, cor, pvalue))
  }
}

# Output correlation results
write.table(file = "corResult.txt", outTab, sep = "\t", quote = FALSE, row.names = FALSE)

# Plot bubble chart for significant correlations
corResult = read.table("corResult.txt", header = TRUE, sep = "\t")

# Extract and factor software/source names for coloring
corResult$Software = sapply(strsplit(corResult[, 1], "_"), '[', 2)
corResult$Software = factor(corResult$Software, levels = as.character(unique(corResult$Software[rev(order(as.character(corResult$Software)))])))

b = corResult[order(corResult$Software), ]
b$immune = factor(b$immune, levels = rev(as.character(b$immune)))

# Define colors for plotting based on software/source
colslabels = rep(hue_pal()(length(levels(b$Software))), table(b$Software))

# Output plot to PDF
pdf(file = "correlation.pdf", width = 9, height = 10)
ggplot(data = b, aes(x = cor, y = immune, color = Software)) +
  labs(x = "Correlation coefficient", y = "Immune cell") +
  geom_point(size = 4.1) +
  theme(
    panel.background = element_rect(fill = "white", size = 1.5, color = "black"),
    panel.grid = element_line(color = "grey75", size = 0.5),
    axis.ticks = element_line(size = 0.5),
    axis.text.y = ggtext::element_markdown(colour = rev(colslabels))
  )
dev.off()





#5.Survival Analysis of Immune Cell Infiltration####
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")

# install.packages('survival')
# install.packages('survminer')

# Load required libraries
library(limma)
library(survival)
library(survminer)

cliFile = "LUAD_time.txt"                  # Clinical survival data file
immFile = "infiltration_estimation_for_tcga.csv"    # Immune cell infiltration file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\23_immuneSur")   # Set working directory

# Read immune cell infiltration data
immune = read.csv(immFile, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
immune = as.matrix(immune)
# Simplify sample names to first three parts separated by hyphens
rownames(immune) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(immune))
immune = avereps(immune)  # Average duplicate samples

# Read clinical survival data
time = read.table(cliFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Merge immune and survival data based on intersecting samples
sameSample = intersect(row.names(immune), row.names(time))
immune = immune[sameSample, ]
time = time[sameSample, ]
rt = cbind(time, immune)

# Convert survival time from days to years
rt$futime = rt$futime / 365

# Survival analysis loop for each immune cell type
outTab = data.frame()
for (i in colnames(rt[, 3:ncol(rt)])) {
  # Skip if near zero variance
  if (sd(rt[, i]) < 0.01) { next }
  
  # Group samples based on median immune cell infiltration level
  group = ifelse(rt[, i] > median(rt[, i]), "high", "low")
  
  # Perform survival difference test (log-rank test)
  diff = survdiff(Surv(futime, fustat) ~ group, data = rt)
  pValue = 1 - pchisq(diff$chisq, df = 1)
  
  # For significant results (p < 0.05), plot survival curves
  if (pValue < 0.05) {
    outVector = cbind(i, pValue)
    outTab = rbind(outTab, outVector)
    
    # Format p-value for plot annotation
    if (pValue < 0.001) {
      pValue = "p<0.001"
    } else {
      pValue = paste0("p=", sprintf("%.03f", pValue))
    }
    
    # Fit Kaplan-Meier survival curves
    fit <- survfit(Surv(futime, fustat) ~ group, data = rt)
    
    # Plot survival curve
    surPlot = ggsurvplot(fit,
                         data = rt,
                         conf.int = FALSE,
                         pval = pValue,
                         pval.size = 6,
                         legend.title = i,
                         legend.labs = c("High", "Low"),
                         palette = c("red", "blue"),
                         xlab = "Time(years)",
                         break.time.by = 1,
                         risk.table = FALSE,
                         risk.table.title = "",
                         risk.table.height = 0.35)
    
    # Save plot to PDF
    pdf(file = paste0("sur.", i, ".pdf"), width = 5, height = 4.5, onefile = FALSE)
    print(surPlot)
    dev.off()
  }
}

# Output immune cell survival p-values to a text file
write.table(outTab, file = "immuneSur.result.txt", sep = "\t", row.names = FALSE, quote = FALSE)





#6.LUAD ssGSEA analysis####
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GSVA")

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("GSEABase")

# Load required libraries
library(GSVA)
library(limma)
library(GSEABase)

inputFile = "impute.txt"          # Expression data file
gmtFile = "immune.gmt"            # Immune gene set file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\24.ssGSEA")    # Set working directory

# Read expression data and process input
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
mat = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
mat = avereps(mat)                 # Average duplicate gene symbols
mat = mat[rowMeans(mat) > 0, ]    # Filter out genes with zero mean expression

# Read gene set file
geneSet = getGmt(gmtFile, geneIdType = SymbolIdentifier())

# Perform single sample GSEA (ssGSEA) analysis
ssgseaScore = gsva(mat, geneSet, method = 'ssgsea', kcdf = 'Gaussian', abs.ranking = TRUE)

# Define function to normalize ssGSEA scores to [0,1]
normalize = function(x){
  return((x - min(x)) / (max(x) - min(x)))
}

# Normalize ssGSEA scores
ssgseaOut = normalize(ssgseaScore)
ssgseaOut = rbind(id = colnames(ssgseaOut), ssgseaOut)

# Write normalized ssGSEA scores to file
write.table(ssgseaOut, file = "ssgseaOut.txt", sep = "\t", quote = FALSE, col.names = FALSE)



#7.ssGSEA difference analysis onlow/high prognostic model risk score####
# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")

# install.packages("ggpubr")
# install.packages("reshape2")

# Load necessary libraries
library(limma)
library(ggpubr)
library(reshape2)

riskFile = "risk.train.txt"        # Risk score file
scoreFile = "ssgseaOut.txt"        # ssGSEA score file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\25.ssGSEAdiff")  # Set working directory

# Read ssGSEA results and preprocess data
data = read.table(scoreFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Remove normal samples
group = sapply(strsplit(colnames(data), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
group = gsub("2", "1", group)
data = data[, group == 0, drop = FALSE]

# Format sample names and average duplicates
colnames(data) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data = avereps(t(data))

# Read risk score file
risk = read.table(riskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Merge data by common samples
sameSample = intersect(row.names(data), row.names(risk))
data = data[sameSample, ]
risk = risk[sameSample, ]
rt = cbind(data, risk[, c("riskScore", "risk")])

# Remove the riskScore column (second last column)
rt = rt[, -(ncol(rt) - 1)]

# Define immune cell types to plot
immCell = c("aDCs", "B_cells", "CD8+_T_cells", "DCs", "iDCs", "Macrophages",
            "Mast_cells", "Neutrophils", "NK_cells", "pDCs", "T_helper_cells",
            "Tfh", "Th1_cells", "Th2_cells", "TIL", "Treg")

# Prepare data for immune cell boxplots
rt1 = rt[, c(immCell, "risk")]
data = melt(rt1, id.vars = c("risk"))
colnames(data) = c("Risk", "Type", "Score")
data$Risk = factor(data$Risk, levels = c("low", "high"))

# Plot immune cell infiltration boxplots grouped by risk
p = ggboxplot(data, x = "Type", y = "Score", fill = "Risk",
              xlab = "", ylab = "Score", add = "none",
              palette = c("#264385", "#b22a2a"))
p = p + rotate_x_text(50)

# Save immune cell boxplot to PDF
pdf(file = "immCell.boxplot.pdf", width = 9, height = 6)
p + stat_compare_means(aes(group = Risk),
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "")),
                       label = "p.signif")
dev.off()

# Define immune-related functions to plot
immFunction = c("APC_co_inhibition", "APC_co_stimulation", "CCR",
                "Check-point", "Cytolytic_activity", "HLA", "Inflammation-promoting",
                "MHC_class_I", "Parainflammation", "T_cell_co-inhibition",
                "T_cell_co-stimulation", "Type_I_IFN_Reponse", "Type_II_IFN_Reponse")

# Prepare data for immune function boxplots
rt1 = rt[, c(immFunction, "risk")]
data = melt(rt1, id.vars = c("risk"))
colnames(data) = c("Risk", "Type", "Score")
data$Risk = factor(data$Risk, levels = c("low", "high"))

# Plot immune function boxplots grouped by risk
p = ggboxplot(data, x = "Type", y = "Score", fill = "Risk",
              xlab = "", ylab = "Score", add = "none",
              palette = c("#264385", "#b22a2a"))
p = p + rotate_x_text(50)

# Save immune function boxplot to PDF
pdf(file = "immFunction.boxplot.pdf", width = 9, height = 6)
p + stat_compare_means(aes(group = Risk),
                       symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                          symbols = c("***", "**", "*", "")),
                       label = "p.signif")
dev.off()





#8.Analysis of immune checkpoint differences between high/low-risk groups in LUAD prognosis model####

# if (!requireNamespace("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("limma")

# install.packages("ggplot2")
# install.packages("ggpubr")

# Load required libraries
library(limma)
library(reshape2)
library(ggplot2)
library(ggpubr)

expFile = "impute.txt"           # Gene expression data file
riskFile = "risk.train.txt"      # Risk group file
geneFile = "79gene.txt"          # Immune checkpoint gene list file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\26.checkpoint")  # Set working directory

# Read gene expression data and preprocess
rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[,1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)  # Average duplicate gene entries

# Read immune checkpoint gene list and extract expression
gene = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
sameGene = intersect(row.names(data), as.vector(gene[,1]))
data = t(data[sameGene, ])
data = log2(data + 1)  # Log2 transform expression data

# Remove normal samples based on barcode annotation
group = sapply(strsplit(row.names(data), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
group = gsub("2", "1", group)
data = data[group == 0, ]
rownames(data) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data = avereps(data)

# Merge with risk group data
risk = read.table(riskFile, sep = "\t", header = TRUE, check.names = FALSE, row.names = 1)
sameSample = intersect(row.names(data), row.names(risk))
rt1 = cbind(data[sameSample, ], risk[sameSample, ])
rt1 = rt1[, c(sameGene, "risk")]

# Loop over genes and select those with significant expression differences between risk groups
sigGene = c()
for(i in colnames(rt1)[1:(ncol(rt1) - 1)]){
  if(sd(rt1[, i]) < 0.001) { next }
  wilcoxTest = wilcox.test(rt1[, i] ~ rt1[, "risk"])
  pvalue = wilcoxTest$p.value
  if(pvalue < 0.05){
    sigGene = c(sigGene, i)
  }
}
sigGene = c(sigGene, "risk")
rt1 = rt1[, sigGene]

# Transform data for ggplot2 plotting
rt1 = melt(rt1, id.vars = c("risk"))
colnames(rt1) = c("risk", "Gene", "Expression")

# Set comparison groups
group = levels(factor(rt1$risk))
rt1$risk = factor(rt1$risk, levels = c("low", "high"))
comp = combn(group, 2)
my_comparisons = list()
for(j in 1:ncol(comp)) { my_comparisons[[j]] = comp[, j] }

# Plot boxplots of gene expression between risk groups with statistical comparison
boxplot = ggboxplot(rt1, x = "Gene", y = "Expression", fill = "risk",
                    xlab = "",
                    ylab = "Gene expression",
                    legend.title = "Risk",
                    width = 0.8,
                    #outlier.shape = NA,
                    palette = c("#264385", "#b22a2a")) +
  rotate_x_text(50) +
  stat_compare_means(aes(group = risk),
                     method = "wilcox.test",
                     symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
                                        symbols = c("***", "**", "*", "ns")),
                     label = "p.signif")

# Save plot to PDF
pdf(file = "checkpoint79-1.diff.pdf", width = 13, height = 4.6)
print(boxplot)
dev.off()




#9.Correlation analysis of immune checkpoints between high/low-risk groups in LUAD prognosis model####

library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(data.table)

expFile = "impute.txt"              # Gene expression data file
geneFile = "interGenes.txt"         # Immune checkpoint gene list file
immFile = "checkpointexp.txt"       # Immune cell infiltration data file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\26.checkpoint\\2.COR_79GENE")  # Set working directory

# Read expression data and preprocess
rt = fread(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)
data = data[rowMeans(data) > 0, ]
data = t(data)

# Remove normal samples based on barcode information
group = sapply(strsplit(row.names(data), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
group = gsub("2", "1", group)
data = data[group == 0, ]
rownames(data) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", rownames(data))
data = avereps(data)

# Read immune checkpoint gene list, select genes from expression data
geneRT = read.table(geneFile, header = FALSE, sep = "\t", check.names = FALSE)
data = data[, as.vector(geneRT[, 1])]

# Read immune cell infiltration data
immune = read.table(immFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Find common samples in both datasets and subset accordingly
sameSample = intersect(row.names(data), row.names(immune))
data = data[sameSample, , drop = FALSE]
immune = immune[sameSample, , drop = FALSE]

# Remove last column from immune data (likely metadata or unwanted column)
immune = immune[, -ncol(immune)]

# Initialize dataframe for correlation results
outTab = data.frame()

# Loop through immune cells and genes to calculate Spearman correlations
for(cell in colnames(immune)){
  if(sd(immune[, cell]) == 0) { next }  # Skip if no variance
  for(gene in colnames(data)){
    x = as.numeric(immune[, cell])
    y = as.numeric(data[, gene])
    corT = cor.test(x, y, method = "spearman")
    cor = corT$estimate
    pvalue = corT$p.value
    text = ifelse(pvalue < 0.001, "***", ifelse(pvalue < 0.01, "**", ifelse(pvalue < 0.05, "*", "")))
    outTab = rbind(outTab, cbind(Gene = gene, Immune = cell, cor, text, pvalue))
  }
}

# Write correlation results to file
write.table(file = "corResult_checkpoint.txt", outTab, sep = "\t", quote = FALSE, row.names = FALSE)

# Plot heatmap of correlation coefficients with significance stars
outTab$cor = as.numeric(outTab$cor)

pdf(file = "cor.pdf", width = 12, height = 4)
ggplot(outTab, aes(Immune, Gene)) + 
  geom_tile(aes(fill = cor), colour = "grey", size = 1) +
  scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
  geom_text(aes(label = text), col = "black", size = 3) +
  theme_minimal() +  # Remove background
  theme(
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 8, face = "bold"),  # Rotate x axis labels
    axis.text.y = element_text(size = 8, face = "bold")
  ) +
  labs(fill = paste0("***  p<0.001", "\n", "**  p<0.01", "\n", "*  p<0.05", "\n", "\n", "Correlation")) +  # Legend title
  scale_x_discrete(position = "bottom")  # Place x axis labels at bottom
dev.off()
