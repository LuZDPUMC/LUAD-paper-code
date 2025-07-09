
#1. Drug sensitivity analysis between high and low risk scores in LUAD prognosis model####

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install(c("limma", "car", "ridge", "preprocessCore", "genefilter", "sva", "biomaRt"))
#BiocManager::install(c("GenomicFeatures", "maftools", "stringr", "org.Hs.eg.db"))
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")

# install.packages("oncoPredict")

# Load packages
library(limma)
library(oncoPredict)
library(parallel)
set.seed(12345)

expFile = "symbol.txt"    # Expression data file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\27.oncoPredict")   # Set working directory

# Read expression data and preprocess
rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)             # Merge duplicate gene symbols
data = data[rowMeans(data) > 0.5, ]  # Filter out genes with low average expression

# Remove normal samples based on barcode information
group = sapply(strsplit(colnames(data), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
group = gsub("2", "1", group)
data = data[, group == 0]
data = t(data)
rownames(data) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-(.*)", "\\1\\-\\2\\-\\3", rownames(data))
data = avereps(data)       # Average replicates again after sample renaming
data = t(data)

# Load GDSC drug sensitivity data
GDSC2_Expr = readRDS(file = 'GDSC2_Expr.rds')   # Expression data from GDSC2
GDSC2_Res = readRDS(file = 'GDSC2_Res.rds')     # Drug response data from GDSC2
GDSC2_Res = exp(GDSC2_Res)                      # Transform IC50 values from log scale

# Predict drug sensitivity in test samples
calcPhenotype(
  trainingExprData = GDSC2_Expr,          # Training expression data
  trainingPtype = GDSC2_Res,              # Training drug response data
  testExprData = data,                    # Test expression data (TCGA-LUAD samples)
  batchCorrect = 'eb',                    # Batch correction method: empirical Bayes
  powerTransformPhenotype = TRUE,
  removeLowVaryingGenes = 0.2,           # Remove genes with low variability (top 20% retained)
  minNumSamples = 10,                    # Minimum number of samples for model fitting
  printOutput = TRUE,                   # Whether to print prediction results
  removeLowVaringGenesFrom = 'rawData'  # Remove low-variance genes from raw data
)





#2.Boxplot visualization of drug sensitivity analysis results####

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")
#install.packages("ggpubr")

# Load packages
library(limma)
library(ggplot2)
library(ggpubr)

riskFile = "risk.train.txt"          # Risk score file
drugFile = "DrugPredictions.csv"     # Drug sensitivity prediction file
pFilter = 0.001                      # p-value threshold for plotting
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\28.boxplot")   # Set working directory

# Read risk score file
risk = read.table(riskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Read drug sensitivity prediction file
senstivity = read.csv(drugFile, header = TRUE, sep = ",", check.names = FALSE, row.names = 1)
colnames(senstivity) = gsub("(.*)\\_(\\d+)", "\\1", colnames(senstivity))   # Remove replicate numbers in column names

# Merge data based on common samples
sameSample = intersect(row.names(risk), row.names(senstivity))
risk = risk[sameSample, "risk", drop = FALSE]
senstivity = senstivity[sameSample, , drop = FALSE]
rt = cbind(risk, senstivity)

# Define comparison groups
rt$risk = factor(rt$risk, levels = c("low", "high"))
type = levels(factor(rt[,"risk"]))
comp = combn(type, 2)
my_comparisons = list()
for(i in 1:ncol(comp)) { my_comparisons[[i]] <- comp[, i] }

# Loop through each drug and plot boxplot if significant
for(drug in colnames(rt)[2:ncol(rt)]) {
  rt1 = rt[, c(drug, "risk")]
  colnames(rt1) = c("Drug", "Risk")
  rt1 = na.omit(rt1)
  rt1$Drug = log2(rt1$Drug + 1)    # Log2 transform for better visualization
  
  # Wilcoxon test for difference between groups
  test = wilcox.test(Drug ~ Risk, data = rt1)
  diffPvalue = test$p.value
  
  # Plot boxplot if p-value < threshold
  if(diffPvalue < pFilter) {
    boxplot = ggboxplot(rt1, x = "Risk", y = "Drug", fill = "Risk",
                        xlab = "Risk",
                        ylab = paste0(drug, " sensitivity"),
                        legend.title = "Risk",
                        #outlier.shape = NA,
                        palette = c("#264385", "#b22a2a")
    ) +
      stat_compare_means(comparisons = my_comparisons)
    
    # Save plot as PDF
    pdf(file = paste0("drugSenstivity.", drug, ".pdf"), width = 3.5, height = 4.5)
    print(boxplot)
    dev.off()
  }
}




#3.MSI analysis between LUAD prognostic model high/low risk groups####

#install.packages("ggplot2")
#install.packages("ggpubr")

# Load packages
library(plyr)
library(ggplot2)
library(ggpubr)

riskFile = "risk.train.txt"     # Risk score file
cliFile = "MSI.txt"             # Microsatellite instability (MSI) file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\28_2.MSI")   # Set working directory

# Read input files
risk = read.table(riskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
cli = read.table(cliFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Merge samples based on common samples
sameSample = intersect(row.names(risk), row.names(cli))
rt = cbind(risk[sameSample, , drop = FALSE], cli[sameSample, , drop = FALSE])
rt$MSI = factor(rt$MSI, levels = c("MSS", "MSI-L", "MSI-H"))
rt$risk = factor(rt$risk, levels = c("low", "high"))

# Define color palette
bioCol = c("#0066FF", "#FF0000", "#FF9900", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
bioCol = bioCol[1:length(unique(rt[,"MSI"]))]

# Count number of patients in each risk group and MSI type
rt1 = rt[, c("MSI", "risk")]
df = as.data.frame(table(rt1))

# Calculate percentage within each risk group
df = ddply(df, .(risk), transform, percent = Freq / sum(Freq) * 100)

# Position for text labels
df = ddply(df, .(risk), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label = paste0(sprintf("%.0f", df$percent), "%")

# Draw stacked barplot
p = ggplot(df, aes(x = factor(risk), y = percent, fill = MSI)) +
  geom_bar(position = position_stack(), stat = "identity", width = 0.7) +
  scale_fill_manual(values = bioCol) +
  xlab("Risk score") + ylab("Percent weight") + guides(fill = guide_legend(title = "MSI")) +
  geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
  theme_bw()

# Save barplot
pdf(file = "barplot.pdf", width = 4, height = 5)
print(p)
dev.off()

# Define comparison groups for boxplot
rt2 = rt[, c("MSI", "riskScore")]
type = levels(factor(rt2[,"MSI"]))
comp = combn(type, 2)
my_comparisons = list()
for(i in 1:ncol(comp)) { my_comparisons[[i]] <- comp[, i] }

# Draw boxplot comparing risk score among MSI groups
boxplot = ggboxplot(rt2, x = "MSI", y = "riskScore", fill = "MSI",
                    xlab = "",
                    ylab = "Risk score",
                    legend.title = "MSI",
                    palette = bioCol
) + 
  stat_compare_means(comparisons = my_comparisons)

# Save boxplot
pdf(file = "boxplot.pdf", width = 4, height = 4.5)
print(boxplot)
dev.off()






#4.TIDE analysis between LUAD prognostic model high/low risk groups####
library(limma)
library(plyr)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(data.table)
library(impute)

setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\28_3.TIDE\\20240523TCGA")   # Set working directory

# Read gene expression data
data1 = fread("TCGA-LUAD.txt", sep = "\t", header = TRUE, check.names = FALSE, quote = "!")
data1 = as.matrix(data1)
rownames(data1) = data1[,1]
geneExp = data1[,2:ncol(data1)]
matrixName = list(rownames(geneExp), colnames(geneExp))
geneExp = matrix(as.numeric(as.matrix(geneExp)), nrow = nrow(geneExp), dimnames = matrixName)

# Keep genes expressed in at least 50% of samples
geneExp = geneExp[apply(geneExp, 1, function(x) sum(x > 0) >= 0.5 * ncol(geneExp)), ]

# Impute missing values
geneExp1 = impute.knn(geneExp)
geneExp = geneExp1$data

# Remove duplicated genes
geneExp = avereps(geneExp)

# Check number of missing values
sum(is.na(geneExp))

# Remove normal samples
group = sapply(strsplit(colnames(geneExp), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
group = gsub("2", "1", group)
geneExp = geneExp[, group == 0]
colnames(geneExp) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1-\\2-\\3", colnames(geneExp))

# Mean centering normalization
Expr <- t(apply(geneExp, 1, function(x){ x - mean(x) }))
colnames(Expr) <- substr(colnames(Expr), 1, 15)

write.table(Expr, "LUAD_TIDE.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)

### Build TIDE-risk merged data
riskdata <- as.data.frame(fread('risk.train.txt'))
m <- as.numeric(ncol(riskdata))
riskdata <- riskdata[, c(1, m)]

TIDEdata <- as.data.frame(fread('tcga_luad_tide.csv'))
TIDE_risk <- merge(riskdata, TIDEdata, by.x = 1, by.y = 1)
TIDE_risk[, ncol(TIDE_risk) + 1] <- TIDE_risk[, 2]
TIDE_risk <- TIDE_risk[, -2]
colnames(TIDE_risk)[ncol(TIDE_risk)] <- "Risk"

write.table(TIDE_risk, "LUAD_TIDE_risk.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Read merged result file
result <- fread('LUAD_TIDE_risk.txt')
colnames(result)[16] <- "Risk"

# Update Risk column labels
result$Risk <- ifelse(result$Risk == 'low', 'Low_Risk', 'High_Risk')
result$Risk <- factor(result$Risk, levels = c('Low_Risk', 'High_Risk'))

# Define comparison groups
my_comparisons <- list(c("Low_Risk", "High_Risk"))

# 1. TIDE violin plot
p1 <- ggviolin(result, x = 'Risk', y = 'TIDE', fill = 'Risk',
               palette = c("#264385", "#b22a2a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.02, method = 't.test')

# 2. Dysfunction violin plot
p2 <- ggviolin(result, x = 'Risk', y = 'Dysfunction', fill = 'Risk',
               palette = c("#264385", "#b22a2a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.02, method = 't.test')

# 3. Exclusion violin plot
p3 <- ggviolin(result, x = 'Risk', y = 'Exclusion', fill = 'Risk',
               palette = c("#264385", "#b22a2a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.02, method = 't.test')

# 4. MSI violin plot
colnames(result)[6] <- 'MSI'
p4 <- ggviolin(result, x = 'Risk', y = 'MSI', fill = 'Risk',
               palette = c("#264385", "#b22a2a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.02, method = 't.test')

# 5. CAF violin plot
colnames(result)[14] <- 'CAF'
p5 <- ggviolin(result, x = 'Risk', y = 'CAF', fill = 'Risk',
               palette = c("#264385", "#b22a2a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.02, method = 't.test')

# 6. IFNG violin plot
p6 <- ggviolin(result, x = 'Risk', y = "IFNG", fill = 'Risk',
               palette = c("#264385", "#b22a2a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.02, method = 't.test')

# 7. Merck18 violin plot
p7 <- ggviolin(result, x = 'Risk', y = "Merck18", fill = 'Risk',
               palette = c("#264385", "#b22a2a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.02, method = 't.test')

# 8. CD274 violin plot
p8 <- ggviolin(result, x = 'Risk', y = "CD274", fill = 'Risk',
               palette = c("#264385", "#b22a2a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.02, method = 't.test')

# 9. MDSC violin plot
p9 <- ggviolin(result, x = 'Risk', y = "MDSC", fill = 'Risk',
               palette = c("#264385", "#b22a2a"),
               add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.02, method = 't.test')

# 10. TAM M2 violin plot
p10 <- ggviolin(result, x = 'Risk', y = "TAM M2", fill = 'Risk',
                palette = c("#264385", "#b22a2a"),
                add = 'boxplot', add.params = list(fill = "white")) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", 
                     bracket.size = 0.5, tip.length = 0.02, method = 't.test')

# Combine plots
p_plot1 <- p1 + p2 + p3 + p4
p_plot2 <- p5 + p6 + p7 + p8 + p9 + p10

p_plot1
p_plot2




# 5.TCIA analysis between LUAD prognostic model high/low risk groups####

#install.packages("ggpubr")

library(ggpubr)       # Load package

riskFile = "risk.train.txt"               # Risk score file
tciaFile = "TCIA-ClinicalData.tsv"       # TCIA immunotherapy score file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\28_4.TCIA")   # Set working directory

# Read TCIA immunotherapy score file
tcia = read.table(tciaFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
tcia = tcia[, c("ips_ctla4_neg_pd1_neg", "ips_ctla4_neg_pd1_pos", 
                "ips_ctla4_pos_pd1_neg", "ips_ctla4_pos_pd1_pos")]

# Read risk file
risk = read.table(riskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
colnames(risk)[ncol(risk)] <- "Risk"

# Merge data
sameSample = intersect(row.names(tcia), row.names(risk))
tcia = tcia[sameSample, , drop = FALSE]
risk = risk[sameSample, "Risk", drop = FALSE]
data = cbind(tcia, risk)

# Set comparison groups
data$Risk = ifelse(data$Risk == "high", "High-risk", "Low-risk")
group = levels(factor(data$Risk))
data$Risk = factor(data$Risk, levels = c("Low-risk", "High-risk"))
group = levels(factor(data$Risk))
comp = combn(group, 2)
my_comparisons = list()
for (i in 1:ncol(comp)) { my_comparisons[[i]] <- comp[, i] }

# Loop through IPS columns and draw violin plots
for (i in colnames(data)[1:(ncol(data) - 1)]) {
  rt = data[, c(i, "Risk")]
  colnames(rt) = c("IPS", "Risk")
  
  gg1 = ggviolin(rt, x = "Risk", y = "IPS", fill = "Risk", 
                 xlab = "", ylab = i,
                 legend.title = "Risk",
                 palette = c("#264385", "#b22a2a"),
                 add = "boxplot", add.params = list(fill = "white")) + 
    stat_compare_means(comparisons = my_comparisons)
  # Optionally, you can use p-value symbol annotation:
  # stat_compare_means(comparisons = my_comparisons, 
  #    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1),
  #                       symbols = c("***", "**", "*", "ns")), 
  #    label = "p.signif")
  
  # Save violin plot
  pdf(file = paste0(i, ".pdf"), width = 4, height = 4.25)
  print(gg1)
  dev.off()
}

# Save merged data table
write.table(data, "data.txt", sep = "\t", row.names = TRUE, col.names = TRUE)








               