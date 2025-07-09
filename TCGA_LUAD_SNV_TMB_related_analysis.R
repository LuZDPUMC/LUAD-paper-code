
#1.SNV DATA download#####

setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\19_snv_data\\1.martix_data")

# This function has two parameters:
# @metadata: sample sheet file downloaded from TCGA
# @path: the directory path where MAF files are stored

merge_maf <- function(metadata, path) {
  # Construct the full file paths by combining path, file_id, and file_name columns
  filenames <- file.path(path, metadata$file_id, metadata$file_name, 
                         fsep = .Platform$file.sep)
  
  message('############### Merging maf data ################\n',
          '### This step may take a few minutes ###\n')
  
  # Read each sample's MAF file using lapply and combine them by row (rbind)
  # colClasses is set to "character" to ensure all columns are read as character type
  mafMatrix <- do.call("rbind", lapply(filenames, function(fl) 
    read.table(gzfile(fl), header = TRUE, sep = "\t", quote = "", fill = TRUE, colClasses = "character")))
  
  return(mafMatrix)
}

# Define a function to remove duplicate samples
FilterDuplicate <- function(metadata) {
  filter <- which(duplicated(metadata[,'sample']))
  if (length(filter) != 0) {
    metadata <- metadata[-filter, ]
  }
  message(paste('Removed', length(filter), 'samples', sep = ' '))
  return(metadata)
}

# Read the sample sheet file for MAF
metaMatrix.maf = read.table("maf_sample_sheet.tsv", sep = "\t", header = TRUE)

# Replace '.' with '_' and convert all column names to lowercase; also replace 'sample_id' with 'sample'
names(metaMatrix.maf) = gsub("sample_id", "sample", gsub("\\.", "_", tolower(names(metaMatrix.maf))))

# Remove spaces in the last column sample_type
metaMatrix.maf$sample_type = gsub(" ", "", metaMatrix.maf$sample_type)

# Remove duplicate samples
metaMatrix.maf <- FilterDuplicate(metaMatrix.maf)

# Merge MAF matrices using merge_maf function
maf_value = merge_maf(metadata = metaMatrix.maf, 
                      path = "maf_data")

# View first 3 rows and first 10 columns
maf_value[1:3, 1:10]

# Save the merged MAF matrix
write.table(file = "combined_maf_value.txt", maf_value, row.names = FALSE, quote = FALSE, sep = "\t")






#2.Calculate TMB (tumor mutational burden) values####
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\19_snv_data\\2.TMB_Score")

# Install maftools package if not already installed
# BiocManager::install("maftools")
library(maftools)

# Read the merged MAF file
laml <- read.maf(maf = "combined_maf_value.txt")

# Calculate TMB (tumor mutational burden) values
tmb_table_wt_log = tmb(maf = laml)

# Save the calculated TMB values
write.table(tmb_table_wt_log, file = "TMB_log.txt", sep = "\t", row.names = FALSE)



#3.Predictive Model riskTMB####

# install.packages("ggpubr")

# Load required packages
library(ggpubr)
library(reshape2)
library(tidyverse)

RiskFile = "risk.train.txt"    # Risk file
tmbFile = "TMB.txt"           # Tumor mutational burden (TMB) file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\19_snv_data\\6.riskTMB")    # Set working directory

# Read input files
rt = read.table(tmbFile, header = TRUE, sep = "\t", check.names = FALSE)      # Read TMB file
Risk = read.table(RiskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)  # Read risk file

# Extract required columns
tmb <- rt[, c(1, 3)]

# Truncate sample barcode to first 12 characters
tmb$Tumor_Sample_Barcode <- substr(tmb$Tumor_Sample_Barcode, 1, 12)

# Remove duplicate rows
tmb <- tmb[!duplicated(tmb$Tumor_Sample_Barcode), ]

# Set Tumor_Sample_Barcode column as row names
rownames(tmb) <- tmb$Tumor_Sample_Barcode

# Convert to matrix and remove Tumor_Sample_Barcode column
m <- as.matrix(tmb[, -1])
rownames(m) <- rownames(tmb)
tmb <- m
colnames(tmb) <- "TMB"

# Merge data
tmb = as.matrix(tmb)
tmb = log2(tmb + 1)

# Find common samples
sameSample = intersect(row.names(tmb), row.names(Risk))
tmb = tmb[sameSample, , drop = FALSE]
Risk = Risk[sameSample, , drop = FALSE]

# Combine into one data frame
data = cbind(Risk, tmb)
data = data[, c("riskScore", "risk", "TMB")]

# Cap riskScore values above 99th percentile
data[,"riskScore"][data[,"riskScore"] > quantile(data[,"riskScore"], 0.99)] = quantile(data[,"riskScore"], 0.99)

# Set group comparison
data$risk = factor(data$risk, levels = c("low", "high"))
Risk = levels(factor(data$risk))
comp = combn(Risk, 2)
my_comparisons = list()
for (i in 1:ncol(comp)) {
  my_comparisons[[i]] <- comp[, i]
}

# Draw boxplot comparing TMB between risk groups
boxplot = ggboxplot(data, x = "risk", y = "TMB", fill = "risk",
                    xlab = "",
                    ylab = "Tumor Mutational Burden",
                    legend.title = "Risk",
                    palette = c("#264385", "#b22a2a")) + 
  stat_compare_means(comparisons = my_comparisons)

pdf(file = "boxplot.pdf", width = 5, height = 4.25)
print(boxplot)
dev.off()

# Correlation plot between risk score and TMB
p1 = ggplot(data, aes(riskScore, TMB)) + 
  xlab("Risk score") + ylab("Tumor Mutational Burden") +
  geom_point(aes(colour = risk)) +
  scale_color_manual(values = c("#264385", "#b22a2a")) +
  geom_smooth(method = "lm", formula = y ~ x) + 
  theme_bw() +
  stat_cor(method = 'spearman', aes(x = riskScore, y = TMB))

pdf(file = "cor.pdf", width = 5.5, height = 4.25)
print(p1)
dev.off()





#4.TCGA-LUAD-SNP-related plot##### 
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\19_snv_data\\TME")

# Load required libraries
library(TCGAbiolinks)
library(wordcloud2)

# Query mutation data from TCGA-LUAD
query <- GDCquery(
  project = "TCGA-LUAD", 
  data.category = "Simple Nucleotide Variation",
  data.type = "Masked Somatic Mutation",
  access = "open"
)

# Download data
GDCdownload(query)

# Prepare and save the data
GDCprepare(query, save = TRUE, save.filename = "TCGA-LUAD_SNP.Rdata")

library(maftools)

# Load the prepared mutation data
load(file = "TCGA-LUAD_SNP.Rdata")

maf.coad <- data

# Save the raw mutation data as a table
write.table(data, "LUAD_TMB.txt", sep = "\t", quote = FALSE, row.names = FALSE)

# Check the class of the data
class(maf.coad)
## [1] "data.frame"

# Check the dimensions of the data
dim(maf.coad)
## [1] 252664    141

# Read the MAF object
maf <- read.maf(maf.coad)
## -Validating
## -Silent variants: 63597 
## -Summarizing
## --Multiple centers found
## BCM;WUGSC;BCM;WUGSC;BCM;BI -- Possible FLAGS among top ten genes:
##   TTN
##   SYNE1
##   MUC16
## -Processing clinical data
## --Missing clinical data
## -Finished

# Summary plot
pdf(file = "summary.pdf", width = 7, height = 6)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
dev.off()

# Oncoplot (waterfall plot)
pdf(file = "waterfall.pdf", width = 7, height = 6)
oncoplot(maf = maf, top = 30, fontSize = 0.8, showTumorSampleBarcodes = FALSE)
dev.off()

# Somatic interaction plot
pdf(file = "interaction.pdf", width = 7, height = 6)
somaticInteractions(maf = maf, top = 20, pvalue = c(0.05, 0.001), showCounts = FALSE, fontSize = 0.6)
dev.off()

# Gene cloud plot
pdf(file = "Genecloud.pdf", width = 7, height = 6)
geneCloud(input = maf, minMut = 30)
dev.off()

# Calculate TMB values
coad.tmb <- tmb(maf, captureSize = 38, logScale = TRUE)
dim(coad.tmb)

head(coad.tmb)

# Save TMB table
write.table(coad.tmb, "TMB.txt", sep = "\t", quote = FALSE, row.names = FALSE)




#5.draw survival curve for combined TMB and risk####
# install.packages("survival")
# install.packages("survminer")

# Load required libraries
library(survival)
library(survminer)

riskFile = "risk.train.txt"    # Risk file
tmbFile = "TMB.txt"           # Tumor mutational burden (TMB) file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\19_snv_data\\7.tmbSur")    # Set working directory

# Read input files
risk = read.table(riskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)  # Read risk file
rt = read.table(tmbFile, header = TRUE, sep = "\t", check.names = FALSE)                   # Read TMB file

# Extract required columns
tmb <- rt[, c(1, 3)]

# Truncate sample barcodes to the first 12 characters
tmb$Tumor_Sample_Barcode <- substr(tmb$Tumor_Sample_Barcode, 1, 12)

# Remove duplicate rows
tmb <- tmb[!duplicated(tmb$Tumor_Sample_Barcode), ]

# Set Tumor_Sample_Barcode as row names
rownames(tmb) <- tmb$Tumor_Sample_Barcode

# Convert to matrix and remove Tumor_Sample_Barcode column
m <- as.matrix(tmb[, -1])
rownames(m) <- rownames(tmb)
tmb <- m
colnames(tmb) <- "TMB"

# Merge data
sameSample = intersect(row.names(tmb), row.names(risk))
tmb = tmb[sameSample, , drop = FALSE]
risk = risk[sameSample, , drop = FALSE]
data = cbind(risk, tmb)

# Determine optimal cutoff value for TMB
res.cut = surv_cutpoint(data, time = "futime", event = "fustat", variables = c("TMB"))
cutoff = as.numeric(res.cut$cutpoint[1])
tmbType = ifelse(data[, "TMB"] <= cutoff, "L-TMB", "H-TMB")
scoreType = ifelse(data$risk == "low", "low risk", "high risk")
mergeType = paste0(tmbType, "+", scoreType)

# Define survival analysis plotting function
bioSurvival = function(surData = NULL, outFile = NULL) {
  diff = survdiff(Surv(futime, fustat) ~ group, data = surData)
  length = length(levels(factor(surData[, "group"])))
  pValue = 1 - pchisq(diff$chisq, df = length - 1)
  if (pValue < 0.001) {
    pValue = "p<0.001"
  } else {
    pValue = paste0("p=", sprintf("%.03f", pValue))
  }
  fit <- survfit(Surv(futime, fustat) ~ group, data = surData)
  
  # Define color palette
  bioCol = c("#b22a2a", "#264385", "#6E568C", "#7CC767", "#223D6C", "#D20A13", "#FFD121", "#088247", "#11AA4D")
  bioCol = bioCol[1:length]
  
  # Draw survival curve
  surPlot = ggsurvplot(fit,
                       data = surData,
                       conf.int = FALSE,
                       pval = pValue,
                       pval.size = 6,
                       legend.title = "",
                       legend.labs = levels(factor(surData[, "group"])),
                       font.legend = 10,
                       legend = c(0.8, 0.8),
                       xlab = "Time (years)",
                       break.time.by = 1,
                       palette = bioCol,
                       risk.table = FALSE,
                       cumevents = FALSE,
                       risk.table.height = 0.25)
  
  # Save plot
  pdf(file = outFile, width = 5.5, height = 4.8, onefile = FALSE)
  print(surPlot)
  dev.off()
}

# Call function to draw survival curve for TMB alone
data$group = tmbType
bioSurvival(surData = data, outFile = "TMB.survival.pdf")

# Call function to draw survival curve for combined TMB and risk
data$group = mergeType
bioSurvival(surData = data, outFile = "TMB-risk.survival.pdf")


