
#1.Merge the TCGA LUAD RNAseq data and clinical survival data.####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("impute")

# Load required packages
library(limma)
library(impute)
library(data.table)

expFile = "TCGA-LUAD.txt"     # Expression data file
cliFile = "LUAD_time.txt"     # Clinical survival data file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\3.1tcgaMergeTime")   # Set working directory

# Read expression data file and organize
rf = fread(expFile, header = TRUE, sep = "\t", check.names = FALSE)
colnames(rf)[1] <- "gene"     # Set first column name to "gene"

exp.pl <- rf

# Set gene names as row names
rownames(exp.pl) = rf$gene

# Remove duplicated genes and keep the first occurrence
exp.pl.1 = distinct(exp.pl, gene, .keep_all = TRUE)

# Set gene names as row names again
rownames(exp.pl.1) = exp.pl.1$gene

# Convert to data frame
rt <- as.data.frame(exp.pl.1)
# Move gene column to row names
rt <- column_to_rownames(rt, "gene")
class(rt)

# # Data imputation for missing values (commented out)
# mat = impute.knn(rt)
# data = mat$data
data <- rt

# Log2 transformation
data = log2(data + 1)

# Output imputed (or preprocessed) data
outData = cbind(id = row.names(data), data)
write.table(outData, file = "impute.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# Remove normal samples
group = sapply(strsplit(colnames(data), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
group = gsub("2", "1", group)
data = data[, group == 0]
colnames(data) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))
data = t(data)
data = avereps(data)

# Read survival data
cli = read.table(cliFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Merge expression data and survival data, and output
sameSample = intersect(row.names(data), row.names(cli))
data = data[sameSample, ]
cli = cli[sameSample, ]
out = cbind(cli, data)
out = cbind(id = row.names(out), out)
write.table(out, file = "expTime.txt", sep = "\t", row.names = FALSE, quote = FALSE)






#2.Merge the GEO LUAD RNAseq data and clinical survival data.####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("impute")

# Load required packages
library(limma)
library(impute)
library(data.table)

expFile = "GSE140343EXPR.txt"     # Expression data file
cliFile = "GSE140343_time.txt"     # Clinical survival data file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\4.1geoMergeTime")   # Set working directory

# Read expression data file and organize
rf = fread(expFile, header = TRUE, sep = "\t", check.names = FALSE)
colnames(rf)[1] <- "gene"      # Set first column name as "gene"

exp.pl <- rf

# Set gene names as row names
rownames(exp.pl) = rf$gene

# Remove duplicated genes and keep the first occurrence
exp.pl.1 = distinct(exp.pl, gene, .keep_all = TRUE)

# Set gene names as row names again
rownames(exp.pl.1) = exp.pl.1$gene

# Convert to data frame
rt <- as.data.frame(exp.pl.1)
# Move gene column to row names
rt <- column_to_rownames(rt, "gene")
class(rt)

# # Data imputation for missing values (commented out)
# mat = impute.knn(rt)
# data = mat$data
data <- rt

# Output imputed (or preprocessed) data
outData = cbind(id = row.names(data), data)
# write.table(outData, file = "impute.txt", sep = "\t", row.names = FALSE, quote = FALSE)

# # Remove normal samples (commented out for this dataset)
# group = sapply(strsplit(colnames(data), "\\-"), "[", 4)
# group = sapply(strsplit(group, ""), "[", 1)
# group = gsub("2", "1", group)
# data = data[, group == 0]
# colnames(data) = gsub("(.*?)\\-(.*?)\\-(.*?)\\-(.*?)\\-.*", "\\1\\-\\2\\-\\3", colnames(data))

# Transpose data so that samples are rows and genes are columns
data = t(data)

# Average replicate gene expression values
data = avereps(data)

# Read clinical survival data
cli = read.table(cliFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Merge expression data and survival data, and output
sameSample = intersect(row.names(data), row.names(cli))
data = data[sameSample, ]
cli = cli[sameSample, ]
out = cbind(cli, data)
out = cbind(id = row.names(out), out)
write.table(out, file = "GSE_expTime.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#3. Batch survival analysis using Kaplan-Meier plots####
# Load required packages
library(survival)
library(survminer)
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\5.K-M\\surivival1")

inputFile = "TCGA_veengenexpr.txt"  # Input data file
outFile = "Survival_result"         # Prefix for output files

# Read data file
rt <- read.table(inputFile, header = TRUE, sep = "\t")

# Replace survival time <= 0 with 1 (to avoid errors)
rt$futime[rt$futime <= 0] = 1

# Convert survival time from days to years
rt$futime <- rt$futime / 365

# Loop through each gene to perform survival analysis
genes <- colnames(rt)[4:ncol(rt)]
for (gene in genes) {
  # Divide samples into high and low expression groups based on median
  rt[[gene]] <- ifelse(rt[[gene]] > median(rt[[gene]]), "high", "low")
  
  # Perform log-rank test
  diff <- survdiff(Surv(futime, fustat) ~ get(gene), data = rt)
  pValue <- 1 - pchisq(diff$chisq, df = 1)
  pValue <- signif(pValue, 4)
  pValue <- format(pValue, scientific = TRUE)
  
  # Fit Kaplan-Meier survival curves
  fit <- survfit(Surv(futime, fustat) ~ get(gene), data = rt)
  
  # Generate survival plot
  surPlot <- ggsurvplot(fit, 
                        data = rt,
                        conf.int = TRUE,                     # Show confidence interval
                        pval = paste0("p=", pValue),        # Add p-value
                        pval.size = 4,
                        risk.table = TRUE,                  # Show risk table
                        legend.labs = c("High_EXP", "Low_EXP"),  # Legend labels
                        legend.title = gene,                # Legend title
                        xlab = "Time (years)",             # X-axis label
                        break.time.by = 1,                 # Interval for x-axis breaks
                        risk.table.title = "",             # Risk table title
                        palette = c("red", "blue"),       # Colors for groups
                        risk.table.height = 0.25)         # Risk table height
  
  # Save plot as PDF
  pdf(file = paste0(outFile, "_", gene, ".pdf"), onefile = FALSE, width = 6.5, height = 5.5)
  print(surPlot)
  dev.off()
}




#4.Univariate Cox regression analysis in batch ####
# Install survival package if not already installed
# install.packages('survival')

# Set working directory
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\6.Unicox")

# Load the survival package
library(survival)

# Read expression and survival data
rt = read.table("TCGA_veengenexpr.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Replace survival time values <= 0 with 1
rt$futime[rt$futime <= 0] = 1

# Convert survival time from days to years
rt$futime = rt$futime / 365

# Initialize an empty data frame to store univariate Cox regression results
outTab = data.frame()

# Perform univariate Cox regression analysis for each gene
for (i in colnames(rt[, 3:ncol(rt)])) {
  cox <- coxph(Surv(futime, fustat) ~ rt[, i], data = rt)
  coxSummary = summary(cox)
  coxP = coxSummary$coefficients[, "Pr(>|z|)"]
  outTab = rbind(outTab,
                 cbind(id = i,
                       HR = coxSummary$conf.int[, "exp(coef)"],
                       HR.95L = coxSummary$conf.int[, "lower .95"],
                       HR.95H = coxSummary$conf.int[, "upper .95"],
                       pvalue = coxSummary$coefficients[, "Pr(>|z|)"]
                 )
  )
}

# Write results to file
write.table(outTab, file = "uniCox.xls", sep = "\t", row.names = FALSE, quote = FALSE)




#5.LASSO Cox Regression Analysis for Feature Selection####

# install.packages("glmnet")
# install.packages("survival")

# Load required packages
library("glmnet")
library("survival")

# Set working directory
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\10.lasso_plot")

# Read input data file
rt = read.table("34uni.SigExp.txt", header = TRUE, sep = "\t", row.names = 1, check.names = FALSE)

# Replace survival time values less than or equal to 0 with 1
rt$futime[rt$futime <= 0] = 1

# Convert survival time from days to years
rt$futime = rt$futime / 365

# Set random seed for reproducibility
set.seed(12345)

# Extract gene expression matrix as numeric matrix
x = as.matrix(rt[, c(3:ncol(rt))])

# Create survival object matrix
y = data.matrix(Surv(rt$futime, rt$fustat))

# Fit LASSO Cox regression model
fit <- glmnet(x, y, family = "cox", maxit = 10000)

# Plot coefficient paths versus log(lambda)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

# Perform cross-validation to determine optimal lambda
cvfit <- cv.glmnet(x, y, family = "cox", maxit = 1000)

# Plot cross-validation curve
pdf("cvfit.pdf")
plot(cvfit)
abline(v = log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty = "dashed")
dev.off()

# Extract coefficients at lambda.min
coef <- coef(fit, s = cvfit$lambda.min)

# Identify non-zero coefficients (selected genes)
index <- which(coef != 0)

# Get actual coefficient values
actCoef <- coef[index]

# Get names of selected genes
lassoGene = row.names(coef)[index]

# Combine survival time, status, and selected genes
lassoGene = c("futime", "fustat", lassoGene)
lassoSigExp = rt[, lassoGene]

# Add sample ID as a separate column
lassoSigExp = cbind(id = row.names(lassoSigExp), lassoSigExp)

# Write selected gene expression data to file
write.table(lassoSigExp, file = "lassoSigExp.txt", sep = "\t", row.names = FALSE, quote = FALSE)



#6.Construction of the LUAD Prognostic Model####
# install.packages("survival")
# install.packages("caret")
# install.packages("glmnet")
# install.packages("survminer")
# install.packages("timeROC")

# Load required packages
library(survival)
library(caret)
library(glmnet)
library(survminer)
library(timeROC)
library(data.table)
library(tidyverse)
set.seed(124)

coxPfilter = 0.01  # p-value cutoff for univariate Cox analysis

# Set working directory
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\11.model")

# Read training dataset
rt = read.table("TCGA_uniSigExp.txt", header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
rt$futime[rt$futime <= 0] = 1/365
# rt$futime = rt$futime / 365
# rt[,3:ncol(rt)] = log2(rt[,3:ncol(rt)] + 1)

# Read external validation dataset
vad = fread("GSE_expTime.txt", header = TRUE, sep = "\t", check.names = FALSE)
vad = as.data.frame(vad)
vad <- column_to_rownames(vad, "id")
vad$futime[vad$futime <= 0] = 1
vad$futime = vad$futime / 365
vad[, 3:ncol(vad)] <- log2(vad[, 3:ncol(vad)] + 1)

############ Function to draw forest plot 
bioForest = function(coxFile = NULL, forestFile = NULL, forestCol = NULL) {
  # Read input Cox result file
  rt <- read.table(coxFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f", rt$"HR")
  hrLow <- sprintf("%.3f", rt$"HR.95L")
  hrLow[hrLow < 0.001] = 0.001
  hrHigh <- sprintf("%.3f", rt$"HR.95H")
  Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")")
  pVal <- ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  # Create forest plot
  pdf(file = forestFile, width = 7, height = 6)
  n <- nrow(rt)
  nRow <- n + 1
  ylim <- c(1, nRow)
  layout(matrix(c(1, 2), nc = 2), width = c(3, 2.5))
  
  # Draw left side (gene info)
  xlim = c(0, 3)
  par(mar = c(4, 2.5, 2, 1))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
  text.cex = 0.8
  text(0, n:1, gene, adj = 0, cex = text.cex)
  text(1.5 - 0.5 * 0.2, n:1, pVal, adj = 1, cex = text.cex)
  text(1.5 - 0.5 * 0.2, n + 1, 'pvalue', cex = text.cex, adj = 1)
  text(3, n:1, Hazard.ratio, adj = 1, cex = text.cex)
  text(3, n + 1, 'Hazard ratio', cex = text.cex, adj = 1)
  
  # Draw forest plot on right side
  par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
  LOGindex = 10
  hrLow = log(as.numeric(hrLow), LOGindex)
  hrHigh = log(as.numeric(hrHigh), LOGindex)
  hr = log(as.numeric(hr), LOGindex)
  xlim = c(floor(min(hrLow, hrHigh)), ceiling(max(hrLow, hrHigh)))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, ylab = "", xaxs = "i", xlab = "Hazard ratio")
  arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
  abline(v = log(1, LOGindex), col = "black", lty = 2, lwd = 2)
  boxcolor = ifelse(as.numeric(hr) > log(1, LOGindex), forestCol[1], forestCol[2])
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.3)
  a1 = axis(1, labels = FALSE, tick = FALSE)
  axis(1, a1, 10^a1)
  dev.off()
}
############ End of forest plot function 


# Data grouping and prognostic model construction
n = 1  # number of random groupings
for (i in 1:n) {
  ############# Data grouping 
  # inTrain <- createDataPartition(y = rt[,2], p = 0.5, list = FALSE)
  train <- rt
  test <- vad
  trainOut = cbind(id = row.names(train), train)
  testOut = cbind(id = row.names(test), test)
  
  # Univariate Cox analysis
  outUniTab = data.frame()
  sigGenes = c("futime", "fustat")
  for (i in colnames(train[, 3:ncol(train)])) {
    cox <- coxph(Surv(futime, fustat) ~ train[, i], data = train)
    coxSummary = summary(cox)
    coxP = coxSummary$coefficients[,"Pr(>|z|)"]
    
    # Select significant genes
    if (coxP < coxPfilter) {
      sigGenes = c(sigGenes, i)
      outUniTab = rbind(outUniTab,
                        cbind(id = i,
                              HR = coxSummary$conf.int[,"exp(coef)"],
                              HR.95L = coxSummary$conf.int[,"lower .95"],
                              HR.95H = coxSummary$conf.int[,"upper .95"],
                              pvalue = coxSummary$coefficients[,"Pr(>|z|)"]))
    }
  }
  uniSigExp = train[, sigGenes]
  if (ncol(uniSigExp) < 5) { next }
  uniSigExpOut = cbind(id = row.names(uniSigExp), uniSigExp)
  
  # LASSO regression
  x = as.matrix(uniSigExp[, c(3:ncol(uniSigExp))])
  y = data.matrix(Surv(uniSigExp$futime, uniSigExp$fustat))
  fit <- glmnet(x, y, family = "cox", maxit = 10000)
  cvfit <- cv.glmnet(x, y, family = "cox", maxit = 10000)
  coef <- coef(fit, s = cvfit$lambda.min)
  index <- which(coef != 0)
  actCoef <- coef[index]
  lassoGene = row.names(coef)[index]
  lassoSigExp = uniSigExp[, c("futime", "fustat", lassoGene)]
  lassoSigExpOut = cbind(id = row.names(lassoSigExp), lassoSigExp)
  geneCoef = cbind(Gene = lassoGene, Coef = actCoef)
  if (nrow(geneCoef) < 2) { next }
  
  ############# Construct multivariate Cox model 
  multiCox <- coxph(Surv(futime, fustat) ~ ., data = lassoSigExp)
  multiCox = step(multiCox, direction = "both")
  multiCoxSum = summary(multiCox)
  
  # Save model coefficients
  outMultiTab = data.frame()
  outMultiTab = cbind(
    coef = multiCoxSum$coefficients[, "coef"],
    HR = multiCoxSum$conf.int[, "exp(coef)"],
    HR.95L = multiCoxSum$conf.int[, "lower .95"],
    HR.95H = multiCoxSum$conf.int[, "upper .95"],
    pvalue = multiCoxSum$coefficients[, "Pr(>|z|)"])
  outMultiTab = cbind(id = row.names(outMultiTab), outMultiTab)
  outMultiTab = outMultiTab[, 1:2]
  
  # Calculate risk score for train set
  riskScore = predict(multiCox, type = "risk", newdata = train)
  coxGene = rownames(multiCoxSum$coefficients)
  coxGene = gsub("`", "", coxGene)
  outCol = c("futime", "fustat", coxGene)
  medianTrainRisk = median(riskScore)
  risk = as.vector(ifelse(riskScore > medianTrainRisk, "high", "low"))
  trainRiskOut = cbind(id = rownames(cbind(train[, outCol], riskScore, risk)),
                       cbind(train[, outCol], riskScore, risk))
  
  # Calculate risk score for test set
  riskScoreTest = predict(multiCox, type = "risk", newdata = test)
  riskTest = as.vector(ifelse(riskScoreTest > medianTrainRisk, "high", "low"))
  testRiskOut = cbind(id = rownames(cbind(test[, outCol], riskScoreTest, riskTest)),
                      cbind(test[, outCol], riskScore = riskScoreTest, risk = riskTest))
  
  # Survival difference p-value
  diff = survdiff(Surv(futime, fustat) ~ risk, data = train)
  pValue = 1 - pchisq(diff$chisq, df = 1)
  diffTest = survdiff(Surv(futime, fustat) ~ riskTest, data = test)
  pValueTest = 1 - pchisq(diffTest$chisq, df = 1)
  
  # ROC AUC
  predictTime = 3
  roc = timeROC(T = train$futime, delta = train$fustat,
                marker = riskScore, cause = 1,
                times = c(predictTime), ROC = TRUE)
  rocTest = timeROC(T = test$futime, delta = test$fustat,
                    marker = riskScoreTest, cause = 1,
                    times = c(predictTime), ROC = TRUE)
  
  if ((pValue < 0.05) & (roc$AUC[2] > 0) & (pValueTest < 0.05) & (rocTest$AUC[2] > 0)) {
    # Output results
    write.table(outUniTab, file = "uni.trainCox.txt", sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(uniSigExpOut, file = "uni.SigExp.txt", sep = "\t", row.names = FALSE, quote = FALSE)
    bioForest(coxFile = "uni.trainCox.txt", forestFile = "uni.foreast.pdf", forestCol = c("red", "green"))
    
    write.table(lassoSigExpOut, file = "lasso.SigExp.txt", sep = "\t", row.names = FALSE, quote = FALSE)
    pdf("lasso.lambda.pdf")
    plot(fit, xvar = "lambda", label = TRUE)
    dev.off()
    pdf("lasso.cvfit.pdf")
    plot(cvfit)
    abline(v = log(c(cvfit$lambda.min, cvfit$lambda.1se)), lty = "dashed")
    dev.off()
    
    write.table(outMultiTab, file = "multiCox.txt", sep = "\t", row.names = FALSE, quote = FALSE)
    write.table(trainRiskOut, file = "risk.train.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    write.table(testRiskOut, file = "risk.test.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    
    allRiskOut = rbind(trainRiskOut, testRiskOut)
    write.table(allRiskOut, file = "risk.all.txt", sep = "\t", quote = FALSE, row.names = FALSE)
    break
  }
}


#7.modelgenes_barplot####

# install.packages("reshape2")
# install.packages("ggpubr")

# Load required packages
library(reshape2)
library(ggpubr)

inputFile = "multiCox.txt"  # Input file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\11_2.modelgenes_barplot")  # Set working directory

# Read input file
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)

# Draw bar plot
pdf(file = "barplot.pdf", width = 5.5, height = 5)
rt$coef = as.numeric(rt$coef)  # Convert coefficient column to numeric
rt$Sig = ifelse(rt$coef > 0, "Positive", "Negative")  # Classify coefficients as positive or negative

# Create bar plot using ggbarplot
gg1 = ggbarplot(rt, x = "id", y = "coef", fill = "Sig", color = "white",
                sort.val = "asc", sort.by.groups = TRUE,
                palette = "jco", 
                rotate = TRUE, title = "", legend = "",
                xlab = "Gene", ylab = "Coefficient", legend.title = "Group", x.text.angle = 60)
print(gg1)
dev.off()





#8.uniSigHeatmap####
#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("pheatmap")

# Load required packages
library(limma)
library(pheatmap)
library(data.table)

expFile = "TCGA-LUAD.txt"           # Expression data file
uniCoxFile = "uni.trainCox.txt"     # Univariate COX analysis results file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\12.uniSigHeatmap")  # Set working directory

# Read expression data file and preprocess the data
rt = fread(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]              # Set first column (gene names) as row names
exp = rt[, 2:ncol(rt)]              # Extract expression matrix (exclude gene names column)
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)                 # Average expression values of duplicated genes

# Read univariate Cox results
uniCox = read.table(uniCoxFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
data = data[row.names(uniCox), ]   # Subset data to genes significant in Cox analysis

# Get the number of normal and tumor samples
group = sapply(strsplit(colnames(data), "\\-"), "[", 4)
group = sapply(strsplit(group, ""), "[", 1)
group = gsub("2", "1", group)      # Replace '2' with '1' to merge sample types (normal=1)
conNum = length(group[group == 1])     # Number of normal samples
treatNum = length(group[group == 0])   # Number of tumor samples
sampleType = c(rep(1, conNum), rep(2, treatNum))  # Sample type vector: 1 = Normal, 2 = Tumor

# Initialize vector to store gene significance labels
sigVec = c()
outTab = data.frame()
for(i in rownames(data)){
  if(sd(data[i, ]) < 0.001) { next }   # Skip genes with near-zero variance
  wilcoxTest = wilcox.test(data[i, ] ~ sampleType)   # Wilcoxon test between Normal and Tumor
  pvalue = wilcoxTest$p.value
  Sig = ifelse(pvalue < 0.001, "***",
               ifelse(pvalue < 0.01, "**",
                      ifelse(pvalue < 0.05, "*", "")))    # Assign significance stars
  sigVec = c(sigVec, paste0(i, Sig))    # Append gene name with significance stars
}

# Log2 transform expression data with a small offset
exp = log2(data + 0.1)
row.names(exp) = sigVec

# Create annotation dataframe for sample types
Type = c(rep("Normal", conNum), rep("Tumor", treatNum))
names(Type) = colnames(data)
Type = as.data.frame(Type)

# Plot heatmap
pdf(file = "heatmap.pdf", width = 8, height = 5)
pheatmap(exp, 
         annotation = Type, 
         color = colorRampPalette(c(rep("blue", 5), "white", rep("red", 5)))(50),
         cluster_cols = FALSE,      # Do not cluster columns (samples)
         cluster_rows = TRUE,       # Cluster rows (genes)
         show_colnames = FALSE,
         show_rownames = TRUE,
         scale = "row",             # Scale gene expression by row (z-score)
         fontsize = 8,
         fontsize_row = 7,
         fontsize_col = 8)
dev.off()





#9.Survival Analysis on Different Risk Groups####
#install.packages("survival")
#install.packages("survminer")

# Load necessary packages
library(survival)
library(survminer)

# Set working directory
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\13.survival")       

# Define Survival Analysis Function

bioSurvival = function(inputFile = NULL, outFile = NULL) {
  # Read input data file
  rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE)
  
  # Compare survival difference between high- and low-risk groups and calculate p-value
  diff = survdiff(Surv(futime, fustat) ~ risk, data = rt)
  pValue = 1 - pchisq(diff$chisq, df = 1)
  if (pValue < 0.001) {
    pValue = "p<0.001"
  } else {
    pValue = paste0("p=", sprintf("%.03f", pValue))
  }
  
  # Fit Kaplan-Meier survival curves
  fit <- survfit(Surv(futime, fustat) ~ risk, data = rt)
  
  # Plot survival curves with confidence intervals and risk table
  surPlot = ggsurvplot(fit, 
                       data = rt,
                       conf.int = TRUE,
                       # conf.int.style = "step",
                       pval = pValue,
                       pval.size = 6,
                       legend.title = "Risk",
                       legend.labs = c("High risk", "Low risk"),
                       xlab = "Time (years)",
                       break.time.by = 1,
                       palette = c("#264385", "#b22a2a"),
                       risk.table = TRUE,
                       risk.table.title = "",
                       risk.table.col = "strata",
                       risk.table.height = 0.25)
  
  # Save the plot as a PDF file
  pdf(file = outFile, width = 9.5, height = 5.5, onefile = FALSE)
  print(surPlot)
  dev.off()
}

# Run Survival Analysis Function on Different Risk Groups

bioSurvival(inputFile = "risk.train.txt", outFile = "surv.train-1.pdf")
bioSurvival(inputFile = "risk.test.txt", outFile = "surv.test-1.pdf")
bioSurvival(inputFile = "risk.all.txt", outFile = "surv.all-1.pdf")




#10.calibration curve on Different Risk Groups####
# install.packages("rms")

library(survival)
library(rms)

dir = "G:\\20240429_ATAC_project\\7.prog_surivival\\14.2_model_Calibration"
setwd(dir)
inputfile = "risk.train.txt"
miRNA <- read.table(inputfile, header = TRUE, sep = "\t", row.names = 1, check.names = FALSE, stringsAsFactors = FALSE)
colnames(miRNA)
# Log2 transform expression columns (columns 3 to 7)
miRNAEXP = log2(miRNA[, 3:7] + 1)
miRNA = cbind(miRNA[, 1:2], miRNAEXP)

# Prepare data distribution object for rms package

ddist <- datadist(miRNA)
options(datadist = 'ddist')

# Specify units for survival time
units(miRNA$futime) <- "Year"

rt <- miRNA

# Plot Calibration Curves

pdf(file = "calibration_1_3_5.pdf", width = 5, height = 5)

# 1-year calibration curve
f2 <- cph(Surv(futime, fustat) ~ TNS4 + RHOV + YWHAZ + CASZ1,
          x = TRUE, y = TRUE, surv = TRUE, data = rt, time.inc = 1)
cal <- calibrate(f2, cmethod = "KM", method = "boot", u = 1, m = (nrow(rt)/3), B = 1000)
plot(cal, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Nomogram-predicted OS (%)", ylab = "Observed OS (%)",
     lwd = 1.5, col = "green", sub = FALSE)

# 3-year calibration curve
f3 <- cph(Surv(futime, fustat) ~ TNS4 + RHOV + YWHAZ + CASZ1,
          x = TRUE, y = TRUE, surv = TRUE, data = rt, time.inc = 3)
cal <- calibrate(f3, cmethod = "KM", method = "boot", u = 3, m = (nrow(rt)/3), B = 1000)
plot(cal, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "", ylab = "", lwd = 1.5, col = "blue", sub = FALSE, add = TRUE)

# 5-year calibration curve
f4 <- cph(Surv(futime, fustat) ~ TNS4 + RHOV + YWHAZ + CASZ1,
          x = TRUE, y = TRUE, surv = TRUE, data = rt, time.inc = 5)
cal <- calibrate(f4, cmethod = "KM", method = "boot", u = 5, m = (nrow(rt)/3), B = 1000)
plot(cal, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "", ylab = "", lwd = 1.5, col = "red", sub = FALSE, add = TRUE)

# Add legend
legend('topleft', c('1-year', '3-year', '5-year'),
       col = c("green", "blue", "red"), lwd = 1.5, bty = 'n')

# Calculate and display C-index with confidence interval

sum.surv = summary(coxph(Surv(futime, fustat) ~ TNS4 + RHOV + YWHAZ + CASZ1, data = rt))
c_index_se = sum.surv$concordance
c_index = c_index_se[1]
c_index.ci_low = sprintf("%.03f", c_index - c_index_se[2])
c_index.ci_high = sprintf("%.03f", c_index + c_index_se[2])
c_index = sprintf("%.03f", c_index)
cindexLabel = paste0(c_index, " (95% CI: ", c_index.ci_low, "-", c_index.ci_high, ")")

# Add text to plot
text(0.5, 0.1, "C-index:")
text(0.7, 0.03, cindexLabel)

dev.off()





#11.Generate risk plots for train, test, and all samples####

# install.packages("pheatmap")

library(pheatmap)

setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\14.riskPlot")


# Define function to plot risk curves, survival status, and heatmap

bioRiskPlot = function(inputFile = NULL, project = NULL) {
  # Read input data file
  rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  
  # Sort samples by patient risk score
  rt = rt[order(rt$riskScore),]
  
  # Plot risk score curve

  riskClass = rt[,"risk"]
  lowLength = length(riskClass[riskClass == "low"])
  highLength = length(riskClass[riskClass == "high"])
  lowMax = max(rt$riskScore[riskClass == "low"])
  line = rt[,"riskScore"]
  line[line > 10] = 10    # Cap risk scores at 10 for plotting
  
  pdf(file = paste0(project, ".riskScore.pdf"), width = 7, height = 4)
  plot(line, type = "p", pch = 20,
       xlab = "Patients (increasing risk score)",
       ylab = "Risk score",
       col = c(rep("#264385", lowLength), rep("#b22a2a", highLength)))
  abline(h = lowMax, v = lowLength, lty = 2)  # Add horizontal and vertical dashed lines
  legend("topleft", c("High risk", "Low Risk"), bty = "n", pch = 19, col = c("#b22a2a", "#264385"), cex = 1.2)
  dev.off()
  

  # Plot survival status

  color = as.vector(rt$fustat)
  color[color == 1] = "#b22a2a"    # Dead samples in red
  color[color == 0] = "#264385"    # Alive samples in blue
  
  pdf(file = paste0(project, ".survStat.pdf"), width = 7, height = 4)
  plot(rt$futime, pch = 19,
       xlab = "Patients (increasing risk score)",
       ylab = "Survival time (years)",
       col = color)
  legend("topleft", c("Dead", "Alive"), bty = "n", pch = 19, col = c("#b22a2a", "#264385"), cex = 1.2)
  abline(v = lowLength, lty = 2)   # Vertical dashed line separating risk groups
  dev.off()
  

  # Define annotation colors for heatmap

  ann_colors = list()
  bioCol = c("#264385", "#b22a2a")
  names(bioCol) = c("low", "high")
  ann_colors[["Risk"]] = bioCol
  
  # Plot heatmap of risk-related gene expression

  rt1 = rt[c(3:(ncol(rt) - 2))]  # Select expression columns only
  rt1 = t(rt1)
  annotation = data.frame(Risk = rt[, ncol(rt)])
  rownames(annotation) = rownames(rt)
  
  pdf(file = paste0(project, ".heatmap.pdf"), width = 6, height = 2.5)
  pheatmap(rt1,
           annotation = annotation,
           annotation_colors = ann_colors,
           cluster_cols = FALSE,
           cluster_rows = FALSE,
           show_colnames = FALSE,
           scale = "row",
           color = colorRampPalette(c(rep("#264385", 3.5), "white", rep("#b22a2a", 3.5)))(50),
           fontsize_col = 3,
           fontsize = 7,
           fontsize_row = 8)
  dev.off()
}

# Generate risk plots for train, test, and all samples

bioRiskPlot(inputFile = "risk.train.txt", project = "train")
bioRiskPlot(inputFile = "risk.test.txt", project = "test")
bioRiskPlot(inputFile = "risk.all.txt", project = "all")







#12.independent prognostic analysis for model####
# install.packages('survival')

library(survival)

setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\15.indep")

# Define function to draw forest plot

bioForest = function(coxFile = NULL, forestFile = NULL, forestCol = NULL){
  # Read input file
  rt <- read.table(coxFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  gene <- rownames(rt)
  hr <- sprintf("%.3f", rt$"HR")
  hrLow  <- sprintf("%.3f", rt$"HR.95L")
  hrHigh <- sprintf("%.3f", rt$"HR.95H")
  Hazard.ratio <- paste0(hr, "(", hrLow, "-", hrHigh, ")")
  pVal <- ifelse(rt$pvalue < 0.001, "<0.001", sprintf("%.3f", rt$pvalue))
  
  # Output plot
  pdf(file = forestFile, width = 6.6, height = 4.5)
  n <- nrow(rt)
  nRow <- n + 1
  ylim <- c(1, nRow)
  layout(matrix(c(1, 2), nc = 2), width = c(3, 2.5))
  
  # Plot clinical info on left side of forest plot
  xlim = c(0, 3)
  par(mar = c(4, 2.5, 2, 1))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, xlab = "", ylab = "")
  text.cex = 0.8
  text(0, n:1, gene, adj = 0, cex = text.cex)
  text(1.5 - 0.5 * 0.2, n:1, pVal, adj = 1, cex = text.cex)
  text(1.5 - 0.5 * 0.2, n + 1, 'pvalue', cex = text.cex, font = 2, adj = 1)
  text(3.1, n:1, Hazard.ratio, adj = 1, cex = text.cex)
  text(3.1, n + 1, 'Hazard ratio', cex = text.cex, font = 2, adj = 1)
  
  # Plot forest plot on right side
  par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
  xlim = c(0, max(as.numeric(hrLow), as.numeric(hrHigh)))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = FALSE, ylab = "", xaxs = "i", xlab = "Hazard ratio")
  arrows(as.numeric(hrLow), n:1, as.numeric(hrHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 2.5)
  abline(v = 1, col = "black", lty = 2, lwd = 2)
  boxcolor = ifelse(as.numeric(hr) > 1, forestCol, forestCol)
  points(as.numeric(hr), n:1, pch = 15, col = boxcolor, cex = 1.5)
  axis(1)
  dev.off()
}

# Define function for independent prognostic analysis

indep = function(riskFile = NULL, cliFile = NULL, project = NULL){
  # Read risk score file
  risk = read.table(riskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  # Read clinical data file
  cli = read.table(cliFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
  
  # Merge data by common samples
  sameSample = intersect(row.names(cli), row.names(risk))
  risk = risk[sameSample, ]
  cli = cli[sameSample, ]
  rt = cbind(futime = risk[, 1], fustat = risk[, 2], cli, riskScore = risk[, (ncol(risk) - 1)])
  
  # Univariate independent prognostic analysis
  uniCoxFile = paste0(project, ".uniCox.txt")
  uniCoxPdf = paste0(project, ".uniCox.pdf")
  uniTab = data.frame()
  for(i in colnames(rt[, 3:ncol(rt)])){
    cox <- coxph(Surv(futime, fustat) ~ rt[, i], data = rt)
    coxSummary = summary(cox)
    uniTab = rbind(uniTab,
                   cbind(id = i,
                         HR = coxSummary$conf.int[, "exp(coef)"],
                         HR.95L = coxSummary$conf.int[, "lower .95"],
                         HR.95H = coxSummary$conf.int[, "upper .95"],
                         pvalue = coxSummary$coefficients[, "Pr(>|z|)"]))
  }
  write.table(uniTab, file = uniCoxFile, sep = "\t", row.names = FALSE, quote = FALSE)
  bioForest(coxFile = uniCoxFile, forestFile = uniCoxPdf, forestCol = "#008B45FF")
  
  # Multivariate independent prognostic analysis
  multiCoxFile = paste0(project, ".multiCox.txt")
  multiCoxPdf = paste0(project, ".multiCox.pdf")
  uniTab = uniTab[as.numeric(uniTab[, "pvalue"]) < 1, ]
  rt1 = rt[, c("futime", "fustat", as.vector(uniTab[, "id"]))]
  multiCox = coxph(Surv(futime, fustat) ~ ., data = rt1)
  multiCoxSum = summary(multiCox)
  multiTab = data.frame()
  multiTab = cbind(
    HR = multiCoxSum$conf.int[, "exp(coef)"],
    HR.95L = multiCoxSum$conf.int[, "lower .95"],
    HR.95H = multiCoxSum$conf.int[, "upper .95"],
    pvalue = multiCoxSum$coefficients[, "Pr(>|z|)"])
  multiTab = cbind(id = row.names(multiTab), multiTab)
  write.table(multiTab, file = multiCoxFile, sep = "\t", row.names = FALSE, quote = FALSE)
  bioForest(coxFile = multiCoxFile, forestFile = multiCoxPdf, forestCol = "#EE0000FF")
}


# Call function to perform independent prognostic analysis

indep(riskFile = "risk.train.txt", cliFile = "clinical.txt", project = "all")




#13.Plot nomogram####
# Load required packages
library(survival)
library(regplot)
library(rms)

# Define input files and set working directory
riskFile = "risk.train.txt"       # Risk score file
cliFile = "clinical.txt"          # Clinical data file
setwd("G:\\20240429_ATAC_project\\7.prog_surivival\\18.Nomo")    # Set working directory

# Read risk score input file
risk = read.table(riskFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
# Replace non-positive survival time with minimum value (1 day)
risk$futime[risk$futime <= 0] = 1 / 365

# Read clinical data file
cli = read.table(cliFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
# Remove samples containing 'unknow' in any clinical field
cli = cli[apply(cli, 1, function(x) !any(is.na(match('unknow', x)))), , drop = FALSE]
cli$Age = as.numeric(cli$Age)
str(cli)

# Merge risk and clinical data by common samples
samSample = intersect(row.names(risk), row.names(cli))
risk1 = risk[samSample, , drop = FALSE]
cli = cli[samSample, , drop = FALSE]
rt = cbind(risk1[, c("futime", "fustat", "risk")], cli)
str(rt)

# Fit Cox proportional hazards model for nomogram
res.cox = coxph(Surv(futime, fustat) ~ risk + Age + Gender + Stage + N + T, data = rt)

# Plot nomogram
nom1 = regplot(res.cox,
               plots = c("density", "boxes"),   # Plot types to display
               clickable = FALSE,
               title = "",                      # Plot title
               points = TRUE,                   # Show points on plot
               droplines = TRUE,                # Show drop lines
               observation = rt[1, ],          # Sample to highlight on plot
               rank = "sd",
               failtime = c(1, 2, 3),          # Time points for survival prediction (years)
               prfail = FALSE)
dev.copy2pdf(file = "Nomo.pdf", width = 8, height = 6, out.type = "pdf")

# Output nomogram risk scores
nomoRisk = predict(res.cox, data = rt, type = "risk")
rt = cbind(risk1, Nomogram = nomoRisk)
outTab = rbind(ID = colnames(rt), rt)
write.table(outTab, file = "nomoRisk.txt", sep = "\t", col.names = FALSE, quote = FALSE)

# Calibration curves
pdf(file = "calibration.pdf", width = 5, height = 5)

# 1-year calibration curve
f <- cph(Surv(futime, fustat) ~ Nomogram, x = TRUE, y = TRUE, surv = TRUE, data = rt, time.inc = 1)
cal <- calibrate(f, cmethod = "KM", method = "boot", u = 1, m = (nrow(rt) / 3), B = 1000)
plot(cal, xlim = c(0, 1), ylim = c(0, 1),
     xlab = "Nomogram-predicted OS (%)", ylab = "Observed OS (%)", lwd = 1.5, col = "green", sub = FALSE)

# 2-year calibration curve
f <- cph(Surv(futime, fustat) ~ Nomogram, x = TRUE, y = TRUE, surv = TRUE, data = rt, time.inc = 2)
cal <- calibrate(f, cmethod = "KM", method = "boot", u = 2, m = (nrow(rt) / 3), B = 1000)
plot(cal, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", lwd = 1.5, col = "blue", sub = FALSE, add = TRUE)

# 3-year calibration curve
f <- cph(Surv(futime, fustat) ~ Nomogram, x = TRUE, y = TRUE, surv = TRUE, data = rt, time.inc = 3)
cal <- calibrate(f, cmethod = "KM", method = "boot", u = 3, m = (nrow(rt) / 3), B = 1000)
plot(cal, xlim = c(0, 1), ylim = c(0, 1), xlab = "", ylab = "", lwd = 1.5, col = "red", sub = FALSE, add = TRUE)

# Add legend for calibration curves
legend('topleft', c('1-year', '2-year', '3-year'),
       col = c("green", "blue", "red"), lwd = 1.5, bty = 'n')

# Calculate concordance index (C-index) for model performance
sum.surv = summary(coxph(Surv(futime, fustat) ~ ., data = rt))
c_index_se = sum.surv$concordance
c_index = c_index_se[1]
c_index.ci_low = sprintf("%.03f", c_index - c_index_se[2])
c_index.ci_high = sprintf("%.03f", c_index + c_index_se[2])
c_index = sprintf("%.03f", c_index)
cindexLabel = paste0(c_index, " (95% CI: ", c_index.ci_low, "-", c_index.ci_high, ")")

# Add C-index text to calibration plot
text(0.5, 0.1, "C-index:")
text(0.7, 0.03, cindexLabel)

dev.off()





