

#1.Feature gene selection and visualization based on random forest analysis####
library(randomForest)
library(ggplot2)  # Load ggplot2 package for plotting

set.seed(123456)
conNum <- 347   # Number of samples in normal group
treatNum <- 515 # Number of samples in LUAD group (or treatment group)
inputFile <- "veengenexpr.txt"       # Input file
setwd("G:\\20240429_ATAC_project\\6.SVM+RF+ANN\\14.RF")      # Set working directory

# Read input file
data1 <- read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
data <- t(data1)  # Transpose data
colnames(data) <- sub("-", "_", colnames(data))
colnames(data) <- sub("-", "_", colnames(data))
group <- c(rep(1, conNum), rep(2, treatNum))  # Define group labels

# Random forest model
rf <- randomForest(as.factor(group) ~ ., data = data, ntree = 500)

# Plot error curve
plot(rf, main = "Random forest", lwd = 2)

# Find the number of trees with the minimum error
optionTrees <- which.min(rf$err.rate[, 1])
optionTrees
rf2 <- randomForest(as.factor(group) ~ ., data = data, ntree = optionTrees)

# Add a red vertical line at the optimal number of trees
abline(v = optionTrees, col = "red", lwd = 2)

dev.off()

# Create a data frame for ggplot visualization of error rates
oob.error.data <- data.frame(
  Trees = rep(1:nrow(rf$err.rate), times = 3),
  Type = rep(c("OOB", "Normal", "LUAD"), each = nrow(rf$err.rate)),
  Error = c(rf$err.rate[, "OOB"],
            rf$err.rate[, "1"],
            rf$err.rate[, "2"]))

# Plot error rates using ggplot
ggplot(oob.error.data, aes(x = Trees, y = Error, color = Type)) +
  geom_line(size = 2) +
  geom_vline(xintercept = optionTrees, linetype = "dashed", color = "red", size = 1.5) +
  labs(title = "Error Rate by Number of Trees", x = "Number of Trees", y = "Error Rate") +
  theme(plot.title = element_text(hjust = 0.5))

# Check gene importance
importance <- importance(x = rf2)

# Plot variable importance (built-in plot)
pdf(file = "geneImportance.pdf", width = 6.2, height = 5.8)
varImpPlot(rf2, main = "")
dev.off()

# Select top genes based on importance
rfGenes <- importance[order(importance[, "MeanDecreaseGini"], decreasing = TRUE), ]
# rfGenes <- names(rfGenes[rfGenes > 5])     # Uncomment if selecting by threshold
rfGenes <- names(rfGenes[1:40])
write.table(rfGenes, file = "rfGenes_40.txt", sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

# Output expression data for selected genes
sigExp <- t(data[, rfGenes])
sigExpOut <- rbind(ID = colnames(sigExp), sigExp)
write.table(sigExpOut, file = "rfGeneExp_40.txt", sep = "\t", quote = FALSE, col.names = FALSE)

# Get feature importance and sort
importance_scores <- importance(rf2)
feature_importance <- data.frame(Gene = rownames(importance_scores), Importance = importance_scores[, "MeanDecreaseGini"])
ordered_features <- feature_importance[order(-feature_importance$Importance), ]

# Get importance based on Gini and accuracy
importance_gini <- importance(rf2, type = 2)
importance_accuracy <- importance(rf2, type = 1)

# Select top 40 important genes
top_genes <- head(ordered_features, 40)

# Plot horizontal barplot for top 40 genes
library(ggplot2)
library(ggthemes)
ggplot(top_genes, aes(x = reorder(Gene, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_few() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) +
  labs(title = "Top 40 Important Genes in Heart Disease Prediction", x = "Gene", y = "Importance") +
  coord_flip()

# Plot with gradient color
color_range <- c("#D8BFD8", "#8B008B")

ggplot(top_genes, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance)) +
  geom_bar(stat = "identity") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) +
  labs(title = "Top 40 Important Genes in Heart Disease Prediction", x = "Gene", y = "Importance") +
  coord_flip()

# Plot with text labels on bars
ggplot(top_genes, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 12), 
        axis.text = element_text(size = 12), 
        axis.title = element_text(size = 15)) +
  labs(title = "Top 40 Important Genes in Heart Disease Prediction", x = "Gene", y = "Importance") +
  coord_flip()

# Circular plot (polar coordinate)
top_genes$Importance <- round(top_genes$Importance, 2)
ggplot(top_genes, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  ylim(-40, 50) +
  geom_text(size = 2.5, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() +
  theme(plot.title = element_text(hjust = 0.5, size = 16), 
        plot.caption = element_text(size = 15), 
        axis.text = element_text(size = 8), 
        axis.title = element_text(size = 15)) +
  labs(title = "Top 40 Important Genes in Heart Disease Prediction", x = "Gene", y = "Importance") +
  coord_polar()

# Step 4: Comprehensive output (not used here because too many genes)

library(gridExtra)

# Retrieve importance scores and sort
importance_scores <- importance(rf_model)
importance_gini <- importance(rf_model, type = 2)
importance_accuracy <- importance(rf_model, type = 1)

feature_importance <- data.frame(Gene = rownames(importance_scores), Importance = importance_scores[, "MeanDecreaseGini"])
ordered_features <- feature_importance[order(-feature_importance$Importance), ]

feature_importance_gini <- data.frame(Gene = rownames(importance_gini), Importance = importance_gini[, "MeanDecreaseGini"])
ordered_features_gini <- feature_importance_gini[order(-feature_importance_gini$Importance), ]

feature_importance_accuracy <- data.frame(Gene = rownames(importance_accuracy), Importance = importance_accuracy[, "MeanDecreaseAccuracy"])
ordered_features_accuracy <- feature_importance_accuracy[order(-feature_importance_accuracy$Importance), ]

color_range <- c("#D8BFD8", "#8B008B")

# Plot importance (Gini)
plot_gini <- ggplot(ordered_features_gini, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.caption = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15)) +
  labs(title = "Feature Importance (Gini Index)", x = "Gene", y = "Importance") +
  coord_flip()

# Plot importance (Accuracy)
plot_accuracy <- ggplot(ordered_features_accuracy, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.caption = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15)) +
  labs(title = "Feature Importance (Accuracy)", x = "Gene", y = "Importance") +
  coord_flip()

# Plot general importance
plot_immitation <- ggplot(ordered_features, aes(x = reorder(Gene, Importance), y = Importance, fill = Importance, label = round(Importance, 2))) +
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5), color = "white") +
  scale_fill_gradient(low = color_range[1], high = color_range[2]) +
  theme_few() +
  theme(plot.title = element_text(hjust = 0.5, size = 16),
        plot.caption = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 15)) +
  labs(title = "Feature Importance", x = "Gene", y = "Importance") +
  coord_flip()

# Arrange plots side by side on one page
grid.arrange(plot_gini, plot_accuracy, plot_immitation, ncol = 3)








#2.Feature gene selection and visualization based on LASSO analysis####
set.seed(123)                        # Set random seed for reproducibility
library(glmnet)                      # Load glmnet package

inputFile = "veengenexpr.txt"       # Input file
setwd("G:\\20240429_ATAC_project\\6.SVM+RF+ANN\\14_1.lasso")  # Set working directory

# Read the input file
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
rt = t(rt)                           # Transpose the expression matrix so samples are rows

conNum <- 347                       # Number of normal samples
treatNum <- 515                    # Number of treatment (disease) samples

group <- c(rep(1, conNum), rep(2, treatNum))  # Define group labels (1 = normal, 2 = disease)

# Build LASSO logistic regression model
x = as.matrix(rt)                  # Convert expression data to matrix format
y = group                          # Define response variable

fit = glmnet(x, y, family = "binomial", alpha = 1)  # Fit LASSO model (alpha = 1 means LASSO)

# Plot coefficient paths along lambda values
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
dev.off()

# Perform cross-validation to determine optimal lambda
cvfit = cv.glmnet(x, y, family = "binomial", alpha = 1, type.measure = 'deviance', nfolds = 10)
pdf(file = "cvfit.pdf", width = 6, height = 5.5)
plot(cvfit)
dev.off()

# Print optimal lambda value
cvfit$lambda.min

# Output selected feature genes
coef = coef(fit, s = cvfit$lambda.min)   # Get coefficients at lambda.min
index = which(coef != 0)                 # Get indices of non-zero coefficients
lassoGene = row.names(coef)[index]       # Get gene names
lassoGene = lassoGene[-1]                # Remove intercept term
write.table(lassoGene, file = "LASSO.gene.txt", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Output expression matrix for selected genes
rt1 = t(rt)                              # Transpose back to gene x sample
lassoexp = rt1[lassoGene, , drop = FALSE]  # Subset expression of selected genes
lassoexp = as.data.frame(lassoexp)
write.table(lassoexp, file = "LASSO.geneExp.txt", sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE)





#3.Construction and Validation of the ANN Predictive Model####
##3.1  Gene Score Calculation####

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)            # Load limma package

expFile = "featuregenexpr.txt"    # Expression data file
diffFile = "diff.txt"             # Differential gene file
setwd("G:\\20240429_ATAC_project\\6.SVM+RF+ANN\\16.ANN\\14.geneScore")    # Set working directory

# Read expression file and preprocess
rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]                 # Use first column as row names
exp = rt[, 2:ncol(rt)]                 # Remove first column (gene ID column)
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)  # Convert to numeric matrix
data = avereps(data)                  # Average expression values if there are duplicate gene names

# Read differential gene file
diffRT = read.table(diffFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
diffRT = diffRT[row.names(data), ]    # Match differential gene info to expression matrix rows

# Gene scoring
dataUp = data[diffRT[, "logFC"] > 0, ]    # Select upregulated genes (logFC > 0)
dataDown = data[diffRT[, "logFC"] < 0, ]  # Select downregulated genes (logFC < 0)

# For upregulated genes: score 1 if expression > median, otherwise 0
dataUp2 = t(apply(dataUp, 1, function(x) ifelse(x > median(x), 1, 0)))
# For downregulated genes: score 1 if expression < median, otherwise 0
dataDown2 = t(apply(dataDown, 1, function(x) ifelse(x > median(x), 0, 1)))

# Output gene score results
outTab = rbind(dataUp2, dataDown2)        # Combine up- and downregulated gene scores
outTab = rbind(id = colnames(outTab), outTab)   # Add column IDs as first row
write.table(outTab, file = "geneScore.txt", sep = "\t", quote = FALSE, col.names = FALSE)



##3.2 Neural Network Construction and Prediction Based on Gene Scores####
# install.packages("neuralnet")
# install.packages("NeuralNetTools")

# Load packages
library(neuralnet)
library(NeuralNetTools)
set.seed(12345678)

conNum = 347    # Number of normal samples
treatNum = 515  # Number of treated samples

inputFile = "geneScore.txt"       # Input file
setwd("G:\\20240429_ATAC_project\\6.SVM+RF+ANN\\16.ANN\\15.neuralNet")      # Set working directory

# Read the input file
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
rt = as.data.frame(t(rt))   # Transpose data to make samples as rows

# Add sample type information (1: normal, 2: treated)
rt$type <- c(rep(1, conNum), rep(2, treatNum))

# Adjust columns to put 'type' as the first column
rt <- rt[, c(ncol(rt), 2:(ncol(rt) - 1))]

data <- rt[, 2:ncol(rt)]

# Extract grouping information
group = rt$type
data$con = ifelse(group == "1", 1, 0)    # Create binary label for normal samples
data$treat = ifelse(group == "2", 1, 0)  # Create binary label for treated samples

# Build neural network model
fit = neuralnet(con + treat ~ ., data, hidden = 5)
fit$result.matrix  # View result matrix
fit$weight         # View weights

# Save weights if needed
# write.table(fit$weight, file="weight.txt")

# Plot network structure
plot(fit)

# Save network structure plot as PDF
pdf(file = "neuralnet.pdf", width = 13, height = 8)
plotnet(fit)
dev.off()

# Predict using the trained model
net.predict = compute(fit, data)$net.result
net.prediction = c("con", "treat")[apply(net.predict, 1, which.max)]
predict.table = table(group, net.prediction)
predict.table

# Calculate accuracy for normal and treated groups
conAccuracy = predict.table[1, 1] / (predict.table[1, 1] + predict.table[1, 2])
treatAccuracy = predict.table[2, 2] / (predict.table[2, 1] + predict.table[2, 2])
paste0("Con accuracy: ", sprintf("%.3f", conAccuracy))
paste0("Treat accuracy: ", sprintf("%.3f", treatAccuracy))

# Output prediction probabilities
colnames(net.predict) = c("con", "treat")
outTab = rbind(id = colnames(net.predict), net.predict)
write.table(outTab, file = "neural.predict.txt", sep = "\t", quote = FALSE, col.names = FALSE)




##3.3 ROC analysis####
# Load the pROC package for ROC curve analysis
library(pROC)

# Input file containing prediction results
inputFile = "neural.predict.txt"
# Input file containing sample group/type information
groupFile = "type2.txt"
# Set working directory
setwd("G:\\20240429_ATAC_project\\6.SVM+RF+ANN\\16.ANN\\16.ROC")

# Read prediction results file, with header, tab-separated, without modifying column names, and using the first column as row names
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
# Read group information file similarly
group = read.table(groupFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
# Extract the 'type' column as the true class labels
y = group$type
# Convert class labels: if label is "1" then 0, else 1 (binary classification labels transformation)
y = ifelse(y == "1", 0, 1)

# Calculate ROC curve object, with true labels y and predicted probabilities from second column of rt (converted to numeric)
roc1 = roc(y, as.numeric(rt[, 2]))
# Calculate 95% confidence interval for AUC using bootstrap method
ci1 = ci.auc(roc1, method = "bootstrap")
# Convert confidence interval result to numeric vector
ciVec = as.numeric(ci1)

# Open PDF device to save plot, with specified size
pdf(file = "ROC.pdf", width = 5, height = 5)
# Plot the ROC curve, print AUC on plot, use red color, and use legacy x-axis for specificity
plot(roc1, print.auc = TRUE, col = "red", legacy.axes = TRUE, main = "Train group")
# Add text on plot to show 95% confidence interval of AUC
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "red")
# Close the PDF device
dev.off()



#4.Test Dataset####
#4.1 Gene Score Calculation for Test Dataset####

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

library(limma)                   # Load limma package for differential expression analysis

expFile = "GSE140343EXPR.txt"          # Expression data file
conFile = "GSE140343_s1.txt"           # Control group sample information file
treatFile = "GSE140343_s2.txt"         # Treatment group sample information file
diffFile = "diff.txt"                  # Differential gene list file
setwd("G:\\20240429_ATAC_project\\6.SVM+RF+ANN\\16.ANN\\17.testGeneScore")  # Set working directory

# Read expression data file
rt = read.table(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]                 # Use first column as row names (gene names)
exp = rt[, 2:ncol(rt)]                 # Extract expression values (exclude gene name column)
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
rt = avereps(data)                     # Average replicate gene expression values

# Log2 transform data if raw values are large
qx = as.numeric(quantile(rt, c(0, 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
LogC = ((qx[5] > 100) || ((qx[6] - qx[1]) > 50 && qx[2] > 0))  # Condition to decide if log transformation is needed
if (LogC) {
  rt[rt < 0] = 0                      # Set any negative values to zero
  rt = log2(rt + 1)                  # Log2 transform with pseudocount 1
}
data = normalizeBetweenArrays(rt)     # Normalize expression data between arrays

# Read sample group information files
con = read.table(conFile, header = FALSE, sep = "\t", check.names = FALSE)
treat = read.table(treatFile, header = FALSE, sep = "\t", check.names = FALSE)
conData = data[, as.vector(con[, 1])]       # Extract control group data columns
treatData = data[, as.vector(treat[, 1])]   # Extract treatment group data columns
data = cbind(conData, treatData)             # Combine control and treatment data
conNum = ncol(conData)                        # Number of control samples
treatNum = ncol(treatData)                    # Number of treatment samples
Type = c(rep("con", conNum), rep("treat", treatNum))  # Sample type labels
colnames(data) = paste0(colnames(data), "_", Type)   # Append type info to column names

# Read differential gene list
diffRT = read.table(diffFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
sameGene = intersect(row.names(data), row.names(diffRT))  # Get genes present in both expression and diff gene list
diffRT = diffRT[sameGene, ]
data = data[sameGene, ]

# Calculate gene scores based on differential expression direction
dataUp = data[diffRT[, "logFC"] > 0, ]       # Genes upregulated in treatment group
dataDown = data[diffRT[, "logFC"] < 0, ]     # Genes downregulated in treatment group
dataUp2 = t(apply(dataUp, 1, function(x) ifelse(x > median(x), 1, 0)))      # For upregulated genes, assign 1 if above median, else 0
dataDown2 = t(apply(dataDown, 1, function(x) ifelse(x > median(x), 0, 1)))  # For downregulated genes, assign 1 if below median, else 0

# Output gene score results
outTab = rbind(dataUp2, dataDown2)
outTab = rbind(id = colnames(outTab), outTab)  # Add column names as first row
write.table(outTab, file = "testGeneScore.txt", sep = "\t", quote = FALSE, col.names = FALSE)  # Write output to file


##4.2.Neural Network Model Training and Prediction on Test Dataset####
#install.packages("neuralnet")
#install.packages("NeuralNetTools")

# Load required packages for neural network modeling and visualization
library(neuralnet)
library(NeuralNetTools)

set.seed(12345678)   # Set random seed for reproducibility

trainFile = "geneScore.txt"        # Input file for training group
testFile = "testGeneScore.txt"     # Input file for test group
inputFile2 = "type2.txt"             # Input file containing group labels
setwd("G:\\20240429_ATAC_project\\6.SVM+RF+ANN\\16.ANN\\18.testNeuralPre")   # Set working directory

# Read training data file, transpose, and convert to data frame
rt = read.table(trainFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
rt = as.data.frame(t(rt))

# Read group labels file
m <- read.table(inputFile2, header = TRUE)
# Add group labels as a new column to training data
rt$type <- m$type
# Reorder columns: put 'type' column first, then the rest of the data columns
rt <- rt[, c(ncol(rt), 2:(ncol(rt) - 1))]

# Extract expression data (excluding the type column)
data <- rt[, 2:ncol(rt)]

# Get group information from 'type' column
group = rt$type
# Create binary indicator columns for each class: 'con' for group "1", 'treat' for group "2"
data$con = ifelse(group == "1", 1, 0)
data$treat = ifelse(group == "2", 1, 0)

# Train neural network model to predict 'con' and 'treat' simultaneously based on input features
fit = neuralnet(con + treat ~ ., data, hidden = 5)

# Read test group gene score data, transpose it to match training data format
data2 = read.table(testFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)
data2 = t(data2)

# Extract group labels from row names of test data (assuming format "sample_group")
group2 = gsub("(.*)\\_(.*)", "\\2", row.names(data2))

# Find common genes/features between training and test data
sameGene = intersect(colnames(data), colnames(data2))
# Subset test data to include only common genes/features
data2 = data2[, sameGene]

# Predict group membership probabilities for test data using trained neural network
net.predict = compute(fit, data2)$net.result

# Convert prediction probabilities to predicted class labels ("con" or "treat") based on max probability
net.prediction = c("con", "treat")[apply(net.predict, 1, which.max)]

# Generate confusion matrix comparing true and predicted groups
predict.table = table(group2, net.prediction)
predict.table

# Calculate prediction accuracy for control group
conAccuracy = predict.table[1, 1] / (predict.table[1, 1] + predict.table[1, 2])
# Calculate prediction accuracy for treatment group
treatAccuracy = predict.table[2, 2] / (predict.table[2, 1] + predict.table[2, 2])
paste0("Con accuracy: ", sprintf("%.3f", conAccuracy))
paste0("Treat accuracy: ", sprintf("%.3f", treatAccuracy))

# Output prediction probabilities for the test group
colnames(net.predict) = c("con", "treat")
outTab = rbind(id = colnames(net.predict), net.predict)
write.table(outTab, file = "test.neuralPredict.txt", sep = "\t", quote = FALSE, col.names = FALSE)



#4.3.neuralDiagnostic_testROC####
#install.packages("pROC")

library(pROC)        # Load the pROC package for ROC curve analysis

inputFile = "test.neuralPredict.txt"    # Input file containing prediction results
setwd("G:\\20240429_ATAC_project\\6.SVM+RF+ANN\\16.ANN\\19.testROC")   # Set working directory

# Read input file, with header, tab-separated, no check on column names, and first column as row names
rt = read.table(inputFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Extract true class labels from row names by splitting string pattern (assumes format like "sample_con" or "sample_treat")
y = gsub("(.*)\\_(.*)", "\\2", row.names(rt))
# Convert class labels: "con" to 0, others to 1 (binary classification)
y = ifelse(y == "con", 0, 1)

# Calculate ROC curve with true labels y and predicted probabilities from second column
roc1 = roc(y, as.numeric(rt[, 2]))
# Compute 95% confidence interval for AUC using bootstrap method
ci1 = ci.auc(roc1, method = "bootstrap")
ciVec = as.numeric(ci1)

# Open PDF device to save ROC plot with specified width and height
pdf(file = "ROC.pdf", width = 5, height = 5)
# Plot ROC curve, print AUC on plot, set color to red, use legacy axes style, and title
plot(roc1, print.auc = TRUE, col = "red", legacy.axes = TRUE, main = "Test group")
# Add text on plot showing 95% confidence interval of AUC at specified coordinates
text(0.39, 0.43, paste0("95% CI: ", sprintf("%.03f", ciVec[1]), "-", sprintf("%.03f", ciVec[3])), col = "red")
# Close the PDF device to write plot to file
dev.off()
