#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("ComplexHeatmap")

#install.packages("dplyr")
#install.packages("NMF")
#install.packages("patchwork")
#install.packages("igraph")
#install.packages("ggplot2")
#install.packages("ggalluvial")
#install.packages("circlize")
#install.packages("svglite")

#install.packages("devtools")
#library(devtools)
#devtools::install_github("sqjin/CellChat")



# Load required libraries
library(limma)
library(NMF)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(CellChat)
library(data.table)

expFile = "20240518_expMatrix.txt"    # Expression data file
annFile = "04.cellAnn.txt"             # Cell annotation file
geneFile = "Genes.csv"                 # Gene list file
setwd("D:\\lu_singlecell\\3.cellchat")  # Set working directory

# Read expression data file and process the data
rt = fread(expFile, header = TRUE, sep = "\t", check.names = FALSE)
rt = as.matrix(rt)
rownames(rt) = rt[, 1]
exp = rt[, 2:ncol(rt)]
dimnames = list(rownames(exp), colnames(exp))
data = matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)
data = avereps(data)                 # Average expression values for duplicated genes
data = data[rowMeans(data) > 0, ]   # Filter genes with average expression > 0

# Read cell annotation file
meta = read.table(annFile, header = TRUE, sep = "\t", check.names = FALSE, row.names = 1)

# Create CellChat object
cellchat <- createCellChat(object = data, meta = meta, group.by = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))    # Number of cells in each cell type

# Load ligand-receptor interaction database
CellChatDB <- CellChatDB.human        # For mouse data, replace with CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")  # Subset secreted signaling interactions
cellchat@DB <- CellChatDB.use

# Visualize types of ligand-receptor pairs in the database
pdf(file = "COMM01.DatabaseCategory.pdf", width = 7, height = 5)
showDatabaseCategory(CellChatDB)
dev.off()

# Preprocess expression data for CellChat analysis
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)          # Identify overexpressed genes in each cell group
cellchat <- identifyOverExpressedInteractions(cellchat)   # Identify overexpressed ligand-receptor interactions
cellchat <- projectData(cellchat, PPI.human)              # Project data onto protein-protein interaction network (human)

# Compute communication probability between cell groups
cellchat <- computeCommunProb(cellchat)
# Filter out cell-cell communications involving fewer than 10 cells
cellchat <- filterCommunication(cellchat, min.cells = 10)
# Extract and save the inferred cell-cell communication network
df.net = subsetCommunication(cellchat)
write.table(file = "COMM02.Comm.network.xls", df.net, sep = "\t", row.names = FALSE, quote = FALSE)

# Infer cell-cell communication at signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

# Aggregate the computed communication network
cellchat <- aggregateNet(cellchat)
# Visualize the number of interactions in the cell communication network
pdf(file = "COMM03.cellNetworkCount.pdf", width = 7, height = 6)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
dev.off()
# Visualize interaction strength in the cell communication network
pdf(file = "COMM04.cellNetworkWeight.pdf", width = 7, height = 6)
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction strength")
dev.off()

# Visualize cell-cell communication networks for each cell type individually
pdf(file = "COMM05.singleCell.pdf", width = 8, height = 6)
weight_mat <- cellchat@net$weight
par(mfrow = c(2,3), mgp = c(0,0,0), xpd = TRUE)
for (cel in unique(cellchat@idents)) {
  cir_mat <- matrix(0, nrow = nrow(weight_mat), ncol = ncol(weight_mat), dimnames = dimnames(weight_mat))
  cir_mat[cel, ] <- weight_mat[cel, ]
  netVisual_circle(cir_mat, vertex.weight = groupSize, weight.scale = TRUE, edge.weight.max = max(weight_mat), vertex.label.cex = 0.8, title.name = cel)
}
dev.off()

# Draw bubble plot showing ligand-receptor pairs and their interactions
pdf(file = paste0("COMM06.bubble.pdf"), width = 8, height = 5.5)
netVisual_bubble(cellchat, remove.isolate = FALSE, angle.x = 45)
dev.off()

# Read gene list file to identify core genes involved in cell communication and their pathways
geneRT = read.csv(geneFile, header = TRUE, sep = ",", check.names = FALSE)
hubGenes = as.vector(geneRT[, 1])
def.hub = df.net[((df.net$ligand %in% hubGenes) | (df.net$receptor %in% hubGenes)), ]
write.table(file = "COMM07.Comm.hubNetwork.xls", def.hub, sep = "\t", row.names = FALSE, quote = FALSE)

# Visualization of signaling pathway level communication
cellchat@netP$pathways     # Show all signaling pathways involved
pathways.show = "SPP1"     # Select the pathway to visualize (can be modified)
# Plot circle diagram of cell communication for selected signaling pathway
pdf(file = paste0("COMM08.", pathways.show , ".circle.pdf"), width = 8, height = 6)
circle = netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle", vertex.size = groupSize)
print(circle)
dev.off()

# Plot hierarchy diagram for the selected signaling pathway
pdf(file = paste0("COMM09.", pathways.show , ".hierarchy.pdf"), width = 12, height = 6)
hierarchy = netVisual_aggregate(cellchat, signaling = pathways.show, layout = "hierarchy", vertex.receiver = seq(1,4), vertex.size = groupSize)
print(hierarchy)
dev.off()

# Plot heatmap showing communication strength among cell groups for selected pathway
pdf(file = paste0("COMM10.", pathways.show , ".heatmap.pdf"), width = 8, height = 6)
heatmap = netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds", measure = 'weight')    
print(heatmap)
dev.off()

# Analyze signaling roles of cells in the selected pathway
pdf(file = paste0("COMM11.", pathways.show , ".netAnalysis.pdf"), width = 6, height = 5)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") 
netAnalysis = netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, width = 8, height = 5, font.size = 12)
print(netAnalysis)
dev.off()

# Identify which ligand-receptor pairs contribute to the selected pathway
pdf(file = paste0("COMM12.", pathways.show , ".contribution.pdf"), width = 8, height = 6)
contribution = netAnalysis_contribution(cellchat, signaling = pathways.show)
print(contribution)
dev.off()

# Plot expression of genes involved in the selected signaling pathway
pdf(file = paste0("COMM13.", pathways.show , ".geneExp.pdf"), width = 8, height = 6)
geneExp = plotGeneExpression(cellchat, signaling = pathways.show)
print(geneExp)
dev.off()

# Visualize individual ligand-receptor pair communication for the selected pathway (circle layout)
pairLR <- extractEnrichedLR(cellchat, signaling = pathways.show, geneLR.return = FALSE)
pdf(file = paste0("COMM14.", pathways.show , ".pairLR.pdf"), width = 9, height = 8)
pairCircos = netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = pairLR[1], layout = "circle")
print(pairCircos)
dev.off()

# Loop through all ligand-receptor pairs to visualize their communication in chord diagrams
for (i in 1:nrow(pairLR)) {
  pdf(file = paste0("COMM15.", pairLR[i,], ".pairLR.pdf"), width = 8, height = 6)
  pairChord = netVisual_individual(cellchat, signaling = pathways.show, pairLR.use = pairLR[i,], layout = "chord")
  print(pairChord)
  dev.off()
}
