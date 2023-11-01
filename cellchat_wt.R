##cellchat for WT

library(NMF)
library (circlize)
library (devtools)
library (BiocGenerics)
library (CellChat)
library (Seurat)
library(Seurat)
library(rhdf5)
library(dplyr)
library(viridis)
library(writexl)
library(ggplot2)
library(metap)
library(cowplot)
library(pheatmap)
library(hdf5r)
library(patchwork)
library(ggalluvial)
library(igraph)
library(ComplexHeatmap)
options(stringsAsFactors = FALSE)

setwd("/Users/steed/Dropbox (University of Michigan)/bi_working/guzman")
getwd()
future::plan("multisession", workers = 4) # do parallel

proj <- readRDS("../proj.rds")
wt <-subset(proj, subset = orig.ident == c("WT_GFP+", "WT_GFP-"), 
                idents = c("m-Schwann Cells_A", "m-Schwann Cells_B",
                           "tSC_A", "tSC_B",
                           "Tenm2-High", "Mesenchymal Stem Cells",
                           "Satellite Cells", "Immune Cells", 'Proliferating'))

DimPlot(wt, reduction = 'umap', split.by = 'orig.ident')

data.input <- GetAssayData(wt, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(wt)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels


cellchat_wt <- createCellChat(object = data.input, meta = meta, group.by = "group")

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat_wt@DB <- CellChatDB.use


# subset the expression data of signaling genes for saving computation cost
cellchat_wt <- subsetData(cellchat_wt) # This step is necessary even if using the whole database


cellchat_wt <- identifyOverExpressedGenes(cellchat_wt)
cellchat_wt <- identifyOverExpressedInteractions(cellchat_wt)
cellchat_wt <- computeCommunProb(cellchat_wt)

saveRDS(cellchat_wt, file = "../cellchat.rds")
cellchat_wt <- readRDS("../cellchat.rds")


cellchat_wt <- filterCommunication(cellchat_wt, min.cells = 10)
cellchat_wt <- computeCommunProbPathway(cellchat_wt)
cellchat_wt <- aggregateNet(cellchat_wt)

groupSize <- as.numeric(table(cellchat_wt@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_wt@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_wt@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


#sod1 dataset
mat <- cellchat_wt@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


#Heiarchy plot
vertex.receiver = seq(1,4) # a numeric vector. 
netVisual_aggregate(cellchat_wt, signaling = pathways.show,  vertex.receiver = vertex.receiver)

sig_pathways <- cellchat_wt@netP$pathways
# Compute the network centrality scores
cellchat_wt <- netAnalysis_computeCentrality(cellchat_wt, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways

