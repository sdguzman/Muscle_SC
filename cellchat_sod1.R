##cellchat analysis for Guzman et al (in review)

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

#setup working directory
setwd("/Users/steed/Dropbox (University of Michigan)/bi_working/guzman")
getwd()
options(future.globals.maxSize= 786432000)

#Read in RDS File
proj <- readRDS("../proj.rds")
sod1 <-subset(proj, subset = orig.ident == c("Sod1KO_GFP+", "Sod1KO_GFP-"), 
              idents = c("m-Schwann Cells_A", "m-Schwann Cells_B",
                         "tSC_A", "tSC_B",
                         "Tenm2-High", "Mesenchymal Stem Cells",
                         "Satellite Cells", "Immune Cells", "Proliferating"))

DimPlot(sod1, reduction = "umap", split.by = "orig.ident")


data.input <- GetAssayData(sod1, assay = "RNA", slot = "data") # normalized data matrix
labels <- Idents(sod1)
meta <- data.frame(group = labels, row.names = names(labels)) # create a dataframe of the cell labels

cellchat_sod1 <- createCellChat(object = data.input, meta = meta, group.by = "group")

#Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat_sod1@DB <- CellChatDB.use

# subset the expression data of signaling genes for saving computation cost
cellchat_sod1 <- subsetData(cellchat_sod1) # This step is necessary even if using the whole database


cellchat_sod1 <- identifyOverExpressedGenes(cellchat_sod1)
cellchat_sod1 <- identifyOverExpressedInteractions(cellchat_sod1)
cellchat_sod1 <- computeCommunProb(cellchat_sod1)

saveRDS(cellchat_sod1, file = "../cellchat_sod1.rds")
cellchat_sod1 <- readRDS("../cellchat_sod1.rds")


cellchat_sod1 <- filterCommunication(cellchat_sod1, min.cells = 10)
cellchat_sod1 <- computeCommunProbPathway(cellchat_sod1)
cellchat_sod1 <- aggregateNet(cellchat_sod1)

groupSize <- as.numeric(table(cellchat_sod1@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat_sod1@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat_sod1@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")


#sod1 dataset
mat <- cellchat_sod1@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}


