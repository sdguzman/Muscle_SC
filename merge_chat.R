##Cell Chat analysis of WT and Sod1, combined.
# Please run cellchat_wt and cellchat_sod1 before running this code.

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


cellchat_wt <- readRDS("../cellchat.rds")
cellchat_sod1 <- readRDS("../cellchat_sod1.rds")
#saveRDS(cellchat, file = "../guzman/cellchat.rds")
cellchat <- readRDS("../guzman/cellchat.rds")


#merging objects
object.list <- list(WT = cellchat_wt, Sod1KO = cellchat_sod1)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))


#Differential number of interactions or interaction strength among different cell populations
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#control of nodes
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Number of interactions - ", names(object.list)[i]))
}

#Compare the major sources and targets in 2D space
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax)
}


#Figure 5a
#Compare the overall information flow of each signaling pathway
pdf("information_flow.pdf", height = 4, width = 6)
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE, color.use =c("#1E4485","#138CA3"))
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE, color.use =c("#1E4485","#138CA3"))
gg1 + gg2
dev.off()

#Compare outgoing (or incoming) signaling associated with each cell population
i = 1
cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WT", "Sod1KO")) # set factor level

#Extended Data Figures 4a and 4b
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
top20 <- pathway_subset <- pathway.union[1:20]

pdf('incoming_secreted.pdf', height = 6, width = 8)
i = 1
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "incoming", signaling = top20, title = names(object.list)[i], width = 5, height = 6, color.heatmap = "GnBu")
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = top20, title = names(object.list)[i+1], width = 5, height = 6, color.heatmap = "GnBu")
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()

pdf('outgoing_secreted.pdf', height = 6, width = 8)
i = 1
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = top20, title = names(object.list)[i], width = 5, height = 6)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = top20, title = names(object.list)[i+1], width = 5, height = 6)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
dev.off()






