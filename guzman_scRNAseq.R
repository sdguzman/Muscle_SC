# Contains the code used for scGEX data analysis for Guzman et al., (Currently in Review)

# Workflow: 
# 1. Initialize variables
# 2. Load and merge data, QC, and Non-linear reduction analysis
# 3. Plots for Figure 4 

library(Seurat) # v4.3
library(celda) # v1.10.0 decontX containing package
library(SingleCellExperiment)
library(rhdf5)
library(dplyr)
library(ggplot2) # v3.4.1
library(hdf5r)
library(SeuratDisk)
library(tidyverse)
library(openxlsx)



setwd("/Users/steed/Dropbox (University of Michigan)/bi_working/guzman")
getwd()
options(future.globals.maxSize= 786432000)

# list of objects to read
files <- list(
  s1 = list(mtx = "WT_GFPpos_matrix.mtx.gz", features = "WT_GFPpos_features.tsv.gz", cells = "WT_GFPpos_barcodes.tsv.gz"),
  s2 = list(mtx = "WT_GFPneg_matrix.mtx.gz", features = "WT_GFPneg_features.tsv.gz", cells = "WT_GFPneg_barcodes.tsv.gz"),
  s3 = list(mtx = "Sod1KO_GFPpos_matrix.mtx.gz", features = "Sod1KO_GFPpos_features.tsv.gz", cells = "Sod1KO_GFPpos_barcodes.tsv.gz"),
  s4 = list(mtx = "Sod1KO_GFPneg_matrix.mtx.gz", features = "Sod1KO_GFPneg_features.tsv.gz", cells = "Sod1KO_GFPneg_barcodes.tsv.gz")
)

# Initialize a named list to hold the results
results <- list(s1 = NULL, s2 = NULL, s3 = NULL, s4 = NULL)

# apply same operations to each object in the list
for (x in names(files)) {
  
  # read mtx file
  s <- ReadMtx(mtx = files[[x]]$mtx, features = files[[x]]$features, cells = files[[x]]$cells)
  
  # create SingleCellExperiment object
  sce <- SingleCellExperiment(list(counts = s))
  
  # apply decontX function
  s <- decontX(sce)
  
  # compute UMAP
  umap <- reducedDim(s, "decontX_UMAP")
  
  # plot dimensions
  plotDimReduceCluster(x = s$decontX_clusters, dim1 = umap[, 1], dim2 = umap[, 2])
  
  # plot decontX contamination
  plotDecontXContamination(s)
  
  # Store the result in the results list
  results[[x]] <- s
}

# assign results back to global environment
list2env(results, envir = .GlobalEnv)


wt_GFPpos <- CreateSeuratObject(round(decontXcounts(s1)), project = "WT_GFP+", min.cells = 3, min.features = 200)
wt_GFPneg <- CreateSeuratObject(round(decontXcounts(s2)), project = "WT_GFP-", min.cells = 3, min.features = 200)
Sod1KO_GFPpos <- CreateSeuratObject(round(decontXcounts(s3)), project = "Sod1KO_GFP+", min.cells = 3, min.features = 200)
Sod1KO_GFPneg <- CreateSeuratObject(round(decontXcounts(s4)), project = "Sod1KO_GFP-", min.cells = 3, min.features = 200)


#QC1
wt_GFPpos[["percent.mt"]] <- PercentageFeatureSet(wt_GFPpos, pattern = "^mt-")
VlnPlot(wt_GFPpos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
wt_GFPpos <- subset(wt_GFPpos, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 55000 & percent.mt < 10)


#QC2
wt_GFPneg[["percent.mt"]] <- PercentageFeatureSet(wt_GFPneg, pattern = "^mt-")
VlnPlot(wt_GFPneg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
wt_GFPneg <- subset(wt_GFPneg, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 50000 & percent.mt < 10)


#QC3
Sod1KO_GFPpos[["percent.mt"]] <- PercentageFeatureSet(Sod1KO_GFPpos, pattern = "^mt-")
VlnPlot(Sod1KO_GFPpos, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Sod1KO_GFPpos <- subset(Sod1KO_GFPpos, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 60000 & percent.mt < 10)

#QC4
Sod1KO_GFPneg[["percent.mt"]] <- PercentageFeatureSet(Sod1KO_GFPneg, pattern = "^mt-")
VlnPlot(Sod1KO_GFPneg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
Sod1KO_GFPneg <- subset(Sod1KO_GFPneg, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & nCount_RNA < 65000 & percent.mt < 10)


#Merging was applied as all samples were collected on the same day. Library preparation and sequneincing was also performed from a single batch.
proj <- merge(wt_GFPpos, y = c(wt_GFPneg, Sod1KO_GFPpos, Sod1KO_GFPneg), add.cell.ids = c("WT_S100GFP+", "WT_S100GFP-", "SOD1KO_S100GFP+", "SOD1KO_S100GFP-"), project = "Sod1")
VlnPlot(proj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

proj<- NormalizeData(proj)
proj <- FindVariableFeatures(proj, selection.method = "vst", nfeatures = 2000)




all.genes <- rownames(proj)
proj <- ScaleData(proj, features = all.genes)
proj <- RunPCA(proj, features = VariableFeatures(object = proj))
proj <- JackStraw(proj, num.replicate = 100)
proj <- ScoreJackStraw(proj, dims = 1:20)
ElbowPlot(proj)
proj <- FindNeighbors(proj, dims = 1:15)
proj <- FindClusters(proj, resolution = 0.5)
proj <- RunUMAP(proj, dims = 1:15)
DimPlot(proj, reduction = "umap")

# ########SAVING AND READING FILES
saveRDS(proj, file = "../proj.rds")
proj <- readRDS("../proj.rds")


# find markers for every cluster compared to all remaining cells, report only the positive. 
proj.markers <- FindAllMarkers(proj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#save results into an excel sheet and determine cell clusters.
openxlsx::write.xlsx(proj.markers, "markers.xlsx", rowNames=TRUE)

##Renaming idents/clusters

cluster.ids <- c("Tenocytes",
                 "Endothelial",
                  "m-Schwann Cells_A",
                 "Smooth Muscle",
                  "Mesenchymal Stem Cells",
                 "FAP_A",
                 "Tenm2-High",
                  "Smooth Muscle",
                 "tSC_B",
                 "Smooth Muscle",
                 "Endothelial",
                 "tSC_A",
                 "m-Schwann Cells_B",
                 "Immune Cells",
                 "Smooth Muscle",
                 "FAP_B",
                 "Mesenchymal Stem Cells",
                 "RBC",
                 "Satellite Cells",
                 "Proliferating",
                 "Endothelial",
                 "Rbc",
                 "Rbc",
                  "Fibroblasts")

names(cluster.ids) <- levels(proj)
proj <- RenameIdents(proj, cluster.ids)

#changing levels for plotting
levels(proj) <- c("m-Schwann Cells_A",  "m-Schwann Cells_B","tSC_A", "tSC_B", "Tenm2-High", "Mesenchymal Stem Cells",  "FAP_A",  "FAP_B", "Proliferating", "Fibroblasts",
                  "Satellite Cells", "Immune Cells",  "Endothelial", "Smooth Muscle", "Tenocytes")

###bioinformatically remove RBCs
proj <- subset(x = proj, idents = 'Rbc', invert = TRUE)


#Change levels of groups and add "group" meta-data
proj$orig.ident <- factor(proj$orig.ident, levels = c("WT_GFP+","WT_GFP-" ,'Sod1KO_GFP+', 'Sod1KO_GFP-'))

#add metadata
proj$group <- ifelse(proj$orig.ident %in% c("WT_GFP+", "WT_GFP-"), "WT", 
                     ifelse(proj$orig.ident %in% c("Sod1KO_GFP+", "Sod1KO_GFP-"), "KO", NA))
proj$group <- factor(proj$group, levels = c('WT','KO'))


#subset schwann cell clusters
sc_obj = subset(x = proj, idents = c("m-Schwann Cells_A", "m-Schwann Cells_B", "tSC_A", "tSC_B"))

#number ofmuscle residentSCs
cell.no <- as.data.frame.matrix(table(sc_obj@active.ident, sc_obj@meta.data$orig.ident))
cell.no
write_xlsx(cell.no, "/Users/steed/Dropbox (University of Michigan)/bi_working/guzman/sc_numbers.xlsx")




#Figure 4B
s100b_data <- data.frame(
  OrigIdent = proj@meta.data$orig.ident,
  Expression = FetchData(proj, vars = "S100b")$S100b
)

# Filtering NAs if any
s100b_data <- s100b_data %>% drop_na()

# Define the colors
colors <- c('#30B276', '#CCCCCC', '#30B276', '#CCCCCC')

pdf("s100b.pdf", height = 3, width = 3.5)
ggplot(s100b_data, aes(x = OrigIdent, y = Expression, fill = OrigIdent)) +
  geom_jitter(
    color = "grey", 
    position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.9), 
    size = 1, alpha = 0.5, 
    stroke = 0
  ) +
  geom_boxplot(alpha = 0.5, outlier.shape = NA) +
  geom_hline(yintercept = 0, color = "black", size = 0.5) +
  scale_fill_manual(values = colors) +
  labs(title = "Expression of S100b", x = "Original Identity", y = "Expression Level") +
  theme_minimal() + 
  theme(
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(), 
    axis.line.y = element_line(size = 0.5, colour = "black"),
    axis.ticks.y = element_line(color = "black", size = 0.25) 
  )
dev.off() 


#Figure 4c
# Define the colors
colors <- c("#D078FF", "#F8766D", "#00C0B5", "#619CFF", "#53B400", "#DA8F00", "#EB8335", "#00ABFD", 
            "#00BDD2", "#86AC00", "#A58AFF", "#FF63B9", "#FF6B96", "#00B6EB", "#FB61D7", "#00BE6D", 
            "#00C094", "#C49A00", "#00BA38", "#A9A400", "#EC69EF")
pdf("dimplot.pdf", width = 3.5, height = 3.5, units = 'in', res = 300)
DimPlot(proj, reduction = 'umap', cols = colors, label = FALSE) + NoLegend() & NoAxes()
dev.off()


#Figure 4d
pdf("sc_markers_dimplot.pdf", width = 8, height = 6)
FeaturePlot(proj, features = c("S100b", "Sox10", "Agrn", 'Mbp', "Mpz", "Ngfr", "Cspg4", "Kcnj10", "Ccnd1"), ncol = 3) & NoAxes()
dev.off()


#Figure 4e
pdf("dimplot_samples.pdf", width = 15, height = 3, units = 'in', res = 600) 
DimPlot(sc_obj, reduction = 'umap', split.by = "orig.ident", cols = c('#003F5C', '#58508D', '#BC5090', '#FF6361', '#FFA600')) + NoLegend() & NoAxes()
dev.off()

#Figure 4f
pdf("tsc_genes.pdf", width = 4, height = 7, units = 'in', res = 300) 
VlnPlot(sc_obj, 
        features = c("Cd44","Pde10a", 'Cadm1', "Runx2","Sox6", "Tnc"),
        cols = c('#003F5C', '#58508D', '#BC5090', '#FF6361', '#FFA600'),
        pt.size = 0, ncol = 2)  &  theme(axis.title.x = element_blank())
dev.off()




