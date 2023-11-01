#Pathway analysis for Guzman et al. (In Review)

library(clusterProfiler)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(enrichplot)
library(ComplexHeatmap)
library(tidyr)


#read muscle schwann cell object
sc_obj <- readRDS("../sc_obj.rds")

###pathway for tscB compared to all other SC subtypes
terminal_diff.markers <- FindMarkers(sc_obj, ident.1 = "tSC_B", ident.2 = c('tSC_A','m-Schwann Cells_A', 'm-Schwann Cells_B'), min.pct = 0.25)
openxlsx::write.xlsx(terminal_diff.markers, "tscb_diff_markers.xlsx", rowNames=TRUE)


# Subset the data frame for rows with positive "avg_log2FC"
positive_diff_genes <- terminal_diff.markers[terminal_diff.markers$avg_log2FC > 1, ]
# Subset the data frame for rows with negative "avg_log2FC"
negative_diff_genes <- terminal_diff.markers[terminal_diff.markers$avg_log2FC < -1, ]

# Get the row names of the subsetted data frame, which are the gene names
genes_pos <- rownames(positive_diff_genes)
genes_neg <- rownames(negative_diff_genes)

# Run bitr function
gene_up <- bitr(genes_pos, fromType="SYMBOL", toType="ENTREZID", OrgDb= "org.Mm.eg.db")
gene_down <- bitr(genes_neg, fromType="SYMBOL", toType="ENTREZID", OrgDb= "org.Mm.eg.db")

enrich_up <- enrichGO(gene         = gene_up$SYMBOL,
                 OrgDb         = org.Mm.eg.db,
                 keyType       = 'SYMBOL',
                 ont           = "BP",
                 pAdjustMethod = "BH",
                 pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05)

enrich_down <- enrichGO(gene         = gene_down$SYMBOL,
                   OrgDb         = org.Mm.eg.db,
                   keyType       = 'SYMBOL',
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.01,
                   qvalueCutoff  = 0.05)


# Fig 4g
pdf("tscb_enrich_up_bp.pdf", width = 5.5, height = 3)
ggplot(enrich_up, showCategory=10, aes(GeneRatio, Description)) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(
    colours = c('#fac484', '#eb7f86', '#ce6693', '#5c53a5'),
    trans = "log10",
    guide = guide_colorbar(reverse = TRUE, order = 1)
  ) +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  xlab("Gene Ratio") +
  labs(title = "Biological Processes")
dev.off()

#Fig 4h
pdf("tscb_enric_dowh_bp.pdf", width = 5.5, height = 3)
ggplot(enrich_down, showCategory=10, aes(GeneRatio, Description)) +
  geom_segment(aes(xend = 0, yend = Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_gradientn(
    colours = c('#fac484', '#eb7f86', '#ce6693', '#5c53a5'),
    trans = "log10",
    guide = guide_colorbar(reverse = TRUE, order = 1)
  ) +
  scale_size_continuous(range = c(2, 10)) +
  theme_minimal() +
  xlab("Gene Ratio") +
  labs(title = "Biological Processes")
dev.off()


### Extended Data Fig 3

#extract the Fold Change
avg_log2FC <- data.frame(
  Gene = rownames(positive_diff_genes),
  avg_log2FC = positive_diff_genes$avg_log2FC)

# Initialize the df dataframe
df <- data.frame(Term = character(), 
                 Gene = character(), 
                 FoldChange = numeric())

# Define the top_terms variable (if not already defined)
top_terms <- head(enrich_up$Description, 10)  # Replace this line with actual top_terms if they are defined elsewhere


# Loop over top_terms to extract the associated genes and their fold changes
for (term in top_terms) {
  # Extract the genes associated with the current term
  genes <- unlist(strsplit(as.character(enrich_up$geneID[enrich_up$Description == term]), "/"))
  
  # Check if any extracted genes are present in avg_log2FC
  present_genes <- avg_log2FC$Gene[avg_log2FC$Gene %in% genes]
  
  # If no genes are present, continue to the next term
  if(length(present_genes) == 0) next
  
  # Extract the fold changes for the present genes
  fold_changes <- avg_log2FC$avg_log2FC[avg_log2FC$Gene %in% present_genes]
  
  # Create a data frame with the current term, present genes, and their fold changes
  temp_df <- data.frame(Term = rep(term, length(present_genes)), 
                        Gene = present_genes, 
                        FoldChange = fold_changes)
  
  # Append to the main data frame
  df <- rbind(df, temp_df)
}

print(df)

# Set the factor levels of Term in df to be in the same order as top_terms
df$Term <- factor(df$Term, levels = top_terms)

# Now, when you use tapply, the order of top_terms will be preserved in the matrix
mat <- with(df, tapply(FoldChange, list(Term, Gene), mean, default = NA))

# Reverse order of the pathway labels
mat <- mat[nrow(mat):1,]
mat <- mat[nrow(mat):1,]
print(mat)

# Identify the range for breaks
min_val <- min(df$FoldChange, na.rm = TRUE)
max_val <- max(df$FoldChange, na.rm = TRUE)



colors <- colorRampPalette(c('#fac484','#f8a07e','#eb7f86','#ce6693','#a059a0','#5c53a5'))(99)

# Make the Heatmap
pdf("heatmap_geneset.pdf", width = 10, height = 3)
Heatmap(mat, name = "Fold Change",
        col = colors,
        na_col = "white", 
        cluster_rows = FALSE, 
        cluster_columns = FALSE, 
        row_names_side = "left") 
dev.off()



