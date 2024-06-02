library(DESeq2)
library(pheatmap)
library(GenomicFeatures)
library(rtracklayer)
library(dplyr)
library(RColorBrewer)
library(biomaRt)
library(dendextend)
library(fgsea)
library(ggrepel)
library(clusterProfiler)
library(enrichplot)
library(EnhancedVolcano)


# Read Data
countData <- read.table("C:/Users/HP/Desktop/Dhaval/DDP/Results/merged_counts_VC.txt", header=TRUE, row.names=1)

# countData1 <- countData[, !(names(countData) %in% c("DMP11",'CFP08'))]

# Create colData with the 'condition' variable
colData <- data.frame(condition = c("healthy","healthy","healthy","healthy","diseased","diseased","diseased","diseased"))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData, colData = colData, design = ~ condition)

ntd <- normTransform(dds)

######################################################################################################################################################################################
vsd <- vst(dds, blind = FALSE)
head(assay(vsd), 3)
sampleDists <- dist(t(assay(vsd)))
sampleDists
sampleDistMatrix <- as.matrix( sampleDists )
rownames(sampleDistMatrix) <-  rownames(t(assay(vsd)))
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix, clustering_distance_rows = sampleDists, 
         clustering_distance_cols = sampleDists, col = colors)


rldl <- rlog(dds, blind = TRUE)

# Plot PCA using ggplot2
pca_data <- plotPCA(rldl, intgroup = c("condition"), returnData = TRUE)

# Extract percentage variance explained by each PC
variance_explained <- attr(pca_data, "percentVar")

# Rename PC1 and PC2 with the percentage variance
PC1_label <- paste0("PC1 (", round(variance_explained[1]*100, 2), "%)")
PC2_label <- paste0("PC2 (", round(variance_explained[2]*100, 2), "%)")

# Plotting
ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = rownames(pca_data))) +
  xlab(PC1_label) +
  ylab(PC2_label)


####################################################################################################
ntd <- normTransform(dds)
dds <- estimateSizeFactors(dds)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:40]
df <- as.data.frame(colData(dds)[,c("condition","sizeFactor")])
pheatmap(assay(ntd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
         cluster_cols=TRUE, annotation_col=df)

#pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#         cluster_cols=TRUE, annotation_col=df)
######################################################################################################################################################################################


# Run DESeq2 analysis
dds <- DESeq(dds)

contrast_condition <- c("condition", "diseased", "healthy")

# Get differential expression results for the specified contrast
results <- results(dds, contrast = contrast_condition)


# Get differential expression results
results <- results(dds)



# Filter based on both log2FC and Padj
filtered_genes <- subset(results, abs(log2FoldChange) > 1.5)

# Filter based on both log2FC and Padj
filtered_genes_padj <- subset(filtered_genes, padj > 0.7)


# Arrange in Descending order
filtered_genes <- filtered_genes[order(filtered_genes$log2FoldChange, decreasing = FALSE), ]

# Display top DE genes
#topGenes <- head(filtered_genes, n = 25)
topGenes <- filtered_genes

topGenes$gene_id <- rownames(topGenes)


# Extract normalized counts for DE genes
counts_DE <- counts(dds, normalized = TRUE)[rownames(topGenes), ]

# Define the column order for the heatmap
column_order <- c("CMVC14","CFVC15","CFVC16","CMVC17","DMVC18","DMVC19","DFVC20","DFVC21")

# Reorder columns in counts_DE based on column_order
counts_DE <- counts_DE[, column_order]

# Convert the matrix to a data frame
counts_DE <- as.data.frame(counts_DE)

# Add a new column named gene_id with values from the row names (index)
counts_DE$gene_id <- rownames(counts_DE)

gene_id <- gsub("\\.\\d+$", "", counts_DE$gene_id)




gtf_file <- "C:/Users/HP/Desktop/Dhaval/DDP/Bioinformatics/GTF file/gencode.v38.annotation.gtf"
gtf <- readGFF(gtf_file)

# Extract unique combinations of gene_id and gene_name
unique_gene_info <- unique(gtf[, c("gene_id", "gene_name")])
filtered_geness <- filtered_genes
filtered_geness$gene_id <- rownames(counts_DE)
deseq_results_df <- as.data.frame(filtered_geness)
gene_log2F_VC <- merge(deseq_results_df, unique_gene_info, by = "gene_id", all.x = TRUE)



merged_counts_DE <- merge(counts_DE, unique_gene_info, by = "gene_id", all.x = TRUE)
unique_gene_names <- make.unique(as.character(merged_counts_DE$gene_name))
rownames(merged_counts_DE) <- unique_gene_names
merged_counts_DE <- merged_counts_DE[, -which(colnames(merged_counts_DE) == "gene_name")]
merged_counts_DE <- merged_counts_DE[, -which(colnames(merged_counts_DE) == "gene_id")]



# Create a heatmap with specified column order
pheatmap(merged_counts_DE, cluster_rows = TRUE, cluster_cols = TRUE, show_rownames = F,treeheight_row = 0,main = "Heatmap of Differentially Expressed Genes", plot = FALSE)


z_scores <- scale(merged_counts_DE, center = TRUE, scale = TRUE)

# Define custom order for columns
custom_column_order <- c("CMVC14","CFVC15","CFVC16","CMVC17","DMVC18","DMVC19","DFVC20","DFVC21")  # Add your sample names in the desired order

# Plot heatmap using pheatmap with custom column order
pheatmap(z_scores, 
              cluster_rows = TRUE,                  # Cluster rows (genes)
              cluster_cols = TRUE,   # Do not cluster columns
              clustering_method = "complete",
              scale = "row",                        # Scale rows (genes)
              main = "Differentially Expressed Genes Heatmap",  # Main title
              color = colorRampPalette(c("blue", "white", "red"))(100),  # Color palette
              show_rownames = FALSE,               # Remove row annotations
              treeheight_row = 0,                  # Remove row tree
              labels_col = custom_column_order     # Specify custom column order
)


#############################################

# Volcano Plot
# Assuming 'results' is your data frame
plot(results$log2FoldChange, -log10(results$pvalue), 
     main="Differentially Expressed Genes for Healthy vs Alzheimer's", 
     xlab="log2(Fold Change)", 
     ylab="-log10(P-value)", 
     pch=20, 
     col=ifelse(results$log2FoldChange < -1.5, "blue",
                ifelse(results$log2FoldChange > 1.5, "red", "black")))

# Add labels for upregulated and downregulated genes at the top
text(0, max(-log10(results$pvalue)) + 1, labels="Downregulated", col="blue", pos=3)
text(0, max(-log10(results$pvalue)) + 1, labels="Upregulated", col="red", pos=2)

# Draw a vertical line at log2FoldChange = 0
abline(v=0, col="gray", lty=2)


# Volcano Plot
png("volcano_plot.png", width = 10, height = 4, units = "in", res = 300) 

# Volcano Plot
# Assuming 'results' is your data frame

# Define the x-axis limits
x_limits <- c(-4, 4)

# Plot the volcano plot with specified x-axis limits
plot(results$log2FoldChange, -log10(results$pvalue), 
     main="Differentially Expressed Genes for Healthy vs Alzheimer's", 
     xlab="log2(Fold Change)", 
     ylab="-log10(P-value)", 
     pch=20, 
     col=ifelse(results$log2FoldChange < -1.5 , "blue",
                ifelse(results$log2FoldChange > 1.5, "red", "black")),
     xlim=x_limits  # Set x-axis limits
)

# Add labels for upregulated and downregulated genes at the top
text(0, max(-log10(results$pvalue)) + 1, labels="Downregulated", col="blue", pos=3)
text(0, max(-log10(results$pvalue)) + 1, labels="Upregulated", col="red", pos=2)

# Draw a vertical line at log2FoldChange = 0
abline(v=0, col="gray", lty=2)

##############################################


# Extract normalized counts for DE genes
#results$gene_id <- rownames(results)

counts_DE_all <- counts(dds, normalized = TRUE)[rownames(results), ]
# Define the column order for the heatmap
column_order <- c("CMVC14","CFVC15","CFVC16","CMVC17","DMVC18","DMVC19","DFVC20","DFVC21")

# Reorder columns in counts_DE based on column_order
counts_DE_all <- counts_DE_all[, column_order]

# Convert the matrix to a data frame
counts_DE_all <- as.data.frame(counts_DE_all)

# Add a new column named gene_id with values from the row names (index)
counts_DE_all$gene_id <- rownames(counts_DE_all)
counts_DE_all$gene_id <- rownames(counts_DE_all)
# Extract unique combinations of gene_id and gene_name
unique_gene_info <- unique(gtf[, c("gene_id", "gene_name","gene_type")])

merged_counts_DE_all <- merge(counts_DE_all, unique_gene_info, by = "gene_id", all.x = TRUE)
duplicates <- duplicated(merged_counts_DE_all$gene_name)
merged_counts_DE_all <- merged_counts_DE_all[!duplicates, ]
merged_counts_DE_all <- merged_counts_DE_all[!is.na(merged_counts_DE_all$gene_name), ]
rownames(merged_counts_DE_all) <- merged_counts_DE_all$gene_name
merged_counts_DE_all <- merged_counts_DE_all[, -which(colnames(merged_counts_DE_all) == "gene_name")]
merged_counts_DE_all <- merged_counts_DE_all[, -which(colnames(merged_counts_DE_all) == "gene_id")]

# Define the entries you want to keep
entries_to_keep <- c("protein_coding")  
merged_counts_DE_alll <- merged_counts_DE_all[merged_counts_DE_all$gene_type %in% entries_to_keep, ]
merged_counts_DE_alll$gene_type <- NULL
z_scores_all <- scale(merged_counts_DE_alll, center = FALSE, scale = TRUE)
# Remove rows with NaN values
# Remove rows with NA values
z_scores_all_complete <- z_scores_all[complete.cases(z_scores_all), ]

# Remove rows with NA values
z_scores_all_complete <- z_scores_all_complete[complete.cases(z_scores_all_complete), ]
z_scores_all_complete <- z_scores_all[complete.cases(t(z_scores_all)), ]

# Plot the heatmap
pheatmap(z_scores_all_complete, 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "complete",
         scale = "row",
         main = "Differentially Expressed Genes Heatmap",
         color = colorRampPalette(c("blue", "white", "red"))(100),
         show_rownames = FALSE,  # Remove row annotations
         treeheight_row = 0
)



############################################

# Extract gene IDs, gene names, and expression values from topGenes
gene_ids <- rownames(topGenes)
gene_names <- rownames(topGenes)
expression_values <- topGenes$log2FoldChange

# Create a dataframe with the extracted information
topGenes_df <- data.frame(gene_id = gene_ids, gene_name = gene_names, log2FoldChange = expression_values)

# Merge topGenes_df with unique_gene_info
merged_counts_genes_brain <- merge(topGenes_df, unique_gene_info, by = "gene_id", all.x = TRUE)


gene_name.y <- c("SIGLEC1", "RPH3A", "CRYM", "MET", "MS4A6A", "TRIP10", 
                 "PARVG", "FBLN7", "LRMDA", "PAK1", "SLC7A11", "TLR5", 
                 "SIRPB2", "NFAM1", "RP11-399B17.1")

# Create the dataframe
gene_names_df <- data.frame(gene_name.y = gene_name.y)
common_genes_brain <- merge(merged_counts_genes_brain, gene_names_df, by = "gene_name.y", all.x = FALSE)
##############################################

corMatrix <- cor(merged_counts_DE)
pheatmap(corMatrix)

# Volcano Plot
# Assuming 'results' is your data frame
plot(results$log2FoldChange, -log10(results$pvalue), 
     main="Differentially Expressed Genes for Healthy vs Alzheimer's", 
     xlab="log2(Fold Change)", 
     ylab="-log10(P-value)", 
     pch=20, 
     col=ifelse(results$padj < 0.05 & results$log2FoldChange < 5, "blue",
                ifelse(results$padj < 0.05 & results$log2FoldChange > 5, "red", "black")))

# Add labels for upregulated and downregulated genes at the top
text(0, max(-log10(results$pvalue)) + 1, labels="Downregulated", col="blue", pos=3)
text(0, max(-log10(results$pvalue)) + 1, labels="Upregulated", col="red", pos=2)

# Draw a vertical line at log2FoldChange = 0
abline(v=0, col="gray", lty=2)
#######################################################################################

rownames(results) <- c("gene_id", rownames(results)[-1])

# Convert the matrix to a data frame
results <- as.data.frame(results)

# Add a new column named gene_id with values from the row names (index)
results$gene_id <- rownames(results)

merged_counts_volcano <- merge(results, unique_gene_info, by = "gene_id", all.x = TRUE)
# Filter out rows with missing gene names
# Subset merged_counts_volcano to remove rows with NA, NaN, or Inf in log2FoldChange and pvalue
merged_counts_volcano_clean <- merged_counts_volcano[complete.cases(merged_counts_volcano[, c('log2FoldChange', 'pvalue')]), ]

# Subset lab to match the rows in merged_counts_volcano_clean
lab_subset <- merged_counts_volcano$gene_name[complete.cases(merged_counts_volcano[, c('log2FoldChange', 'pvalue')])]

# Assuming merged_counts_volcano_clean is your cleaned data frame
# Subset top 5 genes with highest log2FoldChange
top_genes <- merged_counts_volcano_clean[order(merged_counts_volcano_clean$log2FoldChange, decreasing = TRUE), ][1:5, ]

# Subset bottom 5 genes with lowest log2FoldChange
bottom_genes <- merged_counts_volcano_clean[order(merged_counts_volcano_clean$log2FoldChange, decreasing = FALSE), ][1:5, ]


# Plot EnhancedVolcano with the corrected lab argument
# Plot EnhancedVolcano with specific labels for top 5 and bottom 5 genes
EnhancedVolcano(merged_counts_volcano_clean, 
                lab = merged_counts_volcano_clean$gene_name,  # Use gene names as labels
                x = 'log2FoldChange', 
                y = 'pvalue', 
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 3.0,
                legendLabSize = 12,
                legendIconSize = 4.0,
                drawConnectors = TRUE,
                widthConnectors = 0.25,
                title = "Healthy vs Alzheimer's",
                subtitle = "",
                selectLab = c(top_genes$gene_name, bottom_genes$gene_name)  # Only label top and bottom genes
)

#########################################################################################
  
# Gene Set Enrichment Analysis
  
  #if (!require("BiocManager", quietly = TRUE))
  #  install.packages("BiocManager")
  
  #BiocManager::install("clusterProfiler", force = TRUE)
  #BiocManager::install("pathview", force = TRUE)
  #BiocManager::install("enrichplot")
  library(clusterProfiler)
  library(enrichplot)
  # we use ggplot2 to add x axis labels (ex: ridgeplot)
  library(ggplot2)
  
  
  organism <- "org.Hs.eg.db"
  # Install the annotation package for human genome
  if (!requireNamespace(organism, quietly = TRUE)) {
    BiocManager::install(organism, character.only = TRUE)
  }
  
  # Load the annotation package for human genome
  library(org.Hs.eg.db)
  
  # we want the log2 fold change 
  original_gene_list <- topGenes$log2FoldChange
  
  gene_id <- rownames(topGenes)
  
  gene_list_trimmed <- gsub("\\.\\d+$", "", gene_id)
  # name the vector
  names(original_gene_list) <- gene_list_trimmed
  
  # omit any NA values 
  gene_list<-na.omit(original_gene_list)
  
  # sort the list in decreasing order (required for clusterProfiler)
  gene_list = sort(gene_list, decreasing = TRUE)
  
  
  gse <- gseGO(geneList=gene_list, 
               ont ="ALL", 
               keyType = "ENSEMBL", 
               nPerm = 10000, 
               minGSSize = 3, 
               maxGSSize = 800, 
               pvalueCutoff = 0.05, 
               verbose = TRUE, 
               OrgDb = organism, 
               pAdjustMethod = "none")
  
  require(DOSE)
  dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
  
  gse1 <- pairwise_termsim(gse)
  emapplot(gse1, showCategory = 10)
  
  cnetplot(gse, categorySize="pvalue", foldChange=gene_list, showCategory = 5)
  
  ridgeplot(gse) + labs(x = "enrichment distribution")
  
  ########################################################################
  topGenes$gene_id <- rownames(topGenes)
  topGenes <- as.data.frame(topGenes)
  topGenes_df_brain <- merge(topGenes, unique_gene_info, by = "gene_id", all.x = TRUE)
  filtered_df <- topGenes_df_brain[topGenes_df_brain$gene_type == "protein_coding", ]
  write.csv(filtered_df, "topGenes_df_brain_VC.csv", row.names = TRUE)
  
  
  results$gene_id <- rownames(results)
  results <- as.data.frame(results)
  results_df_brain <- merge(results, unique_gene_info, by = "gene_id", all.x = TRUE)
  
  
  # Your original code (assuming 'gse' is your data frame)
  your_plot <- dotplot(gse, showCategory = 8, split = ".sign") +
    facet_grid(.~.sign)
  
  # Modify the y-axis text size
  your_plot + theme(axis.text.y = element_text(size = 7))
  
  #####################################################################################

