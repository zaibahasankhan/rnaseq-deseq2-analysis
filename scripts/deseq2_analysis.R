############################################################
# DESeq2 RNA-seq differential expression analysis
# Author: Dr. Zaiba Hasan Khan
############################################################

# This script includes:
# - Loading count data and metadata
# - DESeq2 end-to-end workflow
# - Visualization (MA plot, volcano plot, PCA, heatmap)
# - Saving results
# - Reproducibility information

############################################################
# STEP 1: Required packages
############################################################

# NOTE:
# Required packages should be installed ONCE before running this script.
# This script assumes the following packages are already installed:
# - BiocManager
# - DESeq2
# - EnhancedVolcano
# - ggplot2
# - pheatmap

############################################################
# STEP 2: Load libraries
############################################################

library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(pheatmap)
############################################################
# STEP 3: Load data
############################################################

Count_data <- read.csv(
  file = "data/counts.csv",
  header = TRUE,
  row.names = 1,
  check.names = FALSE
)

Col_data <- read.csv(
  file = "data/column.csv",
  header = TRUE,
  row.names = 1
)

############################################################
# STEP 4: Data checks
############################################################

dim(Count_data)
dim(Col_data)

all(rownames(Col_data) == colnames(Count_data))

############################################################
# STEP 5: NA handling
############################################################

sum(is.na(Count_data))
which(is.na(Count_data), arr.ind = TRUE)

# Optional NA handling
# Count_data <- na.omit(Count_data)
# Count_data[is.na(Count_data)] <- 0

############################################################
# STEP 6: Filter genes with all zero counts (optional)
############################################################

Count_data <- Count_data[rowSums(Count_data) > 0, ]

############################################################
# STEP 7: Create DESeq2 object
############################################################

dds <- DESeqDataSetFromMatrix(
  countData = Count_data,
  colData = Col_data,
  design = ~ condition
)

############################################################
# STEP 8: Set factor level
############################################################

dds$condition <- relevel(dds$condition, ref = "Control")

############################################################
# STEP 9: Low-count filtering
############################################################

filtered_count <- rowSums(counts(dds)) >= 10
dds <- dds[filtered_count, ]

############################################################
# STEP 10: Run DESeq2
############################################################

dds <- DESeq(dds)
res <- results(dds)
summary(res)

############################################################
# STEP 11: Extract significant genes
############################################################

resSigUp <- subset(res, padj < 0.05 & log2FoldChange > 1)
resSigDown <- subset(res, padj < 0.05 & log2FoldChange < -1)

resSig_combined <- subset(
  res,
  padj < 0.05 & abs(log2FoldChange) > 1
)

############################################################
# STEP 12: Save outputs
############################################################

write.csv(resSigUp, "results/Upregulated_genes.csv")
write.csv(resSigDown, "results/Downregulated_genes.csv")
write.csv(resSig_combined, "results/All_DEGs.csv")

############################################################
# STEP 13: Volcano plot
############################################################

pdf("results/volcano_plot_DESeq2.pdf", width = 8, height = 8)


EnhancedVolcano(res,
                lab = "",
                x = 'log2FoldChange',
                y = 'padj',
                selectLab = "NULL",
                title = 'Volcano plot',
                cutoffLineType = 'twodash',
                cutoffLineWidth = 0.8,
                xlim = c(-8,8),
                xlab = bquote(~Log[2]~ 'fold change'),
                ylim = c(0,12),
                ylab = bquote(~-Log[10]~italic(P)),
                pCutoff = 0.05,
                FCcutoff = 1.0,
                #transcriptPointSize = 0.5,
                #transcriptLabSize = 4.0,
                colAlpha = 1,
                shape = 19,
                subtitle = NULL,
                legendPosition = 'top',
                legendLabSize = 12,
                legendIconSize = 4.0,
                gridlines.major = FALSE,
                gridlines.minor = FALSE,
                drawConnectors = FALSE,
                widthConnectors = 0.2,
                colConnectors = 'grey50',
                border = 'full' )


dev.off()

############################################################
# STEP 14: MA PLOT
############################################################
pdf("results/MA_plot_DESeq2.pdf", width = 6, height = 5)

plotMA(
  res,
  ylim = c(-5, 5),
  alpha = 0.05,
  main = "MA plot"
)

dev.off()

############################################################
# STEP 15: Normalized data (median to ratio normalized counts)
############################################################
norm_counts <- counts(dds, normalized = TRUE)


############################################################
# STEP 16: PCA plot (sample clustering) on normalized count data
############################################################
#VST-transformation
vsd <- vst(dds, blind = FALSE)

pcaData <- plotPCA(vsd, intgroup = "condition", returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

p <- ggplot(pcaData, aes(PC1, PC2, color = condition)) +
  geom_point(size = 4) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal() +
  ggtitle("PCA plot")

pdf("results/PCA_plot.pdf", width = 6, height = 5)

print(p)

dev.off()

############################################################
# STEP 17: Heatmap
############################################################
#vsd <- vst(dds, blind = FALSE) #VST-transformation if not done in PCA step

# Select top 100 DEGs
topGenes <- head(order(res$padj), 100)

# Extract VST normalized expression
mat <- assay(vsd)[topGenes, ]

# Center genes (important for visualization)
mat <- mat - rowMeans(mat)

# Sample annotation
annotation_col <- as.data.frame(colData(dds)[, "condition", drop = FALSE])

pdf("results/Heatmap_top100_DEGs.pdf", width = 7, height = 8)

pheatmap(
  mat,
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 6,
  cluster_cols = TRUE,
  cluster_rows = TRUE,
  scale = "row"
)

dev.off()

############################################################
# STEP 18: Reproducibility
############################################################

sessionInfo()
