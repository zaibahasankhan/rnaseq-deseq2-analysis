############################################################
# DESeq2 RNA-seq differential expression analysis
# Author: Dr. Zaiba Hasan Khan
############################################################

# This script includes:
# - Loading count data
# - DESeq2 end-to-end workflow
# - Visualization (volcano plot)
# - Saving results
# - Reproducibility information

############################################################
# STEP 1: Package installation (run once if not installed)
############################################################

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("DESeq2", "EnhancedVolcano", "ggplot2"))

############################################################
# STEP 2: Load libraries
############################################################

library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
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
pdf("MA_plot_DESeq2.pdf", width = 6, height = 5)

plotMA(
  res,
  ylim = c(-5, 5),
  alpha = 0.05,
  main = "MA plot"
)

dev.off()

############################################################
# STEP 14: Reproducibility
############################################################

sessionInfo()
