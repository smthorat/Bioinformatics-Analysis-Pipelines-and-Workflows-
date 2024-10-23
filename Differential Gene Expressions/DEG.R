# Load necessary libraries
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm"))
library(DESeq2)
library(apeglm)
library(ggplot2)

# Set working directory
setwd("/Users/swaraj/Desktop/")
working_dir <- "/Users/swaraj/Desktop/workingdirectory"

if (!dir.exists(working_dir)) {
  dir.create(working_dir, recursive = TRUE)
}

setwd(working_dir)

# Read in the counts data
counts <- read.csv("/Users/swaraj/Desktop/workingdirectory/GSE183947_gene_counts.csv", row.names = 1)
counts[] <- lapply(counts, as.integer)

# Create colData dataframe
condition <- factor(c(rep("control", 30), rep("treated", 30)))
colData <- data.frame(row.names = colnames(counts), condition = condition)

# Create DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)

# Run the DESeq pipeline
dds <- DESeq(dds)

# Get results and order by p-value
res <- results(dds)
resOrdered <- res[order(res$padj), ]

# Write the ordered results to a CSV file
write.csv(as.data.frame(resOrdered), file = "DEGs_results.csv")

# Optional: Plot a volcano plot
volcano <- ggplot(as.data.frame(resOrdered), aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  ggtitle("Volcano Plot") +
  xlab("Log2 Fold Change") +
  ylab("-Log10 Adjusted P-Value") +
  geom_hline(yintercept = -log10(0.05), col = "red", linetype = "dashed") +
  geom_vline(xintercept = c(-1, 1), col = "pink", linetype = "dashed")
print(volcano)

# MA plot
ma_plot <- ggplot(as.data.frame(resOrdered), aes(x = baseMean, y = log2FoldChange)) +
  geom_point(alpha = 0.4) +
  theme_minimal() +
  ggtitle("MA Plot") +
  xlab("Mean of Normalized Counts") +
  ylab("Log2 Fold Change") +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  scale_y_continuous(limits = c(-5, 5))
print(ma_plot)

# Save MA plot
ggsave("MA_plot_from_DEGs.png", plot = ma_plot)


# Assuming 'dds' is already created as in your previous analysis
# Normalize the data
vsd <- vst(dds, blind=FALSE)

# Perform PCA
pcaData <- plotPCA(vsd, intgroup="condition", returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

# Create PCA plot
pca_plot <- ggplot(pcaData, aes(x=PC1, y=PC2, color=condition)) +
  geom_point(size=3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  ggtitle("PCA Plot of RNA-Seq Data") +
  theme_minimal()

# Print the plot
print(pca_plot)

# Optionally, save the plot
ggsave("PCA_plot.png", plot = pca_plot)




if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")

# Load the DESeq2 library
library(DESeq2)

install.packages("pheatmap")
install.packages("RColorBrewer")
library(RColorBrewer)
library(pheatmap)
brewer.pal.info  # This should list all available color palettes


# Assuming 'dds' is already created as in your previous analysis
# Normalize the data
vsd <- vst(dds, blind=FALSE)

# Select the top 20 most differentially expressed genes (by adjusted p-value)
top_genes <- head(order(resOrdered$padj), 20)
top_genes_vsd <- assay(vsd)[top_genes, ]

# Create a heatmap of the top 20 genes
heatmap_colors <- colorRampPalette(rev(brewer.pal(9, "RdBu")))(255)

pheatmap(top_genes_vsd, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         color = heatmap_colors,
         main = "Heatmap of Top 20 Differentially Expressed Genes")

# Optionally, save the heatmap to a file
pheatmap(top_genes_vsd, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         show_rownames = TRUE, 
         show_colnames = TRUE, 
         color = heatmap_colors,
         main = "Heatmap of Top 20 Differentially Expressed Genes",
         filename = "heatmap_top20_genes.png")