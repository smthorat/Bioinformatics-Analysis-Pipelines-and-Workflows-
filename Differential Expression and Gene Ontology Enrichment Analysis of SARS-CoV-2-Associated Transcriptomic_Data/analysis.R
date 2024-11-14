# Load necessary libraries
library(DESeq2)
library(ggplot2)
library(ggrepel)

# Load the count data
counts <- read.delim("counts.txt", sep = "\t", header = TRUE, comment.char = "#")
head(counts)  # Preview the data

# Rename BAM file columns to sample 1, sample 2 and so on. For better readability and interpretation
bam_columns <- grep("^X", colnames(counts))
colnames(counts)[bam_columns] <- paste0("Sample", 1:length(bam_columns))

# Extract only the sample data
data_samples <- counts[, grepl("^Sample", colnames(counts))]
rownames(data_samples) <- counts$GeneID  # Set rownames to gene IDs

metadata = read.csv("assignment_2_info.csv")
metadata

colnames(metadata)
row.names(metadata)

# Remove unwanted columns
colData <- colData[, !(names(colData) %in% c("X", "X.1", "X.2", "X.3", "X.4", "X.5"))]

all(rownames(metadata) == colnames(data_samples)) 
dim(data_samples)  # Returns number of rows and columns
dim(metadata)       # Returns number of rows and columns

rownames(metadata) <- colnames(data_samples)
rownames(metadata) <- colnames(data_samples)
all(rownames(metadata) == colnames(data_samples)) 


# Pre-processing and Quality Control
# Check data dimensions
dim(data_samples)
dim(metadata)

# Summarize the data
summary(data_samples)

# Remove low-count genes
keep <- rowSums(data_samples >= 10) >= 2  # Keep genes with at least 10 counts in at least 2 samples
filtered_data <- data_samples[keep, ]
dim(filtered_data)  # Reduced dataset dimensions

# Create DESeq Dataset
dds <- DESeqDataSetFromMatrix(countData = filtered_data,
                              colData = metadata,
                              design = ~ Condition + Time.Point)

# Regularized log transformation (for visualization)
rld <- rlog(dds)

# Perform Differential Expression Analysis
# Analyze for condition SARS-CoV-2 vs Mock
dds <- DESeq(dds)
results_condition <- results(dds, contrast = c("Condition", "SARS-CoV-2", "Mock "))

# Ensure Time.Point is a factor
metadata$Time.Point <- factor(metadata$Time.Point, levels = c("24", "72"))

# Create DESeq Dataset with updated metadata
dds <- DESeqDataSetFromMatrix(countData = filtered_data,
                              colData = metadata,
                              design = ~ Condition + Time.Point)

# DESeq analysis
dds <- DESeq(dds)

# Analyze results for time points (72 H vs 24 H)
results_time <- results(dds, contrast = c("Time.Point", "72", "24"))

# Analyze for time points (72 H vs 24 H)
results_time <- results(dds, contrast = c("Time.Point", "72", "24"))

# Extract Significant Genes
# Significant genes for condition
significant_genes_condition <- subset(results_condition, padj < 0.05)
write.csv(as.data.frame(significant_genes_condition), "significant_genes_condition.csv")

# Significant genes for time points
significant_genes_time <- subset(results_time, padj < 0.05)
write.csv(as.data.frame(significant_genes_time), "significant_genes_time.csv")

# Visualization
# MA Plot for condition
plotMA(results_condition, main = "MA Plot for SARS-CoV-2 vs Mock", ylim = c(-5, 5))

# MA Plot for time points
plotMA(results_time, main = "MA Plot for 72 H vs 24 H", ylim = c(-5, 5))

# PCA Plot
plotPCA(rld, intgroup = c("Condition"))
plotPCA(rld, intgroup = c("Time.Point"))

# Volcano Plot for condition
volcano_data <- as.data.frame(results_condition)
volcano_data$Significance <- ifelse(volcano_data$padj < 0.05, "Significant", "Not Significant")
volcano_data$Label <- ifelse(volcano_data$padj < 0.05 & abs(volcano_data$log2FoldChange) > 2, rownames(volcano_data), NA)

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(aes(label = Label), max.overlaps = 10, size = 3) +
  scale_color_manual(values = c("Not Significant" = "gray", "Significant" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot for SARS-CoV-2 vs Mock",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.title = element_blank())

# Volcano Plot for time points
volcano_data_time <- as.data.frame(results_time)
volcano_data_time$Significance <- ifelse(volcano_data_time$padj < 0.05, "Significant", "Not Significant")
volcano_data_time$Label <- ifelse(volcano_data_time$padj < 0.05 & abs(volcano_data_time$log2FoldChange) > 2, rownames(volcano_data_time), NA)

ggplot(volcano_data_time, aes(x = log2FoldChange, y = -log10(padj), color = Significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_text_repel(aes(label = Label), max.overlaps = 10, size = 3) +
  scale_color_manual(values = c("Not Significant" = "black", "Significant" = "orange")) +
  theme_minimal() +
  labs(title = "Volcano Plot for 72 H vs 24 H",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value") +
  theme(legend.title = element_blank())

# Save Results
# Save normalized data
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(normalized_counts), "normalized_counts.csv")

## GO Term Enrichment Analysis

# Install required packages
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"))

# Load necessary libraries
library(clusterProfiler)
library(org.Hs.eg.db)  # Human annotation database
library(DOSE)  # For visualization

significant_genes_condition <- rownames(subset(results_condition, padj < 0.05))
significant_genes_time <- rownames(subset(results_time, padj < 0.05))

# GO Enrichment for Condition
go_condition <- enrichGO(gene = significant_genes_condition,  # Directly use these
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",  # Ensure this matches your identifier type
                         ont = "BP",  # Biological Process
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.05,
                         qvalueCutoff = 0.05)


# GO Enrichment for Time Point
go_time <- enrichGO(gene = significant_genes_time,
                    OrgDb = org.Hs.eg.db,
                    keyType = "ENTREZID", 
                    ont = "BP",  
                    pAdjustMethod = "BH",
                    pvalueCutoff = 0.05,
                    qvalueCutoff = 0.05)


# Barplot for Condition
barplot(go_condition, showCategory = 20, title = "GO Enrichment: Condition (SARS-CoV-2 vs Mock)")
head(go_condition)

length(go_condition)

# Relax the thresholds (Not able to find significant genes)
go_condition <- enrichGO(gene = condition_entrez$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",  
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.1,  
                         qvalueCutoff = 0.2)
# Check the number of enriched terms
length(go_condition)

# Preview the enriched terms
head(as.data.frame(go_condition))

# Barplot for Time Point
barplot(go_time, showCategory = 20, title = "GO Enrichment: Time Point (72 H vs 24 H)")

# Barplot with relaxed thresholds
if (length(go_condition) > 0) {
  barplot(go_condition, showCategory = 20, title = "GO Enrichment: Condition (SARS-CoV-2 vs Mock)")
} else {
  message("No significant GO terms found, even with relaxed thresholds.")
}


background_genes <- rownames(filtered_data)  # Use all genes from the dataset
go_condition <- enrichGO(gene = condition_entrez$ENTREZID,
                         OrgDb = org.Hs.eg.db,
                         keyType = "ENTREZID",
                         ont = "BP",
                         pAdjustMethod = "BH",
                         pvalueCutoff = 0.1,
                         qvalueCutoff = 0.2,
                         universe = background_genes)  # Define custom background

# Dotplot with relaxed thresholds
dotplot(go_condition, showCategory = 20, title = "GO Enrichment: Condition (SARS-CoV-2 vs Mock)")


# Save GO Enrichment Results for Condition
write.csv(as.data.frame(go_condition), "GO_Enrichment_Condition.csv")

# Save GO Enrichment Results for Time Point
write.csv(as.data.frame(go_time), "GO_Enrichment_Time_Point.csv")


