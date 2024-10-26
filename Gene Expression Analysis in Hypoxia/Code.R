library("DESeq2")
library("tidyverse")

counts_data = read.csv("raw_counts.csv")
print(counts_data)
head(counts_data)

# Assuming counts_data is already loaded
rownames(counts_data) <- counts_data$ensembl_id
counts_data <- counts_data[, -which(names(counts_data) == "ensembl_id")]

# Verify the changes
print(colnames(counts_data))
print(head(counts_data))



# Sometimes condition columns are not in order 
counts_data <- counts_data[,sort(colnames(counts_data))]
head(counts_data)

counts_data[,-1]

# Total reads per sample 
colSums(counts_data[,-1]) # need to exclude the non-numeric column(s) when calculating column sums. Assuming the first column contains gene names

# Creating the DESeq2 object 
# To create DESeq2 object we need three things. Counts_data, colData, Design. 

# Construct colData
# Identify the biological replicates or conditions
colnames(counts_data)
condition <- c(rep("LNCAP_Hypoxia", 2), rep("LNCAP_Normoxia", 2), rep("PC3_Hypoxia", 2), rep("PC3_Normoxia", 2))

# Assign replicates to each sample name to construct colData

colData <- as.data.frame(condition)
rownames(colData) <- colnames(counts_data)
colData
                          
# Lets cretae the DESeq2 object
# Convert 'condition' to a factor
colData$condition <- factor(colData$condition)

# Now create the DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_data,
                              colData = colData,
                              design = ~condition)

dds <- DESeq(dds)
dds


# We can also sift through the dds object using @ and $. For example, we can find the raw counts data this way
head(dds@assays@data$counts)

# Another way to get the raw counts from the dds object is by using the function counts(dds, normalized = F). If you want the normalized counts instead, we set normalized = T. For example
normalized_counts <- counts(dds, normalized = T)

head(normalized_counts)

# Annontation file
annotation <- read.csv("GRCh38.p13_annotation.csv", header = T, stringsAsFactors = F)
head(annotation)

normalized_counts <- rownames_to_column(as.data.frame(normalized_counts), var = "ensembl_id")
annotated_data <- right_join(annotation, normalized_counts, by = c("Gene.stable.ID" = "ensembl_id"))
head(annotated_data)

results(lncap, contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia"))



## Generate PCA Plot 

## Generate PCA plot
vsd <- vst(dds, blind = FALSE)
plotPCA(vsd, intgroup = "condition")

## Create a heatmap of sample distances
sampleDists <- dist(t(assay(vsd)))
library(pheatmap)
pheatmap(as.matrix(sampleDists))

## MA plot
plotMA(dds, ylim = c(-2, 2))


## Differential Expression Analysis

#Get results for all pairwise comparisons
results_list <- list()
conditions <- levels(colData$condition)
for (i in 1:(length(conditions) - 1)) {
  for (j in (i + 1):length(conditions)) {
    contrast <- c("condition", conditions[i], conditions[j])
    res <- results(dds, contrast = contrast)
    results_list[[paste(conditions[i], "vs", conditions[j])]] <- res
  }
}

## Example: Get results for LNCAP_Hypoxia vs LNCAP_Normoxia
res_LNCAP <- results(dds, contrast = c("condition", "LNCAP_Hypoxia", "LNCAP_Normoxia"))
summary(res_LNCAP)


## Filter and annotate significant genes
## Function to filter and annotate results
filter_and_annotate <- function(res, annotation, padj_cutoff = 0.05, lfc_cutoff = 1) {
  res_df <- as.data.frame(res) %>%
    rownames_to_column("Gene.stable.ID") %>%
    filter(!is.na(padj), padj < padj_cutoff, abs(log2FoldChange) > lfc_cutoff) %>%
    left_join(annotation, by = "Gene.stable.ID")
  return(res_df)
}

## Apply to all results
annotated_results <- lapply(results_list, filter_and_annotate, annotation = annotation)
# Print the names of all comparisons in the results list
print(names(annotated_results))
# Choose one comparison to examine (for example, the first one)
comparison_name <- names(annotated_results)[1]
# Print the first few rows of the chosen comparison
print(head(annotated_results[[comparison_name]]))
# Get a summary of the data frame
summary(annotated_results[[comparison_name]])
# Check the dimensions of the data frame
print(dim(annotated_results[[comparison_name]]))
# Check the column names
print(colnames(annotated_results[[comparison_name]]))



## Visualize results
library(ggplot2)
library(dplyr)

# Assuming 'annotated_results' is your list of dataframes
# Let's use the first comparison as an example
comparison_name <- names(annotated_results)[1]
df <- annotated_results[[comparison_name]]

# Define thresholds
padj_threshold <- 0.05
lfc_threshold <- 1

# Calculate -log10(padj) for plotting
df$`-log10(padj)` <- -log10(df$padj)

# Define significance
df$Significance <- ifelse(df$padj < padj_threshold & abs(df$log2FoldChange) > lfc_threshold, "Significant", "Not Significant")

# Create the enhanced volcano plot
ggplot(df, aes(x = log2FoldChange, y = `-log10(padj)`, color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("blue", "red")) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), color = "green", linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_threshold), color = "blue", linetype = "dashed") +
  labs(title = paste("Enhanced Volcano Plot -", comparison_name),
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  theme_minimal() +
  theme(legend.title = element_blank())

# Optionally, label the top significant genes
top_genes <- df %>% filter(Significance == "Significant") %>% top_n(10, wt = `-log10(padj)`)
ggplot(df, aes(x = log2FoldChange, y = `-log10(padj)`, color = Significance)) +
  geom_point(alpha = 0.7, size = 2) +
  scale_color_manual(values = c("blue", "red")) +
  geom_vline(xintercept = c(-lfc_threshold, lfc_threshold), color = "green", linetype = "dashed") +
  geom_hline(yintercept = -log10(padj_threshold), color = "black", linetype = "dashed") +
  labs(title = paste("Enhanced Volcano Plot -", comparison_name),
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted p-value") +
  theme_minimal() +
  theme(legend.title = element_blank()) +
  geom_text(data=top_genes, aes(label=Gene.name), vjust=1.5, hjust=1.5, size=3)


## Pathway Enrichment Analysis
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggplot2)
library(gridExtra)

# Function for GO enrichment
perform_GO_enrichment <- function(gene_list, ont = "BP") {
  ego <- enrichGO(gene = gene_list,
                  OrgDb = org.Hs.eg.db,
                  keyType = "ENSEMBL",
                  ont = ont,
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2)
  return(ego)
}

# Apply to each comparison
go_results <- lapply(annotated_results, function(res) {
  perform_GO_enrichment(res$Gene.stable.ID)
})

# Visualize top GO terms with improved clarity
plot_top_GO <- function(ego, title) {
  dotplot(ego, showCategory = 5, title = title) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 16),
      axis.text.y = element_text(size = 10),
      axis.text.x = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "right"
    ) +
    labs(y = "Pathway Names")
}

# Generate and store plots for each comparison
go_plots <- lapply(names(go_results), function(name) {
  plot_top_GO(go_results[[name]], name)
})

# Arrange the plots in a grid
grid.arrange(grobs = go_plots, ncol = 2, nrow = ceiling(length(go_plots) / 2))