library(Signac)
library(Seurat)
library(EnsDb.Hsapiens.v75)
library(tidyverse)
library(biovizBase)


frag.file <- read.delim('10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz', header = F, nrows = 10)
head(frag.file)

# 1. Read in data -----------------

counts <- Read10X_h5('10k_pbmc_ATACv2_nextgem_Chromium_Controller_filtered_peak_bc_matrix.h5')
counts[1:10,1:10]

# Creating chromatics assay, it contains inof on chromosomes, annotations, fragments etc
chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  fragments = "10k_pbmc_ATACv2_nextgem_Chromium_Controller_fragments.tsv.gz",
  min.cells = 10,
  min.features = 200
)

str(chrom_assay)

# Reading metadata 

metadata <- read.csv(file = '10k_pbmc_ATACv2_nextgem_Chromium_Controller_singlecell.csv', header = T, row.names = 1)
View(metadata)

# create a seurat Object
pbmc <- CreateSeuratObject(
  counts = chrom_assay,
  meta.data = metadata,
  assay = 'ATAC'
)

str(pbmc)

# ....Adding Gene Annotation -------------------

pbmc@assays$ATAC@annotation
# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)


# change to UCSC style since the data was mapped to hg19
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))


# add the gene information to the object
Annotation(pbmc) <- annotations
pbmc@assays$ATAC@annotation


# 2. Computing QC ---------------------

# compute nucleosome signal score per cell
pbmc <- NucleosomeSignal(pbmc)

# compute TSS enrichment score per cell
pbmc <- TSSEnrichment(object = pbmc, fast = FALSE)

# add blacklist ratio and fraction of reads in peaks
pbmc$blacklist_ratio <- pbmc$blacklist_region_fragments / pbmc$peak_region_fragments
pbmc$pct_reads_in_peaks <- pbmc$peak_region_fragments / pbmc$passed_filters * 100

View(pbmc@meta.data)
# Save the metadata to a CSV file
write.csv(pbmc@meta.data, file = "pbmc_metadata.csv", row.names = TRUE)


# ....Visualizing QC --------------------
# Save the column names from pbmc@meta.data into a CSV file
# Save the column names from pbmc@meta.data into a CSV file
meta_colnames <- colnames(pbmc@meta.data)
write.csv(meta_colnames, file = "pbmc_meta_colnames.csv", row.names = FALSE)

# Create density scatter plots
a1 <- DensityScatter(pbmc, x = 'nCount_ATAC', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)
a2 <- DensityScatter(pbmc, x = 'nucleosome_signal', y = 'TSS.enrichment', log_x = TRUE, quantiles = TRUE)

# Combine the two density scatter plots
combined_plot <- a1 | a2
print(combined_plot)

# Generate a violin plot for selected features
violin_plot <- VlnPlot(object = pbmc, 
                       features = c('nCount_ATAC', 'nFeature_ATAC', 'TSS.enrichment', 
                                    'nucleosome_signal', 'blacklist_ratio', 'pct_reads_in_peaks'),
                       pt.size = 0.1,
                       ncol = 6)
print(violin_plot)


# ....Filtering poor quality cells --------------------

# 3. Normalization and linear dimensional reduction ------------------
pbmc <- RunTFIDF(pbmc) # normalization
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0') # selecting top features
pbmc <- RunSVD(pbmc) # dimensionality reduction

DepthCor(pbmc)



# Number of cells before filtering
cat("Number of cells before filtering:", ncol(pbmc), "\n")

# Violin plots before filtering
p_before <- VlnPlot(
  pbmc,
  features = c("nCount_ATAC", "pct_reads_in_peaks", "blacklist_ratio",
               "nucleosome_signal", "TSS.enrichment"),
  ncol = 5,
  pt.size = 0.2
) + NoLegend()

# Print the 'before' plot
print(p_before)

# -------------------Filtering poor quality cells --------------------
pbmc_filtered <- subset(
  x = pbmc,
  subset = nCount_ATAC > 3000 &
    nCount_ATAC < 30000 &
    pct_reads_in_peaks > 15 &
    blacklist_ratio < 0.05 &
    nucleosome_signal < 4 &
    TSS.enrichment > 3
)


# Number of cells after filtering
cat("Number of cells after filtering:", ncol(pbmc_filtered), "\n")

# Violin plots after filtering
p_after <- VlnPlot(
  pbmc_filtered,
  features = c("nCount_ATAC", "pct_reads_in_peaks", "blacklist_ratio",
               "nucleosome_signal", "TSS.enrichment"),
  ncol = 5,
  pt.size = 0.2
) + NoLegend()

# Print the 'after' plot
print(p_after)

# 3. Normalization and linear dimensional reduction ------------------
pbmc <- RunTFIDF(pbmc) # normalization
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0') # selecting top features
pbmc <- RunSVD(pbmc) # dimensionality reduction

DepthCor(pbmc)


# 4. Non-linear dimensional reduction and Clustering -------------------

# See what reductions exist
Reductions(pbmc)

# Check the dimensions of the LSI embedding
dim(Embeddings(pbmc, reduction = "lsi"))
# 4. Non-linear dimensional reduction and Clustering -------------------
# Run UMAP with dims 2:15 and set n.neighbors to a value smaller than the cell count (e.g., 5)
pbmc <- RunUMAP(object = pbmc, reduction = 'lsi', dims = 2:15, n.neighbors = 5)

# Find neighbors using dims 2:15
pbmc <- FindNeighbors(object = pbmc, reduction = 'lsi', dims = 2:15)

# Find clusters (using algorithm 3)
pbmc <- FindClusters(object = pbmc, algorithm = 3)

# Plot the clusters with labels and remove the legend
DimPlot(object = pbmc, label = TRUE) + NoLegend()



# Run normalization and dimensionality reduction on the full dataset (pbmc)
pbmc <- RunTFIDF(pbmc)
pbmc <- FindTopFeatures(pbmc, min.cutoff = 'q0')
pbmc <- RunSVD(pbmc)

# Check the dimensions of the LSI embedding
lsi_dims <- dim(Embeddings(pbmc, reduction = "lsi"))
cat("LSI dimensions:", lsi_dims[1], "x", lsi_dims[2], "\n")

# Run UMAP on the unfiltered dataset using the LSI reduction (skip the first component)
pbmc <- RunUMAP(pbmc, reduction = 'lsi', dims = 2:lsi_dims[2], n.neighbors = 30)

# Optionally, run neighbors and clustering on the unfiltered dataset
pbmc <- FindNeighbors(pbmc, reduction = 'lsi', dims = 2:lsi_dims[2])
pbmc <- FindClusters(pbmc, algorithm = 3, resolution = 0.5)

# Plot the UMAP with cluster labels
DimPlot(pbmc, label = TRUE) + NoLegend()
































