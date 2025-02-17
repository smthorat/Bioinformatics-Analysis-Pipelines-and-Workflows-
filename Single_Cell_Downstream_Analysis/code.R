library(Seurat)
library(ggplot2)

# Set data directory
data_dir <- "/Volumes/Jagannath/Projects/single_cell_downstream/Downstream/data/"

# 1. Load data ------------------------------------------------------------
# Load gene expression data
data <- Read10X(paste0(data_dir, "filtered_feature_bc_matrix/"))

# Load spatial metadata
spatial_data <- Read10X_Image(
  image.dir = paste0(data_dir, "spatial/"),
  image.name = "tissue_hires_image.png"
)

# 2. Create Seurat object -------------------------------------------------
seurat_obj <- CreateSeuratObject(counts = data, assay = "Spatial")

# 3. Add spatial data -----------------------------------------------------
seurat_obj[["slice1"]] <- spatial_data
seurat_obj@images$slice1@assay <- "Spatial"  # Link to assay

# 4. Add metadata ---------------------------------------------------------
seurat_obj$nCount_Spatial <- colSums(seurat_obj[["Spatial"]]$counts)

# 5. Plot ------------------------------------------------------------------
SpatialFeaturePlot(seurat_obj, features = "nCount_Spatial") +
  ggtitle("Total UMI Counts per Spot") +
  scale_fill_gradient(low = "blue", high = "red", limits = c(0, 20000))

# Calculate mitochondrial gene percentage
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

head(seurat_obj@meta.data)

# Visualize QC metrics
VlnPlot(seurat_obj, features = c("nCount_Spatial", "nFeature_Spatial", "percent.mt"))

# Copy the original metrics to new columns
seurat_obj$nFeature_Spatial <- seurat_obj$nFeature_RNA
seurat_obj$nCount_Spatial <- seurat_obj$nCount_RNA

# Now subset using the new names
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_Spatial > 200 & 
                       nFeature_Spatial < 6000 & 
                       percent.mt < 20)

# Normalization and Feature Selection -------------------------------------------------------------------------------
# Normalize (SCTransform is preferred for spatial data)
seurat_obj <- SCTransform(seurat_obj, assay = "Spatial", verbose = FALSE)

#Check how data look like
head(GetAssayData(seurat_obj, assay = "SCT", slot = "data"))

# Identify HVGs (automatically done by SCTransform) Highly Variable Genes
VariableFeatures(seurat_obj)  # View top HVGs


# Diamentionility Reduction and clustering------------------------------------------------------------------------------
# PCA
# Run principal component analysis (PCA) on the Seurat object,
# computing 30 principal components (PCs) based on the variable features.
seurat_obj <- RunPCA(seurat_obj, npcs = 30)

# Print a summary of the PCA results:
# - For the first 5 principal components (dims = 1:5)
# - Displaying the top 5 genes (nfeatures = 5) that contribute to each component.
print(seurat_obj[["pca"]], dims = 1:5, nfeatures = 5)

# Visualize the loadings (contribution) of genes for the first 2 PCs:
# This plot shows which genes have the highest influence on the first two principal components.
VizDimLoadings(seurat_obj, dims = 1:2, reduction = "pca")

# Generate a scatter plot of the cells (or spots) in PCA space:
# This DimPlot visualizes how cells are distributed based on their PCA scores,
# which can help reveal underlying clusters or patterns in the data.
DimPlot(seurat_obj, reduction = "pca")

# Clustering
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:20)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)

# View the first few rows of the metadata to see the cluster assignments
head(seurat_obj@meta.data)

# Alternatively, view just the cluster assignments
head(seurat_obj$seurat_clusters)

# Get a summary table of how many cells are in each cluster
table(seurat_obj$seurat_clusters)

# Visualize clusters in PCA space
DimPlot(seurat_obj, reduction = "pca", group.by = "seurat_clusters")

# UMAP/TSNE
# Run UMAP on the Seurat object using the first 20 principal components (PCs)
# This computes a two-dimensional embedding (UMAP coordinates) for each cell.
seurat_obj <- RunUMAP(seurat_obj, dims = 1:20)

# Visualize the UMAP embedding without any additional grouping
# This plot displays each cell based solely on the UMAP coordinates.
DimPlot(seurat_obj, reduction = "umap")

# Visualize the UMAP embedding with cells colored by their cluster assignments
# Here, the cells are grouped by the 'seurat_clusters' metadata, highlighting the clusters identified during clustering.
DimPlot(seurat_obj, reduction = "umap", group.by = "seurat_clusters")


#Spatial-Specific Analyses--------------------------------------------------------------------

spatially_variable_genes <- Seurat::FindSpatiallyVariableFeatures(
  seurat_obj,
  assay = "SCT",
  features = VariableFeatures(seurat_obj),
  selection.method = "moransi"
)

# Plot top spatially variable genes
SpatialFeaturePlot(seurat_obj, features = head(spatially_variable_genes, 3))

#Diagnose----------------------------
n_cells <- ncol(GetAssayData(seurat_obj, assay = "SCT", slot = "data"))
print(n_cells)

n_spots <- nrow(seurat_obj@images$slice1@coordinates)
print(n_spots)

# Extract the x and y coordinates (make sure the order is x then y)
coords <- spatial_coords[, c("X12203", "X4650")]

# Optionally, check the first few rows. (You face issue because you are passing coordinates internally as well as manually)
head(coords)


# Read the CSV file
spatial_coords <- read.csv(
  paste0(data_dir, "spatial/tissue_positions.csv"),
  header = TRUE,
  row.names = 1  # Ensure barcodes are row names
)

# Extract x/y coordinates (adjust column names to match your CSV)
coords <- spatial_coords[, c("X12203", "X4650")]  # x = pxl_col, y = pxl_row
colnames(coords) <- c("x", "y")  # Standardize names


# Subset coordinates to match cells in the Seurat object
cells_in_seurat <- colnames(seurat_obj)
coords_aligned <- coords[cells_in_seurat, ]

# Check alignment
nrow(coords_aligned) == ncol(seurat_obj)  # Should return TRUE


spatially_variable_genes <- Seurat::FindSpatiallyVariableFeatures(
  seurat_obj,
  assay = "SCT",
  features = VariableFeatures(seurat_obj),
  selection.method = "moransi",
  spatial.location = coords_aligned  # Explicitly provide coordinates
)

# Add aligned coordinates to the Seurat object's metadata
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, coords_aligned)

# Add the coordinates to the metadata (if not already added)
seurat_obj <- AddMetaData(seurat_obj, metadata = coords_aligned)

spatially_variable_genes <- Seurat::FindSpatiallyVariableFeatures(
  seurat_obj,
  assay = "SCT",
  features = VariableFeatures(seurat_obj),
  selection.method = "moransi",
  spatial.location = c("X12203", "X4650")
)


#CellType annotation------------------------------------------------
# Find cluster markers
markers <- FindAllMarkers(seurat_obj, assay = "SCT", only.pos = TRUE)

# Example: Assign cell types (customize based on your markers)
seurat_obj <- RenameIdents(seurat_obj, 
                           "0" = "Neurons", 
                           "1" = "Astrocytes", 
                           "2" = "Oligodendrocytes"
)

# Plot annotated clusters
SpatialDimPlot(seurat_obj, label = TRUE)

#Advanced Spatial Analysis---------------------------------------------------------------

## Cell cell interaction
# Spatial UMAP overlay
SpatialDimPlot(seurat_obj, combine = FALSE)

# Co-expression of two genes
# Check the default assay
DefaultAssay(seurat_obj)

# Look at a few row names
head(rownames(seurat_obj[[DefaultAssay(seurat_obj)]]))

# See if "GAD1" and "SLC17A6" are present
"GAD1" %in% rownames(seurat_obj[[DefaultAssay(seurat_obj)]])
"SLC17A6" %in% rownames(seurat_obj[[DefaultAssay(seurat_obj)]])

DefaultAssay(seurat_obj)


"GAD1" %in% rownames(seurat_obj[[DefaultAssay(seurat_obj)]])
# returns TRUE
"SLC17A6" %in% rownames(seurat_obj[[DefaultAssay(seurat_obj)]])
# returns FALSE

grep("SLC17A6", rownames(seurat_obj[[DefaultAssay(seurat_obj)]]), value = TRUE)

"SLC17A7" %in% rownames(seurat_obj[[DefaultAssay(seurat_obj)]])
p1 <- SpatialFeaturePlot(seurat_obj, features = "GAD1")
SpatialFeaturePlot(seurat_obj, features = "GAD1")

p1 <- SpatialFeaturePlot(seurat_obj, features = "GAD1")
SpatialFeaturePlot(seurat_obj, features = "GAD1")