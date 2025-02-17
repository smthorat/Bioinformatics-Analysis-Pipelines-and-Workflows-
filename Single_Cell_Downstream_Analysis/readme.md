
# Spatial Transcriptomics Analysis Report

This report outlines the analysis performed on the 10x Genomics Visium dataset (human breast cancer) using the Seurat pipeline in R. The analysis covers data loading, quality control, normalization, dimensionality reduction, clustering, spatial feature analysis, and cell type annotation.

## Download the data: https://www.10xgenomics.com/datasets/human-breast-cancer-visium-fresh-frozen-whole-transcriptome-1-standard

## 1. Data Loading

- **Gene Expression Data:**  
  The gene expression data was loaded using the `Read10X` function, creating a counts matrix.
  
- **Spatial Metadata:**  
  The spatial image data (high-resolution tissue image) and associated metadata were loaded using `Read10X_Image`, linking spatial coordinates to the expression data.

- **Seurat Object Creation:**  
  A Seurat object was created with the spatial assay, and the spatial data was added to this object. Metadata such as total UMI counts per spot and mitochondrial gene percentages were calculated and added.

## 2. Quality Control (QC)

- **UMI Counts and Feature Detection:**  
  The `SpatialFeaturePlot` and `VlnPlot` functions were used to visualize the distribution of total UMI counts and the number of detected features (genes) per spot.  
  These plots help identify outlier spots and ensure data quality.

- **Mitochondrial Gene Percentage:**  
  The percentage of mitochondrial gene expression was calculated using the `PercentageFeatureSet` function. Spots with more than 20% mitochondrial reads were filtered out, as high mitochondrial content can indicate low-quality or dying cells.

- **Filtering Criteria:**  
  Spots were retained based on the following criteria:
  - More than 200 features detected
  - Fewer than 6000 features detected
  - Less than 20% mitochondrial gene content

## 3. Data Normalization and Feature Selection

- **Normalization with SCTransform:**  
  The SCTransform method was applied to normalize the data. This step also identifies highly variable genes (HVGs) that capture the biological variability in the dataset.
  
- **Identification of Highly Variable Genes (HVGs):**  
  HVGs were automatically identified during the SCTransform step and are used for subsequent analyses.

## 4. Dimensionality Reduction and Clustering

- **Principal Component Analysis (PCA):**  
  PCA was run on the Seurat object using 30 principal components (PCs).  
  - The top genes contributing to the variation in each of the first five PCs were examined using print and visualization commands.
  - Loadings for the first two PCs were visualized to highlight the genes that drive variation.

- **Clustering:**  
  - A nearest-neighbor graph was constructed using the first 20 PCs.
  - Clusters were identified at a resolution of 0.5 using `FindClusters`.
  - Cluster assignments were inspected via metadata and visualized in PCA space.

- **UMAP Embedding:**  
  UMAP was computed using the first 20 PCs to obtain a two-dimensional embedding of the data.  
  - UMAP plots were generated, both uncolored and colored by cluster assignment, to reveal distinct cell populations.

## 5. Spatial Analysis

- **Extraction and Alignment of Spatial Coordinates:**  
  Spatial coordinates were read from a CSV file (`tissue_positions.csv`), aligned with the Seurat object’s cell barcodes, and added as metadata.  
  This ensures that each spot’s gene expression data is correctly mapped to its physical location in the tissue.

- **Identification of Spatially Variable Genes:**  
  Spatially variable genes were identified using the `FindSpatiallyVariableFeatures` function with the “moransi” method.  
  - The top spatially variable genes were visualized with `SpatialFeaturePlot`, highlighting regions of the tissue where these genes exhibit differential expression.

## 6. Cell Type Annotation

- **Marker Gene Identification:**  
  Marker genes for each cluster were identified using `FindAllMarkers` on the normalized data (assay "SCT").

- **Annotation of Clusters:**  
  Based on marker gene expression, clusters were manually annotated:
  - Cluster "0" was labeled as **Neurons**
  - Cluster "1" as **Astrocytes**
  - Cluster "2" as **Oligodendrocytes**
  
- **Visualization of Annotated Clusters:**  
  The annotated clusters were visualized using `SpatialDimPlot`, providing spatial context to the cell type distribution.

## 7. Advanced Spatial Feature Analysis

- **Gene-Specific Spatial Expression:**  
  - The script checked for the presence of key genes (e.g., `GAD1`, `SLC17A6`/`SLC17A7`) in the dataset.
  - Spatial expression of the gene `GAD1` was visualized to further understand its spatial distribution in the tissue.

## 8. Considerations for Data Completeness

The analysis includes all the critical components:
- **High-resolution and Low-resolution Tissue Images:**  
  The high-resolution tissue image is linked to the spatial metadata, and the corresponding low-resolution image is available for alignment and quality control.
- **Scale Factors and Fiducial Markers:**  
  The `scalefactors_json.json` file provides essential scale information, ensuring proper calibration between image pixels and actual spatial dimensions.
- **Tissue Positions:**  
  The CSV file (`tissue_positions.csv`) contains the spatial coordinates required to align the gene expression data with the tissue morphology.
  
**Note:** Ensure consistency in file formats and naming conventions. For example, if your script references a `.png` file for the high-resolution image, verify that the file is available in the correct format.

## 9. Conclusion

The spatial transcriptomics analysis pipeline successfully integrates:
- **QC and Filtering:** Ensuring high-quality spots for analysis.
- **Normalization and Feature Selection:** Highlighting biologically relevant variation.
- **Dimensionality Reduction and Clustering:** Revealing distinct cellular subpopulations.
- **Spatial Analysis and Cell Type Annotation:** Providing insight into tissue architecture and gene expression patterns.

This comprehensive analysis allows for in-depth exploration of the spatial dynamics within the tissue, identifying both global patterns and localized gene expression features that can be linked to specific cell types.

---

*This report serves as documentation for the analysis workflow and the results obtained. Adjustments can be made based on further experimental needs or data quality assessments.*
