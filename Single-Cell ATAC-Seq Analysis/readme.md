# Single-Cell ATAC-Seq Analysis of Human PBMCs

## Introduction
Single-cell Assay for Transposase-Accessible Chromatin using sequencing (**scATAC-seq**) is a technique used to study chromatin accessibility at the single-cell level. This analysis focuses on a **publicly available dataset of human Peripheral Blood Mononuclear Cells (PBMCs)** provided by 10x Genomics.

### **Objectives**
- Process and structure single-cell ATAC-seq data.
- Perform quality control (**QC**) to remove low-quality cells.
- Conduct **dimensionality reduction and clustering** to identify distinct cell populations.
- Interpret **chromatin accessibility profiles** based on clustering results.

---

## **Step 1: Data Processing**
### **Data Sources**
The dataset includes:
1. **Fragments File** – Contains transposase-inserted reads.
2. **Filtered Peak-Barcode Matrix** – Maps chromatin accessibility peaks to cells.
3. **Metadata File** – Provides quality control metrics per cell.

### **Chromatin Assay Creation**
- The data is structured into a **Chromatin Assay**, a Seurat object designed for chromatin accessibility data.
- **Gene annotations** are added using the EnsDb.Hsapiens.v75 database to map peaks to known genes.

---

## **Step 2: Quality Control (QC)**
### **Key QC Metrics**
- **Nucleosome Signal Score** – Measures the presence of nucleosome-bound vs. nucleosome-free fragments.
- **TSS (Transcription Start Site) Enrichment Score** – Indicates the capture efficiency of active TSS.
- **Blacklist Ratio** – Measures the fraction of reads mapping to known **blacklisted genomic regions**.
- **Percentage of Reads in Peaks** – The proportion of reads found in accessible chromatin regions.

### **QC Visualization**
#### **Before Filtering:**
- **Violin Plot:** [Violin plot before filtering](Violin%20plot.png)
- **Density Scatter Plots:** [Quantile scatter density plot](Rplot.png)

---

## **Step 3: Filtering Low-Quality Cells**
Cells were filtered based on:
- `nCount_ATAC > 3000` (remove low-coverage cells).
- `nCount_ATAC < 30000` (remove doublets).
- `pct_reads_in_peaks > 15%`.
- `blacklist_ratio < 0.05`.
- `nucleosome_signal < 4`.
- `TSS.enrichment > 3`.

### **Post-Filtering Quality Check**
- 📌 **Violin Plot After Filtering:** [Filtered QC Metrics](Quality%20control%20metrics.jpeg)

---

## **Step 4: Dimensionality Reduction and Clustering**
### **Normalization**
- **TF-IDF normalization** was applied to standardize data.
- **Feature selection** was performed to retain the most informative chromatin regions.

### **Correlation Between Sequencing Depth and Principal Components**
- **Correlation Plot:** [Correlation with sequencing depth](Correlation%20depth.jpeg)

### **UMAP Embedding and Clustering**
- **UMAP (Uniform Manifold Approximation and Projection)** was used for clustering.
- Cells were grouped into **distinct clusters** based on chromatin accessibility.

- **UMAP Visualization:** [Clustered UMAP Plot](Umap.jpeg)

---

## **Step 5: Results and Interpretation**
- **Distinct cell populations** were identified based on **chromatin accessibility**.
- **High-quality cells were retained** after filtering.
- **UMAP visualization shows well-separated clusters**, indicating different regulatory landscapes in PBMC subtypes.

---

## **Conclusion**
### **Key Findings**
**Effective Quality Control**: Low-quality cells were successfully removed.  
**Identification of Cell Clusters**: PBMCs were grouped based on chromatin accessibility.  
**Potential Insights into Gene Regulation**: Chromatin accessibility correlates with distinct immune cell types.

### **Next Steps**
**Cluster Annotation** – Identify immune cell types using marker peaks.  
**Differential Accessibility Analysis** – Compare chromatin accessibility across clusters.  
**Integration with Gene Expression Data** – Link chromatin accessibility with transcriptomic profiles.

---

## **References**
**Tutorial Reference:** [10x Genomics PBMC Vignette](https://stuartlab.org/signac/articles/pbmc_vignette#non-linear-dimension-reduction-and-clustering)

---
