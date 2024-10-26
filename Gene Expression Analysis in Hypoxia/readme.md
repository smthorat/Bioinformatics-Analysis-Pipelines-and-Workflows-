# Differential Gene Expression and Pathway Analysis

This project performs differential gene expression analysis and pathway enrichment analysis on RNA-seq data, comparing hypoxia and normoxia conditions in two prostate cancer cell lines, LNCaP and PC3. The analysis uses the DESeq2 package for identifying differentially expressed genes (DEGs) and the clusterProfiler package for functional enrichment analysis.

## Table of Contents
- [Introduction](#introduction)
- [Data Description](#data-description)
- [Analysis Pipeline](#analysis-pipeline)
- [Installation](#installation)
- [Usage](#usage)
- [Results](#results)
- [Visualization](#visualization)
- [File Structure](#file-structure)
- [References](#references)

## Introduction
The goal of this project is to explore the gene expression changes under different oxygen conditions (hypoxia vs. normoxia) in two prostate cancer cell lines (LNCaP and PC3). Understanding these changes could help reveal pathways involved in the response to hypoxia, which is relevant to cancer progression and therapy resistance.

## Data Description
The analysis requires two main data files:
1. **Raw Count Data (`raw_counts.csv`)**: Contains raw gene expression counts for each sample. The rows correspond to genes (Ensembl IDs), and the columns represent different samples.
2. **Gene Annotation File (`GRCh38.p13_annotation.csv`)**: Provides annotations for the genes, including Ensembl IDs and gene names.

The data should be structured as follows:
- **Raw Counts Data**: The first column should contain Ensembl gene IDs, followed by columns representing different samples. Each cell contains the raw read counts for a gene in a given sample.
- **Annotation File**: Must include at least two columns: "Gene.stable.ID" (Ensembl ID) and "Gene.name" (gene name).

## Analysis Pipeline
The analysis follows these steps:

1. **Data Preparation**:
   - Import the raw counts and annotation data.
   - Format the data to make it suitable for analysis (set Ensembl IDs as row names, sort columns).

2. **DESeq2 Analysis**:
   - Create a `DESeqDataSet` object using the count data and sample conditions.
   - Run the DESeq2 pipeline to normalize the data and perform differential expression analysis.
   - Obtain lists of differentially expressed genes based on adjusted p-value and log2 fold change criteria.

3. **Visualization**:
   - Generate PCA plots to check for clustering based on conditions.
   - Plot heatmaps showing sample distances.
   - Create MA plots to visualize differential expression.

4. **Pathway Enrichment Analysis**:
   - Perform Gene Ontology (GO) enrichment analysis using `clusterProfiler`.
   - Analyze the top enriched pathways for Biological Processes (BP).

5. **Results Saving**:
   - Save lists of DEGs and enriched pathways as CSV files.
   - Export visualizations (PCA plots, heatmaps, GO enrichment plots).

## Installation
To run this analysis, make sure you have R and the required packages installed. You can install the necessary packages with:

```R
install.packages("BiocManager")
BiocManager::install(c("DESeq2", "clusterProfiler", "org.Hs.eg.db", "pheatmap"))
install.packages(c("tidyverse", "ggplot2", "gridExtra"))
