This tutorial is for Bioinformatics class B students

# Single-Cell Transcriptomic Profiling of PBMCs using Seurat

Hi there! Welcome to my single-cell RNA sequencing project on PBMCs (Peripheral Blood Mononuclear Cells) using the **Seurat** package in R. This project showcases how to analyze scRNA-seq data, from initial data loading all the way to visualizing and interpreting cell clusters.

## Project Overview

The goal of this project is to understand the diversity within PBMCs by identifying different cell types and observing gene expression patterns across individual cells. Single-cell RNA sequencing (scRNA-seq) is super powerful for this because it allows us to see the variability in gene expression at a single-cell level, which is way more detailed than traditional bulk RNA sequencing.

Here’s a quick rundown of the steps we’ll go through:
1. Load the data and prepare it for analysis.
2. Clean the data by filtering out low-quality cells.
3. Normalize and identify genes that show high variability.
4. Cluster the cells based on their gene expression patterns.
5. Identify and visualize the clusters with different cell types.

## Data Requirements
For this project, I used data from 10X Genomics, which comes in the format of Matrix.mtx, genes.tsv, and barcodes.tsv files. I stored my data in a folder named hg19, but you can rename this if you’d like—just update the data.dir argument in the code to match your folder’s name.

To get started, place the data in a directory named data/ at the project root. 
Link: https://cf.10xgenomics.com/samples/cell/pbmc3k/pbmc3k_filtered_gene_bc_matrices.tar.gz

## Workflow Summary

1. Data Loading & Initialization

First, we load the scRNA-seq data and create a Seurat object. The Seurat object is a great structure because it holds all our data, from raw counts to metadata, and allows us to add on things like quality control (QC) metrics later.

2. Quality Control & Filtering

Next, I filter out cells that don’t meet quality standards. I removed cells with low counts of unique genes (probably low-quality cells) and cells with too high mitochondrial gene content (often stressed or dying cells).

3. Normalization & Identification of Variable Genes

Once we have clean data, it’s time to normalize it. I used a log-normalization method and identified the top 2,000 highly variable genes, which will help us see what makes each cell type unique.

4. Dimensionality Reduction & Clustering

With the variable genes identified, we move on to dimensionality reduction using PCA. Then, we use UMAP to visualize clusters based on their expression profiles. Each cluster represents a group of cells with similar gene expression, which often translates to distinct cell types.

5. Differential Expression & Marker Identification

After clustering, we find specific genes that are expressed differently in each cluster. These are known as “marker genes” and are critical for identifying what type of cells we’re looking at in each cluster.

6. Visualization & Results Export

Finally, we visualize our results using UMAP plots, violin plots, and feature plots, and save the final Seurat object so that we can come back to it anytime.


## Results

The analysis produces a few key visualizations:

	•	UMAP Plots: These 2D plots show clusters, which correspond to different cell types or states.
	•	Violin Plots: These display the expression levels of specific genes across clusters.
	•	Feature Plots: Highlight the expression of selected genes on UMAP plots to help with cell-type identification.
	•	Heatmaps: Display top marker genes for each cluster, giving an overview of the distinguishing features of each cell type.
