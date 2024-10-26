{\rtf1\ansi\ansicpg1252\cocoartf2818
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\froman\fcharset0 Times-Roman;}
{\colortbl;\red255\green255\blue255;\red255\green255\blue255;}
{\*\expandedcolortbl;;\cssrgb\c100000\c100000\c100000;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\partightenfactor0

\f0\fs32 \cf0 \expnd0\expndtw0\kerning0
# Gene Expression Analysis Project\
\
## Overview\
\
This project analyzes gene expression data to identify enriched biological pathways under different conditions in **LNCAP** and **PC3** cell lines. The analysis includes differential expression and pathway enrichment using the **DESeq2** and **clusterProfiler** packages in R.\
\
## Objectives\
\
- Analyze gene expression data from **LNCAP** and **PC3** cell lines.\
- Identify differentially expressed genes (DEGs) between **hypoxia** and **normoxia** conditions.\
- Perform pathway enrichment analysis to find significant biological processes.\
\
## Steps\
\
### 1. Data Preparation\
\
- **Load Data**: Import raw count data from CSV files.\
- **Prepare Data**: Organize the data by setting row names and sorting columns.\
\
### 2. Differential Expression Analysis\
\
- **Create DESeq2 Object**: Use `DESeqDataSetFromMatrix` to create a DESeq2 object with count data and condition information.\
- **Run DESeq2**: Perform differential expression analysis to identify DEGs.\
\
### 3. Normalization\
\
- **Normalize Counts**: Obtain normalized counts using the `counts` function with `normalized = TRUE`.\
\
### 4. Pathway Enrichment Analysis\
\
- **Perform Enrichment**: Use `clusterProfiler` to conduct GO enrichment analysis for biological processes (BP).\
- **Visualize Results**: Create dot plots to visualize enriched pathways for each condition comparison.\
\
### 5. Visualization\
\
- **Volcano Plots**: Generate volcano plots to visualize DEGs.\
- **Pathway Plots**: Use dot plots to display top enriched pathways, focusing on gene names for clarity.\
\
## Key Functions and Packages\
\
- **DESeq2**: For differential expression analysis.\
- **clusterProfiler**: For GO enrichment analysis.\
- **org.Hs.eg.db**: Annotation database for human genes.\
- **ggplot2**: For creating visualizations.\
- **gridExtra**: To arrange multiple plots in a grid layout.\
- **tidyverse**: For data manipulation and visualization.\
\
## Requirements\
\
- R version 4.0 or higher\
- R packages: DESeq2, clusterProfiler, org.Hs.eg.db, ggplot2, gridExtra, tidyverse\
\
## Files\
\
- `raw_counts.csv`: Contains raw gene expression counts.\
- `GRCh38.p13_annotation.csv`: Annotation file for gene IDs.\
- `*_GO_enrichment.csv`: Output files for GO enrichment results.\
\
## Visualization Example\
\
The attached image shows dot plots of enriched pathways for different conditions, highlighting the most significant biological processes.\
\
## Conclusion\
\
This project provides insights into the biological processes affected by hypoxia and normoxia in LNCAP and PC3 cell lines, aiding in understanding the underlying mechanisms of these conditions.\
\
Feel free to explore the code and modify it as needed for your specific research questions!}