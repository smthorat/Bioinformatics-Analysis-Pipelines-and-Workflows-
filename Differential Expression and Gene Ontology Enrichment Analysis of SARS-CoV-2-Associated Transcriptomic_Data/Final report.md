# Assignment 2

## Objective

This report outlines the process for analyzing miRNA transcriptomic data to identify differential gene expression between mock-control and SARS-CoV-2-infected samples at two time points: 24 hours and 72 hours post-infection. The analysis was performed using RNA sequencing data from the NCBI BioProject PRJNA901149.

## Data Source

* BioProject ID: PRJNA901149
* Reference Genome: GRCh38.p14
* Files Used:
* Reference genome FASTA: GCF_000001405.40_GRCh38.p14_genomic.fna
* Annotation file GTF: GCF_000001405.40_GRCh38.p14_genomic.gtf

## Tools and Software

* Data Preprocessing: SRA Toolkit, FastQC, Trim Galore
* Alignment: HISAT2
* File Conversion: Samtools
* Quantification: Subread (featureCounts)
* Differential Expression Analysis: R (DESeq2)
* Visualization: ggplot2, ggrepel, R base plotting functions
* Gene Ontology Enrichment: clusterProfiler, org.Hs.eg.db

## Steps in the Analysis

### 1.Data Acquisition

Tools Used: SRA Toolkit, Python
I downloaded the raw sequencing data from the NCBI BioProject repository using the SRA Toolkit. The sample accession numbers were extracted into a text file to enable bulk downloading. A Python script with a loop was used to automate the downloading process, ensuring that all 12 samples were obtained efficiently.

## 2.Reference Genome and Annotation Files

Files Used:

	1.Genome Assembly: GRCh38.p14

	2.Reference Genome File: GCF_000001405.40_GRCh38.p14_genomic.fna

	3.Annotation File (GTF): GCF_000001405.40_GRCh38.p14_genomic.gtf


The human genome reference sequence and its  GTF annotation file were downloaded. These files are essential for alignment and quantification steps, as they provide a reference for mapping the sequencing reads and identifying transcript regions.

## 3.Quality Control of Raw Data

Tools Used: FastQC, Trim Galore
The raw sequencing files were assessed for quality using FastQC. This tool generates reports on read quality, GC content, adapter contamination, etc.

Following quality checks, Trim Galore was used for trimming low-quality bases and removing adapter sequences. This step ensured that the data were clean and suitable for downstream analyses.

## 4.Alignment to the Reference Genome

Tools Used: HISAT2
The trimmed sequencing reads were aligned to the GRCh38.p14 human reference genome using HISAT2.

Indexing: Before alignment, the reference genome was indexed with HISAT2. Indexing creates a data structure that allows for rapid alignment of sequencing reads.

Output: The output from this step was in SAM format.

## 5.Conversion to BAM Format

Tools Used: Samtools
The SAM files generated from HISAT2 were converted to BAM format using Samtools.

## 6.Gene Quantification

Tools Used: Subread
Gene expression levels were quantified using the featureCounts function from the Subread package. This tool calculates the number of reads mapped to each gene based on the annotation file.

## Differential Gene Expression Analysis

Identify differentially expressed genes between SARS-CoV-2-infected and mock-control samples and between two time points (24 hours and 72 hours). Additionally, gene ontology enrichment performed to understand the biological processes involved.

Steps in the Analysis

1.  Data Import
The count data were imported and prepared, with gene IDs as row names and sample metadata loaded into a separate metadata file for annotation.

2.	Data Preprocessing
Low-count genes were filtered out to improve data quality.

3.	Creation of DESeq Dataset
A DESeq dataset was created with filtered counts and metadata using the specified experimental design (Condition + Time Point).

4.	Regularized Log Transformation
The count data were transformed using regularized log for better visualization and downstream PCA analysis.

5.	Differential Expression Analysis
Differential expression analysis was performed for two comparisons: 	   	       SARS-CoV-2-infected versus mock-control and 72 hours versus 24 hours.

6.	Significant Gene Extraction
Significant genes were identified using an adjusted p-value threshold (padj < 0.05) for both conditions and time points.

7.	Visualization
MA plots and volcano plots were generated to visualize the differential expression results for both comparisons. PCA plots were created to observe sample clustering.

8.	Gene Ontology (GO) Enrichment Analysis
GO enrichment performed for significant genes to identify enriched biological processes. Separate analyses done for condition and time point comparisons.

9.	Visualization of GO Terms
Bar plots and dot plots were generated to visualize the enriched GO terms.

10.	Saving Results
Normalized count data, significant genes, and GO enrichment results were saved as CSV files for further interpretation.


Conclusion

This analysis successfully identified differentially expressed genes and enriched biological processes associated with SARS-CoV-2 infection and varying time points, providing insights into host immune responses and viral pathology. However, no significant GO terms were found for the SARS-CoV-2 infection condition, indicating a potential need for further exploration of the data