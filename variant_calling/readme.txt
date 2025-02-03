# E. coli Variant Calling Pipeline


A robust bioinformatics pipeline for calling variants in Escherichia coli sequencing data using Illumina paired-end reads.

## Features

- **End-to-end workflow** from raw reads to final variants
- Quality control with FastQC
- BWA-based alignment to reference genome
- GATK-based variant calling following best practices
- SLURM cluster compatibility
- Comprehensive error checking and validation
- Automated directory structure setup
- Detailed QC reports and statistics

## Prerequisites

### Software Requirements
- BWA (v0.7.17)
- GATK (v4.4.0.0)
- samtools (v1.20)
- bcftools (v1.17)
- FastQC (v0.12.1)
- Java (v11+)

### Hardware Requirements
Basic, I have used Macbook Pro M2 

## Directory Structure
Ecoli_variant_calling/
├── align_reads/ # Alignment intermediate files
├── aligned_reads/ # Processed BAM files
├── data/ # Raw data (optional)
├── reads/ # Input FASTQ files
├── ref/ # Reference genome files
│ └── ecoli_ref.fasta # Reference genome
└── results/ # Final outputs and reports
├── qc_reports/ # FastQC reports
└── intermediate_files/ # Processing intermediates

