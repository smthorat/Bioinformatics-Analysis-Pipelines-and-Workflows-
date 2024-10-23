{\rtf1\ansi\ansicpg1252\cocoartf2818
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/bin/bash\
# RNA-Seq Pipeline Script for Bacillus subtilis\
# Author: Swaraj\
# Date: September 2024\
\
# -------------------------------------------------\
# Pre-requisites: Downloading Data and Required Files\
# -------------------------------------------------\
\
# Load required modules\
# Ensure the following tools are available on your HPC and replace 'module_name' with the appropriate module name on your system\
\
module load sra-toolkit    # For downloading data from SRA\
module load fastqc         # For quality control\
module load trimmomatic    # For trimming reads\
module load hisat2         # For alignment\
module load samtools       # For BAM/SAM file handling\
module load subread        # For featureCounts (quantification)\
\
# Step 0: Download the FASTQ files (using the SRA toolkit)\
# Download FASTQ files using SRA Toolkit (you can replace SRR30615974 with your sample ID)\
echo "Downloading FASTQ data from SRA..."\
fasterq-dump SRR30615974\
\
# Required files:\
# - Genome fasta file: GCF_000009045.1_ASM904v1_genomic.fna\
# - Annotation file (GTF): genomic.gtf\
# - Paired-end FASTQ files: SRR30615974_1.fastq, SRR30615974_2.fastq\
# - Trimmomatic adapter file: TruSeq3-PE-2.fa\
\
# -------------------------------------------------\
# RNA-Seq Pipeline\
# -------------------------------------------------\
\
# Step 1: Quality Control with FastQC\
# Perform quality control on the raw FASTQ files\
# Input: Paired-end FASTQ files\
\
echo "Running FastQC on paired-end FASTQ files..."\
fastqc SRR30615974_1.fastq SRR30615974_2.fastq -o ./fastqc_output/\
\
# Step 2: Trimming Reads with Trimmomatic\
# Trim low-quality reads and adapter sequences\
\
echo "Trimming reads with Trimmomatic..."\
java -jar trimmomatic-0.39.jar PE -phred33 \\\
SRR30615974_1.fastq SRR30615974_2.fastq \\\
SRR30615974_1_paired_trimmed.fastq SRR30615974_1_unpaired.fastq \\\
SRR30615974_2_paired_trimmed.fastq SRR30615974_2_unpaired.fastq \\\
ILLUMINACLIP:TruSeq3-PE-2.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36\
\
# Step 3: Build HISAT2 Index for the Reference Genome\
# Input: Reference genome (FASTA)\
\
echo "Building HISAT2 index for the reference genome..."\
hisat2-build GCF_000009045.1_ASM904v1_genomic.fna genome_index\
\
# Step 4: Align Reads to the Genome Using HISAT2\
# Input: Paired-end trimmed FASTQ files and genome index\
# Condensed into a one-liner\
\
hisat2 -p 8 -x genome_index -1 SRR30615974_1_paired_trimmed.fastq -2 SRR30615974_2_paired_trimmed.fastq -S SRR30615974_aligned.sam\
\
# Step 5: Convert SAM to BAM and Sort BAM\
# Convert the SAM file to BAM format and sort the BAM file using Samtools\
\
echo "Converting SAM to BAM and sorting the BAM file..."\
samtools view -bS SRR30615974_aligned.sam > SRR30615974_aligned.bam\
samtools sort SRR30615974_aligned.bam -o SRR30615974_sorted.bam\
\
# Step 6: Index the Sorted BAM File\
# Indexing the BAM file to prepare for counting\
\
echo "Indexing the sorted BAM file..."\
samtools index SRR30615974_sorted.bam\
\
# Step 7: Quantify Gene Expression with featureCounts\
# Input: Sorted BAM file and GTF annotation file\
\
echo "Quantifying gene expression using featureCounts..."\
featureCounts -T 8 -p -t exon -g gene_id -a genomic.gtf -o SRR30615974_counts.txt SRR30615974_sorted.bam\
\
# Step 8: Transfer File from HPC to Local Machine Using SCP\
# Command to transfer the output file from the HPC to a local machine\
\
echo "Transferring the counts file to the local machine..."\
scp smthorat@quartz.uits.iu.edu:/N/u/smthorat/Quartz/Swaraj/RNA_pipeline/SRR30615974_counts.txt /Users/swaraj/Downloads/\
\
# End of Pipeline\
echo "RNA-seq pipeline completed successfully!"}