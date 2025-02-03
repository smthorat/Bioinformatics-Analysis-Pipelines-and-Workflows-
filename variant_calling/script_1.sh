#!/bin/bash

#SBATCH -J Variant_calling
#SBATCH -p general
#SBATCH -o logs/variant_calling_%j.out  # Store logs in a separate directory
#SBATCH -e logs/variant_calling_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=smthorat@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=95:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750


# E. coli Variant Calling Pipeline
# 
# A complete workflow for calling variants in E. coli sequencing data
# Includes quality control, alignment, duplicate marking, and variant calling
#
# Input Requirements:
# - Paired-end FASTQ files in reads/ directory (SRR519926_1.fastq.gz, SRR519926_2.fastq.gz)
# - Reference genome in ref/ directory (ecoli_ref.fasta)
#
# Output:
# - Final VCF with variants in results/ directory
# - QC reports in results/qc_reports
#
# Dependencies: GATK4, BWA, samtools, FastQC, bcftools


# --------------------------------------
# 1. Load Required Modules
# --------------------------------------
module purge
module load java            # Ensure Java 8/11 for GATK
module load fastqc
module load sra-toolkit
module load samtools        # Always specify versions when possible
module load bwa
module load gatk
module load bcftools

# --------------------------------------
# 2. Set Paths and Create Directories
# --------------------------------------
base_dir="/N/project/Krolab/Swaraj/Ecoli_variant_calling"
ref="$base_dir/ref/ecoli_ref.fasta"
reads="$base_dir/reads"
aligned="$base_dir/aligned_reads"
results="$base_dir/results"
logs="$base_dir/logs"

# Ensure required directories exist
mkdir -p "$aligned" "$results/qc_reports" "$logs"

# --------------------------------------
# 3. Prepare Reference Genome
# --------------------------------------
# Index reference genome (GATK requires both .fai and .dict files)
if [ ! -f "$ref.fai" ]; then
  echo "Indexing reference genome..."
  samtools faidx "$ref" || { echo "ERROR: Failed to index reference genome"; exit 1; }
fi

if [ ! -f "${ref%.fasta}.dict" ]; then
  echo "Creating sequence dictionary..."
  gatk CreateSequenceDictionary -R "$ref" -O "${ref%.fasta}.dict" || { echo "ERROR: Failed to create dictionary"; exit 1; }
fi

# --------------------------------------
# 4. Quality Control (FastQC on raw reads)
# --------------------------------------
fastqc -o "$results/qc_reports" "$reads/SRR519926_1.fastq.gz" "$reads/SRR519926_2.fastq.gz"

# --------------------------------------
# 5. Align Reads with BWA
# --------------------------------------
echo "Aligning reads with BWA..."
bwa mem -t 4 \
  -R "@RG\tID:SRR519926\tPL:ILLUMINA\tSM:SRR519926" \
  "$ref" \
  "$reads/SRR519926_1.fastq.gz" "$reads/SRR519926_2.fastq.gz" \
  > "$aligned/SRR519926.sam"

# Verify SAM file was created
if [ ! -s "$aligned/SRR519926.sam" ]; then
  echo "ERROR: Alignment failed. SAM file is empty."
  exit 1
fi

# --------------------------------------
# 6. Process BAM File
# --------------------------------------
echo "Converting and sorting BAM file..."
samtools view -@4 -Sb "$aligned/SRR519926.sam" | samtools sort -@4 -o "$aligned/SRR519926.sorted.bam"

# Index BAM file
samtools index "$aligned/SRR519926.sorted.bam"

# Remove intermediate SAM file to save space
rm "$aligned/SRR519926.sam"

# Mark duplicates
echo "Marking duplicates..."
gatk MarkDuplicates \
  -I "$aligned/SRR519926.sorted.bam" \
  -O "$aligned/SRR519926.sorted.dedup.bam" \
  -M "$aligned/marked_dup_metrics.txt" || { echo "ERROR: MarkDuplicates failed"; exit 1; }

# Index deduplicated BAM
samtools index "$aligned/SRR519926.sorted.dedup.bam"

# --------------------------------------
# 7. Variant Calling
# --------------------------------------
# Step 7a: HaplotypeCaller (generate GVCF)
echo "Running GATK HaplotypeCaller..."
gatk HaplotypeCaller \
  -R "$ref" \
  -I "$aligned/SRR519926.sorted.dedup.bam" \
  -O "$results/raw_variants.g.vcf" \
  -ERC GVCF \
  -ploidy 1 || { echo "ERROR: HaplotypeCaller failed"; exit 1; }

# Step 7b: GenotypeGVCFs (convert GVCF to VCF)
echo "Running GATK GenotypeGVCFs..."
gatk GenotypeGVCFs \
  -R "$ref" \
  -V "$results/raw_variants.g.vcf" \
  -O "$results/final_variants.vcf" || { echo "ERROR: GenotypeGVCFs failed"; exit 1; }

# --------------------------------------
# 8. Validate Output
# --------------------------------------
if [ -s "$results/final_variants.vcf" ]; then
  echo "Success! VCF file generated."
else
  echo "ERROR: VCF file is empty. Check logs."
  exit 1
fi

# Compress and index VCF
echo "Compressing and indexing VCF..."
bgzip -c "$results/final_variants.vcf" > "$results/final_variants.vcf.gz"
bcftools index "$results/final_variants.vcf.gz"

echo "Pipeline completed successfully."