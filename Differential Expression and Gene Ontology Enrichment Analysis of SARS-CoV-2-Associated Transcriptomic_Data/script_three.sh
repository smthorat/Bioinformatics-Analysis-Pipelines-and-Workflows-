#!/bin/bash

#SBATCH -J Gene_quantification
#SBATCH -p general
#SBATCH -o filename_%j.txt
#SBATCH -e filename_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=smthorat@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=50:00:00
#SBATCH --mem=100GB
#SBATCH -A r00750

# Load required modules
module load conda
module load hisat
module load samtools
module load subread

# Set paths
SAM_DIR="/N/project/Krolab/Swaraj/project/alignment_output"
BAM_DIR="/N/project/Krolab/Swaraj/project/alignment_output/bam"
QUANT_DIR="/N/project/Krolab/Swaraj/project/alignment_output/quantification"
GTF_FILE="/N/project/Krolab/Swaraj/project/genome/GCF_000001405.40_GRCh38.p14_genomic.gtf"

# Create directories if they do not exist
mkdir -p "$BAM_DIR"
mkdir -p "$QUANT_DIR"

# Convert SAM to BAM, sort, and index
echo "Converting SAM files to sorted BAM files..."
for sam_file in "$SAM_DIR"/*.sam; do
    base=$(basename "$sam_file" .sam)
    echo "Processing $sam_file..."
    # Convert SAM to BAM
    samtools view -S -b "$sam_file" > "$BAM_DIR/${base}.bam"
    # Sort BAM
    samtools sort "$BAM_DIR/${base}.bam" -o "$BAM_DIR/${base}_sorted.bam"
    # Index BAM
    samtools index "$BAM_DIR/${base}_sorted.bam"
done

echo "SAM to BAM conversion and sorting completed."

# Perform gene quantification using featureCounts
echo "Running featureCounts for gene quantification..."
featureCounts -T 4 -a "$GTF_FILE" -o "$QUANT_DIR/counts.txt" "$BAM_DIR"/*_sorted.bam

echo "Gene quantification completed. Results are in '$QUANT_DIR/counts.txt'."

# Display summary of quantification results
echo "Preview of quantification output:"
head "$QUANT_DIR/counts.txt"

echo "All steps completed successfully!"