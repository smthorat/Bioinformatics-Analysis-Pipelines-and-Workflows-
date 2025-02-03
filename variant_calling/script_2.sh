#!/bin/bash

#SBATCH -J Filter_Annotate
#SBATCH -p general
#SBATCH -o logs/filter_annotate_%j.out
#SBATCH -e logs/filter_annotate_%j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=smthorat@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=48:00:00
#SBATCH --mem=50GB
#SBATCH -A r00750

# --------------------------------------
# 1. Load Required Modules
# --------------------------------------
module purge
module load java
module load gatk
module load bcftools

# --------------------------------------
# 2. Set Paths and Create Directories
# --------------------------------------
base_dir="/N/project/Krolab/Swaraj/Ecoli_variant_calling"
ref="$base_dir/ref/ecoli_ref.fasta"
results="$base_dir/results"
logs="$base_dir/logs"

mkdir -p "$logs"

echo "Starting variant filtering and annotation pipeline..."

# --------------------------------------
# 3. Filter Variants
# --------------------------------------

echo "Filtering SNPs..."
gatk VariantFiltration \
	-R "$ref" \
	-V "$results/final_variants.vcf" \
	-O "$results/filtered_snps.vcf" \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter" || { echo "ERROR: SNP filtering failed"; exit 1; }

echo "Filtering INDELs..."
gatk VariantFiltration \
	-R "$ref" \
	-V "$results/final_variants.vcf" \
	-O "$results/filtered_indels.vcf" \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter" || { echo "ERROR: INDEL filtering failed"; exit 1; }

# --------------------------------------
# 4. Select High-Quality Variants
# --------------------------------------

echo "Selecting high-quality SNPs..."
gatk SelectVariants \
	--exclude-filtered \
	-V "$results/filtered_snps.vcf" \
	-O "$results/analysis-ready-snps.vcf" || { echo "ERROR: Selecting SNPs failed"; exit 1; }

echo "Selecting high-quality INDELs..."
gatk SelectVariants \
	--exclude-filtered \
	-V "$results/filtered_indels.vcf" \
	-O "$results/analysis-ready-indels.vcf" || { echo "ERROR: Selecting INDELs failed"; exit 1; }

# Remove low-quality genotype-filtered variants
echo "Removing low-quality genotype-filtered SNPs..."
grep -v -E "DP_filter|GQ_filter" "$results/analysis-ready-snps.vcf" > "$results/analysis-ready-snps-filteredGT.vcf"

echo "Removing low-quality genotype-filtered INDELs..."
grep -v -E "DP_filter|GQ_filter" "$results/analysis-ready-indels.vcf" > "$results/analysis-ready-indels-filteredGT.vcf"

# --------------------------------------
# 5. Annotate Variants
# --------------------------------------
# NOTE: Funcotator is designed for human genomes, but we assume you will replace this with a bacterial annotation tool if needed.

echo "Annotating SNPs with Funcotator..."
gatk Funcotator \
	--variant "$results/analysis-ready-snps-filteredGT.vcf" \
	--reference "$ref" \
	--ref-version ecoli \
	--data-sources-path "$base_dir/funcotator_data" \
	--output "$results/analysis-ready-snps-filteredGT-functotated.vcf" \
	--output-file-format VCF || { echo "ERROR: SNP annotation failed"; exit 1; }

echo "Annotating INDELs with Funcotator..."
gatk Funcotator \
	--variant "$results/analysis-ready-indels-filteredGT.vcf" \
	--reference "$ref" \
	--ref-version ecoli \
	--data-sources-path "$base_dir/funcotator_data" \
	--output "$results/analysis-ready-indels-filteredGT-functotated.vcf" \
	--output-file-format VCF || { echo "ERROR: INDEL annotation failed"; exit 1; }

# --------------------------------------
# 6. Extract VCF Information to a Table
# --------------------------------------

echo "Extracting SNP data to table..."
gatk VariantsToTable \
	-V "$results/analysis-ready-snps-filteredGT-functotated.vcf" \
	-F AC -F AN -F DP -F AF -F FUNCOTATION \
	-O "$results/output_snps.table" || { echo "ERROR: Variant-to-table conversion failed"; exit 1; }

echo "Pipeline completed successfully."