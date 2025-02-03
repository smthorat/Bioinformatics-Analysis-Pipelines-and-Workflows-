#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Define Parameters
params.base_dir = "/N/project/Krolab/Swaraj/Ecoli_variant_calling"
params.ref = "${params.base_dir}/ref/ecoli_ref.fasta"
params.reads_dir = "${params.base_dir}/reads"
params.aligned_dir = "${params.base_dir}/aligned_reads"
params.results_dir = "${params.base_dir}/results"
params.logs_dir = "${params.base_dir}/logs"

process index_reference {
    tag "Index Reference"

    input:
    path ref from params.ref

    output:
    path "${params.ref}.fai"
    path "${params.ref.replace('.fasta', '.dict')}"

    script:
    """
    samtools faidx $ref
    gatk CreateSequenceDictionary -R $ref -O ${params.ref.replace('.fasta', '.dict')}
    """
}

process fastqc {
    tag "Quality Control"

    input:
    path reads from file("${params.reads_dir}/SRR519926_*.fastq.gz")

    output:
    path "qc_reports"

    script:
    """
    mkdir -p ${params.results_dir}/qc_reports
    fastqc -o ${params.results_dir}/qc_reports $reads
    """
}

process align_reads {
    tag "Align Reads"

    input:
    path ref from params.ref
    path r1 from file("${params.reads_dir}/SRR519926_1.fastq.gz")
    path r2 from file("${params.reads_dir}/SRR519926_2.fastq.gz")

    output:
    path "${params.aligned_dir}/SRR519926.sorted.dedup.bam"

    script:
    """
    bwa mem -t 4 -R "@RG\\tID:SRR519926\\tPL:ILLUMINA\\tSM:SRR519926" $ref $r1 $r2 > ${params.aligned_dir}/SRR519926.sam
    samtools view -Sb ${params.aligned_dir}/SRR519926.sam | samtools sort -o ${params.aligned_dir}/SRR519926.sorted.bam
    rm ${params.aligned_dir}/SRR519926.sam
    samtools index ${params.aligned_dir}/SRR519926.sorted.bam
    gatk MarkDuplicates -I ${params.aligned_dir}/SRR519926.sorted.bam -O ${params.aligned_dir}/SRR519926.sorted.dedup.bam -M ${params.aligned_dir}/marked_dup_metrics.txt
    samtools index ${params.aligned_dir}/SRR519926.sorted.dedup.bam
    """
}

process call_variants {
    tag "Call Variants"

    input:
    path ref from params.ref
    path bam from file("${params.aligned_dir}/SRR519926.sorted.dedup.bam")

    output:
    path "${params.results_dir}/final_variants.vcf.gz"

    script:
    """
    gatk HaplotypeCaller -R $ref -I $bam -O ${params.results_dir}/raw_variants.g.vcf -ERC GVCF -ploidy 1
    gatk GenotypeGVCFs -R $ref -V ${params.results_dir}/raw_variants.g.vcf -O ${params.results_dir}/final_variants.vcf
    bgzip -c ${params.results_dir}/final_variants.vcf > ${params.results_dir}/final_variants.vcf.gz
    bcftools index ${params.results_dir}/final_variants.vcf.gz
    """
}

process filter_variants {
    tag "Filter Variants"

    input:
    path ref from params.ref
    path vcf from file("${params.results_dir}/final_variants.vcf.gz")

    output:
    path "${params.results_dir}/analysis-ready-snps.vcf"
    path "${params.results_dir}/analysis-ready-indels.vcf"

    script:
    """
    gatk VariantFiltration -R $ref -V $vcf -O ${params.results_dir}/filtered_snps.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 60.0" \
        -filter-name "MQ_filter" -filter "MQ < 40.0" \
        -filter-name "SOR_filter" -filter "SOR > 4.0" \
        -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"

    gatk VariantFiltration -R $ref -V $vcf -O ${params.results_dir}/filtered_indels.vcf \
        -filter-name "QD_filter" -filter "QD < 2.0" \
        -filter-name "FS_filter" -filter "FS > 200.0" \
        -filter-name "SOR_filter" -filter "SOR > 10.0" \
        -genotype-filter-expression "DP < 10" -genotype-filter-name "DP_filter" \
        -genotype-filter-expression "GQ < 10" -genotype-filter-name "GQ_filter"

    gatk SelectVariants --exclude-filtered -V ${params.results_dir}/filtered_snps.vcf -O ${params.results_dir}/analysis-ready-snps.vcf
    gatk SelectVariants --exclude-filtered -V ${params.results_dir}/filtered_indels.vcf -O ${params.results_dir}/analysis-ready-indels.vcf
    """
}

process annotate_variants {
    tag "Annotate Variants"

    input:
    path ref from params.ref
    path snp_vcf from file("${params.results_dir}/analysis-ready-snps.vcf")
    path indel_vcf from file("${params.results_dir}/analysis-ready-indels.vcf")

    output:
    path "${params.results_dir}/analysis-ready-snps-filteredGT-functotated.vcf"
    path "${params.results_dir}/analysis-ready-indels-filteredGT-functotated.vcf"

    script:
    """
    gatk Funcotator --variant $snp_vcf --reference $ref --ref-version ecoli \
        --data-sources-path ${params.base_dir}/funcotator_data \
        --output ${params.results_dir}/analysis-ready-snps-filteredGT-functotated.vcf --output-file-format VCF

    gatk Funcotator --variant $indel_vcf --reference $ref --ref-version ecoli \
        --data-sources-path ${params.base_dir}/funcotator_data \
        --output ${params.results_dir}/analysis-ready-indels-filteredGT-functotated.vcf --output-file-format VCF
    """
}

workflow {
    index_reference()
    fastqc()
    align_reads()
    call_variants()
    filter_variants()
    annotate_variants()
}