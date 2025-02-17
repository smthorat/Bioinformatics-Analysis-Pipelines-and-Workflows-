# Nextflow Workflow Explanation: E. coli Variant Calling Pipeline

This Nextflow workflow is a structured and automated pipeline for detecting genetic variants in *E. coli* sequencing data. It ensures reproducibility, scalability, and efficiency by leveraging parallel execution and optimized data processing techniques. Below is a step-by-step breakdown of each stage of the pipeline.

---

## 1️⃣ Reference Genome Indexing

Before starting variant calling, the reference genome is prepared by generating index files that allow efficient and fast access to specific regions during mapping.

- **FASTA Indexing:**  
  Creates a small index file enabling tools like Samtools, BWA, and GATK to rapidly retrieve sections of the reference genome rather than scanning the entire file.

- **Sequence Dictionary:**  
  GATK requires a sequence dictionary file containing metadata about the reference genome (e.g., chromosome names and lengths).

*Without these index files, subsequent steps like read alignment and variant calling would be slow and inefficient.*

---

## 2️⃣ Quality Control of Raw Reads

Quality control is performed on the raw sequencing reads (FASTQ files) to ensure high-quality data before alignment. Poor-quality reads can lead to false positives in variant detection.

- **Quality Score Distribution:**  
  Evaluates the accuracy of base calls made by the sequencer.

- **GC Content Analysis:**  
  Checks for significant bias in the sequencing process.

- **Overrepresented Sequences:**  
  Identifies sequences that appear more frequently than expected, which might indicate contamination or sequencing artifacts.

*Performing quality control helps determine if trimming or other pre-processing steps are necessary before alignment.*

---

## 3️⃣ Read Alignment to Reference Genome

Sequencing reads are mapped to the reference genome using **BWA-MEM**, a widely used alignment algorithm optimized for next-generation sequencing data.

- **Mapping Reads:**  
  Aligns each short read to the most similar region in the reference genome.

- **Read Groups:**  
  Adds metadata (e.g., sample ID and platform) to each read, assisting in downstream analyses like duplicate marking and variant calling.

- **Handling Unmapped Reads:**  
  Some reads may not align due to sequencing errors or contamination and are typically discarded or analyzed separately.

*The output of this step is a SAM file (Sequence Alignment Map) containing alignment information in a human-readable format.*

---

## 4️⃣ Sorting, Deduplication, and BAM Processing

Aligned reads are processed to improve computational efficiency and accuracy before variant calling.

- **Converting SAM to BAM:**  
  SAM files are converted to BAM files, which are compressed, indexed, and more efficient to process.

- **Sorting Alignments:**  
  Alignments are sorted by genomic position to allow faster processing.

- **Removing PCR Duplicates:**  
  Duplicate reads from library preparation are removed to prevent bias, as they do not represent true biological variation.

*After this step, you obtain a clean, indexed BAM file ready for variant calling.*

---

## 5️⃣ Variant Calling

High-quality, aligned reads are used to identify genetic variants.

- **HaplotypeCaller (GVCF Mode):**  
  Detects single nucleotide polymorphisms (SNPs) and small insertions/deletions (INDELs) and outputs results in Genomic VCF (GVCF) format, which includes information on both variant and non-variant sites.

- **Genotyping Variants:**  
  The GVCF file is processed using GenotypeGVCFs to produce the final VCF file containing all detected variants.

*At this stage, you have an unfiltered variant file that includes both high-confidence and low-confidence calls.*

---

## 6️⃣ Filtering Variants for Quality

Not all detected variants are biologically relevant. Filtering is applied to remove potential sequencing errors, mapping artifacts, or variants from low-quality regions.

**Filters Applied:**

1. **Quality by Depth (QD):**  
   Low QD values indicate uncertain variant calls.

2. **Fisher Strand Bias (FS):**  
   Measures strand bias, which can signal sequencing errors.

3. **Mapping Quality (MQ):**  
   Ensures that variants are in well-mapped regions.

4. **Strand Odds Ratio (SOR):**  
   Detects imbalanced strand representation.

5. **Genotype Quality (GQ) and Read Depth (DP):**  
   Ensures variants have sufficient supporting reads.

*After filtering, the pipeline produces high-confidence SNPs and INDELs by removing false positives.*

---

## 7️⃣ Variant Annotation

With a filtered, high-quality set of variants, the next step is to understand their biological significance.

**During Annotation:**

- Each SNP and INDEL is compared against known databases.
- Functional impact on genes, proteins, and regulatory regions is assessed.
- Information such as amino acid changes, gene names, and mutation effects is added.

*Note:* Since Funcotator (GATK’s annotation tool) is designed for human genomes, for *E. coli*, tools like **SnpEff** or **Prokka** may be more appropriate for annotation.

---

## Final Outputs and Interpretation

At the end of the pipeline, the following key files are generated:

1. **Final Variant Call File (`final_variants.vcf.gz`):**  
   Contains all detected genetic variants in *E. coli*.

2. **Filtered High-Quality SNPs (`analysis-ready-snps.vcf`):**  
   A refined list of high-confidence single nucleotide polymorphisms.

3. **Filtered High-Quality INDELs (`analysis-ready-indels.vcf`):**  
   High-confidence insertions and deletions.

4. **Annotated Variants (`analysis-ready-snps-filteredGT-functotated.vcf`):**  
   Variants with gene names, functional impact, and additional metadata.

*These results are useful for comparative genomics, antibiotic resistance studies, and evolutionary analysis in *E. coli**.

---

## Why Use Nextflow for This Pipeline?

- **Scalability:**  
  Can handle multiple samples in parallel.

- **Reproducibility:**  
  Ensures consistent results regardless of the execution environment.

- **Portability:**  
  Works on local systems, HPC clusters, or cloud environments.

- **Checkpointing:**  
  Resumes processing from the point of interruption, reducing downtime.

*This workflow can be easily extended to include additional analyses such as structural variant calling, phylogenetic analysis, or antimicrobial resistance (AMR) detection based on research goals.*

---

This Nextflow workflow provides a robust, scalable, and reproducible pipeline for *E. coli* variant calling, ensuring high-quality results that are ready for downstream analysis and interpretation.
