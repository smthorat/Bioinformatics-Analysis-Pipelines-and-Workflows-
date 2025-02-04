
# Create a working dir and go to the working directory
## The && operator helps execute two commands using a single line of code.
# mkdir af_xmpl_run && cd af_xmpl_run

# Fetch the example dataset and CB permit list and decompress them
## The pipe operator (|) passes the output of the wget command to the tar command.
## The dash operator (-) after `tar xzf` captures the output of the first command.
## - example dataset
# wget -qO- https://umd.box.com/shared/static/lx2xownlrhz3us8496tyu9c4dgade814.gz | tar xzf - --strip-components=1 -C .
## The fetched folder containing the fastq files is called toy_read_fastq.
fastq_dir="/N/project/Krolab/Swaraj/single_cell/af_xmpl_run/toy_read_fastq"
## The fetched folder containing the human ref files is called toy_human_ref.
ref_dir="/N/project/Krolab/Swaraj/single_cell/af_xmpl_run/toy_human_ref"

# # Fetch CB permit list
# ## the right chevron (>) redirects the STDOUT to a file.
# wget -qO- https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz | gunzip - > 3M-february-2018.txt

# ## Above steps involve downloading fastq and reference genome file ## 

# ### Now index the reference genome ###

# simpleaf needs the environment variable ALEVIN_FRY_HOME to store configuration and data.
# # For example, the paths to the underlying programs it uses and the CB permit list
mkdir alevin_fry_home 

export ALEVIN_FRY_HOME='alevin_fry_home'

# the simpleaf set-paths command finds the path to the required tools and writes a configuration JSON file in the ALEVIN_FRY_HOME folder.
simpleaf set-paths

# # simpleaf index
# # Usage: simpleaf index -o out_dir [-f genome_fasta -g gene_annotation_GTF|--refseq transcriptome_fasta] -r read_length -t number_of_threads
# ## The -r read_lengh is the number of sequencing cycles performed by the sequencer to generate biological reads (read2 in Illumina).
# ## Publicly available datasets usually have the read length in the description. Sometimes they are called the number of cycles.
# simpleaf index \
# -o simpleaf_index \
# -f /N/project/Krolab/Swaraj/single_cell/af_xmpl_run/toy_human_ref/fasta/genome.fa \
# -g /N/project/Krolab/Swaraj/single_cell/af_xmpl_run/toy_human_ref/genes/genes.gtf \
# -r 90 \
# -t 8

## After indexing we want to create gene count matrix

# Collecting sequencing read files
# ## The reads1 and reads2 variables are defined by finding the filenames with the pattern "_R1_" and "_R2_" from the toy_read_fastq directory.
# reads1_pat="selected_R1_reads.fastq"
# reads2_pat="selected_R2_reads.fastq"
# ## The read files must be sorted and separated by a comma.
# ### The find command finds the files in the fastq_dir with the name pattern
# ### The sort command sorts the file names
# ### The awk command and the paste command together convert the file names into a comma-separated string.
# reads1="$(find -L ${fastq_dir} -name "*${reads1_pat}*" -type f | sort | awk -v OFS=, '{$1=$1; print}' | paste -sd,)"
# reads2="$(find -L ${fastq_dir} -name "*${reads2_pat}*" -type f | sort | awk -v OFS=, '{$1=$1; print}' | paste -sd,)"
# # simpleaf quant
# ## Usage: simpleaf quant -c chemistry -t threads -1 reads1 -2 reads2 -i index -u [unspliced permit list] -r resolution -m t2g_3col -o output_dir
# # Now run simpleaf quant with the collected files
# simpleaf quant \
# -c 10xv3 -t 8 \
# -1 $reads1 -2 $reads2 \
# -i simpleaf_index/index \
# -u -r cr-like \
# -m simpleaf_index/index/t2g_3col.tsv \
# -o simpleaf_quant


# Each line in `quants_mat.mtx` represents
# a non-zero entry in the format row column entry
tail -3 simpleaf_quant/af_quant/alevin/quants_mat.mtx


# Each line in `quants_mat_cols.txt` is a splice status
# of a gene in the format (gene name)-(splice status)
tail -3 simpleaf_quant/af_quant/alevin/quants_mat_cols.txt


# Each line in `quants_mat_rows.txt` is a corrected
# (and, potentially, filtered) cell barcode
tail -3 simpleaf_quant/af_quant/alevin/quants_mat_rows.txt
