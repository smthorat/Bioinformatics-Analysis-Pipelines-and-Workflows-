#!/bin/bash

# Convert to FASTQ (if needed)
for sra in $(cat run_list.txt); do
    fasterq-dump $sra
done
[smthorat@h1 assignment_2]$ cat trimming.sh
#!/bin/bash
mkdir -p trimming

for file in *.fastq; do
    trim_galore --phred33 --fastqc "$file" -o trimming
done