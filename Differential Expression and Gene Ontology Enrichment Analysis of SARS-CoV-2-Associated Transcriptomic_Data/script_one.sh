#!/bin/bash

# Convert to FASTQ (if needed)
for sra in $(cat run_list.txt); do
    fasterq-dump $sra
done