#!/bin/bash

# Directory containing the fastq files
fastq_dir="/mnt/Drive-E/seq/Sc_glycerol_t200/quality_files"
output_dir="/mnt/Drive-E/seq/Sc_glycerol_t200/quality_files/output_fastq"

# Run FastQC on all fastq.gz files
for file in "$fastq_dir"/*.fastq.gz
do
    fastqc "$file" -o "$output_dir"
done

# Run MultiQC on the FastQC output
multiqc "$output_dir" -o "$output_dir"

