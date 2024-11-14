#!/bin/bash

# Set variables
IND="/mnt/Drive-D/DATA/Scer_reference_files/S288C_R64-2-1_20150113_renamed.fsa"  # Path to the indexed reference genome
OUTPUT_DIR="/mnt/Drive-E/seq/Sc_glycerol_t200/aln"  # Output directory for alignment files

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Change directory to where the FastQ files are located
cd /mnt/Drive-E/seq/Sc_glycerol_t200

# Loop through all R1 files in the quality_files directory
for R1_file in quality_files/*.trim.R1.fastq.gz; do
    # Extract sample name from R1 file name
    a=${R1_file#quality_files\/}
    sample=${a%.trim.R1.fastq.gz}

    # Construct R2 file name based on sample name
    R2_file="quality_files/${sample}.trim.R2.fastq.gz"

    # Perform alignment with bowtie2
    bowtie2 -p 30 -x "$IND" -1 "$R1_file" -2 "$R2_file" -N 1 -q --local -S "${OUTPUT_DIR}/${sample}.sam" > "${OUTPUT_DIR}/${sample}.bawtie.info"

    # Convert SAM to BAM and sort
    samtools view -b -S -o "${OUTPUT_DIR}/${sample}.bam" "${OUTPUT_DIR}/${sample}.sam"
    samtools sort -n -o "${OUTPUT_DIR}/${sample}.sorted.byname.bam" "${OUTPUT_DIR}/${sample}.bam"

    # Optional: Remove intermediate files
    # rm "${OUTPUT_DIR}/${sample}.sam"
    # rm "${OUTPUT_DIR}/${sample}.bam"
done


