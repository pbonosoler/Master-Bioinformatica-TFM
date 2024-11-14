#!/bin/bash

#index file for fasta (built with bowtie2-build
IND=/mnt/Drive-D/DATA/Scer_reference_files/S288C_R64-2-1_20150113_renamed.fsa
GFF=/mnt/Drive-D/DATA/Scer_reference_files/R64-2-1_BSMedited2noFasta.gff

cd /mnt/Drive-E/seq/Sc_glycerol_t200

for i in quality_files2/*_.fastq.gz
do
# this removes the path and the end of the file name to keep only the sample "L*"
a=${i#quality_files2\/}
sample=${a%_.fastq.gz}

fq=${i%.R1.fastq.gz}

echo $sample

if [ ! -e counts/$sample\_counts.table ]; then
echo "run pipeline on counts/$sample\_counts.table"

#mq_AS remove reads with quality score < 2 and filter reads with more than one alignment. Secondary alignment score > primary alignment score are removed.

samtools view -h aln2/$sample.sorted.byname.bam | perl /home/christina/bin/MQ_AS_check.pl - aln2/$sample-mq_as.sam
rm aln2/$sample.bam
rm aln2/$sample.sam
rm aln2/$sample.sorted.byname.bam
samtools view -h -b -S -o aln2/$sample.mq_as.bam aln2/$sample-mq_as.sam
rm aln2/$sample-mq_as.sam


#htseq-count and assign counts to each strain
#samtools view $sample.mq_as.bam | htseq-count -r name -a 1 --idattr ID -t CDS -m union --stranded=no - $GFF > $sample\_counts.table
htseq-count -f bam -r name -a 1 --idattr Parent -t mRNA -m union --stranded=no aln2/$sample.mq_as.bam $GFF > counts/$sample\_counts.table

fi
done

