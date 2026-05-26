#!/bin/bash

############################################################
# RNA-seq preprocessing pipeline
############################################################

# Quality control
fastqc *.fastq.gz

# Adapter trimming
trimmomatic PE \
GSE_cortex_R1.fastq.gz GSE_cortex_R2.fastq.gz \
GSE_cortex_R1_trimmed.fastq.gz GSE_cortex_R1_unpaired.fastq.gz \
GSE_cortex_R2_trimmed.fastq.gz GSE_cortex_R2_unpaired.fastq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# STAR alignment
STAR \
--runThreadN 8 \
--genomeDir GRCm39_STARindex \
--readFilesIn GSE_cortex_R1_trimmed.fastq.gz GSE_cortex_R2_trimmed.fastq.gz \
--readFilesCommand zcat \
--outSAMtype BAM SortedByCoordinate

# featureCounts
featureCounts \
-a gencode.vM32.annotation.gtf \
-o counts.txt \
*.bam
