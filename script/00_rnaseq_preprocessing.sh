#!/bin/bash

############################################################
# MULTI-DATASET RNA-SEQ RAW UPSTREAM PREPROCESSING PIPELINE
# Script: script/01_rnaseq_preprocessing.sh
############################################################

# Exit immediately if a command exits with a non-zero status
set -e

# Establish repository path configurations
DATA_DIR="data/RNAseq_raw"
RAW_COUNTS_DIR="data/RNAseq_counts"
LOG_DIR="results/logs"

mkdir -p "$RAW_COUNTS_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "figures"

# Define the explicit RNA-seq dataset list mapping directly to your manuscript's Table 3
datasets=("GSE132074" "GSE144957" "GSE222410" "GSE113754" "GSE211088" "GSE211301" "GSE156925" "GSE166831")

echo "======================================================================"
# Step 1: Broad Quality Control Evaluation via FastQC
echo "Starting global FastQC baseline assessment..."
fastqc "$DATA_DIR"/*.fastq.gz -o figures/ 2>> "$LOG_DIR"/fastqc_error.log

# Step 2: Loop through each study to perform automated Trimming, Alignment, and Quantitation
for ds in "${datasets[@]}"; do
    echo "----------------------------------------------------------------------"
    echo "Processing Dataset: $ds"
    echo "----------------------------------------------------------------------"
    
    # Locate sample forward reads for the specific dataset
    # Expects files named with structure: DATA_DIR/GSE132074_sample1_R1.fastq.gz
    forward_reads=($(ls "$DATA_DIR"/${ds}_*_R1.fastq.gz 2>/dev/null))
    
    if [ ${#forward_reads[@]} -eq 0 ]; then
        echo "Warning: No raw sequencing files discovered for cohort $ds in $DATA_DIR. Skipping..."
        continue
    fi

    # Create temporary scratch directories for processing steps
    mkdir -p "${DATA_DIR}/${ds}_trimmed"
    mkdir -p "${DATA_DIR}/${ds}_aligned"

    for r1_file in "${forward_reads[@]}"; do
        # Extract base sample prefix identifier
        sample_base=$(basename "$r1_file" _R1.fastq.gz)
        r2_file="${DATA_DIR}/${sample_base}_R2.fastq.gz"
        
        if [ ! -f "$r2_file" ]; then
            echo "Error: Paired-end reverse read file missing for $sample_base. Aborting loop."
            exit 1
        fi

        echo "Executing adapter trimming via Trimmomatic for: $sample_base"
        trimmomatic PE -threads 8 \
            "$r1_file" "$r2_file" \
            "${DATA_DIR}/${ds}_trimmed/${sample_base}_R1_trimmed.fastq.gz" "${DATA_DIR}/${ds}_trimmed/${sample_base}_R1_unpaired.fastq.gz" \
            "${DATA_DIR}/${ds}_trimmed/${sample_base}_R2_trimmed.fastq.gz" "${DATA_DIR}/${ds}_trimmed/${sample_base}_R2_unpaired.fastq.gz" \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
            LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 >> "$LOG_DIR"/${ds}_trimmomatic.log 2>&1

        echo "Executing splice-aware genomic alignment via STAR for: $sample_base"
        STAR \
            --runThreadN 8 \
            --genomeDir GRCm39_STARindex \
            --readFilesIn "${DATA_DIR}/${ds}_trimmed/${sample_base}_R1_trimmed.fastq.gz" "${DATA_DIR}/${ds}_trimmed/${sample_base}_R2_trimmed.fastq.gz" \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outFileNamePrefix "${DATA_DIR}/${ds}_aligned/${sample_base}_" >> "$LOG_DIR"/${ds}_star_alignment.log 2>&1
    done

    echo "Executing digital expression quantitation via featureCounts for dataset: $ds"
    # Runs featureCounts globally across all generated coordinate-sorted BAM files for this dataset
    featureCounts \
        -T 8 \
        -p \
        -a gencode.vM32.annotation.gtf \
        -o "${RAW_COUNTS_DIR}/${ds}_raw_counts.txt" \
        "${DATA_DIR}/${ds}_aligned"/*_Aligned.sortedByCoord.out.bam >> "$LOG_DIR"/${ds}_featureCounts.log 2>&1

    # Optional: Clean up intermediate huge trimmed FASTQ and heavy BAM alignments to conserve infrastructure storage
    # rm -rf "${DATA_DIR}/${ds}_trimmed" "${DATA_DIR}/${ds}_aligned"
    
    echo "Dataset $ds completed successfully. Raw expression matrices saved to ${RAW_COUNTS_DIR}/${ds}_raw_counts.txt"
done
