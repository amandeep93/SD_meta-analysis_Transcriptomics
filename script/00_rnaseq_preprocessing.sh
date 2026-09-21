#!/bin/bash

############################################################
# MULTI-DATASET RNA-SEQ RAW UPSTREAM PREPROCESSING PIPELINE
# Script: script/01_rnaseq_preprocessing.sh
#
# Workflow:
# FASTQ → FastQC → Trimmomatic → STAR → featureCounts
#
# Reference:
#   Genome: GRCm39
#   Annotation: GENCODE vM32
#
# Expected datasets: 8 RNA-seq GEO datasets
############################################################

# Exit immediately if a command fails
set -e

############################################################
# 1. DIRECTORY CONFIGURATION
############################################################

DATA_DIR="data/RNAseq_raw"
RAW_COUNTS_DIR="data/RNAseq_counts"
LOG_DIR="results/logs"

mkdir -p "$RAW_COUNTS_DIR"
mkdir -p "$LOG_DIR"
mkdir -p "figures"

############################################################
# 2. RNA-SEQ DATASET LIST
############################################################

# These correspond to the 8 RNA-seq datasets included
# in the manuscript.

datasets=(
    "GSE132074"
    "GSE144957"
    "GSE222410"
    "GSE113754"
    "GSE211088"
    "GSE211301"
    "GSE156925"
    "GSE166831"
)

# Keep track of datasets that could not be processed
missing_datasets=()

echo "======================================================================"
echo "RNA-SEQ PREPROCESSING PIPELINE"
echo "======================================================================"
echo "Expected datasets: ${#datasets[@]}"
echo "======================================================================"

############################################################
# 3. GLOBAL FASTQC
############################################################

echo ""
echo "Starting FastQC quality assessment..."

fastqc \
    "$DATA_DIR"/*.fastq.gz \
    -o figures/ \
    2> "$LOG_DIR/fastqc_error.log"

echo "FastQC completed."
echo ""

############################################################
# 4. PROCESS EACH DATASET
############################################################

for ds in "${datasets[@]}"; do

    echo "----------------------------------------------------------------------"
    echo "Processing Dataset: $ds"
    echo "----------------------------------------------------------------------"

    ########################################################
    # Locate R1 files
    ########################################################

    forward_reads=("$DATA_DIR"/${ds}_*_R1.fastq.gz)

    if [ ! -e "${forward_reads[0]}" ]; then

        echo "ERROR: No raw sequencing files found for $ds."
        echo "Expected files such as:"
        echo "  ${DATA_DIR}/${ds}_sample_R1.fastq.gz"
        echo ""

        missing_datasets+=("$ds")
        continue
    fi

    ########################################################
    # Create dataset-specific directories
    ########################################################

    TRIM_DIR="${DATA_DIR}/${ds}_trimmed"
    ALIGN_DIR="${DATA_DIR}/${ds}_aligned"

    mkdir -p "$TRIM_DIR"
    mkdir -p "$ALIGN_DIR"

    ########################################################
    # Process each sample
    ########################################################

    for r1_file in "${forward_reads[@]}"; do

        # Extract sample identifier
        sample_base=$(basename "$r1_file" _R1.fastq.gz)

        r2_file="${DATA_DIR}/${sample_base}_R2.fastq.gz"

        ####################################################
        # Check for paired-end R2 file
        ####################################################

        if [ ! -f "$r2_file" ]; then

            echo "ERROR: Paired-end R2 file missing."
            echo "Sample: $sample_base"
            echo "Expected: $r2_file"
            echo ""

            exit 1
        fi

        echo ""
        echo "Sample: $sample_base"

        ####################################################
        # Step 1: Adapter trimming and quality filtering
        ####################################################

        echo "Running Trimmomatic..."

        trimmomatic PE \
            -threads 8 \
            "$r1_file" \
            "$r2_file" \
            "${TRIM_DIR}/${sample_base}_R1_trimmed.fastq.gz" \
            "${TRIM_DIR}/${sample_base}_R1_unpaired.fastq.gz" \
            "${TRIM_DIR}/${sample_base}_R2_trimmed.fastq.gz" \
            "${TRIM_DIR}/${sample_base}_R2_unpaired.fastq.gz" \
            ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 \
            LEADING:3 \
            TRAILING:3 \
            SLIDINGWINDOW:4:15 \
            MINLEN:36 \
            >> "$LOG_DIR/${ds}_trimmomatic.log" 2>&1

        echo "Trimmomatic completed."

        ####################################################
        # Step 2: Splice-aware alignment using STAR
        ####################################################

        echo "Running STAR alignment..."

        STAR \
            --runThreadN 8 \
            --genomeDir GRCm39_STARindex \
            --readFilesIn \
                "${TRIM_DIR}/${sample_base}_R1_trimmed.fastq.gz" \
                "${TRIM_DIR}/${sample_base}_R2_trimmed.fastq.gz" \
            --readFilesCommand zcat \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMunmapped Within \
            --outFileNamePrefix "${ALIGN_DIR}/${sample_base}_" \
            >> "$LOG_DIR/${ds}_star_alignment.log" 2>&1

        echo "STAR alignment completed."

    done

    ########################################################
    # Step 3: Gene-level quantification using featureCounts
    ########################################################

    echo ""
    echo "Running featureCounts for $ds..."

    featureCounts \
        -T 8 \
        -p \
        -B \
        -C \
        -a gencode.vM32.annotation.gtf \
        -o "${RAW_COUNTS_DIR}/${ds}_raw_counts.txt" \
        "${ALIGN_DIR}"/*_Aligned.sortedByCoord.out.bam \
        >> "$LOG_DIR/${ds}_featureCounts.log" 2>&1

    echo ""
    echo "Dataset $ds completed successfully."
    echo "Raw count matrix:"
    echo "  ${RAW_COUNTS_DIR}/${ds}_raw_counts.txt"
    echo ""

done

############################################################
# 5. FINAL DATASET CHECK
############################################################

echo "======================================================================"
echo "RNA-SEQ PREPROCESSING SUMMARY"
echo "======================================================================"

if [ ${#missing_datasets[@]} -gt 0 ]; then

    echo ""
    echo "ERROR: The following expected datasets were NOT processed:"
    printf '  %s\n' "${missing_datasets[@]}"

    echo ""
    echo "Please check that the corresponding FASTQ files are present in:"
    echo "  $DATA_DIR"

    exit 1

else

    echo ""
    echo "All ${#datasets[@]} expected RNA-seq datasets were processed successfully."

fi

echo ""
echo "Count matrices:"
echo "  $RAW_COUNTS_DIR"

echo ""
echo "Logs:"
echo "  $LOG_DIR"

echo ""
echo "FastQC results:"
echo "  figures/"

echo ""
echo "======================================================================"
echo "RNA-SEQ PREPROCESSING COMPLETED"
echo "======================================================================"
