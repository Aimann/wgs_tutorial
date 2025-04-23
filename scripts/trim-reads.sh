#!/bin/bash

# Script to trim adapter sequences from FASTQ files
# Usage: ./02_trim_reads.sh samples.txt

# Source configuration
source ./config.sh

# Check if samples file was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 samples.txt"
    echo "samples.txt format should be: sample_name"
    echo "Expects files in raw_data directory named: sample_name_R1.fastq.gz and sample_name_R2.fastq.gz"
    exit 1
fi

SAMPLES_FILE=$1

# Check if samples file exists
if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file $SAMPLES_FILE does not exist."
    exit 1
fi

echo "Starting adapter trimming at $(date)"
echo "========================================================"

# Process each sample in the samples file
while read -r SAMPLE; do
    # Skip empty lines and comments
    [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && continue
    
    echo "Processing sample: $SAMPLE"
    
    # Check if input files exist
    R1="${RAW_DATA_DIR}/${SAMPLE}_R1.fastq.gz"
    R2="${RAW_DATA_DIR}/${SAMPLE}_R2.fastq.gz"
    
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        echo "Error: Input files for sample $SAMPLE not found."
        echo "Expected: $R1 and $R2"
        continue
    fi
    
    # Output files
    OUT_R1="${TRIMMED_READS_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    OUT_R1_UNPAIRED="${TRIMMED_READS_DIR}/${SAMPLE}_R1_unpaired.fastq.gz"
    OUT_R2="${TRIMMED_READS_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
    OUT_R2_UNPAIRED="${TRIMMED_READS_DIR}/${SAMPLE}_R2_unpaired.fastq.gz"
    LOG="${LOGS_DIR}/${SAMPLE}_trimmomatic.log"
    
    # Skip if output already exists
    if [ -f "$OUT_R1" ] && [ -f "$OUT_R2" ]; then
        echo "Trimmed files for $SAMPLE already exist, skipping."
        continue
    fi
    
    echo "Trimming adapters for sample: $SAMPLE"
    
    # Run Trimmomatic
    ${TRIMMOMATIC} PE \
        -threads ${THREADS} \
        "$R1" "$R2" \
        "$OUT_R1" "$OUT_R1_UNPAIRED" \
        "$OUT_R2" "$OUT_R2_UNPAIRED" \
        ILLUMINACLIP:${ADAPTER_FILE}:2:30:10 \
        LEADING:3 TRAILING:3 \
        SLIDINGWINDOW:4:15 \
        MINLEN:36 \
        2> "$LOG"
    
    # Check if trimming was successful
    if [ $? -eq 0 ]; then
        echo "Adapter trimming completed for sample: $SAMPLE"
    else
        echo "Error: Adapter trimming failed for sample: $SAMPLE"
        echo "Check log file: $LOG"
    fi
    
done < "$SAMPLES_FILE"

echo "Adapter trimming completed at $(date)"
echo "========================================================"
