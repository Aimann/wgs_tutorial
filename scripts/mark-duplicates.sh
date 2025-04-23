#!/bin/bash

# Script to mark duplicates using Picard
# Usage: ./04_mark_duplicates.sh samples.txt

# Source configuration
source ./config.sh

# Check if samples file was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 samples.txt"
    echo "samples.txt format should be: sample_name"
    exit 1
fi

SAMPLES_FILE=$1

# Check if samples file exists
if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file $SAMPLES_FILE does not exist."
    exit 1
fi

echo "Starting duplicate marking at $(date)"
echo "========================================================"

# Process each sample in the samples file
while read -r SAMPLE; do
    # Skip empty lines and comments
    [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && continue
    
    echo "Processing sample: $SAMPLE"
    
    # Check if sorted BAM exists
    SORTED_BAM="${ALIGNED_READS_DIR}/${SAMPLE}.sorted.bam"
    
    if [ ! -f "$SORTED_BAM" ]; then
        echo "Error: Sorted BAM file for sample $SAMPLE not found."
        echo "Expected: $SORTED_BAM"
        echo "Please run 03_align_reads.sh first."
        continue
    fi
    
    # Output files
    DEDUP_BAM="${PROCESSED_BAMS_DIR}/${SAMPLE}.dedup.bam"
    METRICS="${PROCESSED_BAMS_DIR}/${SAMPLE}.metrics.txt"
    LOG="${LOGS_DIR}/${SAMPLE}_markdup.log"
    
    # Skip if output already exists
    if [ -f "$DEDUP_BAM" ]; then
        echo "Deduplication BAM for $SAMPLE already exists, skipping."
        continue
    fi
    
    echo "Marking duplicates for sample: $SAMPLE"
    
    # Run Picard MarkDuplicates
    java -Xmx${JAVA_MEM} -jar $(which picard).jar MarkDuplicates \
        INPUT=${SORTED_BAM} \
        OUTPUT=${DEDUP_BAM} \
        METRICS_FILE=${METRICS} \
        ASSUME_SORTED=true \
        VALIDATION_STRINGENCY=LENIENT \
        REMOVE_DUPLICATES=false \
        CREATE_INDEX=true \
        OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
        2> ${LOG}
    
    # Check if duplicate marking was successful
    if [ $? -eq 0 ]; then
        echo "Duplicate marking completed for sample: $SAMPLE"
    else
        echo "Error: Duplicate marking failed for sample: $SAMPLE"
        echo "Check log file: $LOG"
    fi
    
done < "$SAMPLES_FILE"

echo "Duplicate marking completed at $(date)"
echo "========================================================"
