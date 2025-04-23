#!/bin/bash

# Master script to run the entire WGS processing pipeline
# Usage: ./run_pipeline.sh samples.txt [normal_sample]

# Check if samples file was provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 samples.txt [normal_sample]"
    echo "samples.txt format should be: sample_name (one per line)"
    echo "If normal_sample is provided, tumor-normal analysis will be performed"
    exit 1
fi

SAMPLES_FILE=$1
NORMAL_SAMPLE=$2

# Check if samples file exists
if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file $SAMPLES_FILE does not exist."
    exit 1
fi

# Make scripts executable
chmod +x config.sh
chmod +x 01_download_reference.sh
chmod +x 02_trim_reads.sh
chmod +x 03_align_reads.sh
chmod +x 04_mark_duplicates.sh
chmod +x 05_variant_calling.sh

# Start time
START_TIME=$(date +%s)
echo "Starting WGS processing pipeline at $(date)"
echo "========================================================"

# Step 1: Download and prepare reference genome
echo "Step 1: Downloading and preparing reference genome"
./01_download_reference.sh
if [ $? -ne 0 ]; then
    echo "Reference genome preparation failed. Exiting."
    exit 1
fi

# Step 2: Trim adapters
echo "Step 2: Trimming adapters"
./02_trim_reads.sh $SAMPLES_FILE
if [ $? -ne 0 ]; then
    echo "Adapter trimming failed. Exiting."
    exit 1
fi

# Step 3: Align reads
echo "Step 3: Aligning reads to reference genome"
./03_align_reads.sh $SAMPLES_FILE
if [ $? -ne 0 ]; then
    echo "Read alignment failed. Exiting."
    exit 1
fi

# Step 4: Mark duplicates
echo "Step 4: Marking duplicates"
./04_mark_duplicates.sh $SAMPLES_FILE
if [ $? -ne 0 ]; then
    echo "Duplicate marking failed. Exiting."
    exit 1
fi

# Step 5: Variant calling
echo "Step 5: Calling variants"
if [ -z "$NORMAL_SAMPLE" ]; then
    ./05_variant_calling.sh $SAMPLES_FILE
else
    ./05_variant_calling.sh $SAMPLES_FILE $NORMAL_SAMPLE
fi

if [ $? -ne 0 ]; then
    echo "Variant calling failed. Exiting."
    exit 1
fi

# End time
END_TIME=$(date +%s)
DURATION=$((END_TIME - START_TIME))
HOURS=$((DURATION / 3600))
MINUTES=$(( (DURATION % 3600) / 60 ))
SECONDS=$((DURATION % 60))

echo "========================================================"
echo "WGS processing pipeline completed at $(date)"
echo "Total runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s"
echo "========================================================"

# Generate a summary report
source ./config.sh
REPORT="${PIPELINE_DIR}/pipeline_summary.txt"

echo "WGS Processing Pipeline Summary" > $REPORT
echo "=======================================================" >> $REPORT
echo "Date: $(date)" >> $REPORT
echo "Runtime: ${HOURS}h ${MINUTES}m ${SECONDS}s" >> $REPORT
echo "Samples processed:" >> $REPORT
cat $SAMPLES_FILE | grep -v '^#' | grep -v '^$' | wc -l >> $REPORT
echo "=======================================================" >> $REPORT
echo "Sample statistics:" >> $REPORT

# Process each sample and add to report
while read -r SAMPLE; do
    # Skip empty lines and comments
    [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && continue
    
    echo "Sample: $SAMPLE" >> $REPORT
    
    # Check if variant file exists
    if [ -z "$NORMAL_SAMPLE" ]; then
        VARIANT_FILE="${VARIANTS_DIR}/${SAMPLE}.filtered.vcf.gz"
    else
        VARIANT_FILE="${VARIANTS_DIR}/${SAMPLE}_vs_${NORMAL_SAMPLE}.filtered.vcf.gz"
    fi
    
    if [ -f "$VARIANT_FILE" ]; then
        echo "  Variant calling: SUCCESS" >> $REPORT
        echo "  Variant count: $(zcat $VARIANT_FILE | grep -v '^#' | wc -l)" >> $REPORT
    else
        echo "  Variant calling: FAILED" >> $REPORT
    fi
    
    # Add metrics from duplicate marking
    METRICS="${PROCESSED_BAMS_DIR}/${SAMPLE}.metrics.txt"
    if [ -f "$METRICS" ]; then
        DUPLICATION_RATE=$(grep -A 1 "LIBRARY" $METRICS | tail -n 1 | cut -f 9)
        echo "  Duplication rate: $DUPLICATION_RATE" >> $REPORT
    fi
    
    echo "" >> $REPORT
    
done < "$SAMPLES_FILE"

echo "Pipeline summary saved to: $REPORT"
