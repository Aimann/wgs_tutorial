#!/bin/bash

# Script to perform variant calling using GATK Mutect2
# Usage: ./05_variant_calling.sh samples.txt [normal_sample]
# If normal_sample is provided, tumor-normal analysis will be performed
# Otherwise, tumor-only mode will be used

# Source configuration
source ./config.sh

# Check if samples file was provided
if [ "$#" -lt 1 ]; then
    echo "Usage: $0 samples.txt [normal_sample]"
    echo "samples.txt format should be: sample_name (one per line)"
    echo "If normal_sample is provided, all other samples will be treated as tumor samples"
    exit 1
fi

SAMPLES_FILE=$1
NORMAL_SAMPLE=$2

# Check if samples file exists
if [ ! -f "$SAMPLES_FILE" ]; then
    echo "Error: Samples file $SAMPLES_FILE does not exist."
    exit 1
fi

# Check if reference genome exists
if [ ! -f "$REF_GENOME" ] || [ ! -f "${REF_GENOME}.fai" ]; then
    echo "Error: Reference genome or index not found."
    echo "Please run 01_download_reference.sh first."
    exit 1
fi

# Check if dbSNP exists
if [ ! -f "$DBSNP" ]; then
    echo "Warning: dbSNP file not found. Variant calling will proceed without germline resource."
    GERMLINE_RESOURCE_ARG=""
else
    GERMLINE_RESOURCE_ARG="--germline-resource ${DBSNP}"
fi

echo "Starting variant calling at $(date)"
echo "========================================================"

# Process each sample in the samples file
while read -r SAMPLE; do
    # Skip empty lines and comments
    [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && continue
    
    # Skip if this is the normal sample
    if [ "$SAMPLE" == "$NORMAL_SAMPLE" ]; then
        echo "Skipping normal sample: $SAMPLE (will be used as reference)"
        continue
    fi
    
    echo "Processing sample: $SAMPLE"
    
    # Check if deduplicated BAM exists
    DEDUP_BAM="${PROCESSED_BAMS_DIR}/${SAMPLE}.dedup.bam"
    
    if [ ! -f "$DEDUP_BAM" ]; then
        echo "Error: Deduplicated BAM file for sample $SAMPLE not found."
        echo "Expected: $DEDUP_BAM"
        echo "Please run 04_mark_duplicates.sh first."
        continue
    fi
    
    # Set output files
    if [ -z "$NORMAL_SAMPLE" ]; then
        # Tumor-only mode
        OUTPUT_VCF="${VARIANTS_DIR}/${SAMPLE}.vcf.gz"
        FILTERED_VCF="${VARIANTS_DIR}/${SAMPLE}.filtered.vcf.gz"
        LOG="${LOGS_DIR}/${SAMPLE}_mutect2.log"
        
        # Skip if output already exists
        if [ -f "$FILTERED_VCF" ]; then
            echo "Filtered VCF for $SAMPLE already exists, skipping."
            continue
        fi
        
        echo "Running Mutect2 in tumor-only mode for sample: $SAMPLE"
        
        # Run Mutect2
        ${GATK} --java-options "-Xmx${JAVA_MEM}" Mutect2 \
            -R ${REF_GENOME} \
            -I ${DEDUP_BAM} \
            ${GERMLINE_RESOURCE_ARG} \
            -O ${OUTPUT_VCF} \
            2> ${LOG}
        
    else
        # Tumor-normal mode
        NORMAL_BAM="${PROCESSED_BAMS_DIR}/${NORMAL_SAMPLE}.dedup.bam"
        
        # Check if normal BAM exists
        if [ ! -f "$NORMAL_BAM" ]; then
            echo "Error: Normal BAM file not found: $NORMAL_BAM"
            echo "Please run 04_mark_duplicates.sh for the normal sample first."
            continue
        fi
        
        OUTPUT_VCF="${VARIANTS_DIR}/${SAMPLE}_vs_${NORMAL_SAMPLE}.vcf.gz"
        FILTERED_VCF="${VARIANTS_DIR}/${SAMPLE}_vs_${NORMAL_SAMPLE}.filtered.vcf.gz"
        LOG="${LOGS_DIR}/${SAMPLE}_vs_${NORMAL_SAMPLE}_mutect2.log"
        
        # Skip if output already exists
        if [ -f "$FILTERED_VCF" ]; then
            echo "Filtered VCF for $SAMPLE vs $NORMAL_SAMPLE already exists, skipping."
            continue
        fi
        
        echo "Running Mutect2 in tumor-normal mode for sample: $SAMPLE vs $NORMAL_SAMPLE"
        
        # Run Mutect2
        ${GATK} --java-options "-Xmx${JAVA_MEM}" Mutect2 \
            -R ${REF_GENOME} \
            -I ${DEDUP_BAM} \
            -I ${NORMAL_BAM} \
            -normal ${NORMAL_SAMPLE} \
            ${GERMLINE_RESOURCE_ARG} \
            -O ${OUTPUT_VCF} \
            2> ${LOG}
    fi
    
    # Check if Mutect2 was successful
    if [ $? -ne 0 ]; then
        echo "Error: Mutect2 failed for sample: $SAMPLE"
        echo "Check log file: $LOG"
        continue
    fi
    
    # Filter variants
    FILTER_LOG="${LOGS_DIR}/${SAMPLE}_filter.log"
    echo "Filtering variants for sample: $SAMPLE"
    
    ${GATK} --java-options "-Xmx${JAVA_MEM}" FilterMutectCalls \
        -R ${REF_GENOME} \
        -V ${OUTPUT_VCF} \
        -O ${FILTERED_VCF} \
        2> ${FILTER_LOG}
    
    # Check if filtering was successful
    if [ $? -eq 0 ]; then
        echo "Variant calling and filtering completed for sample: $SAMPLE"
    else
        echo "Error: Variant filtering failed for sample: $SAMPLE"
        echo "Check log file: $FILTER_LOG"
    fi
    
done < "$SAMPLES_FILE"

echo "Variant calling completed at $(date)"
echo "========================================================"
