#!/bin/bash

# Script to align trimmed reads to reference genome using BWA-MEM
# Usage: ./03_align_reads.sh samples.txt

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

# Check if reference genome index exists
if [ ! -f "${REF_GENOME}.bwt" ]; then
    echo "Error: BWA index for reference genome not found."
    echo "Please run 01_download_reference.sh first."
    exit 1
fi

echo "Starting read alignment at $(date)"
echo "========================================================"

# Process each sample in the samples file
while read -r SAMPLE; do
    # Skip empty lines and comments
    [[ -z "$SAMPLE" || "$SAMPLE" == \#* ]] && continue
    
    echo "Processing sample: $SAMPLE"
    
    # Check if trimmed files exist
    R1_TRIMMED="${TRIMMED_READS_DIR}/${SAMPLE}_R1_trimmed.fastq.gz"
    R2_TRIMMED="${TRIMMED_READS_DIR}/${SAMPLE}_R2_trimmed.fastq.gz"
    
    if [ ! -f "$R1_TRIMMED" ] || [ ! -f "$R2_TRIMMED" ]; then
        echo "Error: Trimmed files for sample $SAMPLE not found."
        echo "Expected: $R1_TRIMMED and $R2_TRIMMED"
        echo "Please run 02_trim_reads.sh first."
        continue
    fi
    
    # Output files
    SAM="${ALIGNED_READS_DIR}/${SAMPLE}.sam"
    BAM="${ALIGNED_READS_DIR}/${SAMPLE}.bam"
    SORTED_BAM="${ALIGNED_READS_DIR}/${SAMPLE}.sorted.bam"
    LOG_BWA="${LOGS_DIR}/${SAMPLE}_bwa.log"
    LOG_SORT="${LOGS_DIR}/${SAMPLE}_sort.log"
    
    # Skip if sorted BAM already exists
    if [ -f "$SORTED_BAM" ]; then
        echo "Sorted BAM file for $SAMPLE already exists, skipping alignment."
        continue
    fi
    
    # Align with BWA-MEM
    echo "Aligning reads for sample: $SAMPLE"
    ${BWA} mem \
        -t ${THREADS} \
        -M \
        -R "@RG\\tID:${SAMPLE}\\tSM:${SAMPLE}\\tPL:ILLUMINA\\tLB:${SAMPLE}\\tPU:unit1" \
        ${REF_GENOME} \
        ${R1_TRIMMED} \
        ${R2_TRIMMED} \
        > ${SAM} \
        2> ${LOG_BWA}
    
    # Check if alignment was successful
    if [ $? -ne 0 ]; then
        echo "Error: BWA alignment failed for sample: $SAMPLE"
        echo "Check log file: $LOG_BWA"
        continue
    fi
    
    # Convert SAM to BAM
    echo "Converting SAM to BAM for sample: $SAMPLE"
    ${SAMTOOLS} view -bS ${SAM} > ${BAM}
    
    # Check if conversion was successful
    if [ $? -ne 0 ]; then
        echo "Error: SAM to BAM conversion failed for sample: $SAMPLE"
        continue
    fi
    
    # Sort BAM file
    echo "Sorting BAM file for sample: $SAMPLE"
    ${SAMTOOLS} sort -@ ${THREADS} -o ${SORTED_BAM} ${BAM} 2> ${LOG_SORT}
    
    # Check if sorting was successful
    if [ $? -ne 0 ]; then
        echo "Error: BAM sorting failed for sample: $SAMPLE"
        echo "Check log file: $LOG_SORT"
        continue
    fi
    
    # Index BAM file
    echo "Indexing BAM file for sample: $SAMPLE"
    ${SAMTOOLS} index ${SORTED_BAM}
    
    # Check if indexing was successful
    if [ $? -ne 0 ]; then
        echo "Error: BAM indexing failed for sample: $SAMPLE"
        continue
    fi
    
    # Remove intermediate files
    echo "Removing intermediate SAM and unsorted BAM files"
    rm -f ${SAM} ${BAM}
    
    echo "Alignment completed for sample: $SAMPLE"
    
done < "$SAMPLES_FILE"

echo "Read alignment completed at $(date)"
echo "========================================================"
