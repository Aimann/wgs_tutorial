#!/bin/bash

# Configuration file for WGS processing pipeline
# Edit these parameters according to your environment

# Directories
export PIPELINE_DIR="$PWD"
export RAW_DATA_DIR="${PIPELINE_DIR}/raw_data"
export TRIMMED_READS_DIR="${PIPELINE_DIR}/trimmed_reads"
export ALIGNED_READS_DIR="${PIPELINE_DIR}/aligned_reads"
export PROCESSED_BAMS_DIR="${PIPELINE_DIR}/processed_bams"
export VARIANTS_DIR="${PIPELINE_DIR}/variants"
export LOGS_DIR="${PIPELINE_DIR}/logs"

# Reference files
export REF_GENOME="${PIPELINE_DIR}/reference/hg38.fa"
export ADAPTER_FILE="/usr/share/trimmomatic/adapters/TruSeq3-PE.fa"
export DBSNP="${PIPELINE_DIR}/reference/Homo_sapiens_assembly38.dbsnp138.vcf"

# Tool paths (modify if not in PATH)
export TRIMMOMATIC="trimmomatic"
export BWA="bwa"
export SAMTOOLS="samtools"
export PICARD="picard"
export GATK="gatk"

# Performance settings
export THREADS=8
export JAVA_MEM="16G"

# Create directories if they don't exist
mkdir -p ${RAW_DATA_DIR}
mkdir -p ${TRIMMED_READS_DIR}
mkdir -p ${ALIGNED_READS_DIR}
mkdir -p ${PROCESSED_BAMS_DIR}
mkdir -p ${VARIANTS_DIR}
mkdir -p ${LOGS_DIR}
mkdir -p ${PIPELINE_DIR}/reference

# Check if tools are available
command -v ${BWA} >/dev/null 2>&1 || { echo "BWA is not installed or not in PATH. Aborting."; exit 1; }
command -v ${SAMTOOLS} >/dev/null 2>&1 || { echo "Samtools is not installed or not in PATH. Aborting."; exit 1; }
command -v ${TRIMMOMATIC} >/dev/null 2>&1 || { echo "Trimmomatic is not installed or not in PATH. Aborting."; exit 1; }
command -v ${PICARD} >/dev/null 2>&1 || { echo "Picard is not installed or not in PATH. Aborting."; exit 1; }
command -v ${GATK} >/dev/null 2>&1 || { echo "GATK is not installed or not in PATH. Aborting."; exit 1; }

echo "Configuration loaded successfully!"
