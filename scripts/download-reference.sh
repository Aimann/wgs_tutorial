#!/bin/bash

# Script to download and prepare reference genome
# Usage: ./01_download_reference.sh

# Source configuration
source ./config.sh

echo "Starting reference genome download and preparation at $(date)"
echo "========================================================"

cd ${PIPELINE_DIR}/reference

# Download reference genome if not present
if [ ! -f "hg38.fa" ]; then
  echo "Downloading human reference genome (GRCh38)..."
  wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
  gunzip hg38.fa.gz
  echo "Reference genome downloaded successfully!"
else
  echo "Reference genome already exists, skipping download."
fi

# Download dbSNP if not present
if [ ! -f "Homo_sapiens_assembly38.dbsnp138.vcf" ]; then
  echo "Downloading dbSNP database..."
  wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
  wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx
  echo "dbSNP database downloaded successfully!"
else
  echo "dbSNP database already exists, skipping download."
fi

# Check if BWA index exists, create if not
if [ ! -f "${REF_GENOME}.bwt" ]; then
  echo "Building BWA index for reference genome..."
  ${BWA} index -a bwtsw ${REF_GENOME} 2> ${LOGS_DIR}/bwa_index.log
  echo "BWA index built successfully!"
else
  echo "BWA index already exists, skipping creation."
fi

# Create FASTA index for reference
if [ ! -f "${REF_GENOME}.fai" ]; then
  echo "Creating FASTA index..."
  ${SAMTOOLS} faidx ${REF_GENOME}
  echo "FASTA index created successfully!"
else
  echo "FASTA index already exists, skipping creation."
fi

# Create sequence dictionary for GATK
if [ ! -f "${REF_GENOME%.fa}.dict" ]; then
  echo "Creating sequence dictionary..."
  ${PICARD} CreateSequenceDictionary \
    R=${REF_GENOME} \
    O=${REF_GENOME%.fa}.dict
  echo "Sequence dictionary created successfully!"
else
  echo "Sequence dictionary already exists, skipping creation."
fi

echo "Reference genome preparation completed at $(date)"
echo "========================================================"
