# WGS Processing Pipeline Tutorial

## Overview
This tutorial provides a step-by-step guide for processing Whole Genome Sequencing (WGS) reads, from quality control to variant calling. The workflow includes adapter trimming, reference genome indexing, read alignment, duplicate removal, and somatic variant calling.

## Prerequisites
- Linux/Unix-based system
- At least 16GB RAM (32GB+ recommended for WGS data)
- At least 500GB free disk space
- The following tools installed:
  - Trimmomatic (v0.39+)
  - BWA (v0.7.17+)
  - Samtools (v1.13+)
  - Picard Tools (v2.25.0+)
  - GATK (v4.2.0+)

## Installation
```bash
# Install required tools (Ubuntu/Debian)
conda create -c bioconda -c conda-forge --name heini bwa samtools picard zlib
```

## Workflow

### 1. Data Preparation
Create a working directory for your analysis:

```bash
mkdir -p wgs_pipeline/{alignments,variants,genome}
cd wgs_pipeline
```

Download or copy your raw FASTQ files to the `raw_data` directory:

```bash
# Example (replace with your actual data paths)
cp /path/to/sample_R1.fastq.gz /path/to/sample_R2.fastq.gz raw_data/
```

### 2. Quality Control and Adapter Trimming
Use Trimmomatic to remove adapter sequences and low-quality bases:

```bash
trimmomatic PE \
  raw_data/sample_R1.fastq.gz \
  raw_data/sample_R2.fastq.gz \
  trimmed_reads/sample_R1_trimmed.fastq.gz \
  trimmed_reads/sample_R1_unpaired.fastq.gz \
  trimmed_reads/sample_R2_trimmed.fastq.gz \
  trimmed_reads/sample_R2_unpaired.fastq.gz \
  ILLUMINACLIP:/path/to/adapters/TruSeq3-PE.fa:2:30:10 \
  LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
```

### 3. Reference Genome Preparation
Download the reference genome and build BWA index:

```bash
# Download human reference genome (GRCh38)
cd genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# Build BWA index
bwa index -a bwtsw GRCh38.primary_assembly.genome.fa
```

### 4. Read Alignment with BWA-MEM
Align the trimmed reads to the reference genome:

```bash
# Align reads
bwa mem -t 8 -M \
  hg38.fa \
  trimmed_reads/sample_R1_trimmed.fastq.gz \
  trimmed_reads/sample_R2_trimmed.fastq.gz \
  > aligned_reads/sample.sam

# Convert SAM to BAM (more compact)
samtools view -bS aligned_reads/sample.sam > aligned_reads/sample.bam

# Sort BAM file
samtools sort -o aligned_reads/sample.sorted.bam aligned_reads/sample.bam

# Index BAM file
samtools index aligned_reads/sample.sorted.bam
```

```bash
# Align reads
bwa mem -t 8 -M \
  -R '@RG\tID:unique_readgroup_id\tSM:your_sample_name\tPL:sequencing_platform\tLB:your_library_id\tPU:unique_platform_unit'
  hg38.fa \
  trimmed_reads/sample_R1_trimmed.fastq.gz \
  trimmed_reads/sample_R2_trimmed.fastq.gz \
  > aligned_reads/sample.sam

# Convert SAM to BAM (more compact)
samtools view -bS aligned_reads/sample.sam > aligned_reads/sample.bam

# Sort BAM file
samtools sort -o aligned_reads/sample.sorted.bam aligned_reads/sample.bam

# Index BAM file
samtools index aligned_reads/sample.sorted.bam
```
Example -R
```bash
-R '@RG\tID:L001\tSM:Sample1\tLB:WGSH\tPL:ILLUMINA'
```


### 5. Mark and Remove Duplicates with Picard
Identify and mark PCR and optical duplicates:

```bash
picard MarkDuplicates \
  INPUT=aligned_reads/sample.sorted.bam \
  OUTPUT=processed_bams/sample.dedup.bam \
  METRICS_FILE=processed_bams/sample.metrics.txt \
  ASSUME_SORTED=true \
  VALIDATION_STRINGENCY=LENIENT \
  REMOVE_DUPLICATES=false

# Index the deduplicated BAM
samtools index processed_bams/sample.dedup.bam
```

### 6. Base Quality Score Recalibration (Optional)
Perform BQSR to correct for systematic technical errors in base quality scores:

```bash
# Download known sites
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf
wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx

# Generate recalibration table
gatk BaseRecalibrator \
  -I processed_bams/sample.dedup.bam \
  -R hg38.fa \
  --known-sites Homo_sapiens_assembly38.dbsnp138.vcf \
  -O processed_bams/sample.recal.table

# Apply recalibration
gatk ApplyBQSR \
  -I processed_bams/sample.dedup.bam \
  -R hg38.fa \
  --bqsr-recal-file processed_bams/sample.recal.table \
  -O processed_bams/sample.recal.bam
```
### 7. Variant Calling with bcftools
Call variants with bcftools:
```bash
bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf
bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ov -o calls.vcf

```

### 8. Filter Variants
Filter variants to improve quality:


bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ob -o calls.bcf

### 7. Variant Calling with Mutect2
Call somatic variants using Mutect2:

```bash
# For tumor-only calling
gatk Mutect2 \
  -R hg38.fa \
  -I processed_bams/sample.recal.bam \
  --germline-resource Homo_sapiens_assembly38.dbsnp138.vcf \
  -O variants/sample.vcf.gz

# For tumor-normal paired analysis (if applicable)
gatk Mutect2 \
  -R hg38.fa \
  -I tumor_sample.recal.bam \
  -I normal_sample.recal.bam \
  -normal normal_sample_name \
  --germline-resource Homo_sapiens_assembly38.dbsnp138.vcf \
  -O variants/tumor_vs_normal.vcf.gz
```

### 8. Filter Variants
Filter variants to improve quality:

```bash
# Calculate contamination
gatk CalculateContamination \
  -I tumor_sample.pileups.table \
  -matched normal_sample.pileups.table \
  -O contamination.table

# Filter variants
gatk FilterMutectCalls \
  -R hg38.fa \
  -V variants/sample.vcf.gz \
  --contamination-table contamination.table \
  -O variants/sample.filtered.vcf.gz
```

## Output Files
- `trimmed_reads/`: Quality-controlled FASTQ files
- `aligned_reads/`: Initial SAM/BAM alignment files
- `processed_bams/`: Deduplicated and recalibrated BAM files
- `variants/`: VCF files containing called variants

## Tips and Troubleshooting
- Allocate sufficient memory for larger genomes (human genome requires 16GB+)
- Use `-t` parameter in BWA to specify number of threads based on your system
- If limited by disk space, remove intermediate files after each step
- Consider using FASTQC before and after trimming to assess read quality
- For production pipelines, consider using workflow managers like Snakemake or Nextflow

## References
- BWA-MEM: [https://github.com/lh3/bwa](https://github.com/lh3/bwa)
- Picard Tools: [https://broadinstitute.github.io/picard/](https://broadinstitute.github.io/picard/)
- GATK Best Practices: [https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
