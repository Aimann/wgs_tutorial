# WGS Processing Pipeline Tutorial

## Overview
This tutorial provides a step-by-step guide for processing Whole Genome Sequencing (WGS) reads, from quality control to variant calling. The workflow includes adapter trimming, reference genome indexing, read alignment, duplicate removal, and somatic variant calling.

## Installation
```bash
conda create -c conda-forge --name heini bwa samtools picard zlib bcftools openssl=1.0
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

### 2. Reference Genome Preparation
Download the reference genome and build BWA index:

```bash
# Download human reference genome (GRCh38)
cd genome
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
gunzip GRCh38.primary_assembly.genome.fa.gz

# Build BWA index
bwa index -a bwtsw GRCh38.primary_assembly.genome.fa
```

### 3. Read Alignment with BWA-MEM
Align the trimmed reads to the reference genome:

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

### 4. Mark and Remove Duplicates with Picard
Identify and mark PCR and optical duplicates:

```bash
picard MarkDuplicates \
    I=${OUTPUT_DIR}/${SAMPLE}.sorted.bam \
    O=${OUTPUT_DIR}/${SAMPLE}.dedup.bam \
    M=${OUTPUT_DIR}/${SAMPLE}.metrics.txt \
    REMOVE_DUPLICATES=true

# Index the deduplicated BAM
samtools index processed_bams/sample.dedup.bam
```

### 5. Variant Calling with bcftools
Call variants with bcftools:
```bash
bcftools mpileup -f reference.fa alignments.bam | bcftools call -mv -Ov -o calls.vcf
```

### 6. Filter Variants
Filter variants to improve quality:
bcftools view -i 'QUAL>=20' calls.vcf > calls.filtered.vcf
```bash
bcftools view -i 'QUAL>=20' calls.vcf > calls.filtered.vcf
```


## Output Files
- `trimmed_reads/`: Quality-controlled FASTQ files
- `aligned_reads/`: Initial SAM/BAM alignment files
- `processed_bams/`: Deduplicated and recalibrated BAM files
- `variants/`: VCF files containing called variants


## References
- BWA-MEM: [https://github.com/lh3/bwa](https://github.com/lh3/bwa)
- Picard Tools: [https://broadinstitute.github.io/picard/](https://broadinstitute.github.io/picard/)
- GATK Best Practices: [https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows](https://gatk.broadinstitute.org/hc/en-us/sections/360007226651-Best-Practices-Workflows)
