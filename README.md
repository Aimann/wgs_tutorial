# WGS Processing Pipeline Tutorial

## Overview
This tutorial provides a step-by-step guide for processing Whole Genome Sequencing (WGS) reads, from quality control to variant calling. The workflow includes adapter trimming, reference genome indexing, read alignment, duplicate removal, and somatic variant calling.

## Installation

### 1. Set up conda environment with required tools

```bash
conda create -c conda-forge -c bioconda --name wgs bwa samtools bedtools picard zlib bcftools openssl=1.0 
```

## Workflow

### 2. Data Preparation
Navigate to a directory of your choice to perform the analysis

```bash
cd /path/to/directory
```

Create a working directory for your analysis:

```bash
mkdir -p wgs/{alignments,genome,fastq}
```

Move into the newly created directory

```bash
cd wgs
```
Note the path to your fastq files with your input reads
```bash
# Example (replace with your actual data paths)
/path/to/sample_R1.fastq.gz /path/to/sample_R2.fastq.gz
```

OR

(optional)
Download or copy your FASTQ files to the `fastq` directory:
```bash
cp /path/to/sample_R1.fastq.gz /path/to/sample_R2.fastq.gz fastq/
```

### 3. Reference Genome Preparation

Navigate to the `genome` directory
```bash
cd genome
```

Download the latest human reference genome from GENCODE (hg38/GRCh38)

```bash
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_47/GRCh38.primary_assembly.genome.fa.gz
```

Unzip the file

```bash
gunzip GRCh38.primary_assembly.genome.fa.gz
```

Build the bwa index for aligning the reads to the human reference genome

```bash
bwa index -a bwtsw GRCh38.primary_assembly.genome.fa
```

### 4. Read Alignment with BWA-MEM

Navigate back to your main `wgs` directory

if you are already in the `genome` directory
```bash
cd ..
```

OR 

if you are anywhere else
```bash
cd /path/to/wgs
```

Align the trimmed reads to the reference genome:

```bash
bwa mem -t 8 -M \
  -R '@RG\tID:unique_readgroup_id\tSM:your_sample_name\tPL:sequencing_platform\tLB:your_library_id\tPU:unique_platform_unit'
  genome/GRCh38.primary_assembly.genome.fa \
  fastq/sample_R1_trimmed.fq.gz \
  fastq/sample_R2_trimmed.fq.gz \
  > alignments/sample.sam
```

Example -R
```bash
-R '@RG\tID:L001\tSM:Sample1\tLB:WGSH\tPL:ILLUMINA'
```

# Convert SAM to BAM
```bash
samtools view -bS alignments/sample.sam > alignments/sample.bam
```

# Sort BAM file
```bash
samtools sort -o alignments/sample.sorted.bam alignments/sample.bam
```

# Index BAM file
```bash
samtools index alignments/sample.sorted.bam
```

### 5. Mark and Remove Duplicates with Picard
Identify and mark PCR and optical duplicates:

```bash
picard MarkDuplicates \
    I=alignments/${SAMPLE}.sorted.bam \
    O=alignments/${SAMPLE}.dedup.bam \
    M=alignments/${SAMPLE}.metrics.txt \
    REMOVE_DUPLICATES=true
```

### 6. Variant Calling with bcftools
Call variants with bcftools:
```bash
bcftools mpileup -f genome/GRCh38.primary_assembly.genome.fa ${SAMPLE}.dedup.bam | bcftools call -mv -Ov -o calls.vcf
```

### 7. Filter Variants
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
