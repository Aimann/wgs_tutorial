# Overview

RiboMarker is a pipeline for mapping and feature counting of RNA-seq data. It is designed to process ribosome profiling data and identify different types of RNA features, such as rRNA, tRNA, miRNA, lncRNA, snRNA, snoRNA, and protein-coding genes. The pipeline uses bowtie or bowtie2 for mapping. It also provides options for trimming the input FASTQ files and collapsing the reads before alignment.

The pipeline flows as follows:

1. Trim the input FASTQ files (optional).
2. Collapse the reads (optional).
3. Map the reads to the reference genome using bowtie or bowtie2.
4. Parse the mapping results and count the reads for each feature.

## Single-End Workflow
![Overview](assets/overview.single.png)
## Paired-End Workflow
![Overview](assets/overview.paired.png)

# Installation

1. Clone this repository:
   ```bash
   git clone https://github.com/Aimann/rsbalign-priv.git
    ```

2. Build the Docker image:
   ```bash
   docker build -t rsbalign-priv .
   ```

3. Verify that the Docker image was created successfully:
   ```bash
   docker images
   ```

# Running the test dataset

To run the test dataset, follow these steps:

1. Navigate to the `test_data` directory:
   ```bash
   cd test_data
   ```

2. Modify the `input_dir` parameter in the `config.yaml` file to specify the path to the `test_data` directory:
   ```yaml
   input_dir: /path/to/rsb-pub-run/test_data/input/
   ```

3. Run the pipeline using the wrapper script:
   ```bash
   ./run_pipeline.sh config.single.yaml
   ```
   This script will run the pipeline using the `config.single.yaml` configuration file.

4. Alternatively, you can run paired-end data using the `config.paired.yaml` configuration file:
   ```bash
   ./run_pipeline.sh config.paired.yaml
   ```

5. The output files will be saved in the specified `output_prefix` directory.

# Output files

The pipeline generates the following output files in the following structure (using `RiboMarker` as an example `output_prefix`):

```
/path/to/your/data/RiboMarker/
├── alignments/
│   ├── sample1.bam
│   ├── sample2.bam
│   ├── ...
├── logs/
│   ├── trim.log (if trimming is 'True')
│   ├── map.log
│   ├── parse.log
│   ├── ...
├── db/
│   ├── species/
│   │   ├── genome/
│   │   ├── rrna/
│   │   ├── mirna/
│   │   ├── trna/
│   │   ├── ...
├── trimmed/ (if trimming is 'True')
│   ├── sample1.trimmed.fq.gz
│   ├── sample2.trimmed.fq.gz
│   ├── ...
├── collapse/ (if collapse_reads is 'True')
│   ├── sample1.collapse.fa
│   ├── sample2.collapse.fa
│   ├── ...
├── RiboMarker-trimstats.csv (if trimming is 'True') # Trimming statistics
├── RiboMarker-mapstats.csv # Mapping statistics
├── RiboMarker-readcounts.csv # Read counts for all features
├── RiboMarker-featuretypes.csv # Feature types
├── RiboMarker-typecounts.csv # RNA type counts
├── RiboMarker-typelengths.csv # Read lengths for each RNA type
├── RiboMarker-filtered-readcounts.csv # Read counts filtered by read minimum and fraction samples minimum and for biotypes
├── ...
```

# Usage

## Required input files

The pipeline requires the following input files:

path_to_data: /path/to/your/data
sample_sheet: samples.csv
species: human (e.g., human, mouse, rat)
fastq_files: raw or trimmed fastq files

The `samples.csv` file should have the following format:
   ```
   sample1,group1,sample1.fastq.gz
   sample2,group1,sample2.fastq.gz
   sample3,group2,sample3.fastq.gz
   sample4,group2,sample4.fastq.gz
   ```

For paired-end reads the `samples.csv` file should have the following format:
   ```
   sample1,group1,sample1.R1.fastq.gz,sample1.R2.fastq.gz
   sample2,group1,sample2.R1.fastq.gz,sample2.R2.fastq.gz
   sample3,group2,sample3.R1.fastq.gz,sample3.R2.fastq.gz
   sample4,group2,sample4.R1.fastq.gz,sample4.R2.fastq.gz
   ```

The sample sheet file should be specified in the `config.yaml` file. Do not include spaces or special characters in the sample IDs or group IDs, only use letters, numbers, and underscores.

## 1. Create a directory for your data and place your input files there.

```bash
mkdir /path/to/your/data
```
This directory should contain the following files:
   - `samples.csv`: Sample sheet file with the following format: sample_id,group_id,fastq1
   - `sample1.fastq.gz`, `sample2.fastq.gz`, ...: FASTQ files containing the raw sequencing reads
```
/path/to/your/data/
   ├── samples.csv
   ├── sample1.fastq.gz
   ├── sample2.fastq.gz
   ├── ...
   ├── pairs.csv (optional)
```

## 2. Modify the `config.yaml` file to specify the paths to your input files and the desired output directory.
Required parameters:
   ```yaml
   input_dir: /path/to/your/data # Path to the directory containing the input files
   output_prefix: RiboMarker     # Prefix to use for the output files
   samples: samples.csv          # Path to the sample sheet file with the following format: sample_id,group_id,fastq1
   species: human                # Species to use for alignment and feature counting (e.g., human, mouse, rat)
   trim_fastq: True/False        # Whether to trim the input FASTQ files, True/False
   collapse_reads: True/False    # Whether to collapse the reads before alignment, True/False
   paired_end: True/False        # Whether input reads are paired-end or single-end
   ```
Optional parameters:
   ```yaml
   pairs: pairs.csv              # Optional: Sample pairs for DESeq2 comparisons (relative to input_dir). Format "Group1,Group2"
   subsample_reads: 0            # Subsample trimmed reads to this amount for each sample. Default 0; include all reads
   dump_unmapped: True/False     # Whether to dump the unmapped reads to a separate bam file, True/False
   trim_only: True/False         # Whether to only trim the input FASTQ files, True/False
   count_ribospike: True/False   # Whether to count RiboSpike reads, True/False
   multimap: 100                 # Maximum number of multimapping reads to keep
   read_minimum: 10              # Minimum # of reads a feature must have to be included in the analysis; combined with fraction_minimum
   fraction_samples_minimum: 0.1 # Minimum fraction of total reads a feature must have to be included in the analysis; combined with read_minimum
   unique_only: True/False       # Whether to only inlcude uniquely mapping reads in the output, True/False
   ```

   - Modify the `species` parameter to specify the species to use for alignment and feature counting (e.g., human, mouse, rat).

   - Modify the `trim_fastq` parameter to specify whether to trim the input FASTQ files. If set to `True`, the pipeline will trim the input FASTQ files using Cutadapt.

   - Modify the `collapse_reads` parameter to specify whether to collapse the reads before alignment. If set to `True`, the pipeline will collapse the reads prior to mapping. This significantly speeds up the mapping process; however it may make downstream bam manipulation more difficult.


