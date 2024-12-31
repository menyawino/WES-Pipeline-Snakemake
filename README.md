# ICC Pipeline

This pipeline is designed for DNA sequencing analysis, specifically tailored for cardiology research. It leverages Snakemake for workflow management and Conda for environment management.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Workflow Overview](#workflow-overview)
- [Rules](#rules)
  - [Quality Control](#quality-control)
  - [Trimming](#trimming)
  - [Alignment](#alignment)
  - [Variant Calling](#variant-calling)
  - [Variant Filtering](#variant-filtering)
  - [Annotation](#annotation)
- [Resource Tracking](#resource-tracking)

## Installation

1. **Clone the repository:**
    ```sh
    git clone https://github.com/yourusername/icc-pipeline.git
    cd icc-pipeline
    ```

2. **Install Conda environments:**
    ```sh
    conda env create -f environment.yml
    ```

3. **Activate the environment:**
    ```sh
    conda activate icc_pipeline
    ```

## Usage

To run the pipeline, use the `cidna.py` script with the required arguments:

```sh
./cidna.py run workflow/config.yml -i /path/to/input -o /path/to/output -- -c88 --printshellcmds --rerun-incomplete
```

## Workflow Overview

The pipeline consists of several steps, each managed by Snakemake rules. The main steps include:

1. **Quality Control:** Initial quality checks on raw sequencing data.
2. **Trimming:** Removing adapters and low-quality bases.
3. **Alignment:** Aligning reads to the reference genome.
4. **Variant Calling:** Identifying variants from the aligned reads.
5. **Variant Filtering:** Filtering the identified variants.
6. **Annotation:** Annotating the filtered variants.

## Rules

### Quality Control

- **Rule:** `raw_fastqc`
- **Description:** Runs FastQC on raw sequencing data.
- **Input:** Raw FASTQ files.
- **Output:** FastQC reports.

### Trimming

- **Rule:** `trimming`
- **Description:** Trims adapters and low-quality bases using Prinseq.
- **Input:** Raw FASTQ files.
- **Output:** Trimmed FASTQ files.

### Alignment

- **Rule:** `bwa_alignment`
- **Description:** Aligns trimmed reads to the reference genome using BWA.
- **Input:** Trimmed FASTQ files.
- **Output:** BAM files.

### Variant Calling

- **Rule:** `haplotypecaller`
- **Description:** Calls variants using GATK HaplotypeCaller.
- **Input:** Realigned BAM files.
- **Output:** VCF files.

### Variant Filtering

- **Rule:** `filter_variants`
- **Description:** Filters variants using GATK VariantFiltration.
- **Input:** VCF files.
- **Output:** Filtered VCF files.

### Annotation

- **Rule:** `annotate_variants`
- **Description:** Annotates variants using VEP.
- **Input:** Filtered VCF files.
- **Output:** Annotated VCF files.

## Resource Tracking

The pipeline tracks resource usage, including CPU, memory, and network usage. A detailed report is generated at the end of the pipeline run.


## Inhouse vs new pipeline comparison

- Removed RealignerTargetCreator because gatk4 HaplotypeCaller realigns the reads on the fly.

