# ICC Pipeline - WES Analysis

DNAseq Analysis Toolkit for Cardiovascular Disease Research. This pipeline is designed for Whole Exome Sequencing (WES) analysis with standardized quality control, alignment, and variant calling workflows.

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Pipeline Architecture](#pipeline-architecture)
- [Workflow Steps](#workflow-steps)
- [Configuration](#configuration)
- [Output Structure](#output-structure)
- [Resource Tracking](#resource-tracking)

## Installation

1. **Clone the repository:**
    ```sh
    git clone https://github.com/yourusername/icc-pipeline.git
    cd icc-pipeline
    ```

2. **Configure reference data:**
    Update `workflow/config.yml` with appropriate paths to:
    - Reference genome (GRCh38 or GRCh37)
    - Target regions (gene panels)
    - Known variant databases (dbSNP, 1000G, etc.)

## Usage

Run the pipeline using the `cidna.py` wrapper script:

```sh
./cidna.py run workflow/config.yml -i /path/to/input -o /path/to/output -- --cores 8
```

**Optional flags:**
- `-i, --inputdir`: Input directory with FASTQ files (required)
- `-o, --outdir`: Output directory for results (required)
- `--verbose`: Enable verbose output
- Additional Snakemake arguments after `--`

## Pipeline Architecture

The pipeline follows a modular sequential design:
1. **Raw Data QC** → Extract sequencing metrics
2. **Adapter Trimming** → Remove low-quality bases
3. **Trimmed QC** → Verify trimming quality
4. **Read Alignment** → Map to reference genome
5. **BAM Processing** → Coordinate sorting, deduplication
6. **BAM QC** → Coverage metrics and flagstat
7. **Variant Calling** → Identify variants with HaplotypeCaller
8. **Variant Filtering** → Apply quality filters
9. **Annotation** → VEP/SnpEff annotation (optional)
10. **Summary** → Generate per-sample and cohort variant reports

## Workflow Steps

| Step | Rule File | Tool | Input | Output |
|------|-----------|------|-------|--------|
| 01 | `001_qc.smk` | FastQC | FASTQ files | HTML/ZIP reports |
| 02 | `002_trimming.smk` | fastp | Raw FASTQ | Trimmed FASTQ |
| 03 | `003_posttrim_qc.smk` | FastQC | Trimmed FASTQ | HTML/ZIP reports |
| 04 | `004_alignment.smk` | BWA-MEM2 + Sambamba | Trimmed FASTQ | Sorted BAM |
| 05 | `005_bam_prep.smk` | Sambamba + GATK | BAM | Processed BAM |
| 06 | `006_bam_qc.smk` | Samtools/GATK | BAM | Coverage/Flagstat |
| 07 | `007_variant_calling.smk` | GATK HaplotypeCaller | BAM | gVCF / VCF |
| 08 | `008_variant_filtering.smk` | GATK VariantFiltration | VCF | Filtered VCF |
| 09 | `009_annotation.smk` | VEP | VCF | Annotated VCF |
| 10 | `010_summary.smk` | Custom scripts | Filtered VCFs | Sample/cohort variant reports |

## Configuration

Edit `workflow/config.yml` to customize:

**Thread allocation:**
```yaml
threads_high: 11
threads_mid: 4
threads_low: 1
```

**Reference genomes (GRCh38 or GRCh37):**
```yaml
reference_genome: "/path/to/grch38.fa"
icc_panel: "/path/to/target_regions.bed"
```

**Tool parameters:**
```yaml
fastp:
  min_read_length: 35
  window_size: 5
gatk:
  HaplotypeCaller:
    dcovg: 1000
```

## Output Structure

```
output_dir/
├── analysis/
│   ├── 001_qc/pretrim/          # Pre-trimming FastQC reports
│   ├── 002_trimming/            # Trimmed FASTQ files
│   ├── 003_qc/posttrim/         # Post-trimming FastQC reports
│   ├── 004_alignment/           # Aligned BAM files
│   ├── 005_bam_prep/            # Processed BAM files
│   ├── 006_qc/bam/              # BAM QC metrics
│   ├── 007_variant_calling/     # gVCF/VCF files
│   ├── 008_variant_filtering/   # Filtered VCF files
│   ├── 009_annotation/          # Annotated variants
│   └── 010_summary/             # Sample and cohort variant reports
├── logs/                        # Execution logs per rule
├── benchmarks/                  # Resource usage per rule
└── results/                     # Final outputs
```

## Resource Tracking

Pipeline execution generates a resource usage report including:
- Runtime duration
- CPU and memory usage
- Network I/O statistics
- Output file sizes

Reports are saved to `benchmarks/resource_usage.txt`

## Variant Reporting

The workflow now includes a local, ClawBio-inspired reporting layer after variant filtering.

- Per-sample outputs in `analysis/010_summary/<sample>/`:
  - `variant_summary.md`
  - `variant_summary.tsv`
  - `variant_summary.json`
- Cohort outputs in `analysis/010_summary/`:
  - `cohort_variant_report.md`
  - `cohort_variant_summary.tsv`
  - `cohort_variant_summary.json`

These summaries are generated from the filtered SNP and indel VCFs and include:

- total and PASS variant counts
- SNP/indel breakdown
- genotype zygosity counts when sample genotypes are present
- transition/transversion summary for SNPs
- chromosome-level burden table
- top PASS variants ranked by QUAL

## Implementation Notes

### GATK Changes from InHouse Pipeline
- Removed `IndelRealigner` and `RealignerTargetCreator`: HaplotypeCaller in GATK4 performs realignment on-the-fly
- Replaced `samtools` with `sambamba`: Equivalent filtering with improved performance
- Kept `DepthOfCoverage`: Provides comprehensive coverage metrics

### Sample Naming Convention
- Input: `SAMPLE_ID_SX_LYYYY_RZ_NNN.fastq.gz`
  - X = sample number
  - Y = lane number  
  - Z = read direction (1/2)
  - N = chunk number
- Output: Organized by sample ID with lane tracking

## Troubleshooting

**Pipeline fails during sample discovery:**
- Verify FASTQ file naming matches expected pattern
- Check `samplesfile` in Snakefile points to correct CSV

**Dry-run before execution:**
```sh
./cidna.py run workflow/config.yml -i input/ -o output/ -- --dry-run
```

**View DAG visualization:**
```sh
snakemake -s workflow/Snakefile --dag | dot -Tpng > dag.png
```


