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

Run the pipeline using the `wes_pipeline.py` CLI wrapper:

```sh
# Local multi-core execution (single node)
./wes_pipeline.py run -i /path/to/input -o /path/to/output --cores 16

# Distributed HPC Slurm execution across all 4 nodes (320 cores)
./wes_pipeline.py run -i ~/project/input -o ~/project/results --slurm --jobs 320

# Or use the cluster runner script (recommended inside tmux):
./run_cluster.sh -i ~/project/input -o ~/project/results
```

**Available CLI subcommands:**
- `run`: Execute the full WES pipeline (supports `--slurm` for cluster-wide scheduling).
- `cluster-status`: Check current Slurm node status, allocated/idle core counts, and running jobs.
- `sync-data`: Fast rsync helper between QNAP and shared `/home` storage.
- `plan`: Preview execution plan and render high-resolution DAG/rulegraph diagrams (`results/dag.png`, `results/rulegraph.png`).
- `validate`: Perform standalone pre-flight configuration and resource checks (automatically downloads GRCh38 if missing).
- `download-ref`: Download and index the GRCh38 reference genome (`resources/ref/grch38/GRCh38.primary_assembly.genome.fa`).
- `download-giab`: Download NIST Genome in a Bottle (GIAB) HG001 / HG002 benchmark truth sets for GRCh38.
- `benchmark`: Benchmark pipeline variant calls (GATK & DeepVariant) against GIAB high-confidence truth standard (GA4GH).
- `report`: Generate a Snakemake HTML execution report.

**Options for `run`:**
- `-i, --inputdir`: Input directory containing sample FASTQ files (required, must reside in `/home`)
- `-o, --outdir`: Target output directory for pipeline results (required, must reside in `/home`)
- `--slurm`: Enable Slurm cluster execution across 4 nodes (320 cores total)
- `--jobs`: Maximum concurrent jobs submitted to Slurm (default: 320)
- `--default-mem`: Default memory fallback in MB (default: 4000)
- `--deepvariant`: Enable dual variant calling (GATK + DeepVariant)
- `--skip-validation`: Bypass pre-flight configuration and resource checks
- `--verbose`: Enable detailed verbose logging

## Pipeline Architecture

The pipeline follows a modular sequential design:
1. **Trimming & QC (fastp)** → Adapter removal, poly-G clipping, and comprehensive pre/post-filter QC reports
2. **Read Alignment** → Map to reference genome with BWA-MEM2
3. **BAM Processing** → Coordinate sorting, Sambamba deduplication, GATK4 BQSR
4. **BAM QC** → Exon-level coverage metrics and flagstat
5. **Variant Calling** → Distributed parallel HaplotypeCaller (and optional DeepVariant)
6. **Variant Filtering** → Apply quality filters (GATK VariantFiltration / bcftools)
7. **Annotation** → Ensembl-VEP annotation & ACMG clinical classification
8. **Summary** → Cohort variant reports and MultiQC dashboard

## Workflow Steps

| Step | Rule File | Tool | Input | Output |
|------|-----------|------|-------|--------|
| 00 | `000_ref.smk` | Python | - | Downloaded & Staged GRCh38 Genome |
| 01 | `001_trimming.smk` | fastp | Raw FASTQ | Trimmed FASTQ + HTML/JSON QC Reports |
| 02 | `002_alignment.smk` | BWA-MEM2 + Samtools | Trimmed FASTQ | Coordinate-sorted BAM |
| 03 | `003_bam_prep.smk` | Sambamba + GATK4 | BAM | Deduplicated & BQSR-recalibrated BAM |
| 04 | `004_bam_qc.smk` | Samtools + Bedtools | BAM | Exon-level Coverage & Flagstat Reports |
| 05 | `005_variant_calling.smk` | GATK4 HaplotypeCaller + DeepVariant | BAM | Distributed gVCF / VCF |
| 06 | `006_variant_filtering.smk` | GATK4 / bcftools | VCF | High-confidence Filtered SNPs & Indels |
| 07 | `007_annotation.smk` | Ensembl-VEP + GeneBe | Filtered VCF | VEP Annotated VCF & ACMG TSV |
| 08 | `008_summary.smk` | Custom Python | VCF / TSV | Markdown, TSV, and JSON Cohort Summaries |
| 09 | `009_multiqc.smk` | MultiQC | fastp JSON / QC Metrics | Aggregated Interactive MultiQC HTML Report |
| 10 | `010_benchmark_giab.smk` | GA4GH Benchmark | Filtered VCF | GIAB NIST Accuracy Metrics & Dashboard |


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
│   ├── 000_ref/                 # Staged reference indices and interval chunks
│   ├── 001_trimming/            # Trimmed FASTQ files and fastp QC reports
│   ├── 002_alignment/           # BWA aligned & merged BAM files
│   ├── 003_bam_prep/            # Deduplicated, BQSR, and filtered BAM files
│   ├── 004_bam_qc/              # Exon coverage & flagstat QC reports
│   ├── 005_variant_calling/     # gVCF/VCF files (GATK & DeepVariant)
│   ├── 006_variant_filtering/   # High-confidence filtered SNPs & Indels
│   ├── 007_annotation/          # VEP annotations & ACMG classification
│   ├── 008_summary/             # Sample & cohort variant reports and dashboards
│   └── 010_benchmark/giab/      # GIAB NIST benchmark verification
├── logs/                        # Execution logs per rule (000..010)
├── benchmarks/                  # Resource usage per rule (000..010)
└── results/                     # Final MultiQC & coverage reports
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

- Per-sample outputs in `analysis/008_summary/<sample>/`:
  - `variant_summary.md`
  - `variant_summary.tsv`
  - `variant_summary.json`
- Cohort outputs in `analysis/008_summary/`:
  - `cohort_variant_report.md`
  - `cohort_variant_summary.tsv`
  - `cohort_variant_summary.json`
  - `cohort_variant_dashboard.html`

These summaries are generated from the filtered SNP and indel VCFs and include:

- total and PASS variant counts
- SNP/indel breakdown
- genotype zygosity counts when sample genotypes are present
- transition/transversion summary for SNPs
- chromosome-level burden table
- top PASS variants ranked by QUAL

When `analysis/007_annotation/<sample>.annotated.vcf` is present, the summary layer also picks up annotation-aware fields without making VEP a hard dependency. The dashboard and markdown reports then include:

- gene-level burden summaries
- impact tiers such as `HIGH` and `MODERATE`
- consequence labels from VEP/SnpEff-style annotations
- clinical significance tags when present in the annotated VCF

## GIAB Benchmarking

Evaluate pipeline variant calling accuracy against the NIST Genome in a Bottle (GIAB) benchmark truth sets (HG001 / NA12878 and HG002 / NA24385 on GRCh38) using GA4GH standards.

1. **Download GIAB Truth Sets:**
   ```sh
   ./wes_pipeline.py download-giab --sample HG001
   ```

2. **Benchmark Variant Callers (GATK vs. DeepVariant):**
   ```sh
   ./wes_pipeline.py benchmark -o /path/to/output --sample HG001 --download-truth
   ```

Outputs generated in `analysis/010_benchmark/giab/`:

- `HG001_giab_benchmark_report.md`: Markdown summary of Recall, Precision, F1-scores, and Ti/Tv ratios.
- `HG001_giab_benchmark_summary.tsv`: Machine-readable TSV matrix of SNP and Indel GA4GH metrics.
- `HG001_giab_benchmark_summary.json`: JSON output for programmatic integration.
- `HG001_giab_benchmark_dashboard.html`: Aesthetic interactive dashboard visualizing caller accuracy.

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
./wes_pipeline.py run workflow/config.yml -i input/ -o output/ -- --dry-run
```

**Generate DAG & Rulegraph visualization:**
```sh
./wes_pipeline.py plan workflow/config.yml -i input/
```


