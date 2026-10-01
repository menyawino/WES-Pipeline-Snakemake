#!/usr/bin/env python3

import sys
import os
import click

from workflow.scripts.utils import (
    build_folders,
    run_snakemake,
    run_snakemake_plan,
    get_snakemake_report,
    GRE,
    NC
)
from workflow.scripts.validate import validate_pipeline_config

def print_banner():
    """Print pipeline CLI header banner."""
    banner = f"""
    ╭─────────────────────────────────────────────────────────────╮
    │                                                             │
    │         ██████  █████  ██████  ██████  ██  ██████           │
    │        ██      ██   ██ ██   ██ ██   ██ ██ ██    ██          │
    │        ██      ███████ ██████  ██   ██ ██ ██    ██          │
    │        ██      ██   ██ ██   ██ ██   ██ ██ ██    ██          │
    │         ██████ ██   ██ ██   ██ ██████  ██  ██████           │
    │                                                             │
    │     ██ ███    ██ ██████  ██████  ██████  ███      ███       │
    │     ██ ████   ██ ██     ██    ██ ██   ██ ████    ████       │
    │     ██ ██ ██  ██ ██████ ██    ██ ██████  ██ ██  ██ ██       │
    │     ██ ██  ██ ██ ██     ██    ██ ██   ██ ██  ████  ██       │
    │     ██ ██   ████ ██      ██████  ██   ██ ██   ██   ██       │
    │                                                             │
    │   ██████  ███    ██  █████       ███████ ███████  ██████    │
    │   ██   ██ ████   ██ ██   ██      ██      ██      ██    ██   │
    │   ██   ██ ██ ██  ██ ███████ ████ ███████ ███████ ██    ██   │
    │   ██   ██ ██  ██ ██ ██   ██           ██ ██      ██    ██   │
    │   ██████  ██   ████ ██   ██      ███████ ███████  ██████▄   │
    │                                                             │
    │{GRE} DNAseq Analysis Toolkit for Cardiovascular Disease Research{NC} │
    │                    Author: {GRE}Omar Ahmed{NC}                       │
    ╰─────────────────────────────────────────────────────────────╯
"""
    print(banner)

@click.group()
def cli():
    """WES Analysis Pipeline CLI - Whole Exome Sequencing workflow automation."""
    pass

@cli.command(context_settings={"ignore_unknown_options": True})
@click.option('-i', '--inputdir', required=True, help="Path to raw FASTQ input directory.")
@click.option('-o', '--outdir', required=True, help="Path to pipeline output directory.")
@click.option('-c', '--configfile', default='workflow/config.yml', show_default=True, help="Path to Snakemake config YAML.")
@click.option('-s', '--samplesfile', default='', show_default=True, help="Optional sample metadata CSV (leave empty for auto-discovery).")
@click.option('--slurm', is_flag=True, default=False, help="Execute workflow across the cluster via Slurm (4 nodes x 80 = 320 cores).")
@click.option('--jobs', default=320, type=int, show_default=True, help="Maximum concurrent Slurm jobs.")
@click.option('--default-mem', default=4000, type=int, show_default=True, help="Default memory in MB for Slurm jobs.")
@click.option('--cores', default=None, type=int, help="Number of CPU cores for local mode (default: from config / 88).")
@click.option('-n', '--dry-run', is_flag=True, help="Perform dry-run without executing jobs.")
@click.option('--deepvariant', is_flag=True, default=False, help="Enable dual variant calling with Google DeepVariant alongside GATK.")
@click.option('--singularity-args', default='', show_default=True, help="Bind mounts for Singularity containers.")
@click.option('--email', default=None, help="Recipient email address(es) for completion/failure notifications.")
@click.option('--verbose', is_flag=True, help="Enable verbose output logging.")
@click.option('--skip-validation', is_flag=True, help="Skip pre-flight configuration and resource validation checks.")
@click.argument('snakemake_args', nargs=-1)
def run(inputdir, outdir, configfile, samplesfile, slurm, jobs, default_mem, cores, dry_run, deepvariant, singularity_args, email, verbose, skip_validation, snakemake_args):
    """Execute the WES analysis pipeline."""
    if not skip_validation:
        print(f"{GRE}[INFO] Running pre-flight pipeline validation...{NC}")
        if not validate_pipeline_config(configfile, inputdir, outdir):
            print("\n[ERROR] Pipeline validation failed. Fix configuration or pass --skip-validation to bypass.")
            sys.exit(1)
        print(f"{GRE}[INFO] Validation passed successfully.{NC}\n")

    build_folders(outdir)

    return_code = run_snakemake(
        configfile=configfile,
        inputdir=inputdir,
        outdir=outdir,
        samplesfile=samplesfile,
        cores=cores,
        slurm=slurm,
        jobs=jobs,
        default_mem=default_mem,
        dry_run=dry_run,
        deepvariant=deepvariant,
        singularity_args=singularity_args,
        verbose=verbose,
        extra_args=snakemake_args,
        email=email
    )

    if return_code == 0 or return_code is None:
        print(f"\n{GRE}[SUCCESS] WES Pipeline execution completed successfully.{NC}")
    else:
        print(f"\n\033[91m[FAILURE] Pipeline terminated early with errors (exit code: {return_code}).{NC}")
        sys.exit(return_code)

@cli.command()
@click.option('-c', '--configfile', default='workflow/config.yml', show_default=True, help="Path to Snakemake config YAML.")
@click.option('-i', '--inputdir', help="Optional path to input directory for sample DAG generation.")
def plan(configfile, inputdir):
    """Preview pipeline execution DAG and render aesthetic topology diagrams."""
    print(f"{GRE}[INFO] Generating aesthetic workflow DAG and rulegraph diagrams...{NC}")
    run_snakemake_plan(configfile, inputdir=inputdir)
    print(f"{GRE}[INFO] Rendered graphs saved to results/dag.png and results/rulegraph.png{NC}")

@cli.command()
@click.option('-i', '--inputdir', required=True, help="Path to raw FASTQ input directory.")
@click.option('-o', '--outdir', required=True, help="Path to pipeline output directory.")
@click.option('-c', '--configfile', default='workflow/config.yml', show_default=True, help="Path to Snakemake config YAML.")
def validate(inputdir, outdir, configfile):
    """Run standalone pre-flight configuration and resource validation checks."""
    print(f"{GRE}[INFO] Running pre-flight pipeline validation...{NC}")
    success = validate_pipeline_config(configfile, inputdir, outdir)
    if success:
        print(f"\n{GRE}[SUCCESS] All pre-flight checks passed! Pipeline is ready to run.{NC}")
    else:
        print("\n\033[91m[ERROR] Pre-flight validation failed. Please address the errors above.{NC}")
        sys.exit(1)

@cli.command()
@click.option('-c', '--configfile', default='workflow/config.yml', show_default=True, help="Path to Snakemake config YAML.")
@click.option('--force', is_flag=True, help="Force re-download even if target file exists.")
def download_ref(configfile, force):
    """Download and index GRCh38 reference genome."""
    import yaml
    from workflow.scripts.download_ref import download_reference_genome
    
    print(f"{GRE}[INFO] Downloading GRCh38 reference genome data...{NC}")
    target_path = "resources/ref/grch38/Homo_sapiens_assembly38.fasta"
    genome = "grch38"
    
    if os.path.exists(configfile):
        try:
            with open(configfile, 'r') as f:
                cfg = yaml.safe_load(f) or {}
                target_path = cfg.get('reference_genome', target_path)
                genome = cfg.get('Genome', genome)
        except Exception as e:
            print(f"Warning: Could not parse config file: {e}")

    success = download_reference_genome(target_path=target_path, genome=genome, force=force)
    if success:
        print(f"\n{GRE}[SUCCESS] Reference genome GRCh38 downloaded and indexed successfully!{NC}")
    else:
        print("\n\033[91m[ERROR] Failed to download or index GRCh38 reference genome.{NC}")
        sys.exit(1)

@cli.command()
@click.option('-c', '--configfile', default='workflow/config.yml', show_default=True, help="Path to Snakemake config YAML.")
def report(configfile):
    """Generate a Snakemake HTML execution report."""
    print(f"{GRE}[INFO] Generating HTML execution report...{NC}")
    get_snakemake_report(configfile)

@cli.command("download-giab")
@click.option('--sample', default='HG001', type=click.Choice(['HG001', 'NA12878', 'HG002', 'NA24385']), show_default=True, help="GIAB sample name.")
@click.option('--outdir', default='resources/giab/HG001', show_default=True, help="Destination directory for truth files.")
@click.option('--force', is_flag=True, help="Force re-download even if files already exist.")
def download_giab(sample, outdir, force):
    """Download NIST Genome in a Bottle (GIAB) benchmark truth sets for GRCh38."""
    from workflow.scripts.download_giab_truth import download_giab_truth_set
    print(f"{GRE}[INFO] Fetching GIAB benchmark truth set for {sample}...{NC}")
    res = download_giab_truth_set(sample=sample, outdir=outdir, force=force)
    if res:
        print(f"\n{GRE}[SUCCESS] GIAB {sample} truth data ready in {outdir}!{NC}")
    else:
        print("\n\033[91m[ERROR] Failed to download GIAB truth data.{NC}")
        sys.exit(1)

@cli.command()
@click.option('-o', '--outdir', required=True, help="Path to pipeline output directory containing variant filtering results.")
@click.option('-s', '--sample', default='HG001', show_default=True, help="Sample name to benchmark.")
@click.option('-c', '--configfile', default='workflow/config.yml', show_default=True, help="Path to Snakemake config YAML.")
@click.option('--giab-reference', default='HG001', type=click.Choice(['HG001', 'NA12878', 'HG002', 'NA24385']), show_default=True, help="GIAB reference truth set to use.")
@click.option('--caller', 'callers', multiple=True, default=['gatk', 'deepvariant'], show_default=True, help="Variant callers to evaluate.")
@click.option('--truth-vcf', default=None, help="Path to GIAB truth VCF.gz (default: auto-detected from resources/giab/<ref>/).")
@click.option('--truth-bed', default=None, help="Path to GIAB high-confidence BED (default: auto-detected).")
@click.option('--target-bed', default=None, help="Path to target capture BED (default: from config).")
@click.option('--download-truth', is_flag=True, help="Automatically download GIAB truth data if missing.")
def benchmark(outdir, sample, configfile, giab_reference, callers, truth_vcf, truth_bed, target_bed, download_truth):
    """Benchmark caller VCFs against GIAB high-confidence truth sets (GA4GH standard)."""
    import yaml
    import subprocess
    from workflow.scripts.download_giab_truth import download_giab_truth_set

    giab_ref = "HG002" if giab_reference.upper() in ["HG002", "NA24385"] or sample.upper() in ["HG002", "NA24385"] else "HG001"
    default_truth_dir = os.path.join("resources", "giab", giab_ref)

    if not truth_vcf:
        truth_vcf = os.path.join(default_truth_dir, f"{giab_ref}_GRCh38_1_22_v4.2.1_benchmark.vcf.gz")
    if not truth_bed:
        truth_bed = os.path.join(default_truth_dir, f"{giab_ref}_GRCh38_1_22_v4.2.1_benchmark.bed")

    if not target_bed and os.path.exists(configfile):
        try:
            with open(configfile, "r") as f:
                cfg = yaml.safe_load(f) or {}
                target_bed = cfg.get("icc_panel")
        except Exception:
            pass

    if (not os.path.exists(truth_vcf) or not os.path.exists(truth_bed)) and download_truth:
        print(f"{GRE}[INFO] Truth files missing. Downloading GIAB truth data for {giab_sample}...{NC}")
        download_giab_truth_set(sample=giab_sample, outdir=default_truth_dir)

    if not os.path.exists(truth_vcf) or not os.path.exists(truth_bed):
        print(f"\n\033[91m[ERROR] Truth files not found at {truth_vcf} or {truth_bed}.\nPass --download-truth or run: ./wes_pipeline.py download-giab --sample {giab_sample}{NC}")
        sys.exit(1)

    print(f"\n{GRE}======================================================={NC}")
    print(f"{GRE} Running GIAB Benchmark for Sample: {sample}{NC}")
    print(f" Output Directory : {outdir}")
    print(f" Evaluated Callers: {', '.join(callers)}")
    print(f" Truth VCF        : {truth_vcf}")
    print(f" Confident BED    : {truth_bed}")
    print(f"{GRE}=======================================================\n{NC}")

    cmd = [
        sys.executable,
        "workflow/scripts/benchmark_giab.py",
        "--outdir", outdir,
        "--sample", sample,
        "--callers", *callers,
        "--truth-vcf", truth_vcf,
        "--truth-bed", truth_bed,
    ]
    if target_bed and os.path.exists(target_bed):
        cmd.extend(["--target-bed", target_bed])

    res = subprocess.run(cmd)
    if res.returncode == 0:
        print(f"\n{GRE}[SUCCESS] GIAB benchmark completed successfully!{NC}")
    else:
        print(f"\n\033[91m[FAILURE] Benchmark evaluation failed with exit code: {res.returncode}{NC}")
        sys.exit(res.returncode)

@cli.command("cluster-status")
def cluster_status():
    """Display current Slurm cluster nodes, core allocations, and queued jobs."""
    import subprocess
    print(f"\n{GRE}=== Cluster Node & Core Allocation (sinfo) ==={NC}")
    try:
        subprocess.run(["sinfo", "-o", "%10P %10N %10c %20C %15m %10e %10T"], check=False)
        print(f"\n{GRE}=== Core Counts (Allocated/Idle/Other/Total) ==={NC}")
        subprocess.run(["sinfo", "-o", "%C"], check=False)
    except FileNotFoundError:
        print("[WARNING] 'sinfo' command not found. Ensure Slurm client utilities are in PATH.")

    print(f"\n{GRE}=== Running / Queued Jobs (squeue) ==={NC}")
    try:
        user = os.environ.get("USER", "omar")
        subprocess.run(["squeue", "-u", user], check=False)
    except FileNotFoundError:
        print("[WARNING] 'squeue' command not found.")
    print()

@cli.command("sync-data")
@click.option('--direction', type=click.Choice(['in', 'out'], case_sensitive=False), required=True, help="'in' (from QNAP to ~/project) or 'out' (from ~/project to QNAP).")
@click.option('--qnap-path', required=True, help="Source/Destination path on QNAP storage.")
@click.option('--local-path', required=True, help="Source/Destination path on local shared /home.")
def sync_data(direction, qnap_path, local_path):
    """High-performance rsync helper to stage data between QNAP and local /home."""
    import subprocess
    if direction.lower() == 'in':
        print(f"{GRE}[INFO] Staging input data from QNAP to shared /home...{NC}")
        cmd = ["rsync", "-ah", "--info=progress2", qnap_path.rstrip('/') + '/', local_path.rstrip('/') + '/']
    else:
        print(f"{GRE}[INFO] Archiving results from shared /home to QNAP...{NC}")
        cmd = ["rsync", "-ah", "--info=progress2", local_path.rstrip('/') + '/', qnap_path.rstrip('/') + '/']
    print(f"Executing: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)

def main():
    """CLI Main Entry Point."""
    print_banner()
    cli()

if __name__ == '__main__':
    main()
