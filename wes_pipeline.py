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
@click.argument('configfile', default='workflow/config.yml')
@click.option('-i', '--inputdir', required=True, help="Path to raw FASTQ input directory.")
@click.option('-o', '--outdir', required=True, help="Path to pipeline output directory.")
@click.option('--verbose', is_flag=True, help="Enable verbose output logging.")
@click.option('--skip-validation', is_flag=True, help="Skip pre-flight configuration and resource validation checks.")
@click.argument('snakemake_args', nargs=-1)
def run(configfile, inputdir, outdir, verbose, skip_validation, snakemake_args):
    """Execute the WES analysis pipeline."""
    if not skip_validation:
        print(f"{GRE}[INFO] Running pre-flight pipeline validation...{NC}")
        if not validate_pipeline_config(configfile, inputdir, outdir):
            print("\n[ERROR] Pipeline validation failed. Fix configuration or pass --skip-validation to bypass.")
            sys.exit(1)
        print(f"{GRE}[INFO] Validation passed successfully.{NC}\n")

    build_folders(outdir)

    return_code = run_snakemake(
        configfile,
        inputdir,
        outdir,
        verbose=verbose,
        extra_args=snakemake_args
    )

    if return_code == 0 or return_code is None:
        print(f"\n{GRE}[SUCCESS] WES Pipeline execution completed successfully.{NC}")
    else:
        print(f"\n\033[91m[FAILURE] Pipeline terminated early with errors (exit code: {return_code}).{NC}")
        sys.exit(return_code)

@cli.command()
@click.argument('configfile', default='workflow/config.yml')
@click.option('-i', '--inputdir', help="Optional path to input directory for sample DAG generation.")
def plan(configfile, inputdir):
    """Preview pipeline execution DAG and render aesthetic topology diagrams."""
    print(f"{GRE}[INFO] Generating aesthetic workflow DAG and rulegraph diagrams...{NC}")
    run_snakemake_plan(configfile, inputdir=inputdir)
    print(f"{GRE}[INFO] Rendered graphs saved to results/dag.png and results/rulegraph.png{NC}")

@cli.command()
@click.argument('configfile', default='workflow/config.yml')
@click.option('-i', '--inputdir', required=True, help="Path to raw FASTQ input directory.")
@click.option('-o', '--outdir', required=True, help="Path to pipeline output directory.")
def validate(configfile, inputdir, outdir):
    """Run standalone pre-flight configuration and resource validation checks."""
    print(f"{GRE}[INFO] Running pre-flight pipeline validation...{NC}")
    success = validate_pipeline_config(configfile, inputdir, outdir)
    if success:
        print(f"\n{GRE}[SUCCESS] All pre-flight checks passed! Pipeline is ready to run.{NC}")
    else:
        print("\n\033[91m[ERROR] Pre-flight validation failed. Please address the errors above.{NC}")
        sys.exit(1)

@cli.command()
@click.argument('configfile', default='workflow/config.yml')
def report(configfile):
    """Generate a Snakemake HTML execution report."""
    print(f"{GRE}[INFO] Generating HTML execution report...{NC}")
    get_snakemake_report(configfile)

def main():
    """CLI Main Entry Point."""
    print_banner()
    cli()

if __name__ == '__main__':
    main()
