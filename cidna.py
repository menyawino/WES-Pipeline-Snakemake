#!/usr/bin/env python

from workflow.scripts.utils import *
from workflow.scripts.validate import validate_pipeline_config

@click.group()
def cli():
    """Define the CLI group."""
    pass

@click.command(context_settings={"ignore_unknown_options": True})
@click.argument('configfile')
@click.option('-i', '--inputdir', required=True, help="Input directory.")
@click.option('-o', '--outdir', required=True, help="Output directory.")
@click.option('--verbose', is_flag=True, help="Enable verbose output.")
@click.option('--skip-validation', is_flag=True, help="Skip pre-flight validation checks.")
@click.argument('snakemake_args', nargs=-1)

def run(configfile, inputdir, outdir, verbose, skip_validation, snakemake_args):
    """Execute workflow (using Snakemake underneath)."""
    
    if not skip_validation:
        if not validate_pipeline_config(configfile, inputdir, outdir):
            print("\nPipeline validation failed. Use --skip-validation to bypass checks (not recommended).")
            sys.exit(1)
        print()
    
    build_folders(outdir)

    return_code = run_snakemake(configfile, inputdir, outdir, verbose=verbose,
                                extra_args=snakemake_args)

    if return_code == 0 or return_code is None:
        print("Pipeline completed successfully.")
    else:
        print("Pipeline terminated early due to errors.")


cli.add_command(run)

def main():
    """Main entry point."""
    print(f"""
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
""")
    cli()

if __name__ == '__main__':
    main()
