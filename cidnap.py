#!/usr/bin/env python

import sys
import os
import subprocess
import click
import time
import psutil
import shutil
from datetime import timedelta
import re
from tqdm import tqdm

# build folders for the pipeline: analysis, benchmarks, results, logs if they don't exist
def build_folders(outdir):
    """Build folders for the pipeline if they don't exist."""
    folders = [os.path.join(outdir, subfolder) for subfolder in ['analysis', 'benchmarks', 'results', 'logs']]
    for folder in folders:
        if not os.path.exists(folder):
            os.makedirs(folder)
       
# Track the resources used by the pipeline     
def track_resources(start_time, net_start, outdir, verbose=False):
    """Track and display resource usage."""
    end_time = time.time()
    elapsed_time = end_time - start_time

    # Get CPU and memory usage
    cpu_usage = psutil.cpu_percent(interval=1)
    memory_info = psutil.virtual_memory()
    memory_usage = memory_info.used / (1024 ** 3)  # Convert to GB
    total_memory = memory_info.total / (1024 ** 3)  # Convert to GB

    # Get network usage
    net_end = psutil.net_io_counters()
    bytes_sent = (net_end.bytes_sent - net_start.bytes_sent) / (1024 ** 2)  # Convert to MB
    bytes_recv = (net_end.bytes_recv - net_start.bytes_recv) / (1024 ** 2)  # Convert to MB

    # Get the size of the files in Gb
    output_folder_size = get_folder_size(os.path.join(outdir, 'results')) / (1024 ** 3) + get_folder_size(os.path.join(outdir, 'analysis')) / (1024 ** 3)


    # Print resource usage stats
    print(f"""
     _____________________________________________________________
    |                                                             |
    |                 Pipeline Resource Usage Report              |
    |_____________________________________________________________|
    |   Runtime:                    {str(timedelta(seconds=elapsed_time))}                |
    |   CPU Usage:                  {cpu_usage}%                          |
    |   Memory Usage:               {memory_usage:.2f} GB / {total_memory:.2f} GB           |
    |   Network Sent:               {bytes_sent:.2f} MB                       |
    |   Network Received:           {bytes_recv:.2f} MB                       |
    |   Output Files Size:         {output_folder_size:.2f} GB                       |
    |_____________________________________________________________|
    """)
    
    # output the resource usage to a file with date and time before runtime
    with open(os.path.join(outdir, 'benchmarks/resource_usage.txt'), 'w') as f:
        f.write(f"""
     _____________________________________________________________
    |                                                             |
    |                 Pipeline Resource Usage Report              |
    |_____________________________________________________________|
    |   Date and Time:              {time.ctime()}                |
    |   Runtime:                    {str(timedelta(seconds=elapsed_time))}                |
    |   CPU Usage:                  {cpu_usage}%                          |
    |   Memory Usage:               {memory_usage:.2f} GB / {total_memory:.2f} GB           |
    |   Network Sent:               {bytes_sent:.2f} MB                       |
    |   Network Received:           {bytes_recv:.2f} MB                       |
    |   Output Files Size:         {output_folder_size:.2f} GB                       |
    |_____________________________________________________________|
    """)
    

#  get the size of a folder in bytes 
def get_folder_size(folder):
    """Return the size of a folder in bytes."""
    total_size = 0
    for dirpath, dirnames, filenames in os.walk(folder):
        for f in filenames:
            fp = os.path.join(dirpath, f)
            total_size += os.path.getsize(fp)
    return total_size


# run snakemake with the specified options and configuration
def run_snakemake(configfile, inputdir, outdir, verbose=False, extra_args=[]):
    """Run Snakemake with the specified options and configuration."""

    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, 'workflow/Snakefile')

    # Basic Snakemake command
    cmd = ["snakemake", "-s", snakefile, "--use-conda", "-k", "--benchmark-extended"]

    # Add additional Snakemake arguments
    cmd += list(extra_args)

    if configfile:
        # Only add the specified config file without defaults and system confs
        cmd += ["--configfile", configfile]

    # Add input and output directories to the command
    cmd += ["--directory", outdir]
    cmd += ["--config", f"inputdir={inputdir}"]

    # Print the final command if verbose with cmd list as a string
    if verbose:
        print('Command executed:', ' '.join(cmd))

    start_time = time.time()
    net_start = psutil.net_io_counters()  # Capture network usage at start

    # Initialize progress bar
    total_jobs = 0
    completed_jobs = 0
    running_jobs = 0
    progress_bar = None

    # Run Snakemake and capture stdout
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    for line in process.stdout:
        # Parse total jobs
        if "total" in line:
            match = re.search(r'total\s+(\d+)', line)
            if match:
                total_jobs = int(match.group(1))
                progress_bar = tqdm(total=total_jobs, desc="Jobs", unit="job")

        # Parse running jobs
        if "Job" in line and "Reason" in line:
            running_jobs += 1
            if progress_bar:
                elapsed_time = time.time() - start_time
                progress_bar.set_postfix_str(f"{completed_jobs}/{total_jobs} done, {running_jobs} running, elapsed time: {str(timedelta(seconds=elapsed_time))}")

        # Parse finished jobs
        if "Finished job" in line:
            completed_jobs += 1
            running_jobs -= 1
            if progress_bar:
                progress_bar.update(1)
                progress_bar.set_postfix_str(f"{completed_jobs}/{total_jobs} done, {running_jobs} running")

    process.wait()

    if progress_bar:
        progress_bar.close()

    # Check for errors
    if process.returncode != 0:
        print(f'Error in Snakemake invocation: {process.returncode}', file=sys.stderr)
        return process.returncode

    # Display resource usage
    track_resources(start_time, net_start, outdir, verbose=verbose)
    print("Pipeline completed successfully.")


# run snakemake plan to preview the pipeline
def run_snakemake_plan(configfile):
    """Preview the Snakemake plan before running the pipeline."""
    
    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, 'workflow/Snakefile')
    
    os.system("snakemake -s " + snakefile + " --use-conda --dag \
        --configfile " + configfile + " --quiet \
        | dot -Tpng > results/dag.png")
    os.system("snakemake -s " + snakefile + " --use-conda --rulegraph \
        --configfile " + configfile + " --quiet \
        | dot -Tpng > results/rulegraph.png")


# generate snakemake report for the pipeline
def get_snakemake_report(configfile):
    """Generate a Snakemake report for the pipeline."""
    
    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, 'workflow/Snakefile')
    
    os.system("snakemake -s " + snakefile + " --use-conda --report " + 
              "--configfile " + configfile + " --quiet")
    

@click.group()
def cli():
    """Define the CLI group."""
    pass


@click.command(context_settings={"ignore_unknown_options": True})
@click.argument('configfile')
@click.option('-i', '--inputdir', required=True, help="Input directory.")
@click.option('-o', '--outdir', required=True, help="Output directory.")
@click.option('--verbose', is_flag=True, help="Enable verbose output.")
@click.argument('snakemake_args', nargs=-1)
def run(configfile, inputdir, outdir, snakemake_args, verbose):
    """Execute workflow (using Snakemake underneath)."""
    build_folders(outdir)
    
    start_time = time.time()
    net_start = psutil.net_io_counters()
    
    return_code = run_snakemake(configfile, inputdir, outdir, verbose=verbose,
                                extra_args=snakemake_args)
    
    if return_code == 0:
        track_resources(start_time, net_start, outdir, verbose=verbose)
    else:
        print("Pipeline terminated early due to errors.")


cli.add_command(run)


# ANSI color codes
GRE = '\033[92m'  # Green color
NC = '\033[0m'    # No color


def main():
    """Main entry point."""
    print(f"""
 _____________________________________________________________
|                                                             |
|         ██████  █████  ██████  ██████  ██  ██████           |
|        ██      ██   ██ ██   ██ ██   ██ ██ ██    ██          |
|        ██      ███████ ██████  ██   ██ ██ ██    ██          |
|        ██      ██   ██ ██   ██ ██   ██ ██ ██    ██          |
|         ██████ ██   ██ ██   ██ ██████  ██  ██████           |
|                                                             |
|     ██ ███    ██ ██████  ██████  ██████  ███      ███       |
|     ██ ████   ██ ██     ██    ██ ██   ██ ████    ████       |
|     ██ ██ ██  ██ ██████ ██    ██ ██████  ██ ██  ██ ██       |
|     ██ ██  ██ ██ ██     ██    ██ ██   ██ ██  ████  ██       |
|     ██ ██   ████ ██      ██████  ██   ██ ██   ██   ██       |
|                                                             |
|   ██████  ███    ██  █████       ███████ ███████  ██████    |
|   ██   ██ ████   ██ ██   ██      ██      ██      ██    ██   |
|   ██   ██ ██ ██  ██ ███████ ████ ███████ ███████ ██    ██   |
|   ██   ██ ██  ██ ██ ██   ██           ██ ██      ██    ██   |
|   ██████  ██   ████ ██   ██      ███████ ███████  ██████▄   |
|                                                             |
| {GRE}      DNAseq Analysis Toolkit for Cardiology Research{NC}       |
|                    Author: {GRE}Omar Ahmed{NC}                       |
|_____________________________________________________________|

""")
    cli()

if __name__ == '__main__':
    main()
