import sys
import os
import subprocess
import click
import time
import psutil
from datetime import timedelta
import yaml

# ANSI color codes
GRE = '\033[92m'  # Green color
NC = '\033[0m'    # No color

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
    output_folder_size = get_folder_size(outdir) / (1024 ** 3)

    report = f"""
    ╭─────────────────────────────────────────────────────────────╮
    │                                                             │
    │               Pipeline Resource Usage Report                │
    │─────────────────────────────────────────────────────────────│
    │   Date and Time:         {time.ctime():<30}     │
    │   Runtime:               {str(timedelta(seconds=elapsed_time)):>20}               │
    │   CPU Usage:             {cpu_usage:>10.1f}%                        │
    │   Memory Usage:          {memory_usage:>7.2f} GB /{total_memory:>7.2f} GB             │
    │   Network Sent:          {bytes_sent:>10.2f} MB                      │
    │   Network Received:      {bytes_recv:>10.2f} MB                      │
    │   Output Files Size:     {output_folder_size:>10.2f} GB                      │
    ╰─────────────────────────────────────────────────────────────╯
    """
    print(report)

    # output the resource usage to a file with date and time before runtime
    with open(os.path.join(outdir, 'benchmarks/resource_usage.txt'), 'w') as f:
        f.write(report)
    

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
    snakefile = os.path.join(thisdir, '../Snakefile')

    # Load snakemake settings from config if available
    snakemake_opts = {}
    if configfile and os.path.exists(configfile):
        try:
            with open(configfile, 'r') as f:
                config = yaml.safe_load(f)
                snakemake_opts = config.get('snakemake', {})
        except Exception as e:
            if verbose:
                print(f"Warning: Could not load snakemake options from config: {e}")

    # Build Snakemake command with configurable options
    cmd = ["snakemake", "-s", snakefile]
    
    # Add configurable flags (with sane defaults)
    if snakemake_opts.get('use_conda', True):
        cmd.append("--use-conda")
    if snakemake_opts.get('keep_going', True):
        cmd.append("-k")
    if snakemake_opts.get('benchmark_extended', True):
        cmd.append("--benchmark-extended")
    if snakemake_opts.get('print_shell_commands', True):
        cmd.append("--printshellcmds")
    if snakemake_opts.get('rerun_incomplete', True):
        cmd.append("--rerun-incomplete")
    if snakemake_opts.get('generate_report', True):
        cmd.append("--report=results/snakemake_report.html")

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

    # run Snakemake
    return_code = 0
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError as e:
        print(f'Error in Snakemake invocation: {e}', file=sys.stderr)
        return_code = e.returncode
    except FileNotFoundError as e:
        print(f'Snakemake not found: {e}', file=sys.stderr)
        return_code = 1

    if return_code == 0:
        # Aggregate benchmarks and logs on success
        try:
            from aggregate import create_execution_report
            create_execution_report(outdir)
        except Exception as e:
            if verbose:
                print(f"Warning: Could not generate execution report: {e}")
        
        # display resource usage
        track_resources(start_time, net_start, outdir, verbose=verbose)
    
    return return_code


# run snakemake plan to preview the pipeline
def run_snakemake_plan(configfile):
    """Preview the Snakemake plan before running the pipeline."""
    
    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, '../Snakefile')  # Updated path
    
    os.system("snakemake -s " + snakefile + " --use-conda --dag \
        --configfile " + configfile + " --quiet \
        │ dot -Tpng > results/dag.png")
    os.system("snakemake -s " + snakefile + " --use-conda --rulegraph \
        --configfile " + configfile + " --quiet \
        │ dot -Tpng > results/rulegraph.png")


# generate snakemake report for the pipeline
def get_snakemake_report(configfile):
    """Generate a Snakemake report for the pipeline."""
    
    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.join(thisdir, '../Snakefile')  # Updated path
    
    os.system("snakemake -s " 
              + snakefile 
              + " --use-conda"
              + " --report" 
              + " --configfile " 
              + configfile 
              + " --quiet")