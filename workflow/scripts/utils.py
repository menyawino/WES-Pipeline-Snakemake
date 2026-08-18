import sys
import os
import subprocess
import re
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
    """Return the size of a folder in bytes using os.scandir for fast traversal."""
    total_size = 0
    try:
        for entry in os.scandir(folder):
            if entry.is_file(follow_symlinks=False):
                total_size += entry.stat().st_size
            elif entry.is_dir(follow_symlinks=False):
                total_size += get_folder_size(entry.path)
    except (PermissionError, FileNotFoundError):
        pass
    return total_size


# run snakemake with the specified options and configuration
def run_snakemake(configfile, inputdir, outdir, verbose=False, extra_args=[], email=None):
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

    # Resolve Snakemake executable
    import shutil
    snakemake_bin = "snakemake"
    sm_bin = "/home/omar/Downloads/miniconda3/envs/sm/bin/snakemake"
    if not shutil.which("snakemake") and os.path.exists(sm_bin):
        snakemake_bin = sm_bin

    # Initialize RAM disk tmp directory if configured
    tmp_dir = snakemake_opts.get('tmpdir', '/dev/shm/wes_pipeline_tmp')
    try:
        os.makedirs(tmp_dir, exist_ok=True)
    except Exception:
        pass

    # Build Snakemake command with configurable options
    cmd = [snakemake_bin, "-s", snakefile]
    
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

    # Add additional Snakemake arguments
    cmd += list(extra_args)

    if configfile:
        # Only add the specified config file without defaults and system confs
        cmd += ["--configfile", configfile]

    # Pass inputdir and outdir to Snakemake configuration
    config_params = [f"inputdir={inputdir}", f"outdir={outdir}"]
    if email:
        config_params.append(f"email_recipient={email}")
    cmd += ["--config"] + config_params

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


def stylize_dot(dot_text, title="WES Pipeline Workflow"):
    """Transform raw Snakemake DOT output into a modern, publication-ready Graphviz graph."""
    if not dot_text or "digraph" not in dot_text:
        return dot_text

    start_idx = dot_text.find("digraph")
    end_idx = dot_text.rfind("}")
    if start_idx != -1 and end_idx != -1:
        dot_text = dot_text[start_idx:end_idx + 1]

    category_styles = {
        "qc": {"fill": "#F0F9FF", "color": "#0284C7", "text": "#0369A1"},
        "trim": {"fill": "#F0FDF4", "color": "#16A34A", "text": "#15803D"},
        "align": {"fill": "#EEF2FF", "color": "#4F46E5", "text": "#3730A3"},
        "bam_prep": {"fill": "#FAF5FF", "color": "#9333EA", "text": "#6B21A8"},
        "bam_qc": {"fill": "#FDF2F8", "color": "#DB2777", "text": "#9D174D"},
        "variant": {"fill": "#FFFBEB", "color": "#D97706", "text": "#92400E"},
        "annot": {"fill": "#FFEDD5", "color": "#EA580C", "text": "#9A3412"},
        "summary": {"fill": "#F8FAFC", "color": "#475569", "text": "#0F172A"},
    }

    def get_style(label):
        lbl = label.lower()
        if "fastqc" in lbl or "qc_report" in lbl:
            return category_styles["qc"]
        if "trim" in lbl or "fastp" in lbl:
            return category_styles["trim"]
        if "bwa" in lbl or "align" in lbl or "merge" in lbl:
            return category_styles["align"]
        if "markdup" in lbl or "bqsr" in lbl or "recal" in lbl or "filter_bam" in lbl:
            return category_styles["bam_prep"]
        if "flagstat" in lbl or "coverage" in lbl or "depth" in lbl or "metrics" in lbl:
            return category_styles["bam_qc"]
        if "haplotype" in lbl or "genotype" in lbl or "split_vcf" in lbl or "filter_snp" in lbl or "filter_indel" in lbl:
            return category_styles["variant"]
        if "annot" in lbl or "vep" in lbl:
            return category_styles["annot"]
        return category_styles["summary"]

    graph_attrs = f"""
    graph [
        rankdir=LR,
        bgcolor="#FFFFFF",
        pad="0.5",
        nodesep="0.45",
        ranksep="0.65",
        fontname="Inter, Helvetica, Arial, sans-serif",
        dpi=300,
        label="{title}\\n ",
        labelloc="t",
        fontsize=14,
        fontcolor="#1E293B",
        compound=true,
        splines=ortho
    ];
    node [
        shape=box,
        style="filled,rounded",
        fontname="Inter, Helvetica, Arial, sans-serif",
        fontsize=10,
        penwidth=1.5,
        margin="0.2,0.12"
    ];
    edge [
        penwidth=1.5,
        color="#94A3B8",
        arrowsize=0.85,
        fontname="Inter, Helvetica, Arial, sans-serif",
        fontsize=9,
        fontcolor="#64748B"
    ];
    """

    lines = dot_text.splitlines()
    new_lines = []

    for line in lines:
        stripped = line.strip()
        if stripped.startswith("digraph"):
            new_lines.append(line)
            new_lines.append(graph_attrs)
            continue
        if stripped.startswith("graph[") or stripped.startswith("node[") or stripped.startswith("edge["):
            continue

        if "label =" in stripped:
            label_match = re.search(r'label\s*=\s*"([^"]+)"', stripped)
            if label_match:
                label_val = label_match.group(1)
                st = get_style(label_val)
                line = re.sub(
                    r'color\s*=\s*"[^"]+"',
                    f'fillcolor="{st["fill"]}", color="{st["color"]}", fontcolor="{st["text"]}"',
                    line
                )
                if 'fillcolor=' not in line:
                    line = line.replace(']', f', fillcolor="{st["fill"]}", color="{st["color"]}", fontcolor="{st["text"]}"]')
        new_lines.append(line)

    return "\n".join(new_lines)


def render_graph_png(cmd, output_png, title):
    """Execute Snakemake graph command, apply modern DOT styling, and render high-res PNG."""
    try:
        res = subprocess.run(cmd, capture_output=True, text=True, check=True)
        raw_dot = res.stdout
        styled_dot = stylize_dot(raw_dot, title=title)
        
        os.makedirs(os.path.dirname(os.path.abspath(output_png)), exist_ok=True)
        
        # Try running dot from environment or system path
        dot_bin = "dot"
        sm_dot = "/home/omar/Downloads/miniconda3/envs/sm/bin/dot"
        if os.path.exists(sm_dot):
            dot_bin = sm_dot
            
        dot_res = subprocess.run(
            [dot_bin, "-Tpng", "-o", output_png],
            input=styled_dot,
            capture_output=True,
            text=True
        )
        if dot_res.returncode == 0:
            print(f"✓ Rendered aesthetic graph: {output_png}")
        else:
            print(f"Warning: dot rendering returned error: {dot_res.stderr}", file=sys.stderr)
    except Exception as e:
        print(f"Warning: Could not render graph {output_png}: {e}", file=sys.stderr)


# run snakemake plan to preview the pipeline
def run_snakemake_plan(configfile, inputdir=None):
    """Preview the Snakemake plan before running the pipeline."""
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.abspath(os.path.join(thisdir, '../Snakefile'))

    snakemake_bin = "snakemake"
    sm_bin = "/home/omar/Downloads/miniconda3/envs/sm/bin/snakemake"
    if os.path.exists(sm_bin):
        snakemake_bin = sm_bin

    extra_config = []
    if inputdir:
        extra_config = ["--config", f"inputdir={inputdir}"]

    # Generate DAG
    dag_cmd = [snakemake_bin, "-s", snakefile, "--dag", "--configfile", configfile, "--quiet"] + extra_config
    render_graph_png(dag_cmd, "results/dag.png", "Snakemake Rule & Sample Execution DAG")

    # Generate Rulegraph
    rule_cmd = [snakemake_bin, "-s", snakefile, "--rulegraph", "--configfile", configfile, "--quiet"] + extra_config
    render_graph_png(rule_cmd, "results/rulegraph.png", "WES Snakemake Rule Graph Topology")


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