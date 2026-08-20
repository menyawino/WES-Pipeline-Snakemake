import sys
import os
import subprocess
import re
import click
import time
import psutil
import collections
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
# run snakemake with the specified options and configuration
def run_snakemake(
    configfile="workflow/config.yml",
    inputdir=None,
    outdir=None,
    samplesfile="",
    cores=None,
    singularity_args="-B /mnt/bucket -B /mnt/qnap-public",
    rerun_triggers="mtime",
    dry_run=False,
    deepvariant=False,
    verbose=False,
    extra_args=[],
    email=None
):
    """Run Snakemake with the specified options and configuration."""

    # Find the Snakefile relative to the package path
    thisdir = os.path.dirname(__file__)
    snakefile = os.path.abspath(os.path.join(thisdir, '../Snakefile'))

    # Load snakemake settings from config if available
    snakemake_opts = {}
    if configfile and os.path.exists(configfile):
        try:
            with open(configfile, 'r') as f:
                config = yaml.safe_load(f) or {}
                snakemake_opts = config.get('snakemake', {})
        except Exception as e:
            if verbose:
                print(f"Warning: Could not load snakemake options from config: {e}")

    # Resolve Snakemake executable and ensure its environment is in PATH
    import shutil
    snakemake_bin = "snakemake"
    sm_env_bin = "/home/omar/Downloads/miniconda3/envs/sm/bin"
    sm_bin = os.path.join(sm_env_bin, "snakemake")
    if not shutil.which("snakemake") and os.path.exists(sm_bin):
        snakemake_bin = sm_bin
        os.environ["PATH"] = sm_env_bin + os.pathsep + os.environ.get("PATH", "")

    # Initialize RAM disk tmp directory if configured
    tmp_dir = snakemake_opts.get('tmpdir', '/dev/shm/wes_pipeline_tmp')
    try:
        os.makedirs(tmp_dir, exist_ok=True)
    except Exception:
        pass

    # Build Snakemake command with standard flags
    cmd = [snakemake_bin, "-s", snakefile]
    
    # Add conda and singularity flags
    if snakemake_opts.get('use_conda', True):
        cmd.append("--use-conda")
    if snakemake_opts.get('use_singularity', True):
        cmd.append("--use-singularity")
        if singularity_args:
            cmd.extend(["--singularity-args", singularity_args])

    if rerun_triggers:
        cmd.extend(["--rerun-triggers", rerun_triggers])

    if snakemake_opts.get('keep_going', True):
        cmd.append("-k")
    if snakemake_opts.get('print_shell_commands', True):
        cmd.append("--printshellcmds")
    if snakemake_opts.get('benchmark_extended', True):
        cmd.append("--benchmark-extended")
    if snakemake_opts.get('rerun_incomplete', True):
        cmd.append("--rerun-incomplete")

    if dry_run:
        cmd.append("-n")

    # Determine CPU cores
    if not any(arg.startswith("--cores") or arg == "-j" or arg == "-c" for arg in extra_args):
        allocated_cores = cores or snakemake_opts.get('cores', 88)
        cmd.extend(["--cores", str(allocated_cores)])

    # Add config file
    if configfile:
        cmd.extend(["--configfile", configfile])

    # Pass config parameters (inputdir, outdir, samplesfile, email, variant_callers)
    config_params = []
    if inputdir:
        config_params.append(f"inputdir={inputdir}")
    if outdir:
        config_params.append(f"outdir={outdir}")
    if samplesfile is not None:
        config_params.append(f"samplesfile={samplesfile}")
    if email:
        config_params.append(f"email_recipient={email}")
    if deepvariant:
        config_params.append("variant_callers=['gatk','deepvariant']")
    else:
        config_params.append("variant_callers=['gatk']")

    if config_params:
        cmd.append("--config")
        cmd.extend(config_params)

    # Append any extra Snakemake arguments
    cmd += list(extra_args)

    # Print the final command if verbose or info
    print(f"\n{GRE}[INFO] Executing Snakemake command:{NC}")
    print(" ".join(f'"{c}"' if " " in c else c for c in cmd))
    print()

    start_time = time.time()
    net_start = psutil.net_io_counters()  # Capture network usage at start

    # Run Snakemake with live tqdm progress bar or verbose output
    return_code = execute_snakemake_with_progress(
        cmd=cmd,
        outdir=outdir,
        verbose=verbose,
        dry_run=dry_run
    )

    if return_code == 0 and not dry_run:
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


class PipelineDashboard:
    """Rich live dashboard manager for WES Snakemake execution."""

    def __init__(self, outdir, configfile="workflow/config.yml"):
        self.outdir = outdir
        self.configfile = configfile
        self.total_jobs = 0
        self.completed_jobs = 0
        self.start_time = time.time()
        self.spinner_frames = ["⠋", "⠙", "⠹", "⠸", "⠼", "⠴", "⠦", "⠧", "⠇", "⠏"]
        self.spinner_idx = 0
        self.active_jobs = collections.deque(maxlen=3)
        self.setup_status = "Building DAG of jobs & preparing environment..."
        self.setup_substatus = ""

        self.stages = collections.OrderedDict([
            ("01_qc_trim", {
                "name": "01. Quality Control & Trimming",
                "rules": ["trimming_fp", "fastqc", "fastqc_pretrim", "fastqc_posttrim"],
                "total": 0, "done": 0, "status": "pending"
            }),
            ("02_alignment", {
                "name": "02. Alignment & BAM Prep",
                "rules": ["bwa_mem", "merge_bams", "mark_duplicates", "index_markdup_bam"],
                "total": 0, "done": 0, "status": "pending"
            }),
            ("03_bqsr", {
                "name": "03. BQSR & Panel BAM Filtering",
                "rules": ["base_recalibrator", "apply_bqsr", "filter_bam_target", "filter_bam_prot_coding", "filter_bam_canon_tran"],
                "total": 0, "done": 0, "status": "pending"
            }),
            ("04_bam_qc", {
                "name": "04. Coverage & BAM Quality QC",
                "rules": [
                    "flagstat_original", "flagstat_target", "flagstat_prot_coding", "flagstat_canon_tran",
                    "coverage_stats", "coverage_stats_target", "coverage_stats_per_base", "coverage_stats_per_base_target",
                    "coverage_hist", "coverage_hist_target", "fast_bam_qc_prot_coding", "fast_bam_qc_target",
                    "mean_coverage_per_exon", "mean_coverage_per_exon_target", "qc_report"
                ],
                "total": 0, "done": 0, "status": "pending"
            }),
            ("05_variant_calling", {
                "name": "05. Variant Calling (GATK & DeepVariant)",
                "rules": ["split_target_intervals", "haplotypecaller_chunk", "gather_gvcfs", "genotype_gvcfs", "deepvariant_call", "split_vcfs"],
                "total": 0, "done": 0, "status": "pending"
            }),
            ("06_filtering_vep", {
                "name": "06. Variant Filtering & VEP",
                "rules": ["filter_snps", "filter_indels", "vep_genebe_annotate_variants"],
                "total": 0, "done": 0, "status": "pending"
            }),
            ("07_summary_multiqc", {
                "name": "07. ACMG, Reporting & MultiQC",
                "rules": ["aggregate_acmg_annotations", "summarize_variants", "aggregate_variant_summaries", "compare_callers", "multiqc_report", "all"],
                "total": 0, "done": 0, "status": "pending"
            }),
        ])

        self.rule_to_stage = {}
        for s_id, s_info in self.stages.items():
            for r in s_info["rules"]:
                self.rule_to_stage[r] = s_id

    def register_job_stats(self, rule_counts, total_count):
        """Populate initial stage job counts from Snakemake Job stats."""
        self.total_jobs = total_count
        self.setup_status = "Executing workflow tasks..."
        self.setup_substatus = ""
        # Reset totals before populating to avoid double-counting if Job stats prints twice
        for s_info in self.stages.values():
            s_info["total"] = 0
        for r_name, count in rule_counts.items():
            s_id = self.rule_to_stage.get(r_name)
            if s_id and s_id in self.stages:
                self.stages[s_id]["total"] += count

    def notify_job_start(self, rule_name, sample_name=None):
        """Update stage and active job state when a rule starts."""
        s_id = self.rule_to_stage.get(rule_name)
        if s_id and s_id in self.stages:
            if self.stages[s_id]["status"] == "pending":
                self.stages[s_id]["status"] = "running"
        desc = f"[bold cyan]{rule_name}[/]"
        if sample_name:
            desc += f" (sample: [yellow]{sample_name}[/])"
        if desc not in self.active_jobs:
            self.active_jobs.append(desc)

    def notify_job_done(self, rule_name):
        """Update stage and completion counter when a rule finishes."""
        s_id = self.rule_to_stage.get(rule_name)
        if s_id and s_id in self.stages:
            self.stages[s_id]["done"] += 1
            if self.stages[s_id]["total"] > 0 and self.stages[s_id]["done"] >= self.stages[s_id]["total"]:
                self.stages[s_id]["status"] = "done"
            else:
                self.stages[s_id]["status"] = "running"

    def render(self):
        """Render the complete rich dashboard UI components."""
        from rich.console import Group
        from rich.panel import Panel
        from rich.table import Table
        from rich.text import Text

        self.spinner_idx = (self.spinner_idx + 1) % len(self.spinner_frames)
        spin_char = self.spinner_frames[self.spinner_idx]

        elapsed = time.time() - self.start_time
        elapsed_str = str(timedelta(seconds=int(elapsed)))

        pct = (self.completed_jobs / self.total_jobs * 100.0) if self.total_jobs > 0 else 0.0
        rate = (self.completed_jobs / elapsed) if elapsed > 1 else 0.0
        remaining_sec = int((self.total_jobs - self.completed_jobs) / rate) if rate > 0.05 else 0
        eta_str = str(timedelta(seconds=remaining_sec)) if remaining_sec > 0 else "--:--:--"

        # 1. Overall Progress Bar / Pre-flight Header
        if self.total_jobs == 0:
            overall_text = (
                f" [bold cyan]{spin_char} Pre-flight Setup:[/] [yellow]{self.setup_status}[/]"
                + (f"\n [dim]{self.setup_substatus}[/]" if self.setup_substatus else "")
                + f" • Elapsed: [white]{elapsed_str}[/]"
            )
        else:
            bar_len = 36
            filled_len = int(bar_len * (pct / 100.0))
            prog_bar = f"[bold green]{'━' * filled_len}[/][dim]{'━' * (bar_len - filled_len)}[/]"
            overall_text = (
                f" [bold white]Overall Progress[/]  {prog_bar}  "
                f"[bold green]{pct:5.1f}%[/] • [bold cyan]{self.completed_jobs}/{self.total_jobs}[/] Jobs • "
                f"Elapsed: [white]{elapsed_str}[/] • ETA: [yellow]{eta_str}[/] • [cyan]{rate:.2f} j/s[/]"
            )

        # 2. Stage Breakdown Table
        stage_table = Table(box=None, expand=True, padding=(0, 1), show_header=True)
        stage_table.add_column("Pipeline Stage", style="bold white", width=34)
        stage_table.add_column("Status", width=14)
        stage_table.add_column("Progress", justify="right", width=10)
        stage_table.add_column("Jobs", justify="right", width=14)

        for s_id, s_info in self.stages.items():
            tot = s_info["total"]
            done = min(s_info["done"], tot) if tot > 0 else s_info["done"]
            st_pct = (done / tot * 100.0) if tot > 0 else (100.0 if s_info["status"] == "done" else 0.0)

            if s_info["status"] == "done" or (tot > 0 and done >= tot):
                status_str = "[bold green]✔ Complete[/]"
                pct_str = "[green]100%[/]"
            elif s_info["status"] == "running":
                status_str = f"[bold cyan]{spin_char} Running[/]"
                pct_str = f"[yellow]{st_pct:4.0f}%[/]"
            else:
                status_str = "[dim]○ Pending[/]"
                pct_str = "[dim]  0%[/]"

            jobs_str = f"{done}/{tot}" if tot > 0 else "--"
            stage_table.add_row(s_info["name"], status_str, pct_str, jobs_str)

        # 3. System Metrics & Active Tasks Footer
        cpu_usage = psutil.cpu_percent(interval=None)
        mem = psutil.virtual_memory()
        mem_used_gb = mem.used / (1024 ** 3)
        mem_tot_gb = mem.total / (1024 ** 3)

        if self.total_jobs == 0:
            active_display = f"[bold cyan]{spin_char} {self.setup_status}[/]"
            if self.setup_substatus:
                active_display += f" • [yellow]{self.setup_substatus}[/]"
        else:
            active_list = list(self.active_jobs)
            active_display = " • ".join(active_list) if active_list else "[dim]Waiting for next task...[/]"

        footer_panel = Panel(
            f" [bold white]Active Jobs:[/] {active_display}\n"
            f" [dim]System:[/] [cyan]CPU: {cpu_usage:.1f}%[/] • [magenta]RAM: {mem_used_gb:.1f}/{mem_tot_gb:.1f} GB[/] • "
            f"[yellow]Logs:[/] {os.path.join(self.outdir, 'logs/snakemake_execution.log')}",
            title="[bold green]Live Execution Monitor[/]",
            border_style="green",
            padding=(0, 1)
        )

        header_panel = Panel(
            f" [bold green]Output Directory:[/] [white]{self.outdir}[/]\n"
            f"{overall_text}",
            title="[bold green]🧬 WES Analysis Pipeline | Snakemake Execution[/]",
            border_style="green",
            padding=(0, 1)
        )

        return Group(header_panel, stage_table, footer_panel)


def execute_snakemake_with_progress(cmd, outdir, verbose=False, dry_run=False):
    """Execute Snakemake command while displaying a sleek Rich live dashboard and logging full output."""
    from rich.console import Console
    from rich.live import Live

    console = Console()
    logs_dir = os.path.join(outdir, "logs")
    os.makedirs(logs_dir, exist_ok=True)
    log_file_path = os.path.join(logs_dir, "snakemake_execution.log")

    console.print(f"{GRE}[INFO] Full execution log:{NC} {log_file_path}")

    dashboard = PipelineDashboard(outdir=outdir)
    recent_lines = collections.deque(maxlen=60)
    return_code = 0
    in_job_stats = False
    rule_counts = {}
    current_sample = None

    try:
        with open(log_file_path, "w") as log_f:
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                universal_newlines=True
            )

            if verbose:
                for line in proc.stdout:
                    log_f.write(line)
                    log_f.flush()
                    recent_lines.append(line)
                    sys.stdout.write(line)
                    sys.stdout.flush()
                proc.wait()
                return_code = proc.returncode
            elif dry_run:
                print(f"\n{GRE}[INFO] Plan / Dry-Run Job Summary:{NC}")
                for line in proc.stdout:
                    log_f.write(line)
                    log_f.flush()
                    if "Job stats:" in line:
                        in_job_stats = True
                    if in_job_stats:
                        if "Reasons:" in line or "Rules with missing" in line:
                            in_job_stats = False
                        else:
                            sys.stdout.write(line)
                            sys.stdout.flush()
                proc.wait()
                return_code = proc.returncode
            else:
                with Live(dashboard.render(), console=console, refresh_per_second=8, transient=False) as live:
                    for line in proc.stdout:
                        log_f.write(line)
                        log_f.flush()
                        recent_lines.append(line)

                        # Parse Pre-execution setup events (Conda environment creation, package downloads, Singularity pulls)
                        if "Creating conda environment" in line or "Creating environment" in line:
                            env_path = line.strip().split()[-1]
                            dashboard.setup_status = f"Creating conda environment: {env_path}"
                            dashboard.setup_substatus = "Solving dependencies..."
                        elif "Downloading and installing remote packages" in line:
                            dashboard.setup_substatus = "Downloading and installing remote conda packages..."
                        elif "Pulling singularity image" in line:
                            dashboard.setup_status = f"Pulling Singularity container..."
                            dashboard.setup_substatus = line.strip()
                        elif "Building DAG of jobs" in line:
                            dashboard.setup_status = "Building DAG of jobs & calculating execution tree..."
                            dashboard.setup_substatus = ""

                        # Parse Job stats table
                        if "Job stats:" in line:
                            in_job_stats = True
                        elif in_job_stats:
                            if "total" in line:
                                in_job_stats = False
                                m_tot = re.search(r'total\s+(\d+)', line)
                                if m_tot:
                                    tot_count = int(m_tot.group(1))
                                    dashboard.register_job_stats(rule_counts, tot_count)
                            else:
                                m_cnt = re.match(r'^\s*([a-zA-Z0-9_]+)\s+(\d+)\s*$', line)
                                if m_cnt:
                                    rule_counts[m_cnt.group(1)] = int(m_cnt.group(2))

                        # Parse active sample from context
                        m_samp = re.search(r'(?:sample[ =]|for sample\s+)([a-zA-Z0-9_\/.-]+)', line)
                        if m_samp:
                            current_sample = m_samp.group(1).split('/')[-1]

                        # Parse active rule starting
                        m_rule = re.search(r'\b(?:Rule|rule):\s*([a-zA-Z0-9_]+)', line)
                        if not m_rule:
                            m_rule = re.search(r'\bStarting\s+job\s+\d+:\s*([a-zA-Z0-9_]+)', line)

                        if m_rule:
                            r_name = m_rule.group(1)
                            if r_name not in ["stats", "all", "reasons", "reasons:"]:
                                dashboard.notify_job_start(r_name, current_sample)

                        # Parse step completion
                        m_steps = re.search(r'(\d+)\s+of\s+(\d+)\s+steps', line)
                        if m_steps:
                            completed = int(m_steps.group(1))
                            total = int(m_steps.group(2))
                            if dashboard.total_jobs == 0:
                                dashboard.total_jobs = total
                            dashboard.completed_jobs = completed

                        # Parse finished rule
                        m_fin = re.search(r'Finished job \d+:\s*([a-zA-Z0-9_]+)', line)
                        if m_fin:
                            dashboard.notify_job_done(m_fin.group(1))

                        live.update(dashboard.render())

                    proc.wait()
                    return_code = proc.returncode
                    if return_code == 0:
                        dashboard.completed_jobs = dashboard.total_jobs
                        for s_info in dashboard.stages.values():
                            s_info["status"] = "done"
                            if s_info["total"] > 0:
                                s_info["done"] = s_info["total"]
                    live.update(dashboard.render())

    except KeyboardInterrupt:
        console.print(f"\n[bold yellow][INTERRUPTED] Pipeline execution paused/aborted by user.{NC}")
        return 130
    except Exception as e:
        console.print(f"\n[bold red][ERROR] Exception while running Snakemake: {e}{NC}", file=sys.stderr)
        return 1

    if return_code != 0 and not dry_run:
        from rich.panel import Panel
        err_text = "".join(list(recent_lines)[-30:])
        console.print(Panel(
            err_text,
            title=f"[bold red]Pipeline Error (Exit Code: {return_code})[/]",
            subtitle=f"Full log: {log_file_path}",
            border_style="red"
        ))

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
              + " --conda-frontend mamba"
              + " --report" 
              + " --configfile " 
              + configfile 
              + " --quiet")