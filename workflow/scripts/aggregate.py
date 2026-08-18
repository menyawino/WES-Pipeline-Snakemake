#!/usr/bin/env python3
"""Pipeline Aggregator and Comprehensive Analysis Report Generator.

Compiles disk space and storage consumption, pipeline execution logs, variant burden statistics,
ACMG clinical classifications, and produces a highly organized, clustered, and interactive
pipeline topology visualizer without runtime benchmark tracking.
"""

from __future__ import annotations

import base64
import glob
import gzip
import html
import json
import logging
import os
import re
import shutil
import subprocess
import sys
import tarfile
from datetime import datetime
from pathlib import Path
from typing import Any

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger("pipeline_aggregator")

PIPELINE_STAGES = {
    "00_ref": {
        "title": "Stage 0: Reference & Indexing",
        "icon": "📦",
        "desc": "Downloads and indexes GRCh38 reference FASTA and interval targets",
        "color": "#38bdf8",
        "rules": ["download_grch38_reference", "split_target_intervals"],
    },
    "01_qc_trim": {
        "title": "Stage 1: Read QC & Adapter Trimming",
        "icon": "🔍",
        "desc": "Pre-trim FastQC, poly-G clipping & quality filtering with fastp, post-trim FastQC",
        "color": "#0284c7",
        "rules": ["raw_fastqc", "trimming_fp", "posttrim_fastqc"],
    },
    "02_alignment": {
        "title": "Stage 2: BWA-MEM2 Alignment & BAM Prep",
        "icon": "🧬",
        "desc": "BWA-MEM2 alignment, coordinate sorting, lane merging, Picard MarkDuplicates, GATK4 BaseRecalibrator & ApplyBQSR",
        "color": "#2563eb",
        "rules": ["bwa_mem", "merge_bams", "mark_duplicates", "index_markdup_bam", "base_recalibrator", "apply_bqsr"],
    },
    "03_bam_qc": {
        "title": "Stage 3: Target BAM Filtering & Coverage QC",
        "icon": "📊",
        "desc": "Target, CDS, and canonical transcript extraction; Picard alignment metrics; exon-level depth of coverage",
        "color": "#0d9488",
        "rules": [
            "filter_bam_target", "filter_bam_prot_coding", "filter_bam_canon_tran",
            "flagstat_original", "flagstat_target", "flagstat_prot_coding", "flagstat_canon_tran",
            "depth_of_coverage", "depth_of_coverage_target", "mean_coverage_per_exon", "mean_coverage_per_exon_target",
            "collect_alignment_summary_metrics", "collect_alignment_summary_metrics_target",
            "coverage_stats", "coverage_stats_target", "coverage_hist", "coverage_hist_target", "qc_report"
        ],
    },
    "04_calling": {
        "title": "Stage 4: GATK4 HaplotypeCaller (Parallel)",
        "icon": "🔬",
        "desc": "Scattered interval chunking, parallel HaplotypeCaller GVCF generation, GatherGVCFs, GenotypeGVCFs",
        "color": "#7c3aed",
        "rules": ["split_bed", "haplotypecaller_chunk", "gather_gvcfs", "genotype_gvcfs", "split_vcfs"],
    },
    "05_filtering": {
        "title": "Stage 5: Hard Variant Quality Filtering",
        "icon": "⚡",
        "desc": "GATK VariantFiltration with standard strict thresholds for SNPs (QD, FS, MQ, SOR) and Indels",
        "color": "#d97706",
        "rules": ["filter_snps", "filter_indels"],
    },
    "06_annotation": {
        "title": "Stage 6: Ensembl VEP & GeneBe ACMG Engine",
        "icon": "🏷️",
        "desc": "Unified single-step Ensembl REST VEP annotation and automated 5-tier ACMG/AMP 2015 evidence grading",
        "color": "#db2777",
        "rules": ["vep_genebe_annotate_variants"],
    },
    "07_reporting": {
        "title": "Stage 7: Cohort Dashboards & MultiQC",
        "icon": "📈",
        "desc": "Sample variant summaries, cohort health reports, ACMG clinical summaries, and full MultiQC aggregation",
        "color": "#059669",
        "rules": ["summarize_variants", "aggregate_acmg_annotations", "aggregate_variant_summaries", "multiqc_report", "all"],
    },
}


def format_bytes(size_bytes: float | int) -> str:
    if size_bytes <= 0:
        return "0 B"
    units = ["B", "KB", "MB", "GB", "TB"]
    i = 0
    size = float(size_bytes)
    while size >= 1024 and i < len(units) - 1:
        size /= 1024.0
        i += 1
    return f"{size:.2f} {units[i]}"


def stylize_dot_graph(dot_str: str, title: str = "WES Pipeline Topology") -> str:
    """Transform raw Snakemake DOT graph into a clustered, modern layout."""
    lines = dot_str.strip().split("\n")
    node_map = {}
    edges = []

    for line in lines:
        m_node = re.search(r'^\s*(\d+)\[label\s*=\s*"([^"]+)"', line)
        if m_node:
            node_id, label = m_node.group(1), m_node.group(2)
            node_map[node_id] = label
        m_edge = re.search(r"^\s*(\d+)\s*->\s*(\d+)", line)
        if m_edge:
            edges.append((m_edge.group(1), m_edge.group(2)))

    new_dot = [
        "digraph WES_Pipeline {",
        f'    graph [rankdir=TB, compound=true, bgcolor="#0b1120", fontname="Helvetica-Bold", fontsize=14, fontcolor="#f8fafc", nodesep=0.55, ranksep=0.7, pad=0.5, label="{title}", labelloc=t];',
        '    node [shape=box, style="filled,rounded", fontname="Helvetica-Bold", fontsize=11, fontcolor="#f8fafc", penwidth=1.8, margin="0.25,0.12"];',
        '    edge [color="#475569", penwidth=1.6, arrowhead=vee, arrowsize=0.85];',
        "",
    ]

    for stage_id, info in PIPELINE_STAGES.items():
        stage_nodes = [nid for nid, lbl in node_map.items() if lbl in info["rules"]]
        if stage_nodes:
            new_dot.append(f"    subgraph cluster_{stage_id} {{")
            new_dot.append(f'        label="{info["title"]}";')
            new_dot.append('        style="rounded,dashed";')
            new_dot.append(f'        color="{info["color"]}";')
            new_dot.append(f'        fontcolor="{info["color"]}";')
            new_dot.append('        bgcolor="#0f172a";')
            new_dot.append("        penwidth=1.5;")
            for nid in stage_nodes:
                lbl = node_map[nid]
                new_dot.append(f'        node_{nid} [label="{lbl}", fillcolor="#1e293b", color="{info["color"]}"];')
            new_dot.append("    }\n")

    all_clustered = {nid for s in PIPELINE_STAGES.values() for nid, lbl in node_map.items() if lbl in s["rules"]}
    for nid, lbl in node_map.items():
        if nid not in all_clustered:
            new_dot.append(f'    node_{nid} [label="{lbl}", fillcolor="#1e293b", color="#38bdf8"];')

    for src, dst in edges:
        new_dot.append(f"    node_{src} -> node_{dst};")

    new_dot.append("}")
    return "\n".join(new_dot)


def export_reproducible_snakefile(outdir: str) -> tuple[str, str]:
    results_dir = os.path.join(outdir, "results")
    os.makedirs(results_dir, exist_ok=True)

    this_dir = os.path.dirname(os.path.abspath(__file__))
    workflow_dir = os.path.abspath(os.path.join(this_dir, ".."))
    snakefile_path = os.path.join(workflow_dir, "Snakefile")

    header_injections = [
        "# =============================================================================\n",
        "# CONSOLIDATED STANDALONE WES SNAKEMAKE PIPELINE\n",
        f"# Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n",
        "# Contains all inlined rules and workflow definitions in a single file\n",
        "# =============================================================================\n\n",
        "import sys, os\n",
        'for p in [os.getcwd(), os.path.join(os.getcwd(), "workflow"), getattr(workflow, "basedir", "")]:\n',
        '    if p and os.path.exists(p) and p not in sys.path:\n',
        "        sys.path.insert(0, p)\n\n",
    ]

    consolidated_lines = list(header_injections)

    if os.path.exists(snakefile_path):
        with open(snakefile_path, "r", encoding="utf-8") as f:
            for line in f:
                if line.strip().startswith("include:"):
                    inc_file = line.strip().split("include:")[1].strip().strip("\"'")
                    full_inc_path = os.path.join(workflow_dir, inc_file)
                    if os.path.exists(full_inc_path):
                        consolidated_lines.append(f"\n# >>>>>>>>>> INLINED RULE: {inc_file} >>>>>>>>>>\n")
                        with open(full_inc_path, "r", encoding="utf-8") as inc_f:
                            consolidated_lines.extend(inc_f.readlines())
                        consolidated_lines.append(f"\n# <<<<<<<<<< END INLINED RULE: {inc_file} <<<<<<<<<<\n\n")
                    else:
                        consolidated_lines.append(line)
                else:
                    consolidated_lines.append(line)

    full_content = "".join(consolidated_lines)

    out_snakefile_gz = os.path.join(results_dir, "Snakefile.gz")
    with gzip.open(out_snakefile_gz, "wt", encoding="utf-8") as out_f:
        out_f.write(full_content)

    tarball_path = os.path.join(results_dir, "reproducible_pipeline.tar.gz")
    with tarfile.open(tarball_path, "w:gz") as tar:
        temp_sf = os.path.join(results_dir, "Snakefile.standalone")
        with open(temp_sf, "w", encoding="utf-8") as tf:
            tf.write(full_content)
        tar.add(temp_sf, arcname="Snakefile")
        os.remove(temp_sf)

        for folder_name in ["scripts", "plugins", "envs"]:
            f_path = os.path.join(workflow_dir, folder_name)
            if os.path.exists(f_path):
                tar.add(f_path, arcname=f"workflow/{folder_name}")
        cfg_path = os.path.join(workflow_dir, "config.yml")
        if os.path.exists(cfg_path):
            tar.add(cfg_path, arcname="workflow/config.yml")

    return out_snakefile_gz, tarball_path


def render_pipeline_graphs(outdir: str, configfile: str | None = None) -> tuple[str | None, str | None]:
    results_dir = os.path.join(outdir, "results")
    os.makedirs(results_dir, exist_ok=True)

    dag_png = os.path.join(results_dir, "dag.png")
    rulegraph_png = os.path.join(results_dir, "rulegraph.png")

    this_dir = os.path.dirname(os.path.abspath(__file__))
    snakefile = os.path.abspath(os.path.join(this_dir, "../Snakefile"))

    snakemake_bin = shutil.which("snakemake") or "/home/omar/Downloads/miniconda3/envs/sm/bin/snakemake"
    dot_bin = shutil.which("dot") or "/home/omar/Downloads/miniconda3/envs/sm/bin/dot"

    cfg_args = []
    if configfile and os.path.exists(configfile):
        cfg_args = ["--configfile", configfile]

    try:
        cmd = [snakemake_bin, "-s", snakefile, "--rulegraph", "--quiet"] + cfg_args
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if res.returncode == 0 and res.stdout and "digraph" in res.stdout:
            styled_dot = stylize_dot_graph(res.stdout, title="WES Pipeline Stage Architecture & Rule Topology")
            subprocess.run([dot_bin, "-Tpng", "-o", rulegraph_png], input=styled_dot, text=True, check=True)
            logger.info("✓ Rendered tidy rulegraph: %s", rulegraph_png)
    except Exception as e:
        logger.warning("Could not render rulegraph: %s", e)

    try:
        cmd = [snakemake_bin, "-s", snakefile, "--dag", "--quiet"] + cfg_args
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=60)
        if res.returncode == 0 and res.stdout and "digraph" in res.stdout:
            styled_dag = stylize_dot_graph(res.stdout, title="Sample Execution DAG Topology")
            subprocess.run([dot_bin, "-Tpng", "-o", dag_png], input=styled_dag, text=True, check=True)
            logger.info("✓ Rendered tidy DAG: %s", dag_png)
    except Exception as e:
        logger.warning("Could not render DAG: %s", e)

    return (
        dag_png if os.path.exists(dag_png) else None,
        rulegraph_png if os.path.exists(rulegraph_png) else None,
    )


def encode_image_base64(path: str | None) -> str | None:
    if not path or not os.path.exists(path):
        return None
    try:
        with open(path, "rb") as f:
            encoded = base64.b64encode(f.read()).decode("utf-8")
            return f"data:image/png;base64,{encoded}"
    except Exception:
        return None


def get_dir_space_breakdown(outdir: str) -> dict[str, Any]:
    breakdown = {"total_bytes": 0, "stages": [], "largest_files": []}
    all_files = []
    analysis_dir = os.path.join(outdir, "analysis")
    stage_dirs = {}

    if os.path.exists(analysis_dir):
        for entry in os.scandir(analysis_dir):
            if entry.is_dir():
                stage_dirs[entry.name] = entry.path

    for top_sub in ["benchmarks", "logs", "results"]:
        p = os.path.join(outdir, top_sub)
        if os.path.exists(p):
            stage_dirs[top_sub] = p

    for stage_name, stage_path in sorted(stage_dirs.items()):
        stage_bytes = 0
        file_count = 0
        for root, _, files in os.walk(stage_path):
            for f in files:
                f_path = os.path.join(root, f)
                try:
                    f_size = os.path.getsize(f_path)
                    stage_bytes += f_size
                    file_count += 1
                    all_files.append((f_path, f_size, stage_name))
                except OSError:
                    pass
        breakdown["stages"].append(
            {
                "stage": stage_name,
                "bytes": stage_bytes,
                "formatted": format_bytes(stage_bytes),
                "file_count": file_count,
            }
        )
        breakdown["total_bytes"] += stage_bytes

    all_files.sort(key=lambda x: x[1], reverse=True)
    for path, sz, stg in all_files[:15]:
        rel_path = os.path.relpath(path, outdir)
        breakdown["largest_files"].append(
            {
                "path": rel_path,
                "size_bytes": sz,
                "formatted": format_bytes(sz),
                "stage": stg,
            }
        )

    return breakdown


def aggregate_logs(log_dir: str) -> dict[str, Any]:
    if not os.path.isdir(log_dir):
        return {"total_logs": 0, "logs_with_errors": 0, "logs_with_warnings": 0, "errors": [], "warnings": []}

    log_files = glob.glob(os.path.join(log_dir, "**/*.log"), recursive=True)
    summary = {
        "total_logs": len(log_files),
        "logs_with_errors": 0,
        "logs_with_warnings": 0,
        "errors": [],
        "warnings": [],
        "timestamp": datetime.now().isoformat(),
    }

    for log_file in log_files:
        try:
            has_err = False
            has_warn = False
            rel = os.path.relpath(log_file, log_dir)
            with open(log_file, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    low = line.lower()
                    if not has_err and any(x in low for x in ("error:", "exception:", "fatal:", "failed:")):
                        has_err = True
                    if not has_warn and "warning:" in low:
                        has_warn = True
                    if has_err and has_warn:
                        break
            if has_err:
                summary["logs_with_errors"] += 1
                summary["errors"].append(rel)
            if has_warn:
                summary["logs_with_warnings"] += 1
                summary["warnings"].append(rel)
        except Exception as e:
            logger.warning("Could not read log file %s: %s", log_file, e)

    return summary


def create_analysis_report(outdir: str, configfile: str | None = None) -> str:
    """Generate comprehensive, clean analysis report HTML without runtime metrics."""
    logger.info("Compiling clean analysis report for %s...", outdir)

    snakefile_gz_path, tarball_path = export_reproducible_snakefile(outdir)
    dag_path, rulegraph_path = render_pipeline_graphs(outdir, configfile)
    dag_b64 = encode_image_base64(dag_path)
    rulegraph_b64 = encode_image_base64(rulegraph_path)

    logs_data = aggregate_logs(os.path.join(outdir, "logs"))
    space_data = get_dir_space_breakdown(outdir)

    var_summary_json = os.path.join(outdir, "analysis/010_summary/cohort_variant_summary.json")
    acmg_summary_json = os.path.join(outdir, "analysis/007_annotation/cohort_acmg_summary.json")

    var_summary = {}
    if os.path.exists(var_summary_json):
        try:
            with open(var_summary_json, "r", encoding="utf-8") as f:
                var_summary = json.load(f)
        except Exception:
            pass

    acmg_summary = {}
    if os.path.exists(acmg_summary_json):
        try:
            with open(acmg_summary_json, "r", encoding="utf-8") as f:
                acmg_summary = json.load(f)
        except Exception:
            pass

    # Build Interactive Stage Cards HTML
    stage_flow_cards = []
    for sid, sinfo in PIPELINE_STAGES.items():
        rules_badges = "".join([f'<span class="rule-badge">{r}</span>' for r in sinfo["rules"]])
        stage_flow_cards.append(
            f"""
            <div class="stage-card" style="border-top: 3px solid {sinfo['color']};">
                <div class="stage-header">
                    <span class="stage-icon">{sinfo['icon']}</span>
                    <div>
                        <h4 style="margin:0;font-size:14px;color:{sinfo['color']};">{sinfo['title']}</h4>
                        <p style="margin:2px 0 0;font-size:12px;color:var(--muted);">{sinfo['desc']}</p>
                    </div>
                </div>
                <div class="stage-rules">
                    {rules_badges}
                </div>
            </div>
            """
        )
    stage_flow_html = "".join(stage_flow_cards)

    # Tables
    stage_rows = []
    tot_bytes = max(space_data["total_bytes"], 1)
    for st in space_data["stages"]:
        pct = (st["bytes"] / tot_bytes) * 100
        stage_rows.append(
            f"""
            <tr>
                <td><strong>{html.escape(st['stage'])}</strong></td>
                <td>{st['file_count']} files</td>
                <td><strong>{st['formatted']}</strong></td>
                <td>
                    <div style="background:#334155;border-radius:4px;overflow:hidden;height:10px;width:100%;">
                        <div style="background:#38bdf8;height:100%;width:{pct:.1f}%;"></div>
                    </div>
                    <small style="color:#94a3b8;">{pct:.1f}%</small>
                </td>
            </tr>
            """
        )

    largest_file_rows = []
    for lf in space_data["largest_files"]:
        largest_file_rows.append(
            f"""
            <tr>
                <td><code>{html.escape(lf['path'])}</code></td>
                <td><span class="badge" style="background:#1e293b;border:1px solid #475569;">{html.escape(lf['stage'])}</span></td>
                <td><strong>{lf['formatted']}</strong></td>
            </tr>
            """
        )

    total_vars = var_summary.get("total_variants", 0)
    pass_vars = var_summary.get("pass_variants", 0)
    annot_vars = var_summary.get("annotated_variants", 0)
    high_impact = var_summary.get("high_impact_variants", 0)
    titv = var_summary.get("mean_titv_ratio", "N/A")

    p_cnt = acmg_summary.get("total_pathogenic", 0)
    lp_cnt = acmg_summary.get("total_likely_pathogenic", 0)
    vus_cnt = acmg_summary.get("total_vus", 0)
    lb_cnt = acmg_summary.get("total_likely_benign", 0)
    b_cnt = acmg_summary.get("total_benign", 0)

    top_gene_items = []
    for g in var_summary.get("top_genes", [])[:10]:
        top_gene_items.append(f"<li><span>{html.escape(g['gene'])}</span><strong>{g['count']}</strong></li>")
    top_genes_html = "".join(top_gene_items) or "<li>No gene hotspots found.</li>"

    sf_size_str = format_bytes(os.path.getsize(snakefile_gz_path)) if os.path.exists(snakefile_gz_path) else "N/A"
    tar_size_str = format_bytes(os.path.getsize(tarball_path)) if os.path.exists(tarball_path) else "N/A"

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>WES Pipeline Execution & Analysis Report</title>
    <style>
        :root {{
            --bg: #0b1120;
            --panel: #1e293b;
            --panel-header: #0f172a;
            --border: #334155;
            --text: #f8fafc;
            --muted: #94a3b8;
            --cyan: #38bdf8;
            --purple: #c084fc;
            --green: #4ade80;
            --yellow: #facc15;
            --red: #f87171;
        }}
        * {{ box-sizing: border-box; }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background: var(--bg);
            color: var(--text);
            margin: 0;
            padding: 0;
            line-height: 1.5;
        }}
        .container {{
            max-width: 1440px;
            margin: 0 auto;
            padding: 36px 24px 64px;
        }}
        .header {{
            background: linear-gradient(135deg, #1e293b 0%, #0f172a 100%);
            border: 1px solid var(--border);
            border-radius: 12px;
            padding: 28px 32px;
            margin-bottom: 28px;
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 20px;
            box-shadow: 0 10px 25px rgba(0,0,0,0.25);
        }}
        .header h1 {{ margin: 0 0 8px 0; font-size: 26px; color: var(--cyan); letter-spacing: -0.5px; }}
        .header p {{ margin: 0; color: var(--muted); font-size: 14px; }}
        .badge-status {{
            background: rgba(74, 222, 128, 0.15);
            border: 1px solid var(--green);
            color: var(--green);
            padding: 6px 14px;
            border-radius: 20px;
            font-weight: 700;
            font-size: 13px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        .grid-4 {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 16px;
            margin-bottom: 28px;
        }}
        .card {{
            background: var(--panel);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 20px;
            position: relative;
            overflow: hidden;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        }}
        .card-label {{ font-size: 12px; font-weight: 600; text-transform: uppercase; color: var(--muted); letter-spacing: 0.8px; margin-bottom: 6px; }}
        .card-val {{ font-size: 30px; font-weight: 800; line-height: 1.1; }}
        .card-sub {{ font-size: 12px; color: var(--muted); margin-top: 6px; }}
        
        .section {{
            background: var(--panel);
            border: 1px solid var(--border);
            border-radius: 10px;
            padding: 24px;
            margin-bottom: 28px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
        }}
        .section h2 {{
            margin: 0 0 18px 0;
            font-size: 18px;
            color: var(--cyan);
            border-bottom: 1px solid var(--border);
            padding-bottom: 12px;
            display: flex;
            justify-content: space-between;
            align-items: center;
        }}
        
        .stage-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fill, minmax(310px, 1fr));
            gap: 16px;
            margin-bottom: 20px;
        }}
        .stage-card {{
            background: #0f172a;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 16px;
            display: flex;
            flex-direction: column;
            justify-content: space-between;
        }}
        .stage-header {{
            display: flex;
            align-items: flex-start;
            gap: 12px;
            margin-bottom: 12px;
        }}
        .stage-icon {{ font-size: 20px; }}
        .stage-rules {{
            display: flex;
            flex-wrap: wrap;
            gap: 6px;
            margin-top: 8px;
        }}
        .rule-badge {{
            background: #1e293b;
            border: 1px solid #334155;
            color: #cbd5e1;
            padding: 3px 8px;
            border-radius: 4px;
            font-size: 11px;
            font-family: monospace;
        }}

        .split-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 24px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            font-size: 13px;
            text-align: left;
        }}
        th, td {{
            padding: 10px 14px;
            border-bottom: 1px solid var(--border);
        }}
        th {{
            background: var(--panel-header);
            color: var(--muted);
            font-size: 11px;
            text-transform: uppercase;
            letter-spacing: 0.5px;
        }}
        tr:hover {{ background: rgba(255,255,255,0.02); }}
        
        .tab-bar {{
            display: flex;
            gap: 10px;
            margin-bottom: 16px;
            border-bottom: 1px solid var(--border);
            padding-bottom: 10px;
        }}
        .tab-btn {{
            background: #0f172a;
            border: 1px solid var(--border);
            color: var(--muted);
            padding: 8px 16px;
            border-radius: 6px;
            cursor: pointer;
            font-weight: 600;
            font-size: 13px;
            transition: all 0.2s;
        }}
        .tab-btn.active {{
            background: #0284c7;
            color: #fff;
            border-color: #0284c7;
        }}
        
        .image-viewer {{
            background: #0f172a;
            border: 1px solid var(--border);
            border-radius: 8px;
            padding: 20px;
            text-align: center;
            overflow: auto;
            max-height: 700px;
        }}
        .image-viewer img {{
            max-width: 100%;
            height: auto;
            border-radius: 4px;
            box-shadow: 0 4px 15px rgba(0,0,0,0.5);
        }}
        .list-box ul {{
            list-style: none;
            padding: 0;
            margin: 0;
        }}
        .list-box li {{
            display: flex;
            justify-content: space-between;
            padding: 9px 0;
            border-bottom: 1px solid var(--border);
            font-size: 13px;
        }}
        .badge {{
            display: inline-block;
            padding: 3px 8px;
            border-radius: 4px;
            font-size: 11px;
            font-weight: 600;
        }}
        .btn-link {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
            background: #0284c7;
            color: #fff;
            padding: 7px 14px;
            border-radius: 6px;
            text-decoration: none;
            font-size: 12px;
            font-weight: 600;
        }}
        .btn-link:hover {{ background: #0369a1; }}
        @media (max-width: 960px) {{
            .split-grid {{ grid-template-columns: 1fr; }}
        }}
    </style>
</head>
<body>
    <div class="container">
        <!-- Header -->
        <header class="header">
            <div>
                <h1>WES Pipeline Execution & Analysis Report</h1>
                <p>Output Directory: <code>{html.escape(outdir)}</code> | Generated: <strong>{datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</strong></p>
            </div>
            <div style="display:flex;gap:12px;align-items:center;">
                <a href="multiqc_report.html" class="btn-link" target="_blank">📊 Open MultiQC Report</a>
                <span class="badge-status">Workflow Complete</span>
            </div>
        </header>

        <!-- KPI Metrics -->
        <div class="grid-4">
            <div class="card">
                <div class="card-label">Cohort & Log Audit</div>
                <div class="card-val" style="color:var(--cyan);">{logs_data['total_logs']} <span style="font-size:16px;color:var(--muted)">Logs</span></div>
                <div class="card-sub">Errors: {logs_data['logs_with_errors']} | Warnings: {logs_data['logs_with_warnings']}</div>
            </div>
            <div class="card">
                <div class="card-label">Total Storage Footprint</div>
                <div class="card-val" style="color:var(--purple);">{format_bytes(space_data['total_bytes'])}</div>
                <div class="card-sub">{len(space_data['stages'])} pipeline directories</div>
            </div>
            <div class="card">
                <div class="card-label">Total Variants (PASS)</div>
                <div class="card-val" style="color:var(--green);">{total_vars} <span style="font-size:16px;color:var(--muted)">({pass_vars} PASS)</span></div>
                <div class="card-sub">Ti/Tv: {titv} | Annotated: {annot_vars}</div>
            </div>
            <div class="card">
                <div class="card-label">ACMG Clinical Signals</div>
                <div class="card-val" style="color:var(--yellow);">{p_cnt + lp_cnt} <span style="font-size:16px;color:var(--muted)">Pathogenic</span></div>
                <div class="card-sub">High Impact: {high_impact} | VUS: {vus_cnt}</div>
            </div>
        </div>

        <!-- Section: Structured Architecture Stages -->
        <div class="section">
            <h2>
                <span>Pipeline Stage Architecture (7 Modular Stages)</span>
                <span style="font-size:12px;color:var(--muted);font-weight:400;">Logical progression from Raw Reads to Clinical Interpretation</span>
            </h2>
            <div class="stage-grid">
                {stage_flow_html}
            </div>
        </div>

        <!-- Section: Tidy Topology Visualizer -->
        <div class="section">
            <h2>Pipeline Topology & Clustered Graphs</h2>
            <div class="tab-bar">
                <button class="tab-btn active" onclick="showGraph('rulegraph')">🏛️ Rule Topology Graph (Clustered)</button>
                <button class="tab-btn" onclick="showGraph('dag')">🔬 Sample Execution DAG</button>
            </div>
            
            <div id="view-rulegraph" class="image-viewer">
                {f'<img src="{rulegraph_b64}" alt="Clustered Rule Graph">' if rulegraph_b64 else '<p style="color:var(--muted);">Rulegraph image not available.</p>'}
            </div>
            <div id="view-dag" class="image-viewer" style="display:none;">
                {f'<img src="{dag_b64}" alt="Sample Execution DAG">' if dag_b64 else '<p style="color:var(--muted);">DAG image not available.</p>'}
            </div>
        </div>

        <!-- Section: Reproducible Pipeline Artifacts -->
        <div class="section">
            <h2>
                <span>Reproducible Pipeline Artifacts</span>
                <span style="font-size:12px;color:var(--muted);font-weight:400;">Single-file executable bundle</span>
            </h2>
            <div class="split-grid">
                <div class="list-box">
                    <ul>
                        <li>
                            <span><strong>Gzipped Standalone Snakefile</strong><br><small style="color:var(--muted)">results/Snakefile.gz (all rules inlined)</small></span>
                            <span class="badge" style="background:#0284c7;color:#fff;">{sf_size_str}</span>
                        </li>
                        <li>
                            <span><strong>Reproducible Pipeline Tarball</strong><br><small style="color:var(--muted)">results/reproducible_pipeline.tar.gz (Snakefile + scripts + plugins)</small></span>
                            <span class="badge" style="background:#7c3aed;color:#fff;">{tar_size_str}</span>
                        </li>
                    </ul>
                </div>
                <div style="background:#0f172a;padding:16px;border-radius:8px;border:1px solid var(--border);font-size:12px;">
                    <p style="margin:0 0 6px;color:var(--cyan);font-weight:600;">How to re-run with standalone Snakefile:</p>
                    <code style="color:#e2e8f0;display:block;white-space:pre-wrap;"># Extract and execute standalone pipeline
gzip -dc results/Snakefile.gz > Snakefile
snakemake --configfile run_config.yml --use-conda -j 16</code>
                </div>
            </div>
        </div>

        <!-- Section: Space & Storage Analysis -->
        <div class="section">
            <h2>Disk Space & Stage Storage Footprint</h2>
            <div class="split-grid">
                <div>
                    <h3 style="font-size:14px;color:var(--muted);margin-top:0;">Storage Consumption by Stage</h3>
                    <table>
                        <thead>
                            <tr>
                                <th>Stage</th>
                                <th>Files</th>
                                <th>Size</th>
                                <th>Share</th>
                            </tr>
                        </thead>
                        <tbody>
                            {''.join(stage_rows) or '<tr><td colspan="4">No storage data.</td></tr>'}
                        </tbody>
                    </table>
                </div>
                <div>
                    <h3 style="font-size:14px;color:var(--muted);margin-top:0;">Largest Generated Files</h3>
                    <table>
                        <thead>
                            <tr>
                                <th>File</th>
                                <th>Stage</th>
                                <th>Size</th>
                            </tr>
                        </thead>
                        <tbody>
                            {''.join(largest_file_rows) or '<tr><td colspan="3">No files found.</td></tr>'}
                        </tbody>
                    </table>
                </div>
            </div>
        </div>

        <!-- Section: Variant Health & Gene Hotspots -->
        <div class="section">
            <h2>Cohort Health & Top Variant Burden</h2>
            <div class="split-grid">
                <div class="list-box">
                    <h3 style="font-size:14px;color:var(--muted);margin-top:0;">Top Variant Burden Genes</h3>
                    <ul>
                        {top_genes_html}
                    </ul>
                </div>
                <div class="list-box">
                    <h3 style="font-size:14px;color:var(--muted);margin-top:0;">ACMG Classification Breakdown</h3>
                    <ul>
                        <li><span>Pathogenic</span><strong style="color:var(--red)">{p_cnt}</strong></li>
                        <li><span>Likely Pathogenic</span><strong style="color:var(--yellow)">{lp_cnt}</strong></li>
                        <li><span>Variant of Uncertain Significance (VUS)</span><strong style="color:var(--purple)">{vus_cnt}</strong></li>
                        <li><span>Likely Benign / Benign</span><strong style="color:var(--green)">{lb_cnt + b_cnt}</strong></li>
                    </ul>
                </div>
            </div>
        </div>
    </div>

    <script>
        function showGraph(type) {{
            const ruleView = document.getElementById('view-rulegraph');
            const dagView = document.getElementById('view-dag');
            const btns = document.querySelectorAll('.tab-btn');
            
            if (type === 'rulegraph') {{
                ruleView.style.display = 'block';
                dagView.style.display = 'none';
                btns[0].classList.add('active');
                btns[1].classList.remove('active');
            }} else {{
                ruleView.style.display = 'none';
                dagView.style.display = 'block';
                btns[0].classList.remove('active');
                btns[1].classList.add('active');
            }}
        }}
    </script>
</body>
</html>
"""

    report_path = os.path.join(outdir, "results", "analysis_report.html")
    with open(report_path, "w", encoding="utf-8") as f:
        f.write(html_content)

    # Clean up duplicate execution_report.html if exists
    exec_path = os.path.join(outdir, "results", "execution_report.html")
    if os.path.exists(exec_path):
        try:
            os.remove(exec_path)
        except OSError:
            pass

    logger.info("✓ Master Analysis Report generated at %s", report_path)
    return report_path


def main():
    if len(sys.argv) < 2:
        print("Usage: aggregate.py <outdir> [configfile]")
        sys.exit(1)

    outdir = sys.argv[1]
    cfg = sys.argv[2] if len(sys.argv) > 2 else None
    create_analysis_report(outdir, cfg)


if __name__ == "__main__":
    main()
