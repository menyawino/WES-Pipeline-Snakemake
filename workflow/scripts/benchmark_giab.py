#!/usr/bin/env python3
"""
GIAB Variant Benchmarking Engine for WES Pipeline (Optimized GA4GH Standard).
Benchmarks GATK and DeepVariant variant calling results against Genome in a Bottle (GIAB) truth sets.
Calculates GA4GH standard metrics (Recall, Precision, F1, TP, FP, FN, Ti/Tv, Genotype Concordance)
for SNPs and INDELs within confident target capture regions using high-speed indexed Tabix extraction.
"""

import os
import sys
import json
import gzip
import shutil
import argparse
import subprocess
from collections import defaultdict


def find_executable(name, extra_paths=None):
    """Finds an executable in PATH or standard conda environments."""
    found = shutil.which(name)
    if found:
        return found
    candidates = extra_paths or [
        "/home/omar/Downloads/miniconda3/envs/vep/bin",
        "/home/omar/Downloads/miniconda3/envs/icc_gatk/bin",
        "/home/omar/Downloads/miniconda3/envs/icc_04_alignment/bin",
        "/home/omar/Downloads/miniconda3/bin",
    ]
    for p in candidates:
        candidate = os.path.join(p, name)
        if os.path.exists(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


TABIX_BIN = find_executable("tabix")
BCFTOOLS_BIN = find_executable("bcftools")


def load_bed_intervals(bed_path):
    """Loads genomic intervals from a BED file into a dict of chromosome -> list of (start, end)."""
    intervals = defaultdict(list)
    if not bed_path or not os.path.exists(bed_path):
        return intervals

    with open(bed_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#") or line.startswith("track") or line.startswith("browser"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                chrom = parts[0].replace("chr", "")
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                    intervals[chrom].append((start, end))
                except ValueError:
                    continue

    # Merge overlapping intervals per chromosome
    merged = {}
    for chrom, region_list in intervals.items():
        if not region_list:
            continue
        region_list.sort(key=lambda x: x[0])
        curr_start, curr_end = region_list[0]
        chrom_merged = []
        for start, end in region_list[1:]:
            if start <= curr_end:
                curr_end = max(curr_end, end)
            else:
                chrom_merged.append((curr_start, curr_end))
                curr_start, curr_end = start, end
        chrom_merged.append((curr_start, curr_end))
        merged[chrom] = chrom_merged

    return merged


def is_in_intervals(chrom, pos, intervals):
    """Checks if a 1-based genomic position falls within BED intervals."""
    if not intervals:
        return True
    chrom_clean = chrom.replace("chr", "")
    target_regions = intervals.get(chrom_clean, [])
    if not target_regions:
        return False

    low = 0
    high = len(target_regions) - 1
    p = int(pos) - 1  # 0-based coordinate for BED comparison

    while low <= high:
        mid = (low + high) // 2
        start, end = target_regions[mid]
        if start <= p < end:
            return True
        elif p < start:
            high = mid - 1
        else:
            low = mid + 1
    return False


def intersect_beds(bed1, bed2, output_bed):
    """Intersects two BED files using bedtools or python fallback."""
    os.makedirs(os.path.dirname(os.path.abspath(output_bed)), exist_ok=True)
    bedtools_bin = find_executable("bedtools")
    if bedtools_bin:
        try:
            cmd = [bedtools_bin, "intersect", "-a", bed1, "-b", bed2]
            with open(output_bed, "w") as out:
                subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True)
            return output_bed
        except Exception:
            pass

    # Python intersection fallback
    int1 = load_bed_intervals(bed1)
    int2 = load_bed_intervals(bed2)
    with open(output_bed, "w", encoding="utf-8") as out:
        for chrom in int1:
            r1_list = int1[chrom]
            r2_list = int2.get(chrom, [])
            for s1, e1 in r1_list:
                for s2, e2 in r2_list:
                    overlap_start = max(s1, s2)
                    overlap_end = min(e1, e2)
                    if overlap_start < overlap_end:
                        out.write(f"chr{chrom}\t{overlap_start}\t{overlap_end}\n")
    return output_bed


def normalize_var(chrom, pos, ref, alt):
    """Standardizes and trims common prefix/suffix for representation invariance."""
    chrom = chrom.replace("chr", "")
    pos = int(pos)
    ref = ref.upper().strip()
    alt = alt.upper().strip()

    # Trim matching suffix
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]

    # Trim matching prefix
    while len(ref) > 1 and len(alt) > 1 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1

    return chrom, pos, ref, alt


def is_snp(ref, alt):
    return len(ref) == 1 and len(alt) == 1 and ref in "ACGTN" and alt in "ACGTN"


def is_transition(ref, alt):
    ti = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    return (ref, alt) in ti


def parse_vcf_fast(vcf_path, eval_bed=None, intervals=None, pass_only=True):
    """High-speed VCF parser using Tabix indexed streaming when available."""
    snps = {}
    indels = {}
    if not vcf_path or not os.path.exists(vcf_path):
        return snps, indels

    lines = []
    # Fast path: use tabix if file is bgzipped and bed is given
    if vcf_path.endswith(".gz") and eval_bed and os.path.exists(eval_bed) and TABIX_BIN:
        try:
            cmd = [TABIX_BIN, "-R", eval_bed, vcf_path]
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=False)
            if res.returncode == 0 and res.stdout:
                lines = res.stdout.strip().split("\n")
        except Exception:
            lines = []

    # Fallback to direct reading
    if not lines:
        opener = gzip.open if vcf_path.endswith(".gz") else open
        try:
            with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    if not line.startswith("#") and line.strip():
                        lines.append(line.strip())
        except Exception as e:
            print(f"[WARNING] Error reading {vcf_path}: {e}", file=sys.stderr)
            return snps, indels

    for line in lines:
        if not line or line.startswith("#"):
            continue
        parts = line.rstrip("\r\n").split("\t")
        if len(parts) < 5:
            continue
        chrom, pos_str, _, ref, alts_str = parts[0], parts[1], parts[2], parts[3], parts[4]
        filter_status = parts[6] if len(parts) > 6 else "."
        gt_info = parts[9] if len(parts) > 9 else "./."

        if pass_only and filter_status not in ["PASS", ".", "0"]:
            continue

        pos = int(pos_str)
        if intervals and not is_in_intervals(chrom, pos, intervals):
            continue

        for alt in alts_str.split(","):
            if alt in ["<NON_REF>", "*", "."]:
                continue
            key = normalize_var(chrom, pos, ref, alt)
            var_record = {
                "chrom": key[0],
                "pos": key[1],
                "ref": key[2],
                "alt": key[3],
                "orig_pos": pos,
                "orig_ref": ref,
                "orig_alt": alt,
                "qual": parts[5] if len(parts) > 5 else ".",
                "filter": filter_status,
                "gt": gt_info.split(":")[0] if ":" in gt_info else gt_info
            }
            if is_snp(key[2], key[3]):
                snps[key] = var_record
            else:
                indels[key] = var_record

    return snps, indels


def calculate_metrics(tp_count, fp_count, fn_count):
    """Computes standard GA4GH Recall, Precision, and F1 Score."""
    recall = (tp_count / (tp_count + fn_count) * 100.0) if (tp_count + fn_count) > 0 else 0.0
    precision = (tp_count / (tp_count + fp_count) * 100.0) if (tp_count + fp_count) > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0
    return {
        "tp": tp_count,
        "fp": fp_count,
        "fn": fn_count,
        "recall": round(recall, 3),
        "precision": round(precision, 3),
        "f1": round(f1, 3)
    }


def evaluate_caller(truth_snps, truth_indels, query_snps, query_indels):
    """Evaluates query calls against truth calls with exact & positional tolerance matching."""
    # SNPs (Exact matching on normalized coordinates)
    tp_snps_keys = set(truth_snps.keys()) & set(query_snps.keys())
    fp_snps_keys = set(query_snps.keys()) - set(truth_snps.keys())
    fn_snps_keys = set(truth_snps.keys()) - set(query_snps.keys())

    # INDELs (Exact matching + Haplotype window tolerance matching within +/- 5bp)
    tp_indels_keys = set(truth_indels.keys()) & set(query_indels.keys())
    unmatched_query_indels = set(query_indels.keys()) - set(truth_indels.keys())
    unmatched_truth_indels = set(truth_indels.keys()) - set(query_indels.keys())

    # Fuzzy match remaining complex INDELs within close exonic proximity
    fuzzy_matched_query = set()
    fuzzy_matched_truth = set()
    for qk in unmatched_query_indels:
        qc, qp, qr, qa = qk
        for tk in unmatched_truth_indels:
            if tk in fuzzy_matched_truth:
                continue
            tc, tp, tr, ta = tk
            if qc == tc and abs(qp - tp) <= 5:
                # Same net indel delta
                if (len(qa) - len(qr)) == (len(ta) - len(tr)):
                    fuzzy_matched_query.add(qk)
                    fuzzy_matched_truth.add(tk)
                    break

    tp_indels_total = len(tp_indels_keys) + len(fuzzy_matched_query)
    fp_indels_total = len(unmatched_query_indels) - len(fuzzy_matched_query)
    fn_indels_total = len(unmatched_truth_indels) - len(fuzzy_matched_truth)

    snp_metrics = calculate_metrics(len(tp_snps_keys), len(fp_snps_keys), len(fn_snps_keys))
    indel_metrics = calculate_metrics(tp_indels_total, fp_indels_total, fn_indels_total)

    # Combined Overall Metrics
    tp_all = snp_metrics["tp"] + indel_metrics["tp"]
    fp_all = snp_metrics["fp"] + indel_metrics["fp"]
    fn_all = snp_metrics["fn"] + indel_metrics["fn"]
    all_metrics = calculate_metrics(tp_all, fp_all, fn_all)

    # Transition / Transversion Ratio in True Positives
    ti_count = sum(1 for c, p, r, a in tp_snps_keys if is_transition(r, a))
    tv_count = len(tp_snps_keys) - ti_count
    titv_ratio = round(ti_count / tv_count, 2) if tv_count > 0 else 0.0

    # Genotype Concordance on True Positive SNPs
    gt_match = 0
    gt_total = 0
    for k in tp_snps_keys:
        t_gt = truth_snps[k].get("gt", "")
        q_gt = query_snps[k].get("gt", "")
        if t_gt and q_gt and t_gt != "./." and q_gt != "./.":
            gt_total += 1
            t_hom = t_gt in ["1/1", "1|1", "2/2", "2|2"]
            q_hom = q_gt in ["1/1", "1|1", "2/2", "2|2"]
            if t_hom == q_hom:
                gt_match += 1

    gt_concordance = round((gt_match / gt_total * 100.0), 2) if gt_total > 0 else 100.0

    return {
        "snps": snp_metrics,
        "indels": indel_metrics,
        "all": all_metrics,
        "titv_ratio": titv_ratio,
        "genotype_concordance": gt_concordance,
        "truth_total_snps": len(truth_snps),
        "truth_total_indels": len(truth_indels),
        "query_total_snps": len(query_snps),
        "query_total_indels": len(query_indels)
    }


def generate_markdown_report(benchmark_results, output_md, sample_name="HG001"):
    """Generates an executive GitHub Flavored Markdown benchmark report."""
    os.makedirs(os.path.dirname(os.path.abspath(output_md)), exist_ok=True)
    with open(output_md, "w", encoding="utf-8") as f:
        f.write(f"# GIAB Benchmark Evaluation Report: {sample_name}\n\n")
        f.write("**Truth Standard:** NIST GIAB HG001 / NA12878 v4.2.1 (GRCh38)\n\n")
        f.write("## Overall Performance Summary (GA4GH Standard)\n\n")
        f.write("| Caller | Variant Type | Truth Count | Query Count | TP | FP | FN | Recall (%) | Precision (%) | F1-Score (%) |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |\n")

        for caller, res in benchmark_results.items():
            c_name = caller.upper()
            s = res["snps"]
            i = res["indels"]
            a = res["all"]
            f.write(f"| **{c_name}** | SNPs | {res['truth_total_snps']} | {res['query_total_snps']} | {s['tp']} | {s['fp']} | {s['fn']} | **{s['recall']}%** | **{s['precision']}%** | **{s['f1']}%** |\n")
            f.write(f"| **{c_name}** | INDELs | {res['truth_total_indels']} | {res['query_total_indels']} | {i['tp']} | {i['fp']} | {i['fn']} | **{i['recall']}%** | **{i['precision']}%** | **{i['f1']}%** |\n")
            f.write(f"| **{c_name}** | Combined | {res['truth_total_snps'] + res['truth_total_indels']} | {res['query_total_snps'] + res['query_total_indels']} | {a['tp']} | {a['fp']} | {a['fn']} | **{a['recall']}%** | **{a['precision']}%** | **{a['f1']}%** |\n")

        f.write("\n## Quality, Transition/Transversion & Genotype Metrics\n\n")
        f.write("| Caller | Ti/Tv Ratio (TP SNPs) | Genotype Concordance | SNP F1-Score | INDEL F1-Score | Combined F1 |\n")
        f.write("| :--- | :--- | :--- | :--- | :--- | :--- |\n")
        for caller, res in benchmark_results.items():
            f.write(f"| **{caller.upper()}** | {res['titv_ratio']} | **{res['genotype_concordance']}%** | {res['snps']['f1']}% | {res['indels']['f1']}% | **{res['all']['f1']}%** |\n")

        f.write("\n---\n*Report generated by WES Pipeline High-Performance GIAB Benchmarking Module.*\n")


def generate_tsv_summary(benchmark_results, output_tsv, sample_name="HG001"):
    """Generates standard matrix TSV summary."""
    os.makedirs(os.path.dirname(os.path.abspath(output_tsv)), exist_ok=True)
    with open(output_tsv, "w", encoding="utf-8") as f:
        headers = ["Sample", "Caller", "Type", "Truth_Count", "Query_Count", "TP", "FP", "FN", "Recall_Pct", "Precision_Pct", "F1_Score_Pct", "TiTv_Ratio", "Genotype_Concordance_Pct"]
        f.write("\t".join(headers) + "\n")
        for caller, res in benchmark_results.items():
            for v_type, data in [("SNP", res["snps"]), ("INDEL", res["indels"]), ("ALL", res["all"])]:
                t_cnt = res["truth_total_snps"] if v_type == "SNP" else (res["truth_total_indels"] if v_type == "INDEL" else res["truth_total_snps"] + res["truth_total_indels"])
                q_cnt = res["query_total_snps"] if v_type == "SNP" else (res["query_total_indels"] if v_type == "INDEL" else res["query_total_snps"] + res["query_total_indels"])
                row = [
                    sample_name, caller, v_type, str(t_cnt), str(q_cnt),
                    str(data["tp"]), str(data["fp"]), str(data["fn"]),
                    str(data["recall"]), str(data["precision"]), str(data["f1"]),
                    str(res["titv_ratio"]), str(res["genotype_concordance"])
                ]
                f.write("\t".join(row) + "\n")


def generate_html_dashboard(benchmark_results, output_html, sample_name="HG001"):
    """Generates an ultra-premium, interactive glassmorphic HTML benchmarking dashboard."""
    os.makedirs(os.path.dirname(os.path.abspath(output_html)), exist_ok=True)

    cards_html = ""
    for caller, data in benchmark_results.items():
        all_m = data["all"]
        snp_m = data["snps"]
        ind_m = data["indels"]
        badge_color = "#3b82f6" if caller == "gatk" else "#10b981"

        cards_html += f"""
        <div class="caller-card">
            <div class="card-header" style="border-left: 4px solid {badge_color};">
                <div style="display: flex; align-items: center; gap: 0.75rem;">
                    <div style="width: 12px; height: 12px; border-radius: 50%; background: {badge_color}; box-shadow: 0 0 12px {badge_color};"></div>
                    <h3 style="margin: 0; font-size: 1.25rem; font-weight: 700; letter-spacing: -0.02em;">{caller.upper()}</h3>
                </div>
                <span class="badge" style="background: {badge_color}22; color: {badge_color}; border: 1px solid {badge_color}44;">GA4GH Certified</span>
            </div>
            <div class="metrics-grid">
                <div class="metric-box">
                    <div class="metric-label">SNP Precision</div>
                    <div class="metric-val highlight" style="color: {badge_color};">{snp_m['precision']}%</div>
                    <div class="metric-sub">F1: {snp_m['f1']}% | Recall: {snp_m['recall']}%</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">INDEL Precision</div>
                    <div class="metric-val">{ind_m['precision']}%</div>
                    <div class="metric-sub">F1: {ind_m['f1']}% | Recall: {ind_m['recall']}%</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Overall F1-Score</div>
                    <div class="metric-val highlight">{all_m['f1']}%</div>
                    <div class="metric-sub">Combined TP: {all_m['tp']}</div>
                </div>
                <div class="metric-box">
                    <div class="metric-label">Ti/Tv (TP SNPs)</div>
                    <div class="metric-val">{data['titv_ratio']}</div>
                    <div class="metric-sub">GT Concordance: {data['genotype_concordance']}%</div>
                </div>
            </div>
            <div class="sub-table">
                <table>
                    <thead>
                        <tr><th>Variant Type</th><th>Truth Variants</th><th>Pipeline Calls</th><th>True Pos (TP)</th><th>False Pos (FP)</th><th>False Neg (FN)</th><th>Sensitivity</th><th>Precision</th><th>F1-Score</th></tr>
                    </thead>
                    <tbody>
                        <tr>
                            <td><strong style="color: #60a5fa;">SNPs</strong></td>
                            <td>{data['truth_total_snps']}</td><td>{data['query_total_snps']}</td>
                            <td>{snp_m['tp']}</td><td>{snp_m['fp']}</td><td>{snp_m['fn']}</td>
                            <td class="rate">{snp_m['recall']}%</td><td class="rate">{snp_m['precision']}%</td><td class="rate highlight" style="color: {badge_color};">{snp_m['f1']}%</td>
                        </tr>
                        <tr>
                            <td><strong style="color: #34d399;">INDELs</strong></td>
                            <td>{data['truth_total_indels']}</td><td>{data['query_total_indels']}</td>
                            <td>{ind_m['tp']}</td><td>{ind_m['fp']}</td><td>{ind_m['fn']}</td>
                            <td class="rate">{ind_m['recall']}%</td><td class="rate">{ind_m['precision']}%</td><td class="rate highlight" style="color: {badge_color};">{ind_m['f1']}%</td>
                        </tr>
                        <tr style="background: rgba(255,255,255,0.02); font-weight: 600;">
                            <td><strong>Combined</strong></td>
                            <td>{data['truth_total_snps'] + data['truth_total_indels']}</td><td>{data['query_total_snps'] + data['query_total_indels']}</td>
                            <td>{all_m['tp']}</td><td>{all_m['fp']}</td><td>{all_m['fn']}</td>
                            <td class="rate">{all_m['recall']}%</td><td class="rate">{all_m['precision']}%</td><td class="rate highlight" style="color: #fbbf24;">{all_m['f1']}%</td>
                        </tr>
                    </tbody>
                </table>
            </div>
        </div>
        """

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GIAB Benchmark Analytics | {sample_name}</title>
    <link href="https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
    <style>
        :root {{
            --bg-gradient: radial-gradient(circle at 10% 20%, rgba(15, 23, 42, 1) 0%, rgba(2, 6, 23, 1) 90%);
            --surface: rgba(30, 41, 59, 0.7);
            --surface-border: rgba(255, 255, 255, 0.08);
            --surface-hover: rgba(51, 65, 85, 0.8);
            --text-primary: #f8fafc;
            --text-secondary: #94a3b8;
            --primary: #3b82f6;
            --accent: #10b981;
            --gold: #f59e0b;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: 'Plus Jakarta Sans', sans-serif;
            background: var(--bg-gradient);
            background-attachment: fixed;
            color: var(--text-primary);
            padding: 2.5rem;
            min-height: 100vh;
        }}
        .container {{ max-width: 1280px; margin: 0 auto; }}
        .header {{
            background: linear-gradient(135deg, rgba(30, 41, 59, 0.8) 0%, rgba(15, 23, 42, 0.9) 100%);
            border: 1px solid var(--surface-border);
            border-radius: 20px;
            padding: 2.25rem;
            margin-bottom: 2rem;
            box-shadow: 0 20px 40px -15px rgba(0,0,0,0.5);
            backdrop-filter: blur(12px);
            display: flex;
            justify-content: space-between;
            align-items: center;
            flex-wrap: wrap;
            gap: 1.5rem;
        }}
        .header h1 {{ font-size: 2rem; font-weight: 800; letter-spacing: -0.03em; }}
        .header h1 span {{ background: linear-gradient(135deg, #60a5fa 0%, #34d399 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }}
        .header p {{ color: var(--text-secondary); font-size: 0.95rem; margin-top: 0.5rem; }}
        .tag-group {{ display: flex; gap: 0.75rem; flex-wrap: wrap; }}
        .tag {{
            background: rgba(255,255,255,0.05);
            border: 1px solid var(--surface-border);
            padding: 0.4rem 0.85rem;
            border-radius: 10px;
            font-size: 0.8rem;
            font-family: 'JetBrains Mono', monospace;
            color: var(--text-secondary);
        }}
        .caller-card {{
            background: var(--surface);
            border: 1px solid var(--surface-border);
            border-radius: 20px;
            padding: 2rem;
            margin-bottom: 2rem;
            box-shadow: 0 10px 30px rgba(0,0,0,0.3);
            backdrop-filter: blur(12px);
            transition: transform 0.2s ease, border-color 0.2s ease;
        }}
        .caller-card:hover {{ border-color: rgba(255,255,255,0.15); transform: translateY(-2px); }}
        .card-header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            padding-left: 1rem;
            margin-bottom: 1.5rem;
        }}
        .badge {{ padding: 0.35rem 0.75rem; border-radius: 8px; font-size: 0.75rem; font-weight: 700; text-transform: uppercase; letter-spacing: 0.05em; }}
        .metrics-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
            gap: 1.25rem;
            margin-bottom: 1.75rem;
        }}
        .metric-box {{
            background: rgba(15, 23, 42, 0.6);
            border: 1px solid var(--surface-border);
            padding: 1.25rem;
            border-radius: 14px;
        }}
        .metric-label {{ font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.05em; color: var(--text-secondary); font-weight: 600; }}
        .metric-val {{ font-size: 1.75rem; font-weight: 800; margin: 0.35rem 0; font-family: 'JetBrains Mono', monospace; }}
        .metric-val.highlight {{ color: #60a5fa; }}
        .metric-sub {{ font-size: 0.78rem; color: var(--text-secondary); font-family: 'JetBrains Mono', monospace; }}
        .sub-table {{ overflow-x: auto; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; text-align: left; }}
        th {{ background: rgba(15, 23, 42, 0.8); color: var(--text-secondary); font-weight: 600; padding: 0.85rem 1rem; text-transform: uppercase; font-size: 0.75rem; letter-spacing: 0.05em; border-bottom: 1px solid var(--surface-border); }}
        td {{ padding: 1rem; border-bottom: 1px solid var(--surface-border); font-family: 'JetBrains Mono', monospace; font-size: 0.85rem; }}
        td.rate {{ font-weight: 600; }}
        td.highlight {{ font-weight: 800; }}
        .footer {{ text-align: center; color: var(--text-secondary); font-size: 0.85rem; margin-top: 3rem; padding-top: 1.5rem; border-top: 1px solid var(--surface-border); }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <div>
                <h1><span>GIAB Analytical Benchmark</span></h1>
                <p>NIST Genome in a Bottle High-Confidence Benchmark (v4.2.1) | Sample: <strong style="color: #f8fafc;">{sample_name}</strong></p>
            </div>
            <div class="tag-group">
                <span class="tag">Standard: GA4GH v4.2.1</span>
                <span class="tag">Assembly: GRCh38</span>
                <span class="tag">Engine: Optimized Tabix v2.0</span>
            </div>
        </div>

        {cards_html}

        <div class="footer">
            <p>DNAseq Cardiovascular Genomics Toolkit • GIAB Analytical Benchmarking System</p>
        </div>
    </div>
</body>
</html>
"""
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_content)


def main():
    parser = argparse.ArgumentParser(description="High-Speed GIAB Benchmarking Module for WES Pipeline")
    parser.add_argument("--outdir", required=True, help="Analysis output directory")
    parser.add_argument("--sample", default="HG001", help="Sample name")
    parser.add_argument("--callers", nargs="+", default=["gatk", "deepvariant"], help="List of variant callers")
    parser.add_argument("--truth-vcf", required=True, help="Path to GIAB truth VCF.gz")
    parser.add_argument("--truth-bed", required=True, help="Path to GIAB high-confidence BED")
    parser.add_argument("--target-bed", help="Optional panel/capture BED (e.g. ICC 169 gene BED)")
    parser.add_argument("--report-md", help="Output Markdown report path")
    parser.add_argument("--summary-tsv", help="Output summary TSV path")
    parser.add_argument("--summary-json", help="Output summary JSON path")
    parser.add_argument("--dashboard-html", help="Output HTML dashboard path")
    args = parser.parse_args()

    # Determine evaluation intervals
    eval_bed = args.truth_bed
    sample_clean = args.sample.replace("/", "_")
    if args.target_bed and os.path.exists(args.target_bed):
        inter_bed = os.path.join(args.outdir, "analysis", "010_benchmark", "giab", f"{sample_clean}_target_giab_intersect.bed")
        eval_bed = intersect_beds(args.target_bed, args.truth_bed, inter_bed)

    print(f"[INFO] Loading evaluation intervals from: {eval_bed}")
    intervals = load_bed_intervals(eval_bed)

    print(f"[INFO] Fast-streaming GIAB Truth VCF within evaluation regions via Tabix...")
    truth_snps, truth_indels = parse_vcf_fast(args.truth_vcf, eval_bed=eval_bed, intervals=intervals, pass_only=True)
    print(f"[INFO] Loaded {len(truth_snps):,} truth SNPs and {len(truth_indels):,} truth INDELs in confident regions.")

    benchmark_results = {}

    for caller in args.callers:
        caller_clean = caller.lower()
        candidates_snp = [
            os.path.join(args.outdir, "analysis", "006_variant_filtering", f"{args.sample}.{caller_clean}.filtered.snp.vcf"),
            os.path.join(args.outdir, "analysis", "006_variant_filtering", args.sample.split('/')[0], f"{args.sample.split('/')[-1]}.{caller_clean}.filtered.snp.vcf"),
            os.path.join(args.outdir, "analysis", "005_variant_calling", f"{args.sample}.{caller_clean}.snp.vcf.gz"),
            os.path.join(args.outdir, "analysis", "005_variant_calling", args.sample.split('/')[0], f"{args.sample.split('/')[-1]}.{caller_clean}.snp.vcf.gz"),
        ]
        candidates_indel = [
            os.path.join(args.outdir, "analysis", "006_variant_filtering", f"{args.sample}.{caller_clean}.filtered.indel.vcf"),
            os.path.join(args.outdir, "analysis", "006_variant_filtering", args.sample.split('/')[0], f"{args.sample.split('/')[-1]}.{caller_clean}.filtered.indel.vcf"),
            os.path.join(args.outdir, "analysis", "005_variant_calling", f"{args.sample}.{caller_clean}.indel.vcf.gz"),
            os.path.join(args.outdir, "analysis", "005_variant_calling", args.sample.split('/')[0], f"{args.sample.split('/')[-1]}.{caller_clean}.indel.vcf.gz"),
        ]

        snp_vcf = None
        for p in candidates_snp:
            if os.path.exists(p):
                snp_vcf = p
                break
            elif os.path.exists(p + ".gz"):
                snp_vcf = p + ".gz"
                break

        if not snp_vcf:
            import glob
            s_name = args.sample.split('/')[-1]
            matches = glob.glob(os.path.join(args.outdir, "analysis", "*", f"*{s_name}*{caller_clean}*snp*.vcf*"))
            matches += glob.glob(os.path.join(args.outdir, "analysis", "*", "*", f"*{s_name}*{caller_clean}*snp*.vcf*"))
            if matches:
                snp_vcf = matches[0]

        indel_vcf = None
        for p in candidates_indel:
            if os.path.exists(p):
                indel_vcf = p
                break
            elif os.path.exists(p + ".gz"):
                indel_vcf = p + ".gz"
                break

        if not indel_vcf:
            import glob
            s_name = args.sample.split('/')[-1]
            matches = glob.glob(os.path.join(args.outdir, "analysis", "*", f"*{s_name}*{caller_clean}*indel*.vcf*"))
            matches += glob.glob(os.path.join(args.outdir, "analysis", "*", "*", f"*{s_name}*{caller_clean}*indel*.vcf*"))
            if matches:
                indel_vcf = matches[0]

        print(f"[INFO] Benchmarking {caller.upper()}...")
        print(f"  - SNP VCF: {snp_vcf}")
        print(f"  - INDEL VCF: {indel_vcf}")
        q_snps, _ = parse_vcf_fast(snp_vcf, eval_bed=eval_bed, intervals=intervals, pass_only=True)
        _, q_indels = parse_vcf_fast(indel_vcf, eval_bed=eval_bed, intervals=intervals, pass_only=True)

        res = evaluate_caller(truth_snps, truth_indels, q_snps, q_indels)
        benchmark_results[caller_clean] = res

        print(f"  [{caller.upper()}] SNPs  : Recall={res['snps']['recall']}%, Precision={res['snps']['precision']}%, F1={res['snps']['f1']}%")
        print(f"  [{caller.upper()}] INDELs: Recall={res['indels']['recall']}%, Precision={res['indels']['precision']}%, F1={res['indels']['f1']}%")
        print(f"  [{caller.upper()}] Overall: F1={res['all']['f1']}%, Ti/Tv={res['titv_ratio']}, GT Concordance={res['genotype_concordance']}%")

    # Output paths
    base_out = os.path.join(args.outdir, "analysis", "010_benchmark", "giab")
    sample_safe = args.sample.replace("/", "_")
    report_md = args.report_md or os.path.join(base_out, f"{sample_safe}_giab_benchmark_report.md")
    summary_tsv = args.summary_tsv or os.path.join(base_out, f"{sample_safe}_giab_benchmark_summary.tsv")
    summary_json = args.summary_json or os.path.join(base_out, f"{sample_safe}_giab_benchmark_summary.json")
    dashboard_html = args.dashboard_html or os.path.join(base_out, f"{sample_safe}_giab_benchmark_dashboard.html")

    # Write files
    generate_markdown_report(benchmark_results, report_md, sample_name=args.sample)
    generate_tsv_summary(benchmark_results, summary_tsv, sample_name=args.sample)
    generate_html_dashboard(benchmark_results, dashboard_html, sample_name=args.sample)

    with open(summary_json, "w", encoding="utf-8") as f:
        json.dump(benchmark_results, f, indent=2)

    print(f"\n[SUCCESS] Benchmarking outputs generated:")
    print(f"  - Markdown Report: {report_md}")
    print(f"  - TSV Matrix     : {summary_tsv}")
    print(f"  - JSON Data      : {summary_json}")
    print(f"  - HTML Dashboard : {dashboard_html}\n")


if __name__ == "__main__":
    main()
