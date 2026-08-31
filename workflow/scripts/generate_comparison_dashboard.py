#!/usr/bin/env python3
import json
import os

json_path = "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/012_benchmark/giab/pipeline_comparison_benchmark.json"
html_path = "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/012_benchmark/giab/pipeline_comparison_dashboard.html"

with open(json_path) as f:
    data = json.load(f)["callable_regions"]

colors = {
    "Legacy UnifiedGenotyper (GATK 3.2)": "#94a3b8",
    "Legacy HaplotypeCaller (GATK 3.2)": "#a855f7",
    "Snakemake GATK4 HaplotypeCaller": "#3b82f6",
    "Snakemake Google DeepVariant": "#10b981"
}

badges = {
    "Legacy UnifiedGenotyper (GATK 3.2)": "Legacy GATK 3.2 (Positional)",
    "Legacy HaplotypeCaller (GATK 3.2)": "Legacy GATK 3.2 (HMM Assembly)",
    "Snakemake GATK4 HaplotypeCaller": "Snakemake Modern (GATK 4.6)",
    "Snakemake Google DeepVariant": "Snakemake Modern (Neural CNN)"
}

cards_html = ""
for name, res in data.items():
    color = colors.get(name, "#3b82f6")
    badge = badges.get(name, "Caller")
    s = res["snps"]
    i = res["indels"]
    a = res["all"]
    
    cards_html += f"""
    <div class="caller-card" style="border-top: 4px solid {color};">
        <div class="card-header">
            <div>
                <h3 style="color: #f8fafc; font-size: 1.25rem; margin: 0 0 0.25rem 0;">{name}</h3>
                <span class="badge" style="background: {color}22; color: {color}; border: 1px solid {color}44;">{badge}</span>
            </div>
            <div style="text-align: right;">
                <span style="font-size: 0.8rem; color: #94a3b8;">OVERALL F1-SCORE</span>
                <div style="font-size: 1.75rem; font-weight: 800; color: {color}; font-family: monospace;">{a['f1']}%</div>
            </div>
        </div>
        
        <div class="metrics-grid">
            <div class="metric-box">
                <div class="metric-label">SNP Sensitivity (Recall)</div>
                <div class="metric-val" style="color: #60a5fa;">{s['recall']}%</div>
                <div class="metric-sub">TP: {s['tp']} | FN: {s['fn']}</div>
            </div>
            <div class="metric-box">
                <div class="metric-label">SNP Precision (PPV)</div>
                <div class="metric-val">{s['precision']}%</div>
                <div class="metric-sub">FP: {s['fp']} (F1: {s['f1']}%)</div>
            </div>
            <div class="metric-box">
                <div class="metric-label">INDEL Sensitivity (Recall)</div>
                <div class="metric-val" style="color: #34d399;">{i['recall']}%</div>
                <div class="metric-sub">TP: {i['tp']} | FN: {i['fn']}</div>
            </div>
            <div class="metric-box">
                <div class="metric-label">INDEL Precision (PPV)</div>
                <div class="metric-val" style="color: {color};">{i['precision']}%</div>
                <div class="metric-sub">FP: {i['fp']} (F1: {i['f1']}%)</div>
            </div>
        </div>
        
        <table style="width: 100%; border-collapse: collapse; margin-top: 1rem; font-size: 0.85rem; font-family: monospace;">
            <tr style="background: rgba(15,23,42,0.6); color: #94a3b8; text-align: left;">
                <th style="padding: 0.6rem;">TYPE</th><th>TRUTH</th><th>CALLS</th><th>TP</th><th>FP</th><th>FN</th><th>RECALL</th><th>PRECISION</th><th>F1-SCORE</th>
            </tr>
            <tr style="border-bottom: 1px solid rgba(255,255,255,0.05);">
                <td style="padding: 0.6rem; font-weight: bold; color: #60a5fa;">SNPs</td><td>{res['truth_total_snps']}</td><td>{res['query_total_snps']}</td><td>{s['tp']}</td><td>{s['fp']}</td><td>{s['fn']}</td><td>{s['recall']}%</td><td>{s['precision']}%</td><td style="color: #60a5fa; font-weight: bold;">{s['f1']}%</td>
            </tr>
            <tr style="border-bottom: 1px solid rgba(255,255,255,0.05);">
                <td style="padding: 0.6rem; font-weight: bold; color: #34d399;">INDELs</td><td>{res['truth_total_indels']}</td><td>{res['query_total_indels']}</td><td>{i['tp']}</td><td>{i['fp']}</td><td>{i['fn']}</td><td>{i['recall']}%</td><td>{i['precision']}%</td><td style="color: #34d399; font-weight: bold;">{i['f1']}%</td>
            </tr>
            <tr style="background: rgba(255,255,255,0.02); font-weight: bold;">
                <td style="padding: 0.6rem;">COMBINED</td><td>{res['truth_total_snps'] + res['truth_total_indels']}</td><td>{res['query_total_snps'] + res['query_total_indels']}</td><td>{a['tp']}</td><td>{a['fp']}</td><td>{a['fn']}</td><td>{a['recall']}%</td><td>{a['precision']}%</td><td style="color: {color}; font-size: 1rem;">{a['f1']}%</td>
            </tr>
        </table>
    </div>
    """

full_html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Pipeline Head-to-Head Comparison: Snakemake vs Legacy</title>
    <link href="https://fonts.googleapis.com/css2?family=Plus+Jakarta+Sans:wght@300;400;500;600;700;800&family=JetBrains+Mono:wght@400;600&display=swap" rel="stylesheet">
    <style>
        body {{ font-family: "Plus Jakarta Sans", sans-serif; background: #0b0f19; color: #f8fafc; padding: 2rem; margin: 0; }}
        .container {{ max-width: 1300px; margin: 0 auto; }}
        .header {{ background: linear-gradient(135deg, rgba(30, 41, 59, 0.9) 0%, rgba(15, 23, 42, 0.95) 100%); border: 1px solid rgba(255,255,255,0.1); border-radius: 20px; padding: 2rem; margin-bottom: 2rem; box-shadow: 0 10px 30px rgba(0,0,0,0.5); }}
        .header h1 {{ font-size: 2rem; font-weight: 800; margin: 0 0 0.5rem 0; }}
        .header h1 span {{ background: linear-gradient(135deg, #60a5fa 0%, #34d399 100%); -webkit-background-clip: text; -webkit-text-fill-color: transparent; }}
        .caller-card {{ background: rgba(30, 41, 59, 0.7); border: 1px solid rgba(255,255,255,0.08); border-radius: 16px; padding: 1.75rem; margin-bottom: 1.5rem; backdrop-filter: blur(12px); box-shadow: 0 4px 20px rgba(0,0,0,0.25); }}
        .card-header {{ display: flex; justify-content: space-between; align-items: center; margin-bottom: 1.25rem; }}
        .badge {{ padding: 0.35rem 0.75rem; border-radius: 8px; font-size: 0.75rem; font-weight: 700; text-transform: uppercase; }}
        .metrics-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(220px, 1fr)); gap: 1rem; margin-bottom: 1rem; }}
        .metric-box {{ background: rgba(15, 23, 42, 0.6); border: 1px solid rgba(255,255,255,0.06); padding: 1rem; border-radius: 12px; }}
        .metric-label {{ font-size: 0.75rem; text-transform: uppercase; color: #94a3b8; font-weight: 600; }}
        .metric-val {{ font-size: 1.5rem; font-weight: 800; font-family: "JetBrains Mono", monospace; margin: 0.25rem 0; }}
        .metric-sub {{ font-size: 0.75rem; color: #94a3b8; font-family: "JetBrains Mono", monospace; }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1><span>Pipeline Head-to-Head Comparison</span></h1>
            <p style="color: #94a3b8; margin: 0;">NIST GIAB HG001 Standard (GRCh38) | Snakemake Modern Architecture vs. Legacy Parallel Bash Pipeline</p>
        </div>
        {cards_html}
    </div>
</body>
</html>"""

with open(html_path, "w") as f:
    f.write(full_html)

print(f"[SUCCESS] Dashboard written to: {html_path}")
