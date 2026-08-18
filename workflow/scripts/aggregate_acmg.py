#!/usr/bin/env python3
"""
Cohort ACMG Summary & Clinical Dashboard Generator
Aggregates per-sample ACMG classification tables into cohort-level summaries and interactive HTML dashboard.
"""

import sys
import os
import argparse
import glob
import pandas as pd
import json

def aggregate_acmg(input_tsvs, out_report_md, out_summary_tsv, out_summary_json, out_dashboard_html):
    records = []
    pathogenic_variants = []
    
    for tsv_path in input_tsvs:
        if not os.path.exists(tsv_path):
            continue
        try:
            df = pd.read_csv(tsv_path, sep='\t')
            if df.empty:
                continue
            sample = os.path.basename(os.path.dirname(tsv_path))
            if not sample or sample == '009_annotation':
                sample = os.path.basename(tsv_path).replace('.acmg_variants.tsv', '')
                
            counts = df['ACMG_CLASS'].value_counts().to_dict()
            rec = {
                'sample': sample,
                'total_variants': len(df),
                'pathogenic': counts.get('Pathogenic', 0),
                'likely_pathogenic': counts.get('Likely Pathogenic', 0),
                'vus': counts.get('Uncertain Significance (VUS)', 0),
                'likely_benign': counts.get('Likely Benign', 0),
                'benign': counts.get('Benign', 0)
            }
            records.append(rec)
            
            # Extract any Pathogenic / Likely Pathogenic / High-impact variants
            actionable = df[df['ACMG_CLASS'].isin(['Pathogenic', 'Likely Pathogenic'])]
            for _, r in actionable.iterrows():
                p_var = r.to_dict()
                p_var['sample'] = sample
                pathogenic_variants.append(p_var)
        except Exception:
            continue
            
    df_summary = pd.DataFrame(records)
    if df_summary.empty:
        df_summary = pd.DataFrame(columns=['sample', 'total_variants', 'pathogenic', 'likely_pathogenic', 'vus', 'likely_benign', 'benign'])
        
    # Write TSV
    os.makedirs(os.path.dirname(os.path.abspath(out_summary_tsv)), exist_ok=True)
    df_summary.to_csv(out_summary_tsv, sep='\t', index=False)
    
    # Write JSON
    summary_data = {
        'cohort_size': len(df_summary),
        'total_variants': int(df_summary['total_variants'].sum()) if not df_summary.empty else 0,
        'total_pathogenic': int(df_summary['pathogenic'].sum()) if not df_summary.empty else 0,
        'total_likely_pathogenic': int(df_summary['likely_pathogenic'].sum()) if not df_summary.empty else 0,
        'total_vus': int(df_summary['vus'].sum()) if not df_summary.empty else 0,
        'samples': df_summary.to_dict(orient='records'),
        'actionable_variants': pathogenic_variants
    }
    with open(out_summary_json, 'w') as f:
        json.dump(summary_data, f, indent=2)
        
    # Write Markdown Report
    with open(out_report_md, 'w') as f:
        f.write("# Cohort ACMG/AMP Clinical Variant Annotation Report\n\n")
        f.write(f"- **Cohort Size**: {summary_data['cohort_size']} samples\n")
        f.write(f"- **Total Variants Classified**: {summary_data['total_variants']}\n")
        f.write(f"- **Pathogenic Variants**: {summary_data['total_pathogenic']}\n")
        f.write(f"- **Likely Pathogenic Variants**: {summary_data['total_likely_pathogenic']}\n")
        f.write(f"- **Variants of Uncertain Significance (VUS)**: {summary_data['total_vus']}\n\n")
        f.write("## Sample-Level ACMG Breakdown\n\n")
        if not df_summary.empty:
            f.write(df_summary.to_markdown(index=False))
        else:
            f.write("No sample records available.\n")
        f.write("\n")
        
    # Write Dashboard HTML
    write_cohort_dashboard(df_summary, pathogenic_variants, out_dashboard_html)

def write_cohort_dashboard(df_summary, pathogenic_variants, out_html):
    os.makedirs(os.path.dirname(os.path.abspath(out_html)), exist_ok=True)
    
    total_samples = len(df_summary)
    tot_vars = df_summary['total_variants'].sum() if not df_summary.empty else 0
    tot_p = df_summary['pathogenic'].sum() if not df_summary.empty else 0
    tot_lp = df_summary['likely_pathogenic'].sum() if not df_summary.empty else 0
    tot_vus = df_summary['vus'].sum() if not df_summary.empty else 0
    
    sample_rows = []
    if not df_summary.empty:
        for _, r in df_summary.iterrows():
            sample_rows.append(f"""
            <tr>
                <td><strong>{r['sample']}</strong></td>
                <td>{r['total_variants']}</td>
                <td><span style="color:#ef4444;font-weight:700">{r['pathogenic']}</span></td>
                <td><span style="color:#f97316;font-weight:700">{r['likely_pathogenic']}</span></td>
                <td><span style="color:#eab308;font-weight:700">{r['vus']}</span></td>
                <td>{r['likely_benign']}</td>
                <td>{r['benign']}</td>
            </tr>
            """)
    sample_table = "\n".join(sample_rows)
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>Cohort ACMG Variant Classification Dashboard</title>
    <style>
        :root {{
            --bg: #0b1120; --card: #1e293b; --text: #f8fafc; --muted: #94a3b8;
            --pathogenic: #ef4444; --lp: #f97316; --vus: #eab308; --border: #334155;
        }}
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 32px; }}
        .header {{ margin-bottom: 24px; border-bottom: 1px solid var(--border); padding-bottom: 16px; }}
        .header h1 {{ margin: 0 0 8px 0; font-size: 26px; color: #38bdf8; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 16px; margin-bottom: 32px; }}
        .stat-card {{ background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 20px; text-align: center; }}
        .stat-num {{ font-size: 32px; font-weight: 700; margin-top: 4px; }}
        .stat-label {{ font-size: 12px; text-transform: uppercase; color: var(--muted); letter-spacing: 0.5px; }}
        .table-container {{ background: var(--card); border: 1px solid var(--border); border-radius: 8px; overflow-x: auto; margin-bottom: 32px; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 14px; text-align: left; }}
        th, td {{ padding: 12px 18px; border-bottom: 1px solid var(--border); }}
        th {{ background: #0f172a; color: var(--muted); font-size: 11px; text-transform: uppercase; }}
        tr:hover {{ background: rgba(255,255,255,0.02); }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Cohort ACMG/AMP Clinical Variant Dashboard</h1>
        <p>Comprehensive 5-Tier ACMG classification across {total_samples} samples</p>
    </div>
    
    <div class="stats-grid">
        <div class="stat-card"><div class="stat-label">Samples</div><div class="stat-num">{total_samples}</div></div>
        <div class="stat-card"><div class="stat-label">Total Variants</div><div class="stat-num">{tot_vars}</div></div>
        <div class="stat-card"><div class="stat-label" style="color:var(--pathogenic)">Pathogenic</div><div class="stat-num" style="color:var(--pathogenic)">{tot_p}</div></div>
        <div class="stat-card"><div class="stat-label" style="color:var(--lp)">Likely Pathogenic</div><div class="stat-num" style="color:var(--lp)">{tot_lp}</div></div>
        <div class="stat-card"><div class="stat-label" style="color:var(--vus)">VUS</div><div class="stat-num" style="color:var(--vus)">{tot_vus}</div></div>
    </div>
    
    <h2>Sample Variant Classification Summary</h2>
    <div class="table-container">
        <table>
            <thead>
                <tr>
                    <th>Sample</th>
                    <th>Total Variants</th>
                    <th>Pathogenic</th>
                    <th>Likely Pathogenic</th>
                    <th>VUS</th>
                    <th>Likely Benign</th>
                    <th>Benign</th>
                </tr>
            </thead>
            <tbody>
                {sample_table}
            </tbody>
        </table>
    </div>
</body>
</html>
"""
    with open(out_html, 'w') as f:
        f.write(html)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Cohort ACMG Aggregator")
    parser.add_argument("--inputs", nargs="*", default=[], help="Input sample TSVs")
    parser.add_argument("--report-md", required=True, help="Output markdown report")
    parser.add_argument("--summary-tsv", required=True, help="Output summary TSV")
    parser.add_argument("--summary-json", required=True, help="Output summary JSON")
    parser.add_argument("--dashboard-html", required=True, help="Output dashboard HTML")
    
    args = parser.parse_args()
    aggregate_acmg(
        input_tsvs=args.inputs,
        out_report_md=args.report_md,
        out_summary_tsv=args.summary_tsv,
        out_summary_json=args.summary_json,
        out_dashboard_html=args.dashboard_html
    )
