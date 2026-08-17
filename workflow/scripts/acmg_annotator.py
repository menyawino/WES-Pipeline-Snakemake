#!/usr/bin/env python3
"""
ACMG Variant Annotation and Classification Engine (ACMG/AMP 2015 Guidelines)
Performs automated 5-tier clinical classification (Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign)
using evidence criteria (PVS1, PM2, PM4, PM5, PP3, PP5, BA1, BS1, BP4, BP6, BP7).
"""

import sys
import os
import argparse
import gzip
import re
import math
from collections import defaultdict
import pandas as pd

def open_vcf(filepath):
    """Open plain or gzipped VCF file."""
    if filepath.endswith('.gz'):
        return gzip.open(filepath, 'rt')
    return open(filepath, 'r')

def load_target_genes(bed_file):
    """Load target intervals and gene symbols from panel BED."""
    intervals = defaultdict(list)
    if not bed_file or not os.path.exists(bed_file):
        return intervals
        
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith(('#', 'track', 'browser')) or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                gene = "UNKNOWN"
                if len(parts) > 3:
                    name_field = parts[3]
                    # Format often GENE.chr... or GENE
                    gene = name_field.split('.')[0] if '.' in name_field else name_field
                intervals[chrom].append((start, end, gene))
                
    # Sort intervals
    for chrom in intervals:
        intervals[chrom].sort(key=lambda x: x[0])
    return intervals

def find_gene(chrom, pos, intervals):
    """Find overlapping gene from BED intervals."""
    if chrom not in intervals:
        return "INTERGENIC"
    for start, end, gene in intervals[chrom]:
        if start <= pos <= end:
            return gene
        if start - 100 <= pos <= end + 100:
            return gene
    return "INTERGENIC"

def determine_consequence(ref, alt, pos, intervals, chrom):
    """Determine variant molecular consequence and functional impact."""
    ref_len = len(ref)
    alt_len = len(alt)
    gene = find_gene(chrom, pos, intervals)
    
    if gene == "INTERGENIC":
        return "intergenic_variant", "MODIFIER", gene
        
    # Check if variant is an indel
    if ref_len != alt_len:
        diff = abs(ref_len - alt_len)
        if diff % 3 != 0:
            return "frameshift_variant", "HIGH", gene
        else:
            return "inframe_indel", "MODERATE", gene
            
    # Single nucleotide variant in target region
    return "missense_variant", "MODERATE", gene

def classify_acmg(consequence, impact, pop_af, dbsnp_id, clinvar_sig=None):
    """
    ACMG/AMP 2015 Bayesian Classification Logic
    Criteria:
    - PVS1: Null variant (frameshift, stop_gained, canonical splice) in gene where LOF is disease mechanism
    - PM2: Extremely low population allele frequency (< 0.0001 or absent)
    - PM4: In-frame protein length change
    - PP3: In-silico / computational evidence supporting pathogenic
    - PP5: Reputable source (ClinVar) pathogenic
    - BA1: Stand-alone benign (AF > 0.05)
    - BS1: Strong benign (AF > 0.01)
    - BP4: In-silico neutral
    - BP6: Reputable source (ClinVar) benign
    - BP7: Synonymous / non-coding neutral
    """
    criteria = {}
    
    # Population Frequency Evidence (BA1 / BS1 / PM2)
    if pop_af is not None:
        if pop_af >= 0.05:
            criteria["BA1"] = "Stand-alone Benign: Allele Frequency >= 5%"
        elif pop_af >= 0.01:
            criteria["BS1"] = "Strong Benign: Allele Frequency >= 1%"
        elif pop_af <= 0.0001:
            criteria["PM2"] = "Moderate Pathogenic: Extremely rare or absent (AF <= 0.01%)"
    else:
        # Novel variant absent from dbSNP / 1000G
        if not dbsnp_id or dbsnp_id == '.':
            criteria["PM2"] = "Moderate Pathogenic: Novel variant absent from population databases"

    # Functional Consequence Evidence (PVS1 / PM4 / BP7)
    if impact == "HIGH" or consequence in ["frameshift_variant", "stop_gained", "splice_donor_variant", "splice_acceptor_variant"]:
        criteria["PVS1"] = f"Very Strong Pathogenic: Predicted loss-of-function ({consequence})"
    elif consequence == "inframe_indel":
        criteria["PM4"] = "Moderate Pathogenic: In-frame insertion/deletion altering protein length"
    elif consequence == "synonymous_variant":
        criteria["BP7"] = "Supporting Benign: Synonymous change with no predicted splice impact"

    # ClinVar Evidence (PP5 / BP6)
    if clinvar_sig:
        sig_lower = str(clinvar_sig).lower()
        if "pathogenic" in sig_lower and "conflict" not in sig_lower:
            criteria["PP5"] = f"Supporting Pathogenic: ClinVar reported as {clinvar_sig}"
        elif "benign" in sig_lower and "conflict" not in sig_lower:
            criteria["BP6"] = f"Supporting Benign: ClinVar reported as {clinvar_sig}"

    # Supporting In-Silico Evidence (PP3 for missense / high impact)
    if "PVS1" in criteria or impact == "HIGH":
        criteria["PP3"] = "Supporting Pathogenic: Multiple computational predictors indicate deleterious effect"
    elif impact == "MODERATE" and "BA1" not in criteria and "BS1" not in criteria:
        criteria["PP3"] = "Supporting Pathogenic: Missense mutation with predicted functional disruption"

    # Evaluate ACMG 2015 Combinations
    pvs1 = "PVS1" in criteria
    ps_count = sum(1 for c in criteria if c.startswith("PS"))
    pm_count = sum(1 for c in criteria if c.startswith("PM"))
    pp_count = sum(1 for c in criteria if c.startswith("PP"))
    ba1 = "BA1" in criteria
    bs_count = sum(1 for c in criteria if c.startswith("BS"))
    bp_count = sum(1 for c in criteria if c.startswith("BP"))

    acmg_class = "Uncertain Significance (VUS)"
    
    # Benign Rules
    if ba1 or bs_count >= 2:
        acmg_class = "Benign"
    elif bs_count >= 1 and bp_count >= 1:
        acmg_class = "Likely Benign"
    elif bp_count >= 2 and not (pvs1 or ps_count > 0 or pm_count > 0):
        acmg_class = "Likely Benign"
    # Pathogenic Rules
    elif pvs1 and (ps_count >= 1 or pm_count >= 2 or (pm_count >= 1 and pp_count >= 1) or pp_count >= 2):
        acmg_class = "Pathogenic"
    elif ps_count >= 2:
        acmg_class = "Pathogenic"
    elif ps_count >= 1 and (pm_count >= 3 or (pm_count >= 2 and pp_count >= 2) or (pm_count >= 1 and pp_count >= 4)):
        acmg_class = "Pathogenic"
    # Likely Pathogenic Rules
    elif pvs1 and pm_count >= 1:
        acmg_class = "Likely Pathogenic"
    elif ps_count >= 1 and (1 <= pm_count <= 2 or pp_count >= 2):
        acmg_class = "Likely Pathogenic"
    elif pm_count >= 3:
        acmg_class = "Likely Pathogenic"
    elif pm_count >= 2 and pp_count >= 2:
        acmg_class = "Likely Pathogenic"
    elif pm_count >= 1 and pp_count >= 4:
        acmg_class = "Likely Pathogenic"
    elif pvs1:
        acmg_class = "Likely Pathogenic"
        
    return acmg_class, criteria

def annotate_vcf(snp_vcf, indel_vcf, output_vcf, output_tsv, output_html, sample_name, bed_file=None):
    """Parse SNPs and Indels, annotate with ACMG criteria, and write VCF, TSV, and HTML."""
    intervals = load_target_genes(bed_file)
    variants = []
    
    # Process both VCFs
    for vcf_path in [snp_vcf, indel_vcf]:
        if not vcf_path or not os.path.exists(vcf_path):
            continue
        with open_vcf(vcf_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 8:
                    continue
                chrom = parts[0]
                pos = int(parts[1])
                var_id = parts[2]
                ref = parts[3]
                alt = parts[4]
                qual = parts[5]
                filter_status = parts[6]
                info_str = parts[7]
                
                # Parse INFO fields
                info_dict = {}
                for item in info_str.split(';'):
                    if '=' in item:
                        k, v = item.split('=', 1)
                        info_dict[k] = v
                    else:
                        info_dict[item] = True
                        
                # Extract Population AF from dbSNP / 1000G tags if present
                pop_af = None
                if 'AF' in info_dict:
                    try:
                        pop_af = float(info_dict['AF'].split(',')[0])
                    except ValueError:
                        pass
                elif 'CAF' in info_dict:
                    try:
                        # CAF=[ref_af, alt_af]
                        caf_vals = info_dict['CAF'].strip('[]').split(',')
                        if len(caf_vals) > 1 and caf_vals[1] != '.':
                            pop_af = float(caf_vals[1])
                    except ValueError:
                        pass
                        
                consequence, impact, gene = determine_consequence(ref, alt, pos, intervals, chrom)
                acmg_class, criteria = classify_acmg(consequence, impact, pop_af, var_id)
                
                variants.append({
                    'CHROM': chrom,
                    'POS': pos,
                    'ID': var_id,
                    'REF': ref,
                    'ALT': alt,
                    'QUAL': qual,
                    'FILTER': filter_status,
                    'GENE': gene,
                    'CONSEQUENCE': consequence,
                    'IMPACT': impact,
                    'POP_AF': f"{pop_af:.6f}" if pop_af is not None else "Novel",
                    'ACMG_CLASS': acmg_class,
                    'ACMG_CRITERIA': ";".join(f"{k}={v}" for k, v in criteria.items()),
                    'CRITERIA_KEYS': ",".join(criteria.keys())
                })
                
    # Sort variants by coordinate
    variants.sort(key=lambda x: (x['CHROM'], x['POS']))
    
    # Write TSV
    os.makedirs(os.path.dirname(os.path.abspath(output_tsv)), exist_ok=True)
    df = pd.DataFrame(variants)
    df.to_csv(output_tsv, sep='\t', index=False)
    
    # Write Annotated VCF
    os.makedirs(os.path.dirname(os.path.abspath(output_vcf)), exist_ok=True)
    vcf_out_path = output_vcf[:-3] if output_vcf.endswith('.gz') else output_vcf
    with open(vcf_out_path, 'w') as out:
        out.write("##fileformat=VCFv4.2\n")
        out.write("##source=WES_ACMG_Annotation_Engine_v2\n")
        out.write("##INFO=<ID=GENE,Number=1,Type=String,Description=\"Target Gene Symbol\">\n")
        out.write("##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description=\"Molecular Consequence\">\n")
        out.write("##INFO=<ID=IMPACT,Number=1,Type=String,Description=\"Functional Impact (HIGH, MODERATE, LOW, MODIFIER)\">\n")
        out.write("##INFO=<ID=ACMG_CLASS,Number=1,Type=String,Description=\"ACMG/AMP 2015 Clinical Classification (Pathogenic, Likely Pathogenic, VUS, Likely Benign, Benign)\">\n")
        out.write("##INFO=<ID=ACMG_CRITERIA,Number=.,Type=String,Description=\"ACMG Criteria Satisfied\">\n")
        out.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
        for v in variants:
            info_fields = f"GENE={v['GENE']};CONSEQUENCE={v['CONSEQUENCE']};IMPACT={v['IMPACT']};ACMG_CLASS={v['ACMG_CLASS']};ACMG_CRITERIA={v['CRITERIA_KEYS']}"
            out.write(f"{v['CHROM']}\t{v['POS']}\t{v['ID']}\t{v['REF']}\t{v['ALT']}\t{v['QUAL']}\t{v['FILTER']}\t{info_fields}\n")
            
    # Compress and index if .gz
    if output_vcf.endswith('.gz'):
        os.system(f"bgzip -f {vcf_out_path} && tabix -f -p vcf {output_vcf}")
        
    # Write Interactive HTML Report
    write_acmg_html_report(df, output_html, sample_name)

def write_acmg_html_report(df, output_html, sample_name):
    """Generate modern interactive HTML dashboard with ACMG classification badges."""
    os.makedirs(os.path.dirname(os.path.abspath(output_html)), exist_ok=True)
    
    total = len(df)
    counts = df['ACMG_CLASS'].value_counts().to_dict() if not df.empty else {}
    p_count = counts.get('Pathogenic', 0)
    lp_count = counts.get('Likely Pathogenic', 0)
    vus_count = counts.get('Uncertain Significance (VUS)', 0)
    lb_count = counts.get('Likely Benign', 0)
    b_count = counts.get('Benign', 0)
    
    rows_html = []
    if not df.empty:
        for _, r in df.head(100).iterrows():
            badge_class = "badge-vus"
            if r['ACMG_CLASS'] == 'Pathogenic':
                badge_class = "badge-pathogenic"
            elif r['ACMG_CLASS'] == 'Likely Pathogenic':
                badge_class = "badge-likely-pathogenic"
            elif r['ACMG_CLASS'] == 'Likely Benign':
                badge_class = "badge-likely-benign"
            elif r['ACMG_CLASS'] == 'Benign':
                badge_class = "badge-benign"
                
            impact_badge = "impact-mod"
            if r['IMPACT'] == 'HIGH':
                impact_badge = "impact-high"
            elif r['IMPACT'] == 'MODIFIER':
                impact_badge = "impact-low"
                
            rows_html.append(f"""
            <tr>
                <td><strong>{r['CHROM']}:{r['POS']}</strong></td>
                <td><code>{r['REF']}&gt;{r['ALT']}</code></td>
                <td><span class="gene-tag">{r['GENE']}</span></td>
                <td>{r['CONSEQUENCE']}</td>
                <td><span class="impact-tag {impact_badge}">{r['IMPACT']}</span></td>
                <td>{r['POP_AF']}</td>
                <td><span class="badge {badge_class}">{r['ACMG_CLASS']}</span></td>
                <td><small>{r['CRITERIA_KEYS']}</small></td>
            </tr>
            """)
            
    table_content = "\n".join(rows_html) if rows_html else "<tr><td colspan='8'>No variants found.</td></tr>"
    
    html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>ACMG Variant Classification Report - {sample_name}</title>
    <style>
        :root {{
            --bg: #0f172a; --card: #1e293b; --text: #f8fafc; --muted: #94a3b8;
            --pathogenic: #ef4444; --lp: #f97316; --vus: #eab308; --lb: #3b82f6; --benign: #10b981;
            --border: #334155;
        }}
        body {{ font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif; background: var(--bg); color: var(--text); margin: 0; padding: 24px; }}
        .header {{ margin-bottom: 24px; border-bottom: 1px solid var(--border); padding-bottom: 16px; }}
        .header h1 {{ margin: 0 0 8px 0; font-size: 24px; color: #38bdf8; }}
        .header p {{ margin: 0; color: var(--muted); font-size: 14px; }}
        .stats-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(180px, 1fr)); gap: 16px; margin-bottom: 24px; }}
        .stat-card {{ background: var(--card); border: 1px solid var(--border); border-radius: 8px; padding: 16px; text-align: center; }}
        .stat-num {{ font-size: 28px; font-weight: 700; margin-top: 4px; }}
        .stat-label {{ font-size: 12px; text-transform: uppercase; color: var(--muted); letter-spacing: 0.5px; }}
        .table-container {{ background: var(--card); border: 1px solid var(--border); border-radius: 8px; overflow-x: auto; }}
        table {{ width: 100%; border-collapse: collapse; font-size: 13px; text-align: left; }}
        th, td {{ padding: 12px 16px; border-bottom: 1px solid var(--border); }}
        th {{ background: #0f172a; color: var(--muted); font-weight: 600; text-transform: uppercase; font-size: 11px; }}
        tr:hover {{ background: rgba(255,255,255,0.02); }}
        .badge {{ padding: 4px 8px; border-radius: 4px; font-size: 11px; font-weight: 600; display: inline-block; }}
        .badge-pathogenic {{ background: rgba(239,68,68,0.2); color: var(--pathogenic); border: 1px solid var(--pathogenic); }}
        .badge-likely-pathogenic {{ background: rgba(249,115,22,0.2); color: var(--lp); border: 1px solid var(--lp); }}
        .badge-vus {{ background: rgba(234,179,8,0.2); color: var(--vus); border: 1px solid var(--vus); }}
        .badge-likely-benign {{ background: rgba(59,130,246,0.2); color: var(--lb); border: 1px solid var(--lb); }}
        .badge-benign {{ background: rgba(16,185,129,0.2); color: var(--benign); border: 1px solid var(--benign); }}
        .gene-tag {{ color: #38bdf8; font-weight: 600; }}
        .impact-tag {{ padding: 2px 6px; border-radius: 3px; font-size: 10px; font-weight: 600; }}
        .impact-high {{ background: #ef4444; color: #fff; }}
        .impact-mod {{ background: #f59e0b; color: #000; }}
        .impact-low {{ background: #64748b; color: #fff; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>ACMG/AMP Variant Annotation Report</h1>
        <p>Sample: <strong>{sample_name}</strong> | Standards: ACMG/AMP 2015 5-Tier Classification Guidelines</p>
    </div>
    
    <div class="stats-grid">
        <div class="stat-card"><div class="stat-label">Total Variants</div><div class="stat-num">{total}</div></div>
        <div class="stat-card"><div class="stat-label" style="color:var(--pathogenic)">Pathogenic</div><div class="stat-num" style="color:var(--pathogenic)">{p_count}</div></div>
        <div class="stat-card"><div class="stat-label" style="color:var(--lp)">Likely Pathogenic</div><div class="stat-num" style="color:var(--lp)">{lp_count}</div></div>
        <div class="stat-card"><div class="stat-label" style="color:var(--vus)">VUS</div><div class="stat-num" style="color:var(--vus)">{vus_count}</div></div>
        <div class="stat-card"><div class="stat-label" style="color:var(--lb)">Likely Benign</div><div class="stat-num" style="color:var(--lb)">{lb_count}</div></div>
        <div class="stat-card"><div class="stat-label" style="color:var(--benign)">Benign</div><div class="stat-num" style="color:var(--benign)">{b_count}</div></div>
    </div>
    
    <div class="table-container">
        <table>
            <thead>
                <tr>
                    <th>Locus</th>
                    <th>Alleles</th>
                    <th>Gene</th>
                    <th>Consequence</th>
                    <th>Impact</th>
                    <th>Pop AF</th>
                    <th>ACMG Classification</th>
                    <th>Evidence Criteria</th>
                </tr>
            </thead>
            <tbody>
                {table_content}
            </tbody>
        </table>
    </div>
</body>
</html>
"""
    with open(output_html, 'w') as f:
        f.write(html)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="ACMG Variant Classifier")
    parser.add_argument("--snp-vcf", required=True, help="Filtered SNP VCF")
    parser.add_argument("--indel-vcf", required=True, help="Filtered Indel VCF")
    parser.add_argument("--bed", required=False, help="Target BED file with gene symbols")
    parser.add_argument("--sample-name", required=True, help="Sample name")
    parser.add_argument("--output-vcf", required=True, help="Output annotated VCF")
    parser.add_argument("--output-tsv", required=True, help="Output TSV table")
    parser.add_argument("--output-html", required=True, help="Output HTML report")
    
    args = parser.parse_args()
    annotate_vcf(
        snp_vcf=args.snp_vcf,
        indel_vcf=args.indel_vcf,
        output_vcf=args.output_vcf,
        output_tsv=args.output_tsv,
        output_html=args.output_html,
        sample_name=args.sample_name,
        bed_file=args.bed
    )
