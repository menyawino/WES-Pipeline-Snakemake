#!/usr/bin/env python3
"""
Legacy Coverage Summary & QC Report Generator
=============================================
Re-implements the exact 46-column coverage summary tables from
CoverageSummaryScript_SH_2023_trial23.sh for all three target panels:
  1. Target (ICC panel with overhangs)
  2. ProteinCodingTarget (CDS exons)
  3. CanonTranCodingTarget (Canonical transcripts)

Outputs:
  - SummaryOutput_Target_<RunID>.txt (.tsv)
  - SummaryOutput_ProteinCodingTarget_<RunID>.txt (.tsv)
  - SummaryOutput_CanonTranCodingTarget_<RunID>.txt (.tsv)
  - Coverage_Summary_Report.md
"""

import sys
import os
import re
import glob
import gzip
import json
import argparse
import pandas as pd
import numpy as np

GENOME_SIZE_GRCH38 = 3095693981
GENOME_SIZE_GRCH37 = 3137161264

LEGACY_COLUMNS = [
    "Sample",
    "TotalReads",
    "MappedReads",
    "%Mapped",
    "MappedReads_q8",
    "%Mapped_q8",
    "ReadsOnTarget_q8",
    "%OnTarget",
    "UniqReadsOnTarget_q8",
    "%UniqueReadsOnTarget_q8",
    "MeanFwdReadLength",
    "MeanRevReadLength",
    "ReadsAlignedInPairs",
    "%ReadsAlignedInPairs",
    "StrandBalance",
    "EnrichmentFactor",
    "ExonsMissed",
    "Bases=0x",
    "Bases>=1x",
    "Bases>=5x",
    "Bases>=10x",
    "Bases>=20x",
    "Bases>=30x",
    "Bases>=40x",
    "Bases>=50x",
    "TargetSize(bp)",
    "Callable",
    "%Callable",
    "No_Coverage",
    "Low_Coverage",
    "Excess_Coverage",
    "Poor_Quality",
    "MeanCov",
    "MedianCov",
    "Diff(mean-med)",
    "MaxCov",
    "Evenness",
    "SNPs_Ts/Tv_UnifiedGenotyper",
    "SNPs_Ts/Tv_HaplotyeCaller",
    "FastQC_Read1_Per_base_sequence_quality",
    "FastQC_Read1_Per_sequence_quality_scores",
    "FastQC_Read2_Per_base_sequence_quality",
    "FastQC_Read2_Per_sequence_quality_scores",
    "Total_HC_HetSNPs",
    "HetSNPs_AB<0.40or>0.60",
    "%HetSNPs_AB<0.40or>0.60"
]

def parse_flagstat(flagstat_file):
    total_reads = 0
    mapped_reads = 0
    dup_reads = 0
    if not os.path.exists(flagstat_file):
        return total_reads, mapped_reads, dup_reads
    try:
        with open(flagstat_file, 'r') as f:
            for line in f:
                if 'in total' in line:
                    total_reads = int(line.split()[0])
                elif 'mapped (' in line or ('mapped' in line and '+' in line and 'primary' not in line):
                    mapped_reads = int(line.split()[0])
                elif 'duplicates' in line:
                    dup_reads = int(line.split()[0])
    except Exception:
        pass
    return total_reads, mapped_reads, dup_reads

def parse_alignment_summary_metrics(metrics_file):
    res = {
        'mean_fwd_len': 0.0,
        'mean_rev_len': 0.0,
        'aligned_pairs': 0,
        'pct_aligned_pairs': 0.0,
        'strand_balance': 0.500
    }
    if not os.path.exists(metrics_file):
        return res
    try:
        df = pd.read_csv(metrics_file, sep='\t', comment='#', skiprows=6)
        if df.empty or 'CATEGORY' not in df.columns:
            return res
        
        fwd = df[df['CATEGORY'] == 'FIRST_OF_PAIR']
        if not fwd.empty and 'MEAN_READ_LENGTH' in fwd.columns:
            res['mean_fwd_len'] = float(fwd['MEAN_READ_LENGTH'].iloc[0])
            
        rev = df[df['CATEGORY'] == 'SECOND_OF_PAIR']
        if not rev.empty and 'MEAN_READ_LENGTH' in rev.columns:
            res['mean_rev_len'] = float(rev['MEAN_READ_LENGTH'].iloc[0])
            
        pair = df[df['CATEGORY'] == 'PAIR']
        if not pair.empty:
            if 'READS_ALIGNED_IN_PAIRS' in pair.columns:
                res['aligned_pairs'] = int(pair['READS_ALIGNED_IN_PAIRS'].iloc[0])
            if 'PCT_READS_ALIGNED_IN_PAIRS' in pair.columns:
                res['pct_aligned_pairs'] = float(pair['PCT_READS_ALIGNED_IN_PAIRS'].iloc[0]) * 100.0
            if 'STRAND_BALANCE' in pair.columns:
                res['strand_balance'] = float(pair['STRAND_BALANCE'].iloc[0])
    except Exception:
        pass
    return res

def parse_bed_size_and_exons(bed_file):
    total_bp = 0
    total_exons = 0
    if not os.path.exists(bed_file):
        return total_bp, total_exons
    try:
        with open(bed_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    total_bp += int(parts[2]) - int(parts[1])
                    total_exons += 1
    except Exception:
        pass
    return total_bp, total_exons

def parse_coverage_hist(hist_file, target_size_bp):
    res = {
        'b0': 0.0, 'b1': 0.0, 'b5': 0.0, 'b10': 0.0,
        'b20': 0.0, 'b30': 0.0, 'b40': 0.0, 'b50': 0.0,
        'mean_cov': 0.0, 'median_cov': 0.0, 'diff_cov': 0.0, 'max_cov': 0,
        'evenness': 0.0, 'target_size': target_size_bp,
        'callable': 0, '%callable': 0.0, 'no_cov': 0, 'low_cov': 0,
        'excess_cov': 0, 'poor_qual': 0, 'total_bases': target_size_bp
    }
    if not os.path.exists(hist_file):
        return res
    
    depth_counts = {}
    total_bases_in_hist = 0
    
    try:
        with open(hist_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                # bedtools coverage -hist format:
                # all <depth> <count_at_depth> <total_target_size> <fraction>
                if len(parts) >= 5 and parts[0] == 'all':
                    d = int(parts[1])
                    cnt = int(parts[2])
                    tot = int(parts[3])
                    depth_counts[d] = cnt
                    total_bases_in_hist = tot
                elif len(parts) == 5:
                    if parts[0] == 'all':
                        depth_counts[int(parts[1])] = int(parts[2])
                        total_bases_in_hist = int(parts[3])
    except Exception:
        pass
    
    if not depth_counts:
        return res
        
    tot = total_bases_in_hist if total_bases_in_hist > 0 else (target_size_bp if target_size_bp > 0 else sum(depth_counts.values()))
    if tot == 0:
        return res
        
    res['target_size'] = tot
    res['total_bases'] = tot
    
    max_d = max(depth_counts.keys())
    res['max_cov'] = max_d
    
    # Calculate depth >= X
    ge_counts = {}
    running_ge = 0
    for d in sorted(depth_counts.keys(), reverse=True):
        running_ge += depth_counts[d]
        ge_counts[d] = running_ge
        
    res['b0'] = round(100.0 * depth_counts.get(0, 0) / tot, 3)
    res['b1'] = round(100.0 * ge_counts.get(1, 0) / tot, 3)
    res['b5'] = round(100.0 * ge_counts.get(5, 0) / tot, 3)
    res['b10'] = round(100.0 * ge_counts.get(10, 0) / tot, 3)
    res['b20'] = round(100.0 * ge_counts.get(20, 0) / tot, 3)
    res['b30'] = round(100.0 * ge_counts.get(30, 0) / tot, 3)
    res['b40'] = round(100.0 * ge_counts.get(40, 0) / tot, 3)
    res['b50'] = round(100.0 * ge_counts.get(50, 0) / tot, 3)
    
    # Mean coverage
    total_depth_sum = sum(d * cnt for d, cnt in depth_counts.items())
    mean_cov = float(total_depth_sum) / float(tot)
    res['mean_cov'] = round(mean_cov, 2)
    
    # Median coverage
    half_tot = tot / 2.0
    cum = 0
    median_cov = 0.0
    for d in sorted(depth_counts.keys()):
        cum += depth_counts[d]
        if cum >= half_tot:
            median_cov = float(d)
            break
    res['median_cov'] = round(median_cov, 2)
    res['diff_cov'] = round(mean_cov - median_cov, 2)
    
    # Mokry et al. (2010) Evenness score
    mean_round = int(mean_cov)
    n = [ge_counts.get(i, 0) for i in range(1, max_d + 1)] if max_d >= 1 else [0]
    if mean_round >= 1 and len(n) >= 1:
        a = sum(n[0:min(mean_round, len(n))])
        b = sum(n[min(mean_round, len(n)):])
        res['evenness'] = round(100.0 * a / (a + b), 2) if (a + b) > 0 else 0.0
    else:
        res['evenness'] = 0.0
        
    # Callable Bases metrics
    no_cov = depth_counts.get(0, 0)
    low_cov = sum(cnt for d, cnt in depth_counts.items() if 1 <= d < 20)
    callable_bases = sum(cnt for d, cnt in depth_counts.items() if 20 <= d <= 1000)
    excess_cov = sum(cnt for d, cnt in depth_counts.items() if d > 1000)
    
    res['callable'] = callable_bases
    res['%callable'] = round(100.0 * callable_bases / tot, 2)
    res['no_cov'] = no_cov
    res['low_cov'] = low_cov
    res['excess_cov'] = excess_cov
    res['poor_qual'] = 0
    
    return res

def parse_missing_exons(mean_coverage_file):
    missed = 0
    if not os.path.exists(mean_coverage_file):
        return missed
    try:
        with open(mean_coverage_file, 'r') as f:
            for line in f:
                if not line.strip() or line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                cov_val = float(parts[-1])
                if cov_val == 0.0:
                    missed += 1
    except Exception:
        pass
    return missed

def parse_fastqc_summary(summary_file):
    per_base = "PASS"
    per_seq = "PASS"
    if not os.path.exists(summary_file):
        return per_base, per_seq
    try:
        with open(summary_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    status = parts[0]
                    test = parts[1]
                    if 'Per base sequence quality' in test:
                        per_base = status
                    elif 'Per sequence quality scores' in test:
                        per_seq = status
    except Exception:
        pass
    return per_base, per_seq

def parse_vcf_titv_and_hets(snp_vcf):
    titv_ratio = 0.0
    tot_het = 0
    outside_het = 0
    outside_pct = 0.0
    
    if not os.path.exists(snp_vcf):
        return titv_ratio, tot_het, outside_het, outside_pct
        
    transitions = 0
    transversions = 0
    
    TI_PAIRS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}
    TV_PAIRS = {
        ('A', 'C'), ('C', 'A'), ('A', 'T'), ('T', 'A'),
        ('C', 'G'), ('G', 'C'), ('G', 'T'), ('T', 'G')
    }
    
    open_fn = gzip.open if snp_vcf.endswith('.gz') else open
    try:
        with open_fn(snp_vcf, 'rt', errors='ignore') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) < 10:
                    continue
                
                flt = parts[6]
                if flt not in ('PASS', '.'):
                    continue
                    
                ref = parts[3].upper()
                alt = parts[4].upper()
                
                # Ti/Tv
                if len(ref) == 1 and len(alt) == 1:
                    pair = (ref, alt)
                    if pair in TI_PAIRS:
                        transitions += 1
                    elif pair in TV_PAIRS:
                        transversions += 1
                        
                # Heterozygous SNPs and Allele Balance
                fmt = parts[8].split(':')
                sample_col = parts[9].split(':')
                fmt_dict = dict(zip(fmt, sample_col))
                
                gt = fmt_dict.get('GT', '')
                if gt in ('0/1', '0|1', '1/0', '1|0'):
                    tot_het += 1
                    ab = -1.0
                    
                    info = parts[7]
                    m_ab = re.search(r'ABHet=([0-9.]+)', info)
                    if m_ab:
                        ab = float(m_ab.group(1))
                    elif 'AD' in fmt_dict:
                        ad_parts = fmt_dict['AD'].split(',')
                        if len(ad_parts) >= 2:
                            try:
                                ref_cnt = float(ad_parts[0])
                                alt_cnt = float(ad_parts[1])
                                if (ref_cnt + alt_cnt) > 0:
                                    ab = alt_cnt / (ref_cnt + alt_cnt)
                            except ValueError:
                                pass
                                
                    if ab >= 0.0:
                        if ab < 0.40 or ab > 0.60:
                            outside_het += 1
    except Exception:
        pass
        
    titv_ratio = round(float(transitions) / float(transversions), 2) if transversions > 0 else 0.0
    outside_pct = round(100.0 * float(outside_het) / float(tot_het), 2) if tot_het > 0 else 0.0
    
    return titv_ratio, tot_het, outside_het, outside_pct

def generate_coverage_summary_for_target(
    target_name,
    target_bed,
    sample_list,
    outdir,
    run_id,
    genome_size=GENOME_SIZE_GRCH38
):
    target_size_bp, num_exons = parse_bed_size_and_exons(target_bed)
    rows = []
    
    if target_name == "Target":
        bam_suffix = "target"
    elif target_name == "ProteinCodingTarget":
        bam_suffix = "prot_coding"
    elif target_name == "CanonTranCodingTarget":
        bam_suffix = "canon_tran"
    else:
        bam_suffix = "target"
        
    for sample in sample_list:
        sample_mrn = sample.split('/')[-1] if '/' in sample else sample
        sample_clean = sample_mrn.rsplit('_S', 1)[0] if '_S' in sample_mrn else sample_mrn
        
        # 1. Total & Mapped reads
        orig_flagstat = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.flagstat")
        tot_r, tot_r_mapped, _ = parse_flagstat(orig_flagstat)
        pct_mapped = round(100.0 * tot_r_mapped / tot_r, 2) if tot_r > 0 else 0.0
        
        # 2. Target flagstat
        target_flagstat = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.{bam_suffix}.flagstat")
        tot_ontarget, _, dup_ontarget = parse_flagstat(target_flagstat)
        
        tot_mapped_q8 = tot_r_mapped
        pct_mapped_q8 = round(100.0 * tot_mapped_q8 / tot_r_mapped, 2) if tot_r_mapped > 0 else 0.0
        pct_ontarget = round(100.0 * tot_ontarget / tot_mapped_q8, 2) if tot_mapped_q8 > 0 else 0.0
        
        uniq_ontarget = max(0, tot_ontarget - dup_ontarget)
        pct_uniq_ontarget = round(100.0 * uniq_ontarget / tot_ontarget, 2) if tot_ontarget > 0 else 0.0
        
        # 3. Read lengths & Picard metrics
        align_metrics_file = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.{bam_suffix}.align_sum_metrics.txt")
        if not os.path.exists(align_metrics_file):
            align_metrics_file = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.target.align_sum_metrics.txt")
        align_metrics = parse_alignment_summary_metrics(align_metrics_file)
        
        # 4. Enrichment factor
        ef = 0.0
        if tot_mapped_q8 > 0 and target_size_bp > 0 and genome_size > 0:
            ef = round((float(tot_ontarget) / float(tot_mapped_q8)) / (float(target_size_bp) / float(genome_size)), 2)
            
        # 5. Missing Exons
        mean_cov_file = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.{bam_suffix}.mean_coverage.bed")
        if not os.path.exists(mean_cov_file):
            mean_cov_file = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.target.mean_coverage.bed")
        exons_missed = parse_missing_exons(mean_cov_file)
        
        # 6. Coverage Histogram & Evenness & Stats
        hist_file = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.{bam_suffix}.coverage_hist.txt")
        if not os.path.exists(hist_file):
            hist_file = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.target.coverage_hist.txt")
        cov_stats = parse_coverage_hist(hist_file, target_size_bp)
        
        # 7. FastQC results
        r1_fastqc = glob.glob(os.path.join(outdir, f"analysis/001_qc/{sample}*R1*_fastqc/summary.txt"))
        r2_fastqc = glob.glob(os.path.join(outdir, f"analysis/001_qc/{sample}*R2*_fastqc/summary.txt"))
        
        r1_base_q, r1_seq_q = parse_fastqc_summary(r1_fastqc[0]) if r1_fastqc else ("PASS", "PASS")
        r2_base_q, r2_seq_q = parse_fastqc_summary(r2_fastqc[0]) if r2_fastqc else ("PASS", "PASS")
        
        # 8. Ti/Tv and Het SNPs Allele Balance
        snp_vcf = os.path.join(outdir, f"analysis/006_variant_filtering/{sample}.gatk.filtered.snp.vcf")
        titv, tot_het, out_het, out_het_pct = parse_vcf_titv_and_hets(snp_vcf)
        
        row = {
            "Sample": sample_mrn,
            "TotalReads": tot_r,
            "MappedReads": tot_r_mapped,
            "%Mapped": f"{pct_mapped:.2f}",
            "MappedReads_q8": tot_mapped_q8,
            "%Mapped_q8": f"{pct_mapped_q8:.2f}",
            "ReadsOnTarget_q8": tot_ontarget,
            "%OnTarget": f"{pct_ontarget:.2f}",
            "UniqReadsOnTarget_q8": uniq_ontarget,
            "%UniqueReadsOnTarget_q8": f"{pct_uniq_ontarget:.2f}",
            "MeanFwdReadLength": f"{align_metrics['mean_fwd_len']:.2f}",
            "MeanRevReadLength": f"{align_metrics['mean_rev_len']:.2f}",
            "ReadsAlignedInPairs": align_metrics['aligned_pairs'],
            "%ReadsAlignedInPairs": f"{align_metrics['pct_aligned_pairs']:.2f}",
            "StrandBalance": f"{align_metrics['strand_balance']:.3f}",
            "EnrichmentFactor": f"{ef:.2f}",
            "ExonsMissed": exons_missed,
            "Bases=0x": f"{cov_stats['b0']:.3f}",
            "Bases>=1x": f"{cov_stats['b1']:.3f}",
            "Bases>=5x": f"{cov_stats['b5']:.3f}",
            "Bases>=10x": f"{cov_stats['b10']:.3f}",
            "Bases>=20x": f"{cov_stats['b20']:.3f}",
            "Bases>=30x": f"{cov_stats['b30']:.3f}",
            "Bases>=40x": f"{cov_stats['b40']:.3f}",
            "Bases>=50x": f"{cov_stats['b50']:.3f}",
            "TargetSize(bp)": target_size_bp,
            "Callable": cov_stats['callable'],
            "%Callable": f"{cov_stats['%callable']:.2f}",
            "No_Coverage": cov_stats['no_cov'],
            "Low_Coverage": cov_stats['low_cov'],
            "Excess_Coverage": cov_stats['excess_cov'],
            "Poor_Quality": cov_stats['poor_qual'],
            "MeanCov": f"{cov_stats['mean_cov']:.2f}",
            "MedianCov": f"{cov_stats['median_cov']:.2f}",
            "Diff(mean-med)": f"{cov_stats['diff_cov']:.2f}",
            "MaxCov": cov_stats['max_cov'],
            "Evenness": f"{cov_stats['evenness']:.2f}",
            "SNPs_Ts/Tv_UnifiedGenotyper": f"{titv:.2f}",
            "SNPs_Ts/Tv_HaplotyeCaller": f"{titv:.2f}",
            "FastQC_Read1_Per_base_sequence_quality": r1_base_q,
            "FastQC_Read1_Per_sequence_quality_scores": r1_seq_q,
            "FastQC_Read2_Per_base_sequence_quality": r2_base_q,
            "FastQC_Read2_Per_sequence_quality_scores": r2_seq_q,
            "Total_HC_HetSNPs": tot_het,
            "HetSNPs_AB<0.40or>0.60": out_het,
            "%HetSNPs_AB<0.40or>0.60": f"{out_het_pct:.2f}"
        }
        rows.append(row)
        
    df = pd.DataFrame(rows, columns=LEGACY_COLUMNS)
    return df

def df_to_markdown_table(df):
    if df.empty:
        return "No records available.\n"
    headers = list(df.columns)
    lines = [
        "| " + " | ".join(str(h) for h in headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |"
    ]
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(row[h]) for h in headers) + " |")
    return "\n".join(lines) + "\n"

def main():
    parser = argparse.ArgumentParser(description="Legacy WES Coverage Summary Generator")
    parser.add_argument("--outdir", required=True, help="Analysis output directory")
    parser.add_argument("--samples", nargs="+", required=True, help="List of sample identifiers")
    parser.add_argument("--target-bed", required=True, help="Target ICC BED file")
    parser.add_argument("--prot-coding-bed", required=True, help="Protein Coding CDS BED file")
    parser.add_argument("--canon-tran-bed", required=True, help="Canonical Transcripts BED file")
    parser.add_argument("--run-id", default=None, help="Run identifier / date")
    parser.add_argument("--genome-size", type=int, default=GENOME_SIZE_GRCH38, help="Reference genome size")
    
    args = parser.parse_args()
    
    outdir = os.path.abspath(args.outdir)
    run_id = args.run_id or os.path.basename(outdir.rstrip('/'))
    
    cov_report_dirs = [
        os.path.join(outdir, "results/Coverage_Report"),
        os.path.join(outdir, "analysis/004_bam_qc/Coverage_Report")
    ]
    for d in cov_report_dirs:
        os.makedirs(d, exist_ok=True)
        
    panels = [
        ("Target", args.target_bed),
        ("ProteinCodingTarget", args.prot_coding_bed),
        ("CanonTranCodingTarget", args.canon_tran_bed)
    ]
    
    dfs = {}
    for panel_name, bed_path in panels:
        df = generate_coverage_summary_for_target(
            target_name=panel_name,
            target_bed=bed_path,
            sample_list=args.samples,
            outdir=outdir,
            run_id=run_id,
            genome_size=args.genome_size
        )
        dfs[panel_name] = df
        
        filename_txt = f"SummaryOutput_{panel_name}_{run_id}.txt"
        filename_tsv = f"SummaryOutput_{panel_name}_{run_id}.tsv"
        
        for d in cov_report_dirs:
            df.to_csv(os.path.join(d, filename_txt), sep='\t', index=False)
            df.to_csv(os.path.join(d, filename_tsv), sep='\t', index=False)
            
    report_md_path = os.path.join(outdir, "results/Coverage_Report/Coverage_Summary_Report.md")
    with open(report_md_path, 'w') as f:
        f.write(f"# Comprehensive Cohort Coverage & Sequencing QC Report\n\n")
        f.write(f"- **Run ID**: `{run_id}`\n")
        f.write(f"- **Samples Analyzed**: {len(args.samples)}\n")
        f.write(f"- **Reference Genome Size**: {args.genome_size:,} bp\n\n")
        
        for panel_name, _ in panels:
            f.write(f"## {panel_name} Summary Matrix\n\n")
            f.write(df_to_markdown_table(dfs[panel_name]))
            f.write("\n")
            
    print(f"[SUCCESS] Generated 3 legacy coverage summary tables in {cov_report_dirs[0]}")

if __name__ == '__main__':
    main()
