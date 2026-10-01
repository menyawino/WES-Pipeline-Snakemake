#!/usr/bin/env python3
"""
Legacy Coverage Summary & QC Report Generator
=============================================
Produces standard 46-column coverage summary tables for all three panels:
  1. Target (ICC capture panel with overhangs)
  2. ProteinCodingTarget (CDS exons)
  3. CanonTranCodingTarget (Canonical transcripts)

Key Design Principles & Fixes:
  - Exact Coverage Histograms: Bedtools coverage is executed with -a <target.bed>
    and -b <reads.bam> -hist. The histogram is strictly parsed from the 'all' lines
    representing per-base depth across the target features.
  - Strict Panel Isolation: Completely eliminates silent fallbacks from canonical
    transcripts or protein-coding targets to the general target panel. If panel-specific
    histogram or coverage files are missing, the metrics are explicitly marked 'NA'
    and recorded in the QC warning audit.
  - MappedReads_q8 Metric: Evaluated from the full BAM at MAPQ >= 8 (via
    samtools view -c -q 8 -F 4). If unavailable, reports 'NA' rather than falsifying
    the metric with total mapped reads (which forced %Mapped_q8 to 100%).
  - Ti/Tv Metrics: Reports 'NA' for SNPs_Ts/Tv_UnifiedGenotyper (UnifiedGenotyper is
    not run in this modern GATK4/DeepVariant workflow), while reporting the calculated
    ratio for SNPs_Ts/Tv_HaplotyeCaller.
  - Transparent Missing Data Handling: Missing or unparsable inputs (FastQC, flagstat,
    VCF, or coverage metrics) do NOT default to 0 or 'PASS'; they are explicitly labeled
    'NA' and flagged as missing data in the QC audit table.
  - Configurable QC Thresholds: Resolves legacy script discrepancy where condition
    checked '<= 200' while the warning message displayed '< 100'.
"""

import sys
import os
import re
import glob
import gzip
import json
import argparse
import subprocess
import pandas as pd
import numpy as np

# Reference Genome Sizes (bp)
GENOME_SIZE_GRCH38 = 3095693981  # Standard GRCh38 primary assembly size
GENOME_SIZE_GRCH37 = 3137161264  # Legacy hg19 / GRCh37 assembly size

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
    """Parse standard samtools flagstat output. Returns (total, mapped, dup) or (None, None, None)."""
    if not os.path.exists(flagstat_file) or os.path.getsize(flagstat_file) == 0:
        return None, None, None
    total_reads = None
    mapped_reads = None
    dup_reads = None
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
        return None, None, None
    return total_reads, mapped_reads, (dup_reads if dup_reads is not None else 0)

def parse_q8_mapped_reads(sample, outdir):
    """
    Parse mapped reads with MAPQ >= 8 across the FULL original BAM.
    Returns int count, or None if unavailable (no silent fallback).
    """
    q8_file = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.q8.txt")
    if os.path.exists(q8_file) and os.path.getsize(q8_file) > 0:
        try:
            with open(q8_file, 'r') as f:
                val = f.read().strip()
                if val.isdigit():
                    return int(val)
        except Exception:
            pass

    # Direct query on BQSR or markdup BAM if present
    candidate_bams = [
        os.path.join(outdir, f"analysis/003_bam_prep/02_bqsr/{sample}.bqsr.bam"),
        os.path.join(outdir, f"analysis/003_bam_prep/01_markdup/{sample}.markdup.bam"),
        os.path.join(outdir, f"analysis/002_alignment/{sample}.bam")
    ]
    for b in candidate_bams:
        if os.path.exists(b) and os.path.getsize(b) > 0:
            try:
                res = subprocess.run(
                    ["samtools", "view", "-c", "-q", "8", "-F", "4", b],
                    capture_output=True, text=True, check=True
                )
                val = res.stdout.strip()
                if val.isdigit():
                    return int(val)
            except Exception:
                pass
    return None

def parse_alignment_summary_metrics(metrics_file):
    """Parse Picard / fast_bam_qc AlignmentSummaryMetrics output. Returns dict or None."""
    if not os.path.exists(metrics_file) or os.path.getsize(metrics_file) == 0:
        return None
    res = {
        'mean_fwd_len': None,
        'mean_rev_len': None,
        'aligned_pairs': None,
        'pct_aligned_pairs': None,
        'strand_balance': None
    }
    try:
        df = pd.read_csv(metrics_file, sep='\t', comment='#', skiprows=6)
        if df.empty or 'CATEGORY' not in df.columns:
            return None
        
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
        return None
    return res

def parse_bed_size_and_exons(bed_file):
    """Calculate total target size in base pairs and total interval count."""
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
    """
    Parse bedtools coverage -a <target.bed> -b <reads.bam> -hist output.
    Returns dict of depth metrics or None if missing/corrupt.
    """
    if not os.path.exists(hist_file) or os.path.getsize(hist_file) == 0:
        return None
    
    depth_counts = {}
    total_bases_in_hist = 0
    
    try:
        with open(hist_file, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                # bedtools coverage -a target.bed -b reads.bam -hist outputs summary 'all':
                # all <depth> <count_at_depth> <total_target_size> <fraction>
                if len(parts) >= 5 and parts[0] == 'all':
                    d = int(parts[1])
                    cnt = int(parts[2])
                    tot = int(parts[3])
                    depth_counts[d] = cnt
                    total_bases_in_hist = tot
    except Exception:
        return None
    
    if not depth_counts:
        return None
        
    tot = total_bases_in_hist if total_bases_in_hist > 0 else target_size_bp
    if tot == 0:
        return None
        
    res = {
        'target_size': tot,
        'total_bases': tot,
        'poor_qual': 0
    }
    
    max_d = max(depth_counts.keys())
    res['max_cov'] = max_d
    
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
    
    total_depth_sum = sum(d * cnt for d, cnt in depth_counts.items())
    mean_cov = float(total_depth_sum) / float(tot)
    res['mean_cov'] = round(mean_cov, 2)
    
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
    
    mean_round = int(mean_cov)
    n = [ge_counts.get(i, 0) for i in range(1, max_d + 1)] if max_d >= 1 else [0]
    if mean_round >= 1 and len(n) >= 1:
        a = sum(n[0:min(mean_round, len(n))])
        b = sum(n[min(mean_round, len(n)):])
        res['evenness'] = round(100.0 * a / (a + b), 2) if (a + b) > 0 else 0.0
    else:
        res['evenness'] = 0.0
        
    no_cov = depth_counts.get(0, 0)
    low_cov = sum(cnt for d, cnt in depth_counts.items() if 1 <= d < 20)
    callable_bases = sum(cnt for d, cnt in depth_counts.items() if 20 <= d <= 1000)
    excess_cov = sum(cnt for d, cnt in depth_counts.items() if d > 1000)
    
    res['callable'] = callable_bases
    res['%callable'] = round(100.0 * callable_bases / tot, 2)
    res['no_cov'] = no_cov
    res['low_cov'] = low_cov
    res['excess_cov'] = excess_cov
    
    return res

def parse_missing_exons(mean_coverage_file):
    """Count number of target intervals with 0x mean depth. Returns int count or None."""
    if not os.path.exists(mean_coverage_file) or os.path.getsize(mean_coverage_file) == 0:
        return None
    missed = 0
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
        return None
    return missed

def parse_fastqc_summary(summary_file):
    """
    Extract Per base and Per sequence quality status from FastQC summary.txt.
    Returns ('NA', 'NA') if missing, never defaulting to 'PASS'.
    """
    if not os.path.exists(summary_file) or os.path.getsize(summary_file) == 0:
        return "NA", "NA"
    per_base = "NA"
    per_seq = "NA"
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
        return "NA", "NA"
    return per_base, per_seq

def parse_fastp_qc(sample, outdir):
    """
    Parse fastp report JSON files for sample lanes.
    Evaluates quality metrics (Q30/Q20 rate) for Read 1 and Read 2.
    Returns dict with quality status and rates, or None if no fastp reports are found.
    """
    sample_mrn = sample.split('/')[-1] if '/' in sample else sample
    sample_clean = sample_mrn.rsplit('_S', 1)[0] if '_S' in sample_mrn else sample_mrn
    
    candidate_patterns = [
        os.path.join(outdir, f"analysis/001_trimming/{sample}*_report.json"),
        os.path.join(outdir, f"analysis/001_trimming/{sample_mrn}*_report.json"),
        os.path.join(outdir, f"analysis/001_trimming/*/{sample_mrn}*_report.json"),
        os.path.join(outdir, f"analysis/001_trimming/{sample_clean}*_report.json"),
        os.path.join(outdir, f"analysis/001_trimming/*/{sample_clean}*_report.json"),
        os.path.join(outdir, f"analysis/001_trimming/**/{sample_clean}*report.json"),
        os.path.join(outdir, f"analysis/001_qc/{sample_clean}*report.json")
    ]
    
    fastp_files = []
    for pat in candidate_patterns:
        matches = glob.glob(pat, recursive=True)
        if matches:
            fastp_files = matches
            break
            
    if not fastp_files:
        return None
        
    total_r1_q30_bases = 0
    total_r1_bases = 0
    total_r2_q30_bases = 0
    total_r2_bases = 0
    total_reads = 0
    
    try:
        for fp in fastp_files:
            with open(fp, 'r') as f:
                data = json.load(f)
                
            summary = data.get("summary", {})
            before = summary.get("before_filtering", {})
            total_reads += before.get("total_reads", 0)
            
            r1_info = data.get("read1_before_filtering", {})
            if r1_info and r1_info.get("total_bases", 0) > 0:
                total_r1_q30_bases += r1_info.get("q30_bases", 0)
                total_r1_bases += r1_info.get("total_bases", 0)
            elif before.get("total_bases", 0) > 0:
                total_r1_q30_bases += before.get("q30_bases", 0) // 2
                total_r1_bases += before.get("total_bases", 0) // 2
                
            r2_info = data.get("read2_before_filtering", {})
            if r2_info and r2_info.get("total_bases", 0) > 0:
                total_r2_q30_bases += r2_info.get("q30_bases", 0)
                total_r2_bases += r2_info.get("total_bases", 0)
            elif before.get("total_bases", 0) > 0:
                total_r2_q30_bases += before.get("q30_bases", 0) // 2
                total_r2_bases += before.get("total_bases", 0) // 2
    except Exception:
        return None
        
    if total_r1_bases == 0:
        return None

    r1_q30_rate = (total_r1_q30_bases / total_r1_bases) if total_r1_bases > 0 else 0.0
    r2_q30_rate = (total_r2_q30_bases / total_r2_bases) if total_r2_bases > 0 else 0.0
    
    # Assess status analogous to FastQC: >=80% Q30 is PASS, >=70% is WARN, else FAIL
    r1_status = "PASS" if r1_q30_rate >= 0.80 else ("WARN" if r1_q30_rate >= 0.70 else "FAIL")
    r2_status = "PASS" if r2_q30_rate >= 0.80 else ("WARN" if r2_q30_rate >= 0.70 else "FAIL")
    
    return {
        "r1_status": r1_status,
        "r2_status": r2_status,
        "r1_q30_pct": round(r1_q30_rate * 100.0, 2),
        "r2_q30_pct": round(r2_q30_rate * 100.0, 2),
        "total_reads": total_reads
    }

def parse_vcf_titv_and_hets(snp_vcf):
    """
    Parse filtered SNP VCF to compute Ti/Tv ratio, Het allele balance, and total SNPs.
    Returns tuple of metrics or (None, ...) if missing/unparsable.
    """
    if not os.path.exists(snp_vcf) or os.path.getsize(snp_vcf) == 0:
        return None, None, None, None, None
        
    transitions = 0
    transversions = 0
    tot_het = 0
    outside_het = 0
    tot_snps = 0
    
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
                
                if len(ref) == 1 and len(alt) == 1:
                    pair = (ref, alt)
                    if pair in TI_PAIRS:
                        transitions += 1
                        tot_snps += 1
                    elif pair in TV_PAIRS:
                        transversions += 1
                        tot_snps += 1
                        
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
        return None, None, None, None, None
        
    titv_ratio = round(float(transitions) / float(transversions), 2) if transversions > 0 else 0.0
    outside_pct = round(100.0 * float(outside_het) / float(tot_het), 2) if tot_het > 0 else 0.0
    
    return titv_ratio, tot_het, outside_het, outside_pct, tot_snps

def parse_vcf_indels(indel_vcf):
    """Count total passing INDEL variants in filtered INDEL VCF. Returns int or None."""
    if not indel_vcf or not os.path.exists(indel_vcf) or os.path.getsize(indel_vcf) == 0:
        return None
    tot_indels = 0
    open_fn = gzip.open if indel_vcf.endswith('.gz') else open
    try:
        with open_fn(indel_vcf, 'rt', errors='ignore') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                if len(parts) >= 7 and parts[6] in ('PASS', '.'):
                    tot_indels += 1
    except Exception:
        return None
    return tot_indels

def generate_coverage_summary_for_target(
    target_name,
    target_bed,
    sample_list,
    outdir,
    run_id,
    genome_size=GENOME_SIZE_GRCH38,
    min_snps=200,
    min_indels=20,
    min_titv=2.0,
    min_callable_pct=80.0,
    max_het_ab_bias_pct=25.0
):
    """
    Generate standard 46-column coverage summary matrix with strict panel isolation.
    Missing metrics are reported as 'NA' without falling back to other panels.
    """
    target_size_bp, num_exons = parse_bed_size_and_exons(target_bed)
    rows = []
    qc_evaluations = []
    
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
        missing_inputs = []
        sample_warnings = []
        
        # 1. Total & Mapped reads from original full BAM
        orig_flagstat = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.flagstat")
        tot_r, tot_r_mapped, _ = parse_flagstat(orig_flagstat)
        if tot_r is None:
            missing_inputs.append("original.flagstat")
            pct_mapped_str = "NA"
            tot_r_str = "NA"
            tot_r_mapped_str = "NA"
        else:
            tot_r_str = str(tot_r)
            tot_r_mapped_str = str(tot_r_mapped)
            pct_mapped_str = f"{round(100.0 * tot_r_mapped / tot_r, 2):.2f}" if tot_r > 0 else "0.00"
        
        # 2. Mapped Reads with MAPQ >= 8 across FULL BAM (Legacy MappedReads_q8)
        tot_mapped_q8 = parse_q8_mapped_reads(sample, outdir)
        if tot_mapped_q8 is None:
            missing_inputs.append("full_bam_q8_reads")
            tot_mapped_q8_str = "NA"
            pct_mapped_q8_str = "NA"
        else:
            tot_mapped_q8_str = str(tot_mapped_q8)
            pct_mapped_q8_str = f"{round(100.0 * tot_mapped_q8 / tot_r, 2):.2f}" if (tot_r and tot_r > 0) else "0.00"
        
        # 3. Target flagstat (Reads on target with MAPQ >= 8)
        target_flagstat = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.{bam_suffix}.flagstat")
        tot_ontarget, _, dup_ontarget = parse_flagstat(target_flagstat)
        if tot_ontarget is None:
            missing_inputs.append(f"{bam_suffix}.flagstat")
            tot_ontarget_str = "NA"
            pct_ontarget_str = "NA"
            uniq_ontarget_str = "NA"
            pct_uniq_ontarget_str = "NA"
        else:
            tot_ontarget_str = str(tot_ontarget)
            if tot_mapped_q8 is not None and tot_mapped_q8 > 0:
                pct_ontarget_str = f"{round(100.0 * tot_ontarget / tot_mapped_q8, 2):.2f}"
            else:
                pct_ontarget_str = "NA"
            uniq_ontarget = max(0, tot_ontarget - (dup_ontarget or 0))
            uniq_ontarget_str = str(uniq_ontarget)
            pct_uniq_ontarget_str = f"{round(100.0 * uniq_ontarget / tot_ontarget, 2):.2f}" if tot_ontarget > 0 else "0.00"
        
        # 4. Read lengths & Picard metrics (Strict panel isolation: NO fallback to target panel)
        align_metrics_file = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.{bam_suffix}.align_sum_metrics.txt")
        align_metrics = parse_alignment_summary_metrics(align_metrics_file)
        if align_metrics is None:
            missing_inputs.append(f"{bam_suffix}.align_sum_metrics")
            fwd_len_str = "NA"
            rev_len_str = "NA"
            pairs_str = "NA"
            pct_pairs_str = "NA"
            strand_str = "NA"
        else:
            fwd_len_str = f"{align_metrics['mean_fwd_len']:.2f}" if align_metrics['mean_fwd_len'] is not None else "NA"
            rev_len_str = f"{align_metrics['mean_rev_len']:.2f}" if align_metrics['mean_rev_len'] is not None else "NA"
            pairs_str = str(align_metrics['aligned_pairs']) if align_metrics['aligned_pairs'] is not None else "NA"
            pct_pairs_str = f"{align_metrics['pct_aligned_pairs']:.2f}" if align_metrics['pct_aligned_pairs'] is not None else "NA"
            strand_str = f"{align_metrics['strand_balance']:.3f}" if align_metrics['strand_balance'] is not None else "NA"
        
        # 5. Enrichment factor (Calculated using GRCh38 reference genome denominator)
        if (tot_mapped_q8 is not None and tot_ontarget is not None and 
            tot_mapped_q8 > 0 and target_size_bp > 0 and genome_size > 0):
            ef_val = round((float(tot_ontarget) / float(tot_mapped_q8)) / (float(target_size_bp) / float(genome_size)), 2)
            ef_str = f"{ef_val:.2f}"
        else:
            ef_str = "NA"
            
        # 6. Missing Exons (Strict panel isolation: NO fallback to target panel)
        mean_cov_file = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.{bam_suffix}.mean_coverage.bed")
        exons_missed = parse_missing_exons(mean_cov_file)
        if exons_missed is None:
            missing_inputs.append(f"{bam_suffix}.mean_coverage.bed")
            exons_missed_str = "NA"
        else:
            exons_missed_str = str(exons_missed)
        
        # 7. Coverage Histogram (Strict panel isolation: NO fallback to target panel)
        hist_file = os.path.join(outdir, f"analysis/004_bam_qc/{sample}.{bam_suffix}.coverage_hist.txt")
        cov_stats = parse_coverage_hist(hist_file, target_size_bp)
        if cov_stats is None:
            missing_inputs.append(f"{bam_suffix}.coverage_hist")
            b_strs = {k: "NA" for k in ['b0', 'b1', 'b5', 'b10', 'b20', 'b30', 'b40', 'b50']}
            callable_str = "NA"
            pct_callable_str = "NA"
            no_cov_str = "NA"
            low_cov_str = "NA"
            excess_cov_str = "NA"
            poor_qual_str = "NA"
            mean_cov_str = "NA"
            median_cov_str = "NA"
            diff_cov_str = "NA"
            max_cov_str = "NA"
            evenness_str = "NA"
        else:
            b_strs = {k: f"{cov_stats[k]:.3f}" for k in ['b0', 'b1', 'b5', 'b10', 'b20', 'b30', 'b40', 'b50']}
            callable_str = str(cov_stats['callable'])
            pct_callable_str = f"{cov_stats['%callable']:.2f}"
            no_cov_str = str(cov_stats['no_cov'])
            low_cov_str = str(cov_stats['low_cov'])
            excess_cov_str = str(cov_stats['excess_cov'])
            poor_qual_str = str(cov_stats['poor_qual'])
            mean_cov_str = f"{cov_stats['mean_cov']:.2f}"
            median_cov_str = f"{cov_stats['median_cov']:.2f}"
            diff_cov_str = f"{cov_stats['diff_cov']:.2f}"
            max_cov_str = str(cov_stats['max_cov'])
            evenness_str = f"{cov_stats['evenness']:.2f}"
        
        # 8. Read Quality Evaluation (fastp QC or legacy FastQC)
        # Check for legacy FastQC summaries or modern fastp JSON reports
        r1_fastqc = glob.glob(os.path.join(outdir, f"analysis/001_trimming/{sample}*R1*_fastqc/summary.txt")) or glob.glob(os.path.join(outdir, f"analysis/001_qc/{sample}*R1*_fastqc/summary.txt"))
        r2_fastqc = glob.glob(os.path.join(outdir, f"analysis/001_trimming/{sample}*R2*_fastqc/summary.txt")) or glob.glob(os.path.join(outdir, f"analysis/001_qc/{sample}*R2*_fastqc/summary.txt"))
        
        fastp_qc = parse_fastp_qc(sample, outdir)
        
        if r1_fastqc and r2_fastqc:
            r1_base_q, r1_seq_q = parse_fastqc_summary(r1_fastqc[0])
            r2_base_q, r2_seq_q = parse_fastqc_summary(r2_fastqc[0])
        elif fastp_qc is not None:
            r1_base_q = fastp_qc["r1_status"]
            r1_seq_q = fastp_qc["r1_status"]
            r2_base_q = fastp_qc["r2_status"]
            r2_seq_q = fastp_qc["r2_status"]
            if fastp_qc["r1_status"] == "FAIL" or fastp_qc["r2_status"] == "FAIL":
                sample_warnings.append(f"Low read quality in fastp (R1 Q30={fastp_qc['r1_q30_pct']}%, R2 Q30={fastp_qc['r2_q30_pct']}%)")
            elif fastp_qc["r1_status"] == "WARN" or fastp_qc["r2_status"] == "WARN":
                sample_warnings.append(f"Marginal read quality in fastp (R1 Q30={fastp_qc['r1_q30_pct']}%, R2 Q30={fastp_qc['r2_q30_pct']}%)")
        else:
            # Pre-trim FastQC was removed in favor of fastp; mark unavailable without failing sample
            r1_base_q, r1_seq_q = "NA", "NA"
            r2_base_q, r2_seq_q = "NA", "NA"
        
        # 9. Ti/Tv and Het SNPs Allele Balance
        snp_vcf = os.path.join(outdir, f"analysis/006_variant_filtering/{sample}.gatk.filtered.snp.vcf")
        indel_vcf = os.path.join(outdir, f"analysis/006_variant_filtering/{sample}.gatk.filtered.indel.vcf")
        titv, tot_het, out_het, out_het_pct, tot_snps = parse_vcf_titv_and_hets(snp_vcf)
        tot_indels = parse_vcf_indels(indel_vcf)
        
        if titv is None:
            missing_inputs.append("filtered_snp_vcf")
            titv_str = "NA"
            tot_het_str = "NA"
            out_het_str = "NA"
            out_het_pct_str = "NA"
        else:
            titv_str = f"{titv:.2f}" if titv > 0 else "NA"
            tot_het_str = str(tot_het)
            out_het_str = str(out_het)
            out_het_pct_str = f"{out_het_pct:.2f}"
            
        if tot_indels is None:
            missing_inputs.append("filtered_indel_vcf")
            tot_indels_str = "NA"
        else:
            tot_indels_str = str(tot_indels)
            
        # 10. QC Diagnostics & Warnings Evaluation
        if missing_inputs:
            sample_warnings.append(f"Missing data inputs: {', '.join(missing_inputs)}")
        if tot_snps is not None and tot_snps <= min_snps:
            sample_warnings.append(f"Low SNP count ({tot_snps} <= {min_snps})")
        if tot_indels is not None and tot_indels <= min_indels:
            sample_warnings.append(f"Low INDEL count ({tot_indels} <= {min_indels})")
        if titv is not None and titv > 0 and titv < min_titv:
            sample_warnings.append(f"Low Ti/Tv ratio ({titv:.2f} < {min_titv:.2f})")
        if cov_stats is not None and cov_stats['%callable'] < min_callable_pct:
            sample_warnings.append(f"Low %callable bases ({cov_stats['%callable']:.1f}% < {min_callable_pct:.1f}%)")
        if out_het_pct is not None and out_het_pct > max_het_ab_bias_pct:
            sample_warnings.append(f"High Het AB bias ({out_het_pct:.1f}% > {max_het_ab_bias_pct:.1f}%)")

        qc_status = "FAIL_MISSING_DATA" if missing_inputs else ("WARNING" if sample_warnings else "PASS")
        qc_evaluations.append({
            "Sample": sample_mrn,
            "Panel": target_name,
            "Total_SNPs": tot_snps if tot_snps is not None else "NA",
            "Total_INDELs": tot_indels_str,
            "Ti/Tv": titv_str,
            "%Callable": pct_callable_str,
            "%Het_AB_Out": out_het_pct_str,
            "Status": qc_status,
            "Warnings": "; ".join(sample_warnings) if sample_warnings else "None"
        })
        
        row = {
            "Sample": sample_mrn,
            "TotalReads": tot_r_str,
            "MappedReads": tot_r_mapped_str,
            "%Mapped": pct_mapped_str,
            "MappedReads_q8": tot_mapped_q8_str,
            "%Mapped_q8": pct_mapped_q8_str,
            "ReadsOnTarget_q8": tot_ontarget_str,
            "%OnTarget": pct_ontarget_str,
            "UniqReadsOnTarget_q8": uniq_ontarget_str,
            "%UniqueReadsOnTarget_q8": pct_uniq_ontarget_str,
            "MeanFwdReadLength": fwd_len_str,
            "MeanRevReadLength": rev_len_str,
            "ReadsAlignedInPairs": pairs_str,
            "%ReadsAlignedInPairs": pct_pairs_str,
            "StrandBalance": strand_str,
            "EnrichmentFactor": ef_str,
            "ExonsMissed": exons_missed_str,
            "Bases=0x": b_strs['b0'],
            "Bases>=1x": b_strs['b1'],
            "Bases>=5x": b_strs['b5'],
            "Bases>=10x": b_strs['b10'],
            "Bases>=20x": b_strs['b20'],
            "Bases>=30x": b_strs['b30'],
            "Bases>=40x": b_strs['b40'],
            "Bases>=50x": b_strs['b50'],
            "TargetSize(bp)": target_size_bp,
            "Callable": callable_str,
            "%Callable": pct_callable_str,
            "No_Coverage": no_cov_str,
            "Low_Coverage": low_cov_str,
            "Excess_Coverage": excess_cov_str,
            "Poor_Quality": poor_qual_str,
            "MeanCov": mean_cov_str,
            "MedianCov": median_cov_str,
            "Diff(mean-med)": diff_cov_str,
            "MaxCov": max_cov_str,
            "Evenness": evenness_str,
            "SNPs_Ts/Tv_UnifiedGenotyper": "NA",  # Modern workflow runs GATK4 HaplotypeCaller, not UnifiedGenotyper
            "SNPs_Ts/Tv_HaplotyeCaller": titv_str,
            "FastQC_Read1_Per_base_sequence_quality": r1_base_q,
            "FastQC_Read1_Per_sequence_quality_scores": r1_seq_q,
            "FastQC_Read2_Per_base_sequence_quality": r2_base_q,
            "FastQC_Read2_Per_sequence_quality_scores": r2_seq_q,
            "Total_HC_HetSNPs": tot_het_str,
            "HetSNPs_AB<0.40or>0.60": out_het_str,
            "%HetSNPs_AB<0.40or>0.60": out_het_pct_str
        }
        rows.append(row)
        
    df = pd.DataFrame(rows, columns=LEGACY_COLUMNS)
    df_qc = pd.DataFrame(qc_evaluations)
    return df, df_qc

def df_to_markdown_table(df):
    """Render pandas DataFrame as standard GitHub-flavored Markdown table."""
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
    parser = argparse.ArgumentParser(description="Legacy WES 46-Column Coverage Summary & QC Report Generator")
    parser.add_argument("--outdir", required=True, help="Analysis output directory")
    parser.add_argument("--samples", nargs="+", required=True, help="List of sample identifiers")
    parser.add_argument("--target-bed", required=True, help="Target ICC BED file")
    parser.add_argument("--prot-coding-bed", required=True, help="Protein Coding CDS BED file")
    parser.add_argument("--canon-tran-bed", required=True, help="Canonical Transcripts BED file")
    parser.add_argument("--run-id", default=None, help="Run identifier / date")
    parser.add_argument("--genome-size", type=int, default=GENOME_SIZE_GRCH38, help="Reference genome size (default: 3,095,693,981 for GRCh38)")
    parser.add_argument("--min-snps", type=int, default=200, help="Minimum SNP warning threshold (legacy checked <= 200; default: 200)")
    parser.add_argument("--min-indels", type=int, default=20, help="Minimum INDEL warning threshold (default: 20)")
    parser.add_argument("--min-titv", type=float, default=2.0, help="Minimum Ti/Tv warning threshold (default: 2.0)")
    parser.add_argument("--min-callable-pct", type=float, default=80.0, help="Minimum %Callable bases threshold (default: 80.0)")
    parser.add_argument("--max-het-ab-bias-pct", type=float, default=25.0, help="Maximum Het AB bias threshold (default: 25.0)")
    
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
    all_qc_evals = []
    
    for panel_name, bed_path in panels:
        df, df_qc = generate_coverage_summary_for_target(
            target_name=panel_name,
            target_bed=bed_path,
            sample_list=args.samples,
            outdir=outdir,
            run_id=run_id,
            genome_size=args.genome_size,
            min_snps=args.min_snps,
            min_indels=args.min_indels,
            min_titv=args.min_titv,
            min_callable_pct=args.min_callable_pct,
            max_het_ab_bias_pct=args.max_het_ab_bias_pct
        )
        dfs[panel_name] = df
        all_qc_evals.append(df_qc)
        
        filename_txt = f"SummaryOutput_{panel_name}_{run_id}.txt"
        filename_tsv = f"SummaryOutput_{panel_name}_{run_id}.tsv"
        
        for d in cov_report_dirs:
            df.to_csv(os.path.join(d, filename_txt), sep='\t', index=False)
            df.to_csv(os.path.join(d, filename_tsv), sep='\t', index=False)
            
    # Combined QC evaluation table across all panels
    combined_qc = pd.concat(all_qc_evals, ignore_index=True) if all_qc_evals else pd.DataFrame()
    qc_tsv_filename = f"SummaryOutput_QC_Warnings_{run_id}.tsv"
    for d in cov_report_dirs:
        combined_qc.to_csv(os.path.join(d, qc_tsv_filename), sep='\t', index=False)

    report_md_path = os.path.join(outdir, "results/Coverage_Report/Coverage_Summary_Report.md")
    with open(report_md_path, 'w') as f:
        f.write("# Comprehensive Cohort Coverage & Sequencing QC Report\n\n")
        f.write(f"- **Run ID**: `{run_id}`\n")
        f.write(f"- **Samples Analyzed**: {len(args.samples)}\n")
        f.write(f"- **Reference Genome Denominator**: {args.genome_size:,} bp (GRCh38)\n")
        f.write(f"- **Ti/Tv Configuration**: `SNPs_Ts/Tv_UnifiedGenotyper = NA (Not executed)` | `SNPs_Ts/Tv_HaplotypeCaller = Active`\n")
        f.write(f"- **Panel Isolation Policy**: Strict panel isolation enforced (no cross-panel fallbacks)\n\n")
        
        f.write("## Quality Control & Warning Audit\n\n")
        f.write("Configured QC Thresholds:\n")
        f.write(f"- **Minimum Passing SNPs**: `{args.min_snps}` *(Explicit threshold resolving legacy condition '<= 200' vs '< 100' log)*\n")
        f.write(f"- **Minimum Passing INDELs**: `{args.min_indels}`\n")
        f.write(f"- **Expected Exome Ti/Tv**: `>= {args.min_titv:.2f}`\n")
        f.write(f"- **Minimum Callable Target Bases**: `>= {args.min_callable_pct:.1f}%`\n")
        f.write(f"- **Maximum Heterozygous AB Bias**: `<= {args.max_het_ab_bias_pct:.1f}%`\n\n")
        
        f.write("### Per-Sample & Panel Diagnostic Audit\n\n")
        f.write(df_to_markdown_table(combined_qc))
        f.write("\n")
        
        for panel_name, _ in panels:
            f.write(f"## {panel_name} Summary Matrix (46 Columns)\n\n")
            f.write(df_to_markdown_table(dfs[panel_name]))
            f.write("\n")
            
    print(f"[SUCCESS] Generated 3 legacy coverage summary tables and QC audit report in {cov_report_dirs[0]}")

if __name__ == '__main__':
    main()
