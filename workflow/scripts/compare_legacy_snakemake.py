#!/usr/bin/env python3
"""
Comprehensive Benchmarking & Head-to-Head Comparison:
Snakemake Modern Pipeline (GATK4 & Google DeepVariant) vs.
Legacy Parallel Pipeline (GATK 3.2 UnifiedGenotyper & GATK 3.2 HaplotypeCaller).
Evaluated on NIST GIAB HG001 / NA12878 standard.
"""

import os
import sys
sys.path.insert(0, os.path.abspath("."))
import json
from workflow.scripts.benchmark_giab import load_bed_intervals, parse_vcf_fast, evaluate_caller

truth_vcf = "resources/giab/HG001/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"
callable_bed = "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/012_benchmark/giab/HG001_callable_min10x_giab_intersect.bed"
panel_bed = "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/012_benchmark/giab/HG001_HG001_S1_target_giab_intersect.bed"

callers_files = {
    "Legacy UnifiedGenotyper (GATK 3.2)": {
        "pipeline": "Legacy Pipeline (GATK 3.2)",
        "caller": "UnifiedGenotyper",
        "snp": "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/legacy/HG001_S1.legacy_ug.filtered.snp.vcf",
        "indel": "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/legacy/HG001_S1.legacy_ug.filtered.indel.vcf"
    },
    "Legacy HaplotypeCaller (GATK 3.2)": {
        "pipeline": "Legacy Pipeline (GATK 3.2)",
        "caller": "HaplotypeCaller v3.2",
        "snp": "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/legacy/HG001_S1.legacy_hc.filtered.snp.vcf",
        "indel": "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/legacy/HG001_S1.legacy_hc.filtered.indel.vcf"
    },
    "Snakemake GATK4 HaplotypeCaller": {
        "pipeline": "Snakemake Modern (GATK 4.6)",
        "caller": "HaplotypeCaller v4.6",
        "snp": "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/006_variant_filtering/HG001/HG001_S1.gatk.filtered.snp.vcf",
        "indel": "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/006_variant_filtering/HG001/HG001_S1.gatk.filtered.indel.vcf"
    },
    "Snakemake Google DeepVariant": {
        "pipeline": "Snakemake Modern (DeepVariant 1.6)",
        "caller": "DeepVariant CNN",
        "snp": "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/006_variant_filtering/HG001/HG001_S1.deepvariant.filtered.snp.vcf",
        "indel": "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/006_variant_filtering/HG001/HG001_S1.deepvariant.filtered.indel.vcf"
    }
}

def evaluate_all(bed_path, title):
    intervals = load_bed_intervals(bed_path)
    truth_snps, truth_indels = parse_vcf_fast(truth_vcf, eval_bed=bed_path, intervals=intervals, pass_only=True)
    
    print("\n" + "=" * 96)
    print(f" {title.upper()}")
    print(f" Truth Standards: {len(truth_snps)} SNPs, {len(truth_indels)} INDELs | Evaluation Region: {bed_path}")
    print("=" * 96)
    print(f"{'Pipeline & Caller Engine':<38} | {'Type':<6} | {'Recall':<8} | {'Precision':<10} | {'F1-Score':<8} | {'TP/FP/FN':<10} | {'Ti/Tv'}")
    print("-" * 96)
    
    results = {}
    for name, info in callers_files.items():
        q_snps, _ = parse_vcf_fast(info["snp"], eval_bed=bed_path, intervals=intervals, pass_only=True)
        _, q_indels = parse_vcf_fast(info["indel"], eval_bed=bed_path, intervals=intervals, pass_only=True)
        res = evaluate_caller(truth_snps, truth_indels, q_snps, q_indels)
        results[name] = res
        
        s = res["snps"]
        i = res["indels"]
        a = res["all"]
        print(f"{name:<38} | {'SNP':<6} | {s['recall']:>6.2f}% | {s['precision']:>8.2f}% | {s['f1']:>6.2f}% | {s['tp']}/{s['fp']}/{s['fn']:<4} | {res['titv_ratio']}")
        print(f"{'':<38} | {'INDEL':<6} | {i['recall']:>6.2f}% | {i['precision']:>8.2f}% | {i['f1']:>6.2f}% | {i['tp']}/{i['fp']}/{i['fn']:<4} | --")
        print(f"{'':<38} | {'ALL':<6} | {a['recall']:>6.2f}% | {a['precision']:>8.2f}% | {a['f1']:>6.2f}% | {a['tp']}/{a['fp']}/{a['fn']:<4} | GT: {res['genotype_concordance']}%")
        print("-" * 96)
    return results

callable_results = evaluate_all(callable_bed, "1. Callable Regions (Depth >= 10X Intersect — 564,620 bp)")
panel_results = evaluate_all(panel_bed, "2. Full 169-Gene Custom Capture Panel Intersect (Includes Low-Depth Overhangs)")

# Export summary JSON and Markdown report
out_json = "/mnt/bucket/VC/Snakemake_Analysis/GIAB_HG001_Full_Benchmark/analysis/012_benchmark/giab/pipeline_comparison_benchmark.json"
with open(out_json, "w") as f:
    json.dump({"callable_regions": callable_results, "full_panel_regions": panel_results}, f, indent=2)

print(f"\n[SUCCESS] Comparison results saved to: {out_json}\n")
