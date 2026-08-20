#!/usr/bin/env python3
"""
Caller Concordance & Comparison Tool for WES Pipeline.
Compares variant calls between GATK HaplotypeCaller and Google DeepVariant.
Generates Venn intersection counts, Jaccard similarity indices, and markdown reports.
"""

import os
import sys
import gzip
import argparse
from collections import defaultdict


def parse_vcf_variants(vcf_path):
    """Parses passed / high-quality variants from a VCF file."""
    variants = set()
    if not os.path.exists(vcf_path):
        return variants
    
    opener = gzip.open if vcf_path.endswith(".gz") else open
    try:
        with opener(vcf_path, "rt", encoding="utf-8", errors="ignore") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\r\n").split("\t")
                if len(parts) < 5:
                    continue
                chrom, pos, ref, alt = parts[0], parts[1], parts[3], parts[4]
                # Filter status
                filter_status = parts[6] if len(parts) > 6 else "."
                if filter_status not in ["PASS", ".", "0"]:
                    continue
                # Split multiallelics
                for a in alt.split(","):
                    variants.add((chrom, pos, ref, a))
    except Exception as e:
        print(f"[WARNING] Error reading {vcf_path}: {e}", file=sys.stderr)
    return variants


def main():
    parser = argparse.ArgumentParser(description="Compare GATK and DeepVariant VCF calls")
    parser.add_argument("--outdir", required=True, help="Analysis output root directory")
    parser.add_argument("--samples", nargs="+", required=True, help="List of sample names")
    parser.add_argument("--report-md", required=True, help="Output Markdown report path")
    parser.add_argument("--summary-tsv", required=True, help="Output TSV matrix path")
    args = parser.parse_args()

    os.makedirs(os.path.dirname(args.report_md), exist_ok=True)
    os.makedirs(os.path.dirname(args.summary_tsv), exist_ok=True)

    results = []
    
    for sample in args.samples:
        gatk_snp = os.path.join(args.outdir, "analysis", "006_variant_filtering", f"{sample}.gatk.filtered.snp.vcf")
        gatk_indel = os.path.join(args.outdir, "analysis", "006_variant_filtering", f"{sample}.gatk.filtered.indel.vcf")
        dv_snp = os.path.join(args.outdir, "analysis", "006_variant_filtering", f"{sample}.deepvariant.filtered.snp.vcf")
        dv_indel = os.path.join(args.outdir, "analysis", "006_variant_filtering", f"{sample}.deepvariant.filtered.indel.vcf")

        gatk_vars = parse_vcf_variants(gatk_snp) | parse_vcf_variants(gatk_indel)
        dv_vars = parse_vcf_variants(dv_snp) | parse_vcf_variants(dv_indel)

        shared = gatk_vars & dv_vars
        gatk_only = gatk_vars - dv_vars
        dv_only = dv_vars - gatk_vars
        union = gatk_vars | dv_vars

        n_gatk = len(gatk_vars)
        n_dv = len(dv_vars)
        n_shared = len(shared)
        n_gatk_only = len(gatk_only)
        n_dv_only = len(dv_only)
        jaccard = (n_shared / len(union) * 100.0) if union else 100.0
        gatk_concordance = (n_shared / n_gatk * 100.0) if n_gatk else 100.0
        dv_concordance = (n_shared / n_dv * 100.0) if n_dv else 100.0

        results.append({
            "sample": sample,
            "gatk_total": n_gatk,
            "dv_total": n_dv,
            "shared": n_shared,
            "gatk_only": n_gatk_only,
            "dv_only": n_dv_only,
            "jaccard_pct": jaccard,
            "gatk_concordance": gatk_concordance,
            "dv_concordance": dv_concordance
        })

    # Write TSV
    with open(args.summary_tsv, "w", encoding="utf-8") as f:
        headers = ["Sample", "GATK_Total", "DeepVariant_Total", "Concordant_Shared", "GATK_Only", "DeepVariant_Only", "Jaccard_Similarity_Pct", "GATK_Concordance_Pct", "DeepVariant_Concordance_Pct"]
        f.write("\t".join(headers) + "\n")
        for r in results:
            f.write(f"{r['sample']}\t{r['gatk_total']}\t{r['dv_total']}\t{r['shared']}\t{r['gatk_only']}\t{r['dv_only']}\t{r['jaccard_pct']:.2f}\t{r['gatk_concordance']:.2f}\t{r['dv_concordance']:.2f}\n")

    # Write Markdown Report
    with open(args.report_md, "w", encoding="utf-8") as f:
        f.write("# Variant Caller Concordance Report: GATK vs. Google DeepVariant\n\n")
        f.write("This report evaluates concordance and discordance between **GATK HaplotypeCaller** and **Google DeepVariant (WES model)** across the cohort.\n\n")
        f.write("## Cohort Summary Table\n\n")
        f.write("| Sample | GATK Calls | DeepVariant Calls | Shared (Concordant) | GATK Only | DeepVariant Only | Jaccard Concordance (%) |\n")
        f.write("| :--- | :---: | :---: | :---: | :---: | :---: | :---: |\n")
        for r in results:
            f.write(f"| **{r['sample']}** | {r['gatk_total']:,} | {r['dv_total']:,} | {r['shared']:,} | {r['gatk_only']:,} | {r['dv_only']:,} | **{r['jaccard_pct']:.1f}%** |\n")
        
        f.write("\n## Clinical Interpretation Guidelines\n")
        f.write("- **Shared Variants**: Highest confidence calls supported by both graph assembly (PairHMM) and deep learning computer vision (CNN).\n")
        f.write("- **DeepVariant-Only Variants**: Often identify challenging Indels in homopolymer runs and repetitive regions where graph assembly may undercall.\n")
        f.write("- **GATK-Only Variants**: Clinically standard calls that may benefit from manual IGV inspection or Sanger confirmation.\n")

    print(f"[INFO] Caller concordance report generated successfully: {args.report_md}")


if __name__ == "__main__":
    main()
