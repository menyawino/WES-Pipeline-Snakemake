#!/usr/bin/env python

"""Create sample-level and cohort-level variant summaries from filtered VCFs."""

from __future__ import annotations

import argparse
import gzip
import json
from collections import Counter
from pathlib import Path

import pandas as pd


def open_maybe_gzip(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")


def ensure_parent(path: str) -> None:
    Path(path).parent.mkdir(parents=True, exist_ok=True)


def normalize_sample_name(sample_name: str) -> str:
    return sample_name.replace("/", "_")


def parse_gt(sample_field: str) -> str | None:
    if not sample_field:
        return None
    genotype = sample_field.split(":", 1)[0].replace("|", "/")
    if genotype in {"0/1", "1/0"}:
        return "het"
    if genotype in {"1/1", "2/2", "3/3"}:
        return "hom_alt"
    if genotype == "0/0":
        return "hom_ref"
    if genotype == "./.":
        return "missing"
    if "/" in genotype:
        alleles = genotype.split("/")
        if "." in alleles:
            return "missing"
        if len(set(alleles)) == 1 and alleles[0] != "0":
            return "hom_alt"
        if "0" in alleles and len(set(alleles)) > 1:
            return "het"
        return "complex"
    return None


def variant_kind(ref: str, alt: str) -> str:
    if len(ref) == 1 and len(alt) == 1:
        return "snp"
    return "indel"


def transition_class(ref: str, alt: str) -> str | None:
    pair = (ref.upper(), alt.upper())
    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}
    transversions = {
        ("A", "C"), ("C", "A"), ("A", "T"), ("T", "A"),
        ("C", "G"), ("G", "C"), ("G", "T"), ("T", "G"),
    }
    if pair in transitions:
        return "transition"
    if pair in transversions:
        return "transversion"
    return None


def read_variants(path: str, source: str) -> list[dict]:
    records = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, variant_id, ref, alt_field, qual, filt, info = fields[:8]
            sample_field = fields[9] if len(fields) > 9 else ""
            for alt in alt_field.split(","):
                alt = alt.strip()
                if not alt:
                    continue
                kind = variant_kind(ref, alt)
                transition = transition_class(ref, alt) if kind == "snp" else None
                qual_value = None
                if qual not in {".", ""}:
                    try:
                        qual_value = float(qual)
                    except ValueError:
                        qual_value = None
                records.append(
                    {
                        "source": source,
                        "chrom": chrom,
                        "pos": int(pos),
                        "id": variant_id if variant_id != "." else "",
                        "ref": ref,
                        "alt": alt,
                        "qual": qual_value,
                        "filter": filt,
                        "pass": filt in {"PASS", "."},
                        "type": kind,
                        "gt_class": parse_gt(sample_field),
                        "transition_class": transition,
                        "variant_key": f"{chrom}:{pos}:{ref}>{alt}",
                        "length_delta": abs(len(ref) - len(alt)),
                    }
                )
    return records


def build_sample_summary(sample_name: str, snp_vcf: str, indel_vcf: str, top_variants: int, top_chromosomes: int) -> tuple[dict, pd.DataFrame]:
    variants = read_variants(snp_vcf, "snp_vcf") + read_variants(indel_vcf, "indel_vcf")
    df = pd.DataFrame(variants)

    if df.empty:
        summary = {
            "sample": sample_name,
            "normalized_sample": normalize_sample_name(sample_name),
            "total_variants": 0,
            "pass_variants": 0,
            "snp_count": 0,
            "indel_count": 0,
            "filter_counts": {},
            "genotype_counts": {},
            "transition_count": 0,
            "transversion_count": 0,
            "titv_ratio": None,
            "top_chromosomes": [],
            "top_pass_variants": [],
        }
        return summary, df

    total_variants = int(len(df))
    pass_df = df[df["pass"]].copy()
    transition_count = int((df["transition_class"] == "transition").sum())
    transversion_count = int((df["transition_class"] == "transversion").sum())
    titv_ratio = round(transition_count / transversion_count, 3) if transversion_count else None

    top_chromosome_df = (
        pass_df.groupby("chrom").size().sort_values(ascending=False).head(top_chromosomes).reset_index(name="pass_variants")
    )

    sortable_pass = pass_df.copy()
    sortable_pass["qual_sort"] = sortable_pass["qual"].fillna(-1)
    top_pass_variants = (
        sortable_pass.sort_values(["qual_sort", "chrom", "pos"], ascending=[False, True, True])
        .head(top_variants)[["variant_key", "type", "qual", "filter", "gt_class"]]
        .to_dict(orient="records")
    )

    summary = {
        "sample": sample_name,
        "normalized_sample": normalize_sample_name(sample_name),
        "total_variants": total_variants,
        "pass_variants": int(len(pass_df)),
        "snp_count": int((df["type"] == "snp").sum()),
        "indel_count": int((df["type"] == "indel").sum()),
        "filter_counts": {str(k): int(v) for k, v in df["filter"].value_counts().to_dict().items()},
        "genotype_counts": {str(k): int(v) for k, v in df["gt_class"].dropna().value_counts().to_dict().items()},
        "transition_count": transition_count,
        "transversion_count": transversion_count,
        "titv_ratio": titv_ratio,
        "top_chromosomes": top_chromosome_df.to_dict(orient="records"),
        "top_pass_variants": top_pass_variants,
    }
    return summary, df


def write_sample_outputs(summary: dict, df: pd.DataFrame, report_md: str, summary_tsv: str, summary_json: str) -> None:
    ensure_parent(report_md)
    ensure_parent(summary_tsv)
    ensure_parent(summary_json)

    output_df = df.copy()
    if not output_df.empty:
        output_df = output_df[["variant_key", "chrom", "pos", "type", "filter", "qual", "gt_class", "source"]]
    output_df.to_csv(summary_tsv, sep="\t", index=False)

    with open(summary_json, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    lines = [
        f"# Variant Summary: {summary['sample']}",
        "",
        "## Overview",
        "",
        f"- Total variants: {summary['total_variants']}",
        f"- PASS variants: {summary['pass_variants']}",
        f"- SNPs: {summary['snp_count']}",
        f"- Indels: {summary['indel_count']}",
        f"- Ti/Tv ratio: {summary['titv_ratio'] if summary['titv_ratio'] is not None else 'NA'}",
        "",
    ]

    if summary["genotype_counts"]:
        lines.extend(["## Genotype Classes", ""])
        for genotype, count in sorted(summary["genotype_counts"].items()):
            lines.append(f"- {genotype}: {count}")
        lines.append("")

    if summary["filter_counts"]:
        lines.extend(["## Filter Outcomes", ""])
        for filter_name, count in sorted(summary["filter_counts"].items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- {filter_name}: {count}")
        lines.append("")

    if summary["top_chromosomes"]:
        lines.extend(["## Top Chromosomes By PASS Burden", "", "| Chromosome | PASS variants |", "|---|---:|"])
        for row in summary["top_chromosomes"]:
            lines.append(f"| {row['chrom']} | {row['pass_variants']} |")
        lines.append("")

    if summary["top_pass_variants"]:
        lines.extend(["## Top PASS Variants", "", "| Variant | Type | QUAL | Genotype |", "|---|---|---:|---|"])
        for row in summary["top_pass_variants"]:
            qual = row["qual"] if row["qual"] is not None else "NA"
            lines.append(f"| {row['variant_key']} | {row['type']} | {qual} | {row.get('gt_class') or 'NA'} |")
        lines.append("")

    with open(report_md, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def build_cohort_summary(input_paths: list[str]) -> tuple[list[dict], pd.DataFrame]:
    summaries = []
    for path in input_paths:
        with open(path, "r", encoding="utf-8") as handle:
            summaries.append(json.load(handle))

    rows = []
    filter_counter = Counter()
    genotype_counter = Counter()
    for summary in summaries:
        rows.append(
            {
                "sample": summary["sample"],
                "total_variants": summary["total_variants"],
                "pass_variants": summary["pass_variants"],
                "snp_count": summary["snp_count"],
                "indel_count": summary["indel_count"],
                "titv_ratio": summary["titv_ratio"],
            }
        )
        filter_counter.update(summary.get("filter_counts", {}))
        genotype_counter.update(summary.get("genotype_counts", {}))

    cohort_df = pd.DataFrame(rows).sort_values(["pass_variants", "total_variants", "sample"], ascending=[False, False, True])
    titv_values = cohort_df["titv_ratio"].dropna() if not cohort_df.empty else pd.Series(dtype=float)
    cohort_summary = {
        "samples": summaries,
        "sample_count": len(summaries),
        "total_variants": int(cohort_df["total_variants"].sum()) if not cohort_df.empty else 0,
        "pass_variants": int(cohort_df["pass_variants"].sum()) if not cohort_df.empty else 0,
        "mean_titv_ratio": round(titv_values.mean(), 3) if not titv_values.empty else None,
        "filter_counts": dict(filter_counter),
        "genotype_counts": dict(genotype_counter),
    }
    return [cohort_summary], cohort_df


def write_cohort_outputs(summary: dict, cohort_df: pd.DataFrame, report_md: str, summary_tsv: str, summary_json: str) -> None:
    ensure_parent(report_md)
    ensure_parent(summary_tsv)
    ensure_parent(summary_json)

    cohort_df.to_csv(summary_tsv, sep="\t", index=False)

    with open(summary_json, "w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)

    lines = [
        "# Cohort Variant Summary",
        "",
        "## Overview",
        "",
        f"- Samples: {summary['sample_count']}",
        f"- Total variants: {summary['total_variants']}",
        f"- PASS variants: {summary['pass_variants']}",
        f"- Mean Ti/Tv ratio: {summary['mean_titv_ratio'] if summary['mean_titv_ratio'] is not None else 'NA'}",
        "",
    ]

    if not cohort_df.empty:
        lines.extend(["## Sample Table", "", "| Sample | Total | PASS | SNPs | Indels | Ti/Tv |", "|---|---:|---:|---:|---:|---:|"])
        for _, row in cohort_df.iterrows():
            titv = row["titv_ratio"] if pd.notna(row["titv_ratio"]) else "NA"
            lines.append(f"| {row['sample']} | {int(row['total_variants'])} | {int(row['pass_variants'])} | {int(row['snp_count'])} | {int(row['indel_count'])} | {titv} |")
        lines.append("")

    if summary["filter_counts"]:
        lines.extend(["## Cohort Filter Outcomes", ""])
        for filter_name, count in sorted(summary["filter_counts"].items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- {filter_name}: {count}")
        lines.append("")

    with open(report_md, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    sample_parser = subparsers.add_parser("sample")
    sample_parser.add_argument("--sample-name", required=True)
    sample_parser.add_argument("--snp-vcf", required=True)
    sample_parser.add_argument("--indel-vcf", required=True)
    sample_parser.add_argument("--report-md", required=True)
    sample_parser.add_argument("--summary-tsv", required=True)
    sample_parser.add_argument("--summary-json", required=True)
    sample_parser.add_argument("--top-variants", type=int, default=10)
    sample_parser.add_argument("--top-chromosomes", type=int, default=10)

    cohort_parser = subparsers.add_parser("cohort")
    cohort_parser.add_argument("--inputs", nargs="+", required=True)
    cohort_parser.add_argument("--report-md", required=True)
    cohort_parser.add_argument("--summary-tsv", required=True)
    cohort_parser.add_argument("--summary-json", required=True)

    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if args.command == "sample":
        summary, df = build_sample_summary(
            sample_name=args.sample_name,
            snp_vcf=args.snp_vcf,
            indel_vcf=args.indel_vcf,
            top_variants=args.top_variants,
            top_chromosomes=args.top_chromosomes,
        )
        write_sample_outputs(summary, df, args.report_md, args.summary_tsv, args.summary_json)
        return

    cohort_summaries, cohort_df = build_cohort_summary(args.inputs)
    write_cohort_outputs(cohort_summaries[0], cohort_df, args.report_md, args.summary_tsv, args.summary_json)


if __name__ == "__main__":
    main()