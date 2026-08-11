#!/usr/bin/env python

"""Create annotation-aware sample and cohort variant summaries from VCF outputs."""

from __future__ import annotations

import argparse
import gzip
import html
import json
import re
from collections import Counter
from pathlib import Path

import pandas as pd

ANNOTATION_COLUMNS = ["gene", "impact", "consequence", "clin_sig", "transcript"]
IMPACT_RANK = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "": 0, None: 0}


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


def parse_info_field(info_field: str) -> dict[str, str]:
    info_map = {}
    for item in info_field.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info_map[key] = value
        else:
            info_map[item] = "True"
    return info_map


def parse_annotation_format(line: str, field_id: str) -> list[str]:
    pattern = rf"ID={field_id}.*?Format: ([^\">]+)"
    match = re.search(pattern, line)
    if not match:
        return []
    return [field.strip() for field in match.group(1).split("|")]


def parse_annotation_entries(raw_value: str, annotation_fields: list[str]) -> list[dict[str, str]]:
    entries = []
    if not raw_value or not annotation_fields:
        return entries
    for block in raw_value.split(","):
        values = block.split("|")
        if len(values) < len(annotation_fields):
            values.extend([""] * (len(annotation_fields) - len(values)))
        entries.append(dict(zip(annotation_fields, values)))
    return entries


def annotation_rank(annotation: dict[str, str]) -> tuple[int, int, int, str]:
    impact = annotation.get("IMPACT") or annotation.get("Annotation_Impact") or ""
    clin_sig = (annotation.get("CLIN_SIG") or annotation.get("CLNSIG") or "").lower()
    consequence = annotation.get("Consequence") or annotation.get("Annotation") or ""
    pathogenic_score = 1 if any(token in clin_sig for token in ["pathogenic", "likely_pathogenic"]) else 0
    transcript_support = 1 if annotation.get("Feature") or annotation.get("Feature_ID") else 0
    return (IMPACT_RANK.get(impact, 0), pathogenic_score, transcript_support, consequence)


def normalize_annotation(annotation: dict[str, str]) -> dict[str, str | None]:
    return {
        "gene": annotation.get("SYMBOL") or annotation.get("Gene_Name") or annotation.get("Gene") or None,
        "impact": annotation.get("IMPACT") or annotation.get("Annotation_Impact") or None,
        "consequence": annotation.get("Consequence") or annotation.get("Annotation") or None,
        "clin_sig": annotation.get("CLIN_SIG") or annotation.get("CLNSIG") or None,
        "transcript": annotation.get("Feature") or annotation.get("Feature_ID") or None,
    }


def load_annotation_map(path: str | None) -> tuple[dict[str, dict[str, str | None]], dict[str, object]]:
    if not path:
        return {}, {"enabled": False, "source": None}

    annotation_path = Path(path)
    if not annotation_path.exists():
        return {}, {"enabled": False, "source": str(annotation_path), "missing": True}

    annotation_fields: list[str] = []
    annotation_id: str | None = None
    annotation_map: dict[str, dict[str, str | None]] = {}
    gene_counter: Counter = Counter()
    impact_counter: Counter = Counter()
    clin_sig_counter: Counter = Counter()
    consequence_counter: Counter = Counter()

    with open_maybe_gzip(str(annotation_path)) as handle:
        for line in handle:
            if line.startswith("##INFO=<ID=CSQ"):
                annotation_fields = parse_annotation_format(line, "CSQ")
                annotation_id = "CSQ"
                continue
            if line.startswith("##INFO=<ID=ANN") and not annotation_fields:
                annotation_fields = parse_annotation_format(line, "ANN")
                annotation_id = "ANN"
                continue
            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8 or not annotation_fields or not annotation_id:
                continue
            chrom, pos, variant_id, ref, alt_field, _, _, info_field = fields[:8]
            info_map = parse_info_field(info_field)
            raw_annotations = info_map.get(annotation_id)
            if not raw_annotations:
                continue

            parsed_entries = parse_annotation_entries(raw_annotations, annotation_fields)
            for alt in alt_field.split(","):
                alt = alt.strip()
                if not alt:
                    continue
                alt_entries = [entry for entry in parsed_entries if not entry or entry.get(annotation_fields[0], "") in {"", alt}]
                if not alt_entries:
                    alt_entries = parsed_entries
                best_annotation = normalize_annotation(max(alt_entries, key=annotation_rank)) if alt_entries else {}
                if variant_id and variant_id != ".":
                    best_annotation["id"] = variant_id
                variant_key = f"{chrom}:{int(pos)}:{ref}>{alt}"
                annotation_map[variant_key] = best_annotation

                if best_annotation.get("gene"):
                    gene_counter.update([best_annotation["gene"]])
                if best_annotation.get("impact"):
                    impact_counter.update([best_annotation["impact"]])
                if best_annotation.get("clin_sig"):
                    for value in str(best_annotation["clin_sig"]).split("&"):
                        if value:
                            clin_sig_counter.update([value])
                if best_annotation.get("consequence"):
                    for value in str(best_annotation["consequence"]).split("&"):
                        if value:
                            consequence_counter.update([value])

    return annotation_map, {
        "enabled": True,
        "source": str(annotation_path),
        "field_id": annotation_id,
        "annotated_variant_count": len(annotation_map),
        "impact_counts": dict(impact_counter),
        "clin_sig_counts": dict(clin_sig_counter),
        "consequence_counts": dict(consequence_counter),
        "top_genes": [{"gene": gene, "count": count} for gene, count in gene_counter.most_common(10)],
    }


def read_variants(path: str, source: str) -> list[dict]:
    records = []
    with open_maybe_gzip(path) as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            chrom, pos, variant_id, ref, alt_field, qual, filt, info_field = fields[:8]
            sample_field = fields[9] if len(fields) > 9 else ""
            gt_class = parse_gt(sample_field)
            
            rs_id = None
            if variant_id == ".":
                info_map = parse_info_field(info_field)
                rs_id = info_map.get("RS") or ""
            else:
                rs_id = variant_id
                
            qual_value = None
            if qual not in {".", ""}:
                try:
                    qual_value = float(qual)
                except ValueError:
                    qual_value = None

            is_pass = filt in {"PASS", "."}
            pos_int = int(pos)

            for alt in alt_field.split(","):
                alt = alt.strip()
                if not alt:
                    continue
                kind = variant_kind(ref, alt)
                transition = transition_class(ref, alt) if kind == "snp" else None
                records.append(
                    {
                        "source": source,
                        "chrom": chrom,
                        "pos": pos_int,
                        "id": rs_id,
                        "ref": ref,
                        "alt": alt,
                        "qual": qual_value,
                        "filter": filt,
                        "pass": is_pass,
                        "type": kind,
                        "gt_class": gt_class,
                        "transition_class": transition,
                        "variant_key": f"{chrom}:{pos}:{ref}>{alt}",
                        "length_delta": abs(len(ref) - len(alt)),
                    }
                )
    return records


def empty_sample_summary(sample_name: str, annotation_metadata: dict[str, object]) -> dict[str, object]:
    return {
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
        "annotation": {
            "enabled": bool(annotation_metadata.get("enabled")),
            "source": annotation_metadata.get("source"),
            "annotated_variant_count": 0,
            "impact_counts": {},
            "clin_sig_counts": {},
            "consequence_counts": {},
            "top_genes": [],
        },
    }


def build_sample_summary(
    sample_name: str,
    snp_vcf: str,
    indel_vcf: str,
    top_variants: int,
    top_chromosomes: int,
    annotated_vcf: str | None = None,
) -> tuple[dict[str, object], pd.DataFrame]:
    variants = read_variants(snp_vcf, "snp_vcf") + read_variants(indel_vcf, "indel_vcf")
    df = pd.DataFrame(variants)
    annotation_map, annotation_metadata = load_annotation_map(annotated_vcf)

    if df.empty:
        return empty_sample_summary(sample_name, annotation_metadata), df

    for column in ANNOTATION_COLUMNS:
        df[column] = None

    if annotation_map:
        annotation_df = pd.DataFrame.from_dict(annotation_map, orient="index")
        annotation_df.index.name = "variant_key"
        annotation_df.reset_index(inplace=True)
        df = df.merge(annotation_df, on="variant_key", how="left", suffixes=("", "_annot"))
        for column in ANNOTATION_COLUMNS:
            merged_column = f"{column}_annot"
            if merged_column in df.columns:
                df[column] = df[merged_column].combine_first(df[column])
                df.drop(columns=[merged_column], inplace=True)

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
    top_columns = ["variant_key", "type", "qual", "filter", "gt_class", "gene", "impact", "consequence", "clin_sig"]
    top_pass_variants = sortable_pass.sort_values(["qual_sort", "chrom", "pos"], ascending=[False, True, True]).head(top_variants)[top_columns].to_dict(orient="records")

    impact_counts = {str(key): int(value) for key, value in df["impact"].dropna().value_counts().to_dict().items()}
    consequence_counts = Counter()
    for value in df["consequence"].dropna():
        for part in str(value).split("&"):
            if part:
                consequence_counts.update([part])
    clin_sig_counts = Counter()
    for value in df["clin_sig"].dropna():
        for part in str(value).split("&"):
            if part:
                clin_sig_counts.update([part])
    gene_counts = df["gene"].dropna().value_counts().head(10)

    summary = {
        "sample": sample_name,
        "normalized_sample": normalize_sample_name(sample_name),
        "total_variants": total_variants,
        "pass_variants": int(len(pass_df)),
        "snp_count": int((df["type"] == "snp").sum()),
        "indel_count": int((df["type"] == "indel").sum()),
        "filter_counts": {str(key): int(value) for key, value in df["filter"].value_counts().to_dict().items()},
        "genotype_counts": {str(key): int(value) for key, value in df["gt_class"].dropna().value_counts().to_dict().items()},
        "transition_count": transition_count,
        "transversion_count": transversion_count,
        "titv_ratio": titv_ratio,
        "top_chromosomes": top_chromosome_df.to_dict(orient="records"),
        "top_pass_variants": top_pass_variants,
        "annotation": {
            "enabled": bool(annotation_metadata.get("enabled")),
            "source": annotation_metadata.get("source"),
            "annotated_variant_count": int(df["gene"].notna().sum() or df["impact"].notna().sum()),
            "impact_counts": impact_counts,
            "consequence_counts": dict(consequence_counts),
            "clin_sig_counts": dict(clin_sig_counts),
            "top_genes": [{"gene": gene, "count": int(count)} for gene, count in gene_counts.items()],
        },
    }
    return summary, df


def write_sample_outputs(summary: dict[str, object], df: pd.DataFrame, report_md: str, summary_tsv: str, summary_json: str) -> None:
    ensure_parent(report_md)
    ensure_parent(summary_tsv)
    ensure_parent(summary_json)

    output_df = df.copy()
    if not output_df.empty:
        output_columns = ["variant_key", "chrom", "pos", "type", "filter", "qual", "gt_class", "gene", "impact", "consequence", "clin_sig", "source"]
        output_df = output_df[output_columns]
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

    annotation_summary = summary["annotation"]
    if annotation_summary["enabled"]:
        lines.extend([
            "## Annotation Snapshot",
            "",
            f"- Annotated variants: {annotation_summary['annotated_variant_count']}",
            f"- Annotation source: {annotation_summary['source']}",
            "",
        ])
        if annotation_summary["impact_counts"]:
            lines.extend(["### Impact Counts", ""])
            for impact, count in sorted(annotation_summary["impact_counts"].items(), key=lambda item: (-IMPACT_RANK.get(item[0], 0), item[0])):
                lines.append(f"- {impact}: {count}")
            lines.append("")
        if annotation_summary["top_genes"]:
            lines.extend(["### Top Genes", "", "| Gene | Variants |", "|---|---:|"])
            for row in annotation_summary["top_genes"]:
                lines.append(f"| {row['gene']} | {row['count']} |")
            lines.append("")

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
        lines.extend(["## Top PASS Variants", "", "| Variant | Gene | Impact | Consequence | QUAL | Genotype |", "|---|---|---|---|---:|---|"])
        for row in summary["top_pass_variants"]:
            qual = row["qual"] if row["qual"] is not None else "NA"
            lines.append(
                f"| {row['variant_key']} | {row.get('gene') or 'NA'} | {row.get('impact') or 'NA'} | {row.get('consequence') or 'NA'} | {qual} | {row.get('gt_class') or 'NA'} |"
            )
        lines.append("")

    with open(report_md, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def build_cohort_summary(input_paths: list[str]) -> tuple[dict[str, object], pd.DataFrame]:
    summaries = []
    for path in input_paths:
        with open(path, "r", encoding="utf-8") as handle:
            summaries.append(json.load(handle))

    rows = []
    filter_counter = Counter()
    genotype_counter = Counter()
    impact_counter = Counter()
    clin_sig_counter = Counter()
    gene_counter = Counter()

    for summary in summaries:
        annotation = summary.get("annotation", {})
        high_impact_count = int(annotation.get("impact_counts", {}).get("HIGH", 0))
        rows.append(
            {
                "sample": summary["sample"],
                "total_variants": summary["total_variants"],
                "pass_variants": summary["pass_variants"],
                "snp_count": summary["snp_count"],
                "indel_count": summary["indel_count"],
                "titv_ratio": summary["titv_ratio"],
                "annotated_variants": int(annotation.get("annotated_variant_count", 0)),
                "high_impact_variants": high_impact_count,
            }
        )
        filter_counter.update(summary.get("filter_counts", {}))
        genotype_counter.update(summary.get("genotype_counts", {}))
        impact_counter.update(annotation.get("impact_counts", {}))
        clin_sig_counter.update(annotation.get("clin_sig_counts", {}))
        gene_counter.update({row["gene"]: row["count"] for row in annotation.get("top_genes", [])})

    if rows:
        cohort_df = pd.DataFrame(rows).sort_values(
            ["pass_variants", "annotated_variants", "sample"], ascending=[False, False, True]
        )
    else:
        cohort_df = pd.DataFrame(columns=["sample", "total_variants", "pass_variants", "snp_count", "indel_count", "titv_ratio", "annotated_variants", "high_impact_variants"])
    titv_values = cohort_df["titv_ratio"].dropna() if not cohort_df.empty else pd.Series(dtype=float)
    cohort_summary = {
        "samples": summaries,
        "sample_count": len(summaries),
        "total_variants": int(cohort_df["total_variants"].sum()) if not cohort_df.empty else 0,
        "pass_variants": int(cohort_df["pass_variants"].sum()) if not cohort_df.empty else 0,
        "annotated_variants": int(cohort_df["annotated_variants"].sum()) if not cohort_df.empty else 0,
        "high_impact_variants": int(cohort_df["high_impact_variants"].sum()) if not cohort_df.empty else 0,
        "mean_titv_ratio": round(titv_values.mean(), 3) if not titv_values.empty else None,
        "filter_counts": dict(filter_counter),
        "genotype_counts": dict(genotype_counter),
        "impact_counts": dict(impact_counter),
        "clin_sig_counts": dict(clin_sig_counter),
        "top_genes": [{"gene": gene, "count": count} for gene, count in gene_counter.most_common(15)],
    }
    return cohort_summary, cohort_df


def format_metric(value: object) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        return f"{value:.3f}".rstrip("0").rstrip(".")
    return str(value)


def render_counter_list(counter_map: dict[str, int], empty_label: str) -> str:
    if not counter_map:
        return f"<li>{html.escape(empty_label)}</li>"
    items = []
    for key, value in sorted(counter_map.items(), key=lambda item: (-item[1], item[0])):
        items.append(f"<li><span>{html.escape(str(key))}</span><strong>{value}</strong></li>")
    return "".join(items)


def render_dashboard(summary: dict[str, object], cohort_df: pd.DataFrame, dashboard_html: str) -> None:
    ensure_parent(dashboard_html)
    sample_rows = []
    for _, row in cohort_df.iterrows():
        sample_rows.append(
            "".join(
                [
                    "<tr>",
                    f"<td>{html.escape(str(row['sample']))}</td>",
                    f"<td>{int(row['total_variants'])}</td>",
                    f"<td>{int(row['pass_variants'])}</td>",
                    f"<td>{int(row['annotated_variants'])}</td>",
                    f"<td>{int(row['high_impact_variants'])}</td>",
                    f"<td>{html.escape(format_metric(row['titv_ratio']) if pd.notna(row['titv_ratio']) else 'NA')}</td>",
                    "</tr>",
                ]
            )
        )

    gene_rows = []
    for row in summary.get("top_genes", []):
        gene_rows.append(f"<tr><td>{html.escape(str(row['gene']))}</td><td>{int(row['count'])}</td></tr>")

    html_output = f"""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>WES Cohort Variant Dashboard</title>
  <style>
    :root {{
      --paper: #f4efe6;
      --ink: #18222c;
      --accent: #a63d2f;
      --accent-2: #2f6b7d;
      --panel: rgba(255, 251, 245, 0.82);
      --line: rgba(24, 34, 44, 0.14);
      --muted: #59626c;
      --shadow: 0 18px 40px rgba(24, 34, 44, 0.09);
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      color: var(--ink);
      background:
        radial-gradient(circle at top left, rgba(166, 61, 47, 0.18), transparent 28%),
        radial-gradient(circle at top right, rgba(47, 107, 125, 0.2), transparent 24%),
        linear-gradient(160deg, #efe7da 0%, var(--paper) 55%, #e7dfd2 100%);
      font-family: Baskerville, "Iowan Old Style", "Palatino Linotype", serif;
      min-height: 100vh;
    }}
    .shell {{
      max-width: 1240px;
      margin: 0 auto;
      padding: 48px 24px 64px;
    }}
    .hero {{
      display: grid;
      grid-template-columns: 1.4fr 0.9fr;
      gap: 24px;
      align-items: end;
      margin-bottom: 28px;
    }}
    .headline {{
      padding: 28px 30px;
      background: linear-gradient(145deg, rgba(255,255,255,0.75), rgba(255,247,238,0.62));
      border: 1px solid var(--line);
      box-shadow: var(--shadow);
      position: relative;
      overflow: hidden;
    }}
    .headline::after {{
      content: "";
      position: absolute;
      right: -72px;
      top: -72px;
      width: 180px;
      height: 180px;
      border-radius: 50%;
      background: rgba(166, 61, 47, 0.09);
    }}
    .eyebrow {{
      font-family: "Avenir Next Condensed", "Franklin Gothic Medium", sans-serif;
      font-size: 0.82rem;
      letter-spacing: 0.18em;
      text-transform: uppercase;
      color: var(--accent-2);
      margin: 0 0 12px;
    }}
    h1 {{
      margin: 0;
      font-size: clamp(2.4rem, 4.2vw, 4.8rem);
      line-height: 0.94;
      letter-spacing: -0.04em;
      max-width: 11ch;
    }}
    .subtitle {{
      margin: 14px 0 0;
      color: var(--muted);
      max-width: 54ch;
      line-height: 1.5;
    }}
    .metrics {{
      display: grid;
      grid-template-columns: repeat(2, minmax(0, 1fr));
      gap: 14px;
    }}
    .metric {{
      background: var(--panel);
      border: 1px solid var(--line);
      box-shadow: var(--shadow);
      padding: 18px 20px;
    }}
    .metric span {{
      display: block;
      font-family: "Avenir Next Condensed", "Franklin Gothic Medium", sans-serif;
      font-size: 0.74rem;
      letter-spacing: 0.18em;
      text-transform: uppercase;
      color: var(--muted);
      margin-bottom: 10px;
    }}
    .metric strong {{
      font-size: 2rem;
      line-height: 1;
    }}
    .grid {{
      display: grid;
      grid-template-columns: 1.2fr 0.8fr;
      gap: 24px;
    }}
    .panel {{
      background: var(--panel);
      border: 1px solid var(--line);
      box-shadow: var(--shadow);
      padding: 24px;
    }}
    .panel h2 {{
      margin: 0 0 18px;
      font-size: 1.3rem;
      letter-spacing: -0.02em;
    }}
    table {{ width: 100%; border-collapse: collapse; }}
    th, td {{ padding: 11px 10px; text-align: left; border-bottom: 1px solid var(--line); }}
    th {{
      font-family: "Avenir Next Condensed", "Franklin Gothic Medium", sans-serif;
      text-transform: uppercase;
      letter-spacing: 0.12em;
      font-size: 0.72rem;
      color: var(--muted);
    }}
    td {{ font-size: 0.98rem; }}
    .lists {{ display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 24px; margin-top: 24px; }}
    .list-box ul {{ list-style: none; padding: 0; margin: 0; }}
    .list-box li {{ display: flex; justify-content: space-between; gap: 16px; padding: 9px 0; border-bottom: 1px solid var(--line); }}
    .ribbon {{
      display: inline-block;
      font-family: "Avenir Next Condensed", "Franklin Gothic Medium", sans-serif;
      letter-spacing: 0.14em;
      text-transform: uppercase;
      font-size: 0.72rem;
      padding: 7px 10px;
      background: rgba(166, 61, 47, 0.1);
      color: var(--accent);
      margin-bottom: 12px;
    }}
    @media (max-width: 960px) {{
      .hero, .grid, .lists {{ grid-template-columns: 1fr; }}
      .metrics {{ grid-template-columns: repeat(2, minmax(0, 1fr)); }}
    }}
    @media (max-width: 640px) {{
      .shell {{ padding: 24px 14px 40px; }}
      .headline, .panel, .metric {{ padding: 18px; }}
      .metrics {{ grid-template-columns: 1fr; }}
      th, td {{ padding-left: 6px; padding-right: 6px; }}
    }}
  </style>
</head>
<body>
  <main class="shell">
    <section class="hero">
      <div class="headline">
        <p class="eyebrow">Whole Exome Cohort Dashboard</p>
        <h1>Variant Burden, Annotation, and Sample Triage</h1>
        <p class="subtitle">A single artifact for filtered variant burden and optional VEP-derived annotations. When annotation files are present, the dashboard surfaces impact tiers, gene hotspots, and clinically tagged findings without making VEP a hard dependency.</p>
      </div>
      <div class="metrics">
        <article class="metric"><span>Samples</span><strong>{summary['sample_count']}</strong></article>
        <article class="metric"><span>Total Variants</span><strong>{summary['total_variants']}</strong></article>
        <article class="metric"><span>Annotated Variants</span><strong>{summary['annotated_variants']}</strong></article>
        <article class="metric"><span>High Impact</span><strong>{summary['high_impact_variants']}</strong></article>
      </div>
    </section>

    <section class="grid">
      <article class="panel">
        <span class="ribbon">Sample Ranking</span>
        <h2>Per-Sample Summary</h2>
        <table>
          <thead>
            <tr>
              <th>Sample</th>
              <th>Total</th>
              <th>PASS</th>
              <th>Annotated</th>
              <th>High Impact</th>
              <th>Ti/Tv</th>
            </tr>
          </thead>
          <tbody>{''.join(sample_rows) or '<tr><td colspan="6">No sample summaries available.</td></tr>'}</tbody>
        </table>
      </article>

      <article class="panel">
        <span class="ribbon">Cohort Health</span>
        <h2>Pipeline Snapshot</h2>
        <div class="lists">
          <div class="list-box">
            <h3>Impact Counts</h3>
            <ul>{render_counter_list(summary.get('impact_counts', {}), 'No annotation impacts found')}</ul>
          </div>
          <div class="list-box">
            <h3>Filter Outcomes</h3>
            <ul>{render_counter_list(summary.get('filter_counts', {}), 'No filter counts found')}</ul>
          </div>
        </div>
      </article>
    </section>

    <section class="grid" style="margin-top: 24px;">
      <article class="panel">
        <span class="ribbon">Gene Hotspots</span>
        <h2>Top Cohort Genes</h2>
        <table>
          <thead><tr><th>Gene</th><th>Variants</th></tr></thead>
          <tbody>{''.join(gene_rows) or '<tr><td colspan="2">No annotated genes found.</td></tr>'}</tbody>
        </table>
      </article>

      <article class="panel">
        <span class="ribbon">Clinical Signals</span>
        <h2>ClinVar-like Labels</h2>
        <ul>{render_counter_list(summary.get('clin_sig_counts', {}), 'No clinical labels captured')}</ul>
        <p class="subtitle" style="margin-top: 16px;">Mean cohort Ti/Tv: <strong>{html.escape(format_metric(summary.get('mean_titv_ratio')))}</strong></p>
      </article>
    </section>
  </main>
</body>
</html>
"""

    with open(dashboard_html, "w", encoding="utf-8") as handle:
        handle.write(html_output)


def write_cohort_outputs(summary: dict[str, object], cohort_df: pd.DataFrame, report_md: str, summary_tsv: str, summary_json: str, dashboard_html: str | None = None) -> None:
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
        f"- Annotated variants: {summary['annotated_variants']}",
        f"- High-impact variants: {summary['high_impact_variants']}",
        f"- Mean Ti/Tv ratio: {summary['mean_titv_ratio'] if summary['mean_titv_ratio'] is not None else 'NA'}",
        "",
    ]

    if not cohort_df.empty:
        lines.extend([
            "## Sample Table",
            "",
            "| Sample | Total | PASS | Annotated | High impact | SNPs | Indels | Ti/Tv |",
            "|---|---:|---:|---:|---:|---:|---:|---:|",
        ])
        for _, row in cohort_df.iterrows():
            titv = row["titv_ratio"] if pd.notna(row["titv_ratio"]) else "NA"
            lines.append(
                f"| {row['sample']} | {int(row['total_variants'])} | {int(row['pass_variants'])} | {int(row['annotated_variants'])} | {int(row['high_impact_variants'])} | {int(row['snp_count'])} | {int(row['indel_count'])} | {titv} |"
            )
        lines.append("")

    if summary["impact_counts"]:
        lines.extend(["## Cohort Impact Counts", ""])
        for impact, count in sorted(summary["impact_counts"].items(), key=lambda item: (-IMPACT_RANK.get(item[0], 0), item[0])):
            lines.append(f"- {impact}: {count}")
        lines.append("")

    if summary["top_genes"]:
        lines.extend(["## Top Genes", "", "| Gene | Variants |", "|---|---:|"])
        for row in summary["top_genes"]:
            lines.append(f"| {row['gene']} | {row['count']} |")
        lines.append("")

    if summary["filter_counts"]:
        lines.extend(["## Cohort Filter Outcomes", ""])
        for filter_name, count in sorted(summary["filter_counts"].items(), key=lambda item: (-item[1], item[0])):
            lines.append(f"- {filter_name}: {count}")
        lines.append("")

    with open(report_md, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))

    if dashboard_html:
        render_dashboard(summary, cohort_df, dashboard_html)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    sample_parser = subparsers.add_parser("sample")
    sample_parser.add_argument("--sample-name", required=True)
    sample_parser.add_argument("--snp-vcf", required=True)
    sample_parser.add_argument("--indel-vcf", required=True)
    sample_parser.add_argument("--annotated-vcf")
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
    cohort_parser.add_argument("--dashboard-html")

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
            annotated_vcf=args.annotated_vcf,
        )
        write_sample_outputs(summary, df, args.report_md, args.summary_tsv, args.summary_json)
        return

    cohort_summary, cohort_df = build_cohort_summary(args.inputs)
    write_cohort_outputs(cohort_summary, cohort_df, args.report_md, args.summary_tsv, args.summary_json, args.dashboard_html)


if __name__ == "__main__":
    main()