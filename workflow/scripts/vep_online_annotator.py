#!/usr/bin/env python3
"""Ensembl VEP Online Annotator with integrated GeneBe ACMG Classification Plugin.

Annotates VCF variants with high-confidence Ensembl transcript annotations,
consequences, functional impacts, HGVS, SIFT, PolyPhen, ClinVar, and automated
ACMG/AMP 2015 5-tier classification using the official GeneBe algorithm in a single unified command.
"""

from __future__ import annotations

import argparse
import gzip
import html
import json
import logging
import math
import os
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Any

import pandas as pd
import requests

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)
logger = logging.getLogger("vep_genebe_annotator")

ENSEMBL_REST_URL = "https://rest.ensembl.org/vep/homo_sapiens/region"
GENEBE_API_URL = "https://api.genebe.net/cloud/api-public/v1/variant"
BATCH_SIZE = 150
MAX_RETRIES = 4
BACKOFF_FACTOR = 2.0

CSQ_FORMAT_FIELDS = [
    "Allele",
    "Consequence",
    "IMPACT",
    "SYMBOL",
    "Gene",
    "Feature_type",
    "Feature",
    "BIOTYPE",
    "EXON",
    "INTRON",
    "HGVSc",
    "HGVSp",
    "cDNA_position",
    "CDS_position",
    "Protein_position",
    "Amino_acids",
    "Codons",
    "Existing_variation",
    "DISTANCE",
    "STRAND",
    "FLAGS",
    "VARIANT_CLASS",
    "SYMBOL_SOURCE",
    "HGNC_ID",
    "CANONICAL",
    "CLIN_SIG",
    "SIFT",
    "PolyPhen",
    "GeneBe_ACMG_score",
    "GeneBe_ACMG_classification",
    "GeneBe_ACMG_criteria",
]

CSQ_HEADER_LINE = (
    '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP with GeneBe ACMG Plugin. Format: '
    + "|".join(CSQ_FORMAT_FIELDS)
    + '">'
)

IMPACT_ORDER = {"HIGH": 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1, "": 0}


def open_file(path: str):
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, "r", encoding="utf-8", errors="ignore")


def vcf_to_region_query(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Convert VCF 1-based coordinates to Ensembl REST region string."""
    clean_chr = chrom.replace("chr", "")
    ref = ref.upper()
    alt = alt.upper()

    if len(ref) == 1 and len(alt) == 1:
        return f"{clean_chr} {pos} {pos} {ref}/{alt} 1"
    elif len(ref) < len(alt):
        ins = alt[len(ref):]
        return f"{clean_chr} {pos + 1} {pos} -/{ins} 1"
    elif len(ref) > len(alt):
        del_len = len(ref) - len(alt)
        del_seq = ref[len(alt):]
        return f"{clean_chr} {pos + 1} {pos + del_len} {del_seq}/- 1"
    else:
        return f"{clean_chr} {pos} {pos + len(ref) - 1} {ref}/{alt} 1"


def query_ensembl_batch(batch_regions: list[str]) -> list[dict[str, Any]]:
    """Query Ensembl REST API for batch of variants."""
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    payload = {
        "variants": batch_regions,
        "hgvs": 1,
        "protein": 1,
        "canonical": 1,
        "numbers": 1,
        "domains": 1,
        "variant_class": 1,
    }

    for attempt in range(1, MAX_RETRIES + 1):
        try:
            response = requests.post(ENSEMBL_REST_URL, headers=headers, json=payload, timeout=60)
            if response.status_code == 200:
                return response.json()
            elif response.status_code == 429:
                wait_time = float(response.headers.get("Retry-After", BACKOFF_FACTOR * attempt))
                logger.warning("Ensembl rate limit (429). Backing off for %.1f seconds...", wait_time)
                time.sleep(wait_time)
            elif response.status_code >= 500:
                wait_time = BACKOFF_FACTOR ** attempt
                logger.warning("Ensembl server error (%d). Retrying in %.1f seconds...", response.status_code, wait_time)
                time.sleep(wait_time)
            else:
                break
        except requests.RequestException as e:
            wait_time = BACKOFF_FACTOR ** attempt
            time.sleep(wait_time)

    # Fallback to single queries
    results = []
    for reg in batch_regions:
        try:
            res = requests.post(ENSEMBL_REST_URL, headers=headers, json={"variants": [reg], "hgvs": 1, "canonical": 1}, timeout=20)
            if res.status_code == 200:
                results.extend(res.json())
            time.sleep(0.05)
        except Exception:
            pass
    return results


def query_single_genebe(var_info: tuple[str, int, str, str]) -> tuple[str, dict[str, Any]]:
    """Query GeneBe ACMG classification for a single variant."""
    chrom, pos, ref, alt = var_info
    clean_chr = chrom.replace("chr", "")
    key = f"{chrom}:{pos}:{ref}>{alt}"
    url = f"{GENEBE_API_URL}?chr={clean_chr}&pos={pos}&ref={ref}&alt={alt}&genome=hg38"
    headers = {"User-Agent": "GeneBe_VEP_plugin"}

    for _ in range(3):
        try:
            r = requests.get(url, headers=headers, timeout=15)
            if r.status_code == 200:
                data = r.json()
                variants = data.get("variants", [])
                if variants:
                    v = variants[0]
                    return key, {
                        "score": str(v.get("acmg_score", "")),
                        "classification": v.get("acmg_classification", "Uncertain Significance (VUS)"),
                        "criteria": v.get("acmg_criteria", ""),
                        "gene": v.get("gene_symbol", ""),
                    }
            elif r.status_code == 429:
                time.sleep(1.0)
        except Exception:
            time.sleep(0.5)

    return key, {"score": "", "classification": "Uncertain Significance (VUS)", "criteria": "", "gene": ""}


def query_genebe_parallel(variant_list: list[tuple[str, int, str, str]], max_workers: int = 12) -> dict[str, dict[str, Any]]:
    """Query GeneBe API for all variants concurrently."""
    genebe_map: dict[str, dict[str, Any]] = {}
    logger.info("Querying GeneBe ACMG plugin for %d variants with %d threads...", len(variant_list), max_workers)
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        results = executor.map(query_single_genebe, variant_list)
        for key, res in results:
            genebe_map[key] = res
    logger.info("Received GeneBe ACMG classifications for %d variants.", len(genebe_map))
    return genebe_map


def format_csq_field(val: Any) -> str:
    if val is None or val is False or val == "":
        return ""
    if isinstance(val, list):
        return "&".join(str(x) for x in val if x is not None)
    return str(val).replace("|", "%7C").replace(";", "%3B").replace(",", "%2C")


def parse_transcript_csq(
    alt: str,
    tc: dict[str, Any],
    colocated: list[dict[str, Any]] | None = None,
    genebe_info: dict[str, Any] | None = None,
) -> str:
    consequences = tc.get("consequence_terms", [])
    consequence_str = "&".join(consequences) if consequences else ""
    impact = tc.get("impact", "")
    symbol = tc.get("gene_symbol", "")
    gene_id = tc.get("gene_id", "")
    feature_type = "Transcript" if tc.get("transcript_id") else ""
    feature = tc.get("transcript_id", "")
    biotype = tc.get("biotype", "")
    exon = tc.get("exon", "")
    intron = tc.get("intron", "")
    hgvsc = tc.get("hgvsc", "")
    hgvsp = tc.get("hgvsp", "")
    cdna_pos = tc.get("cdna_start", "")
    cds_pos = tc.get("cds_start", "")
    protein_pos = tc.get("protein_start", "")
    amino_acids = tc.get("amino_acids", "")
    codons = tc.get("codons", "")
    strand = tc.get("strand", "")
    flags = "&".join(tc.get("flags", [])) if tc.get("flags") else ""
    variant_class = tc.get("variant_allele", alt)
    symbol_source = tc.get("gene_symbol_source", "")
    hgnc_id = tc.get("hgnc_id", "")
    canonical = "YES" if tc.get("canonical") == 1 else ""
    sift = tc.get("sift_prediction", "")
    polyphen = tc.get("polyphen_prediction", "")

    # ClinVar / Existing variation
    existing_vars = []
    clin_sigs = []
    if colocated:
        for cv in colocated:
            if cv.get("id"):
                existing_vars.append(cv["id"])
            if cv.get("clin_sig"):
                for cs in cv["clin_sig"]:
                    if cs not in clin_sigs:
                        clin_sigs.append(cs)

    existing_var_str = "&".join(existing_vars)
    clin_sig_str = "&".join(clin_sigs)

    # GeneBe ACMG fields
    gb_score = genebe_info.get("score", "") if genebe_info else ""
    gb_class = genebe_info.get("classification", "") if genebe_info else ""
    gb_crit = genebe_info.get("criteria", "") if genebe_info else ""

    fields = [
        format_csq_field(alt),
        format_csq_field(consequence_str),
        format_csq_field(impact),
        format_csq_field(symbol),
        format_csq_field(gene_id),
        format_csq_field(feature_type),
        format_csq_field(feature),
        format_csq_field(biotype),
        format_csq_field(exon),
        format_csq_field(intron),
        format_csq_field(hgvsc),
        format_csq_field(hgvsp),
        format_csq_field(cdna_pos),
        format_csq_field(cds_pos),
        format_csq_field(protein_pos),
        format_csq_field(amino_acids),
        format_csq_field(codons),
        format_csq_field(existing_var_str),
        "",  # DISTANCE
        format_csq_field(strand),
        format_csq_field(flags),
        format_csq_field(variant_class),
        format_csq_field(symbol_source),
        format_csq_field(hgnc_id),
        format_csq_field(canonical),
        format_csq_field(clin_sig_str),
        format_csq_field(sift),
        format_csq_field(polyphen),
        format_csq_field(gb_score),
        format_csq_field(gb_class),
        format_csq_field(gb_crit),
    ]
    return "|".join(fields)


def annotate_vcfs_with_vep_genebe(
    snp_vcf: str,
    indel_vcf: str,
    output_vcf: str,
    output_tsv: str | None = None,
    output_html: str | None = None,
    sample_name: str | None = None,
) -> None:
    raw_records = []
    genebe_query_list = []
    seen_keys = set()

    for vcf_path in [snp_vcf, indel_vcf]:
        if not vcf_path or not Path(vcf_path).exists():
            continue
        with open_file(vcf_path) as handle:
            for line in handle:
                if line.startswith("#"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 8:
                    continue
                chrom, pos, vid, ref, alts, qual, filt, info = fields[:8]
                pos_int = int(pos)
                for alt in alts.split(","):
                    alt = alt.strip()
                    if not alt:
                        continue
                    vkey = f"{chrom}:{pos_int}:{ref}>{alt}"
                    if vkey in seen_keys:
                        continue
                    seen_keys.add(vkey)
                    region = vcf_to_region_query(chrom, pos_int, ref, alt)
                    raw_records.append(
                        {
                            "key": vkey,
                            "chrom": chrom,
                            "pos": pos_int,
                            "id": vid,
                            "ref": ref,
                            "alt": alt,
                            "qual": qual,
                            "filter": filt,
                            "info": info,
                            "region": region,
                        }
                    )
                    genebe_query_list.append((chrom, pos_int, ref, alt))

    logger.info("Loaded %d variants across input VCFs.", len(raw_records))

    # 1. Ensembl VEP REST Query
    annotation_by_region: dict[str, dict[str, Any]] = {}
    total_batches = math.ceil(len(raw_records) / BATCH_SIZE) if raw_records else 0

    for i in range(0, len(raw_records), BATCH_SIZE):
        batch = raw_records[i : i + BATCH_SIZE]
        regions = [r["region"] for r in batch]
        batch_idx = (i // BATCH_SIZE) + 1
        logger.info("Querying Ensembl VEP REST batch %d/%d (%d variants)...", batch_idx, total_batches, len(regions))
        res_list = query_ensembl_batch(regions)
        for item in res_list:
            inp = item.get("input")
            if inp:
                annotation_by_region[inp] = item
        time.sleep(0.1)

    # 2. GeneBe ACMG Plugin Query
    genebe_results = query_genebe_parallel(genebe_query_list)

    # Write Output VCF
    Path(output_vcf).parent.mkdir(parents=True, exist_ok=True)
    tsv_rows = []

    with open(output_vcf, "w", encoding="utf-8") as out_f:
        out_f.write("##fileformat=VCFv4.2\n")
        out_f.write("##source=Ensembl_VEP_with_GeneBe_Plugin\n")
        out_f.write(f"{CSQ_HEADER_LINE}\n")
        out_f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Most severe gene symbol">\n')
        out_f.write('##INFO=<ID=CONSEQUENCE,Number=1,Type=String,Description="Most severe consequence">\n')
        out_f.write('##INFO=<ID=IMPACT,Number=1,Type=String,Description="Most severe impact">\n')
        out_f.write('##INFO=<ID=ACMG_CLASS,Number=1,Type=String,Description="ACMG/AMP 2015 Classification from GeneBe Plugin">\n')
        out_f.write('##INFO=<ID=ACMG_CRITERIA,Number=.,Type=String,Description="ACMG Criteria from GeneBe Plugin">\n')
        out_f.write('##INFO=<ID=ACMG_SCORE,Number=1,Type=String,Description="ACMG Score from GeneBe Plugin">\n')
        out_f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

        for rec in raw_records:
            annot = annotation_by_region.get(rec["region"])
            gb = genebe_results.get(rec["key"], {})
            csq_entries = []
            most_severe = annot.get("most_severe_consequence", "") if annot else ""
            top_gene = gb.get("gene", "")
            top_impact = ""
            top_clin_sig = ""

            if annot:
                colocated = annot.get("colocated_variants", [])
                tcs = annot.get("transcript_consequences", [])

                def tc_rank(tc):
                    imp = IMPACT_ORDER.get(tc.get("impact", ""), 0)
                    canon = 1 if tc.get("canonical") == 1 else 0
                    return (imp, canon)

                sorted_tcs = sorted(tcs, key=tc_rank, reverse=True)
                for tc in sorted_tcs:
                    csq_entries.append(parse_transcript_csq(rec["alt"], tc, colocated, gb))

                if sorted_tcs:
                    best_tc = sorted_tcs[0]
                    if not top_gene:
                        top_gene = best_tc.get("gene_symbol", "")
                    top_impact = best_tc.get("impact", "")

                if colocated:
                    for cv in colocated:
                        if cv.get("clin_sig"):
                            top_clin_sig = "&".join(cv["clin_sig"])
                            break

            acmg_class = gb.get("classification", "Uncertain Significance (VUS)")
            acmg_crit = gb.get("criteria", "")
            acmg_score = gb.get("score", "")

            info_parts = []
            if rec["info"] and rec["info"] != ".":
                info_parts.append(rec["info"])
            if top_gene:
                info_parts.append(f"GENE={top_gene}")
            if most_severe:
                info_parts.append(f"CONSEQUENCE={most_severe}")
            if top_impact:
                info_parts.append(f"IMPACT={top_impact}")
            if acmg_class:
                info_parts.append(f"ACMG_CLASS={acmg_class}")
            if acmg_crit:
                info_parts.append(f"ACMG_CRITERIA={acmg_crit}")
            if acmg_score != "":
                info_parts.append(f"ACMG_SCORE={acmg_score}")
            if csq_entries:
                info_parts.append(f"CSQ={','.join(csq_entries)}")

            info_str = ";".join(info_parts) if info_parts else "."
            vid_out = rec["id"] if rec["id"] else "."
            out_f.write(
                f"{rec['chrom']}\t{rec['pos']}\t{vid_out}\t{rec['ref']}\t{rec['alt']}\t{rec['qual']}\t{rec['filter']}\t{info_str}\n"
            )

            tsv_rows.append(
                {
                    "CHROM": rec["chrom"],
                    "POS": rec["pos"],
                    "ID": rec["id"],
                    "REF": rec["ref"],
                    "ALT": rec["alt"],
                    "FILTER": rec["filter"],
                    "GENE": top_gene,
                    "CONSEQUENCE": most_severe,
                    "IMPACT": top_impact,
                    "CLINVAR": top_clin_sig,
                    "ACMG_CLASS": acmg_class,
                    "ACMG_SCORE": acmg_score,
                    "ACMG_CRITERIA": acmg_crit,
                }
            )

    logger.info("Successfully wrote VEP+GeneBe annotated VCF to %s", output_vcf)

    df = pd.DataFrame(tsv_rows)
    if output_tsv:
        Path(output_tsv).parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(output_tsv, sep="\t", index=False)
        logger.info("Wrote TSV to %s", output_tsv)

    if output_html:
        Path(output_html).parent.mkdir(parents=True, exist_ok=True)
        write_dashboard_html(df, output_html, sample_name or "Sample")


def write_dashboard_html(df: pd.DataFrame, output_html: str, sample_name: str) -> None:
    total = len(df)
    counts = df["ACMG_CLASS"].value_counts().to_dict() if not df.empty else {}
    p_count = counts.get("Pathogenic", 0)
    lp_count = counts.get("Likely Pathogenic", 0) + counts.get("Likely pathogenic", 0)
    vus_count = counts.get("Uncertain Significance (VUS)", 0) + counts.get("Uncertain significance", 0)
    lb_count = counts.get("Likely Benign", 0) + counts.get("Likely benign", 0)
    b_count = counts.get("Benign", 0)

    rows_html = []
    if not df.empty:
        for _, r in df.head(200).iterrows():
            acmg_val = str(r["ACMG_CLASS"])
            badge_class = "badge-vus"
            if "Pathogenic" in acmg_val and "Likely" not in acmg_val:
                badge_class = "badge-pathogenic"
            elif "Likely Pathogenic" in acmg_val or "Likely pathogenic" in acmg_val:
                badge_class = "badge-likely-pathogenic"
            elif "Likely Benign" in acmg_val or "Likely benign" in acmg_val:
                badge_class = "badge-likely-benign"
            elif "Benign" in acmg_val:
                badge_class = "badge-benign"

            impact_badge = "impact-mod"
            if r["IMPACT"] == "HIGH":
                impact_badge = "impact-high"
            elif r["IMPACT"] == "LOW":
                impact_badge = "impact-low"
            elif r["IMPACT"] == "MODIFIER":
                impact_badge = "impact-modf"

            rows_html.append(
                f"""
            <tr>
                <td><strong>{html.escape(str(r['CHROM']))}:{r['POS']}</strong></td>
                <td><code>{html.escape(str(r['REF']))}&gt;{html.escape(str(r['ALT']))}</code></td>
                <td><span class="gene-tag">{html.escape(str(r['GENE']))}</span></td>
                <td>{html.escape(str(r['CONSEQUENCE']))}</td>
                <td><span class="impact-tag {impact_badge}">{html.escape(str(r['IMPACT']))}</span></td>
                <td><span class="badge {badge_class}">{html.escape(str(r['ACMG_CLASS']))}</span></td>
                <td><strong>{html.escape(str(r['ACMG_SCORE']))}</strong></td>
                <td><small>{html.escape(str(r['ACMG_CRITERIA']))}</small></td>
            </tr>
            """
            )

    table_content = "\n".join(rows_html) if rows_html else "<tr><td colspan='8'>No variants found.</td></tr>"

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <title>VEP + GeneBe ACMG Classification Report - {html.escape(sample_name)}</title>
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
        .impact-low {{ background: #38bdf8; color: #000; }}
        .impact-modf {{ background: #64748b; color: #fff; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>Ensembl VEP + GeneBe ACMG Classification Report</h1>
        <p>Sample: <strong>{html.escape(sample_name)}</strong> | Standards: Ensembl VEP Consequence Prediction + Official GeneBe ACMG 2015 Engine</p>
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
                    <th>ACMG Classification</th>
                    <th>ACMG Score</th>
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
    with open(output_html, "w", encoding="utf-8") as f:
        f.write(html_content)


def main():
    parser = argparse.ArgumentParser(description="Ensembl VEP Annotator with GeneBe ACMG Plugin")
    parser.add_argument("--snp-vcf", required=True, help="Input filtered SNP VCF")
    parser.add_argument("--indel-vcf", required=True, help="Input filtered Indel VCF")
    parser.add_argument("--sample-name", default=None, help="Sample name")
    parser.add_argument("--output-vcf", required=True, help="Output annotated VCF with CSQ & GeneBe ACMG")
    parser.add_argument("--output-tsv", default=None, help="Optional output TSV")
    parser.add_argument("--output-html", default=None, help="Optional output HTML")

    args = parser.parse_args()
    annotate_vcfs_with_vep_genebe(
        snp_vcf=args.snp_vcf,
        indel_vcf=args.indel_vcf,
        output_vcf=args.output_vcf,
        output_tsv=args.output_tsv,
        output_html=args.output_html,
        sample_name=args.sample_name,
    )


if __name__ == "__main__":
    main()
