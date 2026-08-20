#!/usr/bin/env python3
"""
Liftover BED files from GRCh37/hg19 to GRCh38/hg38 using CrossMap and UCSC chain file.
Sorts output BED intervals strictly according to the reference genome FAI order and start coordinates.
"""

import os
import sys
import shutil
import subprocess
import logging

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("liftover_panels")

def get_fai_chrom_order(fai_path):
    """Load chromosome order from reference FASTA index (.fai)."""
    chrom_order = {}
    if not os.path.exists(fai_path):
        logger.warning(f"FAI file not found at {fai_path}, using default ordering")
        return chrom_order
    with open(fai_path, "r") as f:
        for idx, line in enumerate(f):
            parts = line.strip().split("\t")
            if parts:
                chrom_order[parts[0]] = idx
    return chrom_order

def liftover_and_sort(chain_path, in_bed, out_bed, fai_path, unmap_file=None):
    """Perform CrossMap liftover, sort output by FAI chromosome order and coordinate, and write clean BED."""
    logger.info(f"Lifting over {in_bed} -> {out_bed}")
    os.makedirs(os.path.dirname(os.path.abspath(out_bed)), exist_ok=True)
    if unmap_file:
        os.makedirs(os.path.dirname(os.path.abspath(unmap_file)), exist_ok=True)
    else:
        unmap_file = out_bed + ".unmapped.txt"

    tmp_out = out_bed + ".tmp.bed"
    cmd = [
        "CrossMap", "bed",
        "--chromid", "l",
        "--unmap-file", unmap_file,
        chain_path,
        in_bed,
        tmp_out
    ]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        logger.error(f"CrossMap failed for {in_bed}:\n{res.stderr}")
        raise RuntimeError(f"CrossMap failed for {in_bed}")

    chrom_order = get_fai_chrom_order(fai_path)

    # Read and sort records
    records = []
    with open(tmp_out, "r") as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith(("#", "track", "browser")):
                continue
            parts = line_str.split("\t")
            if len(parts) >= 3:
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                rest = parts[3:]
                order_idx = chrom_order.get(chrom, 999999)
                records.append((order_idx, chrom, start, end, rest))

    records.sort(key=lambda x: (x[0], x[2], x[3]))

    with open(out_bed, "w") as f_out:
        for _, chrom, start, end, rest in records:
            if rest:
                f_out.write(f"{chrom}\t{start}\t{end}\t" + "\t".join(rest) + "\n")
            else:
                f_out.write(f"{chrom}\t{start}\t{end}\n")

    if os.path.exists(tmp_out):
        os.remove(tmp_out)

    with open(in_bed) as f:
        in_count = sum(1 for _ in f if _.strip())
    logger.info(f"✓ Completed {os.path.basename(out_bed)}: {in_count} input lines -> {len(records)} lifted records.")
    return len(records)

def main():
    base_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
    mother_dir = os.path.abspath(os.path.join(base_dir, ".."))
    
    chain_file = os.path.join(mother_dir, "hg19ToHg38.over.chain.gz")
    if not os.path.exists(chain_file):
        chain_file = os.path.join(base_dir, "resources/ref/hg19ToHg38.over.chain.gz")
    if not os.path.exists(chain_file):
        logger.error(f"Chain file not found at {chain_file}")
        sys.exit(1)

    ref_fai = os.path.join(base_dir, "resources/ref/grch38/GRCh38.primary_assembly.genome.fa.fai")

    panels = [
        {
            "name": "TargetFile (Overhangs exon +/- 40bp)",
            "in_bed": os.path.join(mother_dir, "ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.mergeBed.bed"),
            "out_grch38": os.path.join(base_dir, "resources/ref/grch38/ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.hg38.mergeBed.bed"),
            "out_grch37": os.path.join(base_dir, "resources/ref/grch37/ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.mergeBed.bed"),
        },
        {
            "name": "CDSFile (Protein Coding Region)",
            "in_bed": os.path.join(mother_dir, "ICC_169Genes_Nextera_V4_ProteinCodingExons.mergeBed.bed"),
            "out_grch38": os.path.join(base_dir, "resources/ref/grch38/ICC_169Genes_Nextera_V4_ProteinCodingExons.hg38.mergeBed.bed"),
            "out_grch37": os.path.join(base_dir, "resources/ref/grch37/ICC_169Genes_Nextera_V4_ProteinCodingExons.mergeBed.bed"),
        },
        {
            "name": "CanonTranFile (Canonical Transcript)",
            "in_bed": os.path.join(mother_dir, "ICC_169Genes_Nextera_V4_ProteinCoding_CanonicalTrans.mergeBed.bed"),
            "out_grch38": os.path.join(base_dir, "resources/ref/grch38/ICC_169Genes_Nextera_V4_ProteinCoding_CanonicalTrans.hg38.mergeBed.bed"),
            "out_grch37": os.path.join(base_dir, "resources/ref/grch37/ICC_169Genes_Nextera_V4_ProteinCoding_CanonicalTrans.mergeBed.bed"),
        },
    ]

    for p in panels:
        if not os.path.exists(p["in_bed"]):
            logger.error(f"Input file not found: {p['in_bed']}")
            sys.exit(1)
        
        # Liftover to GRCh38
        liftover_and_sort(chain_file, p["in_bed"], p["out_grch38"], ref_fai)
        
        # Copy GRCh37 original
        os.makedirs(os.path.dirname(p["out_grch37"]), exist_ok=True)
        shutil.copyfile(p["in_bed"], p["out_grch37"])
        logger.info(f"✓ Copied GRCh37 original to {p['out_grch37']}")

    logger.info("=" * 60)
    logger.info("ALL PANELS SUCCESSFULLY LIFTED OVER AND ORGANIZED")
    logger.info("=" * 60)

if __name__ == "__main__":
    main()
