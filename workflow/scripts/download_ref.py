#!/usr/bin/env python3
"""
Reference genome download, indexing, and validation utility for GRCh38/hg38.
Downloads official Broad Institute GATK Best Practices reference bundle from fast Google Cloud mirrors.
"""

import os
import sys
import shutil
import urllib.request
import subprocess
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

GCS_BROAD_HG38_BASE = "https://storage.googleapis.com/gcp-public-data--broad-references/hg38/v0"

REF_URLS = {
    "grch38": f"{GCS_BROAD_HG38_BASE}/Homo_sapiens_assembly38.fasta",
    "hg38": f"{GCS_BROAD_HG38_BASE}/Homo_sapiens_assembly38.fasta",
    "fai": f"{GCS_BROAD_HG38_BASE}/Homo_sapiens_assembly38.fasta.fai",
    "dict": f"{GCS_BROAD_HG38_BASE}/Homo_sapiens_assembly38.dict",
    "dbsnp": f"{GCS_BROAD_HG38_BASE}/Homo_sapiens_assembly38.dbsnp138.vcf",
    "dbsnp_idx": f"{GCS_BROAD_HG38_BASE}/Homo_sapiens_assembly38.dbsnp138.vcf.idx",
}

def fast_download(url, target_file):
    """Download a file with curl / wget or urllib."""
    os.makedirs(os.path.dirname(os.path.abspath(target_file)), exist_ok=True)
    if os.path.exists(target_file) and os.path.getsize(target_file) > 0:
        logger.info(f"✓ Already present: {target_file}")
        return True
    logger.info(f"Downloading from {url} -> {target_file}")
    try:
        if shutil.which("curl"):
            cmd = ["curl", "-f", "-L", "-C", "-", "-o", target_file, url]
            res = subprocess.run(cmd)
            if res.returncode == 0 and os.path.getsize(target_file) > 0:
                return True
        if shutil.which("wget"):
            cmd = ["wget", "-c", "-O", target_file, url]
            res = subprocess.run(cmd)
            if res.returncode == 0 and os.path.getsize(target_file) > 0:
                return True
        urllib.request.urlretrieve(url, target_file)
        return os.path.exists(target_file) and os.path.getsize(target_file) > 0
    except Exception as e:
        logger.error(f"Download failed for {url}: {e}")
        return False

def validate_bed_compatibility(fasta_path, bed_path):
    """Validate chromosome compatibility between reference FASTA index and target BED."""
    logger.info(f"Validating target BED compatibility: {bed_path}")
    if not os.path.exists(bed_path):
        logger.warning(f"Target BED file does not exist at: {bed_path}")
        return False
        
    fai_path = f"{fasta_path}.fai"
    if not os.path.exists(fai_path):
        logger.error(f"Cannot validate BED: Missing FASTA index {fai_path}")
        return False

    with open(fai_path, 'r') as f:
        fasta_chroms = set(line.split('\t')[0] for line in f if line.strip())

    bed_chroms = set()
    with open(bed_path, 'r') as f:
        for line in f:
            if line.strip() and not line.startswith(('#', 'track', 'browser')):
                parts = line.strip().split('\t')
                if parts:
                    bed_chroms.add(parts[0])

    missing_chroms = bed_chroms - fasta_chroms
    if missing_chroms:
        logger.error(f"✗ Chromosome naming mismatch! BED contains {len(missing_chroms)} chromosomes not found in reference FASTA: {sorted(list(missing_chroms))[:5]}")
        return False

    logger.info(f"✓ BED compatibility verified! All {len(bed_chroms)} chromosomes in BED match the reference FASTA ({fasta_path}).")
    return True

def index_reference(fasta_path):
    """Ensure .fai, .dict, and bwa-mem2 indexes exist."""
    fai_path = f"{fasta_path}.fai"
    dict_path = os.path.splitext(fasta_path)[0] + ".dict"
    
    # 1. Download or generate .fai
    if not os.path.exists(fai_path) or os.path.getsize(fai_path) == 0:
        if not fast_download(REF_URLS["fai"], fai_path):
            samtools_bin = "/home/omar/Downloads/miniconda3/envs/icc_04_alignment/bin/samtools"
            if not os.path.exists(samtools_bin):
                samtools_bin = shutil.which("samtools")
            if samtools_bin:
                logger.info("Indexing FASTA with samtools faidx...")
                subprocess.run([samtools_bin, "faidx", fasta_path])

    # 2. Download or generate .dict
    if not os.path.exists(dict_path) or os.path.getsize(dict_path) == 0:
        fast_download(REF_URLS["dict"], dict_path)

    # 3. Check / build bwa-mem2 index
    bwa_bwt = f"{fasta_path}.bwt.2bit.64"
    if not os.path.exists(bwa_bwt):
        bwa_bin = "/home/omar/Downloads/miniconda3/envs/icc_04_alignment/bin/bwa-mem2"
        if not os.path.exists(bwa_bin):
            bwa_bin = shutil.which("bwa-mem2")
        if bwa_bin:
            logger.info("Generating bwa-mem2 index (this runs once)...")
            subprocess.run([bwa_bin, "index", fasta_path])

    return os.path.exists(fai_path)

def download_reference_genome(target_path="resources/ref/grch38/GRCh38.primary_assembly.genome.fa", genome="grch38", force=False):
    """Download GRCh38 reference genome from fast Broad GCS mirror."""
    target_path = os.path.abspath(target_path)
    os.makedirs(os.path.dirname(target_path), exist_ok=True)

    if not os.path.exists(target_path) or os.path.getsize(target_path) == 0 or force:
        url = REF_URLS.get(genome.lower(), REF_URLS["grch38"])
        logger.info(f"Downloading reference genome from GCS Broad mirror:\n  {url}")
        success = fast_download(url, target_path)
        if not success:
            logger.error("Failed to download reference genome")
            return False
            
    # Also download dbSNP if configured in resources
    dbsnp_target = "resources/ref/grch38/dbsnp/Homo_sapiens_assembly38.dbsnp138.vcf"
    fast_download(REF_URLS["dbsnp"], dbsnp_target)
    fast_download(REF_URLS["dbsnp_idx"], f"{dbsnp_target}.idx")

    return index_reference(target_path)

if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "resources/ref/grch38/GRCh38.primary_assembly.genome.fa"
    force = "--force" in sys.argv
    success = download_reference_genome(target, force=force)
    sys.exit(0 if success else 1)
