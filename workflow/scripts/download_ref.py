#!/usr/bin/env python3
"""
Reference genome download, indexing, and validation utility for GRCh38.
Ensures full compatibility with UCSC/GATK chr-prefixed target BED files
(e.g., TruSight_Cardio_TargetedRegions_v1.0.hg38.bed).
"""

import os
import sys
import shutil
import urllib.request
import gzip
import subprocess
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

# Standard reference genome download URLs
REF_URLS = {
    "grch38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    "hg38": "https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz",
    "ensembl_grch38": "https://ftp.ensembl.org/pub/release-110/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"
}

def create_fai_index_python(fasta_path, fai_path):
    """Fallback pure-Python FASTA indexer (.fai generator)."""
    logger.info(f"Generating FASTA index ({fai_path}) via Python indexer...")
    with open(fasta_path, 'rb') as f_in, open(fai_path, 'w') as f_out:
        seq_name = None
        seq_len = 0
        offset = 0
        line_blen = 0
        line_clen = 0
        
        curr_offset = 0
        line = f_in.readline()
        
        while line:
            next_offset = curr_offset + len(line)
            if line.startswith(b'>'):
                if seq_name:
                    f_out.write(f"{seq_name}\t{seq_len}\t{offset}\t{line_clen}\t{line_blen}\n")
                seq_name = line[1:].split()[0].decode('utf-8')
                seq_len = 0
                offset = next_offset
                line_blen = 0
                line_clen = 0
            else:
                line_str = line.rstrip(b'\r\n')
                if line_clen == 0 and len(line_str) > 0:
                    line_clen = len(line_str)
                    line_blen = len(line)
                seq_len += len(line_str)
            curr_offset = next_offset
            line = f_in.readline()
            
        if seq_name:
            f_out.write(f"{seq_name}\t{seq_len}\t{offset}\t{line_clen}\t{line_blen}\n")

def create_dict_index_python(fasta_path, dict_path, fai_path):
    """Generate GATK/Picard sequence dictionary (.dict) file."""
    logger.info(f"Generating Sequence Dictionary ({dict_path})...")
    with open(fai_path, 'r') as f_fai, open(dict_path, 'w') as f_dict:
        f_dict.write("@HD\tVN:1.6\tSO:unsorted\n")
        for line in f_fai:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                chrom, length = parts[0], parts[1]
                f_dict.write(f"@SQ\tSN:{chrom}\tLN:{length}\n")

def index_reference(fasta_path):
    """Index FASTA file using samtools/gatk or Python fallbacks."""
    fai_path = f"{fasta_path}.fai"
    dict_path = os.path.splitext(fasta_path)[0] + ".dict"
    
    # 1. FASTA .fai index
    if not os.path.exists(fai_path):
        if shutil.which("samtools"):
            logger.info("Indexing FASTA with samtools faidx...")
            subprocess.run(["samtools", "faidx", fasta_path])
        if not os.path.exists(fai_path):
            create_fai_index_python(fasta_path, fai_path)
    else:
        logger.info(f"✓ FASTA index exists: {fai_path}")

    # 2. Sequence dictionary .dict
    if not os.path.exists(dict_path):
        if shutil.which("gatk"):
            logger.info("Creating sequence dictionary with GATK...")
            subprocess.run(["gatk", "CreateSequenceDictionary", "-R", fasta_path, "-O", dict_path], stderr=subprocess.DEVNULL)
        elif shutil.which("picard"):
            logger.info("Creating sequence dictionary with Picard...")
            subprocess.run(["picard", "CreateSequenceDictionary", f"R={fasta_path}", f"O={dict_path}"], stderr=subprocess.DEVNULL)
        if not os.path.exists(dict_path):
            create_dict_index_python(fasta_path, dict_path, fai_path)
    else:
        logger.info(f"✓ Sequence dictionary exists: {dict_path}")

    return os.path.exists(fai_path) and os.path.exists(dict_path)

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

def download_reference_genome(target_path="resources/ref/grch38/GRCh38.primary_assembly.genome.fa", genome="grch38", force=False):
    """Download, decompress, index, and validate GRCh38 reference genome."""
    target_path = os.path.abspath(target_path)
    fai_path = f"{target_path}.fai"
    dict_path = os.path.splitext(target_path)[0] + ".dict"
    
    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    gz_path = f"{target_path}.gz" if not target_path.endswith(".gz") else target_path

    if not os.path.exists(target_path) or force:
        url = REF_URLS.get(genome.lower(), REF_URLS["grch38"])
        logger.info(f"Downloading {genome.upper()} reference genome FASTA archive from:\n  {url}")
        logger.info(f"Target location: {target_path}")

        try:
            if shutil.which("wget"):
                subprocess.check_call(["wget", "-c", "-O", gz_path, url])
            elif shutil.which("curl"):
                subprocess.check_call(["curl", "-L", "-C", "-", "-o", gz_path, url])
            else:
                urllib.request.urlretrieve(url, gz_path)

            if gz_path != target_path:
                logger.info("Extracting gzipped FASTA archive and normalizing chromosome names...")
                with gzip.open(gz_path, 'rt') as f_in, open(target_path, 'w') as f_out:
                    for line in f_in:
                        if line.startswith('>'):
                            # Ensure chr prefix for compatibility if needed
                            parts = line[1:].split()
                            chrom = parts[0]
                            if not chrom.startswith('chr'):
                                if chrom == 'MT':
                                    chrom = 'chrM'
                                elif chrom.isalnum() and len(chrom) <= 5:
                                    chrom = f'chr{chrom}'
                                line = f">{chrom} {' '.join(parts[1:])}\n"
                        f_out.write(line)
                if os.path.exists(gz_path):
                    os.remove(gz_path)

            logger.info(f"✓ Downloaded and extracted reference genome to {target_path}")

        except Exception as e:
            logger.error(f"Failed to download reference genome: {e}")
            return False

    indexed = index_reference(target_path)
    return indexed

if __name__ == "__main__":
    target = sys.argv[1] if len(sys.argv) > 1 else "resources/ref/grch38/GRCh38.primary_assembly.genome.fa"
    bed = sys.argv[2] if len(sys.argv) > 2 else "resources/ref/grch38/TruSight_Cardio_TargetedRegions_v1.0.hg38.bed"
    success = download_reference_genome(target)
    if success and os.path.exists(bed):
        validate_bed_compatibility(target, bed)
    sys.exit(0 if success else 1)
