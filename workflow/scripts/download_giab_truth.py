#!/usr/bin/env python3
"""
GIAB Truth Data Downloader for WES Pipeline.
Downloads and verifies NIST / GIAB high-confidence benchmark truth sets (HG001 / HG002) for GRCh38.
"""

import os
import sys
import argparse
import urllib.request
import urllib.error
import subprocess
import time

GIAB_URLS = {
    "HG001": {
        "vcf": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        "tbi": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
        "bed": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/NA12878_HG001/latest/GRCh38/HG001_GRCh38_1_22_v4.2.1_benchmark.bed",
    },
    "HG002": {
        "vcf": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz",
        "tbi": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi",
        "bed": "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/NISTv4.2.1/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.bed",
    }
}


def download_file(url, destination, force=False):
    """Downloads a file with progress tracking."""
    if os.path.exists(destination) and not force:
        if os.path.getsize(destination) > 0:
            print(f"[INFO] File already exists and is not empty: {destination}")
            return True

    os.makedirs(os.path.dirname(os.path.abspath(destination)), exist_ok=True)
    temp_dest = f"{destination}.tmp.{int(time.time())}"
    print(f"[INFO] Downloading {url} -> {destination} ...")

    # Try curl / wget if available for speed and reliability, otherwise fallback to urllib
    try:
        cmd = ["curl", "-L", "-f", "-o", temp_dest, "--retry", "3", "--retry-delay", "2", url]
        res = subprocess.run(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        if res.returncode == 0 and os.path.exists(temp_dest) and os.path.getsize(temp_dest) > 0:
            os.replace(temp_dest, destination)
            print(f"[SUCCESS] Downloaded: {destination} ({os.path.getsize(destination):,} bytes)")
            return True
    except Exception:
        pass

    try:
        req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
        with urllib.request.urlopen(req, timeout=120) as response, open(temp_dest, 'wb') as out_file:
            total_size = int(response.info().get('Content-Length', -1))
            downloaded = 0
            block_size = 1024 * 1024  # 1 MB chunks
            while True:
                buffer = response.read(block_size)
                if not buffer:
                    break
                downloaded += len(buffer)
                out_file.write(buffer)
                if total_size > 0:
                    pct = (downloaded / total_size) * 100
                    sys.stdout.write(f"\r  Downloaded {downloaded / (1024*1024):.1f} MB / {total_size / (1024*1024):.1f} MB ({pct:.1f}%)")
                    sys.stdout.flush()
            print()
        if os.path.exists(temp_dest) and os.path.getsize(temp_dest) > 0:
            os.replace(temp_dest, destination)
            print(f"[SUCCESS] Downloaded: {destination} ({os.path.getsize(destination):,} bytes)")
            return True
    except Exception as e:
        if os.path.exists(temp_dest):
            os.remove(temp_dest)
        print(f"[ERROR] Failed to download {url}: {e}", file=sys.stderr)
        return False
    return False


def download_giab_truth_set(sample="HG001", outdir="resources/giab/HG001", force=False):
    """Downloads all truth resources for a given GIAB sample."""
    sample = sample.upper()
    if sample in ["NA12878", "HG001"]:
        sample_key = "HG001"
    elif sample in ["NA24385", "HG002"]:
        sample_key = "HG002"
    else:
        print(f"[ERROR] Unsupported GIAB sample '{sample}'. Supported samples: HG001 (NA12878), HG002 (NA24385)", file=sys.stderr)
        return False

    targets = GIAB_URLS[sample_key]
    os.makedirs(outdir, exist_ok=True)

    vcf_filename = os.path.basename(targets["vcf"])
    tbi_filename = os.path.basename(targets["tbi"])
    bed_filename = os.path.basename(targets["bed"])

    vcf_dest = os.path.join(outdir, vcf_filename)
    tbi_dest = os.path.join(outdir, tbi_filename)
    bed_dest = os.path.join(outdir, bed_filename)

    print(f"\n=======================================================")
    print(f" Downloading GIAB Truth Set for {sample_key} (GRCh38)")
    print(f" Target Directory: {outdir}")
    print(f"=======================================================\n")

    ok_vcf = download_file(targets["vcf"], vcf_dest, force=force)
    ok_tbi = download_file(targets["tbi"], tbi_dest, force=force)
    ok_bed = download_file(targets["bed"], bed_dest, force=force)

    # If tabix index failed to download, build it locally if tabix is available
    if ok_vcf and not ok_tbi:
        print("[INFO] Attempting to index truth VCF with tabix locally...")
        try:
            subprocess.run(["tabix", "-f", "-p", "vcf", vcf_dest], check=True)
            ok_tbi = True
            print(f"[SUCCESS] Indexed truth VCF with tabix: {tbi_dest}")
        except Exception as e:
            print(f"[WARNING] Could not create tabix index: {e}")

    if ok_vcf and ok_bed:
        print(f"\n[SUCCESS] All GIAB {sample_key} truth files are ready in {outdir}\n")
        return {
            "vcf": vcf_dest,
            "tbi": tbi_dest if ok_tbi else None,
            "bed": bed_dest
        }
    else:
        print(f"\n[ERROR] One or more truth files failed to download.", file=sys.stderr)
        return False


def main():
    parser = argparse.ArgumentParser(description="Download GIAB high-confidence benchmark truth sets.")
    parser.add_argument("--sample", default="HG001", choices=["HG001", "NA12878", "HG002", "NA24385"], help="GIAB sample name (default: HG001)")
    parser.add_argument("--outdir", default="resources/giab/HG001", help="Target output directory")
    parser.add_argument("--force", action="store_true", help="Force re-download even if files already exist")
    args = parser.parse_args()

    result = download_giab_truth_set(sample=args.sample, outdir=args.outdir, force=args.force)
    if not result:
        sys.exit(1)


if __name__ == "__main__":
    main()
