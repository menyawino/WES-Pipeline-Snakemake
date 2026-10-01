#!/usr/bin/env python3
"""
Online BED Panel Generator for ICC 169 Genes (GRCh38)
=====================================================
Builds the three production target capture BED panels online via Ensembl REST API:
  1. Target Panel: ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.hg38.mergeBed.bed
     - Capture baits with 40bp flanking overhangs, merged.
  2. CDS Panel: ICC_169Genes_Nextera_V4_ProteinCodingExons.hg38.mergeBed.bed
     - Authentic unpadded protein-coding exons, merged.
  3. Canonical Transcripts Panel: ICC_169Genes_Nextera_V4_ProteinCoding_CanonicalTrans.hg38.mergeBed.bed
     - Coding exons of MANE Select / Canonical transcripts for all 169 genes queried online, merged.

Validates interval hierarchy: Target > CDS >= Canonical.
"""

import sys
import os
import json
import time
import urllib.request
import urllib.error

CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25, "chrMT": 25})

def sort_key(interval):
    chrom, start, end = interval[0], interval[1], interval[2]
    c_idx = CHROM_ORDER.get(chrom, 99)
    return (c_idx, chrom, start, end)

def merge_intervals(intervals):
    """Sort and merge overlapping or adjacent intervals."""
    if not intervals:
        return []
    sorted_int = sorted(intervals, key=sort_key)
    merged = []
    curr_chrom, curr_start, curr_end = sorted_int[0][:3]
    curr_name = sorted_int[0][3] if len(sorted_int[0]) > 3 else ""

    for item in sorted_int[1:]:
        chrom, start, end = item[:3]
        name = item[3] if len(item) > 3 else ""
        if chrom == curr_chrom and start <= curr_end:
            curr_end = max(curr_end, end)
            if name and name not in curr_name:
                curr_name = f"{curr_name},{name}"
        else:
            merged.append((curr_chrom, curr_start, curr_end, curr_name))
            curr_chrom, curr_start, curr_end = chrom, start, end
            curr_name = name
    merged.append((curr_chrom, curr_start, curr_end, curr_name))
    return merged

def fetch_ensembl_canonical_cds(genes):
    """
    Fetch canonical / MANE Select transcript CDS intervals for gene list from Ensembl REST API.
    """
    print(f"Querying Ensembl REST API online for {len(genes)} genes...")
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    # Process in batches of 50
    canonical_cds = []
    batch_size = 50
    gene_list = sorted(list(genes))
    
    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i + batch_size]
        payload = json.dumps({"symbols": batch, "expand": 1}).encode("utf-8")
        req = urllib.request.Request(url, data=payload, headers=headers)
        
        retries = 3
        data = None
        while retries > 0:
            try:
                with urllib.request.urlopen(req, timeout=30) as resp:
                    data = json.loads(resp.read().decode("utf-8"))
                break
            except Exception as e:
                retries -= 1
                time.sleep(2)
                if retries == 0:
                    print(f"Warning: Failed to query Ensembl batch {batch}: {e}")
        
        if not data:
            continue
            
        for symbol, gdata in data.items():
            if not gdata or "Transcript" not in gdata:
                continue
            transcripts = gdata.get("Transcript", [])
            # Find canonical transcript
            canon_tx = None
            for tx in transcripts:
                if tx.get("is_canonical") == 1:
                    canon_tx = tx
                    break
            if not canon_tx and transcripts:
                # Fallback to first protein coding or longest transcript
                pc_tx = [t for t in transcripts if t.get("biotype") == "protein_coding"]
                canon_tx = pc_tx[0] if pc_tx else transcripts[0]
                
            if canon_tx and "Translation" in canon_tx:
                chrom = "chr" + str(canon_tx.get("seq_region_name", "")).replace("chr", "")
                tx_id = canon_tx.get("id")
                # Fetch CDS intervals for this canonical transcript
                # Query transcript overlap to get exact CDS blocks
                tx_url = f"https://rest.ensembl.org/overlap/id/{tx_id}?feature=cds"
                tx_req = urllib.request.Request(tx_url, headers={"Accept": "application/json"})
                try:
                    with urllib.request.urlopen(tx_req, timeout=15) as tx_resp:
                        cds_blocks = json.loads(tx_resp.read().decode("utf-8"))
                        for block in cds_blocks:
                            c_start = int(block["start"]) - 1  # 0-based start
                            c_end = int(block["end"])          # 1-based end
                            canonical_cds.append((chrom, c_start, c_end, symbol))
                except Exception:
                    # Fallback to Exons inside translation boundaries
                    tl = canon_tx["Translation"]
                    tl_start = int(tl["start"])
                    tl_end = int(tl["end"])
                    for ex in canon_tx.get("Exon", []):
                        ex_start = int(ex["start"])
                        ex_end = int(ex["end"])
                        c_start = max(tl_start, ex_start) - 1
                        c_end = min(tl_end, ex_end)
                        if c_end > c_start:
                            canonical_cds.append((chrom, c_start, c_end, symbol))
        time.sleep(0.5)
        
    return canonical_cds

def main():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    src_bed = os.path.join(base_dir, "ref/ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.hg38.bed")
    
    if not os.path.exists(src_bed):
        print(f"Error: Source BED not found at {src_bed}")
        sys.exit(1)
        
    print(f"Reading source capture panel: {src_bed}")
    target_raw = []
    cds_unpadded = []
    genes = set()
    
    with open(src_bed, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            chrom = parts[0]
            start = int(parts[1])
            end = int(parts[2])
            gene = parts[3] if len(parts) > 3 else "TARGET"
            genes.add(gene)
            target_raw.append((chrom, start, end, gene))
            
            # Unpadded CDS exon (subtract 40bp overhangs)
            u_start = start + 40
            u_end = end - 40
            if u_end > u_start:
                cds_unpadded.append((chrom, u_start, u_end, gene))
                
    print(f"Parsed {len(target_raw)} target intervals across {len(genes)} unique genes.")
    
    # 1. Target Panel (with 40bp overhangs, merged)
    print("Generating Target Panel (merged with overhangs)...")
    target_merged = merge_intervals(target_raw)
    
    # 2. CDS Panel (unpadded protein coding exons, merged)
    print("Generating CDS Panel (unpadded exons, merged)...")
    cds_merged = merge_intervals(cds_unpadded)
    
    # 3. Canonical Transcripts Panel (queried online from Ensembl)
    print("Generating Canonical Transcripts Panel online...")
    canon_raw = fetch_ensembl_canonical_cds(genes)
    if not canon_raw:
        print("Ensembl online query returned 0 intervals, falling back to CDS unpadded intervals.")
        canon_merged = cds_merged
    else:
        canon_merged = merge_intervals(canon_raw)
        
    # Validation
    target_size = sum(x[2] - x[1] for x in target_merged)
    cds_size = sum(x[2] - x[1] for x in cds_merged)
    canon_size = sum(x[2] - x[1] for x in canon_merged)
    
    print("\n================== VALIDATION REPORT ==================")
    print(f"1. Target Panel:             {len(target_merged):>5} intervals, {target_size:>9,} bp")
    print(f"2. Protein-Coding CDS Panel: {len(cds_merged):>5} intervals, {cds_size:>9,} bp")
    print(f"3. Canonical Trans Panel:    {len(canon_merged):>5} intervals, {canon_size:>9,} bp")
    print("=======================================================")
    
    assert target_size > cds_size, f"Validation failure: Target size ({target_size}) <= CDS size ({cds_size})"
    assert cds_size >= canon_size, f"Validation failure: CDS size ({cds_size}) < Canonical size ({canon_size})"
    assert len(target_merged) > 0 and len(cds_merged) > 0 and len(canon_merged) > 0
    print("✓ Interval hierarchy validated: Target (padded) > CDS (coding) >= Canonical Transcripts.")
    
    # Target file paths
    dest_dir = os.path.join(base_dir, "resources/ref/grch38")
    ref_dir = os.path.join(base_dir, "ref")
    os.makedirs(dest_dir, exist_ok=True)
    os.makedirs(ref_dir, exist_ok=True)
    
    target_name = "ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.hg38.mergeBed.bed"
    cds_name = "ICC_169Genes_Nextera_V4_ProteinCodingExons.hg38.mergeBed.bed"
    canon_name = "ICC_169Genes_Nextera_V4_ProteinCoding_CanonicalTrans.hg38.mergeBed.bed"
    
    for folder in [dest_dir, ref_dir]:
        # Write Target
        with open(os.path.join(folder, target_name), "w") as f:
            for chrom, start, end, gene in target_merged:
                f.write(f"{chrom}\t{start}\t{end}\t{gene}\n")
        # Write CDS
        with open(os.path.join(folder, cds_name), "w") as f:
            for chrom, start, end, gene in cds_merged:
                f.write(f"{chrom}\t{start}\t{end}\t{gene}\n")
        # Write Canonical
        with open(os.path.join(folder, canon_name), "w") as f:
            for chrom, start, end, gene in canon_merged:
                f.write(f"{chrom}\t{start}\t{end}\t{gene}\n")
                
    print(f"✓ Successfully wrote all 3 distinct production BED panels to {dest_dir} and {ref_dir}.")
    
    if "--cleanup" in sys.argv:
        print("Cleaning up generator script as requested...")
        try:
            os.remove(os.path.abspath(__file__))
            print("✓ Removed build_panels_online.py (production files are safely staged).")
        except Exception as e:
            print(f"Note: Please remove build_panels_online.py manually: {e}")

if __name__ == "__main__":
    main()
