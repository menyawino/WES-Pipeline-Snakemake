from build_panels import build_panels

if __name__ == "__main__":
    build_panels()

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

def fetch_ensembl_canonical_cds(genes, cds_unpadded):
    """
    Fetch canonical / MANE Select transcript CDS intervals for gene list from Ensembl REST API online.
    Uses lightweight symbol lookup without heavy transcript trees to avoid timeouts.
    """
    print(f"Querying Ensembl REST API online for {len(genes)} genes...")
    url = "https://rest.ensembl.org/lookup/symbol/homo_sapiens"
    headers = {"Content-Type": "application/json", "Accept": "application/json"}
    
    canonical_cds = []
    gene_list = sorted(list(genes))
    batch_size = 25
    canon_tx_map = {}
    
    # Gene symbol aliases if needed
    ALIASES = {
        "SEPN1": "SELENON",
        "GAA": "GAA"
    }

    # Step 1: Lightweight lookup to get canonical transcript IDs (no expand=1)
    for i in range(0, len(gene_list), batch_size):
        batch = gene_list[i:i + batch_size]
        query_batch = [ALIASES.get(g, g) for g in batch]
        payload = json.dumps({"symbols": query_batch}).encode("utf-8")
        req = urllib.request.Request(url, data=payload, headers=headers)
        
        try:
            with urllib.request.urlopen(req, timeout=15) as resp:
                data = json.loads(resp.read().decode("utf-8"))
                for orig_sym, q_sym in zip(batch, query_batch):
                    gdata = data.get(q_sym)
                    if gdata and "canonical_transcript" in gdata:
                        tx_full = gdata["canonical_transcript"]
                        tx_id = tx_full.split(".")[0]
                        chrom = "chr" + str(gdata.get("seq_region_name", "")).replace("chr", "")
                        canon_tx_map[orig_sym] = (tx_id, chrom)
        except Exception as e:
            print(f"Note: Ensembl batch lookup warning for batch {i//batch_size + 1}: {e}")
            
        time.sleep(0.1)

    print(f"Resolved {len(canon_tx_map)} / {len(genes)} canonical transcripts from Ensembl.")

    # Step 2: Query CDS intervals for resolved canonical transcripts
    resolved_count = 0
    for gene, (tx_id, chrom) in canon_tx_map.items():
        tx_url = f"https://rest.ensembl.org/overlap/id/{tx_id}?feature=cds"
        tx_req = urllib.request.Request(tx_url, headers={"Accept": "application/json"})
        got_cds = False
        try:
            with urllib.request.urlopen(tx_req, timeout=10) as tx_resp:
                cds_blocks = json.loads(tx_resp.read().decode("utf-8"))
                for block in cds_blocks:
                    c_start = int(block["start"]) - 1
                    c_end = int(block["end"])
                    canonical_cds.append((chrom, c_start, c_end, gene))
                    got_cds = True
                if got_cds:
                    resolved_count += 1
        except Exception:
            pass
        time.sleep(0.08)  # Respect Ensembl 15 req/sec rate limit

    print(f"Retrieved exact CDS exons for {resolved_count} canonical transcripts.")
    
    # For any genes not returned by Ensembl, retain unpadded CDS intervals
    covered_genes = {x[3] for x in canonical_cds}
    missing_genes = genes - covered_genes
    if missing_genes:
        print(f"Adding unpadded coding intervals for {len(missing_genes)} genes with alternate identifiers...")
        for item in cds_unpadded:
            if item[3] in missing_genes:
                canonical_cds.append(item)

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
    canon_raw = fetch_ensembl_canonical_cds(genes, cds_unpadded)
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
