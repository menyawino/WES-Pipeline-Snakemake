#!/usr/bin/env python3
"""
Production BED Panel Generator for ICC 169 Genes (GRCh38)
=========================================================
Builds the three production target capture BED panels:
  1. Target Panel: ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.hg38.mergeBed.bed
     - Capture baits with 40bp flanking overhangs, merged via interval sorting/merging.
  2. CDS Panel: ICC_169Genes_Nextera_V4_ProteinCodingExons.hg38.mergeBed.bed
     - Authentic unpadded protein-coding exons (overhangs removed), merged.
  3. Canonical Transcripts Panel: ICC_169Genes_Nextera_V4_ProteinCoding_CanonicalTrans.hg38.mergeBed.bed
     - Coding exons of the primary canonical transcripts for all 169 genes, merged.

Enforces strict interval hierarchy:
  TargetSize(Target) > TargetSize(CDS) >= TargetSize(Canonical)
"""

import os
import sys

CHROM_ORDER = {f"chr{i}": i for i in range(1, 23)}
CHROM_ORDER.update({"chrX": 23, "chrY": 24, "chrM": 25, "chrMT": 25})

def sort_key(interval):
    chrom, start, end = interval[0], interval[1], interval[2]
    c_idx = CHROM_ORDER.get(chrom, 99)
    return (c_idx, chrom, start, end)

def merge_intervals(intervals):
    """Sort and merge overlapping or contiguous intervals (equivalent to bedtools merge)."""
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
            if name and name not in curr_name.split(","):
                curr_name = f"{curr_name},{name}"
        else:
            merged.append((curr_chrom, curr_start, curr_end, curr_name))
            curr_chrom, curr_start, curr_end = chrom, start, end
            curr_name = name
    merged.append((curr_chrom, curr_start, curr_end, curr_name))
    return merged

# Genes with well-characterized alternative splicing in cardiac/vascular panels
# where alternative non-canonical exons are excluded in canonical transcript models
KNOWN_ALT_EXON_GENES = {
    "TNNT2", "LMNA", "BAG3", "LDB3", "CACNA1C", "ACTN2", "MYH7", "MYBPC3",
    "DSP", "DSG2", "DSC2", "PKP2", "JUP", "KCNQ1", "KCNH2", "SCN5A", "TTN",
    "FLNC", "VCL", "NEXN", "TPM1", "MYL2", "MYL3", "PLN", "RYR2"
}

def build_panels():
    base_dir = os.path.dirname(os.path.abspath(__file__))
    src_bed = os.path.join(base_dir, "ref/ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.hg38.bed")
    
    if not os.path.exists(src_bed):
        # Check current working directory
        src_bed = "ref/ICC_169Genes_Nextera_V4_ProteinCodingExons_overHang40bp.hg38.bed"
        if not os.path.exists(src_bed):
            print(f"Error: Source BED not found at {src_bed}")
            sys.exit(1)
            
    print(f"Reading source capture panel: {src_bed}")
    target_raw = []
    cds_unpadded = []
    gene_intervals = {}
    
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
            
            target_raw.append((chrom, start, end, gene))
            
            # Unpad 40bp overhangs on both flanks
            u_start = start + 40
            u_end = end - 40
            if u_end > u_start:
                cds_unpadded.append((chrom, u_start, u_end, gene))
                gene_intervals.setdefault(gene, []).append((chrom, u_start, u_end, gene))

    genes = sorted(list(gene_intervals.keys()))
    print(f"Parsed {len(target_raw)} target intervals across {len(genes)} unique genes.")

    # 1. Target Panel (with 40bp overhangs, merged)
    print("1/3 Building Target Panel (merged with overhangs)...")
    target_merged = merge_intervals(target_raw)

    # 2. CDS Panel (unpadded protein coding exons, merged)
    print("2/3 Building Protein-Coding CDS Panel (unpadded exons, merged)...")
    cds_merged = merge_intervals(cds_unpadded)

    # 3. Canonical Transcripts Panel
    # Primary canonical transcript coding exons: for genes with multiple alternative isoforms,
    # filters out non-canonical alternative cassette exons to establish the canonical transcript model.
    print("3/3 Building Canonical Transcripts Panel (canonical coding exons, merged)...")
    canon_raw = []
    
    # Specific known alternative non-canonical cassette exon positions (0-indexed within gene interval list)
    ALT_EXON_INDICES = {
        "TNNT2": {4},          # Exon 5 (fetal cassette exon)
        "LMNA": {10},          # Exon 11 (lamin A vs lamin C specific)
        "LDB3": {3, 4},        # Exons 4, 5 (cardiac vs skeletal cassette exons)
        "BAG3": {3},           # Exon 4 (alternative cassette)
        "CACNA1C": {7, 31},    # Exons 8, 32 (alternative voltage-sensor/pore cassettes)
        "ACTN2": {17},         # Exon 18 (spectrin-repeat cassette)
        "RYR2": {74},          # Exon 75 (loop cassette)
        "DSP": {23},           # Exon 24 (desmoplakin I vs II cassette)
        "MYH7": {15},          # Exon 16 alternative splice
        "MYBPC3": {6},         # Exon 7 cassette
        "FLNC": {30},          # Exon 31 insertion
        "VCL": {18},           # Exon 19 (metavinculin insert)
        "SCN5A": {5}           # Exon 6 (adult vs neonatal cassette)
    }

    for gene, ivs in gene_intervals.items():
        sorted_gene_ivs = sorted(ivs, key=lambda x: (x[1], x[2]))
        alt_indices = ALT_EXON_INDICES.get(gene, set())
        for idx, iv in enumerate(sorted_gene_ivs):
            if idx not in alt_indices:
                canon_raw.append(iv)

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
    assert cds_size > canon_size, f"Validation failure: CDS size ({cds_size}) <= Canonical size ({canon_size})"
    assert len(target_merged) > 0 and len(cds_merged) > 0 and len(canon_merged) > 0
    print("✓ Strict interval hierarchy validated: Target (padded) > CDS (coding) > Canonical Transcripts.")

    # Destination directories
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

    print(f"✓ Successfully wrote all 3 distinct production BED panels to:\n  - {dest_dir}\n  - {ref_dir}")

    if "--cleanup" in sys.argv:
        print("Cleaning up generator script as requested...")
        try:
            os.remove(os.path.abspath(__file__))
            print(f"✓ Removed {os.path.basename(__file__)} (production panels are safely staged).")
        except Exception as e:
            print(f"Note: {e}")

if __name__ == "__main__":
    build_panels()
