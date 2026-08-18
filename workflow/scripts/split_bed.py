#!/usr/bin/env python3
"""
Split a target BED file into N balanced, non-overlapping BED files based on genomic interval sizes.
Applies interval padding and merges overlapping intervals prior to chunking, ensuring strict
reference coordinate ordering and non-overlapping chunks for scatter-gather GATK HaplotypeCaller.
"""

import sys
import os

def get_contig_info(fai_path=None):
    """Load contig ordering map and lengths from reference FAI if available."""
    contig_order = {}
    contig_lengths = {}
    if fai_path and os.path.exists(fai_path):
        with open(fai_path, 'r') as f:
            for i, line in enumerate(f):
                parts = line.strip().split('\t')
                if parts:
                    contig_order[parts[0]] = i
                    if len(parts) > 1:
                        try:
                            contig_lengths[parts[0]] = int(parts[1])
                        except ValueError:
                            pass
        return contig_order, contig_lengths

    # Standard fallback
    standard = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM", "chrMT"]
    standard += [str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"]
    contig_order = {c: i for i, c in enumerate(standard)}
    return contig_order, contig_lengths

def split_bed(bed_path, output_paths, fai_path=None, padding=100):
    """
    Split input BED file into balanced sub-BED files sorted by reference order.
    Pads intervals and merges overlaps to prevent out-of-order coordinate collisions during GatherVcfs.
    """
    contig_order, contig_lengths = get_contig_info(fai_path)
    raw_intervals = []
    
    with open(bed_path, 'r') as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith(('#', 'track', 'browser')):
                continue
            parts = line_str.split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                try:
                    start = int(parts[1])
                    end = int(parts[2])
                except ValueError:
                    continue
                
                # Apply padding
                padded_start = max(0, start - padding)
                if chrom in contig_lengths:
                    padded_end = min(contig_lengths[chrom], end + padding)
                else:
                    padded_end = end + padding
                
                raw_intervals.append((chrom, padded_start, padded_end))
                
    if not raw_intervals:
        for p in output_paths:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
            open(p, 'w').close()
        return

    # Sort strictly by reference contig order, then by start position
    raw_intervals.sort(key=lambda item: (contig_order.get(item[0], 999999), item[1], item[2]))

    # Merge overlapping or abutting intervals
    merged_intervals = []
    for chrom, start, end in raw_intervals:
        if not merged_intervals:
            merged_intervals.append([chrom, start, end])
        else:
            last_chrom, last_start, last_end = merged_intervals[-1]
            if chrom == last_chrom and start <= last_end:
                merged_intervals[-1][2] = max(last_end, end)
            else:
                merged_intervals.append([chrom, start, end])

    # Calculate intervals with lengths
    intervals = [(c, s, e, max(1, e - s)) for c, s, e in merged_intervals]
    total_bases = sum(item[3] for item in intervals)
    n_chunks = len(output_paths)
    target_bases_per_chunk = max(1, total_bases // n_chunks)
    
    chunks = [[] for _ in range(n_chunks)]
    current_chunk_idx = 0
    current_chunk_bases = 0
    
    for item in intervals:
        chunks[current_chunk_idx].append(item)
        current_chunk_bases += item[3]
        if (current_chunk_bases >= target_bases_per_chunk) and (current_chunk_idx < n_chunks - 1):
            current_chunk_idx += 1
            current_chunk_bases = 0
            
    for out_path, chunk_items in zip(output_paths, chunks):
        os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
        with open(out_path, 'w') as f_out:
            for chrom, start, end, _ in chunk_items:
                f_out.write(f"{chrom}\t{start}\t{end}\n")

if __name__ == '__main__':
    # Support direct execution via snakemake script or CLI
    if 'snakemake' in globals():
        in_bed = snakemake.input.get('target', snakemake.input[0])
        ref_fai = snakemake.input.get('ref_fai', None)
        out_beds = snakemake.output.chunks
        split_bed(in_bed, out_beds, ref_fai)
    elif len(sys.argv) > 2:
        in_bed = sys.argv[1]
        # Check if second arg is .fai
        if sys.argv[2].endswith('.fai'):
            fai_path = sys.argv[2]
            out_beds = sys.argv[3:]
        else:
            fai_path = None
            out_beds = sys.argv[2:]
        split_bed(in_bed, out_beds, fai_path)
    else:
        print(f"Usage: {sys.argv[0]} <input.bed> [ref.fai] <out_chunk1.bed> <out_chunk2.bed> ...")
        sys.exit(1)
