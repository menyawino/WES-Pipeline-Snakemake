#!/usr/bin/env python3
"""
Split a target BED file into N balanced BED files based on genomic interval sizes.
Ensures reference coordinate ordering and uniform load balancing for scatter-gather GATK HaplotypeCaller.
"""

import sys
import os

def get_contig_order(fai_path=None):
    """Load contig ordering map from reference FAI if available, or fallback to standard human order."""
    if fai_path and os.path.exists(fai_path):
        contig_order = {}
        with open(fai_path, 'r') as f:
            for i, line in enumerate(f):
                parts = line.strip().split('\t')
                if parts:
                    contig_order[parts[0]] = i
        return contig_order

    # Standard fallback
    standard = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY", "chrM", "chrMT"]
    standard += [str(i) for i in range(1, 23)] + ["X", "Y", "MT", "M"]
    return {c: i for i, c in enumerate(standard)}

def split_bed(bed_path, output_paths, fai_path=None):
    """Split input BED file into balanced sub-BED files sorted by reference order."""
    contig_order = get_contig_order(fai_path)
    intervals = []
    
    with open(bed_path, 'r') as f:
        for line in f:
            line_str = line.strip()
            if not line_str or line_str.startswith(('#', 'track', 'browser')):
                continue
            parts = line_str.split('\t')
            if len(parts) >= 3:
                chrom = parts[0]
                start = int(parts[1])
                end = int(parts[2])
                length = max(1, end - start)
                rest = parts[3:] if len(parts) > 3 else []
                intervals.append((chrom, start, end, length, rest))
                
    if not intervals:
        for p in output_paths:
            os.makedirs(os.path.dirname(os.path.abspath(p)), exist_ok=True)
            open(p, 'w').close()
        return

    # Sort strictly by reference contig order, then by start position
    intervals.sort(key=lambda item: (contig_order.get(item[0], 999999), item[1], item[2]))

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
            for chrom, start, end, _, rest in chunk_items:
                if rest:
                    f_out.write(f"{chrom}\t{start}\t{end}\t" + "\t".join(rest) + "\n")
                else:
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
