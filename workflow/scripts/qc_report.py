import pandas as pd
import re
from snakemake.shell import shell

def parse_flagstat(file_path):
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        total_match = re.search(r'(\d+) \+ \d+ in total', content)
        mapped_match = re.search(r'(\d+) \+ \d+ mapped', content)
        total_reads = int(total_match.group(1)) if total_match else 0
        mapped_reads = int(mapped_match.group(1)) if mapped_match else 0
        pct_mapped = round(mapped_reads / total_reads * 100, 2) if total_reads > 0 else 0.0
        return {
            'TotalReads': total_reads,
            'MappedReads': mapped_reads,
            '%Mapped': pct_mapped
        }
    except Exception:
        return {'TotalReads': 0, 'MappedReads': 0, '%Mapped': 0.0}

def parse_coverage_stats(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')
        if df.empty:
            return {'ReadsOnTarget_q8': 0, '%OnTarget': 0.0, 'UniqReadsOnTarget_q8': 0, '%UniqueReadsOnTarget_q8': 0.0, 'TargetSize(bp)': 0}
        return {
            'ReadsOnTarget_q8': df['reads_on_target'].iloc[0],
            '%OnTarget': round(df['pct_on_target'].iloc[0] * 100, 2),
            'UniqReadsOnTarget_q8': df['uniq_reads_on_target'].iloc[0],
            '%UniqueReadsOnTarget_q8': round(df['pct_uniq_on_target'].iloc[0] * 100, 2),
            'TargetSize(bp)': df['target_size'].iloc[0]
        }
    except Exception:
        return {'ReadsOnTarget_q8': 0, '%OnTarget': 0.0, 'UniqReadsOnTarget_q8': 0, '%UniqueReadsOnTarget_q8': 0.0, 'TargetSize(bp)': 0}

def parse_coverage_hist(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t')
        if df.empty or 'depth' not in df.columns or 'bases' not in df.columns:
            return {('Bases=0x' if d == 0 else f'Bases>={d}x'): 0 for d in [0, 1, 5, 10, 20, 30, 50]}
        depths = df['depth'].values
        bases = df['bases'].values
        coverage_metrics = {}
        for depth in [0, 1, 5, 10, 20, 30, 50]:
            key = 'Bases=0x' if depth == 0 else f'Bases>={depth}x'
            coverage_metrics[key] = int(bases[depths >= depth].sum())
        return coverage_metrics
    except Exception:
        return {('Bases=0x' if d == 0 else f'Bases>={d}x'): 0 for d in [0, 1, 5, 10, 20, 30, 50]}

def parse_depth_of_coverage(file_path):
    try:
        with open(file_path, 'r') as f:
            content = f.readlines()
        if len(content) < 2:
            return {'MeanCov': 0.0, 'MedianCov': 0.0, 'PCT_TARGET_BASES_1X (updated)': 0.0, 'PCT_TARGET_BASES_10X (updated)': 0.0, 'PCT_TARGET_BASES_20X (updated)': 0.0, 'PCT_TARGET_BASES_30X (updated)': 0.0, 'GC_DROPOUT (updated)': 0.0, 'AT_DROPOUT (updated)': 0.0}
        sample_line = content[1].split('\t')
        return {
            'MeanCov': float(sample_line[2]) if len(sample_line) > 2 else 0.0,
            'MedianCov': float(sample_line[3]) if len(sample_line) > 3 else 0.0,
            'PCT_TARGET_BASES_1X (updated)': float(sample_line[11]) if len(sample_line) > 11 else 0.0,
            'PCT_TARGET_BASES_10X (updated)': float(sample_line[13]) if len(sample_line) > 13 else 0.0,
            'PCT_TARGET_BASES_20X (updated)': float(sample_line[15]) if len(sample_line) > 15 else 0.0,
            'PCT_TARGET_BASES_30X (updated)': float(sample_line[17]) if len(sample_line) > 17 else 0.0,
            'GC_DROPOUT (updated)': float(sample_line[-2]) if len(sample_line) > 2 else 0.0,
            'AT_DROPOUT (updated)': float(sample_line[-1]) if len(sample_line) > 1 else 0.0
        }
    except Exception:
        return {'MeanCov': 0.0, 'MedianCov': 0.0, 'PCT_TARGET_BASES_1X (updated)': 0.0, 'PCT_TARGET_BASES_10X (updated)': 0.0, 'PCT_TARGET_BASES_20X (updated)': 0.0, 'PCT_TARGET_BASES_30X (updated)': 0.0, 'GC_DROPOUT (updated)': 0.0, 'AT_DROPOUT (updated)': 0.0}

def parse_alignment_summary_metrics(file_path):
    try:
        df = pd.read_csv(file_path, sep='\t', comment='#', skiprows=6)
        first_pair = df[df['CATEGORY'] == 'FIRST_OF_PAIR'].iloc[0] if not df[df['CATEGORY'] == 'FIRST_OF_PAIR'].empty else {}
        second_pair = df[df['CATEGORY'] == 'SECOND_OF_PAIR'].iloc[0] if not df[df['CATEGORY'] == 'SECOND_OF_PAIR'].empty else {}
        pair = df[df['CATEGORY'] == 'PAIR'].iloc[0] if not df[df['CATEGORY'] == 'PAIR'].empty else {}
        return {
            'MeanFwdReadLength': first_pair.get('MEAN_READ_LENGTH', 0),
            'MeanRevReadLength': second_pair.get('MEAN_READ_LENGTH', 0),
            'ReadsAlignedInPairs': pair.get('READS_ALIGNED_IN_PAIRS', 0),
            '%ReadsAlignedInPairs': round(pair.get('PCT_READS_ALIGNED_IN_PAIRS', 0) * 100, 2),
            'StrandBalance': pair.get('STRAND_BALANCE', 0),
            'PCT_PF_READS_ALIGNED (updated)': round(pair.get('PCT_PF_READS_ALIGNED', 0) * 100, 2),
            'PCT_CHIMERAS (updated)': round(pair.get('PCT_CHIMERAS', 0) * 100, 2),
            'PCT_ADAPTER (updated)': round(pair.get('PCT_ADAPTER', 0) * 100, 2)
        }
    except Exception:
        return {'MeanFwdReadLength': 0, 'MeanRevReadLength': 0, 'ReadsAlignedInPairs': 0, '%ReadsAlignedInPairs': 0.0, 'StrandBalance': 0, 'PCT_PF_READS_ALIGNED (updated)': 0.0, 'PCT_CHIMERAS (updated)': 0.0, 'PCT_ADAPTER (updated)': 0.0}

def parse_mean_coverage(file_path):
    df = pd.read_csv(file_path, sep='\t', header=None, names=['chrom', 'start', 'end', 'coverage'])
    return {
        'MEAN_COVERAGE (updated)': df['coverage'].mean()
    }

def main(snakemake):
    metrics = {}
    
    # Parse flagstat files
    metrics.update(parse_flagstat(snakemake.input.flagstat_original))
    metrics.update({'Target_' + k: v for k, v in parse_flagstat(snakemake.input.flagstat_target).items()})
    
    # Parse coverage stats files
    metrics.update(parse_coverage_stats(snakemake.input.coverage_stats))
    metrics.update({'Target_' + k: v for k, v in parse_coverage_stats(snakemake.input.coverage_stats_target).items()})
    
    # Parse coverage hist files
    metrics.update(parse_coverage_hist(snakemake.input.coverage_hist))
    metrics.update({'Target_' + k: v for k, v in parse_coverage_hist(snakemake.input.coverage_hist_target).items()})
    
    # Parse depth of coverage files
    metrics.update(parse_depth_of_coverage(snakemake.input.depth_of_coverage))
    metrics.update({'Target_' + k: v for k, v in parse_depth_of_coverage(snakemake.input.depth_of_coverage_target).items()})
    
    # Parse alignment summary metrics files
    metrics.update(parse_alignment_summary_metrics(snakemake.input.alignment_summary_metrics))
    metrics.update({'Target_' + k: v for k, v in parse_alignment_summary_metrics(snakemake.input.alignment_summary_metrics_target).items()})
    
    # Parse mean coverage files
    metrics.update(parse_mean_coverage(snakemake.input.mean_coverage))
    metrics.update({'Target_' + k: v for k, v in parse_mean_coverage(snakemake.input.mean_coverage_target).items()})
    
    # Create DataFrame and save to TSV
    df = pd.DataFrame([metrics])
    df.to_csv(snakemake.output.qc_metrics, sep='\t', index=False)

if __name__ == "__main__":
    main(snakemake)
