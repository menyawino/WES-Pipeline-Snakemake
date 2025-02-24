import pandas as pd
import re
from snakemake.shell import shell

def parse_flagstat(file_path):
    with open(file_path, 'r') as f:
        content = f.read()
    total_reads = int(re.search(r'(\d+) \+ \d+ in total', content).group(1))
    mapped_reads = int(re.search(r'(\d+) \+ \d+ mapped', content).group(1))
    return {
        'TotalReads': total_reads,
        'MappedReads': mapped_reads,
        '%Mapped': round(mapped_reads / total_reads * 100, 2)
    }

def parse_coverage_stats(file_path):
    df = pd.read_csv(file_path, sep='\t')
    return {
        'ReadsOnTarget_q8': df['reads_on_target'].iloc[0],
        '%OnTarget': round(df['pct_on_target'].iloc[0] * 100, 2),
        'UniqReadsOnTarget_q8': df['uniq_reads_on_target'].iloc[0],
        '%UniqueReadsOnTarget_q8': round(df['pct_uniq_on_target'].iloc[0] * 100, 2),
        'TargetSize(bp)': df['target_size'].iloc[0]
    }

def parse_coverage_hist(file_path):
    df = pd.read_csv(file_path, sep='\t')
    coverage_metrics = {}
    for depth in [0, 1, 5, 10, 20, 30, 50]:
        key = f'Bases>={depth}x'
        if depth == 0:
            key = 'Bases=0x'
        coverage_metrics[key] = df[df['depth'] >= depth]['bases'].sum()
    return coverage_metrics

def parse_depth_of_coverage(file_path):
    with open(file_path, 'r') as f:
        content = f.readlines()
    sample_line = content[1].split('\t')
    return {
        'MeanCov': float(sample_line[2]),
        'MedianCov': float(sample_line[3]),
        'PCT_TARGET_BASES_1X (updated)': float(sample_line[11]),
        'PCT_TARGET_BASES_10X (updated)': float(sample_line[13]),
        'PCT_TARGET_BASES_20X (updated)': float(sample_line[15]),
        'PCT_TARGET_BASES_30X (updated)': float(sample_line[17]),
        'GC_DROPOUT (updated)': float(sample_line[-2]),
        'AT_DROPOUT (updated)': float(sample_line[-1])
    }

def parse_alignment_summary_metrics(file_path):
    df = pd.read_csv(file_path, sep='\t', comment='#', skiprows=6)
    first_pair = df[df['CATEGORY'] == 'FIRST_OF_PAIR'].iloc[0]
    second_pair = df[df['CATEGORY'] == 'SECOND_OF_PAIR'].iloc[0]
    pair = df[df['CATEGORY'] == 'PAIR'].iloc[0]
    return {
        'MeanFwdReadLength': first_pair['MEAN_READ_LENGTH'],
        'MeanRevReadLength': second_pair['MEAN_READ_LENGTH'],
        'ReadsAlignedInPairs': pair['READS_ALIGNED_IN_PAIRS'],
        '%ReadsAlignedInPairs': round(pair['PCT_READS_ALIGNED_IN_PAIRS'] * 100, 2),
        'StrandBalance': pair['STRAND_BALANCE'],
        'PCT_PF_READS_ALIGNED (updated)': round(pair['PCT_PF_READS_ALIGNED'] * 100, 2),
        'PCT_CHIMERAS (updated)': round(pair['PCT_CHIMERAS'] * 100, 2),
        'PCT_ADAPTER (updated)': round(pair['PCT_ADAPTER'] * 100, 2)
    }

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
