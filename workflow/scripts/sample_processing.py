#!/usr/bin/env python3

'''Functions to process and validate sample data for the pipeline (supports NextSeq, MiSeq, NovaSeq).'''

import pandas as pd
import os
import re
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def resolve_effective_input_dir(input_dir):
    """Detect if input_dir is a root Illumina run folder with Data/Intensities/BaseCalls."""
    basecalls_dir = os.path.join(input_dir, "Data", "Intensities", "BaseCalls")
    if os.path.isdir(basecalls_dir):
        logger.info(f"Detected Illumina run folder. Using BaseCalls directory: {basecalls_dir}")
        return basecalls_dir
    return input_dir

def get_sample_data(csv_file, input_dir, lanes=[1, 2, 3, 4], fail_on_missing=False):
    """
    Get sample metadata and FASTQ file mapping without modifying the input directory.
    Supports nested sample directories, flat BaseCalls directories, and automated sample discovery.
    """
    effective_dir = resolve_effective_input_dir(input_dir)
    
    # 1. Fast scan for all .fastq.gz in effective_dir (flat + immediate subdirectories)
    all_fastqs = []
    try:
        for entry in os.scandir(effective_dir):
            if entry.is_file(follow_symlinks=False) and entry.name.endswith(".fastq.gz"):
                all_fastqs.append(entry.path)
            elif entry.is_dir(follow_symlinks=False):
                try:
                    for sub in os.scandir(entry.path):
                        if sub.is_file(follow_symlinks=False) and sub.name.endswith(".fastq.gz"):
                            all_fastqs.append(sub.path)
                except (PermissionError, FileNotFoundError):
                    pass
    except (PermissionError, FileNotFoundError) as e:
        logger.warning(f"Error scanning input directory {effective_dir}: {e}")

    # Map discovered FASTQ files by prefix
    discovered_samples = {}
    for f in all_fastqs:
        fname = os.path.basename(f)
        if fname.startswith("Undetermined"):
            continue
        # Standard Illumina naming: <SampleID>_S<SampleNum>_L00<Lane>_R<Read>_001.fastq.gz
        m = re.search(r'^(.+?)_S(\d+)_L00(\d)_R([12])_(\d+)\.fastq\.gz$', fname)
        if m:
            s_name = m.group(1)
            s_num = m.group(2)
            s_lane = f"L00{m.group(3)}"
            s_read = f"R{m.group(4)}"
            
            if s_name not in discovered_samples:
                discovered_samples[s_name] = []
            discovered_samples[s_name].append((s_num, s_lane, s_read, f))

    # 2. Check if a valid sample CSV was provided
    sample_names = []
    df = None
    if csv_file and os.path.exists(csv_file):
        try:
            df = pd.read_csv(csv_file)
            if 'sample' in df.columns:
                cand_names = df['sample'].astype(str).tolist()
                # Check if at least one candidate matches discovered fastqs
                if any(c in discovered_samples for c in cand_names):
                    sample_names = cand_names
        except Exception as e:
            logger.warning(f"Could not parse CSV {csv_file}: {e}")

    # 3. If no matching sample_names from CSV, use all discovered FASTQ sample prefixes
    if not sample_names:
        sample_names = sorted(list(discovered_samples.keys()))
        df = pd.DataFrame({'sample': sample_names, 'condition': ['unknown'] * len(sample_names)})
        logger.info(f"Using {len(sample_names)} sample prefixes discovered directly from FASTQ files.")

    sample_fastq_files = {}
    for sample in sample_names:
        matched_records = discovered_samples.get(sample, [])
        if matched_records:
            sample_fastq_files[sample] = matched_records

    logger.info(f"Found {len(sample_fastq_files)} matching samples with FASTQ files.")

    if not sample_fastq_files and fail_on_missing:
        logger.error(f"No FASTQ files matched in {effective_dir}")
        sys.exit(1)

    # 4. Construct result DataFrame
    metadata_columns = [col for col in (df.columns.tolist() if df is not None else ['condition']) if col != 'sample']
    output_rows = []
    
    for index, row in df.iterrows():
        sample = str(row['sample'])
        if sample not in sample_fastq_files:
            continue
        metadata = [row[col] for col in metadata_columns]
        for s_num, lane, read, file in sample_fastq_files[sample]:
            merged_sample = f"{sample}_S{s_num}"
            output_rows.append([merged_sample, lane, read, file] + metadata)

    header = ['sample', 'lane', 'read', 'file'] + metadata_columns
    result_df = pd.DataFrame(output_rows, columns=header)
    logger.info(f"Processed {len(result_df)} sample FASTQ records.")
    return result_df

if __name__ == '__main__':
    if len(sys.argv) > 2:
        res = get_sample_data(sys.argv[1], sys.argv[2])
        print(res.head())