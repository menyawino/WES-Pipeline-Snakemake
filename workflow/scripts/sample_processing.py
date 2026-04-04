#!/usr/bin/env python

'''Functions to process and validate sample data for the pipeline.'''

import pandas as pd
import glob
import os
import re
import sys
import logging

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def get_sample_data(csv_file, input_dir, lanes=[1, 2, 3, 4], fail_on_missing=False):
    """
    Get the sample data needed for the pipeline.
    
    Args:
        csv_file: Path to sample metadata CSV
        input_dir: Path to input directory containing FASTQ files
        lanes: List of lane numbers to expect (default: [1,2,3,4])
        fail_on_missing: If True, exit on missing files; if False, continue with warning
    
    Returns:
        DataFrame with sample metadata and file paths
    """
    df = pd.read_csv(csv_file)
    sample_names = df['sample'].astype(str).tolist()

    sample_fastq_files = {}
    all_samples_in_dir = set(os.listdir(input_dir))

    missing_files = []
    not_found_samples = []

    for sample in sample_names:
        if sample not in all_samples_in_dir:
            not_found_samples.append(sample)
            continue
        
        sample_fastq_files[sample] = []
        
        for lane in lanes:
            # Build glob patterns for R1 and R2 files
            r1_pattern = os.path.join(input_dir, sample, "{}_S*_L00{}_R1_*.fastq.gz".format(sample, lane))
            r2_pattern = os.path.join(input_dir, sample, "{}_S*_L00{}_R2_*.fastq.gz".format(sample, lane))
            
            r1_files = glob.glob(r1_pattern)
            r2_files = glob.glob(r2_pattern)

            if not r1_files:
                missing_files.append(f"Missing R1 file for sample {sample}, lane {lane}")
            if not r2_files:
                missing_files.append(f"Missing R2 file for sample {sample}, lane {lane}")

            for file in r1_files:
                match = re.search(r'_S(\d+)_L00(\d)_R1_', file)
                if match:
                    sample_number = match.group(1)
                    matched_lane = match.group(2)
                    sample_fastq_files[sample].append((sample_number, f'L00{matched_lane}', 'R1', file))
            
            for file in r2_files:
                match = re.search(r'_S(\d+)_L00(\d)_R2_', file)
                if match:
                    sample_number = match.group(1)
                    matched_lane = match.group(2)
                    sample_fastq_files[sample].append((sample_number, f'L00{matched_lane}', 'R2', file))

    # Log validation results
    logger.info(f"Found {len(sample_fastq_files)} samples with FASTQ files.")
    
    extra_samples = all_samples_in_dir - set(sample_names)
    if extra_samples:
        logger.warning(f"{len(extra_samples)} extra samples in directory not in sample list: {', '.join(extra_samples)}")

    if not_found_samples:
        msg = f"{len(not_found_samples)} samples not found: {', '.join(not_found_samples)}"
        if fail_on_missing:
            logger.error(msg)
            sys.exit(1)
        else:
            logger.warning(msg)

    if missing_files:
        msg = f"{len(missing_files)} missing files detected"
        if fail_on_missing:
            for f in missing_files[:5]:  # Show first 5
                logger.error(f)
            if len(missing_files) > 5:
                logger.error(f"... and {len(missing_files) - 5} more")
            sys.exit(1)
        else:
            logger.warning(msg)
            for f in missing_files[:5]:  # Show first 5
                logger.warning(f)
            if len(missing_files) > 5:
                logger.warning(f"... and {len(missing_files) - 5} more")

    # Merge metadata with fastq file information
    metadata_columns = df.columns.tolist()
    if 'sample' in metadata_columns:
        metadata_columns.remove('sample')

    output_rows = []
    for index, row in df.iterrows():
        sample = str(row['sample'])
        if sample not in sample_fastq_files:
            continue
        
        metadata = [row[col] for col in metadata_columns]
        for sample_number, lane, read, file in sample_fastq_files[sample]:
            merged_sample = f"{sample}_S{sample_number}"
            output_rows.append([merged_sample, lane, read, file] + metadata)

    # Save results to CSV
    output_file = 'sample_data.csv'
    with open(output_file, 'w') as f:
        header = ['sample', 'lane', 'read', 'file'] + metadata_columns
        f.write(','.join(header) + '\n')
        for row in output_rows:
            f.write(','.join(map(str, row)) + '\n')
    
    logger.info(f"Sample data saved to {output_file}")
    return pd.read_csv(output_file)