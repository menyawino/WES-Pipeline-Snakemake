#!/usr/bin/env python3

'''Functions to process and validate sample data for the pipeline (supports NextSeq, MiSeq, NovaSeq).'''

import pandas as pd
import os
import re
import sys
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)

def find_all_fastqs(input_dir, max_depth=5):
    """Recursively search for all .fastq.gz files in input directory up to max_depth."""
    all_fastqs = []
    if not os.path.exists(input_dir):
        return all_fastqs
    base_depth = input_dir.rstrip(os.path.sep).count(os.path.sep)
    for root, dirs, files in os.walk(input_dir, followlinks=True):
        cur_depth = root.rstrip(os.path.sep).count(os.path.sep) - base_depth
        if cur_depth >= max_depth:
            dirs.clear()
            continue
        dirs[:] = [d for d in dirs if not d.startswith('.')]
        for f in files:
            if f.endswith(".fastq.gz") and not f.startswith("._"):
                all_fastqs.append(os.path.join(root, f))
    return all_fastqs

def discover_samples_from_analysis(outdir):
    """Scan existing analysis directory to infer sample naming and _S numbers."""
    discovered = {}
    if not outdir or not os.path.exists(outdir):
        return discovered
    for check_dir in ["analysis/001_QC/pretrim", "analysis/002_trimming", "analysis/003_alignment"]:
        full_dir = os.path.join(outdir, check_dir)
        if os.path.exists(full_dir):
            for root, _, files in os.walk(full_dir):
                for f in files:
                    m = re.search(r'^(\d+|[A-Za-z0-9_-]+)_S(\d+)_L00(\d)_', f)
                    if m:
                        s_name = m.group(1)
                        s_num = m.group(2)
                        s_lane = f"L00{m.group(3)}"
                        if s_name not in discovered:
                            discovered[s_name] = {}
                        discovered[s_name][s_lane] = s_num
    return discovered

def get_sample_data(csv_file, input_dir, lanes=[1, 2, 3, 4], fail_on_missing=False, outdir=None):
    """
    Get sample metadata and FASTQ file mapping without modifying the input directory.
    Supports nested sample directories, NextSeq/MiSeq run folders, and automated sample discovery.
    """
    all_fastqs = find_all_fastqs(input_dir)
    logger.info(f"Discovered {len(all_fastqs)} total FASTQ files in {input_dir}")

    # Map discovered FASTQ files by prefix
    discovered_samples = {}
    for f in all_fastqs:
        fname = os.path.basename(f)
        if fname.startswith("Undetermined"):
            continue
        m = re.search(r'^(.+?)_S(\d+)_L00(\d)_R([12])_(\d+)\.fastq\.gz$', fname)
        if m:
            s_name = m.group(1)
            s_num = m.group(2)
            s_lane = f"L00{m.group(3)}"
            s_read = f"R{m.group(4)}"
            
            if s_name not in discovered_samples:
                discovered_samples[s_name] = []
            discovered_samples[s_name].append((s_num, s_lane, s_read, f))

    # Also discover from analysis outputs if raw fastqs not present
    analysis_samples = discover_samples_from_analysis(outdir or os.path.join(os.path.dirname(csv_file or ""), "..", "output_38"))

    # 2. Check if a valid sample CSV was provided
    sample_names = []
    df = None
    if csv_file and os.path.exists(csv_file):
        try:
            df = pd.read_csv(csv_file)
            if 'sample' in df.columns:
                sample_names = [str(x).strip() for x in df['sample'].dropna().tolist() if str(x).strip()]
        except Exception as e:
            logger.warning(f"Could not parse CSV {csv_file}: {e}")

    # 3. If no sample CSV, discover from FASTQ files
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
        logger.error(f"No FASTQ files matched in {input_dir}")
        sys.exit(1)

    # 4. Construct result DataFrame
    metadata_columns = [col for col in (df.columns.tolist() if df is not None else ['condition']) if col != 'sample']
    output_rows = []
    
    # Format lanes list
    formatted_lanes = []
    for l in lanes:
        str_l = str(l)
        if not str_l.startswith("L00"):
            str_l = f"L00{str_l}"
        formatted_lanes.append(str_l)

    for index, row in (df.iterrows() if df is not None else []):
        sample = str(row['sample']).strip()
        metadata = [row[col] for col in metadata_columns]
        if sample in sample_fastq_files:
            for s_num, lane, read, file in sample_fastq_files[sample]:
                merged_sample = f"{sample}_S{s_num}"
                output_rows.append([merged_sample, lane, read, file] + metadata)
        elif sample in analysis_samples:
            lane_map = analysis_samples[sample]
            for lane in formatted_lanes:
                s_num = lane_map.get(lane, list(lane_map.values())[0] if lane_map else str(index + 1))
                merged_sample = f"{sample}_S{s_num}" if not sample.endswith(f"_S{s_num}") else sample
                for read in ["R1", "R2"]:
                    f_name = f"{sample}_S{s_num}_{lane}_{read}_001.fastq.gz"
                    dummy_file = os.path.join(input_dir, f_name)
                    output_rows.append([merged_sample, lane, read, dummy_file] + metadata)
        else:
            s_num = index + 1
            merged_sample = f"{sample}_S{s_num}" if not sample.endswith(f"_S{s_num}") else sample
            for lane in formatted_lanes:
                for read in ["R1", "R2"]:
                    f_name = f"{sample}_S{s_num}_{lane}_{read}_001.fastq.gz"
                    dummy_file = os.path.join(input_dir, f_name)
                    output_rows.append([merged_sample, lane, read, dummy_file] + metadata)

    header = ['sample', 'lane', 'read', 'file'] + metadata_columns
    result_df = pd.DataFrame(output_rows, columns=header)
    logger.info(f"Processed {len(result_df)} sample FASTQ records.")
    return result_df

if __name__ == '__main__':
    if len(sys.argv) > 2:
        res = get_sample_data(sys.argv[1], sys.argv[2])
        print(res.head())