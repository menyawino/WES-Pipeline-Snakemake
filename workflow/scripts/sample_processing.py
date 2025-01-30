#!/usr/bin/env python

'''A function to get the sample data needed for the pipeline.'''

import pandas as pd
import glob
import os
import re
import sys

def get_sample_data(csv_file, input_dir):
    """
    A function to get the sample data needed for the pipeline.
    """
    df = pd.read_csv(csv_file)

    # Assuming the sample names are in a column named 'sample'
    sample_names = df['sample'].astype(str).tolist()

    # Dictionary to hold the fastq file paths for each sample
    sample_fastq_files = {}
    all_samples_in_dir = set(os.listdir(input_dir))  # Get all samples present in the input directory

    missing_files = []
    not_found_samples = []

    for sample in sample_names:
        if sample not in all_samples_in_dir:
            not_found_samples.append(sample)
            continue
        sample_fastq_files[sample] = []
        for lane in range(1, 5):  # Loop through lanes 1 to 4 (1, 5), for deployment change to (1, 2)
            # Update the pattern to include sample directories
            r1_pattern = os.path.join(input_dir, sample, "{}_S*_L00{}_R1_*.fastq.gz".format(sample, lane))
            r2_pattern = os.path.join(input_dir, sample, "{}_S*_L00{}_R2_*.fastq.gz".format(sample, lane))
            
            r1_files = glob.glob(r1_pattern)
            r2_files = glob.glob(r2_pattern)

            if not r1_files:
                missing_files.append("Missing R1 file for sample {}, lane {}".format(sample, lane))
            if not r2_files:
                missing_files.append("Missing R2 file for sample {}, lane {}".format(sample, lane))

            for file in r1_files:
                match = re.search(r'_S(\d+)_L00(\d)_R1_', file)
                if match:
                    sample_number = match.group(1)
                    matched_lane = match.group(2)  # Use a different variable name for the matched lane
                    sample_fastq_files[sample].append((sample_number, 'L00{}'.format(matched_lane), 'R1', file))
            
            for file in r2_files:
                match = re.search(r'_S(\d+)_L00(\d)_R2_', file)
                if match:
                    sample_number = match.group(1)
                    matched_lane = match.group(2)  # Use a different variable name for the matched lane
                    sample_fastq_files[sample].append((sample_number, 'L00{}'.format(matched_lane), 'R2', file))

    # Check for missing files
    if missing_files:
        for missing_file in missing_files:
            print(missing_file)
            
        # give option to continue without the missing files: Do you want to continue without the missing files? (y/n)
        print("Do you want to continue without the missing files? (y/n)")
        response = input()
        if response.lower() == 'n':
            print("Pipeline terminated.")
            sys.exit(1)
        elif response.lower() == 'y':
            print("Continuing without the missing files.")
        else:
            raise ValueError("Invalid response. Please enter 'y' or 'n'.")
        
    # Print the number of found samples and samples not mentioned in the input list but present in the input directory
    found_samples = len(sample_fastq_files)
    extra_samples = all_samples_in_dir - set(sample_names)
    print(f"Found {found_samples} samples.")
    print(f"{len(extra_samples)} samples are not mentioned in the input list but present in the input directory: {', '.join(extra_samples)}")

    # Check for samples in the input list but not found in the input directory
    if not_found_samples:
        print(f"{len(not_found_samples)} samples are mentioned in the input list but not found in the input directory: {', '.join(not_found_samples)}")
        print("Do you want to start the analysis anyway? (y/n)")
        response = input()
        if response.lower() == 'n':
            print("Pipeline terminated.")
            sys.exit(1)
        elif response.lower() == 'y':
            print("Continuing without the missing samples.")
        else:
            raise ValueError("Invalid response. Please enter 'y' or 'n'.")

    # Merge the metadata with the fastq file information
    metadata_columns = df.columns.tolist()
    metadata_columns.remove('sample')

    output_rows = []

    for index, row in df.iterrows():
        sample = str(row['sample'])
        metadata = [row[col] for col in metadata_columns]
        
        for sample_number, lane, read, file in sample_fastq_files[sample]:
            merged_sample = "{}_S{}".format(sample, sample_number)
            output_rows.append([merged_sample, lane, read, file] + metadata)

    # Save the results to a CSV file
    output_file = 'sample_data.csv'
    with open(output_file, 'w') as f:
        header = ['sample', 'lane', 'read', 'file'] + metadata_columns
        f.write(','.join(header) + '\n')
        
        for row in output_rows:
            f.write(','.join(map(str, row)) + '\n')
    
    # print("Results have been saved to {}".format(output_file))
    
    # Return the CSV file as a dataframe
    return pd.read_csv(output_file)