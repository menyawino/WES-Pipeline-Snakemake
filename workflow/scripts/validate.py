#!/usr/bin/env python

'''Input validation utilities for the pipeline.'''

import os
import sys
import logging
import pandas as pd
from pathlib import Path

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def validate_reference_files(config):
    """
    Validate that all required reference files exist and are properly indexed.
    If the reference genome is missing, automatically download and index GRCh38.
    
    Args:
        config: Dictionary containing reference file paths
    
    Returns:
        bool: True if all validations pass
    """
    logger.info("Validating reference files...")
    errors = []
    
    ref_path = config.get('reference_genome')
    if not ref_path:
        errors.append("Missing configuration key: reference_genome")
    elif not os.path.exists(ref_path):
        logger.info(f"Reference genome not found at {ref_path}. Attempting automatic download of GRCh38...")
        try:
            from workflow.scripts.download_ref import download_reference_genome
        except ImportError:
            try:
                from download_ref import download_reference_genome
            except ImportError:
                sys.path.append(os.path.dirname(__file__))
                from download_ref import download_reference_genome
                
        success = download_reference_genome(target_path=ref_path, genome=config.get('Genome', 'grch38'))
        if not success or not os.path.exists(ref_path):
            errors.append(f"Reference genome not found and download failed: {ref_path}")
        else:
            logger.info(f"✓ Reference genome downloaded and verified: {ref_path}")

    if ref_path and os.path.exists(ref_path):
        logger.info(f"✓ Found Reference genome: {ref_path}")
        fasta_index = f"{ref_path}.fai"
        if not os.path.exists(fasta_index):
            logger.info(f"Index missing for {ref_path}. Generating index...")
            try:
                from workflow.scripts.download_ref import index_reference
            except ImportError:
                from download_ref import index_reference
            index_reference(ref_path)

        icc_panel = config.get('icc_panel')
        if icc_panel and os.path.exists(icc_panel):
            try:
                from workflow.scripts.download_ref import validate_bed_compatibility
            except ImportError:
                from download_ref import validate_bed_compatibility
            validate_bed_compatibility(ref_path, icc_panel)

    other_refs = {
        'icc_panel': 'ICC panel BED file',
        'dbsnp': 'dbSNP variants',
    }
    
    for ref_key, ref_name in other_refs.items():
        path = config.get(ref_key)
        if not path:
            errors.append(f"Missing configuration key: {ref_key}")
        elif not os.path.exists(path):
            logger.warning(f"  Optional/Panel reference dataset not found at: {path} (will be created or passed if rules run)")
        else:
            logger.info(f"✓ Found {ref_name}: {path}")
    
    if errors:
        logger.error("Reference file validation failed:")
        for error in errors:
            logger.error(f"  ✗ {error}")
        return False
    
    logger.info("✓ All required reference files validated")
    return True


def validate_directories(inputdir, outdir):
    """
    Validate input and output directories.
    
    Args:
        inputdir: Path to input directory
        outdir: Path to output directory
    
    Returns:
        bool: True if directories are valid and writable
    """
    logger.info("Validating directories...")
    errors = []
    
    # Check input directory
    if not os.path.isdir(inputdir):
        errors.append(f"Input directory not found: {inputdir}")
    elif not os.access(inputdir, os.R_OK):
        errors.append(f"Input directory not readable: {inputdir}")
    else:
        logger.info(f"✓ Input directory accessible: {inputdir}")
    
    # Check/create output directory
    try:
        os.makedirs(outdir, exist_ok=True)
        if not os.access(outdir, os.W_OK):
            errors.append(f"Output directory not writable: {outdir}")
        else:
            logger.info(f"✓ Output directory writable: {outdir}")
    except Exception as e:
        errors.append(f"Cannot create output directory: {outdir} ({e})")
    
    if errors:
        logger.error("Directory validation failed:")
        for error in errors:
            logger.error(f"  ✗ {error}")
        return False
    
    logger.info("✓ All directories validated")
    return True


def validate_sample_file(sample_file):
    """
    Validate sample metadata CSV file.
    
    Args:
        sample_file: Path to sample CSV file
    
    Returns:
        bool: True if sample file is valid
    """
    logger.info("Validating sample file...")
    
    if not os.path.exists(sample_file):
        logger.error(f"Sample file not found: {sample_file}")
        return False
    
    try:
        df = pd.read_csv(sample_file)
        
        if 'sample' not in df.columns:
            logger.error("Sample CSV must have 'sample' column")
            return False
        
        if df.empty:
            logger.error("Sample CSV is empty")
            return False
        
        logger.info(f"✓ Sample file contains {len(df)} samples")
        logger.info(f"✓ Columns: {', '.join(df.columns)}")
        
        return True
    
    except Exception as e:
        logger.error(f"Error reading sample file: {e}")
        return False


def validate_fastq_files(inputdir, sample_df, lanes=[1, 2, 3, 4]):
    """
    Validate that expected FASTQ files exist for all samples.
    
    Args:
        inputdir: Path to input directory
        sample_df: DataFrame with sample information
        lanes: List of lane numbers to check
    
    Returns:
        tuple: (bool, list) - validation result and list of missing files
    """
    logger.info("Validating FASTQ files...")
    import glob
    
    missing = []
    found_count = 0
    
    for idx, row in sample_df.iterrows():
        sample = str(row['sample'])
        sample_path = os.path.join(inputdir, sample)
        
        if not os.path.isdir(sample_path):
            missing.append(f"Sample directory not found: {sample_path}")
            continue
        
        for lane in lanes:
            for read in ['R1', 'R2']:
                pattern = os.path.join(sample_path, f"{sample}_S*_L00{lane}_{read}_*.fastq.gz")
                files = glob.glob(pattern)
                
                if not files:
                    missing.append(f"No {read} files found for {sample} lane {lane}")
                else:
                    found_count += 1
    
    if missing:
        logger.warning(f"Found {len(missing)} missing/problematic FASTQ files:")
        for m in missing[:5]:
            logger.warning(f"  ✗ {m}")
        if len(missing) > 5:
            logger.warning(f"  ... and {len(missing) - 5} more")
        return False, missing
    
    logger.info(f"✓ Found {found_count} FASTQ file pairs")
    return True, []


def validate_pipeline_config(configfile, inputdir, outdir):
    """
    Comprehensive pipeline configuration validation.
    
    Args:
        configfile: Path to Snakemake config YAML
        inputdir: Path to input directory
        outdir: Path to output directory
    
    Returns:
        bool: True if all validations pass
    """
    import yaml
    
    logger.info("=" * 60)
    logger.info("PIPELINE PRE-FLIGHT VALIDATION")
    logger.info("=" * 60)
    
    # Load config
    try:
        with open(configfile, 'r') as f:
            config = yaml.safe_load(f)
    except Exception as e:
        logger.error(f"Cannot read config file: {e}")
        return False
    
    # Run validations
    checks = [
        ("Reference files", lambda: validate_reference_files(config)),
        ("Directories", lambda: validate_directories(inputdir, outdir)),
        ("Sample file", lambda: validate_sample_file(config.get('samplesfile', 'samples.csv'))),
    ]
    
    results = []
    for check_name, check_func in checks:
        try:
            result = check_func()
            results.append(result)
        except Exception as e:
            logger.error(f"Error during {check_name} check: {e}")
            results.append(False)
    
    logger.info("=" * 60)
    if all(results):
        logger.info("✓ ALL VALIDATIONS PASSED")
        logger.info("=" * 60)
        return True
    else:
        logger.error("✗ VALIDATION FAILED")
        logger.error("=" * 60)
        return False


if __name__ == '__main__':
    if len(sys.argv) < 4:
        print("Usage: validate.py <config.yml> <inputdir> <outputdir>")
        sys.exit(1)
    
    configfile = sys.argv[1]
    inputdir = sys.argv[2]
    outdir = sys.argv[3]
    
    valid = validate_pipeline_config(configfile, inputdir, outdir)
    sys.exit(0 if valid else 1)
