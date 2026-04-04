#!/usr/bin/env python

'''Post-processing utilities to aggregate pipeline benchmarks and logs.'''

import os
import glob
import json
import pandas as pd
import logging
from pathlib import Path
from datetime import datetime

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')
logger = logging.getLogger(__name__)


def aggregate_benchmarks(benchmark_dir, output_file='benchmark_summary.csv'):
    """
    Aggregate benchmark files from all rules into a single CSV.
    
    Args:
        benchmark_dir: Path to benchmarks directory
        output_file: Name of output summary file
    
    Returns:
        pandas.DataFrame: Aggregated benchmark data
    """
    logger.info("Aggregating benchmark files...")
    
    if not os.path.isdir(benchmark_dir):
        logger.warning(f"Benchmark directory not found: {benchmark_dir}")
        return None
    
    benchmark_files = glob.glob(os.path.join(benchmark_dir, '**/*.txt'), recursive=True)
    
    if not benchmark_files:
        logger.warning("No benchmark files found")
        return None
    
    data = []
    
    for bf in benchmark_files:
        try:
            with open(bf, 'r') as f:
                lines = f.readlines()
                
            # Parse Snakemake benchmark format (tab-separated)
            if len(lines) >= 2:
                header = lines[0].strip().split('\t')
                values = lines[1].strip().split('\t')
                
                row = dict(zip(header, values))
                row['rule_file'] = os.path.basename(bf)
                row['relative_path'] = os.path.relpath(bf, benchmark_dir)
                
                data.append(row)
        except Exception as e:
            logger.warning(f"Could not parse benchmark file {bf}: {e}")
    
    if not data:
        logger.warning("No valid benchmark data found")
        return None
    
    df = pd.DataFrame(data)
    
    # Convert numeric columns
    numeric_cols = ['s', 'h:m:s', 'max_rss', 'start_time', 'end_time']
    for col in numeric_cols:
        if col in df.columns:
            try:
                df[col] = pd.to_numeric(df[col], errors='coerce')
            except:
                pass
    
    # Save aggregate
    output_path = os.path.join(benchmark_dir, '..', output_file)
    df.to_csv(output_path, index=False)
    logger.info(f"✓ Aggregated {len(df)} benchmark entries to {output_path}")
    
    return df


def aggregate_logs(log_dir, output_file='log_summary.txt'):
    """
    Aggregate log files with error/warning detection.
    
    Args:
        log_dir: Path to logs directory
        output_file: Name of output summary file
    
    Returns:
        dict: Summary statistics
    """
    logger.info("Analyzing log files...")
    
    if not os.path.isdir(log_dir):
        logger.warning(f"Log directory not found: {log_dir}")
        return None
    
    log_files = glob.glob(os.path.join(log_dir, '**/*.log'), recursive=True)
    
    if not log_files:
        logger.warning("No log files found")
        return None
    
    summary = {
        'total_logs': len(log_files),
        'logs_with_errors': 0,
        'logs_with_warnings': 0,
        'errors': [],
        'warnings': [],
        'timestamp': datetime.now().isoformat()
    }
    
    for log_file in log_files:
        try:
            with open(log_file, 'r', encoding='utf-8', errors='ignore') as f:
                content = f.read().lower()
                
            rule_name = os.path.relpath(log_file, log_dir)
            
            # Search for error patterns
            if any(x in content for x in ['error', 'exception', 'failed', 'fatal']):
                summary['logs_with_errors'] += 1
                summary['errors'].append(rule_name)
            
            # Search for warning patterns
            if 'warning' in content:
                summary['logs_with_warnings'] += 1
                summary['warnings'].append(rule_name)
        
        except Exception as e:
            logger.warning(f"Could not read log file {log_file}: {e}")
    
    # Write summary report
    output_path = os.path.join(log_dir, '..', output_file)
    with open(output_path, 'w') as f:
        f.write("=" * 70 + "\n")
        f.write("PIPELINE LOG SUMMARY\n")
        f.write("=" * 70 + "\n\n")
        
        f.write(f"Generated: {summary['timestamp']}\n")
        f.write(f"Total log files: {summary['total_logs']}\n")
        f.write(f"Logs with errors: {summary['logs_with_errors']}\n")
        f.write(f"Logs with warnings: {summary['logs_with_warnings']}\n\n")
        
        if summary['errors']:
            f.write("ERRORS DETECTED:\n")
            f.write("-" * 70 + "\n")
            for error in summary['errors'][:10]:
                f.write(f"  • {error}\n")
            if len(summary['errors']) > 10:
                f.write(f"  ... and {len(summary['errors']) - 10} more\n")
            f.write("\n")
        
        if summary['warnings']:
            f.write("WARNINGS DETECTED:\n")
            f.write("-" * 70 + "\n")
            for warning in summary['warnings'][:10]:
                f.write(f"  • {warning}\n")
            if len(summary['warnings']) > 10:
                f.write(f"  ... and {len(summary['warnings']) - 10} more\n")
    
    logger.info(f"✓ Log analysis saved to {output_path}")
    logger.info(f"  - {summary['logs_with_errors']} logs with errors")
    logger.info(f"  - {summary['logs_with_warnings']} logs with warnings")
    
    return summary


def create_execution_report(outdir):
    """
    Create a comprehensive execution report combining benchmarks, logs, and stats.
    
    Args:
        outdir: Output directory containing benchmarks and logs subdirectories
    
    Returns:
        str: Path to generated report
    """
    logger.info("Creating execution report...")
    
    benchmark_dir = os.path.join(outdir, 'benchmarks')
    log_dir = os.path.join(outdir, 'logs')
    
    # Aggregate data
    benchmarks_df = aggregate_benchmarks(benchmark_dir)
    log_summary = aggregate_logs(log_dir)
    
    # Create HTML report
    report_path = os.path.join(outdir, 'results', 'execution_report.html')
    os.makedirs(os.path.dirname(report_path), exist_ok=True)
    
    html_content = """
    <html>
    <head>
        <title>Pipeline Execution Report</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }
            .header { background: #2c3e50; color: white; padding: 20px; border-radius: 5px; }
            .section { background: white; margin: 20px 0; padding: 20px; border-radius: 5px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
            table { width: 100%; border-collapse: collapse; }
            th, td { padding: 10px; text-align: left; border-bottom: 1px solid #ddd; }
            th { background: #34495e; color: white; }
            tr:hover { background: #f5f5f5; }
            .error { color: #e74c3c; }
            .warning { color: #f39c12; }
            .success { color: #27ae60; }
            .metric { display: inline-block; background: #ecf0f1; padding: 15px; margin: 10px; border-radius: 5px; }
        </style>
    </head>
    <body>
        <div class="header">
            <h1>Pipeline Execution Report</h1>
            <p><strong>Generated:</strong> """ + datetime.now().strftime("%Y-%m-%d %H:%M:%S") + """</p>
        </div>
    """
    
    # Add log summary
    if log_summary:
        html_content += """
        <div class="section">
            <h2>Log Summary</h2>
            <div class="metric"><strong>Total Logs:</strong> """ + str(log_summary['total_logs']) + """</div>
            <div class="metric error"><strong>Errors:</strong> """ + str(log_summary['logs_with_errors']) + """</div>
            <div class="metric warning"><strong>Warnings:</strong> """ + str(log_summary['logs_with_warnings']) + """</div>
        </div>
        """
    
    # Add benchmarks table
    if benchmarks_df is not None and not benchmarks_df.empty:
        html_content += """
        <div class="section">
            <h2>Rule Benchmarks</h2>
            <table>
                <tr>
                    <th>Rule</th>
                    <th>Duration (s)</th>
                    <th>Memory (MB)</th>
                    <th>CPU %</th>
                </tr>
        """
        for _, row in benchmarks_df.head(20).iterrows():
            html_content += f"""
                <tr>
                    <td>{row.get('rule_file', 'unknown')}</td>
                    <td>{row.get('s', 'N/A')}</td>
                    <td>{row.get('max_rss', 'N/A')}</td>
                    <td>N/A</td>
                </tr>
            """
        html_content += "</table></div>"
    
    html_content += """
    </body>
    </html>
    """
    
    with open(report_path, 'w') as f:
        f.write(html_content)
    
    logger.info(f"✓ Execution report generated: {report_path}")
    return report_path


if __name__ == '__main__':
    import sys
    
    if len(sys.argv) < 2:
        print("Usage: aggregate.py <outdir>")
        sys.exit(1)
    
    outdir = sys.argv[1]
    create_execution_report(outdir)
