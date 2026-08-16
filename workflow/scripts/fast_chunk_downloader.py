#!/usr/bin/env python3
"""High-speed multi-threaded chunk downloader with file concatenation."""

import urllib.request
import concurrent.futures
import os
import sys
import time

def download_part(url, start_byte, end_byte, part_path):
    req = urllib.request.Request(url, headers={
        'User-Agent': 'Mozilla/5.0',
        'Range': f'bytes={start_byte}-{end_byte}'
    })
    with urllib.request.urlopen(req) as resp, open(part_path, 'wb') as f:
        while True:
            chunk = resp.read(1024 * 1024)
            if not chunk:
                break
            f.write(chunk)

def download_file(url, output_path, num_workers=16):
    req = urllib.request.Request(url, method='HEAD', headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req) as resp:
        total_size = int(resp.headers['Content-Length'])
    
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    chunk_size = total_size // num_workers
    part_files = []
    chunks = []
    
    for i in range(num_workers):
        start = i * chunk_size
        end = total_size - 1 if i == num_workers - 1 else (i + 1) * chunk_size - 1
        part_path = f"{output_path}.part_{i}"
        part_files.append(part_path)
        chunks.append((start, end, part_path))
        
    start_time = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(download_part, url, s, e, p) for s, e, p in chunks]
        for f in concurrent.futures.as_completed(futures):
            if f.exception():
                raise f.exception()
                
    # Concatenate parts
    with open(output_path, 'wb') as out_f:
        for p in part_files:
            with open(p, 'rb') as in_f:
                while True:
                    b = in_f.read(8 * 1024 * 1024)
                    if not b:
                        break
                    out_f.write(b)
            os.remove(p)
            
    elapsed = time.time() - start_time
    speed_mb = (total_size / (1024**2)) / max(elapsed, 0.001)
    print(f"✓ Downloaded {total_size/(1024**2):.1f} MB in {elapsed:.1f}s ({speed_mb:.1f} MB/s) -> {output_path}")

if __name__ == '__main__':
    if len(sys.argv) > 2:
        url = sys.argv[1]
        out = sys.argv[2]
        workers = int(sys.argv[3]) if len(sys.argv) > 3 else 16
        download_file(url, out, num_workers=workers)
