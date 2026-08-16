#!/usr/bin/env python3
"""High-speed multi-connection HTTP range downloader."""

import urllib.request
import concurrent.futures
import os
import sys
import time

def download_chunk(url, start_byte, end_byte, file_path):
    req = urllib.request.Request(url, headers={'Range': f'bytes={start_byte}-{end_byte}'})
    with urllib.request.urlopen(req) as resp:
        data = resp.read()
    with open(file_path, 'r+b') as f:
        f.seek(start_byte)
        f.write(data)

def parallel_download(url, output_path, num_workers=32):
    req = urllib.request.Request(url, method='HEAD')
    with urllib.request.urlopen(req) as resp:
        total_size = int(resp.headers['Content-Length'])
    
    os.makedirs(os.path.dirname(os.path.abspath(output_path)), exist_ok=True)
    with open(output_path, 'wb') as f:
        f.seek(total_size - 1)
        f.write(b'\0')
        
    chunk_size = total_size // num_workers
    chunks = []
    for i in range(num_workers):
        start = i * chunk_size
        end = total_size - 1 if i == num_workers - 1 else (i + 1) * chunk_size - 1
        chunks.append((start, end))
        
    start_time = time.time()
    with concurrent.futures.ThreadPoolExecutor(max_workers=num_workers) as executor:
        futures = [executor.submit(download_chunk, url, s, e, output_path) for s, e in chunks]
        for f in concurrent.futures.as_completed(futures):
            if f.exception():
                raise f.exception()
        
    elapsed = time.time() - start_time
    speed_mb = (total_size / (1024**2)) / max(elapsed, 0.001)
    print(f"✓ Downloaded {total_size/(1024**3):.2f} GB in {elapsed:.1f}s ({speed_mb:.1f} MB/s) -> {output_path}")

if __name__ == '__main__':
    if len(sys.argv) > 2:
        url = sys.argv[1]
        out = sys.argv[2]
        workers = int(sys.argv[3]) if len(sys.argv) > 3 else 32
        parallel_download(url, out, num_workers=workers)
