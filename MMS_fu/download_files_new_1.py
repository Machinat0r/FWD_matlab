#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan  7 19:34:42 2026

------- written by Wending Fu in Beijing ------
"""
import sys
import time
import threading
import requests
from concurrent.futures import ThreadPoolExecutor

_progress_lock = threading.Lock()


def download_chunk(args):
    """
    下载指定 Range，并在连接中断时自动续传补齐。
    """
    url, start_byte, end_byte, local_file_path, file_name, progress_callback = args

    timeout = (5, 30)

    current_pos = start_byte
    max_retries = 8
    attempt = 0

    base_headers = {"Accept-Encoding": "identity"}

    while current_pos <= end_byte:
        headers = dict(base_headers)
        headers["Range"] = f"bytes={current_pos}-{end_byte}"

        try:
            with requests.get(url, headers=headers, stream=True, timeout=timeout, allow_redirects=True) as response:

                if response.status_code != 206:
                    raise RuntimeError(f"Expected 206 for range request, got {response.status_code}")

                
                content_range = response.headers.get("Content-Range", "")
                
                if not content_range.startswith("bytes"):
                    raise RuntimeError(f"Missing/invalid Content-Range: {content_range}")

                
                try:
                    cr_range = content_range.split(" ")[1].split("/")[0]  # 0-1023
                    cr_start = int(cr_range.split("-")[0])
                except Exception:
                    raise RuntimeError(f"Unparseable Content-Range: {content_range}")

                if cr_start != current_pos:
                    raise RuntimeError(f"Content-Range start mismatch: {cr_start} != {current_pos}")

                with open(local_file_path, "r+b") as file:
                    file.seek(current_pos)

                    for chunk in response.iter_content(chunk_size=1024 * 256):
                        if not chunk:
                            continue
                        file.write(chunk)
                        current_pos += len(chunk)
                        progress_callback(len(chunk))

        except (requests.exceptions.ChunkedEncodingError,
                requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout,
                RuntimeError) as e:
            attempt += 1
            if attempt > max_retries:
                raise RuntimeError(
                    f"Failed downloading range {start_byte}-{end_byte} after {max_retries} retries. "
                    f"Last error: {e}"
                ) from e

            
            backoff = min(2.0, 0.25 * (2 ** (attempt - 1)))
            time.sleep(backoff)
            
            continue


def write_files(file_url, local_file_path, file_name, num_threads):
    
    response = requests.head(file_url, allow_redirects=True, timeout=(5, 30))
    response.raise_for_status()

    if "content-length" not in response.headers:
        raise RuntimeError("HEAD response missing Content-Length; cannot safely multi-thread range download.")

    content_size = int(response.headers["content-length"])
    if content_size <= 0:
        raise RuntimeError(f"Invalid Content-Length: {content_size}")

    chunk_size = content_size // num_threads
    if chunk_size <= 0:
        
        num_threads = 1
        chunk_size = content_size

    
    with open(local_file_path, "wb") as file:
        file.truncate(content_size)

    start_time = time.time()
    current_size = 0

    def progress_callback(nbytes):
        nonlocal current_size
        with _progress_lock:
            current_size += nbytes
            elapsed_time = time.time() - start_time
            pct = (current_size / content_size) * 100
            print(
                f"\r Downloaded: [{current_size / 1024 / 1024:.2f} MB]"
                f"[{pct:.2f}%]"
                f"Elapsed Time: {elapsed_time:.2f} seconds",
                end="",
                flush=True,
            )

    futures = []

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for i in range(num_threads):
            start_byte = i * chunk_size
            end_byte = start_byte + chunk_size - 1 if i < num_threads - 1 else content_size - 1
            args = (file_url, start_byte, end_byte, local_file_path, file_name, progress_callback)
            futures.append(executor.submit(download_chunk, args))

        
        for future in futures:
            future.result()

    print()


if __name__ == "__main__":
    file_url = sys.argv[1]
    local_file_path = sys.argv[2]
    file_name = sys.argv[3]
    num_threads = int(sys.argv[4])

    write_files(file_url, local_file_path, file_name, num_threads)
