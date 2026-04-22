#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 11 10:43:31 2026

------- written by Wending Fu in Beijing ------
"""

import sys
import time
import random
import threading
import requests
from concurrent.futures import ThreadPoolExecutor
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

_progress_lock = threading.Lock()

try:
    import warnings
    from urllib3.exceptions import NotOpenSSLWarning
    warnings.filterwarnings("ignore", category=NotOpenSSLWarning)
except Exception:
    pass


def _build_session(trust_env: bool) -> requests.Session:
    s = requests.Session()
    s.trust_env = trust_env

    retry = Retry(
        total=6,
        connect=6,
        read=6,
        status=6,
        backoff_factor=0.5,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["HEAD", "GET"]),
        raise_on_status=False,
        respect_retry_after_header=True,
    )

    adapter = HTTPAdapter(
        max_retries=retry,
        pool_connections=64,
        pool_maxsize=64,
    )

    s.mount("http://", adapter)
    s.mount("https://", adapter)

    s.headers.update({
        "User-Agent": "Mozilla/5.0 (Macintosh; Intel Mac OS X) "
                      "AppleWebKit/537.36 (KHTML, like Gecko) "
                      "Chrome/120.0 Safari/537.36",
        "Accept": "*/*",
    })
    return s


_SESSION_ENV = _build_session(trust_env=True)
_SESSION_NOENV = _build_session(trust_env=False)


def _sleep_backoff(attempt: int, cap: float = 8.0):
    base = min(cap, 0.5 * (2 ** (attempt - 1)))
    jitter = random.uniform(0, 0.25)
    time.sleep(base + jitter)


def _head_with_retries(session: requests.Session, url: str, timeout, attempts: int = 6):
    last_err = None
    for i in range(1, attempts + 1):
        try:
            with session.head(url, allow_redirects=True, timeout=timeout) as r:
                r.raise_for_status()
                return r.url, dict(r.headers)
        except (requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout,
                requests.exceptions.SSLError,
                requests.exceptions.HTTPError) as e:
            last_err = e
            if i == attempts:
                break
            _sleep_backoff(i)
    raise last_err


def _probe_size_by_range(session: requests.Session, url: str, timeout, attempts: int = 6):
    last_err = None
    headers = {"Range": "bytes=0-0", "Accept-Encoding": "identity"}
    for i in range(1, attempts + 1):
        try:
            with session.get(url, headers=headers, stream=True, allow_redirects=True, timeout=timeout) as r:
                if r.status_code != 206:
                    raise RuntimeError(f"Probe range expected 206, got {r.status_code}")
                cr = r.headers.get("Content-Range", "")
                if "/" not in cr:
                    raise RuntimeError(f"Missing/invalid Content-Range: {cr}")
                total = int(cr.split("/")[-1])
                return r.url, total
        except (requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout,
                requests.exceptions.SSLError,
                RuntimeError,
                ValueError) as e:
            last_err = e
            if i == attempts:
                break
            _sleep_backoff(i)
    raise last_err


def _resolve_url_and_size(file_url: str, timeout=(5, 30)):
    sessions = (_SESSION_ENV, _SESSION_NOENV)

    for sess in sessions:
        try:
            final_url, headers = _head_with_retries(sess, file_url, timeout=timeout, attempts=6)
            cl = headers.get("Content-Length") or headers.get("content-length")
            if cl is not None:
                content_size = int(cl)
                if content_size > 0:
                    accept_ranges = (headers.get("Accept-Ranges") or headers.get("accept-ranges") or "").lower()
                    return final_url, content_size, sess, (accept_ranges == "bytes")
        except Exception:
            continue

    for sess in sessions:
        try:
            final_url, total = _probe_size_by_range(sess, file_url, timeout=timeout, attempts=6)
            if total > 0:
                return final_url, total, sess, True
        except Exception:
            continue

    raise RuntimeError("Failed to resolve Content-Length via HEAD or Range probe (proxy/no-proxy both failed).")


def download_chunk(args):
    session, url, start_byte, end_byte, local_file_path, file_name, progress_callback = args
    timeout = (5, 30)

    current_pos = start_byte
    max_retries = 10
    attempt = 0

    base_headers = {
        "Accept-Encoding": "identity",
    }

    while current_pos <= end_byte:
        headers = dict(base_headers)
        headers["Range"] = f"bytes={current_pos}-{end_byte}"

        try:
            with session.get(url, headers=headers, stream=True, timeout=timeout, allow_redirects=True) as response:
                if response.status_code != 206:
                    raise RuntimeError(f"Expected 206 for range request, got {response.status_code}")

                content_range = response.headers.get("Content-Range", "")
                if not content_range.startswith("bytes"):
                    raise RuntimeError(f"Missing/invalid Content-Range: {content_range}")

                try:
                    cr_range = content_range.split(" ")[1].split("/")[0]  # "X-Y"
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
                requests.exceptions.SSLError,
                RuntimeError) as e:
            attempt += 1
            if attempt > max_retries:
                raise RuntimeError(
                    f"Failed downloading range {start_byte}-{end_byte} after {max_retries} retries. "
                    f"Last error: {e}"
                ) from e

            _sleep_backoff(attempt, cap=6.0)
            continue


def write_files(file_url, local_file_path, file_name, num_threads):
    final_url, content_size, session, range_hint = _resolve_url_and_size(file_url, timeout=(5, 30))

    range_ok = False
    try:
        _final2, _total2 = _probe_size_by_range(session, final_url, timeout=(5, 30), attempts=3)
        range_ok = True
        if _total2 and _total2 > 0 and _total2 != content_size:
            content_size = _total2
            final_url = _final2
    except Exception:
        range_ok = False

    if not range_ok or content_size <= 0:
        timeout = (5, 30)
        headers = {"Accept-Encoding": "identity"}
        with session.get(final_url, headers=headers, stream=True, allow_redirects=True, timeout=timeout) as r:
            r.raise_for_status()
            with open(local_file_path, "wb") as f:
                start_time = time.time()
                current_size = 0

                def progress_callback(nbytes):
                    nonlocal current_size
                    with _progress_lock:
                        current_size += nbytes
                        elapsed = time.time() - start_time
                        print(
                            f"\r Downloaded: [{current_size / 1024 / 1024:.2f} MB]"
                            f"Elapsed Time: {elapsed:.2f} seconds",
                            end="",
                            flush=True,
                        )

                for chunk in r.iter_content(chunk_size=1024 * 256):
                    if chunk:
                        f.write(chunk)
                        progress_callback(len(chunk))
        print()
        return

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
            args = (session, final_url, start_byte, end_byte, local_file_path, file_name, progress_callback)
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
