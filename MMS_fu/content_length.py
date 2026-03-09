#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 14:54:31 2024

------- written by Wending Fu in Beijing ------
"""

import sys
import warnings

try:
    from urllib3.exceptions import NotOpenSSLWarning
    warnings.filterwarnings("ignore", category=NotOpenSSLWarning)
except Exception:
    pass

import time
import random
import math
import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

def _build_session() -> requests.Session:
    s = requests.Session()

    s.trust_env = True

    retry = Retry(
        total=3,
        connect=3,
        read=3,
        status=3,
        backoff_factor=0.3,
        status_forcelist=(429, 500, 502, 503, 504),
        allowed_methods=frozenset(["HEAD", "GET"]),
        raise_on_status=False,
        respect_retry_after_header=True,
    )
    adapter = HTTPAdapter(max_retries=retry, pool_connections=16, pool_maxsize=16)
    s.mount("http://", adapter)
    s.mount("https://", adapter)

    s.headers.update({
        "User-Agent": "Mozilla/5.0",
        "Accept": "*/*",
        "Accept-Encoding": "identity",
    })
    return s


def _sleep_backoff(attempt: int, cap: float = 10.0):
    base = min(cap, 0.6 * (2 ** (attempt - 1)))
    time.sleep(base + random.uniform(0, 0.4))


def _parse_length(val):
    if val is None:
        return None
    s = str(val).strip()
    if not s:
        return None
    try:
        f = float(s)
        if math.isnan(f) or f <= 0:
            return None
        if not f.is_integer():
            return None
        return int(f)
    except Exception:
        return None


def get_file_size(url: str, timeout=(8, 30), attempts: int = 8) -> int:
    session = _build_session()
    last_err = None

    for i in range(1, attempts + 1):
        try:
            with session.head(url, allow_redirects=True, timeout=timeout) as r:
                r.raise_for_status()
                length = _parse_length(r.headers.get("Content-Length"))
                if length:
                    return length
                final_url = r.url

            headers = {"Range": "bytes=0-0", "Accept-Encoding": "identity"}
            with session.get(final_url, headers=headers, stream=True, allow_redirects=True, timeout=timeout) as r2:
                r2.raise_for_status()
                cr = r2.headers.get("Content-Range", "")
                if "/" in cr:
                    total = _parse_length(cr.split("/")[-1])
                    if total:
                        return total

            raise RuntimeError("No valid Content-Length and no valid Content-Range total.")

        except (requests.exceptions.ConnectTimeout,
                requests.exceptions.ReadTimeout,
                requests.exceptions.ConnectionError,
                requests.exceptions.SSLError,
                requests.exceptions.HTTPError,
                RuntimeError) as e:
            last_err = e
            if i == attempts:
                break
            _sleep_backoff(i)

    raise RuntimeError(f"Failed to get file size after {attempts} attempts. Last error: {last_err}")


if __name__ == "__main__":
    file_url = sys.argv[1]
    try:
        size = get_file_size(file_url)
        sys.stdout.write(f"{size}\n")
        sys.exit(0)
    except Exception as e:        
        print(str(e), file=sys.stderr)
        sys.exit(1)
