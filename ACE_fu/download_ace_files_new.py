#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACE CDAWeb crawler and downloader.

This script follows the same split used by the MMS downloader:
1. crawl CDAWeb directory listings and choose the newest CDF version per day;
2. skip complete local files;
3. download missing/incomplete files with optional multi-thread Range requests.
"""

from __future__ import annotations

import argparse
import datetime as dt
import html
import json
import os
import random
import re
import sys
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple
from urllib.parse import urljoin

import requests
from requests.adapters import HTTPAdapter
from urllib3.util.retry import Retry

try:
    import warnings
    from urllib3.exceptions import NotOpenSSLWarning

    warnings.filterwarnings("ignore", category=NotOpenSSLWarning)
except Exception:
    pass


CDAWEB_ACE_ROOT = "https://cdaweb.gsfc.nasa.gov/pub/data/ace/"


@dataclass(frozen=True)
class Product:
    key: str
    path: str
    note: str = ""

    @property
    def base_url(self) -> str:
        return urljoin(CDAWEB_ACE_ROOT, self.path.strip("/") + "/")


PRODUCTS: Dict[str, Product] = {
    "mfi_h0": Product("mfi_h0", "mag/level_2_cdaweb/mfi_h0", "MAG/MFI H0"),
    "mfi_h1": Product("mfi_h1", "mag/level_2_cdaweb/mfi_h1", "MAG/MFI H1"),
    "mfi_h2": Product("mfi_h2", "mag/level_2_cdaweb/mfi_h2", "MAG/MFI H2"),
    "mfi_h3": Product("mfi_h3", "mag/level_2_cdaweb/mfi_h3", "MAG/MFI H3"),
    "mfi_k0": Product("mfi_k0", "mag/level_2_cdaweb/mfi_k0", "MAG/MFI key parameter"),
    "mfi_k1": Product("mfi_k1", "mag/level_2_cdaweb/mfi_k1", "MAG/MFI key parameter"),
    "mfi_k2": Product("mfi_k2", "mag/level_2_cdaweb/mfi_k2", "MAG/MFI key parameter"),
    "swe_h0": Product("swe_h0", "swepam/level_2_cdaweb/swe_h0", "SWEPAM H0"),
    "swe_h2": Product("swe_h2", "swepam/level_2_cdaweb/swe_h2", "SWEPAM H2"),
    "swe_k0": Product("swe_k0", "swepam/level_2_cdaweb/swe_k0", "SWEPAM key parameter"),
    "swe_k1": Product("swe_k1", "swepam/level_2_cdaweb/swe_k1", "SWEPAM key parameter"),
    "epm_h1": Product("epm_h1", "epam/level_2_cdaweb/epm_h1", "EPAM H1"),
    "epm_h2": Product("epm_h2", "epam/level_2_cdaweb/epm_h2", "EPAM H2"),
    "epm_h3": Product("epm_h3", "epam/level_2_cdaweb/epm_h3", "EPAM H3"),
    "epm_k0": Product("epm_k0", "epam/level_2_cdaweb/epm_k0", "EPAM key parameter"),
    "epm_k1": Product("epm_k1", "epam/level_2_cdaweb/epm_k1", "EPAM key parameter"),
    "sis_h1": Product("sis_h1", "sis/level_2_cdaweb/sis_h1", "SIS H1"),
    "sis_h2": Product("sis_h2", "sis/level_2_cdaweb/sis_h2", "SIS H2"),
    "sis_k0": Product("sis_k0", "sis/level_2_cdaweb/sis_k0", "SIS key parameter"),
}

ALIASES = {
    "mag": "mfi_h0",
    "mfi": "mfi_h0",
    "mfi_h": "mfi_h0",
    "fgm": "mfi_h0",
    "swe": "swe_h0",
    "swepam": "swe_h0",
    "plasma": "swe_h0",
    "epam": "epm_h1",
    "epm": "epm_h1",
    "sis": "sis_h1",
}

_progress_lock = threading.Lock()


def _build_session(trust_env: bool = True) -> requests.Session:
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
    adapter = HTTPAdapter(max_retries=retry, pool_connections=64, pool_maxsize=64)
    s.mount("http://", adapter)
    s.mount("https://", adapter)
    s.headers.update(
        {
            "User-Agent": (
                "Mozilla/5.0 (Macintosh; Intel Mac OS X) "
                "AppleWebKit/537.36 (KHTML, like Gecko) "
                "Chrome/120.0 Safari/537.36"
            ),
            "Accept": "*/*",
        }
    )
    return s


_SESSION_ENV = _build_session(True)
_SESSION_NOENV = _build_session(False)


def _sleep_backoff(attempt: int, cap: float = 8.0) -> None:
    base = min(cap, 0.5 * (2 ** (attempt - 1)))
    time.sleep(base + random.uniform(0, 0.25))


def parse_bool(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y", "on"}


def parse_date(value: str) -> dt.date:
    value = str(value).strip()
    if "T" in value:
        value = value.split("T", 1)[0]
    value = value.replace(".", "-")
    if re.fullmatch(r"\d{8}", value):
        return dt.datetime.strptime(value, "%Y%m%d").date()
    return dt.datetime.strptime(value, "%Y-%m-%d").date()


def parse_date_range(value: str) -> Tuple[dt.date, dt.date]:
    parts = re.split(r"/|,", str(value).strip())
    parts = [p for p in parts if p]
    if not parts:
        raise ValueError("date is empty")
    start = parse_date(parts[0])
    end = parse_date(parts[-1]) if len(parts) > 1 else start
    if end < start:
        start, end = end, start
    return start, end


def iter_years(start: dt.date, end: dt.date) -> Iterable[int]:
    for year in range(start.year, end.year + 1):
        yield year


def normalize_product(value: str) -> str:
    key = str(value).strip().lower().replace("-", "_")
    key = ALIASES.get(key, key)
    if key not in PRODUCTS:
        raise KeyError(
            f"Unknown ACE product '{value}'. Use --show-products to list supported products."
        )
    return key


def parse_product_list(value: str) -> List[str]:
    products: List[str] = []
    for item in re.split(r"[,;\s]+", str(value).strip()):
        if not item:
            continue
        products.append(normalize_product(item))
    if not products:
        products.append("mfi_h0")
    return list(dict.fromkeys(products))


def fetch_html(url: str, timeout: Tuple[float, float] = (8, 45)) -> str:
    last_err: Optional[BaseException] = None
    for sess in (_SESSION_ENV, _SESSION_NOENV):
        for attempt in range(1, 7):
            try:
                with sess.get(url, timeout=timeout, allow_redirects=True) as r:
                    if r.status_code == 404:
                        return ""
                    r.raise_for_status()
                    return r.text
            except (
                requests.exceptions.ConnectionError,
                requests.exceptions.ReadTimeout,
                requests.exceptions.SSLError,
                requests.exceptions.HTTPError,
            ) as exc:
                last_err = exc
                if attempt == 6:
                    break
                _sleep_backoff(attempt)
    raise RuntimeError(f"Failed to read directory: {url}. Last error: {last_err}")


def extract_cdf_links(index_html: str) -> List[str]:
    links = re.findall(r'href=["\']([^"\']+\.cdf)["\']', index_html, flags=re.IGNORECASE)
    names = [os.path.basename(html.unescape(link)) for link in links]
    return sorted(dict.fromkeys(names))


def filename_date_version(filename: str) -> Optional[Tuple[dt.date, int]]:
    m = re.search(r"_(\d{8})_v(\d+)\.cdf$", filename, flags=re.IGNORECASE)
    if not m:
        return None
    return dt.datetime.strptime(m.group(1), "%Y%m%d").date(), int(m.group(2))


def crawl_product(product_key: str, start: dt.date, end: dt.date) -> List[dict]:
    product = PRODUCTS[product_key]
    by_date: Dict[dt.date, dict] = {}

    for year in iter_years(start, end):
        year_url = urljoin(product.base_url, f"{year}/")
        index_html = fetch_html(year_url)
        if not index_html:
            print(f"Warning: remote year directory not found: {year_url}", file=sys.stderr)
            continue

        for filename in extract_cdf_links(index_html):
            parsed = filename_date_version(filename)
            if parsed is None:
                continue
            day, version = parsed
            if day < start or day > end:
                continue
            old = by_date.get(day)
            if old is None or version > old["version"]:
                by_date[day] = {
                    "product": product.key,
                    "date": day.isoformat(),
                    "year": year,
                    "version": version,
                    "filename": filename,
                    "url": urljoin(year_url, filename),
                }

    return [by_date[k] for k in sorted(by_date)]


def crawl_products(product_keys: Sequence[str], start: dt.date, end: dt.date) -> List[dict]:
    files: List[dict] = []
    for product_key in product_keys:
        files.extend(crawl_product(product_key, start, end))
    return sorted(files, key=lambda item: (item["product"], item["date"], item["filename"]))


def _head_with_retries(
    session: requests.Session, url: str, timeout: Tuple[float, float], attempts: int = 6
) -> Tuple[str, dict]:
    last_err: Optional[BaseException] = None
    for i in range(1, attempts + 1):
        try:
            with session.head(url, allow_redirects=True, timeout=timeout) as r:
                r.raise_for_status()
                return r.url, dict(r.headers)
        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.ReadTimeout,
            requests.exceptions.SSLError,
            requests.exceptions.HTTPError,
        ) as exc:
            last_err = exc
            if i == attempts:
                break
            _sleep_backoff(i)
    raise RuntimeError(f"HEAD failed for {url}. Last error: {last_err}")


def _probe_size_by_range(
    session: requests.Session, url: str, timeout: Tuple[float, float], attempts: int = 6
) -> Tuple[str, int]:
    last_err: Optional[BaseException] = None
    headers = {"Range": "bytes=0-0", "Accept-Encoding": "identity"}
    for i in range(1, attempts + 1):
        try:
            with session.get(
                url, headers=headers, stream=True, allow_redirects=True, timeout=timeout
            ) as r:
                if r.status_code != 206:
                    raise RuntimeError(f"Range probe expected 206, got {r.status_code}")
                cr = r.headers.get("Content-Range", "")
                if "/" not in cr:
                    raise RuntimeError(f"Missing/invalid Content-Range: {cr}")
                total = int(cr.split("/")[-1])
                return r.url, total
        except (
            requests.exceptions.ConnectionError,
            requests.exceptions.ReadTimeout,
            requests.exceptions.SSLError,
            RuntimeError,
            ValueError,
        ) as exc:
            last_err = exc
            if i == attempts:
                break
            _sleep_backoff(i)
    raise RuntimeError(f"Range probe failed for {url}. Last error: {last_err}")


def resolve_url_and_size(url: str, timeout: Tuple[float, float] = (8, 45)) -> Tuple[str, int, requests.Session, bool]:
    for sess in (_SESSION_ENV, _SESSION_NOENV):
        try:
            final_url, headers = _head_with_retries(sess, url, timeout=timeout, attempts=6)
            cl = headers.get("Content-Length") or headers.get("content-length")
            if cl is not None:
                content_size = int(cl)
                if content_size > 0:
                    accept_ranges = (
                        headers.get("Accept-Ranges")
                        or headers.get("accept-ranges")
                        or ""
                    ).lower()
                    return final_url, content_size, sess, accept_ranges == "bytes"
        except Exception:
            continue

    for sess in (_SESSION_ENV, _SESSION_NOENV):
        try:
            final_url, total = _probe_size_by_range(sess, url, timeout=timeout, attempts=6)
            if total > 0:
                return final_url, total, sess, True
        except Exception:
            continue

    raise RuntimeError(f"Failed to resolve Content-Length for {url}")


def download_chunk(args) -> None:
    session, url, start_byte, end_byte, local_file_path, progress_callback = args
    timeout = (8, 45)
    current_pos = start_byte
    max_retries = 10
    attempt = 0
    base_headers = {"Accept-Encoding": "identity"}

    while current_pos <= end_byte:
        headers = dict(base_headers)
        headers["Range"] = f"bytes={current_pos}-{end_byte}"
        try:
            with session.get(
                url, headers=headers, stream=True, timeout=timeout, allow_redirects=True
            ) as response:
                if response.status_code != 206:
                    raise RuntimeError(
                        f"Expected 206 for range request, got {response.status_code}"
                    )
                content_range = response.headers.get("Content-Range", "")
                if not content_range.startswith("bytes"):
                    raise RuntimeError(f"Missing/invalid Content-Range: {content_range}")
                try:
                    cr_range = content_range.split(" ")[1].split("/")[0]
                    cr_start = int(cr_range.split("-")[0])
                except Exception as exc:
                    raise RuntimeError(f"Unparseable Content-Range: {content_range}") from exc
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
        except (
            requests.exceptions.ChunkedEncodingError,
            requests.exceptions.ConnectionError,
            requests.exceptions.ReadTimeout,
            requests.exceptions.SSLError,
            RuntimeError,
        ) as exc:
            attempt += 1
            if attempt > max_retries:
                raise RuntimeError(
                    f"Failed downloading range {start_byte}-{end_byte} after "
                    f"{max_retries} retries. Last error: {exc}"
                ) from exc
            _sleep_backoff(attempt, cap=6.0)


def stream_download(
    session: requests.Session, final_url: str, local_file_path: Path, expected_size: Optional[int]
) -> None:
    headers = {"Accept-Encoding": "identity"}
    start_time = time.time()
    current_size = 0
    with session.get(
        final_url, headers=headers, stream=True, allow_redirects=True, timeout=(8, 45)
    ) as r:
        r.raise_for_status()
        with open(local_file_path, "wb") as f:
            for chunk in r.iter_content(chunk_size=1024 * 256):
                if not chunk:
                    continue
                f.write(chunk)
                current_size += len(chunk)
                elapsed = max(time.time() - start_time, 0.001)
                speed = current_size / elapsed / 1024 / 1024
                if expected_size:
                    pct = current_size / expected_size * 100
                    msg = (
                        f"\rDownloaded: {current_size / 1024 / 1024:.2f} MB "
                        f"[{pct:6.2f}%] {speed:.2f} MB/s"
                    )
                else:
                    msg = (
                        f"\rDownloaded: {current_size / 1024 / 1024:.2f} MB "
                        f"{speed:.2f} MB/s"
                    )
                print(msg, end="", flush=True)
    print()


def range_download(
    session: requests.Session,
    final_url: str,
    local_file_path: Path,
    content_size: int,
    num_threads: int,
) -> None:
    num_threads = max(1, int(num_threads))
    chunk_size = content_size // num_threads
    if chunk_size <= 0:
        num_threads = 1
        chunk_size = content_size

    local_file_path.parent.mkdir(parents=True, exist_ok=True)
    with open(local_file_path, "wb") as file:
        file.truncate(content_size)

    start_time = time.time()
    current_size = 0

    def progress_callback(nbytes: int) -> None:
        nonlocal current_size
        with _progress_lock:
            current_size += nbytes
            elapsed_time = max(time.time() - start_time, 0.001)
            pct = (current_size / content_size) * 100
            speed = current_size / elapsed_time / 1024 / 1024
            print(
                f"\rDownloaded: {current_size / 1024 / 1024:.2f} MB "
                f"[{pct:6.2f}%] {speed:.2f} MB/s",
                end="",
                flush=True,
            )

    futures = []
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        for i in range(num_threads):
            start_byte = i * chunk_size
            end_byte = start_byte + chunk_size - 1 if i < num_threads - 1 else content_size - 1
            args = (session, final_url, start_byte, end_byte, str(local_file_path), progress_callback)
            futures.append(executor.submit(download_chunk, args))
        for future in futures:
            future.result()
    print()


def target_path(out_dir: Path, item: dict, keep_tree: bool) -> Path:
    if keep_tree:
        return out_dir / item["product"] / str(item["year"]) / item["filename"]
    return out_dir / item["filename"]


def is_complete(path: Path, size: Optional[int], check_size: bool) -> bool:
    if not path.is_file():
        return False
    if not check_size:
        return True
    if size is None or size <= 0:
        return False
    return path.stat().st_size == size


def download_one(item: dict, out_dir: Path, threads: int, check_size: bool, keep_tree: bool) -> str:
    local_path = target_path(out_dir, item, keep_tree)
    local_path.parent.mkdir(parents=True, exist_ok=True)

    final_url, content_size, session, range_hint = resolve_url_and_size(item["url"])
    item["remote_bytes"] = content_size
    item["local_path"] = str(local_path)

    if is_complete(local_path, content_size, check_size):
        return "skip"

    if local_path.exists() and check_size:
        print(f"Incomplete old file, re-downloading: {local_path.name}")

    range_ok = False
    if threads > 1:
        try:
            _final2, total2 = _probe_size_by_range(session, final_url, timeout=(8, 45), attempts=3)
            if total2 > 0:
                final_url = _final2
                content_size = total2
                range_ok = True
        except Exception:
            range_ok = False

    if threads > 1 and range_ok and range_hint:
        range_download(session, final_url, local_path, content_size, threads)
    elif threads > 1 and range_ok:
        range_download(session, final_url, local_path, content_size, threads)
    else:
        stream_download(session, final_url, local_path, content_size)

    if check_size and local_path.stat().st_size != content_size:
        raise RuntimeError(
            f"Size check failed for {local_path.name}: "
            f"{local_path.stat().st_size} != {content_size}"
        )
    return "download"


def show_products() -> None:
    for key in sorted(PRODUCTS):
        product = PRODUCTS[key]
        print(f"{key:8s}  {product.path:42s}  {product.note}")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Download ACE CDF data from NASA CDAWeb public directories."
    )
    parser.add_argument("--date", required=False, help="YYYY-MM-DD/YYYY-MM-DD or YYYYMMDD/YYYYMMDD")
    parser.add_argument(
        "--product",
        default="mfi_h0",
        help="Comma-separated products, e.g. mfi_h0,swe_h0,epm_h1. Default: mfi_h0",
    )
    parser.add_argument("--out", default=".", help="Output directory")
    parser.add_argument("--threads", type=int, default=8, help="Download threads per file")
    parser.add_argument("--check-size", default="1", help="1/0, verify local size with Content-Length")
    parser.add_argument("--keep-tree", default="0", help="1/0, save as product/year/filename")
    parser.add_argument("--list-only", action="store_true", help="Only crawl and print matched files")
    parser.add_argument("--json", action="store_true", help="Print list-only result as JSON")
    parser.add_argument("--show-products", action="store_true", help="Show supported products and exit")
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.show_products:
        show_products()
        return 0

    if not args.date:
        parser.error("--date is required unless --show-products is used")

    start, end = parse_date_range(args.date)
    product_keys = parse_product_list(args.product)
    out_dir = Path(args.out).expanduser()
    keep_tree = parse_bool(args.keep_tree)
    check_size = parse_bool(args.check_size)

    files = crawl_products(product_keys, start, end)

    if args.list_only:
        if args.json:
            sys.stdout.write(json.dumps(files, ensure_ascii=False, indent=2) + "\n")
        else:
            for item in files:
                print(f"{item['product']} {item['date']} {item['filename']} {item['url']}")
            print(f"Matched {len(files)} file(s).")
        return 0

    if not files:
        print("No ACE files matched this request.")
        return 0

    out_dir.mkdir(parents=True, exist_ok=True)
    print("========== ACE download started ==========")
    print(f"Date range: {start.isoformat()} / {end.isoformat()}")
    print(f"Products: {', '.join(product_keys)}")
    print(f"Output: {out_dir}")
    print(f"Files: {len(files)}")

    downloaded = 0
    skipped = 0
    for idx, item in enumerate(files, start=1):
        print(f"[{idx}/{len(files)}] {item['filename']}")
        t0 = time.time()
        status = download_one(item, out_dir, args.threads, check_size, keep_tree)
        dt_s = max(time.time() - t0, 0.001)
        local = Path(item.get("local_path", target_path(out_dir, item, keep_tree)))
        mb = local.stat().st_size / 1024 / 1024 if local.exists() else 0.0
        if status == "skip":
            skipped += 1
            print(f"Already exists: {local}")
        else:
            downloaded += 1
            print(f"Saved: {local} ({mb:.2f} MB, {mb / dt_s:.2f} MB/s)")

    print("========== ACE download finished ==========")
    print(f"Downloaded: {downloaded}; skipped: {skipped}; total: {len(files)}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        raise SystemExit(130)
